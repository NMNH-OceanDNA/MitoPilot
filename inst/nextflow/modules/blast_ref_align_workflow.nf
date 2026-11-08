include {blast_ref_align} from './blast_ref_align.nf'

// SQL fragment: assemble_notes with any segment starting at `tag` stripped (and a
// preceding '; '). Mirrors the helpers in the blast_* workflows; duplicated here so
// this module is self-contained.
def stripTagSql(String tag) {
    def lit = tag.replace("'", "''")
    return "RTRIM(" +
        "CASE WHEN INSTR(COALESCE(assemble_notes,''), '${lit}') > 0 " +
            "THEN SUBSTR(COALESCE(assemble_notes,''), 1, INSTR(COALESCE(assemble_notes,''), '${lit}') - 1) " +
            "ELSE COALESCE(assemble_notes,'') END" +
    ", '; ')"
}

def appendTaggedNoteSql(String tag, String msg) {
    def stripped = stripTagSql(tag)
    def tagged = (tag + ' ' + msg).replace("'", "''")
    return "CASE WHEN ${stripped} = '' THEN '${tagged}' ELSE ${stripped} || '; ${tagged}' END"
}

// Reference sequence, topology and gene table per candidate reference for
// newly-curated samples (gated on curate_out). One row per (ID, candidate
// accession). rank=1 is the top hit. The rotation itself is computed in the
// align process (MitoPilot::ref_rotation_offset) from the reference gene table
// and the unit's own annotations: a circular reference is rotated so the
// sample's first annotated gene anchors both sequences; against a linear GenBank
// reference (native coordinates, rotation 0) a circular sample's alignment copy
// is rotated to the reference's first gene instead (sample_rotation; the stored
// sequence is untouched).
// Gene tables travel as compact 'gene|pos1|pos2|direction;...' strings.
// Assembly sequence comes from the curate output FASTA (already rotated to start_gene).
def featsSql(String tbl, String where) {
    return "COALESCE((SELECT GROUP_CONCAT(x.gene || '|' || x.pos1 || '|' || x.pos2 || '|' || x.direction, ';') " +
           "FROM ${tbl} x WHERE ${where}), '')"
}
params.sqlReadRef =
    'SELECT c.ID, c.accession, c.rank, s.sequence, COALESCE(s.topology, \'linear\'), ' +
        featsSql('blast_ref_annotations', 'x.accession = c.accession') + ' ' +
    'FROM (SELECT ID, accession, MIN(rank) AS rank FROM blast_ref_candidates ' +
         'GROUP BY ID, accession) c ' +
    'JOIN blast_ref_sequences s ON c.accession = s.accession'

// Backfill: samples already validated in a prior run (assemblies.sequence is the
// rotated post-curate sequence) that are missing an alignment row. annotate_switch = 2
// is the post-VALIDATE state and is the only safe signal that the stored sequence
// has been rotated to start_gene - for userAsmb projects, assemblies.sequence is
// pre-populated with the unrotated user FASTA in WF1_userAsmb and only becomes
// rotated after CURATE runs in WF2. The unit's annotations are in the DB here, so
// the sample gene table comes from SQL too.
params.sqlBackfill =
    'SELECT b.ID, cd.accession, cd.rank, a.sequence, s.sequence, COALESCE(s.topology, \'linear\'), ' +
        featsSql('blast_ref_annotations', 'x.accession = cd.accession') + ', ' +
        featsSql('annotations', 'x.ID = a.ID AND x.path = a.path AND x.scaffold = a.scaffold') + ', ' +
    'a.path, a.scaffold, COALESCE(a.topology, \'linear\') ' +
    'FROM assemble b ' +
    'JOIN (SELECT ID, accession, MIN(rank) AS rank FROM blast_ref_candidates ' +
         'GROUP BY ID, accession) cd ON cd.ID = b.ID ' +
    'JOIN assemblies a   ON b.ID = a.ID AND a.ignore = 0 ' +
    'JOIN annotate c     ON c.ID = a.ID AND c.path = a.path AND c.scaffold = a.scaffold ' +
    'JOIN blast_ref_sequences s ON cd.accession = s.accession ' +
    'WHERE a.sequence IS NOT NULL ' +
      'AND c.annotate_switch = 2 ' +
      'AND NOT EXISTS (SELECT 1 FROM blast_ref_alignment al WHERE al.ID = b.ID AND al.path = a.path AND al.scaffold = a.scaffold AND al.accession = cd.accession)'

// Sample gene table for one scaffold from the curate annotations CSV (per path,
// 'contig' = ID.path.scaffold). Same 'gene|pos1|pos2|direction;...' format as featsSql.
// Quote-aware CSV split. Not splitCsv: the curate CSV's refHits column is JSON
// with RFC-4180 doubled quotes, which Nextflow's splitCsv rejects ("Invalid CSV
// value"). Verified 2026-11-07 on circ_fwd_annotations_1.1.csv.
def csvSplit(String l) {
    return l.split(',(?=(?:[^"]*"[^"]*")*[^"]*$)', -1).collect { it.replaceAll('^"|"$', '') }
}
def sampleFeats(annFile, scaffold) {
    def lines = annFile.readLines()
    if (lines.size() < 2) return ''
    def ix = [:]; csvSplit(lines[0]).eachWithIndex { c, i -> ix[c.trim()] = i }
    def need = ['contig', 'gene', 'pos1', 'pos2', 'direction']
    if (!need.every { ix.containsKey(it) }) return ''
    def maxIx = need.collect { ix[it] }.max()
    lines.drop(1).collect { csvSplit(it) }
        .findAll { f -> f.size() > maxIx && f[ix['contig']].endsWith(".${scaffold}") && f[ix['pos1']].isInteger() }
        .collect { f -> [f[ix['gene']], f[ix['pos1']], f[ix['pos2']], f[ix['direction']]].join('|').replace("'", '') }
        .join(';')
}
// One scaffold's record from the per-path curate FASTA ('>ID.path.scaffold topology').
// Returns [sequence, topology]; falls back to every record joined when no header
// names the scaffold (single-record files written before headers carried it).
def scaffoldRecord(fastaFile, scaffold) {
    def recs = fastaFile.splitFasta(record: [id: true, desc: true, seqString: true])
    def r = recs.find { it.id.endsWith(".${scaffold}") }
    if (r) return [r.seqString, r.desc ?: 'linear']
    return [recs.collect { it.seqString }.join(''), 'linear']
}

params.sqlWriteAlignment = '''INSERT OR REPLACE INTO blast_ref_alignment
    (ID, path, scaffold, accession, aligned_sample, aligned_ref, rotation, ref_length, ref_start, strand, sample_rotation, time_stamp)
    VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)'''

// Mark the 'BLAST Ref Align' field (poor_blast_ref) = 'failed' and add an [align]
// note when blast_ref_align fails for a sample (errorStrategy 'ignore' suppresses
// task failure; we detect missing output downstream). This is a WF2 step that runs
// after the sample is already complete/locked in WF1, so the assemble_switch state
// is intentionally left unchanged: the assembly and annotation are unaffected, only
// the whole-genome alignment view is unavailable.
params.refAlignFailedMsg = "Whole-genome BLAST reference alignment failed. Annotation and curation are unaffected; the alignment view in the Assemble details will be unavailable."
params.sqlWriteBlastRefAlignFailed = "UPDATE assemble SET " +
    "poor_blast_ref = 'failed', " +
    "assemble_notes = ${appendTaggedNoteSql('[align]', params.refAlignFailedMsg)} " +
    "WHERE ID = ?"

workflow BLAST_REF_ALIGN {
    take:
        validated   // (id, path, scaffold) - gates timing: fires after VALIDATE completes
        curate_out  // (id, path, scaffold, annotations, assembly_fasta, coverage, work_dir) [per unit]

    main:
        // Reference data (sequence, topology, gene table) available at WF2 startup from WF1.
        // combine(by:0) below buffers this finite fromQuery side, so no deadlock.
        channel.fromQuery(params.sqlReadRef, db: 'sqlite')
            .map { row -> tuple(row[0], row[1], (row[2] as Integer) == 1, row[3], row[4].toString(), row[5].toString()) }
            .filter { id, accession, is_top, ref_seq, topology, ref_feats -> ref_seq }
            .set { ref_ch }  // (id, accession, is_top, ref_seq, topology, ref_feats)

        // Newly-validated units: gate on `validated` so alignment runs after
        // VALIDATE completes, then align EACH scaffold's own curated assembly to
        // every candidate reference. blast_ref_alignment is keyed
        // (ID, path, scaffold, accession).
        curate_out
            .map { tuple(it[0], it[1], it[2], it[3], it[4]) }   // (id, path, scaffold, annotations, assembly_fasta)
            .set { curate_assembly }
        validated
            .map { tuple(it[0], it[1], it[2]) }          // (id, path, scaffold)
            .unique()
            .join(curate_assembly, by: [0, 1, 2])
            .combine(ref_ch, by: 0)                       // ref_ch keyed by id -> broadcast to every scaffold
            .map { id, path, scaffold, annotations, assembly_fasta, accession, is_top, ref_seq, topology, ref_feats ->
                def rec = scaffoldRecord(assembly_fasta, scaffold)
                tuple(id, path, scaffold, accession, is_top, rec[0], ref_seq, topology, ref_feats,
                      sampleFeats(annotations, scaffold), rec[1])
            }
            .set { new_align_in }

        // Backfill: already-validated units missing an alignment row for a candidate.
        channel.fromQuery(params.sqlBackfill, db: 'sqlite')
            .map { row ->
                def accession = row[1]?.toString()
                def is_top    = (row[2] as Integer) == 1
                def asm_seq   = row[3]?.toString()
                def ref_seq   = row[4]?.toString()
                if (!accession || !asm_seq || !ref_seq) return null
                if (!ref_seq.matches('[ACGTNacgtnRYSWKMBDHVryswkmbdhv]+')) return null
                tuple(row[0], row[8], row[9], accession, is_top, asm_seq, ref_seq,
                      row[5].toString(), row[6].toString(), row[7].toString().replace("'", ''), row[10].toString())
            }
            .filter { it != null }
            .set { backfill_ch }

        new_align_in.mix(backfill_ch).set { align_in }

        // Track top-hit (id, path, scaffold, accession) entering align so a failed
        // TOP-ref alignment can flag the sample. Scaffold is in the key so one
        // scaffold's success does not mask a sibling's top-ref failure.
        align_in
            .filter { it[4] }
            .map    { tuple("${it[0]}|${it[1]}|${it[2]}|${it[3]}".toString(), it[0]) }
            .set { top_align_keys }

        blast_ref_align(
            align_in.map { it[0..3] + it[5..-1] }   // drop is_top
        ).set { align_out }

        // Parse the one-row CSV and write to blast_ref_alignment table
        align_out
            .map { id, path, scaffold, accession, csv_file ->
                def rows = csv_file.splitCsv(header: true, quote: '"')
                if (!rows) return null
                def r = rows[0]
                tuple(id, path, scaffold, accession, r.aligned_sample, r.aligned_ref, r.rotation.toLong(),
                      r.ref_length.toLong(), r.ref_start.toLong(), r.strand, r.sample_rotation.toLong(),
                      java.time.Instant.now().getEpochSecond())
            }
            .filter { it != null }
            .sqlInsert(statement: params.sqlWriteAlignment, db: 'sqlite')

        // Detect TOP-ref align failures: top (id, path, scaffold, accession) that
        // entered but produced no CSV output.
        align_out
            .map { id, path, scaffold, accession, csv_file -> tuple("${id}|${path}|${scaffold}|${accession}".toString(), true) }
            .set { succeeded_align_keys }

        top_align_keys
            .join(succeeded_align_keys, remainder: true)
            .filter { key, id, success_flag -> success_flag == null }
            .map    { key, id, success_flag -> tuple(id) }
            .sqlInsert(statement: params.sqlWriteBlastRefAlignFailed, db: 'sqlite')
}
