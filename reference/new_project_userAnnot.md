# Initialize new MitoPilot Project with user-provided assemblies and annotations

Each sample brings one single-sequence mitogenome FASTA and one
annotation file (GFF3, BED, or 5-column feature table). The annotation
tools are not run: the pipeline curates, validates, and exports the
imported features. Gene names found in the files are guessed at setup
and must be reviewed in the app before the pipeline runs.

## Usage

``` r
new_project_userAnnot(
  path = ".",
  mapping_fn = NULL,
  mapping_id = "ID",
  mapping_geome = "GEOME_BCID",
  fetch_geome = TRUE,
  mapping_gbif = "GBIF_ID",
  fetch_gbif = TRUE,
  mapping_ncbi = "BioSample",
  fetch_ncbi = TRUE,
  link_sources = FALSE,
  data_path = NULL,
  no_raw_data = FALSE,
  assembly_path = "NA",
  annotation_path = "NA",
  genetic_code = NULL,
  executor = c("local", "awsbatch", "slurm", "sge", "pbs", "lsf", "NMNH_Hydra",
    "NOAA_SEDNA"),
  container = paste0("macguigand/mitopilot:", utils::packageVersion("MitoPilot")),
  config = NULL,
  profile_dir = mitopilot_config_dir(),
  ncbi_api_key = NULL,
  Rproj = FALSE,
  force = FALSE,
  ...
)
```

## Arguments

- path:

  Path to the project directory (default = current working directory)

- mapping_fn:

  Path to a mapping file. Should be a csv that minimally includes an
  \`ID\` column with a unique identifier for each sample, a \`Taxon\`
  column containing taxonomic information for each sample, and columns
  \`R1\` and \`R2\` specifying the names of the raw paired read inputs
  (unless \`no_raw_data = TRUE\`), an \`Assembly\` column containing
  names of assembly fasta files (one sequence each), and an
  \`Annotation\` column containing names of annotation files (\`.gff\`,
  \`.gff3\`, \`.bed\`, or \`.tbl\`). An optional \`Topology\` column
  ("circular" or "linear", blank = linear) declares each assembly's
  topology; a feature spanning the origin needs "circular". May include
  additional columns with other sample metadata.

- mapping_id:

  The name of the column in the mapping file that contains the unique
  sample identifiers (default = "ID").

- mapping_geome:

  Name of the mapping-file column holding GEOME BCIDs (optional). Stored
  as \`GEOME_BCID\`. See \`vignette("Specimen-Metadata")\`. Passed to
  \`new_db_userAsmb()\`.

- fetch_geome:

  Fetch GEOME metadata for samples with a BCID during setup (default
  TRUE). Set FALSE when offline and run \[fetch_geome()\] later. Passed
  to \`new_db_userAsmb()\`.

- mapping_gbif:

  Name of the mapping-file column holding GBIF occurrence IDs
  (optional). Stored as \`GBIF_ID\`. See
  \`vignette("Specimen-Metadata")\`. Passed to \`new_db_userAsmb()\`.

- fetch_gbif:

  Fetch GBIF metadata for samples with a GBIF ID during setup (default
  TRUE). Set FALSE when offline and run \[fetch_gbif()\] later.

- mapping_ncbi:

  Name of the mapping-file column holding NCBI BioSample or SRA
  accessions (optional). Must not be the sample ID column. Stored as
  \`BioSample\`. See \`vignette("Specimen-Metadata")\`.

- fetch_ncbi:

  Fetch NCBI metadata for samples with a BioSample value during setup
  (default TRUE). Set FALSE when offline and run \[fetch_ncbi()\] later.

- link_sources:

  Follow links between GEOME, GBIF, and NCBI records to fill in IDs a
  sample does not have yet (default FALSE). Saved as the project
  setting. See \`vignette("Specimen-Metadata")\`. Passed to
  \`new_db_userAsmb()\`.

- data_path:

  Path to the directory where the raw data is located. Can be a AWS s3
  bucket even if not using AWS for pipeline execution. Not required when
  \`no_raw_data = TRUE\`.

- no_raw_data:

  (logical) Run without raw sequence data (default = FALSE). When TRUE,
  the read-mapping coverage step is skipped: \`data_path\` is ignored,
  coverage/depth statistics are left empty, and annotation coverage
  trimming is disabled. Use this to annotate an assembly you already
  have.

- assembly_path:

  Path to the directory where the mitogenome assemblies are located. Can
  be a AWS s3 bucket even if not using AWS for pipeline execution.

- annotation_path:

  Path to the directory where the annotation files are located. Each
  file must name one seqid equal to the first word of its FASTA header.
  Files are copied to \`user_annotations/raw/\` in the project.

- genetic_code:

  Optional NCBI translation table override. Default \`NULL\`
  auto-selects from each sample's curation ruleset; a number sets a
  project-wide override.
  https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi

- executor:

  The executor to use for running the nextflow pipeline. May be a
  built-in template ("local" (default), "awsbatch", "slurm", "sge",
  "pbs", "lsf", "NMNH_Hydra", "NOAA_SEDNA") or the name of a saved
  cluster profile created with \[generate_config()\]. See
  \[list_configs()\] for available names.

- container:

  The docker container to use for pipeline execution.

- config:

  (optional) provide a path to an existing custom nextflow config file.
  If not provided a config file template will be created based on the
  specified executor.

- profile_dir:

  Directory searched for saved cluster profiles when resolving
  \`executor\` (default \[mitopilot_config_dir()\]).

- ncbi_api_key:

  Optional NCBI API key string. Used to raise NCBI request rate limits
  for the remote BLAST + GenBank fetch steps. See
  \<https://www.ncbi.nlm.nih.gov/datasets/docs/v2/api/api-keys/\>. May
  be left empty and edited later in \`.config\`
  (\`params.ncbi_api_key\`).

- Rproj:

  (logical) Initialize and open an RStudio project in the project
  directory (default = FALSE). This option has no effect if not running
  interactively in RStudio.

- force:

  (logical) Force recreating of existing project database and config
  files (default = FALSE).

- ...:

  Additional arguments passed as default processing parameters to
  \[new_db_userAsmb()\] through \[new_db_userAnnot()\]. Mitogenome
  search, circularization, and scaffold joining are always off.
