ua_gff <- function(rows, seqid = "mito_contig") {
  c("##gff-version 3", paste(seqid, "src", rows, sep = "\t"))
}

# ann: named list of annotation file lines, one per sample (UA_LINEAR, then UA_CIRCULAR)
ua_project <- function(td, ann) {
  ids <- c("UA_LINEAR", "UA_CIRCULAR")[seq_along(ann)]
  asm <- file.path(td, "asm")
  ann_dir <- file.path(td, "ann")
  dir.create(asm)
  dir.create(ann_dir)
  file.copy(file.path(app_sys("test_data/assemblies"), paste0(ids, ".fasta")), asm)
  for (fn in names(ann)) writeLines(ann[[fn]], file.path(ann_dir, fn))
  m <- data.frame(ID = ids, Taxon = "Danio rerio", Assembly = paste0(ids, ".fasta"),
                  Topology = c("linear", "circular")[seq_along(ids)], Annotation = names(ann))
  utils::write.csv(m, file.path(td, "mapping.csv"), row.names = FALSE)
  proj <- file.path(td, "proj")
  suppressMessages(new_project_userAnnot(
    path = proj, mapping_fn = file.path(td, "mapping.csv"), assembly_path = asm,
    annotation_path = ann_dir, no_raw_data = TRUE, Rproj = FALSE, executor = "local"))
  proj
}
