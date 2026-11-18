# Initialize a new user-annotation project database

Builds the user-assembly schema with \[new_db_userAsmb()\] (mitogenome
search, circularization, and scaffold joining off), then parses each
sample's annotation file into the \`user_annotations\` table and the
gene names found into \`gene_name_map\`. The annotation tools are
switched off in the default \`annotate_opts\` set.

## Usage

``` r
new_db_userAnnot(
  db_path = "./.sqlite",
  mapping_fn = NULL,
  mapping_id = "ID",
  assembly_path = NULL,
  annotation_path = NULL,
  ...
)
```

## Arguments

- db_path:

  Path to the new database file

- mapping_fn:

  Path to the mapping CSV file. Must contain columns "ID", "Taxon",
  "Assembly", and "Annotation", plus "R1" and "R2" unless \`no_raw_data
  = TRUE\`.

- mapping_id:

  Column name of the mapping file to use as the primary key

- assembly_path:

  Directory holding the assembly FASTA files

- annotation_path:

  Directory holding the annotation files

- ...:

  Additional arguments passed to \[new_db_userAsmb()\].

## Value

The import summary (invisibly).
