# Which side of a synteny alignment to rotate, and by how much

A circular reference is rotated to the sample's first gene. When the
reference is linear but the sample is circular, the sample's alignment
copy is rotated to the reference's first gene instead (same rule, roles
swapped); the stored sample sequence is untouched and the synteny viewer
maps coordinates back. Two linear sequences: neither moves.

## Usage

``` r
alignment_rotations(
  ref_feats,
  sample_feats,
  ref_topology = "linear",
  sample_topology = "linear"
)
```

## Arguments

- ref_feats, sample_feats:

  data frames with \`gene\`, \`pos1\`, \`pos2\`, \`direction\`, or the
  compact \`gene\|pos1\|pos2\|direction;...\` string the Nextflow
  workflow builds from those columns.

- ref_topology:

  "circular" or anything else (no rotation).

- sample_topology:

  "circular" or anything else.

## Value

list(rotation, sample_rotation), both 0-based offsets.
