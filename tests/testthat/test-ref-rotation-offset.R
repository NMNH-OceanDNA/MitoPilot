# ref_rotation_offset(): circular reference rotates to the sample's first
# annotated gene; strand decides start vs end anchoring.

ref <- data.frame(
  gene = c("trnF", "rrnS", "cox1", "nad5"),
  pos1 = c(1L, 70L, 5000L, 12000L), pos2 = c(69L, 1000L, 6500L, 13800L),
  direction = c("+", "+", "+", "-"), stringsAsFactors = FALSE
)

test_that("linear reference is never rotated", {
  smp <- data.frame(gene = "cox1", pos1 = 1L, pos2 = 1500L, direction = "+")
  expect_identical(ref_rotation_offset(ref, smp, "linear"), 0L)
  expect_identical(ref_rotation_offset(ref, smp, NA_character_), 0L)
})

test_that("same strand anchors the first sample gene at the reference start", {
  smp <- data.frame(gene = c("cox1", "nad5"), pos1 = c(1L, 3000L), pos2 = c(1500L, 4800L),
                    direction = c("+", "-"))
  expect_identical(ref_rotation_offset(ref, smp, "circular"), 4999L)
})

test_that("opposite strand anchors the first sample gene at the reference end", {
  smp <- data.frame(gene = "cox1", pos1 = 1L, pos2 = 1500L, direction = "-")
  expect_identical(ref_rotation_offset(ref, smp, "circular"), 6500L)
  smp <- data.frame(gene = "nad5", pos1 = 10L, pos2 = 1800L, direction = "-")
  expect_identical(ref_rotation_offset(ref, smp, "circular"), 11999L)
})

test_that("a rotated circular sample still anchors on start_gene", {
  smp <- data.frame(gene = c("trnF", "rrnS"), pos1 = c(1L, 70L), pos2 = c(69L, 1000L),
                    direction = "+")
  expect_identical(ref_rotation_offset(ref, smp, "circular"), 0L)
})

test_that("genes missing from the reference are skipped; nothing shared gives 0", {
  smp <- data.frame(gene = c("orf314", "nad5"), pos1 = c(1L, 400L), pos2 = c(300L, 2000L),
                    direction = c("+", "-"))
  expect_identical(ref_rotation_offset(ref, smp, "circular"), 11999L)
  smp <- data.frame(gene = "orf314", pos1 = 1L, pos2 = 300L, direction = "+")
  expect_identical(ref_rotation_offset(ref, smp, "circular"), 0L)
  expect_identical(ref_rotation_offset(ref, smp[0, ], "circular"), 0L)
})

test_that("compact string form matches the data frame form", {
  r <- "trnF|1|69|+;rrnS|70|1000|+;cox1|5000|6500|+;nad5|12000|13800|-"
  expect_identical(ref_rotation_offset(r, "cox1|1|1500|-;nad5|3000|4800|+", "circular"), 6500L)
  expect_identical(ref_rotation_offset(r, "", "circular"), 0L)
  expect_identical(ref_rotation_offset(r, character(), "circular"), 0L)
  expect_identical(ref_rotation_offset(r, "junk;cox1|x|y|+", "circular"), 0L)
})

test_that("alignment_rotations rotates the circular side only", {
  smp <- data.frame(gene = c("trnF", "cox1"), pos1 = c(1L, 3000L), pos2 = c(69L, 4500L),
                    direction = "+", stringsAsFactors = FALSE)
  # circular ref: ref moves to the sample's first gene
  expect_equal(alignment_rotations(ref, smp, "circular", "circular"),
               list(rotation = 0L, sample_rotation = 0L))
  smp2 <- data.frame(gene = "cox1", pos1 = 1L, pos2 = 1500L, direction = "+")
  expect_equal(alignment_rotations(ref, smp2, "circular", "linear"),
               list(rotation = 4999L, sample_rotation = 0L))
  # linear ref, circular sample: the sample copy moves to the ref's first gene (trnF)
  smp3 <- data.frame(gene = c("cox1", "trnF"), pos1 = c(1L, 3000L), pos2 = c(1500L, 3068L),
                     direction = c("+", "+"))
  expect_equal(alignment_rotations(ref, smp3, "linear", "circular"),
               list(rotation = 0L, sample_rotation = 2999L))
  # ref's first gene on the opposite strand in the sample: anchor at the sample end
  smp4 <- data.frame(gene = c("cox1", "trnF"), pos1 = c(1L, 3000L), pos2 = c(1500L, 3068L),
                     direction = c("+", "-"))
  expect_equal(alignment_rotations(ref, smp4, "linear", "circular")$sample_rotation, 3068L)
  # both linear, or topology unknown: nothing moves
  expect_equal(alignment_rotations(ref, smp3, "linear", "linear"),
               list(rotation = 0L, sample_rotation = 0L))
  expect_equal(alignment_rotations(ref, smp3, NA_character_, NA_character_),
               list(rotation = 0L, sample_rotation = 0L))
})

test_that("compute_blast_ref_alignment rotates the sample copy and reports it", {
  set.seed(7)
  smp <- paste(sample(c("A", "C", "G", "T"), 300, TRUE), collapse = "")
  rotated <- paste0(substr(smp, 101, 300), substr(smp, 1, 100))
  f <- tempfile(fileext = ".csv")
  compute_blast_ref_alignment(smp, rotated, rotation = 0L, output_file = f, sample_rotation = 100L)
  res <- utils::read.csv(f, stringsAsFactors = FALSE)
  expect_equal(res$sample_rotation, 100L)
  expect_equal(gsub("-", "", res$aligned_sample), rotated)
  expect_equal(res$strand, "+")
})
