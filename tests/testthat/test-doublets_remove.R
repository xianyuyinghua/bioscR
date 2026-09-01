test_that("doublets_remove requires a Seurat object", {
  expect_error(
    doublets_remove(matrix(1, nrow = 2, ncol = 2)),
    "'seu' must be a Seurat object.",
    fixed = TRUE
  )
})

test_that("basic Seurat input arguments are validated", {
  skip_if_not_installed("SeuratObject")

  counts <- matrix(
    c(1, 0, 2, 1, 3, 0),
    nrow = 3,
    dimnames = list(
      paste0("gene", seq_len(3)),
      paste0("cell", seq_len(2))
    )
  )
  seu <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_error(
    doublets_remove(seu, assay = "missing"),
    "Assay 'missing' is absent from the Seurat object.",
    fixed = TRUE
  )
  expect_error(
    doublets_remove(seu, remove_doublets = NA),
    "'remove_doublets' must be TRUE or FALSE.",
    fixed = TRUE
  )
  expect_error(
    doublets_remove(seu, remove_doublets = c(TRUE, FALSE)),
    "'remove_doublets' must be TRUE or FALSE.",
    fixed = TRUE
  )
})

test_that("sample metadata is validated", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scDblFinder")

  counts <- matrix(
    c(1, 0, 2, 1, 3, 0),
    nrow = 3,
    dimnames = list(
      paste0("gene", seq_len(3)),
      paste0("cell", seq_len(2))
    )
  )
  seu <- SeuratObject::CreateSeuratObject(counts = counts)

  expect_error(
    doublets_remove(seu, sample_col = "missing", verbose = FALSE),
    "Metadata column 'missing' is absent from the Seurat object.",
    fixed = TRUE
  )

  seu$sample <- c("sample1", NA_character_)
  expect_error(
    doublets_remove(seu, sample_col = "sample", verbose = FALSE),
    "Metadata column 'sample' contains missing/empty values.",
    fixed = TRUE
  )
})

test_that("predictions are added and cells are optionally filtered", {
  skip_if_not_installed("SeuratObject")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scDblFinder")

  counts <- matrix(
    rep(c(1, 0, 2, 3), 3),
    nrow = 4,
    dimnames = list(
      paste0("gene", seq_len(4)),
      paste0("cell", seq_len(3))
    )
  )
  seu <- SeuratObject::CreateSeuratObject(counts = counts)
  seu$sample <- c("sample1", "sample1", "sample2")

  fake_scDblFinder <- function(sce, samples = NULL, verbose = TRUE, ...) {
    SummarizedExperiment::colData(sce)$scDblFinder.class <- c(
      "singlet", "doublet", "singlet"
    )
    SummarizedExperiment::colData(sce)$scDblFinder.score <- c(0.1, 0.9, 0.2)
    sce
  }

  testthat::local_mocked_bindings(
    scDblFinder = fake_scDblFinder,
    .package = "scDblFinder"
  )

  annotated <- doublets_remove(
    seu,
    remove_doublets = FALSE,
    verbose = FALSE
  )
  expect_equal(colnames(annotated), colnames(seu))
  expect_equal(
    as.character(annotated$scDblFinder.class),
    c("singlet", "doublet", "singlet")
  )
  expect_equal(as.logical(annotated$is_doublet), c(FALSE, TRUE, FALSE))
  expect_equal(as.logical(annotated$is_singlet), c(TRUE, FALSE, TRUE))

  filtered <- doublets_remove(
    seu,
    remove_doublets = TRUE,
    verbose = FALSE
  )
  expect_equal(colnames(filtered), c("cell1", "cell3"))
  expect_true(all(filtered$is_singlet))
})

test_that("multiple Seurat v5 count layers are joined", {
  skip_if_not_installed("SeuratObject", minimum_version = "5.0.0")
  skip_if_not_installed("SingleCellExperiment")
  skip_if_not_installed("SummarizedExperiment")
  skip_if_not_installed("scDblFinder")

  counts_a <- matrix(
    c(1, 0, 2, 1),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("cell1", "cell2"))
  )
  counts_b <- matrix(
    c(3, 1, 0, 2),
    nrow = 2,
    dimnames = list(c("gene1", "gene2"), c("cell3", "cell4"))
  )
  seu <- SeuratObject::CreateSeuratObject(
    counts = list(sample1 = counts_a, sample2 = counts_b)
  )
  seu$sample <- c("sample1", "sample1", "sample2", "sample2")

  expect_length(
    grep(
      "^counts($|\\.)",
      SeuratObject::Layers(seu[["RNA"]], search = NA),
      value = TRUE
    ),
    2L
  )

  fake_scDblFinder <- function(sce, samples = NULL, verbose = TRUE, ...) {
    expect_identical(colnames(sce), colnames(seu))
    expect_identical(as.character(samples), as.character(seu$sample))
    SummarizedExperiment::colData(sce)$scDblFinder.class <- rep(
      "singlet", ncol(sce)
    )
    sce
  }

  testthat::local_mocked_bindings(
    scDblFinder = fake_scDblFinder,
    .package = "scDblFinder"
  )

  result <- expect_message(
    doublets_remove(seu, remove_doublets = FALSE, verbose = TRUE),
    "Joining 2 count layers"
  )

  expect_identical(
    SeuratObject::Layers(result[["RNA"]], search = NA),
    "counts"
  )
  expect_equal(ncol(result), 4L)
  expect_true(all(result$is_singlet))
})
