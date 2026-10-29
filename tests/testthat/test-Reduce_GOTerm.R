# Pruebas de processGO() y Reduce_GOTerm()

skip_if_no_go_pkgs <- function() {
  for (p in c("rrvgo", "org.Hs.eg.db", "GO.db", "treemap", "dplyr", "ggplot2")) {
    skip_if_not_installed(p)
  }
}

# Guarda una tabla de enriquecimiento con las columnas que espera processGO()
save_go_file <- function(dir, file, n = 60) {
  d <- go_celltype_results[go_celltype_results$query == "SLE_UP_moDC" &
                             go_celltype_results$source == "GO:BP", ]
  d <- utils::head(d[order(d$p_value), ], n)
  go_table <- data.frame(go_id = d$term_id, term = d$term_name, p.adjust = d$p_value)
  path <- file.path(dir, file)
  save(go_table, file = path)
  path
}

local_dirs <- function(env = parent.frame()) {
  root <- withr::local_tempdir(.local_envir = env)
  dirs <- file.path(root, c("in", "out", "fig"))
  for (d in dirs) dir.create(d)
  stats::setNames(as.list(dirs), c("indir", "outdir", "figdir"))
}

test_that("processGO reduces GO terms and saves results and plots", {
  skip_if_no_go_pkgs()
  dirs <- local_dirs()
  f <- save_go_file(dirs$indir, "UP_GO_SLE_moDC.RData")

  res <- suppressMessages(processGO(f, "UP", dirs$outdir, dirs$figdir))
  expect_null(res)
  expect_setequal(list.files(dirs$outdir),
                  c("reducedTerms_simMatrix_UP_SLE_moDC.RData",
                    "geneListAndGOterms_UP_SLE_moDC.RData"))
  expect_setequal(list.files(dirs$figdir),
                  c("GOTerm_UP_Treemap_SLE_moDC.png", "GOTerm_UP_ReducedBar_SLE_moDC.png"))

  env <- new.env()
  load(file.path(dirs$outdir, "reducedTerms_simMatrix_UP_SLE_moDC.RData"), envir = env)
  expect_true(is.matrix(env$simMatrix))
  expect_true(all(env$reducedTerms$size > 0))
  load(file.path(dirs$outdir, "geneListAndGOterms_UP_SLE_moDC.RData"), envir = env)
  expect_false(anyNA(env$joined_data$parentTerm))
  expect_true(all(c("go_id", "term", "p.adjust", "parentTerm") %in% names(env$joined_data)))
})

test_that("processGO validates direction and the file contents", {
  skip_if_no_go_pkgs()
  dirs <- local_dirs()
  f <- save_go_file(dirs$indir, "UP_GO_x.RData")
  expect_error(processGO(f, "LEFT", dirs$outdir, dirs$figdir), "should be one of")

  a <- data.frame(go_id = "GO:0006955")
  b <- a
  two <- file.path(dirs$indir, "UP_GO_two.RData")
  save(a, b, file = two)
  expect_error(suppressMessages(processGO(two, "UP", dirs$outdir, dirs$figdir)),
               "exactamente un objeto")
})

test_that("processGO skips files without go_id or with few GO terms", {
  skip_if_no_go_pkgs()
  dirs <- local_dirs()
  no_id <- data.frame(term = "immune response", p.adjust = 0.01)
  f_no_id <- file.path(dirs$indir, "UP_GO_noid.RData")
  save(no_id, file = f_no_id)
  msgs <- capture_messages(res <- processGO(f_no_id, "UP", dirs$outdir, dirs$figdir))
  expect_true(any(grepl("No 'go_id' column found", msgs, fixed = TRUE)))
  expect_null(res)

  f_few <- save_go_file(dirs$indir, "DOWN_GO_few.RData", n = 10)
  msgs <- capture_messages(processGO(f_few, "DOWN", dirs$outdir, dirs$figdir))
  expect_true(any(grepl("fewer than 40 unique GO terms", msgs, fixed = TRUE)))
  expect_length(list.files(dirs$outdir), 0)
  expect_length(list.files(dirs$figdir), 0)
})

test_that("processGO skips GO IDs that are not in org.Hs.eg.db", {
  skip_if_no_go_pkgs()
  dirs <- local_dirs()
  bad <- data.frame(go_id = sprintf("GO:99%05d", 1:45), term = paste("term", 1:45),
                    p.adjust = 0.01)
  f <- file.path(dirs$indir, "UP_GO_bad.RData")
  save(bad, file = f)
  msgs <- capture_messages(suppressWarnings(processGO(f, "UP", dirs$outdir, dirs$figdir)))
  expect_true(any(grepl("Invalid simMatrix structure", msgs, fixed = TRUE)))
  expect_length(list.files(dirs$outdir), 0)
})

test_that("Reduce_GOTerm processes every UP_GO_ and DOWN_GO_ file in a folder", {
  skip_if_no_go_pkgs()
  dirs <- local_dirs()
  save_go_file(dirs$indir, "UP_GO_first.RData", n = 10)
  save_go_file(dirs$indir, "DOWN_GO_second.RData", n = 10)
  save_go_file(dirs$indir, "other_file.RData", n = 10)

  msgs <- capture_messages(res <- Reduce_GOTerm(dirs$indir, dirs$outdir, dirs$figdir))
  expect_null(res)
  expect_true(any(grepl("Dataset: first", msgs, fixed = TRUE)))
  expect_true(any(grepl("Dataset: second", msgs, fixed = TRUE)))
  expect_false(any(grepl("other_file", msgs, fixed = TRUE)))
})
