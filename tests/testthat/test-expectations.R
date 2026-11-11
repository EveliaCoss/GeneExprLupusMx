# Pruebas de las expectativas propias (helper-expectations.R)

test_that("expect_de_table passes or fails once, with a clear message", {
  d <- de_results$monocyte
  genes <- d$gene
  expect_success(expect_de_table(d, genes))
  expect_failure(expect_de_table(as.list(d), genes), "to be a data frame")
  expect_failure(expect_de_table(d[-1], genes), "to have columns")
  expect_failure(expect_de_table(d, rev(genes)), "genes of filtered_txi")
  bad <- d
  bad$pvalue[1] <- 2
  expect_failure(expect_de_table(bad, genes), "between 0 and 1")
  bad <- d
  bad$padj[1] <- 0
  bad$pvalue[1] <- 0.5
  expect_failure(expect_de_table(bad, genes), "at least `pvalue`")
  bad <- d
  bad$baseMean[1] <- -1
  expect_failure(expect_de_table(bad, genes), "non-negative")
})

test_that("expect_go_table passes or fails once, with a clear message", {
  expect_success(expect_go_table(go_results))
  expect_failure(expect_go_table(as.list(go_results)), "to be a data frame")
  expect_failure(expect_go_table(go_results[-1]), "to have columns")
  bad <- go_results
  bad$term_id[1] <- NA
  expect_failure(expect_go_table(bad), "no missing values")
  bad <- go_results
  bad$p_value[1] <- 0
  expect_failure(expect_go_table(bad), "in \\(0, 1\\]")
  bad <- go_results
  bad$intersection_size[1] <- bad$term_size[1] + 1
  expect_failure(expect_go_table(bad), "at most")
})
