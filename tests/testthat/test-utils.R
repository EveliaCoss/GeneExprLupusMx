# Pruebas de check_pkgs(): mensaje claro cuando falta un paquete

test_that("check_pkgs passes when the packages are installed", {
  expect_true(check_pkgs("stats", "f"))
  expect_identical(pkgs_missing(c("stats", "notAPackage123")), "notAPackage123")
})

test_that("check_pkgs says which package is missing and how to install it", {
  expect_error(check_pkgs(c("stats", "notAPackage123"), "my_fun"),
               "my_fun\\(\\) needs the package\\(s\\) 'notAPackage123'")
  expect_error(check_pkgs("notAPackage123", "my_fun"),
               'BiocManager::install\\(c\\("notAPackage123"\\)\\)')
})

test_that("each function asks for the packages it needs", {
  local_mocked_bindings(pkgs_missing = function(pkgs) pkgs)
  expect_error(create_heatmap(matrix(1), "black", NULL),
               "create_heatmap\\(\\) needs the package\\(s\\) 'ComplexHeatmap'")
  expect_error(get_heatmap_mat(NULL),
               "get_heatmap_mat\\(\\) needs the package\\(s\\) 'DESeq2', 'SummarizedExperiment'")
  expect_error(process_go("UP_GO_x.RData", "UP", tempdir(), tempdir()),
               "process_go\\(\\) needs the package\\(s\\) 'rrvgo'")
})
