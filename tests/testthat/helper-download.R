# filtered_txi, dds y zscore_mat no vienen en el paquete: se descargan
# (ver ?download_data). Las pruebas que los usan se saltan sin internet y en CRAN
# (skip_if_offline() ya incluye skip_on_cran()).
local_download <- function(name) {
  skip_if_offline()
  download_data(name)
}
