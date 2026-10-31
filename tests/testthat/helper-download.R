# filtered_txi, dds y zscore_mat no vienen en el paquete: se descargan
# (ver ?download_data). Las pruebas que los usan se saltan sin internet y en CRAN.
local_download <- function(name) {
  skip_on_cran()
  skip_if_offline()
  download_data(name)
}
