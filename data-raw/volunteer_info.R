# Información clínica por voluntaria (una fila por persona).
# Fuente: tabla de reclutamiento del estudio, transcrita en data-raw/volunteer_info.csv
# Ejecutar desde la raíz del paquete: source("data-raw/volunteer_info.R")

raw <- read.csv("data-raw/volunteer_info.csv", na.strings = c("", "NA"))
is_sle <- raw$Group == "SLE"

# En pacientes, una celda vacía de tratamiento = no lo recibe (FALSE);
# en controles no aplica (NA)
took <- function(x) ifelse(is_sle, !is.na(x), NA)

# Meses en inglés con month.name para no depender del idioma del sistema
month <- match(sub(" .*", "", raw$Time_of_recruitment), month.name)
year  <- sub(".* ", "", raw$Time_of_recruitment)

volunteer_info <- data.frame(
  Volunteer_ID           = raw$Volunteer_ID,
  Volunteer              = raw$Volunteer,
  Group                  = factor(raw$Group, levels = c("Ctrl", "SLE")),
  Age_years              = raw$Age_years,
  Years_with_diagnosis   = raw$Years_with_diagnosis,
  Age_at_diagnosis_years = raw$Age_at_diagnosis_years,
  SLEDAI                 = raw$SLEDAI,
  Nephritis              = factor(raw$Nephritis, levels = c("No", "Yes")),
  Treatment_GCC          = took(raw$Treatment_GCC),
  GCC_drug               = sub("^(PD[NS]).*", "\\1", raw$Treatment_GCC),
  GCC_dose_mg            = as.numeric(ifelse(grepl("[0-9]", raw$Treatment_GCC),
                                             gsub("[^0-9.]", "", raw$Treatment_GCC), NA)),
  Treatment_AM           = took(raw$Treatment_AM),
  Treatment_MMF          = took(raw$Treatment_MMF),
  Treatment_AZT          = took(raw$Treatment_AZT),
  Treatment_MTT          = took(raw$Treatment_MTT),
  Time_of_recruitment    = as.Date(sprintf("%s-%02d-01", year, month))
)
volunteer_info <- volunteer_info[order(volunteer_info$Volunteer), ]
rownames(volunteer_info) <- NULL

# Verificaciones
load("data/metadata_cleaned.rda")
stopifnot(
  nrow(volunteer_info) == 33,
  !anyDuplicated(volunteer_info$Volunteer_ID),
  !anyNA(volunteer_info$Time_of_recruitment),
  setequal(volunteer_info$Volunteer_ID, sub("_.*", "", metadata_cleaned$sample_ID)),
  with(subset(volunteer_info, Group == "SLE"),
       all(Age_years - Years_with_diagnosis == Age_at_diagnosis_years))
)

usethis::use_data(volunteer_info, overwrite = TRUE, compress = "xz")
