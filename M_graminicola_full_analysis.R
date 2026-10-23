###############################################################################
# Meloidogyne graminicola species distribution model - complete analysis
# One script: present-day model, revision analyses, future projection.
#
# HOW TO RUN
#   1. Save this file anywhere (for example in submission_code).
#   2. Check ROOT_DIR below. It is already set to your folder.
#   3. Run everything:   source("M_graminicola_full_analysis.R", echo = TRUE)
#      or run part by part in RStudio (each PART starts with a banner; use
#      the document outline, Ctrl+Shift+O, to jump between sections).
#
# PARTS
#   PART 0  Setup. Sets the project folder (ROOT_DIR/submission_code) and, on
#           the first run, copies the files of the original analysis that are
#           already on this computer into data/ and env_data/ under the names
#           the script expects (EPPO database and tables, geocoded EPPO records,
#           GBIF records, WorldClim and cropland rasters, CMIP6 rasters).
#           Nothing is downloaded again if it is already there.
#   PART 1  Present-day model (GBIF + EPPO, Maxent, cross-validation, figures).
#           Writes outputs/analysis_handoff.rds, read by Parts 2 and 3.
#   PART 2  Revision analyses for the editor's comments (location weights,
#           withheld validation, block sizes, MESS, rice production systems).
#           Writes outputs/revision/ and outputs/revision_handoff.rds.
#           BEFORE THE FULL RUN: check outputs/revision/eppo_unit_crosswalk_DRAFT.csv
#           and put corrections in data/eppo_unit_crosswalk_manual.csv.
#   PART 3  Future projection (5 GCMs x 2 SSPs x 2 periods), projecting the
#           revised model by default (FUT$model_source = "revised").
#
# To skip a part that has already been run, set its switch to FALSE.
# Each part reloads what it needs from outputs/, so later parts work alone
# once earlier parts have been run once.
#
# Requirements: R >= 4.3; packages are checked at the start of each part.
# Downloads still needed: MapSPAM rice and GADM boundaries (Part 2, ~1 GB) and
# the two CMIP6 models / 2061-2080 period not already on disk (Part 3).
# Chunks marked "(not run ...)" are commented out on purpose, as in the
# original .qmd files.
###############################################################################

# Folder that holds submission_code/, Supplementary files/, wc_data/ etc.
ROOT_DIR <- "C:/Users/davit/Documents/M. graminicola"

RUN_PART <- c(present_day = FALSE, revision = TRUE, future = TRUE)   # Part 1 already done; set TRUE to redo


###############################################################################
# PART 0 - SETUP: PROJECT FOLDER AND FILES FROM THE ORIGINAL RUN ##############
###############################################################################

PROJ_DIR <- file.path(ROOT_DIR, "submission_code")
SUPP_DIR <- file.path(ROOT_DIR, "Supplementary files")
if (!dir.exists(PROJ_DIR)) stop("Project folder not found: ", PROJ_DIR,
                               "\nCheck ROOT_DIR at the top of this script.", call. = FALSE)
setwd(PROJ_DIR)
message("Working directory: ", getwd())
for (d in c("data", "outputs", "figures", "env_data/climate/wc2.1_2.5m", "env_data/landuse"))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)

# Copy a file only if the target is missing, so nothing is ever overwritten.
setup_copy <- function(from, to) {
  if (file.exists(to)) return(invisible("present"))
  if (!file.exists(from)) { message("  [setup] not found, skipped: ", from); return(invisible("missing")) }
  ok <- file.copy(from, to)
  message(sprintf("  [setup] %s %s", if (ok) "copied" else "FAILED to copy", basename(to)))
  invisible(if (ok) "copied" else "failed")
}

# 1. EPPO database and distribution table from the original run
setup_copy(file.path(ROOT_DIR, "sqlite_all.zip"), "data/sqlite_all.zip")
setup_copy(file.path(SUPP_DIR, "eppocodes_all.sqlite"), "data/eppocodes_all.sqlite")

# 2. EPPO records. The EPPO tables saved in "Supplementary files" come from a
#    later query (146 rows, Pacific islands, different absences) and are NOT the
#    71-row set used for the submitted model. The exact training set of the
#    submitted model is data/M_graminicola_occurrence_points_FROZEN.csv, so the
#    EPPO rows are rebuilt here by matching the geocoded table to the frozen
#    points (same Nominatim coordinates), plus the two archipelago national rows
#    that fell in sea cells and the confirmed absence (Slovakia).
f_geo_target  <- "data/Meloidogyne_graminicola_EPPO_geocoded.csv"
f_dist_target <- "data/Meloidogyne_graminicola_EPPO_distribution.csv"
f_frozen      <- "data/M_graminicola_occurrence_points_FROZEN.csv"
f_pres <- file.path(SUPP_DIR, "Meloidogyne_graminicola_EPPO_presence_geocoded.csv")
f_abs  <- file.path(SUPP_DIR, "Meloidogyne_graminicola_EPPO_absence_geocoded.csv")
f_dist <- file.path(SUPP_DIR, "Meloidogyne_graminicola_EPPO_distribution.csv")

std_cols <- function(df) {
  names(df) <- make.names(names(df))
  rn <- function(df, want, candidates) {
    if (want %in% names(df)) return(df)
    hit <- names(df)[tolower(names(df)) %in% tolower(candidates)]
    if (length(hit)) names(df)[names(df) == hit[1]] <- want
    df
  }
  df <- rn(df, "eppo_lat", c("lat", "latitude", "decimalLatitude", "y"))
  df <- rn(df, "eppo_lon", c("lon", "long", "lng", "longitude", "decimalLongitude", "x"))
  df <- rn(df, "Status",   c("status"))
  df <- rn(df, "country",  c("country.name", "countryname"))
  df <- rn(df, "state",    c("state.name", "statename", "subdivision"))
  if (!"state" %in% names(df)) df$state <- NA_character_
  df$state[!is.na(df$state) & !nzchar(trimws(df$state))] <- NA_character_
  df
}

SETUP_VERSION <- 3
rebuild_eppo <- !file.exists(f_geo_target) || {
  g0 <- utils::read.csv(f_geo_target)
  !("setup_version" %in% names(g0)) || !all(g0$setup_version == SETUP_VERSION)
}
if (rebuild_eppo && file.exists(f_pres) && file.exists(f_frozen)) {
  geo_all <- std_cols(utils::read.csv(f_pres, stringsAsFactors = FALSE))
  if ("eppocode" %in% names(geo_all) && !all(geo_all$eppocode %in% c("MELGGC", NA))) {
    stop("The EPPO table in 'Supplementary files' is for EPPO code ",
         paste(unique(geo_all$eppocode), collapse = ", "),
         " (ROTYRE = Rotylenchulus reniformis), not M. graminicola (MELGGC).\n",
         "Put the corrected Meloidogyne_graminicola_EPPO_geocoded.csv and ",
         "Meloidogyne_graminicola_EPPO_distribution.csv in submission_code/data/ ",
         "and run again.", call. = FALSE)
  }
  need <- c("country", "state", "Status", "eppo_lat", "eppo_lon")
  if (!all(need %in% names(geo_all)))
    stop("Geocoded EPPO file lacks columns. Found: ", paste(names(geo_all), collapse = ", "), call. = FALSE)
  geo_all <- geo_all[!grepl("^Absent", geo_all$Status), ]
  geo_all$level <- ifelse(is.na(geo_all$state), 0L, 1L)

  frozen <- utils::read.csv(f_frozen, stringsAsFactors = FALSE)
  fe     <- frozen[grepl("EPPO", frozen$Data_Source, ignore.case = TRUE), ]
  fe$level <- ifelse(grepl("^state", fe$Coordinate_Precision, ignore.case = TRUE), 1L, 0L)

  km <- function(lon1, lat1, lon2, lat2) {
    r <- pi / 180
    a <- sin((lat2 - lat1) * r / 2)^2 + cos(lat1 * r) * cos(lat2 * r) * sin((lon2 - lon1) * r / 2)^2
    2 * 6371 * asin(pmin(sqrt(a), 1))
  }
  # Coordinates are taken from the frozen file, so the training data are exactly
  # those of the submitted model. Country/state names (needed only by the
  # revision analysis) come from the nearest row of the same level in the saved
  # geocoded table; the revision crosswalk (Part 2) checks them again.
  rows <- lapply(seq_len(nrow(fe)), function(i) {
    cand <- geo_all[geo_all$level == fe$level[i], ]
    if (!nrow(cand)) cand <- geo_all
    d <- km(fe$Longitude[i], fe$Latitude[i], cand$eppo_lon, cand$eppo_lat)
    j <- which.min(d)
    r <- cand[j, , drop = FALSE]
    r$name_match_km <- round(d[j], 1)
    if (d[j] > 300) {                     # no close same-level row: keep country, mark state
      jj <- which.min(km(fe$Longitude[i], fe$Latitude[i], geo_all$eppo_lon, geo_all$eppo_lat))
      r[, c("country", intersect("country.code", names(r)))] <-
        geo_all[jj, c("country", intersect("country.code", names(geo_all)))]
      r$state <- if (fe$level[i] == 1L) "unmatched subdivision" else NA_character_
      r$Status <- "Present, no details"
    }
    r$eppo_lon <- fe$Longitude[i]; r$eppo_lat <- fe$Latitude[i]
    r
  })
  m <- do.call(rbind, rows)
  m$level <- NULL
  message(sprintf("  [setup] EPPO records rebuilt from the frozen training set: %d (names matched within 50 km: %d, 50-300 km: %d, >300 km: %d)",
                  nrow(m), sum(m$name_match_km <= 50), sum(m$name_match_km > 50 & m$name_match_km <= 300),
                  sum(m$name_match_km > 300)))
  if (any(m$name_match_km > 50))
    print(m[m$name_match_km > 50, c("country", "state", "eppo_lon", "eppo_lat", "name_match_km")])

  # the two national centroids dropped in the submitted run for falling in sea cells
  arch <- geo_all[is.na(geo_all$state) & geo_all$country %in% c("Indonesia", "Philippines"), ]
  arch <- arch[!duplicated(arch$country), ]
  arch$level <- NULL
  if (nrow(arch)) arch$name_match_km <- 0
  # the confirmed absence (Slovakia), the independent negative of the submitted run
  sk <- m[0, ][1, ]
  sk[1, c("country", "Status", "eppo_lat", "eppo_lon")] <-
    list("Slovakia", "Absent, confirmed by survey", 48.7411522, 19.4528646)
  if ("country.code" %in% names(sk)) sk$country.code <- "SK"
  if ("continent" %in% names(sk)) sk$continent <- "Europe"
  sk$state <- NA_character_
  geo <- rbind(m, arch[, names(m)], sk[, names(m)])
  geo$geocode_query <- ifelse(!is.na(geo$state), paste(geo$state, geo$country, sep = ", "), geo$country)
  geo$setup_version <- SETUP_VERSION
  utils::write.csv(geo, f_geo_target, row.names = FALSE)
  message(sprintf("  [setup] built %s: %d presence + %d absence rows (submitted run: 70 + 1)",
                  basename(f_geo_target), sum(!grepl("^Absent", geo$Status)),
                  sum(grepl("^Absent", geo$Status))))

  # Distribution table = the same rows (used for counts only)
  keep_cols <- intersect(c("continent", "country", "state", "country.code", "state.code", "Status"), names(geo))
  utils::write.csv(geo[, keep_cols], f_dist_target, row.names = FALSE)

  # Everything computed from earlier, wrong EPPO sets must be recomputed
  stale <- c(list.files("outputs", pattern = "^M_graminicola_", full.names = TRUE),
             list.files("outputs", pattern = "^(analysis_handoff|revision_handoff)\\.rds$", full.names = TRUE),
             list.files("outputs", pattern = "^Table_", full.names = TRUE))
  stale <- stale[!dir.exists(stale)]
  if (length(stale)) {
    file.remove(stale)
    message("  [setup] removed ", length(stale), " outputs computed from the earlier EPPO set")
  }
} else if (!file.exists(f_geo_target)) {
  message("  [setup] geocoded EPPO table or frozen points not found; Part 1 will geocode again (Nominatim).")
}

# 3. GBIF records from the original run, so GBIF is not queried again (a fresh
#    query now returns more records and would change the results).
setup_copy("data/Meloidogyne_graminicola_GBIF.csv", "data/gbif_raw.csv")
if (!file.exists("data/gbif_raw.csv"))
  setup_copy(file.path(SUPP_DIR, "Meloidogyne_graminicola_GBIF.csv"), "data/gbif_raw.csv")

# 4. Rasters already downloaded during the original run (several GB)
for (src in c(file.path(ROOT_DIR, "wc_data/climate/wc2.1_2.5m"),
              file.path(ROOT_DIR, "wc_future/climate/wc2.1_2.5m"))) {
  for (f in list.files(src, pattern = "\\.tif$", full.names = TRUE))
    setup_copy(f, file.path("env_data/climate/wc2.1_2.5m", basename(f)))
}
setup_copy(file.path(ROOT_DIR, "wc_data/landuse/WorldCover_cropland_30s.tif"),
           "env_data/landuse/WorldCover_cropland_30s.tif")

# 5. Report
setup_check <- c(
  `GBIF records (data/gbif_raw.csv)`           = file.exists("data/gbif_raw.csv"),
  `EPPO distribution table`                    = file.exists("data/Meloidogyne_graminicola_EPPO_distribution.csv"),
  `EPPO geocoded records`                      = file.exists(f_geo_target),
  `WorldClim BIO1 raster`                      = file.exists("env_data/climate/wc2.1_2.5m/wc2.1_2.5m_bio_1.tif"),
  `Cropland raster`                            = file.exists("env_data/landuse/WorldCover_cropland_30s.tif"))
print(setup_check)
if (!all(setup_check[1:3])) {
  stop("Input data from the original run are missing (see FALSE above). ",
       "Fix this before running Part 1, or Part 1 will query GBIF/EPPO afresh.", call. = FALSE)
}



###############################################################################
# PART 1 - PRESENT-DAY SDM ANALYSIS (M_graminicola_SDM_analysis.qmd) ##########
###############################################################################

if (RUN_PART[["present_day"]]) {



# Scope and how to read this document =====================================


# Three deviations between the published run and this script --------------


# Live services this analysis depends on ----------------------------------


# Setup ===================================================================


# Packages ----------------------------------------------------------------

## [packages]
required <- c(
  # data acquisition
  "rgbif",            # GBIF taxonomy and occurrence search
  "pestr",            # EPPO Global Database interface
  "DBI",              # database connection closed after the EPPO query
  "CoordinateCleaner",# automated coordinate quality flags
  "tmaptools",        # geocode_OSM(): OSM/Nominatim geocoding
  # spatial
  "terra",            # raster and vector operations
  "geodata",          # WorldClim 2.1 and ESA WorldCover downloads
  "sf",               # simple features, required by blockCV
  "predicts",         # mess(): multivariate environmental similarity surface
  # modelling
  "maxnet",           # Maxent via glmnet
  "ENMeval",          # hyperparameter tuning and cross-validation
  "blockCV",          # spatial block cross-validation
  "usdm",             # vifstep(): variance inflation screening
  "spThin",           # spatial thinning
  "pROC",             # AUC for the permutation importance routine
  # data manipulation and graphics
  "dplyr",
  "ggplot2"
)

missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing packages: ", paste(missing, collapse = ", "),
       "\nInstall them, or restore the pinned environment with renv::restore().",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(maxnet)
})

# terra and dplyr both define extract(); be explicit throughout.
select <- dplyr::select


# Paths, constants and seeds ----------------------------------------------

## [config]
# ---- directories -----------------------------------------------------------
PROJ      <- normalizePath(".", winslash = "/", mustWork = FALSE)
DIR_DATA  <- file.path(PROJ, "data")        # archived / cached inputs
DIR_ENV   <- file.path(PROJ, "env_data")    # large downloaded rasters
DIR_OUT   <- file.path(PROJ, "outputs")     # rasters, models, result tables
DIR_FIG   <- file.path(PROJ, "figures")     # publication figures
for (d in c(DIR_DATA, DIR_ENV, DIR_OUT, DIR_FIG)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

# ---- study species ---------------------------------------------------------
SPECIES      <- "Meloidogyne graminicola"
SPECIES_SLUG <- "M_graminicola"
EPPO_CODE    <- "MELGGC"      # confirmed via eppo_names_tables(); codeid 14038

# ---- analysis constants ----------------------------------------------------
CFG <- list(
  gbif_limit        = 5000,       # occ_search() record cap
  dedup_round_dp    = 2,          # coordinate rounding for near-duplicate removal (about 1 km)
  thin_km           = 10,         # spThin minimum nearest-neighbour distance
  thin_reps         = 100,        # spThin iterations; the most-retained rep is kept
  bioclim_res       = 2.5,        # arc-minutes
  bioclim_vars      = c("wc2.1_2.5m_bio_1",  # BIO1  annual mean temperature
                        "wc2.1_2.5m_bio_4",  # BIO4  temperature seasonality
                        "wc2.1_2.5m_bio_12", # BIO12 annual precipitation
                        "wc2.1_2.5m_bio_15"),# BIO15 precipitation seasonality
  extent_buffer_deg = 10,         # buffer around occurrences defining the study extent
  vif_threshold     = 5,          # usdm::vifstep th
  n_background      = 1000,       # background / pseudo-absence points
  feature_classes   = c("L", "LQ"),          # ENMeval tune.args fc
  reg_multipliers   = c(1, 1.5, 2, 3),       # ENMeval tune.args rm
  kfolds            = 5,
  block_size_m      = 500000,     # 500 km spatial blocks (the published scheme)
  threshold_quantile= 0.10,       # 10th percentile of training-presence suitability
  perm_reps         = 1,          # permutations per variable; 1 reproduces the published run
  weights_in_model  = FALSE       # FALSE reproduces the published run; see @sec-weighting
)

# ---- seeds -----------------------------------------------------------------
# The published run seeded only the two background draws and the two cv_spatial()
# calls, all with 42. Seeds are added here for every other stochastic stage.
SEED <- list(
  clean        = 42,
  vif          = 42,
  thin         = 42,
  background   = 42,   # as published
  enmeval      = 42,
  block        = 42,   # as published
  block_autocor= 42,
  permutation  = 42
)

# ---- re-run switches -------------------------------------------------------
# TRUE forces a stage to hit the network even when a cached file exists.
REFRESH <- list(
  gbif     = FALSE,
  eppo     = FALSE,
  geocode  = FALSE   # leave FALSE: Nominatim is mutable and the extent depends on it
)


# Cache-aware IO helpers --------------------------------------------------

## [cache-helpers]
cache_csv <- function(path, expr, refresh = FALSE) {
  if (file.exists(path) && !refresh) {
    message("[cache] reading ", basename(path))
    return(utils::read.csv(path, stringsAsFactors = FALSE, check.names = TRUE))
  }
  message("[build] computing ", basename(path))
  obj <- force(expr)
  utils::write.csv(obj, path, row.names = FALSE)
  obj
}

cache_rast <- function(path, expr, refresh = FALSE) {
  if (file.exists(path) && !refresh) {
    message("[cache] reading ", basename(path))
    return(terra::rast(path))
  }
  message("[build] computing ", basename(path))
  r <- force(expr)
  terra::writeRaster(r, path, overwrite = TRUE)
  r
}

cache_rds <- function(path, expr, refresh = FALSE) {
  if (file.exists(path) && !refresh) {
    message("[cache] reading ", basename(path))
    return(readRDS(path))
  }
  message("[build] computing ", basename(path))
  obj <- force(expr)
  saveRDS(obj, path)
  obj
}

# Convenience: absolute path inside a project directory.
p_data <- function(...) file.path(DIR_DATA, ...)
p_out  <- function(...) file.path(DIR_OUT,  ...)
p_fig  <- function(...) file.path(DIR_FIG,  ...)


# Occurrence data =========================================================


# GBIF records ------------------------------------------------------------

## [gbif-taxonomy]
if (!file.exists(file.path(DIR_DATA, "gbif_raw.csv")) || REFRESH$gbif) {
  taxon_key <- rgbif::name_backbone(name = SPECIES)
  print(taxon_key[, c("usageKey", "scientificName", "rank", "status", "matchType")])

  # Record count available before committing to a download.
  # NB the published run used the deprecated `georeferenced=` argument; the modern
  # equivalent below returns the same quantity.
  n_available <- rgbif::occ_count(
    taxonKey           = taxon_key$usageKey,
    hasCoordinate      = TRUE,
    hasGeospatialIssue = FALSE
  )
  message("Georeferenced records available on GBIF: ", n_available)
}

## [gbif-download]
gbif_raw <- cache_csv(p_data("gbif_raw.csv"), {
  res <- rgbif::occ_search(
    scientificName     = SPECIES,
    hasCoordinate      = TRUE,
    hasGeospatialIssue = FALSE,
    limit              = CFG$gbif_limit
  )
  keep <- c("scientificName", "decimalLatitude", "decimalLongitude",
            "country", "year", "basisOfRecord", "coordinateUncertaintyInMeters")
  as.data.frame(res$data)[, intersect(keep, names(res$data)), drop = FALSE]
}, refresh = REFRESH$gbif)

message("GBIF raw rows: ", nrow(gbif_raw))

## [gbif-clean]
gbif_clean <- gbif_raw |>
  filter(!is.na(decimalLatitude), !is.na(decimalLongitude)) |>
  distinct(decimalLatitude, decimalLongitude, .keep_all = TRUE)

message("After removing missing and duplicate coordinates: ", nrow(gbif_clean))

# Archive the cleaned set under the name used in the published run.
utils::write.csv(gbif_clean,
                 p_data(paste0(gsub(" ", "_", SPECIES), "_GBIF.csv")),
                 row.names = FALSE)


# Automated coordinate quality screening ----------------------------------

## [coordinate-cleaning]
set.seed(SEED$clean)
cc_flags <- CoordinateCleaner::clean_coordinates(
  x       = gbif_clean,
  lon     = "decimalLongitude",
  lat     = "decimalLatitude",
  species = "scientificName"
)

gbif_final <- gbif_clean[cc_flags$.summary, ]

message(sprintf("CoordinateCleaner retained %d of %d records (%d flagged).",
                nrow(gbif_final), nrow(gbif_clean),
                nrow(gbif_clean) - nrow(gbif_final)))

# Per-test flag counts, for Table 1.
cc_summary <- vapply(
  setdiff(names(cc_flags)[startsWith(names(cc_flags), ".")], ".summary"),
  function(nm) sum(!cc_flags[[nm]]),
  integer(1)
)
print(cc_summary)


# EPPO Global Database records --------------------------------------------

## [eppo-token]
eppo_token_str <- Sys.getenv("EPPO_TOKEN", unset = NA_character_)
have_token     <- !is.na(eppo_token_str) && nzchar(eppo_token_str)

if (!have_token) {
  message("EPPO_TOKEN is not set. This is fine when the cached EPPO tables are ",
          "present; set it in .Renviron only if you need to re-query EPPO.")
}

## [eppo-distribution]
f_eppo_dist <- p_data(paste0(gsub(" ", "_", SPECIES), "_EPPO_distribution.csv"))

eppo_distribution <- cache_csv(f_eppo_dist, {

  sqlite_zip <- p_data("sqlite_all.zip")
  sqlite_db  <- p_data("eppocodes_all.sqlite")

  if (!file.exists(sqlite_db)) {
    if (!file.exists(sqlite_zip)) {
      stop("EPPO SQLite database not found.\n",
           "pestr::eppo_database_download() is non-functional upstream. Log in at\n",
           "  https://data.eppo.int/  ->  Data services  ->  SQLite full database\n",
           "and save the archive to: ", sqlite_zip, call. = FALSE)
    }
    utils::unzip(sqlite_zip, exdir = DIR_DATA)
  }
  stopifnot(file.exists(sqlite_db))

  con <- pestr::eppo_database_connect(filepath = dirname(sqlite_db),
                                      filename = basename(sqlite_db))
  on.exit(try(DBI::dbDisconnect(con), silent = TRUE), add = TRUE)

  names_table <- pestr::eppo_names_tables(SPECIES, con)
  found_code  <- unique(unlist(lapply(names_table, function(x)
    if (is.data.frame(x) && "eppocode" %in% names(x)) x$eppocode)))
  message("EPPO code(s) resolved for ", SPECIES, ": ",
          paste(found_code, collapse = ", "))
  if (length(found_code) && !EPPO_CODE %in% found_code) {
    warning("Expected EPPO code ", EPPO_CODE, " was not returned.")
  }

  distri <- pestr::eppo_tabletools_distri(names_table)$long_table

  # Column names arrive with spaces ("country code"); make them syntactic before
  # anything downstream relies on them, then canonicalise the three the rest of
  # the script reads, since pestr has renamed them between versions.
  names(distri) <- make.names(names(distri))

  canon <- function(df, want, pattern) {
    hit <- names(df)[grepl(pattern, names(df), ignore.case = TRUE)]
    hit <- setdiff(hit, want)
    if (length(hit) && !want %in% names(df)) names(df)[names(df) == hit[1]] <- want
    df
  }
  distri <- canon(distri, "country", "^country(\\.?name)?$")
  distri <- canon(distri, "state",   "^(state|subdiv|subdivision)(\\.?name)?$")
  distri <- canon(distri, "Status",  "^status$")
  if (!"state" %in% names(distri)) distri$state <- NA_character_
  stopifnot(all(c("country", "state", "Status") %in% names(distri)))

  distri
}, refresh = REFRESH$eppo)

message("EPPO distribution rows: ", nrow(eppo_distribution))
print(table(sub(",.*$", "", eppo_distribution$Status)))


# Geocoding administrative units ------------------------------------------

## [eppo-geocode]
f_geo <- p_data(paste0(gsub(" ", "_", SPECIES), "_EPPO_geocoded.csv"))

eppo_geocoded <- cache_csv(f_geo, {
  to_geocode <- eppo_distribution |>
    mutate(geocode_query = ifelse(!is.na(state) & nzchar(state),
                                  paste(state, country, sep = ", "),
                                  country))

  coords <- vapply(to_geocode$geocode_query, function(loc) {
    Sys.sleep(1)   # Nominatim usage policy: at most one request per second
    hit <- tryCatch(tmaptools::geocode_OSM(loc, as.data.frame = TRUE),
                    error = function(e) NULL)
    if (is.null(hit) || !nrow(hit)) return(c(lat = NA_real_, lon = NA_real_))
    c(lat = hit$lat[1], lon = hit$lon[1])
  }, numeric(2))

  to_geocode$eppo_lat <- coords["lat", ]
  to_geocode$eppo_lon <- coords["lon", ]
  to_geocode
}, refresh = REFRESH$geocode)

message(sprintf("Geocoded %d of %d rows (%d failures).",
                sum(!is.na(eppo_geocoded$eppo_lat)), nrow(eppo_geocoded),
                sum(is.na(eppo_geocoded$eppo_lat))))

## [eppo-split]
# "Absent, confirmed by survey" rows are held back entirely from model fitting
# and used only as independent validation.
eppo_presence <- eppo_geocoded |>
  filter(!grepl("^Absent", Status), !is.na(eppo_lat), !is.na(eppo_lon)) |>
  mutate(
    coordinatePrecision = ifelse(!is.na(state) & nzchar(state),
                                 "state-centroid", "country-centroid"),
    basisOfRecord = "EPPO_DISTRIBUTION"
  )

eppo_absence <- eppo_geocoded |>
  filter(grepl("^Absent", Status), !is.na(eppo_lat), !is.na(eppo_lon))

message("EPPO presence: ", nrow(eppo_presence),
        " | EPPO absence: ", nrow(eppo_absence))
print(eppo_absence[, c("country", "Status")])


# Precision weighting -----------------------------------------------------

## [precision-weights]
WEIGHTS <- c(exact = 1.0, `state-centroid` = 0.5, `country-centroid` = 0.25)

gbif_pts <- gbif_final |>
  transmute(
    lon = decimalLongitude,
    lat = decimalLatitude,
    coordinatePrecision = "exact",
    source = "GBIF",
    weight = WEIGHTS[["exact"]]
  )

eppo_pts <- eppo_presence |>
  transmute(
    lon = eppo_lon,
    lat = eppo_lat,
    coordinatePrecision = coordinatePrecision,
    source = "EPPO",
    weight = unname(WEIGHTS[coordinatePrecision])
  )

stopifnot(!any(is.na(eppo_pts$weight)))
print(table(c(gbif_pts$coordinatePrecision, eppo_pts$coordinatePrecision)))


# What these weights did, and did not, do in the published run ------------


# Combining and deduplication ---------------------------------------------

## [combine-dedup]
occ_all <- bind_rows(gbif_pts, eppo_pts)
message("Combined records: ", nrow(occ_all))
print(table(occ_all$source, occ_all$coordinatePrecision))

# Rounding to 2 dp, about 1 km, catches near-duplicates, not only exact ties.
# arrange(desc(weight)) first, so the most precise version of a duplicated
# location is the one retained.
occ_combined <- occ_all |>
  mutate(lon_round = round(lon, CFG$dedup_round_dp),
         lat_round = round(lat, CFG$dedup_round_dp)) |>
  arrange(desc(weight)) |>
  distinct(lon_round, lat_round, .keep_all = TRUE) |>
  select(-lon_round, -lat_round)

message(sprintf("After deduplication: %d (%d removed)",
                nrow(occ_combined), nrow(occ_all) - nrow(occ_combined)))


# Environmental predictors ================================================


# WorldClim bioclimatic variables -----------------------------------------

## [worldclim]
bioclim_all <- geodata::worldclim_global(var = "bio", res = CFG$bioclim_res,
                                         path = DIR_ENV)
bioclim_sel <- bioclim_all[[CFG$bioclim_vars]]
print(bioclim_sel)


# Study extent ------------------------------------------------------------

## [study-extent]
b <- CFG$extent_buffer_deg
study_ext <- terra::ext(
  min(occ_combined$lon) - b, max(occ_combined$lon) + b,
  min(occ_combined$lat) - b, max(occ_combined$lat) + b
)
print(study_ext)

env_clim <- terra::crop(bioclim_sel, study_ext)
print(env_clim)

# The published rasters have this snapped geometry; assert it so that any drift
# in the geocoded inputs is caught loudly rather than propagating silently.
PUBLISHED_GEOM <- list(ncol = 5837, nrow = 2219,
                       xmin = -110.4583333333, xmax = 132.75,
                       ymin =  -39.8333333333, ymax =  52.625)
geom_now <- list(ncol = terra::ncol(env_clim), nrow = terra::nrow(env_clim),
                 xmin = terra::xmin(env_clim), xmax = terra::xmax(env_clim),
                 ymin = terra::ymin(env_clim), ymax = terra::ymax(env_clim))
geom_match <- isTRUE(all.equal(PUBLISHED_GEOM, geom_now, tolerance = 1e-6))
if (!geom_match) {
  warning("Study extent differs from the published run:\n",
          paste(capture.output(str(geom_now)), collapse = "\n"),
          "\nDownstream numbers will not match the manuscript exactly.")
} else {
  message("Study extent matches the published run exactly.")
}


# Cropland fraction -------------------------------------------------------

## [cropland]
cropland_raw <- geodata::landcover(var = "cropland", path = DIR_ENV)

env_full <- cache_rast(p_out(paste0(SPECIES_SLUG, "_predictors.tif")), {
  cl <- terra::crop(cropland_raw, study_ext)
  cl <- terra::resample(cl, env_clim[[1]], method = "bilinear")
  s  <- c(env_clim, cl)
  names(s) <- c(CFG$bioclim_vars, "cropland_fraction")
  s
})
names(env_full) <- c(CFG$bioclim_vars, "cropland_fraction")
print(env_full)


# Multicollinearity screening ---------------------------------------------

## [vif]
set.seed(SEED$vif)
vif_clim <- usdm::vifstep(env_clim, th = CFG$vif_threshold)
print(vif_clim)

set.seed(SEED$vif)
vif_full <- usdm::vifstep(env_full, th = CFG$vif_threshold)
print(vif_full)


# Screening occurrences against predictor coverage ------------------------

## [na-screen]
occ_vect  <- terra::vect(occ_combined, geom = c("lon", "lat"), crs = "EPSG:4326")
env_at_pt <- terra::extract(env_clim, occ_vect)

na_rows <- which(rowSums(is.na(env_at_pt[, -1, drop = FALSE])) > 0)
if (length(na_rows)) {
  message("Records falling outside predictor coverage (dropped):")
  print(occ_combined[na_rows, c("lon", "lat", "source", "coordinatePrecision")])
}

occ_valid <- occ_combined[setdiff(seq_len(nrow(occ_combined)), na_rows), ]
message("Records with complete predictor coverage: ", nrow(occ_valid))

env_valid <- terra::extract(
  env_clim,
  terra::vect(occ_valid, geom = c("lon", "lat"), crs = "EPSG:4326")
)
stopifnot(!anyNA(env_valid[, -1, drop = FALSE]))

occ_with_env <- cbind(occ_valid, env_valid[, -1, drop = FALSE])
utils::write.csv(occ_with_env, p_out(paste0(SPECIES_SLUG, "_occurrences_with_env.csv")),
                 row.names = FALSE)


# Spatial thinning ========================================================

## [thinning]
set.seed(SEED$thin)
thin_out <- spThin::thin(
  loc.data                 = occ_with_env |> mutate(species = SPECIES_SLUG),
  lat.col                  = "lat",
  long.col                 = "lon",
  spec.col                 = "species",
  thin.par                 = CFG$thin_km,
  reps                     = CFG$thin_reps,
  locs.thinned.list.return = TRUE,
  write.files              = FALSE,
  write.log.file           = FALSE
)

retained  <- vapply(thin_out, nrow, integer(1))
best_thin <- thin_out[[which.max(retained)]]
message(sprintf("Thinning: max retained %d of %d across %d reps (%d reps hit the max).",
                max(retained), nrow(occ_with_env), CFG$thin_reps,
                sum(retained == max(retained))))

occ_thinned <- occ_with_env |>
  inner_join(best_thin, by = c("lat" = "Latitude", "lon" = "Longitude"))

# The join is on floating-point coordinates, which is safe only because spThin
# returns the original values unaltered. Verify rather than assume.
if (nrow(occ_thinned) != nrow(best_thin)) {
  stop(sprintf("Thinning join lost rows: %d joined vs %d retained by spThin.",
               nrow(occ_thinned), nrow(best_thin)), call. = FALSE)
}

message("Thinned occurrences: ", nrow(occ_thinned))
print(table(occ_thinned$source, occ_thinned$coordinatePrecision))


# Background points =======================================================

## [background]
set.seed(SEED$background)
bg_points <- terra::spatSample(env_full, size = CFG$n_background,
                               method = "random", na.rm = TRUE, xy = TRUE)
message("Background points drawn: ", nrow(bg_points))

bg_data <- bg_points |>
  rename(lon = x, lat = y) |>
  mutate(presence = 0, weight = 1)   # background is not precision-weighted

## [model-frame]
PRED_VARS <- c(CFG$bioclim_vars, "cropland_fraction")

# Re-extract at the thinned presences from the five-layer stack, so that presence
# and background rows are guaranteed to come from the identical raster. The
# climate-only columns carried over from the coverage screen are replaced, not
# duplicated.
occ_thin_vect <- terra::vect(occ_thinned, geom = c("lon", "lat"), crs = "EPSG:4326")
env_thin      <- terra::extract(env_full, occ_thin_vect)
stopifnot(!anyNA(env_thin[, -1, drop = FALSE]))

occ_thinned <- occ_thinned |>
  select(-any_of(PRED_VARS)) |>
  bind_cols(as.data.frame(env_thin)[, -1, drop = FALSE])
stopifnot(all(PRED_VARS %in% names(occ_thinned)))

pres_data <- occ_thinned |>
  select(lon, lat, weight, all_of(PRED_VARS)) |>
  mutate(presence = 1)

model_data <- bind_rows(pres_data, bg_data)
print(table(model_data$presence))

predictors <- model_data[, PRED_VARS, drop = FALSE]
response   <- model_data$presence


# Hyperparameter tuning ===================================================

## [enmeval-tuning]
# ENMevaluate wrapper. Uses the predictor raster (as in the submitted run) and,
# if that runs out of memory (std::bad_alloc on machines with < 8 GB RAM),
# repeats the same evaluation in "samples with data" mode: predictor values at
# the occurrence and background points are passed directly and no raster is
# processed. Folds, metrics and model selection are unaffected; AICc is then
# computed over the background instead of the full raster.
enm_eval <- function(occ_df, mode = c("auto", "swd"), ...) {
  mode <- match.arg(mode)
  occs_xy <- occ_df[, c("lon", "lat")]
  bg_xy   <- stats::setNames(bg_points[, c("x", "y")], c("lon", "lat"))
  swd <- function() {
    message("  [ENMeval] samples-with-data mode (no raster in memory)")
    ENMeval::ENMevaluate(
      occs = cbind(occs_xy, as.data.frame(occ_df)[, PRED_VARS]),
      bg   = cbind(bg_xy, as.data.frame(bg_points)[, PRED_VARS]),
      envs = NULL, algorithm = "maxnet", ...)
  }
  if (mode == "swd") return(swd())
  terra::terraOptions(memfrac = 0.5)
  gc()
  tryCatch(
    ENMeval::ENMevaluate(occs = occs_xy, envs = env_full, bg = bg_xy,
                         algorithm = "maxnet", ...),
    error = function(e) {
      message("  [ENMeval] raster mode failed (", conditionMessage(e), "); retrying")
      gc(); swd()
    })
}

set.seed(SEED$enmeval)
eval_random <- cache_rds(p_out(paste0(SPECIES_SLUG, "_ENMeval_random.rds")), {
  enm_eval(occ_thinned,
    partitions = "randomkfold",
    partition.settings = list(kfolds = CFG$kfolds),
    tune.args = list(fc = CFG$feature_classes, rm = CFG$reg_multipliers))
})

res_random <- eval_random@results
print(res_random[, c("fc", "rm", "auc.train", "cbi.train", "auc.val.avg",
                     "auc.val.sd", "cbi.val.avg", "or.10p.avg", "AICc", "delta.AICc",
                     "ncoef")])

## [enmeval-select]
best_random <- res_random |> filter(delta.AICc == 0)
print(best_random)

BEST_FC <- as.character(best_random$fc[1])
BEST_RM <- as.numeric(as.character(best_random$rm[1]))
message(sprintf("Selected by AICc: fc = %s, rm = %s (AICc = %.3f, ncoef = %d)",
                BEST_FC, BEST_RM, best_random$AICc[1], best_random$ncoef[1]))


# Final model =============================================================

## [final-fit]
mx_final <- cache_rds(p_out(paste0(SPECIES_SLUG, "_FINAL_model.rds")), {
  f <- maxnet::maxnet.formula(response, predictors,
                              classes = tolower(BEST_FC))
  maxnet::maxnet(p = response, data = predictors, f = f, regmult = BEST_RM)
})

print(mx_final)
message("Non-zero coefficients: ", sum(mx_final$betas != 0))


# Optional refit with precision weights applied ---------------------------

## [weighted-refit]  (not run: kept for reference, as in the .qmd)
# # Capability check first: maxnet's signature has changed across versions, so test
# # for the argument rather than assuming it exists.
# if (!"weights" %in% names(formals(maxnet::maxnet))) {
#   message("This maxnet build does not accept `weights`. ",
#           "Weighting can instead be applied by replicating rows in proportion ",
#           "to weight, or by moving to a glmnet fit called directly.")
# } else {
#   f_w <- maxnet::maxnet.formula(response, predictors, classes = tolower(BEST_FC))
#   mx_weighted <- maxnet::maxnet(p = response, data = predictors, f = f_w,
#                                 regmult = BEST_RM,
#                                 weights = model_data$weight)
#   saveRDS(mx_weighted, p_out(paste0(SPECIES_SLUG, "_WEIGHTED_model.rds")))
#
#   pred_unw <- predict(mx_final,    predictors, type = "cloglog", clamp = TRUE)
#   pred_w   <- predict(mx_weighted, predictors, type = "cloglog", clamp = TRUE)
#
#   message(sprintf("Correlation between weighted and unweighted predictions: %.4f",
#                   cor(as.vector(pred_unw), as.vector(pred_w))))
#   message(sprintf("Mean absolute difference in suitability: %.4f",
#                   mean(abs(as.vector(pred_unw) - as.vector(pred_w)))))
# }


# Prediction and thresholding =============================================

## [predict-raster]
suitability <- cache_rast(p_out(paste0(SPECIES_SLUG, "_FINAL_suitability.tif")), {
  terra::predict(env_full, mx_final, type = "cloglog", na.rm = TRUE)
})
names(suitability) <- "suitability"
print(suitability)
print(terra::global(suitability, c("min", "max", "mean"), na.rm = TRUE))

## [threshold]
train_suit <- terra::extract(
  suitability,
  terra::vect(occ_thinned, geom = c("lon", "lat"), crs = "EPSG:4326")
)[, 2]

threshold_10p <- stats::quantile(train_suit, probs = CFG$threshold_quantile,
                                 na.rm = TRUE)
message(sprintf("10th-percentile training-presence threshold (5-predictor model): %.6f",
                threshold_10p))

binary_suit <- suitability >= threshold_10p
names(binary_suit) <- "suitable"
terra::writeRaster(binary_suit, p_out(paste0(SPECIES_SLUG, "_FINAL_binary.tif")),
                   overwrite = TRUE)

pct_suitable <- 100 * unlist(terra::global(binary_suit, "mean", na.rm = TRUE))
message(sprintf("Study extent classified suitable: %.1f%%", pct_suitable))


# Independent validation against EPPO status ==============================

## [validate-presence]
eppo_pres_valid <- occ_valid |> filter(source == "EPPO")

eppo_suit <- terra::extract(
  suitability,
  terra::vect(eppo_pres_valid, geom = c("lon", "lat"), crs = "EPSG:4326")
)[, 2]

n_eppo_cov <- sum(!is.na(eppo_suit))
n_hit      <- sum(eppo_suit >= threshold_10p, na.rm = TRUE)

message(sprintf("EPPO presence records above threshold: %d of %d (%.1f%%)",
                n_hit, n_eppo_cov, 100 * n_hit / n_eppo_cov))

eppo_missed <- eppo_pres_valid[which(eppo_suit < threshold_10p), ]
eppo_missed$suitability <- eppo_suit[which(eppo_suit < threshold_10p)]
print(eppo_missed[, c("lon", "lat", "coordinatePrecision", "suitability")])

## [validate-absence]
slovakia_suit <- terra::extract(
  suitability,
  terra::vect(eppo_absence, geom = c("eppo_lon", "eppo_lat"), crs = "EPSG:4326")
)[, 2]

message(sprintf("Slovakia (absence, confirmed by survey): suitability = %.7f, threshold = %.6f -> %s",
                slovakia_suit, threshold_10p,
                ifelse(slovakia_suit < threshold_10p,
                       "correctly below threshold", "INCORRECTLY above threshold")))


# Variable importance =====================================================


# Permutation importance --------------------------------------------------

## [permutation-importance]
calc_permutation_importance <- function(model, data, response, reps = 1, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  base_pred <- as.vector(predict(model, data, type = "cloglog", clamp = FALSE))
  base_auc  <- as.numeric(pROC::auc(pROC::roc(response, base_pred, quiet = TRUE)))

  drops <- vapply(names(data), function(v) {
    per_rep <- vapply(seq_len(reps), function(i) {
      d <- data
      d[[v]] <- sample(d[[v]])
      pr <- as.vector(predict(model, d, type = "cloglog", clamp = FALSE))
      base_auc - as.numeric(pROC::auc(pROC::roc(response, pr, quiet = TRUE)))
    }, numeric(1))
    mean(per_rep)
  }, numeric(1))

  drops <- pmax(drops, 0)
  data.frame(
    variable       = names(data),
    auc_drop       = as.numeric(drops),
    importance_pct = 100 * as.numeric(drops) / sum(drops),
    row.names      = NULL
  ) |> arrange(desc(importance_pct))
}

perm_imp <- calc_permutation_importance(
  mx_final, predictors, response,
  reps = CFG$perm_reps, seed = SEED$permutation
)
print(perm_imp)

utils::write.csv(perm_imp, p_out(paste0(SPECIES_SLUG, "_variable_importance.csv")),
                 row.names = FALSE)


# Jackknife on regularised training gain ----------------------------------

## [jackknife]
calc_gain <- function(model, data, response) {
  preds <- predict(model, data, type = "cloglog", clamp = FALSE)
  preds <- pmin(pmax(as.vector(preds), 1e-10), 1 - 1e-10)
  mean(log(preds[response == 1]))
}

full_gain <- calc_gain(mx_final, predictors, response)
message(sprintf("Full-model gain: %.7f", full_gain))

jack <- lapply(PRED_VARS, function(v) {
  # without-one model
  d_wo <- predictors[, setdiff(PRED_VARS, v), drop = FALSE]
  m_wo <- tryCatch(
    maxnet::maxnet(p = response, data = d_wo,
                   f = maxnet::maxnet.formula(response, d_wo,
                                              classes = tolower(BEST_FC)),
                   regmult = BEST_RM),
    error = function(e) NULL)
  g_wo <- if (is.null(m_wo)) NA_real_ else calc_gain(m_wo, d_wo, response)

  # with-only-one model
  d_only <- predictors[, v, drop = FALSE]
  m_only <- tryCatch(
    maxnet::maxnet(p = response, data = d_only,
                   f = maxnet::maxnet.formula(response, d_only,
                                              classes = tolower(BEST_FC)),
                   regmult = BEST_RM),
    error = function(e) NULL)
  g_only <- if (is.null(m_only)) NA_real_ else calc_gain(m_only, d_only, response)

  data.frame(variable = v, gain_without = g_wo, gain_only = g_only,
             delta_gain = full_gain - g_wo)
}) |> bind_rows() |> arrange(desc(delta_gain))

print(jack)
utils::write.csv(jack, p_out(paste0(SPECIES_SLUG, "_jackknife_results.csv")),
                 row.names = FALSE)


# Extrapolation assessment ================================================

## [mess]
ref_env <- as.data.frame(env_thin)[, -1, drop = FALSE]

mess_r <- cache_rast(p_out(paste0(SPECIES_SLUG, "_MESS.tif")), {
  predicts::mess(env_full, ref_env)
})
names(mess_r) <- "mess"
print(mess_r)

novel_pct <- 100 * mean(terra::values(mess_r) < 0, na.rm = TRUE)
message(sprintf("Study extent requiring extrapolation (MESS < 0): %.1f%%", novel_pct))

## [mess-at-records]
# MESS at the New World and African records, which are the ones the manuscript
# discusses as surveillance priorities.
occ_thinned$mess <- terra::extract(
  mess_r, terra::vect(occ_thinned, geom = c("lon", "lat"), crs = "EPSG:4326")
)[, 2]

amer_afr <- occ_thinned |>
  filter(lon < -25 | (lon > -20 & lon < 52 & lat < 20))

message(sprintf("Americas/Africa records: n = %d, MESS range %.6f to %.6f, all positive: %s",
                nrow(amer_afr), min(amer_afr$mess), max(amer_afr$mess),
                all(amer_afr$mess > 0)))
print(amer_afr[order(amer_afr$mess), c("lon", "lat", "source", "mess")])


# Spatially blocked cross-validation ======================================

## [sac-range]
# blockCV builds its blocks with sf. With the s2 spherical engine on (the sf
# default since 1.0), blocks that span the whole longitude range of this extent
# fail with "Loop 0 is not valid: Edge 1 crosses edge 3". Planar geometry is
# what blockCV assumed when the submitted run was made.
sf::sf_use_s2(FALSE)
occ_sf <- sf::st_as_sf(occ_thinned, coords = c("lon", "lat"), crs = 4326)
bg_sf  <- sf::st_as_sf(bg_points,  coords = c("x", "y"),     crs = 4326)

set.seed(SEED$block_autocor)
sac <- blockCV::cv_spatial_autocor(r = env_full, num_sample = 5000, progress = FALSE)
message(sprintf("Median spatial autocorrelation range across predictors: %.0f m (%.0f km)",
                sac$range, sac$range / 1000))

## [blockcv-autocorr-range]
set.seed(SEED$block)
blocks_sac <- blockCV::cv_spatial(
  x         = occ_sf,
  column    = NULL,
  r         = env_full,
  size      = sac$range,
  k         = CFG$kfolds,
  selection = "random",
  iteration = 100,
  progress  = FALSE,
  report    = TRUE,
  plot      = FALSE
)
print(blocks_sac$records)

## [blockcv-500km]
set.seed(SEED$block)
blocks_500 <- blockCV::cv_spatial(
  x         = occ_sf,
  column    = NULL,
  r         = env_full,
  size      = CFG$block_size_m,
  k         = CFG$kfolds,
  selection = "systematic",
  progress  = FALSE,
  report    = TRUE,
  plot      = FALSE
)
print(blocks_500$records)

set.seed(SEED$block)
bg_blocks_500 <- blockCV::cv_spatial(
  x         = bg_sf,
  column    = NULL,
  r         = env_full,
  size      = CFG$block_size_m,
  k         = CFG$kfolds,
  selection = "systematic",
  progress  = FALSE,
  report    = FALSE,
  plot      = FALSE
)

## [enmeval-blocked]
set.seed(SEED$enmeval)
eval_blocked <- cache_rds(p_out(paste0(SPECIES_SLUG, "_ENMeval_blocked500km.rds")), {
  # Only validation metrics are used from this run, so samples-with-data mode
  # is used directly (it avoids the memory error seen with the full raster).
  enm_eval(occ_thinned, mode = "swd",
    partitions = "user",
    user.grp   = list(occs.grp = blocks_500$folds_ids,
                      bg.grp   = bg_blocks_500$folds_ids),
    tune.args  = list(fc = CFG$feature_classes, rm = CFG$reg_multipliers))
})

res_blocked <- eval_blocked@results
best_blocked <- res_blocked |> filter(fc == BEST_FC, rm == BEST_RM)
print(best_blocked[, c("fc", "rm", "auc.val.avg", "auc.val.sd", "cbi.val.avg",
                       "auc.diff.avg", "or.10p.avg")])


# Sensitivity to occurrence data composition ==============================

## [sensitivity-fit]
occ_gbif_only <- occ_thinned |> filter(source == "GBIF")
message("GBIF-only occurrences after thinning: ", nrow(occ_gbif_only))

set.seed(SEED$enmeval)
eval_gbif <- cache_rds(p_out(paste0(SPECIES_SLUG, "_ENMeval_GBIFonly.rds")), {
  enm_eval(occ_gbif_only,
    partitions = "randomkfold",
    partition.settings = list(kfolds = CFG$kfolds),
    tune.args = list(fc = CFG$feature_classes, rm = CFG$reg_multipliers))
})

res_gbif  <- eval_gbif@results
best_gbif <- res_gbif |> filter(delta.AICc == 0)
print(best_gbif[, c("fc", "rm", "auc.train", "auc.val.avg", "auc.val.sd",
                    "cbi.val.avg", "cbi.val.sd", "auc.diff.avg", "or.10p.avg", "AICc")])

## [sensitivity-predict]
pres_gbif <- occ_gbif_only |> select(lon, lat, all_of(PRED_VARS)) |> mutate(presence = 1)
md_gbif   <- bind_rows(pres_gbif, bg_data |> select(lon, lat, all_of(PRED_VARS), presence))

pred_gbif <- md_gbif[, PRED_VARS, drop = FALSE]
resp_gbif <- md_gbif$presence

mx_gbif <- cache_rds(p_out(paste0(SPECIES_SLUG, "_GBIFonly_model.rds")), {
  maxnet::maxnet(p = resp_gbif, data = pred_gbif,
                 f = maxnet::maxnet.formula(resp_gbif, pred_gbif,
                                            classes = tolower(as.character(best_gbif$fc[1]))),
                 regmult = as.numeric(as.character(best_gbif$rm[1])))
})

suit_gbif <- cache_rast(p_out(paste0(SPECIES_SLUG, "_GBIFonly_suitability.tif")), {
  terra::predict(env_full, mx_gbif, type = "cloglog", na.rm = TRUE)
})
names(suit_gbif) <- "suitability"

## [sensitivity-validate]
train_suit_gbif <- terra::extract(
  suit_gbif, terra::vect(occ_gbif_only, geom = c("lon", "lat"), crs = "EPSG:4326")
)[, 2]
thr_gbif <- stats::quantile(train_suit_gbif, probs = CFG$threshold_quantile, na.rm = TRUE)

eppo_suit_gbif <- terra::extract(
  suit_gbif, terra::vect(eppo_pres_valid, geom = c("lon", "lat"), crs = "EPSG:4326")
)[, 2]
hit_gbif <- 100 * mean(eppo_suit_gbif >= thr_gbif, na.rm = TRUE)

slovakia_gbif <- terra::extract(
  suit_gbif, terra::vect(eppo_absence, geom = c("eppo_lon", "eppo_lat"), crs = "EPSG:4326")
)[, 2]

message(sprintf("GBIF-only model: threshold %.6f | EPPO hit rate %.1f%% | Slovakia %.7f",
                thr_gbif, hit_gbif, slovakia_gbif))

## [sensitivity-table]
sensitivity <- data.frame(
  Metric = c("Occurrence records (thinned)",
             "Training AUC",
             "Validation AUC (mean)",
             "Validation AUC (SD)",
             "Validation CBI (mean)",
             "Validation CBI (SD)",
             "AUC difference (train - val)",
             "Omission rate (10th pct)",
             "AICc",
             "EPPO presence hit rate (%)",
             "Slovakia suitability"),
  Integrated = c(nrow(occ_thinned),
                 best_random$auc.train[1], best_random$auc.val.avg[1],
                 best_random$auc.val.sd[1], best_random$cbi.val.avg[1],
                 best_random$cbi.val.sd[1], best_random$auc.diff.avg[1],
                 best_random$or.10p.avg[1], best_random$AICc[1],
                 100 * n_hit / n_eppo_cov, slovakia_suit),
  GBIF_only  = c(nrow(occ_gbif_only),
                 best_gbif$auc.train[1], best_gbif$auc.val.avg[1],
                 best_gbif$auc.val.sd[1], best_gbif$cbi.val.avg[1],
                 best_gbif$cbi.val.sd[1], best_gbif$auc.diff.avg[1],
                 best_gbif$or.10p.avg[1], best_gbif$AICc[1],
                 hit_gbif, slovakia_gbif),
  check.names = FALSE
)
print(sensitivity, digits = 4)

utils::write.csv(sensitivity, p_out(paste0(SPECIES_SLUG, "_sensitivity_comparison.csv")),
                 row.names = FALSE)


# Tables ==================================================================

## [table-1]
table1 <- data.frame(
  Stage = c("GBIF records returned",
            "After removing missing/duplicate coordinates",
            "Flagged by CoordinateCleaner",
            "GBIF records retained",
            "EPPO distribution rows",
            "EPPO presence rows geocoded",
            "EPPO absence rows (independent test)",
            "Combined presence records",
            "Removed by coordinate deduplication",
            "Removed for missing predictor coverage",
            "Records entering thinning",
            "Records retained after 10 km thinning"),
  n = c(nrow(gbif_raw),
        nrow(gbif_clean),
        nrow(gbif_clean) - nrow(gbif_final),
        nrow(gbif_final),
        nrow(eppo_distribution),
        nrow(eppo_presence),
        nrow(eppo_absence),
        nrow(occ_all),
        nrow(occ_all) - nrow(occ_combined),
        nrow(occ_combined) - nrow(occ_valid),
        nrow(occ_with_env),
        nrow(occ_thinned)),
  check.names = FALSE
)
utils::write.csv(table1, p_out("Table_1_data_assembly.csv"), row.names = FALSE)
knitr::kable(table1, caption = "Table 1. Occurrence data assembly and filtering.")

## [table-2]
table2 <- occ_thinned |>
  count(source, coordinatePrecision, name = "n") |>
  mutate(weight = unname(WEIGHTS[coordinatePrecision]),
         percent = round(100 * n / sum(n), 1)) |>
  arrange(desc(weight))
utils::write.csv(table2, p_out("Table_2_record_composition.csv"), row.names = FALSE)
knitr::kable(table2, caption = "Table 2. Composition of the modelling data set by source and locational precision.")

## [table-3]
vif_tab <- vif_full@results
names(vif_tab) <- c("Variable", "VIF")
table3 <- vif_tab |>
  mutate(Description = c(
    "wc2.1_2.5m_bio_1"   = "BIO1: annual mean temperature",
    "wc2.1_2.5m_bio_4"   = "BIO4: temperature seasonality",
    "wc2.1_2.5m_bio_12"  = "BIO12: annual precipitation",
    "wc2.1_2.5m_bio_15"  = "BIO15: precipitation seasonality",
    "cropland_fraction"  = "ESA WorldCover cropland fraction"
  )[Variable]) |>
  select(Variable, Description, VIF)
utils::write.csv(table3, p_out("Table_3_predictors_vif.csv"), row.names = FALSE)
knitr::kable(table3, digits = 3,
             caption = "Table 3. Environmental predictors and variance inflation factors.")

## [table-4]
table4 <- res_random |>
  select(fc, rm, auc.train, cbi.train, auc.val.avg, auc.val.sd,
         cbi.val.avg, auc.diff.avg, or.10p.avg, AICc, delta.AICc, ncoef) |>
  arrange(AICc)
utils::write.csv(table4, p_out("Table_4_enmeval_tuning.csv"), row.names = FALSE)
knitr::kable(table4, digits = 4,
             caption = "Table 4. Model tuning and cross-validation results under random 5-fold partitioning.")

## [table-5]
table5 <- perm_imp |>
  rename(Variable = variable,
         `AUC drop` = auc_drop,
         `Permutation importance (%)` = importance_pct) |>
  left_join(jack |> select(variable, delta_gain, gain_without),
            by = c("Variable" = "variable")) |>
  rename(`Delta gain (jackknife)` = delta_gain,
         `Gain without variable` = gain_without)
utils::write.csv(table5, p_out("Table_5_variable_importance.csv"), row.names = FALSE)
knitr::kable(table5, digits = 4,
             caption = "Table 5. Permutation importance and jackknife contribution of each predictor.")

## [table-6]
table6 <- sensitivity
utils::write.csv(table6, p_out("Table_6_sensitivity.csv"), row.names = FALSE)
knitr::kable(table6, digits = 4,
             caption = "Table 6. Sensitivity of model performance and independent validation to occurrence data composition.")

## [table-s1]
# Supplementary: the full occurrence table used for fitting.
table_s1 <- occ_thinned |>
  transmute(species = SPECIES,
            longitude = round(lon, 5), latitude = round(lat, 5),
            source, coordinatePrecision, weight,
            bio_1 = .data[[CFG$bioclim_vars[1]]],
            bio_4 = .data[[CFG$bioclim_vars[2]]],
            bio_12 = .data[[CFG$bioclim_vars[3]]],
            bio_15 = .data[[CFG$bioclim_vars[4]]],
            cropland_fraction, mess) |>
  arrange(desc(weight), longitude)
utils::write.csv(table_s1, p_out(paste0(SPECIES_SLUG, "_occurrence_points_TABLE.csv")),
                 row.names = FALSE)
message("Supplementary occurrence table rows: ", nrow(table_s1))


# Figures =================================================================

## [fig-theme]
var_labels <- c(
  "wc2.1_2.5m_bio_1"  = "BIO1 annual mean temperature",
  "wc2.1_2.5m_bio_4"  = "BIO4 temperature seasonality",
  "wc2.1_2.5m_bio_12" = "BIO12 annual precipitation",
  "wc2.1_2.5m_bio_15" = "BIO15 precipitation seasonality",
  "cropland_fraction" = "Cropland fraction"
)

theme_sdm <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
      plot.title    = element_text(face = "bold", size = base + 2),
      plot.subtitle = element_text(colour = "grey30", size = base - 1),
      legend.position = "bottom",
      legend.key.width = grid::unit(2.2, "cm")
    )
}

world <- geodata::world(path = DIR_ENV)
world_crop <- terra::crop(world, study_ext)
world_df   <- suppressWarnings(as.data.frame(terra::geom(world_crop)))

save_fig <- function(plot, file, width = 10, height = 6) {
  ggsave(p_fig(file), plot = plot, width = width, height = height,
         dpi = 300, bg = "white")
  message("[figure] ", file)
  invisible(p_fig(file))
}

# Convert a SpatRaster to a plotting data frame, downsampled for speed.
rast_df <- function(r, max_cells = 6e5) {
  rr <- if (terra::ncell(r) > max_cells) {
    terra::aggregate(r, fact = ceiling(sqrt(terra::ncell(r) / max_cells)),
                     fun = "mean", na.rm = TRUE)
  } else r
  d <- terra::as.data.frame(rr, xy = TRUE, na.rm = TRUE)
  names(d) <- c("x", "y", "value")
  d
}

## [figure-1]
fig1 <- ggplot() +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = "grey95", colour = "grey70", linewidth = 0.2) +
  geom_point(data = occ_thinned,
             aes(lon, lat, colour = source, size = weight, shape = coordinatePrecision),
             alpha = 0.85, stroke = 0.6) +
  scale_colour_manual(values = c(GBIF = "#1B7837", EPPO = "#762A83"), name = "Source") +
  scale_size_continuous(range = c(1.4, 3.4), name = "Locational weight") +
  scale_shape_manual(values = c(exact = 16, `state-centroid` = 17,
                                `country-centroid` = 15),
                     name = "Precision") +
  coord_equal(xlim = c(terra::xmin(study_ext), terra::xmax(study_ext)),
              ylim = c(terra::ymin(study_ext), terra::ymax(study_ext)),
              expand = FALSE) +
  labs(title = expression(paste("Occurrence records of ", italic("Meloidogyne graminicola"))),
       subtitle = sprintf("n = %d after 10 km spatial thinning; %d GBIF, %d EPPO",
                          nrow(occ_thinned), sum(occ_thinned$source == "GBIF"),
                          sum(occ_thinned$source == "EPPO")),
       x = "Longitude", y = "Latitude") +
  theme_sdm()

save_fig(fig1, "Figure_1_occurrence_records.png")
fig1

## [figure-2]
suit_df <- rast_df(suitability)

fig2 <- ggplot() +
  geom_raster(data = suit_df, aes(x, y, fill = value)) +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = NA, colour = "grey25", linewidth = 0.18) +
  geom_point(data = occ_thinned, aes(lon, lat), colour = "red",
             size = 0.7, alpha = 0.8) +
  scale_fill_viridis_c(option = "D", name = "Suitability (cloglog)",
                       limits = c(0, 1)) +
  coord_equal(xlim = c(terra::xmin(study_ext), terra::xmax(study_ext)),
              ylim = c(terra::ymin(study_ext), terra::ymax(study_ext)),
              expand = FALSE) +
  labs(title = "Predicted environmental suitability",
       subtitle = sprintf("Maxent (%s, rm = %s) | validation AUC = %.3f | threshold = %.3f | occurrences in red",
                          BEST_FC, BEST_RM, best_random$auc.val.avg[1], threshold_10p),
       x = "Longitude", y = "Latitude") +
  theme_sdm()

save_fig(fig2, "Figure_2_suitability_map.png")
fig2

## [figure-3]
imp_long <- perm_imp |>
  select(variable, value = importance_pct) |>
  mutate(metric = "Permutation importance (%)") |>
  bind_rows(jack |> select(variable, value = delta_gain) |>
              mutate(metric = "Jackknife contribution (delta gain)"))

imp_long$label <- unname(var_labels[imp_long$variable])

fig3 <- ggplot(imp_long,
               aes(x = reorder(label, value), y = value,
                   fill = value > 0)) +
  geom_col(width = 0.68) +
  geom_hline(yintercept = 0, linewidth = 0.3, colour = "grey40") +
  coord_flip() +
  facet_wrap(~ metric, scales = "free_x") +
  scale_fill_manual(values = c(`TRUE` = "#2166AC", `FALSE` = "#B2182B"),
                    guide = "none") +
  labs(title = "Predictor contributions under two measures",
       subtitle = "Permutation importance reflects reliance as fitted; jackknife delta gain reflects unique information",
       x = NULL, y = NULL) +
  theme_sdm()

save_fig(fig3, "Figure_3_variable_importance.png", width = 11, height = 5.5)
fig3

## [figure-4]
# Marginal response: vary one predictor across its observed range, holding the
# others at their mean across the modelling data.
response_curve <- function(model, data, var, n = 200) {
  rng  <- range(data[[var]], na.rm = TRUE)
  grid <- data[rep(1, n), , drop = FALSE]
  for (v in names(data)) grid[[v]] <- mean(data[[v]], na.rm = TRUE)
  grid[[var]] <- seq(rng[1], rng[2], length.out = n)
  data.frame(variable = var,
             x = grid[[var]],
             suitability = as.vector(predict(model, grid, type = "cloglog",
                                             clamp = TRUE)))
}

curves <- bind_rows(lapply(PRED_VARS, function(v)
  response_curve(mx_final, predictors, v)))
curves$label <- unname(var_labels[curves$variable])

fig4 <- ggplot(curves, aes(x, suitability)) +
  geom_line(linewidth = 0.9, colour = "#1B7837") +
  facet_wrap(~ label, scales = "free_x", ncol = 3) +
  labs(title = "Marginal response curves",
       subtitle = "Each predictor varied across its observed range with the others held at their mean",
       x = NULL, y = "Suitability (cloglog)") +
  theme_sdm()

save_fig(fig4, "Figure_4_response_curves.png", width = 11, height = 6.5)

# Single-panel version of the cropland curve, as used in the manuscript.
fig4b <- ggplot(curves |> filter(variable == "cropland_fraction"),
                aes(x, suitability)) +
  geom_line(linewidth = 1.1, colour = "#1B7837") +
  scale_y_continuous(limits = c(0, 1)) +
  labs(title = "Response to cropland fraction",
       subtitle = "Non-zero suitability at zero cropland, saturating above roughly 60 to 70 per cent cover",
       x = "Cropland fraction", y = "Suitability (cloglog)") +
  theme_sdm()

save_fig(fig4b, "Figure_4_response_curve_cropland.png", width = 8, height = 5.5)
fig4

## [figure-5]
mess_df <- rast_df(mess_r)
mess_df$class <- ifelse(mess_df$value < 0, "Novel (extrapolation)", "Analogous")

fig5 <- ggplot() +
  geom_raster(data = mess_df, aes(x, y, fill = pmax(pmin(value, 100), -100))) +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = NA, colour = "grey25", linewidth = 0.18) +
  geom_point(data = occ_thinned, aes(lon, lat), colour = "black",
             size = 0.6, alpha = 0.7) +
  scale_fill_gradient2(low = "#B2182B", mid = "grey96", high = "#2166AC",
                       midpoint = 0, name = "MESS") +
  coord_equal(xlim = c(terra::xmin(study_ext), terra::xmax(study_ext)),
              ylim = c(terra::ymin(study_ext), terra::ymax(study_ext)),
              expand = FALSE) +
  labs(title = "Multivariate environmental similarity surface",
       subtitle = sprintf("Negative values indicate extrapolation beyond the training envelope; %.1f%% of the extent is novel",
                          novel_pct),
       x = "Longitude", y = "Latitude") +
  theme_sdm()

save_fig(fig5, "Figure_5_MESS.png")
fig5

## [figure-6]
comp <- data.frame(
  Metric = rep(c("Validation AUC", "Validation CBI",
                 "AUC difference", "EPPO hit rate (proportion)"), each = 2),
  Model  = rep(c("Integrated", "GBIF only"), times = 4),
  Value  = c(best_random$auc.val.avg[1], best_gbif$auc.val.avg[1],
             best_random$cbi.val.avg[1], best_gbif$cbi.val.avg[1],
             best_random$auc.diff.avg[1], best_gbif$auc.diff.avg[1],
             n_hit / n_eppo_cov,          hit_gbif / 100),
  SD     = c(best_random$auc.val.sd[1],  best_gbif$auc.val.sd[1],
             best_random$cbi.val.sd[1],  best_gbif$cbi.val.sd[1],
             NA, NA, NA, NA)
)

fig6 <- ggplot(comp, aes(Model, Value, fill = Model)) +
  geom_col(width = 0.62) +
  geom_errorbar(aes(ymin = pmax(Value - SD, 0), ymax = Value + SD),
                width = 0.16, na.rm = TRUE, linewidth = 0.4) +
  facet_wrap(~Metric, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c(Integrated = "#762A83", `GBIF only` = "#1B7837"),
                    guide = "none") +
  labs(title = "Consequences of occurrence data composition",
       subtitle = "Error bars are cross-validation standard deviations; the GBIF-only model scores higher on AUC yet misses far more independent records",
       x = NULL, y = NULL) +
  theme_sdm()

save_fig(fig6, "Figure_6_sensitivity_comparison.png", width = 11, height = 5)
fig6


# Verification against the published run ==================================

## [verify-occurrences]
f_frozen <- p_data("M_graminicola_occurrence_points_FROZEN.csv")

if (!file.exists(f_frozen)) {
  message("Frozen occurrence set not present; skipping verification.")
} else {
  frozen <- utils::read.csv(f_frozen, stringsAsFactors = FALSE)

  message(sprintf("Frozen set: %d records | this run: %d records",
                  nrow(frozen), nrow(occ_thinned)))

  # Match on coordinates rounded to 4 dp, about 11 m, which is well inside a 2.5
  # arc-minute cell and tolerant of CSV round-tripping.
  key <- function(lon, lat) paste(round(lon, 4), round(lat, 4))
  k_frozen <- key(frozen$Longitude, frozen$Latitude)
  k_now    <- key(occ_thinned$lon,  occ_thinned$lat)

  message(sprintf("Shared records: %d | only in frozen set: %d | only in this run: %d",
                  length(intersect(k_frozen, k_now)),
                  length(setdiff(k_frozen, k_now)),
                  length(setdiff(k_now, k_frozen))))

  if (!setequal(k_frozen, k_now)) {
    message("Point sets differ. This is expected if spatial thinning selected a ",
            "different replicate, or if GBIF has added records since the ",
            "reported run. Downstream metrics will differ correspondingly.")
  } else {
    message("Point sets are identical to the published run.")
  }

  # Composition should hold even when the specific points differ.
  print(table(frozen$Data_Source, frozen$Coordinate_Precision))
}

## [verify-reference-results]
compare_to_reference <- function(file, join_by, this_run, cols, tol = 5e-3) {
  path <- file.path(PROJ, "reference_results", file)
  if (!file.exists(path)) {
    message("Reference file absent, skipping: ", file)
    return(invisible(NULL))
  }
  ref <- utils::read.csv(path, stringsAsFactors = FALSE)
  cmp <- merge(ref, this_run, by = join_by, suffixes = c("_published", "_now"))
  for (cc in cols) {
    a <- cmp[[paste0(cc, "_published")]]
    bb <- cmp[[paste0(cc, "_now")]]
    if (is.null(a) || is.null(bb)) next
    worst <- max(abs(a - bb), na.rm = TRUE)
    verdict <- if (worst <= tol) "agrees" else "DIFFERS"
    message(sprintf("  %-14s largest absolute difference %.6g  [%s]", cc, worst, verdict))
  }
  invisible(cmp)
}

message("Jackknife delta gain, published against this run:")
compare_to_reference(
  "M_graminicola_jackknife_results.csv",
  join_by  = "variable",
  this_run = jack |> select(variable, gain_without, delta_gain),
  cols     = c("gain_without")
)

# The published jackknife file names the delta column `gain_drop`.
ref_jack <- file.path(PROJ, "reference_results", "M_graminicola_jackknife_results.csv")
if (file.exists(ref_jack)) {
  rj <- utils::read.csv(ref_jack, stringsAsFactors = FALSE)
  cj <- merge(rj[, c("variable", "gain_drop")], jack[, c("variable", "delta_gain")],
              by = "variable")
  cj$difference <- cj$gain_drop - cj$delta_gain
  print(cj)
  message(sprintf("Published full-model gain %.7f | this run %.7f",
                  unique(rj$gain_full)[1], full_gain))
}

## [verify-sensitivity]
ref_sens <- file.path(PROJ, "reference_results", "M_graminicola_sensitivity_comparison.csv")
if (file.exists(ref_sens)) {
  rs <- utils::read.csv(ref_sens, stringsAsFactors = FALSE)

  published <- data.frame(
    Metric = c("Training AUC", "Validation AUC (mean)", "Validation CBI (mean)",
               "AUC difference (train - val)", "Omission rate (10th pct)",
               "EPPO presence hit rate (%)", "Slovakia suitability"),
    Integrated_published = c(rs$AUC_train[1], rs$AUC_val[1], rs$CBI_val[1],
                             rs$auc_diff[1], rs$or_10p[1],
                             rs$EPPO_hit_rate_pct[1], rs$Slovakia_suitability[1]),
    GBIFonly_published   = c(rs$AUC_train[2], rs$AUC_val[2], rs$CBI_val[2],
                             rs$auc_diff[2], rs$or_10p[2],
                             rs$EPPO_hit_rate_pct[2], rs$Slovakia_suitability[2])
  )

  side_by_side <- merge(published, sensitivity, by = "Metric", all.x = TRUE)
  side_by_side$Integrated_delta <- side_by_side$Integrated - side_by_side$Integrated_published
  side_by_side$GBIFonly_delta   <- side_by_side$GBIF_only  - side_by_side$GBIFonly_published
  print(side_by_side, digits = 4)

  utils::write.csv(side_by_side, p_out("VERIFICATION_vs_published.csv"), row.names = FALSE)
}


# Fitting the reported model without EPPO access --------------------------

## [offline-reproduction]  (not run: kept for reference, as in the .qmd)
# frozen <- utils::read.csv(p_data("M_graminicola_occurrence_points_FROZEN.csv"),
#                           stringsAsFactors = FALSE)
#
# occ_f <- frozen |>
#   transmute(
#     lon = Longitude,
#     lat = Latitude,
#     source = ifelse(grepl("EPPO", Data_Source), "EPPO", "GBIF"),
#     coordinatePrecision = dplyr::case_when(
#       grepl("^State",   Coordinate_Precision) ~ "state-centroid",
#       grepl("^Country", Coordinate_Precision) ~ "country-centroid",
#       TRUE                                    ~ "exact"),
#     weight = unname(WEIGHTS[coordinatePrecision])
#   )
#
# env_f <- terra::extract(
#   env_full, terra::vect(occ_f, geom = c("lon", "lat"), crs = "EPSG:4326"))
# pres_f <- bind_cols(occ_f, as.data.frame(env_f)[, -1, drop = FALSE]) |>
#   mutate(presence = 1)
#
# md_f   <- bind_rows(pres_f |> select(lon, lat, weight, all_of(PRED_VARS), presence),
#                     bg_data |> select(lon, lat, weight, all_of(PRED_VARS), presence))
# pr_f   <- md_f[, PRED_VARS, drop = FALSE]
# rs_f   <- md_f$presence
#
# mx_f <- maxnet::maxnet(p = rs_f, data = pr_f,
#                        f = maxnet::maxnet.formula(rs_f, pr_f, classes = "lq"),
#                        regmult = 1)
# suit_f <- terra::predict(env_full, mx_f, type = "cloglog", na.rm = TRUE)
#
# train_f <- terra::extract(
#   suit_f, terra::vect(occ_f, geom = c("lon", "lat"), crs = "EPSG:4326"))[, 2]
# message(sprintf("Threshold from frozen points: %.6f",
#                 stats::quantile(train_f, probs = 0.10, na.rm = TRUE)))


# Reproducibility =========================================================


# Which numbers are exactly reproducible ----------------------------------


# Session information -----------------------------------------------------

## [session-info]
sessionInfo()


# Handoff to the future-projection document -------------------------------

## [save-handoff]
handoff <- list(
  species          = SPECIES,
  species_slug     = SPECIES_SLUG,
  cfg              = CFG,
  seeds            = SEED,
  pred_vars        = PRED_VARS,
  bioclim_vars     = CFG$bioclim_vars,
  best_fc          = BEST_FC,
  best_rm          = BEST_RM,
  threshold_10p    = as.numeric(threshold_10p),
  threshold_source = "10th percentile of cloglog suitability at the thinned training presences, 5-predictor model",
  pct_suitable     = as.numeric(pct_suitable),
  novel_pct        = as.numeric(novel_pct),
  n_occ            = nrow(occ_thinned),
  study_ext        = as.vector(terra::ext(env_full)),
  grid             = list(ncol = terra::ncol(env_full),
                          nrow = terra::nrow(env_full),
                          res  = terra::res(env_full),
                          crs  = terra::crs(env_full, proj = TRUE)),
  ref_env          = ref_env,
  bg_points        = bg_points,   # read by the revision analysis
  occ_coords       = occ_thinned[, c("lon", "lat", "source",
                                     "coordinatePrecision", "weight")],
  r_version        = paste(R.version$major, R.version$minor, sep = "."),
  built            = format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")
)

saveRDS(handoff, p_out("analysis_handoff.rds"))
message("Handoff written: ", p_out("analysis_handoff.rds"))

## [write-manifest]
manifest <- data.frame(
  file  = list.files(c(DIR_OUT, DIR_FIG), full.names = FALSE),
  path  = list.files(c(DIR_OUT, DIR_FIG), full.names = TRUE)
)
manifest$size_kb  <- round(file.size(manifest$path) / 1024, 1)
manifest$modified <- format(file.mtime(manifest$path), "%Y-%m-%d %H:%M:%S")
utils::write.csv(manifest[, c("file", "size_kb", "modified")],
                 p_out("OUTPUT_MANIFEST.csv"), row.names = FALSE)
knitr::kable(manifest[, c("file", "size_kb")], caption = "Generated outputs.")


# Appendix: approaches tried and not used =================================

## [appendix-spocc]  (not run: kept for reference, as in the .qmd)
# # Aggregated multi-provider search. Returned 51 records, 43 after deduplication,
# # essentially all of them the same GBIF records already retrieved directly, and
# # without the provenance detail needed for precision weighting. Abandoned.
# sp <- spocc::occ(query = SPECIES,
#                  from  = c("gbif", "inat", "bison", "vertnet", "idigbio"),
#                  limit = 500, has_coords = TRUE)
# spocc::occ2df(sp)
#
# # iDigBio queried directly returned zero records: this species is not
# # represented in digitised natural history collections.
# ridigbio::idig_search_records(rq = list(scientificname = SPECIES))

## [appendix-eppo-download]  (not run: kept for reference, as in the .qmd)
# # Non-functional upstream: EPPO changed its data service and the packaged
# # downloader no longer resolves the SQLite archive. Manual download is required,
# # as handled in Section on EPPO records.
# pestr::eppo_database_download()

## [appendix-glad-cropland]  (not run: kept for reference, as in the .qmd)
# # GLAD cropland was downloaded during exploration but not used; ESA WorldCover
# # was preferred for its consistent global 30 arc-second grid.
# geodata::cropland(source = "glad", year = 2015, path = DIR_ENV)

## [appendix-earlier-models]  (not run: kept for reference, as in the .qmd)
# # Two earlier fits, superseded by the five-predictor model. Retained here because
# # their AICc values are the basis for the statement that adding cropland fraction
# # improved AICc by 15.3.
#
# # Four predictors, unthinned occurrences, rm = 1.5. AICc 2979.737.
# maxnet::maxnet(p = resp_unthinned, data = pred_unthinned,
#                f = maxnet::maxnet.formula(resp_unthinned, pred_unthinned,
#                                           classes = "lq"),
#                regmult = 1.5)
#
# # Four predictors, thinned occurrences, rm = 1. AICc 2900.279; threshold 0.196,
# # which is the value quoted in the manuscript text.
# maxnet::maxnet(p = resp_thinned4, data = pred_thinned4,
#                f = maxnet::maxnet.formula(resp_thinned4, pred_thinned4,
#                                           classes = "lq"),
#                regmult = 1)

## [appendix-jackknife-singlevar]  (not run: kept for reference, as in the .qmd)
# # Single-variable fits did not converge for any predictor, including with
# # linear-only features, which is why gain_only is NA throughout the jackknife
# # table. Recorded so the NAs are not mistaken for a coding error.
# for (v in PRED_VARS) {
#   d <- predictors[, v, drop = FALSE]
#   try(maxnet::maxnet(p = response, data = d,
#                      f = maxnet::maxnet.formula(response, d, classes = "l"),
#                      regmult = BEST_RM))
# }



}  # end of PART 1
###############################################################################
# PART 2 - REVISION ANALYSES (M_graminicola_revision_analyses.qmd) ############
###############################################################################

if (RUN_PART[["revision"]]) {



# Scope and how to read this document =====================================


# Two facts about the submitted model that the revision must state --------


# Setup ===================================================================

## [rev-packages]
required_rev <- c("terra", "geodata", "predicts", "maxnet", "glmnet",
                  "dplyr", "ggplot2", "knitr")
missing_rev <- required_rev[!vapply(required_rev, requireNamespace, logical(1),
                                    quietly = TRUE)]
if (length(missing_rev)) {
  stop("Missing packages: ", paste(missing_rev, collapse = ", "), call. = FALSE)
}
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(maxnet)
})
select <- dplyr::select


# Paths and settings ------------------------------------------------------

## [rev-config]
PROJ     <- normalizePath(".", winslash = "/", mustWork = FALSE)
DIR_DATA <- file.path(PROJ, "data")
DIR_ENV  <- file.path(PROJ, "env_data")
DIR_OUT  <- file.path(PROJ, "outputs")
DIR_FIG  <- file.path(PROJ, "figures")
DIR_REV  <- file.path(DIR_OUT, "revision")
for (d in c(DIR_DATA, DIR_ENV, DIR_OUT, DIR_FIG, DIR_REV)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}
p_data <- function(...) file.path(DIR_DATA, ...)
p_out  <- function(...) file.path(DIR_OUT, ...)
p_rev  <- function(...) file.path(DIR_REV, ...)
p_fig  <- function(...) file.path(DIR_FIG, ...)

REV <- list(
  # ---- location uncertainty (Comment 4) ------------------------------------
  # A cell is treated as potential host land if it carries rice (MapSPAM physical
  # rice area >= rice_min_ha); units without rice fall back to cropland fraction
  # >= cropland_min, and units with neither fall back to all land cells.
  rice_min_ha       = 1,
  cropland_min      = 0.05,
  pts_per_unit      = 20,     # points drawn per administrative unit
  n_redraws         = 5,      # independent redraws used to check stability
  r0_km             = 50,     # tolerance radius of the area-based weight
  max_cells_eval    = 5000,   # host cells sampled per unit for unit-level tests
  max_cells_unit    = 50000,  # cells kept per unit (area-weighted subsample above this)
  fallback_buffer_km= 100,    # polygon used when an EPPO unit cannot be matched

  # Pre-specified primary specification of the revised model. It is fixed here,
  # before the sensitivity results are seen, and reported whatever they show.
  primary = list(location = "poly_host", weighting = "area", drop_nested = TRUE),

  # ---- modelling -------------------------------------------------------------
  feature_classes   = c("L", "LQ"),
  reg_multipliers   = c(1, 1.5, 2, 3),
  n_bg              = 1000,   # as in the submitted model
  n_bg_eval         = 2000,   # background drawn inside each withheld area
  aicc_cells        = 1e5,    # cells used to normalise raw output for AICc
  threshold_q       = 0.10,

  # ---- withheld validation (Comments 2 and 3) ------------------------------
  group_folds       = 5,      # country-grouped folds
  group_reps        = 10,     # repetitions with different random grouping
  min_gbif_train    = 10,     # GBIF-only model is not fitted below this

  # ---- spatial partitions (Comment 5) --------------------------------------
  block_sizes_km    = c(250, 500, 1000, 2000, 4482),  # last = autocorrelation range
  block_reps        = 20,     # random fold assignments per block size
  buffer_km         = c(0, 250, 500, 1000, 2000),

  # ---- rice production systems (Comment 1) ---------------------------------
  spam_var          = "phys_area",   # MapSPAM 2010 v2r0 physical area (ha)
  adoption          = c(0.10, 0.25, 0.50),  # share of irrigated rice moving to AWD / DSR
  # Yield-loss fractions for the production-at-risk scenario. Left NA on purpose:
  # set them from the field loss studies cited in the manuscript before running.
  loss_rates        = c(low = NA_real_, mid = NA_real_, high = NA_real_),

  # ---- optional local inputs (used instead of downloads when present) ------
  admin0_file       = p_data("admin0.gpkg"),
  admin1_file       = p_data("admin1.gpkg"),
  rice_dir          = p_data("spam"),              # folder with SPAM rice GeoTIFFs
  crosswalk_manual  = p_data("eppo_unit_crosswalk_manual.csv"),
  independent_file  = p_data("independent_literature_records.csv"),
  # Optional country-specific adoption shares of AWD/DSR (ISO3, share)
  adoption_file     = p_data("water_saving_adoption.csv")
)

SEED <- 2026


# Inputs from the main analysis -------------------------------------------

## [rev-inputs]
HO <- readRDS(p_out("analysis_handoff.rds"))
SPECIES      <- HO$species
SPECIES_SLUG <- HO$species_slug
PRED_VARS    <- HO$pred_vars
BIOCLIM_VARS <- HO$bioclim_vars

env_full <- terra::rast(p_out(paste0(SPECIES_SLUG, "_predictors.tif")))
names(env_full) <- PRED_VARS
mx_published   <- readRDS(p_out(paste0(SPECIES_SLUG, "_FINAL_model.rds")))
suit_published <- terra::rast(p_out(paste0(SPECIES_SLUG, "_FINAL_suitability.tif")))
names(suit_published) <- "suitability"
THR_PUBLISHED  <- HO$threshold_10p

# Pre-thinning records with complete predictor coverage (GBIF and EPPO), and the
# thinned set actually used for fitting.
occ_with_env <- utils::read.csv(p_out(paste0(SPECIES_SLUG, "_occurrences_with_env.csv")))
occ_thinned  <- HO$occ_coords

eppo_geo <- utils::read.csv(p_data(paste0(gsub(" ", "_", SPECIES), "_EPPO_geocoded.csv")))
names(eppo_geo) <- make.names(names(eppo_geo))
eppo_geo$eppo_id <- sprintf("E%03d", seq_len(nrow(eppo_geo)))

message(sprintf("Thinned records: %d (%d GBIF, %d EPPO); EPPO rows: %d",
                nrow(occ_thinned), sum(occ_thinned$source == "GBIF"),
                sum(occ_thinned$source == "EPPO"), nrow(eppo_geo)))

# Background: the same draw as the submitted model (seed 42, 1,000 cells).
if (!is.null(HO$bg_points)) {
  bg_points <- HO$bg_points
} else {
  set.seed(42)
  bg_points <- terra::spatSample(env_full, size = REV$n_bg, method = "random",
                                 na.rm = TRUE, xy = TRUE)
}
bg_df <- bg_points |> rename(lon = x, lat = y)
stopifnot(all(PRED_VARS %in% names(bg_df)))

# Land mask, cell areas and one large fixed cell sample used for AICc
# normalisation and for area-weighted prevalence.
LAND      <- !is.na(sum(env_full))
CELL_AREA <- terra::cellSize(env_full[[1]], unit = "km", mask = TRUE)
set.seed(SEED)
aicc_samp <- terra::spatSample(c(env_full, CELL_AREA), size = REV$aicc_cells,
                               method = "random", na.rm = TRUE, xy = TRUE,
                               cells = TRUE)
names(aicc_samp)[names(aicc_samp) == "area"] <- "area_km2"
N_LAND_CELLS <- as.numeric(terra::global(LAND, "sum", na.rm = TRUE)[1, 1])


# Modelling engine with case weights ======================================

## [rev-engine]
wsd <- function(x, w) {
  if (sum(w) <= 0 || length(x) < 2) return(0)
  mu <- sum(w * x) / sum(w)
  v  <- sum(w * (x - mu)^2) / (sum(w) - sum(w^2) / sum(w))
  sqrt(max(v, 0))
}

reg_w <- function(p, m, w, n_eff) {
  isproduct   <- function(x) grepl(":", x) & !grepl("\\(", x)
  isquadratic <- function(x) grepl("^I\\(.*\\^2\\)", x)
  lregtable <- list(c(0, 10, 30, 100), c(1, 1, 0.2, 0.05))
  qregtable <- list(c(0, 10, 17, 30, 100), c(1.3, 0.8, 0.5, 0.25, 0.05))
  pregtable <- list(c(0, 10, 17, 30, 100), c(2.6, 1.6, 0.9, 0.55, 0.05))
  mm <- m[p == 1, , drop = FALSE]
  wp <- w[p == 1]
  tab <- lregtable
  if (sum(isquadratic(colnames(mm)))) tab <- qregtable
  if (sum(isproduct(colnames(mm))))   tab <- pregtable
  classreg <- vapply(colnames(mm), function(n)
    stats::approx(tab[[1]], tab[[2]], n_eff, rule = 2)$y, numeric(1)) / sqrt(n_eff)
  sds <- vapply(seq_len(ncol(mm)), function(j) wsd(mm[, j], wp), numeric(1))
  pmax(0.001 * (apply(m, 2, max) - apply(m, 2, min)), sds * classreg)
}

maxnet_w <- function(p, data, w = NULL, f, regmult = 1, n_eff = NULL) {
  if (anyNA(data)) stop("NA values in data table.")
  data <- as.data.frame(data)
  if (is.null(w)) w <- rep(1, length(p))
  w[p == 0] <- 1
  if (is.null(n_eff)) n_eff <- sum(p == 1)

  # Add to the background every presence row not already there (as maxnet does),
  # carrying the presence's weight.
  key   <- do.call(paste, c(data, sep = "\r"))
  pres  <- which(p == 1)
  inbg  <- key[pres] %in% key[p == 0]
  addk  <- unique(key[pres][!inbg])
  if (length(addk)) {
    idx  <- pres[match(addk, key[pres])]
    wadd <- vapply(addk, function(k) mean(w[pres][key[pres] == k]), numeric(1))
    data <- rbind(data, data[idx, , drop = FALSE])
    p    <- c(p, rep(0, length(idx)))
    w    <- c(w, unname(wadd))
  }

  mm  <- stats::model.matrix(f, data)
  reg <- reg_w(p, mm, w, n_eff) * regmult
  ow  <- ifelse(p == 1, w, 100 * w)
  glmnet::glmnet.control(pmin = 1.0e-8, fdev = 0)
  lam <- 10^(seq(4, 0, length.out = 200)) * sum(reg) / length(reg) *
    sum(ow[p == 1]) / sum(ow)
  model <- glmnet::glmnet(x = mm, y = as.factor(p), family = "binomial",
                          standardize = FALSE, penalty.factor = reg,
                          lambda = lam, weights = ow)
  class(model) <- c("maxnet", class(model))
  if (length(model$lambda) < 200) stop("glmnet failed to complete the regularisation path.")
  bb <- model$beta[, 200]
  model$betas <- bb[bb != 0]
  model$alpha <- 0
  bgi <- p == 0
  rr  <- as.vector(predict(model, data[bgi, , drop = FALSE], type = "exponential",
                           clamp = FALSE))
  wb  <- w[bgi]
  raw <- rr / sum(wb * rr)
  model$entropy <- -sum(wb * raw * log(raw))
  model$alpha   <- -log(sum(wb * rr))
  model$penalty.factor <- reg
  model$featuremins <- apply(mm, 2, min)
  model$featuremaxs <- apply(mm, 2, max)
  model$varmin <- apply(data, 2, min)
  model$varmax <- apply(data, 2, max)
  model$samplemeans <- colSums(data[p == 1, , drop = FALSE] * w[p == 1]) / sum(w[p == 1])
  model$levels <- lapply(data, levels)
  model$n_eff  <- n_eff
  model
}

# predict.maxnet builds its design matrix from the model's own terms, so extra
# columns in `df` are ignored.
pred_cll <- function(m, df) as.vector(predict(m, as.data.frame(df), type = "cloglog",
                                              clamp = TRUE))
pred_raw <- function(m, df) as.vector(predict(m, as.data.frame(df), type = "exponential",
                                              clamp = TRUE))

# Threshold: the q-th quantile of suitability over training *records*, each GBIF
# point and each administrative unit counting once. A unit's value is the
# weighted mean suitability over the points drawn for it. Location weights decide
# how much a record shapes the fitted surface, not how much it counts when the
# surface is cut into suitable and unsuitable.
rec_threshold <- function(m, pres, q) {
  s  <- pred_cll(m, pres)
  rv <- tapply(s * pres$w, pres$rec_id, sum) / tapply(pres$w, pres$rec_id, sum)
  as.numeric(stats::quantile(rv, q, type = 7))
}

wquantile <- function(x, w, q) {
  o <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  x[which(cw >= q)[1]]
}

## [rev-engine-check]
# Equivalence check against maxnet() on the submitted model's own data.
pres_chk <- occ_thinned
env_chk  <- terra::extract(env_full, terra::vect(pres_chk, geom = c("lon", "lat"),
                                                 crs = "EPSG:4326"))[, PRED_VARS]
dat_chk  <- rbind(env_chk, bg_df[, PRED_VARS])
p_chk    <- c(rep(1, nrow(env_chk)), rep(0, nrow(bg_df)))
f_chk    <- maxnet::maxnet.formula(p_chk, dat_chk, classes = "lq")
m_ref    <- maxnet::maxnet(p_chk, dat_chk, f_chk, regmult = 1)
m_new    <- maxnet_w(p_chk, dat_chk, w = NULL, f = f_chk, regmult = 1)
common   <- union(names(m_ref$betas), names(m_new$betas))
b_ref    <- setNames(rep(0, length(common)), common); b_ref[names(m_ref$betas)] <- m_ref$betas
b_new    <- setNames(rep(0, length(common)), common); b_new[names(m_new$betas)] <- m_new$betas
message(sprintf("maxnet() vs maxnet_w() with unit weights: max |beta difference| = %.2e",
                max(abs(b_ref - b_new))))
stopifnot(max(abs(b_ref - b_new)) < 1e-6)


# Tuning, AICc and evaluation metrics -------------------------------------

## [rev-metrics]
auc_w <- function(pres, bg, w = NULL) {
  if (!length(pres) || !length(bg)) return(NA_real_)
  if (is.null(w)) w <- rep(1, length(pres))
  bg <- sort(bg)
  lo <- findInterval(pres, bg, left.open = TRUE)   # bg strictly below
  le <- findInterval(pres, bg)                     # bg below or equal
  f  <- (lo + 0.5 * (le - lo)) / length(bg)
  sum(w * f) / sum(w)
}

boyce_w <- function(pres, bg, w = NULL, res = 100) {
  if (length(pres) < 3) return(NA_real_)
  if (is.null(w)) w <- rep(1, length(pres))
  lo <- min(c(pres, bg)); hi <- max(c(pres, bg))
  ww <- (hi - lo) / 10
  mv <- seq(lo, hi - ww, length.out = res + 1)
  pe <- vapply(mv, function(a) {
    pi <- sum(w[pres >= a & pres <= a + ww]) / sum(w)
    ei <- mean(bg >= a & bg <= a + ww)
    if (ei == 0) NA_real_ else pi / ei
  }, numeric(1))
  ok <- !is.na(pe)
  if (sum(ok) < 3) return(NA_real_)
  suppressWarnings(stats::cor(pe[ok], mv[ok] + ww / 2, method = "spearman"))
}

aicc_w <- function(m, pres_df, w, n_eff, samp = aicc_samp) {
  raw_p <- pred_raw(m, pres_df)
  raw_s <- pred_raw(m, samp)
  log_norm <- log(mean(raw_s)) + log(N_LAND_CELLS)
  ll <- sum(w * (log(raw_p) - log_norm))
  k  <- length(m$betas)
  n  <- n_eff
  if (k >= n - 1) return(NA_real_)
  2 * k - 2 * ll + (2 * k * (k + 1)) / (n - k - 1)
}

# Area-weighted share of the land in the extent predicted suitable.
prevalence <- function(m, thr, samp = aicc_samp) {
  s <- pred_cll(m, samp)
  sum(samp$area_km2 * (s >= thr)) / sum(samp$area_km2)
}

## [rev-fit-tuned]
# Fit every fc x rm candidate on one training set and keep the lowest AICc.
# pres: data frame with PRED_VARS, w (weight) and rec_id (one id per independent
# record, shared by all points drawn for one administrative unit).
# With `fixed` (a list of fc and rm) no tuning is done and, because the fit is
# only needed for cross-validation, AICc and prevalence are not computed.
fit_tuned <- function(pres, bg = bg_df, fcs = REV$feature_classes,
                      rms = REV$reg_multipliers, fixed = NULL,
                      vars = PRED_VARS, samp = aicc_samp) {
  light <- !is.null(fixed)
  dat <- rbind(pres[, vars], bg[, vars])
  p   <- c(rep(1, nrow(pres)), rep(0, nrow(bg)))
  w   <- c(pres$w, rep(1, nrow(bg)))
  n_eff <- length(unique(pres$rec_id))
  grid <- if (is.null(fixed)) expand.grid(fc = fcs, rm = rms, stringsAsFactors = FALSE)
          else data.frame(fc = fixed$fc, rm = fixed$rm, stringsAsFactors = FALSE)
  fits <- lapply(seq_len(nrow(grid)), function(i) {
    f <- maxnet::maxnet.formula(p, dat, classes = tolower(grid$fc[i]))
    m <- tryCatch(maxnet_w(p, dat, w = w, f = f, regmult = grid$rm[i], n_eff = n_eff),
                  error = function(e) NULL)
    if (is.null(m)) return(NULL)
    list(m = m, aicc = if (light) NA_real_ else aicc_w(m, pres, pres$w, n_eff, samp),
         fc = grid$fc[i], rm = grid$rm[i])
  })
  fits <- Filter(Negate(is.null), fits)
  if (!length(fits)) return(NULL)
  aic  <- vapply(fits, `[[`, numeric(1), "aicc")
  best <- if (all(is.na(aic))) fits[[1]] else fits[[which.min(aic)]]
  best$thr   <- rec_threshold(best$m, pres, REV$threshold_q)
  best$prev  <- if (light) NA_real_ else prevalence(best$m, best$thr, samp)
  best$n_eff <- n_eff
  best$tuning <- data.frame(fc = vapply(fits, `[[`, "", "fc"),
                            rm = vapply(fits, `[[`, numeric(1), "rm"),
                            AICc = aic, ncoef = vapply(fits, function(z) length(z$m$betas), 1L))
  best
}


# Administrative units, duplication and overlap ===========================

## [rev-admin-polygons]
cc <- geodata::country_codes()
eppo_geo$iso2 <- if ("country.code" %in% names(eppo_geo)) toupper(eppo_geo$country.code) else NA
eppo_geo$iso3 <- cc$ISO3[match(eppo_geo$iso2, cc$ISO2)]
miss_iso <- is.na(eppo_geo$iso3)
if (any(miss_iso)) {
  eppo_geo$iso3[miss_iso] <- cc$ISO3[match(tolower(eppo_geo$country[miss_iso]),
                                          tolower(cc$NAME))]
}
if (anyNA(eppo_geo$iso3)) {
  warning("EPPO rows without an ISO3 code: ",
          paste(eppo_geo$country[is.na(eppo_geo$iso3)], collapse = ", "))
}
eppo_geo$level <- ifelse(!is.na(eppo_geo$state) & nzchar(eppo_geo$state), 1L, 0L)

ext_env <- terra::ext(env_full)

if (file.exists(REV$admin0_file)) {
  adm0 <- terra::vect(REV$admin0_file)
} else {
  adm0 <- geodata::world(resolution = 2, level = 0, path = DIR_ENV)
}
adm0 <- terra::crop(adm0, ext_env)
names(adm0)[names(adm0) %in% c("GID_0", "ISO3", "iso3")][1] <- "GID_0"

iso_l1 <- sort(unique(eppo_geo$iso3[eppo_geo$level == 1L & !is.na(eppo_geo$iso3)]))
if (file.exists(REV$admin1_file)) {
  adm1 <- terra::vect(REV$admin1_file)
} else {
  adm1 <- do.call(rbind, lapply(iso_l1, function(i) {
    g <- tryCatch(geodata::gadm(country = i, level = 1, path = DIR_ENV, resolution = 2),
                  error = function(e) NULL)
    if (is.null(g)) message("  [gadm] level 1 not available for ", i)
    g
  }))
}
stopifnot(all(c("GID_0", "GID_1", "NAME_1") %in% names(adm1)))
if (!"VARNAME_1" %in% names(adm1)) adm1$VARNAME_1 <- NA_character_
message(sprintf("Admin polygons: %d countries, %d subdivisions", nrow(adm0), nrow(adm1)))

## [rev-crosswalk]
norm_nm <- function(x) {
  x <- iconv(as.character(x), to = "ASCII//TRANSLIT")
  x <- tolower(gsub("[^A-Za-z]", "", x))
  x[is.na(x)] <- ""
  x
}

# Island groups used by EPPO for Indonesia, mapped to GADM provinces by pattern.
IDN_GROUPS <- c(
  java        = "^Jawa|Banten|Jakarta|Yogyakarta",
  sumatra     = "Sumatera|Aceh|^Riau$|Jambi|Bengkulu|Lampung|Bangka",
  sulawesi    = "Sulawesi|Gorontalo",
  kalimantan  = "Kalimantan",
  irianjaya   = "Papua",
  nusatenggara= "Nusa Tenggara|^Bali$",
  moluccas    = "Maluku"
)

adm1_df <- as.data.frame(adm1)[, c("GID_0", "GID_1", "NAME_1", "VARNAME_1")]
eppo_pt <- terra::vect(eppo_geo[!is.na(eppo_geo$eppo_lat), ],
                       geom = c("eppo_lon", "eppo_lat"), crs = "EPSG:4326")
spatial_hit <- terra::extract(adm1[, c("GID_1")], eppo_pt)
spatial_gid <- setNames(rep(NA_character_, nrow(eppo_geo)), eppo_geo$eppo_id)
spatial_gid[eppo_pt$eppo_id[spatial_hit$id.y]] <- spatial_hit$GID_1

xw <- lapply(seq_len(nrow(eppo_geo)), function(i) {
  r <- eppo_geo[i, ]
  if (r$level == 0L) {
    return(data.frame(eppo_id = r$eppo_id, gid = r$iso3, method = "country",
                      flag = ifelse(r$iso3 %in% adm0$GID_0, "", "country polygon not found")))
  }
  a  <- adm1_df[adm1_df$GID_0 == r$iso3, , drop = FALSE]
  st <- norm_nm(r$state)
  if (r$iso3 == "IDN" && st %in% names(IDN_GROUPS)) {
    g <- a$GID_1[grepl(IDN_GROUPS[[st]], a$NAME_1)]
    return(data.frame(eppo_id = r$eppo_id, gid = paste(g, collapse = ";"),
                      method = "island group table",
                      flag = ifelse(length(g), "", "no province matched")))
  }
  nm  <- norm_nm(a$NAME_1)
  vnm <- strsplit(as.character(a$VARNAME_1), "\\|")
  exact <- which(nm == st | vapply(vnm, function(v) st %in% norm_nm(v), logical(1)))
  method <- "name"
  if (!length(exact)) {
    exact <- agrep(st, nm, max.distance = 0.15)
    method <- "fuzzy name"
  }
  sp <- spatial_gid[[r$eppo_id]]
  if (length(exact) == 1) {
    g <- a$GID_1[exact]
    flag <- if (!is.na(sp) && sp != g) paste("centroid falls in", sp) else ""
  } else if (!is.na(sp)) {
    g <- sp; method <- "centroid location"
    flag <- if (length(exact) > 1) "several name matches" else "no name match"
  } else {
    g <- ""; method <- "unmatched"; flag <- "no match: buffer used"
  }
  data.frame(eppo_id = r$eppo_id, gid = g, method = method, flag = flag)
})
xw <- dplyr::bind_rows(xw)
xw <- merge(eppo_geo[, c("eppo_id", "country", "state", "Status", "iso3", "level")],
            xw, by = "eppo_id")
utils::write.csv(xw, p_rev("eppo_unit_crosswalk_DRAFT.csv"), row.names = FALSE)

if (file.exists(REV$crosswalk_manual)) {
  man <- utils::read.csv(REV$crosswalk_manual, stringsAsFactors = FALSE)
  hit <- match(man$eppo_id, xw$eppo_id)
  xw$gid[hit]    <- man$gid
  xw$method[hit] <- "manual"
  xw$flag[hit]   <- ifelse(is.na(man$note), "", man$note)
  message("Manual crosswalk applied to ", sum(!is.na(hit)), " rows")
}
print(table(xw$method))
if (any(nzchar(xw$flag))) print(xw[nzchar(xw$flag), c("eppo_id", "country", "state", "gid", "method", "flag")])

## [rev-unit-geometry]
# One polygon per EPPO record; area, equivalent radius, and whether the geocoded
# centroid used in the submitted model lies inside it.
unit_poly <- function(i) {
  r <- xw[i, ]
  g <- strsplit(r$gid, ";")[[1]]
  v <- if (r$level == 0L) adm0[adm0$GID_0 %in% g] else adm1[adm1$GID_1 %in% g]
  if (!nrow(v)) {
    e <- eppo_geo[eppo_geo$eppo_id == r$eppo_id, ]
    v <- terra::buffer(terra::vect(e, geom = c("eppo_lon", "eppo_lat"), crs = "EPSG:4326"),
                       width = REV$fallback_buffer_km * 1000)
  }
  v <- terra::aggregate(v)
  v$eppo_id <- r$eppo_id
  v[, "eppo_id"]
}
units <- do.call(rbind, lapply(seq_len(nrow(xw)), unit_poly))
units <- merge(units, xw, by = "eppo_id")
units$area_km2 <- terra::expanse(units, unit = "km")
units$r_eq_km  <- sqrt(units$area_km2 / pi)
inside <- terra::relate(terra::vect(eppo_geo[match(units$eppo_id, eppo_geo$eppo_id), ],
                                    geom = c("eppo_lon", "eppo_lat"), crs = "EPSG:4326"),
                        units, "intersects", pairs = FALSE)
units$centroid_inside <- diag(as.matrix(inside))
terra::writeVector(units, p_rev("eppo_units.gpkg"), overwrite = TRUE)


# Host mask from MapSPAM rice area ----------------------------------------

## [rev-rice-layers]
# MapSPAM 2010 v2r0, physical rice area (ha per 5 arc-minute cell), by production
# system. Brought onto the 2.5 arc-minute grid by splitting each cell into four
# and dividing by four, which conserves area.
spam_local <- list.files(REV$rice_dir, pattern = "\\.tif$", full.names = TRUE)
if (length(spam_local)) {
  spam <- terra::rast(spam_local)
} else {
  spam <- geodata::crop_spam(crop = "rice", var = REV$spam_var, path = DIR_ENV)
}
sys_of <- function(nm) {
  dplyr::case_when(grepl("all$|_A$", nm) ~ "all",
                   grepl("irrigated$|_I$", nm) ~ "irrigated",
                   grepl("highinput$|_H$", nm) ~ "rainfed_high",
                   grepl("lowinput$|_L$", nm)  ~ "rainfed_low",
                   grepl("subsistence$|_S$", nm) ~ "rainfed_subsistence",
                   grepl("rainfed$|_R$", nm) ~ "rainfed",
                   TRUE ~ NA_character_)
}
names(spam) <- sys_of(names(spam))
spam <- spam[[c("all", "irrigated")]]
spam <- terra::crop(spam, terra::extend(ext_env, 1))
spam <- terra::disagg(spam, fact = 2) / 4
spam <- terra::resample(spam, env_full[[1]], method = "near")
spam <- terra::mask(terra::classify(spam, cbind(NA, 0)), env_full[[1]])
RICE <- c(spam[["all"]], spam[["irrigated"]], spam[["all"]] - spam[["irrigated"]])
names(RICE) <- c("rice_all_ha", "rice_irrigated_ha", "rice_rainfed_ha")
RICE <- terra::ifel(RICE < 0, 0, RICE)
terra::writeRaster(RICE, p_rev("rice_phys_area_systems_2p5min.tif"), overwrite = TRUE)

message(sprintf("Rice physical area in extent: %.1f Mha (irrigated %.1f, rainfed %.1f)",
                terra::global(RICE[[1]], "sum", na.rm = TRUE)[1, 1] / 1e6,
                terra::global(RICE[[2]], "sum", na.rm = TRUE)[1, 1] / 1e6,
                terra::global(RICE[[3]], "sum", na.rm = TRUE)[1, 1] / 1e6))

HOST_RICE <- RICE[["rice_all_ha"]] >= REV$rice_min_ha
HOST_CROP <- env_full[["cropland_fraction"]] >= REV$cropland_min

## [rev-unit-cells]
# All candidate cells of each unit with predictor values, area and host class.
stack_all <- c(env_full, CELL_AREA, HOST_RICE, HOST_CROP, RICE)
names(stack_all) <- c(PRED_VARS, "area_km2", "host_rice", "host_crop", names(RICE))

unit_cells <- lapply(seq_len(nrow(units)), function(i) {
  d <- terra::extract(stack_all, units[i], cells = TRUE, xy = TRUE, ID = FALSE)
  d <- d[stats::complete.cases(d[, PRED_VARS]), , drop = FALSE]
  if (!nrow(d)) return(NULL)
  if (nrow(d) > REV$max_cells_unit) {          # very large countries: area-weighted subsample
    set.seed(SEED)
    d <- d[sample.int(nrow(d), REV$max_cells_unit, prob = d$area_km2), , drop = FALSE]
  }
  d$host_class <- dplyr::case_when(any(d$host_rice == 1) ~ "rice",
                                   any(d$host_crop == 1) ~ "cropland",
                                   TRUE ~ "land")
  d$is_host <- switch(d$host_class[1],
                      rice = d$host_rice == 1, cropland = d$host_crop == 1,
                      land = rep(TRUE, nrow(d)))
  d$eppo_id <- units$eppo_id[i]
  d
})
names(unit_cells) <- units$eppo_id
units$n_cells   <- vapply(unit_cells, function(d) if (is.null(d)) 0L else nrow(d), 1L)
units$host_class <- vapply(unit_cells, function(d) if (is.null(d)) NA_character_ else d$host_class[1], "")
units$host_km2   <- vapply(unit_cells, function(d) if (is.null(d)) 0 else sum(d$area_km2[d$is_host]), 1)
print(table(units$host_class, units$level, useNA = "ifany"))


# Duplication and overlap audit -------------------------------------------

## [rev-overlap-audit]
gbif_pts <- occ_with_env |> filter(source == "GBIF")
gv <- terra::vect(gbif_pts, geom = c("lon", "lat"), crs = "EPSG:4326")
gin <- terra::relate(gv, units, "intersects")           # GBIF x units
units$n_gbif_inside <- colSums(gin)

pres_rows <- !grepl("^Absent", units$Status)
countries_with_l1 <- unique(units$iso3[units$level == 1L & pres_rows])
units$nested_national <- units$level == 0L & units$iso3 %in% countries_with_l1 & pres_rows
dup_key <- paste(units$level, units$gid)
units$duplicate_of <- ifelse(duplicated(dup_key), units$eppo_id[match(dup_key, dup_key)], "")

overlap_tab <- as.data.frame(units)[, c("eppo_id", "country", "state", "Status", "level", "gid",
                                        "method", "flag", "area_km2", "r_eq_km",
                                        "centroid_inside", "host_class", "host_km2",
                                        "n_gbif_inside", "nested_national", "duplicate_of")]
utils::write.csv(overlap_tab, p_rev("Table_R1_overlap_audit.csv"), row.names = FALSE)

audit <- c(
  `EPPO presence records`                          = sum(pres_rows),
  `  national`                                     = sum(pres_rows & units$level == 0L),
  `  subnational`                                  = sum(pres_rows & units$level == 1L),
  `National records nested within subnational ones`= sum(units$nested_national),
  `Duplicated units`                               = sum(nzchar(units$duplicate_of)),
  `EPPO units containing >= 1 GBIF record`         = sum(pres_rows & units$n_gbif_inside > 0),
  `GBIF records inside an EPPO presence unit`      = sum(rowSums(gin[, pres_rows, drop = FALSE]) > 0),
  `Geocoded centroids outside their own unit`      = sum(pres_rows & !units$centroid_inside),
  `Units without rice; cropland host mask used`    = sum(pres_rows & units$host_class == "cropland", na.rm = TRUE),
  `Units with neither; all land used`              = sum(pres_rows & units$host_class == "land", na.rm = TRUE)
)
audit_df <- data.frame(Item = names(audit), n = as.integer(audit), check.names = FALSE)
utils::write.csv(audit_df, p_rev("Table_R1b_overlap_summary.csv"), row.names = FALSE)
knitr::kable(audit_df, caption = "Table R1. Duplication and overlap between and within the two sources.")


# Location assignment and weighting variants ==============================

## [rev-variant-builder]
gbif_train <- occ_thinned |> filter(source == "GBIF")
gbif_env   <- terra::extract(env_full, terra::vect(gbif_train, geom = c("lon", "lat"),
                                                   crs = "EPSG:4326"))[, PRED_VARS]
gbif_train <- cbind(gbif_train[, c("lon", "lat")], gbif_env)
gbif_train$rec_id <- paste0("G", seq_len(nrow(gbif_train)))
gbif_train$w      <- 1
gbif_train$iso3   <- terra::extract(adm0[, "GID_0"],
                                    terra::vect(gbif_train, geom = c("lon", "lat"),
                                                crs = "EPSG:4326"))$GID_0
gbif_train$kind   <- "GBIF"

unit_weight <- function(u, scheme) {
  switch(scheme,
         none      = rep(1, nrow(u)),
         published = ifelse(u$level == 1L, 0.5, 0.25),
         steep     = ifelse(u$level == 1L, 0.25, 0.10),
         area      = pmin(1, REV$r0_km / u$r_eq_km))
}

eppo_units_pres <- function(drop_nested = TRUE) {
  u <- as.data.frame(units)
  u <- u[!grepl("^Absent", u$Status) & u$n_cells > 0 & !nzchar(u$duplicate_of), ]
  if (drop_nested) u <- u[!u$nested_national, ]
  u
}

# Build the EPPO part of a training set.
eppo_points <- function(u, location, scheme, draw = 1) {
  wts <- unit_weight(u, scheme)
  if (location == "centroid") {
    e <- eppo_geo[match(u$eppo_id, eppo_geo$eppo_id), ]
    d <- data.frame(lon = e$eppo_lon, lat = e$eppo_lat)
    ev <- terra::extract(env_full, terra::vect(d, geom = c("lon", "lat"), crs = "EPSG:4326"))
    d <- cbind(d, ev[, PRED_VARS])
    d$rec_id <- u$eppo_id; d$w <- wts; d$iso3 <- u$iso3; d$kind <- "EPPO"
    return(d[stats::complete.cases(d[, PRED_VARS]), ])
  }
  set.seed(SEED + draw)
  dplyr::bind_rows(lapply(seq_len(nrow(u)), function(i) {
    cells <- unit_cells[[u$eppo_id[i]]]
    pool  <- if (location == "poly_host") cells[cells$is_host, , drop = FALSE] else cells
    if (!nrow(pool)) pool <- cells
    k   <- min(REV$pts_per_unit, nrow(pool))
    pick <- pool[sample.int(nrow(pool), k, prob = pool$area_km2), , drop = FALSE]
    data.frame(lon = pick$x, lat = pick$y, pick[, PRED_VARS],
               rec_id = u$eppo_id[i], w = wts[i] / k, iso3 = u$iso3[i], kind = "EPPO")
  }))
}

training_set <- function(location, scheme, drop_nested = TRUE, draw = 1,
                         include_eppo = TRUE, exclude_iso3 = character(0)) {
  g <- gbif_train[!(gbif_train$iso3 %in% exclude_iso3), ]
  if (!include_eppo) return(g)
  u <- eppo_units_pres(drop_nested)
  u <- u[!(u$iso3 %in% exclude_iso3), ]
  rbind(g[, c("lon", "lat", PRED_VARS, "rec_id", "w", "iso3", "kind")],
        eppo_points(u, location, scheme, draw)[, c("lon", "lat", PRED_VARS,
                                                    "rec_id", "w", "iso3", "kind")])
}

## [rev-variants-fit]
variant_grid <- rbind(
  expand.grid(location = c("centroid", "poly_host", "poly_land"),
              weighting = c("none", "published", "steep", "area"),
              drop_nested = TRUE, stringsAsFactors = FALSE),
  data.frame(location = c("centroid", REV$primary$location),
             weighting = c("none", REV$primary$weighting),
             drop_nested = FALSE)
)
variant_grid$id <- with(variant_grid, paste(location, weighting,
                                            ifelse(drop_nested, "dedup", "all"), sep = "_"))
PRIMARY_ID <- with(REV$primary, paste(location, weighting,
                                      ifelse(drop_nested, "dedup", "all"), sep = "_"))
SUBMITTED_ID <- "centroid_none_all"   # what the submitted model actually was

# Fixed evaluation background spanning the extent, for correlation and AUC summaries.
set.seed(SEED + 1)
eval_cells <- aicc_samp[sample.int(nrow(aicc_samp), min(20000, nrow(aicc_samp))), ]

variant_fits <- lapply(seq_len(nrow(variant_grid)), function(i) {
  v  <- variant_grid[i, ]
  tr <- training_set(v$location, v$weighting, v$drop_nested)
  ft <- fit_tuned(tr)
  ft$train <- tr
  ft$id <- v$id
  message(sprintf("  %-28s n_rows = %4d, n_eff = %3d, fc = %s, rm = %s, thr = %.3f",
                  v$id, nrow(tr), ft$n_eff, ft$fc, ft$rm, ft$thr))
  ft
})
names(variant_fits) <- variant_grid$id
saveRDS(variant_fits, p_rev("variant_fits.rds"))

## [rev-variants-table]
s_sub <- pred_cll(variant_fits[[SUBMITTED_ID]]$m, eval_cells)
variant_tab <- dplyr::bind_rows(lapply(names(variant_fits), function(id) {
  ft <- variant_fits[[id]]
  s  <- pred_cll(ft$m, eval_cells)
  tr <- ft$train
  # importance of cropland: drop in training AUC when cropland is permuted
  set.seed(SEED)
  pb <- pred_cll(ft$m, bg_df)
  base_auc <- auc_w(pred_cll(ft$m, tr), pb, tr$w)
  trp <- tr; trp$cropland_fraction <- sample(trp$cropland_fraction)
  bgp <- bg_df; bgp$cropland_fraction <- sample(bgp$cropland_fraction)
  perm_auc <- auc_w(pred_cll(ft$m, trp), pred_cll(ft$m, bgp), tr$w)
  data.frame(
    variant = id,
    location = variant_grid$location[variant_grid$id == id],
    weighting = variant_grid$weighting[variant_grid$id == id],
    drop_nested = variant_grid$drop_nested[variant_grid$id == id],
    rows = nrow(tr), n_eff = ft$n_eff,
    eppo_weight_share = round(sum(tr$w[tr$kind == "EPPO"]) / sum(tr$w), 3),
    fc = ft$fc, rm = ft$rm, ncoef = length(ft$m$betas),
    threshold = round(ft$thr, 4),
    suitable_pct_area = round(100 * ft$prev, 1),
    train_auc = round(base_auc, 3),
    cropland_auc_drop = round(base_auc - perm_auc, 3),
    cor_with_submitted = round(stats::cor(s, s_sub), 3)
  )
}))
utils::write.csv(variant_tab, p_rev("Table_R2_variants.csv"), row.names = FALSE)
knitr::kable(variant_tab, caption = "Table R2. Training-data variants: location assignment, weighting and duplication handling.")

## [rev-redraw-stability]
# Polygon-sampled variants depend on the random draw. Refit the primary
# specification on independent redraws and compare the surfaces.
redraws <- lapply(seq_len(REV$n_redraws), function(d) {
  tr <- training_set(REV$primary$location, REV$primary$weighting,
                     REV$primary$drop_nested, draw = d + 100)
  fit_tuned(tr)
})
S_red <- sapply(redraws, function(ft) pred_cll(ft$m, eval_cells))
cor_red <- stats::cor(S_red)
redraw_summary <- data.frame(
  redraws = REV$n_redraws,
  min_pairwise_cor = round(min(cor_red[upper.tri(cor_red)]), 3),
  suitable_pct_range = paste(round(100 * range(vapply(redraws, `[[`, 1, "prev")), 1),
                             collapse = " - "),
  thresholds = paste(round(vapply(redraws, `[[`, 1, "thr"), 3), collapse = ", ")
)
utils::write.csv(redraw_summary, p_rev("Table_R2b_redraw_stability.csv"), row.names = FALSE)
knitr::kable(redraw_summary)

## [rev-primary-model]
PRIMARY <- variant_fits[[PRIMARY_ID]]
suit_primary <- terra::predict(env_full, PRIMARY$m, type = "cloglog", clamp = TRUE,
                               na.rm = TRUE)
names(suit_primary) <- "suitability"
terra::writeRaster(suit_primary, p_out(paste0(SPECIES_SLUG, "_REVISED_suitability.tif")),
                   overwrite = TRUE)
saveRDS(PRIMARY$m, p_out(paste0(SPECIES_SLUG, "_REVISED_model.rds")))
THR_PRIMARY <- PRIMARY$thr
bin_primary <- suit_primary >= THR_PRIMARY

ref_env_primary <- PRIMARY$train[, PRED_VARS]
saveRDS(list(model_file = p_out(paste0(SPECIES_SLUG, "_REVISED_model.rds")),
             suitability_file = p_out(paste0(SPECIES_SLUG, "_REVISED_suitability.tif")),
             threshold_10p = THR_PRIMARY,
             threshold_source = "Weighted 10th percentile of cloglog suitability at the training presences of the revised primary model",
             ref_env = ref_env_primary,
             spec = REV$primary, fc = PRIMARY$fc, rm = PRIMARY$rm),
        p_out("revision_handoff.rds"))
message(sprintf("Revised primary model %s: fc = %s, rm = %s, threshold = %.4f, suitable = %.1f%% of land",
                PRIMARY_ID, PRIMARY$fc, PRIMARY$rm, THR_PRIMARY, 100 * PRIMARY$prev))
message(sprintf("Submitted model threshold (five predictors): %.4f", THR_PUBLISHED))


# Withheld validation =====================================================

## [rev-regions]
region_of <- function(iso3) {
  r1 <- cc$UNREGION1[match(iso3, cc$ISO3)]
  co <- cc$continent[match(iso3, cc$ISO3)]
  dplyr::case_when(r1 == "Southern Asia" ~ "South Asia",
                   r1 == "Eastern Asia" ~ "East Asia",
                   r1 == "South-Eastern Asia" ~ "South-East Asia",
                   co %in% c("North America", "South America") ~ "Americas",
                   co == "Africa" ~ "Africa",
                   co == "Europe" ~ "Europe",
                   TRUE ~ "Other")
}
units$region      <- region_of(units$iso3)
gbif_train$region <- region_of(gbif_train$iso3)
print(table(units$region[!grepl("^Absent", units$Status)]))
print(table(gbif_train$region))

# Numeric country raster plus lookup table (avoids version-dependent handling of
# categorical rasters).
adm0$zone   <- seq_len(nrow(adm0))
COUNTRY_LUT <- data.frame(zone = adm0$zone, iso3 = adm0$GID_0)
COUNTRY_R   <- terra::rasterize(adm0, env_full[[1]], field = "zone")
iso_at <- function(lon, lat) {
  z <- terra::extract(COUNTRY_R, cbind(lon, lat))[, 1]
  COUNTRY_LUT$iso3[match(z, COUNTRY_LUT$zone)]
}

## [rev-evaluate-heldout]
# Background inside a set of countries.
bg_in <- function(iso3s, n = REV$n_bg_eval, seed = SEED) {
  v <- adm0[adm0$GID_0 %in% iso3s]
  if (!nrow(v)) return(NULL)
  m <- terra::mask(terra::crop(env_full, v), v)
  set.seed(seed)
  d <- terra::spatSample(m, size = n, method = "random", na.rm = TRUE, xy = TRUE)
  names(d)[1:2] <- c("lon", "lat")
  d
}

# Unit-level host cells, subsampled for speed.
unit_eval_cells <- function(id) {
  d <- unit_cells[[id]]
  d <- d[d$is_host, , drop = FALSE]
  if (nrow(d) > REV$max_cells_eval) {
    set.seed(SEED)
    d <- d[sample.int(nrow(d), REV$max_cells_eval, prob = d$area_km2), ]
  }
  d
}

# Score one fitted model on one withheld set.
score_model <- function(ft, test_units, test_gbif, test_iso3, bgE) {
  if (is.null(ft)) return(NULL)
  m <- ft$m; thr <- ft$thr
  # unit-level
  su <- vapply(test_units$eppo_id, function(id) {
    d <- unit_eval_cells(id)
    s <- pred_cll(m, d)
    sum(d$area_km2 * (s >= thr)) / sum(d$area_km2)
  }, numeric(1))
  smean <- vapply(test_units$eppo_id, function(id) {
    d <- unit_eval_cells(id); s <- pred_cll(m, d)
    sum(d$area_km2 * s) / sum(d$area_km2)
  }, numeric(1))
  # centroid hits
  ec <- eppo_geo[match(test_units$eppo_id, eppo_geo$eppo_id), ]
  ce <- terra::extract(env_full, terra::vect(ec, geom = c("eppo_lon", "eppo_lat"),
                                             crs = "EPSG:4326"))[, PRED_VARS]
  ok <- stats::complete.cases(ce)
  sc <- if (any(ok)) pred_cll(m, ce[ok, , drop = FALSE]) else numeric(0)
  # local prevalence within the withheld countries
  sb <- pred_cll(m, bgE)
  prev_local <- mean(sb >= thr)
  # point AUC / CBI: GBIF test points + host-cell samples (unit-weighted)
  pts <- list(); wts <- c()
  if (nrow(test_gbif)) { pts[[1]] <- pred_cll(m, test_gbif); wts <- c(wts, rep(1, nrow(test_gbif))) }
  for (id in test_units$eppo_id) {
    d <- unit_eval_cells(id)
    if (nrow(d) > 50) d <- d[seq(1, nrow(d), length.out = 50), ]
    pts[[length(pts) + 1]] <- pred_cll(m, d)
    wts <- c(wts, rep(1 / nrow(d), nrow(d)))
  }
  ps <- unlist(pts)
  k  <- sum(sc >= thr); n <- length(sc)
  list(
    unit = data.frame(eppo_id = test_units$eppo_id, S_u = su, mean_suit = smean),
    summary = data.frame(
      n_units = nrow(test_units), n_gbif = nrow(test_gbif),
      auc = if (length(ps)) auc_w(ps, sb, wts) else NA_real_,
      cbi = if (length(ps)) boyce_w(ps, sb, wts) else NA_real_,
      unit_hit_rate = mean(su >= 0.5), unit_mean_S = mean(su),
      gbif_hit_rate = if (nrow(test_gbif)) mean(pred_cll(m, test_gbif) >= thr) else NA_real_,
      centroid_hit_rate = if (n) k / n else NA_real_,
      local_prevalence = prev_local,
      binom_p = if (n) stats::pbinom(k - 1, n, prev_local, lower.tail = FALSE) else NA_real_,
      fc = ft$fc, rm = ft$rm, threshold = thr, n_eff_train = ft$n_eff
    )
  )
}

run_heldout <- function(test_iso3, label, seed = SEED) {
  tu <- eppo_units_pres(TRUE)
  tu <- tu[tu$iso3 %in% test_iso3, , drop = FALSE]
  tg <- gbif_train[gbif_train$iso3 %in% test_iso3, , drop = FALSE]
  if (nrow(tu) + nrow(tg) < 3) return(NULL)
  bgE <- bg_in(test_iso3, seed = seed)
  if (is.null(bgE) || nrow(bgE) < 100) return(NULL)
  bg_train <- bg_df[!(iso_at(bg_df$lon, bg_df$lat) %in% test_iso3), ]

  tr_primary <- training_set(REV$primary$location, REV$primary$weighting, TRUE,
                             exclude_iso3 = test_iso3)
  tr_submit  <- training_set("centroid", "none", FALSE, exclude_iso3 = test_iso3)
  tr_gbif    <- training_set("centroid", "none", TRUE, include_eppo = FALSE,
                             exclude_iso3 = test_iso3)
  fits <- list(
    revised   = fit_tuned(tr_primary, bg = bg_train),
    submitted = fit_tuned(tr_submit,  bg = bg_train),
    gbif_only = if (nrow(tr_gbif) >= REV$min_gbif_train) fit_tuned(tr_gbif, bg = bg_train) else NULL
  )
  res <- lapply(names(fits), function(nm) {
    s <- score_model(fits[[nm]], tu, tg, test_iso3, bgE)
    if (is.null(s)) return(NULL)
    s$summary$model <- nm
    s$unit$model    <- nm
    s
  })
  res <- Filter(Negate(is.null), res)
  list(label = label, test_iso3 = test_iso3,
       summary = dplyr::bind_rows(lapply(res, `[[`, "summary")) |> mutate(design = label),
       unit = dplyr::bind_rows(lapply(res, `[[`, "unit")) |> mutate(design = label),
       train_env = list(revised = tr_primary[, PRED_VARS], gbif_only = tr_gbif[, PRED_VARS]),
       test_units = tu, test_gbif = tg)
}

## [rev-design-A]
all_iso <- unique(c(units$iso3, gbif_train$iso3))
non_asia <- all_iso[!(cc$continent[match(all_iso, cc$ISO3)] %in% "Asia")]
design_A <- run_heldout(non_asia, "A: outside Asia")
knitr::kable(design_A$summary, digits = 3,
             caption = "Design A. Records outside Asia withheld from every model.")

## [rev-design-B]
REGIONS <- c("South Asia", "East Asia", "South-East Asia", "Americas", "Africa", "Europe")
design_B <- lapply(REGIONS, function(r) {
  iso <- all_iso[region_of(all_iso) == r]
  message("Withholding ", r, " (", length(iso), " countries)")
  run_heldout(iso, paste("B:", r))
})
design_B <- Filter(Negate(is.null), design_B)
B_sum <- dplyr::bind_rows(lapply(design_B, `[[`, "summary"))
knitr::kable(B_sum, digits = 3, caption = "Design B. Leave one region out.")

## [rev-design-C]
rec_iso <- unique(c(units$iso3[!grepl("^Absent", units$Status)], gbif_train$iso3))
rec_iso <- rec_iso[!is.na(rec_iso)]
design_C <- list()
for (rp in seq_len(REV$group_reps)) {
  set.seed(SEED + rp)
  fold <- sample(rep(seq_len(REV$group_folds), length.out = length(rec_iso)))
  for (k in seq_len(REV$group_folds)) {
    r <- run_heldout(rec_iso[fold == k], sprintf("C: rep %d fold %d", rp, k), seed = SEED + rp * 10 + k)
    if (!is.null(r)) { r$summary$rep <- rp; r$unit$rep <- rp; design_C[[length(design_C) + 1]] <- r }
  }
}
C_unit <- dplyr::bind_rows(lapply(design_C, `[[`, "unit"))
C_sum  <- dplyr::bind_rows(lapply(design_C, `[[`, "summary"))

## [rev-withheld-tables]
# Pooled unit-level results per design and model, with the paired test.
# Each unit enters once: S_u is first averaged over the repetitions in which the
# unit was withheld, so repeated folds do not inflate the sample size.
paired <- function(u, a, b) {
  u <- u[u$model %in% c(a, b), ] |> group_by(eppo_id, model) |>
    summarise(S_u = mean(S_u), .groups = "drop") |> as.data.frame()
  w <- reshape(u, idvar = "eppo_id", timevar = "model", direction = "wide")
  ca <- paste0("S_u.", a); cb <- paste0("S_u.", b)
  w <- w[stats::complete.cases(w[, c(ca, cb)]), ]
  if (nrow(w) < 3) return(c(n = nrow(w), diff = NA, p = NA))
  c(n = nrow(w), diff = mean(w[[ca]] - w[[cb]]),
    p = suppressWarnings(stats::wilcox.test(w[[ca]], w[[cb]], paired = TRUE)$p.value))
}

summ_design <- function(S, U, name) {
  S <- S |> group_by(model) |>
    summarise(test_sets = n(), units = sum(n_units), gbif_points = sum(n_gbif),
              auc_sd = stats::sd(auc, na.rm = TRUE), auc = mean(auc, na.rm = TRUE),
              cbi = mean(cbi, na.rm = TRUE),
              unit_hit = weighted.mean(unit_hit_rate, n_units, na.rm = TRUE),
              unit_S = weighted.mean(unit_mean_S, n_units, na.rm = TRUE),
              centroid_hit = weighted.mean(centroid_hit_rate, n_units, na.rm = TRUE),
              local_prev = weighted.mean(local_prevalence, n_units, na.rm = TRUE),
              .groups = "drop") |>
    mutate(design = name)
  if (is.null(U$rep)) U$rep <- 1L
  pr <- rbind(revised_vs_gbif = paired(U, "revised", "gbif_only"),
              submitted_vs_gbif = paired(U, "submitted", "gbif_only"),
              revised_vs_submitted = paired(U, "revised", "submitted"))
  list(summary = S, paired = data.frame(design = name, comparison = rownames(pr), pr,
                                        row.names = NULL))
}

wA <- summ_design(design_A$summary, design_A$unit, "A: outside Asia")
wB <- summ_design(B_sum, dplyr::bind_rows(lapply(design_B, `[[`, "unit")), "B: leave one region out")
wC <- summ_design(C_sum, C_unit, "C: country-grouped folds")
table_R3 <- dplyr::bind_rows(wA$summary, wB$summary, wC$summary)
table_R3b <- dplyr::bind_rows(wA$paired, wB$paired, wC$paired)
utils::write.csv(table_R3,  p_rev("Table_R3_withheld_validation.csv"), row.names = FALSE)
utils::write.csv(table_R3b, p_rev("Table_R3b_paired_comparison.csv"), row.names = FALSE)
utils::write.csv(B_sum, p_rev("Table_R3c_leave_one_region_out_detail.csv"), row.names = FALSE)
knitr::kable(table_R3, digits = 3, caption = "Table R3. Withheld validation: the same test records for every model.")
knitr::kable(table_R3b, digits = 4, caption = "Table R3b. Paired comparison of unit suitable share (Wilcoxon signed-rank).")

## [rev-absence-and-independent]
# The single confirmed absence, scored by every model.
absent <- eppo_geo[grepl("^Absent", eppo_geo$Status), ]
abs_tab <- NULL
if (nrow(absent)) {
  ae <- terra::extract(env_full, terra::vect(absent, geom = c("eppo_lon", "eppo_lat"),
                                             crs = "EPSG:4326"))[, PRED_VARS]
  au <- units[units$eppo_id %in% absent$eppo_id]
  abs_tab <- dplyr::bind_rows(lapply(c(PRIMARY_ID, SUBMITTED_ID), function(id) {
    ft <- variant_fits[[id]]
    d  <- unit_eval_cells(au$eppo_id[1])
    s  <- pred_cll(ft$m, d)
    data.frame(model = id, country = absent$country,
               centroid_suitability = round(pred_cll(ft$m, ae), 4),
               threshold = round(ft$thr, 4),
               share_of_host_area_suitable = round(sum(d$area_km2 * (s >= ft$thr)) / sum(d$area_km2), 3))
  }))
  utils::write.csv(abs_tab, p_rev("Table_R3d_confirmed_absence.csv"), row.names = FALSE)
  knitr::kable(abs_tab)
}

# Optional: georeferenced records from the literature that are in neither GBIF
# nor EPPO. Columns: lon, lat, country_iso3, year, reference, coord_uncertainty_km.
indep_tab <- NULL
if (file.exists(REV$independent_file)) {
  ind <- utils::read.csv(REV$independent_file, stringsAsFactors = FALSE)
  iv  <- terra::vect(ind, geom = c("lon", "lat"), crs = "EPSG:4326")
  # Independence check: distance to the nearest training record of any kind.
  tr_all <- PRIMARY$train
  dmin <- apply(terra::distance(iv, terra::vect(tr_all, geom = c("lon", "lat"),
                                                crs = "EPSG:4326")), 1, min) / 1000
  in_gbif_cell <- terra::cellFromXY(env_full, as.matrix(ind[, c("lon", "lat")])) %in%
    terra::cellFromXY(env_full, as.matrix(gbif_train[, c("lon", "lat")]))
  ind$dist_nearest_training_km <- round(dmin, 1)
  ind$same_cell_as_gbif <- in_gbif_cell
  ind$independent <- dmin > 10 & !in_gbif_cell
  ie <- terra::extract(env_full, iv)[, PRED_VARS]
  for (id in c(PRIMARY_ID, SUBMITTED_ID, "gbif_only")) {
    ft <- if (id == "gbif_only") fit_tuned(training_set("centroid", "none", TRUE, include_eppo = FALSE)) else variant_fits[[id]]
    ind[[paste0("suit_", id)]] <- round(pred_cll(ft$m, ie), 4)
    ind[[paste0("hit_", id)]]  <- pred_cll(ft$m, ie) >= ft$thr
  }
  utils::write.csv(ind, p_rev("Table_R3e_independent_literature_records.csv"), row.names = FALSE)
  indep_tab <- ind
  message(sprintf("Independent literature records: %d (%d pass the independence check)",
                  nrow(ind), sum(ind$independent)))
}


# Spatial partitions and separation =======================================

## [rev-blocks]
moll <- "+proj=moll +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs"
to_moll <- function(lon, lat) terra::project(cbind(lon, lat), from = "EPSG:4326", to = moll)

hav_km <- function(lon1, lat1, lon2, lat2) {
  rad <- pi / 180
  dlat <- outer(lat1, lat2, "-") * rad
  dlon <- outer(lon1, lon2, "-") * rad
  a <- sin(dlat / 2)^2 + outer(cos(lat1 * rad), cos(lat2 * rad)) * sin(dlon / 2)^2
  2 * 6371 * asin(pmin(sqrt(a), 1))
}

assign_blocks <- function(rep_xy, bg_xy, size_km, k, rule, seed) {
  sz <- size_km * 1000
  bid <- function(xy) paste(floor(xy[, 1] / sz), floor(xy[, 2] / sz))
  pb <- bid(rep_xy); bb <- bid(bg_xy)
  blocks <- unique(c(pb, bb))
  bxy <- do.call(rbind, strsplit(blocks, " "))
  ord <- order(as.numeric(bxy[, 2]), as.numeric(bxy[, 1]))
  set.seed(seed)
  one <- function() {
    if (rule == "systematic") { f <- integer(length(blocks)); f[ord] <- rep_len(seq_len(k), length(blocks)); f }
    else sample(rep_len(seq_len(k), length(blocks)))
  }
  if (rule == "balanced") {
    cand <- replicate(100, one(), simplify = FALSE)
    sdv  <- vapply(cand, function(f) {
      tab <- tabulate(f[match(pb, blocks)], k); if (any(tab == 0)) Inf else stats::sd(tab)
    }, numeric(1))
    fb <- cand[[which.min(sdv)]]
  } else fb <- one()
  list(pres = fb[match(pb, blocks)], bg = fb[match(bb, blocks)])
}

cv_partition <- function(train, size_km, rule, seed, fixed) {
  rep_pts <- train |> group_by(rec_id) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
  rxy <- to_moll(rep_pts$lon, rep_pts$lat)
  bxy <- to_moll(bg_df$lon, bg_df$lat)
  k <- 5
  a <- assign_blocks(rxy, bxy, size_km, k, rule, seed)
  fold_p <- a$pres[match(train$rec_id, rep_pts$rec_id)]
  res <- lapply(seq_len(k), function(j) {
    trn <- train[fold_p != j, ]; tst <- train[fold_p == j, ]
    bgt <- bg_df[a$bg != j, ];   bgv <- bg_df[a$bg == j, ]
    if (!nrow(tst) || nrow(bgv) < 20 || length(unique(trn$rec_id)) < 10) return(NULL)
    ft <- fit_tuned(trn, bg = bgt, fixed = fixed)
    if (is.null(ft)) return(NULL)
    pt <- pred_cll(ft$m, tst); pv <- pred_cll(ft$m, bgv)
    rt <- trn |> group_by(rec_id) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
    rv <- tst |> group_by(rec_id) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
    sep <- apply(hav_km(rv$lon, rv$lat, rt$lon, rt$lat), 1, min)
    data.frame(fold = j, n_train = nrow(rt), n_test = nrow(rv),
               auc = auc_w(pt, pv, tst$w), cbi = boyce_w(pt, pv, tst$w),
               or10 = sum(tst$w[pt < ft$thr]) / sum(tst$w),
               auc_train = auc_w(pred_cll(ft$m, trn), pred_cll(ft$m, bgt), trn$w),
               sep_median_km = stats::median(sep))
  })
  dplyr::bind_rows(res)
}

block_runs <- list()
for (spec in c(PRIMARY_ID, SUBMITTED_ID)) {
  ft0 <- variant_fits[[spec]]
  fixed <- list(fc = ft0$fc, rm = ft0$rm)
  for (sz in REV$block_sizes_km) for (rule in c("systematic", "random", "balanced")) {
    reps <- if (rule == "systematic") 1 else REV$block_reps
    for (rp in seq_len(reps)) {
      r <- cv_partition(ft0$train, sz, rule, SEED + rp, fixed)
      if (nrow(r)) block_runs[[length(block_runs) + 1]] <-
          cbind(spec = spec, block_km = sz, rule = rule, rep = rp, r)
    }
  }
}
block_folds <- dplyr::bind_rows(block_runs)
utils::write.csv(block_folds, p_rev("Table_R4_block_cv_folds.csv"), row.names = FALSE)

table_R4 <- block_folds |>
  group_by(spec, block_km, rule) |>
  summarise(partitions = n_distinct(rep), folds = n(),
            auc_sd = stats::sd(auc, na.rm = TRUE), auc = mean(auc, na.rm = TRUE),
            cbi = mean(cbi, na.rm = TRUE), or10 = mean(or10, na.rm = TRUE),
            auc_diff = mean(auc_train - auc, na.rm = TRUE),
            min_test = min(n_test), max_test = max(n_test),
            sep_median_km = stats::median(sep_median_km), .groups = "drop")
utils::write.csv(table_R4, p_rev("Table_R4_block_size_sensitivity.csv"), row.names = FALSE)
knitr::kable(table_R4, digits = 3, caption = "Table R4. Sensitivity to block size and fold assignment.")

# Does imbalance explain poor performance at large block sizes?
imb <- block_folds |> filter(spec == PRIMARY_ID, block_km == max(REV$block_sizes_km))
imb_fit <- if (nrow(imb) > 10) stats::lm(auc ~ n_train + sep_median_km, data = imb) else NULL
if (!is.null(imb_fit)) print(summary(imb_fit)$coefficients)

## [rev-buffered-loo]
# Spatial jackknife: each record is predicted by a model fitted without any record
# within d km of it (and without background within d km). AUC is pooled against
# the fixed evaluation background.
buffered_loo <- function(train, d_km, fixed) {
  rep_pts <- train |> group_by(rec_id) |> summarise(lon = mean(lon), lat = mean(lat), .groups = "drop")
  D  <- hav_km(rep_pts$lon, rep_pts$lat, rep_pts$lon, rep_pts$lat)
  Db <- hav_km(rep_pts$lon, rep_pts$lat, bg_df$lon, bg_df$lat)
  out <- lapply(seq_len(nrow(rep_pts)), function(i) {
    keep_rec <- rep_pts$rec_id[D[i, ] > d_km | (d_km == 0 & seq_len(nrow(rep_pts)) != i)]
    keep_rec <- setdiff(keep_rec, rep_pts$rec_id[i])
    trn <- train[train$rec_id %in% keep_rec, ]
    if (length(unique(trn$rec_id)) < 10) return(NULL)
    bgt <- bg_df[Db[i, ] > d_km, ]
    ft  <- fit_tuned(trn, bg = bgt, fixed = fixed)
    if (is.null(ft)) return(NULL)
    tst <- train[train$rec_id == rep_pts$rec_id[i], ]
    pt  <- pred_cll(ft$m, tst)
    data.frame(rec_id = rep_pts$rec_id[i], n_train = length(keep_rec),
               suit = sum(tst$w * pt) / sum(tst$w),
               omitted = sum(tst$w[pt < ft$thr]) / sum(tst$w),
               auc_point = auc_w(pt, pred_cll(ft$m, eval_cells[seq_len(5000), ]), tst$w))
  })
  dplyr::bind_rows(out)
}

loo_runs <- list()
for (spec in c(PRIMARY_ID, SUBMITTED_ID, "gbif_only")) {
  ft0 <- if (spec == "gbif_only") fit_tuned(training_set("centroid", "none", TRUE, include_eppo = FALSE)) else variant_fits[[spec]]
  tr0 <- if (spec == "gbif_only") training_set("centroid", "none", TRUE, include_eppo = FALSE) else ft0$train
  for (d in REV$buffer_km) {
    r <- buffered_loo(tr0, d, list(fc = ft0$fc, rm = ft0$rm))
    if (nrow(r)) loo_runs[[length(loo_runs) + 1]] <- cbind(spec = spec, buffer_km = d, r)
  }
}
loo_all <- dplyr::bind_rows(loo_runs)
table_R5 <- loo_all |> group_by(spec, buffer_km) |>
  summarise(records = n(), median_train = stats::median(n_train),
            omission = mean(omitted), pooled_auc = mean(auc_point), .groups = "drop")
utils::write.csv(loo_all, p_rev("Table_R5_buffered_loo_records.csv"), row.names = FALSE)
utils::write.csv(table_R5, p_rev("Table_R5_buffered_loo.csv"), row.names = FALSE)
knitr::kable(table_R5, digits = 3, caption = "Table R5. Buffered leave-one-out (spatial jackknife).")


# Environmental novelty of withheld records ===============================

## [rev-mess-withheld]
# Univariate MESS, as in predicts:::.messi (Elith et al. 2010), on plain vectors.
messi <- function(p, v) {
  v <- sort(v)
  f <- 100 * findInterval(p, v) / length(v)
  minv <- v[1]; maxv <- v[length(v)]
  res <- 2 * f
  i <- f > 50 & f < 100; res[i] <- 200 - res[i]
  i <- f == 0;   res[i] <- 100 * (p[i] - minv) / (maxv - minv)
  i <- f == 100; res[i] <- 100 * (maxv - p[i]) / (maxv - minv)
  res
}
mess_df <- function(test_env, ref_env) {
  per_var <- sapply(PRED_VARS, function(v) messi(test_env[[v]], ref_env[[v]]))
  per_var <- matrix(per_var, ncol = length(PRED_VARS), dimnames = list(NULL, PRED_VARS))
  data.frame(mess = apply(per_var, 1, min), mod = PRED_VARS[apply(per_var, 1, which.min)])
}

mess_rows <- list()
for (ds in c(list(design_A), design_B)) {
  ref <- ds$train_env$revised
  te  <- dplyr::bind_rows(
    if (nrow(ds$test_gbif)) data.frame(ds$test_gbif[, PRED_VARS], kind = "GBIF point") else NULL,
    dplyr::bind_rows(lapply(ds$test_units$eppo_id, function(id) {
      d <- unit_eval_cells(id); data.frame(d[, PRED_VARS], kind = "EPPO host cell", area = d$area_km2, unit = id)
    })))
  if (!nrow(te)) next
  ms <- mess_df(te[, PRED_VARS], ref)
  te <- cbind(te, ms)
  if (is.null(te$area)) te$area <- 1
  te$area[is.na(te$area)] <- 1
  mess_rows[[length(mess_rows) + 1]] <- data.frame(
    design = ds$label,
    gbif_points_novel = if (any(te$kind == "GBIF point")) sprintf("%d of %d", sum(te$mess[te$kind == "GBIF point"] < 0), sum(te$kind == "GBIF point")) else "-",
    eppo_units = length(unique(te$unit[te$kind == "EPPO host cell"])),
    eppo_host_area_novel_pct = round(100 * sum(te$area[te$kind == "EPPO host cell" & te$mess < 0]) /
                                       max(1e-9, sum(te$area[te$kind == "EPPO host cell"])), 1),
    most_dissimilar_variable = if (any(te$mess < 0))
      names(sort(table(te$mod[te$mess < 0]), decreasing = TRUE))[1] else "-"
  )
}
table_R6 <- dplyr::bind_rows(mess_rows)
utils::write.csv(table_R6, p_rev("Table_R6_mess_withheld.csv"), row.names = FALSE)
knitr::kable(table_R6, caption = "Table R6. Environmental novelty of withheld records relative to the training data that excluded them.")


# Rice production systems, water management and surveillance ==============

## [rev-rice-exposure]
COUNTRY_ZONES <- COUNTRY_R
exposure_by_country <- function(bin, tag) {
  lyr <- c(RICE * bin, RICE)
  names(lyr) <- c("suit_all", "suit_irr", "suit_rf", "all", "irr", "rf")
  z <- terra::zonal(lyr, COUNTRY_ZONES, fun = "sum", na.rm = TRUE)
  names(z)[1] <- "zone"
  z$iso3 <- COUNTRY_LUT$iso3[match(z$zone, COUNTRY_LUT$zone)]
  z <- dplyr::group_by(z, iso3) |> dplyr::summarise(dplyr::across(-zone, sum), .groups = "drop")
  z$model <- tag
  z
}

thr5 <- rec_threshold(PRIMARY$m, PRIMARY$train, 0.05)
bins <- list(
  revised_10p   = suit_primary >= THR_PRIMARY,
  revised_5p    = suit_primary >= thr5,
  submitted_10p = suit_published >= THR_PUBLISHED
)
expo <- dplyr::bind_rows(lapply(names(bins), function(b) exposure_by_country(bins[[b]], b)))

# Variant envelope: share of rice area suitable under each of the 12 variants.
var_share <- vapply(variant_fits, function(ft) {
  s <- terra::predict(env_full, ft$m, type = "cloglog", clamp = TRUE, na.rm = TRUE)
  as.numeric(terra::global(RICE[["rice_all_ha"]] * (s >= ft$thr), "sum", na.rm = TRUE)[1, 1]) /
    as.numeric(terra::global(RICE[["rice_all_ha"]], "sum", na.rm = TRUE)[1, 1])
}, numeric(1))

glob <- expo |> group_by(model) |>
  summarise(rice_Mha = sum(all) / 1e6, suit_Mha = sum(suit_all) / 1e6,
            suit_irr_Mha = sum(suit_irr) / 1e6, suit_rf_Mha = sum(suit_rf) / 1e6,
            .groups = "drop") |>
  mutate(share_suitable_pct = round(100 * suit_Mha / rice_Mha, 1),
         irr_share_of_suitable_pct = round(100 * suit_irr_Mha / suit_Mha, 1))
utils::write.csv(glob, p_rev("Table_R7_rice_exposure_global.csv"), row.names = FALSE)
utils::write.csv(data.frame(variant = names(var_share), share_rice_suitable_pct = round(100 * var_share, 1)),
                 p_rev("Table_R7b_rice_exposure_variants.csv"), row.names = FALSE)
knitr::kable(glob, digits = 2, caption = "Table R7. Rice physical area in climatically suitable cells, by production system.")
message(sprintf("Share of rice area suitable across the %d variants: %.1f%% to %.1f%%",
                length(var_share), 100 * min(var_share), 100 * max(var_share)))

## [rev-rice-country]
recorded_gbif <- unique(gbif_train$iso3)
recorded_eppo <- unique(units$iso3[!grepl("^Absent", units$Status)])
absent_eppo   <- unique(units$iso3[grepl("^Absent", units$Status)])

ctab <- expo |> filter(model == "revised_10p") |>
  transmute(iso3, country = cc$NAME[match(iso3, cc$ISO3)],
            rice_kha = all / 1e3, irrigated_share_pct = 100 * irr / pmax(all, 1),
            suitable_rice_kha = suit_all / 1e3,
            suitable_share_pct = 100 * suit_all / pmax(all, 1),
            suitable_irrigated_kha = suit_irr / 1e3, suitable_rainfed_kha = suit_rf / 1e3) |>
  left_join(expo |> filter(model == "revised_5p") |> transmute(iso3, suitable_share_5p_pct = 100 * suit_all / pmax(all, 1)), by = "iso3") |>
  left_join(expo |> filter(model == "submitted_10p") |> transmute(iso3, suitable_share_submitted_pct = 100 * suit_all / pmax(all, 1)), by = "iso3") |>
  mutate(record_status = dplyr::case_when(iso3 %in% absent_eppo ~ "EPPO: absent",
                                          iso3 %in% recorded_eppo & iso3 %in% recorded_gbif ~ "GBIF + EPPO",
                                          iso3 %in% recorded_eppo ~ "EPPO",
                                          iso3 %in% recorded_gbif ~ "GBIF",
                                          TRUE ~ "no record")) |>
  arrange(desc(rice_kha))
utils::write.csv(ctab, p_rev("Table_R8_rice_exposure_by_country_full.csv"), row.names = FALSE)
knitr::kable(utils::head(ctab, 25), digits = 1,
             caption = "Table R8. The 25 largest rice-growing countries in the extent: rice area in suitable climate and record status.")

## [rev-system-association]
# Exact GBIF records in rice cells: irrigated share at the record versus the
# rice-area-weighted distribution of irrigated share across rice cells of the same
# countries. Uses only exact coordinates, so it is free of any host-mask circularity.
irr_share <- RICE[["rice_irrigated_ha"]] / RICE[["rice_all_ha"]]
names(irr_share) <- "irr_share"
gpts <- terra::vect(gbif_pts, geom = c("lon", "lat"), crs = "EPSG:4326")
gb   <- terra::extract(c(irr_share, RICE[["rice_all_ha"]]), gpts, ID = FALSE)
names(gb) <- c("irr_share", "rice_ha")
gb$iso3 <- iso_at(gbif_pts$lon, gbif_pts$lat)
gb <- gb[!is.na(gb$irr_share) & gb$rice_ha >= REV$rice_min_ha, ]
assoc <- NULL
if (nrow(gb) >= 5) {
  ref <- dplyr::bind_rows(lapply(unique(gb$iso3), function(i) {
    m <- terra::mask(c(irr_share, RICE[["rice_all_ha"]]), terra::ifel(COUNTRY_R == COUNTRY_LUT$zone[COUNTRY_LUT$iso3 == i][1], 1, NA))
    d <- terra::as.data.frame(m, na.rm = TRUE)
    names(d) <- c("irr_share", "rice_ha")
    d <- d[d$rice_ha >= REV$rice_min_ha, ]
    set.seed(SEED)
    d[sample.int(nrow(d), 5000, replace = TRUE, prob = d$rice_ha), , drop = FALSE] |> mutate(iso3 = i)
  }))
  wt <- stats::wilcox.test(gb$irr_share, ref$irr_share)
  assoc <- data.frame(n_records_in_rice = nrow(gb),
                      median_irr_share_records = round(stats::median(gb$irr_share), 3),
                      median_irr_share_rice_landscape = round(stats::median(ref$irr_share), 3),
                      wilcoxon_p = signif(wt$p.value, 3))
  utils::write.csv(assoc, p_rev("Table_R9_system_association.csv"), row.names = FALSE)
  knitr::kable(assoc)
} else message("Fewer than five exact records fall in rice cells; association test skipped.")

## [rev-water-scenarios]
# Conducive rice area = rainfed rice in suitable cells (already aerobic or
# intermittently flooded) + a share a of irrigated rice in suitable cells, for the
# share assumed to move from continuous flooding to AWD or direct seeding.
adopt <- REV$adoption
adopt_cty <- NULL
if (file.exists(REV$adoption_file)) {
  adopt_cty <- utils::read.csv(REV$adoption_file, stringsAsFactors = FALSE)  # iso3, share
}
scen <- ctab |>
  select(iso3, country, rice_kha, suitable_irrigated_kha, suitable_rainfed_kha, record_status)
for (a in adopt) {
  scen[[sprintf("conducive_kha_adopt_%02d", round(100 * a))]] <-
    scen$suitable_rainfed_kha + a * scen$suitable_irrigated_kha
}
if (!is.null(adopt_cty)) {
  scen$adopt_country <- adopt_cty$share[match(scen$iso3, adopt_cty$iso3)]
  scen$conducive_kha_country_adoption <- scen$suitable_rainfed_kha +
    ifelse(is.na(scen$adopt_country), 0, scen$adopt_country) * scen$suitable_irrigated_kha
}
scen_tot <- data.frame(
  scenario = c("Rainfed only (no transition)", sprintf("+ %d%% of irrigated to AWD/DSR", round(100 * adopt))),
  conducive_Mha = c(sum(scen$suitable_rainfed_kha), sum(scen$suitable_rainfed_kha) + adopt * sum(scen$suitable_irrigated_kha)) / 1e3
)
utils::write.csv(scen, p_rev("Table_R10_water_management_by_country.csv"), row.names = FALSE)
utils::write.csv(scen_tot, p_rev("Table_R10b_water_management_totals.csv"), row.names = FALSE)
knitr::kable(scen_tot, digits = 2, caption = "Table R10. Rice area in suitable climate under conducive water regimes.")

# Production at risk, only if loss rates have been set from the literature.
if (all(!is.na(REV$loss_rates))) {
  prod <- geodata::crop_spam(crop = "rice", var = "prod", path = DIR_ENV)
  names(prod) <- sys_of(names(prod))
  prod <- terra::resample(terra::disagg(terra::crop(prod[[c("all", "irrigated")]],
                                                    terra::extend(ext_env, 1)), 2) / 4,
                          env_full[[1]], method = "near")
  prod_rf  <- prod[["all"]] - prod[["irrigated"]]
  suit_bin <- bins$revised_10p
  par <- vapply(adopt, function(a) {
    as.numeric(terra::global(suit_bin * (prod_rf + a * prod[["irrigated"]]), "sum", na.rm = TRUE)[1, 1])
  }, numeric(1))
  risk <- expand.grid(adoption = adopt, loss = names(REV$loss_rates))
  risk$production_at_risk_Mt <- par[match(risk$adoption, adopt)] *
    REV$loss_rates[risk$loss] / 1e6
  utils::write.csv(risk, p_rev("Table_R10c_production_at_risk.csv"), row.names = FALSE)
  knitr::kable(risk, digits = 2)
} else {
  message("REV$loss_rates not set: production-at-risk scenario skipped.")
}

## [rev-surveillance]
prio <- ctab |>
  filter(rice_kha >= 10) |>
  mutate(priority = dplyr::case_when(
    record_status == "no record" & suitable_share_pct >= 50 ~ "1: survey (suitable, rice-growing, unrecorded)",
    record_status == "no record" & suitable_share_pct >= 10 ~ "2: survey where suitable",
    record_status != "no record" & suitable_rainfed_kha + 0.25 * suitable_irrigated_kha >= 100 ~ "3: manage (recorded; large conducive area)",
    TRUE ~ "4: low"),
    robust = abs(suitable_share_pct - suitable_share_5p_pct) < 20) |>
  arrange(priority, desc(suitable_rice_kha))
utils::write.csv(prio, p_rev("Table_R11_surveillance_priorities.csv"), row.names = FALSE)
knitr::kable(prio |> filter(grepl("^1|^2", priority)) |>
               select(country, rice_kha, irrigated_share_pct, suitable_share_pct,
                      suitable_share_5p_pct, suitable_rainfed_kha, priority),
             digits = 1, caption = "Table R11. Unrecorded rice-growing countries with climatically suitable rice area.")

## [rev-host-model]
# Exploratory: replace generic cropland with irrigated and rainfed rice fractions.
rice_frac <- c(RICE[["rice_irrigated_ha"]], RICE[["rice_rainfed_ha"]]) / (CELL_AREA * 100)
names(rice_frac) <- c("rice_irrigated_frac", "rice_rainfed_frac")
env_rice <- c(env_full[[BIOCLIM_VARS]], rice_frac)
RICE_VARS <- names(env_rice)

tr_r <- PRIMARY$train
tr_r <- cbind(tr_r, terra::extract(rice_frac, as.matrix(tr_r[, c("lon", "lat")])))
bg_r <- cbind(bg_df, terra::extract(rice_frac, as.matrix(bg_df[, c("lon", "lat")])))
ok_t <- stats::complete.cases(tr_r[, RICE_VARS]); ok_b <- stats::complete.cases(bg_r[, RICE_VARS])
samp_r <- cbind(aicc_samp, terra::extract(rice_frac, as.matrix(aicc_samp[, c("x", "y")])))
samp_r <- samp_r[stats::complete.cases(samp_r[, RICE_VARS]), ]
host_fit <- fit_tuned(tr_r[ok_t, ], bg = bg_r[ok_b, ], vars = RICE_VARS, samp = samp_r)
host_tab <- data.frame(model = c("climate + cropland (revised primary)", "climate + irrigated/rainfed rice"),
                       AICc = c(min(PRIMARY$tuning$AICc, na.rm = TRUE), min(host_fit$tuning$AICc, na.rm = TRUE)),
                       fc = c(PRIMARY$fc, host_fit$fc), rm = c(PRIMARY$rm, host_fit$rm))
utils::write.csv(host_tab, p_rev("Table_R12_host_model.csv"), row.names = FALSE)
knitr::kable(host_tab, digits = 1)
print(host_fit$m$betas)


# Figures =================================================================

## [rev-fig-helpers]
theme_rev <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          plot.title = element_text(face = "bold"),
          legend.position = "bottom")
}
save_fig <- function(p, file, w = 10, h = 6) {
  ggsave(p_fig(file), p, width = w, height = h, dpi = 300, bg = "white")
  invisible(p_fig(file))
}
world_df <- suppressWarnings(as.data.frame(terra::geom(adm0)))

## [fig-R1-validation]
long <- dplyr::bind_rows(lapply(c("auc", "unit_hit", "centroid_hit"), function(v)
  data.frame(design = table_R3$design, model = table_R3$model, metric = v, value = table_R3[[v]])))
long$metric <- factor(long$metric, c("auc", "unit_hit", "centroid_hit"),
                      c("AUC vs local background", "Unit hit rate (S_u >= 0.5)", "Centroid hit rate"))
figR1 <- ggplot(long, aes(model, value, fill = model)) +
  geom_col(width = 0.65) +
  facet_grid(metric ~ design, scales = "free_y") +
  scale_fill_manual(values = c(revised = "#762A83", submitted = "#C2A5CF", gbif_only = "#1B7837"),
                    guide = "none") +
  labs(title = "Withheld validation: identical test records for every model", x = NULL, y = NULL) +
  theme_rev()
save_fig(figR1, "Figure_R1_withheld_validation.png", 11, 7)
figR1

## [fig-R2-blocks]
figR2 <- ggplot(table_R4, aes(block_km, auc, colour = rule)) +
  geom_line() + geom_point(aes(size = min_test)) +
  facet_wrap(~ spec) +
  scale_x_log10() +
  labs(title = "Validation AUC by block size and fold assignment",
       subtitle = "Point size: smallest test fold (records)", x = "Block size (km, log scale)",
       y = "Mean validation AUC", size = "Smallest test fold") +
  theme_rev()
save_fig(figR2, "Figure_R2_block_sensitivity.png", 10, 5.5)
figR2

## [fig-R3-rice]
agg <- terra::aggregate(c(RICE[["rice_irrigated_ha"]] * bins$revised_10p,
                          RICE[["rice_rainfed_ha"]] * bins$revised_10p),
                        fact = 24, fun = "sum", na.rm = TRUE)   # 1 degree
names(agg) <- c("Irrigated rice in suitable climate", "Rainfed rice in suitable climate")
ad <- terra::as.data.frame(agg, xy = TRUE, na.rm = TRUE)
ad <- dplyr::bind_rows(lapply(names(agg), function(n) data.frame(x = ad$x, y = ad$y, system = n, kha = ad[[n]] / 1e3)))
ad <- ad[ad$kha > 0.5, ]
figR3 <- ggplot() +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = "grey96", colour = "grey75", linewidth = 0.15) +
  geom_tile(data = ad, aes(x, y, fill = kha)) +
  scale_fill_viridis_c(trans = "log10", name = "Rice area per 1° cell (thousand ha)") +
  facet_wrap(~ system, ncol = 1) +
  coord_equal(expand = FALSE) +
  labs(title = "Rice area in climatically suitable cells, by production system", x = NULL, y = NULL) +
  theme_rev()
save_fig(figR3, "Figure_R3_rice_exposure_map.png", 10, 9)
figR3

## [fig-R4-countries]
top <- utils::head(ctab[order(-ctab$suitable_rice_kha), ], 20)
bd <- dplyr::bind_rows(
  data.frame(country = top$country, part = "Rainfed (conducive now)", kha = top$suitable_rainfed_kha),
  data.frame(country = top$country, part = "Irrigated (conducive if moved to AWD/DSR)", kha = top$suitable_irrigated_kha))
bd$country <- factor(bd$country, rev(top$country))
figR4 <- ggplot(bd, aes(kha / 1e3, country, fill = part)) +
  geom_col(width = 0.7) +
  geom_text(data = transform(top, country = factor(country, rev(top$country))),
            aes(x = suitable_rice_kha / 1e3, y = country, label = record_status),
            inherit.aes = FALSE, hjust = -0.1, size = 3) +
  scale_fill_manual(values = c("#E08214", "#2166AC"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
  labs(title = "Rice area in suitable climate, by system", x = "Million ha", y = NULL) +
  theme_rev()
save_fig(figR4, "Figure_R4_rice_exposure_countries.png", 10, 7)
figR4


# Results digest ==========================================================

## [rev-digest]
f1 <- function(x, d = 1) formatC(x, format = "f", digits = d, big.mark = ",")
g <- function(tab, des, mod, col) tab[[col]][tab$design == des & tab$model == mod][1]
rA <- "A: outside Asia"; rC <- "C: country-grouped folds"
lines <- c(
  "# Revision results digest", "",
  sprintf("Generated %s.", format(Sys.time(), "%Y-%m-%d %H:%M")), "",
  "## Data audit (Comment 4)",
  paste0("- ", audit_df$Item, ": ", audit_df$n), "",
  "## Revised primary model",
  sprintf("- Specification %s; fc = %s, rm = %s; n_eff = %d; threshold = %.3f (submitted five-predictor model: %.3f).",
          PRIMARY_ID, PRIMARY$fc, PRIMARY$rm, PRIMARY$n_eff, THR_PRIMARY, THR_PUBLISHED),
  sprintf("- Suitable share of land in the extent: %s%% (submitted model: %s%%).",
          f1(100 * PRIMARY$prev), f1(100 * variant_fits[[SUBMITTED_ID]]$prev)),
  sprintf("- Correlation of the %d variant surfaces with the submitted surface: %s to %s.",
          nrow(variant_tab), f1(min(variant_tab$cor_with_submitted), 2), f1(max(variant_tab$cor_with_submitted), 2)),
  sprintf("- Redraw stability: minimum pairwise correlation %s across %d redraws.",
          redraw_summary$min_pairwise_cor, REV$n_redraws), "",
  "## Withheld validation (Comments 2 and 3)",
  sprintf("- Design A (outside Asia withheld): %d units, %d GBIF points. AUC revised %s, submitted %s, GBIF-only %s. Unit hit rate %s / %s / %s. Local prevalence %s / %s / %s.",
          g(table_R3, rA, "revised", "units"), g(table_R3, rA, "revised", "gbif_points"),
          f1(g(table_R3, rA, "revised", "auc"), 3), f1(g(table_R3, rA, "submitted", "auc"), 3), f1(g(table_R3, rA, "gbif_only", "auc"), 3),
          f1(g(table_R3, rA, "revised", "unit_hit"), 2), f1(g(table_R3, rA, "submitted", "unit_hit"), 2), f1(g(table_R3, rA, "gbif_only", "unit_hit"), 2),
          f1(g(table_R3, rA, "revised", "local_prev"), 2), f1(g(table_R3, rA, "submitted", "local_prev"), 2), f1(g(table_R3, rA, "gbif_only", "local_prev"), 2)),
  sprintf("- Design C (country-grouped folds): AUC revised %s, submitted %s, GBIF-only %s.",
          f1(g(table_R3, rC, "revised", "auc"), 3), f1(g(table_R3, rC, "submitted", "auc"), 3), f1(g(table_R3, rC, "gbif_only", "auc"), 3)),
  paste0("- Paired tests: ", paste(sprintf("%s %s: mean dS_u = %s (n = %s, p = %s)",
                                           table_R3b$design, table_R3b$comparison,
                                           f1(table_R3b$diff, 3), table_R3b$n,
                                           signif(table_R3b$p, 3)), collapse = "; "), "."), "",
  "## Spatial partitions (Comment 5)",
  paste0("- ", sprintf("%s, %s km, %s: AUC %s (SD %s), smallest test fold %d, median separation %s km",
                       table_R4$spec, table_R4$block_km, table_R4$rule, f1(table_R4$auc, 3),
                       f1(table_R4$auc_sd, 3), table_R4$min_test, f1(table_R4$sep_median_km, 0))),
  paste0("- Buffered LOO: ", paste(sprintf("%s %s km: omission %s, pooled AUC %s",
                                           table_R5$spec, table_R5$buffer_km, f1(table_R5$omission, 2),
                                           f1(table_R5$pooled_auc, 3)), collapse = "; "), "."), "",
  "## MESS of withheld records (Comment 6)",
  paste0("- ", table_R6$design, ": host area novel ", table_R6$eppo_host_area_novel_pct,
         "%; GBIF points novel ", table_R6$gbif_points_novel), "",
  "## Rice production systems (Comment 1)",
  sprintf("- Rice area in suitable climate (revised, 10th pct): %s Mha of %s Mha (%s%%); irrigated %s Mha, rainfed %s Mha. 5th pct threshold: %s%%. Across variants: %s-%s%%.",
          f1(glob$suit_Mha[glob$model == "revised_10p"], 1), f1(glob$rice_Mha[glob$model == "revised_10p"], 1),
          f1(glob$share_suitable_pct[glob$model == "revised_10p"]), f1(glob$suit_irr_Mha[glob$model == "revised_10p"], 1),
          f1(glob$suit_rf_Mha[glob$model == "revised_10p"], 1), f1(glob$share_suitable_pct[glob$model == "revised_5p"]),
          f1(100 * min(var_share)), f1(100 * max(var_share))),
  paste0("- Conducive area: ", paste(sprintf("%s = %s Mha", scen_tot$scenario, f1(scen_tot$conducive_Mha, 1)), collapse = "; "), "."),
  if (!is.null(assoc)) sprintf("- System association: median irrigated share %s at %d records vs %s in the rice landscape of the same countries (Wilcoxon p = %s).",
                               assoc$median_irr_share_records, assoc$n_records_in_rice,
                               assoc$median_irr_share_rice_landscape, assoc$wilcoxon_p) else "- System association: not run.",
  paste0("- Priority 1 survey countries: ", paste(prio$country[grepl("^1", prio$priority)], collapse = ", "), "."),
  paste0("- Priority 2 survey countries: ", paste(prio$country[grepl("^2", prio$priority)], collapse = ", "), "."),
  sprintf("- Host model AICc: cropland %s vs rice systems %s.", f1(host_tab$AICc[1]), f1(host_tab$AICc[2]))
)
writeLines(lines, p_rev("REVISION_results_digest.md"))
cat(lines, sep = "\n")

## [rev-session]
sessionInfo()



}  # end of PART 2
###############################################################################
# PART 3 - FUTURE PROJECTION (M_graminicola_future_projection.qmd) ############
###############################################################################

if (RUN_PART[["future"]]) {



# Scope and prerequisites =================================================


# What is projected, and what is held fixed -------------------------------


# The cropland layer is held constant, and this is a real limitation ------


# Clamping, and why it is on by default -----------------------------------


# Model availability varies -----------------------------------------------


# Setup ===================================================================


# Packages ----------------------------------------------------------------

## [fut-packages]
required_fut <- c(
  "terra",      # raster operations
  "geodata",    # cmip6_world(): downscaled CMIP6 bioclimatic surfaces
  "predicts",   # mess(): multivariate environmental similarity surface
  "maxnet",     # the fitted model object and its predict method
  "dplyr",
  "ggplot2",
  "knitr"       # kable() for the tables
)

missing_fut <- required_fut[
  !vapply(required_fut, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_fut)) {
  stop("Missing packages: ", paste(missing_fut, collapse = ", "),
       "\nInstall them, or restore the pinned environment with renv::restore().",
       call. = FALSE)
}

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(terra)
  library(maxnet)
})

# terra and dplyr both define extract(); be explicit throughout.
select <- dplyr::select


# Paths, and the handoff from the present-day analysis --------------------

## [fut-paths]
PROJ      <- normalizePath(".", winslash = "/", mustWork = FALSE)
DIR_DATA  <- file.path(PROJ, "data")
DIR_ENV   <- file.path(PROJ, "env_data")
DIR_OUT   <- file.path(PROJ, "outputs")
DIR_FIG   <- file.path(PROJ, "figures")
DIR_FUT   <- file.path(DIR_OUT, "future")

for (d in c(DIR_DATA, DIR_ENV, DIR_OUT, DIR_FIG, DIR_FUT)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

p_out <- function(...) file.path(DIR_OUT, ...)
p_fig <- function(...) file.path(DIR_FIG, ...)
p_fut <- function(...) file.path(DIR_FUT, ...)

## [fut-handoff]
handoff_path <- p_out("analysis_handoff.rds")
if (!file.exists(handoff_path)) {
  stop("Handoff object not found at ", handoff_path,
       "\nRender M_graminicola_SDM_analysis.qmd first: this document depends on it.",
       call. = FALSE)
}
HO <- readRDS(handoff_path)

SPECIES      <- HO$species
SPECIES_SLUG <- HO$species_slug
CFG          <- HO$cfg
PRED_VARS    <- HO$pred_vars
BIOCLIM_VARS <- HO$bioclim_vars
THR          <- HO$threshold_10p
REF_ENV      <- HO$ref_env

message(sprintf("Handoff built %s under R %s", HO$built, HO$r_version))
message(sprintf("Predictors: %s", paste(PRED_VARS, collapse = ", ")))
message(sprintf("Present-day threshold: %.6f", THR))
message(sprintf("Occurrences used: %d; study extent classified suitable: %.1f%%",
                HO$n_occ, HO$pct_suitable))

stopifnot(length(PRED_VARS) == 5L,
          identical(PRED_VARS[1:4], BIOCLIM_VARS),
          identical(PRED_VARS[5], "cropland_fraction"),
          is.finite(THR), THR > 0, THR < 1,
          all(PRED_VARS %in% names(REF_ENV)))


# Projection settings -----------------------------------------------------

## [fut-config]
FUT <- list(
  # Five models spanning different modelling centres and different equilibrium
  # climate sensitivities. Edit freely; unavailable combinations are skipped.
  gcms    = c("ACCESS-CM2", "CNRM-ESM2-1", "IPSL-CM6A-LR",
              "MIROC6", "MPI-ESM1-2-HR"),
  ssps    = c("245", "585"),
  periods = c("2041-2060", "2061-2080"),

  min_gcms = 3L,        # drop a scenario contributed to by fewer models than this

  # Which present-day model to project. "revised" uses the model written by
  # M_graminicola_revision_analyses.qmd (weights applied, polygon-based placement
  # of administrative records); "submitted" reproduces the originally submitted run.
  model_source = "revised",
  clamp    = TRUE,      # hold predictors at the training range when predicting

  # Persist the 5-layer future predictor stack for every model? About 100 MB
  # each, 20 of them. Off by default; the scenario ensemble-mean climate is
  # always kept because the novelty analysis needs it.
  keep_gcm_stacks = FALSE,

  # A cell counts as agricultural in the crop-area summaries if its present-day
  # cropland fraction is at least this value.
  cropland_min = 0.05,

  # Optional: path to a rice-specific harvested-area raster, for example a
  # MapSPAM or IRRI product. If supplied it replaces the generic cropland mask
  # in the crop-area summaries. See Section @sec-fut-country.
  rice_raster = NA_character_,

  lat_band  = 0.5,      # degrees, for the latitudinal profile
  gdal_opts = c("COMPRESS=DEFLATE", "PREDICTOR=2", "TILED=YES")
)

# BIO indices to pull out of the 19-layer CMIP6 stack, read off the training
# names: wc2.1_2.5m_bio_1 -> 1, wc2.1_2.5m_bio_4 -> 4, and so on.
BIO_INDEX <- suppressWarnings(as.integer(sub("^.*_bio_?", "", BIOCLIM_VARS)))
stopifnot(!anyNA(BIO_INDEX), all(BIO_INDEX >= 1L), all(BIO_INDEX <= 19L))
message("Bioclim indices requested from CMIP6: ",
        paste(BIO_INDEX, collapse = ", "))

# One row per scenario, ordered pathway then horizon.
scen_grid <- expand.grid(period = FUT$periods, ssp = FUT$ssps,
                         stringsAsFactors = FALSE)
scen_grid <- scen_grid[order(scen_grid$ssp, scen_grid$period), , drop = FALSE]
scen_grid$key <- paste0("ssp", scen_grid$ssp, "_", scen_grid$period)
scen_grid$label <- sprintf(
  "SSP%s-%s.%s, %s",
  substr(scen_grid$ssp, 1, 1), substr(scen_grid$ssp, 2, 2),
  substr(scen_grid$ssp, 3, 3), scen_grid$period
)
rownames(scen_grid) <- NULL
SCEN_LABEL <- stats::setNames(scen_grid$label, scen_grid$key)
print(scen_grid)

CHANGE_LEVELS <- c("Stable unsuitable", "Loss", "Gain", "Stable suitable")
CHANGE_COL <- c("Stable unsuitable" = "grey88", "Loss" = "#c0392b",
                "Gain" = "#1f6f8b", "Stable suitable" = "#4d9e5a")


# Cache helpers -----------------------------------------------------------

## [fut-cache-helpers]
write_rast <- function(r, path) {
  terra::writeRaster(r, path, overwrite = TRUE, gdal = FUT$gdal_opts)
  invisible(path)
}

cache_rast <- function(path, expr, refresh = FALSE) {
  if (file.exists(path) && !refresh) {
    message("[cache] reading ", basename(path))
    return(terra::rast(path))
  }
  message("[build] computing ", basename(path))
  r <- force(expr)
  write_rast(r, path)
  r
}


# Present-day baseline ----------------------------------------------------

## [fut-baseline]
mx_final    <- readRDS(p_out(paste0(SPECIES_SLUG, "_FINAL_model.rds")))
env_present <- terra::rast(p_out(paste0(SPECIES_SLUG, "_predictors.tif")))
names(env_present) <- PRED_VARS

suit_present <- terra::rast(p_out(paste0(SPECIES_SLUG, "_FINAL_suitability.tif")))
names(suit_present) <- "suitability"

if (identical(FUT$model_source, "revised")) {
  rho <- p_out("revision_handoff.rds")
  if (!file.exists(rho)) stop("FUT$model_source is 'revised' but ", rho,
                              " is missing. Render M_graminicola_revision_analyses.qmd first.",
                              call. = FALSE)
  RH <- readRDS(rho)
  mx_final     <- readRDS(RH$model_file)
  suit_present <- terra::rast(RH$suitability_file)
  names(suit_present) <- "suitability"
  THR     <- RH$threshold_10p
  REF_ENV <- RH$ref_env
  # The revised model is fitted with its own (weights-aware) maxnet copy; its
  # predict method is maxnet's, so nothing else changes downstream.
  message(sprintf("Projecting the REVISED model (%s, %s, fc = %s, rm = %s); threshold %.4f",
                  RH$spec$location, RH$spec$weighting, RH$fc, RH$rm, THR))
} else {
  message("Projecting the SUBMITTED model; threshold ", round(THR, 4))
}

# Recompute the binary layer from the handed-over threshold rather than reading
# the stored one, so that present and future are thresholded by demonstrably the
# same number.
bin_present <- suit_present >= THR
names(bin_present) <- "suitable"

# One reference geometry and NA mask for everything that follows. Every future
# raster is masked to this, which is what makes the area accounting in Section
# @sec-fut-change exact: gain, loss and the two stable classes then partition one
# common set of cells.
TEMPLATE  <- suit_present
CELL_AREA <- terra::cellSize(TEMPLATE, unit = "km", mask = TRUE)
names(CELL_AREA) <- "area_km2"

# Cell-wise area sum under a 0/1 mask. Every area figure in this document goes
# through this one function, so they are all computed the same way.
area_of <- function(mask_layer) {
  as.numeric(terra::global(CELL_AREA * mask_layer, "sum", na.rm = TRUE)[1, 1])
}

CROPLAND_CONST <- env_present[["cropland_fraction"]]

print(env_present)
message(sprintf("Grid: %d by %d cells, resolution %.7f degrees",
                terra::ncol(TEMPLATE), terra::nrow(TEMPLATE),
                terra::res(TEMPLATE)[1]))

area_total   <- as.numeric(terra::global(CELL_AREA, "sum", na.rm = TRUE)[1, 1])
area_present <- area_of(bin_present)
message(sprintf("Present-day suitable area: %s km2 (%.1f%% of the %s km2 of land in extent)",
                format(round(area_present), big.mark = ","),
                100 * area_present / area_total,
                format(round(area_total), big.mark = ",")))

# Guard against a stale handoff paired with a re-run analysis.
stopifnot(terra::nlyr(env_present) == length(PRED_VARS))
if (!isTRUE(all.equal(as.vector(terra::ext(env_present)),
                      HO$study_ext, tolerance = 1e-6))) {
  warning("Predictor stack extent does not match the handoff. ",
          "Re-render M_graminicola_SDM_analysis.qmd so the two agree.")
}


# Assembling the future predictor stacks ==================================

## [fut-stack-builder]
# The 4 requested bioclim layers for one model-pathway-horizon, renamed to the
# training names and snapped to the baseline grid. NULL if unavailable.
future_climate <- function(gcm, ssp, period) {
  bioc <- tryCatch(
    geodata::cmip6_world(model = gcm, ssp = ssp, time = period,
                         var = "bioc", res = CFG$bioclim_res, path = DIR_ENV),
    error = function(e) {
      message("  [skip] ", gcm, " ssp", ssp, " ", period, ": ",
              conditionMessage(e))
      NULL
    }
  )
  if (is.null(bioc)) return(NULL)

  if (terra::nlyr(bioc) != 19L) {
    warning("Expected 19 bioclim layers for ", gcm, " ssp", ssp, " ", period,
            ", got ", terra::nlyr(bioc), ". Skipping.")
    return(NULL)
  }

  # Match by trailing index rather than trusting layer order.
  idx <- suppressWarnings(
    as.integer(sub("^.*[^0-9]([0-9]+)$", "\\1", names(bioc)))
  )
  if (anyNA(idx) || !identical(sort(idx), 1:19)) {
    warning("Layer names for ", gcm, " ssp", ssp, " ", period,
            " carry no usable bioclim index; falling back to layer order.")
    idx <- seq_len(19L)
  }
  pick <- match(BIO_INDEX, idx)
  stopifnot(!anyNA(pick))

  sel <- bioc[[pick]]
  names(sel) <- BIOCLIM_VARS

  sel <- terra::crop(sel, terra::ext(TEMPLATE))
  if (!isTRUE(terra::compareGeom(sel, TEMPLATE, stopOnError = FALSE,
                                 messages = FALSE))) {
    message("  [align] resampling ", gcm, " ssp", ssp, " ", period,
            " onto the baseline grid")
    sel <- terra::resample(sel, TEMPLATE, method = "bilinear")
  }
  sel
}

# The full 5-predictor stack the model expects.
future_stack <- function(gcm, ssp, period) {
  clim <- future_climate(gcm, ssp, period)
  if (is.null(clim)) return(NULL)
  s <- c(clim, CROPLAND_CONST)          # land use held constant, by design
  names(s) <- PRED_VARS                 # exact training names, exact order
  s <- terra::mask(s, TEMPLATE)         # one common NA footprint
  stopifnot(identical(names(s), PRED_VARS))
  s
}


# Projecting the model ====================================================

## [fut-project]
proj_log <- vector("list", nrow(scen_grid))

for (i in seq_len(nrow(scen_grid))) {
  key    <- scen_grid$key[i]
  ssp    <- scen_grid$ssp[i]
  period <- scen_grid$period[i]
  message("=== ", SCEN_LABEL[[key]], " ===")

  clim_path <- p_fut(paste0("climate_", key, ".tif"))
  clim_have <- file.exists(clim_path)

  clim_sum <- NULL
  n_sum    <- 0L
  used     <- character(0)

  for (gcm in FUT$gcms) {
    tag       <- paste0(gcm, "_", key)
    suit_path <- p_fut(paste0(SPECIES_SLUG, "_suitability_", tag, ".tif"))

    # Nothing to build: both the projection and the scenario climate are cached,
    # so the future stack is never assembled and no download is triggered.
    if (clim_have && file.exists(suit_path)) {
      message("  [cache] ", basename(suit_path))
      used <- c(used, gcm)
      next
    }

    fs <- future_stack(gcm, ssp, period)
    if (is.null(fs)) next

    if (isTRUE(FUT$keep_gcm_stacks)) {
      write_rast(fs, p_fut(paste0("predictors_", tag, ".tif")))
    }

    if (file.exists(suit_path)) {
      message("  [cache] ", basename(suit_path))
    } else {
      message("  [build] ", basename(suit_path))
      sp <- terra::predict(fs, mx_final, type = "cloglog",
                           clamp = FUT$clamp, na.rm = TRUE)
      names(sp) <- "suitability"
      write_rast(sp, suit_path)
      rm(sp)
    }

    if (!clim_have) {
      clim_sum <- if (is.null(clim_sum)) fs[[BIOCLIM_VARS]] else
        clim_sum + fs[[BIOCLIM_VARS]]
      n_sum <- n_sum + 1L
    }
    used <- c(used, gcm)
    rm(fs)
  }

  # Ensemble-mean climate for the scenario, needed by the novelty analysis. The
  # divisor is the count actually accumulated, not the count of models used.
  if (!clim_have && n_sum > 0L) {
    cm <- clim_sum / n_sum
    names(cm) <- BIOCLIM_VARS
    write_rast(cm, clim_path)
    rm(cm)
  }
  rm(clim_sum)

  proj_log[[i]] <- data.frame(
    scenario = key,
    label    = unname(SCEN_LABEL[[key]]),
    n_gcm    = length(used),
    gcms     = paste(used, collapse = "; "),
    stringsAsFactors = FALSE
  )
  gc()
}

projection_log <- dplyr::bind_rows(proj_log)
utils::write.csv(projection_log, p_out("Table_S2_projection_log.csv"),
                 row.names = FALSE)
knitr::kable(projection_log[, c("label", "n_gcm", "gcms")],
             caption = "Table S2. Models contributing to each scenario.")

## [fut-scenario-screen]
kept    <- projection_log$scenario[projection_log$n_gcm >= FUT$min_gcms]
dropped <- setdiff(projection_log$scenario, kept)

if (length(dropped)) {
  warning("Scenarios dropped for having fewer than ", FUT$min_gcms,
          " contributing models: ", paste(dropped, collapse = ", "))
}
if (!length(kept)) {
  stop("No scenario met the minimum model count. Nothing downloaded successfully; ",
       "check network access and the model names in FUT$gcms.", call. = FALSE)
}
SCEN <- kept
message("Scenarios carried forward: ", paste(SCEN, collapse = ", "))


# Ensemble mean and model agreement =======================================

## [fut-ensemble]
ens_mean  <- list()
ens_agree <- list()
ens_n     <- integer(0)

for (key in SCEN) {
  # Restrict to the models currently configured, by construction rather than by
  # pattern matching, so that suitability rasters left behind by a run with a
  # different FUT$gcms cannot silently enter the ensemble.
  wanted <- p_fut(paste0(SPECIES_SLUG, "_suitability_", FUT$gcms, "_", key, ".tif"))
  files  <- wanted[file.exists(wanted)]
  stopifnot(length(files) >= FUT$min_gcms)

  s <- terra::rast(files)
  ens_n[key] <- terra::nlyr(s)

  em <- cache_rast(p_fut(paste0(SPECIES_SLUG, "_ensmean_", key, ".tif")), {
    terra::mask(terra::app(s, fun = "mean", na.rm = TRUE), TEMPLATE)
  })
  names(em) <- "suitability"
  ens_mean[[key]] <- em

  ea <- cache_rast(p_fut(paste0(SPECIES_SLUG, "_agreement_", key, ".tif")), {
    terra::mask(terra::app(s >= THR, fun = "sum", na.rm = TRUE), TEMPLATE)
  })
  names(ea) <- "n_models_suitable"
  ens_agree[[key]] <- ea

  message(sprintf("%s: ensemble from %d models", SCEN_LABEL[[key]], ens_n[[key]]))
  rm(s)
}

ens_bin <- lapply(ens_mean, function(r) {
  b <- r >= THR
  names(b) <- "suitable"
  b
})

## [fut-agreement-summary]
agree_tab <- dplyr::bind_rows(lapply(SCEN, function(key) {
  n     <- ens_n[[key]]
  any_m <- area_of(ens_agree[[key]] >= 1L)
  maj   <- area_of(ens_agree[[key]] >= ceiling(n / 2))
  full  <- area_of(ens_agree[[key]] == n)
  data.frame(
    Scenario = unname(SCEN_LABEL[[key]]),
    Models   = n,
    `Suitable in at least one model (km2)` = round(any_m),
    `Suitable in a majority (km2)`         = round(maj),
    `Suitable in all models (km2)`         = round(full),
    `Unanimous share of any-model area (%)` = round(100 * full / any_m, 1),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))

utils::write.csv(agree_tab, p_out("Table_9_model_agreement.csv"), row.names = FALSE)
knitr::kable(agree_tab, format.args = list(big.mark = ","),
             caption = "Table 9. Model agreement on suitability, by scenario. The unanimous share indicates how much of the projected range every contributing model supports.")


# Range change and area statistics ========================================

## [fut-change-raster]
change_rast <- list()
for (key in SCEN) {
  cr <- cache_rast(p_fut(paste0(SPECIES_SLUG, "_change_", key, ".tif")), {
    r <- 1 * bin_present + 2 * ens_bin[[key]]
    names(r) <- "change"
    r
  })
  names(cr) <- "change"
  change_rast[[key]] <- cr
}

# The four classes must partition the land area exactly. Check it rather than
# assume it: a mismatch would mean the NA masks came apart somewhere upstream.
for (key in SCEN) {
  parts <- vapply(0:3, function(k) area_of(change_rast[[key]] == k), numeric(1))
  if (!isTRUE(all.equal(sum(parts), area_total, tolerance = 1e-6))) {
    warning(SCEN_LABEL[[key]], ": change classes sum to ",
            format(round(sum(parts)), big.mark = ","), " km2 but the land area is ",
            format(round(area_total), big.mark = ","), " km2.")
  }
}

## [fut-change-areas]
change_areas <- dplyr::bind_rows(lapply(SCEN, function(key) {
  cr <- change_rast[[key]]
  a  <- vapply(0:3, function(k) area_of(cr == k), numeric(1))
  data.frame(scenario = key, label = unname(SCEN_LABEL[[key]]),
             stable_unsuitable_km2 = a[1], loss_km2 = a[2],
             gain_km2 = a[3], stable_suitable_km2 = a[4],
             stringsAsFactors = FALSE)
}))

change_areas <- change_areas |>
  mutate(present_km2 = stable_suitable_km2 + loss_km2,
         future_km2  = stable_suitable_km2 + gain_km2,
         net_km2     = future_km2 - present_km2,
         net_pct     = 100 * net_km2 / present_km2,
         stable_pct_of_present = 100 * stable_suitable_km2 / present_km2,
         gain_pct_of_present = 100 * gain_km2 / present_km2,
         loss_pct_of_present = 100 * loss_km2 / present_km2,
         turnover_pct        = 100 * (gain_km2 + loss_km2) / future_km2)

# The present-day column is recovered from the change classes, so it is a second
# independent route to the same number computed in the baseline chunk.
stopifnot(all(abs(change_areas$present_km2 - area_present) < 1))

table_7 <- data.frame(
  Scenario                     = change_areas$label,
  `Present suitable (km2)`     = round(change_areas$present_km2),
  `Projected suitable (km2)`   = round(change_areas$future_km2),
  `Net change (km2)`           = round(change_areas$net_km2),
  `Stable (% of present)`      = round(change_areas$stable_pct_of_present, 1),
  `Net change (%)`             = round(change_areas$net_pct, 1),
  `Gain (% of present)`        = round(change_areas$gain_pct_of_present, 1),
  `Loss (% of present)`        = round(change_areas$loss_pct_of_present, 1),
  `Range turnover (%)`         = round(change_areas$turnover_pct, 1),
  check.names = FALSE, stringsAsFactors = FALSE
)
utils::write.csv(table_7, p_out("Table_7_projected_area_change.csv"),
                 row.names = FALSE)
knitr::kable(table_7, format.args = list(big.mark = ","),
             caption = "Table 7. Projected change in climatically suitable area, ensemble mean per scenario, thresholded at the present-day 10th-percentile training-presence value. All percentages except turnover use present-day suitable area as the denominator: stable + loss = 100%, and net change = gain - loss. Range turnover is gain plus loss as a percentage of projected suitable area. Areas are sums of true cell areas (terra::cellSize), not cell counts.")

## [fut-change-by-gcm]
# The same accounting for every individual model, so that the spread behind each
# ensemble mean is visible (one row per model and scenario).
gcm_change <- dplyr::bind_rows(lapply(SCEN, function(key) {
  dplyr::bind_rows(lapply(FUT$gcms, function(gcm) {
    f <- p_fut(paste0(SPECIES_SLUG, "_suitability_", gcm, "_", key, ".tif"))
    if (!file.exists(f)) return(NULL)
    fb <- terra::mask(terra::rast(f) >= THR, TEMPLATE)
    st <- area_of(bin_present * fb); ls <- area_of(bin_present * (1 - fb))
    gn <- area_of((1 - bin_present) * fb)
    data.frame(Scenario = unname(SCEN_LABEL[[key]]), GCM = gcm,
               `Stable (% of present)` = round(100 * st / area_present, 1),
               `Loss (% of present)`   = round(100 * ls / area_present, 1),
               `Gain (% of present)`   = round(100 * gn / area_present, 1),
               `Net change (% of present)` = round(100 * (gn - ls) / area_present, 1),
               check.names = FALSE)
  }))
}))
utils::write.csv(gcm_change, p_out("Table_7b_change_by_gcm.csv"), row.names = FALSE)
knitr::kable(gcm_change, caption = "Table 7b. Range change by individual model. Denominator: present-day suitable area.")

## [fut-change-classes-table]
table_8 <- dplyr::bind_rows(lapply(seq_len(nrow(change_areas)), function(i) {
  r <- change_areas[i, ]
  data.frame(
    Scenario = rep(r$label, 4),
    Class    = CHANGE_LEVELS,
    `Area (km2)` = round(c(r$stable_unsuitable_km2, r$loss_km2,
                           r$gain_km2, r$stable_suitable_km2)),
    `Share of land in extent (%)` = round(
      100 * c(r$stable_unsuitable_km2, r$loss_km2,
              r$gain_km2, r$stable_suitable_km2) / area_total, 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))
utils::write.csv(table_8, p_out("Table_8_range_change_classes.csv"),
                 row.names = FALSE)
knitr::kable(table_8, format.args = list(big.mark = ","),
             caption = "Table 8. Range change classes by scenario, in absolute area and as a share of the land within the study extent.")


# Future extrapolation and novelty ========================================

## [fut-mess]
mess_rast <- list()
for (key in SCEN) {
  clim <- terra::rast(p_fut(paste0("climate_", key, ".tif")))
  names(clim) <- BIOCLIM_VARS

  fstack <- c(clim, CROPLAND_CONST)
  names(fstack) <- PRED_VARS
  fstack <- terra::mask(fstack, TEMPLATE)

  mr <- cache_rast(p_fut(paste0(SPECIES_SLUG, "_MESS_", key, ".tif")), {
    terra::mask(predicts::mess(fstack, REF_ENV[, PRED_VARS, drop = FALSE]),
                TEMPLATE)
  })
  names(mr) <- "mess"
  mess_rast[[key]] <- mr
  rm(clim, fstack)
}

mess_present <- terra::rast(p_out(paste0(SPECIES_SLUG, "_MESS.tif")))
names(mess_present) <- "mess"

## [fut-mess-table]
novelty_row <- function(label, mess_r, bin_r) {
  novel <- mess_r < 0
  data.frame(
    Scenario                              = label,
    `Novel area in extent (km2)`          = round(area_of(novel)),
    `Novel share of extent (%)`           = round(100 * area_of(novel) / area_total, 1),
    `Suitable area that is novel (km2)`   = round(area_of(novel * bin_r)),
    `Novel share of suitable area (%)`    = round(
      100 * area_of(novel * bin_r) / area_of(bin_r), 1),
    `Minimum MESS`                        = round(
      as.numeric(terra::global(mess_r, "min", na.rm = TRUE)[1, 1]), 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}

table_10 <- dplyr::bind_rows(
  novelty_row("Present day (baseline)", mess_present, bin_present),
  dplyr::bind_rows(lapply(SCEN, function(key)
    novelty_row(unname(SCEN_LABEL[[key]]), mess_rast[[key]], ens_bin[[key]])))
)
utils::write.csv(table_10, p_out("Table_10_future_novelty.csv"), row.names = FALSE)
knitr::kable(table_10, format.args = list(big.mark = ","),
             caption = "Table 10. Multivariate environmental similarity of the ensemble-mean climate to the training envelope. Negative MESS marks cells where at least one predictor falls outside the range spanned by the training presences.")

## [fut-mess-gain-overlap]
# The sharpest version of the question: how much of the *newly gained* area is
# climatically novel? A gain that is largely novel is a projection into
# conditions the model has never seen, and should be reported as a hypothesis
# rather than a prediction.
gain_novelty <- dplyr::bind_rows(lapply(SCEN, function(key) {
  gain  <- change_rast[[key]] == 2
  novel <- mess_rast[[key]] < 0
  g  <- area_of(gain)
  gn <- area_of(gain * novel)
  data.frame(
    Scenario                       = unname(SCEN_LABEL[[key]]),
    `Gained area (km2)`            = round(g),
    `Gained and novel (km2)`       = round(gn),
    `Novel share of gain (%)`      = round(100 * gn / g, 1),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))
utils::write.csv(gain_novelty, p_out("Table_10b_gain_novelty.csv"),
                 row.names = FALSE)
knitr::kable(gain_novelty, format.args = list(big.mark = ","),
             caption = "Table 10b. How much of the projected gain requires extrapolation beyond the training envelope.")


# Latitudinal and centroid shift ==========================================


# Both hemispheres are in the extent --------------------------------------

## [fut-shift-helpers]
LON_R <- terra::init(TEMPLATE, "x")
LAT_R <- terra::init(TEMPLATE, "y")
names(LON_R) <- "lon"
names(LAT_R) <- "lat"

LAT_BAND <- round(LAT_R / FUT$lat_band) * FUT$lat_band
names(LAT_BAND) <- "lat_band"

# Area-weighted centroid of the suitable cells.
suit_centroid <- function(bin) {
  w  <- CELL_AREA * bin
  sw <- as.numeric(terra::global(w, "sum", na.rm = TRUE)[1, 1])
  if (!is.finite(sw) || sw <= 0) {
    return(c(lon = NA_real_, lat = NA_real_, area_km2 = 0))
  }
  c(lon = as.numeric(terra::global(w * LON_R, "sum", na.rm = TRUE)[1, 1]) / sw,
    lat = as.numeric(terra::global(w * LAT_R, "sum", na.rm = TRUE)[1, 1]) / sw,
    area_km2 = sw)
}

# Suitable area per latitude band.
lat_profile <- function(bin) {
  w   <- CELL_AREA * bin
  tab <- terra::zonal(w, LAT_BAND, fun = "sum", na.rm = TRUE)
  names(tab) <- c("lat_band", "area_km2")
  tab <- tab[is.finite(tab$lat_band) & is.finite(tab$area_km2), , drop = FALSE]
  tab[order(tab$lat_band), , drop = FALSE]
}

# Area-weighted quantiles of a latitude profile.
wquant <- function(x, w, probs) {
  keep <- is.finite(x) & is.finite(w) & w > 0
  x <- x[keep]; w <- w[keep]
  if (!length(x)) return(rep(NA_real_, length(probs)))
  o  <- order(x); x <- x[o]; w <- w[o]
  cw <- cumsum(w) / sum(w)
  vapply(probs, function(p) x[which(cw >= p)[1]], numeric(1))
}

# Haversine distance, so no extra package dependency is introduced for one line.
gc_dist_km <- function(lon1, lat1, lon2, lat2) {
  R   <- 6371.0088
  d2r <- pi / 180
  a <- sin((lat2 - lat1) * d2r / 2)^2 +
    cos(lat1 * d2r) * cos(lat2 * d2r) * sin((lon2 - lon1) * d2r / 2)^2
  2 * R * asin(pmin(1, sqrt(a)))
}

## [fut-shift-compute]
# Bands holding less than this are treated as noise when locating a range margin,
# so that a handful of isolated cells cannot set the reported limit.
LIMIT_MIN_AREA <- 500

margin_limits <- function(prof) {
  solid <- prof[prof$area_km2 >= LIMIT_MIN_AREA, , drop = FALSE]
  if (!nrow(solid)) return(c(north = NA_real_, south = NA_real_))
  c(north = max(solid$lat_band), south = min(solid$lat_band))
}

profiles <- c(
  list(Present = lat_profile(bin_present)),
  stats::setNames(lapply(SCEN, function(key) lat_profile(ens_bin[[key]])), SCEN)
)

cent_present <- suit_centroid(bin_present)
lim_present  <- margin_limits(profiles$Present)
q_present    <- wquant(profiles$Present$lat_band, profiles$Present$area_km2,
                       c(0.05, 0.5, 0.95))

shift_tab <- dplyr::bind_rows(lapply(SCEN, function(key) {
  cent <- suit_centroid(ens_bin[[key]])
  prof <- profiles[[key]]
  lim  <- margin_limits(prof)
  q    <- wquant(prof$lat_band, prof$area_km2, c(0.05, 0.5, 0.95))
  data.frame(
    Scenario              = unname(SCEN_LABEL[[key]]),
    `Centroid lon`        = round(cent[["lon"]], 3),
    `Centroid lat`        = round(cent[["lat"]], 3),
    `Northward (deg)`     = round(cent[["lat"]] - cent_present[["lat"]], 3),
    `Eastward (deg)`      = round(cent[["lon"]] - cent_present[["lon"]], 3),
    `Displacement (km)`   = round(
      gc_dist_km(cent_present[["lon"]], cent_present[["lat"]],
                 cent[["lon"]], cent[["lat"]]), 1),
    `Median lat`          = q[2],
    `5th pct lat`         = q[1],
    `95th pct lat`        = q[3],
    `North margin`        = lim[["north"]],
    `North shift (deg)`   = round(lim[["north"]] - lim_present[["north"]], 2),
    `South margin`        = lim[["south"]],
    `South shift (deg)`   = round(lim_present[["south"]] - lim[["south"]], 2),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))

table_11 <- dplyr::bind_rows(
  data.frame(
    Scenario              = "Present day (baseline)",
    `Centroid lon`        = round(cent_present[["lon"]], 3),
    `Centroid lat`        = round(cent_present[["lat"]], 3),
    `Northward (deg)`     = 0,
    `Eastward (deg)`      = 0,
    `Displacement (km)`   = 0,
    `Median lat`          = q_present[2],
    `5th pct lat`         = q_present[1],
    `95th pct lat`        = q_present[3],
    `North margin`        = lim_present[["north"]],
    `North shift (deg)`   = 0,
    `South margin`        = lim_present[["south"]],
    `South shift (deg)`   = 0,
    check.names = FALSE, stringsAsFactors = FALSE
  ),
  shift_tab
)

utils::write.csv(table_11, p_out("Table_11_centroid_and_latitude_shift.csv"),
                 row.names = FALSE)
knitr::kable(table_11,
             caption = "Table 11. Area-weighted centroid of suitable habitat and the latitudinal structure of the range. Margin shifts are signed so that a positive value means the margin moves away from the equator; margins are the outermost 0.5 degree bands holding at least 500 km2 of suitable area.")

## [fut-profile-export]
profile_long <- dplyr::bind_rows(lapply(names(profiles), function(nm) {
  p <- profiles[[nm]]
  data.frame(series = if (nm == "Present") "Present day" else unname(SCEN_LABEL[[nm]]),
             lat_band = p$lat_band, area_km2 = p$area_km2,
             stringsAsFactors = FALSE)
}))
utils::write.csv(profile_long, p_out("future_latitudinal_profile.csv"),
                 row.names = FALSE)
message("Latitude bands with suitable area: ",
        sum(profile_long$area_km2 > 0), " rows across ",
        length(unique(profile_long$series)), " series.")


# Country and crop-area summaries =========================================


# These are cropland summaries, not rice summaries ------------------------

## [fut-crop-mask]
if (!is.na(FUT$rice_raster) && nzchar(FUT$rice_raster) &&
    file.exists(FUT$rice_raster)) {
  rr <- terra::rast(FUT$rice_raster)[[1]]
  if (!isTRUE(terra::compareGeom(rr, TEMPLATE, stopOnError = FALSE,
                                 messages = FALSE))) {
    message("[align] resampling the supplied rice raster onto the analysis grid")
    rr <- terra::resample(rr, TEMPLATE, method = "bilinear")
  }
  CROP_MASK  <- terra::mask(rr > 0, TEMPLATE)
  CROP_LABEL <- "rice"
  message("Crop mask: user-supplied rice raster ", basename(FUT$rice_raster))
} else {
  CROP_MASK  <- terra::mask(CROPLAND_CONST >= FUT$cropland_min, TEMPLATE)
  CROP_LABEL <- "cropland"
  message(sprintf("Crop mask: ESA WorldCover cropland fraction >= %.2f",
                  FUT$cropland_min))
}
names(CROP_MASK) <- "crop"

crop_area_total   <- area_of(CROP_MASK)
crop_suit_present <- area_of(bin_present * CROP_MASK)
message(sprintf("Area under the crop mask: %s km2 (%.1f%% of the extent)",
                format(round(crop_area_total), big.mark = ","),
                100 * crop_area_total / area_total))
message(sprintf("Of which suitable today: %s km2 (%.1f%%)",
                format(round(crop_suit_present), big.mark = ","),
                100 * crop_suit_present / crop_area_total))

## [fut-country-zones]
wv <- geodata::world(path = DIR_ENV)
wv <- terra::crop(wv, terra::ext(TEMPLATE))

wv_fields <- names(wv)
name_field <- if ("NAME_0" %in% wv_fields) "NAME_0" else wv_fields[1]
message("Country name field: ", name_field)

wv$zone <- seq_len(nrow(wv))
wv_df   <- as.data.frame(wv)
country_lut <- data.frame(zone = wv_df$zone,
                          country = as.character(wv_df[[name_field]]),
                          stringsAsFactors = FALSE)
utils::write.csv(country_lut, p_out("future_country_lookup.csv"), row.names = FALSE)

country_r <- cache_rast(p_fut("country_zones.tif"), {
  terra::mask(terra::rasterize(wv, TEMPLATE, field = "zone"), TEMPLATE)
})
names(country_r) <- "zone"
message("Countries intersecting the study extent: ", nrow(country_lut))

## [fut-country-summary]
# One zonal pass per scenario over an 8-layer area stack, rather than eight
# passes over one layer each. On this grid that is the difference between a
# couple of minutes and most of an hour.
country_by_scen <- dplyr::bind_rows(lapply(SCEN, function(key) {
  cr   <- change_rast[[key]]
  fbin <- ens_bin[[key]]

  lyrs <- c(CELL_AREA * bin_present,
            CELL_AREA * fbin,
            CELL_AREA * (cr == 2),
            CELL_AREA * (cr == 1),
            CELL_AREA * bin_present * CROP_MASK,
            CELL_AREA * fbin * CROP_MASK,
            CELL_AREA * (cr == 2) * CROP_MASK,
            CELL_AREA * (cr == 1) * CROP_MASK)
  names(lyrs) <- c("present_km2", "future_km2", "gain_km2", "loss_km2",
                   "present_crop_km2", "future_crop_km2",
                   "gain_crop_km2", "loss_crop_km2")

  z <- terra::zonal(lyrs, country_r, fun = "sum", na.rm = TRUE)
  names(z)[1] <- "zone"
  z <- merge(country_lut, z, by = "zone", all.y = TRUE)
  z$scenario <- key
  z$label    <- unname(SCEN_LABEL[[key]])
  z$net_km2  <- z$future_km2 - z$present_km2
  z$net_pct  <- ifelse(z$present_km2 > 0, 100 * z$net_km2 / z$present_km2, NA_real_)
  z
}))

utils::write.csv(country_by_scen, p_out("Table_S3_country_change_full.csv"),
                 row.names = FALSE)
message("Country by scenario rows written: ", nrow(country_by_scen))

## [fut-country-table]
# Headline table: the countries with the largest absolute change under the most
# severe scenario, shown across all scenarios so the gradient is visible.
focus_key <- SCEN[length(SCEN)]
focus <- country_by_scen |>
  filter(scenario == focus_key) |>
  arrange(desc(abs(net_km2)))
top_countries <- utils::head(focus$country[!is.na(focus$country)], 15)

table_12 <- country_by_scen |>
  filter(country %in% top_countries) |>
  mutate(country = factor(country, levels = top_countries)) |>
  arrange(country, scenario) |>
  transmute(Country              = as.character(country),
            Scenario             = label,
            `Present (km2)`      = round(present_km2),
            `Projected (km2)`    = round(future_km2),
            `Gain (km2)`         = round(gain_km2),
            `Loss (km2)`         = round(loss_km2),
            `Net (km2)`          = round(net_km2),
            `Net (%)`            = round(net_pct, 1),
            `Projected on crop land (km2)` = round(future_crop_km2),
            `Crop-land gain (km2)`         = round(gain_crop_km2))

utils::write.csv(table_12, p_out("Table_12_country_change_top.csv"),
                 row.names = FALSE)
knitr::kable(table_12, format.args = list(big.mark = ","),
             caption = paste0("Table 12. Projected change in suitable area for the fifteen countries with the largest absolute change under ",
                              SCEN_LABEL[[focus_key]],
                              ", shown across all scenarios. Crop-land columns are restricted to the ",
                              CROP_LABEL, " mask."))

## [fut-crop-table]
table_13 <- dplyr::bind_rows(lapply(SCEN, function(key) {
  cr   <- change_rast[[key]]
  fbin <- ens_bin[[key]]
  fc <- area_of(fbin * CROP_MASK)
  data.frame(
    Scenario                                  = unname(SCEN_LABEL[[key]]),
    `Present suitable crop land (km2)`         = round(crop_suit_present),
    `Projected suitable crop land (km2)`       = round(fc),
    `Net change (km2)`                         = round(fc - crop_suit_present),
    `Net change (%)`                           = round(
      100 * (fc - crop_suit_present) / crop_suit_present, 1),
    `Gain on crop land (km2)`                  = round(area_of((cr == 2) * CROP_MASK)),
    `Loss on crop land (km2)`                  = round(area_of((cr == 1) * CROP_MASK)),
    `Share of crop land suitable (%)`          = round(100 * fc / crop_area_total, 1),
    check.names = FALSE, stringsAsFactors = FALSE
  )
}))

utils::write.csv(table_13, p_out("Table_13_crop_area_change.csv"), row.names = FALSE)
knitr::kable(table_13, format.args = list(big.mark = ","),
             caption = paste0("Table 13. Projected change in suitable area restricted to the ",
                              CROP_LABEL, " mask. Because the ", CROP_LABEL,
                              " layer is held constant, every change here is driven by climate alone."))


# Rice area by production system ------------------------------------------

## [fut-rice-systems]
rice_f <- file.path(DIR_OUT, "revision", "rice_phys_area_systems_2p5min.tif")
if (file.exists(rice_f)) {
  RICE <- terra::rast(rice_f)
  rice_tab <- dplyr::bind_rows(lapply(SCEN, function(key) {
    cr <- change_rast[[key]]
    dplyr::bind_rows(lapply(c("rice_irrigated_ha", "rice_rainfed_ha"), function(v) {
      r <- RICE[[v]]
      a <- vapply(c(3, 1, 2), function(k)
        as.numeric(terra::global(r * (cr == k), "sum", na.rm = TRUE)[1, 1]), numeric(1))
      pres <- a[1] + a[2]
      data.frame(Scenario = unname(SCEN_LABEL[[key]]),
                 System = sub("rice_(.*)_ha", "\\1", v),
                 `Suitable now (Mha)` = round(pres / 1e6, 2),
                 `Stable (Mha)` = round(a[1] / 1e6, 2), `Lost (Mha)` = round(a[2] / 1e6, 2),
                 `Gained (Mha)` = round(a[3] / 1e6, 2),
                 `Net change (% of suitable now)` = round(100 * (a[3] - a[2]) / pres, 1),
                 check.names = FALSE)
    }))
  }))
  utils::write.csv(rice_tab, p_out("Table_13b_rice_system_change.csv"), row.names = FALSE)
  knitr::kable(rice_tab, caption = "Table 13b. Projected change in climatically suitable rice area, by production system (MapSPAM 2010 physical area).")
} else {
  message("Revision rice layers not found; rice-system table skipped.")
}


# Figures =================================================================

## [fut-fig-theme]
theme_sdm <- function(base = 11) {
  theme_minimal(base_size = base) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(linewidth = 0.2, colour = "grey88"),
      plot.title    = element_text(face = "bold", size = base + 2),
      plot.subtitle = element_text(colour = "grey30", size = base - 1),
      legend.position = "bottom",
      legend.key.width = grid::unit(2.0, "cm")
    )
}

world_df <- suppressWarnings(as.data.frame(terra::geom(wv)))

XLIM <- c(terra::xmin(TEMPLATE), terra::xmax(TEMPLATE))
YLIM <- c(terra::ymin(TEMPLATE), terra::ymax(TEMPLATE))

save_fig <- function(plot, file, width = 12, height = 8) {
  ggsave(p_fig(file), plot = plot, width = width, height = height,
         dpi = 300, bg = "white")
  message("[figure] ", file)
  invisible(p_fig(file))
}

# Long data frame from a multi-layer raster, downsampled for plotting. `fun`
# must be "modal" for categorical layers: averaging class codes would invent
# classes that do not exist.
rast_df_multi <- function(r, max_cells = 3.5e5, fun = "mean") {
  rr <- if (terra::ncell(r) > max_cells) {
    terra::aggregate(r, fact = ceiling(sqrt(terra::ncell(r) / max_cells)),
                     fun = fun, na.rm = TRUE)
  } else r
  d <- terra::as.data.frame(rr, xy = TRUE, na.rm = FALSE)
  names(d) <- c("x", "y", names(rr))
  out <- do.call(rbind, lapply(names(rr), function(nm) {
    data.frame(x = d$x, y = d$y, layer = nm, value = d[[nm]],
               stringsAsFactors = FALSE)
  }))
  out[!is.na(out$value), , drop = FALSE]
}

# Combine a per-scenario list of single-layer rasters into one stack whose layer
# names are the scenario keys, then attach ordered display labels.
stack_scen <- function(lst) {
  s <- terra::rast(unname(lst[SCEN]))
  names(s) <- SCEN
  s
}
add_scen_label <- function(df) {
  df$label <- factor(unname(SCEN_LABEL[df$layer]),
                     levels = unname(SCEN_LABEL[SCEN]))
  df
}

occ_pts <- HO$occ_coords

## [figure-7]
df7 <- add_scen_label(rast_df_multi(stack_scen(ens_mean)))

fig7 <- ggplot() +
  geom_raster(data = df7, aes(x, y, fill = value)) +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = NA, colour = "grey25", linewidth = 0.15) +
  geom_point(data = occ_pts, aes(lon, lat), colour = "red",
             size = 0.35, alpha = 0.7) +
  scale_fill_viridis_c(option = "D", name = "Suitability (cloglog)",
                       limits = c(0, 1)) +
  coord_equal(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  facet_wrap(~ label, ncol = 2) +
  labs(title = "Projected suitability under CMIP6 scenarios",
       subtitle = sprintf("Ensemble mean of up to %d general circulation models | cropland held at present day | occurrences in red",
                          max(ens_n)),
       x = "Longitude", y = "Latitude") +
  theme_sdm()

save_fig(fig7, "Figure_7_future_suitability.png", width = 13, height = 8.5)
fig7

## [figure-8]
df8 <- add_scen_label(rast_df_multi(stack_scen(change_rast), fun = "modal"))
df8$class <- factor(CHANGE_LEVELS[df8$value + 1L], levels = CHANGE_LEVELS)

fig8 <- ggplot() +
  geom_raster(data = df8, aes(x, y, fill = class)) +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = NA, colour = "grey35", linewidth = 0.15) +
  scale_fill_manual(values = CHANGE_COL, name = NULL, drop = FALSE) +
  coord_equal(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  facet_wrap(~ label, ncol = 2) +
  labs(title = "Projected range change",
       subtitle = sprintf("Present and future both thresholded at the present-day 10th-percentile value (%.3f)",
                          THR),
       x = "Longitude", y = "Latitude") +
  theme_sdm() +
  guides(fill = guide_legend(override.aes = list(size = 4), nrow = 1))

save_fig(fig8, "Figure_8_range_change.png", width = 13, height = 8.5)
fig8

## [figure-9]
df9 <- add_scen_label(rast_df_multi(stack_scen(ens_agree), fun = "modal"))
df9$n <- factor(as.integer(round(df9$value)),
                levels = 0:max(ens_n))

fig9 <- ggplot() +
  geom_raster(data = df9, aes(x, y, fill = n)) +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = NA, colour = "grey35", linewidth = 0.15) +
  scale_fill_viridis_d(option = "C", name = "Models projecting suitable",
                       drop = FALSE) +
  coord_equal(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  facet_wrap(~ label, ncol = 2) +
  labs(title = "Agreement among general circulation models",
       subtitle = "Number of models placing the cell at or above the present-day threshold",
       x = "Longitude", y = "Latitude") +
  theme_sdm() +
  guides(fill = guide_legend(nrow = 1))

save_fig(fig9, "Figure_9_model_agreement.png", width = 13, height = 8.5)
fig9

## [figure-10]
bar_df <- rbind(
  data.frame(label = change_areas$label, class = "Gain",
             area = change_areas$gain_km2, stringsAsFactors = FALSE),
  data.frame(label = change_areas$label, class = "Loss",
             area = -change_areas$loss_km2, stringsAsFactors = FALSE)
)
bar_df$label <- factor(bar_df$label, levels = unname(SCEN_LABEL[SCEN]))
net_df <- change_areas
net_df$label <- factor(net_df$label, levels = unname(SCEN_LABEL[SCEN]))

fig10 <- ggplot() +
  geom_col(data = bar_df, aes(label, area / 1e6, fill = class), width = 0.6) +
  geom_hline(yintercept = 0, colour = "grey25", linewidth = 0.4) +
  geom_point(data = net_df, aes(label, net_km2 / 1e6),
             shape = 23, size = 3.4, fill = "white", colour = "black",
             stroke = 0.9) +
  geom_text(data = net_df,
            aes(label, net_km2 / 1e6,
                label = sprintf("net %+.2f Mkm2 (%+.1f%%)",
                                net_km2 / 1e6, net_pct)),
            vjust = -1.4, size = 3.1) +
  scale_fill_manual(values = c(Gain = "#1f6f8b", Loss = "#c0392b"), name = NULL) +
  labs(title = "Change in climatically suitable area",
       subtitle = sprintf("Present-day suitable area %.2f million km2 | diamonds mark net change",
                          area_present / 1e6),
       x = NULL, y = "Area (million km2)") +
  theme_sdm() +
  theme(axis.text.x = element_text(angle = 12, hjust = 1))

save_fig(fig10, "Figure_10_area_change.png", width = 10, height = 6)
fig10

## [figure-11]
df11 <- add_scen_label(rast_df_multi(stack_scen(mess_rast)))
mess_cap <- stats::quantile(abs(df11$value), 0.99, na.rm = TRUE)
df11$value <- pmax(pmin(df11$value, mess_cap), -mess_cap)

fig11 <- ggplot() +
  geom_raster(data = df11, aes(x, y, fill = value)) +
  geom_polygon(data = world_df, aes(x, y, group = interaction(geom, part)),
               fill = NA, colour = "grey25", linewidth = 0.15) +
  scale_fill_gradient2(low = "#8c2d04", mid = "grey95", high = "#08519c",
                       midpoint = 0, name = "MESS",
                       limits = c(-mess_cap, mess_cap)) +
  coord_equal(xlim = XLIM, ylim = YLIM, expand = FALSE) +
  facet_wrap(~ label, ncol = 2) +
  labs(title = "Environmental novelty of the projected climate",
       subtitle = "Negative values mark cells where at least one predictor falls outside the training range; colour scale clipped at the 99th percentile of absolute MESS",
       x = "Longitude", y = "Latitude") +
  theme_sdm()

save_fig(fig11, "Figure_11_future_MESS.png", width = 13, height = 8.5)
fig11

## [figure-12]
prof_df <- profile_long
prof_df$series <- factor(prof_df$series,
                         levels = c("Present day", unname(SCEN_LABEL[SCEN])))

cent_df <- data.frame(
  series = factor(c("Present day", unname(SCEN_LABEL[SCEN])),
                  levels = levels(prof_df$series)),
  lat = c(cent_present[["lat"]], shift_tab[["Centroid lat"]]),
  stringsAsFactors = FALSE
)

series_cols <- stats::setNames(
  c("grey20", grDevices::hcl.colors(length(SCEN), "Dark 3")),
  levels(prof_df$series)
)

fig12 <- ggplot(prof_df, aes(lat_band, area_km2 / 1e3, colour = series)) +
  geom_line(linewidth = 0.7) +
  geom_vline(data = cent_df, aes(xintercept = lat, colour = series),
             linetype = "22", linewidth = 0.5, show.legend = FALSE) +
  scale_colour_manual(values = series_cols, name = NULL) +
  scale_x_continuous(breaks = seq(-40, 60, by = 10)) +
  coord_flip() +
  labs(title = "Latitudinal distribution of suitable area",
       subtitle = sprintf("Suitable area per %.1f degree latitude band | dashed lines mark the area-weighted centroid latitude",
                          FUT$lat_band),
       x = "Latitude", y = "Suitable area per band (thousand km2)") +
  theme_sdm()

save_fig(fig12, "Figure_12_latitudinal_profile.png", width = 9, height = 7)
fig12

## [figure-13]
ct_df <- country_by_scen |>
  filter(country %in% top_countries) |>
  mutate(country = factor(country, levels = rev(top_countries)),
         label   = factor(label, levels = unname(SCEN_LABEL[SCEN])))

fig13 <- ggplot(ct_df, aes(net_km2 / 1e3, country, fill = label)) +
  geom_col(position = position_dodge(width = 0.78), width = 0.72) +
  geom_vline(xintercept = 0, colour = "grey25", linewidth = 0.4) +
  scale_fill_viridis_d(option = "D", end = 0.88, name = NULL) +
  labs(title = "Net change in suitable area by country",
       subtitle = sprintf("The fifteen countries with the largest absolute change under %s",
                          SCEN_LABEL[[focus_key]]),
       x = "Net change in suitable area (thousand km2)", y = NULL) +
  theme_sdm() +
  guides(fill = guide_legend(nrow = 2))

save_fig(fig13, "Figure_13_country_change.png", width = 10, height = 8)
fig13


# How to read and report these results ====================================


# Reproducibility =========================================================


# What is deterministic here ----------------------------------------------

## [fut-session-info]
sessionInfo()

## [fut-write-manifest]
fut_files <- list.files(c(DIR_OUT, DIR_FUT, DIR_FIG), full.names = TRUE)
fut_files <- fut_files[!dir.exists(fut_files)]

fut_manifest <- data.frame(
  file     = basename(fut_files),
  dir      = basename(dirname(fut_files)),
  size_kb  = round(file.size(fut_files) / 1024, 1),
  modified = format(file.mtime(fut_files), "%Y-%m-%d %H:%M:%S"),
  stringsAsFactors = FALSE
)
utils::write.csv(fut_manifest, p_out("FUTURE_OUTPUT_MANIFEST.csv"),
                 row.names = FALSE)

total_mb <- round(sum(fut_manifest$size_kb) / 1024, 1)
message(sprintf("Outputs on disk: %d files, %s MB",
                nrow(fut_manifest), format(total_mb, big.mark = ",")))

knitr::kable(
  fut_manifest[grepl("^future$", fut_manifest$dir),
               c("file", "size_kb")],
  format.args = list(big.mark = ","),
  caption = "Projection rasters written to outputs/future/."
)


# Appendix: variations not used ===========================================

## [fut-appendix-mean-climate-then-predict]  (not run: kept for reference, as in the .qmd)
# # Predicting once from the ensemble-mean climate, rather than averaging the
# # per-model predictions. Cheaper by a factor of five, and it is what the
# # ensemble-mean climate raster written in Section @sec-fut-project would support
# # directly. Not used, because the model is non-linear: averaging climate first
# # blends incompatible projections into a state no model produced, and the result
# # systematically differs from the mean of the predictions. Kept here because the
# # comparison is a useful diagnostic of how much the non-linearity matters.
# for (key in SCEN) {
#   clim <- terra::rast(p_fut(paste0("climate_", key, ".tif")))
#   names(clim) <- BIOCLIM_VARS
#   fs <- c(clim, CROPLAND_CONST)
#   names(fs) <- PRED_VARS
#   alt <- terra::predict(terra::mask(fs, TEMPLATE), mx_final,
#                         type = "cloglog", clamp = FUT$clamp, na.rm = TRUE)
#   d <- terra::global(abs(alt - ens_mean[[key]]), c("mean", "max"), na.rm = TRUE)
#   message(sprintf("%s: mean |difference| %.4f, max %.4f",
#                   SCEN_LABEL[[key]], d[1, 1], d[1, 2]))
# }

## [fut-appendix-unclamped]  (not run: kept for reference, as in the .qmd)
# # Unclamped projection, for comparison with the clamped default. Run this to
# # quantify how much of the projected change depends on holding predictors at the
# # training limits. Expect the largest divergence in the hot and dry interiors,
# # where BIO12 falls below anything in the training set.
# FUT$clamp <- FALSE
# # ... then re-run Section @sec-fut-project with a different output prefix, and
# # difference the two ensemble means. Writing to the same paths would silently
# # overwrite the clamped results, so change SPECIES_SLUG or DIR_FUT first.

## [fut-appendix-scenario-thresholds]  (not run: kept for reference, as in the .qmd)
# # Re-deriving a 10th-percentile threshold from each future surface, instead of
# # carrying the present-day threshold across. Rejected: the training presences
# # are present-day observations, so a "training-presence percentile" of a future
# # surface has no clear meaning, and gain and loss areas computed this way are
# # partly an artefact of recalibration rather than a statement about climate.
# for (key in SCEN) {
#   v <- terra::extract(ens_mean[[key]],
#                       terra::vect(occ_pts, geom = c("lon", "lat"),
#                                   crs = "EPSG:4326"))[, 2]
#   message(sprintf("%s: future-surface 10th percentile at present occurrences %.4f (present-day threshold %.4f)",
#                   SCEN_LABEL[[key]],
#                   stats::quantile(v, CFG$threshold_quantile, na.rm = TRUE), THR))
# }

## [fut-appendix-luh2]  (not run: kept for reference, as in the .qmd)
# # Land-use-aware projection sketch. LUH2 supplies fractional crop area per
# # scenario on a 0.25 degree grid, which would need downscaling to this 2.5
# # arc-minute grid and reconciling with the ESA WorldCover definition the model
# # was fitted on. The reconciliation is the hard part: a bilinear downscale of a
# # coarse fraction is not the same quantity as a WorldCover-derived fraction, and
# # feeding one into a model calibrated on the other would introduce a bias of
# # unknown sign. Left undone deliberately rather than done badly.
# luh2_crop <- terra::rast("path/to/LUH2_cropland_ssp245_2050.tif")
# luh2_crop <- terra::resample(luh2_crop, TEMPLATE, method = "bilinear")
# names(luh2_crop) <- "cropland_fraction"


}  # end of PART 3
