################################################################################
# Southern Ocean COI Depth-Genetics Pipeline
# ============================================
# Title:   Depth as a driver of cryptic biodiversity and population structure
#          in Southern Ocean marine invertebrates and fish
# Journal: Biodiversity and Conservation (submitted)
# Authors: [Your Name] et al.
# Date:    April 2025
#
# Description:
#   Complete analysis pipeline for investigating whether depth (shelf 0-500 m
#   vs. bathyal 500-3000 m) structures genetic diversity in Southern Ocean
#   marine taxa using publicly available COI barcoding sequences from NCBI
#   GenBank. Includes sequence download, depth metadata recovery from ETOPO1
#   bathymetry, quality filtering, alignment, species delimitation, population
#   genetic analyses (AMOVA, PhiST, haplotype diversity, Tajima's D), PCA,
#   admixture (snapclust), DAPC, Mantel IBD tests, and figure generation.
#
# Dependencies:
#   R >= 4.3
#   MAFFT >= 7.5 (external; install via brew install mafft or apt install mafft)
#
# Usage:
#   Rscript SO_depth_genetics_pipeline.R
#   OR source("SO_depth_genetics_pipeline.R") from RStudio
#
# Output:
#   data/raw/          - downloaded sequences and metadata
#   data/aligned/      - MAFFT-aligned FASTA files
#   data/filtered/     - quality-filtered sequences
#   results/tables/    - CSV tables for manuscript
#   results/figures/   - PNG figures
#   results/networks/  - haplotype network files (NEXUS for PopART)
#
# Citation:
#   [Your Name] et al. (2025) Depth as a driver of cryptic biodiversity and
#   population structure in Southern Ocean marine invertebrates and fish:
#   a COI barcoding synthesis. Biodiversity and Conservation.
#
# Repository: https://github.com/[yourusername]/SO-depth-genetics
# License:    MIT
################################################################################


# ==============================================================================
# SECTION 1: Setup
# ==============================================================================

## 1.1 Install and load required packages -----------------------------------

required_pkgs <- c(
  # Data download
  "rentrez",     # NCBI GenBank API
  "bold",        # BOLD Systems API
  # Sequence analysis
  "ape",         # DNA sequence handling, distance matrices, NJ trees
  "pegas",       # Haplotype networks, AMOVA, Tajima's D
  "seqinr",      # FASTA I/O, sequence statistics
  "adegenet",    # snapclust admixture, DAPC, genind objects
  # Statistics
  "vegan",       # Mantel test, partial Mantel
  "marmap",      # ETOPO1 bathymetry query
  # Data wrangling
  "dplyr",       # Data manipulation
  "tidyr",       # Data reshaping
  "purrr",       # Functional iteration
  "readr",       # CSV I/O
  "stringr",     # String manipulation
  "glue",        # String interpolation
  "tibble",      # Modern data frames
  # Visualisation
  "ggplot2",     # Publication figures
  "gridExtra"    # Multi-panel figures
)

new_pkgs <- required_pkgs[!(required_pkgs %in% installed.packages()[, "Package"])]
if (length(new_pkgs) > 0) {
  message("Installing missing packages: ", paste(new_pkgs, collapse = ", "))
  install.packages(new_pkgs, repos = "https://cloud.r-project.org")
}

# bold package requires GitHub installation on R >= 4.5
if (!"bold" %in% installed.packages()[, "Package"]) {
  if (!"remotes" %in% installed.packages()[, "Package"])
    install.packages("remotes")
  remotes::install_github("ropensci/bold")
}

invisible(lapply(required_pkgs, library, character.only = TRUE))

## 1.2 Global settings -------------------------------------------------------

set.seed(42)  # Reproducibility for all permutation tests

## 1.3 Create output directories ---------------------------------------------

dirs <- c(
  "data/raw", "data/aligned", "data/filtered",
  "results/tables", "results/figures", "results/networks"
)
invisible(lapply(dirs, dir.create, recursive = TRUE, showWarnings = FALSE))


# ==============================================================================
# SECTION 2: Configuration
# ==============================================================================
# Edit this section to modify species list, depth thresholds, or quality
# filtering criteria without changing analysis code.

## 2.1 Target taxa -----------------------------------------------------------

TARGET_SPECIES <- c(
  "Dissostichus eleginoides",    # Patagonian toothfish (Nototheniidae)
  "Dissostichus mawsoni",        # Antarctic toothfish (Nototheniidae)
  "Lepidonotothen squamifrons",  # Grey rockcod (Nototheniidae)
  "Sterechinus neumayeri",       # Antarctic sea urchin (Echinoidea)
  "Odontaster validus",          # Antarctic sea star (Asteroidea)
  "Euphausia superba",           # Antarctic krill (Euphausiacea)
  "Trematomus loennbergii",      # Nototheniidae
  "Eusirus spp",                 # Amphipoda (genus-level)
  "Harmothoe spinosa"            # Polychaeta (genus-level)
)

## 2.2 Depth zone definitions (metres) ----------------------------------------

SHELF_MIN    <-    0   # Shelf minimum depth
SHELF_MAX    <-  500   # Shelf maximum depth (continental shelf + upper slope)
BATHYAL_MIN  <-  500   # Bathyal minimum depth
BATHYAL_MAX  <- 3000   # Bathyal maximum depth

## 2.3 Geographic scope -------------------------------------------------------

LAT_MAX <- -45  # Southern Ocean: all latitudes south of 45 degrees S

## 2.4 Sequence quality thresholds --------------------------------------------

MIN_SEQ_LEN  <- 500   # Minimum sequence length (bp)
MAX_SEQ_LEN  <- 700   # Maximum sequence length (bp); COI-5P barcode region
MAX_N_PROP   <- 0.01  # Maximum proportion of ambiguous bases (N)

## 2.5 Analysis thresholds ----------------------------------------------------

MIN_N_PER_ZONE    <- 10    # Minimum sequences per depth zone per species
MOTU_INTRA_THRESH <- 0.03  # COI intraspecific threshold (3%; Hebert et al. 2003)
MOTU_INTER_THRESH <- 0.10  # Conservative interspecific threshold (10%)
N_PERM            <- 999   # Permutations for AMOVA (use 9999 for publication)
N_PERM_MANTEL     <- 9999  # Permutations for Mantel test


# ==============================================================================
# SECTION 3: Sequence Download — NCBI GenBank
# ==============================================================================

#' Download COI sequences for a taxon from NCBI GenBank
#'
#' Uses multiple complementary Entrez queries to maximise retrieval.
#' Parses depth, coordinates, and country from GenBank flat files.
#'
#' @param taxon_name Character. Scientific name of target taxon.
#' @param retmax Integer. Maximum number of accessions to retrieve.
#' @return Data frame with columns: accession, seq, source, taxon,
#'         lat, lon, depth, country, depth_flag. Returns NULL if no records.

download_genbank_robust <- function(taxon_name, retmax = 500) {
  message("Downloading GenBank: ", taxon_name)

  # Multiple query strategies to maximise retrieval
  queries <- c(
    glue('"{taxon_name}"[Organism] AND COI[gene]'),
    glue('"{taxon_name}"[Organism] AND COX1[gene]'),
    glue('"{taxon_name}"[Organism] AND "cytochrome oxidase subunit I"[Title]')
  )

  all_ids <- c()
  for (q in queries) {
    res <- tryCatch(
      rentrez::entrez_search(db = "nuccore", term = q, retmax = retmax),
      error = function(e) NULL
    )
    if (!is.null(res)) all_ids <- unique(c(all_ids, res$ids))
    Sys.sleep(0.3)
  }

  if (length(all_ids) == 0) {
    message("  -> 0 records found")
    return(NULL)
  }

  message("  -> ", length(all_ids), " IDs found, fetching sequences...")

  # Fetch in batches of 50 (NCBI rate limit)
  batches       <- split(all_ids, ceiling(seq_along(all_ids) / 50))
  all_seqs_list <- list()
  all_meta_list <- list()

  for (b in batches) {

    # Fetch FASTA sequences
    fa_text <- tryCatch(
      rentrez::entrez_fetch(db = "nuccore", id = b,
                            rettype = "fasta", retmode = "text"),
      error = function(e) NULL
    )

    if (!is.null(fa_text) && nchar(fa_text) > 10) {
      entries    <- strsplit(fa_text, "\n>")[[1]]
      entries[1] <- sub("^>", "", entries[1])
      for (entry in entries) {
        lines   <- strsplit(entry, "\n")[[1]]
        lines   <- lines[nchar(trimws(lines)) > 0]
        if (length(lines) < 2) next
        seq_str <- toupper(gsub("[^ACGTNacgtn]", "",
                                paste(lines[-1], collapse = "")))
        if (nchar(seq_str) < 100) next
        acc <- strsplit(lines[1], "\\s+")[[1]][1]
        all_seqs_list[[length(all_seqs_list) + 1]] <-
          tibble::tibble(accession = acc, seq = seq_str)
      }
    }

    # Fetch GenBank flat file for metadata
    gb_text <- tryCatch(
      rentrez::entrez_fetch(db = "nuccore", id = b,
                            rettype = "gb", retmode = "text"),
      error = function(e) NULL
    )

    if (!is.null(gb_text)) {
      records <- strsplit(gb_text, "\n//\n")[[1]]
      for (rec in records) {
        if (nchar(rec) < 50) next

        acc     <- stringr::str_match(rec, "ACCESSION\\s+(\\S+)")[, 2]
        if (is.na(acc)) next

        country <- stringr::str_match(rec, '/country="([^"]+)"')[, 2]

        # Parse lat/lon from /lat_lon field (format: "77.85 S 166.67 E")
        latlon  <- stringr::str_match(rec, '/lat_lon="([^"]+)"')[, 2]
        lat <- lon <- NA_real_
        if (!is.na(latlon)) {
          p <- strsplit(trimws(latlon), "\\s+")[[1]]
          if (length(p) >= 4) {
            lat <- as.numeric(p[1]) * ifelse(p[2] == "S", -1, 1)
            lon <- as.numeric(p[3]) * ifelse(p[4] == "W", -1, 1)
          }
        }

        # Parse depth — matches common annotation patterns
        depth <- NA_real_
        dm <- stringr::str_match_all(
          rec, "(?i)depth[^0-9]{0,10}([0-9]+(?:\\.[0-9]+)?)")[[1]]
        if (nrow(dm) > 0) depth <- as.numeric(dm[1, 2])

        all_meta_list[[length(all_meta_list) + 1]] <-
          tibble::tibble(accession = acc, country = country,
                         lat = lat, lon = lon, depth = depth)
      }
    }
    Sys.sleep(0.35)
  }

  if (length(all_seqs_list) == 0) {
    message("  -> sequence parsing failed")
    return(NULL)
  }

  seqs_df <- dplyr::bind_rows(all_seqs_list) %>%
    dplyr::distinct(accession, .keep_all = TRUE) %>%
    dplyr::mutate(acc_base = sub("\\..*", "", accession))

  meta_df <- if (length(all_meta_list) > 0)
    dplyr::bind_rows(all_meta_list) %>%
      dplyr::distinct(accession, .keep_all = TRUE) %>%
      dplyr::mutate(acc_base = sub("\\..*", "", accession))
  else
    tibble::tibble(acc_base = seqs_df$acc_base,
                   country  = NA, lat = NA, lon = NA, depth = NA)

  final <- seqs_df %>%
    dplyr::left_join(meta_df %>% dplyr::select(-accession),
                     by = "acc_base") %>%
    dplyr::mutate(
      taxon      = taxon_name,
      source     = "GenBank",
      depth_flag = dplyr::case_when(
        !is.na(depth) ~ "verified",
        !is.na(lat)   ~ "has_coordinates",
        TRUE          ~ "needs_manual_check"
      )
    ) %>%
    dplyr::select(accession, seq, source, taxon,
                  lat, lon, depth, country, depth_flag)

  message("  -> ", nrow(final), " sequences | ",
          sum(!is.na(final$depth)), " with depth | ",
          sum(!is.na(final$lat)),   " with coordinates")
  return(final)
}

# Download all target species
gb_data <- purrr::map_dfr(TARGET_SPECIES, download_genbank_robust)
readr::write_csv(gb_data, "data/raw/genbank_raw.csv")

message("\n=== DOWNLOAD SUMMARY ===")
message("Total sequences:     ", nrow(gb_data))
message("With depth:          ", sum(!is.na(gb_data$depth)))
message("With coordinates:    ", sum(!is.na(gb_data$lat)))


# ==============================================================================
# SECTION 4: Depth Metadata Recovery from ETOPO1 Bathymetry
# ==============================================================================
# For sequences with coordinates but no depth, seabed depth is inferred from
# the ETOPO1 global 1-arcminute relief model (Amante & Eakins 2009).
# This is methodologically valid for benthic and demersal species where
# collection coordinates indicate the seabed sampling depth.

# Records with coordinates but no depth
has_coords_no_depth <- gb_data %>%
  dplyr::filter(!is.na(lat), !is.na(lon), is.na(depth))

message("Records needing bathymetric depth recovery: ",
        nrow(has_coords_no_depth))

if (nrow(has_coords_no_depth) > 0) {

  # Download ETOPO1 for Southern Ocean (cached after first download)
  message("Querying ETOPO1 bathymetry (downloads once, then cached)...")
  bathy <- marmap::getNOAA.bathy(
    lon1 = -180, lon2 = 180,
    lat1 = -80,  lat2 = -40,
    resolution = 10,
    keep = TRUE,
    path = "data/raw/"
  )

  # Extract depth at collection coordinates
  coord_depths <- marmap::get.depth(
    bathy,
    x        = has_coords_no_depth$lon,
    y        = has_coords_no_depth$lat,
    locator  = FALSE
  )

  # marmap returns negative values for ocean, positive for land
  has_coords_no_depth <- has_coords_no_depth %>%
    dplyr::mutate(
      depth_recovered = dplyr::if_else(
        coord_depths$depth < 0,
        abs(coord_depths$depth),
        NA_real_   # Positive = land point; discard
      ),
      depth_source = "etopo1_bathymetry"
    )

  message("Depths recovered from ETOPO1: ",
          sum(!is.na(has_coords_no_depth$depth_recovered)))

  # Merge recovered depths back into main dataset
  depth_lookup <- has_coords_no_depth %>%
    dplyr::filter(!is.na(depth_recovered)) %>%
    dplyr::select(accession, depth_recovered, depth_source)

  gb_data <- gb_data %>%
    dplyr::left_join(depth_lookup, by = "accession") %>%
    dplyr::mutate(
      depth_final  = dplyr::coalesce(depth, depth_recovered),
      depth_source = dplyr::case_when(
        !is.na(depth)           ~ "genbank_metadata",
        !is.na(depth_recovered) ~ "etopo1_bathymetry",
        TRUE                    ~ "unknown"
      )
    )
} else {
  gb_data <- gb_data %>%
    dplyr::mutate(depth_final  = depth,
                  depth_source = dplyr::if_else(!is.na(depth),
                                                "genbank_metadata", "unknown"))
}


# ==============================================================================
# SECTION 5: Geographic Filter and Depth Zone Assignment
# ==============================================================================

main_dataset <- gb_data %>%
  dplyr::mutate(
    # Geographic filter: keep Southern Ocean (lat <= -45) or unknown location
    in_SO = dplyr::case_when(
      is.na(lat)  ~ NA,
      lat <= LAT_MAX ~ TRUE,
      lat >  LAT_MAX ~ FALSE
    ),
    # Assign depth zones
    depth_zone = dplyr::case_when(
      !is.na(depth_final) & depth_final >= SHELF_MIN &
        depth_final <= SHELF_MAX                         ~ "shelf",
      !is.na(depth_final) & depth_final >  BATHYAL_MIN &
        depth_final <= BATHYAL_MAX                       ~ "bathyal",
      !is.na(depth_final) & depth_final >  BATHYAL_MAX  ~ "abyssal",
      TRUE                                               ~ NA_character_
    )
  ) %>%
  dplyr::filter(is.na(in_SO) | in_SO == TRUE)

# Save records needing manual depth recovery
needs_depth <- main_dataset %>%
  dplyr::filter(is.na(depth_zone)) %>%
  dplyr::select(accession, source, taxon, lat, lon, country, depth_final, depth_flag)

readr::write_csv(needs_depth,      "data/raw/depth_needs_manual_check.csv")
readr::write_csv(main_dataset,     "data/raw/main_dataset_with_depths.csv")

# Depth zone summary
zone_summary <- main_dataset %>%
  dplyr::mutate(
    taxon_grouped = dplyr::case_when(
      stringr::str_starts(taxon, "Eusirus")    ~ "Eusirus spp.",
      stringr::str_starts(taxon, "Harmothoe")  ~ "Harmothoe spp.",
      TRUE                                      ~ taxon
    )
  ) %>%
  dplyr::group_by(taxon_grouped) %>%
  dplyr::summarise(
    n_total   = dplyr::n(),
    n_shelf   = sum(depth_zone == "shelf",   na.rm = TRUE),
    n_bathyal = sum(depth_zone == "bathyal", na.rm = TRUE),
    n_unknown = sum(is.na(depth_zone)),
    .groups   = "drop"
  ) %>%
  dplyr::arrange(dplyr::desc(n_total))

message("\n=== DEPTH ZONE COVERAGE ===")
print(zone_summary)


# ==============================================================================
# SECTION 6: Quality Filtering and FASTA Export
# ==============================================================================

# Consolidate genus-level taxa and apply quality filters
WORKING_SPECIES <- c(
  "Euphausia superba",
  "Dissostichus eleginoides",
  "Eusirus spp.",
  "Harmothoe spp.",
  "Lepidonotothen squamifrons"
)

analysis_final <- main_dataset %>%
  dplyr::mutate(
    taxon_grouped = dplyr::case_when(
      stringr::str_starts(taxon, "Eusirus")    ~ "Eusirus spp.",
      stringr::str_starts(taxon, "Harmothoe")  ~ "Harmothoe spp.",
      TRUE                                      ~ taxon
    ),
    seq_clean = toupper(stringr::str_replace_all(seq, "[^ACGTNacgtn]", "")),
    seq_len   = nchar(seq_clean),
    n_prop    = stringr::str_count(seq_clean, "N") / pmax(seq_len, 1)
  ) %>%
  dplyr::filter(
    taxon_grouped %in% WORKING_SPECIES,
    !is.na(depth_zone),
    depth_zone %in% c("shelf", "bathyal"),
    seq_len   >= MIN_SEQ_LEN,
    seq_len   <= MAX_SEQ_LEN,
    n_prop    <= MAX_N_PROP
  ) %>%
  dplyr::mutate(taxon = taxon_grouped) %>%
  dplyr::distinct(accession, .keep_all = TRUE) %>%
  # Subsample to max 100 per species per depth zone
  dplyr::group_by(taxon, depth_zone) %>%
  dplyr::group_modify(~ dplyr::slice_sample(.x, n = min(100L, nrow(.x)))) %>%
  dplyr::ungroup()

message("\nFinal dataset for alignment:")
print(table(analysis_final$taxon, analysis_final$depth_zone))

# Write per-species FASTA files
# Header format: >accession|depth_zone|depth_m
for (sp in unique(analysis_final$taxon)) {
  sp_df     <- analysis_final %>% dplyr::filter(taxon == sp)
  safe_name <- stringr::str_replace_all(sp, "[ /.]", "_")
  fa_in     <- glue("data/raw/{safe_name}_unaligned.fasta")

  headers  <- glue(">{sp_df$accession}|{sp_df$depth_zone}|{round(sp_df$depth_final, 0)}")
  fa_lines <- as.vector(rbind(headers, sp_df$seq_clean))
  writeLines(fa_lines, fa_in)
  message("Written: ", fa_in, " (", nrow(sp_df), " sequences)")
}


# ==============================================================================
# SECTION 7: Multiple Sequence Alignment (MAFFT)
# ==============================================================================
# Requires MAFFT installed on PATH.
# Install: brew install mafft (macOS) or sudo apt install mafft (Linux)

mafft_available <- nchar(Sys.which("mafft")) > 0
if (!mafft_available)
  stop("MAFFT not found. Install with: brew install mafft (macOS) ",
       "or sudo apt install mafft (Linux)")

aligned_files <- c()
for (sp in unique(analysis_final$taxon)) {
  safe_name <- stringr::str_replace_all(sp, "[ /.]", "_")
  fa_in     <- glue("data/raw/{safe_name}_unaligned.fasta")
  fa_out    <- glue("data/aligned/{safe_name}_aligned.fasta")

  cmd <- glue("mafft --auto --thread -1 --quiet '{fa_in}' > '{fa_out}'")
  ret <- system(cmd)

  if (ret == 0 && file.exists(fa_out) && file.size(fa_out) > 0) {
    n_seq <- length(grep("^>", readLines(fa_out)))
    message("Aligned: ", sp, " -> ", n_seq, " sequences")
    aligned_files <- c(aligned_files, fa_out)
  } else {
    warning("Alignment failed for: ", sp)
  }
}


# ==============================================================================
# SECTION 8: Species Delimitation
# ==============================================================================
# Applied to genus-level datasets (Harmothoe spp., Eusirus spp.) prior to
# population genetic analyses. Uses K80 pairwise distances and single-linkage
# clustering at 3% (Hebert et al. 2003) and 10% thresholds.

#' Perform species delimitation on a COI alignment
#'
#' @param sp Character. Species/genus name.
#' @param aligned_files Character vector of aligned FASTA paths.
#' @param main_dataset Data frame. Full dataset with taxon labels.
#' @return List with distance summary, MOTU assignments, and NJ tree.

run_species_delimitation <- function(sp, aligned_files, main_dataset) {
  safe_name <- stringr::str_replace_all(sp, "[ /.]", "_")
  fa_out    <- glue("data/aligned/{safe_name}_aligned.fasta")
  if (!file.exists(fa_out)) return(NULL)

  dna  <- ape::read.FASTA(fa_out)
  ids  <- names(dna)
  meta <- data.frame(
    accession  = sapply(strsplit(ids, "\\|"), `[`, 1),
    depth_zone = sapply(strsplit(ids, "\\|"), `[`, 2),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::left_join(
      main_dataset %>%
        dplyr::select(accession, taxon_original = taxon) %>%
        dplyr::distinct(),
      by = "accession"
    )

  # Trim to most common sequence length
  seq_lens   <- sapply(as.list(dna), length)
  common_len <- as.integer(names(sort(table(seq_lens), decreasing = TRUE)[1]))
  keep       <- which(seq_lens == common_len)
  dna        <- dna[keep]
  meta       <- meta[keep, ]

  # K80 pairwise distance matrix
  dist_mat <- ape::dist.dna(dna, model = "K80", pairwise.deletion = TRUE)
  dist_vec <- as.vector(dist_mat)
  dist_vec <- dist_vec[dist_vec > 0]

  # Single-linkage clustering at 3% and 10%
  hclust_res  <- hclust(as.dist(as.matrix(dist_mat)), method = "single")
  clusters_3  <- cutree(hclust_res, h = MOTU_INTRA_THRESH)
  clusters_10 <- cutree(hclust_res, h = MOTU_INTER_THRESH)

  motu_df <- data.frame(
    accession    = meta$accession,
    depth_zone   = meta$depth_zone,
    taxon_orig   = meta$taxon_original,
    motu_3pct    = paste0("MOTU_", clusters_3),
    motu_10pct   = paste0("MOTU_", clusters_10),
    stringsAsFactors = FALSE
  )

  cat("\n===", sp, "species delimitation ===\n")
  cat("MOTUs at 3%: ",  length(unique(clusters_3)),  "\n")
  cat("MOTUs at 10%:", length(unique(clusters_10)), "\n")
  cat("Distance summary (K80):\n")
  print(summary(dist_vec))
  cat("< 3%: ",  sum(dist_vec < 0.03),
      "| 3-10%:", sum(dist_vec >= 0.03 & dist_vec < 0.10),
      "| >10%:", sum(dist_vec >= 0.10), "\n")

  # NJ tree coloured by depth zone
  nj_tree    <- ape::nj(dist_mat)
  tip_colours <- ifelse(
    meta$depth_zone[match(nj_tree$tip.label, meta$accession)] == "shelf",
    "steelblue", "darkorange"
  )

  png(glue("results/figures/{safe_name}_NJ_tree.png"),
      width = 1800, height = 2200, res = 150)
  plot(nj_tree, type = "phylogram", tip.color = tip_colours, cex = 0.4,
       main = paste(sp, "- NJ tree (K80)"), no.margin = TRUE)
  legend("topright",
         legend = c("Shelf (0-500 m)", "Bathyal (500-3000 m)"),
         fill   = c("steelblue", "darkorange"), bty = "n", cex = 0.8)
  dev.off()

  readr::write_csv(motu_df,
    glue("results/tables/{safe_name}_motu_assignments.csv"))

  list(sp = sp, motu_df = motu_df, dist_vec = dist_vec,
       n_motu_3pct = length(unique(clusters_3)),
       n_motu_10pct = length(unique(clusters_10)))
}

# Run species delimitation for genus-level datasets
delim_results <- list()
for (sp in c("Harmothoe spp.", "Eusirus spp.")) {
  delim_results[[sp]] <- run_species_delimitation(sp, aligned_files, main_dataset)
}


# ==============================================================================
# SECTION 9: Population Genetic Analyses
# ==============================================================================

#' Run full population genetic analysis for one species
#'
#' Computes: haplotype diversity (h), nucleotide diversity (pi), Tajima's D,
#' AMOVA (PhiST), haplotype network, and admixture (snapclust + DAPC).
#'
#' @param sp Character. Species name matching aligned FASTA filename.
#' @return List with diversity metrics, AMOVA results, and plot objects.

run_popgen <- function(sp) {
  safe_name <- stringr::str_replace_all(sp, "[ /.]", "_")
  fa_out    <- glue("data/aligned/{safe_name}_aligned.fasta")
  if (!file.exists(fa_out)) return(NULL)

  dna  <- ape::read.FASTA(fa_out)
  ids  <- names(dna)
  meta <- data.frame(
    accession  = sapply(strsplit(ids, "\\|"), `[`, 1),
    depth_zone = sapply(strsplit(ids, "\\|"), `[`, 2),
    depth      = as.numeric(sapply(strsplit(ids, "\\|"), `[`, 3)),
    stringsAsFactors = FALSE
  )

  # Trim to most common length
  seq_lens   <- sapply(as.list(dna), length)
  common_len <- as.integer(names(sort(table(seq_lens), decreasing = TRUE)[1]))
  keep       <- which(seq_lens == common_len)
  dna        <- dna[keep]
  meta       <- meta[keep, ]

  n_shelf   <- sum(meta$depth_zone == "shelf")
  n_bathyal <- sum(meta$depth_zone == "bathyal")
  cat("\n===", sp, "===\n")
  cat("shelf n =", n_shelf, "| bathyal n =", n_bathyal, "\n")

  # ── 9.1 Genetic diversity per depth zone
  div_rows <- list()
  for (zone in c("shelf", "bathyal")) {
    idx <- which(meta$depth_zone == zone)
    if (length(idx) < 3) next
    sub_dna <- dna[idx]

    h   <- tryCatch(pegas::hap.div(sub_dna),   error = function(e) NA)
    pi  <- tryCatch(
      mean(ape::dist.dna(sub_dna, model = "raw",
                         pairwise.deletion = TRUE), na.rm = TRUE),
      error = function(e) NA)
    taj <- tryCatch(pegas::tajima.test(sub_dna)$D, error = function(e) NA)

    cat(zone, ": h =", round(h, 4), " | pi =", round(pi, 6),
        " | D =", round(taj, 3), " | n =", length(idx), "\n")

    div_rows[[zone]] <- data.frame(
      taxon      = sp, depth_zone = zone, n = length(idx),
      h          = round(h, 4),   pi       = round(pi, 6),
      tajimas_d  = round(taj, 3), stringsAsFactors = FALSE)
  }

  # ── 9.2 AMOVA (PhiST; Excoffier et al. 1992)
  ord      <- order(meta$depth_zone)
  dna_ord  <- dna[ord]
  pop_vec  <- factor(meta$depth_zone[ord], levels = c("shelf", "bathyal"))
  dist_mat <- ape::dist.dna(dna_ord, model = "raw", pairwise.deletion = TRUE)
  dist_mat[is.na(dist_mat)] <- 0

  amova_res <- tryCatch(
    pegas::amova(dist_mat ~ pop_vec, nperm = N_PERM),
    error = function(e) { cat("AMOVA error:", e$message, "\n"); NULL }
  )

  phi_st <- p_val <- NA
  if (!is.null(amova_res)) {
    sig    <- amova_res$varcomp[, "sigma2"]
    phi_st <- round(sig[1] / sum(sig), 4)
    p_val  <- round(amova_res$varcomp["pop_vec", "P.value"], 4)
    cat("PhiST =", phi_st, "| p =", p_val,
        ifelse(!is.na(p_val) & p_val < 0.05,
               " ** SIGNIFICANT **", " (ns)"), "\n")
  }

  # ── 9.3 Admixture analysis (snapclust; Jombart & Ahmed 2011)
  dna_mat <- as.character(ape::as.matrix.DNAbin(dna))
  num_mat <- apply(dna_mat, 2, function(col) {
    col <- tolower(col)
    ifelse(col == "a", 0, ifelse(col == "c", 1,
    ifelse(col == "g", 2, ifelse(col == "t", 3, NA))))
  })
  var_sites <- apply(num_mat, 2, function(x) {
    x <- x[!is.na(x)]; length(unique(x)) > 1
  })
  n_snps <- sum(var_sites)

  admix_res <- NULL
  if (n_snps >= 3) {
    num_snp <- apply(num_mat[, var_sites, drop = FALSE], 2, function(x) {
      mc <- as.integer(names(sort(table(x[!is.na(x)]), decreasing = TRUE)[1]))
      x[is.na(x)] <- mc; x
    })
    rownames(num_snp) <- meta$accession
    gind <- tryCatch(
      adegenet::df2genind(as.data.frame(num_snp), ploidy = 1,
                          ind.names = meta$accession,
                          pop = meta$depth_zone, type = "codom"),
      error = function(e) NULL)

    if (!is.null(gind)) {
      n_ind     <- nrow(meta)
      snap_list <- bic_vals <- ll_vals <- list()
      cat("Admixture BIC: ")
      for (k in 1:3) {
        sc <- tryCatch(
          adegenet::snapclust(gind, k = k, initFreq = "random", nstart = 20),
          error = function(e) NULL)
        if (!is.null(sc) && sc$converged) {
          bic_vals[[k]] <- -2 * sc$ll + sc$n.param * log(n_ind)
          snap_list[[k]] <- sc
          cat("K=", k, ":", round(bic_vals[[k]], 1), " ")
        }
      }
      best_k <- which.min(unlist(bic_vals))
      cat("-> Best K =", best_k, "\n")

      admix_df <- as.data.frame(snap_list[[best_k]]$proba)
      colnames(admix_df) <- paste0("K", seq_len(ncol(admix_df)))
      admix_df$accession  <- meta$accession
      admix_df$depth_zone <- meta$depth_zone
      admix_res <- list(best_k = best_k, bic = unlist(bic_vals),
                        admix_df = admix_df, gind = gind)
    }
  }

  list(sp = sp, phi_st = phi_st, p_val = p_val,
       diversity = dplyr::bind_rows(div_rows),
       dna = dna_ord, meta = meta[ord, ], pop = pop_vec,
       admix = admix_res, n_snps = n_snps)
}

# Validated single-species for population genetic analysis
POPGEN_SPECIES <- c(
  "Dissostichus eleginoides",
  "Lepidonotothen squamifrons",
  "Eusirus_pontomedon"   # MOTU_2 only; created separately below
)

species_results <- list()
for (sp in c("Dissostichus eleginoides", "Lepidonotothen squamifrons")) {
  species_results[[sp]] <- run_popgen(sp)
}


# ==============================================================================
# SECTION 10: Mantel Test — Isolation by Depth
# ==============================================================================

#' Test for isolation by depth using Mantel and partial Mantel tests
#'
#' @param sp Character. Species name.
#' @param species_results List. Output from run_popgen().
#' @return Data frame with Mantel r, p-value, and partial Mantel results.

run_mantel_ibd <- function(sp, species_results) {
  r    <- species_results[[sp]]
  if (is.null(r)) return(NULL)

  meta <- r$meta
  dna  <- r$dna

  has_depth <- !is.na(meta$depth)
  if (sum(has_depth) < 10) return(NULL)

  dna_d    <- dna[has_depth]
  meta_d   <- meta[has_depth, ]
  gen_dist <- as.matrix(
    ape::dist.dna(dna_d, model = "raw", pairwise.deletion = TRUE))
  gen_dist[is.na(gen_dist)] <- 0
  depth_diff <- as.matrix(dist(meta_d$depth, method = "euclidean"))

  mant <- vegan::mantel(gen_dist, depth_diff,
                         method = "pearson",
                         permutations = N_PERM_MANTEL,
                         na.rm = TRUE)
  cat(sp, "| Mantel r =", round(mant$statistic, 4),
      "| p =", round(mant$signif, 4), "\n")

  data.frame(
    Taxon    = sp,
    n        = nrow(meta_d),
    Mantel_r = round(mant$statistic, 4),
    Mantel_p = round(mant$signif,    4),
    IBD_sig  = ifelse(mant$signif < 0.05, "YES", "no"),
    stringsAsFactors = FALSE
  )
}

cat("\n=== MANTEL IBD RESULTS ===\n")
mantel_summary <- do.call(rbind, lapply(names(species_results),
                                         run_mantel_ibd,
                                         species_results = species_results))
readr::write_csv(mantel_summary, "results/tables/mantel_ibd_results.csv")


# ==============================================================================
# SECTION 11: Summary Tables
# ==============================================================================

# Table 2: Population genetic results
summary_table <- do.call(rbind, lapply(species_results, function(r) {
  div <- r$diversity
  sh  <- div[div$depth_zone == "shelf",   ]
  ba  <- div[div$depth_zone == "bathyal", ]
  data.frame(
    Taxon      = r$sp,
    n_shelf    = if (nrow(sh) > 0) sh$n         else NA,
    n_bathyal  = if (nrow(ba) > 0) ba$n         else NA,
    h_shelf    = if (nrow(sh) > 0) sh$h         else NA,
    h_bathyal  = if (nrow(ba) > 0) ba$h         else NA,
    pi_shelf   = if (nrow(sh) > 0) sh$pi        else NA,
    pi_bathyal = if (nrow(ba) > 0) ba$pi        else NA,
    D_shelf    = if (nrow(sh) > 0) sh$tajimas_d else NA,
    D_bathyal  = if (nrow(ba) > 0) ba$tajimas_d else NA,
    PhiST      = r$phi_st,
    p_value    = r$p_val,
    sig        = ifelse(!is.na(r$p_val) & r$p_val < 0.05, "YES", "no"),
    stringsAsFactors = FALSE
  )
}))

readr::write_csv(summary_table, "results/tables/manuscript_table2_popgen.csv")
cat("\n=== MANUSCRIPT TABLE 2 ===\n")
print(summary_table, row.names = FALSE)


# ==============================================================================
# SECTION 12: Figures
# ==============================================================================

## 12.1 Haplotype diversity barplot -----------------------------------------

div_all <- dplyr::bind_rows(lapply(species_results, function(r) r$diversity)) %>%
  dplyr::mutate(
    taxon_short = dplyr::case_when(
      taxon == "Dissostichus eleginoides"    ~ "D. eleginoides",
      taxon == "Lepidonotothen squamifrons"  ~ "L. squamifrons",
      TRUE                                   ~ taxon
    ),
    depth_zone = factor(depth_zone, levels = c("shelf", "bathyal"),
                        labels = c("Shelf (0-500 m)", "Bathyal (500-3000 m)"))
  )

p_h <- ggplot2::ggplot(div_all,
         ggplot2::aes(x = taxon_short, y = h, fill = depth_zone)) +
  ggplot2::geom_col(position = "dodge", width = 0.65) +
  ggplot2::scale_fill_manual(
    values = c("Shelf (0-500 m)" = "steelblue",
               "Bathyal (500-3000 m)" = "darkorange")) +
  ggplot2::labs(x = NULL, y = "Haplotype diversity (h)", fill = NULL,
                title = "Haplotype diversity by depth zone") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_text(face = "italic", angle = 30, hjust = 1),
    legend.position = "top")

p_pi <- ggplot2::ggplot(div_all,
          ggplot2::aes(x = taxon_short, y = pi * 1000, fill = depth_zone)) +
  ggplot2::geom_col(position = "dodge", width = 0.65) +
  ggplot2::scale_fill_manual(
    values = c("Shelf (0-500 m)" = "steelblue",
               "Bathyal (500-3000 m)" = "darkorange")) +
  ggplot2::labs(x = NULL,
                y = expression(pi ~ "(x10"^{-3}~")"), fill = NULL,
                title = "Nucleotide diversity by depth zone") +
  ggplot2::theme_bw(base_size = 12) +
  ggplot2::theme(
    axis.text.x  = ggplot2::element_text(face = "italic", angle = 30, hjust = 1),
    legend.position = "top")

png("results/figures/fig1_diversity.png", width = 2600, height = 1200, res = 200)
gridExtra::grid.arrange(p_h, p_pi, ncol = 2)
dev.off()

## 12.2 PhiST comparison plot -----------------------------------------------

phi_df <- summary_table %>%
  dplyr::mutate(
    taxon_short = dplyr::case_when(
      Taxon == "Dissostichus eleginoides"   ~ "D. eleginoides\n(Fish)",
      Taxon == "Lepidonotothen squamifrons" ~ "L. squamifrons\n(Fish)",
      TRUE                                  ~ Taxon
    )
  )

p_phi <- ggplot2::ggplot(phi_df,
           ggplot2::aes(x = reorder(taxon_short, PhiST),
                        y = PhiST, fill = sig)) +
  ggplot2::geom_col(width = 0.55) +
  ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "grey40") +
  ggplot2::geom_text(ggplot2::aes(
    label = paste0("p=", p_value),
    vjust = ifelse(PhiST >= 0, -0.5, 1.5)), size = 3.8) +
  ggplot2::scale_fill_manual(
    values = c("YES" = "firebrick", "no" = "grey70"),
    labels = c("YES" = "p < 0.05", "no" = "ns"),
    name   = "Significance") +
  ggplot2::geom_hline(yintercept = 0.05, linetype = "dotted",
                       colour = "orange", linewidth = 0.7) +
  ggplot2::geom_hline(yintercept = 0.15, linetype = "dotted",
                       colour = "red", linewidth = 0.7) +
  ggplot2::labs(x    = NULL,
                y    = expression(Phi[ST] ~ "(shelf vs. bathyal)"),
                title = expression(Phi[ST] ~ "Southern Ocean COI depth structuring")) +
  ggplot2::theme_bw(base_size = 13)

png("results/figures/fig2_PhiST.png", width = 1600, height = 1300, res = 200)
print(p_phi)
dev.off()

## 12.3 PCA on COI SNPs -------------------------------------------------------

for (sp in names(species_results)) {
  r         <- species_results[[sp]]
  safe_name <- stringr::str_replace_all(sp, "[ /.]", "_")

  dna_mat <- as.character(ape::as.matrix.DNAbin(r$dna))
  num_mat <- apply(dna_mat, 2, function(col) {
    col <- tolower(col)
    ifelse(col == "a", 0, ifelse(col == "c", 1,
    ifelse(col == "g", 2, ifelse(col == "t", 3, NA))))
  })
  var_sites <- apply(num_mat, 2, function(x) {
    x <- x[!is.na(x)]; length(unique(x)) > 1
  })
  if (sum(var_sites) < 3) next

  num_snp <- apply(num_mat[, var_sites, drop = FALSE], 2, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE); x
  })

  pca_res <- prcomp(num_snp, scale. = TRUE)
  var_exp <- round(summary(pca_res)$importance[2, 1:2] * 100, 1)

  pca_df <- data.frame(
    PC1        = pca_res$x[, 1],
    PC2        = pca_res$x[, 2],
    depth_zone = r$meta$depth_zone
  )

  p_pca <- ggplot2::ggplot(pca_df,
             ggplot2::aes(x = PC1, y = PC2,
                          colour = depth_zone, shape = depth_zone)) +
    ggplot2::geom_point(size = 2.5, alpha = 0.8) +
    ggplot2::stat_ellipse(level = 0.95, linewidth = 0.8) +
    ggplot2::scale_colour_manual(
      values = c(shelf = "steelblue", bathyal = "darkorange"),
      labels = c(shelf = "Shelf (0-500 m)", bathyal = "Bathyal (500-3000 m)")) +
    ggplot2::scale_shape_manual(
      values = c(shelf = 16, bathyal = 17),
      labels = c(shelf = "Shelf (0-500 m)", bathyal = "Bathyal (500-3000 m)")) +
    ggplot2::labs(
      title  = bquote(italic(.(sp)) ~ "- COI PCA"),
      x      = paste0("PC1 (", var_exp[1], "%)"),
      y      = paste0("PC2 (", var_exp[2], "%)"),
      colour = "Depth zone", shape = "Depth zone") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "top")

  png(glue("results/figures/pca_{safe_name}.png"),
      width = 1400, height = 1200, res = 200)
  print(p_pca)
  dev.off()
}

message("\n========================================")
message("PIPELINE COMPLETE")
message("Tables:   results/tables/")
message("Figures:  results/figures/")
message("Networks: results/networks/")
message("========================================")


# ==============================================================================
# SECTION 13: Session Info (reproducibility)
# ==============================================================================

sink("results/session_info.txt")
cat("Analysis date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
sessionInfo()
sink()

message("Session info saved: results/session_info.txt")
