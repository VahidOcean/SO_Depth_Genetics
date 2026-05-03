##############################################################################
##  ARLEQUIN INPUT GENERATOR
##  Generates .arp files for Arlequin 3.5 from filtered FASTA + metadata
##  Run AFTER SO_depth_genetics_pipeline.R has produced data/filtered/
##
##  Arlequin 3.5: http://cmpg.unibe.ch/software/arlequin35/
##  On Linux/Mac: run via Wine or ArpView
##  The .arp files produced here can also be submitted to the Arlequin web
##  server if desktop installation is not possible.
##############################################################################

library(ape)
library(dplyr)
library(glue)
library(readr)
library(stringr)

# ── Config (must match main pipeline)
MIN_N_PER_ZONE <- 10
TARGET_SPECIES <- c(
  "Dissostichus eleginoides",
  "Dissostichus mawsoni",
  "Lepidonotothen squamifrons",
  "Sterechinus neumayeri",
  "Odontaster validus",
  "Euphausia superba"
)

# ── Load filtered metadata (written by main pipeline)
# If you added manually recovered depths, re-run main pipeline first.

generate_arlequin <- function(sp_name) {
  safe   <- str_replace_all(sp_name, " ", "_")
  fa     <- glue("data/filtered/{safe}_filtered.fasta")

  if (!file.exists(fa)) {
    message("Filtered FASTA not found: ", fa, " — run main pipeline first")
    return(invisible(NULL))
  }

  dna  <- ape::read.FASTA(fa)
  seqs <- as.character(ape::as.character.DNAbin(dna))
  ids  <- names(dna)

  # Parse metadata from FASTA header: accession|depth_zone|depth|lat|lon
  meta <- tibble::tibble(id = ids) %>%
    tidyr::separate(id,
                    into  = c("accession","depth_zone","depth","lat","lon"),
                    sep   = "\\|", extra = "merge", fill = "right") %>%
    mutate(depth_zone = if_else(is.na(depth_zone), "unknown", depth_zone))

  # Split into populations
  shelf_idx   <- which(meta$depth_zone == "shelf")
  bathyal_idx <- which(meta$depth_zone == "bathyal")

  if (length(shelf_idx) < MIN_N_PER_ZONE || length(bathyal_idx) < MIN_N_PER_ZONE) {
    message(sp_name, ": skipping Arlequin — insufficient n")
    return(invisible(NULL))
  }

  n_shelf   <- length(shelf_idx)
  n_bathyal <- length(bathyal_idx)
  seq_len   <- ncol(seqs)

  # Convert sequences to strings
  seq_strings <- apply(seqs, 1, paste, collapse = "")

  # ── Build .arp content
  arp_lines <- c(
    glue("[Profile]"),
    glue('  Title="{sp_name} COI shelf vs bathyal"'),
    glue("  NbSamples=2"),
    glue("  DataType=DNA"),
    glue("  GenotypicData=0"),
    glue("  LocusSeparator=NONE"),
    glue("  MissingData='?'"),
    glue("  CompDistMatrix=0"),
    "",
    "[Data]",
    "[[Samples]]",
    "",
    # Sample 1: Shelf
    glue('SampleName="Shelf_0_500m"'),
    glue("SampleSize={n_shelf}"),
    "SampleData={",
    {
      idx <- shelf_idx
      lines <- character(length(idx))
      for (i in seq_along(idx)) {
        acc <- meta$accession[idx[i]]
        sq  <- seq_strings[idx[i]]
        lines[i] <- glue("  {acc}  1  {sq}")
      }
      lines
    },
    "}",
    "",
    # Sample 2: Bathyal
    glue('SampleName="Bathyal_500_3000m"'),
    glue("SampleSize={n_bathyal}"),
    "SampleData={",
    {
      idx <- bathyal_idx
      lines <- character(length(idx))
      for (i in seq_along(idx)) {
        acc <- meta$accession[idx[i]]
        sq  <- seq_strings[idx[i]]
        lines[i] <- glue("  {acc}  1  {sq}")
      }
      lines
    },
    "}",
    "",
    "[[Structure]]",
    glue('StructureName="Shelf_vs_Bathyal"'),
    "NbGroups=1",
    "Group={",
    '  "Shelf_0_500m"',
    '  "Bathyal_500_3000m"',
    "}"
  )

  out_path <- glue("results/arlequin/{safe}.arp")
  dir.create("results/arlequin", showWarnings = FALSE, recursive = TRUE)
  writeLines(arp_lines, out_path)
  message("Written: ", out_path)

  # ── Also write matching .arl settings file
  arl_lines <- c(
    "[Setting_From_Arlequin]",
    "TaskDNA=110100",           # AMOVA + pairwise distances
    "FreqConsComp=0",
    "NbPermutations=10000",
    "HardyWeinbergTestType=0",
    "SignificanceLevelHW=0.050000",
    "MultipleTestCorrection=1",
    "ComputeStandardErrors=1",
    "DistanceMethod=6",         # TrN model
    "GammaValue=0.000000",
    "NbHaplotypesPerm=1000",
    "NbDifferentiationTests=10000",
    "PrecisionMissingData=0.01",
    "DataRepeatedMeasures=0"
  )

  arl_path <- glue("results/arlequin/{safe}.ars")
  writeLines(arl_lines, arl_path)
  message("Settings: ", arl_path)

  invisible(list(arp = out_path, ars = arl_path))
}

# Generate for all target species
purrr::walk(TARGET_SPECIES, generate_arlequin)

message("\nArlequin files in: results/arlequin/")
message("Open each .arp in Arlequin 3.5 GUI with matching .ars settings,")
message("or run headlessly with: arlequin35s <file.arp>")
message("\nKey outputs from Arlequin to extract for manuscript:")
message("  - Table 2: AMOVA results (% variance among groups)")
message("  - Table 3: Pairwise FST and PhiST with p-values")
message("  - Supplement: Full distance matrices")
