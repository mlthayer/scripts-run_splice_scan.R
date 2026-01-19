#!/usr/bin/env Rscript

# WT vs MUT splice-motif scan (GT/AG) with simple heuristic scores
# Reproducible script: reads WT and MUT sequences, sanitizes, finds first diff,
# scans +/- window for donor/acceptor motifs, scores, builds comparison table,
# prints summary, and writes CSV output.
#
# Usage:
#   Rscript scripts/run_splice_scan.R --wt "ACGT..." --mut "ACGT..." --radius 200 --out outputs/results.csv
#
# Or run with defaults (edit wt_seq/mut_seq below):
#   Rscript scripts/run_splice_scan.R

suppressWarnings(suppressMessages({
  # No required packages; uses base R only
}))

# ---------- Helpers: argument parsing (base R) ----------
get_arg <- function(flag, default = NULL) {
  args <- commandArgs(trailingOnly = TRUE)
  idx <- match(flag, args)
  if (!is.na(idx) && idx < length(args)) return(args[[idx + 1]])
  default
}

flag_present <- function(flag) {
  flag %in% commandArgs(trailingOnly = TRUE)
}

if (flag_present("--help") || flag_present("-h")) {
  cat(
"Spn42Dd splice motif scan (base R)

Flags:
  --wt      WT sequence string (DNA)
  --mut     MUT sequence string (DNA)
  --radius  Window radius around first difference (default: 200)
  --out     Output CSV path (default: outputs/xxxxx_splice_sites_pm200_WT_vs_MUT_simpleScores.csv)

Example:
  Rscript scripts/run_splice_scan.R --wt \"ACGT...\" --mut \"ACGT...\" --radius 200 --out outputs/out.csv
")
  quit(status = 0)
}

# ---------- 1. Input sequences (defaults you can edit) ----------
wt_seq_default  <- "PASTE_WT_SEQUENCE_HERE"
mut_seq_default <- "PASTE_MUT_SEQUENCE_HERE"

wt_seq  <- get_arg("--wt", wt_seq_default)
mut_seq <- get_arg("--mut", mut_seq_default)
window_radius <- as.integer(get_arg("--radius", "200"))

out_default <- sprintf("outputs/xxxxx_splice_sites_pm%d_WT_vs_MUT_simpleScores.csv", window_radius)
out_file <- get_arg("--out", out_default)

if (wt_seq == "PASTE_WT_SEQUENCE_HERE" || mut_seq == "PASTE_MUT_SEQUENCE_HERE") {
  cat("NOTE: Using placeholder sequences. Provide --wt and --mut, or paste sequences into the script.\n")
}

# ---------- 2. Clean sequences to A/C/G/T only ----------
sanitize <- function(s) toupper(gsub("[^ACGT]", "", s))
wt  <- sanitize(wt_seq)
mut <- sanitize(mut_seq)

if (nchar(wt) < 2 || nchar(mut) < 2) {
  stop("WT or MUT sequence is too short after sanitization. Check input (non-ACGT characters are removed).")
}

# ---------- 3. Find the first nucleotide difference ----------
find_first_diff <- function(a, b) {
  L <- min(nchar(a), nchar(b))
  for (i in seq_len(L)) {
    if (substr(a, i, i) != substr(b, i, i)) return(i)
  }
  if (nchar(a) != nchar(b)) return(L + 1)
  NA_integer_
}

mut_nt <- find_first_diff(wt, mut)
cat("First NT difference at position:", mut_nt, "\n")

if (is.na(mut_nt)) {
  cat("No differences detected between WT and MUT after sanitization.\n")
}

# ---------- 4. Define +/- window around the mutation ----------
window_start_wt  <- max(1, mut_nt - window_radius)
window_end_wt    <- min(nchar(wt),  mut_nt + window_radius)
window_start_mut <- max(1, mut_nt - window_radius)
window_end_mut   <- min(nchar(mut), mut_nt + window_radius)

sub_wt  <- substr(wt,  window_start_wt,  window_end_wt)
sub_mut <- substr(mut, window_start_mut, window_end_mut)

cat("WT window:  ", window_start_wt,  "to", window_end_wt,  "\n")
cat("MUT window: ", window_start_mut, "to", window_end_mut, "\n")

# ---------- 5. Find donor (GT) and acceptor (AG) motifs and build windows ----------
# Donor (5'): 9-mer window: 3 nt upstream + "GT" + 4 nt downstream
# Acceptor (3'): 23-mer window: 20 nt upstream + "AG" + 1 nt downstream (heuristic anchor)

find_motifs <- function(seq, offset_start, allele, mut_nt_global) {
  n <- nchar(seq)
  res <- list()

  for (i in 1:(n - 1)) {
    din <- substr(seq, i, i + 1)

    # Donor motif: "GT"
    if (din == "GT") {
      ws <- i - 3
      we <- i + 5
      if (ws >= 1 && we <= n) {
        full_pos <- offset_start + i - 1
        win_seq  <- substr(seq, ws, we)
        res[[length(res) + 1]] <- data.frame(
          allele      = allele,
          site_type   = "donor",
          motif       = "GT",
          nt_pos      = full_pos,
          nt_pos_rel  = full_pos - mut_nt_global,
          window_seq  = win_seq,
          stringsAsFactors = FALSE
        )
      }
    }

    # Acceptor motif: "AG"
    if (din == "AG") {
      # NOTE: this uses a heuristic anchor; if you want exon start AFTER AG, use pos_exon <- i + 2
      pos_exon <- i + 1
      ws <- pos_exon - 20
      we <- pos_exon + 2
      if (ws >= 1 && we <= n) {
        full_pos <- offset_start + pos_exon - 1
        win_seq  <- substr(seq, ws, we)
        res[[length(res) + 1]] <- data.frame(
          allele      = allele,
          site_type   = "acceptor",
          motif       = "AG",
          nt_pos      = full_pos,
          nt_pos_rel  = full_pos - mut_nt_global,
          window_seq  = win_seq,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  if (length(res) == 0) {
    return(data.frame(
      allele=character(), site_type=character(), motif=character(),
      nt_pos=integer(), nt_pos_rel=integer(), window_seq=character(),
      stringsAsFactors=FALSE
    ))
  }
  do.call(rbind, res)
}

wt_sites  <- find_motifs(sub_wt,  window_start_wt,  "WT",  mut_nt_global = mut_nt)
mut_sites <- find_motifs(sub_mut, window_start_mut, "MUT", mut_nt_global = mut_nt)

sites_long <- rbind(wt_sites, mut_sites)

cat("\nNumber of candidate sites in +/-", window_radius, "nt window:\n")
print(table(sites_long$allele, sites_long$site_type))

# ---------- 6. Simple R-only splice-site strength scores ----------
# Donor score: # matches to consensus 9-mer "CAGGTAAGT" (0–9)
# Acceptor score: (# of C/T in positions 1–20) + bonus if AG near end (heuristic)

score_donor <- function(seq9) {
  if (nchar(seq9) != 9) return(NA_real_)
  consensus <- strsplit("CAGGTAAGT", "")[[1]]
  bases     <- strsplit(seq9, "")[[1]]
  sum(bases == consensus)
}

score_acceptor <- function(seq23) {
  if (nchar(seq23) != 23) return(NA_real_)
  bases <- strsplit(seq23, "")[[1]]
  poly_py <- sum(bases[1:20] %in% c("C", "T"))
  ag_end  <- ifelse(paste0(bases[20], bases[21]) == "AG", 5, 0)
  poly_py + ag_end
}

sites_long$score_simple <- NA_real_

sites_long$score_simple[sites_long$site_type == "donor"] <-
  vapply(sites_long$window_seq[sites_long$site_type == "donor"],
         score_donor, numeric(1))

sites_long$score_simple[sites_long$site_type == "acceptor"] <-
  vapply(sites_long$window_seq[sites_long$site_type == "acceptor"],
         score_acceptor, numeric(1))

# ---------- 7. Build WT vs MUT comparison table (wide) ----------
wt_tab  <- subset(sites_long, allele == "WT",
                  select = c(site_type, nt_pos_rel, nt_pos, window_seq, score_simple))
mut_tab <- subset(sites_long, allele == "MUT",
                  select = c(site_type, nt_pos_rel, nt_pos, window_seq, score_simple))

colnames(wt_tab)  <- c("site_type", "nt_pos_rel", "nt_pos_WT",  "window_seq_WT",  "score_WT")
colnames(mut_tab) <- c("site_type", "nt_pos_rel", "nt_pos_MUT", "window_seq_MUT", "score_MUT")

sites_wide <- merge(wt_tab, mut_tab, by = c("site_type", "nt_pos_rel"), all = TRUE)
sites_wide <- sites_wide[order(sites_wide$site_type, sites_wide$nt_pos_rel), ]

if ("score_WT" %in% names(sites_wide) && "score_MUT" %in% names(sites_wide)) {
  sites_wide$delta_score <- sites_wide$score_MUT - sites_wide$score_WT
}

# ---------- 8. Inspect and save ----------
cat("\nHead of WT vs MUT splice-site table (simple scores):\n")
print(utils::head(sites_wide, 20))

dir.create(dirname(out_file), showWarnings = FALSE, recursive = TRUE)
write.csv(sites_wide, out_file, row.names = FALSE)

cat("\nSaved WT vs MUT splice-site table to:\n", normalizePath(out_file), "\n")

# Optional: record session info for reproducibility
session_file <- sub("\\.csv$", "_sessionInfo.txt", out_file)
writeLines(c(capture.output(sessionInfo())), con = session_file)
cat("Saved sessionInfo() to:\n", normalizePath(session_file), "\n")
