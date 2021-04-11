# Functions to count number of reads and calculate readthrough efficiency
library(data.table)
library(dplyr)

# Functions ----------------------------------------------------------

## Count read in the CDS and extension region for each transcript
# dt: data table (a sample table in reads_psite_list format)
# annotation: annotation file containing information about lengths of UTRs, CDS, and position of the next in-frame stop codon from the canonical stop
# cds_m5: the number of nucleotides from the start codon (5' end of CDS) to EXCLUDE from the CDS count
# cds_m3: the number of nucleotides from the stop codon (3' end of CDS) to EXCLUDE from the CDS count
# frame_: reading frame to include in the count. Options are "all", 0, 1, 2, or a vector of pairs. Default is "all"

read_count <- function(dt, annotation, cds_m5 = 15, cds_m3 = 33, frame_ = "all") {
  dt <- data.table(dt)
  
  # Subset data
  message("Subsetting data")
  if (frame_ == "all") {
    message("\tFrame: 0, 1, 2")
    cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3, ]
    utr3_sub <- dt[psite_region == "3utr", ]
  } else if (frame_ < 0 | frame_ > 2) {
    message("\tInvalid frame number. Please put either 0, 1, 2, or all")
  } else {
    message(paste0("\tFrame: ", frame_))
    cds_sub <- dt[psite_region == "cds" & psite_from_start > cds_m5 & psite_from_stop < -cds_m3 & frame %in% frame_, ]
    utr3_sub <- dt[psite_region == "3utr" & frame %in% frame_, ]
  }
  utr3_sub <- merge(utr3_sub, annotation[, c("transcript", "nis_pos")], by = "transcript", all.x = TRUE) # Add nis information
  ext_sub <-  utr3_sub[psite_from_stop < nis_pos, ]   # extension
  distal_sub <- utr3_sub[psite_from_stop >= nis_pos, ]  # distal 3'-UTR

  # Count reads
  message("Counting reads")
  cds_tab <- cds_sub[, list(c_cds = .N), by = list(transcript)]
  utr3_tab <- utr3_sub[, list(c_utr3 = .N), by = list(transcript)]
  ext_tab <- ext_sub[, list(c_ext = .N), by = list(transcript)]
  distal_tab <- distal_sub[, list(c_dist = .N), by = list(transcript)]
  
  # Combine data
  message("Combining read count tables of different regions")
  cu_tab <- Reduce(function(df1, df2) merge(df1, df2, by = "transcript", all.x = TRUE), list(annotation, cds_tab, utr3_tab, ext_tab, distal_tab))
  message("Replacing NA count with 0")
  cu_tab[is.na(cu_tab)] <- 0
  
  # Calculate length of extension and distal 3'UTR
  cu_tab$l_ext <- pmin(cu_tab$nis_pos-1, cu_tab$l_utr3)
  cu_tab$l_dist <- cu_tab$l_utr3 - cu_tab$l_ext
  cu_tab$nis_pos <- NULL
  
  cu_tab <- data.table(cu_tab %>% mutate_if(is.character, as.factor))
  message("Done\n")
  return(cu_tab)
}

## Pool replicates read count
# data_list: list of data tables that are replicates to combine

composite <- function(data_list) {
  data_comp <- rbindlist(data_list, use.names = TRUE)[,lapply(.SD, sum), by = list(transcript, l_tr, l_utr5, l_cds, l_utr3, l_ext, l_dist)]
  return(data_comp)
}

## Calculate readthrough efficiency and RPKM for each region
# dt: read count data table
# cds_m5: the number of nucleotides from the start codon (5' end of CDS) to EXCLUDE from the CDS count
# cds_m3: the number of nucleotides from the stop codon (3' end of CDS) to EXCLUDE from the CDS count

rt_efficiency <- function(dt, cds_m5 = 15, cds_m3 = 33) {
  # Calculate RPKM
  lib_size <- (sum(dt$c_cds) + sum(dt$c_utr3))/10^6
  dt$rpkm_cds <- dt$c_cds/(lib_size * ((dt$l_cds - cds_m5 - cds_m3)/10^3))
  dt$rpkm_utr3 <- dt$c_utr3/(lib_size * ((dt$l_utr3)/10^3))
  dt$rpkm_ext <- dt$c_ext/(lib_size * ((dt$l_ext)/10^3))
  dt$rpkm_dist <- dt$c_dist/(lib_size * ((dt$l_dist)/10^3))
  
  # Calculate readthrough efficiency
  dt$rte_utr3 <- dt$rpkm_utr3/dt$rpkm_cds
  dt$rte_ext <- dt$rpkm_ext/dt$rpkm_cds
  return(dt)
}

# --------------------------------------------------------------------

# Run functions ------------------------------------------------------

# Read count
# This step takes a lot of memory and some time to load data. 
# One set of 4-5 samples may be done at a time and then removed from R environment to clear space.
nis <- read.table("../RData/next_inframe_stop.txt", header = TRUE)

sup45wt_data <- list()
sup45wt_data$WA25 <- read.table("../RData/sup45ts/sup45_wt_25C_A_reads_psite_list.txt.gz", header = TRUE)
sup45wt_data$WB25 <- read.table("../RData/sup45ts/sup45_wt_25C_B_reads_psite_list.txt.gz", header = TRUE)
sup45wt_data$WA37 <- read.table("../RData/sup45ts/sup45_wt_37C_A_reads_psite_list.txt.gz", header = TRUE)
sup45wt_data$WB37 <- read.table("../RData/sup45ts/sup45_wt_37C_B_reads_psite_list.txt.gz", header = TRUE)
f0_sup45wt_rc <- lapply(sup45wt_data, read_count, annotation = nis, cds_m5 = 15, cds_m3 = 33, frame_ = 0)
rm(sup45wt_data)  # remove to clear space

sup45ts_data <- list()
sup45ts_data$MA25 <- read.table("../RData/sup45ts/sup45_ts_25C_A_reads_psite_list.txt.gz", header = TRUE)
sup45ts_data$MB25 <- read.table("../RData/sup45ts/sup45_ts_25C_B_reads_psite_list.txt.gz", header = TRUE)
sup45ts_data$MA37 <- read.table("../RData/sup45ts/sup45_ts_37C_A_reads_psite_list.txt.gz", header = TRUE)
sup45ts_data$MB37 <- read.table("../RData/sup45ts/sup45_ts_37C_B_reads_psite_list.txt.gz", header = TRUE)
f0_sup45ts_rc <- lapply(sup45ts_data, read_count, annotation = nis, cds_m5 = 15, cds_m3 = 33, frame_ = 0)
rm(sup45ts_data)   # remove to clear space

sup45d_data <- list()
sup45d_data$W903 <- read.table("../RData/sup45d/sup45_wt_903_reads_psite_list.txt.gz", header = TRUE)
sup45d_data$W904 <- read.table("../RData/sup45d/sup45_wt_904_reads_psite_list.txt.gz", header = TRUE)
sup45d_data$d908 <- read.table("../RData/sup45d/sup45_d_908_reads_psite_list.txt.gz", header = TRUE)
sup45d_data$d909 <- read.table("../RData/sup45d/sup45_d_909_reads_psite_list.txt.gz", header = TRUE)
f0_sup45d_rc <- lapply(sup45d_data, read_count, annotation = nis, cds_m5 = 15, cds_m3 = 33, frame_ = 0)
rm(sup45d_data)   # remove to clear space

rli1_data <- list()
rli1_data$W309 <- read.table("../RData/rli1/rli1_wt_309_reads_psite_list.txt.gz", header = TRUE)
rli1_data$W310 <- read.table("../RData/rli1/rli1_wt_310_reads_psite_list.txt.gz", header = TRUE)
rli1_data$d311 <- read.table("../RData/rli1/rli1_d_311_reads_psite_list.txt.gz", header = TRUE)
rli1_data$d312 <- read.table("../RData/rli1/rli1_d_312_reads_psite_list.txt.gz", header = TRUE)
rli1_data$d319 <- read.table("../RData/rli1/rli1_d_319_reads_psite_list.txt.gz", header = TRUE)
f0_rli1_rc <- lapply(rli1_data, read_count, annotation = nis, cds_m5 = 15, cds_m3 = 33, frame_ = 0)
rm(rli1_data)   # remove to clear space

# Composite
f0_comp <- list() # Initialize empty list to store data
f0_comp$sup45_wt_25C <- composite(f0_sup45wt_rc[c("WA25", "WB25")])
f0_comp$sup45_wt_37C <- composite(f0_sup45wt_rc[c("WA37", "WB37")])
f0_comp$sup45_ts_25C <- composite(f0_sup45ts_rc[c("MA25", "MB25")])
f0_comp$sup45_ts_37C <- composite(f0_sup45ts_rc[c("MA37", "MB37")])
f0_comp$sup45_wt <- composite(f0_sup45d_rc[c("W903", "W904")])
f0_comp$sup45_d <- composite(f0_sup45d_rc[c("d908", "d909")])
f0_comp$rli1_wt <- composite(f0_rli1_rc[c("W309", "W310")])
f0_comp$rli1_d <- composite(f0_rli1_rc[c("d311", "d312", "d319")])

# RPKM and readthrough efficiency
f0_comp_re <- lapply(f0_comp, rt_efficiency, cds_m5 = 15, cds_m3 = 33)

save(f0_comp_re, file = "rte_f0_cds_m3=33.Rdata")
