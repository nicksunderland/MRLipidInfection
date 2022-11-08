
<!-- README.md is generated from README.Rmd. Please edit that file -->

# MRLipidInfection

<!-- badges: start -->
<!-- badges: end -->

## Installation

You can install the development version of MRLipidInfection like so:

``` r
# Download from GitHub, open package, load 'devtools', install()
```

## Example

``` r
library(MRLipidInfection)
library(devtools)
library(dplyr)
library(tidyr)
library(ieugwasr)
library(tabulizer)
library(TwoSampleMR)
load_all()

## Download the data from the FINNGEN repository
# data <- MRLipidInfection::getFinnGen("finngen_R7_AB1_SEPSIS.gz")

# List available GWASs
ao <- available_outcomes()

# Get instruments (LDL GWAS)
exposure_ldl_dat <- extract_instruments(c("ieu-a-300", "ieu-a-301", "ieu-a-781", "ieu-b-4845", "ieu-b-4846", "ieu-b-110"))

## OR

# Extract just the lipid medication instruments that Yarmolinski et al. used
Yarmolinski_pdf <- MRLipidInfection::getWebFile("https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7042851/bin/jama-323-646-s001.pdf")
snp_col     <- extract_tables(Yarmolinski_pdf, pages=8) %>% as.data.frame() %>% .[,1]
hmg_snps    <- snp_col[(grep("HMG-coA", snp_col, ignore.case=T)+1):(grep("NPC1L1", snp_col, ignore.case=T)-1)]
pcsk9_snps  <- snp_col[-(1:grep("PCSK9", snp_col, ignore.case=T))]
npc1l1_snps <- snp_col[(grep("NPC1L1", snp_col, ignore.case=T)+1):(grep("PCSK9", snp_col, ignore.case=T)-1)]
meds_snps   <- stack(list("HMG-CoA-R" = hmg_snps,
                          "NPC1L1"    = npc1l1_snps,
                          "PCSK9"     = pcsk9_snps)) %>% rename(c("SNP"="values", "med"="ind"))
exposure_meds_dat <- exposure_ldl_dat[exposure_ldl_dat$SNP %in% meds_snps$SNP, ]

# Outcome data
outcomes_studies <- ao[grepl("sepsis", ao$trait, ignore.case=T),]

# Get effects of LDL & med instruments on outcome
outcome_ldl_dat <- extract_outcome_data(snps=exposure_ldl_dat$SNP, outcomes=outcomes_studies$id)
outcome_med_dat <- extract_outcome_data(snps=exposure_meds_dat$SNP, outcomes=outcomes_studies$id)

# Harmonise the exposure and outcome data
dat_ldl <- harmonise_data(exposure_ldl_dat, outcome_ldl_dat)
dat_med <- harmonise_data(exposure_meds_dat, outcome_med_dat)

# Perform MR
res_ldl <- mr(dat_ldl)
res_med <- mr(dat_med)

# View result
res_ldl %>%
  as_tibble() %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(pval)

res_med %>%
  as_tibble() %>%
  filter(method == "Inverse variance weighted") %>%
  arrange(pval)
```
