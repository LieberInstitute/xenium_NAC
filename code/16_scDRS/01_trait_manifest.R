# Build the curated GWAS trait manifest for the scDRS pipeline.
#
# Traits are selected from /dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain
# (see trait_names_keys in that directory). Scope: brain + substance use
# traits relevant to NAc, plus non-brain negative controls to check
# specificity of domain / cell type associations.
#
# Usage: module load conda_R/4.5; Rscript 01_trait_manifest.R

library(here)

gwas_dir <- "/dcs04/lieber/shared/statsgen/LDSC/base/gwas_brain"

manifest <- tibble::tribble(
  ~trait_key,      ~trait_name,                    ~category,        ~file,
  # -- Psychiatric --------------------------------------------------------
  "SCZ",           "Schizophrenia (PGC3)",         "psychiatric",    "scz_PGC3_ldscore.gz",
  "MDD",           "Major depression (ex-23andMe)","psychiatric",    "mdd_ex23andMe_ldscore.gz",
  "BIP",           "Bipolar disorder (PGC3)",      "psychiatric",    "bp3_ldscore.gz",
  "ASD",           "Autism spectrum disorder",     "psychiatric",    "autism_ldscore.gz",
  "ADHD",          "ADHD",                         "psychiatric",    "adhd.gz",
  "PTSD",          "PTSD",                         "psychiatric",    "PTSD.gz",
  "AN",            "Anorexia nervosa",             "psychiatric",    "Anorexia_ldscore.gz",
  "INSOMNIA",      "Insomnia",                     "psychiatric",    "insomnia.gz",
  "NEUROT",        "Neuroticism (UKB)",            "psychiatric",    "UKB_460K.mental_NEUROTICISM.sumstats.gz",
  # -- Cognitive ----------------------------------------------------------
  "IQ",            "Intelligence",                 "cognitive",      "intelligence.gz",
  "EDU",           "Educational attainment (UKB)", "cognitive",      "UKB_460K.cov_EDU_YEARS.sumstats.gz",
  # -- Substance use ------------------------------------------------------
  "OUD",           "Opioid use disorder (MVP1+MVP2+YP+SAGE)", "substance_use", "OUD_EA_MVP1_MVP2_YP_SAGE_Mar12.gz",
  "AUD",           "Alcohol use disorder (EUR meta 2023)",    "substance_use", "AUD.EUR_META.NatMed2023.gz",
  "AGE_SMK",       "Age of smoking initiation (GSCAN)",       "substance_use", "GSCAN_AgeSmk_2022_ldscore.gz",
  "CIG_DAY",       "Cigarettes per day (GSCAN)",              "substance_use", "GSCAN_CigDay_2022_ldscore.gz",
  "DRINKS_WK",     "Drinks per week (GSCAN)",                 "substance_use", "GSCAN_DrnkWk_2022_ldscore.gz",
  "SMK_CES",       "Smoking cessation (GSCAN)",               "substance_use", "GSCAN_SmkCes_2022_ldscore.gz",
  "SMK_INIT",      "Smoking initiation (GSCAN)",              "substance_use", "GSCAN_SmkInit_2022_ldscore.gz",
  # -- Neurological -------------------------------------------------------
  "PD",            "Parkinson disease",            "neurological",   "PD_ldscore.gz",
  "AD",            "Alzheimer disease",            "neurological",   "Alzheimer_ldscore.gz",
  "EPI_ALL",       "Epilepsy (all)",               "neurological",   "epilepsyAll.gz",
  "EPI_FOCAL",     "Epilepsy (focal)",             "neurological",   "epilepsyFocal.gz",
  "EPI_GGE",       "Epilepsy (GGE)",               "neurological",   "epilepsyGGE.gz",
  # -- Negative controls (non-brain) ---------------------------------------
  "HEIGHT",        "Height",                       "negative_control", "height_ldscore.gz",
  "BMI",           "Body mass index",              "negative_control", "bmi_ldscore.gz",
  "T2D",           "Type 2 diabetes",              "negative_control", "PASS_Type_2_Diabetes.sumstats.gz"
)

manifest$file_path <- file.path(gwas_dir, manifest$file)

missing <- manifest$file_path[!file.exists(manifest$file_path)]
if (length(missing) > 0) {
  stop("Missing sumstats files:\n", paste(missing, collapse = "\n"))
}

out <- here("processed-data", "16_scDRS", "trait_manifest.tsv")
write.table(manifest, out, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote ", nrow(manifest), " traits to ", out)

sessioninfo::session_info()
