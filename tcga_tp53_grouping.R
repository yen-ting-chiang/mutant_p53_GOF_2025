setwd("C:/Users/danny/Documents/R_project/mutant_p53_GOF_2025")
getwd()

# 載入必要套件
library(TCGAbiolinks)
library(dplyr)
library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(ggfortify)
library(data.table)
library(shiny)
library(scales)
library(utils)
library(tidyr)

# --- Shiny app for cancer type selection ---
ui <- fluidPage(
  titlePanel("TCGA Cancer Type Selector"),
  sidebarLayout(
    sidebarPanel(
      selectInput("cancer_type", "Select TCGA Project:", 
                  choices = grep("^TCGA", getGDCprojects()$project_id, value = TRUE),
                  selected = "TCGA-BRCA"),
      checkboxInput("download_tpm", "Also download TPM expression data", value = TRUE)
    ),
    mainPanel(
      textOutput("selected")
    )
  )
)

server <- function(input, output, session) {
  output$selected <- renderText({
    paste("Selected TCGA Project:", input$cancer_type)
  })
}

shinyApp(ui, server)

# --- 分析腳本開始 ---
cancer_type <- "TCGA-GBM"

# 下載 mutation data (MAF)
maf_query <- GDCquery(project = cancer_type,
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query, method = "api")
maf <- GDCprepare(maf_query)

# 自動下載 UCSC Xena 提供的 gene-level GISTIC2 CNV 檔案
xena_cnv_url <- "https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"
cnv_file_path <- "GISTIC2_all_data_by_genes"
if (!file.exists(cnv_file_path)) {
  download.file(xena_cnv_url, destfile = cnv_file_path, mode = "wb")
}

# 若尚未解壓縮，先解壓縮檔案（Windows 無法用 zcat）
unzipped_path <- sub(".gz$", "", cnv_file_path)
if (!file.exists(unzipped_path)) {
  R.utils::gunzip(cnv_file_path, destname = unzipped_path, overwrite = TRUE)
}


# 載入 CNV 檔案（以行轉置，取得 TP53 CNV）
gistic_cnv <- fread(cnv_file_path, data.table = FALSE)
rownames(gistic_cnv) <- gistic_cnv$Sample
gistic_cnv <- gistic_cnv[, -1]
gistic_cnv_t <- as.data.frame(t(gistic_cnv))
colnames(gistic_cnv_t) <- make.names(colnames(gistic_cnv_t))
tp53_cnv_long <- data.frame(
  sample_id = rownames(gistic_cnv_t),
  gistic_score = gistic_cnv_t[["TP53"]],
  stringsAsFactors = FALSE
)
tp53_cnv_long$sample_id <- gsub("\\.", "-", tp53_cnv_long$sample_id)
saveRDS(tp53_cnv_long, file = "tp53_cnv_long.rds")
tp53_cnv_long <- readRDS("tp53_cnv_long.rds")

# 安裝 UCSCXenaShiny 套件（如尚未安裝）
if (!requireNamespace("UCSCXenaShiny", quietly = TRUE)) {
  install.packages("UCSCXenaShiny")
}

library(UCSCXenaShiny)
library(dplyr)

# 載入內建的 TCGA purity 資料
data("tcga_purity")

# 假設你要篩選 TCGA-GBM 的資料
purity_df <- tcga_purity %>%
  filter(cancer_type == "GBM") %>%
  select(sample, CPE) %>%
  rename(sample_id = sample, purity = CPE)

head(purity_df)


# 依 cancer_type 過濾樣本（保留當前 cancer 專案樣本）
sample_prefix <- gsub("TCGA-", "", cancer_type)
valid_samples <- maf$Tumor_Sample_Barcode[grep(sample_prefix, maf$Tumor_Sample_Barcode)]
tp53_cnv_long <- tp53_cnv_long %>% filter(sample_id %in% valid_samples)
purity_df <- purity_df %>% filter(sample_id %in% valid_samples)

# 整合 mutation 和 CNV 並考慮 purity 進行分組，分類邏輯依據雙 allele 組合表
core_gof <- c("p.R175H", "p.R248Q", "p.R273C", "p.R273H", "p.R248W", "p.R282W",
              "p.Y220C", "p.G245S", "p.V157F", "p.Y163C", "p.R249S", "p.C176F",
              "p.R158L", "p.Y234C", "p.R280K", "p.G245D", "p.P151S", "p.R280T",
              "p.R158H", "p.L194R")

classify_allele <- Vectorize(function(variant_class, protein_change, cnv_score) {
  if (!is.na(cnv_score) && cnv_score <= -1) {
    return("CNV loss")
  } else if (variant_class %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site")) {
    return("LOF")
  } else if (protein_change %in% core_gof) {
    return("Core GOF")
  } else if (variant_class == "Missense_Mutation") {
    return("Extended GOF")
  } else {
    return("WT")
  }
})

group_samples <- function(maf_df, cnv_df, purity_df) {
  sample_groups <- list()
  allele_info <- list()
  all_samples <- intersect(unique(maf_df$Tumor_Sample_Barcode), purity_df$sample_id)
  
  for (sid in all_samples) {
    purity <- purity_df %>% filter(sample_id == sid) %>% pull(purity)
    if (length(purity) == 0 || is.na(purity) || purity < 0.5) next
    
    muts <- maf_df %>% filter(Tumor_Sample_Barcode == sid, Hugo_Symbol == "TP53")
    cnv_score <- cnv_df %>% filter(sample_id == sid) %>% pull(gistic_score)
    if (length(cnv_score) == 0) cnv_score <- NA
    
    allele_states <- character()
    if (nrow(muts) > 0) {
      allele_states <- classify_allele(muts$Variant_Classification, muts$HGVSp_Short, NA)
    }
    if (length(allele_states) < 2) {
      allele_states <- c(allele_states, classify_allele("", "", cnv_score))
    }
    allele_states <- allele_states[1:2]
    
    combo <- paste(sort(allele_states), collapse = "|")
    group <- switch(combo,
                    "WT|WT" = "WTP53",
                    "WT|LOF" = "LOFP53",
                    "WT|CNV loss" = "LOFP53",
                    "WT|Core GOF" = "Core_GOFP53",
                    "WT|Extended GOF" = "Extended_GOFP53",
                    "LOF|LOF" = "LOFP53",
                    "LOF|CNV loss" = "LOFP53",
                    "LOF|Core GOF" = "Core_GOFP53",
                    "LOF|Extended GOF" = "Extended_GOFP53",
                    "Core GOF|Core GOF" = "Core_GOFP53",
                    "Core GOF|Extended GOF" = "Core_GOFP53",
                    "Core GOF|CNV loss" = "Core_GOFP53",
                    "Extended GOF|Extended GOF" = "Extended_GOFP53",
                    "Extended GOF|CNV loss" = "Extended_GOFP53",
                    "CNV loss|CNV loss" = "LOFP53",
                    "CNV loss|WT" = "LOFP53",
                    "CNV loss|Core GOF" = "Core_GOFP53",
                    "CNV loss|Extended GOF" = "Extended_GOFP53",
                    "Unclassified")
    
    sample_groups[[sid]] <- group
    allele_info[[sid]] <- paste(allele_states, collapse = "+")
  }
  
  return(data.frame(sample_id = names(sample_groups), 
                    group = unlist(sample_groups),
                    allele_states = unlist(allele_info)))
}

sample_group_df <- group_samples(maf, tp53_cnv_long, purity_df)

# 顯示每組 sample 數量
sample_group_counts <- sample_group_df %>% count(group)
print(sample_group_counts)

# 以下 RNA-seq / RPPA 等分析照常執行…（略）