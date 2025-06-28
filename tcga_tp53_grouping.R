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
cancer_type <- "TCGA-BRCA"

# 下載 mutation data (MAF)
maf_query <- GDCquery(project = cancer_type,
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query, method = "api")
maf <- GDCprepare(maf_query)

# 自動下載 UCSC Xena 提供的 gene-level GISTIC2 CNV 檔案
xena_cnv_url <- "https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"
cnv_file_path <- "GISTIC2_all_data_by_genes.gz"
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


# 整合 mutation 和 CNV 以進行 sample 分組（考慮多突變、VAF與 purity）
group_samples <- function(maf_df, cnv_df, vaf_threshold = 0.3) {
  sample_groups <- list()
  all_samples <- unique(c(maf_df$Tumor_Sample_Barcode, cnv_df$sample_id))
  
  core_gof <- c("p.R175H", "p.R248Q", "p.R273C", "p.R273H", "p.R248W", "p.R282W",
                "p.Y220C", "p.G245S", "p.V157F", "p.Y163C", "p.R249S", "p.C176F",
                "p.R158L", "p.Y234C", "p.R280K", "p.G245D", "p.P151S", "p.R280T",
                "p.R158H", "p.L194R")
  
  classify_allele <- Vectorize(function(variant_class, protein_change, cnv_score) {
    if (!is.na(cnv_score) && cnv_score <= -1) {
      return("lof")
    } else if (variant_class %in% c("Frame_Shift_Del", "Frame_Shift_Ins", "Nonsense_Mutation", "Splice_Site")) {
      return("lof")
    } else if (protein_change %in% core_gof) {
      return("core_gof")
    } else if (variant_class == "Missense_Mutation") {
      return("extended_gof")
    } else {
      return("wt")
    }
  })
  
  for (sid in all_samples) {
    muts <- maf_df %>% 
      filter(Tumor_Sample_Barcode == sid, Hugo_Symbol == "TP53") %>% 
      filter(!is.na(t_alt_count) & !is.na(t_ref_count)) %>%
      mutate(vaf = t_alt_count / (t_alt_count + t_ref_count)) %>%
      filter(vaf >= vaf_threshold)
    
    cnv_score <- cnv_df %>% filter(sample_id == sid) %>% pull(gistic_score)
    if (length(cnv_score) == 0) cnv_score <- NA
    
    allele_states <- c()
    if (nrow(muts) >= 1) {
      allele_states <- classify_allele(muts$Variant_Classification, muts$HGVSp_Short, NA)
      allele_states <- allele_states[1:min(2, length(allele_states))]  # 最多兩個
    }
    
    if (length(allele_states) < 2) {
      allele_states <- c(allele_states, classify_allele("", "", cnv_score))
    }
    allele_states <- allele_states[1:2]
    
    if ("core_gof" %in% allele_states) {
      sample_groups[[sid]] <- "Core_GOFP53"
    } else if ("extended_gof" %in% allele_states) {
      sample_groups[[sid]] <- "Extended_GOFP53"
    } else if (all(allele_states %in% c("lof"))) {
      sample_groups[[sid]] <- "LOFP53"
    } else if ("lof" %in% allele_states) {
      sample_groups[[sid]] <- "LOFP53"
    } else {
      sample_groups[[sid]] <- "WTP53"
    }
  }
  
  return(data.frame(sample_id = names(sample_groups), group = unlist(sample_groups)))
}

sample_group_df <- group_samples(maf, tp53_cnv_long)

# 下載 RNA-seq expression 資料（Counts）
exp_query_counts <- GDCquery(project = cancer_type,
                             data.category = "Transcriptome Profiling",
                             data.type = "Gene Expression Quantification",
                             workflow.type = "HTSeq - Counts")
GDCdownload(exp_query_counts, method = "api")
exp_data <- GDCprepare(exp_query_counts)

# 若選擇，同時下載 TPM 資料供其他分析用途
if (TRUE) {
  exp_query_tpm <- GDCquery(project = cancer_type,
                            data.category = "Transcriptome Profiling",
                            data.type = "Gene Expression Quantification",
                            workflow.type = "HTSeq - TPM")
  GDCdownload(exp_query_tpm, method = "api")
  exp_data_tpm <- GDCprepare(exp_query_tpm)
}

# 下載 protein expression (RPPA) 資料
protein_query <- GDCquery(project = cancer_type,
                          data.category = "Protein Expression",
                          data.type = "Protein Expression Quantification")
GDCdownload(protein_query, method = "api")
protein_data <- GDCprepare(protein_query)

# TP53 protein Z-score 計算與整合
protein_df <- assay(protein_data) %>% as.data.frame()
protein_df <- t(scale(t(protein_df)))
protein_df <- as.data.frame(protein_df)
protein_df$Gene <- rownames(protein_df)
protein_tp53 <- protein_df %>% filter(Gene == "TP53") %>% select(-Gene)
protein_tp53_long <- data.frame(sample_id = colnames(protein_tp53), TP53_protein_z = as.numeric(protein_tp53[1,]))

# 整合 protein z-score 資訊回 sample 分組表
sample_group_df <- sample_group_df %>% left_join(protein_tp53_long, by = "sample_id")
