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
library(tibble)
library(UCSCXenaShiny)

# --- 分析腳本開始 ---
cancer_type <- "TCGA-GBM"

# 下載 mutation data (MAF)
maf_query <- GDCquery(project = cancer_type,
                      data.category = "Simple Nucleotide Variation",
                      data.type = "Masked Somatic Mutation",
                      workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking")
GDCdownload(maf_query, method = "api")
maf <- GDCprepare(maf_query)
saveRDS(maf, file = "maf_TCGA_GBM.rds")
maf <- readRDS("maf_TCGA_GBM.rds")


# 下載並解壓 GISTIC2 CNV
xena_cnv_url <- "https://tcga.xenahubs.net/download/TCGA.PANCAN.sampleMap/Gistic2_CopyNumber_Gistic2_all_data_by_genes.gz"
cnv_file_path <- "GISTIC2_all_data_by_genes.gz"
if (!file.exists(cnv_file_path)) download.file(xena_cnv_url, destfile = cnv_file_path, mode = "wb")

unzipped_path <- sub(".gz$", "", cnv_file_path)
if (!file.exists(unzipped_path)) R.utils::gunzip(cnv_file_path, destname = unzipped_path, overwrite = TRUE)

# 讀入 CNV 資料
gistic_cnv <- fread(unzipped_path, data.table = FALSE)
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

# 載入 purity 資料
data("tcga_purity")
purity_df <- tcga_purity %>%
  dplyr::filter(cancer_type == "GBM") %>%
  dplyr::select(sample, CPE) %>%
  dplyr::rename(sample_id = sample, purity = CPE)
purity_df <- purity_df %>%
  filter(!is.na(purity))

# 對 maf 建立簡化 sample ID 欄位
maf$sample_id <- substr(maf$Tumor_Sample_Barcode, 1, 15)

# 確認 CNV 和 purity 的 sample_id 也都是 15 碼
table(nchar(tp53_cnv_long$sample_id))  # 應為 15
table(nchar(purity_df$sample_id))      # 應為 15

valid_samples <- Reduce(intersect, list(
  unique(maf$sample_id),
  unique(tp53_cnv_long$sample_id),
  unique(purity_df$sample_id)
))

maf_filtered <- maf %>% dplyr::filter(sample_id %in% valid_samples)
cnv_filtered <- tp53_cnv_long %>% dplyr::filter(sample_id %in% valid_samples)
purity_filtered <- purity_df %>% dplyr::filter(sample_id %in% valid_samples)

library(stringr)
# 定義 TP53 突變分類
core_gof <- c(
  "p.R175H", "p.R248Q", "p.R273C", "p.R273H", "p.R248W", "p.R282W",
  "p.Y220C", "p.G245S", "p.V157F", "p.Y163C", "p.R249S", "p.C176F",
  "p.R158L", "p.Y234C", "p.R280K", "p.G245D", "p.P151S", "p.R280T",
  "p.R158H", "p.L194R"
)
lof_class <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site", "In_Frame_Del")

# 提取 TP53 突變
tp53_mut <- maf %>%
  dplyr::filter(Hugo_Symbol == "TP53") %>%
  dplyr::mutate(
    sample_id = str_sub(Tumor_Sample_Barcode, 1, 15),
    mutation_class = case_when(
      HGVSp_Short %in% core_gof ~ "Core GOF",
      Variant_Classification %in% lof_class ~ "LOF",
      Variant_Classification == "Missense_Mutation" ~ "Extended GOF",
      TRUE ~ "Other"
    )
  ) %>%
  dplyr::select(sample_id, mutation_class, HGVSp_Short)

# 計算 allele frequency 並保留必要欄位
tp53_mut <- maf %>%
  filter(Hugo_Symbol == "TP53") %>%
  mutate(
    sample_id = substr(Tumor_Sample_Barcode, 1, 15),
    AF = t_alt_count / (t_alt_count + t_ref_count),
    mutation_class = case_when(
      HGVSp_Short %in% core_gof ~ "Core GOF",
      Variant_Classification %in% lof_class ~ "LOF",
      Variant_Classification == "Missense_Mutation" ~ "Extended GOF",
      TRUE ~ "Other"
    )
  ) %>%
  filter(!is.na(HGVSp_Short)) # 避免 NA 影響 group_by


# 設定 AF threshold
AF_threshold <- 0.1  # 可調整為 0.1 或更嚴格

tp53_mut_filtered <- tp53_mut %>%
  filter(AF >= AF_threshold)

# Step 1: 每個 sample_id + HGVSp_Short 保留 AF 最大的一筆
tp53_dedup <- tp53_mut_filtered %>%
  group_by(sample_id, HGVSp_Short) %>%
  slice_max(order_by = AF, n = 1) %>%
  ungroup()

# Step 2: 每個 sample_id 取 AF 前兩名
tp53_summary <- tp53_dedup %>%
  arrange(sample_id, desc(AF)) %>%
  group_by(sample_id) %>%
  summarise(
    allele_1 = mutation_class[1],
    allele_2 = mutation_class[2],
    mutation_1 = HGVSp_Short[1],
    mutation_2 = HGVSp_Short[2],
    af_1 = AF[1],
    af_2 = AF[2],
    .groups = "drop"
  )



# 所有 sample_id（保留 WT）
all_samples <- data.frame(sample_id = unique(str_sub(maf$Tumor_Sample_Barcode, 1, 15)))

# 合併突變分類
sample_mut_df <- all_samples %>%
  left_join(tp53_summary, by = "sample_id") %>%
  mutate(
    allele_1 = ifelse(is.na(allele_1), "WT", allele_1),
    allele_2 = ifelse(is.na(allele_2), NA, allele_2),
    mutation_1 = ifelse(is.na(mutation_1), NA, mutation_1),
    mutation_2 = ifelse(is.na(mutation_2), NA, mutation_2)
  )

# 分類優先順序（Core GOF > Extended GOF > LOF > CNV loss > WT）
priority_rank <- c("Core GOF" = 1, "Extended GOF" = 2, "LOF" = 3, "Other" = 4, "WT" = 5)

# 依據 priority 判定 class
sample_mut_df <- sample_mut_df %>%
  dplyr::mutate(
    rank_1 = priority_rank[allele_1],
    rank_2 = priority_rank[allele_2],
    best_rank = pmin(rank_1, rank_2, na.rm = TRUE),
    class = names(priority_rank)[match(best_rank, priority_rank)]
  ) %>%
  dplyr::select(-rank_1, -rank_2, -best_rank)

# 加入 CNV 資料
sample_mut_df <- sample_mut_df %>%
  left_join(tp53_cnv_long, by = "sample_id") %>%
  dplyr::mutate(
    class = ifelse(class == "WT" & gistic_score <= -0.5, "CNV loss", class)
  )

# 加入 purity 資料
purity_df_unique <- purity_df %>%
  group_by(sample_id) %>%
  summarise(purity = mean(purity, na.rm = TRUE), .groups = "drop")

sample_class_df <- sample_mut_df %>%
  left_join(purity_df_unique, by = "sample_id")

# 設定 purity threshold ≥ 0.6
sample_class_df_filtered <- sample_class_df %>%
  filter(!is.na(purity) & purity >= 0.6)

# 統計每個 TP53 類別的樣本數量
class_counts <- sample_class_df_filtered %>%
  group_by(class) %>%
  summarise(sample_count = n()) %>%
  arrange(desc(sample_count))

# 顯示結果
print(class_counts)


# 下載 RNA-seq expression data
rna_query <- GDCquery(project = cancer_type,
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification",
                      workflow.type = "STAR - Counts")

GDCdownload(rna_query, method = "api")
rna_data <- GDCprepare(rna_query)
saveRDS(rna_data, "rna_data_TCGA_GBM.rds")
rna_data <- readRDS("rna_data_TCGA_GBM.rds")

# 提取 raw counts matrix
counts_matrix <- assay(rna_data, "unstranded")
colnames(counts_matrix) <- substr(colnames(counts_matrix), 1, 15)  # simplify sample ID

# 篩選出與 sample_class_df_filtered 有對應的 sample
valid_sample_ids <- sample_class_df_filtered$sample_id
counts_filtered <- counts_matrix[, colnames(counts_matrix) %in% valid_sample_ids]

# 過濾掉沒有對應 RNA-seq 的樣本
sample_class_matched <- sample_class_df_filtered %>%
  filter(sample_id %in% colnames(counts_filtered)) %>%
  filter(!is.na(class)) %>%  
  distinct(sample_id, class)

# 確保順序一致
counts_filtered <- counts_filtered[, sample_class_matched$sample_id]
stopifnot(all(colnames(counts_filtered) == sample_class_matched$sample_id))

# 準備 DESeq2 對象
coldata <- data.frame(row.names = sample_class_matched$sample_id,
                      condition = sample_class_matched$class)
coldata$condition <- gsub(" ", "_", coldata$condition)
coldata$condition <- factor(coldata$condition, levels = c("WT", "CNV_loss", "LOF", "Core_GOF", "Extended_GOF", "Other"))

dds <- DESeqDataSetFromMatrix(countData = counts_filtered,
                              colData = coldata,
                              design = ~ condition)

# 過濾低表現基因
dds <- dds[rowSums(counts(dds)) > 10, ]

# variance stabilizing transformation (VST)
vsd <- vst(dds, blind = TRUE)

# 繪製 PCA 圖
autoplot(prcomp(t(assay(vsd))), data = coldata, colour = "condition") +
  ggtitle("PCA of TCGA GBM Samples by TP53 Mutation Class") +
  theme_minimal()


coldata_filtered <- coldata %>%
  filter(condition %in% c("WT", "LOF", "Core_GOF", "Extended_GOF"))

# 同步過濾 counts
counts_filtered_subset <- counts_filtered[, rownames(coldata_filtered)]

vsd_subset <- vst(DESeqDataSetFromMatrix(counts_filtered_subset, coldata_filtered, ~ condition))
autoplot(prcomp(t(assay(vsd_subset))), data = coldata_filtered, colour = "condition") +
  theme_minimal()

install.packages("Rtsne")
install.packages("uwot")

library(Rtsne)

# 取 VST 結果並轉置成 sample x gene 的矩陣
tsne_input <- t(assay(vsd_subset))

# 建立 t-SNE
set.seed(42)
tsne_result <- Rtsne(tsne_input, perplexity = 10, verbose = TRUE)

# 整理 t-SNE 結果
tsne_df <- as.data.frame(tsne_result$Y)
tsne_df$sample_id <- rownames(tsne_input)
tsne_df$condition <- coldata_filtered[tsne_df$sample_id, "condition"]

# 繪圖
ggplot(tsne_df, aes(x = V1, y = V2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE of TCGA GBM by TP53 Class") +
  theme_minimal()

library(uwot)

# 建立 UMAP
set.seed(42)
umap_result <- umap(tsne_input)

# 整理結果
umap_df <- as.data.frame(umap_result)
umap_df$sample_id <- rownames(tsne_input)
umap_df$condition <- coldata_filtered[umap_df$sample_id, "condition"]

# 繪圖
ggplot(umap_df, aes(x = V1, y = V2, color = condition)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMAP of TCGA GBM by TP53 Class") +
  theme_minimal()

# 從 MAF 中找出有 IDH1 mutation 的樣本
idh1_status <- maf %>%
  filter(Hugo_Symbol == "IDH1") %>%
  mutate(sample_id = substr(Tumor_Sample_Barcode, 1, 15)) %>%
  distinct(sample_id) %>%
  mutate(idh1_class = "IDH1_MUT")

# 合併進 t-SNE 資料
tsne_df <- tsne_df %>%
  left_join(idh1_status, by = "sample_id") %>%
  mutate(idh1_class = ifelse(is.na(idh1_class), "IDH1_WT", idh1_class))

# 合併進 UMAP 資料
umap_df <- umap_df %>%
  left_join(idh1_status, by = "sample_id") %>%
  mutate(idh1_class = ifelse(is.na(idh1_class), "IDH1_WT", idh1_class))


# t-SNE 上色
ggplot(tsne_df, aes(x = V1, y = V2, color = idh1_class)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "t-SNE of TCGA GBM by IDH1 Mutation Status") +
  theme_minimal()

# UMAP 上色
ggplot(umap_df, aes(x = V1, y = V2, color = idh1_class)) +
  geom_point(size = 3, alpha = 0.8) +
  labs(title = "UMAP of TCGA GBM by IDH1 Mutation Status") +
  theme_minimal()


# 使用 TCGAbiolinks 內建功能抓 GBM 分型資訊
gbm_subtype <- TCGAquery_subtype("GBM")

# 處理 sample ID
gbm_subtype$sample_id <- substr(gbm_subtype$patient, 1, 15)

# 加入到 t-SNE dataframe
tsne_df <- tsne_df %>%
  left_join(gbm_subtype[, c("sample_id", "gene_expression_subtype")], by = "sample_id") %>%
  rename(subtype = gene_expression_subtype) %>%
  mutate(subtype = ifelse(is.na(subtype), "Unknown", subtype))

# 加入到 UMAP dataframe
umap_df <- umap_df %>%
  left_join(gbm_subtype[, c("sample_id", "subtype")], by = "sample_id") %>%
  mutate(subtype = ifelse(is.na(subtype), "Unknown", subtype))








# 加入 IDH1 突變標註後進行 stratified PCA

# --- 步驟 1：從 MAF 中提取 IDH1 mutation status ---
idh1_status <- maf %>%
  filter(Hugo_Symbol == "IDH1") %>%
  mutate(sample_id = substr(Tumor_Sample_Barcode, 1, 15)) %>%
  distinct(sample_id, Variant_Classification) %>%
  mutate(idh1_status = ifelse(Variant_Classification == "Missense_Mutation", "IDH1_MUT", "IDH1_WT")) %>%
  dplyr::select(sample_id, idh1_status)

# --- 步驟 2：合併 IDH1 資訊到 sample_class_df_filtered 中 ---
sample_class_with_idh1 <- sample_class_df_filtered %>%
  left_join(idh1_status, by = "sample_id") %>%
  dplyr::mutate(idh1_status = ifelse(is.na(idh1_status), "IDH1_WT", idh1_status))

# --- 步驟 3：對每個 IDH1 分群單獨做 PCA ---
for (idh_group in c("IDH1_WT", "IDH1_MUT")) {
  cat("\n\n===== Analyzing", idh_group, "group =====\n")
  
  # 篩選目前 group
  sample_subset <- sample_class_with_idh1 %>%
    filter(idh1_status == idh_group,
           class %in% c("WT", "LOF", "Core GOF", "Extended GOF")) %>%
    distinct(sample_id, class)
  
  # 取出 RNA-seq matrix 並同步 sample 順序
  counts_sub <- counts_matrix[, colnames(counts_matrix) %in% sample_subset$sample_id]
  sample_subset <- sample_subset %>% filter(sample_id %in% colnames(counts_sub))
  counts_sub <- counts_sub[, sample_subset$sample_id]
  stopifnot(all(colnames(counts_sub) == sample_subset$sample_id))
  
  # 建立 coldata 並 factor 化
  coldata_sub <- data.frame(row.names = sample_subset$sample_id,
                            condition = sample_subset$class)
  coldata_sub$condition <- gsub(" ", "_", coldata_sub$condition)
  coldata_sub$condition <- factor(coldata_sub$condition,
                                  levels = c("WT", "LOF", "Core_GOF", "Extended_GOF"))
  
  # 建立 DESeq2 dataset 並 vst
  dds_sub <- DESeqDataSetFromMatrix(countData = counts_sub,
                                    colData = coldata_sub,
                                    design = ~ condition)
  dds_sub <- dds_sub[rowSums(counts(dds_sub)) > 10, ]
  vsd_sub <- vst(dds_sub, blind = TRUE)
  
  # 畫圖
  p <- autoplot(prcomp(t(assay(vsd_sub))), data = coldata_sub, colour = "condition") +
    ggtitle(paste("PCA of", idh_group, "Samples by TP53 Class")) +
    theme_minimal()
  print(p)
}
