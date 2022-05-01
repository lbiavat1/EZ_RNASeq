rm(list = ls())


# tidyverse-friendly
library(tidyverse)
library(tidybulk)
library(tidyHeatmap)
library(ComplexHeatmap)
library(ggrepel)
library(plotly)
library(gg3D)
library(GGally)

# other useful libraries
library(DESeq2)
library(limma)
library(edgeR)
library(fgsea)
library(IHW)


library(annotables)

# Set PrimaryDirectory where this script is located
dirname(rstudioapi::getActiveDocumentContext()$path)  
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
getwd()
PrimaryDirectory <- getwd()
PrimaryDirectory

# create working directory
workingDir <- "WorkingDirectory_D0"
dirPath <- file.path(PrimaryDirectory, workingDir)
dir.create(dirPath)
setwd(dirPath)

saveDir <- file.path(dirPath, "results")
savePath <- saveDir
dir.create(saveDir)

filesPath <- file.path(dirPath, "files")
dir.create(filesPath)

plotDir <- file.path(saveDir, "plots")
dir.create(plotDir)

raw_cts <- file.path(PrimaryDirectory, "RNASeq_rawcounts")
files <- list.files(raw_cts)

path.to.files <- file.path(PrimaryDirectory, "RNASeq_rawcounts", files)
list.of.tibbles <- map(path.to.files, read_csv)
myTibble <- list.of.tibbles %>% purrr::reduce(., inner_join)

mtb <- map(file.path(PrimaryDirectory, "RNASeq_rawcounts", files), read_csv) %>% 
  purrr::reduce(., inner_join)
mtb

nfiles <- length(list.files(raw_cts))

t <- read_csv(file.path(PrimaryDirectory, "RNASeq_rawcounts",
                        files[1]))

myTibble <- t
for(i in 2:nfiles){
  myTibble <- inner_join(myTibble, 
                        read_csv(file.path(PrimaryDirectory, "RNASeq_rawcounts",
                                           files[i])))
}
all.equal(mtb, myTibble)
if(all.equal(mtb, myTibble)){
  rm(myTibble)
}

mtb

names(mtb)[1] <- "feature"

grcm38
grcm38_tx2gene

counts_tt <- pivot_longer(mtb, cols = c(2:18),
                          names_to = "sample", values_to = "counts") %>%
  dplyr::inner_join(grcm38, by = c("feature" = "ensgene")) %>%
  dplyr::select(symbol, sample, counts) %>%
  dplyr::rename(., feature = symbol)
counts_tt

  
"Tcf7" %in% counts_tt$feature
"Gzmb" %in% counts_tt$feature
"Elane" %in% counts_tt$feature
"Fcer1g" %in% counts_tt$feature
counts_tt

########################## prep for data analysis #############################



# myInfo <- unique(counts$sample)
# getwd()
# write_csv(as.data.frame(myInfo), "experiment_info.csv")
# myInfo <- read_csv("experiment_info.csv")
# myInfo <- myInfo %>%
#   mutate(mouse = paste("#", mouse, sep = ""))
# write_csv(myInfo, "experiment_info.csv")

myInfo <- read_csv("experiment_info.csv")
myInfo <- myInfo %>%
  mutate(mouse = gsub("#", "", mouse))
myInfo

counts_tt <- counts_tt %>% left_join(myInfo)
counts_tt
counts <- counts_tt %>%
  tidybulk(.sample = sample, .transcript = feature, .abundance = counts) %>%
  aggregate_duplicates()
counts
counts
ggplot(counts_tt, aes(x = sample, weight = counts, fill = sampleName)) +
  geom_bar() +
  theme_bw()

counts %>% distinct(paste(feature, sample, sep = "_"), .keep_all = TRUE) %>%
  dplyr::select(feature, sample, counts, tissue, cell.type, timepoint, mouse, sampleName)

counts_scaled <- counts %>% distinct(paste(feature, sample, sep = "_"), .keep_all = TRUE) %>%
  dplyr::select(feature, sample, counts, tissue, cell.type, timepoint, mouse, sampleName) %>%
  identify_abundant(factor_of_interest = sampleName, minimum_counts = 100, minimum_proportion = 0.50) %>%
  scale_abundance(method = "TMM")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = abundance + 1, color = sample)) +
  geom_density() +
  facet_wrap(~source) +
  scale_x_log10() +
  theme_bw()
# ggsave(file.path(plotDir, "counts_scaled.pdf"), device = "pdf")

counts_scaled %>%
  filter(.abundant) %>%
  pivot_longer(cols = c("counts", "counts_scaled"), names_to = "source", values_to = "abundance") %>%
  ggplot(aes(x = sample, y = abundance + 1, fill = sample)) +
  geom_boxplot() +
  geom_hline(aes(yintercept = median(abundance + 1)), colour="red") +
  facet_wrap(~source) +
  scale_y_log10() +
  theme_bw()
# ggsave(file.path(plotDir, "counts_scaled_boxplot.pdf"), device = "pdf")

counts_scaled %>% group_by(sample) %>%
  summarise(total_reads=sum(counts))

ggplot(counts_scaled, mapping = aes(x = sample, weight = counts, fill = sample)) +
  geom_bar() +
  theme(axis.text.x = element_blank())
# ggsave(file.path(plotDir, "count_reads_per_sample.pdf"), device = "pdf")

################ subset dataset for cell.type - tissue comparison #############

TPEX <- counts_scaled %>% 
  dplyr::filter(cell.type == "TPEX") %>%
  dplyr::filter(timepoint == "late")
TPEX

TEFF <- counts_scaled %>%
  dplyr::filter(cell.type != "TPEX") %>%
  dplyr:: filter(timepoint == "late")
TEFF

################## PCA ###################################
TPEX_PCA <- TPEX %>%
  reduce_dimensions(method = "PCA", top = 500)
attr(TPEX_PCA, "internals")$PCA

TPEX_PCA %>%
  # mutate(tis_cel = paste(tissue, cell.type, sep = "_")) %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = tissue)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  # stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
# ggsave(file.path(plotDir, "PCA_top500_TPEX.pdf"), device = "pdf")

TEFF_PCA <- TEFF %>%
  reduce_dimensions(method = "PCA", top = 500)
attr(TPEX_PCA, "internals")$PCA

TEFF_PCA %>%
  # mutate(tis_cel = paste(tissue, cell.type, sep = "_")) %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = tissue)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  # stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
# ggsave(file.path(plotDir, "PCA_top500_TEFF.pdf"), device = "pdf")

###################### heatmap ################################################

hm <- TPEX %>%

  # filter lowly abundant
  filter(.abundant) %>%

  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%

  as_tibble() %>%

  mutate(genes = feature) %>%

  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    # row_names_gp = gpar(fontsize = 4),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 2,
    column_km_repeats = 100,
    row_km = 2,
    row_km_repeats = 500,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A",
                                       "#33A02C", "#FB9A99"),
                              font = c(1,2,3))
  ) %>%
  add_tile(c(tissue))
hm
# pdf(file = file.path(plotDir, "heatmap_top500_TPEX-late.pdf"))
# hm
# dev.off()

TEFF

hm <- TEFF %>%
  
  filter(cell.type == "TEFF") %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    # row_names_gp = gpar(fontsize = 4),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 2,
    column_km_repeats = 100,
    row_km = 4,
    row_km_repeats = 500,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A",
                                       "#33A02C", "#FB9A99"), 
                              font = c(1,2,3))
  ) %>%
  add_tile(c(tissue))
hm
# pdf(file = file.path(plotDir, "heatmap_top500_TEFF-late.pdf"))
# hm
# dev.off()


######################### DEG testing ########################################
# DESeq2
TPEX_DESeq2 <- TPEX %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = list(c("tissue", "BM", "TUM")),
    method = "DESeq2",
    omit_contrast_in_colnames = TRUE
  )

# edgeR
TPEX_de <- TPEX %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = c("tissueBM - tissueTUM"),
    method = "edgeR_quasi_likelihood",
    omit_contrast_in_colnames = TRUE
  )


deseq2 <- TPEX_DESeq2 %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% filter(padj < 0.05) %>% pull(feature)
length(deseq2)

edgeR <- TPEX_de %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% filter(FDR < 0.05) %>% pull(feature)
length(edgeR)

sum(deseq2 %in% edgeR)/length(deseq2)
sum(edgeR %in% deseq2)/length(edgeR)

TPEX_de %>% filter(.abundant) %>% 
  filter(FDR < 0.05) %>% 
  arrange(desc(logFC)) %>% 
  write_csv(file.path(saveDir, "TPEX_de.csv"))

TPEX_de %>% filter(.abundant) %>%
  filter(FDR < 0.02) %>%
  filter(logFC > 2) %>%
  arrange(desc(logFC)) %>%
  select(feature, logFC) %>%
  distinct() %>%
  write_csv(file.path(saveDir, "TPEX_BM_signature.csv"))

TPEX_DESeq2 %>% filter(.abundant) %>%
  filter(padj < 0.05) %>%
  arrange(desc(stat)) %>%
  write_csv(file.path(saveDir, "TPEX_DESeq2.csv"))
#################### Strip chart graph #######################################
topgenes_symbols <- c("Havcr2", "Entpd1", "Klrb1c", "Fcer1g", "Sell",
                      "Lef1", "Gzmb", "Cd7", "Id3", "Slamf6", "Cx3cr1", "Id2")
topgenes_symbols <- c("Havcr2", "Entpd1", "Sell", "Slamf6",
                      "Cx3cr1", "Tcf7", "Tbx21", "Tox", "Klrb1c", "S1pr5",
                      "S1pr1", "Myb", "Cd101", "Zfp683")

strip_chart <-
  counts_scaled %>%
  
  dplyr::filter(timepoint == "late") %>%
  
  mutate(tis_cel = paste0(cell.type, "_", tissue)) %>%
  
  # extract counts for top differentially expressed genes
  filter(feature %in% topgenes_symbols) %>%
  
  # make stripchart
  ggplot(aes(x = tis_cel, y = counts_scaled + 1, fill = tis_cel, label = "")) +
  geom_boxplot() +
  geom_jitter() +
  facet_wrap(~ feature) +
  scale_y_continuous(trans = "log2") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

strip_chart

grcm38[grep("Def", grcm38$symbol),] %>% select(ensgene, symbol) %>% as.data.frame()

# DESeq2
TEFF_DESeq2 <- TEFF %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = list(c("tissue", "BM", "TUM")),
    method = "DESeq2",
    omit_contrast_in_colnames = TRUE
  )

# edgeR
TEFF_de <- TEFF %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = c("tissueBM - tissueTUM"),
    method = "edgeR_quasi_likelihood",
    omit_contrast_in_colnames = TRUE
  )


deseq2 <- TEFF_DESeq2 %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% filter(padj < 0.05) %>% pull(feature)
length(deseq2)

edgeR <- TEFF_de %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% filter(FDR < 0.05) %>% pull(feature)
length(edgeR)

sum(deseq2 %in% edgeR)/length(deseq2)
sum(edgeR %in% deseq2)/length(edgeR)

TEFF_de %>% filter(.abundant) %>% 
  filter(FDR < 0.05) %>% 
  arrange(desc(logFC)) %>% 
  write_csv(file.path(saveDir, "TEFF_de.csv"))

TEFF_DESeq2 %>% filter(.abundant) %>%
  filter(padj < 0.05) %>%
  arrange(desc(stat)) %>%
  write_csv(file.path(saveDir, "TEFF_DESeq2.csv"))

unique(TEFF_de$feature) %in% unique(TPEX_de$feature)
unique(TPEX_de$feature)[unique(TEFF_de$feature) %in% unique(TPEX_de$feature)]

################# BM signature #################################################
tpex_sig <- TPEX_de %>% filter(.abundant) %>%
  filter(FDR < 0.02) %>%
  filter(logFC > 2) %>%
  arrange(desc(logFC)) %>%
  select(feature, logFC) %>%
  distinct()

teff_sig <- TEFF_de %>% filter(.abundant) %>%
  filter(FDR < 0.02) %>%
  filter(logFC > 2) %>%
  arrange(desc(logFC)) %>%
  select(feature, logFC) %>%
  distinct()

teff_sig %>%
  write_csv(file.path(saveDir, "TEFF_BM_signature.csv"))

BM_sig <- inner_join(tpex_sig, teff_sig, by = c("feature" = "feature"))
grep("Klrb1c", BM_sig$feature)
grep("Cd7", BM_sig$feature)
BM_sig %>% write_csv(file = file.path(saveDir, "BM_sig.csv"))
######################## volcano plots ########################################

topgenes <-
  TPEX_de %>%
  pivot_transcript() %>%
  dplyr::filter(FDR < 0.05) %>% 
  arrange(desc(logFC)) %>%
  head(100)
as.data.frame(topgenes %>% select(feature, logFC, FDR))

topgenes_symbols <- topgenes %>% pull(feature)
topgenes_symbols <- c("Lag3", "Pdcd1", "Rgs16",
                      "Sell", "Lef1", "S1pr1")

TPEX_de %>% filter(feature == c("Klrb1c", "Fcer1g")) %>% dplyr::select(feature, FDR, logFC)

volcano <- TPEX_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  # geom_text_repel() +
  geom_label_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano
ggsave(filename = file.path(plotDir, "Volcano_TPEX.pdf"))

topgenes <-
  TEFF_de %>%
  pivot_transcript() %>%
  dplyr::filter(FDR < 0.05) %>% 
  arrange(desc(logFC)) %>%
  head(100)
as.data.frame(topgenes %>% select(feature, logFC, FDR))

TEFF_de %>% filter(feature == c("Klrb1c", "Havcr2")) %>% dplyr::select(feature, FDR, logFC)

topgenes_symbols <- c("Lag3", "Pdcd1", "Havcr2", "Gzma", "S1pr5", "Lef1")


volcano <- TEFF_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  # geom_text_repel() +
  geom_label_repel(max.overlaps = getOption("ggrepel.max.overlaps", default = 20)) +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano
ggsave(filename = file.path(plotDir, "Volcano_TEFF.pdf"))


counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500, .dims = 3)

attr(counts_scal_PCA, "internals")$PCA

counts_scal_PCA %>%
  mutate(tis_cel = paste(tissue, cell.type, sep = "_")) %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, z = PC3, colour = cell.type, shape = tissue)) +
  theme_void() +
  axes_3D() +
  stat_3D() +
  geom_point(size = 4) +
  geom_text_repel(aes(label = timepoint), show.legend = FALSE)
  # stat_ellipse(type = "norm", level = 0.7)
# ggsave(file.path(plotDir, "PCA_top100_late_noLables.pdf"), device = "pdf")

counts_scal_PCA %>%
  mutate(tis_cel = paste(tissue, cell.type, sep = "_")) %>%
  plot_ly(x = .$PC1, y =.$ PC2, z = .$PC3, type = "scatter3d", mode = "markers", color = .$tis_cel)

counts_scaled <- counts_scaled %>%
  dplyr::filter(timepoint == "late") %>%
  dplyr::filter(cell.type == "TPEX")

# Reduce data dimensionality with arbitrary number of dimensions
tt_mds <- counts_scaled %>% reduce_dimensions(method = "MDS", .dims = 6, top = 500)

tt_mds %>%
  pivot_sample() %>%
  ggplot(aes(x = Dim1, y = Dim2, colour = sampleName, shape = cell.type)) +
  geom_point(size = 4) +
  # stat_ellipse(level = 0.7, type = "norm") +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  theme_bw()
# ggsave(file.path(plotDir, "MDS_top500.pdf"), device = "pdf")

hm <- counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature) %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    # row_names_gp = gpar(fontsize = 4),
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    column_km = 2,
    column_km_repeats = 100,
    row_km = 2,
    row_km_repeats = 500,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A",
                                       "#33A02C", "#FB9A99"), 
                              font = c(1,2,3))
  ) %>%
  add_tile(c(tissue))
hm
pdf(file = file.path(plotDir, "heatmap_top500_TPEX-late.pdf"))
hm
dev.off()

#################### Heatmap - specific genes #################################
row_labels <- c("Cx3cr1", "Fcer1g")

# original matrix

matrix <- counts_scaled %>%
  
  # filter lowly abundant
  filter(.abundant) %>%
  
  # extract most variable genes
  keep_variable( .abundance = counts_scaled, top = 500) %>%
  
  as_tibble() %>%
  
  mutate(genes = feature)

name_list <- matrix %>% 
  pull(feature) 
name_list <- unique(name_list)
row_labels %in% name_list

name_list <- name_list %>% 
  as_tibble() %>%
  rownames_to_column() %>%
  mutate(position = ifelse(value %in% row_labels, rowname, 0)) %>%
  filter(position > 0) %>%
  dplyr::select(value, position)

name_list


ha <- rowAnnotation(foo = anno_mark(at = name_list %>% pull(position) %>% as.integer(),
                                    labels = name_list %>% pull(value),
                                    which = "row", side = "right"))

ha




hm <-  matrix %>%
  
  # create heatmap
  heatmap(
    .column = sample,
    .row = genes,
    .value = counts_scaled,
    transform = log1p,
    palette_value = c("blue", "white", "red"),
    show_column_names = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    column_km = 3,
    column_km_repeats = 500,
    cluster_rows = TRUE,
    row_km = 5,
    row_km_repeats = 500,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A",
                                       "#33A02C", "#FB9A99"), 
                              font = c(1,2,3))
  ) %>%
  add_tile(c(tissue, timepoint, cell.type))

hm
pdf(file = file.path(plotDir, "heatmap_top500_ALL_NAMES.pdf"))
hm
dev.off()

# DESeq2
counts_de_DESeq2 <- counts_scaled %>%
  dplyr::filter(timepoint == "late") %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = list(c("tissue", "BM", "TUM")),
    method = "DESeq2",
    omit_contrast_in_colnames = TRUE
  )

# edgeR
counts_de <- counts_scaled %>%
  dplyr::filter(timepoint == "late") %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = c("tissueBM - tissueTUM"),
    method = "edgeR_quasi_likelihood",
    omit_contrast_in_colnames = TRUE
  )


deseq2 <- counts_de_DESeq2 %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% filter(padj < 0.05) %>% pull(feature)
length(deseq2)

edgeR <- counts_de %>% pivot_transcript(.transcript = feature) %>% 
  filter(.abundant) %>% filter(FDR < 0.05) %>% pull(feature)
length(edgeR)

sum(deseq2 %in% edgeR)/length(deseq2)
sum(edgeR %in% deseq2)/length(edgeR)

counts_de %>% filter(.abundant) %>% 
  filter(FDR < 0.05) %>% 
  arrange(desc(logFC)) %>%
  dplyr::select(feature, logFC, FDR) %>%
  distinct() %>%
  write_csv(file.path(saveDir, "BMvsTUM_de.csv"))

topgenes <-
  counts_de %>%
  pivot_transcript() %>%
  arrange(FDR) %>%
  head(50)

topgenes_symbols <- topgenes %>% pull(feature)

counts_de %>% filter(feature == "") %>% dplyr::select(feature, FDR, logFC)

volcano <- counts_de %>%
  pivot_transcript() %>%
  
  # Subset data
  filter(.abundant) %>%
  mutate(significant = FDR < 0.01 & abs(logFC) >= 2) %>%
  mutate(feature = ifelse(feature %in% topgenes_symbols, as.character(feature), "")) %>%
  
  # Plot
  ggplot(aes(x = logFC, y = PValue, label = feature)) +
  geom_point(aes(color = significant, size = significant, alpha = significant)) +
  geom_text_repel() +
  
  # Custom scales
  scale_y_continuous(trans = "log10_reverse") +
  scale_color_manual(values = c("black", "#e11f28")) +
  scale_size_discrete(range = c(0, 2)) +
  theme_bw()
volcano


