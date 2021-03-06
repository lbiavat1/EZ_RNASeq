rm(list = ls())


# tidyverse-friendly
library(tidyverse)
library(tidybulk)
library(tidyHeatmap)
library(ComplexHeatmap)
library(ggrepel)
library(plotly)
library(GGally)

# other useful libraries
library(DESeq2)
library(limma)
library(edgeR)
library(fgsea)
library(hciR)
library(IHW)


library(biomaRt)
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

names(myTibble)[1] <- "feature"
myTibble
grcm38
grcm38_tx2gene

pivot_longer(myTibble, cols = c(2:18),
             names_to = "sample", values_to = "counts")
counts_tt <- pivot_longer(myTibble, cols = c(2:18),
                          names_to = "sample", values_to = "counts") %>%
  dplyr::inner_join(grcm38, by = c("feature" = "ensgene")) %>%
  dplyr::select(symbol, sample, counts)
counts_tt
names(counts_tt)[1] <- "feature"
  
"Tcf7" %in% counts_tt$feature
"Gzmb" %in% counts_tt$feature
"Elane" %in% counts_tt$feature
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

counts_tt <- counts_tt %>% left_join(myInfo)
counts_tt
counts <- counts_tt %>%
  tidybulk(.sample = sample, .transcript = feature, .abundance = counts)
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

counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 500)
counts_scal_PCA <-
  counts_scaled %>%
  reduce_dimensions(method = "PCA", top = 100)
attr(counts_scal_PCA, "internals")$PCA

counts_scal_PCA %>%
  mutate(tis_cel = paste(tissue, cell.type, sep = "_")) %>%
  pivot_sample() %>%
  ggplot(aes(x = PC1, y = PC2, colour = cell.type, shape = tissue)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = ""), show.legend = FALSE) +
  # stat_ellipse(type = "norm", level = 0.7) +
  theme_bw()
# ggsave(file.path(plotDir, "PCA_top100.pdf"), device = "pdf")

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
    column_km = 3,
    column_km_repeats = 100,
    row_km = 5,
    row_km_repeats = 500,
    row_title = "%s",
    row_title_gp = grid::gpar(fill = c("#A6CEE3", "#1F78B4", "#B2DF8A",
                                       "#33A02C", "#FB9A99"), 
                              font = c(1,2,3))
  ) %>%
  add_tile(c(tissue, timepoint, cell.type))
hm
pdf(file = file.path(plotDir, "heatmap_top500_ALL_RowClusters_0315.pdf"))
hm
dev.off()

#################### Heatmap - specific genes #################################
row_labels <- c("Cx3cr1", "S1pr5")

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


pdf(file = file.path(plotDir, "heatmap_top500_ALL_NAMES.pdf"))
hm
dev.off()

# DESeq2
counts_de_DESeq2 <- counts_scaled %>%
  test_differential_abundance(
    .formula = ~ 0 + tissue + mouse,
    .contrasts = list(c("tissue", "BM", "TUM")),
    method = "DESeq2",
    omit_contrast_in_colnames = TRUE
  )

# edgeR
counts_de <- counts_scaled %>%
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
