rm(list = ls())
library(Seurat)
library(ggplot2)
library(DropSeq.util)



##CB
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Cerebellum_ALT.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- rownames(dge)
cnt <- dge


##FC
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##PC
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##ENT
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60EntoPeduncular.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##GP
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60GlobusPallidus.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##Hippo
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Hippocampus.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##STR
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Striatum.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##SN
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60SubstantiaNigra.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])


##TH
dge.path <- "/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Thalamus.raw.dge.txt.gz"
dge <- loadSparseDge(dge.path) 
genes <- intersect(genes, rownames(dge))
cnt <- cbind(cnt[genes, ], dge[genes, ])



sce <- cnt
save(sce, file="/import/home/share/Portal-reproduce/data_MB3/Drop-seq_all_withNAanno.RData")



load("/import/home/share/Portal-reproduce/data_MB3/Drop-seq_meta_all.RData")
sce <- sce[, as.vector(meta_all$cellid)]
print(dim(sce))
save(sce, file="/import/home/share/Portal-reproduce/data_MB3/Drop-seq_all.RData")




