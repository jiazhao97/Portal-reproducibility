rm(list = ls())
library(Seurat)
library(DropSeq.util)

##CB
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Cerebellum_ALT.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:2, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 3, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Choroid_plexus"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1:2, sep = ""))] <- "Endothelial_tip"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 3:5, sep = ""))] <- "Mural"

meta_CB_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##FC
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Cortex_noRep5_FRONTALonly.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 2, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 3, sep = ""))] <- "Macrophage"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 4, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("12-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("13-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("14-", 1:27, sep = ""))] <- "Endothelial_tip"

meta_FC_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))




##PC
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Cortex_noRep5_POSTERIORonly.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", c(1,3,4), sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 2, sep = ""))] <- "Macrophage"
cluster_assign_vec[which(cluster_assign_vec %in% paste("12-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("13-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("14-", 1:27, sep = ""))] <- "Endothelial_tip"

meta_PC_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##ENT
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60EntoPeduncular.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec
sort(unique(subcluster_assign_vec))

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", c(1,3,4,5,6), sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 2, sep = ""))] <- "Endothelial_tip"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Oligodendrocyte"

meta_ENT_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##GP
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60GlobusPallidus.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec
sort(unique(subcluster_assign_vec))

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", c(1,2,3,4,5,7), sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 6, sep = ""))] <- "Mitotic"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Ependymal"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Endothelial_tip"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 2, sep = ""))] <- "Macrophage"

meta_GP_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##Hippo
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Hippocampus.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec
sort(unique(subcluster_assign_vec))

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 2, sep = ""))] <- "Macrophage"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1:27, sep = ""))] <- "Ependymal"
cluster_assign_vec[which(cluster_assign_vec %in% paste("12-", 1:27, sep = ""))] <- "Choroid_plexus"
cluster_assign_vec[which(cluster_assign_vec %in% paste("13-", 1:27, sep = ""))] <- "Neurogenesis"
cluster_assign_vec[which(cluster_assign_vec %in% paste("14-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("15-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("16-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("17-", 1:27, sep = ""))] <- "Endothelial_tip"

meta_Hippo_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##STR
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Striatum.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec
sort(unique(subcluster_assign_vec))

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Ependymal"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neurogenesis"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1, sep = ""))] <- "Macrophage"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 2, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Endothelial_tip"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("12-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("13-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("14-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("15-", 1:27, sep = ""))] <- "Neuron"

meta_STR_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##SN
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60SubstantiaNigra.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec
sort(unique(subcluster_assign_vec))

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Polydendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Ependymal"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1, sep = ""))] <- "Macrophage"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 2, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("12-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("13-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("14-", 1:27, sep = ""))] <- "Endothelial_tip"

meta_SN_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



##TH
subcluster_assign <- readRDS("/import/home/share/Portal-reproduce/Raw_data_Saunders_full/F_GRCm38.81.P60Thalamus.subcluster.assign.RDS")
cluster_assign_vec <- as.vector(as.factor(subcluster_assign))
names(cluster_assign_vec) <- names(subcluster_assign)
subcluster_assign_vec <- cluster_assign_vec
sort(unique(subcluster_assign_vec))

cluster_assign_vec[which(cluster_assign_vec %in% paste("1-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("2-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("3-", 1:27, sep = ""))] <- "Neuron"
cluster_assign_vec[which(cluster_assign_vec %in% paste("4-", 1:27, sep = ""))] <- "Endothelial_stalk"
cluster_assign_vec[which(cluster_assign_vec %in% paste("5-", 1:27, sep = ""))] <- "Mural"
cluster_assign_vec[which(cluster_assign_vec %in% paste("6-", 1:27, sep = ""))] <- "Endothelial_tip"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 1, sep = ""))] <- "Macrophage"
cluster_assign_vec[which(cluster_assign_vec %in% paste("7-", 2, sep = ""))] <- "Microglia"
cluster_assign_vec[which(cluster_assign_vec %in% paste("8-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("9-", 1:27, sep = ""))] <- "Oligodendrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("10-", 1:27, sep = ""))] <- "Ependymal"
cluster_assign_vec[which(cluster_assign_vec %in% paste("11-", 1:27, sep = ""))] <- "Astrocyte"
cluster_assign_vec[which(cluster_assign_vec %in% paste("12-", 1:27, sep = ""))] <- "Polydendrocyte"

meta_TH_new <- data.frame(cellid = names(cluster_assign_vec), subcluster = as.vector(subcluster_assign), class = as.vector(cluster_assign_vec))



meta_all_new <- rbind(meta_CB_new, meta_FC_new)
meta_all_new <- rbind(meta_all_new, meta_PC_new)
meta_all_new <- rbind(meta_all_new, meta_ENT_new)
meta_all_new <- rbind(meta_all_new, meta_GP_new)
meta_all_new <- rbind(meta_all_new, meta_Hippo_new)
meta_all_new <- rbind(meta_all_new, meta_STR_new)
meta_all_new <- rbind(meta_all_new, meta_SN_new)
meta_all_new <- rbind(meta_all_new, meta_TH_new)

meta_all_new$class <- as.vector(meta_all_new$class)
meta_all_new$class[(as.vector(meta_all_new$class) == "Neuron")] <- "NEURON"
meta_all_new$class[as.vector(meta_all_new$class) == "Oligodendrocyte"] <- "OLIGODENDROCYTE"
meta_all_new$class[as.vector(meta_all_new$class) == "Polydendrocyte"] <- "POLYDENDROCYTE"
meta_all_new$class[as.vector(meta_all_new$class) == "Astrocyte"] <- "ASTROCYTE"
meta_all_new$class[as.vector(meta_all_new$class) == "Endothelial_stalk"] <- "ENDOTHELIAL_STALK"
meta_all_new$class[as.vector(meta_all_new$class) == "Mural"] <- "MURAL"
meta_all_new$class[as.vector(meta_all_new$class) == "Endothelial_tip"] <- "ENDOTHELIAL_TIP"
meta_all_new$class[as.vector(meta_all_new$class) == "Microglia"] <- "MICROGLIA"
meta_all_new$class[as.vector(meta_all_new$class) == "Choroid_plexus"] <- "CHOROID_PLEXUS"
meta_all_new$class[as.vector(meta_all_new$class) == "Macrophage"] <- "MACROPHAGE"
meta_all_new$class[as.vector(meta_all_new$class) == "Mitotic"] <- "MITOTIC"
meta_all_new$class[as.vector(meta_all_new$class) == "Ependymal"] <- "EPENDYMAL"
meta_all_new$class[as.vector(meta_all_new$class) == "Neurogenesis"] <- "NEUROGENESIS"



meta_all <- meta_all_new
save(meta_all, file="/import/home/share/Portal-reproduce/data_MB3/Drop-seq_meta_all.RData")




