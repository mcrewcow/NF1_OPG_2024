library(Seurat)
library(SeuratDisk)

ZH11 <- readRDS('/media/baranov_lab/New Volume/NF_processed2024/ZH11.rds')
ZH12 <- readRDS('/media/baranov_lab/New Volume/NF_processed2024/ZH12.rds')
ZH13 <- readRDS('/media/baranov_lab/New Volume/NF_processed2024/ZH13.rds')
ZH14 <- readRDS('/media/baranov_lab/New Volume/NF_processed2024/ZH14.rds')
ZH15 <- readRDS('/media/baranov_lab/New Volume/NF_processed2024/ZH15.rds')
ZH16 <- readRDS('/media/baranov_lab/New Volume/NF_processed2024/ZH16.rds')

data1 <- LoadH5Seurat('/media/baranov_lab/Seagate/Microglia/GSE123758.h5Seurat')
data2 <- LoadH5Seurat('/media/baranov_lab/Seagate/Microglia/GSE199317_eyecup_CtC57RPECD45R1.h5Seurat')
head(data2)
data2
data3 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57AllOtherR1.rds')
data4 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57AllOtherR2.rds')
data5 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57CD45CD90R1.rds')
data6 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57CD45CD140aR1.rds')
data7 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57CD45CD140aR2.rds')
data8 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57GlastP1.rds')
data9 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE199317_ONC-retina_CtC57GlastR1.rds')
data11 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE243413 (2023) CD90.1.rds')
data12 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE243413 (2023) Chen_CD73.rds')
data13 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSE243413 (2023) WT_CD73.rds')
data10 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM3580725.rds')
data14 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM3580727.rds')
data15 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM3854512.rds')
data16 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM3854514.rds')
data17 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM3854516.rds')
data18 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM3854518.rds')
data19 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4078527.rds')
data20 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4078528.rds')
data21 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4078529.rds')
data22 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4078530.rds')
data23 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4089149.rds')
data24 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4649091.rds')
data25 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4649092.rds')
data26 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4649093.rds')
data27 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4649094.rds')

data28 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4649095.rds')
data29 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4649096.rds')
data30 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4995561.rds')
data31 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM4995562.rds')
data32 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5191739.rds')
data33 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5560840.rds')
data34 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5359947.rds')
data35 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5614917.rds')
data36 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5742457.rds')
data37 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5855990.rds')
data38 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM5855991.rds')
data39 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM6199212.rds')
data40 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM6205478.rds')
data41 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM6807948.rds')
data42 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM6807950.rds')
data43 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM6807951.rds')
data44 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7336928.rds')
data45 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7529034.rds')
data46 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7663480.rds')
data47 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7663481.rds')
data48 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7663482.rds')
data49 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7779359.rds')
data50 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7779360.rds')
data51 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM7986396.rds')
data52 <- readRDS('/media/baranov_lab/Seagate/Microglia/GSM8200462.rds')
gc()


integration_list <- list(data1,
                         data2,
                         data5,
                         data6,
                         data7,
                         data8,
                         data9,
                         data10,
                         data11,
                         data12,
                         data13,
                         data14,
                         data15,
                         data16,
                         data17,
                         data18,
                         data19,
                         data20,
                         data22,
                         data23,
                         data24,
                         data25,
                         data26,
                         data27,
                         data28,
                         data29,
                         data30,
                         data31,
                         data32,
                         data33,
                         data34,
                         data35,
                         data39,
                         data40,
                         data41,
                         data42, data43, data44, data45,
                         data46, data47,
                         data48,
                         data49,data50, data51, data52, ZH11, ZH12, ZH13, ZH14, ZH15, ZH16)

features <- SelectIntegrationFeatures(object.list = integration_list, nfeatures = 2000)
data.anchors <- FindIntegrationAnchors(object.list = integration_list, anchor.features = features)
data.anchors
MOUSE_INTretinaNF1<- IntegrateData(anchorset = data.anchors)

ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"), verbose = T) #, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score")
  data.integrated <- RunPCA(data.integrated, npcs = 30, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:30)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:30)
}
gc()


SaveH5Seurat(MOUSE_INTretinaNF1, '/home/baranov_lab/10X/MOUSE_INTretinaNF1.h5Seurat')



MOUSE_INTretinaNF1 <- ProcessInt(MOUSE_INTretinaNF1)
SaveH5Seurat(MOUSE_INTretinaNF1, '/home/baranov_lab/10X/MOUSE_INTretinaNF1_30.h5Seurat')
gc()
ProcessInt <- function(data.integrated){
  data.integrated <- ScaleData(data.integrated, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score"), verbose = T) #, vars.to.regress = c('percent.mt',"percent.rb","S.Score","G2M.Score")
  data.integrated <- RunPCA(data.integrated, npcs = 150, verbose = T)
  data.integrated <- FindNeighbors(data.integrated, dims = 1:150)
  data.integrated <- FindClusters(data.integrated, resolution = 1)
  data.integrated <- RunUMAP(data.integrated, reduction = "pca", dims = 1:150)
}
MOUSE_INTretinaNF1_150 <- ProcessInt(MOUSE_INTretinaNF1)
SaveH5Seurat(MOUSE_INTretinaNF1_150, '/home/baranov_lab/10X/MOUSE_INTretinaNF1_150.h5Seurat')

head(MOUSE_INTretinaNF1_150)

NF1ret150 <- LoadH5Seurat('/media/baranov_lab/New Volume/NF_processed2024/NF1ret_150npcs_anno.h5Seurat')
NF1_ON150 <- LoadH5Seurat('/media/baranov_lab/New Volume/NF_processed2024/NF1ON150_annotated.h5Seurat', overwrite =  TRUE)

MOUSE_INTconservnovars30 <- LoadH5Seurat('/home/baranov_lab/10X/mouse_atlas/conserv_process30var_anno_cxg.h5Seurat')

MOUSE_INTretinaNF1_150$EK_PB_anno2024 <- 'Adult'
NF1ret150$EK_PB_anno2024 <- as.character(NF1ret150$EK_PB_anno2024)
Idents(MOUSE_INTretinaNF1_150) <- "cluster"


table(mt$EK_PB_anno2024_2)
table(NF1ret150$orig.ident)



NF1ret150 <- SetIdent(NF1ret150, value = 'orig.ident')
NF1ret150 <- RenameIdents(NF1ret150,
                                       'ZH11' = 'NF1 fl/fl',
                                       'ZH12' = 'NF1 fl/fl',
                                       'ZH13' = 'NF1 wt/null',
                                       'ZH14' = 'NF1 fl/null OPG',
                                       'ZH15' = 'NF1 fl/null OPG',
                                       'ZH16' = 'NF1 fl/null OPG')
NF1ret150$condition <- NF1ret150@active.ident
NF1ret150$EK_PB_anno2024_2 <- NF1ret150$EK_PB_anno2024
MOUSE_INTconservnovars30$condition <- 'Atlas Control'


mt <- merge(x = NF1ret150, y = MOUSE_INTconservnovars30)
MOUSE_INTretinaNF1_150$EK_PB_anno2024 <- mt$EK_PB_anno2024_2
MOUSE_INTretinaNF1_150$EK_PB_anno2024 <- AddMetaData(MOUSE_INTretinaNF1_150, metadata = mt$condition, col.name = 'condition')
mt

metadata_NF1ret150 <- NF1ret150@meta.data %>%
  select(EK_PB_anno2024_2, condition)

# Extract metadata columns from MOUSE_INTconservnovars30
metadata_MOUSE_INTconservnovars30 <- MOUSE_INTconservnovars30@meta.data %>%
  select(EK_PB_anno2024_2, condition)

# Combine the extracted metadata
combined_metadata <- bind_rows(metadata_NF1ret150, metadata_MOUSE_INTconservnovars30)

# Ensure rownames in the combined metadata match the integrated object
combined_metadata <- combined_metadata[rownames(MOUSE_INTretinaNF1_150@meta.data), ]

# Assign the combined metadata to the integrated Seurat object
MOUSE_INTretinaNF1_150@meta.data$EK_PB_anno2024_2 <- combined_metadata$EK_PB_anno2024_2
MOUSE_INTretinaNF1_150@meta.data$condition <- combined_metadata$condition



table(NF1ret150$EK_PB_anno2024_2)
RGC1 <- subset(NF1ret150, subset = EK_PB_anno2024_2 == 'RGC')
table(RGC1$condition)
table(MOUSE_INTconservnovars30$EK_PB_anno2024_2)
RGC2 <- subset(MOUSE_INTconservnovars30, subset = EK_PB_anno2024_2 == 'RGC')
RGCS <- merge(x = RGC1, y = RGC2)
SaveH5Seurat(RGCS, '/home/baranov_lab/10X/MOUSE_INTretinaNF1_150_RGCsonly.h5Seurat')


