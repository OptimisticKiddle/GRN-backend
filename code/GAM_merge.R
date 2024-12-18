library(hdf5r)
library(Seurat) # scRNA-seq 主要存储数据的对象类型
library(Matrix)
library(EnsDb.Mmusculus.v79) # 针对此数据集 要用小鼠基因组注释文件的话
library(EnsDb.Hsapiens.v75) # 针对此数据集 要用人类基因组注释文件的话
library(patchwork)
library(Signac)# 处理scATAC-seq数据的包
library(ggplot2)# 绘图
library(spatstat.geom)
library(biovizBase)
library(AnnotationHub) # 主要用于基因注释
library(irlba)
library(dplyr)
library(future)
args<-commandArgs(T)
args[1]
args[2]
# 设置下输出图形大小
options(repr.plot.width = 15, repr.plot.height =9)
plan("multicore", workers = 32)
options(future.globals.maxSize = 50000 * 1024^2) # for 50 Gb RAM
basePath = sprintf("../static/GSE%s/GSM%s",args[1],args[2])
cluster_k = -1 # -1表示后续使用默认参数cluster k分辨率
# 加载上一步生成的 baseGRN 的行名
load(sprintf("%s/GRN_row_names.rda",basePath)) 
# 加载 WT和KO的GAM
wt_GAM <- readRDS(sprintf("%s/GAM.rds",basePath))
ko_GAM <- readRDS(sprintf("%s/KO/GAM_KO.rds",basePath))
# 依row_names进行过滤
wt_GAM = wt_GAM[rownames(wt_GAM) %in% row_names, ]
ko_GAM = ko_GAM[rownames(ko_GAM) %in% row_names, ]
# 调试 查看维度
dim(wt_GAM)
dim(ko_GAM)
# 分别创建wt和ko seurat对象  数据别名为RNA
wt_seurat <- CreateSeuratObject(counts = wt_GAM,assay = "RNA")
ko_seurat <- CreateSeuratObject(counts = ko_GAM,assay = "RNA")
# 添加信息以识别原始数据集
wt_seurat$dataset <- 'WT'
ko_seurat$dataset <- 'KO'
# 创建合并
comb_seurat <- merge(wt_seurat, ko_seurat, add.cell.ids = c("WT", "KO"))
comb_seurat$RNA$data.1 = comb_seurat$RNA$counts.1
comb_seurat$RNA$data.2 = comb_seurat$RNA$counts.2
# 调试 查看信息
comb_seurat
# Visualize QC metrics as a violin plot
# 将QC指标可视化为小提琴图
Vp = VlnPlot(comb_seurat, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
ggsave(sprintf("%s/refine/violin_plot.png",basePath), Vp,width = 15,height = 9)
comb_seurat <- FindVariableFeatures(comb_seurat)
# Identify the 10 most highly variable genes
# 找出10个最易变的基因
top10 <- head(VariableFeatures(comb_seurat), 10)

# plot variable features with and without labels
# 绘制有和没有标签的变量特征
plot1 <- VariableFeaturePlot(comb_seurat)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
pp = plot1 + plot2
ggsave(sprintf('%s/refine/VariableFeature.png',basePath),pp,width = 15,height = 9)
comb_seurat <- ScaleData(comb_seurat)
comb_seurat <- RunPCA(comb_seurat, features = VariableFeatures(object = comb_seurat),verbose = FALSE)

# 集成
comb_seurat <- IntegrateLayers(
    object = comb_seurat, 
    method = CCAIntegration, 
    orig.reduction = "pca", 
    new.reduction = "integrated.cca",
    verbose = FALSE
)
# re-join layers after integration
# 集成后重新连接层
comb_seurat[["RNA"]] <- JoinLayers(comb_seurat[["RNA"]])
comb_seurat <- FindNeighbors(comb_seurat, reduction = "integrated.cca", dims = 1:30)
# 二分resolution来达到想要的cluster k效果。 若没有特定指定cluster k，则进行默认聚类
left = 0
right = cluster_k
now_k = -1
mx_cnt = 0 # 防止死循环
if (cluster_k == -1){
    comb_seurat <- FindClusters(comb_seurat, verbose = FALSE) #  verbose = FALSE
}else {
    # 二分resolution来达到想要的cluster k效果
    while (now_k != cluster_k && mx_cnt < 100) {
        mid = (left + right)/2
        comb_seurat <- FindClusters(comb_seurat, resolution = mid, verbose = FALSE) #  verbose = FALSE
        now_k = max(as.numeric(comb_seurat@meta.data$seurat_clusters), na.rm = TRUE)
        if (now_k > cluster_k) {
            right = mid
        }else {
            left = mid
        }
        mx_cnt = mx_cnt + 1
    }
}
# 使用dplyr包删除列名以"RNA_snn_res."为前缀的所有列
comb_seurat@meta.data <- comb_seurat@meta.data %>% select(-matches("^RNA_snn_res\\..*"))
# 调试 查看信息
head(Idents(comb_seurat), 5) # 此时已经有cluster map关系了，后续的UMAP只是为了可视化 ！！！ 
comb_seurat <- RunUMAP(comb_seurat, dims = 1:30, reduction = "integrated.cca",verbose = F)
# 为了并排可视化这两个条件，我们可以使用参数来显示按聚类着色的每个条件
dimp = DimPlot(comb_seurat, reduction = "umap", split.by = "dataset")
ggsave(sprintf('%s/refine/umap.png',basePath),dimp,width = 15,height = 9)
# 调试 查看信息
comb_seurat
# 抽取信息 查看各cell对应的cluster
wt_df = comb_seurat@meta.data[comb_seurat@meta.data$dataset == "WT", "seurat_clusters", drop = FALSE]
ko_df = comb_seurat@meta.data[comb_seurat@meta.data$dataset == "KO", "seurat_clusters", drop = FALSE]
# 抽取cell barcode信息
wt_ncol = nrow(wt_df)
ko_ncol = nrow(ko_df)
# 获取seurat_clusters的唯一值并排序
cluster_levels <- sort(unique(comb_seurat@meta.data$seurat_clusters))
# wt
# 创建一个空的结果矩阵（ dense 数值型矩阵，初始值为0）
wt_avg_expr_by_cluster <- matrix(0, nrow = nrow(wt_GAM), ncol = length(cluster_levels),
                              dimnames = list(rownames(wt_GAM), cluster_levels))

# 计算每个cluster中各基因的平均表达值
for (i in seq_along(cluster_levels)) {
  current_cluster_cells <- which(wt_df$seurat_clusters == cluster_levels[i])
  
  # 提取当前cluster对应的子矩阵
  sub_mat <- wt_seurat$RNA$count[, current_cluster_cells, drop = FALSE]  # 保持矩阵形式，避免drop掉单列
  
  # 将子矩阵转化为普通矩阵并计算平均表达值
  gene_avg_expr <- apply(as.matrix(sub_mat), 1, mean, na.rm = TRUE)
  
  # 将平均表达值填充到结果矩阵中
  wt_avg_expr_by_cluster[, i] <- gene_avg_expr
}

# 现在avg_expr_by_cluster就是一个n x k维度的矩阵，其中n为基因数，k为cluster数
# 行名是基因名，列名按照cluster编号顺序排列
# ko
# 创建一个空的结果矩阵（ dense 数值型矩阵，初始值为0）
ko_avg_expr_by_cluster <- matrix(0, nrow = nrow(ko_GAM), ncol = length(cluster_levels),
                              dimnames = list(rownames(ko_GAM), cluster_levels))

# 计算每个cluster中各基因的平均表达值
for (i in seq_along(cluster_levels)) {
  current_cluster_cells <- which(ko_df$seurat_clusters == cluster_levels[i])
  
  # 提取当前cluster对应的子矩阵
  sub_mat <- ko_seurat$RNA$count[, current_cluster_cells, drop = FALSE]  # 保持矩阵形式，避免drop掉单列
  
  # 将子矩阵转化为普通矩阵并计算平均表达值
  gene_avg_expr <- apply(as.matrix(sub_mat), 1, mean, na.rm = TRUE)
  
  # 将平均表达值填充到结果矩阵中
  ko_avg_expr_by_cluster[, i] <- gene_avg_expr
}

# 现在avg_expr_by_cluster就是一个n x k维度的矩阵，其中n为基因数，k为cluster数
# 行名是基因名，列名按照cluster编号顺序排列

# 行名按字典序排序
# 对两个矩阵进行按行名字典序排序
rownames_WT <- rownames(wt_avg_expr_by_cluster)
rownames_WT_ordered <- rownames_WT[order(rownames_WT)]
wt_avg_expr_by_cluster <- wt_avg_expr_by_cluster[rownames_WT_ordered, ]

rownames_KO <- rownames(ko_avg_expr_by_cluster)
rownames_KO_ordered <- rownames_KO[order(rownames_KO)]
ko_avg_expr_by_cluster <- ko_avg_expr_by_cluster[rownames_KO_ordered, ]

# 保存gene by cluster格式的GAM 供后续refine使用
write.csv(wt_avg_expr_by_cluster, sprintf("%s/refine/wt_avg_expr_by_cluster.csv",basePath), row.names = TRUE)
write.csv(ko_avg_expr_by_cluster, sprintf("%s/refine/ko_avg_expr_by_cluster.csv",basePath), row.names = TRUE)