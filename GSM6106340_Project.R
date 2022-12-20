###Final Project info-B-573



#Installing required packages
install.packages("Seurat")
install.packages("SCtransform")


#Loading the packages
library(Seurat)
library(sctransform)



#Reading the file and Creating Count Matrix
countmatrix=Read10X("C:/Users/Lenovo/Desktop/project/GSM6106340/")




dim(countmatrix)





#Creating Seurat Object
srobj=CreateSeuratObject(countmatrix, project = "GSM6106340", min.cells = 3, min.features = 200)

srobj










srobj[["MTpercent"]]= PercentageFeatureSet(srobj,pattern = "^MT-")
srobj= SCTransform(srobj,vars.to.regress = "MTpercent", verbose = F)
srobj


srobj=RunPCA(srobj, features = VariableFeatures(srobj))



#UMAP
srobj=RunUMAP(srobj, dims = 1:30)
DimPlot(srobj, reduction = "tsne")


#tsne
srobj= RunTSNE(srobj,dims = 1:30)
DimPlot(srobj, reduction = "tsne")






srobj= FindNeighbors(srobj, dims = 1:30)
srobj= FindClusters(srobj,resolution = 0.5)
DimPlot(srobj, reduction = "umap", label =T)

#tsne
srobj= RunTSNE(srobj,dims = 1:30)
DimPlot(srobj, reduction = "tsne")

#markers
markers= FindAllMarkers(srobj, min.pct = 0.25, logfc.threshold = 1, only.pos = T)
dim(markers)

head(srobj@meta.data)
View(markers) 

cl=markers[which(markers$cluster==0),]

VlnPlot(srobj, cl$gene)

RidgePlot(srobj, cl$gene)

FeaturePlot(srobj, cl$gene)

DoHeatmap(srobj, features=cl$gene)


mean(srobj)

#celltypes
cellClusterIDs=0:21
cellTypes= c("Lake et al.Science.Ex2", "Neuronal cell-Cancer", "Neuron", "Michroglial Cell", "Oligodendrocyte", "Oligodendrocyte", "Astrocyte", "Retinal ganglion cells", "Neuron", "Cajal-Retzius cells", "Neuroendocrine", "Pericyte", "Inhibitory Neurons", "Lake et al.Science.Ex6", "Cajal-Retzius cells","Lake et al.Science.Ex8", "Epidermal-like cell(Cancer)", "Stromal cell", "Microglia", "Chandlier cell", "Lake et al.Science.Ex8", "Lake et al.Science.Ex8")
names(x=cellTypes)=levels(x=srobj)
srobj2= RenameIdents(object = srobj, cellTypes)
DimPlot(srobj2, reduction = "umap", label = T)       






View(markers)
summary(markers$avg_log2FC)
min(markers$avg_log2FC[which(markers$avg_log2FC >0)])
dim(markers[which(markers$avg_log2FC>=1),])
View(markers[which(markers$avg_log2FC>=1),])


markers= markers[which(markers$avg_log2FC>=1),]
clusters=markers[which(markers$cluster==1),]
View(clusters)

#DoHeatmap
DoHeatmap(srobject, features=clusters$gene)


#Featureplot
FeaturePlot(srobject, features = clusters$gene)


#Violine plot
VlnPlot(srobject, features = clusters$gene)


#Ridgeplot
RidgePlot(srobject, features = clusters$gene)
