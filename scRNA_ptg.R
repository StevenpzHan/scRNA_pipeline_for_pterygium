library(Seurat)
library(monocle)
library(ggplot2)
library(ggsci)
library(ggpubr)
library(reshape2)
library(EnhancedVolcano)
library(GSVA)
library(clusterProfiler)
library(limma)
library(MAST)
library(dorothea)
library(CellChat)
library(RColorBrewer)
library(enrichplot)
library(harmony)
library(dplyr)
library(plyr)
library(COSG)
library(pheatmap)
library(msigdbr)
library(GEOquery)
library(AnnoProbe)
#Data input
setwd("E:/scrna_ptg/")
ptg1=Read10X("case1/filtered_feature_bc_matrix/")
T1=CreateSeuratObject(ptg1,project='Pterygium1',min.features = 200,min.cells = 3)
ptg2=Read10X("case2/filtered_feature_bc_matrix/")
T2=CreateSeuratObject(ptg2,project='Pterygium1',min.features = 200,min.cells = 3)
ptg3=Read10X("case3/filtered_feature_bc_matrix/")
T3=CreateSeuratObject(ptg3,project='Pterygium1',min.features = 200,min.cells = 3)
con1=Read10X("control1//filtered_feature_bc_matrix/")
C1=CreateSeuratObject(con1,project='Normal1',min.features = 200,min.cells = 3)
con2=Read10X("control2/filtered_feature_bc_matrix/")
C2=CreateSeuratObject(con2,project='Normal2',min.features = 200,min.cells = 3)
con3=Read10X("control3/filtered_feature_bc_matrix/")
C3=CreateSeuratObject(con3,project='Normal3',min.features = 200,min.cells = 3)
#Doublet removing
T1
T1 = NormalizeData(T1)
T1 = ScaleData(T1)
T1 = FindVariableFeatures(T1, selection.method = "vst", nfeatures = 10000)
T1 = RunPCA(T1, features = T1@assays$RNA@var.features)
T1 =  FindNeighbors(T1, reduction = "pca", dims = 1:50)
T1 = FindClusters(T1, resolution = 0.6)

sweep.data = paramSweep_v3(T1, PCs = 1:50)
sweep.stats = summarizeSweep(sweep.data, GT = FALSE)
bcmvn= find.pK(sweep.stats)

homotypic.prop=modelHomotypic(T1@meta.data$RNA_snn_res.0.6)         
nExp_poi=round(0.075*length(T1$orig.ident))  
nExp_poi.adj=round(nExp_poi*(1-homotypic.prop))
T1=doubletFinder_v3(T1, PCs = 1:50, pN = 0.25, pK = as.numeric(as.character(bcmvn$pK[which.max(bcmvn$BCmetric)])), 
                    nExp = nExp_poi.adj, reuse.pANN = FALSE)
doubletsID = T1@meta.data[,8]
T1 = SubsetData(T1, cells = colnames(T1)[which(doubletsID == "Singlet")])

#Data integration & First-round analysis
T1=RenameCells(T1,add.cell.id = "PTG1")
T2=RenameCells(T2,add.cell.id = "PTG2")
T3=RenameCells(T3,add.cell.id = "PTG3")
C1=RenameCells(C1,add.cell.id = "NC1")
C2=RenameCells(C2,add.cell.id = "NC2")
C3=RenameCells(C3,add.cell.id = "NC3")
allptg=merge(T1,list(T2,T3,C1,C2,C3))
allptg[["percent.mt"]] = PercentageFeatureSet(object = allptg, pattern = "^MT-")
violin1 <- VlnPlot(scRNA,
features = c("nFeature_RNA", "nCount_RNA", "percent.mt","percent.HB"),
cols =rainbow(col.num))
mytotal= subset(allptg, subset =  nFeature_RNA > 500 & nFeature_RNA < 6000 & percent.mt < 10 & nCount_RNA >500)
mytotal=NormalizeData(namedtotal)
mytotal=FindVariableFeatures(mytotal,nfeatures = 3000)
mytotal=ScaleData(mytotal,vars.to.regress = c('percent.mt'))
mytotal=RunPCA(mytotal)
ElbowPlot(mytotal)
mytotal=RunHarmony(mytotal,group.by.vars = "orig.ident")
mytotal=FindNeighbors(mytotal,dims = 1:15,reduction = "harmony")
mytotal=FindClusters(mytotal,resolution = 0.6)
mytotal=RunUMAP(mytotal,dims = 1:15,reduction = 'harmony')
DimPlot(mytotal)
markertab1=FindAllMarkers(mytotal,only.pos = T)
markertab2=cosg(mytotal)
ptgid=read.csv('ptgid.csv',header = F)
ptgid1=ptgid$V2
names(ptgid1)=levels(mytotal)
namedtotal=RenameIdents(mytotal,ptgid1)
namedtotal$FirstLabel=Idents(namedtotal)
DimPlot(namedtotal)+scale_color_npg()
markertab3=cosg(namedtotal)
sampleid=unique(namedtotal$orig.ident)
myid=c(rep('Pterygium',3),rep('NC',3))
namedtotal$Tissue=mapvalues(namedtotal$orig.ident,from = sampleid,to=myid)
namedtotal@meta.data%>%
ggplot(aes(x = Tissue, fill = FirstLabel)) +
geom_bar(position = position_fill()) +
theme_bw() +
theme(panel.grid = element_blank()) +
labs(y = 'Proportion',fill="Cell Type",x='') +
scale_y_continuous(labels = scales::percent)+theme_classic()+theme(legend.position = 'top') #Barplot of percentage
tab=prop.table(table(namedtotal$orig.ident,namedtotal$FirstLabel))*100
tab$ID=c(rep('Normal',3),rep("Pterygium",3))
tab2=melt(tab,id.vars = 'ID')
ggplot(tab2,aes(x=variable,y=value,fill=ID))+geom_boxplot()+scale_fill_jco()+stat_compare_means(label = "p.signif")+theme_bw()+scale_x_discrete(guide = guide_axis(angle = 45))+theme(legend.position = 'top')+ylab('Percentage(%)')

#Second-round analysis (Detailed analysis for each major cell lineages; Taking fibroblast as an example )
myfib=subset(namedtotal,idents='Fibroblasts')
myfib=NormalizeData(myfib)
myfib=FindVariableFeatures(myfib,nfeatures = 3000)
myfib=ScaleData(myfib)
myfib=RunPCA(myfib)
ElbowPlot(myfib)
myfib1=RunHarmony(myfib,'orig.ident')
myfib1=FindNeighbors(myfib1,reduction = 'harmony',dims=1:20)
myfib1=FindClusters(myfib1,resolution = 0.2)
myfib1=RunUMAP(myfib1,reduction = 'harmony',dims = 1:20)
markertabfib=cosg(myfib1)
markertabfib1=FindAllMarkers(myfib1)
fibid=read.csv('fibid.csv',header = F)
fibid1=fibid$V2
names(fibid1)=levels(myfib1)
namedfib=RenameIdents(myfib1,fibid1)
namedfib$SecondLabel=Idents(namedfib)

#GSVA heatmap
exp=AverageExpression(namedfib )
exp=exp[["RNA"]]
m_df<- msigdbr(species = "human",  category = 'H' )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
GSVA_hall <- gsva(expr=as.matrix(exp),
                  gset.idx.list=geneset,
                  mx.diff=T, 
                  kcdf="Poisson", 
                  parallel.sz=1) 
pheatmap(GSVA_hall[selectedways,],cluster_cols = F,cluster_rows = F,show_rownames = T,show_colnames = T,scale='row')

#monocle analysis
onlyfib=subset(namedfib,idents=c('immuno-Fib','myo-Fib'))
HSMM=onlyfib
data <- as(as.matrix(onlyfib@assays$RNA@counts), 'sparseMatrix')
pd <- new('AnnotatedDataFrame', data = onlyfib@meta.data)
fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
fd <- new('AnnotatedDataFrame', data = fData)
HSMM=newCellDataSet(data,phenoData = pd,featureData = fd)
HSMM <- estimateSizeFactors(HSMM)
HSMM <- estimateDispersions(HSMM)
HSMM <- detectGenes(HSMM, min_expr = 0.1)
disp_table <- dispersionTable(HSMM)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
expressed_genes <- row.names(subset(fData(HSMM), num_cells_expressed >= 25))
diff_test_res <- differentialGeneTest(HSMM[expressed_genes,],
fullModelFormulaStr = "~Tissue")
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))
diff_test_res1=subset(diff_test_res,diff_test_res$num_cells_expressed>100&diff_test_res$qval<0.01)
ordering_genes1=diff_test_res1$gene_short_name
HSMM <- setOrderingFilter(HSMM,ordering_genes = ordering_genes1)
plot_ordering_genes(HSMM)
HSMM <- reduceDimension(HSMM, max_components = 2,
residualModelFormulaStr = '~orig.ident')
HSMM <- orderCells(HSMM)
plot_cell_trajectory(HSMM)
plot_cell_trajectory(HSMM, cell_size = 1.5,color_by = "PI16")+scale_color_gsea()
HSMM <- orderCells(HSMM,root_state = 2)
BEAM_res <- BEAM(HSMM, branch_point = 1, cores = 1)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
df=plot_genes_branched_heatmap(HSMM[row.names(subset(BEAM_res,
                                                  qval < 1e-8)),],
                            branch_point = 1,
                            num_clusters = 3,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
write.csv(df$long.res,'dflongres.csv')

#Cell cycle scoring
onlyfib$State=paste0('State',HSMM$State)
onlyfib$pseudotime=onlyfib$Pseudotime
g2m.genes <- cc.genes.2019$g2m.genes
s.genes <- cc.genes.updated.2019$s.genes
onlyfib <- CellCycleScoring(onlyfib, s.features = s.genes, g2m.features = g2m.genes)
ggboxplot(onlyfib@meta.data,x='State',y='G2M.Score')+ylim(0,0.09)+stat_compare_means(ref.group = '2')

#VolcanoPlot
Idents(onlyfib)=onlyfib$State
state2=subset(onlyfib,idents='State2')
pbmc=state2
pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize", scale.factor = 10000)
pbmc <- FindVariableFeatures(pbmc, selection.method = "vst", nfeatures = 2000)
pbmc <- ScaleData(pbmc)
pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
pbmc <- FindNeighbors(pbmc, dims = 1:15)
pbmc <- FindClusters(pbmc, resolution = 0.2)
Idents(pbmc)=pbmc$Tissue
deg=FindMarkers(pbmc,test.use = 'MAST',ident.1 = 'Pterygium',ident.2 = 'NC',min.pct = 0.01,logfc.threshold = 0.01)
deg$diff.pct=abs(deg$pct.1-deg$pct.2)*100
source('Volcano.R')
keyvals <- ifelse(
deg$avg_log2FC < -0.8& deg$p_val_adj<0.05&deg$diff.pct>0.1, 'darkred',
ifelse(deg$avg_log2FC > 0.8& deg$p_val_adj<0.05&deg$diff.pct>0.1, 'blue','black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'blue'] <- 'upregulated '
names(keyvals)[keyvals == 'darkred'] <- 'downregulated '
EnhancedVolcanov1(deg,x = 'avg_log2FC',y='diff.pct',FCcutoff = 0.8,pCutoff=-1,lab = rownames(deg),ylim = 0:1,colCustom = keyvals,drawConnectors = T,cutoffLineType = 'dotted',selectLab=selcvals)+geom_hline(yintercept = 0.1,linetype='dotted')

#Diffrential GSVA analysis
exp=AverageExpression(state3,group.by =  "orig.ident")
exp=exp[["RNA"]]
exp=AverageExpression(scRNA_harmony )
exp=exp[["RNA"]]
counts2=exp[,c(1:6)]
m_df<- msigdbr(species = "human",  category = "H" )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
GSVA_hall <- gsva(expr=as.matrix(counts2),
                  gset.idx.list=geneSets,
                  kcdf="Poisson", 
                  parallel.sz=1)
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
compare <- makeContrasts(Ptg - NC, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
Diff1=subset(Diff,Diff$P.Value<0.05)
m_df<- msigdbr(species = "human",  category = "C5",subcategory = 'GO:BP' )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
GSVA_hall <- gsva(expr=as.matrix(counts2),
                  gset.idx.list=geneSets,
                  kcdf="Poisson", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=1)
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
compare <- makeContrasts(Ptg - NC, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
Diff3=subset(Diff,Diff$P.Value<0.05)
m_df<- msigdbr(species = "human",  category = "C2",subcategory = 'CP:KEGG' )
geneSets <- m_df %>% split(x = .$gene_symbol, f = .$gs_name)
GSVA_hall <- gsva(expr=as.matrix(counts2),
                  gset.idx.list=geneSets,
                  kcdf="Poisson", #CPM, RPKM, TPM数据就用默认值"Gaussian"， read count数据则为"Poisson"，
                  parallel.sz=1)
design <- model.matrix(~0+group)
colnames(design) = levels(factor(group))
rownames(design) = colnames(GSVA_hall)
compare <- makeContrasts(Ptg - NC, levels=design)
fit <- lmFit(GSVA_hall, design)
fit2 <- contrasts.fit(fit, compare)
fit3 <- eBayes(fit2)
Diff <- topTable(fit3, coef=1, number=200)
Diff4=subset(Diff,Diff$P.Value<0.05)
dat_plot=rbind(Diff1,Diff2,Diff4)%>%arrange(t)
dat_plot1$threshold = factor(ifelse(dat_plot1$t  >-1, ifelse(dat_plot1$t >= 1 ,'Up','NoSignifi'),'Down'),levels=c('Up','Down','NoSignifi'))
write.csv(dat_plot1,'tvalue gsva.csv')
myt=read.csv('tvalue gsva.csv',row.names=1)
dat_plot=myt
low1 <- dat_plot %>% filter(t < -1) %>% nrow()
low0 <- dat_plot %>% filter( t < 0) %>% nrow()
high0 <- dat_plot %>% filter(t < 1) %>% nrow()
high1 <- nrow(dat_plot)
p <- ggplot(data = dat_plot,aes(x = id,y = t,fill = threshold)) +
geom_col()+
coord_flip() +
scale_fill_manual(values = c('Up'= '#36638a','NoSignifi'='#cccccc','Down'='#7bcd7b')) +
geom_hline(yintercept = c(-2,2),color = 'white',size = 0.5,lty='dashed') +
xlab('') +
ylab('t value of GSVA score, Pterigium Vs NC') + 
guides(fill=F)+ 
theme_prism(border = T) +
theme(
axis.text.y = element_blank(),
axis.ticks.y = element_blank()
)
p <- p + geom_text(data = dat_plot[1:low1,],aes(x = id,y = 0.1,label = id),
hjust = 0,color = 'black') + 
geom_text(data = dat_plot[(high0 +1):high1,],aes(x = id,y = -0.1,label = id),
hjust = 1,color = 'black')

#Microarray data analysis
source('CIBERSORT.R')
gset=getGEO('GSE2513',getGPL = F,destdir = '.')
expr=gset[["GSE2513_series_matrix.txt.gz"]]@assayData[['exprs']]
ids=idmap(gpl = 'GPL96')
dat=filterEM(expr,ids)
write.csv(dat,'newdat.csv')
dat1=dat[4:11,]
write.table(dat1,'DATA.txt')
result1 <- CIBERSORT('LM22.txt','DATA.txt', perm = 1000, QN = T)
write.csv(result1,'cib_res.csv')
cib1=read.csv('cib_resv2.csv')
ggscatter(cib1, x = "ACKR1", y = 'Total.Macrophage',
add = "reg.line", conf.int = TRUE,
add.params = list(fill = "lightgray")
)+
stat_cor(method = "spearman",)
mydat1=mydat[,4:11]
group_list=c('high','low','high','low','low','high','low','high')
design=model.matrix(~factor(group_list))
fit=lmFit(mydat1,design )
fit=eBayes(fit)
deg=topTable(fit,coef = 2,n=Tnf)
deg$gene=rownames(deg)
write.csv(deg,'hiackr1vslow.csv')


#Cell-cell interaction analysis
ccptg2seurat=readRDS('ccptgseurat.rds')
future::plan("multicore", workers = 3)
data.input  <- ccptg2seurat@assays$RNA@data
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = ccptg2seurat@meta.data)
cellchat <- setIdent(cellchat, ident.use = "SecondLabel")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ccptg=cellchat
ccnc2seurat=readRDS('ccncseurat.rds')
data.input  <- ccnc2seurat@assays$RNA@data
cellchat <- createCellChat(object = data.input)
cellchat <- addMeta(cellchat, meta = ccnc2seurat@meta.data)
cellchat <- setIdent(cellchat, ident.use = "SecondLabel")
levels(cellchat@idents)
CellChatDB <- CellChatDB.human
cellchat@DB <- CellChatDB
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 3)
cellchat <- computeCommunProbPathway(cellchat)
df.net <- subsetCommunication(cellchat)
cellchat <- aggregateNet(cellchat)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ccnc=cellchat
object.list <- list(Ptg = ccptg, NC = ccnc)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")
par(mfrow = c(1,2), xpd=TRUE)
weight.max <- getMaxWeight(object.list, attribute = c("idents","count"))
for (i in 1:length(object.list)) {
netVisual_circle(object.list[[i]]@net$count, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Weight of interactions - ", names(object.list)[i]))
}
netVisual_bubble(cellchat, sources.use = 1, targets.use = 2:19,  comparison = c(1, 2), angle.x = 45,signaling = mysignals,max.dataset = 1)
pathways.show <- c('VEGF')
weight.max <- getMaxWeight(object.list, slot.name = c("netP"), attribute = pathways.show)
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
netVisual_aggregate(object.list[[i]], signaling = 'VEGF', layout = "circle", edge.weight.max = weight.max[1], edge.width.max = 10, signaling.name = paste('FGF', names(object.list)[i]),vertex.label.cex = 1)
}




















