## 8a03a29901b31176e32928321b1349e6
cat("Gift for HMU wei group. --- Lianhao Song. If any questions, please wechat 18746004617.","\n")
library(plyr)
library(dplyr)
library(Matrix)
## 8a03a29901b31176e32928321b1349e6
scRNA_anlysis <- function(path1 = getwd(),path2 = getwd(),pm = 20,Data_name = "temp",Reso = 0.6,detail = T,nGene_R = c(200,Inf),mito_R = c(-Inf,0.4),PC_M = 7,seed = 233){
  library(Seurat)
  cat(" ","Hello!","Now we focus on:",path1,"\n",file = stderr())
  if(detail){
    PBMC <- Read10X(path1)
    PBMC <- CreateSeuratObject(raw.data = PBMC,min.cells = 3,min.genes = 200,project = "PBMC")
    mito_genes <- grep("^MT\\.",x = rownames(PBMC@data),value = T)
    precent_mito <- colSums(PBMC@raw.data[mito_genes,])/colSums(PBMC@raw.data)
    PBMC <- AddMetaData(object = PBMC, metadata = precent_mito, col.name = "percent_mito")
    par(mfrow = c(1, 2))
    GenePlot(object = PBMC, gene1 = "nUMI", gene2 = "percent_mito", cex.use = 0.2)
    GenePlot(object = PBMC, gene1 = "nUMI", gene2 = "nGene", cex.use = 0.2)
    cat(" ","Now let us cut: \n",file = stderr())
    cat(" ","Please input the low & high thresholds for nGene. If none, input '-Inf' . \n",file = stderr())
    nGene_thre <- scan(sep = ";")
    cat(" ","Please input the low & high thresholds for mito. If none, input '-Inf' . \n",file = stderr())
    mito_thre <- scan(sep = ";")
    dev.off()
    PBMC <- FilterCells(object = PBMC, subset.names = c("nGene", "percent_mito"), low.thresholds = c(nGene_thre[1], mito_thre[1]), high.thresholds = c(nGene_thre[2], mito_thre[2]))
    print(summary(PBMC@raw.data[,1]))
    PBMC <- NormalizeData(object = PBMC, normalization.method = "LogNormalize",scale.factor = 10000)
    print(summary(PBMC@data[,1]))
    PBMC <- FindVariableGenes(object = PBMC, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    gc()
    gc()
    PBMC <- ScaleData(object = PBMC, vars.to.regress = c("nUMI", "percent_mito"))
    print(summary(PBMC@scale.data[,1]))
    gc()
    PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes, do.print = F)
    PBMC <- JackStraw(object = PBMC, num.replicate = 100)
    gc()
    cat(" ","Are you ready ? If ok, input 1",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    JackStrawPlot(object = PBMC, PCs = 1:pm)
    cat(" ","Please save your figure. If ok, input 1",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    print(PCElbowPlot(PBMC))
    cat(" ","Please save your figure. If ok, input 1",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    PCHeatmap(object = PBMC, pc.use = 1:pm, cells.use = 500, do.balanced = TRUE, label.columns = FALSE)
    cat(" ","Please input the highest PC well to use.",file = stderr())
    PCmax <- scan()
    dev.off()
    PBMC <- FindClusters(object = PBMC, reduction.type = "pca", dims.use = 1:PCmax, resolution = Reso, print.output = 0, save.SNN = TRUE)
    gc()
    gc()
    PBMC <- RunTSNE(object = PBMC, dims.use = 1:PCmax)
    TSNEPlot(object = PBMC,do.label = T)
    cat(" ","Please save your figure. If ok, input 1",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    cat(" ","Please confirm your cluster in other R. \n",file = stderr())
    rm(mito_genes,precent_mito,nGene_thre,mito_thre,PCmax)
    gc()
    save(PBMC,file = paste0(path2,Data_name,"_backup.rData"))
    return(PBMC)}
  else{
    PBMC <- Read10X(path1)
    PBMC <- CreateSeuratObject(raw.data = PBMC,min.cells = 3,min.genes = 200,project = "PBMC")
    mito_genes <- grep("^MT\\.",x = rownames(PBMC@data),value = T)
    precent_mito <- colSums(PBMC@raw.data[mito_genes,])/colSums(PBMC@raw.data)
    PBMC <- AddMetaData(object = PBMC, metadata = precent_mito, col.name = "percent_mito")
    PBMC <- FilterCells(object = PBMC, subset.names = c("nGene", "percent_mito"), low.thresholds = c(nGene_R[1], mito_R[1]), high.thresholds = c(nGene_R[2], mito_R[2]))
    print(summary(PBMC@raw.data[,1]))
    PBMC <- NormalizeData(object = PBMC, normalization.method = "LogNormalize",scale.factor = 10000)
    print(summary(PBMC@data[,1]))
    PBMC <- FindVariableGenes(object = PBMC, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    gc()
    PBMC <- ScaleData(object = PBMC, vars.to.regress = c("nUMI", "percent_mito"))
    gc()
    print(summary(PBMC@scale.data[,1]))
    PBMC <- RunPCA(object = PBMC, pc.genes = PBMC@var.genes, do.print = F, pcs.print = 1:5, genes.print = 5)
    PBMC <- FindClusters(object = PBMC, reduction.type = "pca", dims.use = 1:PC_M, resolution = Reso, print.output = 0, save.SNN = TRUE)
    gc()
    PBMC <- RunTSNE(object = PBMC, dims.use = 1:PC_M, seed.use = seed)
    Plot_F <- TSNEPlot(object = PBMC,do.label = T)
    HNSC <- list(PBMC,Plot_F)
    rm(mito_genes,precent_mito,nGene_R,mito_R,PC_M,PBMC,Plot_F)
    gc()
    return(HNSC)}}
cat(" ","scRNA_anlysis --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
Enrich <- function(x,dir = NULL,IDname = dir,Cut = 0.01,Go = T,ReactPA = T,Kegg = T,Keggmap = T,save = T,Gomap = T,wid = 8, h = 8){
  library(clusterProfiler)
  library(ReactomePA)
  path <- getwd()
  dir.create(dir)
  setwd(dir)
  cat(" ","Now working in ",getwd(),"\n",file = stderr())
  GeneID <- bitr(x,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
  ac<-NULL
  bc<-NULL
  PA<-NULL
  if(save){write.csv(GeneID,paste(IDname,"GeneID.csv",sep = "_"))}
  gc()
  if(Go){
    ab <- enrichGO(gene =as.character(GeneID[,2]),OrgDb="org.Hs.eg.db",ont = "ALL",pAdjustMethod = "BH",minGSSize = 1,pvalueCutoff = Cut,readable = TRUE)
    aa <- data.frame(ab)
    rm(ab)
    gc()
    if(save){write.csv(aa,paste(IDname,"GO_enrichment.csv",sep = "_"))}
    ac <- aa[c(which(aa$ONTOLOGY=="BP")[1:10], which(aa$ONTOLOGY=="CC")[1:10], which(aa$ONTOLOGY=="MF")[1:10]),]
    rm(aa)
    gc()
    ac <- na.omit(ac)
    if(Gomap){
      library(ggplot2)
      library(stringr)
      dev.set()
      ggplot(data=ac,aes(ac$Count/length(as.character(GeneID[,2])),ac$Description))+geom_point(aes(size=ac$Count,color=-1*log10(ac$qvalue),shape=ac$ONTOLOGY))+scale_colour_gradient(low="blue",high="red")+labs(color=expression(-log[10](Qvalue)),size="Gene number",shape="Ontology",x="GeneRatio",y=NULL,title="GO enrichment")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+scale_y_discrete(labels = function(x) str_wrap(x,width = 50))
      if(save){ggsave(paste(IDname,"tiff",sep = "."),device = "tiff",width = wid,height = h)}
      dev.off()}
    gc()}
  if(Kegg){
    bb <- enrichKEGG(gene = as.character(GeneID[,2]),organism = "human",pvalueCutoff = Cut)
    bc <- data.frame(bb)
    rm(bb)
    gc()
    if(save)
    {num <- ncol(bc)
    for (i in 1:nrow(bc)) 
    {tryCatch(bc[i,num+1]<-paste(as.character(GeneID[,1])[GeneID[,2] %in% strsplit(bc$geneID,"/")[[i]]],collapse = "/"),error = function(e){write.csv(c("NothingForKEGG"),"ErroKEGG.csv",row.names = F)})}
    rm(num)
    write.csv(bc,paste(IDname,"KEGG_enrichment.csv",sep = "_"))}
    if(Keggmap){
      library(pathview)
      tryCatch(pathview(gene.data = as.character(GeneID[,2]),pathway.id = bc$ID,species = "hsa"),error = function(e){write.csv(c("NothingForPath"),"ErroPATH.csv",row.names = F)})}}
  if(ReactPA){
    PA <- matrix(1:2)
    tryCatch(PA <- enrichPathway(gene = as.character(GeneID[,2]),pvalueCutoff = Cut,organism = "human")@result,error = function(e){print("No Reactome")})
    if(save)
    {num <- ncol(PA)
    for (i in 1:nrow(PA)) 
    {tryCatch(PA[i,num+1]<-paste(as.character(GeneID[,1])[GeneID[,2] %in% strsplit(PA$geneID,"/")[[i]]],collapse = "/"),error = function(e){write.csv(c("NothingForReactPA"),"ErroREACTPA.csv",row.names = F)})}
    rm(num)
    write.csv(bc,paste(IDname,"Reactome_enrichment.csv",sep = "_"))}}
  rm(GeneID)
  aa<-list(ac,bc,PA)
  rm(ac,bc,PA)
  gc()
  gc()
  setwd(path)
  return(aa)}
cat(" ","Enrich --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
remRow <- function(x,Rem=0.1,raito = T){
  x[nrow(x)+1,] <- 0
  a <- as.matrix(x)
  a <- t(apply(a,1,as.numeric))
  r <-as.numeric(apply(a,1,function(i) sum(i == 0)))
  if(raito){remove <- which(r>dim(a)[2]*Rem)}
  else{remove <- which(r>Rem)}
  rm(r)
  Frame <- x[-remove,]
  rm(remove)
  gc()
  return(Frame)}                       
## 8a03a29901b31176e32928321b1349e6
WGCNA_TOMmap <- function(x,nCPU = 5,Cutsample = T,nGene = 10,mGene = 12,minMD = 30,Map = F,custom = F){
  library(WGCNA)
  enableWGCNAThreads(nThreads = nCPU)
  if(!custom){
    if(is.numeric(nGene))
    {FPKM<-t(x[[2]][order(apply(x[[2]],1,mad),decreasing = T)[1:(nGene*1000)],])}
    if(!is.numeric(nGene))
    {FPKM<-t(x[[2]][apply(x[[2]],1,mad)>0,])}}
  if(custom){FPKM <- t(x)}
  gsg<-goodSamplesGenes(FPKM, verbose = 3)
  FPKM<-FPKM[gsg$goodSamples, gsg$goodGenes]
  rm(gsg)
  gc()
  if(Cutsample){
    sampleTree<-hclust(dist(FPKM), method = "average")
    par(cex = 0.6)
    par(mar = c(0,4,2,0))
    plot(sampleTree,main = "Sample clustering to detect outliers",sub="",xlab="",cex.lab = 1.5,cex.axis = 1.5,cex.main = 2)
    print("Welcome to Hatsune world. Now let`s cut the samples")
    okline<-scan()
    abline(h = okline, col = "red")
    print("Can it appropriate?")
    ok<-scan(what = "character")
    while(!ok=="yes")
    {print("Please try again")
      okline<-scan()
      abline(h = okline, col = "red")
      print("Can it appropriate?")
      ok<-scan(what = "character")}
    clust<-cutreeStatic(sampleTree, cutHeight = okline, minSize = 10)
    keepSamples<-(clust==1)
    FPKM<-FPKM[keepSamples,]
    print(table(clust))
    rm(sampleTree,okline,ok,clust,keepSamples)
    gc()}
  powers<-c(c(1:10), seq(from = 12, to=20, by=2))
  sft<-pickSoftThreshold(FPKM, powerVector = powers, verbose = 5,blockSize = mGene*1000)
  par(mfrow = c(1,2))
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],labels=powers,cex=0.9,col="red")
  abline(h=0.90,col="red")
  plot(sft$fitIndices[,1], sft$fitIndices[,5],xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=0.9,col="red")
  okpower<-sft$powerEstimate
  print(paste("The power result is",okpower,". Need to modify ?"))
  judg<-scan(what = "character")
  if(judg=="yes"){
    print("Please input your okpower")
    okpower<-scan()
    print(paste("Now the power result is",okpower))}
  rm(sft,powers,judg)
  gc()
  net<-blockwiseModules(FPKM, power = okpower,TOMType = "unsigned", minModuleSize = minMD, reassignThreshold = 0,mergeCutHeight = 0.25,numericLabels = TRUE, pamRespectsDendro = FALSE,verbose = 3,maxBlockSize=mGene*1000*4/3)
  moduleColors<-labels2colors(net$colors)
  geneTree<-net$dendrograms[[1]]
  MEs<-orderMEs(moduleEigengenes(FPKM, moduleColors)$eigengenes)
  rm(net)
  gc()
  gc()
  if(Map){
    dissTOM<-1-TOMsimilarityFromExpr(FPKM, power = okpower)
    plotTOM<-dissTOM^7
    diag(plotTOM) = NA
    rm(dissTOM)
    gc()
    sizeGrWindow(9,9)
    TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap of all modules")
    rm(plotTOM)}
  rm(geneTree)
  collectGarbage()
  gc()
  print("Now we can check some interesting target")
  Targetgene<-scan(what = "character",sep = ",")
  print(colnames(FPKM)[which(colnames(FPKM) %in% Targetgene)])
  print(moduleColors[which(colnames(FPKM) %in% Targetgene)])
  write.csv(data.frame(Gene = Targetgene, Color = moduleColors[which(colnames(FPKM) %in% Targetgene)]),paste(Targetgene[1],"colorFrame.csv",sep = "_"),row.names = F)
  aa<-list(MEs,moduleColors,okpower,FPKM)
  rm(MEs,moduleColors,okpower,FPKM)
  gc()
  return(aa)}
## 8a03a29901b31176e32928321b1349e6
WGCNA_Detail <- function(x,y,Cliorder=NULL,Trait=NULL,Color=NULL,name = "temp",Mapcolor=Color,MMCutgene=0.20,GSCutgene=0.20,coefficient = 0.02,Cys = T,OnlyCys = F){
  Oripath <- getwd()
  dir.create(name)
  setwd(name)
  if(!OnlyCys){
    weight<-as.data.frame(x[[1]][,Cliorder])
    names(weight)<-Trait
    Names<-substring(names(y[[1]]), 3)
    MM<-as.data.frame(cor(y[[4]], y[[1]], use = "p"))
    MMPvalue<-as.data.frame(corPvalueStudent(as.matrix(MM), nrow(y[[4]])))
    names(MM)<-paste("MM", Names, sep="");
    names(MMPvalue)<-paste("p.MM", Names, sep="")
    GS<-as.data.frame(cor(y[[4]], weight, use = "p"));
    GSPvalue<-as.data.frame(corPvalueStudent(as.matrix(GS), nrow(y[[4]])))
    names(GS)<-paste("GS.", names(weight), sep="");
    names(GSPvalue)<-paste("p.GS.", names(weight), sep="")
    column<-match(Color, Names)
    moduleGenes<-y[[2]]==Color
    sizeGrWindow(7, 7)
    par(mfrow = c(1,1))
    verboseScatterplot(abs(MM[moduleGenes, column]), abs(GS[moduleGenes, 1]),xlab = paste("Module Membership in", Color, "module"), ylab = paste("Gene significance for",Trait),abline = 1,abline.lty = 1,abline.color = Mapcolor,main = paste("Module membership vs. Gene significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Mapcolor)
    abline(h = GSCutgene, v = MMCutgene, col = "red")
    aa<-data.frame(MM=MM[moduleGenes, column],GS=GS[moduleGenes, 1],absMM=abs(MM[moduleGenes, column]),absGS=abs(GS[moduleGenes, 1]),row.names =rownames(GS)[moduleGenes])
    print("Now please enter the Dataclass name. Such as: LGG_FPKMUQ_mRNA")
    name<-scan(what = "character")
    write.csv(aa,paste(name,Trait,Color,"MMandGS.csv",sep = "_"))
    aa<-aa[which(aa[,3] > MMCutgene & aa[,4] > GSCutgene),]
    ac<-rownames(aa)
    rm(aa)
    ab<-rownames(GS)[moduleGenes]
    rm(Names,MM,MMPvalue,GS,GSPvalue,column,moduleGenes)
    gc()
    if(Cys){
      TOM<-TOMsimilarityFromExpr(y[[4]], power = y[[3]])
      inModule<-is.finite(match(y[[2]], Color))
      modGenes<-colnames(y[[4]])[inModule]
      modTOM<-TOM[inModule, inModule]
      dimnames(modTOM) = list(modGenes,modGenes)
      cyt<-exportNetworkToCytoscape(modTOM,edgeFile = paste(name,"edges", paste(Color), "Cys.txt", sep="_"),nodeFile = paste(name,"nodes", paste(Color), "Cys.txt", sep="_"),weighted = TRUE,threshold = coefficient, nodeNames = modGenes, nodeAttr = y[[2]][inModule])
      rm(TOM,inModule,modGenes,modTOM,cyt)
      gc()}
    rm(name)
    aa<-list(ac,ab,weight)
    rm(ac,ab,weight)
    print(aa[[1]])
    gc()
    setwd(Oripath)
    rm(Oripath)
    return(aa)}
  else{
    TOM<-TOMsimilarityFromExpr(y[[4]], power = y[[3]])
    inModule<-is.finite(match(y[[2]], Color))
    modGenes<-colnames(y[[4]])[inModule]
    modTOM<-TOM[inModule, inModule]
    dimnames(modTOM) = list(modGenes,modGenes)
    cyt<-exportNetworkToCytoscape(modTOM,edgeFile = paste(name,"edges", paste(Color), "Cys.txt", sep="_"),nodeFile = paste(name,"nodes", paste(Color), "Cys.txt", sep="_"),weighted = TRUE,threshold = coefficient, nodeNames = modGenes, nodeAttr = y[[2]][inModule])
    nodeGenes <- cyt$nodeData$nodeName
    rm(TOM,inModule,modGenes,modTOM,cyt)
    gc()
    setwd(Oripath)
    rm(Oripath)
    return(nodeGenes)}}
cat(" ","WGNCA --- done.","\n",file = stderr())                       
## 8a03a29901b31176e32928321b1349e6
Pesuo <- function(x,gene = x@var.genes){
  library(monocle)
  HNSC_Ps <- newCellDataSet(as.matrix(x@data))
  HNSC_Ps <- estimateSizeFactors(HNSC_Ps)
  HNSC_Ps <- setOrderingFilter(HNSC_Ps,gene)
  HNSC_Ps <- reduceDimension(HNSC_Ps, max_components=2)
  HNSC_Ps <- orderCells(HNSC_Ps)
  return(HNSC_Ps)}
## 8a03a29901b31176e32928321b1349e6
canFil <- function(x,save = F){
  if(sum(duplicated(substr(colnames(x),1,12)))>0)
  {x<-x[,colnames(x) %in% sort(colnames(x))[!duplicated(substr(sort(colnames(x)),1,12))]]}
  Can<-x[,which(substr(colnames(x),14,14)=="0")]
  Nor<-x[,which(substr(colnames(x),14,14)=="1")]
  colnames(Can) <- gsub("\\.","-",substr(toupper(colnames(Can)),1,12))
  if(save){
    print("Please input the name you want to save. Such as LGG_FPKMUQ_mRNA")
    name<-scan(what = "character")
    write.csv(Can,paste(name,"Can.csv",sep = "_"))
    write.csv(Nor,paste(name,"Nor.csv",sep = "_"))}
  rm(Nor)
  gc()
  return(Can)}
## 8a03a29901b31176e32928321b1349e6                       
DESeq2 <- function(countMatrix, pData){
  dds<-DESeqDataSetFromMatrix(countData = countMatrix, colData = pData, design = ~ phenotype)
  dds<-DESeq(dds)
  dds<-replaceOutliersWithTrimmedMean(dds)
  res<-results(dds, cooksCutoff=FALSE)
  rm(dds)
  res<-res[order(res$padj),]
  return(res)}
## 8a03a29901b31176e32928321b1349e6
DErun <- function(x,y,pvalue = 0.01,log2FC = 2,run = T,save = F){
  library(DESeq2)
  if(run){
    positive_ReadCount<-x
    negative_ReadCount<-y
    colnames(negative_ReadCount)<-paste("N",colnames(negative_ReadCount),sep="-")
    colnames(positive_ReadCount)<-paste("P",colnames(positive_ReadCount),sep="-")
    ReadCount<-cbind(negative_ReadCount,positive_ReadCount)
    ReadCount<-round(ReadCount)
    feature<-c(rep("Neg",ncol(negative_ReadCount)),rep("Pos",ncol(positive_ReadCount)))
    rm(positive_ReadCount,negative_ReadCount)
    gc()
    pData<-data.frame(phenotype = factor(feature,levels=c("Neg", "Pos")))
    rownames(pData)<-colnames(ReadCount)
    diffExp<-DESeq2(ReadCount,pData)
    if(save){
      print("Please enter the name you want to save. Such as: LGG_FPKMUQ_lncRNA_Grade")
      name<-scan(what = "character")
      write.csv(diffExp,paste(name,"padj.csv",sep = "_"))
      rm(name)}
    DEgene<-as.data.frame(diffExp)
    DEname<-rownames(DEgene[DEgene$padj<=pvalue & abs(DEgene$log2FoldChange)>=log2FC,])
    print(sort(DEname))
    rm(diffExp,ReadCount,pData)
    gc()}
  if(!run){
    print("Please enter the name you have saved.Such as: LGG_FPKMUQ_lncRNA")
    named<-scan(what = "character")
    DEgene<-read.csv(paste(named,"padj.csv",sep = "_"),header = T,row.names = 1)
    DEname<-rownames(DEgene[DEgene$padj<pvalue & abs(DEgene$log2FoldChange)>log2FC,])
    print(sort(DEname))
    rm(named)
    gc()}
  aa<-list(DEgene,DEname)
  rm(DEgene,DEname)
  gc()
  return(aa)}
cat(" ","DESeq2 --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6                       
Ttest <- function(x,y,Need_T = F,pval = 0.05){
  if(Need_T){
    x <- t(x)
    y <- t(y)}
  Pval <- c()
  St <- c()
  for (i in 1:nrow(x))
  {Pval[i] <- as.numeric(t.test(as.numeric(x[i,]),as.numeric(y[i,]))$p.value)
   St[i] <- as.numeric(t.test(as.numeric(x[i,]),as.numeric(y[i,]))$statistic)}
  RE <- data.frame(rownames(x),St,Pval,Sig = "No")
  rm(Pval,St)
  RE <- RE[order(RE$Pval),]
  RE$Sig[RE$Pval <= pval] <- "Yes"
  return(RE)}                       
## 8a03a29901b31176e32928321b1349e6
CorTest <- function(x,y,method = "pearson",p_cut = 0.01,adj = T,row = T,name = "Main"){
  if(!row){y <- t(y)}
  if(method == "all"){
    method <- "pearson"
    method2 <- "spearman"
    Corlist <- matrix(0, nrow(y),8)
    for(i in 1:nrow(y)){
      cor <- cor.test(as.numeric(x),as.numeric(y[i,]), method = method, alternative="two.sided")
      p_value <- cor$p.value
      pearson_value <- as.numeric(cor[4])
      Corlist[i,1] <- name
      Corlist[i,2] <- rownames(y)[i]
      Corlist[i,3] <- pearson_value
      Corlist[i,4] <- p_value
      cor <- cor.test(as.numeric(x),as.numeric(y[i,]), method = method2, alternative="two.sided")
      p_value <- cor$p.value
      pearson_value <- as.numeric(cor[4])
      Corlist[i,6] <- pearson_value
      Corlist[i,7] <- p_value}
    Corlist[,5] <- p.adjust(Corlist[,4],method="BH")
    Corlist[,8] <- p.adjust(Corlist[,7],method="BH")
    Corlist <- data.frame(Corlist)
    colnames(Corlist) <- c("Mainname","Corname","Cor","P_value","P_adj","Cor_s","P_value_s","P_adj_s")
    Corlist <- Corlist[order(Corlist$P_adj),]
    Corlist$Sig <- "No"
    if(adj){Corlist$Sig[as.numeric(Corlist$P_adj) < p_cut] <- "Yes"}
    else{Corlist$Sig[as.numeric(Corlist$P_value) < p_cut] <- "Yes"}
    Corlist$State <- "None"
    Corlist$State[as.numeric(Corlist$Cor) > 0] <- "pos"
    Corlist$State[as.numeric(Corlist$Cor) < 0] <- "neg"
    return(Corlist)}
  else{
  Corlist <- matrix(0, nrow(y),5)
  for(i in 1:nrow(y)){
    cor <- cor.test(as.numeric(x),as.numeric(y[i,]), method = method, alternative="two.sided")
    p_value <- cor$p.value
    pearson_value <- as.numeric(cor[4])
    Corlist[i,1] <- name
    Corlist[i,2] <- rownames(y)[i]
    Corlist[i,3] <- pearson_value
    Corlist[i,4] <- p_value}
  Corlist[,5] <- p.adjust(Corlist[,4],method="BH")
  Corlist <- data.frame(Corlist)
  colnames(Corlist) <- c("Mainname","Corname","Cor","P_value","P_adj")
  Corlist <- Corlist[order(Corlist$P_adj),]
  Corlist$Sig <- "No"
  if(adj){Corlist$Sig[as.numeric(Corlist$P_adj) < p_cut] <- "Yes"}
  else{Corlist$Sig[as.numeric(Corlist$P_value) < p_cut] <- "Yes"}
  Corlist$State <- "None"
  Corlist$State[as.numeric(Corlist$Cor) > 0] <- "pos"
  Corlist$State[as.numeric(Corlist$Cor) < 0] <- "neg"
  return(Corlist)}}
cat(" ","Test --- done.","\n",file = stderr()) 
## 8a03a29901b31176e32928321b1349e6
scRNA_3 <- function(x,ori = F,pmax = 20,PCmax = NULL,Reso = 0.5,name = "temp",Dim = 2,detail = T,UMap = F,nVar = 2.5,nICA = 35){
  library(Seurat)
  if(ori){HNSC <- Read10X(x)
  cat(" ","Hello!","Now we focus on:",x,"\n",file = stderr())}
  cat(" ","Hello!","Now we locate at:",getwd(),"\n",file = stderr())
  if(detail){
    HNSC <- CreateSeuratObject(x, name, min.cells = 3, min.features = 200)
    HNSC[["percent.mt"]] <- PercentageFeatureSet(object = HNSC, pattern = "^MT\\.")
    print(VlnPlot(HNSC, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.2))
    print(CombinePlots(list(FeatureScatter(HNSC,"nCount_RNA","percent.mt"),FeatureScatter(HNSC,"nCount_RNA","nFeature_RNA"))))
    cat(" ","Now let us cut: \n",file = stderr())
    cat(" ","Please input the low & high thresholds for nFeature and Mito. If none, input '-Inf' . Such as 200;Inf;-Inf;40 \n",file = stderr())
    HNSC <- subset(HNSC,nFeature_RNA >= scan() & nFeature_RNA <= scan() & percent.mt >= scan() & percent.mt <= scan())
    print(CombinePlots(list(FeatureScatter(HNSC,"nCount_RNA","percent.mt"),FeatureScatter(HNSC,"nCount_RNA","nFeature_RNA"))))
    HNSC <- NormalizeData(HNSC)
    HNSC <- FindVariableFeatures(HNSC, selection.method = "vst", nfeatures = 1000*nVar)
    print(LabelPoints(VariableFeaturePlot(HNSC), points = head(VariableFeatures(HNSC), 10), repel = T))
    HNSC <- ScaleData(HNSC, features = rownames(HNSC))
    HNSC <- RunPCA(HNSC, features = VariableFeatures(HNSC),verbose = F)
    DimHeatmap(HNSC, dims = 1:pmax, cells = 500, balanced = TRUE)
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    gc()
    HNSC <- JackStraw(HNSC, num.replicate = 100)
    HNSC <- ScoreJackStraw(HNSC, dims = 1:(pmax))
    print(JackStrawPlot(HNSC, dims = 1:pmax, xmax = 0.1, ymax = 0.5))
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    print(ElbowPlot(HNSC))
    cat(" ","Please input the highest PC well to use. \n",file = stderr())
    PCmax <- scan()
    gc()
    HNSC <- FindNeighbors(HNSC, dims = 1:PCmax)
    HNSC <- FindClusters(HNSC, resolution = Reso)
    if(UMap){
      HNSC <- RunUMAP(HNSC, dims = 1:PCmax)
      DimPlot(HNSC, reduction = "umap")
      cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
      tem <- scan(what = "character")
      if(!is.null(tem)){cat("well done.\n",file = stderr())}
      rm(tem)}
    HNSC <- RunTSNE(HNSC,dim.embed = Dim,reduction = "pca",dims.use = 1:PCmax)
    HNSC <- RunICA(HNSC,nics = nICA,features = HNSC@active.ident)
    print(DimPlot(HNSC,dims = c(1,2),label = T,reduction = "tsne"))
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    dir.create("backup")
    saveRDS(HNSC,paste0("backup/",name,".rds"))
    gc()
    cat(" ","Please confirm your cluster in other R. \n",file = stderr())
    return(HNSC)}
  else{
    HNSC <- CreateSeuratObject(x, name, min.cells = 3, min.features = 200)
    HNSC[["percent.mt"]] <- PercentageFeatureSet(object = HNSC, pattern = "^MT\\.")
    cat(" ","Please input the low & high thresholds for nFeature and Mito. If none, input '-Inf' . Such as 200;Inf;-Inf;40 \n",file = stderr())
    HNSC <- subset(HNSC,nFeature_RNA >= scan() & nFeature_RNA <= scan() & percent.mt >= scan() & percent.mt <= scan())
    HNSC <- NormalizeData(HNSC,verbose = F)
    HNSC <- FindVariableFeatures(HNSC, selection.method = "vst", nfeatures = 1000*nVar,verbose = F)
    HNSC <- ScaleData(HNSC, features = rownames(HNSC),verbose = F)
    HNSC <- RunPCA(HNSC, features = VariableFeatures(HNSC),verbose = F)
    HNSC <- FindNeighbors(HNSC, dims = 1:PCmax,verbose = F)
    HNSC <- FindClusters(HNSC, resolution = Reso,verbose = F)
    HNSC <- RunTSNE(HNSC,dim.embed = Dim,reduction = "pca",dims.use = 1:PCmax)
    print(DimPlot(HNSC,label = T,reduction = "tsne",pt.size = 0.9))
    Plot <- DimPlot(HNSC,label = T,reduction = "tsne",pt.size = 0.9)
    Fast <- list(HNSC,Plot)
    rm(HNSC,Plot)
    gc()
    return(Fast)}}                       
## 8a03a29901b31176e32928321b1349e6
cat(" ","scRNA_3 --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
Lima <- function(x,y,filt = F,log2FC = 2,padj = 0.01,pval = 0.01){
  Data <- cbind(x,y)
  Group <- data.frame(row.names = colnames(Data), Group1 = c(rep(1,ncol(x)),rep(0,ncol(y))), Group2 = c(rep(0,ncol(x)),rep(1,ncol(y))))
  library(limma)
  Sig <- makeContrasts("Group1-Group2", levels = Group)
  fit <- eBayes(contrasts.fit(lmFit(Data,Group),Sig))
  output <- na.omit(topTable(fit,number = nrow(Data),coef = 1,adjust = "BH"))
  rm(Data,Group,Sig,fit)
  gc()
  if(filt)
  {output <- subset(output, adj.P.Val <= padj & abs(logFC) >= log2FC & P.Value <= pval)}
  output$Sig <- "NO"
  output$Sig[output$P.Value < 0.05] <- "Yes"
  output$State <- "None"
  output$State[output$logFC > 0] <- "Up"
  output$State[output$logFC < 0] <- "Down"
  return(output)}
cat(" ","Lima --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
cat(" ","Ready up. Latest update: 2019-05-16. If any questions, please wechat 18746004617. Email: songlianhao233@gmail.com","\n",file = stderr())
