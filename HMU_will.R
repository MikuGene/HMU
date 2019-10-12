## 8a03a29901b31176e32928321b1349e6
cat("Gift for HMU Wei group, 2019-04-09. --- Lianhao Song. If any questions, please wechat 18746004617.","\n")
if(sum(.packages(all.available=T) %in% "plyr") == 0){install.packages("plyr")}
if(sum(.packages(all.available=T) %in% "dplyr") == 0){install.packages("dplyr")}
if(sum(.packages(all.available=T) %in% "Matrix") == 0){install.packages("Matrix")}
if(sum(.packages(all.available=T) %in% "ggplot2") == 0){install.packages("ggplot2")}
if(sum(.packages(all.available=T) %in% "Seurat") == 0){install.packages("Seurat")}
library(plyr)
library(dplyr)
library(Matrix)
library(ggplot2)
library(Seurat)
## 8a03a29901b31176e32928321b1349e6
ggGene <- function(exp,Target,Iden,l_clor = "#00FFF0",h_clor = "#F600FF",lab_clo = "Median",lab_siz = "Pct",Tle = "Markers",Theme = "NULL",Bline = T){
  Gplot <- list()
  num = 1
  for (i in unique(Iden)){
    clu <- exp[Target,Iden == i]
    pct <- apply(clu, 1, function(i) sum(!i == 0))/ncol(clu)
    Gplot[[num]] <- data.frame(Gene = rownames(clu), Median = apply(clu, 1, median), Mean = apply(clu, 1, mean), Pct = pct, Iden = i)
    num = num+1}
  GP <- do.call(rbind,Gplot)
  GP$Gene <- factor(GP$Gene,levels = unique(GP$Gene))
  Gene <- ggplot(GP,aes(Gene,Iden))+geom_point(aes(size = Pct, color = Median))+scale_colour_gradient(low=l_clor,high=h_clor)
  if(Theme == "light"){Gene <- Gene + theme_light()}
  if(Theme == "NULL"){Gene <- Gene + theme_minimal()}
  if(Theme == "frame"){Gene <- Gene + theme_bw()}
  Gene <- Gene+labs(color= lab_clo,size= lab_siz,x=NULL,y=NULL,title=Tle)+theme(plot.title = element_text(hjust = 0.5))+scale_y_discrete(labels = function(x) str_wrap(x,width = 50))
  if(!Bline){Gene <- Gene + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())}
  print(Gene)}
## 8a03a29901b31176e32928321b1349e6
scRNA_2 <- function(path1 = getwd(),path2 = getwd(),mito_name = "^MT\\.",pm = 20,Data_name = "temp",Reso = 0.6,detail = T,nGene_R = c(200,Inf),mito_R = c(-Inf,0.4),PC_M = 7,seed = 233){
  library(Seurat)
  cat(" ","Hello!","Now we focus on:",path1,"\n",file = stderr())
  if(detail){
    PBMC <- Read10X(path1)
    PBMC <- CreateSeuratObject(raw.data = PBMC,min.cells = 3,min.genes = 200,project = "PBMC")
    mito_genes <- grep(mito_name,x = rownames(PBMC@data),value = T)
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
    mito_genes <- grep(mito_name,x = rownames(PBMC@data),value = T)
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
Enrich <- function(x,dir = "temp",IDname = dir,Cut = 0.01,Go = T,ReactPA = T,Kegg = T,Keggmap = T,save = T,Gomap = T,wid = 8, h = 8){
  if(sum(.packages(all.available=T) %in% "clusterProfiler") == 0){install.packages("clusterProfiler")}
  if(sum(.packages(all.available=T) %in% "ReactomePA") == 0){install.packages("ReactomePA")}
  library(clusterProfiler)
  library(ReactomePA)
  path <- getwd()
  dir.create(dir)
  setwd(dir)
  cat(" ","Now working in ",getwd(),"\n",file = stderr())
  GeneT <- gsub("^MT\\.","MT-",x)
  GeneID <- bitr(GeneT,fromType = "SYMBOL",toType = "ENTREZID",OrgDb = "org.Hs.eg.db")
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
      if(sum(.packages(all.available=T) %in% "stringr") == 0){install.packages("stringr")}
      library(stringr)
      print(ggplot(data=ac,aes(ac$Count/length(as.character(GeneID[,2])),ac$Description))+geom_point(aes(size=ac$Count,color=-1*log10(ac$qvalue),shape=ac$ONTOLOGY))+scale_colour_gradient(low="blue",high="red")+labs(color=expression(-log[10](Qvalue)),size="Gene number",shape="Ontology",x="GeneRatio",y=NULL,title="GO enrichment")+theme_bw()+theme(plot.title = element_text(hjust = 0.5))+scale_y_discrete(labels = function(x) str_wrap(x,width = 50)))
      if(save){ggsave(paste(IDname,"tiff",sep = "."),device = "tiff",width = wid,height = h)}}
    gc()}
  if(Kegg){
    bb <- enrichKEGG(gene = as.character(GeneID[,2]),organism = "hsa",pvalueCutoff = Cut)
    bc <- data.frame(bb)
    rm(bb)
    gc()
    if(save)
    {num <- ncol(bc)
    for (i in 1:nrow(bc)) 
    {tryCatch(bc[i,num+1]<-paste(as.character(GeneID[,1])[GeneID[,2] %in% strsplit(bc$geneID,"/")[[i]]],collapse = "/"),error = function(e){write.csv(c("NothingForKEGG"),"ErroKEGG.csv",row.names = F)})}
    rm(num)
    write.csv(bc,paste0(IDname,"_KEGG_enrichment.csv"))}
    if(Keggmap){
      if(sum(.packages(all.available=T) %in% "pathview") == 0){install.packages("pathview")}
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
    write.csv(PA,paste0(IDname,"_Reactome_enrichment.csv"))}}
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
WGCNA_CliCustom <- function(x,y,IDname = "A0_Samples",Need_t = T,save = F,TCGA = F){
  if(Need_t){x <- data.frame(t(x))}
  if(TCGA){
  colnames(y)<-gsub("\\.","-",substr(toupper(colnames(y)),1,12))
  x[IDname,]<-gsub("\\.","-",substr(toupper(x[IDname,]),1,12))}
  ComID<-intersect(colnames(y),x[IDname,])
  FPKMRows<-match(ComID,colnames(y))
  FPKMSim<-y[,FPKMRows]
  traitRows<-match(ComID,x[IDname,])
  simpleCli<-x[-which(rownames(x)==IDname),traitRows]
  colnames(simpleCli) = x[which(rownames(x)==IDname),traitRows]
  rownames(simpleCli) = rownames(x)[-which(rownames(x)==IDname)]
  rm(ComID,FPKMRows,traitRows)
  if(save){
  print("Now please write WGCNA name. Such as: LGG_FPKMUQ_mRNA")
  name<-scan(what = "character")
  write.csv(simpleCli,paste(name,"WGCNA_clinical.csv",sep = "_"))
  write.csv(FPKMSim,paste(name,"WGCNA_FPKM.csv",sep = "_"))
  rm(name)}
  gc()
  a<-list(simpleCli,FPKMSim)
  rm(simpleCli,FPKMSim)
  gc()
  return(a)}
## 8a03a29901b31176e32928321b1349e6
WGCNA_TOMmap <- function(x,nCPU = 5,Cutsample = T,nGene = 10,mGene = 12,minMD = 30,Map = F,custom = F){
  if(sum(.packages(all.available=T) %in% "WGCNA") == 0){install.packages("WGCNA")}
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
  tomP <- NULL
  if(Map){
    dissTOM<-1-TOMsimilarityFromExpr(FPKM, power = okpower)
    plotTOM<-dissTOM^7
    diag(plotTOM) = NA
    rm(dissTOM)
    gc()
    sizeGrWindow(9,9)
    tomP <- TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap of all modules")
    rm(plotTOM)}
  gc()
  print("Now we can check some interesting target")
  Targetgene<-scan(what = "character",sep = ",")
  print(colnames(FPKM)[which(colnames(FPKM) %in% Targetgene)])
  print(moduleColors[which(colnames(FPKM) %in% Targetgene)])
  write.csv(data.frame(Gene = Targetgene, Color = moduleColors[which(colnames(FPKM) %in% Targetgene)]),paste(Targetgene[1],"colorFrame.csv",sep = "_"),row.names = F)
  aa<-list(MEs,moduleColors,okpower,FPKM,tomP)
  rm(MEs,moduleColors,okpower,FPKM,tomP,geneTree)
  gc()
  return(aa)}
## 8a03a29901b31176e32928321b1349e6
WGCNA_CliLink <- function(x,y,xais = T,yais = T,plot = T){
  if(sum(.packages(all.available=T) %in% "WGCNA") == 0){install.packages("WGCNA")}
  library(WGCNA)
  print("Now please enter one class name. Such as: Grade")
  class <- scan(what = "character")
  aa<-x[[1]][match(grep(class,rownames(x[[1]]),value = T),rownames(x[[1]])),match(rownames(y[[4]]),colnames(x[[1]]))]
  aa<-apply(aa, 1, as.numeric)
  TraitCor<-cor(y[[1]], aa, use = "p")
  TraitPvalue<-corPvalueStudent(TraitCor, nrow(y[[1]]))
  textMatrix<-paste(signif(TraitCor, 2), "\n(",signif(TraitPvalue, 1), ")", sep = "")
  dim(textMatrix) = dim(TraitCor)
  if(xais){name = grep(class,rownames(x[[1]]),value = T)}else{name = NULL}
  if(yais){name2 = colnames(y[[1]])}else{name2 = NULL}
  if(plot){labeledHeatmap(Matrix=TraitCor,xLabels = name,yLabels = name2,colorLabels = F,colors = blueWhiteRed(50),textMatrix = textMatrix,setStdMargins = FALSE,cex.text = 0.5,zlim = c(-1,1),main = paste("Module-trait relationships"))}
  rm(TraitCor,TraitPvalue,textMatrix)
  a<-list(aa,substring(names(y[[1]]),3),x[[1]])
  rm(aa)
  print(a[[2]])
  gc()
  return(a)}
## 8a03a29901b31176e32928321b1349e6
WGCNA_Detail <- function(x,y,Cliorder=NULL,plot = T,custom = F,Trait="temp",Color=NULL,name = "temp",Mapcolor=Color,MMCutgene=0.20,GSCutgene=0.20,coefficient = 0.02,Cys = T,OnlyCys = F){
  Oripath <- getwd()
  dir.create(name)
  setwd(name)
  if(!OnlyCys){
    if(!custom){weight<-as.data.frame(x[[1]][,Cliorder])}
    else{weight<-x}
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
    if(plot){verboseScatterplot(abs(MM[moduleGenes, column]), abs(GS[moduleGenes, 1]),xlab = paste("Module Membership in", Color, "module"), ylab = paste("Gene significance for",Trait),abline = 1,abline.lty = 1,abline.color = Mapcolor,main = paste("Module membership vs. Gene significance\n"),cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = Mapcolor)
    abline(h = GSCutgene, v = MMCutgene, col = "red")}
    aa<-data.frame(MM=MM[moduleGenes, column],GS=GS[moduleGenes, 1],absMM=abs(MM[moduleGenes, column]),absGS=abs(GS[moduleGenes, 1]),row.names =rownames(GS)[moduleGenes])
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
Pesuo <- function(x,gene = x@var.genes,Re = F){
  if(sum(.packages(all.available=T) %in% "monocle") == 0){install.packages("monocle")}
  library(monocle)
  HNSC_Ps <- newCellDataSet(as.matrix(x@data))
  HNSC_Ps <- estimateSizeFactors(HNSC_Ps)
  HNSC_Ps <- setOrderingFilter(HNSC_Ps,gene)
  HNSC_Ps <- reduceDimension(HNSC_Ps, max_components=2)
  HNSC_Ps <- orderCells(HNSC_Ps, reverse = Re)
  return(HNSC_Ps)}
## 8a03a29901b31176e32928321b1349e6
canFil <- function(x,can = T,save = F){
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
  gc()
  if(!can){rm(Can)
  return(Nor)}
  else{rm(Nor)
  return(Can)}}
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
DErun <- function(x,y,pvalue = 0.01,log2FC = 2,run = T,save = T,name = "temp"){
  if(sum(.packages(all.available=T) %in% "DESeq2") == 0){install.packages("DESeq2")}
  library(DESeq2)
  if(run){
    pos<-x
    neg<-y
    colnames(neg)<-paste("N",colnames(neg),sep="-")
    colnames(pos)<-paste("P",colnames(pos),sep="-")
    ReadCount<-round(cbind(neg,pos))
    feature<-c(rep("Neg",ncol(neg)),rep("Pos",ncol(pos)))
    rm(pos,neg)
    gc()
    pData<-data.frame(phenotype = factor(feature,levels=c("Neg", "Pos")))
    rownames(pData)<-colnames(ReadCount)
    diffExp<-DESeq2(ReadCount,pData)
    if(save){write.csv(diffExp,paste(name,"padj.csv",sep = "_"))}
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
ggpie <- function(x,Title = "Cell composition",Size = 3.5,Theme = "light"){
  Iden <- data.frame(table(x))
  Iden$Var2 <- factor(c(1:nrow(Iden)))
  Iden_label <- paste0(as.vector(Iden[,1]), " (", round(Iden$Freq / sum(Iden$Freq) * 100, 2), "%)")
  pie <- ggplot(Iden, aes(x = "", y = Iden$Freq, fill = Iden$Var2))+geom_bar(stat = "identity")+coord_polar(theta = "y")+labs(x = "", y = "", title = Title)
  if(Theme == "light"){pie <- pie + theme_light()}
  if(Theme == "NULL"){pie <- pie + theme_minimal()}
  if(Theme == "frame"){pie <- pie + theme_bw()}
  pie <- pie + theme(axis.ticks = element_blank(),axis.text.x = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())+scale_fill_discrete(breaks = Iden$Var2,labels = Iden_label)+theme(legend.title = element_blank(),plot.title = element_text(hjust = 0.5))+geom_text(aes(x = 1.5,y = rev(c(0,cumsum(rev(Iden$Freq))[-length(Iden$Freq)])+rev(Iden$Freq/2)),label = as.vector(Iden$Var2)),size = Size)
  print(pie)}
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
CorTest <- function(x,y,method = "pearson",p_cut = 0.01,adj = T,row = T,name = "Main",Order = F){
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
    if(Order){Corlist <- Corlist[order(Corlist$P_adj),]}
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
## 8a03a29901b31176e32928321b1349e6
CrossCor <- function(x,row = T){
  Re_Cor_p <- matrix(rep(NA,nrow(x)^2),nrow = nrow(x))
  Re_Pval_p <- matrix(rep(NA,nrow(x)^2),nrow = nrow(x))
  Re_Padj_p <- matrix(rep(NA,nrow(x)^2),nrow = nrow(x))
  Re_Cor_s <- matrix(rep(NA,nrow(x)^2),nrow = nrow(x))
  Re_Pval_s <- matrix(rep(NA,nrow(x)^2),nrow = nrow(x))
  Re_Padj_s <- matrix(rep(NA,nrow(x)^2),nrow = nrow(x))
  Re_G <- data.frame()
  if(!row){x <- t(x)}
  for (i in 1:nrow(x)) {
      cat("Cor","Row",i,"\n")
      Re <- CorTest(x[i,],x,name = rownames(x[i,]),method = "all")
      Re_Cor_p[i,] <- Re$Cor
      Re_Pval_p[i,] <- Re$P_value
      Re_Padj_p[i,] <- Re$P_adj
      Re_Cor_s[i,] <- Re$Cor_s
      Re_Pval_s[i,] <- Re$P_value_s
      Re_Padj_s[i,] <- Re$P_adj_s
      Re_G <- rbind(Re_G,Re)
      rm(Re)
      gc()
  cat("Cor","Row",i,"done.","\n",file = stderr())}
  Relist = list(Re_Cor_p,Re_Padj_p,Re_Pval_p,Re_Cor_s,Re_Padj_s,Re_Pval_s,Re_G)
  return(Relist)}
## 8a03a29901b31176e32928321b1349e6
cat(" ","Test --- done.","\n",file = stderr()) 
## 8a03a29901b31176e32928321b1349e6
scRNA_3 <- function(x,y = NULL,Anti = F,if_two = F,if_plot = T,name1 = "temp1_sc",name2 = "temp2_sc",ori = F,Mito = c("^MT\\.","^MT-"),pmax = 20,PCmax = NULL,Reso = 0.6,name = "temp",Dim = 2,detail = T,UMap = F,nVar = 2.5,all_Anc = F,if_var = F,Vars = c("nFeature_RNA","percent.mt"),PCjk = T){
  cat(" ","Hello!","Now we locate at:",getwd(),"\n",file = stderr())
  if(ori){
    cat(" ","Hello!","Now we focus on:",x,"\n",file = stderr())
    x <- Read10X(x)
    if(Anti){
      write.csv(x[[2]],paste0(name,"_Antibody_capture.csv"))
      x <- x[[1]]
      gc()}}
  if(detail){
    if(!if_two){
    HNSC <- CreateSeuratObject(x, name, min.cells = 3, min.features = 200)
    HNSC[["percent.mt"]] <- PercentageFeatureSet(object = HNSC, pattern = Mito[1])
    if(if_plot){print(VlnPlot(HNSC, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.2))}
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    gc()
    if(if_plot){print(CombinePlots(list(FeatureScatter(HNSC,"nCount_RNA","percent.mt"),FeatureScatter(HNSC,"nCount_RNA","nFeature_RNA"))))}
    cat(" ","Now let us cut: \n",file = stderr())
    cat(" ","Please input the low & high thresholds for nFeature and Mito. If none, input '-Inf' . Such as 200;Inf;-Inf;40 \n",file = stderr())
    HNSC <- subset(HNSC,nFeature_RNA >= scan() & nFeature_RNA <= scan() & percent.mt >= scan() & percent.mt <= scan())
    if(if_plot){print(CombinePlots(list(FeatureScatter(HNSC,"nCount_RNA","percent.mt"),FeatureScatter(HNSC,"nCount_RNA","nFeature_RNA"))))}
    HNSC <- NormalizeData(HNSC)
    HNSC <- FindVariableFeatures(HNSC, selection.method = "vst", nfeatures = 1000*nVar)
    if(if_plot){print(LabelPoints(VariableFeaturePlot(HNSC), points = head(VariableFeatures(HNSC), 10), repel = T))}}
    else{
      HNSC1 <- CreateSeuratObject(x, name1, min.cells = 3, min.features = 200)
      HNSC1$Group <- name1
      HNSC1[["percent.mt"]] <- PercentageFeatureSet(object = HNSC1, pattern = Mito[1])
      if(if_plot){print(VlnPlot(HNSC1, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.2))}
      cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
      tem <- scan(what = "character")
      if(!is.null(tem)){cat("well done.\n",file = stderr())}
      rm(tem)
      gc()
      if(if_plot){print(CombinePlots(list(FeatureScatter(HNSC1,"nCount_RNA","percent.mt"),FeatureScatter(HNSC1,"nCount_RNA","nFeature_RNA"))))}
      cat(" ","Now let us cut: \n",file = stderr())
      cat(" ","Please input the low & high thresholds for nFeature and Mito. If none, input '-Inf' . Such as 200;Inf;-Inf;40 \n",file = stderr())
      HNSC1 <- subset(HNSC1,nFeature_RNA >= scan() & nFeature_RNA <= scan() & percent.mt >= scan() & percent.mt <= scan())
      if(if_plot){print(CombinePlots(list(FeatureScatter(HNSC1,"nCount_RNA","percent.mt"),FeatureScatter(HNSC1,"nCount_RNA","nFeature_RNA"))))}
      HNSC1 <- NormalizeData(HNSC1)
      HNSC1 <- FindVariableFeatures(HNSC1, selection.method = "vst", nfeatures = 1000*nVar)
      if(if_plot){print(LabelPoints(VariableFeaturePlot(HNSC1), points = head(VariableFeatures(HNSC1), 10), repel = T))}
      cat(" ","Now First one done. \n",file = stderr())
      cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
      tem <- scan(what = "character")
      if(!is.null(tem)){cat("well done.\n",file = stderr())}
      rm(tem)
      gc()
      SCC090 <- CreateSeuratObject(y, name2, min.cells = 3, min.features = 200)
      SCC090$Group <- name2
      SCC090[["percent.mt"]] <- PercentageFeatureSet(object = SCC090, pattern = Mito[2])
      if(if_plot){print(VlnPlot(SCC090, c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,pt.size = 0.2))}
      cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
      tem <- scan(what = "character")
      if(!is.null(tem)){cat("well done.\n",file = stderr())}
      rm(tem)
      gc()
      if(if_plot){print(CombinePlots(list(FeatureScatter(SCC090,"nCount_RNA","percent.mt"),FeatureScatter(SCC090,"nCount_RNA","nFeature_RNA"))))}
      cat(" ","Now let us cut: \n",file = stderr())
      cat(" ","Please input the low & high thresholds for nFeature and Mito. If none, input '-Inf' . Such as 200;Inf;-Inf;40 \n",file = stderr())
      SCC090 <- subset(SCC090,nFeature_RNA >= scan() & nFeature_RNA <= scan() & percent.mt >= scan() & percent.mt <= scan())
      if(if_plot){print(CombinePlots(list(FeatureScatter(SCC090,"nCount_RNA","percent.mt"),FeatureScatter(SCC090,"nCount_RNA","nFeature_RNA"))))}
      SCC090 <- NormalizeData(SCC090)
      SCC090 <- FindVariableFeatures(SCC090, selection.method = "vst", nfeatures = 1000*nVar)
      if(if_plot){print(LabelPoints(VariableFeaturePlot(SCC090), points = head(VariableFeatures(SCC090), 10), repel = T))}
      cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
      tem <- scan(what = "character")
      if(!is.null(tem)){cat("well done.\n",file = stderr())}
      rm(tem)
      gc()
      if(all_Anc){Anchors <- FindIntegrationAnchors(list(HNSC1, SCC090),anchor.features = intersect(rownames(x),rownames(y)))}
      else{Anchors <- FindIntegrationAnchors(list(HNSC1, SCC090))}
      HNSC <- IntegrateData(anchorset = Anchors)
      DefaultAssay(HNSC) <- "integrated"}
    if(if_var){HNSC <- ScaleData(HNSC, features = rownames(HNSC), vars.to.regress = Vars)}
    else{HNSC <- ScaleData(HNSC, features = rownames(HNSC))}
    HNSC <- RunPCA(HNSC, features = VariableFeatures(HNSC),verbose = F)
    if(if_plot){DimHeatmap(HNSC, dims = 1:pmax, cells = 500, balanced = TRUE)}
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    gc()
    if(PCjk){
    HNSC <- JackStraw(HNSC, num.replicate = 100)
    HNSC <- ScoreJackStraw(HNSC, dims = 1:(pmax))
    if(if_plot){print(JackStrawPlot(HNSC, dims = 1:pmax, xmax = 0.1, ymax = 0.5))}
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)}
    if(if_plot){print(ElbowPlot(HNSC))}
    cat(" ","Please input the highest PC well to use. \n",file = stderr())
    PCmax <- scan()
    gc()
    HNSC <- FindNeighbors(HNSC, dims = 1:PCmax)
    HNSC <- FindClusters(HNSC, resolution = Reso)
    if(UMap){
      HNSC <- RunUMAP(HNSC, dims = 1:PCmax)
      if(if_plot){DimPlot(HNSC, reduction = "umap")}
      cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
      tem <- scan(what = "character")
      if(!is.null(tem)){cat("well done.\n",file = stderr())}
      rm(tem)}
    HNSC <- RunTSNE(HNSC,dim.embed = Dim,reduction = "pca",dims.use = 1:PCmax)
    if(if_plot){print(DimPlot(HNSC,dims = c(1,2),label = T,reduction = "tsne"))}
    cat(" ","Please save your figure. If ok, input 1 \n",file = stderr())
    tem <- scan(what = "character")
    if(!is.null(tem)){cat("well done.\n",file = stderr())}
    rm(tem)
    dir.create("backup")
    saveRDS(HNSC,paste0("backup/",name,"_backup.rds"))
    gc()
    cat(" ","Please confirm your cluster in other R. \n",file = stderr())
    return(HNSC)}
  else{
    HNSC <- CreateSeuratObject(x, name, min.cells = 3, min.features = 200)
    HNSC[["percent.mt"]] <- PercentageFeatureSet(object = HNSC, pattern = Mito[1])
    cat(" ","Please input the low & high thresholds for nFeature and Mito. If none, input '-Inf' . Such as 200;Inf;-Inf;40 \n",file = stderr())
    HNSC <- subset(HNSC,nFeature_RNA >= scan() & nFeature_RNA <= scan() & percent.mt >= scan() & percent.mt <= scan())
    HNSC <- NormalizeData(HNSC,verbose = F)
    HNSC <- FindVariableFeatures(HNSC, selection.method = "vst", nfeatures = 1000*nVar,verbose = F)
    if(if_var){HNSC <- ScaleData(HNSC, features = rownames(HNSC), vars.to.regress = Vars,verbose = F)}
    else{HNSC <- ScaleData(HNSC, features = rownames(HNSC),verbose = F)}
    HNSC <- RunPCA(HNSC, features = VariableFeatures(HNSC),verbose = F)
    HNSC <- FindNeighbors(HNSC, dims = 1:PCmax,verbose = F)
    HNSC <- FindClusters(HNSC, resolution = Reso,verbose = F)
    HNSC <- RunTSNE(HNSC,dim.embed = Dim,reduction = "pca",dims.use = 1:PCmax)
    if(if_plot){print(DimPlot(HNSC,label = T,reduction = "tsne",pt.size = 0.9))}
    Plot <- DimPlot(HNSC,label = T,reduction = "tsne",pt.size = 0.9)
    Fast <- list(HNSC,Plot)
    rm(HNSC,Plot)
    gc()
    return(Fast)}}                     
## 8a03a29901b31176e32928321b1349e6
cat(" ","scRNA_3 --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
Lima <- function(x,y,filt = F,log2FC = 2,padj = 0.01,pval = 0.01,save = T,Order = T,name = "temp"){
  Data <- cbind(x,y)
  Group <- data.frame(row.names = colnames(Data), Group1 = c(rep(1,ncol(x)),rep(0,ncol(y))), Group2 = c(rep(0,ncol(x)),rep(1,ncol(y))))
  if(sum(.packages(all.available=T) %in% "limma") == 0){install.packages("limma")}
  library(limma)
  Sig <- makeContrasts("Group1-Group2", levels = Group)
  fit <- eBayes(contrasts.fit(lmFit(Data,Group),Sig))
  OP <- na.omit(topTable(fit,number = nrow(Data),coef = 1,adjust = "BH"))
  rm(Data,Group,Sig,fit)
  gc()
  colnames(OP) <- c("LogFC","AveExpr","T","P_val","P_adj","B")
  OP$MeanP <- apply(x,1,mean)
  OP$MeanN <- apply(y,1,mean)
  OP$MedianP <- apply(x,1,median)
  OP$MedianN <- apply(y,1,median)
  if(filt){OP <- subset(OP, P_adj <= padj & abs(LogFC) >= log2FC & P_val <= pval)}
  OP$Sig <- "NO"
  OP$Sig[OP$P_val < 0.05] <- "Yes"
  OP$State <- "None"
  OP$State[OP$LogFC > 0] <- "Up"
  OP$State[OP$LogFC < 0] <- "Down"
  if(Order){OP <- OP[order(abs(OP$LogFC)),]}
  if(save){write.csv(OP,paste(name,"lima.csv",sep = "_"))}
  return(OP)}
cat(" ","Lima --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
ggpoint <- function(Data,x,y,size = x,clor = y,l_clor = "grey",h_clor = "red",lab_clor = "Pct",lab_siz = "LogFC",Tle = "Markers",sort = "x",Decr = F,Theme = "NULL",Bline = T){
  if(sum(.packages(all.available=T) %in% "stringr") == 0){install.packages("stringr")}
  library(stringr)
  if(sort == "x"){Data <- Data[order(Data[,x],decreasing = Decr),]}
  if(sort == "y"){Data <- Data[order(Data[,y],decreasing = Decr),]}
  Data[,y] <- factor(Data[,y],levels = unique(Data[,y]))
  point <- ggplot(Data,aes(Data[,x],Data[,y]))+geom_point(aes(size = Data[,size], color = Data[,clor]))+scale_colour_gradient(low=l_clor,high=h_clor)
  if(Theme == "light"){point <- point + theme_light()}
  if(Theme == "NULL"){point <- point + theme_minimal()}
  if(Theme == "frame"){point <- point + theme_bw()}
  point <- point+labs(color=lab_clor,size=lab_siz,x=NULL,y=NULL,title=Tle)+theme(plot.title = element_text(hjust = 0.5))+scale_y_discrete(labels = function(x) str_wrap(x,width = 50))
  if(!Bline){point <- point + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank())}
  print(point)}
cat(" ","ggplot --- done.","\n",file = stderr())
## 8a03a29901b31176e32928321b1349e6
Exct <- function(x,exct = "\\^",filed = 1){
  print(paste0("Before: ",x[1]))
  x <- sapply(strsplit(x,exct),function(i) i[[filed]])
  x <- unique(x[nchar(x) > 0])
  print(paste0("After: ",x[1]))
  return(x)}
## 8a03a29901b31176e32928321b1349e6
gene_line <- function(x,dig = 2){
x <- round(x,digits = dig)
DTB <- data.frame(table(x))
colnames(DTB) <- c("Val","Freq")
DTB$Val <- as.numeric(as.character(DTB$Val))
DTB$Freq <- as.numeric(DTB$Freq)
plot(DTB$Val,DTB$Freq)+lines(DTB$Val,predict(loess(DTB$Freq ~ as.numeric(DTB$Val))))}
## 8a03a29901b31176e32928321b1349e6
DEplot <- function(x, pvalue = 0.01, log2FC = 2, plimit = 30, log2limit = 5, color = 3,Lima = F,adj = T){
  if(Lima){
    if(adj){colnames(x) <- c("log2FoldChange","AveExpr","t","p","padj","B")}
    else{colnames(x) <- c("log2FoldChange","AveExpr","t","padj","p","B")}}
  x$Legend<-as.factor(ifelse(x$padj<pvalue & abs(x$log2FoldChange)>=log2FC, ifelse(x$log2FoldChange>log2FC,'Up','Down'),'Not'))
  print("Please enter your Title name")
  Title<-scan(what = "character",sep = ",")
  if(color == 3){colornum <- c("blue", "black", "red")}
  if(color == 2){colornum <- c("black", "red")}
  print(ggplot(data=x,aes(x=log2FoldChange, y=-log10(padj),colour=Legend))+ggtitle(Title)+xlab("log2 Foldchange")+ylab("-log10 Padj")+geom_vline(xintercept=c(-log2FC,log2FC),lty=6,col="grey",lwd=0.5)+geom_hline(yintercept = -log10(pvalue),lty=4,col="grey",lwd=0.5)+scale_color_manual(values = colornum)+theme(legend.position="right")+theme_bw()+theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),legend.title = element_blank())+xlim(-log2limit,log2limit) + ylim(0,plimit)+theme(plot.title = element_text(hjust = 0.5))+geom_point(alpha=0.4, size=1.2))}
## 8a03a29901b31176e32928321b1349e6
cat(" ","Ready up. Latest update: 2019-10-12-10:31 --- Lianhao Song.","\n","","---If any questions, please wechat 18746004617. Email: songlianhao233@gmail.com","\n",file = stderr())
