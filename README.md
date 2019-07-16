# HMU
I am just curious to try to make R package. If you violate the rights, please contact me in time. songlianhao233@gmail.com and 2743623823@qq.com

Make data analysis simple

TCGA analysis:
RNAdir <- "E:/LGG_HTseq_Counts done/"
RNAname <- "Your_RNA_name.txt"
RNA <- read.table(paste0(RNAdir,RNAname),sep = "\t",header = T,row.names = 1)
rownames(RNA) <- gsub("\\..*","",rownames(RNA))
mRNAID_git <- "https://raw.githubusercontent.com/Citruswalker/HMU/master/mRNA_ID.csv"
mRNAID <- read.csv(mRNAID_git,header = T,row.names = 1)
mRNAComID <- intersect(rownames(RNA),mRNAID$Ensembl_ID)
RNAmRNA <- RNA[mRNAComID,]
rownames(RNAmRNA) <- paste0(1:ncol(RNAmRNA),"__",mRNAID$Symbol[match(mRNAComID,mRNAID$Ensembl_ID)])
RNAmRNA <- remRow(RNAmRNA,Rem = 0.1)
source("https://raw.githubusercontent.com/Citruswalker/HMU/master/HMU_will.R")
mRNACan <- canFil(RNAmRNA,save = F)
mRNANor <- canFil(RNAmRNA,can = F,save = F)
DE_mRNA <- DErun(mRNACan,mRNANor,name = "Name")
DEplot(DE_mRNA[[1]],pvalue = 0.01,log2FC = 2,log2limit = 8,plimit = 120,color = 3)
