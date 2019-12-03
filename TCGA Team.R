## TCGA R 2019.10.30 ##
rm(list = ls())
gc()
##

options(stringsAsFactors = F)
setwd("F:/TCGA/")

RNA <- read.table("HTSeq - FPKM-UQ.merge.txt", header = T, row.names = 1,sep = "\t")

source("https://raw.githubusercontent.com/Citruswalker/HMU/master/HMU_will.R")

## First step Rem row 30% ##
RNA <- remRow(RNA,Rem = 0.3)
rownames(RNA) <- gsub("\\..*","",rownames(RNA))

## Second step match genes. ##
lncRNAID <- read.csv("lncRNA_ID.csv",header = T,row.names = 1)
RNAcomID <- intersect(lncRNAID$Ensembl_ID,rownames(RNA))
RNAcom <- RNA[RNAcomID,]
rm(RNA)
rownames(RNAcom) <- paste0(1:nrow(RNAcom),"_",lncRNAID$Symbol[match(RNAcomID,lncRNAID$Ensembl_ID)])

## Third step extrat Cancer samples. ##
RNAcan <- canFil(RNAcom)
rm(RNAcom,lncRNAID,RNAcomID)

## Fourth step clinical. ##
Cli <- read.table("Clinical BCR XML.merge.txt",header = T,sep = "\t")
table(Cli$A4_N)
N0 <- Cli$A0_Samples[Cli$A4_N == "N0"]
NX <- Cli$A0_Samples[!Cli$A4_N == "N0"]

## Fifth step N0 and NX RNA. ##
N0_RNA <- RNAcan[,colnames(RNAcan) %in% N0]
NX_RNA <- RNAcan[,colnames(RNAcan) %in% NX]

## Sixth step DE genes. ##
# 1. Lima #
DE_NXvsN0 <- Lima(NX_RNA,N0_RNA,name = "NXvsN0")
rm(DE_NXvsN0)
gc()
# 2. DEseq2 #
DE_NXvsN0 <- DErun(NX_RNA,N0_RNA,name = "NXvsN0",pvalue = 0.05,log2FC = 0)
DEplot(DE_NXvsN0[[1]],pvalue = 0.05,log2FC = 0,plimit = 10,log2limit = 2.5)
# Nx vs N0
rm(DE_NXvsN0)
gc()

## Sixth 6.5 step extract DE lncRNA. ##
DE_NXvsN0 <- read.csv("NXvsN0_padj.csv",header = T,row.names = 1)
DE_lnc <- rownames(DE_NXvsN0)[which(DE_NXvsN0$padj < 0.05)]

## Seventh step Survival. ##
library(survival)
library(survminer)

lnc_can <- RNAcan[DE_lnc,]

## Use median. ##
lnc_Med <- apply(lnc_can,1,median)

## Build Surfrm ##
Surfrm <- data.frame(Time = Cli$A1_OS, Event = Cli$A2_Event, Sampl = Cli$A0_Samples)
table(Surfrm$Event)
Surfrm <- Surfrm[!Surfrm$Event == "",]
Surfrm$Event <- ifelse(Surfrm$Event == "Dead",1,0)
sort(Surfrm$Time)

## Combine lnc and sur. ##
Sur_Re <- data.frame(Lnc = NA, Pval = NA)
for (i in 1:length(DE_lnc)) {
  message("Now we process ",DE_lnc[i])
  lnc_can[i,] <- ifelse(lnc_can[i,] < lnc_Med[i],0,1)
  lnc_sur <- Surfrm
  lnc_sur$Exp <- as.character(lnc_can[i,])[match(Surfrm$Sampl,colnames(lnc_can))]
  lnc_sur <- na.omit(lnc_sur)
  lnc_sur <- as.data.frame(apply(t(lnc_sur[,-3]), 1, as.numeric))
  diff <- survdiff(Surv(lnc_sur$Time,lnc_sur$Event) ~ lnc_sur$Exp, data = lnc_sur)
  pval <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
  Sur_Re[i,1] <- DE_lnc[i]
  Sur_Re[i,2] <- pval
  rm(diff,pval,lnc_sur)
  gc()}

Sur_Re <- read.csv("Sur_Re.csv",header = T,row.names = 1)
DE_Sur_lnc <- rownames(Sur_Re)[which(Sur_Re$Pval < 0.05)]

## These for plot ##
p <- 1
for (i in DE_Sur_lnc) {
  i <- as.numeric(i)
  print(paste(i,p))
  lnc_can <- RNAcan[DE_lnc,]
  lnc_can[i,] <- ifelse(lnc_can[i,] < lnc_Med[i],0,1)
  lnc_sur <- Surfrm
  lnc_sur$Exp <- as.character(lnc_can[i,])[match(Surfrm$Sampl,colnames(lnc_can))]
  lnc_sur <- na.omit(lnc_sur)
  lnc_sur <- as.data.frame(apply(t(lnc_sur[,-3]), 1, as.numeric))
  fit <- survfit(Surv(lnc_sur$Time,lnc_sur$Event) ~ lnc_sur$Exp, data = lnc_sur)
  print(ggsurvplot(fit,pval = TRUE,xlab = "Time in days",
             legend.labs = c("Low", "High"),xlim = c(0,3000),
             break.time.by = 600,legend.title = element_blank(),
             risk.table = TRUE, surv.median.line = "hv",
             palette = c("blue", "red"),
             ggtheme = theme_light())+ggtitle(gsub(".*_","",DE_lnc[i])))
  scan()
  p <- p+1}
write.csv(Sur_Re,"Sur_Re.csv")
## Plot done. ##

## Eighth step WGCNA. DE_Sur_lncRNA = marks ##

lnc_can <- RNAcan[DE_lnc,]
Marks <- lnc_can[as.numeric(DE_Sur_lnc),]
rm(lnc_can,RNAcan,fit,DE_NXvsN0,lnc_sur,N0_RNA,NX_RNA,Sur_Re,Surfrm)

RNA <- read.table("HTSeq - FPKM-UQ.merge.txt", header = T, row.names = 1,sep = "\t")
RNA <- remRow(RNA,Rem = 0.3)
rownames(RNA) <- gsub("\\..*","",rownames(RNA))

mRNAID <- read.csv("mRNA_ID.csv",header = T,row.names = 1)
RNAcomID <- intersect(mRNAID$Ensembl_ID,rownames(RNA))
RNAcom <- RNA[RNAcomID,]
rm(RNA)
rownames(RNAcom) <- paste0(1:nrow(RNAcom),"_",mRNAID$Symbol[match(RNAcomID,mRNAID$Ensembl_ID)])

RNAcan <- canFil(RNAcom)
rm(RNAcom,mRNAID,RNAcomID)

mRNA_mark <- RNAcan[order(apply(RNAcan, 1, mad),decreasing = T)[1:(5000-11)],] 
rm(RNAcan)
gc()

RNA <- rbind(Marks,mRNA_mark)
write.csv(RNA,"RNA.csv")
rm(list = ls())

## WGCNA ##
Cli <- read.table("Clinical BCR XML.merge.txt",header = T,sep = "\t")
RNA <- read.csv("RNA.csv",header = T,row.names = 1)
WG_data <- WGCNA_CliCustom(Cli,RNA,TCGA = T)
rm(Cli,RNA)
WG_tom <- WGCNA_TOMmap(WG_data,nCPU = 4,nGene = 5,Map = T,Cutsample = F)
# Input: no #

## Process cli. ##
Cli <- WG_data[[1]]
#Cli[6, Cli[6,] == "T1"] <- 1
#Cli[6, Cli[6,] == "T2"] <- 2
#Cli[6, Cli[6,] == "T3"] <- 3
#Cli[6, Cli[6,] %in% c("T4","T4a","T4b")] <- 4
#Cli[6, Cli[6,] == "TX"] <- NA
Cli[7, Cli[7,] == "N0"] <- 0
Cli[7, Cli[7,] == "N1"] <- 1
Cli[7, Cli[7,] %in% c("N2","N2a","N2b","N2c")] <- 2
Cli[7, Cli[7,] == "N3"] <- 3
Cli[7, Cli[7,] == "NX"] <- NA
#Cli[9, Cli[9,] == "Stage I"] <- 1
#Cli[9, Cli[9,] == "Stage II"] <- 2
#Cli[9, Cli[9,] == "Stage III"] <- 3
#Cli[9, Cli[9,] %in% c("Stage IVA","Stage IVB","Stage IVC")] <- 4
#rownames(Cli)[c(6,7,9)] <- paste0(rownames(Cli)[c(6,7,9)]," WG")

WG_data[[1]] <- Cli
##
WG_cli <- WGCNA_CliLink(WG_data,WG_tom)
# Input: WG or only A4_N #
WG_red <- WGCNA_Detail(WG_cli,WG_tom,Cliorder = 1,Color = "red",Trait = "N",name = "Red",MMCutgene = NULL,GSCutgene = NULL,Cys = T)

## Process Cys. ##
cys <- read.table("Red/Red_edges_red_Cys.txt",sep = "\t",header = T)
cys$fromNode <- gsub(".*_","",cys$fromNode)
cys$toNode <- gsub(".*_","",cys$toNode)
write.csv(cys,"Cys_red.csv")
##

## Enrich genes. ##
Gene <- read.table("Red/Red_nodes_red_Cys.txt",header = T,sep = "\t")
rm(list = ls()[!ls() == "Gene"])
Gene <- gsub(".*_","",WG_red[[2]])
source("https://raw.githubusercontent.com/Citruswalker/HMU/master/HMU_will.R")
En_gene <- Enrich(Gene,dir = "Red-ge")
