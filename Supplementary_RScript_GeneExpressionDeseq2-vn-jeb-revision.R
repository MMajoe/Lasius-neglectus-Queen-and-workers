#'Lasius neglectus Updated analysis final R script:
#'Gene expression data:
#'#"representative" genome/ proteome used for htseq Counts
#Thus only one amino acid/CDS region/isoform associated with every gene
#
#PCAs plotted used in Fig 2 A, C
#Deseq 2 Results compiled in Supplementary table 6 and for further analyses
#posterior log2FC analyses then plotted as Supplementary Fig 2.
#Vg's plotted for Figure 3
#
#
#
#'Packages
library(boot)
library(cluster)
library(MASS)
library(mgcv)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(ggplot2)
library(jcolors)
library(dplyr)
library(VennDiagram)
library(pheatmap)
library(forcats)
library(readxl)
library(apeglm)


#files to use

Lane_htseqcounts<- read.csv("Lasius_neglectus_basedCDSCounts.csv", header = T, row.names = 1)
Lane_htseqcounts<-as.matrix(Lane_htseqcounts,row.names=1,header=T)
colnames(Lane_htseqcounts)
head(Lane_htseqcounts,2) # display first two row in  the table plus header
dim(Lane_htseqcounts) # 14059 genes, 47 samples

#Metadatatable created:
#Lane_metadata_forWGCNA
meta_data_table_Lane <- read.csv("Lane_metadata.csv", header=TRUE, row.names=1)
summary(meta_data_table_Lane)
all(rownames(meta_data_table_Lane) %in% colnames(Lane_htseqcounts))





#############################looking at just queens
LaneQueens_metadata<-meta_data_table_Lane[meta_data_table_Lane$Caste=="Queen",]
Lane_Queens_counts<-Lane_htseqcounts[,LaneQueens_metadata$Sample]


dim(Lane_Queens_counts) # 14059 genes
# create DESeq2 object
deseq_Lane_Q_ALL <- DESeqDataSetFromMatrix(countData = Lane_Queens_counts,
                                           colData= LaneQueens_metadata, 
                                           design=~Colony+Age)
# filtering
deseq_Lane_Q_ALL_filtered <- deseq_Lane_Q_ALL[rowSums(counts(deseq_Lane_Q_ALL)>=10)>=8,]
summary(deseq_Lane_Q_ALL_filtered)
#10022 genes



# differential gene expression (just queens)
deseq_Lane_Q_ALL_filtered$Age<-relevel(deseq_Lane_Q_ALL_filtered$Age,ref="Young")

deseq_Lane_Q_ALL_filtered_LRT_age <- DESeq(deseq_Lane_Q_ALL_filtered,
                                           test="LRT", reduced=~Colony)
resultsNames(deseq_Lane_Q_ALL_filtered_LRT_age)
res_deseq_Lane_Q_ALL_filtered_LRT_age<-results(deseq_Lane_Q_ALL_filtered_LRT_age,                                                                                        name= "Age_Old_vs_Young")

BCV <- sqrt(mcols(deseq_Lane_Q_ALL_filtered_LRT_age)$dispersion)#dispersion
hist(BCV)#dispersion
summary(BCV)#dispersion

plotMA(res_deseq_Lane_Q_ALL_filtered_LRT_age)
#How many genes in queen fatbody change with age?
sum(res_deseq_Lane_Q_ALL_filtered_LRT_age$padj < 0.05, na.rm=TRUE)
#165 genes different between old and young L neglectus queens
res_deseq_Lane_Q_ALL_filtered_LRT_age_sig <- subset(res_deseq_Lane_Q_ALL_filtered_LRT_age, padj < 0.05)
#nrow(res_deseq_Lane_Q_ALL_filtered_LRT_age)
logsQ_Age <- res_deseq_Lane_Q_ALL_filtered_LRT_age_sig$log2FoldChange
pos<-res_deseq_Lane_Q_ALL_filtered_LRT_age_sig[res_deseq_Lane_Q_ALL_filtered_LRT_age_sig$log2FoldChange>0,]
dim(pos) # 85 up in old Queens compared to young queens
neg<-res_deseq_Lane_Q_ALL_filtered_LRT_age_sig[res_deseq_Lane_Q_ALL_filtered_LRT_age_sig$log2FoldChange<0,]
dim(neg) # 80 up in young queens compared to old Queens
plotMA(res_deseq_Lane_Q_ALL_filtered_LRT_age_sig)
nrow(res_deseq_Lane_Q_ALL_filtered_LRT_age_sig)

#################################Workers::

LaneWorkers_metadata<-meta_data_table_Lane[meta_data_table_Lane$Caste=="Worker",]
Lane_Workers_counts<-Lane_htseqcounts[,LaneWorkers_metadata$Sample]




# create DESeq2 object
deseq_Lane_W_ALL_Categories <- 
  DESeqDataSetFromMatrix(countData = Lane_Workers_counts,
                         colData= LaneWorkers_metadata, 
                         design=~Colony+Age+QueenPresence+Age:QueenPresence)
# filtering
deseq_Lane_W_ALL_Categories_filtered <- deseq_Lane_W_ALL_Categories[rowSums(counts
                                                                            (deseq_Lane_W_ALL_Categories)
                                                                            >=10)>=5,]
dim(deseq_Lane_W_ALL_Categories_filtered) #10691


#relevel for easier interpretation of deseq results
deseq_Lane_W_ALL_Categories_filtered$Age<-relevel(deseq_Lane_W_ALL_Categories_filtered$Age,ref="Young")
deseq_Lane_W_ALL_Categories_filtered$QueenPresence<-relevel(deseq_Lane_W_ALL_Categories_filtered$QueenPresence,ref="Queenless")


#LRT INTERACTION Queen presence and Position of worker
deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ <- DESeq(deseq_Lane_W_ALL_Categories_filtered,
                                                      test="LRT", reduced=~Colony+Age+QueenPresence)



###########BCV for revision
BCV <- sqrt(mcols(deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ)$dispersion)#dispersion
hist(BCV)#dispersion
summary(BCV)#dispersion
#####


resultsNames(deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ)
res_deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ<-
  results(deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ,name= "AgeOld.QueenPresenceQueenright")

library(apeglm)
winteraction<-lfcShrink(deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ,coef="AgeOld.QueenPresenceQueenright")
#write.csv(winteraction,file="deseq2-worker-interaction.csv")

plotMA(res_deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ)
sum(res_deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ$padj < 0.05, na.rm=TRUE)
# no interaction, #Also not with 0.1 padj, all genes fail BH FDR correction

##Reducing full model to remove interaction since it could not significantly explain the variation in any gene:
deseq_Lane_W_ALL_Categories1 <- 
  DESeqDataSetFromMatrix(countData = Lane_Workers_counts,
                         colData= LaneWorkers_metadata, 
                         design=~Colony+QueenPresence+Age)
# filtering
deseq_Lane_W_ALL_Categories1_filtered <- deseq_Lane_W_ALL_Categories1[rowSums(counts
                                                                              (deseq_Lane_W_ALL_Categories1)
                                                                              >=10)>=5,]
dim(deseq_Lane_W_ALL_Categories1_filtered)

deseq_Lane_W_ALL_Categories1_filtered$Age<-
  relevel(deseq_Lane_W_ALL_Categories1_filtered$Age,ref="Young")
deseq_Lane_W_ALL_Categories1_filtered$QueenPresence<-
  relevel(deseq_Lane_W_ALL_Categories1_filtered$QueenPresence,ref="Queenless")

#LRT EFFECT of Queen presence on gene expression
deseq_Lane_W_ALL_Categories1_filtered_LRT_Q <- DESeq(deseq_Lane_W_ALL_Categories1_filtered,
                                                     test="LRT", reduced=~Colony+Age)
resultsNames(deseq_Lane_W_ALL_Categories1_filtered_LRT_Q )
res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Q<-
  results(deseq_Lane_W_ALL_Categories1_filtered_LRT_Q,name= "QueenPresence_Queenright_vs_Queenless")
wqueen<-lfcShrink(deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ,coef="QueenPresence_Queenright_vs_Queenless")
#write.csv(wqueen,file="deseq2-worker-queenpresence.csv")
#plotMA(res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Q)
sum(res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Q$padj < 0.05, na.rm=TRUE)
#Effect of queen=0

# Using same additive model, now controlling for any variation explained by colony or queen 
#LRT
deseq_Lane_W_ALL_Categories1_filtered_LRT_Age <- DESeq(deseq_Lane_W_ALL_Categories1_filtered,
                                                       test="LRT", reduced=~Colony+QueenPresence)
resultsNames(deseq_Lane_W_ALL_Categories1_filtered_LRT_Age )
res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age<-
  results(deseq_Lane_W_ALL_Categories1_filtered_LRT_Age,name= "Age_Old_vs_Young")

wposition<-lfcShrink(deseq_Lane_W_ALL_Categories_filtered_LRT_AxQ,coef="Age_Old_vs_Young")
#write.csv(wposition,file="deseq2-worker-pposition.csv")
plotMA(res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age)
sum(res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age$padj < 0.05, na.rm=TRUE)
#2743
res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age_sig<-
  subset(res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age, padj < 0.05)
posW<-res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age_sig[res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age_sig$log2FoldChange>0,]
dim(posW) # 1229 up in old workers compared to young workers
negW<-res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age_sig[res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age_sig$log2FoldChange<0,]
dim(negW) # 1444 up in young workers compared to old workers



######## glmpca  (revision for jeb)

library(glmpca)
fitq <- glmpca(
  Lane_Queens_counts[rowSums(counts(deseq_Lane_Q_ALL)>=10)>=8,],#filtered read counts
  family = "negative.binomial",  # or "negative.binomial" for overdispersed data
  L = 2  # Number of principal components to extract
)

PCAPLOTq=data.frame(fitq$factors,Age=LaneQueens_metadata$Age)

fitw <- glmpca(
  Lane_Workers_counts[rowSums(counts(deseq_Lane_Q_ALL)>=10)>=5,]+1,#filtered read counts
  family = "negative.binomial",  # or "negative.binomial" for overdispersed data
  L = 2  # Number of principal components to extract
)
PCAPLOTw=data.frame(fitw$factors,Age=LaneWorkers_metadata$Age,QueenPresence=LaneWorkers_metadata$QueenPresence)

A<-ggplot(PCAPLOTq, aes(dim1, dim2, color = Age)) +
  geom_point(size = 5, shape = 19) +
  scale_color_manual(values = c("#969696", "#D73027")) +
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_classic()+
  theme(legend.position = "none")+
theme(axis.text = element_text(size = 14))+
  theme(aspect.ratio=1)


C<-ggplot(PCAPLOTw, aes(dim1, dim2, color=Age, shape=QueenPresence)) +
  geom_point(size = 5) +
  scale_color_manual(values = c("mediumseagreen","mediumpurple4")) +
  scale_shape_manual(values=c(17,15))+
  xlab("") +
  ylab("") +
  coord_fixed() +
  theme_classic()+
  theme(legend.position = "none")+
  theme(axis.text = element_text(size = 14))+
  


library(egg)
ggarrange(C, A, ncol=1)

#################



###########comparing posterior lfc: -> Suppl. Fig. 2
par(mfrow=c(3,1))
breaking=c(-5.5,-4.5,-2.5,-2.5,-1.5,-0.5,0.5,1.5,2.5,3.5,4.5,5.5)
hist(wposition$log2FoldChange,breaks=500,xlim=c(-6,6),main="worker position",xlab="")
hist(wqueen$log2FoldChange,xlim=c(-6,6),breaks=200,main="queen presence",xlab="")
hist(winteraction$log2FoldChange,xlim=c(-6,6),breaks=60,xlab="Posterior log2 fold change",main="interaction worker position : queen presence")

summary(abs(wqueen$log2FoldChange))
summary(abs(winteraction$log2FoldChange))
summary(abs(wposition$log2FoldChange))
wpositionn=na.omit(wposition)
min(abs(na.omit(wposition$log2FoldChange[res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age$padj<0.05])))


dat=data.frame(padj=res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age$padj,lfc=wposition$log2FoldChange)
pcol=rep("black",nrow(dat))
pcol[dat$padj<0.05]<-"red"
plot(dat$padj~dat$lfc,data=dat,col=pcol,pch=19)


dat=data.frame(padj=res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Q$padj,lfc=wqueen$log2FoldChange)
pcol=rep("black",nrow(dat))
pcol[dat$padj<0.05]<-"red"
plot(dat$padj~dat$lfc,data=dat,col=pcol,pch=19)

dat=data.frame(padj=res_deseq_Lane_W_ALL_Categories1_filtered_LRT_Age$padj,lfc=wposition$log2FoldChange)
pp=rep("black",nrow(dat))
pcol[dat$padj<0.05]<-"red"
plot(dat$padj~dat$lfc,data=dat,col=pcol,pch=19)

#all padj 1 for queen presence



##########VITELLOGENINS -> FIG 3

LnegWorkersVgConventional<- plotCounts(deseq_Lane_W_ALL_Categories1_filtered_LRT_Age, 
                                       gene="Lneg_g11491", normalized=T,
                                       transform=T,
                                       intgroup=c("Age"),
                                       returnData=TRUE)


LnegWorkersLneg_Vitellogening002<- plotCounts(deseq_Lane_W_ALL_Categories1_filtered_LRT_Age, 
                                              gene="Lneg_Vitellogening002", normalized=T,
                                              transform=T,
                                              intgroup=c("Age"),
                                              returnData=TRUE)

LnegWorkers_g03782<- plotCounts(deseq_Lane_W_ALL_Categories1_filtered_LRT_Age, 
                                gene="Lneg_g03782", normalized=T,
                                transform=T,
                                intgroup=c("Age"),
                                returnData=TRUE)
LnegWorkers_g15702<- plotCounts(deseq_Lane_W_ALL_Categories1_filtered_LRT_Age, 
                                gene="Lneg_g15702", normalized=T,
                                transform=T,
                                intgroup=c("Age"),
                                returnData=TRUE)

library(beeswarm)
par(mfrow=c(2,2))
COL=c("#5c478a","#3cb371")        
par(mar=c(3,5,3,3))
beeswarm(count/1000~Age,LnegWorkersLneg_Vitellogening002
         ,las=1,ylab="",main="A) Vg-like A (Lneg_Vitellogening002)",pch=19,col=COL,xlab="",xaxt="n",spacing=0.3,method="hex")
mtext(side=2,text="K reads",line=4)

beeswarm(count/1000~Age,LnegWorkers_g15702
         ,las=1,ylab="",main="B) Vg-like A (Lneg_g15702)",pch=19,col=COL,xlab="",xaxt="n",spacing=0.3,method="hex")

beeswarm(count/1000~Age,LnegWorkers_g03782
         ,las=1,ylab="",main="C) Vg_like C (Lneg_g03782)",pch=19,method="hex",col=COL,xlab="",xaxt="n",spacing=0.3)
mtext(side=1,at=1:2,text=c("Inside","Outside"),line=1.5)
mtext(side=2,text="K reads",line=4)

beeswarm(count/1000~Age,LnegWorkersVgConventional
         ,las=1,ylab="",main="D) Conventional Vg (Lneg_g11491)",pch=19,col=COL,xlab="",xaxt="n",spacing=0.3,method="hex")
mtext(side=1,at=1:2,text=c("Inside","Outside"),line=1.5)










##################FOR REVISION JEB: BCV PER WORKER CATEGORY
mIN<-meta_data_table_Lane[meta_data_table_Lane$Caste=="Worker"&meta_data_table_Lane$Age=="Young",]
IN<-Lane_htseqcounts[,mIN$Sample]
dIN <-   DESeqDataSetFromMatrix(countData = IN,
                         colData= mIN, 
                         design=~QueenPresence)
fdIN <-dIN[rowSums(counts(dIN)>=10)>=5]# filtering
dim(fdIN) 
dsIN <- DESeq(fdIN)
BCV <- sqrt(mcols(dsIN)$dispersion)#dispersion
hist(BCV)#dispersion
summary(BCV)#dispersion

mOUT<-meta_data_table_Lane[meta_data_table_Lane$Caste=="Worker"&meta_data_table_Lane$Age=="Old",]
OUT<-Lane_htseqcounts[,mOUT$Sample]
dOUT <-   DESeqDataSetFromMatrix(countData = OUT,
                                colData= mOUT, 
                                design=~QueenPresence)
fdOUT <-dOUT[rowSums(counts(dOUT)>=10)>=5]# filtering
dim(fdOUT) 
dsOUT <- DESeq(fdOUT)
BCV <- sqrt(mcols(dsOUT)$dispersion)#dispersion
hist(BCV)#dispersion
summary(BCV)#dispersion



mQ<-meta_data_table_Lane[meta_data_table_Lane$Caste=="Worker"&meta_data_table_Lane$QueenPresence=="Queenright",]
Q<-Lane_htseqcounts[,mQ$Sample]
dQ <-   DESeqDataSetFromMatrix(countData = Q,
                                 colData= mQ, 
                                 design=~Age)
fdQ <-dQ[rowSums(counts(dQ)>=10)>=5]# filtering
dsQ <- DESeq(fdQ)
BCV <- sqrt(mcols(dsQ)$dispersion)#
hist(BCV)#
summary(BCV)#

mNQ<-meta_data_table_Lane[meta_data_table_Lane$Caste=="Worker"&meta_data_table_Lane$QueenPresence=="Queenless",]
NQ<-Lane_htseqcounts[,mNQ$Sample]
dNQ <-   DESeqDataSetFromMatrix(countData = NQ,
                               colData= mNQ, 
                               design=~Age)
fdNQ <-dNQ[rowSums(counts(dNQ)>=10)>=5]# filtering
dsNQ <- DESeq(fdNQ)
BCV <- sqrt(mcols(dsNQ)$dispersion)#
hist(BCV)#
summary(BCV)#
