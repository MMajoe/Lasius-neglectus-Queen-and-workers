#Supplementary TopGO and BH correction 

## Data tables associated are attached:
#LasiusNeglectus_proteome_reference_Unique8255lines.txt
#Up_OldQueens_proteins.txt
#Up_YoungQueens_proteins.txt
#FINAL_Up_OldWorkers_ALLproteins_Sorted.txt
#FINAL_UpyoungWorkers_ALLproteins_Sorted.txt


#Final TopGO Fisher test result values- applied BH correction
#
#Generating sheets in Supplementary table 6. GO terms and adjusted p-values
#Barplots for supplementary fig 3.
#

#Lasius neglectus: TopGo:
#AGE
library(ggplot2)
library(ggwordcloud)
library(enrichplot)
#BiocManager::install('enrichplot',force=TRUE)
#install.packages(c("devtools", "data.table", "enrichR"))

# From Bioconductor
# if(!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
#BiocManager::install(c("org.Mm.eg.db", "clusterProfiler", "enrichplot"))

# From GitHub
# topGO with enrichment_barplot()
#devtools::install_github("ycl6/topGO-feat", ref = "v2.41.0-barplot",force=TRUE)
library(topGO)
library(dplyr)
#you can order.by='Score' as well as 'Ratio'. Which is just gene counts or ratio 
#
#Adapted from following sources: https://datacatz.wordpress.com/2018/01/19/gene-set-enrichment-analysis-with-topgo-part-1/
#http://avrilomics.blogspot.com/2015/07/using-topgo-to-test-for-go-term.html
#colurs: https://ggplot2.tidyverse.org/reference/scale_viridis.html
#https://www.researchgate.net/figure/Bar-chart-of-enrichment-ratios-for-GO-and-KEGG-categories-in-the-three-gene-lists-54_fig1_51096836

#UNIVERSE
Lneguniverse<-readMappings("LasiusNeglectus_proteome_reference_Unique8255lines.txt")
head(Lneguniverse)
LnegUniverse_genelist <- names(Lneguniverse)
summary(LnegUniverse_genelist)

#Old queens-Age
OldQueen_AgeUp<-"Up_OldQueens_proteins.txt"
listOldQueen_AgeUp<-readLines(OldQueen_AgeUp)
head(listOldQueen_AgeUp)

gene_listOldQueen_AgeUp <- factor(as.integer(LnegUniverse_genelist %in% listOldQueen_AgeUp))
head(gene_listOldQueen_AgeUp)


names(gene_listOldQueen_AgeUp) <- LnegUniverse_genelist
GO_data_OldQueen_Up<- new("topGOdata", description="GO_Brain", ontology="BP", allGenes=gene_listOldQueen_AgeUp,  annot = annFUN.gene2GO, 
                          gene2GO=Lneguniverse)
result_topGO_wo1OldQueen_Up<- runTest(GO_data_OldQueen_Up, algorithm = "weight01", statistic = "fisher")

result_topGO_wo1OldQueen_Up

#create a table containing all necessary information and save to file
result_table_w01_Age_wo1OldQueen_Up<- GenTable(GO_data_OldQueen_Up, Fisher = result_topGO_wo1OldQueen_Up, orderBy = "Fisher", 
                                               ranksOf = "Fisher", topNodes = numSigGenes(GO_data_OldQueen_Up),numChar=1000) 
showSigOfNodes(GO_data_OldQueen_Up, score(result_topGO_wo1OldQueen_Up), firstSigNodes = 5, useInfo ='all',useFullNames = TRUE,
               .NO.CHAR=20)

GOterms_OldQueen<-result_table_w01_Age_wo1OldQueen_Up$GO.ID[as.numeric(result_table_w01_Age_wo1OldQueen_Up$Fisher)<0.05]

#What genes are related to the GO terms:

GOterms<-result_table_w01_Age_wo1OldQueen_Up$GO.ID[as.numeric(result_table_w01_Age_wo1OldQueen_Up$Fisher)<0.05]
# get genes corresponding to these terms
genes = genesInTerm(GO_data_OldQueen_Up, GOterms)
GO_genes=c()
GO_genes[1]<-c("GOterm \t Genes")
# iterate over GO terms and write into table with tab-sep between GO and Genes and a Comma between the different Genes (if multiple)
for (j in 1:length(GOterms))
{
  GOterm<-GOterms[j]
  # extract genes from specific GO term and write into format to be written into table
  genesInGO<-genes[GOterm][[1]]
  genesInGO<-paste(genesInGO, collapse=',')
  GO_genes[j+1]=paste(GOterm,"\t",genesInGO, sep="")
}

#View(GO_genes)
#write.table(GO_genes,"LnegOldQueens_TopGOGenes.csv",sep="\t",col.names = T)


####
#Young Queens-Age
YoungQueen_AgeUp<-"Up_YoungQueens_proteins.txt"
listYoungQueen_AgeUp<-readLines(YoungQueen_AgeUp)
head(listYoungQueen_AgeUp)

gene_listYoungQueen_AgeUp <- factor(as.integer(LnegUniverse_genelist %in% listYoungQueen_AgeUp))
head(gene_listYoungQueen_AgeUp)


names(gene_listYoungQueen_AgeUp) <- LnegUniverse_genelist
GO_data_YoungQueen_Up<- new("topGOdata", description="GO_Brain", ontology="BP", allGenes=gene_listYoungQueen_AgeUp,  annot = annFUN.gene2GO, 
                            gene2GO=Lneguniverse)
result_topGO_wo1YoungQueen_Up<- runTest(GO_data_YoungQueen_Up, algorithm = "weight01", statistic = "fisher")
result_topGO_wo1YoungQueen_Up

#create a table containing all necessary information and save to file
result_table_w01_Age_wo1YoungQueen_Up<- GenTable(GO_data_YoungQueen_Up, Fisher = result_topGO_wo1YoungQueen_Up, orderBy = "Fisher", 
                                                 ranksOf = "Fisher", topNodes = numSigGenes(GO_data_YoungQueen_Up),numChar=1000) 
showSigOfNodes(GO_data_YoungQueen_Up, score(result_topGO_wo1YoungQueen_Up), firstSigNodes = 5, useInfo ='all',useFullNames = TRUE,
               .NO.CHAR=20)

#What genes?
GOterms<-result_table_w01_Age_wo1YoungQueen_Up$GO.ID[as.numeric(result_table_w01_Age_wo1YoungQueen_Up$Fisher)<0.05]
# get genes corresponding to these terms
genes = genesInTerm(GO_data_YoungQueen_Up, GOterms)
GO_genes=c()
GO_genes[1]<-c("GOterm \t Genes")
# iterate over GO terms and write into table with tab-sep between GO and Genes and a Comma between the different Genes (if multiple)
for (j in 1:length(GOterms))
{
  GOterm<-GOterms[j]
  # extract genes from specific GO term and write into format to be written into table
  genesInGO<-genes[GOterm][[1]]
  genesInGO<-paste(genesInGO, collapse=',')
  GO_genes[j+1]=paste(GOterm,"\t",genesInGO, sep="")
}

# View(GO_genes)
# write.table(GO_genes,"LnegYoungQueens_TopGOGenes.csv",sep="\t",col.names = T)

###Workers
OldWorkers_AgeUp<-"FINAL_Up_OldWorkers_ALLproteins_Sorted.txt"
listOldWorkers_AgeUp<-readLines(OldWorkers_AgeUp)
head(listOldWorkers_AgeUp)

gene_listOldWorkers_AgeUp <- factor(as.integer(LnegUniverse_genelist %in% listOldWorkers_AgeUp))
head(gene_listOldWorkers_AgeUp)

names(gene_listOldWorkers_AgeUp) <- LnegUniverse_genelist
GO_data_OldWorkers_Up<- new("topGOdata", description="GO_Brain", ontology="BP", allGenes=gene_listOldWorkers_AgeUp,  annot = annFUN.gene2GO, 
                            gene2GO=Lneguniverse)
result_topGO_wo1OldWorkers_Up<- runTest(GO_data_OldWorkers_Up, algorithm = "weight01", statistic = "fisher")
result_topGO_wo1OldWorkers_Up

#create a table containing all necessary information and save to file
result_table_w01_Age_wo1OldWorkers_Up<- GenTable(GO_data_OldWorkers_Up, Fisher = result_topGO_wo1OldWorkers_Up, orderBy = "Fisher", 
                                                 ranksOf = "Fisher", topNodes = numSigGenes(GO_data_OldWorkers_Up),numChar=1000) 
showSigOfNodes(GO_data_OldWorkers_Up, score(result_topGO_wo1OldWorkers_Up), firstSigNodes = 5, useInfo ='all',useFullNames = TRUE,
               .NO.CHAR=20)


#What genes?
GOterms<-result_table_w01_Age_wo1OldWorkers_Up$GO.ID[as.numeric(result_table_w01_Age_wo1OldWorkers_Up$Fisher)<0.05]
# get genes corresponding to these terms
genes = genesInTerm(GO_data_OldWorkers_Up, GOterms)
GO_genes=c()
GO_genes[1]<-c("GOterm \t Genes")
# iterate over GO terms and write into table with tab-sep between GO and Genes and a Comma between the different Genes (if multiple)
for (j in 1:length(GOterms))
{
  GOterm<-GOterms[j]
  # extract genes from specific GO term and write into format to be written into table
  genesInGO<-genes[GOterm][[1]]
  genesInGO<-paste(genesInGO, collapse=',')
  GO_genes[j+1]=paste(GOterm,"\t",genesInGO, sep="")
}

# View(GO_genes)
# write.table(GO_genes,"LnegOutsideWorkers_TopGOGenes.csv",sep="\t",col.names = T)


####
#Up in Young Workers:
YoungWorkers_AgeUp<-"FINAL_UpyoungWorkers_ALLproteins_Sorted.txt"
listYoungWorkers_AgeUp<-readLines(YoungWorkers_AgeUp)
head(listYoungWorkers_AgeUp)

gene_listYoungWorkers_AgeUp <- factor(as.integer(LnegUniverse_genelist %in% listYoungWorkers_AgeUp))
head(gene_listYoungWorkers_AgeUp)


names(gene_listYoungWorkers_AgeUp) <- LnegUniverse_genelist
GO_data_YoungWorkers_Up<- new("topGOdata", description="GO_Brain", ontology="BP", allGenes=gene_listYoungWorkers_AgeUp,  annot = annFUN.gene2GO, 
                              gene2GO=Lneguniverse)
result_topGO_wo1YoungWorkers_Up<- runTest(GO_data_YoungWorkers_Up, algorithm = "weight01", statistic = "fisher")
result_topGO_wo1YoungWorkers_Up

#create a table containing all necessary information and save to file
result_table_w01_Age_wo1YoungWorkers_Up<- GenTable(GO_data_YoungWorkers_Up, Fisher = result_topGO_wo1YoungWorkers_Up, orderBy = "Fisher", 
                                                   ranksOf = "Fisher", topNodes = numSigGenes(GO_data_YoungWorkers_Up),numChar=1000) 
showSigOfNodes(GO_data_YoungWorkers_Up, score(result_topGO_wo1YoungWorkers_Up), firstSigNodes = 5, useInfo ='all',useFullNames = TRUE,
               .NO.CHAR=20)

#What genes?
GOterms<-result_table_w01_Age_wo1YoungWorkers_Up$GO.ID[as.numeric(result_table_w01_Age_wo1YoungWorkers_Up$Fisher)<0.05]
# get genes corresponding to these terms
genes = genesInTerm(GO_data_YoungWorkers_Up, GOterms)
GO_genes=c()
GO_genes[1]<-c("GOterm \t Genes")
# iterate over GO terms and write into table with tab-sep between GO and Genes and a Comma between the different Genes (if multiple)
for (j in 1:length(GOterms))
{
  GOterm<-GOterms[j]
  # extract genes from specific GO term and write into format to be written into table
  genesInGO<-genes[GOterm][[1]]
  genesInGO<-paste(genesInGO, collapse=',')
  GO_genes[j+1]=paste(GOterm,"\t",genesInGO, sep="")
}

# View(GO_genes)
# write.table(GO_genes,"LnegIntsideWorkers_TopGOGenes.csv",sep="\t",col.names = T)


library(openxlsx)
dataset_list_topGO <- list('InsideWorkers'= result_table_w01_Age_wo1YoungWorkers_Up ,  
                           'OutsideWorkers' = result_table_w01_Age_wo1OldWorkers_Up, 
                           'YoungQueens' = result_table_w01_Age_wo1YoungQueen_Up,
                           'OldQueens'= result_table_w01_Age_wo1OldQueen_Up)
#write.xlsx(dataset_list_topGO, file = 'All_Results_TopGO_forpadjust.xlsx')

####padjust of fisher'test for weight01'
library(readxl)
WorkersIn<-read_excel('All_Results_TopGO_forpadjust.xlsx',sheet='InsideWorkers')
WorkersOut<-read_excel('All_Results_TopGO_forpadjust.xlsx',sheet='OutsideWorkers')
QueensY<-read_excel('All_Results_TopGO_forpadjust.xlsx',sheet='YoungQueens')
QueensO<-read_excel('All_Results_TopGO_forpadjust.xlsx',sheet='OldQueens')

WorkerInpadj<-p.adjust(WorkersIn$Fisher,method = 'BH',n=length(WorkersIn$Fisher))
WorkersIn$WorkersInpadj<-WorkerInpadj

WorkerOutpadj<-p.adjust(WorkersOut$Fisher,method = 'BH',n=length(WorkersOut$Fisher))
WorkersOut$WorkersOutpaj<-WorkerOutpadj

QueensYpadj<-p.adjust(QueensY$Fisher,method='BH',n=length(QueensY$Fisher))
QueensY$QueensYpadj<-QueensYpadj #all more than cut off of adj p >0.05

QueensOpadj<-p.adjust(QueensO$Fisher,method='BH',n=length(QueensO$Fisher))
QueensO$QueensOpadj<-QueensOpadj

library(openxlsx)
NewTopGOResults_BHadjust<- list('PadjworkersIn' = WorkersIn, 'PadjworkersOut' = WorkersOut, 'PadjQueensYoung'= QueensY , 'PadjQueensOld'= QueensO )

#above used as sheet names in Suppl. table 6

Theme<-theme(plot.title = element_text(face="bold" ,size=(24),colour="black"),
             axis.text=element_text(size=(22),colour="black"),
             axis.title=element_text(face="bold" ,size=(24),colour="black"))


PadjInW<-WorkersIn

PadjOutW<-WorkersOut

PadjYoungQueen<-QueensY

PadjOldQueen<-QueensO

InsideWorkersPlot_Padj<-ggplot(PadjInW[1:5,], aes(x=reorder(Term,-WorkersInpadj),y=-log10(WorkersInpadj)))+
  geom_bar(stat="identity",fill="mediumpurple4")+theme_classic()+coord_flip()+scale_fill_discrete()+
  geom_text(aes(label= Significant),
            #position=position_dodge(width=0.9),
            vjust=0.25,hjust=-0.25,size=(12))+
  Theme+labs(y = "-log10(adjusted Fisher's p-value)", x = "Enriched GO terms")
#+
# Theme+labs(y = "No. of genes", x = "GO terms")

InsideWorkersPlot_Padj#+Theme+labs(y = "-log10(adjusted Fisher's p-value)", x = "GO terms")
#+geom_text(aes(label= as.ratio(PadjInW[1:5,]$Significant/PadjInW[1:5,]$Annotated)), 
#                                  position=position_dodge(width=0.9), vjust=-0.25)
##
#ggsave("InsideWorkersPlot_Padj.pdf", InsideWorkersPlot_Padj, width=22, height=8, dpi=300)

par(mar = c(0.5,0.5,0.5,1))
OutsideWorkersPlot_Padj<-ggplot(PadjOutW[1:5,], aes(x=reorder(Term,-WorkersOutpaj),y=-log10(WorkersOutpaj)))+
  geom_bar(stat="identity",fill="mediumseagreen")+theme_classic()+coord_flip()+scale_fill_discrete()+
  geom_text(aes(label= Significant),#/Annotated),
            #position=position_dodge(width=0.9),
            vjust=0.25,hjust=-0.25,size=(12))+
  Theme+labs(y = "-log10(adjusted Fisher's p-value)", x = "Enriched GO terms")
#+
# Theme+labs(y = "No. of genes", x = "GO terms")
#par(mar = c(0.5,0.5,0.5,1))

OutsideWorkersPlot_Padj
#ggsave("OutsideWorkersPlot_Padj.pdf", OutsideWorkersPlot_Padj, width=22, height=8, dpi=300)


# OutsideWorkersPlot_Padj1<-ggplot(OutW5, aes(x=Term,y=Significant))+
#   geom_bar(stat="identity",fill="mediumseagreen")+theme_bw()+coord_flip()+scale_fill_discrete()+Theme
#   # geom_text(aes(label= round(WorkersOutpaj,3),
#   #               #position=position_dodge(width=0.9), 
#   #               vjust=-0.25,size=(20)))+
#   # Theme+labs(y = "No. of genes", x = "GO terms")
# 
# OutsideWorkersPlot_Padj1+labs(y = "No. of genes", x = "GO terms")

#
PadjYoungQueenPlot_Padj<-ggplot(PadjYoungQueen[1:5,], aes(x=reorder(Term,-QueensYpadj),y=-log10(QueensYpadj)))+
  geom_bar(stat="identity",fill='firebrick')+theme_classic()+coord_flip()+scale_fill_discrete()+
  geom_text(aes(label= Significant),
            #position=position_dodge(width=0.9),
            vjust=0.25,hjust=-0.25,size=(12))+
  Theme+labs(y = "-log10(adjusted Fisher's p-value)", x = "Enriched GO terms")

PadjYoungQueenPlot_Padj
#ggsave("PadjYoungQueenPlot_Padj.pdf", PadjYoungQueenPlot_Padj, width=22, height=8, dpi=300)
###

PadjOldQueenPlot_Padj<-ggplot(PadjOldQueen[1:5,], aes(x=reorder(Term,-QueensOpadj),y=-log10(QueensOpadj)))+
  geom_bar(stat="identity",fill='azure4')+theme_classic()+coord_flip()+scale_fill_discrete()+
  geom_text(aes(label= Significant),
            #position=position_dodge(width=0.9),
            vjust=0.25,hjust=-0.25,size=(12))+
  Theme+labs(y = "-log10(adjusted Fisher's p-value)", x = "Enriched GO terms")

PadjOldQueenPlot_Padj#+Theme+labs(y = "-log10(adjusted Fisher's p-value)", x = "GO terms")
#ggsave("PadjOldQueenPlot_Padj.pdf", PadjOldQueenPlot_Padj, width=22, height=8, dpi=300)

