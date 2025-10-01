#'Lasius neglectus Updated analysis final R script:
#'
###RESULTING TABLE- neglectus-SupplementaryTable7_Lneg_TIJ_OxStress_Vgs_July.xlsx
#
#
#Used orthofinder to find proteins of L. neglectus and associated D. mellanogaster orthogroups (can be referenced in Supp Table 7 as sheet = "TIJOGsFromFullproteomes". 
#Based on available List of TIJ genes and Ox. stress genes- extracted list to create the two following files
# 
#LnegTIJOGs.xlsx
##LnegOxStressOgs.xlsx
#
#Then compared our DEGs against the above list and found what other genes they BLAST against using the following two files
#neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx
#neglectus-SupplementaryTable4_Allproteins_Blastbesthits.xlsx
#


#
library(readxl)
#Inside workers
LnegTIJInWorkers<-read_excel("LnegTIJOGs.xlsx", sheet='InsideWorkers_DEG')
LnegDEGInWorkers<-read_excel("neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx",sheet='InsideWorkers_DEGs')
InWorkers<-merge(LnegDEGInWorkers,LnegTIJInWorkers,by.x="Gene Name",sort=T,no.dups=T)
#Outside
LnegTIJOutWorkers<-read_excel("LnegTIJOGs.xlsx", sheet='OutsideWorkers_DEG')
LnegDEGOutWorkers<-read_excel("neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx",sheet='OutsideWorkers_DEGs')
OutWorkers<-merge(LnegDEGOutWorkers,LnegTIJOutWorkers,by.x="Gene Name",sort=T,no.dups=T)
OutWorkers[,c(1,9)]
#YoungQueens
LnegTIJYoungQueens<-read_excel("LnegTIJOGs.xlsx", sheet='YoungQueens_DEG')
LnegDEGYoungQueens<-read_excel("neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx",sheet='YoungQueens_DEGs')
YoungQueens<-merge(LnegDEGYoungQueens,LnegTIJYoungQueens,by.x="Gene Name",sort=T,no.dups=T)
#OldQueens
LnegTIJOldQueens<-read_excel("LnegTIJOGs.xlsx", sheet='OldQueens_DEG')
LnegDEGOldQueens<-read_excel("neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx",sheet='OldQueens_DEGs')
OldQueens<-merge(LnegDEGOldQueens,LnegTIJOldQueens,by.x="Gene Name",sort=T,no.dups=T)

# library(openxlsx)
# dataset_names_genes <- list('InsideWorkers_TIJ' = InWorkers, 
#                             'OutsideWorkers_TIJ' = OutWorkers,
#                             'YoungQueens_TIJ'=YoungQueens,
#                             'OldQueens_TIJ'= OldQueens)
#write.xlsx(dataset_names_genes, file = 'Supplementary_Lneg_TIJ_OG.xlsx')

TIJIW<-chisq.test(matrix(c(7,(67-7),1444,(14059-1444)),nrow=2))

TIJOW<-chisq.test(matrix(c(11,(67-11),1229,(14059-1229)),nrow=2))

TIJYQ<-chisq.test(matrix(c(1,(67-1),80,(14059-80)),nrow=2))

TIJOQ<-chisq.test(matrix(c(2,(67-2),85,(14059-85)),nrow=2))

#Oxidative stress
#Inside workers
LnegOxstressInWorkers<-read_excel("LnegOxStressOgs.xlsx", sheet='InsideWorkers_Oxstress')
LnegDEGInWorkers<-read_excel("neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx",sheet='InsideWorkers_DEGs')
InWorkers<-merge(LnegDEGInWorkers,LnegOxstressInWorkers,by.x="Gene Name",sort=T,no.dups=T)
#Outside
LnegOxstressOutWorkers<-read_excel("LnegOxStressOgs.xlsx", sheet='OutsideWorkers_Oxstress')
LnegDEGOutWorkers<-read_excel("neglectus-SupplementaryTable6_DEGs_PadjTopGO.xlsx",sheet='OutsideWorkers_DEGs')
OutWorkers<-merge(LnegDEGOutWorkers,LnegOxstressOutWorkers,by.x="Gene Name",sort=T,no.dups=T)
OutWorkers[,c(1,9)]

OxStressIW<-chisq.test(matrix(c(3,(31-3),1444,(14059-1444)),nrow=2))

OxStressOW<-chisq.test(matrix(c(7,(31-7),1229,(14059-1229)),nrow=2))


Lneg_allProteinBLASTp<-as.data.frame(read_excel("neglectus-SupplementaryTable4_Allproteins_Blastbesthits.xlsx",
  sheet=2))
colnames(Lneg_allProteinBLASTp)
Lneg_BlastP_ALL<-Lneg_allProteinBLASTp[,c(1,4:5)]
head(Lneg_BlastP_ALL)

OutWorkersOxStress_BLASTp<-merge(OutWorkers,Lneg_BlastP_ALL,by.x="Lneg_protein",sort=T,no.dups=T)
InWorkersOxStress_BLASTp<-merge(InWorkers,Lneg_BlastP_ALL,by.x="Lneg_protein",sort=T,no.dups=T)


#Make 
# library(openxlsx)
# dataset_namesOXstress_genes <- list('InsideWorkers_OxStress' = InWorkers, 
#                             'OutsideWorkers_OxStress' = OutWorkers,
#                             'YoungQueens_OxStress'=YoungQueens,
#                             'OldQueens_OxStress'= OldQueens)
#write.xlsx(dataset_names_genes, file = 'Supplementary_Lneg_TIJ_OG.xlsx')


##TIJ-

colnames(Lneg_allProteinBLASTp)
Lneg_BlastP_ALL_TIJ<-Lneg_allProteinBLASTp[,c(2,4:5)]
colnames(Lneg_BlastP_ALL_TIJ)[1] <- 'Lneg_protein'

InworkersTIJ<-as.data.frame(read_excel('Supplementary_Lneg_TIJ_OG.xlsx',sheet="InsideWorkers_TIJ"))
OutworkersTIJ<-as.data.frame(read_excel('Supplementary_Lneg_TIJ_OG.xlsx',sheet="OutsideWorkers_TIJ"))
YoungQueensTIJ<-as.data.frame(read_excel('Supplementary_Lneg_TIJ_OG.xlsx',sheet="YoungQueens_TIJ"))
OldQueensTIJ<-as.data.frame(read_excel('Supplementary_Lneg_TIJ_OG.xlsx',sheet="OldQueens_TIJ"))

head(Lneg_BlastP_ALL_TIJ$Lneg_Protein)
InworkersTIJ_BLASTp<-merge(InworkersTIJ,Lneg_BlastP_ALL_TIJ,by.x="Lneg_protein",sort=T,no.dups=T)
OutworkersTIJ_BLASTp<-merge(OutworkersTIJ,Lneg_BlastP_ALL_TIJ,by.x="Lneg_protein",sort=T,no.dups=T)
YoungQueensTIJ_BLASTp<-merge(YoungQueensTIJ,Lneg_BlastP_ALL_TIJ,by.x="Lneg_protein",sort=T,no.dups=T)
OldQueensTIJ_BLASTp<-merge(OldQueensTIJ,Lneg_BlastP_ALL_TIJ,by.x="Lneg_protein",sort=T,no.dups=T)

#'OutWorkersOxStress_BLASTp
#InWorkersOxStress_BLASTp
TIJandOXstress_genes <- list('InsideWorkers_TIJ' = InworkersTIJ_BLASTp, 
                            'OutsideWorkers_TIJ' = OutworkersTIJ_BLASTp,
                            'YoungQueens_TIJ'=YoungQueensTIJ_BLASTp,
                            'OldQueens_TIJ'= OldQueensTIJ_BLASTp,
                            'InsideWorkers_OxStress' = InWorkersOxStress_BLASTp,
                            'OutsideWorkers_OxStress' = OutWorkersOxStress_BLASTp)
#write.xlsx(TIJandOXstress_genes, file = 'Supplementary_Lneg_TIJ_OxStress_OG.xlsx')
