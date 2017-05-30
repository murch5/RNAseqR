library(dplyr)
library(tidyr)

expMatrix <- read.table("Monocyte.txt",sep="\t")

expMatrixMonocytes <- expMatrix[,grep("Mono", colnames(expMatrix) )]

filterList <- read.table("ToFilter.txt")

t <- which(colnames(expMatrixMonocytes) %in% filterList)

expMatrixMonocytes <- expMatrixMonocytes[,which(!(colnames(expMatrixMonocytes) %in% filterList[,1]))]

#expMatrixClassical <- expMatrixMonocytes[,grep("nonclassical{0}", colnames(expMatrixMonocytes))]

expMatrixMonocytesNoInter <-  expMatrixMonocytes[,grep("classical", colnames(expMatrixMonocytes))]
temp_index_nonclassical <- grep("nonclassical",colnames(expMatrixMonocytesNoInter))

expMatrixNonClassical <- expMatrixMonocytesNoInter[temp_index_nonclassical]

expMatrixClassical <- expMatrixMonocytesNoInter[-temp_index_nonclassical]

classicalID <- as.data.frame(colnames(expMatrixClassical))

classicalID <- classicalID %>%
              mutate(group="Classical")  #%>%
             # mutate(group_color="red")

colnames(classicalID) <- c("Sample","Group")

nonclassicalID <- as.data.frame(colnames(expMatrixNonClassical))

nonclassicalID <- nonclassicalID %>%
            mutate(group="Non-classical")  #%>%
           # mutate(group_color="blue")


colnames(nonclassicalID) <- c("Sample","Group")

IDlist <- bind_rows(classicalID,nonclassicalID)

expComparison <- cbind(expMatrixClassical,expMatrixNonClassical)

expClassicalSum <- rowMeans(expMatrixClassical)

expNonclassicalSum <- rowMeans(expMatrixNonClassical)

expDiff <- expClassicalSum - expNonclassicalSum

exp <- as.data.frame(expDiff)

exp <- tibble::rownames_to_column(exp)

expClassical <- exp[which(exp[,2]>300),]
expNonClassical <- exp[which(exp[,2]<(-300)),]

expClassicalNames <- expClassical[,1]
expNonClassicalNames <- expNonClassical[,1]


expFilteredSet <- exp %>%
                filter(expDiff > 500 | expDiff < -500)

expMatrixFiltered <- expComparison[which(rownames(expMatrix) %in% expFilteredSet[,1]),]

print(colnames(expMatrixFiltered))

targetsNur77HomoSapiens <- read.table("Nur77Targets.txt",sep=",",stringsAsFactors=FALSE)
targetsNur77MusMusculus <- read.table("Nur77GenesMusMusculus.txt",sep=",",stringsAsFactors=FALSE)

targetsGenes <- targetsNur77HomoSapiens[-1,15]
targetsGeneMouse <- targetsNur77MusMusculus[-1,15]

expNur77Target <- expComparison[which(rownames(expComparison) %in% targetsGenes),]

targetsComparison <- cbind(targetsNur77HomoSapiens[,15],targetsNur77MusMusculus[,15])

#targetsComparison <- apply(targetsComparison, c(1,2), function(y)
#  {
#  return(toupper(y))
#})

targetsComparison[1,1] <- "Homo Sapiens"
targetsComparison[1,2] <- "Mus Musculus"

overlapGenes <- read.table("overlap.csv",sep=",",quote="",header=TRUE)

expOverlap <- expComparison[which(rownames(expMatrix) %in% overlapGenes[,1]),]


expOverlapSet <- exp %>%
  filter(expDiff > 10 | expDiff < -10)

expOverlapSetFiltered <- expOverlap[which(rownames(expOverlap) %in% expOverlapSet[,1]),]

mouseRNAseq <- read.table("GSE80036_STAR_ensGenev73.txt",sep="\t",stringsAsFactors=FALSE)

mouseRNAseq <- mouseRNAseq[-1,]
colnames(mouseRNAseq) <- mouseRNAseq[1,]
mouseRNAseqData <- mouseRNAseq[-1,c(1,9:13)]

mouseRNAseqData[,c(2:6)] <- sapply(mouseRNAseqData[,c(2:6)],as.numeric)


colnames(mouseRNAseqData) <- c("EnsemblID","Ly6CNeg_1","Ly6Neg_2","Ly6CPos_1","Ly6CPos_2","Ly6CPos_3")

con<-file("mouseRNAseq.csv",encoding="UTF-8")
write.table(mouseRNAseqData,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)

mouseRNAseqGeneName <- read.table("RNAfixedGeneName.csv",sep=",",stringsAsFactors=FALSE, header=TRUE)


classicalMus <- mouseRNAseqGeneName[,c(1,2:3)]
nonClassicalMus <- mouseRNAseqGeneName[,c(1,4:6)]


classicalMusSum <- rowMeans(classicalMus[,-1])

nonClassicalMusSum <- rowMeans(nonClassicalMus[,-1])

expDiffMus <- classicalMusSum - nonClassicalMusSum

expDiffMus <- as.data.frame(expDiffMus)

expDiffMus <- cbind(mouseRNAseqData[,1],expDiffMus)

rownames(expDiffMus) <- expDiffMus[,1]

expDiffMusFilteredSet <- expDiffMus %>%
  filter(expDiffMus > 100000 | expDiffMus < -100000)

expMusSetFiltered <- mouseRNAseqData[which(rownames(expDiffMus) %in% expDiffMusFilteredSet[,1]),]


con<-file("expMusToFixNames.csv",encoding="UTF-8")
write.table(expMusSetFiltered,file=con,row.names = TRUE,sep=",",quote=FALSE)


expNur77MusTarget <- mouseRNAseqGeneName[which(rownames(expDiffMus) %in% targetsGeneMouse),]


con<-file("expMusMusculusFilt100000.csv",encoding="UTF-8")
write.table(expMusSetFiltered,file=con,row.names = FALSE,sep=",",quote=FALSE)

con<-file("expMatrix.csv",encoding="UTF-8")
write.table(expMatrixFiltered,file=con,row.names = TRUE,sep=",",quote=FALSE)


con<-file("MonocyteComparison.csv",encoding="UTF-8")
write.csv(expFilteredSet,file=con,row.names = TRUE)

con<-file("ClassicalSet.csv",encoding="UTF-8")
write.csv(expClassicalNames,file=con,row.names = FALSE,quote=FALSE)

con<-file("NonclassicalSet.csv",encoding="UTF-8")
write.csv(expNonClassicalNames,file=con,row.names = FALSE, quote=FALSE)

con<-file("Nur77SetGene.csv",encoding="UTF-8")
write.csv(targetsGenes,file=con,row.names = FALSE, quote=FALSE)

con<-file("Nur77Set.csv",encoding="UTF-8")
write.csv(expNur77Target,file=con,row.names = TRUE)

con<-file("Nur77SpeciesOverlapSet.csv",encoding="UTF-8")
write.table(expOverlap,file=con,sep=",",quote=FALSE,row.names = TRUE,col.names=TRUE)

con<-file("SingleCellRNAgroups.csv",encoding="UTF-8")
write.table(IDlist,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)

con<-file("Nur77SpeciesOverlapFiltered10.csv",encoding="UTF-8")
wwrite.table(expOverlapSetFiltered,file=con,sep=",",quote=FALSE,row.names = TRUE,col.names=TRUE)

con<-file("Nur77GeneSpeciesComparison.csv",encoding="UTF-8")
write.table(targetsComparison,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=FALSE)

#Export overarching gene list

genesOverall <- rownames(expMatrixFiltered)
print(genesOverall)

expFilteredPantherAnnot <- read.table("pantherGeneList.txt",sep="\t")

mapGeneNametoID <- expFilteredPantherAnnot[,c(2,1)]
colnames(mapGeneNametoID) <- c("GeneName","ID")
expFilterdPantherPathway <- read.table("pantherGeneListCategories.txt", sep="\t")

geneNameToPathway <- expFilterdPantherPathway[,c(3,2)]

geneNameToPathwayHash <- tidyr::separate_rows(geneNameToPathway,V2,sep=",")

geneNameToPathwayHash <- geneNameToPathwayHash[,c(2,1)]
colnames(geneNameToPathwayHash) <- c("ID","Pathway")

PathwayToColor <- as.data.frame(unique(geneNameToPathwayHash[,2]),stringsAsFactors=FALSE)


                         
con<-file("GeneList.csv",encoding="UTF-8")
write.table(genesOverall,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=FALSE)

con<-file("geneToID.csv",encoding="UTF-8")
write.table(mapGeneNametoID,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)
con<-file("IDtoPathway.csv",encoding="UTF-8")
write.table(geneNameToPathwayHash,file=con,sep=",",quote=FALSE,row.names = FALSE,col.names=TRUE)
