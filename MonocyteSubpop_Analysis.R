library(dplyr)


expMatrix <- read.table("Monocyte.txt",sep="\t")

expMatrixMonocytes <- expMatrix[,grep("Mono", colnames(expMatrix) )]

expMatrixMonocytesSubset <- expMatrixMonocytes[c(1:2000),]

expMatrixMonocytesSubset2 <- expMatrixMonocytes[c(1:2000),]

filterList <- read.table("ToFilter.txt")

t <- which(colnames(expMatrixMonocytes) %in% filterList)

expMatrixMonocytes <- expMatrixMonocytes[,which(!(colnames(expMatrixMonocytes) %in% filterList[,1]))]

expMatrixClassical <- expMatrixMonocytes[,grep("classical", colnames(expMatrixMonocytes))]

expMatrixNonClassical <- expMatrixMonocytes[,grep("nonclassical", colnames(expMatrixMonocytes))]

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


targets <- read.table("Nur77Targets.txt",sep=",",stringsAsFactors=FALSE)


targetsGenes <- targets[-1,15]

expNur77Target <- expComparison[which(rownames(expComparison) %in% targetsGenes),]


con<-file("expMatrix.csv",encoding="UTF-8")
write.csv(expMatrixFiltered,file=con,row.names = TRUE)

con2<-file("MonocyteComparison.csv",encoding="UTF-8")
write.csv(expFiltered,file=con,row.names = TRUE)

con3<-file("ClassicalSet.csv",encoding="UTF-8")
write.csv(expClassicalNames,file=con,row.names = FALSE,quote=FALSE)

con3<-file("NonclassicalSet.csv",encoding="UTF-8")
write.csv(expNonClassicalNames,file=con,row.names = FALSE, quote=FALSE)

con3<-file("Nur77SetGene.csv",encoding="UTF-8")
write.csv(targetsGenes,file=con,row.names = FALSE, quote=FALSE)

con4<-file("Nur77Set.csv",encoding="UTF-8")
write.csv(expNur77Target,file=con,row.names = TRUE)