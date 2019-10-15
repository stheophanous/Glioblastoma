rm(list = ls())
library(data.table)
# Import data
# 1. FPKM counts
primary.FPKM <- as.data.frame(fread(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedFPKM_ALLsamples.txt'))
fpkmlist <- as.list(as.character(colnames(primary.FPKM)))[-1]
# 2. Primary tumour rawcounts
primary.raw <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ForStelios\\PrimaryGBM_rawcounts.txt', sep ="\t", header=TRUE)
row.names(primary.raw) <- primary.raw[,1]
#primary.raw$EnsID <- sub("\\..*", "", primary.raw$EnsID)
#primary.clean <- subset(primary.raw, EnsID %in% fpkmlist)
primary.clean <- primary.raw[,-1:-3]

# NEW Primary tumour rawcounts
primary.rawNEW <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ForStelios\\GSEA\\NewK_STDLOCAL_rawcounts.txt', sep ="\t", header=TRUE)
rownames(primary.rawNEW) <- primary.rawNEW$EnsID
primary.rawNEW <- primary.rawNEW[,-c(1:3)]
primary.rawNEW <- primary.rawNEW[,c(1:8)]

# NEW non standard treatment Primary tumour rawcounts
primary.rawNEWnonstd <- read.table(file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\ForStelios\\GSEA\\NewK_LOCALnSTD_rawcounts.txt', sep ="\t", header=TRUE)
rownames(primary.rawNEWnonstd) <- primary.rawNEWnonstd$EnsID
primary.rawNEWnonstd <- primary.rawNEWnonstd[,-c(1:3)]
primary.rawNEWnonstd <- primary.rawNEWnonstd[,c(1:10)]


#Merge datasets
primary.merged <- merge(primary.clean, primary.rawNEW, by = "row.names", all = TRUE)
rownames(primary.merged) <- primary.merged$Row.names
primary.merged <- primary.merged[,-1]
primary.merged <- merge(primary.merged, primary.rawNEWnonstd, by = "row.names", all = TRUE)
primary.merged$Row.names <- sub("\\..*", "", primary.merged$Row.names)
primary.merged <- subset(primary.merged, Row.names %in% fpkmlist)
rownames(primary.merged) <- primary.merged$Row.names
primary.merged <- primary.merged[,-1]
primary.merged[is.na(primary.merged)] <- 0
write.table(primary.merged, file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\mergedRAW_ALLsamples.txt', sep = "\t")
