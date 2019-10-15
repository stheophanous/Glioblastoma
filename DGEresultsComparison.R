rm(list = ls())

myanalysis <- read.table(file='C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\EdgeDESeq05.txt', sep= ",")
myanalysis <- as.list(as.character(myanalysis$V2))
myanalysis <- myanalysis[-1]

deseq2tool <- read.table(file='C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\OnlineTool\\ComparisonResultsNEW\\deseq2_only.txt')
deseq2tool <- as.list(rownames(deseq2tool))
edgertool <- read.table(file='C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\OnlineTool\\ComparisonResultsNEW\\edger_only.txt')
edgertool <- as.list(rownames(edgertool))

toolanalysis <- read.table(file='C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\OnlineTool\\ComparisonResultsNEW\\edger_deseq2_overlap.txt')
toolanalysis <- as.list(rownames(toolanalysis))

commongenes <- as.data.frame(intersect(myanalysis, toolanalysis))
write.table(commongenes, file = 'C:\\Users\\medsthe\\Desktop\\Glioblastoma\\DifferentialGeneExpressionAnalysis\\ComparisonToolMeNEW.txt', sep = '\t')

library(VennDiagram)
draw.quad.venn(area1=230, 
               area2=54, 
               area3=206, 
               area4=22, 
               n12=45, 
               n34=22, 
               n1234=21, 
               n123=43, 
               n124=21, 
               n134= 22, 
               n234=21, 
               n13=88, 
               n14=22, 
               n23=50, 
               n24=21, category = c('EdgeR-R', 'DESeq2-R', 'EdgeR-Gliovis', 'DESeq2-Gliovis'), fill = c("yellow", "red", "green", "blue"))
