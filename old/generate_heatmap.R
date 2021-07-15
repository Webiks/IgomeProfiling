library("RColorBrewer")
library("gplots")
dropbox_dir="/Users/Oren/Dropbox/"
#dropbox_dir="D:/Dropbox/"
setwd(paste(dropbox_dir, "Projects/gershoni/Experiments/Exp22/analyses/hiv_vs_hcv/PSSM_score_Peptide_Motifs_Pvalues/", sep=""))
min_middle_max = "max"
DataName=paste(min_middle_max, "_merged_features_df", sep="")
DataName="HCV_HIV_selected_motifs.HITS"
DataFile=paste(DataName, ".csv", sep="")
Data=read.table(file=DataFile,sep=",",header=T)
row.names(Data) <- Data$sample_name
Data= Data [,-(dim(Data)[2])]  # remove sample namew
#Data= Data [,-(dim(Data)[2])]  # remove label
#Data= Data [,-1]  # remove index
names(Data)

#for color bar above figure
ClassesSize=paste(min_middle_max, "_merged_features_df_classes_size", sep="")

ClassesSizeFile=paste(ClassesSize, ".csv", sep="")
Sizes=as.numeric(scan(ClassesSizeFile, sep=','))
colors=c("green", "orange", "black", "red", "blue", "yellow", "purple", "grey")
ColSideColors = c()
for (i in 1:length(Sizes)) { 
  ColSideColors = c(ColSideColors, rep(colors[i], Sizes[i]))
}

#extract biological condition names
ClassesNames=paste(min_middle_max, "_merged_features_df_order", sep="")
ClassesNamesFile=paste(ClassesNames, ".csv", sep="")
Names=scan(ClassesNamesFile, sep=',', what=character())

legend_labels = c()
for (i in 1:length(Sizes)) { 
  legend_labels = c(legend_labels, paste("(",Sizes[i],") ",Names[i], sep =''))
}

NumOfMotifs=length(Data)
pdf(paste(DataName, ".pdf", sep=""), width = 17, height = 20)
# res=heatmap.2(as.matrix(log(Data+1)),col=brewer.pal(8,"Blues"),trace="none",dendrogram='row', Rowv=TRUE, Colv=FALSE,cexRow=1.1,cexCol=1.1,labCol = "",xlab=paste(NumOfMotifs,"motifs"),margins=c(3,10),key=TRUE, symkey=FALSE, density.info="none",key.xlab="Log(Number of hits+1)", ColSideColors=ColSideColors)
res=heatmap.2(as.matrix(log(Data+1)),col=brewer.pal(8,"Blues"),trace="none",dendrogram='row', Rowv=TRUE, Colv=TRUE,cexRow=1.1,cexCol=1.1,labCol = "",xlab=paste(NumOfMotifs,"motifs"),margins=c(3,10),key=TRUE, symkey=FALSE, density.info="none",key.xlab="-Log2(pvalue)", ColSideColors=ColSideColors)
#res=heatmap.2(as.matrix(-log(Data)),col=brewer.pal(8,"Blues"),trace="none",dendrogram='none', Rowv=FALSE, Colv=FALSE,cexRow=1.1,cexCol=1.1,labCol = "",xlab=paste(NumOfMotifs,"motifs"),margins=c(3,10),key=TRUE, symkey=FALSE, density.info="none",key.xlab="-Log2(pvalue)", ColSideColors=ColSideColors)
par(lend = 20)          
legend("topright", legend = legend_labels, col = colors, lty= 1, lwd = 10)
dev.off()
