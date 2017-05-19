#' Counts and information from and RNA-seq experiment
#'
#' RNA-seq counts of 1000 genes and 87 individuals
#'
#' @format A matrix (counts) and a dataframe (info)
#' \describe{
#'   \item{counts}{RNA-seq read counts}
#'   \item{info}{id, genotype, treatment and sample order of counts columns}
#' }

# rm(list = ls())
# counts<-read.delim("/Users/John/Desktop/ph2016/rawCounts/Phallii_gsnap.JGI_ID.Genome.All.counts", header=T)
# counts<-counts[,-grep("FH", colnames(counts))]
# libsize<-colSums(counts)
# hist(libsize, breaks=30)
# badlibs<-names(libsize)[libsize<=1e7]
# counts<-counts[,-which(colnames(counts) %in% badlibs)] #toss a small library
#
# info<-read.csv("JGI_Samples_3_6_14.csv", stringsAsFactors=F, header=T)
# info$libname<-with(info, paste(ID.initial, N, Plant.Number, sep="_"))
# info$libname[info$libname == "FIL2_H328_383"]<-"FIL2_H328_83"
# info<-info[,c("libname","id", "Treatment","Order_July")]
# info<-info[info$id %in% c("F1","HAL2","FIL2"),]
# colnames(info)<-c("id","geno","trt","order")
# info<-info[order(info$geno, info$trt),]
# rownames(info)<-info$id
#
# ids<-intersect(colnames(counts),info$id)
# info<-info[ids,]
# info$geno<-factor(info$geno, levels = c("FIL2","F1","HAL2"))
# info$trt<-factor(info$trt, levels = c("Wet","Dry"))
# counts<-counts[,ids]
# counts<-counts[order(rownames(counts)),]
# counts<-as.matrix(counts)
# rownames(counts)<-gsub(".v2.0","", rownames(counts), fixed=T)
#
# set.seed(42)
# counts<-counts[sample(1:nrow(counts),1000),]
# save(counts, file = "/Users/John/Desktop/deTools/deTools/data/counts.rda")
# save(info, file = "/Users/John/Desktop/deTools/deTools/data/info.rda")
