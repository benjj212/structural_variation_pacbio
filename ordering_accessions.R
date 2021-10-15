### function
gffRead <- function(gffFile, nrows = -1) {
  cat("Reading ", gffFile, ": ", sep="")
  gff = read.table(gffFile, sep="\t", as.is=TRUE, quote="",
                   header=FALSE, comment.char="#", nrows = nrows,
                   colClasses=c("character", "character", "character", "integer",  
                                "integer",
                                "character", "character", "character", "character"))
  colnames(gff) = c("seqname", "source", "feature", "start", "end",
                    "score", "strand", "frame", "attributes")
  cat("found", nrow(gff), "rows with classes:",
      paste(sapply(gff, class), collapse=", "), "\n")
  stopifnot(!any(is.na(gff$start)), !any(is.na(gff$end)))
  return(gff)
}
TAIR10 <- gffRead("/groups/nordborg/projects/the1001genomesplus/tools/synteny_tool/data_TAIR/Araport11_GFF3_genes_transposons.201606.gff")
TAIR10_gene <- TAIR10[which(TAIR10$feature=="gene" | TAIR10$feature=="transposable_element_gene"),]
TAIR10_gene$attributes <- substr(TAIR10_gene$attributes,4,12)
args <- commandArgs(trailingOnly = TRUE)
#gene_AT <- as.character(args[1])
bed_file <- c()
file_name <- c()
for(i in seq(-number_of_side_genes, number_of_side_genes,1)){
  bed_file <- rbind(bed_file,c(TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+i,c(1,4,5)], TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+i,9]))
  file_name <- c(file_name, TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+i,9])
}



gene_AT <- "AT1G20400"
#path_for_file <- as.character(args[2])
path_for_file <- "/groups/nordborg/projects/tumv_resistance/005_results/018_structural_rearrangements/"
output_path <- "/groups/nordborg/projects/tumv_resistance/005_results/018_structural_rearrangements/"

dir.create(paste(path_for_file,gene_AT, "/order/", sep=""))
setwd("/groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/")
PACBIO <- dir()
PACBIO <- PACBIO[-grep(".nhr",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nin",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nog",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nsd",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nsi",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nsq",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".fai",PACBIO, fixed=T)]



setwd(paste(path_for_file, gene_AT, sep=""))

pacbio_range_all <- c()
for(i in c(1:length(PACBIO))){
  coordinate_min_max <- c()
for(k in c(1,7)){
    gene_check <- file_name[k]
TAIR_gene_length <- TAIR10_gene[which(TAIR10_gene$attributes==gene_check),"end"]-TAIR10_gene[which(TAIR10_gene$attributes==gene_check),"start"]
coordinate_gene <- read.table(paste(PACBIO[i], gene_check, ".bed.fasta.70.txt", sep=""))
coordinate_gene <- coordinate_gene[which((coordinate_gene[,6]-coordinate_gene[,5])>(TAIR_gene_length*0.4)),]
print(length(coordinate_gene[,1]))
if(length(coordinate_gene[,1])>0){
  coordinate_gene_2 <- coordinate_gene[,c(9,7,8)]
  coordinate_gene_1 <- coordinate_gene_2
for(k in c(1:length(coordinate_gene_2[,1]))){
print(k)

coordinate_gene_1[k,2] <- min(coordinate_gene_2[k,c(2,3)])
coordinate_gene_1[k,3] <- max(coordinate_gene_2[k,c(2,3)])
}
coordinate_gene_1 <- cbind(coordinate_gene_1, PACBIO[i])


}
coordinate_min_max <- rbind(coordinate_min_max, coordinate_gene_1[1,])

}
pacbio_range <- c(strsplit(PACBIO[i],split = ".", fixed = T)[[1]][1],min(coordinate_min_max[1,c(2,3)]), max(coordinate_min_max[2,c(2,3)]), coordinate_min_max[1,1])
pacbio_range_all <-  rbind(pacbio_range_all, pacbio_range)
}
pacbio_range_all[order(as.numeric(pacbio_range_all[,3])-as.numeric(pacbio_range_all[,2]), decreasing = F),1]

write.table(pacbio_range_all[order(as.numeric(pacbio_range_all[,3])-as.numeric(pacbio_range_all[,2]), decreasing = F),1], file=paste(output_path,gene_AT,"/order/", gene_AT,".txt", sep=""), col.names = F, row.names = F, quote = F, sep="\t")


