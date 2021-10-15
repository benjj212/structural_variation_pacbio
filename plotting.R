#### extracting bed file of neigbouring genes
### packages

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
TAIR10_TE <- TAIR10[grep("transposa", TAIR10$feature),]
TAIR10_TE <- TAIR10_TE[-which(TAIR10_TE$feature=="transposable_element_gene"),]
TAIR10_TE$attributes <- substr(TAIR10_TE$attributes,4,13)

args <- commandArgs(trailingOnly = TRUE)
#gene_AT <- as.character(args[1])
gene_AT <- "AT1G31390"
#number_of_side_genes <- as.numeric(args[2])
number_of_side_genes <- 4
#path_for_file <- args[3]
path_for_file <- "/groups/nordborg/projects/tumv_resistance/005_results/018_structural_rearrangements/"
###defining parameters
xlim_plot=c(-100000, 100000)
CHR <- as.numeric(substr(x = gene_AT, 3,3))
centering_1 <- number_of_side_genes-2
centering_2 <- number_of_side_genes+2
gene_with_link <- seq(1,(number_of_side_genes*2)+1,1)

bed_file <- c()
file_name <- c()
for(i in seq(-number_of_side_genes, number_of_side_genes,1)){
bed_file <- rbind(bed_file,c(TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+i,c(1,4,5)], TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+i,9]))
file_name <- c(file_name, TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+i,9])
}
start_bed <- bed_file[1,2]
end_bed <- bed_file[length(bed_file[,1]),3]
chr_te <- bed_file[1,1]
te_bed_file <- TAIR10_TE[which(TAIR10_TE$start>start_bed & TAIR10_TE$end<end_bed & TAIR10_TE$seqname==chr_te),c(1,4,5,9)]
file_name <- c(file_name,te_bed_file[,4])

setwd("/groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/")
PACBIO <- dir()
PACBIO <- PACBIO[-grep(".nhr",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nin",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nog",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nsd",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nsi",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".nsq",PACBIO, fixed=T)]
PACBIO <- PACBIO[-grep(".fai",PACBIO, fixed=T)]

#### setting up the order of the accessions based on multiple alignement



#### setting up the color
setwd(paste(path_for_file,"/",gene_AT,"/", sep=""))
if(number_of_side_genes==3){
color_genes <- c("#b2182b", "#ef8a62", "#fddbc7", "#4d4d4d", "#d1e5f0", "#67a9cf", "#2166ac")
color_genes_trans <- c("#b2182b1A", "#ef8a621A", "#fddbc71A", "#4d4d4d1A", "#d1e5f01A", "#67a9cf1A", "#2166ac1A")
#cbind(seq(1,length(file_name)-length(color_genes),1), "orange")[,2]
if(length(file_name)>length(color_genes)){
color_genes <- cbind(file_name, c(color_genes, cbind(seq(1,length(file_name)-length(color_genes),1), "orange")[,2]))[,2]
}
}

if(number_of_side_genes==4){
  color_genes <- c("#b2182b", "#ef8a62", "#fddbc7", "#4d4d4d", "#4d4d4d","#4d4d4d","#d1e5f0", "#67a9cf", "#2166ac")
  color_genes_trans <- c("#b2182b1A", "#ef8a621A", "#fddbc71A", "#f7f7f71A", "#4d4d4d1A","#f7f7f71A","#d1e5f01A", "#67a9cf1A", "#2166ac1A")
  
}

if(number_of_side_genes==6){
  color_genes <- c("#b2182b", "#ef8a62", "#fddbc7", "#f7f7f7", "#f7f7f7", "#4d4d4d","#f7f7f7","#f7f7f7","#d1e5f0", "#67a9cf", "#2166ac")
  color_genes_trans <- c("#b2182b1A", "#ef8a621A", "#fddbc71A", "#f7f7f71A","#f7f7f71A","#4d4d4d1A", "#f7f7f71A","#f7f7f71A","#d1e5f01A", "#67a9cf1A", "#2166ac1A")
}


pdf(paste(path_for_file,"/",gene_AT,"/",gene_AT,".plot.pdf",sep=""), height=17, width = 15)

par(cex=1.5)
plot(0, type="n", xlim=xlim_plot, ylim=c(0, (length(PACBIO)*10)+40), xlab="relative coordinates in bp", axes = F, ylab="", main=gene_AT)
axis(1)
axis(2)

### plotting gene structure
setwd(paste(path_for_file,"/",gene_AT,"/", sep=""))
#start_gene_D <- (read.table(paste("6909.fasta", file_name[centering_1], ".bed.fasta.70.txt", sep=""))[1,7]+read.table(paste("6909.fasta", file_name[centering_2], ".bed.fasta.70.txt", sep=""))[1,8])/2
start_gene_D <- (TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)-3,"start"]+TAIR10_gene[which(TAIR10_gene$attributes==gene_AT)+3,"end"])/2
TAIR10_sub <- TAIR10
TAIR10_sub$attributes <- substr(TAIR10$attributes,4,13)
#position_AT_gene[k,2]
for(k in c(1:length(file_name))){
  structure_gene <- TAIR10_sub[grep(file_name[k],TAIR10_sub$attributes),]
  if(length(which(structure_gene$feature=="exon"))>1){
    structure_gene <- structure_gene[which(structure_gene$feature=="exon"),]
  }
  structure_gene[,4] <- structure_gene[,4]-start_gene_D
  structure_gene[,5] <- structure_gene[,5]-start_gene_D
  rect(xleft = min(structure_gene[,4]), xright = max(structure_gene[,5]), ybottom = (length(PACBIO)*10)+11.5, ytop = (length(PACBIO)*10)+12.5, col = color_genes[k])
  rect(xleft = structure_gene[,4], xright = structure_gene[,5], ybottom = (length(PACBIO)*10)+10, ytop = (length(PACBIO)*10)+14, col = color_genes[k])
  text(mean(structure_gene[,4]),(length(PACBIO)*10)+28, adj =0.5, labels= file_name[k], cex=0.4, srt = -45)
}



#for columbia 
start_gene_D <- (read.table(paste("6909.fasta", file_name[centering_1], ".bed.fasta.70.txt", sep=""))[1,7]+read.table(paste("6909.fasta", file_name[centering_2], ".bed.fasta.70.txt", sep=""))[1,8])/2

for(j in c(1:length(file_name))){
  if(file.size(paste("6909.fasta", file_name[j], ".bed.fasta.70.txt", sep=""))>280){
  gene_coord <- read.table(paste("6909.fasta", file_name[j], ".bed.fasta.70.txt", sep=""))
  gene_coord[,7] <- gene_coord[,7]-start_gene_D
  gene_coord[,8] <- gene_coord[,8]-start_gene_D
  rect(xleft = gene_coord[,7], xright = gene_coord[,8], ytop = (length(PACBIO)*10)+5, ybottom = (length(PACBIO)*10), col=color_genes[j])
  }
}

#for PACBIO
PACBIO_order <- as.character(read.table(paste(path_for_file,gene_AT, "/order/", gene_AT, ".txt", sep=""))[,1])

for(i in c(1:length(PACBIO_order))){
  start_gene_D <- (read.table(paste(PACBIO_order[i],".fasta", file_name[centering_1], ".bed.fasta.70.txt", sep=""))[1,7]+read.table(paste(PACBIO_order[i],".fasta", file_name[centering_2], ".bed.fasta.70.txt", sep=""))[1,8])/2
  
  for(j in c(1:length(file_name))){
    if(file.size(paste(PACBIO_order[i],".fasta", file_name[j], ".bed.fasta.70.txt", sep=""))>280){
    gene_coord <- read.table(paste(PACBIO_order[i],".fasta", file_name[j], ".bed.fasta.70.txt", sep=""))
    gene_coord[,7] <- gene_coord[,7]-start_gene_D
    gene_coord[,8] <- gene_coord[,8]-start_gene_D
    #gene_coord <- gene_coord[which(gene_coord[,7]>-25000 & gene_coord[,8]<50000),]
    #if(length(gene_coord[,1])>=1){
    rect(xleft = gene_coord[,7], xright = gene_coord[,8], ytop = ((i-1)*10)+5, ybottom = (i-1)*10, col = color_genes[j])
    #}
    }
  }
  text(x = xlim_plot[2]-3000,y = ((i-1)*10)+2.5, labels =PACBIO_order[i], cex=0.7)
}

### adding the link between the genes 

for(i in c(1:(length(PACBIO_order)-1))){
  gene_with_link_n <- c()
  for(l in gene_with_link){
  if(file.size(paste(PACBIO_order[i],".fasta", file_name[l], ".bed.fasta.70.txt", sep=""))>280){
  gene_with_link_n <- c(gene_with_link_n, l)
  }
  }
  gene_with_link_ne <- c()
  for(l in gene_with_link){
    if(file.size(paste(PACBIO_order[i+1],".fasta", file_name[l], ".bed.fasta.70.txt", sep=""))>280){
      gene_with_link_ne <- c(gene_with_link_ne, l)
    }
  }
  if(length(gene_with_link_n)>length(gene_with_link_ne)){
    gene_with_link_n <-gene_with_link_ne 
  }else{
    gene_with_link_n <- gene_with_link_n
  }
  
  link1 <- c()
  for(l in gene_with_link_n){
    if(file.size(paste(PACBIO_order[i],".fasta", file_name[l], ".bed.fasta.70.txt", sep=""))>280){
    gene_coord <- read.table(paste(PACBIO_order[i],".fasta", file_name[l], ".bed.fasta.70.txt", sep=""))
    link1 <- rbind(link1,c(gene_coord[1,7]-start_gene_D, gene_coord[1,8]-start_gene_D))
    }
  }
  start_gene_D <- (read.table(paste(PACBIO_order[i+1],".fasta", file_name[centering_1], ".bed.fasta.70.txt", sep=""))[1,7]+read.table(paste(PACBIO_order[i+1],".fasta", file_name[centering_2], ".bed.fasta.70.txt", sep=""))[1,8])/2
  
  link1_1 <- c()
  for(l in gene_with_link_n){
    gene_coord <- read.table(paste(PACBIO_order[i+1],".fasta", file_name[l], ".bed.fasta.70.txt", sep=""))
    link1_1 <- rbind(link1_1,c(gene_coord[1,7]-start_gene_D, gene_coord[1,8]-start_gene_D))
  }
  link1_2 <- link1_1[which(link1_1[,1]>-100000 & link1_1[,1]<100000),]
  link2 <- link1[which(link1_1[,1]>-100000 & link1_1[,1]<100000),]
  link1_2 <- link1_2[which(link2[,1]>-100000 & link2[,1]<100000),]
  link2 <- link2[which(link2[,1]>-100000 & link2[,1]<100000),]
  if(length(link2[,1]>=1) & length(link1_2[,1]>=1) & length(link2[,1])==length(link1_2[,1])){
  for(m in c(1:length(link2[,1]))){
    polygon(c(link2[m,1],link2[m,2], link1_2[m,2], link1_2[m,1]), c((i*10)-5,(i*10)-5,i*10, i*10), col = color_genes_trans[gene_with_link_n[m]])
  }
  }
}

text(xlim_plot[2]-5000, (length(PACBIO_order)*10)+7.5)

dev.off()
