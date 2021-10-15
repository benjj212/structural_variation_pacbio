
### loading modules
ml r/3.5.1-foss-2018b
ml bedtools/2.27.1-foss-2018b
### generating bed files from the coordinate in TAIR10

AT_number=$1
number_gene_side=5
path_to_output=/groups/nordborg/projects/tumv_resistance/005_results/018_structural_rearrangements/
Rscript /groups/nordborg/projects/tumv_resistance/005_results/018_structural_rearrangements/scripts/extraction_bed.R $AT_number $number_gene_side $path_to_output
path_to_pacbio=/groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/

### extracting fasta files

for bed_file in $path_to_output$AT_number/*.bed
do
bedtools getfasta -fi /groups/nordborg/projects/the1001genomesplus/tools/synteny_tool/data_TAIR/TAIR10_all.fa -bed $bed_file -name -fo $bed_file.fasta
done

echo "extraction of fasta for blast done"

### generating the database forblast only once
#for pacbio in /groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/*.fasta
#do
  #name_pacbio=$(basename $pacbio)
  #/groups/nordborg/projects/the1001genomesplus/tools/synteny_tool/ncbi-blast-2.7.1+/bin/makeblastdb -in /groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/$name_pacbio -parse_seqids -dbtype nucl
#done

#### blasting against all pacbio genomes

cd $path_to_output$AT_number/

for pacbio in /groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/*.fasta
do
name_pacbio=$(basename $pacbio)
for file in $path_to_output$AT_number/*.fasta
do
GENE_FASTA=$file
name_file=$(basename $GENE_FASTA)
OUT_TXT=$name_pacbio$name_file.70.txt
/groups/nordborg/projects/the1001genomesplus/tools/synteny_tool/ncbi-blast-2.7.1+/bin/blastn -query $GENE_FASTA \
-db /groups/nordborg/projects/the1001genomesplus/002_fasta_file/Chr_proper/$name_pacbio \
-task blastn \
-evalue 1 \
-out $OUT_TXT \
-perc_identity 70 \
-outfmt "7 qseqid qacc sacc evalue qstart qend sstart send sseqid" \
-penalty -2 \
-reward 1 \
-word_size 28 \
-gapopen 1 \
-gapextend 1
done
done
echo $AT_number
echo "Blast done"

### extracting sequences around the target gene
#Rscript /groups/nordborg/projects/the1001genomesplus/tools/synteny_tool/scripts/ordering_accessions.R $AT_number $path_to_output

### transforming bed to fasta
#cd $path_to_output$AT_number/
#for pacbio in $path_to_output$AT_number/*.fasta.bed
#do
#pacbio_n=$(basename $pacbio)
#accession_pacbio=${pacbio_n%.*}
#bedtools getfasta -fi $path_to_pacbio$accession_pacbio -bed $pacbio_n -name -fo $accession_pacbio.region.fasta
#done

#echo "extraction of the fasta from region done"

#cat *region.fasta > $AT_number.all_sequence.fasta

Rscript /groups/nordborg/projects/tumv_resistance/005_results/018_structural_rearrangements/scripts/plotting.R $AT_number $number_gene_side $path_to_output
echo "ploting done"
