# Launch genscan on species chromosomess

library(seqinr)
library(data.table)

Dsim <-'chromosomes/dsim-all-chromosome-r2.01.fasta'
Dmoj <- 'chromosomes/dmoj-all-chromosome-r1.04.fasta'
Dsec <- 'chromosomes/dsec-all-chromosome-r1.3.fasta'
Dere <- 'chromosomes/dere-all-chromosome-r1.05.fasta'

files=list(Dsim=Dsim,Dmoj=Dmoj,Dsec=Dsec,Dere=Dere)

genscan <- '/local/data/public/genome_informatics_2018/programs/genscan/'
setwd('/local/data/public/g1/')
system()

for (j in 1:length(files)){
  print("----------------------")
  f <- files[[j]]
  species <- names(files)[[j]]
  print(f)
  print(species)
  # We only take the first sequence, 
  # which is much longer than all the others (I don't know what they are)
  read <-read.fasta(f)
  chr <- read[1][[1]]
  # We split this long sequences in short ones
  seq_length <- 3e5
  seq_overlap <- 3e4
  n_slices <- as.integer(length(chr) / (seq_length-seq_overlap))
  print(paste(length(chr),"bases cut in",n_slices,"slices"))
  # For each sequence, we write a fasta file
  # We then run genscan on it, and save the results
  for (i in 1:n_slices){
    write.fasta(chr[((i-1)*seq_length - (i-1)*seq_overlap):
                      ((i*seq_length-1)-(i-1)*seq_overlap)], 
                paste((i-1)*seq_length - (i-1)*seq_overlap,":",i*seq_length-(i-1)*seq_overlap,sep=""), 
                file.out=paste('chromosomes/',species,'/slices/slice_',i,'.fa',sep=""),
                open="w",as.string=FALSE)
    outex<-system(paste(genscan,'genscan',' ',
                        genscan,'lib/HumanIso.smat',' ',
                        'chromosomes/',species,'/slices/slice_',i,'.fa',' -cds',sep="")
                  ,intern=TRUE)
    if (i %% 10 == 0) print(paste("ran genscan on slice",i,"over",n_slices))
  
    out_file <- file(paste("Genscan/",species,"/slice_",i,".txt",sep=""),open="w")
    write(outex,out_file)
    close(out_file)
  }
  
  # launch perl script
  # which parses the results and creates one csv containing all the lines
  print("-----------")
  print("Launching perl script")
  print(paste("perl parse_genscan_1.pl",
               seq_length,
               seq_overlap,
               n_slices,
               paste0('/local/data/public/g1/Genscan/',species,'_genes_temp.csv')))
  perl<-system(paste("perl parse_genscan_1.pl",
                     seq_length,
                     seq_overlap,
                     3,#n_slices,
                     paste0('/local/data/public/g1/Genscan/',species,'_genes_temp.csv'))
               ,intern=TRUE)
  print("Perl is done")
  
  # Then, load this table and remove duplicates or genes cut in the slicing process 
  # (only leave their full version)
  seqs <- data.table(read.csv(paste0('/local/data/public/g1/Genscan/',species,'_genes_temp.csv')))
  print(paste(nrow(seqs),"sequences"))
  genes <- seqs[,list(start=min(gene_start),end=max(gene_end),strand=max(strand)),by=gene_id]
  genes <- genes[rank(genes$end)]
  genes <- genes[,list(id=first(gene_id),num=length(gene_id)),by=list(start,strand)]
  print(paste("remove",sum(genes$num > 1),"overlaps"))
  seqs <- seqs[seqs$gene_id %in% genes$id]
  # FInally, save the table for analysis
  write.csv(seqs,file=paste0('/local/data/public/g1/Genscan/',species,'_genes.csv'))
  print(paste("wrote file containing",nrow(seqs),"sequences for",nrow(genes),"genes"))
}