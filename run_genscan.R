# Launch genscan on species chromosomess

library(seqinr)

Dsim <-'chromosomes/dsim-all-chromosome-r2.01.fasta'
genscan <- '/local/data/public/genome_informatics_2018/programs/genscan/'
setwd('/local/data/public/g1')

# We only take the first sequence, 
# which is much longer than all the others (I don't know what they are)
chr<-read.fasta(Dsim)[[0]]
# We split this long sequences in short ones
seq_length <- 1e5
seq_overlap <- 1e4
n_slices <- as.integer(length(chr) / (seq_length-seq_overlap))
for (i in 1:n_slices){
  write.fasta(chr[((i-1)*seq_length - (i-1)*seq_overlap):
                    (i*seq_length-(i-1)*seq_overlap)], 
              paste((i-1)*seq_length - (i-1)*seq_overlap,":",i*seq_length-(i-1)*seq_overlap,sep=""), 
              file.out=paste('chromosomes/Dsim/slices/slice_',i,'.fa',sep=""),
              open="w",as.string=FALSE)
  outex<-system(paste(genscan,'genscan',' ',
                      genscan,'lib/HumanIso.smat',' ',
                      'chromosomes/Dsim/slices/slice_',i,'.fa',' -cds',sep="")
                ,intern=TRUE)
  print(paste("ran genscan on slice",i,"over",n_slices))
  if (i == 1){
    print(outex)
  }
  out_file <- file(paste("Genscan/Dsim/slice_",i,".txt",sep=""),open="w")
  write(outex,out_file)
}
