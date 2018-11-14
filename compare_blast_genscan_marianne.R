require(seqinr)


##Function to compute the number of genes found by blast
number.genes.blast<-function(blastpath){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  return(length(out_blast))
}

## Number of genes by genscan
number.genes.genscan<-function(genscanpath){
  out_gen<-read.csv(genscanpath)
  return(nrow(out_gen))
}

## Function to have overlaps between gene regions
## takes as argument one blast output and one genscan output
## Note : query is melanogaster
genes.overlaps<-function(genomepath,blastpath, genscanpath){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  out_gen<-read.csv(genscanpath)
  
  genome_file<-read.fasta(genomepath)
  # seq_length<-rep(0, length(genome_file))
  # for(i in 1:length(genome_file)){
  #   seq_length[i]<-length(genome_file[[i]])
  # }
  # ind<- which.max(seq_length)
  # print(seq_length)
  # print(ind)
  ind<-1
  name<-names(genome_file)[ind]
  
  genome_blast<-rep(0,length(genome_file[[ind]]))
  genome_gen<-rep(0,length(genome_file[[ind]]))
  
  for(g in 1:nrow(out_gen)){
    pos<-seq(from=out_gen[g,3],to=out_gen[g,4],by=1)
    genome_gen[pos]<-rep(1,length(pos))
  }
  print(genome_gen)
  print(genome_gen==1)
  genome_blast <- rep(0,length(genome_file[[ind]]))
  
  indexes_blast<-which(out_blast[,2]==name)
  for(i in indexes_blast){
    if(out_blast[i,9]<=out_blast[i,10]){
      genome_blast[seq(from=out_blast[i,9],to=out_blast[i,10],by=1)]<-rep(1,out_blast[i,10]-out_blast[i,9]+1)
    }
  }
  return(data.frame(overlap.number=sum(genome_gen+genome_blast==2), overlap.over.genscan=sum(genome_gen+genome_blast==2)/sum(genome_gen==1), overlap.over.blast=sum(genome_gen+genome_blast==2)/sum(genome_blast==1)))

}

blast_path<-'/local/data/public/g1/out_ere_n'
genome_path<-'/local/data/public/g1/chromosomes/dere-all-chromosome-r1.05.fasta'
genscan_path<-'/local/data/public/g1/Genscan/genes.csv'

#print(sprintf("number of genes blast %d",number.genes.blast(blast_path)))
#print(sprintf("number of genes genescan %d",number.genes.genscan(genscan_path)))
print(genes.overlaps(genome_path,blast_path,genscan_path))