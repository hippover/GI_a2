require(seqinr)
library(data.table)


##Function to compute the number of genes found by blast
number.genes.blast<-function(blastpath){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  return(length(out_blast))
}

## Number of genes by genscan
number.genes.genscan<-function(genscanpath){
  out_gen<-read.csv(genscanpath)
  genes<-data.table(out_gen)
  genes<-genes[,list(id=first(gene_id),gene_start=min(gene_start), gene_end=max(gene_end)),by=genes$gene_id]
  return(nrow(genes))
}

## Function to have overlaps between gene regions
## takes as argument one blast output and one genscan output
## Note : query is melanogaster
genes.overlaps<-function(species,genomepath,blastpath, genscanpath, cdnapath){
  genome_file<-read.fasta(genomepath)
  
  genome_blast<-vector.match.blast(species,blastpath,genome_file)
  genome_gen<-vector.match.genscan(species,genscanpath,genome_file)
  genome_cdna<-vector.cdna(species,cdnapath,genome_flie)
  
  return(data.frame(overlap.blast=sum(genome_cdna+genome_blast==2),overlap.genscan=sum(genome_cdna+genome_genscan==2), 
                    overlap.blast.percent=sum(genome_cdna+genome_blast==2)/sum(genome_cdna==1), overlap.genscan.percent=sum(genome_cdna+genome_genscan==2)/sum(genome_cdna==1)))

}

## Function that returns the vector of matching genes for blast
vector.match.blast<-function(specie,blastpath, genome_file){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  
  seq_length<-rep(0, length(genome_file))
  for(i in 1:length(genome_file)){
    seq_length[i]<-length(genome_file[[i]])
  }
  ind<- which.max(seq_length)
  print(seq_length)
  print(ind)
  ind<-1
  name<-names(genome_file)[ind]
  
  genome_blast<-rep(0,length(genome_file[[ind]]))
  indexes_blast<-which(out_blast[,2]==name)
  for(i in indexes_blast){
    if(out_blast[i,9]<=out_blast[i,10]){
      pos<-seq(from=out_blast[i,9],to=out_blast[i,10],by=1)
      genome_blast[pos]<-rep(1,length(pos))
    }
  }
  file_name<-paste(specie,'_blast_vector.txt')
  write(genome_blast,file_name)
  return(genome_blast)
}

## Create the same vector for the cdna file
vector.match.cdna<-function(specie,cdnapath,genome_file){
  cda_file<-read.fasta(cdnapath)
  
  file_name<-paste(specie,'_cdna_vector.txt')
  write(genome_cdna,file_name)
  
}


## Function returning the genscan output as a vector of matching positions
vector.match.genscan<-function(specie,genscanpath,genome_file){
  out_gen<-read.csv(genscanpath)
  genes<-data.table(out_gen)
  
  #genes<-genes[,list(id=first(gene_id),gene_start=min(gene_start), gene_end=max(gene_end)),by=(genes$gene_id)]
  scf_list<-genes[,scf=first(scf),by=(genes$scf)]
  print(scf_list)
  genome_gen<-vector(mode="list",length=length(names(genome_file)))
  names(genome_gen)<-names(genome_file)
  for(i in 1:length(names(genome_file))){
    genome_gen[[i]]<-rep(0,length(genome_file[[i]]))
  }
  for(i in 1:length(names(genome_file))){
    scf<-names(genome_file)[i]
    print(scf)
    genes_scf<-genes[genes$scf==scf]
    print(genes_scf)
    genes_scf<-genes_scf[,list(id=first(gene_id),gene_start=min(gene_start), gene_end=max(gene_end)),by=(genes_scf$gene_id)]
    
    for(g in 1:nrow(genes_scf)){
      pos<-seq(from=genes_scf[g]$gene_start,to=genes_scf[g]$gene_end,by=1)
      genome_gen[[i]][pos]<-rep(1,length(pos))
    }
  }
  file_name<-paste(specie,'_genscan_vector.txt')
  write(genome_gen,file_name)
 
  return(genome_gen)
}

blast_path<-'/local/data/public/g1/out_ere_n'
genome_path<-'/local/data/public/g1/chromosomes/dere-all-chromosome-r1.05.fasta'
genscan_path<-'/local/data/public/g1/Genscan/genes.csv'
cdna_path<-'/local/data/public/g1/chromosomes/dere-all-transcript-r1.3.fasta'

#print(sprintf("number of genes blast %d",number.genes.blast(blast_path)))
#print(sprintf("number of genes genescan %d",number.genes.genscan(genscan_path)))
print(genes.overlaps('dere',genome_path,blast_path,genscan_path))

#out_gen<-read.csv(genscan_path)
#genes<-data.table(out_gen)
#genes<-genes[,list(id=first(gene_id),gene_start=min(gene_start), gene_end=max(gene_end)),by=genes$gene_id]