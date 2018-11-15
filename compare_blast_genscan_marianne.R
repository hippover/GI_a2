library(seqinr)
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
genes.overlaps<-function(blastvect, genscanvect, cdnablastvect){
  # genome_file<-read.fasta(genomepath)
  # 
  # genome_blast<-vector.match.blast(species,blastpath,genome_file)
  # genome_gen<-vector.match.genscan(species,genscanpath,genome_file)
  # genome_cdna<-vector.match.blast(species,cdnablastpath,genome_flie)
  
  return(data.frame(overlap.blast=sum(cdnablastvect+blastvect==2),overlap.genscan=sum(cdnablastvect+genscanvect==2), 
                    overlap.blast.percent=sum(cdnablastvect+blastvect==2)/sum(cdnablastvect==1), overlap.genscan.percent=sum(cdnablastvect+genscanvect==2)/sum(genome_cdna==1)))

}

## Function that returns the vector of matching genes for blast
vector.match.blast<-function(specie,blastpath, genome){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  
  out<-unique(out_blast[,c(2,3,9,10)])
  annotations <- vector(mode="list", length=length(genome))
  names(annotations)<-names(genome)
  sum<-rep(0,length(genome))
  for(i in 1:length(genome)){
    annotations[[i]]<-rep(F,length(genome[[i]]))
  }
  for (l in 1:nrow(out)){
    if(out[l,1] %in% names(annotations) & out[l,2]>80){
      i<-which(names(annotations)==out[l,1])
      annotations[[i]][seq(min(out[l,3],out[l,4]),max(out[l,3],out[l,4]),1)]<-T
    }
  }
  
  file_name<-paste(specie,'_blast_vector')
  save(annotations, file= file_name)
  #write.list(annotations,file_name)
  return(annotations)
}

## Create the same vector for the cdna file
# vector.match.cdna<-function(specie,cdnapath,genome){
#   cdna_file<-read.fasta(cdnapath)
#   names_cdna<-names(cdna_file)
#   genome_cdna<- vector(mode="list", length=length(genome))
#   names(genome_cdna)<-names(genome)
#   sum<-rep(0,length(genome))
#   for(i in 1:length(genome)){
#     annotations[[i]]<-rep(F,length(genome[[i]]))
#   }
#   for(i in 1:length(names(genome))){
#     genome_cdna[[i]]<-rep(F,length(genome[[i]]))
#     n<- names(genome)[i]
#     if(n %in% names_cdna){
#       index<-which(names_cdna %in% c(n))
#       pos<- seq
#     }
#   }
#   
#   
#   file_name<-paste(specie,'_cdna_vector.txt')
#   write(genome_cdna,file_name)
#   
# }


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
  #file_name<-paste(specie,'_genscan_vector.txt')
  #write(genome_gen,file_name)
  
 
  return(genome_gen)
}

## Function that returns the vector of matching genes for blast
vector.match.genscan2<-function(specie,blastpath, genome){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  
  out<-unique(out_blast[,c(2,3,9,10)])
  annotations <- vector(mode="list", length=length(genome))
  names(annotations)<-names(genome)
  sum<-rep(0,length(genome))
  for(i in 1:length(genome)){
    annotations[[i]]<-rep(F,length(genome[[i]]))
  }
  for (l in 1:nrow(out)){
    if(out[l,1] %in% names(annotations) & out[l,2]>80){
      i<-which(names(annotations)==out[l,1])
      annotations[[i]][seq(min(out[l,3],out[l,4]),max(out[l,3],out[l,4]),1)]<-T
    }
  }
  
  file_name<-paste(specie,'_blast_vector')
  save(annotations, file= file_name)
  #write.list(annotations,file_name)
  return(annotations)
}

blast_path<-'/local/data/public/g1/out_ere_n'
genome_path<-'/local/data/public/g1/chromosomes/dere-all-chromosome-r1.05.fasta'
genscan_path<-'/local/data/public/g1/Genscan/output_test_marianne.csv'
cdna_path<-'/local/data/public/g1/chromosomes/dere-all-transcript-r1.3.fasta'

#print(sprintf("number of genes blast %d",number.genes.blast(blast_path)))
#print(sprintf("number of genes genescan %d",number.genes.genscan(genscan_path)))
#print(genes.overlaps('dere',genome_path,blast_path,genscan_path))


Dsim<-'/local/data/public/g1/chromosomes/dsim-all-chromosome-r2.01.fasta'
Dere<-'/local/data/public/g1/chromosomes/dere-all-chromosome-r1.05.fasta'
Dmoj<-'/local/data/public/g1/chromosomes/dmoj-all-chromosome-r1.04.fasta'
Dsec<-'/local/data/public/g1/chromosomes/dsec-all-chromosome-r1.3.fasta'
names<-c('sim','moj','ere','sec')
genomes<-c(Dsim,Dmoj,Dere,Dsec)
for(i in 1:length(names)){
  genome<-read.fasta(genomes[i])
  n<-names[i]
  blast_path<-paste(paste('/local/data/public/g1/out_',n,sep=""),'_c',sep="")
  vector.match.blast(paste(n,'_cdna',sep=""), blast_path, genome)
}
