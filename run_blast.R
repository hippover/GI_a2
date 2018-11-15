setwd('/local/data/public/g1/')
Dmel_n<-'chromosomes/Drosophila_melanogaster.BDGP6.cdna.all.fa'
Dmel_p<-'chromosomes/Drosophila_melanogaster.BDGP6.pep.all.fa'
Dsim<-'chromosomes/dsim-all-chromosome-r2.01.fasta'
Dere<-'chromosomes/dere-all-chromosome-r1.05.fasta'
Dmoj<-'chromosomes/dmoj-all-chromosome-r1.04.fasta'
Dsec<-'chromosomes/dsec-all-chromosome-r1.3.fasta'

bl_dr<-'/local/data/public/genome_informatics_2018/programs/ncbi-blast-2.5.0+/bin/'

##make databases
#out_db_sim<-system(paste0(bl_dr,'makeblastdb -in ',Dsim,' -dbtype nucl -out dsim_db'),intern=T)
#out_db_ere<-system(paste0(bl_dr,'makeblastdb -in ',Dere,' -dbtype nucl -out dere_db'),intern=T)
#out_db_moj<-system(paste0(bl_dr,'makeblastdb -in ',Dmoj,' -dbtype nucl -out dmoj_db'),intern=T)
#out_db_sec<-system(paste0(bl_dr,'makeblastdb -in ',Dsec,' -dbtype nucl -out dsec_db'),intern=T)
#print("Databases created")

##align the genome onto the databases
##Dsim
#out_sim_n<-system(paste0(bl_dr,'blastn -db dsim_db -query ',Dmel_n,' -out out_sim_n -outfmt 6'),intern=T)
#print("dsim n alignment done")
#out_sim_p<-system(paste0(bl_dr,'tblastn -db dsim_db -query ',Dmel_p,' -out out_sim_p -outfmt 6'),intern=T)
#print("dsim p alignment done")

##Dere
#out_ere_n<-system(paste0(bl_dr,'blastn -db dere_db -query ',Dmel_n,' -out out_ere_n -outfmt 6'),intern=T)
#print("dere n alignment done")
#out_ere_p<-system(paste0(bl_dr,'tblastn -db dere_db -query ',Dmel_p,' -out out_ere_p -outfmt 6'),intern=T)
#print("dere p alignment done")

##Dmoj
#out_moj_n<-system(paste0(bl_dr,'blastn -db dmoj_db -query ',Dmel_n,' -out out_moj_n -outfmt 6'),intern=T)
#print("dmoj n alignment done")
#out_moj_p<-system(paste0(bl_dr,'tblastn -db dmoj_db -query ',Dmel_p,' -out out_moj_p -outfmt 6'),intern=T)
#print("dmoj p alignment done")

##Dsec
#out_sec_n<-system(paste0(bl_dr,'blastn -db dsec_db -query ',Dmel_n,' -out out_sec_n -outfmt 6'),intern=T)
#print("dsec n alignment done")
#out_sec_p<-system(paste0(bl_dr,'tblastn -db dsec_db -query ',Dmel_p,' -out out_sec_p -outfmt 6'),intern=T)
#print("dsec p alignment done")

##read outputs
require(seqinr)
names<-c('sim','moj','ere','sec')
mel_genes<-read.fasta(Dmel_n)#length=30828
mel_prot<-read.fasta(Dmel_p)#length=30493
genome_sec<-read.fasta(Dsec)
genome_sim<-read.fasta(Dsim)
genome_ere<-read.fasta(Dere)
genome_moj<-read.fasta(Dmoj)
colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
out_sec_n<-read.csv('out_sec_n',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_sec_p<-read.csv('out_sec_p',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_moj_n<-read.csv('out_moj_n',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_moj_p<-read.csv('out_moj_p',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_ere_n<-read.csv('out_ere_n',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_ere_p<-read.csv('out_ere_p',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_sim_n<-read.csv('out_sim_n',col.names=colnames,as.is=c(1,2),sep='\t',header=F)
out_sim_p<-read.csv('out_sim_p',col.names=colnames,as.is=c(1,2),sep='\t',header=F)

find_annotation_blast<-function(genome,out_complete){
  ret<-c()
  out<-unique(out_complete[,c(2,3,9,10)])
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
  for(i in 1:length(genome)){
    sum[i]<-sum(annotations[[i]])
  }
  ret$percent<-mean(out[which(out[,3]>80),3])
  ret$annotations<-annotations
  ret$genes<-unique(out[which(out[,3]>80),1])
  ret$sum<-sum
  return(ret)
}

r_sim_n<-find_annotation_blast(genome_sim,out_sim_n)
r_moj_n<-find_annotation_blast(genome_moj,out_moj_n)
r_ere_n<-find_annotation_blast(genome_ere,out_ere_n)
r_sec_n<-find_annotation_blast(genome_sec,out_sec_n)

r_sim_p<-find_annotation_blast(genome_sim,out_sim_p)
r_moj_p<-find_annotation_blast(genome_moj,out_moj_p)
r_ere_p<-find_annotation_blast(genome_ere,out_ere_p)
r_sec_p<-find_annotation_blast(genome_sec,out_sec_p)


