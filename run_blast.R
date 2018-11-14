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
out_sim_n<-system(paste0(bl_dr,'blastn -db dsim_db -query ',Dmel_n,' -out out_sim_n -outfmt 6'),intern=T)
print("dsim n alignment done")
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

find_annotation_blast<-function(genome_file,blast_file){
  ret<-c()
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out<-read.csv(blast_file,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  genome<-read.fasta(genome_file)
  annotations <- vector(mode="list", length=length(genome))
  names(annotations)<-names(genome)
  for(i in 1:length(names(sim_genome))){
    annotations[[i]]<-rep(0,length(genome[[i]]))
  }
  for (l in 1:length(out)){
    if(out[l,2] %in% names(annotations) & out[l,3]){
      i<-which(names(annotations)==out[l,2])
      I<-c(1:length(annotations[[i]]))
      annotations[[i]][which(I>=out[l,9] & I<=out[l,10])]<-1
    }
  }
  ret$annotations<-annotations
  ret$genes<-unique(out[,1])
  return(ret)
}
annotations_sim_n<-find_annotation_blast(Dsim,'out_sim_n')

percent<-vector(mode="list",length=length(names(sim_genome_file)))
names(percent)<-names(annotations_sim_n)
for (i in 1:length(annotations_sim_n)){
  percent[[i]]<-mean(annotations_sim_n[[i]])
}
