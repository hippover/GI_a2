library(seqinr)

Dmel_n<-'dmel-all-gene-r6.24.fasta'
Dmel_p<-'dmel-all-translation-r6.24.fasta'

Dsim<-'dsim-all-gene-r2.02.fasta'
Dere<-'dere-all-gene-r1.05.fasta'
Dmoj<-'dmoj-all-gene-r1.04.fasta'
Dsec<-'dsec-all-gene-r1.3.fasta'

genes<-c(Dsim,Dere,Dmoj,Dsec)

genescan.path<-'/local/data/public/genome_informatics_2018/programs/genscan/'
home.path<-'/home/md799/GI/A2/'

for(g in genes){
  outfile<-paste(paste('genscan_',substr(g,1,4),sep=""),'.txt',sep="")
  write('',outfile,append=FALSE)
  sequences<-read.fasta(g)
  for(s in sequences){
    write.fasta(getSequence(s), attr(s,"name"), file.out='temp.fa',open="w",as.string=FALSE)
    outex<-system(paste(paste(paste(paste(genescan.path,'genscan',sep=""),
      paste(genescan.path,'/lib/HumanIso.smat',sep="")), paste(home.path,'temp.fa',sep="") ),'-cds'),intern=TRUE)
    
    write(outex,outfile, append=TRUE)
    write("\t\t\t\t\t",outfile, append=TRUE)
  }
}
