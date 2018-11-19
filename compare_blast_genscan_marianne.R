library(seqinr)
library(data.table)

## Function to have overlaps between gene regions
## takes as argument one blast output and one genscan output
## Note : query is melanogaster
genes.overlaps<-function(blastvect, genscanvect, cdnablastvect){
  
  sum_cdna<-sum(cdnablastvect$sum)
  sum_blast<-sum(blastvect$sum)
  sum_genscan<-sum(genscanvect$sum)
  sum_blast_cdna<-0
  sum_genscan_cdna<-0
  sum_blast_genscan<-0
  sum_blast_genscan_cdna<-0

  for(i in 1:length(names(cdnablastvect$annotations))){
    name<- names(cdnablastvect$annotations)[i]
    # sum_cdna<-sum_cdna +cdnablastvect$sum[i]
    # sum_blast<-sum_blast+blastvect$sum[i]
    # sum_genscan<-sum_genscan+genscanvect$sum[i]
    sum_blast_cdna<- sum_blast_cdna + length(intersect(which(cdnablastvect$annotations[[i]],T),which(blastvect$annotations[[i]],T)))
    sum_genscan_cdna<- sum_genscan_cdna + length(intersect(which(cdnablastvect$annotations[[i]],T),which(genscanvect$annotations[[i]],T)))
    sum_blast_genscan<- sum_blast_genscan + length(intersect(which(blastvect$annotations[[i]],T),which(genscanvect$annotations[[i]],T)))
    sum_blast_genscan_cdna<-sum_blast_genscan_cdna + length(intersect(intersect(which(blastvect$annotations[[i]],T),which(cdnablastvect$annotations[[i]],T)),which(genscanvect$annotations[[i]],T)))
    }
  
  
  
  return(data.frame(cdna=sum_cdna, blast = sum_blast, genscan= sum_genscan,overlap.blast=sum_blast_cdna,overlap.genscan=sum_genscan_cdna,
                  overlap.blast.genscan=sum_blast_genscan,overlap.blast.genscan.cdna=sum_blast_genscan_cdna,overlap.blast.percent=sum_blast_cdna/sum_cdna, overlap.genscan.percent=sum_genscan_cdna/sum_cdna,
                   overlap.blast.over.blast=sum_blast_cdna/sum_blast, overlap.genscan.over.genscan=sum_genscan_cdna/sum_genscan))

}

## Function that returns the vector of matching genes for blast
vector.match.blast<-function(specie,blastpath, genome){
  colnames<-c('qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'eval', 'bitscore')
  out_blast<-read.csv(blastpath,col.names=colnames,as.is=c(1,2),sep='\t',header=F)
  ret<-c()
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
  
  for(i in 1:length(genome)){
    sum[i]<-sum(annotations[[i]])
  }
  #print(annotations)
  #file_name<-paste(specie,'_blast_vector.RData',sep="")
  #save(annotations, file= file_name)
  #write.list(annotations,file_name)
  #return(annotations)
  ret$percent<-mean(out[which(out[,3]>80),3])
  ret$annotations<-annotations
  ret$genes<-unique(out[which(out[,3]>80),1])
  ret$sum<-sum
  return(ret)
}


## Function that returns the vector of matching genes for blast
vector.match.genscan2<-function(genscanpath, genome){
  #colnames<-c('gene_id','scf','tag','gene_start','gene_end','strand','start','end','score')
  out_blast<-read.csv(genscanpath,as.is=c(1,2),sep=',',header=T)
  out<-unique(out_blast[,c(1,2,3,4,5,6)])
  annotations <- vector(mode="list", length=length(genome))
  
  names(annotations)<-names(genome)
  sum<-rep(0,length(genome))
  ret<-c()
  for(i in 1:length(genome)){
    annotations[[i]]<-rep(F,length(genome[[i]]))
  }
  for (l in 1:nrow(out)){
    if(out[l,3] %in% names(annotations)){
      i<-which(names(annotations)==out[l,3])
      print(out[l,5])
      print(out[l,6])
      annotations[[i]][seq(out[l,5],out[l,6],1)]<-T
    }
  }
  sum<-0
  for(i in 1:length(annotations)){
    sum<-sum+sum(annotations[[i]])
  }
  #ret$percent<-mean(out[which(out[,3]>80),3])
  ret$annotations<-annotations
  #ret$genes<-unique(out[which(out[,3]>80),1])
  ret$sum<-sum
  print(ret$sum)
  
  return(ret)
}


# Dsim<-'/local/data/public/g1/chromosomes/dsim-all-chromosome-r2.01.fasta'
# Dere<-'/local/data/public/g1/chromosomes/dere-all-chromosome-r1.05.fasta'
# Dmoj<-'/local/data/public/g1/chromosomes/dmoj-all-chromosome-r1.04.fasta'
#Dsec<-'/local/data/public/g1/chromosomes/dsec-all-chromosome-r1.3.fasta'

# Dsim_gen<-'/local/data/public/g1/Genscan/Dsim_genes.csv'
# Dere_gen<-'/local/data/public/g1/Genscan/Dere_genes.csv'
#Dsec_gen<-'/local/data/public/g1/Genscan/Dsec_genes.csv'
#Dmoj_gen<-'/local/data/public/g1/Genscan/Dmoj_genes.csv'

# Dsim_genome<-read.fasta(Dsim)
# Dere_genome<-read.fasta(Dere)
# Dmoj_genome<-read.fasta(Dmoj)
#Dsec_genome<-read.fasta(Dsec)

# print('ere')
# genscan_ere_vector<-vector.match.genscan2(Dere_gen,Dere_genome)
# save(genscan_ere_vector,file='genscan_ere_vector')
# print('sim')
# genscan_sim_vector<-vector.match.genscan2(Dsim_gen,Dsim_genome)
# save(genscan_sim_vector,file='genscan_sim_vector')
# 
# print('sec')
# genscan_sec_vector<-vector.match.genscan2(Dsec_gen,Dsec_genome)
# save(genscan_sec_vector,file='genscan_sec_vector')

# print('moj')
# genscan_moj_vector<-vector.match.genscan2(Dmoj_gen,Dmoj_genome)
# save(genscan_moj_vector,file='genscan_moj_vector')
# 
# load('/local/data/public/g1/r_ere_n')
# load('/local/data/public/g1/r_sim_n')
# load('/local/data/public/g1/r_sec_n')
# load('/local/data/public/g1/r_moj_n')
# 
# load('/local/data/public/g1/r_ere_c')
# load('/local/data/public/g1/r_sim_c')
# load('/local/data/public/g1/r_sec_c')
# load('/local/data/public/g1/r_moj_c')
# 
# load('/local/data/public/g1/genscan_ere_vector')
# load('/local/data/public/g1/genscan_sim_vector')
# load('/local/data/public/g1/genscan_sec_vector')
# load('/local/data/public/g1/genscan_moj_vector')
# 
# #names(results)<-c('sim','ere','sec','moj')
# results<-list(data.frame(),data.frame(),data.frame(),data.frame())
# results[[1]]<- genes.overlaps(r_sim_n,genscan_sim_vector,r_sim_c)
# results[[2]]<- genes.overlaps(r_ere_n,genscan_ere_vector,r_ere_c)
# results[[3]]<- genes.overlaps(r_sec_n,genscan_sec_vector,r_sec_c)
# results[[4]]<- genes.overlaps(r_moj_n,genscan_moj_vector,r_moj_c)
# 
# print(results)
# 
# compare_blast_genscan<-results
# save(compare_blast_genscan,file='compare_blast_genscan')
# 
# 
# load('compare_blast_genscan')
# blast_percent<-NULL
# genscan_percent<-NULL
# for(i in 1:4){
#   blast_percent<-c(blast_percent,compare_blast_genscan[[i]]$overlap.blast.percent)
#   genscan_percent<-c(genscan_percent,compare_blast_genscan[[i]]$overlap.genscan.percent)
# }
# print(blast_percent)
# print(genscan_percent)
# pdf('compare_blast_genscan.pdf')
# mat<-t(matrix(c(blast_percent[1],genscan_percent[1],blast_percent[2],genscan_percent[2],blast_percent[3],genscan_percent[3],blast_percent[4],genscan_percent[4]),
#               ncol=2, byrow=TRUE,
#               dimnames=list(c('sim','ere','sec','moj'),c('Blast','Genscan'))))
# barplot(mat,
#           col=c('darkblue','brown'), space=rep(c(0.7,0.1),4),beside =TRUE, ylab= 'overlap percent with cDNA',
#         xlim=c(0,2*ncol(mat)+7),legend.text = TRUE,args.legend = list(x=2* ncol(mat)+7, y=max(mat),bty='n'))
# #lines(genscan_percent,col='blue')
# dev.off()


load('compare_blast_genscan')

#pdf('intersect.pdf')
#mat<-t(matrix(c(blast_percent[1],genscan_percent[1],blast_percent[2],genscan_percent[2],blast_percent[3],genscan_percent[3],blast_percent[4],genscan_percent[4]),
#              ncol=2, byrow=TRUE,
#              dimnames=list(c('sim','ere','sec','moj'),c('Blast','Genscan'))))
counts<- matrix(ncol=8,cbind(c(compare_blast_genscan[[1]]$blast,compare_blast_genscan[[1]]$genscan,compare_blast_genscan[[2]]$blast,compare_blast_genscan[[2]]$genscan,
                               compare_blast_genscan[[3]]$blast,compare_blast_genscan[[3]]$genscan,compare_blast_genscan[[4]]$blast,compare_blast_genscan[[4]]$genscan),
                             
                             c(compare_blast_genscan[[1]]$overlap.blast,compare_blast_genscan[[1]]$overlap.genscan,compare_blast_genscan[[2]]$overlap.blast,compare_blast_genscan[[2]]$overlap.genscan,
                               compare_blast_genscan[[3]]$overlap.blast,compare_blast_genscan[[3]]$overlap.genscan,compare_blast_genscan[[4]]$overlap.blast,compare_blast_genscan[[4]]$overlap.genscan),
                             
               c(compare_blast_genscan[[1]]$overlap.blast.genscan.cdna, compare_blast_genscan[[1]]$overlap.blast.genscan.cdna,compare_blast_genscan[[2]]$overlap.blast.genscan.cdna, compare_blast_genscan[[2]]$overlap.blast.genscan.cdna,
                 compare_blast_genscan[[3]]$overlap.blast.genscan.cdna, compare_blast_genscan[[3]]$overlap.blast.genscan.cdna,compare_blast_genscan[[4]]$overlap.blast.genscan.cdna, compare_blast_genscan[[4]]$overlap.blast.genscan.cdna)),
               dimnames=list(c('0','1','2'),c('sim blast','sim genscan', 'ere blast', 'ere genscan','sec blast', 'sec genscan','moj blast', 'moj genscan')), byrow=TRUE)


print(counts)

# counts1<- matrix(ncol=4,cbind(c(compare_blast_genscan[[1]]$blast,compare_blast_genscan[[2]]$blast,
#                                compare_blast_genscan[[3]]$blast,compare_blast_genscan[[4]]$blast),
#                              c(compare_blast_genscan[[1]]$overlap.blast,compare_blast_genscan[[2]]$overlap.blast,
#                                compare_blast_genscan[[3]]$overlap.blast,compare_blast_genscan[[4]]$overlap.blast),
#                              c(compare_blast_genscan[[1]]$overlap.blast.genscan.cdna,compare_blast_genscan[[2]]$overlap.blast.genscan.cdna,
#                                compare_blast_genscan[[3]]$overlap.blast.genscan.cdna, compare_blast_genscan[[4]]$overlap.blast.genscan.cdna)),
#                 dimnames=list(c('0','1','2'),c('sim1','ere1','sec1','moj1')), byrow=TRUE)
# 
# counts2<- matrix(ncol=4,cbind(c(compare_blast_genscan[[1]]$genscan,compare_blast_genscan[[2]]$genscan,
#                                compare_blast_genscan[[3]]$genscan,compare_blast_genscan[[4]]$genscan),
#                              c(compare_blast_genscan[[1]]$overlap.genscan,compare_blast_genscan[[2]]$overlap.genscan,
#                                compare_blast_genscan[[3]]$overlap.genscan,compare_blast_genscan[[4]]$overlap.genscan),
#                              c(compare_blast_genscan[[1]]$overlap.blast.genscan.cdna,compare_blast_genscan[[2]]$overlap.blast.genscan.cdna,
#                                compare_blast_genscan[[3]]$overlap.blast.genscan.cdna,compare_blast_genscan[[4]]$overlap.blast.genscan.cdna)),
#                 dimnames=list(c('0','1','2'),c('sim','ere','sec','moj')), byrow=TRUE)

#print(counts1)
#print(counts2)
barplot(counts, col=c('darksalmon','darkseagreen1','cyan4','red','green','blue'), horiz=TRUE, 
        space=c(0.1,0.1,0.7,0.1,0.7,0.1,0.7,0.1,0.1))
#barplot(cbind(counts1[,1], counts2[,1],counts1[,2], counts2[,2],counts1[,3], counts2[,3],counts1[,4], counts2[,4]))

#lines(genscan_percent,col='blue')
#dev.off()
