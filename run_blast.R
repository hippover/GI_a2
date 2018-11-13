setwd('/local/data/public/g1/')
Dmel_n<-'Drosophila_melanogaster.BDGP6.cdna.all.fa'
Dmel_p<-'Drosophila_melanogaster.BDGP6.pep.all.fa'
Dsim<-'dsim-all-chromosome-r2.02.fasta'
Dere<-'dere-all-gene-r1.05.fasta'
Dmoj<-'dmoj-all-gene-r1.04.fasta'
Dsec<-'dsec-all-gene-r1.3.fasta'

bl_dr<-'/local/data/public/genome_informatics_2018/programs/ncbi-blast-2.5.0+/bin/'

##make databases
out_db_sim<-system(paste0(bl_dr,'makeblastdb -in ',Dsim,' -dbtype nucl -out dsim_db'),intern=T)
out_db_ere<-system(paste0(bl_dr,'makeblastdb -in ',Dere,' -dbtype nucl -out dere_db'),intern=T)
out_db_moj<-system(paste0(bl_dr,'makeblastdb -in ',Dmoj,' -dbtype nucl -out dmij_db'),intern=T)
out_db_sec<-system(paste0(bl_dr,'makeblastdb -in ',Dsec,' -dbtype nucl -out dsec_db'),intern=T)
print("Databases created")
##align the genome onto the databases
##Dsim



##Dsim
out_sim_n<-system(paste0(bl_dr,'blastn -db dsim_db -query ',Dmel_n,' -out out_sim_n -outfmt 6'),intern=T)
print("dsim n alignment done")
out_sim_p<-system(paste0(bl_dr,'tblastn -db dsim_db -query ',Dmel_p,' -out out_sim_p -outfmt 6'),intern=T)
print("dsim p alignment done")

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