setwd('/home/sld83/GI/a2')
Dmel_n<-'dmel-all-gene-r6.24.fasta'
Dmel_p<-'dmel-all-translation-r6.24.fasta'
Dsim<-'dsim-all-gene-r2.02.fasta'
Dere<-'dere-all-gene-r1.05.fasta'
Dmoj<-'dmoj-all-gene-r1.04.fasta'
Dsec<-'dsec-all-gene-r1.3.fasta'

bl_dr<-'/local/data/public/genome_informatics_2018/programs/ncbi-blast-2.5.0+/bi
n/'

##make databases
#out_db_n<-system(paste0(bl_dr,'makeblastdb -in ',Dmel_n,' -dbtype nucl -out dme
l_db_n'),intern=T)
#out_db_p<-system(paste0(bl_dr,'makeblastdb -in ',Dmel_p,' -dbtype prot -out dme
l_db_p'),intern=T)

##align the genome onto the two databases
#out_sim_n<-system(paste0(bl_dr,'blastn -db dmel_db_n -query ',Dsim,' -out al_si
m_n -outfmt 6'),intern=T)
#out_sim_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dsim,' -out al_si
m_p -outfmt 6'),intern=T)
out_sim_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dsim,' -out al_sim
_p2 -outfmt 6 -evalue 1e-5 -max_target_seqs 2'),intern=T)
print("dsim alignment done")
#write(out_sim_n,file='out_sim_n')
#write(out_sim_p,file='out_sim_p')
#out_ere_n<-system(paste0(bl_dr,'blastn -db dmel_db_n -query ',Dere,' -out al_er
e_n -outfmt 6'),intern=T)
#out_ere_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dere,' -out al_er
e_p -outfmt 6 -evalue 1e-5'),intern=T)
out_ere_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dere,' -out al_ere
_p2 -outfmt 6 -evalue 1e-5 -max_target_seqs 2'),intern=T)
print("dere alignment done")
#write(out_ere_n,file='out_ere_n')
#write(out_ere_p,file='out_ere_p')
#out_moj_n<-system(paste0(bl_dr,'blastn -db dmel_db_n -query ',Dmoj,' -out al_mo
j_n -outfmt 6'),intern=T)
#out_moj_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dmoj,' -out al_mo
j_p -outfmt 6 -evalue 1e-5'),intern=T)
out_moj_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dmoj,' -out al_moj
_p2 -outfmt 6 -evalue 1e-5 -max_target_seqs 2'),intern=T)
print("dmoj alignment done")
#write(out_moj_n,file='out_moj_n')
#write(out_moj_p,file='out_moj_p')
#out_sec_n<-system(paste0(bl_dr,'blastn -db dmel_db_n -query ',Dsec,' -out al_se
c_n -outfmt 6'),intern=T)
#out_sec_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dsec,' -out al_se
c_p -outfmt 6 -evalue 1e-5'),intern=T)
out_sec_p<-system(paste0(bl_dr,'blastx -db dmel_db_p -query ',Dsec,' -out al_sec
_p2 -outfmt 6 -evalue 1e-5 -max_target_seqs 2'),intern=T)
print("dsec alignment done")
#write(out_sec_n,file='out_sec_n')
#write(out_sec_p,file='out_sec_p')