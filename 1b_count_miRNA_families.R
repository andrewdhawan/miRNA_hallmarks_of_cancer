#This is the code to check the overlap of the signficant miRNA and the miRNA families they come from

load('rank_prod_tables_out_pre_filtered.rda')
library(targetscan.Hs.eg.db)
library(reshape2)
list_of_all_final_tables <- list();
sig_name <- 'Hallmark: Angiogenesis'#'Hallmark: TNFa Signaling via NFKB'#'Hallmark: p53 Pathway'#'Hallmark: TGF Beta Signaling'#'Hallmark: IL6 JAK STAT3 Signaling'#'Hallmark: G2M Checkpoint'#'Hallmark: Inflammatory Response'#'Hallmark: Angiogenesis'#'Hallmark: KRAS Signaling Up'#'Hypoxia, Buffa 2010'#'Hallmark: Epithelial Mesenchymal Transition'#'Hallmark: Apoptosis'#'Invasiveness, Marsan 2014'#'Hallmark: DNA Repair'
all_sigs <- c('Hallmark: Epithelial Mesenchymal Transition','Invasiveness, Marsan 2014',
	'Hallmark: Oxidative Phosphorylation','Hallmark: Reactive Oxygen Species Pathway',
	'Hallmark: G2M Checkpoint',
	'Hallmark: PI3K AKT MTOR Signaling','Hallmark: Xenobiotic Metabolism',
	'Hallmark: DNA Repair','Hallmark: p53 Pathway',
	'Hypoxia, Buffa 2010','Hallmark: Angiogenesis','Hallmark: Hypoxia','Angiogenesis, Desmedt 2008','Angiogenesis, Masiero 2013',
	'Hallmark: Apoptosis','Apoptosis, Desmedt 2008',
	'Proliferation, Desmedt 2008','Hallmark: KRAS Signaling Up',
	'Hallmark: Inflammatory Response','Hallmark: IL2 STAT5 Signaling','Hallmark: IL6 JAK STAT3 Signaling','Hallmark: TGF Beta Signaling','Hallmark: TNFa Signaling via NFKB','Immune, Desmedt 2008')



for(sig_name in all_sigs){
check_fam <- rownames(rank_prod_tables[[sig_name]]$Table1)
ans_list <- c()
for (i in 1:length(check_fam)){
    tryCatch(ans_list<<-c(ans_list,mget(check_fam[i],targetscan.Hs.egMIRBASE2FAMILY)),
    	error=function(e){
    		tryCatch({
    			print(gsub('-3p','',check_fam[i]))
				ans_list<<-c(ans_list,mget(gsub('-3p','',check_fam[i]),targetscan.Hs.egMIRBASE2FAMILY))}, error=function(e1){

					tryCatch(
						{print(gsub('-5p','',check_fam[i]))
						ans_list<<-c(ans_list,mget(gsub('-5p','',check_fam[i]),targetscan.Hs.egMIRBASE2FAMILY))},error=function(e2){
							
							tryCatch({print(gsub('[-]*[0-9]*-3p','',check_fam[i]))
								ans_list<<-c(ans_list,mget(gsub('[-]*[0-9]*-3p','',check_fam[i]),targetscan.Hs.egMIRBASE2FAMILY))},error=function(e3){

								tryCatch(
									{print(gsub('[-]*[0-9]*-5p','',check_fam[i]))
									ans_list<<-c(ans_list,mget(gsub('[-]*[0-9]*-5p','',check_fam[i]),targetscan.Hs.egMIRBASE2FAMILY))},error=function(e4){

									tryCatch(
									{print(gsub('[a-z]*[-]*[0-9]*-3p','',check_fam[i]))
									ans_list<<-c(ans_list,mget(gsub('[a-z]*[-]*[0-9]*-3p','',check_fam[i]),targetscan.Hs.egMIRBASE2FAMILY))},error=function(e5){

										tryCatch(
									{print(gsub('[a-z]*[-]*[0-9]*-5p','',check_fam[i]))
									ans_list<<-c(ans_list,mget(gsub('[a-z]*[-]*[0-9]*-5p','',check_fam[i]),targetscan.Hs.egMIRBASE2FAMILY))},error=function(e6){
											tryCatch(
									{print(substring(text = check_fam[i],first = 1,last = (nchar(check_fam[i])-1)))#gsub('[a-z]*[-]*[0-9]*-5p','',check_fam[i]))
									ans_list<<-c(ans_list,mget(substring(text = check_fam[i],first = 1,last = (nchar(check_fam[i])-1)),targetscan.Hs.egMIRBASE2FAMILY))},error=function(e7){
										print(e7)
								})
							})
						})
					})
				})
			})
		})
    })
}




ans_mat <- as.matrix(ans_list)
frequency_mat <- matrix(0,nrow=length(unique(ans_mat)),ncol=1)

row.names(frequency_mat) <- unique(ans_mat)

for (i in 1:length(rownames(ans_mat))){
	frequency_mat[ans_mat[[i,1]],1] <- frequency_mat[ans_mat[[i,1]],1] + 1
}
frequency_mat <- frequency_mat[order(-frequency_mat),]
print(frequency_mat)

length_of_list <- length(rownames(ans_mat))
human_mirs_mapped <- ls(targetscan.Hs.egMIRBASE2FAMILY)[which(grepl('hsa',ls(targetscan.Hs.egMIRBASE2FAMILY))==T)]
N_reps <- 1000

frequency_mat_resample <- matrix(0,nrow=length(names(frequency_mat)),ncol=N_reps)
row.names(frequency_mat_resample) <- names(frequency_mat)
for (i in 1:N_reps){
	temp <- sample(human_mirs_mapped,length_of_list,replace=F)
	resampled_fams <-mget(temp,targetscan.Hs.egMIRBASE2FAMILY)
	resampled_fams <- as.matrix(resampled_fams)[,1]
	num_in_common <- intersect(rownames(frequency_mat_resample),resampled_fams)
	if(length(num_in_common) > 0){
		for(j in 1:length(num_in_common)){
			frequency_mat_resample[num_in_common[[j]],i] <- length(which(resampled_fams==num_in_common[[j]]))
		}
	}
}
p_values <- c()
for(i in 1:length(names(frequency_mat))){
	p_values<- c(p_values, ((sum(frequency_mat_resample[i,]>frequency_mat[i]))/N_reps))
}
final_answer_mat <- cbind(frequency_mat,p_values)
final_answer_mat <- final_answer_mat[order(final_answer_mat[,'p_values']),]
final_answer_mat <- final_answer_mat[final_answer_mat[,'p_values'] < (0.05/length(final_answer_mat[,'p_values'])),]
list_of_all_final_tables[[sig_name]] <- final_answer_mat
}
#print(final_answer_mat)

all_fam_names <- c()
for (sig_name in all_sigs){
all_fam_names <- c(all_fam_names, rownames(list_of_all_final_tables[[sig_name]]))
all_fam_names <- unique(all_fam_names)

}
heatmap_final_tables <- matrix(0,nrow=length(all_fam_names),ncol=length(all_sigs))
row.names(heatmap_final_tables) <- all_fam_names
colnames(heatmap_final_tables) <- all_sigs
count <-1
for (sig_name in all_sigs){
	heatmap_final_tables[rownames(list_of_all_final_tables[[sig_name]]),count] <- as.matrix(list_of_all_final_tables[[sig_name]][,'frequency_mat'])
	count <- count+1
}

gplots::heatmap.2( heatmap_final_tables,
                   col = gplots::colorpanel(100,"white","red"),#gplots::redgreen(100),#gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
                   trace = "none",
                   xlab = "Gene ID",
                   ylab="Gene ID",
                   na.color="grey",
                   #labRow=rownames(autocors),
                   #labCol=colnames(autocors),#gene_sig,
                   #main = paste0("\n\nAutocorrelation\n", names_datasets[[dataset_ind]] ,' ',names_sigs[[sig_ind]]),
                   dendrogram = "both",
                   #symbreaks = T,
                   Rowv = T,Colv=T ,key.xlab='Rho',key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.25,cexCol=0.35)


