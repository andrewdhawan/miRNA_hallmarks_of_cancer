#linear model of miRNAs predicting various signature scores


library(RankProd)
library(reshape2)
library(penalized)
cancer_types_list <- list();
cancer_types_list[[1]] <- c('BRCA','UCEC','HNSC')
cancer_types_list[[2]] <- c('KIRC','LUAD','THCA')
cancer_types_list[[3]] <- c('PRAD','LUSC','OV')
cancer_types_list[[4]] <- c('STAD','BLCA','COAD')
cancer_types_list[[5]] <- c('LIHC','CESC','KIRP')

all_cancer_types <- melt(cancer_types_list)$value

#load the signatures

sig_fnames_list <- list();
sig_names_list <- list();
categories_of_sigs <- c('invasion','energetics','immortality','growth_suppressors','genome_instability','angiogenesis','apoptosis','proliferation','inflammation')

sig_fnames_list[['invasion']] <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt','invasiveness_gene_sig_entrez_marsan2014.txt')
sig_names_list[['invasion']] <- c('Hallmark: Epithelial Mesenchymal Transition','Invasiveness, Marsan 2014')

sig_fnames_list[['energetics']] <- c('HALLMARK_OXIDATIVE_PHOSPHORYLATION.txt','HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY.txt')
sig_names_list[['energetics']] <- c('Hallmark: Oxidative Phosphorylation','Hallmark: Reactive Oxygen Species Pathway')

sig_fnames_list[['immortality']] <- c('HALLMARK_G2M_CHECKPOINT.txt')
sig_names_list[['immortality']] <- c('Hallmark: G2M Checkpoint')

sig_fnames_list[['growth_suppressors']] <- c('HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt','HALLMARK_XENOBIOTIC_METABOLISM.txt')
sig_names_list[['growth_suppressors']] <- c('Hallmark: PI3K AKT MTOR Signaling','Hallmark: Xenobiotic Metabolism')

sig_fnames_list[['genome_instability']] <- c('HALLMARK_DNA_REPAIR.txt','HALLMARK_P53_PATHWAY.txt')
sig_names_list[['genome_instability']] <- c('Hallmark: DNA Repair','Hallmark: p53 Pathway')

sig_fnames_list[['angiogenesis']] <- c('hypoxia_gene_sig_entrez_probes.txt','HALLMARK_ANGIOGENESIS.txt','HALLMARK_HYPOXIA.txt','angiogenesis_gene_sig_entrez_desmedt2008_pos.txt','Masiero2013angiogenesisENTREZ.txt')
sig_names_list[['angiogenesis']] <- c('Hypoxia, Buffa 2010','Hallmark: Angiogenesis','Hallmark: Hypoxia','Angiogenesis, Desmedt 2008','Angiogenesis, Masiero 2013')

sig_fnames_list[['apoptosis']] <- c('HALLMARK_APOPTOSIS.txt','apoptosis_gene_sig_entrez_desmedt2008_pos.txt')
sig_names_list[['apoptosis']] <- c('Hallmark: Apoptosis','Apoptosis, Desmedt 2008')

sig_fnames_list[['proliferation']] <-c('proliferation_gene_sig_entrez_desmedt2008_pos.txt','HALLMARK_KRAS_SIGNALING_UP.txt')
sig_names_list[['proliferation']] <-c('Proliferation, Desmedt 2008','Hallmark: KRAS Signaling Up')

sig_fnames_list[['inflammation']] <- c('HALLMARK_INFLAMMATORY_RESPONSE.txt','HALLMARK_IL2_STAT5_SIGNALING.txt','HALLMARK_IL6_JAK_STAT3_SIGNALING.txt','HALLMARK_TGF_BETA_SIGNALING.txt','HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt','immune_gene_sig_entrez_desmedt2008_pos.txt')
sig_names_list[['inflammation']] <- c('Hallmark: Inflammatory Response','Hallmark: IL2 STAT5 Signaling','Hallmark: IL6 JAK STAT3 Signaling','Hallmark: TGF Beta Signaling','Hallmark: TNFa Signaling via NFKB','Immune, Desmedt 2008')

sigs_list_by_cat <- list();
for(sig_category in categories_of_sigs){
	sigs_list_by_cat[[sig_category]] <- list();
	for(i in 1:length(sig_fnames_list[[sig_category]])){
		fname <- sig_fnames_list[[sig_category]][i]
		genes = read.csv(paste0('gene_signatures/',fname), header=F, stringsAsFactors=F, colClasses = "character")
		# print(genes)
		sigs_list_by_cat[[sig_category]][[sig_names_list[[sig_category]][i]]]<- genes
	} 
}

#load the datasets

all_mRNA_datasets <- list();
for (cancer_type in all_cancer_types){
	print(cancer_type)
	if(cancer_type!='BRCA'){
		fname_mrna <- paste0('../Reprocessed GDAC data/',cancer_type,'/mRNA/tumour/cleaned_mRNA.txt')
	}else{
		fname_mrna <- paste0('../Reprocessed GDAC data/',cancer_type,'/mRNA/tumour/cleaned_mRNA_ductal.txt')
	}
	all_mRNA_datasets[[cancer_type]] <- read.table(fname_mrna, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(all_mRNA_datasets[[cancer_type]]) <- gsub('[.]','-',colnames(all_mRNA_datasets[[cancer_type]]))
	# want log2 data
	all_mRNA_datasets[[cancer_type]] <- log2(all_mRNA_datasets[[cancer_type]]+1)
	all_mRNA_datasets[[cancer_type]][!is.finite(as.matrix(all_mRNA_datasets[[cancer_type]]))] <- NA
}


#load the miRNA

all_miRNA_datasets <- list();
all_miRNA <- c()
for (cancer_type in all_cancer_types){
	fname_miRNA <- paste0('../Reprocessed GDAC data/',cancer_type,'/miRNA/tumour/cleaned_miRNA_mature.txt')
	all_miRNA_datasets[[cancer_type]] <- read.table(fname_miRNA, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(all_miRNA_datasets[[cancer_type]]) <- gsub('[.]','-',colnames(all_miRNA_datasets[[cancer_type]]))
	all_miRNA <- unique(c(all_miRNA,rownames(all_miRNA_datasets[[cancer_type]])))
}

all_coeffs <- list();
all_rank_product_matrices <- list();
for (category in categories_of_sigs){
	count <- 1
	for (gene_sig in sigs_list_by_cat[[category]]){

		sig_name <- sig_names_list[[category]][count]
		print(sig_name)

		all_coeffs[[sig_name]] <- matrix(0,nrow=length(all_miRNA),ncol=length(all_cancer_types))
		row.names(all_coeffs[[sig_name]]) <- all_miRNA
		colnames(all_coeffs[[sig_name]]) <- all_cancer_types
 		for (cancer_type in all_cancer_types){
 			print(cancer_type)
			genes_present <- intersect(rownames(all_mRNA_datasets[[cancer_type]]),gene_sig$V1)
			#compute and score the scores

			scores <- apply(all_mRNA_datasets[[cancer_type]][genes_present,], 2, function(x) median(x,na.rm=T))

			#cross-validated linear model
			coeffs <- get_coefficients_pre_filter(cancer_type,scores)
			#store the miRNA results
			all_coeffs[[sig_name]][names(coeffs),cancer_type] <- coeffs
		}

		# #compute the rank-product matrix
		# all_rank_product_matrices[[sig_name]] <- make_rank_prod_matrix(all_coeffs[[sig_name]])
		count <- count + 1
	}
}

#the above code is run on a server for each signature in parallel, and the files are saved into a folder
#called 'server_data.' Using this, we re-load everythng in R and then compute the overall rank prod matrices

#the following is the code to load in from all signatures the files from the code running on server
all_signatures <- melt(sig_names_list)$value
rank_prod_tables <- list();
RP_out_values <- list();
all_coeffs_tmp <- list();
for (sig_name in all_signatures){
	load(paste0('server_data/all_coeffs_',sig_name,'.rda'))
	all_coeffs_tmp[[sig_name]] <- all_coeffs[[sig_name]]
}
all_coeffs <- all_coeffs_tmp


#library(rankProd)
rank_prod_tables <- list();
RP_out_values <- list();
# for (category in categories_of_sigs){
# 	count <- 1
# 	for (gene_sig in sigs_list_by_cat[[category]]){
# 		sig_name <- sig_names_list[[category]][count]
#here we need to do the rankprod
for (sig_name in all_signatures){
	all_coeffs[[sig_name]] <- all_coeffs[[sig_name]][which(rowSums(all_coeffs[[sig_name]]==0) < length(colnames(all_coeffs[[sig_name]]))),]
	print(dim(all_coeffs[[sig_name]]))
	RP.out <- RP(all_coeffs[[sig_name]],rep(1,15))
	RP_out_values[[sig_name]] <- RP.out
	rank_prod_tables[[sig_name]] <- topGene(RP.out,cutoff = 0.05,method="pfp",gene.names=rownames(all_coeffs[[sig_name]]))
		# count <- count + 1
}
# }

#then save the outputs

save(file='rank_prod_output_pre_filtered.rda',RP_out_values)
save(file='rank_prod_tables_out_pre_filtered.rda',rank_prod_tables)



for (sig_name in all_signatures){
#for each sig let's save the heatmap of the miRNA coefficients to see whether the cancers act the same
	all_coeffs_tmp_mod <- all_coeffs[[sig_name]][which(rowSums(all_coeffs[[sig_name]]!=0)!=0),]

	gplots::heatmap.2( all_coeffs_tmp_mod,
	                   col = gplots::colorpanel(100,"blue","white","red"),#gplots::colorpanel(100,"white","red"),#gplots::redgreen(100),#gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
	                   trace = "none",
	                   xlab = "Gene ID",
	                   ylab="Gene ID",
	                   na.color="grey",
	                   #labRow=rownames(autocors),
	                   #labCol=colnames(autocors),#gene_sig,
	                   main = paste0("\n\n", sig_name),
	                   dendrogram = "both",
	                   #symbreaks = T,
	                   Rowv = T,Colv=T ,key.xlab='Rho',key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.15,cexCol=0.45)

	dev.copy(pdf,paste0('miRNA_hmap_preFiltered_',sig_name,'.pdf'),width=12,height=12)
  dev.off()
  

}



# #----------the following is to make a  heatmap but for the miRNA that recur among cancer types for each signature themselves, not the families:

# all_sigs_miRNA_list <- list()
# for (sig_name in all_signatures){
# 	cur_miRNAs_list <- c()
# 	for (cancer_type in all_cancer_types){
# 		cur_miRNAs_list <- c(cur_miRNAs_list,rownames(all_coeffs[[sig_name]])[which(all_coeffs[[sig_name]][,cancer_type] < 0)])
# 	}
# 	print(table(cur_miRNAs_list))
# 	all_sigs_miRNA_list[[sig_name]] <- table(cur_miRNAs_list)
# } #counts frequency of the miRNA occurring across all cancer types as significant
# heatmap_matrix <- matrix(0, nrow=length(unique(melt(all_sigs_miRNA_list)[,1])),ncol=length(all_signatures))
# row.names(heatmap_matrix) <- unique(melt(all_sigs_miRNA_list)[,1])
# colnames(heatmap_matrix) <- all_signatures
# for (sig_name in all_signatures){
# 	heatmap_matrix[names(all_sigs_miRNA_list[[sig_name]]),sig_name] <- as.numeric(all_sigs_miRNA_list[[sig_name]])
# }


# gplots::heatmap.2( heatmap_matrix,
# 	                   col = gplots::colorpanel(100,"white","red"),#gplots::colorpanel(100,"white","red"),#gplots::redgreen(100),#gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
# 	                   trace = "none",
# 	                #   xlab = "Gene ID",
# 	                #   ylab="Gene ID",
# 	                   na.color="grey",
# 	                   #labRow=rownames(autocors),
# 	                   #labCol=colnames(autocors),#gene_sig,
# 	                   main = paste0("\n", "miRNA down freq \nof occurrence"),
# 	                   dendrogram = "both",
# 	                   #symbreaks = T,
# 	                   Rowv = T,Colv=T ,key.xlab='Rho',key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.11,cexCol=0.35)


# 	dev.copy(pdf,paste0('miRNA_freq_preFiltered_DOWN_all_sigs.pdf'),width=12,height=12)
#   dev.off()

# #--------------------------------------------------------------------------------------------------

get_coefficients <- function(cancer_type,scores){

	#load in the miRNA data
	#fname_miRNA <- paste0('../Reprocessed GDAC data/',cancer_type,'/miRNA/tumour/cleaned_miRNA_mature_log2.txt')
	miRNA_data <- all_miRNA_datasets[[cancer_type]]#read.table(miRNA_fName, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	#colnames(miRNA_data) <- gsub('[.]','-',colnames(miRNA_data))

	#take only common subset of miRNA and scores
	common_colNames <- intersect(colnames(miRNA_data),names(scores))

	#take just the common pieces
	miRNA_data <- miRNA_data[,common_colNames]
	scores <- scores[common_colNames]

	#z-transform the scores
	scores <- as.numeric(scores) - mean(as.numeric(scores))/sd(as.numeric(scores))
	print(sum(is.na(scores)))

	#expression filter for miRNA
	expression_threshold <- 0.80 # means that at least 10% of samples must have a nonzero value of the mRNA
	miRNA_data <-miRNA_data[which((rowSums(miRNA_data==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]

	#remove NA values from miRNA data
#	expression_threshold <- 0.5
#	miRNA_data <-miRNA_data[which((rowSums(is.na(miRNA_data)==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]
# miRNA_data <-miRNA_data[which(rowSums(is.na(miRNA_data))==0),]
	miRNA_data <- as.matrix(log2(miRNA_data))
	miRNA_data[!(is.finite(miRNA_data))] <- NA
	#z-transform the miRNA data
	for (j in 1:length(rownames(miRNA_data))){
		miRNA_data[j,] <- (as.numeric(miRNA_data[j,]) - mean(as.numeric(miRNA_data[j,])))/sd(as.numeric(miRNA_data[j,]))
	}
	print(paste0("mirna " , sum(is.na(miRNA_data))))

	#penalised linear regression
	new_df <- na.omit(t(rbind(scores,miRNA_data)))
	colnames(new_df) <- c('scores',rownames(miRNA_data))
	print(new_df[1:4,1:4])
	lambda_2_values <- c(0, 0.01, 0.1,1,10,100)
	max_likelihood <- -9999999999
	for (lambda2_val in lambda_2_values){
		cross_val_model <- optL1(response = new_df[,1],penalized = new_df[,2:length(colnames(new_df))], lambda2 = lambda2_val,data=as.data.frame(new_df),model="linear",fold=10,trace=F)#,trace=F,maxiter=1000,tol=.Machine$double.eps^0.23)
		#	cross_val_model <- optL2(response = all_sig_scores[,1],penalized = all_sig_scores[,2:length(colnames(all_sig_scores))], minlambda2 = 0,maxlambda2=100,data=all_sig_scores[,2:length(colnames(all_sig_scores))],model="linear",fold=10)#lambda2 = lambda2_val,data=all_sig_scores,model="linear",fold=10)
		if ((cross_val_model$fullfit)@loglik > max_likelihood){
			best_model <<- cross_val_model
			best_lambda <- lambda2_val
		}
	}

	miRNA_names_reported <- intersect(names(coef(best_model$fullfit)), rownames(miRNA_data))
	#best_coef_matrix[rownames(miRNA_matrix)[i],mRNA_names_reported] <- coef(best_model$fullfit)[mRNA_names_reported]

	#return the coefficients
	coef(best_model$fullfit)[miRNA_names_reported]
}





get_coefficients_pre_filter <- function(cancer_type,scores){

	#load in the miRNA data
	#fname_miRNA <- paste0('../Reprocessed GDAC data/',cancer_type,'/miRNA/tumour/cleaned_miRNA_mature_log2.txt')
	miRNA_data <- all_miRNA_datasets[[cancer_type]]#read.table(miRNA_fName, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	#colnames(miRNA_data) <- gsub('[.]','-',colnames(miRNA_data))

	#take only common subset of miRNA and scores
	common_colNames <- intersect(colnames(miRNA_data),names(scores))

	#take just the common pieces
	miRNA_data <- miRNA_data[,common_colNames]
	scores <- scores[common_colNames]

	#z-transform the scores
	scores <- as.numeric(scores) - mean(as.numeric(scores))/sd(as.numeric(scores))
	print(sum(is.na(scores)))

	#expression filter for miRNA
	expression_threshold <- 0.80 # means that at least 10% of samples must have a nonzero value of the mRNA
	miRNA_data <-miRNA_data[which((rowSums(miRNA_data==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]

	#remove NA values from miRNA data
#	expression_threshold <- 0.5
#	miRNA_data <-miRNA_data[which((rowSums(is.na(miRNA_data)==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]
# miRNA_data <-miRNA_data[which(rowSums(is.na(miRNA_data))==0),]
	miRNA_data <- as.matrix(log2(miRNA_data))
	miRNA_data[!(is.finite(miRNA_data))] <- NA
	#z-transform the miRNA data
	for (j in 1:length(rownames(miRNA_data))){
		miRNA_data[j,] <- (as.numeric(miRNA_data[j,]) - mean(as.numeric(miRNA_data[j,])))/sd(as.numeric(miRNA_data[j,]))
	}
	print(paste0("mirna " , sum(is.na(miRNA_data))))



	#first we need to subset the data into folds


	new_df <- na.omit(t(rbind(scores,miRNA_data)))
	colnames(new_df) <- c('scores',rownames(miRNA_data))

	folds <- 10
	nrows_combined_df <- 1:dim(new_df)[1]
	best_overall_error <- 99999999
	for (i in 0:(folds-1)){
		new_df_subset <- as.data.frame(new_df[!(nrows_combined_df%%folds==i),]) #takes out the 1/nth row of the data set
		#train the univaraite model
		#put these as inputs to the penalized model

		linear_models_miRNA <- matrix(,nrow=length(rownames(miRNA_data)),ncol=1)
		row.names(linear_models_miRNA) <- rownames(miRNA_data)
		for (j in 1:length(rownames(miRNA_data))){
			univariate_data <- as.data.frame(cbind(new_df_subset[,1],new_df_subset[,(j+1)]))
			colnames(univariate_data) <- c('sig_score','miRNA')
			# print(univariate_data)

			univariate_model <- lm(formula = sig_score ~ miRNA,data = univariate_data)
			# print(summary(univariate_model))
	 		linear_models_miRNA[j] <- (summary(univariate_model)$coefficients)[2,4]
			#tmp_model <- coef(summary(coxph(Surv(as.numeric(combined_df_subsetted$times),as.numeric(combined_df_subsetted$events)) ~ combined_df_subsetted[,j+2])))
			#cox_models_circRNA[j] <- tmp_model[5] #c(tmp_model[2],tmp_model[5])
		}

		#significant miRNAs are those w p < 0.2:
		significant_miRNAs <- rownames(linear_models_miRNA)[which(linear_models_miRNA < 0.2 & !is.nan(linear_models_miRNA))]
		# print("sig MiRNA")
		# print(significant_miRNAs)
		#penalised linear regression
		# print(new_df_subset[1:4,1:4])
		lambda_2_values <- c(0, 0.01, 0.1,1,10,100)
		max_likelihood <- -9999999999
		for (lambda2_val in lambda_2_values){

			cross_val_model <- optL1(response = new_df_subset[,1],penalized = new_df_subset[,significant_miRNAs], lambda2 = lambda2_val,data=as.data.frame(new_df_subset),model="linear",fold=10,trace=F)#,trace=F,maxiter=1000,tol=.Machine$double.eps^0.23)
			#	cross_val_model <- optL2(response = all_sig_scores[,1],penalized = all_sig_scores[,2:length(colnames(all_sig_scores))], minlambda2 = 0,maxlambda2=100,data=all_sig_scores[,2:length(colnames(all_sig_scores))],model="linear",fold=10)#lambda2 = lambda2_val,data=all_sig_scores,model="linear",fold=10)
			if ((cross_val_model$fullfit)@loglik > max_likelihood){
				best_model <<- cross_val_model
				best_lambda <- lambda2_val
			}
		}


		#now that we know the best model, let's test it on the other 1/n of the data, and record the error
		unused_df <- as.data.frame(new_df[(nrows_combined_df%%folds==i),])
		current_predictions <- predict(best_model$fullfit, penalized=unused_df[,significant_miRNAs],data=unused_df)


		cur_error <- norm((as.numeric(unused_df[,1]) - as.numeric(current_predictions)),type="2")
		# print(cur_error)
		if (cur_error < best_overall_error){
			best_overall_error <- cur_error
			best_overall_model <- best_model
			best_overall_lambda <- best_lambda
		}
	}

	miRNA_names_reported <- intersect(names(coef(best_overall_model$fullfit)), rownames(miRNA_data))
	#best_coef_matrix[rownames(miRNA_matrix)[i],mRNA_names_reported] <- coef(best_model$fullfit)[mRNA_names_reported]

	#return the coefficients
	coef(best_overall_model$fullfit)[miRNA_names_reported]
}



