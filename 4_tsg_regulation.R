# this code studies the regulation of the TSG

library(org.Hs.eg.db)
library(miRNAtap)
library(ppcor)
library(RankProd)
library(gplots)
library(reshape2)
library(penalized)
library(ggplot2)

#must have all data loaded
cancer_types <- c('BRCA','UCEC','HNSC','KIRC','LUAD','THCA','PRAD','LUSC','OV','STAD','BLCA','COAD','LIHC','CESC','KIRP')

#load the mRNA data:


all_mRNA_datasets <- list();
for (cancer_type in cancer_types){
	#print(cancer_type)
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

#load the miRNA data:

all_miRNA_datasets <- list();
all_miRNA <- c()
for (cancer_type in cancer_types){
	fname_miRNA <- paste0('../Reprocessed GDAC data/',cancer_type,'/miRNA/tumour/cleaned_miRNA_mature.txt')
	all_miRNA_datasets[[cancer_type]] <- read.table(fname_miRNA, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(all_miRNA_datasets[[cancer_type]]) <- gsub('[.]','-',colnames(all_miRNA_datasets[[cancer_type]]))
	all_miRNA <- unique(c(all_miRNA,rownames(all_miRNA_datasets[[cancer_type]])))
}

# load the mutation data:
meth_gene_names <- c('PTEN','ACVR2A','ARHGEF12','CDK12','DNMT3A','FAT4','SFRP4','TGFBR2')
entrez_gene_names <- list()
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
for(gene_name in meth_gene_names){
	entrez_gene_names[[gene_name]] <- xx[[gene_name]]
}

mut_data <- list()
mut_data_overall <- list()
for(gene_name in meth_gene_names){
	print(gene_name)
	mut_data[[gene_name]] <- list()
	mut_data_overall[[gene_name]] <- list()
	for(cancer_type in cancer_types){
		#load mutation data

		fName_mut <- paste0('../Reprocessed GDAC data/',cancer_type,'/mutation/mutations.txt')
		mut_data_ca_type <- read.table(fName_mut, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(mut_data_ca_type) <- gsub('[.]','-',colnames(mut_data_ca_type))
		
		final_mat <- matrix(0,nrow=18,ncol=length(colnames(mut_data_ca_type)))
		colnames(final_mat) <- colnames(mut_data_ca_type)
		row.names(final_mat) <- c('Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site',
			'Translation_Start_Site','Nonstop_Mutation','3\'UTR','5\'UTR','3\'Flank',
			'5\'Flank','RNA','Targeted_Region','In_Frame_Del','In_Frame_Ins',
			'IGR','Intron','Silent')
		for(i in 1:18){
			final_mat[i,which(mut_data_ca_type[gene_name,] == i)] <- 1 
		}
		if(sum(rowSums(final_mat)!=0)==1){
			rname <- rownames(final_mat)[which(rowSums(final_mat)!=0)]
			final_mat <- rbind(final_mat[which(rowSums(final_mat)!=0),])
			row.names(final_mat) <- rname
		}else if(sum(rowSums(final_mat)!=0)==0){
				final_mat <- NULL

		}else{
			final_mat <- final_mat[which(rowSums(final_mat)!=0),]
		}
		mut_data[[gene_name]][[cancer_type]] <- final_mat
		mut_data_overall[[gene_name]][[cancer_type]] <- ((mut_data_ca_type[gene_name,] > 0) & (mut_data_ca_type[gene_name,] < 16)) * 1
	}
}

for(gene_name in meth_gene_names){
	print(gene_name)
	#es(mut_data_ca_type) <- gsub('[.]','-',colnames(mut_data_ca_type))
	for(cancer_type in cancer_types){
		print(dim(mut_data[[gene_name]][[cancer_type]]))
	}
}
# meth_gene_names <- c('PTEN','ACVR2A','ARHGEF12','CYLD','DNMT3A','FAT4','SFRP4','TGFBR2')
#load in methylation data:
methylation_data <- list()
for(gene_name in meth_gene_names){
	methylation_data[[gene_name]] <- list()
	for(cancer_type in cancer_types){
		#load mutation data
		fName_meth <- paste0('../Reprocessed GDAC data/',cancer_type,'/methylation/tumour/',gene_name,'/cleaned_methylation.txt')
		methylation_data[[gene_name]][[cancer_type]] <- read.table(fName_meth, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(methylation_data[[gene_name]][[cancer_type]]) <- gsub('[.]','-',colnames(methylation_data[[gene_name]][[cancer_type]]))

	}
}

#load in CNV data:
cnv_data <- list()
for(cancer_type in cancer_types){
	#load mutation data
	fName_cnv <- paste0('../Reprocessed GDAC data/',cancer_type,'/cnv/tumour/cleaned_cnv.txt')
	cnv_data[[cancer_type]] <- read.table(fName_cnv, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(cnv_data[[cancer_type]]) <- gsub('[.]','-',colnames(cnv_data[[cancer_type]]))

}

for(gene_name in meth_gene_names){
	dir.create(paste0(gene_name,'-analysis'))
	common_samples <- list()
	for(cancer_type in cancer_types){

		print(cancer_type)

		fName_mut <- paste0('../Reprocessed GDAC data/',cancer_type,'/mutation/mutations.txt')
		mut_data_ca_type <- read.table(fName_mut, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
		colnames(mut_data_ca_type) <- gsub('[.]','-',colnames(mut_data_ca_type))

		print(length(intersect(colnames(all_mRNA_datasets[[cancer_type]]), 
			intersect(colnames(all_miRNA_datasets[[cancer_type]]),
				intersect(colnames(methylation_data[[gene_name]][[cancer_type]]),
					intersect(colnames(cnv_data[[cancer_type]]),
						colnames(mut_data_ca_type)))))))

		common_samples[[cancer_type]] <- intersect(colnames(all_mRNA_datasets[[cancer_type]]), 
			intersect(colnames(all_miRNA_datasets[[cancer_type]]),
				intersect(colnames(methylation_data[[gene_name]][[cancer_type]]),
					intersect(colnames(cnv_data[[cancer_type]]),
						colnames(mut_data_ca_type)))))
	}
	# dir.create('PTEN_analysis')
	all_coefficients <- list()
	for(cancer_type in cancer_types) {
		if(length(common_samples[[cancer_type]]) > 9){
			#create the linear model and do penalised regression on it
		pten_expr <- all_mRNA_datasets[[cancer_type]][entrez_gene_names[[gene_name]],common_samples[[cancer_type]]]
		miRNA_data <- all_miRNA_datasets[[cancer_type]][,common_samples[[cancer_type]]]
		meth_data <- methylation_data[[gene_name]][[cancer_type]][,common_samples[[cancer_type]]]
		cnv_pten <- cnv_data[[cancer_type]][entrez_gene_names[[gene_name]],common_samples[[cancer_type]]]

		if(!is.null(mut_data[[gene_name]][[cancer_type]])){
			mutation_data <- rbind(mut_data[[gene_name]][[cancer_type]][,common_samples[[cancer_type]]])
			if(is.null(rownames(mutation_data))){
				# row.names(mutation_data) <- rownames(mut_data[[gene_name]][[cancer_type]])
				mutation_data <- rbind(mutation_data[which(rowSums(mutation_data)!=0),])
				if(dim(mutation_data)[1]==0){
					mutation_data <- NULL
				}else{	
					row.names(mutation_data) <- rownames(mut_data[[gene_name]][[cancer_type]])
				}

			}else{
				mutation_data <- mutation_data[which(rowSums(mutation_data)!=0),]
			}
			if(is.null(dim(mutation_data))){
				mutation_data <- rbind(mutation_data)
			}else if(dim(mutation_data)[1]==0){
				mutation_data <- NULL
			}
		}else{
			mutation_data <- NULL
		}

		#expression filter

		#expression filter for miRNA
		expression_threshold <- 0.80 # means that at least 10% of samples must have a nonzero value of the mRNA
		miRNA_data <-miRNA_data[which((rowSums(miRNA_data==0)) < ((1-expression_threshold) * length(colnames(miRNA_data)))),]
		miRNA_data <- as.matrix(log2(miRNA_data))


		expression_threshold <- 0.50 # means that at least 10% of samples must have a nonzero value of the mRNA
		meth_data <-meth_data[which((rowSums(is.na(meth_data) )) < ((1-expression_threshold) * length(colnames(meth_data)))),]
		

		miRNA_data[!(is.finite(miRNA_data))] <- NA
		miRNA_data <-miRNA_data[which((rowSums(is.na(miRNA_data))) < (0.5 * length(colnames(miRNA_data)))),]
		
		#row transforms for the mutation data
		if(!is.null(mutation_data)){
			mutation_data <-t(apply(as.matrix(mutation_data),1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))
			# (as.numeric(mutation_data) - mean(as.numeric(mutation_data),na.rm=T))/sd(as.numeric(mutation_data),na.rm=T)
		}

		# #z-transform the scores
		pten_expr <- (as.numeric(pten_expr) - mean(as.numeric(pten_expr),na.rm=T))/sd(as.numeric(pten_expr),na.rm=T)
		cnv_pten <- (as.numeric(cnv_pten) - mean(as.numeric(cnv_pten),na.rm=T))/sd(as.numeric(cnv_pten),na.rm=T)
		miRNA_data <- t(apply(miRNA_data,1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))
		meth_data <- t(apply(meth_data,1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))

		#first we need to subset the data into folds
		new_df <- na.omit(t(rbind(pten_expr,cnv_pten,mutation_data,miRNA_data,meth_data)))

		folds <- 10
		nrows_combined_df <- 1:dim(new_df)[1]
		best_overall_error <- 99999999
		for (i in 0:(folds-1)){
			new_df_subset <- as.data.frame(new_df[!(nrows_combined_df%%folds==i),]) #takes out the 1/nth row of the data set
			#train the univaraite model
			#put these as inputs to the penalized model

			linear_models_miRNA_meth <- matrix(0,nrow=(length(rownames(miRNA_data)) + length(rownames(meth_data)) + length(rownames(mutation_data))),ncol=1)
			row.names(linear_models_miRNA_meth) <- c(rownames(mutation_data),rownames(miRNA_data),rownames(meth_data))
			for (j in 1:(length(rownames(miRNA_data)) + length(rownames(meth_data)) + length(rownames(mutation_data)))) {
				univariate_data <- as.data.frame(cbind(new_df_subset[,1],new_df_subset[,(j+2)]))
				colnames(univariate_data) <- c('sig_score','miRNA_or_meth')
				# print(univariate_data)

				univariate_model <- lm(formula = sig_score ~ miRNA_or_meth,data = univariate_data)
				# print(summary(univariate_model))
				if(dim(summary(univariate_model)$coefficients)[1] > 1){
		 			linear_models_miRNA_meth[j] <- (summary(univariate_model)$coefficients)[2,4]
		 		}else{
		 			linear_models_miRNA_meth[j] <- 1
		 		}

			}

			#significant miRNAs are those w p < 0.2:
			significant_miRNAs <- rownames(linear_models_miRNA_meth)[which(linear_models_miRNA_meth < 0.2 & !is.nan(linear_models_miRNA_meth))]
			
			#penalised linear regression

			lambda_2_values <- c(0, 0.01, 0.1,1,10,100)
			max_likelihood <- -9999999999
			for (lambda2_val in lambda_2_values){

				cross_val_model <- optL1(response = new_df_subset[,1],penalized = new_df_subset[,c('cnv_pten',significant_miRNAs)], lambda2 = lambda2_val,data=as.data.frame(new_df_subset),model="linear",fold=10,trace=T)#,trace=F,maxiter=1000,tol=.Machine$double.eps^0.23)
				if ((cross_val_model$fullfit)@loglik > max_likelihood){
					best_model <<- cross_val_model
					best_lambda <- lambda2_val
				}
			}

			#now that we know the best model, let's test it on the other 1/n of the data, and record the error
			unused_df <- as.data.frame(new_df[(nrows_combined_df%%folds==i),])
			current_predictions <- predict(best_model$fullfit, penalized=unused_df[,c('cnv_pten',significant_miRNAs)],data=unused_df)

			cur_error <- norm((as.numeric(unused_df[,1]) - as.numeric(current_predictions)),type="2")
			# print(cur_error)
			if (cur_error < best_overall_error){
				best_overall_error <- cur_error
				best_overall_model <- best_model
				best_overall_lambda <- best_lambda
			}
		}

		# miRNA_names_reported <- intersect(names(coef(best_overall_model$fullfit)), rownames(miRNA_data))
		#best_coef_matrix[rownames(miRNA_matrix)[i],mRNA_names_reported] <- coef(best_model$fullfit)[mRNA_names_reported]

		#return the coefficients
		# coef(best_overall_model$fullfit)[miRNA_names_reported]

		all_coefficients[[cancer_type]] <- coef(best_overall_model$fullfit)
		print(sort(coef(best_overall_model$fullfit)))

		}
	}



	all_combined_rownames <- c()
	for(cancer_type in cancer_types){
		all_combined_rownames <- c(all_combined_rownames, names(all_coefficients[[cancer_type]]))
		all_combined_rownames <- unique(all_combined_rownames)
	}

	overall_mat <- matrix(0,nrow=length(all_combined_rownames),ncol=length(names(all_coefficients)))
	row.names(overall_mat) <- all_combined_rownames
	colnames(overall_mat) <- names(all_coefficients)
	for(rName in all_combined_rownames){
		for(cancer_type in names(all_coefficients)){
			overall_mat[names(all_coefficients[[cancer_type]]),cancer_type] <- all_coefficients[[cancer_type]]
		}
	}
	library(RankProd)

	intercept_row <- which(rownames(overall_mat) == '(Intercept)')
	overall_mat <- overall_mat[-intercept_row,]

	ranked_cors <- RP(overall_mat,cl=rep(1,length(colnames(overall_mat))))
	ranked_cors_table <- topGene(ranked_cors,method='pfp',cutoff=0.05,gene.names = rownames(overall_mat))
	# ranked_cors_table <- ranked_cors_table$Table1

	save(overall_mat,all_coefficients,ranked_cors_table,file=paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))#overall_coeff_mat_M_values_z_transf.rda')
}



#first i want to see all of the overlal rank prod matrices

for(gene_name in meth_gene_names){
	print(gene_name)
	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))
	write.table(ranked_cors_table$Table1,sep='\t',file=paste0(gene_name,'-analysis/overall_rank_prod_neg_pred.txt'),quote=F)
	write.table(ranked_cors_table$Table2,sep='\t',file=paste0(gene_name,'-analysis/overall_rank_prod_pos_pred.txt'),quote=F)

}




#this is just for plotting of the heatmaps that I want
for(gene_name in meth_gene_names){
	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))
	p <- list()
	dir.create(paste0(gene_name,'-analysis/neg_reg_autocors'))
	dir.create(paste0(gene_name,'-analysis/neg_reg_autocors/density'))

	for(cancer_type in colnames(overall_mat)){
		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		neg_meth <- negative_predictors[which(grepl('cg', negative_predictors))]
		mut_types <- c('Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site',
			'Translation_Start_Site','Nonstop_Mutation','3\'UTR','5\'UTR','3\'Flank',
			'5\'Flank','RNA','Targeted_Region','In_Frame_Del','In_Frame_Ins',
			'IGR','Intron','Silent')
		mut_names <- negative_predictors[which(negative_predictors %in% mut_types)]
		if(!is.null(mut_names)){
			miR_meth_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))
			miR_meth_mat_all_values <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))

			row.names(miR_meth_mat) <- c(neg_miRs,neg_meth,mut_names)
			colnames(miR_meth_mat) <- c(neg_miRs,neg_meth,mut_names)

			row.names(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth,mut_names)
			colnames(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth,mut_names)
		}else{
			miR_meth_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))
			miR_meth_mat_all_values <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))

			row.names(miR_meth_mat) <- c(neg_miRs,neg_meth)
			colnames(miR_meth_mat) <- c(neg_miRs,neg_meth)

			row.names(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth)
			colnames(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth)
		}
		
		all_meth_mat_vals <- c()
	 
		for(miR_name in neg_miRs){
			for(meth_name in neg_meth){
				tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[miR_name,meth_name] <- tmp$estimate
					miR_meth_mat[meth_name,miR_name] <- tmp$estimate
				}
				miR_meth_mat_all_values[miR_name,meth_name] <- tmp$estimate
				miR_meth_mat_all_values[meth_name,miR_name] <- tmp$estimate
			}
		}

		for(miR_name in neg_miRs){
			for(miR_name2 in neg_miRs){
				tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(all_miRNA_datasets[[cancer_type]][miR_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[miR_name,miR_name2] <- tmp$estimate
				}
				miR_meth_mat_all_values[miR_name,miR_name2] <- tmp$estimate

			}
			all_meth_mat_vals <- rbind(all_meth_mat_vals,all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]])
		}


		for(meth_name in neg_meth){
			for(meth_name2 in neg_meth){
				tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[meth_name,meth_name2] <- tmp$estimate
				}
				miR_meth_mat_all_values[meth_name,meth_name2] <- tmp$estimate

			}
			all_meth_mat_vals <- rbind(all_meth_mat_vals,methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]])
		}
		if(!is.null(mut_names)){
			for(mut_name in mut_names){
				for(miR_name in neg_miRs){
					tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[miR_name,mut_name] <- tmp$estimate
						miR_meth_mat[mut_name,miR_name] <- tmp$estimate
					}
					miR_meth_mat_all_values[miR_name,mut_name] <- tmp$estimate
					miR_meth_mat_all_values[mut_name,miR_name] <- tmp$estimate
					
				}
				for(meth_name in neg_meth){
					tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[meth_name,mut_name] <- tmp$estimate
						miR_meth_mat[mut_name,meth_name] <- tmp$estimate
					}
					miR_meth_mat_all_values[meth_name,mut_name] <- tmp$estimate
					miR_meth_mat_all_values[mut_name,meth_name] <- tmp$estimate
				}

				for(mut_name2 in mut_names){
					tmp <- cor.test(as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[mut_name,mut_name2] <- tmp$estimate
					}
					miR_meth_mat_all_values[mut_name,mut_name2] <- tmp$estimate

				}

			}
		}
		# if('mutation_data' %in% negative_predictors){
		# 	#compute miRNA-mutatin
		# 	for(miR_name in neg_miRs){
		# 		tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
		# 		if(tmp$p.value < 0.05){
		# 			miR_meth_mat[miR_name,'mutation_data'] <- tmp$estimate
		# 			miR_meth_mat['mutation_data',miR_name] <- tmp$estimate
		# 		}
		# 	}
		# 	#compute meth-mutation
		# 	for(meth_name in neg_meth){
		# 		tmp <- cor.test(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
		# 		if(tmp$p.value < 0.05){
		# 			miR_meth_mat[meth_name,'mutation_data'] <- tmp$estimate
		# 			miR_meth_mat['mutation_data',meth_name] <- tmp$estimate
		# 		}
		# 	}
		# 	#compute muatiosn-mutation (this is 1)
		# 	miR_meth_mat['mutation_data','mutation_data'] <- 1
		# }

		row.names(all_meth_mat_vals) <-  gsub('[-]','.',row.names(all_meth_mat_vals) )
		if(dim(miR_meth_mat)[1]> 1 & dim(miR_meth_mat)[2]>1){
			dev.new()

			p[[cancer_type]] <- gplots::heatmap.2( miR_meth_mat ,
			                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
			                   trace = "none",
			                #   xlab = "Gene ID",
			                #   ylab="Gene ID",
			                   na.color="grey",
			                   labRow=rownames(miR_meth_mat),#rownames(tsg_sig_score_cor_by_gene[[tsg_name]]),#converted_symbols_map[rownames(sig_mut_props_hmap)],
			                   #labCol=colnames(autocors),#gene_sig,
			                   main = paste0('Autocorrelation-',gene_name, ' neg. regulators\n',cancer_type),#paste0("Partial to PTEN mut"),
			                   dendrogram = "both",
			                   breaks = seq(-1,1,length=101),
			                   #symbreaks = T,
			                   Rowv = T,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.8,cexCol=0.8)
			dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/',gene_name,'_neg_reg_autocor',cancer_type,'.pdf'),width=12,height=12)
		  	dev.off()
	  }
	}
	graphics.off()
}

#-------

#now let's test for mutual exlcusivity of the phenomena

all_miR_miR_cor_percentile <- list()
all_miR_meth_cor_percentile <- list()
all_miR_mut_cor_percentile <- list()
all_meth_meth_cor_percentile <- list()
all_meth_mut_cor_percentile <- list()
all_miR_miR_minus_meth_meth_cor_percentile <- list()
all_miR_miR_minus_miR_meth_cor_percentile <- list()
all_meth_meth_minus_miR_meth_cor_percentile <- list()

#first we can check the correlations between the individual variables that are sig in each cancer type
for(gene_name in meth_gene_names){

	all_miR_miR_cor_percentile[[gene_name]] <- list()
	all_miR_meth_cor_percentile[[gene_name]] <- list()
	all_miR_mut_cor_percentile[[gene_name]] <- list()
	all_meth_meth_cor_percentile[[gene_name]] <- list()
	all_meth_mut_cor_percentile[[gene_name]] <- list()
	all_miR_miR_minus_meth_meth_cor_percentile[[gene_name]] <- list()
	all_miR_miR_minus_miR_meth_cor_percentile[[gene_name]] <- list()
	all_meth_meth_minus_miR_meth_cor_percentile[[gene_name]] <- list()


	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))
	p <- list()
	dir.create(paste0(gene_name,'-analysis/neg_reg_autocors'))
	dir.create(paste0(gene_name,'-analysis/neg_reg_autocors/density'))
	for(cancer_type in colnames(overall_mat)){
		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		neg_meth <- negative_predictors[which(grepl('cg', negative_predictors))]
		mut_types <- c('Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site',
			'Translation_Start_Site','Nonstop_Mutation','3\'UTR','5\'UTR','3\'Flank',
			'5\'Flank','RNA','Targeted_Region','In_Frame_Del','In_Frame_Ins',
			'IGR','Intron','Silent')
		mut_names <- negative_predictors[which(negative_predictors %in% mut_types)]
		if(!is.null(mut_names)){
			miR_meth_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))
			miR_meth_mat_all_values <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))

			row.names(miR_meth_mat) <- c(neg_miRs,neg_meth,mut_names)
			colnames(miR_meth_mat) <- c(neg_miRs,neg_meth,mut_names)

			row.names(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth,mut_names)
			colnames(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth,mut_names)
		}else{
			miR_meth_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))
			miR_meth_mat_all_values <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))

			row.names(miR_meth_mat) <- c(neg_miRs,neg_meth)
			colnames(miR_meth_mat) <- c(neg_miRs,neg_meth)

			row.names(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth)
			colnames(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth)
		}
		
		all_meth_mat_vals <- c()
	 
		for(miR_name in neg_miRs){
			for(meth_name in neg_meth){
				tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[miR_name,meth_name] <- tmp$estimate
					miR_meth_mat[meth_name,miR_name] <- tmp$estimate
				}
				miR_meth_mat_all_values[miR_name,meth_name] <- tmp$estimate
				miR_meth_mat_all_values[meth_name,miR_name] <- tmp$estimate
			}
		}

		for(miR_name in neg_miRs){
			for(miR_name2 in neg_miRs){
				tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(all_miRNA_datasets[[cancer_type]][miR_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[miR_name,miR_name2] <- tmp$estimate
				}
				miR_meth_mat_all_values[miR_name,miR_name2] <- tmp$estimate

			}
			all_meth_mat_vals <- rbind(all_meth_mat_vals,all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]])
		}


		for(meth_name in neg_meth){
			for(meth_name2 in neg_meth){
				tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[meth_name,meth_name2] <- tmp$estimate
				}
				miR_meth_mat_all_values[meth_name,meth_name2] <- tmp$estimate

			}
			all_meth_mat_vals <- rbind(all_meth_mat_vals,methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]])
		}
		if(!is.null(mut_names)){
			for(mut_name in mut_names){
				for(miR_name in neg_miRs){
					tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[miR_name,mut_name] <- tmp$estimate
						miR_meth_mat[mut_name,miR_name] <- tmp$estimate
					}
					miR_meth_mat_all_values[miR_name,mut_name] <- tmp$estimate
					miR_meth_mat_all_values[mut_name,miR_name] <- tmp$estimate
					
				}
				for(meth_name in neg_meth){
					tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[meth_name,mut_name] <- tmp$estimate
						miR_meth_mat[mut_name,meth_name] <- tmp$estimate
					}
					miR_meth_mat_all_values[meth_name,mut_name] <- tmp$estimate
					miR_meth_mat_all_values[mut_name,meth_name] <- tmp$estimate
				}

				for(mut_name2 in mut_names){
					tmp <- cor.test(as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[mut_name,mut_name2] <- tmp$estimate
					}
					miR_meth_mat_all_values[mut_name,mut_name2] <- tmp$estimate

				}

			}
		}
		# if('mutation_data' %in% negative_predictors){
		# 	#compute miRNA-mutatin
		# 	for(miR_name in neg_miRs){
		# 		tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
		# 		if(tmp$p.value < 0.05){
		# 			miR_meth_mat[miR_name,'mutation_data'] <- tmp$estimate
		# 			miR_meth_mat['mutation_data',miR_name] <- tmp$estimate
		# 		}
		# 	}
		# 	#compute meth-mutation
		# 	for(meth_name in neg_meth){
		# 		tmp <- cor.test(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
		# 		if(tmp$p.value < 0.05){
		# 			miR_meth_mat[meth_name,'mutation_data'] <- tmp$estimate
		# 			miR_meth_mat['mutation_data',meth_name] <- tmp$estimate
		# 		}
		# 	}
		# 	#compute muatiosn-mutation (this is 1)
		# 	miR_meth_mat['mutation_data','mutation_data'] <- 1
		# }

		row.names(all_meth_mat_vals) <-  gsub('[-]','.',row.names(all_meth_mat_vals) )
		# ggpairs(as.data.frame(t(all_meth_mat_vals)))

		# ###########################PLOTTING ###########################
		dev.new()
		p[[cancer_type]] <- gplots::heatmap.2( miR_meth_mat ,
		                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
		                   trace = "none",
		                #   xlab = "Gene ID",
		                #   ylab="Gene ID",
		                   na.color="grey",
		                   labRow=rownames(miR_meth_mat),#rownames(tsg_sig_score_cor_by_gene[[tsg_name]]),#converted_symbols_map[rownames(sig_mut_props_hmap)],
		                   #labCol=colnames(autocors),#gene_sig,
		                   main = paste0('Autocorrelation-',gene_name, ' neg. regulators\n',cancer_type),#paste0("Partial to PTEN mut"),
		                   dendrogram = "both",
		                   breaks = seq(-1,1,length=101),
		                   #symbreaks = T,
		                   Rowv = T,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.8,cexCol=0.8)
		dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/',gene_name,'_neg_reg_autocor',cancer_type,'.pdf'),width=12,height=12)
	  	dev.off()
	  
# 		print(negative_predictors)
# 		if(length(neg_meth) > 0){
# 			dev.new()
# 			par(mfrow=c(ceiling(sqrt(length(neg_meth))), ceiling(length(neg_meth)/ceiling(sqrt(length(neg_meth))))),mar=c(2,2,2,2))
# 			for(meth_name in neg_meth){
# 				plot(density(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T),main=paste0(meth_name,'\n',cancer_type))
# 			}
# 			dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/density/',gene_name,'_meth_dens',cancer_type,'.pdf'),width=12,height=12)
# 		  	dev.off()
# 	  	}
# 	  	if(length(neg_miRs)> 0){
# 			dev.new()
# 			par(mfrow=c(ceiling(sqrt(length(neg_miRs))), ceiling(length(neg_miRs)/ceiling(sqrt(length(neg_miRs))))),mar=c(2,2,2,2))
# 			for(miR_name in neg_miRs){
# 				plot(density(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),na.rm=T),main=paste0(miR_name,'\n',cancer_type))
# 			}
# 			dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/density/',gene_name,'_mir_dens',cancer_type,'.pdf'),width=12,height=12)
# 		  	dev.off()
# 		}
# 	# }

# 	# for(cancer_type in colnames(overall_mat)){
# 		#also make the density plots for PTEN expression, CNV and mutation status
# 		dev.new()
# 		par(mfrow=c(1,2+length(mut_names)),mar=c(2,2,2,2))
# 		plot(density(as.numeric(all_mRNA_datasets[[cancer_type]][entrez_gene_names[[gene_name]],common_samples[[cancer_type]]]),na.rm=T),main=paste0(gene_name, ' mRNA Expression\n',cancer_type))
# 		plot(density(as.numeric(cnv_data[[cancer_type]][entrez_gene_names[[gene_name]],common_samples[[cancer_type]]]),na.rm=T),main=paste0(gene_name, ' CNV\n',cancer_type))
# 		for(mut_name in mut_names){
# 			plot(density(as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),na.rm=T),main=paste0(gene_name,' mutation status\n',mut_name,' ',cancer_type))
# 		}
# 		dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/density/',gene_name,'_expr_cnv_mut_dens',cancer_type,'.pdf'),width=12,height=4)
# 		dev.off()
# 		graphics.off()
###########################################################################################################################################################################


		#####THIS CODE IS JUST TO COUNT UP SIGNIFICANT CORS IN THE MATRIX TO SEE IF THINGS ARE OCCURRING TOGETHER####


		tot_prob_co_miRs <- 0
		tot_prob_co_miRs_all_values <- 0 # the all values thing repeats the calculation for the Full matrix not just the sig valuess

		for(i in 2:length(neg_miRs)){
			mir_1 <- neg_miRs[i]
			for(j in 1:(i-1)){
				mir_2 <- neg_miRs[j]
				tot_prob_co_miRs <- tot_prob_co_miRs + (miR_meth_mat[mir_1,mir_2] > 0)
				tot_prob_co_miRs_all_values <- tot_prob_co_miRs_all_values + (miR_meth_mat_all_values[mir_1,mir_2] > 0)

			}
		}
		tot_nums <- length(neg_miRs) * (length(neg_miRs) - 1) / 2
		tot_prob_co_miRs <- tot_prob_co_miRs / tot_nums
		tot_prob_co_miRs_all_values <- tot_prob_co_miRs_all_values / tot_nums

		tot_prob_co_miR_meth <- 0
		tot_prob_co_miR_meth_all_values <- 0

		for(miR_name in neg_miRs){
			for(meth_name in neg_meth){
				tot_prob_co_miR_meth <- tot_prob_co_miR_meth + (miR_meth_mat[miR_name,meth_name] > 0)
				tot_prob_co_miR_meth_all_values <- tot_prob_co_miR_meth_all_values + (miR_meth_mat_all_values[miR_name,meth_name] > 0)

			}
		}
		tot_nums <- length(neg_miRs) * length(neg_meth)
		tot_prob_co_miR_meth <- tot_prob_co_miR_meth / tot_nums
		tot_prob_co_miR_meth_all_values <- tot_prob_co_miR_meth_all_values / tot_nums

		tot_prob_co_meth <- 0
		tot_prob_co_meth_all_values <- 0

		if(length(neg_meth) > 1){
			for(i in 2:length(neg_meth)){
				meth_name_1 <- neg_meth[i]
				for(j in 1:(i-1)){
					meth_name_2 <- neg_meth[j]
					tot_prob_co_meth <- tot_prob_co_meth + (miR_meth_mat[meth_name_1,meth_name_2] > 0)
					tot_prob_co_meth_all_values <- tot_prob_co_meth_all_values + (miR_meth_mat_all_values[meth_name_1,meth_name_2] > 0)

				}
			}
			tot_nums <- length(neg_meth) * (length(neg_meth) - 1) / 2
			tot_prob_co_meth <- tot_prob_co_meth / tot_nums
			tot_prob_co_meth_all_values <- tot_prob_co_meth_all_values / tot_nums
		}

		tot_prob_meth_mut <- list()
		tot_prob_miR_mut <- list()
		tot_prob_meth_mut_all_values <- list()
		tot_prob_miR_mut_all_values <- list()
		if(!is.null(mut_names)){

			for(mut_name in mut_names){
				tot_prob_meth_mut[[mut_name]] <- 0
				tot_prob_meth_mut[[mut_name]] <- sum((miR_meth_mat[mut_name,neg_meth] > 0)) / length(neg_meth)

				tot_prob_meth_mut_all_values[[mut_name]] <- 0
				tot_prob_meth_mut_all_values[[mut_name]] <- sum((miR_meth_mat_all_values[mut_name,neg_meth] > 0)) / length(neg_meth)

				tot_prob_miR_mut[[mut_name]] <- 0
				tot_prob_miR_mut[[mut_name]] <- sum((miR_meth_mat[mut_name,neg_miRs] > 0)) / length(neg_miRs)
 
				tot_prob_miR_mut_all_values[[mut_name]] <- 0
				tot_prob_miR_mut_all_values[[mut_name]] <- sum((miR_meth_mat_all_values[mut_name,neg_miRs] > 0)) / length(neg_miRs)


				# for(miR_name in neg_miRs){
				# 	tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				# 	if(tmp$p.value < 0.05){
				# 		miR_meth_mat[miR_name,mut_name] <- tmp$estimate
				# 		miR_meth_mat[mut_name,miR_name] <- tmp$estimate
				# 	}
				# }
				# for(meth_name in neg_meth){
				# 	tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				# 	if(tmp$p.value < 0.05){
				# 		miR_meth_mat[meth_name,mut_name] <- tmp$estimate
				# 		miR_meth_mat[mut_name,meth_name] <- tmp$estimate
				# 	}
				# }

				# for(mut_name2 in mut_names){
				# 	tmp <- cor.test(as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				# 	if(tmp$p.value < 0.05){
				# 		miR_meth_mat[mut_name,mut_name2] <- tmp$estimate
				# 	}
				# }
			}

		}
		null_distr_values <- get_null_distr_p_value_mir_meth_mat(length(neg_meth),length(neg_miRs),mut_names,all_expr_mat_cors[[cancer_type]][[gene_name]])
		null_distr_values_all_values <- get_null_distr_p_value_mir_meth_mat(length(neg_meth),length(neg_miRs),mut_names,all_expr_mat_cors_all_values[[cancer_type]][[gene_name]])


		all_miR_miR_cor_percentile[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[1]])(tot_prob_co_miRs)
		if(length(neg_meth) > 0){
			all_miR_meth_cor_percentile[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[2]])(tot_prob_co_miR_meth)
			all_meth_meth_cor_percentile[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[3]])(tot_prob_co_meth)
		}
		all_miR_mut_cor_percentile[[gene_name]][[cancer_type]] <- list()
		all_meth_mut_cor_percentile[[gene_name]][[cancer_type]] <- list()

		if(!is.null(mut_names)){
			for(mut_name in mut_names){
					all_miR_mut_cor_percentile[[gene_name]][[cancer_type]][[mut_name]] <- ecdf(null_distr_values[[4]][[mut_name]])(tot_prob_miR_mut[[mut_name]])
				
				if(length(neg_meth) > 0){
					all_meth_mut_cor_percentile[[gene_name]][[cancer_type]][[mut_name]] <- ecdf(null_distr_values[[5]][[mut_name]])(tot_prob_meth_mut[[mut_name]])
				}
			}
		}
		if(length(neg_meth) > 0){

			all_miR_miR_minus_meth_meth_cor_percentile[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[1]] - null_distr_values[[3]])(tot_prob_co_miRs - tot_prob_co_meth)
			all_miR_miR_minus_miR_meth_cor_percentile[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[1]] - null_distr_values[[2]])(tot_prob_co_miRs - tot_prob_co_miR_meth)
			all_meth_meth_minus_miR_meth_cor_percentile[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[3]] - null_distr_values[[2]])(tot_prob_co_meth - tot_prob_co_miR_meth)
			print(all_miR_miR_minus_miR_meth_cor_percentile[[gene_name]][[cancer_type]] )#ecdf(null_distr_values[[1]])(tot_prob_co_miRs))

		}
		# print(tot_prob_co_miR_meth)
		# print(tot_prob_co_meth)
		# print(tot_prob_meth_mut)
		# print(tot_prob_miR_mut)

		#######		#######		#######		#######		#######		#######		#######		#######		#######		#######




	}
}




#PLOTTING COMMANDS FOR THE EXLCUSIVITY ANALYSIS
dev.new()
melted_plotting <- melt(all_miR_miR_minus_miR_meth_cor_percentile)
ggplot(as.data.frame(melted_plotting),aes(x=reorder(L1, value,  FUN=mean),y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme_minimal()+  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
labs(color=' ' , x=' ' , y='Empiric CDF percentile',title=expression(paste('eCDF: ',Pi[rho['miR-miR' > 0]] - Pi[rho['miR-methyl' > 0]]))) #miR-miR proportion pos. cor. \nminus miR-meth proportion pos. cor.'))
dev.copy(pdf,'all_miR_miR_minus_miR_meth_cor_percentile.pdf',width=4,height=4)
dev.off()

dev.new()
melted_plotting <- melt(all_miR_miR_minus_meth_meth_cor_percentile)
ggplot(as.data.frame(melted_plotting),aes(x=reorder(L1, value,  FUN=mean),y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme_minimal()+ theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color=' ',x=' ',y='Empiric CDF percentile',title=expression(paste('eCDF: ',Pi[rho['miR-miR' > 0]] - Pi[rho['methyl-methyl' > 0]])))
dev.copy(pdf,'all_miR_miR_minus_meth_meth_cor_percentile.pdf',width=4,height=4)
dev.off()

dev.new()
melted_plotting <- melt(all_meth_meth_minus_miR_meth_cor_percentile)
ggplot(as.data.frame(melted_plotting),aes(x=reorder(L1, value,  FUN=mean),y=value)) +  geom_boxplot() + geom_point(aes(color=L2))+ theme_minimal() + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color=' ',x=' ',y='Empiric CDF percentile',title=expression(paste('eCDF: ', Pi[rho['methyl-methyl' > 0]]- Pi[rho['miR-methyl' > 0]])))
dev.copy(pdf,'all_meth_meth_minus_miR_meth_cor_percentile.pdf',width=4,height=4)
dev.off()

melted_plotting <- melt(all_miR_miR_cor_percentile)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, miR-miR proportion pos. cor. ')


melted_plotting <- melt(all_miR_meth_cor_percentile)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, miR-meth proportion pos. cor. ')

melted_plotting <- melt(all_meth_meth_cor_percentile)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, meth-meth proportion pos. cor. ')

#next task is to make plots of the copy number and types of mutation and expression of the various miRs/methylation probes

#probably easier to show z-score to make rows relatively comparable?
library(ComplexHeatmap)
for(gene_name in meth_gene_names){
	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))

	for(cancer_type in colnames(overall_mat)){
		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		neg_meth <- negative_predictors[which(grepl('cg', negative_predictors))]
		#get the data we need
		miRNA_data <- all_miRNA_datasets[[cancer_type]][neg_miRs,common_samples[[cancer_type]]]
		meth_data <- methylation_data[[gene_name]][[cancer_type]][neg_meth,common_samples[[cancer_type]]]

		miRNA_data <- t(apply(miRNA_data,1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))
		meth_data <- t(apply(meth_data,1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))

		cnv_pten <- cnv_data[[cancer_type]][entrez_gene_names[[gene_name]],common_samples[[cancer_type]]]

		#first we need to subset the samples for the distinct regions of the hmap
		#first we will arrange the heatmap my copy number from lowest to highest and add annotations for the mut types?
		
		if(dim(meth_data)[2] > 0 && dim(miRNA_data)[2] > 0){
			plotting_mat <- rbind(miRNA_data,meth_data)
		}else if(dim(meth_data)[2] == 0) {
				plotting_mat <- miRNA_data
		}else{
			plotting_mat <- meth_data

		}
		colnames(plotting_mat) <-common_samples[[cancer_type]]
		plotting_mat <- plotting_mat[,order(as.numeric(cnv_pten))]

		#load in the mutation data

		if(!is.null(mut_data[[gene_name]][[cancer_type]])){
			mutation_data <- rbind(mut_data[[gene_name]][[cancer_type]][,common_samples[[cancer_type]]])
			if(is.null(rownames(mutation_data))){
				# row.names(mutation_data) <- rownames(mut_data[[gene_name]][[cancer_type]])
				mutation_data <- rbind(mutation_data[which(rowSums(mutation_data)!=0),])
				if(dim(mutation_data)[1]==0){
					mutation_data <- NULL
				}else{	
					row.names(mutation_data) <- rownames(mut_data[[gene_name]][[cancer_type]])
				}

			}else{
				mutation_data <- mutation_data[which(rowSums(mutation_data)!=0),]
			}
			if(is.null(dim(mutation_data))){
				mutation_data <- rbind(mutation_data)
			}else if(dim(mutation_data)[1]==0){
				mutation_data <- NULL
			}
		}else{
			mutation_data <- NULL
		}

		mutation_annotation <- rep('Non-mut',times=length(as.numeric(cnv_pten)))
		
		if(!is.null(mutation_data)){
			for(mut_rName in rownames(mutation_data)){
				mutation_annotation[which(mutation_data[mut_rName,]==1)] <-mut_rName
			}
		}
		mutation_annotation <- mutation_annotation[order(cnv_pten)]
		ha <- HeatmapAnnotation(data.frame(cnv=as.numeric(cnv_pten[order(as.numeric(cnv_pten))]),mutation=mutation_annotation))
		draw(Heatmap(plotting_mat,top_annotation= ha,,cluster_columns = F,show_row_names=T,show_column_names = F,column_title=paste0('Expression of negative predictors for \n ',gene_name,' in ',cancer_type),name='z-transformed expr') )
		dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/density/',gene_name,'_expr_meth_mut_',cancer_type,'.pdf'),width=12,height=12)
		  	dev.off()
	}
}


############ THIS IS THE CODE TO DETERMINE WHETHER THE MIRNA / METHYLATION IS DIFFERENT BETWEEN MUTATED AND UNMUTATED SAMPLES OF THE TSG
wilcox_test_results <- list()
wilcox_test_results_overall_mut <- list()
for(gene_name in meth_gene_names){
	wilcox_test_results[[gene_name]] <- list()
	wilcox_test_results_overall_mut[[gene_name]] <- list()

	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))
	wilcox_test_results_overall_mut[[gene_name]] <- matrix(NA,nrow=dim(overall_mat)[1],ncol=dim(overall_mat)[2])#[[cancer_type]] <- list()
	row.names(wilcox_test_results_overall_mut[[gene_name]] ) <- rownames(overall_mat)
	colnames(wilcox_test_results_overall_mut[[gene_name]]) <- colnames(overall_mat)
	for(cancer_type in colnames(overall_mat)){
		wilcox_test_results[[gene_name]][[cancer_type]] <- list()

		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		neg_meth <- negative_predictors[which(grepl('cg', negative_predictors))]
		#get the data we need
		miRNA_data <- all_miRNA_datasets[[cancer_type]][neg_miRs,common_samples[[cancer_type]]]
		meth_data <- methylation_data[[gene_name]][[cancer_type]][neg_meth,common_samples[[cancer_type]]]

		# miRNA_data <- t(apply(miRNA_data,1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))
		# meth_data <- t(apply(meth_data,1,function(x){(as.numeric(x)-mean(as.numeric(x),na.rm=T))/sd(as.numeric(x),na.rm=T)}))

		cnv_pten <- cnv_data[[cancer_type]][entrez_gene_names[[gene_name]],common_samples[[cancer_type]]]

		#load in the mutation data

		if(!is.null(mut_data[[gene_name]][[cancer_type]])){
			mutation_data <- rbind(mut_data[[gene_name]][[cancer_type]][,common_samples[[cancer_type]]])
			if(is.null(rownames(mutation_data))){
				# row.names(mutation_data) <- rownames(mut_data[[gene_name]][[cancer_type]])
				mutation_data <- rbind(mutation_data[which(rowSums(mutation_data)!=0),])
				if(dim(mutation_data)[1]==0){
					mutation_data <- NULL
				}else{	
					row.names(mutation_data) <- rownames(mut_data[[gene_name]][[cancer_type]])
				}
			}else{
				mutation_data <- mutation_data[which(rowSums(mutation_data)!=0),]
			}
			if(is.null(dim(mutation_data))){
				mutation_data <- rbind(mutation_data)
			}else if(dim(mutation_data)[1]==0){
				mutation_data <- NULL
			}
		}else{
			mutation_data <- NULL
		}


		#now we have appropriately processed the mutation data matrix:
		if(!is.null(mutation_data)){
			#here, we have mutation data and some samples are actually mutated
			#lets pick out the mutated samples

			for(mutation_type in rownames(mutation_data)){
				wilcox_test_results[[gene_name]][[cancer_type]][[mutation_type]] <- list()
				mut_samples <- which(mutation_data[mutation_type,] == 1)
				non_mut_samples <- which(mutation_data[mutation_type,] == 0)
				#loop over all of the miRNA and methylation probes and ask if there is differential expression:

				for(neg_miR_name in neg_miRs){
					if((sum(is.na(as.numeric(miRNA_data[neg_miR_name,mut_samples]))) < length(mut_samples)) && (sum(is.na(as.numeric(miRNA_data[neg_miR_name,non_mut_samples]))) < length(non_mut_samples))){
						wcox_test <- wilcox.test(as.numeric(miRNA_data[neg_miR_name,mut_samples]),as.numeric(miRNA_data[neg_miR_name,non_mut_samples]),alternative='less',paired=F) #tests if in UNMUTATED samples it is higher than in MUTATED samples (ie. when there is mutation, we need less of it because already inactivated)
						wilcox_test_results[[gene_name]][[cancer_type]][[mutation_type]][[neg_miR_name]] <- wcox_test$p.value
					}
				}
				for(neg_meth_name in neg_meth){
					if((sum(is.na(as.numeric(meth_data[neg_meth_name,mut_samples]))) < length(mut_samples)) && (sum(is.na(as.numeric(meth_data[neg_meth_name,non_mut_samples]))) < length(non_mut_samples))){
						wcox_test <- wilcox.test(as.numeric(meth_data[neg_meth_name,mut_samples]),as.numeric(meth_data[neg_meth_name,non_mut_samples]),alternative='less',paired=F) 
						wilcox_test_results[[gene_name]][[cancer_type]][[mutation_type]][[neg_meth_name]] <- wcox_test$p.value
					}
				}
			}
		}

		#NOW WE CAN REPEAT THE ANALYSIS WITH THE overall mutation data as well

		mutation_data <- mut_data_overall[[gene_name]][[cancer_type]][,common_samples[[cancer_type]]]
		mut_samples <- which(mutation_data  > 0)
		non_mut_samples <- which(mutation_data == 0)
		for(neg_miR_name in neg_miRs){
			if((sum(is.na(as.numeric(miRNA_data[neg_miR_name,mut_samples]))) < length(mut_samples)) && (sum(is.na(as.numeric(miRNA_data[neg_miR_name,non_mut_samples]))) < length(non_mut_samples))){
				wcox_test <- wilcox.test(na.omit(as.numeric(miRNA_data[neg_miR_name,mut_samples])),na.omit(as.numeric(miRNA_data[neg_miR_name,non_mut_samples])),alternative='less',paired=F) #tests if in UNMUTATED samples it is higher than in MUTATED samples (ie. when there is mutation, we need less of it because already inactivated)
				wilcox_test_results_overall_mut[[gene_name]][neg_miR_name,cancer_type] <- wcox_test$p.value
				#wilcox_test_results_overall_mut[[gene_name]][[cancer_type]][[neg_miR_name]] <- wcox_test$p.value
				if(wcox_test$p.value < 0.05){
					print(gene_name)
					print(cancer_type)
					print(neg_miR_name)
					print(mean(na.omit(as.numeric(miRNA_data[neg_miR_name,mut_samples]))))
					print(mean(na.omit(as.numeric(miRNA_data[neg_miR_name,non_mut_samples]))))
					print(wilcox_test_results_overall_mut[[gene_name]][neg_miR_name,cancer_type])
				}
			}
		}
		for(neg_meth_name in neg_meth){
			if((sum(is.na(as.numeric(meth_data[neg_meth_name,mut_samples]))) < length(mut_samples)) && (sum(is.na(as.numeric(meth_data[neg_meth_name,non_mut_samples]))) < length(non_mut_samples))){
				wcox_test <- wilcox.test(na.omit(as.numeric(meth_data[neg_meth_name,mut_samples])),na.omit(as.numeric(meth_data[neg_meth_name,non_mut_samples])),alternative='less',paired=F) 
				wilcox_test_results_overall_mut[[gene_name]][neg_meth_name,cancer_type] <- wcox_test$p.value

				# wilcox_test_results_overall_mut[[gene_name]][[cancer_type]][[neg_meth_name]] <- wcox_test$p.value
				if(wcox_test$p.value < 0.05){
					print(gene_name)
					print(cancer_type)
					print(neg_meth_name)
					print(mean(na.omit(as.numeric(meth_data[neg_meth_name,mut_samples]))))
					print(mean(na.omit(as.numeric(meth_data[neg_meth_name,non_mut_samples]))))
					print(wilcox_test_results_overall_mut[[gene_name]][neg_meth_name,cancer_type])
				}
			}
		}
	}
}

##NOW let's reformat this data into something we WANT

for(gene_name in meth_gene_names){

	wilcox_test_results_overall_mut_refined <- wilcox_test_results_overall_mut[[gene_name]][which(rowSums(!is.na(wilcox_test_results_overall_mut[[gene_name]]))>0),]
		print(gene_name)
	wilcox_test_results_overall_mut_refined <- wilcox_test_results_overall_mut_refined[order(rowMeans(wilcox_test_results_overall_mut_refined,na.rm=T)),]
	print(wilcox_test_results_overall_mut_refined)
	write.table(wilcox_test_results_overall_mut_refined,quote=F,file=paste0(gene_name,'-analysis/wilcoxon_test_overall_mut.txt'),sep='\t')
	ranked_cors <- RP(wilcox_test_results_overall_mut_refined,cl=rep(1,length(colnames(wilcox_test_results_overall_mut_refined))))
	ranked_cors_table <- topGene(ranked_cors,method='pfp',cutoff=0.05,gene.names = rownames(wilcox_test_results_overall_mut_refined))
	write.table(ranked_cors_table$Table1,quote=F,file=paste0(gene_name,'-analysis/wilcoxon_test_overall_mut_rankProd_down.txt'),sep='\t')
	write.table(ranked_cors_table$Table2,quote=F,file=paste0(gene_name,'-analysis/wilcoxon_test_overall_mut_rankProd_up.txt'),sep='\t')

	# print(ranked_cors_table)
	# for(cancer_type in colnames(overall_mat)){
	# 	tmp <- melt(wilcox_test_results_overall_mut[[gene_name]][[cancer_type]])
	# 	if(!is.null(dim(tmp))){
	# 		print(tmp[order(tmp[,1]),])
	# 	}
		
	# }
}




###############################################################################################
#############LET'S ALSO FIGURE OUT HOW MANY OF THE NEGATIVELY ASSOICATED MIRS HAVE EACH TSG AS A TARGET

load('all_miRNA_targets.rda')
good_mirs <- list()
for(gene_name in meth_gene_names){
	good_mirs[[gene_name]] <- c()
	for(miR_name in names(targets_list)){
		print(miR_name)
		if(entrez_gene_names[[gene_name]] %in% targets_list[[miR_name]]){
			good_mirs[[gene_name]] <- c(good_mirs[[gene_name]],miR_name)
		}
	}
}

#now do the fisher test of overlap:
fisher_p_values <- c()
proportion_of_overlap <- c()
summary_mat_rnames <- c()
length_neg_mirs <- c()
length_targeting_mirs<- c()
for(gene_name in meth_gene_names){
	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))
	for(cancer_type in colnames(overall_mat)){

		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		summary_mat_rnames <- c(summary_mat_rnames,paste0(gene_name,' ',cancer_type))
		# print(length(intersect(neg_miRs,good_mirs[[gene_name]])) / length(neg_miRs))
		length_neg_mirs <- c(length_neg_mirs,length(neg_miRs))
		length_targeting_mirs <- c(length_targeting_mirs,length(good_mirs[[gene_name]]))
		cur_row_labels <- row_labels[grepl(tsg_name,row_labels)]
		contingency_table <- matrix(0,nrow=2,ncol=2)
		contingency_table[1,1] <- length(intersect(neg_miRs,good_mirs[[gene_name]]))#length(cur_row_labels)
		contingency_table[1,2] <- length(neg_miRs) - contingency_table[1,1]
		contingency_table[2,1] <- length(good_mirs[[gene_name]]) - contingency_table[1,1]
		contingency_table[2,2] <-  length(names(targets_list)) - contingency_table[1,1] - contingency_table[1,2] - contingency_table[2,1]
		fisher_test <- fisher.test(contingency_table)
		print(fisher_test)
		fisher_p_values <- c(fisher_p_values,fisher_test$p.value)
		proportion_of_overlap <- c(proportion_of_overlap,length(intersect(neg_miRs,good_mirs[[gene_name]])) / length(neg_miRs))
		# if(fisher_test$p.value < 0.05){
		# 	print(tsg_name)
		# 	print(contingency_table)
		# 	print(fisher_test)

		# }



		
	#check what the overlap is
	}
}
summary_mat <- cbind(fisher_p_values,proportion_of_overlap,length_neg_mirs,length_targeting_mirs)
row.names(summary_mat) <- summary_mat_rnames
write.table(summary_mat,quote=F,file=paste0('tsg_neg_miRs_targeting_overlap.txt'),sep='\t')

##########################################################################


#----now let's try and use the median cutoff/fisher exact test on the samples to figure out mutual exclusivity
for(gene_name in meth_gene_names){
	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))

	for(cancer_type in colnames(overall_mat)){

		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		neg_meth <- negative_predictors[which(grepl('cg', negative_predictors))]
		mut_types <- c('Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site',
				'Translation_Start_Site','Nonstop_Mutation','3\'UTR','5\'UTR','3\'Flank',
				'5\'Flank','RNA','Targeted_Region','In_Frame_Del','In_Frame_Ins',
				'IGR','Intron','Silent')
			mut_names <- negative_predictors[which(negative_predictors %in% mut_types)]
			if(!is.null(mut_names)){
				miR_meth_fisher_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))
				row.names(miR_meth_fisher_mat) <- c(neg_miRs,neg_meth,mut_names)
				colnames(miR_meth_fisher_mat) <- c(neg_miRs,neg_meth,mut_names)
			}else{
				miR_meth_fisher_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))
				row.names(miR_meth_fisher_mat) <- c(neg_miRs,neg_meth)
				colnames(miR_meth_fisher_mat) <- c(neg_miRs,neg_meth)
			}

		for(miR_name in neg_miRs){
			for(meth_name in neg_meth){
				tmp_mir <- as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]) > median(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),na.rm=T)
				tmp_meth <- as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)

				overlap_lst <- (tmp_mir + tmp_meth) %% 2 
				tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))
				# if(tmp$p.value < 0.05){

				miR_meth_fisher_mat[miR_name,meth_name]<- tmp$p.value
				miR_meth_fisher_mat[meth_name,miR_name]<- tmp$p.value
				# }

				# tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				# if(tmp$p.value < 0.05){
				# 	miR_meth_mat[miR_name,meth_name] <- tmp$estimate
				# 	miR_meth_mat[meth_name,miR_name] <- tmp$estimate

				# }
			}
		}

		for(miR_name in neg_miRs){
			for(miR_name2 in neg_miRs){
				tmp_mir <- as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]) > median(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),na.rm=T)
				tmp_mir2 <- as.numeric(all_miRNA_datasets[[cancer_type]][miR_name2,common_samples[[cancer_type]]]) > median(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name2,common_samples[[cancer_type]]]),na.rm=T)

				overlap_lst <- (tmp_mir + tmp_mir2) %% 2 
				tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))
				# if(tmp$p.value < 0.05){

				miR_meth_fisher_mat[miR_name,miR_name2]<- tmp$p.value
				# }
				# miR_meth_fisher_mat[meth_name,miR_name]<- p_val
			}
			# all_meth_mat_vals <- rbind(all_meth_mat_vals,all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]])
		}


		for(meth_name in neg_meth){
			for(meth_name2 in neg_meth){
				tmp_meth <- as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)
				tmp_meth2 <- as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name2,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name2,common_samples[[cancer_type]]]),na.rm=T)

				overlap_lst <- (tmp_meth + tmp_meth2) %% 2 
				tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))
				# if(tmp$p.value < 0.05){
					miR_meth_fisher_mat[meth_name,meth_name2]<- tmp$p.value
				# }
			}
			# all_meth_mat_vals <- rbind(all_meth_mat_vals,methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]])
		}
		if(!is.null(mut_names)){
			for(mut_name in mut_names){
				for(miR_name in neg_miRs){

					tmp_mir <- as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]) > median(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),na.rm=T)
					tmp_mut <- as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]])#as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)

					overlap_lst <- (tmp_mir + tmp_mut) %% 2 
					tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))

					miR_meth_fisher_mat[miR_name,mut_name]<- tmp$p.value
					miR_meth_fisher_mat[mut_name,miR_name]<- tmp$p.value

				}
				for(meth_name in neg_meth){
					
					tmp_meth <- as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)
					tmp_mut <- as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]])#as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)

					overlap_lst <- (tmp_meth + tmp_mut) %% 2 
					tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))

					miR_meth_fisher_mat[meth_name,mut_name]<- tmp$p.value
					miR_meth_fisher_mat[mut_name,meth_name]<- tmp$p.value
				}

				for(mut_name2 in mut_names){
					tmp_mut <- as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]])#as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)
					tmp_mut2 <- as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name2,common_samples[[cancer_type]]])#as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)

					overlap_lst <- (tmp_mut + tmp_mut2) %% 2 
					tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))
					# if(tmp$p.value < 0.05){
					miR_meth_fisher_mat[mut_name,mut_name2]<- tmp$p.value
					# }
				}

			}
		}

		# if('mutation_data' %in% negative_predictors){
		# 	#compute miRNA-mutatin
		# 	for(miR_name in neg_miRs){
		# 		tmp_mir <- as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]) > median(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),na.rm=T)
		# 		overlap_lst <- (tmp_mir + as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]])) %% 2 			
		# 		tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))
		# 		# if(tmp$p.value < 0.05){
		# 			miR_meth_fisher_mat[miR_name,'mutation_data']<- tmp$p.value
		# 			miR_meth_fisher_mat['mutation_data',miR_name]<-  tmp$p.value
		# 		# }
		# 	}
		# 	#compute meth-mutation
		# 	for(meth_name in neg_meth){
		# 		tmp_meth <- as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]) > median(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),na.rm=T)
		# 		overlap_lst <- (tmp_meth + as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]])) %% 2 
		# 		tmp <- (binom.test(sum(overlap_lst,na.rm=T), length(overlap_lst), p = 0.5,alternative='g'))
		# 		# if(tmp$p.value < 0.05){
		# 			miR_meth_fisher_mat[meth_name,'mutation_data']<- tmp$p.value
		# 			miR_meth_fisher_mat['mutation_data',meth_name]<- tmp$p.value
		# 		# }
		# 	}
		# 	#compute muatiosn-mutation (this is 1)
		# 	miR_meth_fisher_mat['mutation_data','mutation_data'] <- 1
		# }
			p[[cancer_type]] <- gplots::heatmap.2( miR_meth_fisher_mat ,
		                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
		                   trace = "none",
		                #   xlab = "Gene ID",
		                #   ylab="Gene ID",
		                   na.color="grey",
		                   labRow=rownames(miR_meth_fisher_mat),#rownames(tsg_sig_score_cor_by_gene[[tsg_name]]),#converted_symbols_map[rownames(sig_mut_props_hmap)],
		                   #labCol=colnames(autocors),#gene_sig,
		                   main = paste0('Autocorrelation-,',gene_name,' neg. regulators\n',cancer_type),#paste0("Partial to PTEN mut"),
		                   dendrogram = "both",
		                   breaks = seq(0,1,length=101),
		                   #symbreaks = T,
		                   Rowv = T,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.8,cexCol=0.8)
		dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/pten_neg_reg_mut_excl_pVal_',cancer_type,'.pdf'),width=12,height=12)
	  	dev.off()


		
	}
}



pten_entrez <- '5728'

miR_of_interest <- c()

for(cancer_type in cancer_types){
	# plot(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',]))
	print(cancer_type)
	if(length(colnames(all_mRNA_datasets[[cancer_type]]))>9){
		print(sum(mut_data[[cancer_type]]['PTEN',]))
		unmut_cases <- which(mut_data[[cancer_type]]['PTEN',]==0)
		if(length(unmut_cases) > 9){
			print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,unmut_cases]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',unmut_cases]),method='spearman'))
			print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,unmut_cases]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-3p',unmut_cases]),method='spearman'))
			print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,unmut_cases]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-5p',unmut_cases]),method='spearman'))
			print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,unmut_cases]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-3p',unmut_cases]),method='spearman'))
		}
		print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',]),method='spearman'))
		print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-3p',]),method='spearman'))
		print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-5p',]),method='spearman'))
		print(cor.test(as.numeric(all_mRNA_datasets[[cancer_type]][pten_entrez,]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-3p',]),method='spearman'))
		tmp_mat <- cbind(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-3p',]),
			cbind(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',]),
				cbind(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-3p',]),
					as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-5p',]))))
		colnames(tmp_mat) <- c('miR_17_3p','miR_17_5p','miR_21_3p','miR_21_5p')
		p <- ggpairs(as.data.frame(tmp_mat))
		print(p)
		# print(cor.test(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-3p',]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',]),method='spearman'))
		# print(cor.test(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-3p',]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-3p',]),method='spearman'))
		# print(cor.test(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-3p',]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-5p',]),method='spearman'))
		# print(cor.test(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-3p',]),method='spearman'))
		# print(cor.test(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-17-5p',]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-5p',]),method='spearman'))
		# print(cor.test(as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-3p',]),as.numeric(all_miRNA_datasets[[cancer_type]]['hsa-miR-21-5p',]),method='spearman'))


	}
}

#########LET'S TRY AND GET THE DISTRIBUTIONS FOR THE CORRELATION COEFFICIENTS OF MIR-MIR AND MIR-METH AND METH-METH ETC


all_miR_miR_cor_percentile_means <- list()
all_miR_meth_cor_percentile_means <- list()
all_miR_mut_cor_percentile_means <- list()
all_meth_meth_cor_percentile_means <- list()
all_meth_mut_cor_percentile_means <- list()
all_miR_miR_minus_meth_meth_cor_percentile_means <- list()
all_miR_miR_minus_miR_meth_cor_percentile_means <- list()
all_meth_meth_minus_miR_meth_cor_percentile_means <- list()

#first we can check the correlations between the individual variables that are sig in each cancer type
for(gene_name in meth_gene_names){

	all_miR_miR_cor_percentile_means[[gene_name]] <- list()
	all_miR_meth_cor_percentile_means[[gene_name]] <- list()
	all_miR_mut_cor_percentile_means[[gene_name]] <- list()
	all_meth_meth_cor_percentile_means[[gene_name]] <- list()
	all_meth_mut_cor_percentile_means[[gene_name]] <- list()
	all_miR_miR_minus_meth_meth_cor_percentile_means[[gene_name]] <- list()
	all_miR_miR_minus_miR_meth_cor_percentile_means[[gene_name]] <- list()
	all_meth_meth_minus_miR_meth_cor_percentile_means[[gene_name]] <- list()


	load(paste0(gene_name,'-analysis/overall_coeff_mat_beta_values_z_transf_sep_muts.rda'))
	p <- list()
	dir.create(paste0(gene_name,'-analysis/neg_reg_autocors'))
	dir.create(paste0(gene_name,'-analysis/neg_reg_autocors/density'))
	dev.new()
	# par(mfrow=c(3,4)) #because we have 12 cancer types with enough samples
	for(cancer_type in colnames(overall_mat)){

		negative_predictors <- rownames(overall_mat)[which(overall_mat[,cancer_type] < 0)]
		neg_miRs <- negative_predictors[which(grepl('miR', negative_predictors))]
		neg_meth <- negative_predictors[which(grepl('cg', negative_predictors))]
		mut_types <- c('Missense_Mutation','Nonsense_Mutation','Frame_Shift_Del','Frame_Shift_Ins','Splice_Site',
			'Translation_Start_Site','Nonstop_Mutation','3\'UTR','5\'UTR','3\'Flank',
			'5\'Flank','RNA','Targeted_Region','In_Frame_Del','In_Frame_Ins',
			'IGR','Intron','Silent')
		mut_names <- negative_predictors[which(negative_predictors %in% mut_types)]
		if(!is.null(mut_names)){
			miR_meth_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))
			miR_meth_mat_all_values <- matrix(0,nrow=length(neg_miRs)+length(neg_meth)+length(mut_names),ncol=length(neg_miRs)+length(neg_meth)+length(mut_names))

			row.names(miR_meth_mat) <- c(neg_miRs,neg_meth,mut_names)
			colnames(miR_meth_mat) <- c(neg_miRs,neg_meth,mut_names)

			row.names(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth,mut_names)
			colnames(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth,mut_names)
		}else{
			miR_meth_mat <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))
			miR_meth_mat_all_values <- matrix(0,nrow=length(neg_miRs)+length(neg_meth),ncol=length(neg_miRs)+length(neg_meth))

			row.names(miR_meth_mat) <- c(neg_miRs,neg_meth)
			colnames(miR_meth_mat) <- c(neg_miRs,neg_meth)

			row.names(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth)
			colnames(miR_meth_mat_all_values) <- c(neg_miRs,neg_meth)
		}
		
		all_meth_mat_vals <- c()
	 
		for(miR_name in neg_miRs){
			for(meth_name in neg_meth){
				tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[miR_name,meth_name] <- tmp$estimate
					miR_meth_mat[meth_name,miR_name] <- tmp$estimate
				}
				miR_meth_mat_all_values[miR_name,meth_name] <- tmp$estimate
				miR_meth_mat_all_values[meth_name,miR_name] <- tmp$estimate
			}
		}

		for(miR_name in neg_miRs){
			for(miR_name2 in neg_miRs){
				tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(all_miRNA_datasets[[cancer_type]][miR_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[miR_name,miR_name2] <- tmp$estimate
				}
				miR_meth_mat_all_values[miR_name,miR_name2] <- tmp$estimate

			}
			all_meth_mat_vals <- rbind(all_meth_mat_vals,all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]])
		}


		for(meth_name in neg_meth){
			for(meth_name2 in neg_meth){
				tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
				if(tmp$p.value < 0.05){
					miR_meth_mat[meth_name,meth_name2] <- tmp$estimate
				}
				miR_meth_mat_all_values[meth_name,meth_name2] <- tmp$estimate

			}
			all_meth_mat_vals <- rbind(all_meth_mat_vals,methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]])
		}
		if(!is.null(mut_names)){
			for(mut_name in mut_names){
				for(miR_name in neg_miRs){
					tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[miR_name,mut_name] <- tmp$estimate
						miR_meth_mat[mut_name,miR_name] <- tmp$estimate
					}
					miR_meth_mat_all_values[miR_name,mut_name] <- tmp$estimate
					miR_meth_mat_all_values[mut_name,miR_name] <- tmp$estimate
					
				}
				for(meth_name in neg_meth){
					tmp <- cor.test(as.numeric(methylation_data[[gene_name]][[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[meth_name,mut_name] <- tmp$estimate
						miR_meth_mat[mut_name,meth_name] <- tmp$estimate
					}
					miR_meth_mat_all_values[meth_name,mut_name] <- tmp$estimate
					miR_meth_mat_all_values[mut_name,meth_name] <- tmp$estimate
				}

				for(mut_name2 in mut_names){
					tmp <- cor.test(as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[gene_name]][[cancer_type]][mut_name2,common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
					if(tmp$p.value < 0.05){
						miR_meth_mat[mut_name,mut_name2] <- tmp$estimate
					}
					miR_meth_mat_all_values[mut_name,mut_name2] <- tmp$estimate

				}

			}
		}


		#compare the  miR-miR density to the null

		# plot(density(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T))
		all_miR_names <-  rownames(expr_mat_cors_list_all_values[[gene_name]])[which(grepl('miR', rownames(expr_mat_cors_list_all_values[[gene_name]])))]
		###lines(density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_miR_names],na.rm=T),col='red')		
		max_dens_y <- max(density(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T)$y,
			density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_miR_names],na.rm=T)$y)
		max_dens_x <- max(density(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T)$x,
			density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_miR_names],na.rm=T)$x)	
		min_dens_x <- min(density(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T)$x,
			density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_miR_names],na.rm=T)$x)		


		# plot(density(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T),xlim=c(min_dens_x,max_dens_x),ylim=c(0,max_dens_y),lwd=2,
		# 	main=paste0('miRNA-miRNA correlation plots for\n', gene_name, ' in ', cancer_type))
		# lines(density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_miR_names],na.rm=T),,xlim=c(min_dens_x,max_dens_x),ylim=c(0,max_dens_y),lwd=2,lty=2)		

		# 	#compare meth-meth density to the null
		if(length(neg_meth)>1){
		all_meth_names <-  rownames(expr_mat_cors_list_all_values[[gene_name]])[which(grepl('cg', rownames(expr_mat_cors_list_all_values[[gene_name]])))]
		max_dens_y <- max(density(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T)$y,
			density(expr_mat_cors_list_all_values[[gene_name]][all_meth_names,all_meth_names],na.rm=T)$y)
		max_dens_x <- max(density(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T)$x,
			density(expr_mat_cors_list_all_values[[gene_name]][all_meth_names,all_meth_names],na.rm=T)$x)	
		min_dens_x <- min(density(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T)$x,
			density(expr_mat_cors_list_all_values[[gene_name]][all_meth_names,all_meth_names],na.rm=T)$x)		
		# plot(density(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T),xlim=c(min_dens_x,max_dens_x),ylim=c(0,max_dens_y),lwd=2,
		# 	main=paste0('meth-meth correlation plots for\n', gene_name, ' in ', cancer_type))
		# lines(density(expr_mat_cors_list_all_values[[gene_name]][all_meth_names,all_meth_names],na.rm=T),,xlim=c(min_dens_x,max_dens_x),ylim=c(0,max_dens_y),lwd=2,lty=2)		
		}else{
			# plot(1,1)
		}
		# # #compare mir-meth density to the null
		if(length(neg_meth)>1){
		all_meth_names <-  rownames(expr_mat_cors_list_all_values[[gene_name]])[which(grepl('cg', rownames(expr_mat_cors_list_all_values[[gene_name]])))]

		max_dens_y <- max(density(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T)$y,
			density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_meth_names],na.rm=T)$y)
		max_dens_x <- max(density(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T)$x,
			density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_meth_names],na.rm=T)$x)	
		min_dens_x <- min(density(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T)$x,
			density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_meth_names],na.rm=T)$x)		
		# plot(density(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T),xlim=c(min_dens_x,max_dens_x),ylim=c(0,max_dens_y),lwd=2,
		# 	main=paste0('miRNA-meth correlation plots for\n', gene_name, ' in ', cancer_type))
		# lines(density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_meth_names],na.rm=T),,xlim=c(min_dens_x,max_dens_x),ylim=c(0,max_dens_y),lwd=2,lty=2)		
		}else{
			# plot(1,1)
		}

		#now that we have made all of the density plots that we are concerned with, let's 
		#calculate the p-val for the mean autocor being different than that what we observe, given the underling distro
		#let's randomly sample the same length of miRs and see what the means are of the autocors
		#then check of these, how many are above the one we see....


		null_distr_values_comp_means <- get_null_distr_p_value_compare_means(length(neg_meth),length(neg_miRs),mut_names,all_expr_mat_cors_all_values[[cancer_type]][[gene_name]])
		# null_distr_values_all_values <- get_null_distr_p_value_mir_meth_mat(length(neg_meth),length(neg_miRs),mut_names,all_expr_mat_cors_all_values[[cancer_type]][[gene_name]])


		all_miR_miR_cor_percentile_means[[gene_name]][[cancer_type]] <- ecdf(null_distr_values_comp_means[[1]])(mean(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T))
		if(length(neg_meth) > 0){
			all_miR_meth_cor_percentile_means[[gene_name]][[cancer_type]] <- ecdf(null_distr_values_comp_means[[2]])(mean(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T))
			all_meth_meth_cor_percentile_means[[gene_name]][[cancer_type]] <- ecdf(null_distr_values_comp_means[[3]])(mean(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T))
		}
		# all_miR_mut_cor_percentile_means[[gene_name]][[cancer_type]] <- list()
		# all_meth_mut_cor_percentile_means[[gene_name]][[cancer_type]] <- list()

		# if(!is.null(mut_names)){
		# 	for(mut_name in mut_names){
		# 			all_miR_mut_cor_percentile_means[[gene_name]][[cancer_type]][[mut_name]] <- ecdf(null_distr_values[[4]][[mut_name]])(mean(miR_meth_mat_all_values[mut_name,neg_miRs],na.rm=T))
				
		# 		if(length(neg_meth) > 0){
		# 			all_meth_mut_cor_percentile_means[[gene_name]][[cancer_type]][[mut_name]] <- ecdf(null_distr_values[[5]][[mut_name]])(mean(miR_meth_mat_all_values[mut_name,neg_meth],na.rm=T))
		# 		}
		# 	}
		# }
		if(length(neg_meth) > 0){

			all_miR_miR_minus_meth_meth_cor_percentile_means[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[1]] - null_distr_values[[3]])(mean(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T) - mean(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T))
			all_miR_miR_minus_miR_meth_cor_percentile_means[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[1]] - null_distr_values[[2]])(mean(miR_meth_mat_all_values[neg_miRs,neg_miRs],na.rm=T) - mean(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T))
			all_meth_meth_minus_miR_meth_cor_percentile_means[[gene_name]][[cancer_type]] <- ecdf(null_distr_values[[3]] - null_distr_values[[2]])(mean(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T) - mean(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T))
			print(all_miR_miR_minus_miR_meth_cor_percentile_means[[gene_name]][[cancer_type]] )#ecdf(null_distr_values[[1]])(tot_prob_co_miRs))

		}




		

		# plot(density(miR_meth_mat_all_values[neg_meth,neg_meth],na.rm=T))
		# lines(density(expr_mat_cors_list_all_values[[gene_name]][all_meth_names,all_meth_names],na.rm=T),col='red')	


		# plot(density(miR_meth_mat_all_values[neg_miRs,neg_meth],na.rm=T))
		# lines(density(expr_mat_cors_list_all_values[[gene_name]][all_miR_names,all_meth_names],na.rm=T),col='red')	

		# if('mutation_data' %in% negative_predictors){
		# 	#compute miRNA-mutatin
		# 	for(miR_name in neg_miRs){
		# 		tmp <- cor.test(as.numeric(all_miRNA_datasets[[cancer_type]][miR_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
		# 		if(tmp$p.value < 0.05){
		# 			miR_meth_mat[miR_name,'mutation_data'] <- tmp$estimate
		# 			miR_meth_mat['mutation_data',miR_name] <- tmp$estimate
		# 		}
		# 	}
		# 	#compute meth-mutation
		# 	for(meth_name in neg_meth){
		# 		tmp <- cor.test(as.numeric(methylation_data[[cancer_type]][meth_name,common_samples[[cancer_type]]]),as.numeric(mut_data[[cancer_type]]['PTEN',common_samples[[cancer_type]]]),method='spearman',use='pairwise.complete.obs')
		# 		if(tmp$p.value < 0.05){
		# 			miR_meth_mat[meth_name,'mutation_data'] <- tmp$estimate
		# 			miR_meth_mat['mutation_data',meth_name] <- tmp$estimate
		# 		}
		# 	}
		# 	#compute muatiosn-mutation (this is 1)
		# 	miR_meth_mat['mutation_data','mutation_data'] <- 1
		# }

		# row.names(all_meth_mat_vals) <-  gsub('[-]','.',row.names(all_meth_mat_vals) )

		
	}
	# dev.copy(pdf,paste0(gene_name,'-analysis/neg_reg_autocors/',gene_name,'_density_plots_vs_null_distro_miR_meth.pdf'),width=12,height=9)
	#   	dev.off()

}





#PLOTTING COMMANDS FOR THE EXLCUSIVITY ANALYSIS
dev.new()
melted_plotting <- melt(all_miR_miR_minus_miR_meth_cor_percentile_means)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, miR-miR proportion cor. mean \nminus miR-meth proportion cor. mean')
dev.copy(pdf,paste0('eCDF_comparison_miR_miR_minus_miR_meth_mean.pdf'),width=12,height=9)
dev.off()

dev.new()
melted_plotting <- melt(all_miR_miR_minus_meth_meth_cor_percentile_means)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, miR-miR proportion cor. mean \nminus meth-meth proportion cor. mean')
dev.copy(pdf,paste0('eCDF_comparison_miR_miR_minus_meth_meth_mean.pdf'),width=12,height=9)
dev.off()

dev.new()
melted_plotting <- melt(all_meth_meth_minus_miR_meth_cor_percentile_means)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, meth-meth proportion cor. mean \nminus miR-meth proportion cor. mean')
dev.copy(pdf,paste0('eCDF_comparison_meth_meth_minus_miR_meth_mean.pdf'),width=12,height=9)
dev.off()

dev.new()
melted_plotting <- melt(all_miR_miR_cor_percentile_means)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, miR-miR proportion cor. mean')
dev.copy(pdf,paste0('eCDF_comparison_miR_miR_mean.pdf'),width=12,height=9)
dev.off()

dev.new()
melted_plotting <- melt(all_miR_meth_cor_percentile_means)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, miR-meth proportion cor. mean')
dev.copy(pdf,paste0('eCDF_comparison_miR_meth_mean.pdf'),width=12,height=9)
dev.off()

dev.new()
melted_plotting <- melt(all_meth_meth_cor_percentile_means)
ggplot(as.data.frame(melted_plotting),aes(x=L1,y=value)) +  geom_boxplot() + geom_point(aes(color=L2)) + theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
labs(color='Cancer type',x='Gene name',y='Empiric CDF percentile',title='eCDF, meth-meth proportion cor. mean')
dev.copy(pdf,paste0('eCDF_comparison_meth_meth_mean.pdf'),width=12,height=9)
dev.off()

#_------- density plots of genes of interest
genes_of_interest <- c('PTEN','ACVR2A','ARHGEF12','CDK12','DNMT3A','FAT4','SFRP4','TGFBR2')
entrez_gene_names <- list()
xx <- as.list(org.Hs.egALIAS2EG)
# Remove pathway identifiers that do not map to any entrez gene id
xx <- xx[!is.na(xx)]
for(gene_name in genes_of_interest){
	entrez_gene_names[[gene_name]] <- xx[[gene_name]]
}


for (cancer_type in cancer_types){

	for(gene in genes_of_interest){
		dir.create(paste0('../miRNA_hallmarks/',gene,'-analysis/expr_plots'))
		dev.new()
		plot(density(as.numeric(all_mRNA_datasets[[cancer_type]][entrez_gene_names[[gene]],]),na.rm=T),main=paste0(gene, ' expression in ', cancer_type))
		dev.copy(pdf,paste0('../miRNA_hallmarks/',gene,'-analysis/expr_plots/',cancer_type,'_',gene,'.pdf'),width=6,height=6)
		dev.off()
	}
	graphics.off()
}
