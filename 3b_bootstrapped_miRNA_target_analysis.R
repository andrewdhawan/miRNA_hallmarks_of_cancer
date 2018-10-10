#bootstrapped miRNA target analysis


bootstrapped_targets_list <- list()

all_miRNA_up <- c()
# targets_miRNA_up <- list()
# targets_miRNA_down <- list()
for (sig_name in all_sigs){
	print(sig_name)
	# targets_miRNA_down[[sig_name]] <- list()
	# for(miRNA_name in sig_miRNA_down[[sig_name]]){
	# 	targets_miRNA_down[[sig_name]][[miRNA_name]] <- tryCatch(rownames(getPredictedTargets(miRNA_name,species='hsa', method='geom',min_src=2)),error=function(e){NA})
	# 	#print(targets_miRNA_down[[sig_name]][[miRNA_name]])
	# }
	# targets_miRNA_up[[sig_name]] <- list()
	all_miRNA_up <- c(all_miRNA_up, sig_miRNA_up[[sig_name]])
	# for(miRNA_name in sig_miRNA_up[[sig_name]]){
	# 		targets_miRNA_up[[sig_name]][[miRNA_name]]  <- tryCatch(rownames(getPredictedTargets(miRNA_name,species='hsa', method='geom',min_src=2)),error=function(e){NA})
	# }
}
all_miRNA_up <- unique(all_miRNA_up)

length_bootstrapping <- length(all_miRNA_up)


all_miRNA <- rownames(all_miRNA_datasets$BRCA)
all_miRNA_targets <- list()
count <- 1
for(miR_name in all_miRNA){
	all_miRNA_targets[[miR_name]]  <- tryCatch(rownames(getPredictedTargets(miR_name,
		species='hsa', method='geom',min_src=2)),error=function(e){NA})
	count <- count + 1
	print(count)
}
save(file='all_miRNA_targets.rda',all_miRNA_targets)

#---load the data----

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
mut_data <- list()
for(cancer_type in cancer_types){
	#load mutation data
	fName_mut <- paste0('../Reprocessed GDAC data/',cancer_type,'/mutation/mutations.txt')
	mut_data[[cancer_type]] <- read.table(fName_mut, sep='\t',stringsAsFactors = FALSE, header=TRUE,quote="")
	colnames(mut_data[[cancer_type]]) <- gsub('[.]','-',colnames(mut_data[[cancer_type]]))
	mut_data[[cancer_type]] <- ((mut_data[[cancer_type]] > 0) & (mut_data[[cancer_type]] < 16)) * 1

}

#we only want the common ones of all 3 types of data:

for(cancer_type in cancer_types){
	common_samples <- intersect(colnames(mut_data[[cancer_type]]),intersect(colnames(all_miRNA_datasets[[cancer_type]]),colnames(all_mRNA_datasets[[cancer_type]])))
	mut_data[[cancer_type]] <- mut_data[[cancer_type]][,common_samples]
	all_miRNA_datasets[[cancer_type]] <- all_miRNA_datasets[[cancer_type]][,common_samples]
	all_mRNA_datasets[[cancer_type]] <- all_mRNA_datasets[[cancer_type]][,common_samples]
	print(cancer_type)
	print(length(common_samples))
}

#--------


# save(file='miRNA_up_down_targets.rda',targets_miRNA_up,targets_miRNA_down)

all_targets <- c()
for(miR_name in rownames(all_miRNA_datasets[['BRCA']])){
	all_targets<- c(all_targets,all_miRNA_targets[[miR_name]])
}
all_targets <- unique(all_targets)

#to make things more efficient, let's pre-compute the correlation coefficients of EVERY miRNA with EVERY mRNA target across EVERY cancer type
all_miRNA_mRNA_cor_vals_sig_subset_all_genes <- list()

for (cancer_type in cancer_types){
	all_miRNA_mRNA_cor_vals_sig_subset_all_genes[[cancer_type]] <- matrix(NA, nrow=length(rownames(all_miRNA_datasets[[cancer_type]])),ncol=length(rownames(all_mRNA_datasets[[cancer_type]])))
	row.names(all_miRNA_mRNA_cor_vals_sig_subset_all_genes[[cancer_type]]) <- rownames(all_miRNA_datasets[[cancer_type]])
	colnames(all_miRNA_mRNA_cor_vals_sig_subset_all_genes[[cancer_type]]) <- rownames(all_mRNA_datasets[[cancer_type]])
	print(cancer_type)
	for(miRNA_name in rownames(all_miRNA_datasets[[cancer_type]])){
		print(miRNA_name)
		for(gene_name in all_miRNA_targets[[miRNA_name]]){#rownames(all_mRNA_datasets[[cancer_type]])){
			# print(gene_name)
			if(!is.null(converted_symbols[[gene_name]])){
				# start_time <- Sys.time()

				if(converted_symbols[[gene_name]] %in% rownames(mut_data[[cancer_type]])){	
					# end_time <- Sys.time()
					# print(paste0('time1 : ', (end_time - start_time)))	
					com_non_na_vals <- intersect(which(!is.na(all_miRNA_datasets[[cancer_type]][miRNA_name,])),intersect(which(!is.na(all_mRNA_datasets[[cancer_type]][gene_name,])), which(!is.na(mut_data[[cancer_type]][converted_symbols[[gene_name]],]))))
					#ok now we need to know the miRNA-and-all genes correlation distribution partial to the mutation status
						#now we also need to know the miRNA and all targets correlation distribution partial to the mutation status but we have that from above
					if(length(com_non_na_vals) > 10){

						tryCatch({x <- pcor.test(x=as.numeric(all_miRNA_datasets[[cancer_type]][miRNA_name,com_non_na_vals]),
							y=as.numeric(all_mRNA_datasets[[cancer_type]][gene_name,com_non_na_vals]),
							z=as.numeric(mut_data[[cancer_type]][converted_symbols[[gene_name]],com_non_na_vals]),method='spearman')
					
							if(!is.na(x)){
								if(x$p.val < 0.05){
									all_miRNA_mRNA_cor_vals_sig_subset_all_genes[[cancer_type]][miRNA_name,gene_name] <- x$estimate
								}else{
									all_miRNA_mRNA_cor_vals_sig_subset_all_genes[[cancer_type]][miRNA_name,gene_name] <- 0
								}
							}

					# miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]] <- c(x$estimate,x$p.val)
						},error=function(err){print(err)})
						# end_time <- Sys.time()
						# print(paste0('time2 : ', (end_time - start_time)))	
					}
				}
			}
		}
	}
}



#load in list of tumour suppressors
cosmic_tsg <- read.table('cosmic_data/cosmic_tsg.txt',stringsAsFactors=F,quote="",header=F)


N_bootstrap_runs <- 1000
library(RankProd)
tsg_in_sig_list <- rep(0,times=N_bootstrap_runs)

for(i in 1:N_bootstrap_runs){

	all_up_miRNA <- sample(names(all_miRNA_targets),length_bootstrapping,replace=F) #these are the miRNA we are using in the resampling procedure.
	miR_targets_names <- c()

	for(miRNA_name in all_up_miRNA){
		miR_targets_names <-c(miR_targets_names, paste0(miRNA_name,'/',all_miRNA_targets[[miRNA_name]]))
	}

	miR_mRNA_cor_vals_bootstrapped <- matrix(NA,nrow=length(miR_targets_names),ncol=length(cancer_types))
	row.names(miR_mRNA_cor_vals_bootstrapped) <- miR_targets_names
	colnames(miR_mRNA_cor_vals_bootstrapped) <- cancer_types

	for(cancer_type in cancer_types){
		for(miRNA_name in all_up_miRNA){
			for(gene_name in all_miRNA_targets[[miRNA_name]]){
				miR_mRNA_cor_vals_bootstrapped[paste0(miRNA_name,'/',gene_name),cancer_type] <- all_miRNA_mRNA_cor_vals_sig_subset_all_genes[[cancer_type]][miRNA_name,gene_name]
			}
		}
	}

	#now we're going to have to remove any rows or columns that are all na values

	good_cols <- c()
	for(cancer_type in cancer_types){
		if(sum(is.na(miR_mRNA_cor_vals_bootstrapped[,cancer_type])) != length(miR_mRNA_cor_vals_bootstrapped[,cancer_type])){
			good_cols <- c(good_cols,cancer_type)
		}
	}

	miR_mRNA_cor_vals_bootstrapped <- miR_mRNA_cor_vals_bootstrapped[,good_cols]


	good_rows <- which(rowSums(is.na(miR_mRNA_cor_vals_bootstrapped)) <= (length(colnames(miR_mRNA_cor_vals_bootstrapped)) - 5))
	miR_mRNA_cor_vals_bootstrapped <- miR_mRNA_cor_vals_bootstrapped[good_rows,]

	ranked_cors <- RP(miR_mRNA_cor_vals_bootstrapped,cl=rep(1,length(colnames(miR_mRNA_cor_vals_bootstrapped))))
	ranked_cors_table <- topGene(ranked_cors,method='pfp',cutoff=0.05,gene.names = rownames(miR_mRNA_cor_vals_bootstrapped))
	ranked_cors_table <- ranked_cors_table$Table1

	pair_names <- strsplit(rownames(ranked_cors_table),split = '/')
	pair_names <- melt(pair_names)
	pair_names <- cbind(as.character(pair_names$value[seq(from=1,to=length(pair_names$value),by=2)]),
		as.character(pair_names$value[seq(from=2,to=length(pair_names$value),by=2)]))

	# cosmic_oncogenes <- read.table('cosmic_data/cosmic_oncogenes.txt',stringsAsFactors=F,quote="",header=F)
	#object to convert gene entrez id to gene names
	# converted_symbols_map <- melt(converted_symbols)
	# tmp <- converted_symbols_map$L1
	# converted_symbols_map <- converted_symbols_map[,1]
	# names(converted_symbols_map) <- tmp
	#convert the entrez IDs to symbols
	converted_gene_names <- converted_symbols_map[pair_names[,2]]
	common_tsg <- intersect(unique(converted_gene_names),cosmic_tsg[,1])
	tsg_in_sig_list[i] <- common_tsg
}

