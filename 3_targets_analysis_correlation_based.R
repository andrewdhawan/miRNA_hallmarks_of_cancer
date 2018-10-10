#this is the code to do analysis of the targets of the miRNA that we identified as significant earlier.

library(org.Hs.eg.db)
library(miRNAtap)
library(ppcor)
library(RankProd)
library(gplots)
library(reshape2)

#at some point we are going to have to change to the Gene symbol for identification because this is how mutatoin data is recorded
## Bimap interface:
x <- org.Hs.egSYMBOL
# Get the gene symbol that are mapped to an entrez gene identifiers
mapped_genes <- mappedkeys(x)
# Convert to a list
converted_symbols <- as.list(x[mapped_genes])

#check the expression of miRNA and mRNA targets

#first we need to load in the data
cancer_types <- c('BRCA','UCEC','HNSC','KIRC','LUAD','THCA','PRAD','LUSC','OV','STAD','BLCA','COAD','LIHC','CESC','KIRP')


all_sigs <- c('Hallmark: Epithelial Mesenchymal Transition','Invasiveness, Marsan 2014',
	'Hallmark: Oxidative Phosphorylation','Hallmark: Reactive Oxygen Species Pathway',
	'Hallmark: G2M Checkpoint',
	'Hallmark: PI3K AKT MTOR Signaling','Hallmark: Xenobiotic Metabolism',
	'Hallmark: DNA Repair','Hallmark: p53 Pathway',
	'Hypoxia, Buffa 2010','Hallmark: Angiogenesis','Hallmark: Hypoxia','Angiogenesis, Desmedt 2008','Angiogenesis, Masiero 2013',
	'Hallmark: Apoptosis','Apoptosis, Desmedt 2008',
	'Proliferation, Desmedt 2008','Hallmark: KRAS Signaling Up',
	'Hallmark: Inflammatory Response','Hallmark: IL2 STAT5 Signaling','Hallmark: IL6 JAK STAT3 Signaling','Hallmark: TGF Beta Signaling','Hallmark: TNFa Signaling via NFKB','Immune, Desmedt 2008')


#we need to load in all of the signature names here:
sig_fnames_list <- c('HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION.txt','invasiveness_gene_sig_entrez_marsan2014.txt',
	'HALLMARK_OXIDATIVE_PHOSPHORYLATION.txt','HALLMARK_REACTIVE_OXIGEN_SPECIES_PATHWAY.txt',
	'HALLMARK_G2M_CHECKPOINT.txt',
	'HALLMARK_PI3K_AKT_MTOR_SIGNALING.txt','HALLMARK_XENOBIOTIC_METABOLISM.txt',
	'HALLMARK_DNA_REPAIR.txt','HALLMARK_P53_PATHWAY.txt',
	'hypoxia_gene_sig_entrez_probes.txt','HALLMARK_ANGIOGENESIS.txt','HALLMARK_HYPOXIA.txt','angiogenesis_gene_sig_entrez_desmedt2008_pos.txt','Masiero2013angiogenesisENTREZ.txt',
	'HALLMARK_APOPTOSIS.txt','apoptosis_gene_sig_entrez_desmedt2008_pos.txt',
	'proliferation_gene_sig_entrez_desmedt2008_pos.txt','HALLMARK_KRAS_SIGNALING_UP.txt',
	'HALLMARK_INFLAMMATORY_RESPONSE.txt','HALLMARK_IL2_STAT5_SIGNALING.txt','HALLMARK_IL6_JAK_STAT3_SIGNALING.txt','HALLMARK_TGF_BETA_SIGNALING.txt','HALLMARK_TNFA_SIGNALING_VIA_NFKB.txt','immune_gene_sig_entrez_desmedt2008_pos.txt')
names(sig_fnames_list) <- all_sigs

sigs_list_by_name <- list();
for(sig_name in all_sigs){
	fname <- sig_fnames_list[sig_name]
	genes = read.csv(paste0('gene_signatures/',fname), header=F, stringsAsFactors=F, colClasses = "character")
	# print(genes)
	sigs_list_by_name[[sig_name]]<- genes

}

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

load('rank_prod_tables_out_pre_filtered.rda')
sig_miRNA_up <- list()
sig_miRNA_down <- list()
for (sig_name in all_sigs){
	sig_miRNA_up[[sig_name]] <- rownames(rank_prod_tables[[sig_name]]$Table2)
	sig_miRNA_down[[sig_name]] <- rownames(rank_prod_tables[[sig_name]]$Table1)
}
all_up_miRNA <- unique(melt(sig_miRNA_up)$value)
all_down_miRNA <- unique(melt(sig_miRNA_down)$value)

load('miRNA_up_down_targets.rda') #load in the  database of signatures, up miRs and targets

melted_targets_up <- melt(targets_miRNA_up)

miR_targets_cor_and_pval <- list()


for(miRNA_name in all_up_miRNA){
	print(miRNA_name)
	miR_targets_cor_and_pval[[miRNA_name]] <- list()
	targets_list <- unique(melted_targets_up[which(melted_targets_up[,2] == miRNA_name),1])
	for(gene_name in targets_list){
		miR_targets_cor_and_pval[[miRNA_name]][[gene_name]] <- list()
		for(cancer_type in cancer_types){
		#	miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]] <- list()
				if(!is.null(converted_symbols[[gene_name]])){
					if(converted_symbols[[gene_name]] %in% rownames(mut_data[[cancer_type]])){
					com_non_na_vals <- intersect(which(!is.na(all_miRNA_datasets[[cancer_type]][miRNA_name,])),intersect(which(!is.na(all_mRNA_datasets[[cancer_type]][gene_name,])), which(!is.na(mut_data[[cancer_type]][converted_symbols[[gene_name]],]))))
						#now we need to know the miRNA-and-all genes correlation distribution partial to the mutation status
						#now we also need to know the miRNA and all targets correlation distribution partial to the mutation status but we have that from above
					if(length(com_non_na_vals) > 10){
						tryCatch({x <- pcor.test(x=as.numeric(all_miRNA_datasets[[cancer_type]][miRNA_name,com_non_na_vals]),
							y=as.numeric(all_mRNA_datasets[[cancer_type]][gene_name,com_non_na_vals]),
							z=as.numeric(mut_data[[cancer_type]][converted_symbols[[gene_name]],com_non_na_vals]),method='spearman')
					miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]] <- c(x$estimate,x$p.val)
					},error=function(err){print(err)})
					}
				}
			}
		}
	}
}

save(file='miR_target_cor_vals.rda',miR_targets_cor_and_pval)

#take the significnat interactions, then do rank product on them across cancer types?
#should get the most conserved miRNA_mRNA repression across cancer types

#let's make a matrix of the miRNA/mRNA targets as rows and the cancer types as columns

all_miRNA_mRNA_pairs <- c()

for(miRNA_name in all_up_miRNA){
	print(miRNA_name)
	targets_names <- names(miR_targets_cor_and_pval[[miRNA_name]]) 
	all_miRNA_mRNA_pairs <- c(all_miRNA_mRNA_pairs,paste0(miRNA_name,'/',targets_names))
}
all_miRNA_mRNA_pairs <- unique(all_miRNA_mRNA_pairs)

all_miRNA_mRNA_cor_vals <- matrix(NA,nrow=length(all_miRNA_mRNA_pairs),ncol=length(cancer_types))
row.names(all_miRNA_mRNA_cor_vals) <- all_miRNA_mRNA_pairs
colnames(all_miRNA_mRNA_cor_vals) <- cancer_types
for(miRNA_name in all_up_miRNA){
	for(gene_name in names(miR_targets_cor_and_pval[[miRNA_name]])){
		for(cancer_type in cancer_types){
			if(length(miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]]) > 0) {
				all_miRNA_mRNA_cor_vals[paste0(miRNA_name,'/',gene_name),cancer_type] <- miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]][1]
			}
		}
	}
}

save(file='all_miRNA_mRNA_cor_vals.rda',all_miRNA_mRNA_cor_vals)


#now we're going to have to remove any rows or columns that are all na values

good_cols <- c()
for(cancer_type in cancer_types){
	if(sum(is.na(all_miRNA_mRNA_cor_vals[,cancer_type])) != length(all_miRNA_mRNA_cor_vals[,cancer_type])){
		good_cols <- c(good_cols,cancer_type)
	}
}

all_miRNA_mRNA_cor_vals <- all_miRNA_mRNA_cor_vals[,good_cols]


good_rows <- c()
for(miR_mRNA_name in rownames(all_miRNA_mRNA_cor_vals)){
	if(sum(is.na(all_miRNA_mRNA_cor_vals[miR_mRNA_name,])) != length(all_miRNA_mRNA_cor_vals[miR_mRNA_name,])){
		good_rows <- c(good_rows,miR_mRNA_name)
	}
}

all_miRNA_mRNA_cor_vals <- all_miRNA_mRNA_cor_vals[good_rows,]

#-----now retry the above but with a p-value cutoff as well, so that things that are non-sig are omitted from RP calc
all_p_vals <- matrix(NA,nrow=2000000,ncol=1)
count <- 1
for(miRNA_name in all_up_miRNA){
	print(miRNA_name)
	for(gene_name in names(miR_targets_cor_and_pval[[miRNA_name]])){
		for(cancer_type in cancer_types){
			if(length(miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]]) > 0) {
				all_p_vals[count,1] <- miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]][2]
				count <- count + 1
				#all_p_vals <- c(all_p_vals,miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]][2] )
				
			}
		}
	}
}




all_miRNA_mRNA_cor_vals_sig_subset <- matrix(NA,nrow=length(all_miRNA_mRNA_pairs),ncol=length(cancer_types))
row.names(all_miRNA_mRNA_cor_vals_sig_subset) <- all_miRNA_mRNA_pairs
colnames(all_miRNA_mRNA_cor_vals_sig_subset) <- cancer_types


for(miRNA_name in all_up_miRNA){
	print(miRNA_name)
	for(gene_name in names(miR_targets_cor_and_pval[[miRNA_name]])){
		for(cancer_type in cancer_types){
			if(length(miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]]) > 0) {
				if(!is.nan(miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]][2])){
					if(miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]][2] < 0.05){
						all_miRNA_mRNA_cor_vals_sig_subset[paste0(miRNA_name,'/',gene_name),cancer_type] <- miR_targets_cor_and_pval[[miRNA_name]][[gene_name]][[cancer_type]][1]
					}
				}
			}
		}
	}
}

#now we're going to have to remove any rows or columns that are all na values

good_cols <- c()
for(cancer_type in cancer_types){
	if(sum(is.na(all_miRNA_mRNA_cor_vals_sig_subset[,cancer_type])) != length(all_miRNA_mRNA_cor_vals_sig_subset[,cancer_type])){
		good_cols <- c(good_cols,cancer_type)
	}
}

all_miRNA_mRNA_cor_vals_sig_subset <- all_miRNA_mRNA_cor_vals_sig_subset[,good_cols]


good_rows <- which(rowSums(is.na(all_miRNA_mRNA_cor_vals_sig_subset)) <= (length(colnames(all_miRNA_mRNA_cor_vals_sig_subset)) - 5))
all_miRNA_mRNA_cor_vals_sig_subset <- all_miRNA_mRNA_cor_vals_sig_subset[good_rows,]


save(file='all_miRNA_mRNA_cor_vals_sig_subset.rda',all_miRNA_mRNA_cor_vals_sig_subset)
load('all_miRNA_mRNA_cor_vals_sig_subset.rda')

library(RankProd)
ranked_cors <- RP(all_miRNA_mRNA_cor_vals_sig_subset,cl=rep(1,length(colnames(all_miRNA_mRNA_cor_vals_sig_subset))))
ranked_cors_table <- topGene(ranked_cors,method='pfp',cutoff=0.05,gene.names = rownames(all_miRNA_mRNA_cor_vals_sig_subset))
ranked_cors_table <- ranked_cors_table$Table1
save(file='ranked_miRNA_mRNA_cors_table.rda',ranked_cors_table)

pair_names <- strsplit(rownames(ranked_cors_table),split = '/')
pair_names <- melt(pair_names)
pair_names <- cbind(as.character(pair_names$value[seq(from=1,to=length(pair_names$value),by=2)]),
	as.character(pair_names$value[seq(from=2,to=length(pair_names$value),by=2)]))

#load in list of tumour suppressors
cosmic_tsg <- read.table('cosmic_data/cosmic_tsg.txt',stringsAsFactors=F,quote="",header=F)
cosmic_oncogenes <- read.table('cosmic_data/cosmic_oncogenes.txt',stringsAsFactors=F,quote="",header=F)
#object to convert gene entrez id to gene names
converted_symbols_map <- melt(converted_symbols)
tmp <- converted_symbols_map$L1
converted_symbols_map <- converted_symbols_map[,1]
names(converted_symbols_map) <- tmp
#convert the entrez IDs to symbols
converted_gene_names <- converted_symbols_map[pair_names[,2]]
common_tsg <- intersect(unique(converted_gene_names),cosmic_tsg[,1])

#now that we have this list, we should check a number of things:
common_tsg_heatmap <- matrix(0,nrow=1000,ncol=length(colnames(all_miRNA_mRNA_cor_vals_sig_subset)))
colnames(common_tsg_heatmap) <- colnames(all_miRNA_mRNA_cor_vals_sig_subset)
row_count <-1
row_labels <- c()
for(tsg_name in common_tsg){
	row_labels <- c(row_labels," ")
	row_count <- row_count + 1
	num_rows <- length(which(converted_gene_names==tsg_name))
	common_tsg_heatmap[row_count:(row_count + num_rows - 1),] <- all_miRNA_mRNA_cor_vals_sig_subset[paste0(pair_names[which(converted_gene_names==tsg_name),1],'/',pair_names[which(converted_gene_names==tsg_name),2]),]
	row_labels <- c(row_labels,paste0(pair_names[which(converted_gene_names==tsg_name),1],'/',tsg_name))
	row_count <- row_count + num_rows

	# row_labels <- c(row_labels," ")

	# row_count <- row_count + 1

}

common_tsg_heatmap <- common_tsg_heatmap[1:row_count,]
common_tsg_heatmap[is.na(common_tsg_heatmap)] <- 0

p <- gplots::heatmap.2( common_tsg_heatmap ,
	                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
	                   trace = "none",
	                #   xlab = "Gene ID",
	                #   ylab="Gene ID",
	                   na.color="grey",
	                   labRow=row_labels,#rownames(tsg_sig_score_cor_by_gene[[tsg_name]]),#converted_symbols_map[rownames(sig_mut_props_hmap)],
	                   #labCol=colnames(autocors),#gene_sig,
	                   main = 'miRNA/mRNA part. cor. for TSG',#paste0("Partial to PTEN mut"),
	                   dendrogram = "col",
	                   breaks = seq(-1,1,length=101),
	                   #symbreaks = T,
	                   Rowv = F,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.4,cexCol=0.9)

dev.copy(pdf,paste0('miRNA_mRNA_part_cor_for_tsg.pdf'),width=12,height=12)
  dev.off()


#here, let's check the overrepresentation of particular pairs of miRNA./TSG in the matrix!


all_miRNA_mRNA_pairs <- rownames(all_miRNA_mRNA_cor_vals)
pairs <- strsplit(all_miRNA_mRNA_pairs,split='/')

#let's restrict the list pairs to just the TSG to see how many of those there are...

good_pairs <- c()
for(i in 1:length(pairs)){
	if(converted_symbols_map[pairs[[i]][2]] %in% cosmic_tsg[,1]){
		good_pairs <- c(good_pairs,paste0(pairs[[i]][1],'/',converted_symbols_map[pairs[[i]][2]]))
	}
}

significant_TSG <- c()
significant_TSG_pval <- c()

for(tsg_name in common_tsg){

	#tsg_name <- "FAT4"
	cur_row_labels <- row_labels[grepl(tsg_name,row_labels)]
	contingency_table <- matrix(0,nrow=2,ncol=2)
	contingency_table[1,1] <- length(cur_row_labels)
	contingency_table[1,2] <- length(which(nchar(row_labels) > 3)) - length(cur_row_labels)
	contingency_table[2,1] <- length(which(grepl(tsg_name,good_pairs))) - length(cur_row_labels)
	contingency_table[2,2] <-  length(good_pairs) - contingency_table[1,1] - contingency_table[1,2] - contingency_table[2,1]
	fisher_test <- fisher.test(contingency_table)
	if(fisher_test$p.value < 0.05){
		print(tsg_name)
		print(contingency_table)
		print(fisher_test)
		significant_TSG <- c(significant_TSG,tsg_name)
		significant_TSG_pval <- c(significant_TSG_pval,fisher_test$p.value)

	}
#	fisher.test()
}

names(significant_TSG_pval) <- significant_TSG


good_row_labels <- c()
for(tsg_name in significant_TSG){
	good_row_labels <- c(good_row_labels,row_labels[grepl(tsg_name,row_labels)])
}
save(good_row_labels,file='highly_sig_tsg_miRNA_rownames.rda')

splt_names <- melt(strsplit(good_row_labels,split='/'))
splt_names_rows <- splt_names[seq(from=1,to=length(splt_names[,1]),by=2),1]
splt_names_cols <- splt_names[seq(from=2,to=length(splt_names[,1]),by=2),1]

chord_mat <- matrix(0, nrow=length(unique(splt_names_rows)),ncol=length(unique(splt_names_cols)))
row.names(chord_mat) <- unique(splt_names_rows)
colnames(chord_mat) <- unique(splt_names_cols)
for(i in 1:length(splt_names_rows)){
	chord_mat[as.character(splt_names_rows[i]),as.character(splt_names_cols[i])] <- 1
}


library(GOplot)
GOChord(chord_mat, space = 0.02,  gene.space = 0.25, gene.size = 2,border.size = 0,nlfc = 0,process.label=10)
dev.copy(pdf,paste0('sig_TSG_miR_mRNA_chordPlot.pdf'),width=11,height=12)
dev.off()






# dir.create('sig_TSG_miR_mRNA_cor')
# #can check sig by sig what things are there,

dir.create('sig_all_genes_miR_mRNA_cor')
#can check sig by sig what things are there,
for(sig_name in all_sigs){
	print(sig_name)
#print(pair_names[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]]),])
	com_genes <- intersect(unique(converted_gene_names[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]])]),cosmic_tsg[,1])
	all_combos <- c()
	for(miR_name in sig_miRNA_up[[sig_name]]){
		all_combos <- c(all_combos,paste0(miR_name,'/',names(converted_gene_names)[which(converted_gene_names %in% com_genes)]))
	}
	rows_to_show <- intersect(all_combos,rownames(all_miRNA_mRNA_cor_vals_sig_subset))
	splt_rows <- melt(strsplit(rows_to_show,split="/"))
	splt_rows_genes <- converted_gene_names[splt_rows[seq(from=2,to=length(splt_rows[,1]),by=2),1]]
	splt_rows_mirs <- splt_rows[seq(from=1,to=length(splt_rows[,1]),by=2),1]

	plotting_mat <- all_miRNA_mRNA_cor_vals_sig_subset[rows_to_show,]
	plotting_mat[is.na(plotting_mat)] <- 0
	p <- gplots::heatmap.2(  plotting_mat,
	                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
	                   trace = "none",
	                #   xlab = "Gene ID",
	                #   ylab="Gene ID",
	                   na.color="grey",
	                   labRow=paste0(splt_rows_mirs,'/',splt_rows_genes),#rownames(tsg_sig_score_cor_by_gene[[tsg_name]]),#converted_symbols_map[rownames(sig_mut_props_hmap)],
	                   #labCol=colnames(autocors),#gene_sig,
	                   main = paste0('miRNA/mRNA part. cor. for TSG\n',sig_name),#paste0("Partial to PTEN mut"),
	                   dendrogram = "both",
	                   breaks = seq(-1,1,length=101),
	                   #symbreaks = T,
	                   Rowv = T,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.4,cexCol=0.9)
dev.copy(pdf,paste0('sig_TSG_miR_mRNA_cor/',sig_name,'_miRNA_mRNA_part_cor_for_tsg.pdf'),width=12,height=12)
  dev.off()

}


#let's repeat the above calculation but not restrict ourselves to the TSG only this time:

dir.create('sig_all_genes_miR_mRNA_cor')
#can check sig by sig what things are there,

for(sig_name in all_sigs){
	print(sig_name)
#print(pair_names[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]]),])
	pairs_to_plot <- rownames(ranked_cors_table)[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]])]
	
	splt_rows_genes <- converted_symbols_map[pair_names[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]]),2]]
	splt_rows_mirs <- pair_names[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]]),1]#splt_rows[seq(from=1,to=length(splt_rows[,1]),by=2),1]

	plotting_mat <- all_miRNA_mRNA_cor_vals_sig_subset[pairs_to_plot,]
	plotting_mat[is.na(plotting_mat)] <- 0
	p <- gplots::heatmap.2(  plotting_mat,
	                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
	                   trace = "none",
	                #   xlab = "Gene ID",
	                #   ylab="Gene ID",
	                   na.color="grey",
	                   labRow=paste0(splt_rows_mirs,'/',splt_rows_genes),#rownames(tsg_sig_score_cor_by_gene[[tsg_name]]),#converted_symbols_map[rownames(sig_mut_props_hmap)],
	                   #labCol=colnames(autocors),#gene_sig,
	                   main = paste0('miRNA/mRNA part. cor.\n',sig_name),#paste0("Partial to PTEN mut"),
	                   dendrogram = "both",
	                   breaks = seq(-1,1,length=101),
	                   #symbreaks = T,
	                   Rowv = T,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.2,cexCol=0.9)
	dev.copy(pdf,paste0('sig_all_genes_miR_mRNA_cor/',sig_name,'_miRNA_mRNA_part_cor.pdf'),width=12,height=12)
  	dev.off()
}


#-----------check the correlation of all miRNA that target a TSG with tsg expression partial to mutation-----

load('all_miRNA_targets.rda')
miRs_of_interest <- list()
for(tsg_name in cosmic_tsg[,1]){	
	good_miRs <- c()
	tsg_entrez <- names(converted_symbols_map)[which(converted_symbols_map == tsg_name)]
	for(i in 1:length(targets_list)){
		if(tsg_entrez %in% targets_list[[i]]){
			good_miRs<- c(good_miRs,i)
		}
	}
	miRs_of_interest[[tsg_name]] <- names(targets_list)[good_miRs]
}

dir.create('tsg_all_miR_cor_hmaps')
mir_cor_hmap_partial <- list()
for(tsg_name in cosmic_tsg[,1]){
	tsg_entrez <- names(converted_symbols_map)[which(converted_symbols_map == tsg_name)]

	mir_cor_hmap_partial[[tsg_name]] <- matrix(0,nrow=length(miRs_of_interest[[tsg_name]]) + 1,ncol=length(cancer_types))
	row.names(mir_cor_hmap_partial[[tsg_name]]) <- c(miRs_of_interest[[tsg_name]],'SUM')
	colnames(mir_cor_hmap_partial[[tsg_name]]) <- cancer_types
	for(i in 1:length(cancer_types)){
		if(dim(all_mRNA_datasets[[cancer_types[i]]])[2] > 0){
			if(tsg_name %in% rownames(mut_data[[cancer_types[i]]])){
				# mutated_samples <- colnames(mut_data[[cancer_types[i]]])[which(mut_data[[cancer_types[i]]][tsg_name,] > 0)]
				# non_mutated_samples <- colnames(mut_data[[cancer_types[i]]])[which(mut_data[[cancer_types[i]]][tsg_name,] == 0)]
				mRNA_vals <- all_mRNA_datasets[[cancer_types[i]]][tsg_entrez,]
				for(j in 1:length(miRs_of_interest[[tsg_name]])){
					miR_vals <- all_miRNA_datasets[[cancer_types[i]]][miRs_of_interest[[tsg_name]][j],]
					#now we do the partial correlation
					#take union of any non-na samples first
					com_non_na_vals <- intersect(which(!is.na(miR_vals)),intersect(which(!is.na(mRNA_vals)), which(!is.na(mut_data[[cancer_types[i]]][tsg_name,]))))
					if(length(com_non_na_vals) > 9){
						cor_val <- pcor.test(x=as.numeric(miR_vals[com_non_na_vals]),y=as.numeric(mRNA_vals[com_non_na_vals]),z=as.numeric(mut_data[[cancer_types[i]]][tsg_name,com_non_na_vals] > 0),method='spearman')
						if(!is.nan(cor_val$p.val)){
						if(cor_val$p.val < 0.05){
							mir_cor_hmap_partial[[tsg_name]][j,i] <- cor_val$estimate
						}else{
							mir_cor_hmap_partial[[tsg_name]][j,i] <- 0
						}
					}
					}
				}
				com_non_na_vals <- intersect(which(!is.na(mRNA_vals)), which(!is.na(mut_data[[cancer_types[i]]][tsg_name,])))
				if(length(com_non_na_vals) > 9){
					cor_val <- pcor.test(x=as.numeric(colSums(all_miRNA_datasets[[cancer_types[i]]][miRs_of_interest[[tsg_name]],com_non_na_vals],na.rm=T)),y=as.numeric(mRNA_vals[com_non_na_vals]),z=as.numeric(mut_data[[cancer_types[i]]][tsg_name,com_non_na_vals] > 0),method='spearman')
					if(!is.nan(cor_val$p.val)){

					if(cor_val$p.val < 0.05){
						mir_cor_hmap_partial[[tsg_name]][length(miRs_of_interest[[tsg_name]])+1,i] <- cor_val$estimate
					}else{
						mir_cor_hmap_partial[[tsg_name]][length(miRs_of_interest[[tsg_name]])+1,i] <- 0
					}
				}
				}
			}
		}
	}

	#now let's do rankproduct on this!
	if(!is.null(dim(mir_cor_hmap_partial[[tsg_name]]))){
	if((dim(mir_cor_hmap_partial[[tsg_name]])[1] > 1) & (dim(mir_cor_hmap_partial[[tsg_name]])[2] > 1)){

		mir_cor_hmap_partial[[tsg_name]] <- mir_cor_hmap_partial[[tsg_name]][,which(colSums(mir_cor_hmap_partial[[tsg_name]]==0) < length(mir_cor_hmap_partial[[tsg_name]][,1]))]
	}
}
	if(!is.null(dim(mir_cor_hmap_partial[[tsg_name]]))){
	if((dim(mir_cor_hmap_partial[[tsg_name]])[1] > 1) & (dim(mir_cor_hmap_partial[[tsg_name]])[2] > 1)){

		mir_cor_hmap_partial[[tsg_name]] <- mir_cor_hmap_partial[[tsg_name]][which(rowSums(mir_cor_hmap_partial[[tsg_name]] ==  0) <= (length(colnames(mir_cor_hmap_partial[[tsg_name]])) - 5)),]
	}
}
	if(!is.null(dim(mir_cor_hmap_partial[[tsg_name]]))){
	if((dim(mir_cor_hmap_partial[[tsg_name]])[1] > 1) & (dim(mir_cor_hmap_partial[[tsg_name]])[2] > 1)){

		gplots::heatmap.2( mir_cor_hmap_partial[[tsg_name]],
	                   col = gplots::colorpanel(100,"blue","white","red"), #redgreen(100),#colorpanel(100,"red","yellow","green"),
	                   trace = "none",
	                #   xlab = "Gene ID",
	                #   ylab="Gene ID",
	                   na.color="grey",
	                   labRow=rownames(mir_cor_hmap_partial),#converted_symbols_map[rownames(sig_mut_props_hmap)],
	                   #labCol=colnames(autocors),#gene_sig,
	                   main = paste0("Partial to " ,tsg_name," mut"),
	                   dendrogram = "both",
	                   breaks = seq(-1,1,length=101),
	                   #symbreaks = T,
	                   Rowv = T,Colv=T ,key.ylab=NA,  key.title=NA,margins=c(7,7),cexRow=0.5,cexCol=0.9)

		dev.copy(pdf,paste0('tsg_all_miR_cor_hmaps/',tsg_name,'.pdf'),width=12,height=12)
		dev.off()
		ranked_cors <- RP(mir_cor_hmap_partial[[tsg_name]],cl=rep(1,length(colnames(mir_cor_hmap_partial[[tsg_name]]))))
		ranked_cors_table <- topGene(ranked_cors,method='pfp',cutoff=0.05,gene.names = rownames(mir_cor_hmap_partial[[tsg_name]]))
		print(ranked_cors_table$Table1)
		write.table(ranked_cors_table$Table1,file=paste0('tsg_all_miR_cor_hmaps/',tsg_name,'_table_down.txt'),sep='\t',quote=F)
		write.table(ranked_cors_table$Table2,file=paste0('tsg_all_miR_cor_hmaps/',tsg_name,'_table_up.txt'),sep='\t',quote=F)
	}
}
}


#--------------figure out how many signature genes are significantly downregulated by a miR that is significantly pos associated with the signature itself----
load('all_miRNA_mRNA_cor_vals_sig_subset.rda')
prop_down <- list()
for(sig_name in all_sigs){
	print(sig_name)
	#first make a list of sig gene/miRNA pairs for all sig miRs up regulated with the signature
	all_combos <- c()
	# for(gene_name in sigs_list_by_name[[sig_name]][,1]){
		for(miR_name in sig_miRNA_up[[sig_name]]){
			all_combos <- c(all_combos,paste0(miR_name,'/',sigs_list_by_name[[sig_name]][,1]))
		}
	# }
	rows_to_show <- intersect(all_combos,rownames(all_miRNA_mRNA_cor_vals_sig_subset))
	if(length(rows_to_show)>=1){
	splt_rows <- melt(strsplit(rows_to_show,split="/"))
	splt_rows_genes <- converted_gene_names[splt_rows[seq(from=2,to=length(splt_rows[,1]),by=2),1]]
	splt_rows_mirs <- splt_rows[seq(from=1,to=length(splt_rows[,1]),by=2),1]
	print(rows_to_show)
	print(length(unique(splt_rows_genes)))
	print(length(unique(splt_rows_genes))/length(sigs_list_by_name[[sig_name]][,1]))
	prop_down[[sig_name]] <- length(unique(splt_rows_genes))/length(sigs_list_by_name[[sig_name]][,1])
	print(length(unique(splt_rows_mirs)))
	print(length(unique(splt_rows_mirs))/length(sig_miRNA_up[[sig_name]]))
	}else{
		prop_down[[sig_name]] <- 0
	}
	# then check the intersection with the rows of all the signficant interactions


		# repeat for all the miRs sig negatively correlated w the signature?

	# com_genes <- intersect(unique(converted_gene_names[which(pair_names[,1] %in% sig_miRNA_up[[sig_name]])]),cosmic_tsg[,1])
	# all_combos <- c()
	# for(miR_name in sig_miRNA_up[[sig_name]]){
	# 	all_combos <- c(all_combos,paste0(miR_name,'/',names(converted_gene_names)[which(converted_gene_names %in% com_genes)]))
	# }
	# rows_to_show <- intersect(all_combos,rownames(all_miRNA_mRNA_cor_vals_sig_subset))
	# splt_rows <- melt(strsplit(rows_to_show,split="/"))
	# splt_rows_genes <- converted_gene_names[splt_rows[seq(from=2,to=length(splt_rows[,1]),by=2),1]]
	# splt_rows_mirs <- splt_rows[seq(from=1,to=length(splt_rows[,1]),by=2),1]

	# plotting_mat <- all_miRNA_mRNA_cor_vals_sig_subset[rows_to_show,]

}
p <- ggplot(melt(prop_down),aes(x= reorder(L1, value,  FUN=mean) , y= value)) + geom_bar(stat="identity", position=position_dodge()) +
			labs(title = 'Prop of downregulated targets by sig',x='Signature',y='Proportion of downregulated targets') + 
			theme_minimal()+
			theme(axis.text.x = element_text(angle = 60, hjust = 1)) 
			# theme_minimal()
	print(p)


#########
#The following is to create the chord plot of the significantly regulated TSG-miRNA pairs
overall_tsg_down <- read.table('tsg_all_miR_cor_hmaps/allTSG_table_down.txt',header=T,quote="",sep='\t')
overall_tsg_down <- rownames(overall_tsg_down)
hallmarks_tsg_down <- rownames(ranked_cors_table)

pair_names <- strsplit(hallmarks_tsg_down,split = '/')
pair_names <- melt(pair_names)
pair_names <- cbind(as.character(pair_names$value[seq(from=1,to=length(pair_names$value),by=2)]),
	as.character(pair_names$value[seq(from=2,to=length(pair_names$value),by=2)]))

converted_gene_names <- converted_symbols_map[pair_names[,2]]
hallmarks_tsg_down <- paste0(pair_names[,1],'/',converted_gene_names)
sig_tsg_down <- intersect(overall_tsg_down,hallmarks_tsg_down)
print(sig_tsg_down)
	
splt_names <- melt(strsplit(sig_tsg_down,split='/'))
splt_names_rows <- splt_names[seq(from=1,to=length(splt_names[,1]),by=2),1]
splt_names_cols <- splt_names[seq(from=2,to=length(splt_names[,1]),by=2),1]

chord_mat <- matrix(0, nrow=length(unique(splt_names_rows)),ncol=length(unique(splt_names_cols)))
row.names(chord_mat) <- unique(splt_names_rows)
colnames(chord_mat) <- unique(splt_names_cols)
for(i in 1:length(splt_names_rows)){
	chord_mat[as.character(splt_names_rows[i]),as.character(splt_names_cols[i])] <- 1
}


library(GOplot)
dev.new()
GOChord(chord_mat, space = 0.02,  gene.space = 0.25, gene.size = 2,border.size = 0,nlfc = 0,process.label=10)
dev.copy(pdf,paste0('sig_TSG_miR_mRNA_chordPlot_for_paper.pdf'),width=10,height=11)
dev.off()



