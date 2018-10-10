#code to run sigQC on similar signatures over all datasets

library(sigQC)
library(reshape2)

cancer_types_list <- list();
cancer_types_list[[1]] <- c('BRCA','UCEC','HNSC')
cancer_types_list[[2]] <- c('KIRC','LUAD','THCA')
cancer_types_list[[3]] <- c('PRAD','LUSC','OV')
cancer_types_list[[4]] <- c('STAD','BLCA','COAD')
cancer_types_list[[5]] <- c('LIHC','CESC','KIRP')

all_cancer_types <- melt(cancer_types_list)$value

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

#load in all of the signatures

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

#load in all datasets
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


for(sig_category in categories_of_sigs){
	#sig_category <- 'angiogenesis'
	dir.create(sig_category)
	setwd(sig_category)
	for (cancer_types in cancer_types_list){
		print(sig_category)
		print(cancer_types)

		make_all_plots(sigs_list_by_cat[[sig_category]], all_mRNA_datasets[cancer_types], sig_names_list[[sig_category]],cancer_types, covariates=NULL, out_dir = toString(cancer_types),showResults=F)
	}
	setwd("../")
}