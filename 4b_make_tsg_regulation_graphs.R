#make TSG regualtion graphs
library(ggplot2)

tsg_names <- c('ACVR2A','ARHGEF12','CDK12','DNMT3A','FAT4','PTEN','SFRP4','TGFBR2')

#load in the positive and negative regulation tables
pos_reg_tables <- list()
neg_reg_tables <- list()

for(tsg in tsg_names){
	pos_reg_tables[[tsg]] <- read.table(paste0(tsg,'-analysis/overall_rank_prod_pos_pred.txt'),quote="",header=T)
	neg_reg_tables[[tsg]] <- read.table(paste0(tsg,'-analysis/overall_rank_prod_neg_pred.txt'),quote="",header=T)
}

#make the matrix for plotting

pos_overall_toplot <- c()
neg_overall_toplot <- c()

for(tsg in tsg_names){
	pos_overall_toplot <- rbind(pos_overall_toplot, cbind(rep(tsg,times=dim(pos_reg_tables[[tsg]])[1]), rownames(pos_reg_tables[[tsg]]),pos_reg_tables[[tsg]][,'P.value']))
	neg_overall_toplot <- rbind(neg_overall_toplot, cbind(rep(tsg,times=dim(neg_reg_tables[[tsg]])[1]),rownames(neg_reg_tables[[tsg]]),neg_reg_tables[[tsg]][,'P.value']))
}

pos_overall_toplot <- rbind(pos_overall_toplot,c('ACVR2A','Methyl.',1))
pos_overall_toplot <- rbind(pos_overall_toplot,c('SFRP4','Methyl.',1))
pos_overall_toplot <- rbind(pos_overall_toplot,c('TGFBR2','Methyl.',1))

neg_overall_toplot <- rbind(neg_overall_toplot,c('SFRP4','Methyl.',1))



pos_overall_toplot  <- as.data.frame(pos_overall_toplot,stringsAsFactors=FALSE)
neg_overall_toplot  <- as.data.frame(neg_overall_toplot,stringsAsFactors=FALSE)



for(i in 1:dim(pos_overall_toplot)[1]){
	if(grepl('cnv',pos_overall_toplot[i,2])){
		pos_overall_toplot[i,2] <- 'CNV'
	}

	if(grepl('hsa',pos_overall_toplot[i,2])){
		pos_overall_toplot[i,2] <- 'miRNA'
	}

	if(grepl('cg',pos_overall_toplot[i,2])){
		pos_overall_toplot[i,2] <- 'Methyl.'
	}
}

for(i in 1:dim(neg_overall_toplot)[1]){
		if(grepl('cnv',neg_overall_toplot[i,2])){
		neg_overall_toplot[i,2] <- 'CNV'
	}

	if(grepl('hsa',neg_overall_toplot[i,2])){
		neg_overall_toplot[i,2] <- 'miRNA'
	}

	if(grepl('cg',neg_overall_toplot[i,2]) || grepl('ch.11',neg_overall_toplot[i,2])){
		neg_overall_toplot[i,2] <- 'Methyl.'
	}
}

pos_overall_toplot[,3] <- -log10(as.numeric(pos_overall_toplot[,3]))
neg_overall_toplot[,3] <- -log10(as.numeric(neg_overall_toplot[,3]))
colnames(pos_overall_toplot) <- c('gene','pred','value')
colnames(neg_overall_toplot) <- c('gene','pred','value')



pos_overall_toplot$gene <- factor(pos_overall_toplot$gene, levels=rev(sort(unique(pos_overall_toplot$gene))))

neg_overall_toplot$gene <- factor(neg_overall_toplot$gene, levels=rev(sort(unique(neg_overall_toplot$gene))))


group.colors <- c(miRNA = "#67a9cf", Methyl. = "#ef8a62", CNV ="#b2182b", Nonsense_Mutation = "#af8dc3", Frame_Shift_Del = "black")
median.colors <- c(miRNA = "black", Methyl. = "black")

#print(pos_overall_toplot)
dev.new() #geom_boxplot(outlier.shape=NA,lwd=1) #reorder(gene, value, mean)  width=0.3,height=0.3, data=pos_overall_toplot[pos_overall_toplot$pred %in% c('Methyl.', 'miRNA'),] ##+  geom_boxplot(data = subset(pos_overall_toplot,(pred %in% c('Methyl.', 'miRNA'))),aes(col=pred),outlier.shape=NA,width=0.1) #scale_colour_manual(values=median.colors) +
p <- ggplot(pos_overall_toplot, aes(x=gene,y=value,col=pred))  +geom_jitter(alpha=1,size=0.6,position=position_jitterdodge(jitter.width = 0.15))  +theme_minimal() + #+ stat_summary(data = subset(pos_overall_toplot,(pred %in% c('Methyl.', 'miRNA'))),aes(fill=factor(pred)),fun.y = "median",  size = 1, geom = "point",position=position_dodge(width=0.5),shape=2,color='black')
  scale_colour_manual(values=group.colors)  +    coord_flip(ylim = c(2, 30))
# + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x='Cancer subtype',y='log2 expression',title='DICER1 expression across cancer subtypes')
print(p)
dev.copy(pdf,'regs_POS_significance_values.pdf',height=3,width=5)
dev.off()



dev.new() #geom_boxplot(outlier.shape=NA,lwd=1) #reorder(gene, value, mean) width=0.3,height=0.3 data=neg_overall_toplot[neg_overall_toplot$pred %in% c('Methyl.', 'miRNA'),], ## +  geom_boxplot(data = subset(neg_overall_toplot,(pred %in% c('Methyl.', 'miRNA'))),aes(col=pred),outlier.shape=NA,width=0.1) 
p <- ggplot(neg_overall_toplot, aes(x=gene,y=value,col=pred))+ geom_jitter(alpha=1,size=0.6,position=position_jitterdodge(jitter.width = 0.15)) + theme_minimal()+ #+ stat_summary(data = subset(neg_overall_toplot,(pred %in% c('Methyl.', 'miRNA'))),aes(fill=factor(pred)),fun.y = "median", size = 1, geom = "point",position=position_dodge(width=0.5),shape=2,color='black')
  scale_colour_manual(values=group.colors)  + coord_flip(ylim = c(2, 18))
# + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + labs(x='Cancer subtype',y='log2 expression',title='DICER1 expression across cancer subtypes')
print(p)
dev.copy(pdf,'regs_NEG_significance_values.pdf',height=3,width=5)
dev.off()

graphics.off()
