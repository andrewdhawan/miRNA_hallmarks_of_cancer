#make network graph
#this is to make the network graph for the miRs and the signatures to better represent the data
library(reshape2)
library(GOplot)
library(qgraph)
library(igraph)

dev.new()

load('rank_prod_tables_out_pre_filtered.rda')
sig_miRNA_up <- list()
sig_miRNA_down <- list()
for (sig_name in all_sigs){
	sig_miRNA_up[[sig_name]] <- rownames(rank_prod_tables[[sig_name]]$Table2)
	sig_miRNA_down[[sig_name]] <- rownames(rank_prod_tables[[sig_name]]$Table1)
}

all_up_miRNA <- unique(melt(sig_miRNA_up)$value)
all_down_miRNA <- unique(melt(sig_miRNA_down)$value)

all_sigs_modified <- c()
for(sig_name in all_sigs){
	
	if(grepl('Hallmark: ',sig_name)){
		sig_name <- gsub('Hallmark: ','',sig_name)
		sig_name <- paste0(sig_name,', MSigDB')
	}
	all_sigs_modified <- c(all_sigs_modified,sig_name)
}


all_names <- c(as.character(all_up_miRNA),all_sigs_modified)

all_names <-gsub('miR-','', gsub('hsa-','',all_names))
adjacency_mat <- matrix(0,nrow=length(all_names),ncol=length(all_names))
row.names(adjacency_mat) <- all_names
colnames(adjacency_mat) <- all_names

for(sig_name in all_sigs_modified){
	for(mir_name in sig_miRNA_up[[sig_name]]){
		mir_name <- gsub('miR-','',gsub('hsa-','',mir_name))
		adjacency_mat[mir_name,sig_name] <- 1
		# adjacency_mat[sig_name,mir_name] <- 1

	}
}

g1 <- graph_from_adjacency_matrix(adjacency_mat)	
minC <- rep(-Inf, vcount(g1))
maxC <- rep(Inf, vcount(g1))
minC[1] <- maxC[1] <- 0

l <-layout_nicely(g1)#,weight.node.dist=0)#layout_with_graphopt(g1,spring.length = 50, charge = 0.1) #layout_nicely(g1) #layout_with_lgl(g1)#,spring.length	= 4,charge=0.00001)#,niter=500,area=vcount(g1)^2.3,repulserad=vcount(g1)^2.8)
	cols <- c(rep('lightgrey',times=length(all_up_miRNA)),rainbow(length(all_sigs)))
plot.igraph(g1,vertex.size=4.5,vertex.label.cex=0.35, edge.width=0.5,vertex.label.dist=0,vertex.label.color='black',vertex.label=c(gsub('miR-','',gsub('hsa-','',all_up_miRNA)),rep('',times=length(all_sigs))), layout=l,
	vertex.color=cols,vertex.frame.color='darkgrey',edge.curved=T,edge.arrow.size=0,edge.color='black',edge.label.family = 'Helvetica',edge.label.font =2)
plot.new()
legend('topright',all_sigs_modified,fill=cols[(length(cols) - 23):length(cols)],cex=0.3)
dev.copy(pdf,paste0('network_graph_hallmarks-legend','.pdf'),width=14,height=14)
dev.off()

#-----This is a repeat of the above code for the DOWNregulated miRNA

#make network graph
#this is to make the network graph for the miRs and the signatures to better represent the data
library(reshape2)
library(GOplot)
library(qgraph)
library(igraph)
dev.new()
load('rank_prod_tables_out_pre_filtered.rda')
sig_miRNA_up <- list()
sig_miRNA_down <- list()
for (sig_name in all_sigs){
	sig_miRNA_up[[sig_name]] <- rownames(rank_prod_tables[[sig_name]]$Table2)
	sig_miRNA_down[[sig_name]] <- rownames(rank_prod_tables[[sig_name]]$Table1)
}
all_up_miRNA <- unique(melt(sig_miRNA_up)$value)
all_down_miRNA <- unique(melt(sig_miRNA_down)$value)




all_names <- c(as.character(all_down_miRNA),all_sigs)

all_names <-gsub('miR-','', gsub('hsa-','',all_names))
adjacency_mat <- matrix(0,nrow=length(all_names),ncol=length(all_names))
row.names(adjacency_mat) <- all_names
colnames(adjacency_mat) <- all_names

for(sig_name in all_sigs){
	for(mir_name in sig_miRNA_down[[sig_name]]){
		mir_name <- gsub('miR-','',gsub('hsa-','',mir_name))
		adjacency_mat[mir_name,sig_name] <- 1
		# adjacency_mat[sig_name,mir_name] <- 1

	}
}


g1 <- graph_from_adjacency_matrix(adjacency_mat)	
minC <- rep(-Inf, vcount(g1))
maxC <- rep(Inf, vcount(g1))
minC[1] <- maxC[1] <- 0

l <-layout_nicely(g1)#,weight.node.dist=0)#layout_with_graphopt(g1,spring.length = 50, charge = 0.1) #layout_nicely(g1) #layout_with_lgl(g1)#,spring.length	= 4,charge=0.00001)#,niter=500,area=vcount(g1)^2.3,repulserad=vcount(g1)^2.8)
	cols <- c(rep('lightgrey',times=length(all_down_miRNA)),rainbow(length(all_sigs)))
dev.new()
plot.igraph(g1,vertex.size=4.5,vertex.label.cex=0.35, edge.width=0.5,vertex.label.dist=0,vertex.label.color='black',vertex.label=c(gsub('miR-','',gsub('hsa-','',all_down_miRNA)),rep('',times=length(all_sigs))), layout=l,
	vertex.color=cols,vertex.frame.color='darkgrey',edge.curved=T,edge.arrow.size=0,edge.color='black',edge.label.family = 'Helvetica',edge.label.font =2)
#plot.new()
#legend('topright',all_sigs,fill=cols[(length(cols) - 23):length(cols)],cex=0.3)
dev.copy(pdf,paste0('network_graph_hallmarks_DOWN_miRNA','.pdf'),width=14,height=14)
dev.off()