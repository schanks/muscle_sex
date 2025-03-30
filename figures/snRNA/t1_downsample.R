library(data.table)

celltypes=c("Type_1","Type_2a","Type_2x","Mesenchymal_Stem_Cell","Satellite_Cell","Neuromuscular_junction","Endothelial","Smooth_Muscle","Macrophage","Neuronal")

cell_map=as.data.frame(celltypes)
colnames(cell_map)="cell"
cell_map$list=c(10,11,12,4,7,5,2,8,3,6)
cell_map$mean=rep(0, nrow(cell_map))

eds.orig<-readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/data/eds.orig.Rds")
pds.orig<-readRDS("/net/snowwhite/home/aujackso/sn_muscle_2023/data/pds.orig.Rds")


pds_t1=pds.orig[[10]]
pds_t1=pds_t1[which(pds_t1$tot_nuc>=10),]
t1.drop10=pds_t1$labelcode
eds_t1=eds.orig[[10]]
eds_t1=eds_t1[,c("gene",t1.drop10)]
count.t1=eds_t1[,-1]
nsamp_t1=ncol(count.t1)
ct.t1.non0=rowSums(count.t1!=0)
gene_t1=as.data.frame(eds_t1$gene)
colnames(gene_t1)="gene"
gene_t1$prop.non0=ct.t1.non0/nsamp_t1
gene_t1=gene_t1[which(gene_t1$prop.non0>=0.25),]
eds_t1=eds_t1[which(is.element(eds_t1$gene, gene_t1$gene)),]

for (i in 1:nrow(cell_map)){
	if (i!=1){

	index=cell_map$list[i]
	pds=pds.orig[[index]]
	pds=pds[which(pds$tot_nuc>=10),]
	samps.drop10=pds$labelcode

	#Downsample people
	eds_down=eds_t1[,c("gene",samps.drop10)]

	eds=eds.orig[[index]]
	gene=as.data.frame(eds$gene)
	colnames(gene)=c("gene")
	count.tmp=eds[,-1]
	nsamp=ncol(count.tmp)
	ct.non0=rowSums(count.tmp!=0)
	gene$prop.non0=ct.non0/nsamp	
	gene=gene[which(gene$prop.non0>=0.25),]
	eds=eds[which(is.element(eds$gene, gene$gene)),]

	countsum=sum(colSums(count.tmp))
	countprop=countsum/sum(colSums(eds_down[,-1]))

	#Actually downsample counts
	count.t1.down=round(eds_down[,-1]*countprop, digits=0)
	
	print("Cell counts:")
	print(countsum)
	print("Downsample counts:")
	print(sum(colSums(count.t1.down)))

	count.t1.down$gene=eds_down$gene

	#Write out downsampled dataset
	write.table(count.t1.down, paste("Type_1_",cell_map$cell[i],".downsample.tab",sep=""), quote=FALSE, row.names=FALSE, sep="\t")

	}}


