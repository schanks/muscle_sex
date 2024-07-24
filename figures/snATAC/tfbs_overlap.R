library(data.table)
library(optparse)

optionList <- list(
  make_option(c("-c", "--cell"), type="character", help="celltype"),
  make_option(c("-t", "--tfbs"), type="character", help="TFBS bedfile path")
)

parser <- OptionParser(
  usage="%prog -c cell -t tfbs",
  option_list=optionList
)

args <- parse_args(parser, positional_arguments = 0)
opt <- args$options
print(opt)

cell=opt$cell

dat=fread(paste("/net/snowwhite/home/aujackso/sn_muscle_2023/output/DESeq.ATAC/final_drop10nuc/results/",cell,".SEX.M.results.tab",sep=""))

dat$chrom=unlist(strsplit(dat$peak, ":"))[seq(1, nrow(dat)*3, by=3)]
dat$start=unlist(strsplit(dat$peak, ":"))[seq(2, nrow(dat)*3, by=3)]
dat$end=unlist(strsplit(dat$peak, ":"))[seq(3, nrow(dat)*3, by=3)]

tfbs=fread(opt$tfbs, header=FALSE)
colnames(tfbs)=c("chrom","start","stop","tf","score","strand")
tf=tfbs$tf[1]
dat$overlap=rep(0, nrow(dat))
dat$start=as.numeric(dat$start)
dat$end=as.numeric(dat$end)

for (i in 1:nrow(dat)){
	t=tfbs[which(tfbs$chrom==dat$chrom[i] & ((tfbs$start>dat$start[i] & tfbs$start<dat$end[i]) | (tfbs$stop>dat$start[i] & tfbs$stop<dat$end[i]))),]
	if (nrow(t)>0){dat$overlap[i]=1}
	print(i)
}

colnames(dat)[15]=tf
write.table(dat, tf, sep="\t", quote=FALSE, row.names=FALSE)



