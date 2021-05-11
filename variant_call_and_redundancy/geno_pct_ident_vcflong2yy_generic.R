#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}

library(data.table)
fread(paste0('zcat ', args[1])) -> g
setnames(g, c("chr", "pos", "ref", "alt", "qual", "sample.id", "gt", "dp", "dv"))
g[gt == '0/0', gt := paste(0)]
g[gt == '0/1', gt := paste(1)]
g[gt == '1/1', gt := paste(2)]
g[,length(unique(pos))] -> tot.num.snps
fread('~/owwc/data/200701.rosetta.melt.nextera.id.wgs.id.txt') -> mxz

yy=list()
mxz[(!is.na(nextera.id.correct) | !is.na(ksu_id)) & !is.na(wgs.id), unique(wgs.id)] -> WGS.IDs
for (SEL.BW.Line in WGS.IDs){
  mxz[wgs.id==SEL.BW.Line,unique(nextera.id.correct)]-> sel1
  mxz[wgs.id==SEL.BW.Line,c(unique(ksu_id))]-> sel2
  setdiff(c(sel1,sel2), NA) -> sel
  g[sample.id %in% c(SEL.BW.Line, sel)]->x
  if(length(unique(x$sample.id))<2) {next}
  x[sample.id  %in% SEL.BW.Line, .(chr, pos, sample1=sample.id, gt1=gt)]->a
  x[sample.id  %in% sel, .(chr, pos, sample2=sample.id, gt2=gt)]->b
  if(nrow(a)<1){next}
  a[b, on=c("chr", "pos"), allow.cartesian=T][, d := abs(as.integer(gt1)-as.integer(gt2))][]->ab
  dcast(ab[, .N, key=.(sample1, sample2, d)], sample1 + sample2 ~ d, value.var="N", fill=0) -> y
  if(!("0" %in% colnames(y))){y[,N.ident:=0]}
  if(!("1" %in% colnames(y))){y[,N.het:=0]}
  if(!("2" %in% colnames(y))){y[,N.diff:=0]}
  setnames(y, c("0","1", "2"), c("N.ident", "N.het", "N.diff"), skip_absent=T)[, n := N.ident + N.diff]
  y[, pdiff := N.diff/n]
  yy[[SEL.BW.Line]]=y[,.(sample1,sample2,N.ident, N.het, N.diff, n, pdiff)]
}
rbindlist(yy) -> yy
chr=sub('.+/.+_(CM.+\\.\\d):(\\d+)-(\\d+).+.tsv.gz','\\1',args[1])
start=sub('.+/.+_(CM.+\\.\\d):(\\d+)-(\\d+).+.tsv.gz','\\2',args[1])
end=sub('.+/.+_(CM.+\\.\\d):(\\d+)-(\\d+).+.tsv.gz','\\3',args[1])
yy[!is.na(sample1), .(sample1, sample2, N.ident, N.het, N.diff, n, pdiff, chr=chr, start=start, end=end, tot.num.snps=tot.num.snps)] -> yy
write.table(yy, file=args[2], sep='\t', quote=F, row.names=F)

