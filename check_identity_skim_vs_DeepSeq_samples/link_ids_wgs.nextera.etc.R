library(data.table)
library(openxlsx)
library(zoo)

fread('~/owwc/archive/200520.bw_id.ksu_id.lineage.wcm.txt') -> z
### prepare nextera samples
x=fread('find /bulk/jpoland/raw_data/psomagen/1909UNHS-0017/NEX0009/demultiplex/ -type f -name *R1.fq ', head=F)
x[,id:=sub('.+/(.+_idx_\\d+)-R\\d.fq', '\\1', V1)]
fread('/homes/lianggao/owwc/data/191121_KEY_NEX0009-BW-OWWC-samples.txt') -> d
x[,table(substr(id,1,3))]
# Bla BW_  NA T0W TOW 
# 41 788  79  12 136 
x[,id2:=id][,id2:=sub('T0W','TOW', id2)][,id2:=sub('TOWWC_01','TOWWC',id2)][,id2:=sub('TOWWC_','TOWWC',id2)][,id3:=id2][,id3:=sub('-[ABC]','',id3)][,id3:=sub('_idx_\\d+','',id3)]
cnt.sum= fread('/bulk/jpoland/raw_data/psomagen/1909UNHS-0017/NEX0009/demultiplex/summary_report-NEX0009_R1.fastq.gz-200520_key_nextera_v1.txt.txt', col.names=c('id','cnt'))
cnt.sum[x,on='id'] -> x
x[,summary(cnt/1e6)]
x[cnt>100e3]
x[cnt<50e3]
rosetta = data.table(read.xlsx('~/owwc/data/Rosetta_tauschii.xlsx'))
setnames(rosetta, c('source.id','bw15','bw16','gru','bw16.v2','bw17','bw18'))
rosetta[,wgs.id:=na.locf(bw15)][,source.id:=na.locf(source.id)][,wgs.id2:=sub('\\*','', wgs.id)][] -> r
melt(r, id.vars=c('wgs.id2','source.id'), measure.vars=c('bw16','gru','bw16.v2','bw17','bw18')) -> m
m[!(is.na(value)), value]
cnt.cut=25e3
setdiff(m[!(is.na(value)), value], x[cnt>cnt.cut,id3])
intersect(m[!(is.na(value)), value], x[cnt>cnt.cut,id3])
setdiff(x$id3, m[!(is.na(value)), value])
x[,id2:=sub('-','_',id2)]
intersect(x[cnt>cnt.cut, substr(unique(id2),1,8)],m[!(is.na(value)), value] )
m[x[cnt>cnt.cut,], on=c(value='id3')] -> mx

merge(mx,z, by.x='wgs.id2', by.y='id', all=T) -> mxz
setnames(mxz, c('wgs.id', 'source.id','year.bw', 'rosetta.id', 'nextera.id', 'cnt', 'loc', 'nextera.id.correct', 'ksu_id', 'lineage','lineage.pc', 'wcm.cat', 'wcm.qdr','batch' ))
# write.table(mxz, file='~/owwc/data/200701.rosetta.melt.nextera.id.wgs.id.txt', sep='\t', quote=F, row.names=F)

