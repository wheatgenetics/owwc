# load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load('data.table','ape','phyclust','rrBLUP','ggplot2',
               'ggtree','tidytree')

# read the pheno file
pop <-read.table(file="data/selected_lines_interval.txt", header = T, as.is = T, check.names = F)

# samples info
wheat.lines <- pop$line[pop$specie == 'wheat']
tauschii.lines <- pop$line[pop$specie == 'tauschii']

# read hap file numeric
hap01 <- fread("data/hmp_selected_lines.txt", header = T, check.names = F, data.table = F)
geno <- t(as.matrix(hap01[, 16:ncol(hap01)])) # only genotypes

######################
## Ae. tauschii only
######################
geno <- geno[rownames(geno) %in% tauschii.lines, ]
geno[1:15,1:5]

###################################
## Hierarchical Cluster Analysis 
###################################
# compute genetic distance
ifelse(test = file.exists('data/distMat.RData'),
       yes = load('data/distMat.RData'), no = distMat <- dist(geno))
# saving distance matrix because it's a computationally intensive step
if (!file.exists('data/distMat.RData')) save(distMat, file = "data/distMat.RData")

# cluster accessions and convert to phylo object for ape
hc2 <- as.phylo(hclust(distMat))

# cluster coloring
edgecols <- cbind('line'=NA, 1:nrow(hc2$edge), color='black') # create data frame
edgecols[,1] <- hc2$tip.label[hc2$edge[,2]] # get labels
edgecols <- as.matrix(merge(edgecols, pop, by = 'line', all.x = T)) # get samples info
edgecols <- edgecols[order(as.numeric(edgecols[,2])), ] # get samples in original order
edgecols[,5] <- as.numeric(edgecols[,5])

# coloring R/S accessions with different colors
edgecols[,3][edgecols[,5] == 0] = "grey70"
edgecols[,3][edgecols[,5] == 1] = "purple"

# tips coloring
tipcols <- as.matrix(merge(hc2$tip.label, edgecols, by = 1, sort = F))
tipcols[,7] <- as.numeric(tipcols[,7])

# coloring diff lineages with different colors
tipcols[,3][tipcols[,7] == 1] = "blue"
tipcols[,3][tipcols[,7] == 2] = "red"

# plotting tree
pdf("Fig.wcm_chr6D_cluster.pdf", width = 10, height = 10)
plotnj(unrooted.tree=hc2, type='p',
       show.tip.label=T, lab4ut="axial", label.offset=0.2, cex=0.4,
       edge.color=edgecols[,3], edge.width=1, tip.color=tipcols[,3], rotate.tree=-20)
legend(-0.5, 80, lty=1, lwd=2, cex=0.75,
       legend = c('Lineage 1 (L1)', 'Lineage 2 (L2)', 'WCM Resistant', 'WCM Susceptible'),
       text.col = 'black', col = c("blue", "red", "purple", "grey70"))
dev.off()

