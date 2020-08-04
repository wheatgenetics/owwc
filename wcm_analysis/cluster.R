# load required libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load('data.table','ape','phyclust','tidytree')

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

# coloring R/S accessions with different colors
edgecols[,3][edgecols[,7] == 'S'] = "grey70"
edgecols[,3][edgecols[,7] == 'R'] = "purple"

# tips coloring
tipcols <- as.matrix(merge(hc2$tip.label, edgecols, by = 1, sort = F))

# coloring diff lineages with different colors
tipcols[,3][tipcols[,5] == 'L1'] = "blue"
tipcols[,3][tipcols[,5] == 'L2'] = "red"

# plotting tree
pdf("Fig.S1_wcm.pdf", width = 10, height = 10)
plot.phylo(hc2, type='p', show.tip.label=T, label.offset=0.2, cex=0.4,
       edge.color=edgecols[,3], edge.width = 2, tip.color=tipcols[,3])
legend(0, 85, lty = 1, lwd = 3, cex = 1, box.lwd = 1,
       legend = c('Lineage 1 (L1)','Lineage 2 (L2)','WCM Resistant','WCM Susceptible'),
       text.col = c("blue","red","black","black"), col = c("black","black","purple","grey70"))
axisPhylo(1, cex.axis = 0.8)
dev.off()


########################
## Ae. tauschii + wheat
########################
geno <- t(as.matrix(hap01[, 16:ncol(hap01)])) # only genotypes

###################################
## Hierarchical Cluster Analysis 
###################################
# compute genetic distance
ifelse(test = file.exists('data/distMatall.RData'),
       yes = load('data/distMatall.RData'), no = distMatall <- dist(geno))
# saving distance matrix because it's a computationally intensive step
if (!file.exists('data/distMatall.RData')) save(distMatall, file = "data/distMatall.RData")

# cluster accessions and convert to phylo object for ape
hc2 <- as.phylo(hclust(distMatall))

# cluster coloring
edgecols <- cbind('line'=NA, 1:nrow(hc2$edge), color='black') # create data frame
edgecols[,1] <- hc2$tip.label[hc2$edge[,2]] # get labels
edgecols <- as.matrix(merge(edgecols, pop, by = 'line', all.x = T)) # get samples info
edgecols <- edgecols[order(as.numeric(edgecols[,2])), ] # get samples in original order

# coloring R/S accessions with different colors
edgecols[,3][edgecols[,7] == 'S'] = "grey70"
edgecols[,3][edgecols[,7] == 'R'] = "purple"

# tips coloring
tipcols <- as.matrix(merge(hc2$tip.label, edgecols, by = 1, sort = F))

# coloring diff lineages with different colors
tipcols[,3][tipcols[,5] == 'L1'] = "blue"
tipcols[,3][tipcols[,5] == 'L2'] = "red"
tipcols[,3][tipcols[,5] == 'W'] = "green2"

# plotting tree
pdf("Fig.S3_wcm.pdf", width = 10, height = 10)
plot.phylo(hc2, type = 'p', show.tip.label = T, label.offset = 0.2, cex = 0.4,
       edge.color = edgecols[,3], edge.width = 2, tip.color = tipcols[,3]) 
legend(0, 120, lty = 1, lwd = 3, cex = 1, box.lwd = 1,
       legend = c('Lineage 1 (L1)','Lineage 2 (L2)','Wheat (W)','WCM Resistant','WCM Susceptible'),
       text.col = c("blue","red","green2","black","black"), col = c("black","black","black","purple","grey70"))
axisPhylo(1, cex.axis = 0.8)
dev.off()


