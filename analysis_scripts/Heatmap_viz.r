library(plotly)
library(cluster)
library(dplyr)
options(scipen = 999)
heatmap_file = read.csv('herit_Jaccard_total_array.out', header = FALSE, sep = "\t") 
names_vec = as.vector(read.csv('herit_annotated_names_for_map.csv', header = FALSE, sep = "\t"))
names_vec = as.vector(sapply(names_vec, as.character, simplify = TRUE))
colnames(heatmap_file) = names_vec
rownames(heatmap_file) = names_vec
heatmap_file = 1-heatmap_file
clusters <- hclust(as.dist(as.matrix(heatmap_file)), method = 'average')
plot(clusters, labels=FALSE)
clusterCut <- cutree(clusters, 307)

s=NULL
for(i in 10:500){
  k=cutree(clusters,i)
  s[i]=summary(silhouette(k,as.dist(heatmap_file)))$si.summary[3]
}
which.max(s)


names_frame <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
sorted_heatmap <- heatmap_file[as.vector(names_frame$keyName), as.vector(names_frame$keyName)]

write.table(names_frame, file = "herit_for_clustering_jaccard.tsv", quote=FALSE, sep='\t')
write.table(heatmap_file, file = "Herit_jaccard.tsv", quote=FALSE, sep='\t')

colfunc <- colorRampPalette(c("ForestGreen",'lightGreen' , 'white'))
heatmap_matrix = matrix(sorted_heatmap)
p <- plot_ly(z = heatmap_matrix, x = colnames(sorted_heatmap), y = rownames(sorted_heatmap),colors = colfunc(50), type = "heatmap",  reversescale =F)
htmlwidgets::saveWidget(as_widget(p), file = "herit_jaccar.html")

dend <- as.dendrogram(clusters)
dend  = reorder(dend, wts = order(match(names_frame$keyName, rownames(heatmap_file))))
ord_dend <- rev(reorder(dend, agglo.FUN=sum, 543:1))
labels(ord_dend)
plot(ord_dend,labels=FALSE)

library(ape)

plot(as.phylo(as.hclust(as.dendrogram(ord_dend))), type = "phylogram", show.tip.label = F,
     edge.color = "black", edge.width = 0.071, edge.lty = 1,
     tip.color = "black")

library(lattice)
drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    # panel.text(x, y, e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet,
                   at=seq(0, max(df), length.out=100),
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1)),
                   scales=list(y=element_blank(),x=element_blank(),cex=0,tck = c(0,0)), xlab=list(label=''),
                   ylab=list(label=''), panel=myPanel_a))
}


mypal = colorRampPalette(c('#FF4343','#FAFAFA'))
jet = mypal(100)


drawCoolHM(as.matrix(check))
check <- heatmap_file[as.vector(labels(ord_dend)),as.vector(labels(ord_dend))]

occur =vector()
for (cluster in unique(names_frame$value)){
  occur <- c(occur,sum(names_frame$value==cluster))
}



cluster_names_plot <- data.frame(cluster = unique(names_frame$value),occur =occur,fac=1)

names_frame35 <- data.frame(keyName=names(clusterCut), value=clusterCut, row.names=NULL) %>% arrange(value)
sorted_heatmap35 <- jac35[c("I80",  "20002_1094", "6152_5", "I26", "20002_1093","6152_7"), c("I80",  "20002_1094", "6152_5", "I26", "20002_1093","6152_7")]
jac35 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==35]),as.vector(names_frame$keyName[names_frame$value==35])]
clusters35 <- hclust(as.dist(as.matrix(jac35)), method = 'average')
clusterCut <- cutree(clusters35, 4)
names(clusters35$order)
plot(clusters35, hang = -1, axes = FALSE)
drawCoolHM(as.matrix(jac36))

jac21 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==21]),as.vector(names_frame$keyName[names_frame$value==21])]
jac33 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==33]),as.vector(names_frame$keyName[names_frame$value==33])]
jac73 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==73]),as.vector(names_frame$keyName[names_frame$value==73])]
jac34 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==34]),as.vector(names_frame$keyName[names_frame$value==34])]
jac24 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==24]),as.vector(names_frame$keyName[names_frame$value==24])]
rownames(jac24)[7] <- 'Pack years of smoking'
colnames(jac24)[7] <- 'Pack years of smoking'
rownames(jac24)[8] <- 'Pack years adult smoking as proportion of life span'
colnames(jac24)[8] <- 'Pack years adult smoking as proportion of life span'
jac36 <- sorted_heatmap[as.vector(names_frame$keyName[names_frame$value==36]),as.vector(names_frame$keyName[names_frame$value==36])]
rownames(jac36)[9] <- 'None of the above'
colnames(jac36)[9] <- 'None of the above'
rownames(jac36)[10] <- 'Asthma diagnosed by doctor'
colnames(jac36)[10] <- 'Asthma diagnosed by doctor'
rownames(jac36)[11] <- 'Hayfever, allergic rhinitis or eczema'
colnames(jac36)[11] <- 'Hayfever, allergic rhinitis or eczema'



mypal = colorRampPalette(c('#FF4343','#FAFAFA'))
jet = mypal(100)

ggplot(cluster_names_plot, aes(x="", y=occur, fill=as.factor(cluster)))+ coord_flip()+
  geom_bar(width = 1, stat = "identity") + scale_fill_manual(values=c(rep(c("black", "white"),153),"black"))+
  theme(legend.position="none")+scale_y_reverse()+
  theme(axis.line=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks=element_blank(),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.border=element_blank(),
        panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),
        plot.background=element_blank())
