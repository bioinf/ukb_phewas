library(dplyr)
library(ggplot2)
library(gridExtra)

cluster_file = read.csv('herit_clust_diags_jaccard.tsv', header = FALSE, sep = "\t")
str(cluster_file)
cluster_file %>% arrange(V2) %>% write.table(file = "herit_clust_diags_sort_jaccard.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)


cluster_diagram = read.csv('herit_Annotated_clusters_diagram_jaccard.tsv', header = FALSE, sep = "\t")%>% filter(V4>=3)
plot_names <-unique(as.numeric(unlist(as.vector(cluster_diagram %>% filter(V4>=3) %>% select(V1)))))
plot_nums <-unique(as.numeric(unlist(as.vector(cluster_diagram %>% filter(V4>=3) %>% select(V4)))))

cluster_names_plot <- data.frame(name = colnames(dat),occur=1)
ggplot(cluster_names_plot, aes(x="", y=occur, fill=as.factor(name)))+ coord_flip()+
  geom_bar(width = 1, stat = "identity") + scale_fill_manual(values=alpha(c("white","black","white","black","white","black"),0.995))+
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



to_plot <-cluster_diagram %>% arrange(V1) %>% filter(V1==plot_names[8]) 
pie <- ggplot(to_plot, aes(x = "", y=V3, fill = factor(V2))) + 
  geom_bar(width = 1, stat = "identity") +
  theme(axis.line = element_blank(), 
        plot.title = element_text(hjust=0.5)) + 
  labs(fill="class", 
       x=NULL, 
       y=NULL, 
       title="Pie Chart of class", 
       caption="Source: mpg")+ coord_polar(theta = "y", start=0)


library(ggplot2) 
library(treemapify)

ggplot(cluster_diagram, aes(area = V3, fill = V1, label = V2,
                subgroup = V1)) +
  geom_treemap() +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour =
                               "black", fontface = "italic", min.size = 0) +
  geom_treemap_text(colour = "white", place = "topleft", reflow = T)

library(treemap)
pdf("tree3.pdf",20,20)
treemap(cluster_diagram,
        index=c('V1', 'V2'),
        vColor='V3',
        type="index",
        vSize='V3',
        palette = matlab.like(14))

dev.off()
library(colorRamps)
colorRampPalette(matlab.like(20))
for(j in 1:length(plot_names)){
  to_plot <-cluster_diagram %>% arrange(V1) %>% filter(V1==plot_names[j]) 
  assign(paste0("plot", plot_names[j]), ggplot(to_plot, aes(x = reorder(V2, -V3), y = V3, fill =V4)) + 
           theme_bw() + 
           guides(fill=FALSE)+
           geom_bar(stat = "identity", fill ='#6597A7', alpha = 0.6)+labs(x = "Phenotype", y = "Proportion") +
           facet_wrap(~V1, scales = "free_x")+
           theme( axis.text.x = element_text(angle = 90, hjust = 1),
                  panel.background = element_blank(), 
                  axis.line = element_line(colour = "black"),
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank()))
  }


#testing stuff  
signific <- read.csv('herit_valuable_SNP_jaccard.tsv', header = FALSE, sep = "\t")
cluster_file[which(cluster_file$V1 %in% signific$V1),] %>% select(V1,V2)%>% write.table(file = "herit_for_graph_jaccard.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
test = read.csv('herit_for_graph_jaccard.tsv', header = FALSE, sep = "\t")


#graph
graph_file = read.csv('herit_clump_array_for_graph_jaccard.out', header = FALSE, sep = "\t") 
names_vec = as.vector(read.csv('herit_graph_names_jaccard.tsv', header = FALSE, sep = "\t"))
names_vec = as.vector(sapply(names_vec, as.character, simplify = TRUE))
colnames(graph_file) = names_vec
rownames(graph_file) = names_vec
diag(graph_file) <- 0

max(m) #44 max
library(igraph)
m=as.matrix(graph_file)
net=graph.adjacency(m,mode="undirected",weighted=TRUE,diag=FALSE)

E(net)$color<- adjustcolor("#5B0493",alpha.f = 0.21)
for ( i in 1:length (E(net)[[]])){
  num <- paste0('0.',gsub('\\.','',toString(ifelse(E(net)[[i]]$weight/72*100*55/100>2, E(net)[[i]]$weight/72*100*55/100/10, round(E(net)[[i]]$weight/72*100*55/100,3)))))
  E(net)[[i]]$color <- adjustcolor("darkgreen",alpha.f = as.numeric(num))
}
#E(net)[[6]] 72/72*100*21/100

#ifelse(E(net)[[6]]$weight/72*100*21/100>2, E(net)[[6]]$weight/72*100*55/100/10, round(E(net)[[6]]$weight/72*100*55/100,3))
full_g <- plot.igraph(net,vertex.label=V(net)$name,
            vertex.color="white",
            layout=layout.circle,
            edge.width=E(net)$weight/15,
            vertex.size=5.6,vertex.label.cex = 0.34)
cliques <- read.csv('clumped_cliques_analysis.tsv', header = T, sep = "\t")



for (i in 1:length(cliques$Clique)){
  print(i)
  clique = cliques[i,1]
  sub_g <-induced_subgraph(net, unlist(strsplit(toString(clique), "[_]")))
  for ( i in 1:length (E(net)[[]])){
    num <- paste0('0.',gsub('\\.','',toString(ifelse(E(net)[[i]]$weight/72*100*34/100>2, E(net)[[i]]$weight/72*100*34/100/10, round(E(net)[[i]]$weight/72*100*34/100,3)))))
    E(net)[[i]]$color <- adjustcolor("white",alpha.f = 0)
  }
  V(net)$color <- "white"
  E(net)[which(V(net) %in% V(net)[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))])%--%which(V(net) %in% V(net)[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))])]$weight <- min(E(net)[which(V(net) %in% V(net)[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))])%--%which(V(net) %in% V(net)[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))])]$weight)
  E(net)[which(V(net) %in% V(net)[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))])%--%which(V(net) %in% V(net)[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))])]$color <- "#17B4BB"
  
  V(net)$color[V(net)$name %in% unlist(strsplit(toString(clique), "[_]"))] <- adjustcolor("#17B4BB",alpha.f = .6)
  
  
  
  pdf(paste0(toString(clique), "graph.pdf"),4.5,4.5)
  plot.igraph(net,vertex.label=V(net)$name,
              layout=layout.circle,
              edge.width=E(net)$weight/9.5,
              vertex.size=6.5,vertex.label.cex = 0.24)
  dev.off()
}
  
#FA8E14

m %>% write.table(file = "Herit_cluster_graph_matrix__no_names_jaccard.tsv", quote = FALSE, sep = "\t",row.names = FALSE, col.names = FALSE)
test = read.csv('Herit_cluster_graph_matrix_jaccard.tsv', sep = "\t", header = TRUE)
plot(g3)
g3 <-induced_subgraph(net, c('11','36','153','163'))
E(net)$color<- "#5B0493"
V(net)$color <- "white"
E(net)[[c(1,5,6)]] <- "green"
plot.igraph(net,vertex.label=V(net)$name,
            layout=layout.circle,
            edge.width=E(net)$weight/4.5,
            vertex.size=7)
nums=c()
for (i in 1:length(E(net))){
  if (E(net)[[1]] %in% E(g3)[[]]){
    print(E(net)[[i]])
    nums = c(nums,i)
  }
}

which(V(net) %in% V(net)[V(net)$name %in% c('11','36','153','163')])

E(net)[which(V(net) %in% V(net)[V(net)$name %in% c('11','36','153','163')])%--%which(V(net) %in% V(net)[V(net)$name %in% c('11','36','153','163')])]$color <- "green"


induced.subgraph(graph=net,vids=unlist(neighborhood(graph=g2,order=1,nodes=DAT1)))
  
library(qgraph)
qgraph(m,line = 0.5)

groups <- list(A = c('65','21','57','179','195'),B <-colnames(m)[!(colnames(m) %in% c('65','21','57','179','195'))])
# Factor:
groups <- c(rep("A",length(c('65','21','57','179','195'))),rep("B",length(colnames(m)[!(colnames(m) %in% c('65','21','57','179','195'))])))
qgraph(m, posCol = "#5B0493",layout="circular",vsize = 6.5,groups=groups)



groups <- list(A = c(1,2,3,4,5),
               B = c(6,7,8,9,10))
# Factor:
groups <- c("A","A","A","A","A",
            "B","B","B","B","B")
# Result:
qgraph(matrix(1,10,10),layout="circular",groups=groups)



library(plotly)
dat <- data.matrix(read.table("filtered_heatmap_clump_CP.tsv", header = TRUE, row.names = 1,
                            sep = "\t"))
dat[is.na(dat)] <- Inf
dat = -log10(dat)
dat[dat>500] <- max(dat[is.finite(dat)])
dat[!is.finite(dat)] <- 0
non_zero_names <- colnames(dat[,which(colMeans(dat)!=0)])
zero_names <- colnames(dat[,which(colMeans(dat)==0)])
new_dat <- data.matrix(as.data.frame(dat)[,c(non_zero_names,zero_names)])
colfunc <- colorRampPalette(c("DarkRed",'LightCoral' , 'white'))
p <- plot_ly(z = dat, x = colnames(dat), y = rownames(dat), type = "heatmap", colors = colfunc(5), reversescale =T)
htmlwidgets::saveWidget(as_widget(p), file = "filtered_heatmap_clump_set.html")
rownames(dat)[21]='REACTOME_TRANSPORT_OF_GLUCOSE'
rownames_old <- rownames(dat)
rownames(dat)[31]='REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE'
rownames(dat)[4]='REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING'
rownames(dat)[44]='REACTOME_NGF_SIGNALLING'
rownames(dat)[20]='REACTOME_FATTY_ACID_METABOLISM'
rownames(dat)[43]='REACTOME_TRAF6_MEDIATED_INDUCTION_OF_NFKB'
rownames(dat)[34]='REACTOME_GASTRIN_CREB_SIGNALLING'
rownames(dat)[46]='REACTOME_PLATELET_ACTIVATION_SIGNALING'
rownames(dat)[25]='REACTOME_LIPID_DIGESTION'
rownames(dat)[33]='REACTOME_TRANSMEMBRANE_TRANSPORT'
rownames(dat)[5]='REACTOME_METABOLISM_OF_LIPID'
rownames(dat)[32]='REACTOME_SLC_MEDIATED_TRANSPORT'
rownames(dat)[8]='REACTOME_CHYLOMICRON_TRANSPORT'
rownames(dat)[29]='REACTOME_ABACAVIR_TRANSPORT'
rownames(dat)[11]='REACTOME_ORGANIC_ION_TRANSPORT'
rownames(dat)[13]='KEGG_T_CELL_RECEPTOR_PATHWAY'
rownames(dat)[34] ="REACTOME_INTEGRATION_OF_ENERGY"


rownames(dat)
dat <- dat[order(rownames(dat)),]
drawCoolHM(t(dat))

dat_l <- data.matrix(read.table("filtered_heatmap_clump_CP_loci.tsv", header = TRUE, row.names = 1,
                              sep = "\t"))
dat_l[is.na(dat_l)] <- 0
rownames(dat_l)[21]='REACTOME_TRANSPORT_OF_GLUCOSE'
dat_l <- dat_l[order(rownames(dat_l)),]


rectheat = as.vector(dat)
recs= as.vector(dat_l)
dataf$rowv = factor(rownames(dat))
dataf$columnv = factor(colnames(dat))                     
dataf <- data.frame (rowv = rep (rownames(dat), length(colnames(dat))), columnv = rep(colnames(dat), each = length(rownames(dat))),
                     rectheat, circlesize = recs,
                     circlefill =  rectheat)
dataf$columnv <- factor(dataf$columnv,levels(dataf$columnv)[fac_ord])

require(ggplot2)
ggplot(dataf, aes(y = rowv,
                  x = columnv)) +## global aes+ 
  geom_tile(fill = 'lightblue') + 
  ## to get the rect filled
  geom_point(aes(colour = circlefill, 
                 size =ifelse(is.finite(log10(circlesize)*10),log10(circlesize)*10,0)))  +    ## geom_point for circle illusion
  scale_color_gradient(low = "white",  
                       high = "#d50f0dff")+       ## color of the corresponding aes
  scale_size(range = c(0.01, 4))+             ## to tune the size of circles
  theme_bw()+
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  labs(x = NULL,y=NULL)+
  guides(
    colour = guide_colourbar(order = 1, barheight = 53,title.position = 'right',title	
='-log10 (GSEA FDR q-value)',title.theme=element_text(angle = 90,hjust = 0.5))
  )

nums = c()
for (x in colnames(dat)){
  nums = c(nums, substr(x, 2,5))
}
names = c()
for (num in sort(as.numeric(nums))){
  names = c(names, paste0("X",num))
}


fac_ord <- c()
for (x in names){
  fac_ord <- c(fac_ord, which(levels(factor(names))==x))
}


col.order <- c(names )
dat = dat[,col.order]
dat <- dat[order(rownames(dat)),]





library(lattice)
drawCoolHM = function(df){
  e = round(df, digits=3)
  myPanel_a <- function(x, y, z, ...) {
    panel.levelplot(x,y,z,...)
    # panel.text(x, y, e[cbind(x,y)]) ## use handy matrix indexing
  }
  return(levelplot(df, col.regions=jet,
                   at=seq(0, max(df), length.out=100),
                   aspect='fill', colorkey=list(width=3, labels=list(cex=1.0)),
                   scales=list(x=list(rot=45),tck = c(1,0)), xlab=list(label=''),
                   ylab=list(label=''), panel=myPanel_a))
}


mypal = colorRampPalette(c('white', '#ecad2f'))

mypal = colorRampPalette(c('#FFFB74', '#e000ba'))
mypal = colorRampPalette(matlab.like(100))
mypal = colorRampPalette(c('#FFFB74', '#8674FB'))
mypal = colorRampPalette(c('white', '#8674FB'))
jet =  colorRampPalette(blue2red(100))

mypal = colorRampPalette(c('white', '#d50f0dff'))
jet = mypal(100)
drawCoolHM(t(dat))






drawCoolHM(t(mixed_dat))
single <- dat
multi <- dat

single <- single[order(rownames(single)), ]
multi <- multi[order(rownames(multi )), ]

colnames(mixed_dat) <- colnames(dat) 
mixed_dat <- matrix(nrow=length(rownames(single))*2,ncol=length(colnames(single)))
i1=1
i2=1
for (i in 1:(length(rownames(single))*2)){
  if ((i %% 2) == 0){
    mixed_dat[i,]=single[i2,]
    i2=i2+1
    } else {
     mixed_dat[i,]=multi[i1,]
     i1=i1+1
   }
}
rownames(mixed_dat) <-rep(rownames(dat),2)

i1=1
i2=1
for (i in 1:length(rownames(mixed_dat))){
  if ((i %% 2) == 0){
    rownames(mixed_dat)[i] <- paste0(rownames(single)[i2], "_single_tissue")
    i2=i2+1
  } else {
    rownames(mixed_dat)[i] <- paste0(rownames(multi)[i1], "_multi_tissue")
    i1=i1+1
  }
}

pdf("SNP_tree2.pdf",20,20)
treemap(SNP_data,
        index=c('V1', 'V3'),
        vColor='V3',
        type="index",
        vSize='V2')
dev.off()

SNP_data <- read.csv('SNP_diagram_jaccard.tsv', header = FALSE, sep = "\t")%>% filter(V4>=3)
new_plot_names <-unique(as.numeric(unlist(as.vector(SNP_data %>% filter(V4>=3) %>% select(V1)))))
ggplot(SNP_data, aes(area = V3, fill = V1, label = V2,
                            subgroup = V1)) +
  geom_treemap() +
  geom_treemap_subgroup_border() +
  geom_treemap_subgroup_text(place = "centre", grow = T, alpha = 0.5, colour =
                               "black", fontface = "italic", min.size = 0) +
  geom_treemap_text(colour = "white", place = "topleft", reflow = T)


for(j in 1:length(new_plot_names)){
  out_name <- paste0('C:/Users/anton/Desktop/Project/plots',"/herit_cluster_SNP_proportions_",toString(new_plot_names[j]), ".png")
  to_plot <-SNP_data %>% arrange(V1) %>% filter(V1==new_plot_names [j]) 
  g <- ggplot(to_plot, aes(x = V3, y = V2, fill =V1)) + 
    theme_bw() + 
    guides(fill=FALSE)+
    geom_bar(stat = "identity")+labs(x = "Number of phenotype associacions", y = "Proportion") +
    facet_wrap(~V1, scales = "free_x")+
    theme( axis.text.x = element_text(angle = 90, hjust = 1),
           panel.background = element_blank(), 
           axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank())
  ggsave(filename=out_name,width=9, height=9, dpi=200)
}

eqtl_signif <- read.csv('herit_eqtls_for_signif_pairs.tsv', header = FALSE, sep = "\t") #mean=28 median =33 
summary(eqtl_signif$V6)
length(eqtl_signif[eqtl_signif$V6>0,]$V6) 
#1953/2191 - 0.8913738
#more than 50% - 1609/2191 - 0.7343679

eqtl_signif$fac <- 'pleiotropic'
eqtl_s$fac<- 'signifficant'
eqtl_all_s$fac<- 'all'


df_sub_eq = rbind(eqtl_signif[,c(6,7)],eqtl_s[,c(6,7)], eqtl_all_s[,c(6,7)])

ggplot(df_sub_eq, aes(x=fac, y=V6, fill=fac)) + 
  theme_bw()+
  scale_fill_manual(values=c("#F8D66C","#1DAEC5","#AF77BD"))+
  geom_violin(scale = "area", trim = FALSE, adjust=1.5, alpha=0.7)+
  xlab('group') + 
  ylab('Number of eqtls')+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())





library("ggpubr")
ggscatter(eqtl_signif, x = "V6", y = "V3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Number of eqtl enrichments", ylab = "Number of associations")+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

fisher.test(rbind(c(1953,247551-2191),c(2191-1953,469013-247551-2191)), alternative="greater")$p.value
fisher.test(rbind(c(1609,53655-2191),c(2191-1609,469013-53655-2191)), alternative="greater")$p.value
library(dplyr)
eqtl_s <- sample_n(eqtl,2192,replace = FALSE)
eqtl<- read.csv('herit_eqtls_for_clusters_pairs.tsv', header = FALSE, sep = "\t") #mean=7.227 median =1.000
summary(eqtl$V6)
length(eqtl[eqtl$V6>24,]$V6) #247551/469013 - 0.5278127
#more than 50% 53655/469013 - 0.1143998

ggscatter(eqtl, x = "V6", y = "V3", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Number of eqtl enrichments", ylab = "Number of associations")+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

head(eqtl_all)
eqtl_all<- read.table('herit_eqtls_for_all_pairs.tsv', header = FALSE, sep = "\t", 
           col.names = paste0("V",seq_len(6)), fill = TRUE)  #mean=1.672 median = 0
eqtl_all_s <- sample_n(eqtl_all,2192,replace = FALSE)

summary(eqtl_all$V6)
length(eqtl_all[eqtl_all$V6>0,]$V6) 
length(eqtl_all[eqtl_all$V6>24,]$V6) #2462298/10894596 - 0.226011
#more than 50% 210548/10894596 -  0.01932591

eqtl_all$fac <- as.factor(ifelse(eqtl_all$V6>24, ">50%", "<=50%"))
eqtl_signif$fac <- as.factor(ifelse(eqtl_signif$V6>24, ">50%", "<=50%"))
eqtl$fac <- as.factor(ifelse(eqtl$V6>24, ">50%", "<=50%"))

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}


aqv <- ggplot(eqtl_all, aes(x=fac, y=V3, fill=fac)) + 
  theme_bw()+
  ggtitle("All SNPs")+
  geom_violin(scale = "area", trim = FALSE, adjust=75)+
  scale_fill_manual(values=c("#1DAEC5","#AF77BD"))+
  xlab('Eqtl enrichment proportion') + 
  ylab('Number of associations')+ 
  guides(fill=F)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+ 
  stat_summary(fun.data=data_summary, alpha=0.7)



ssqv <- ggplot(eqtl_signif, aes(x=fac, y=V3, fill=fac)) + 
  theme_bw()+
  ggtitle("Pleiotropic signifficant SNPs")+
  scale_fill_manual(values=c("#1DAEC5","#AF77BD"))+
  geom_violin(scale = "area", trim = FALSE, adjust=1.45, alpha=0.7)+
  xlab('Eqtl enrichment proportion') + 
  ylab('Number of associations')+ 
  guides(fill=F)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+ 
  stat_summary(fun.data=data_summary, alpha=0.7)

sqv <- ggplot(eqtl, aes(x=fac, y=V3, fill=fac)) + 
  theme_bw()+
  ggtitle("All signifficant SNPs")+
  geom_violin(scale = "area", trim = FALSE, adjust=20, alpha=0.8)+
  scale_fill_manual(values=c("#1DAEC5","#AF77BD"))+
  xlab('Eqtl enrichment proportion') + 
  ylab('Number of associations')+ 
  guides(fill=F)+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank()) + 
  stat_summary(fun.data=data_summary, alpha=0.7)


lay <- rbind(c(1,1,1,2,3),
             c(1,1,1,4,5),
             c(6,7,8,9,9))
grid.arrange(grobs = gs, layout_matrix = lay)

grid.arrange(
  grobs = gl,
  widths = c(1, 1, 1),
  layout_matrix = rbind(c(aqv, NA, sqv),
                        c(NA, ssqv, NA))
)

grid.arrange(aqv, sqv, ssqv, nrow = 2)

ggplot(eqtl,aes(V3))+ 
  geom_histogram(data=subset(eqtl[eqtl$V3>0,],fac=='<=50%'),aes(fill=fac),alpha=0.5,stat = "count")+
  geom_histogram(data=subset(eqtl[eqtl$V3>0,],fac=='>50%'),aes(fill=fac),alpha=0.6,stat = "count")+
  scale_fill_manual(name="group", values=c("#1DAEC5","#AF77BD"),labels=c("<=50%",">50%"))+
  labs(x = "Number of associated clusters", y = "Occurrence frequency") +theme_bw()+
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"))

ggplot(eqtl_signif,aes(V3))+
  geom_point(data=subset(eqtl_signif[eqtl_signif$V3>0,],fac=='<=50%'),aes(col=fac),alpha=0.5,stat = "count")+
  geom_point(data=subset(eqtl_signif[eqtl_signif$V3>0,],fac=='>50%'),aes(col=fac),alpha=0.5,stat = "count")+
  scale_fill_manual(name="group", values=c("#1DAEC5","#AF77BD"),labels=c("<=50%",">50%"))+theme_bw()+
  labs(x = "Number of associated clusters", y = "Occurrence frequency")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(eqtl,aes(V3))+
  geom_density(data=subset(eqtl[eqtl$V3>0,],fac=='<=50%'),adjust = 15,aes(fill=fac),alpha=0.5)+
  geom_density(data=subset(eqtl[eqtl$V3>0,],fac=='>50%'),adjust = 15,aes(fill=fac),alpha=0.5)+
  scale_fill_manual(name="group", values=c("#1DAEC5","#AF77BD"),labels=c("<=50%",">50%"))+theme_bw()+
  labs(x = "Number of associated clusters")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

Graph_anal <- read.csv('clump_herit_graph_analysis.tsv', header = FALSE, sep = "\t")
Graph_anal <- Graph_anal %>% arrange(desc(V2))
Graph_anal$V1 <- as.factor(Graph_anal$V1)

dd_sub = Graph_anal[,c(1, 3,2,4)]
library(reshape2)
dd = melt(dd_sub, id=c('V1'))



ggplot(Graph_anal) + 
  theme_bw() +
  geom_bar(aes(x = reorder(V1, -V2), y = V2), stat = "identity",alpha=0.8, fill = "#AF77BD")+
  geom_bar(aes(x = reorder(V1, -V4), y = V4), stat = "identity",alpha=0.5,fill = "#1DAEC5")+
  geom_bar(aes(x = reorder(V1, -V3), y = V3), stat = "identity",alpha=0.1,fill = "green")+
  labs(x = "Cluster", y = "Number of edges/Mean weigth of edges") +
  scale_fill_manual("legend", values=c("grey", "darkgreen", "green"),labels=c("Number of edges", "Mean", "Treatment 2")) +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

ggplot(dd,aes(x = reorder(V1, -value), y = value, fill=variable)) + 
  theme_bw() +
  stat_summary(stat='identity', geom="bar",alpha=0.3)+
  labs(x = "Cluster", y = "Node properties") +
  scale_fill_manual("legend", values=c("grey", "darkgreen", "green"),labels=c("Sum of edge weights", "Number of edges", "Mean edge weight")) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())




Graph_anal$num_per <- Graph_anal$V2/62*100
Graph_anal$num_w <- Graph_anal$V3/sum(Graph_anal$V3)*100
dd_sub_2 = Graph_anal[,c(1, 6,7)]
library(reshape2)
dd2 = melt(dd_sub_2, id=c('V1'))

ggplot(dd2,aes(x = reorder(V1, -value), y = value, fill=variable)) + 
  theme_bw() +
  stat_summary(stat='identity', geom="bar",alpha=0.7)+
  labs(x = "Cluster", y = "Node properties") +
  scale_fill_manual("legend", values=c("grey", "darkgreen"),labels=c("Number of edge percent", "Total weight percent")) +
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())



all_vars <- read.csv('annotated_clustered_jaccard.tsv', header = FALSE, sep = "\t")
var_stats <- as.data.frame(table(all_vars$V1))

ggplot(data=var_stats, aes(x=Var1))+ 
  theme_bw() + 
  guides(fill=FALSE)+ 
  geom_bar(aes(x = reorder(Var1, -Freq), y = Freq), stat = "identity",alpha=0.7, fill = "#AF77BD")+
  labs(x = "Variant type", y = "Frequency") +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

pl_vars <- read.csv('pleio_all_stats.tsv', header = FALSE, sep = "\t")
pl_var_stats <- as.data.frame(table(pl_vars$V1))

ggplot(data=pl_var_stats, aes(x=Var1))+ 
  theme_bw() + 
  guides(fill=FALSE)+ 
  geom_bar(aes(x = reorder(Var1, -Freq), y = Freq), stat = "identity",alpha=0.7, fill = "#AF77BD")+
  labs(x = "Variant type", y = "Frequency") +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

#Missense: 3073/469013*100 46/2191*100

fisher.test(rbind(c(46,3073),c(2191-46,469013-3073)), alternative="greater")$p.value
fisher.test(rbind(c(4,633),c(2191-4,469013-633)), alternative="greater")$p.value
 
#splice - 4/2191 633/469013
#stop gained - 0/2191 33/469013


df <- merge(var_stats, pl_var_stats, by="Var1", all = T)
library(dplyr)
df <- df %>%
  mutate(Freq.y = replace(Freq.y,is.na(Freq.y),0))
df$Freq.x_percent <- df$Freq.x/sum(df$Freq.x)*100
df$Freq.y_percent <- df$Freq.y/sum(df$Freq.y)*100

chisq.test(df[, c(2,3)]) 

dd_sub_3 = df[,c(1, 4,5)]
library(reshape2)
dd3 = melt(dd_sub_3, id=c('Var1'))

ggplot(data=dd3, aes(x = reorder(Var1, - value), y=value, fill=variable))+ 
  theme_bw()+
  stat_summary(stat='identity', geom="bar",alpha=0.45)+
  scale_fill_manual("legend", values=c("#AF77BD", "#1DAEC5"),labels=c("All signifficant SNP", "pleiotropic SNP")) +
  labs(x = "Variant type", y = "Percentage") +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())

cliques <- read.csv('clumped_cliques_analysis.tsv', header = T, sep = "\t")
colnames(cliques)
ggplot(cliques, aes(x=Num_of_nodes, y=Num_of_SNPs)) +
  geom_point(size=2, shape=23,  fill ="#1DAEC5")+ 
  theme_bw() + 
  guides(fill=FALSE)+
  labs(x = "Number of nodes in clique", y = "Number of SNPs") +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+geom_smooth()+ylim(c(0,max(cliques$Num_of_SNPs)))
library(ggpubr)
ggscatter(cliques, x = 'Num_of_nodes', y = "Num_of_clumps", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "Number of edges", ylab = "Number of genes")+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())
  
cor.test(cliques$V4, cliques$V2,  method = "pearson", use = "complete.obs")$p.value


pleoi <- read.csv('Enrich_MIR.tsv.csv', header = FALSE, sep = "\t", stringsAsFactors = F)
pleoi[7,1] <- 'REACTOME_FACTORS_INVOLVED_IN_MEGAKARYOCYTE_DEVELOPMENT'
pleoi[12,1] <- 'REACTOME_NEUROTRANSMITTER_RECEPTOR_BINDING'
pleio$V1 <- as.factor(pleio$V1)
ggplot(data=pleoi, aes(x = reorder(V1, -log10(V2)), y=-log10(V2))) +
  geom_bar(stat="identity",fill="blue", alpha=0.5)+ coord_flip()+ 
  theme_bw() +
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  labs( y = "-log10(qvalue)")



library("xlsx")
library(dplyr)
hs<- read.csv('ces_52.metrics', header = T, sep = "\t", stringsAsFactors = F)
hs_sub <- hs %>% select(TOTAL_READS, PCT_PF_UQ_READS,MEAN_BAIT_COVERAGE,PCT_TARGET_BASES_10X,PCT_TARGET_BASES_20X,PCT_TARGET_BASES_30X,PCT_TARGET_BASES_40X,PCT_TARGET_BASES_50X,PCT_TARGET_BASES_100X)



#hs_sub$PCT_PF_UQ_READS <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_PF_UQ_READS)))
#hs_sub$MEAN_BAIT_COVERAGE <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$MEAN_BAIT_COVERAGE)))
#hs_sub$PCT_TARGET_BASES_10X <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_TARGET_BASES_10X)))
#hs_sub$PCT_TARGET_BASES_20X <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_TARGET_BASES_20X)))
#hs_sub$PCT_TARGET_BASES_30X <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_TARGET_BASES_30X)))
#hs_sub$PCT_TARGET_BASES_40X <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_TARGET_BASES_40X)))
#hs_sub$PCT_TARGET_BASES_50X <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_TARGET_BASES_50X)))
#hs_sub$PCT_TARGET_BASES_100X <- as.numeric(gsub(",", ".", gsub("\\.", "", hs_sub$PCT_TARGET_BASES_100X)))


hs_sub$Dups <- 1-as.numeric(hs_sub$PCT_PF_UQ_READS)
colnames(hs_sub) <- c("number_of_reads", "unique_reads_percent", 
                      "Mean_coverage", "PCT_TARGET_BASES_10X", 
                      "PCT_TARGET_BASES_20X",  "PCT_TARGET_BASES_30X",
                     "PCT_TARGET_BASES_40X",  
                     "PCT_TARGET_BASES_50X", 
                      "PCT_TARGET_BASES_100X", "duplicates_percent")
hs_sub <-hs_sub[,c("number_of_reads", "unique_reads_percent", "duplicates_percent", 
                   "Mean_coverage", "PCT_TARGET_BASES_10X", 
                   "PCT_TARGET_BASES_20X",  "PCT_TARGET_BASES_30X",
                   "PCT_TARGET_BASES_40X",  
                   "PCT_TARGET_BASES_50X", 
                   "PCT_TARGET_BASES_100X")]
rownames(hs_sub) <- gsub("ces_52.sample_", "", hs$SAMPLE)
#rownames(hs_sub) <- hs$SAMPLE
write.xlsx(hs_sub, 'ces_52.metrics.xlsx', sheetName = "Sheet1", 
           col.names = TRUE, row.names = TRUE, append = FALSE)

#summary(hs_sub$Mean_coverage)

to_xls = read.csv('ces_52.Final.csv', header = T, sep = ",") 
#xls2 <- to_xls[to_xls$EffType != 'INTER',]
#xls3 <- xls2[xls2$EffType != 'INTRN',]
library("xlsx")
rownames(to_xls) <- NULL
write.xlsx(to_xls, file ='ces_52.Final.xls', sheetName = "Sheet1", 
           col.names = T, row.names = F, append = FALSE)



for (i in 1:6){
  print(paste0('ces_50.',i,'.csv'))
  to_xls = read.csv(paste0('ces_50.',i,'.csv'), header = F, sep = ",")
  rownames(to_xls) <- NULL
  write.xlsx(to_xls, file =paste0('ces_50.',i,'.xls'), sheetName = "Sheet1", 
             col.names = F, row.names = F, append = FALSE)
}


