library(ggplot2)
file = read.csv('herit_Manhattan_09_no_0.tsv', header = FALSE, sep = "\t")
file2 <- read.csv('herit_Clustered_mahattan_jaccard.tsv', header = FALSE, sep = "\t")
file3 <- read.csv('for_all_SNP_stats.tsv', header = FALSE, sep = "\t")
colnames(file)[5] <- c("number_of_associations")
colnames(file2)[5] <- c("number_of_associations")
colnames(file3)[5] <- c("number_of_associations")
file$group <- 'Non-clustered data'
file2$group <- 'Clustered data'


sum(file$number_of_associations>1) #230296
sum(file2$number_of_associations>1) #92224

Merged_data_set <- rbind(file[,c(5,6)],file2[,c(5,6)])
data_set_for_clust_and_phen <- cbind(file[,c(1:5)],file2[,c(4,5)])
colnames(data_set_for_clust_and_phen) <- c("CHR", "BP", "SNP", "Phens","NP", "clust", "NC")

ggplot(data_set_for_clust_and_phen , aes(NC, NP)) +
  geom_count(alpha = 0.6, col = '#1DAEC5',show.legend=F) +
  theme_bw()+ 
  stat_smooth(method = 'lm', col = '#AF77BD',alpha = 0.8)+
  labs(x = "Number of associated clusters") + 
  labs(y = "Number of associated phenotypes") +  
  ylim(min(data_set_for_clust_and_phen$NP), max(data_set_for_clust_and_phen$NP))+
  xlim(min(data_set_for_clust_and_phen$NC), max(data_set_for_clust_and_phen$NC))+
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())


#single data plots
##Density plots
ggplot(file, aes(number_of_associations)) +
  geom_density(adjust = 15, alpha = 0.3, fill = '#1DAEC5')+
  labs(x = "Number of SNP associations") + 
  theme_bw()+
  scale_x_continuous(limits = c(1,52))+
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

ggplot(file2, aes(number_of_associations)) +
  geom_density(adjust = 15, alpha = 0.3, fill = 'blue')+
  labs(x = "Number of associated clusters")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

##Geom points
ggplot(file, aes(number_of_associations)) +
  geom_point(alpha = 0.6, col = '#1DAEC5', stat = "count") + 
  theme_bw()+
  scale_x_continuous(limits = c(1,52))+
  labs(x = "Number of SNP associations", y = "Occurrence frequency")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggplot(file2, aes(number_of_associations)) +
  geom_point(alpha = 0.6, col = 'blue', stat = "count")+
  labs(x = "Number of associated clusters", y = "Occurrence frequency")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
##Histograms
ggplot(file, aes(number_of_associations)) + 
  theme_bw()+
  scale_x_continuous(limits = c(0,52))+
  geom_histogram(alpha = 0.5, fill = '#1DAEC5', stat = "count")+
  labs(x = "Number of SNP associations", y = "Occurrence frequency")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggplot(file2, aes(number_of_associations)) +
  geom_histogram(alpha = 0.5, fill = 'blue', stat = "count")+
  labs(x = "Number of associated clusters", y = "Occurrence frequency")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())






#multi data plots
##histogram
ggplot(Merged_data_set[Merged_data_set$number_of_associations>0,],aes(number_of_associations))+  
  theme_bw()+
  scale_x_continuous(limits = c(3,52))+
  geom_histogram(data=subset(Merged_data_set[Merged_data_set$number_of_associations>3,],group=='Non-clustered data'),aes(fill=group),alpha=0.5,stat = "count")+
  geom_histogram(data=subset(Merged_data_set[Merged_data_set$number_of_associations>3,],group=='Clustered data'),aes(fill=group),alpha=0.5,stat = "count")+
  scale_fill_manual(name="group", values=c("#1DAEC5","#AF77BD"),labels=c("Clustered","Non-clustered"))+
  labs(x = "Number of associated phenotypes/clusters", y = "Occurrence frequency") +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black"), legend.position = 'top')
##geom_point
ggplot(Merged_data_set, aes(number_of_associations))+  
  theme_bw()+
  scale_x_continuous(limits = c(3,52))+
  geom_point(data=subset(Merged_data_set[Merged_data_set$number_of_associations>3,],group=='Non-clustered data'),aes(col=group),alpha=0.5,stat = "count")+
  geom_point(data=subset(Merged_data_set[Merged_data_set$number_of_associations>3,],group=='Clustered data'),aes(col=group),alpha=0.5,stat = "count")+
  scale_fill_manual(name="group", values=c("#1DAEC5","#AF77BD"),labels=c("Clustered","Non-clustered"))+
  labs(x = "Number of associated phenotypes/clusters", y = "Occurrence frequency")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = 'top')
##density plot
ggplot(Merged_data_set, aes(number_of_associations))+  
  theme_bw()+
  scale_x_continuous(limits = c(3,52))+
  geom_density(data=subset(Merged_data_set[Merged_data_set$number_of_associations>3,],group=='Non-clustered data'), adjust = 15, aes(fill=group),alpha=0.5)+
  geom_density(data=subset(Merged_data_set[Merged_data_set$number_of_associations>3,],group=='Clustered data'), adjust = 15, aes(fill=group),alpha=0.615)+
  scale_fill_manual(name="group", values=c("#1DAEC5","#AF77BD"),labels=c("Clustered","Non-clustered"))+
  labs(x = "Number of associated phenotypes/clusters")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), legend.position = 'top')


ggplot(Merged_data_set, aes(x=group, y=number_of_associations, fill=group)) + 
  theme_bw()+
  scale_fill_manual(values=c("#1DAEC5","#AF77BD"))+
  geom_violin(scale = "area", trim = FALSE, adjust=8.45, alpha=0.7)+
  xlab('Data type') + 
  ylab('Number of associations')+ 
  guides(fill=F)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())



Clust_file <- read.csv('herit_for_clustering_jaccard.tsv', header = TRUE, sep = "\t")

ggplot(Clust_file,aes(value)) + 
  theme_bw()+ 
  geom_histogram( fill = "#8E9323", alpha = 0.46,stat = "count")+ 
  labs(x = "Cluster", y = "Number of phenotypes")+ 
  theme(panel.background = element_blank(), 
        axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank()) 

sum(file$V5>1) # 230296/468013
sum(file$V5>30) #7597/468013

mean(file2$V5)
sum(file2$V5>1) # 92224/468013
sum(file2$V5>10)  # 64/468013*100



enrich_dict_CP <- read.csv('enrichment_stats_CP.tsv', header = TRUE, sep = "\t")
colnames(enrich_dict_CP) <- c("pathway","Phenotypes","Clusters","Unions_of_enrichments","Enrichments_of_unions")
library(tidyr)
data_long <- gather(enrich_dict_CP,pathway, Phenotypes:Enrichments_of_unions, factor_key=TRUE)
colnames(data_long) <- c('group', 'value')

ggplot(data_long, aes(x=group, y=value, fill=group)) + 
  theme_bw()+
  geom_violin(scale = "area", trim = FALSE, adjust=1.5, alpha=0.7)+
  scale_fill_manual(values=c("#00ffff","#bf00ff","#ffff00", "#ffbf00"))+
  guides(fill=F)+
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+ 
  ylab('Number of significant clusters/phenotypes') +  geom_violin(data =data_long_t, aes(x=group, y=value, fill=group),scale = "area", trim = FALSE, adjust=2.5, alpha=0.1)

enrich_dict_tot <- read.csv('enrichment_stats_total.tsv', header = TRUE, sep = "\t")
colnames(enrich_dict_tot) <- c("pathway","Phenotypes","Clusters","Unions_of_enrichments","Enrichments_of_unions")
data_long_t <- gather(enrich_dict_tot,pathway,Phenotypes:Enrichments_of_unions, factor_key=TRUE)
colnames(data_long_t) <- c('group', 'value')
data_long_t$fac2 <- ifelse(data_long_t$group %in% c("Phenotypes","Clusters"),1,2)
ggplot(data_long_t, aes(x=value, fill=group)) + 
  theme_bw()+
  geom_density(scale = "area", alpha=0.7)+
  scale_fill_manual(values=c("#00ffff","#bf00ff","#ffff00", "#ffbf00"))+
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+ 
  xlab('Number of significant clusters/phenotypes') +
  facet_grid(~fac2, scales = "free_x")


require(ggplot2)


ggplot(data_long_t, aes(x=value, fill=group)) + 
  theme_bw()+
  geom_bar(stat="count")+
  scale_fill_manual(values=c("#00ffff","#bf00ff","#ffff00", "#ffbf00"))+
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+ 
  xlab('Number of significant clusters/phenotypes') +
  ylab('Frequency (counts)') +
  facet_wrap(~fac2,   scales = "free")


ggplot(data_long_t, aes(x=value, fill=group)) + 
  theme_bw()+
  geom_histogram(aes(y = ..density..),binwidth = 1)+
  scale_fill_manual(values=c("#00ffff","#bf00ff","#ffff00", "#ffbf00"))+
  theme( axis.text.x = element_text(angle = 45, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+ 
  xlab('Number of significant clusters/phenotypes') +
  ylab('Frequency') +
  facet_wrap(~fac2,   scales = "free")


sum(data_long_t[data_long_t$group=='Phenotypes',2]) #403713
sum(data_long_t[data_long_t$group=='Clusters',2]) #98695
sum(data_long_t[data_long_t$group=='Unions_of_enrichments',2]) #38588
sum(data_long_t[data_long_t$group=='Enrichments_of_unions',2]) #29159
