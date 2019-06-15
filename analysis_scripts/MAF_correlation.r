library("ggpubr")
library(ggplot2)
library(dplyr)
coefficient_file = read.csv('herit_Var_MAF_jaccard_all_SNP.tsv', header = FALSE, sep = "\t")
colnames(coefficient_file) = c('rsid','MAF','Num_of_associations')
coefficient_file$MAF <- ifelse(coefficient_file$MAF<0.5, coefficient_file$MAF, 1-coefficient_file$MAF)


#scatter plot
ggscatter(coefficient_file, x = "MAF", y = "Num_of_associations", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = "MAF", ylab = "Number of associations")+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank()) 
str(coefficient_file)

#correlation test
res <- cor.test(coefficient_file$Num_of_associations, coefficient_file$MAF, method = "pearson")
res$estimate

#geom violins
coefficient_file$fac <- as.factor(ifelse(coefficient_file$Num_of_associations<15, "<15", ">=15"))
ggplot(coefficient_file, aes(x=fac, y=MAF, fill=fac)) + 
  geom_violin(scale = "area", trim = FALSE, adjust=1.8)+
  xlab('Number of associations') + 
  ylab('MAF')+ 
  guides(fill=F)+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())
  
max(coefficient_file$MAF)
min(coefficient_file$MAF)

wilcox.test(MAF ~ fac, data=coefficient_file) #W = 8.8255e+10, p-value < 2.2e-16
coefficient_file$Allele_frequency <- ifelse(coefficient_file$MAF<0.005, "MAF<0.005", ifelse(coefficient_file$MAF>0.05, "MAF>0.05","0.005<=MAF<=0.05"))

#Facet graphs with MAF
ggplot(coefficient_file, aes(Num_of_associations,fill = Allele_frequency)) +
  geom_density(adjust = 15, alpha = 0.5)+
  labs(x = "Number of associations")+ 
  guides(fill=guide_legend(title="Allele frequencies"))+
  facet_wrap(~Allele_frequency, nrow=2, ncol=2, scales = "free_x")+
    theme( axis.text.x = element_text(angle = 90, hjust = 1),
           panel.background = element_blank(), 
           axis.line = element_line(colour = "black"),
           panel.grid.major = element_blank(), 
           panel.grid.minor = element_blank())


#lm Graph
g1 <- ggplot(coefficient_file, 
       aes(y = Num_of_associations, x = -log10(MAF))) + geom_point(color ='grey', alpha = 0.5)+geom_smooth(method = "lm")+
  labs(y = "Number of cluster associations")


# New stuff

coefficient_file = read.csv('herit_Clust_Var_MAF_jaccard_all_SNP.tsv', header = FALSE, sep = "\t")
colnames(coefficient_file) = c('rsid','MAF','Num_of_associations')
coefficient_file$MAF <- ifelse(coefficient_file$MAF<0.5, coefficient_file$MAF, 1 - coefficient_file$MAF)
coefficient_file$maf_bin = floor(coefficient_file$MAF / 0.001)
a <- by(data = coefficient_file$Num_of_associations, INDICES = coefficient_file$maf_bin, mean, simplify = T)

new_df = data.frame(mean_k = as.numeric(a), maf_bin = 1:500)
ggplot(new_df, aes(x = maf_bin, y = mean_k)) + 
  geom_line()+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  labs(x = "MAF bin", y ="Mean k")

coefficient_file$kexp = new_df[coefficient_file$maf_bin, 'mean_k']
coefficient_file$pval = ppois(coefficient_file$Num_of_associations-1, lambda = coefficient_file$kexp, lower.tail = F)

ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlab=expression(Expected~~-log[10](italic(p))),
       ylab=expression(Observed~~-log[10](italic(p))),
       xlim=c(0,max(e)), ylim=c(0,max(o)))
  lines(e,e,col="red")
}


png('new_qqplot.png')
ggd.qqplot(coefficient_file$pval)
dev.off()

signif = coefficient_file[coefficient_file$pval < 5e-9, ]
signif %>% select(rsid, Num_of_associations, pval) %>% arrange(pval) %>% write.table(file = "herit_valuable_SNP_jaccard.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
coefficient_file %>% filter(Num_of_associations>0) %>% write.table(file = "herit_pval_SNP_jaccard.tsv", quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)




coefficient_file$greater = ifelse(coefficient_file$Num_of_associations > 0, 1, 0)
b <- by(data = coefficient_file$greater, INDICES = coefficient_file$maf_bin, mean, simplify = T)
new_df = data.frame(mean_k = as.numeric(b), maf_bin = 1:500)

ggplot(new_df, aes(x = maf_bin, y = mean_k)) + geom_line()+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  labs(x = "MAF bin", y ="Mean k")


new_df$obs = as.numeric(by(data = coefficient_file$greater, 
                           INDICES = coefficient_file$maf_bin, mean, simplify = T))
new_df$exp = sapply(new_df$mean_k, function(x) mean(rpois(10000, lambda=x) > 0))

library(car)
qqPlot(-log10(new_df$exp))
library(reshape2)

tpl = melt(new_df[, c('obs', 'exp', 'maf_bin')], id.vars='maf_bin')
ggplot(tpl, aes(x = maf_bin, y = value, col=variable)) + 
  geom_line(lwd=1.2)+theme_bw()+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  labs(x = "MAF bin", y ="Mean k")

new_df$error = new_df$obs - new_df$exp

er<-ggplot(new_df, aes(x=maf_bin, y=error)) + geom_line() + scale_y_continuous(limits=c(-0.05, 0.10))+
  theme( axis.text.x = element_text(angle = 90, hjust = 1),
         panel.background = element_blank(), 
         axis.line = element_line(colour = "black"),
         panel.grid.major = element_blank(), 
         panel.grid.minor = element_blank())+
  labs(x = "MAF bin", y ="Error")

new_df$z = (new_df$obs - mean(new_df$obs))^2
new_df$resid = new_df$error ^ 2

r_squared = 1 - sum(new_df$resid)/sum(new_df$z) # R=0.9818536