
#箱线图
rm(list=ls())
library(lme4)
library(lmerTest)

data<- read.table("all.txt", header = TRUE) 
library(ggnewscale)
library(ggplot2)
library(ggbreak)
library(gg.gap)
l<-list()##循环时，建新表
positions <- c("Seeding","Streching","Flowering")
#这里的i in ***你自己看一下要做哪个指标，可以随意改，例如 i in c(5, 8, 9, 10)
for (i in 2){
data1 <- data[,c(i,2)]
p = ggplot(data, aes(x=treatment, y= data[[i]], fill=treatment)) +  
  geom_boxplot(alpha=0.6, outlier.size=0, size=0.7, width=0.8) +
  scale_x_discrete(limits = positions)+
  scale_fill_manual(values=c("Seeding"="#36BED9","Streching"="#FBAD01","Flowering"="#FF6B70"))+
  #添加点 
  new_scale_fill() +#很重要，加另外的图,颜色重置
  geom_jitter(position=position_jitter(0.17),shape = 21,size=2.5, alpha=1,colour = "black",aes(fill = treatment))+                             
  scale_fill_manual(values=c("Seeding"="#36BED9","Streching"="#FBAD01","Flowering"="#FF6B70"))+  
  theme_bw() +
  xlab("")+
  ylab(colnames(data1)[1])+
  theme(legend.position="none")+
  theme(#panel.grid.major = element_blank(),
    #panel.grid.minor = element_blank(),
    strip.text = element_text(size=8,face="plain",color="black"),
    text=element_text(size=8,face="plain",color="black"),
    axis.title=element_text(size=12,face="plain",color="black"),
    axis.text = element_text(size=12,face="plain",color="black")
  )
p

l[[i-1]] <- ggplot_build(p)
}

library(cowplot)#组图
q <- plot_grid(plotlist = lapply(l, ggplot_gtable))
q
#q
ggsave("box_plot2.pdf",q,width = 5,height = 5)

#柱状图
#非参检验
##Kruskal-Wallis和wilcox test检验（整体差异分析和两两组间差异分析）
##非参数组间差异的Kruskal-Wallis检验
boxplerk <-  function(X,
                      Y,
                      n_permutations = 1000,
                      main = NULL,
                      xlab = NULL,
                      ylab = NULL,
                      bcol = "bisque",
                      p.adj = "none",
                      cexy = 1,
                      varwidth = TRUE,
                      las = 1,
                      paired = FALSE)
{
  aa <- levels(as.factor(Y))
  an <- as.character(c(1:length(aa)))
  tt1 <- matrix(nrow = length(aa), ncol = 7)    
  for (i in 1:length(aa))
  {
    temp <- X[Y == aa[i]]
    tt1[i, 1] <- mean(temp, na.rm = TRUE)
    tt1[i, 2] <- sd(temp, na.rm = TRUE) / sqrt(length(temp))
    tt1[i, 3] <- sd(temp, na.rm = TRUE)
    tt1[i, 4] <- min(temp, na.rm = TRUE)
    tt1[i, 5] <- max(temp, na.rm = TRUE)
    tt1[i, 6] <- median(temp, na.rm = TRUE)
    tt1[i, 7] <- length(temp)
  }
  
  tt1 <- as.data.frame(tt1)
  row.names(tt1) <- aa
  colnames(tt1) <- c("mean", "se", "sd", "min", "max", "median", "n")
  
  boxplot(
    X ~ Y,
    main = main,
    xlab = xlab,
    ylab = ylab,
    las = las,
    col = bcol,
    cex.axis = cexy,
    cex.lab = cexy,
    varwidth = varwidth
  )    
  require(agricolae)
  Yn <- factor(Y, labels = an)
  comp <- kruskal(X, Yn, p.adj = p.adj)
  sig <- "ns"
  
  if (paired == TRUE & length(aa) == 2)
  {
    coms <- wilcox.test(X ~ Yn, paired = TRUE)
    pp <- coms$p.value
  }    else
  {
    pp <- comp$statistics$p.chisq
  }    
  if(pp <= 0.1)
    sig <- "."
  if(pp <= 0.05)
    sig <- "*"
  if(pp <= 0.01)
    sig <- "**"
  if(pp <= 0.001)
    sig <- "***"
  
  gror <- comp$groups[order(rownames(comp$groups)), ]
  tt1$rank <- gror$X
  tt1$group <- gror$groups
  mtext(
    sig,
    side = 3,
    line = 0.5,
    adj = 0,
    cex = 2,
    font = 1
  )   
  if(pp <= 0.1)
    mtext(
      tt1$group,
      side = 3,
      at = c(1:length(aa)),
      line = 0.5,
      cex = 1,
      font = 4
    )
  
  list(comparison = tt1, p.value = pp)
  
}

#上面不要动，除非你会改

l<-list()##循环时，建新表
for (i in 2)
 {a=boxplerk(data[[i]],data$treatment, 
           ylab = colnames(data)[i], xlab = "Compartment",
           bcol = "bisque",p.adj = "fdr",las = 1) ##对p值进行“fdr”校准 
##用wilcox test标两两显著性
stat = a$comparison
p_value = a$p.value
stat[,c("p_value")] =p_value
stat[,c("index")] = colnames(data)[i]
stat$treatment=rownames(stat)
l[[i]]<-stat
}
f=do.call(rbind,l)#合并生成的新表

positions <- c("Seeding","Streching","Flowering")
library(ggplot2)
p = ggplot(f, aes(x=treatment,y=mean)) + #
  geom_bar(stat = "identity",position = "dodge",width = 0.8,color = "black",aes(fill=treatment),alpha=0.6) +
  scale_x_discrete(limits = positions)+
  facet_wrap(~index,scales= "free",nrow=4) + 
  scale_fill_manual(values=rev(c("Seeding"="#36BED9","Streching"="#FBAD01","Flowering"="#FF6B70")))+
  geom_errorbar(aes(ymin=mean-se, ymax=mean +se),position=position_dodge(.8), width=.2,cex=0.5)+
  geom_text(aes(label=group),color="black",size=5,vjust=0,hjust=0.5)+
  labs(x="", y="mean") +
  theme_bw()+
  theme(panel.grid.major=element_blank(),
        panel.grid.minor=element_blank())+
  theme(axis.text=element_text(colour='black',size=8))
p 
ggsave("bar_plot2.pdf",p,width = 15,height = 12)

