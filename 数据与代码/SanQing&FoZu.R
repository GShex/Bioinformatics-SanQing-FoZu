GSE44077_expre <-read.csv('GSE44077_series_matrix.csv', ,header=TRUE,row.names=1) #或者采用下面的命令读入亦可
show(GSE44077_expre) #可以查看读入的数据结果究竟是怎么样的。后面的每个命令操作过程，均可用这个命令进行查验。这也是初学者必备的素质。
#GSE44077_expre<-read.table('GSE44077_matrix.csv',header=TRUE,row.names=1) #需要注意
#GSE44077_log2 <- log2(GSE44077_expre) #将表达值进行对数转换
boxplot(data.frame(GSE44077_expre),col="red")

library('Biobase')
library('limma')
library('ggplot2')
GSE44077_expr_N <- GSE44077_expre

GSE44077_expr_N <- normalizeBetweenArrays(GSE44077_expre)
#GSE44077_N_log2 <- log2(GSE44077_expr_N)
boxplot(data.frame(GSE44077_expr_N),col="red")
GSE44077_phenodata <- read.table('GSE44077_type.csv', sep = ',', header=TRUE, row.names=1, check.names = FALSE)

GSE44077_group <- GSE44077_phenodata [,1] #获取第二列的表型标签
GSE44077_design <- model.matrix(~0+factor(GSE44077_group)) #将表型标签进行数字0和1的转换
head(GSE44077_design) #查看表型标签进行数字转换后，整个表型表变换的结果，它和show()的结果是异曲同工的
colnames(GSE44077_design) <- levels(factor(GSE44077_group)) #将表型标签作为列
rownames(GSE44077_design) <- colnames(GSE44077_expr_N) #将对数表达值作为行
GSE44077_cont.matrix <- makeContrasts(paste0(unique(GSE44077_group),collapse="-"),levels=GSE44077_design)


GSE44077_fit <- lmFit(GSE44077_expr_N,GSE44077_design) #给定系列数据，为每个基因进行线性模型拟合）
GSE44077_fit2=contrasts.fit(GSE44077_fit,GSE44077_cont.matrix) #基于拟合线性模型进行系数估计和标准误差计算
GSE44077_fit3 <- eBayes(GSE44077_fit2)#基于经验贝叶斯进行进行t统计量、F统计量和差异对数表达计算
GSE44077_output = topTable(GSE44077_fit3,coef=1,n=Inf,adjust="BH") #从线性模型拟合结果中，进行排序，并提取排序靠前的基因列表
GSE44077_nrDEG = na.omit(GSE44077_output) #删除基因排序列表中出现的缺失值，并将结果返回
head(GSE44077_nrDEG)
write.table(GSE44077_nrDEG, "GSE44077_DE.csv", quote=F, sep=",") #该命令加入了相关文档数字分隔的字符定义
GSE33532_alldiff<-read.csv("GSE33532_DE.csv",header=TRUE,row.names=1,check.names=FALSE)
GSE33532_alldiff$threshold[GSE33532_alldiff$P.Value < 0.05 & GSE33532_alldiff$logFC > 1] = "up"
GSE33532_alldiff$threshold[GSE33532_alldiff$P.Value < 0.05 & GSE33532_alldiff$logFC < -1] = "down"
GSE33532_alldiff$threshold[GSE33532_alldiff$P.Value > 0.05 | (GSE33532_alldiff$logFC >= -1 & GSE33532_alldiff$logFC <= 1)] = "non"
ggplot(GSE33532_alldiff,aes(x=logFC,y=-log10(P.Value),colour=threshold))+ xlab("log2 Fold Change")+ ylab("-log10 P.Value")+ geom_point()+ scale_color_manual(values =c("red","black","green"))+ geom_hline(yintercept=1.30103,linetype=3)+ geom_vline(xintercept=c(-1,1),linetype=3)+ theme(legend.title = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank())
