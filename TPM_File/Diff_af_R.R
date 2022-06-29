library(ggplot2)
library(limma)
library(pheatmap)
library(ggsci)
library(dplyr)

rt=read.table("matrix.txt",sep="\t",header=T,check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
rt=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
rt=avereps(rt)

#seq生成步长为2的等差序列，将正常组按列提取到前面，肿瘤组提取到后面
rt=rt[,
      c(
        seq(1,115,2) ,seq(2,116,2)
        )
      ]

#判断原始数据是否去了log
max(rt)
if(max(rt)>30) rt=log2(rt+1)     #rt最大值大于30则取log

#使用normalizeBetweenArrays进行矫正，矫正后赋值为rt1
rt1=normalizeBetweenArrays(as.matrix(rt))

#未标准化
cols=rainbow(ncol(rt)) ###针对24个样本，设置颜色，整体呈现彩虹色
par(cex = 0.7)
if(ncol(rt)>40) par(cex = 0.5)   ###设置字体大小
#pdf(file = "raw.pdf",width=5,height = 4)
boxplot(rt,las=2,col =cols ) ###绘图
#dev.off()

#标准化
cols=rainbow(ncol(rt1)) ###针对24个样本，设置颜色，整体呈现彩虹色
par(cex = 0.7)
if(ncol(rt1)>40) par(cex = 0.5)   ###设置字体大小
#pdf(file = "nor.pdf",width=5,height = 4.5)
boxplot(rt1,las=2,col =cols ) ###绘图
#dev.off()

#保存标准化后结果
rt2=rbind(ID=colnames(rt1),rt1)
write.table(rt2,file="norexp.txt",sep="\t",quote=F,col.names = F)

####差异分析
data=rt1
#data=rt

#控制组的数量，因为前面已经把正常组提到了最前面几列（一定要根据自己的数据集情况操作），所以直接提就行了
afcon=58
conData=data[,as.vector(colnames(data)[1:afcon])]
aftreat=afcon+1
treatData=data[,as.vector(colnames(data)[aftreat:ncol(data)])]
rt=cbind(conData,treatData)
conNum=ncol(conData)
treatNum=ncol(treatData)

#limma差异标准流程
Type=c(rep("con",conNum),rep("treat",treatNum))
design <- model.matrix(~0+factor(Type))
colnames(design) <- c("con","treat")
fit <- lmFit(rt,design)
cont.matrix<-makeContrasts(treat-con,levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

Diff=topTable(fit2,adjust='fdr',number=length(rownames(data)))
#保存所有基因的差异结果
DIFFOUT=rbind(id=colnames(Diff),Diff)
write.table(DIFFOUT,file="DIFF_all.xls",sep="\t",quote=F,col.names=F)

#热图展示差异最大的前50个基因
Diff=Diff[order(as.numeric(as.vector(Diff$logFC))),]
diffGene=as.vector(rownames(Diff))
diffLength=length(diffGene)
afGene=c()
if(diffLength>(100)){
  afGene=diffGene[c(1:50,(diffLength-50+1):diffLength)]
}else{
  afGene=diffGene
}
afExp=rt[afGene,]
#分组标签
Type=c(rep("N",conNum),rep("T",treatNum))
names(Type)=colnames(rt)
Type=as.data.frame(Type)
#分组标签的注释颜色
anncolor=list(Type=c(T=pal_npg()(1),N=pal_npg()(2)[2]))

##pdf(file="DIFF_heatmap.pdf",height=7,width=8)
pheatmap(afExp,                                                                      #热图数据
         annotation=Type,                                                            #分组
         color = colorRampPalette(c(pal_npg()(2)[2],"white", pal_npg()(1)))(50),     #热图颜色
         cluster_cols =F,                                                           #不添加列聚类树
         show_colnames = F,                                                         #展示列名
         scale="row", 
         fontsize = 8,
         fontsize_row=6,
         fontsize_col=8,
         annotation_colors=anncolor
)
#dev.off()

#火山图差异标准设置
adjP=0.05
aflogFC=1
Significant=ifelse((Diff$P.Value<adjP & abs(Diff$logFC)>aflogFC), ifelse(Diff$logFC>aflogFC,"Up","Down"), "Not")
#开始绘制
p = ggplot(Diff, aes(logFC, -log10(P.Value)))+
  geom_point(aes(col=Significant),size=3)+
  scale_color_manual(values=c(pal_npg()(2)[2], "#838B8B", pal_npg()(1)))+
  labs(title = " ")+
  theme(plot.title = element_text(size = 16, hjust = 0.5, face = "bold"))+
  geom_hline(aes(yintercept=-log10(adjP)), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=aflogFC), colour="gray", linetype="twodash",size=1)+
  geom_vline(aes(xintercept=-aflogFC), colour="gray", linetype="twodash",size=1)
#添加标记，按照
point.Pvalue=0.0001
point.logFc=3.8
#继续绘制
Diff$symbol=rownames(Diff)
#pdf("DIFF_vol.pdf",width=5.5,height=5)
p=p+theme_bw()
for_label <- Diff %>% 
  filter(abs(logFC) >point.logFc & P.Value< point.Pvalue )
p+geom_point(size = 2, shape = 1, data = for_label) +
  ggrepel::geom_label_repel(
    aes(label = symbol),
    data = for_label,
    color="black",
    label.size =0.1
  )
#dev.off()