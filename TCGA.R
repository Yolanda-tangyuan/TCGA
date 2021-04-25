#step1 从临床样本信息中筛选出配对样本
#read.delim用于读取默认分隔符为TAB的文件
#stringsAsFactors= F意味着，“在读入数据是，遇到字符串之后，不将其转换为factors,仍然保留为字符串格式”
rm(list = ls())
getwd()
setwd("/home/LYX/drafts/CRC/TCGA450k/")
getwd()
pd.all<-read.delim("COAD_clinicalMatrix",header = T,stringsAsFactors = F)
#查看一下变量维度，得到行列个数
dim(pd.all)
#总共样本数量是551个，有一些样本来自于同一个人。这就需要看COAD_clinicalMatrix文件中的“_PATIENT”列，如果相同，表示来自于同一个人。不过不被R读入的时候，列名_前面加了X，变成了“X_PATIENT”
#文件总共133列，只有三列是我们后续分析用到的，把这三列提取出来“sampleID”，“X_PATIENT”,"sample_type","pathologic_stage","pathologic_M","pathologic_N","pathologic_T"
head(pd.all)#查看数据结构
colnames(pd.all)#查看列名
pd<-pd.all[,c("sampleID","X_PATIENT","gender","sample_type","pathologic_stage","pathologic_M","pathologic_N","pathologic_T")]
#样本类型包括了，原位肿瘤，复发肿瘤，转移肿瘤和Normal,我们这里只取原位肿瘤和Normal的样本类型做下一类分析
#肿瘤分期包括了“Stage I”，“Stage IA”，“Stage II”，“Stage IIA”，“Stage IIB”，“Stage IIC”，“Stage III”，“Stage IIIA”，“Stage IIIB”，“Stage IIIC”，“Stage IV”，“Stage IVA”，“Stage IVB”，我们这里只选取分期“Stage I”，“Stage II”

table(pd$sample_type)
table(pd$pathologic_stage)

tissue=c("Primary Tumor","Solid Tissue Normal")
pd_tissue<-pd[pd$sample_type %in% tissue,]
stage=c("Stage I","Stage IA","Stage II","Stage IIA","Stage IIB","Stage IIC")
pd_stage<-pd_tissue[pd_tissue$pathologic_stage %in% stage,]
table(pd_stage$sample_type)
table(pd_stage$pathologic_stage)
#然后剔除非配对样本数据#一定要剔除么？
#先定义一个函数，目的是去掉X_PATIENT列值出现过一次所在的行。其中，下面函数可以根据自己需要再具体修改
#来自于：https://blog.csdn.net/dingchenxixi/article/details/50865277
deleteuniquelines<-function(x){#x为输入的数据框
  stand.col<-x$X_PATIENT
  count<-table(stand.col) #table函数可以得到每个上述列每个数所出现的频数
  if(all(count<2)) stop("no repeated records")
  else{
    ind<-sapply(stand.col, function(t) ifelse(count[as.character(t)]>1,TRUE,FALSE))
  }
  return(x[ind,])
}

#应用这个函数
pd_stage_filtered=deleteuniquelines(pd_stage)
#查看一下得到的过滤后的结果，现在样本类型只剩下Primary Tumor Solid和Tissue Normal两种，并且各为49个，且pd_tissue_filtered的行数为98
table(pd_stage_filtered$sample_type)
dim(pd_stage_filtered)

save(pd_stage_filtered,file = "pd_stage_filtered")
#现在，所需要的pair样本提取完毕
#将数据保存,确定文件保存路径
#quote=F是去掉“”，row.names=F不要行名称，需要列名，以\t分割。默认空格分隔
getwd()
setwd("/home/LYX/drafts/CRC/TCGA450k")
write.table(pd_stage_filtered,file = "TCGA_coad_pair_sampleid_from_COAD_clinicalMatrix",quote = F,row.names = F,sep = "\t")

#step2 从芯片数据中提取配对样本的甲基化信息
getwd()

rm(list = ls())
library("data.table")
methy_data<-data.table::fread(file = "/home/LYX/drafts/CRC/TCGA450k/HumanMethylation450.gz",data.table = F)
methy_data[1:4,1:4]
dim(methy_data)
#从上一步得到的pd_stage_filtered配对样本中分别选出normal和tumor样本
load("pd_stage_filtered")
normal=pd_stage_filtered[pd_stage_filtered$sample_type=="Solid Tissue Normal",]
tumor=pd_stage_filtered[pd_stage_filtered$sample_type=="Primary Tumor",]
colnames(methy_data)
normal$sampleID
#从芯片数据中分别挑出normal和tumor样本
methy_normal=methy_data[,colnames(methy_data) %in% normal$sampleID]
head(methy_data)
methy_tumor=methy_data[,colnames(methy_data) %in% tumor$sampleID]
#长这样
methy_normal[1:4,1:4]
methy_tumor[1:4,1:4]
#这里发现从芯片里挑出来的normal样本和tumor样本数量并不一致，应该去掉7个tumor样本
dim(methy_normal)
dim(methy_tumor)
#合并normal和tumor数据
methy_combine=cbind(methy_normal,methy_tumor)
dim(methy_combine)
methy_combine[1:4,1:4]
#需要删除掉只出现过一次的PATIENT。如果注意观察就会发现，这时候的列名的前1-12字段就是PATIENT.
#继续用到上一步的函数
deleteuniquelines<-function(x){ # x为输入的数据框
  m=substr(colnames(x),1,12)
  stand.col<-m
  count<-table(stand.col) #table函数可以得到每个上述列每个数所出现的频数
  if(all(count<2)) stop("no repeated records")
  else{
    ind<-sapply(stand.col, function(t) ifelse(count[as.character(t)]>1,TRUE,FALSE))
  }
  return(x[,ind])
}
methy_paired=deleteuniquelines(methy_combine)
dim(methy_paired)
#此时行名还没更改好
methy_paired[1:4,1:4]
#更改行名
methy_data$sample
rownames(methy_paired)=methy_data$sample
methy_paired[1:4,1:4]
#将数据保存
save(methy_paired,file = "methy_paired.Rdata")
write.csv(methy_paired,file = "TCGA_coad_pair_from_HumanMethylation450.csv")


#step3 ChAMP过滤甲基化芯片数据
rm(list =ls())
library(ChAMP)
library(dplyr)
library(tibble)
#读入之前整理好的肿瘤配对样本的甲基化数据
methy=data.table::fread(file = "TCGA_coad_pair_from_HumanMethylation450.csv",data.table = F)
pd<-read.delim(file = "/home/LYX/drafts/CRC/TCGA450k/TCGA_coad_pair_sampleid_from_COAD_clinicalMatrix",header = T,stringsAsFactors = F)
#长这样：
pd[1:3,1:3]
methy[1:4,1:4]
#先把列名去掉
a=column_to_rownames(methy,"V1")
colnames(a)
pd$sampleID
#从样本信息中选出承兑的甲基化数据样本
pd_info<-pd[pd$sampleID %in% colnames(a),]
#对methy的列进行排序，如果这里不排序，后面分析的组别就会乱
methy_sort<-a[,pd_info$sampleID]
methy_sort[1:4,1:4]
pd_info[1:4,]
beta_value=as.matrix(methy_sort)
#这一步需要一些时间
beta=impute.knn(beta_value)
sum(is.na(beta))
beta=beta$data
head(beta)
beta=beta+0.0000001
pd_1<-as.data.frame(colnames(beta))
#筛选出成对的样本信息
pd_info<-pd[pd$sampleID %in% pd_1$`colnames(beta)`,]
dim(pd_info)
dim(pd)
head(pd_1)
#改一下列名
colnames(pd_1)="sampleID"
pdmerge<-merge(pd_1,pd_info,by="sampleID",all.x=TRUE)
head(pdmerge)
#完成，下面开始过滤
myload=champ.filter(beta = beta,pd=pdmerge)
dim(myload$beta)
dim(myload$pd)
dim(beta)
myload$beta[1:4,1:4]
myload$pd[1:7,1:7]
#保存为CHAMP格式，用于后面的质控
save(myload,file = "TCGA_coad_pair_CHAMPfilter.Rdata")

#step4 ChAMP对读入rdata数据进行QC、标准化
rm(list = ls())
getwd()
setwd("/home/LYX/drafts/CRC/TCGA450k")
library(ChAMP)
library(stringr)
load("TCGA_coad_pair_CHAMPfilter.Rdata")
QC=champ.QC(beta = myload$beta,pheno = myload$pd$sample_type)

#norm
# library("parallel")
# detectCores(logical = F)
# cl <- makeCluster(getOption("cl.cores", 15))
norm_file="TCGA_coad_pair_CHAMPfilter_Norm.Rdata"
if(!file.exists(norm_file)){
  myNorm<-champ.norm(beta = myload$beta,arraytype = "450K",cores = 8)
  save(myNorm,file = norm_file)
}
# stopCluster(cl)

#去掉含有NA值的样本
load("TCGA_coad_pair_CHAMPfilter_Norm.Rdata")
load(norm_file)
num.na<-apply(myNorm,2,function(x)(sum(is.na(x))))
table(num.na)
names(num.na)=colnames(myNorm)
dt=names(num.na[num.na>0])
dn=str_replace(dt,"-01","-11")
keep=setdiff(colnames(myNorm),c(dt,dn))
myNorm=myNorm[,keep]
pd=myload$pd
dim(pd)
pd=pd[pd$sampleID %in% keep,]
dim(pd)
save(pd,file = "pd.Rdata")
save(myNorm,file = "myNorm.Rdata")
load("pd.Rdata")
load("myNorm.Rdata")
#画图

#聚类热图
install.packages("pheatmap")
library("pheatmap")
cg=names(tail(sort(apply(myNorm,1,sd)),1000))
library(pheatmap)
group_list<-pd$sample_type
ac=data.frame(group=group_list)
rownames(ac)=colnames(myNorm)
pdf("TCGA_coad_pair_CHAMPfilter_Norm_88_pheatmap.pdf")
pheatmap(myNorm[cg,],show_colnames=F,show_rownames=F,annotation_col=ac)
dev.off()

#相关性图
pdf("TCGA_coad_pair_CHAMPfilter_Norm_88_correlationpheatmap.pdf")
pheatmap(cor(myNorm[cg,]),show_colnames=F,show_rownames=F,annotation_col=ac)
dev.off()

#SVD画图
#根据myNorm（矫正后的探针Beta）和pd(样本信息)来计算变异程度
#champ.SVD()分析时会把协变量打印在屏幕上，结果是热图，保存为SVDsummary.pdf文件。黑色表示最显著的P值。如果发现技术因素有影响，就需要用ComBat等方法重新标准化数据
#dev.new()
champ.SVD(beta=myNorm,pd=pd)

#Step5差异分析
#用myDMP做差异分析，其中分组名不能用空格，需要替换一下
rm(list = ls())
x<-pd$sample_type
x[which(x=="Primary Tumor")]<-"Tumor"
x[which(x=="Solid Tissue Normal")]<-"Normal"
myDMP<-champ.DMP(beta = myNorm,pheno = x)

head(myDMP$Tumor_to_Normal)
save(myDMP,file = "myDMP.Rdata")
write.csv(myDMP$Tumor_to_Normal,file = "TCGA_coadread_pair_myDMP.csv")
load("myDMP.Rdata")
#标出UP和down
df_DMP<-myDMP$Tumor_to_Normal
df_DMP<-df_DMP[df_DMP$gene!="",]
head(df_DMP)
# logFC_t<-0.45
# P.Value_t<-10^-15
deltaBeta_t<-0.2
P.Value_t<-0.05
df_DMP$change<-ifelse(df_DMP$adj.P.Val<P.Value_t &abs(df_DMP$deltaBeta)>deltaBeta_t,
                      ifelse(df_DMP$deltaBeta>deltaBeta_t,"UP","DOWN"),
                      "NOT")
head(df_DMP)
save(df_DMP,file = "TCGA_coad_pair_CHAMPfilter_Norm_DMP.Rdata")
write.csv(df_DMP,file = "TCGA_coad_pair_myDMP_chang.csv")

#做火山图
getwd()
load("TCGA_coad_pair_CHAMPfilter_Norm_DMP.Rdata")
library(dplyr)
library(ggplot2)
library(tibble)
dat=rownames_to_column(df_DMP)
dat[1:4,1:4]
pdf("Volcano.pdf")
ggplot(data = dat,aes(x=deltaBeta,y=-log10(adj.P.Val)))+
  geom_point(alpha=0.4,size=3.5,aes(color=df_DMP$change))+
  ylab("-log10(Pvalue)")+
  scale_color_manual(values = c("blue","grey","red"))+
  geom_vline(xintercept = c(-deltaBeta_t,deltaBeta_t),lty=4,col="black",lwd=0.8)+
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8)+
  theme_bw()
dev.new()
dev.off()

#差异基因热图
load("TCGA_coad_pair_CHAMPfilter_Norm.Rdata")
num.na<-apply(myNorm,2,function(x)(sum(is.na(x))))
names(num.na)=colnames(myNorm)
dt=names(num.na[num.na>0])
library(stringr)
dn=str_replace(dt,"-01","-11")
keep=setdiff(colnames(myNorm),c(dt,dn))
myNorm=myNorm[,keep]
load("TCGA_coad_pair_CHAMPfilter.Rdata")
pd=myload$pd
pd=pd[pd$sampleID %in% keep,]
x<-pd$sample_type
x[which(x=="Primary Tumor")]<-"Tumor"
x[which(x=="Solid Tissue Normal")]<-"Normal"
pd$sample_type=x
group_list<-pd$sample_type
cg<-rownames(df_DMP[df_DMP$change!="NOT",])
plot_matrix<-myNorm[cg,]
annotation_col<-data.frame(Sample=pd$sample_type)
rownames(annotation_col)<-colnames(plot_matrix)
ann_colors=list(Sample=c(Normal="#4DAF4A",Tumor="#E41A1C"))
library(pheatmap)
pdf("DEG.pheatmap.pdf")
pheatmap(plot_matrix,show_colnames = T,annotation_col = annotation_col,
         border_color = NA,color = colorRampPalette(colors = c("white","navy"))(50),annotation_colors = ann_colors,show_rownames = F,cellwidth = 15, cellheight = 12)
dev.off()
dev.new()

#甲基化差异基因的整理
#将Tumor中高甲基化基因和低甲基化基因挑选出来，保存，用于后续分析
downgene<-df_DMP[df_DMP$change=="DOWN",]$gene
upgene<-df_DMP[df_DMP$change=="UP",]$gene
write.table(upgene,file = "Tumor_UP_gene.xls",quote = F,row.names = F,sep = "\t",col.names = F)
write.table(downgene,file = "Tumor_DOWN_gene.xls",quote = F,row.names = F,sep = "\t",col.names = F)
write.csv(upgene,file = "Tumor_UP_gene.csv")
write.csv(downgene,file = "Tumor_DOWN_gene.csv")
