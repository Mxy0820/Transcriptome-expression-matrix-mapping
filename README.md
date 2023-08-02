# Transcriptome-expression-matrix-mapping

#差异表达基因的富集分析

library(clusterProfiler)
library(tidyverse)
library(AnnotationForge)
library(dplyr)
library(readr)

gene_exp <- read_delim("D:/R/R_practice/1.transcriptome/data/rsem_matrix.TMM.EXPR.matrix.txt", 
                       delim = "\t", escape_double = FALSE, 
                       trim_ws = TRUE)


de_result <- read_delim("D:/R/R_practice/1.transcriptome/data/DE_results.txt", 
                        delim = "\t", escape_double = FALSE, 
                        trim_ws = TRUE)


temp <- select(de_result, id, log2FoldChange, pvalue, padj) %>%
    mutate(direction = if_else(abs(log2FoldChange) < 1 | padj > 0.05, 'ns',
                               if_else(log2FoldChange >= 1, 'up', 'down')))

de_result <- left_join(temp, gene_exp, by='id')

group_by(de_result,direction)%>%
    summarise(count=n())

#准备输入文件
##差异表达表
gene <- filter(de_result,
               abs(log2FoldChange) >= 1 & padj < 0.05) %>%
    pull(id)

##基因注释信息表
emapper <- read_delim("D:/R/R_practice/1.transcriptome/data/out-2.emapper.annotations",
                      "\t", escape_double = FALSE, col_names = FALSE, comment = "#",
                      trim_ws = TRUE) %>%
    dplyr::select(GID = X1,
                  Gene_Symbol = X9,
                  GO = X10,
                  KO = X12,
                  Pathway = X13,
                  OG = X7,
                  Gene_Name = X8)

##基因差异向量
geneList <- de_result$log2FoldChange
names(geneList) <- de_result$id
geneList <- sort(geneList, decreasing = T)


##构建OrgDB
###提取基因信息
gene_info_1 <- dplyr::select(emapper, GID, Gene_Name) %>%
    dplyr::filter(!is.na(Gene_Name))

gene_info <- dplyr::filter(gene_info_1, Gene_Name != "-")


###提取GO信息
gene2go_1 <- dplyr::select(emapper,GID, GO)%>%
    separate_rows(GO, sep = ',', convert = F)%>%
    filter(!is.na(GO))%>%
    mutate(EVIDENCE = 'IEA')#把基因与go的关系表示出来

gene2go <- dplyr::filter(gene2go_1, GO != "-")

###构建OrgDB
AnnotationForge::makeOrgPackage(gene_info=gene_info,
                                go=gene2go,
                                maintainer='Xinyi_Ma <Xinyi_Ma@genek.tv>',
                                author = 'Tengzhuo_Ma',
                                outputDir = "./",
                                tax_id = 0000,
                                genus = 'A',
                                species = 'chi',
                                goTable = 'go',
                                version = '1.0')

###打包
pkgbuild::build('./org.Achi.eg.db', dest_path = "./Enrichment-clusterprofiler")


###创建一个文件夹
dir.create('R_Library', recursive = T)#在当前目录下创建一个文件夹

###将包安装在该文件夹下
install.packages('Enrichment-clusterprofiler/org.Achi.eg.db_1.0.tar.gz',
                 repos =NULL,#从本地安装
                 lib = 'R_Library')#安装在文件夹





#绘图
##准备输入文件
top_de_149 <- filter(de_result,
                     abs(log2FoldChange) >= 1 & padj < 0.05) %>%
    dplyr::select(-log2FoldChange, -pvalue, -padj, -direction) %>%
    column_to_rownames(var = 'id')

##Fig.4a
library(EnhancedVolcano)

ed <- de_result
ed
ed$log2FoldChange[(ed$log2FoldChange) < -5] <- -5
ed$log2FoldChange[(ed$log2FoldChange) > 5] <- 5

ed$padj[(ed$padj) < 1e-10] <- 1e-10

vals <- c('Cchi04059')

group_color_1 <-ifelse(
    ed$log2FoldChange<(-1)&ed$padj<0.05,'green',
    ifelse(ed$log2FoldChange>(1)&ed$padj<0.05,'red','#b5b5b5'))

group_color_1[is.na(group_color_1)]<-'#b5b5b5'

names(group_color_1)[group_color_1=='red']<-'Up'

names(group_color_1)[group_color_1=='#b5b5b5']<-'NS'

names(group_color_1)[group_color_1=='green']<-'Down'

EnhancedVolcano(
    ed,
    lab = ed$id,
    x = 'log2FoldChange',
    y = 'padj',
    xlim=c(-5, 5),
    ylim=c(0,11),
    FCcutoff = 1,
    pCutoff = 0.05,
    selectLab=c(vals),#使用selectLab参数选定所关注的标签
    xlab=bquote(~Log[2]~'fold change'),#将内容传递给xlab
    labCol='black',#标签颜色
    labFace='bold',#标签字体
    boxedLabels=TRUE,#是否在框中绘制标签
    drawConnectors=TRUE,#是否通过连线将标签连接到对应的点上
    widthConnectors=0.4,#连线的宽度
    endsConnectors= "last",#连线绘制箭头的方向，可选first、both、last
    colConnectors='black',#连线的颜色
    colCustom=group_color_1,#用group覆盖默认配色方案
    caption=bquote(~Log[2]~ "fold change cutoff,1;padj cutoff,0.05"),#注释说明
    cutoffLineType='longdash',#阈值线类型，可选“blank”、“solid”、“dashed”、“dotted”、“dotdash”、“longdash”和“twodash”
    cutoffLineCol='pink',#阈值线颜色
    cutoffLineWidth=0.1,#阈值线粗细
    pointSize=c(ifelse(de_result$id %in% vals,5,3))#通过ifelse条件判断函数指定散点大小
)


##Supplymentary-Fig.S3_热图绘制——ggplot2
library (ggplot2)
library (reshape2)#数据转换
require(scales)#数据缩放
library(ggtree)#聚类
library(aplot)#拼图



M <- log2(top_de_149 +1)

B <- log2(top_de_149 +1)
B <- t(B)
B %>% scale() %>% as.data.frame()#缩放数据
C <- B %>% scale(center = T) %>% as.data.frame()#以行缩放
C <- C%>% mutate(B=row.names(.)) %>% melt()#转化为ggplot画图需要的长列表


group <- colnames(M) %>% as.data.frame() %>% 
    mutate(group=c(rep("D1_star",1),rep("D2_star",1),rep("D3_star",1),rep("N1_star",1),rep("N2_star",1),rep("N3_star",1))) %>%
    #mutate(group=c(rep("D1_star","D2_star","D3_star",1),rep("N1_star","N2_star","N3_star",1))) %>%
    mutate(p="group") %>%
    ggplot(aes(.,y=p,fill=group))+
    geom_tile() + 
    scale_y_discrete(position="right") +
    theme_minimal()+xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_blank())+
    labs(fill = "Group")


###然后用ggtree做聚类。

p <- M %>% scale(center = T) %>% as.data.frame()
phr <- hclust(dist(p)) %>% 
    ggtree(layout="rectangular",branch.length="none")
phc <- hclust(dist(t(p))) %>% 
    ggtree() + layout_dendrogram()

###画热图并将以上信息添加进去：

pH1<-ggplot(C,aes(x=B,y=variable,fill=value)) #热图绘制
pH2 <- p1+geom_raster()+scale_fill_gradient2(low="green", high="red", mid="white")+
    geom_tile()+theme_minimal()+
    theme(axis.text.x =element_text(angle =90,hjust =0.1,vjust = 0.1))+
    xlab(NULL) + ylab(NULL) +
    #theme(axis.text.y = element_text(size = 6))+
    #theme(axis.text.x = element_text(size = 0))+
    theme(axis.text.x = element_blank())+
    theme(axis.text.y = element_blank())
pH2 %>%
    insert_top(group, height = .03)%>%
    #insert_left(Pathway, width = .05)%>%
    insert_left(phr,width=.2)


##GO富集
library(enrichplot)
library(ggnewscale)
library(ggplot2)
library(org.Achi.eg.db, lib = 'D:/R/R_practice/1.transcriptome/R_Library')#从某文件夹下加载包，即该软件包没有装在公共目录下（平时软件默认安装位置），防止构建的db文件被替换掉

de_ego_0612 <- enrichGO(gene = gene,
                        OrgDb = org.Achi.eg.db,
                        keyType = 'GID',
                        ont = 'all',
                        pvalueCutoff = 0.5,
                        qvalueCutoff = 0.5)
###导出文件到指定位置  
write.table(de_ego_0612,
            file = "D:/R/R_practice/1.transcriptome/de_ego_0612.txt",
            sep = "\t",
            quote = F,
            row.names = F) 

###由于有些图需要用到Enrichfaxtor绘图，自行计算
eGo = read.table("de_ego_0612.txt",
                 header = TRUE,
                 sep = "\t",
                 quote = "") 
####劈分GeneRatio为两列    
eGo <- separate(data = eGo,
                col = GeneRatio,
                into = c("GR1", "GR2"),
                sep = "/")   
####劈分BgRatio为两列
eGo <- separate(data = eGo,
                col = BgRatio,
                into = c("BR1", "BR2"),
                sep = "/")  
####计算EnrichFactor
eGo <- mutate(eGo,
              enrichment_factor = (as.numeric(GR1)/as.numeric(GR2))/(as.numeric(BR1)/as.numeric(BR2)))

##Supplymentary-Fig.S4
eGoBP <- eGo %>%
    filter(ONTOLOGY == "BP") %>%
    filter(row_number() >= 1,
           row_number() <= 10)

eGoCC <- eGo %>%
    filter(ONTOLOGY == "CC") %>%
    filter(row_number() >= 1,
           row_number() <= 10)

eGoMF <- eGo %>%
    filter(ONTOLOGY == "MF") %>%
    filter(row_number() >= 1,
           row_number() <= 10)
eGo10 <- rbind(eGoBP, eGoCC, eGoMF)

library(ggplot2)
pE <- ggplot(eGo10,
            aes(enrichment_factor,
                fct_reorder(factor(Description),enrichment_factor))) + 
    geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) + 
    scale_color_gradient(low="green",high ="red") + 
    labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology", x="Enrichment Factor",y="GO term",title="GO enrichment") + 
    theme(axis.text.y = element_text(angle = -45)) + 
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    theme_bw() 
pE + facet_wrap( ~ ONTOLOGY,ncol= 1,scale='free')


##Fig.4b
eGoBP_20 <- eGo %>%
    filter(ONTOLOGY == "BP") %>%
    filter(row_number() >= 1,
           row_number() <= 20)

pE_20 <- ggplot(eGoBP_20,
             aes(enrichment_factor,
                 fct_reorder(factor(Description),enrichment_factor))) + 
    geom_point(aes(size=Count,color=-1*log10(pvalue),shape=ONTOLOGY)) + 
    scale_color_gradient(low="green",high ="red") + 
    labs(color=expression(-log[10](p_value)),size="Count", shape="Ontology", x="Enrichment Factor",y="GO term",title="GO enrichment") + 
    theme(axis.text.y = element_text(angle = -45)) + #size = 3, family = "myFont", color = "black", face = "bold", vjust = 0.5, hjust = 0.5,
    #scale_y_continuous(limits = c(-5,10))+
    theme(axis.title.x = element_blank())+
    theme(axis.title.y = element_blank())+
    theme_bw() 
pE_20


###通路树状图绘制
#de_ego_tree <- enrichGO(gene = gene,
#                        OrgDb = org.Achi.eg.db,
#                        keyType = 'GID',
#                        ont = 'BP',
#                        pvalueCutoff = 0.5,
#                        qvalueCutoff = 0.5)

write.table(de_ego_tree,
            file = "D:/R/R_practice/1.transcriptome/de_ego_tree.txt",
            sep = "\t",
            quote = F,
            row.names = F) 

plotGOgraph(de_ego_tree,
            firstSigNodes = 10,
            useInfo = "all",
            sigForAll = TRUE,
            useFullNames = TRUE)


##Fig.5a
y <- c("response to hormone","lipid catabolic process","lipid homeostasis","positive regulation of lipid transport")
p_cnteplot <- cnetplot(de_ego,
         foldChange = geneList,
         showCategory = y,
         node_lable = 'all',
         circular = TRUE,
         colorEdge = TRUE)

p_cnteplot + scale_color_gradientn(colours = c("green","white", "red"))

##Fig.5b_热图绘制——ggplot2

fpath_r_genes <- `4pathway.r.genes`
colnames(fpath_r_genes)[1] <- c("id")
temp_20230613 <- select(de_result, id, D1_star, D2_star, D3_star, N1_star, N2_star, N3_star)

fpath_r_genes <- left_join(fpath_r_genes, temp_20230613, by='id')%>%
    column_to_rownames(var = 'id')
#fpath_r_genes <- merge(fpath_r_genes, top_de_149, by = "row.names", all = F)

fpath_r_genes <- log2(fpath_r_genes +1)
fpath_r_genes_B <- t(fpath_r_genes)
fpath_r_genes_B %>% scale() %>% as.data.frame()#缩放数据
fpath_r_genes_C <- fpath_r_genes_B %>% scale(center = T) %>% as.data.frame()#以行缩放
fpath_r_genes_C <- fpath_r_genes_C%>% mutate(fpath_r_genes_B=row.names(.)) %>% melt()#转化为ggplot画图需要的长列表


group_4 <- colnames(fpath_r_genes) %>% as.data.frame() %>% 
    mutate(group = c(rep("D1_star",1),rep("D2_star",1),rep("D3_star",1),rep("N1_star",1),rep("N2_star",1),rep("N3_star",1))) %>%
    mutate(p="group") %>%
    ggplot(aes(.,y=p,fill=group))+
    geom_tile() + 
    scale_y_discrete(position="right") +
    theme_minimal()+xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_blank())+
    labs(fill = "Group")

#然后用ggtree做聚类。

p_4 <- fpath_r_genes %>% scale(center = T) %>% as.data.frame()
phr_4 <- hclust(dist(p_4)) %>% 
    ggtree(layout="rectangular",branch.length="none")
phc_4 <- hclust(dist(t(p_4))) %>% 
    ggtree() + layout_dendrogram()

#画热图并将以上信息添加进去：
pf1<-ggplot(fpath_r_genes_C,aes(x=fpath_r_genes_B,y=variable,fill=value)) #热图绘制
pf2 <- pf1+geom_raster()+scale_fill_gradient2(low="green", high="red", mid="white")+
    geom_tile()+theme_minimal()+
    theme(axis.text.x =element_text(angle =90,hjust =0.1,vjust = 0.1))+
    xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_text(size = 10))
pf2 %>%
    insert_top(group_4, height = .03)%>%
    #insert_left(Pathway, width = .05)%>%
    insert_left(phr_4,width=.2)


##Fig.6a
Fig_6a <- qPCR_ts
colnames(Fig_6a) <- Fig_6a[1,]
Fig_6a <- Fig_6a[-(1),]
rownames(Fig_6a) <- Fig_6a[,1]
Fig_6a <- Fig_6a[,-1]

Fig_6a <- as.data.frame(lapply(Fig_6a,
                               function(x) as.numeric(as.character(x))))
row.names(Fig_6a) <- colnames(Fig_6a_B)

Fig_6a_B <- t(Fig_6a)
Fig_6a_B %>% scale() %>% as.data.frame()#缩放数据
Fig_6a_C <- Fig_6a_B %>% scale(center = T) %>% as.data.frame()#以行缩放
Fig_6a_C <- Fig_6a_C%>% mutate(Fig_6a_B=row.names(.)) %>% melt()#转化为ggplot画图需要的长列表


group_Fig_6a <- colnames(Fig_6a) %>% as.data.frame() %>% 
    mutate(group = c(rep("D1",1),rep("D2",1),rep("D3",1),rep("N1",1),rep("N2",1),rep("N3",1))) %>%
    mutate(p="group") %>%
    ggplot(aes(.,y=p,fill=group))+
    geom_tile() + 
    scale_y_discrete(position="right") +
    theme_minimal()+xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_blank())+
    labs(fill = "Group")

Subfamily_Fig_6a_A <- rownames(Fig_6a) %>% as.data.frame() %>%
  mutate(group=rep(c("AcHsp90","AcHsp70","AcHsp60","AcsHsp"),
                   times=c(2,14,4,5))) %>%
  mutate(p="")%>%
  ggplot(aes(p,.,fill=group))+
  geom_tile() + 
  scale_y_discrete(position="right") +
  theme_minimal()+xlab(NULL) + ylab(NULL) +
  theme(axis.text.y = element_blank(),
        axis.text.x =element_text(
          angle =90,hjust =0.5,vjust = 0.5))+
  labs(fill = "Subfamily")

#用ggtree做聚类。

Fig_6a_p <- Fig_6a %>% scale(center = T) %>% as.data.frame()
Fig_6a_phr <- hclust(dist(Fig_6a_p)) %>% 
    ggtree(layout="rectangular",branch.length="none")
Fig_6a_phc <- hclust(dist(t(Fig_6a_p))) %>% 
    ggtree() + layout_dendrogram()

#画热图并将以上信息添加进去：
Fig_6a_p1 <- ggplot(Fig_6a_C,aes(x=Fig_6a_B,y=variable,fill=value)) #热图绘制

Fig_6a_p2 <- Fig_6a_p1+geom_raster()+scale_fill_gradient2(low="green", high="red", mid="white")+
    geom_tile()+theme_minimal()+
    theme(axis.text.x =element_text(angle =90,hjust =0.1,vjust = 0.1))+
    xlab(NULL) + ylab(NULL) +
    theme(axis.text.x = element_blank()) +
    theme(axis.text.y = element_text(size = 10))
Fig_6a_p2 %>%
    insert_top(group_Fig_6a, height = .03)%>%
    insert_left(Subfamily_Fig_6a_A, width = .2)%>%
    insert_left(Fig_6a_phr,width=.2)

