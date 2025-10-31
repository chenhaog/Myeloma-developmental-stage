pkgs<- c("futile.logger", "configr", "stringr", "ggpubr", "ggthemes","glue", "ggsci",
"patchwork", "tidyverse", "dplyr",'Seurat','GEOquery','SeuratData','ggplot2','dplyr',
'patchwork','AnnoProbe','rjags','infercnv','factoextra','FactoMineR','scater',
'scRNAseq','ROCR','M3Drop','monocle3','dbscan','reticulate','devtools','flexclust',
'mcclust','stringr','clusterProfiler','org.Hs.eg.db','phylogram','gridExtra','grid',
'dendextend','miscTools','SeuratDisk','SeuratWrappers')

for(pkg in pkgs){
  suppressPackageStartupMessages(library(pkg, character.only = T))
}
project<- 'project'
dataset<- 'dataset'
species<- 'species'
workdir<- glue::glue('C:\/Projects/{project}/')
setwd(workdir)
sce<- LoadH5Seurat('F://HD_plasma/sce.h5seurat')
c<- sce@meta.data
e<- read.csv('./test1.csv',row.names = 1,header = T)
value<- read.csv('./test1_proba.csv')
value<- value[,-1]
rownames(value)<- rownames(e)
e<- cbind(e,value)
colnames(e)<- c('prediction','blood_cycling_plasmacell','blood_end_plasmacell',
                'blood_highCXCR4_plasmacell','blood_lowCXCR4_plasmacell',
                'blood_plasmablast','tonsil_cycling_plasmacell','tonsil_highCD74_plasmacell',
                'tonsil_lowCD74_plasmacell')
e$celltype<- sce@meta.data[rownames(e),'celltype']
blood_cycling_plasmacell<- e
blood_cycling_plasmacell$prediction[which(blood_cycling_plasmacell$prediction != 'blood_cycling_plasmacell')]<- 'No'
blood_cycling_plasmacell$celltype[which(blood_cycling_plasmacell$celltype != 'blood_cycling_plasmacell')]<- 'No'
# install.packages('pROC')
library(pROC)
dfroc1<- roc(blood_cycling_plasmacell$celltype, blood_cycling_plasmacell$blood_cycling_plasmacell)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#216b9a")+
  labs(x= '',y= '',title= 'blood_cycling_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.35,y=0.75,label='AUC = 0.9997061',size=5)

blood_end_plasmacell<- e
blood_end_plasmacell$prediction[which(blood_end_plasmacell$prediction != 'blood_end_plasmacell')]<- 'No'
blood_end_plasmacell$celltype[which(blood_end_plasmacell$celltype != 'blood_end_plasmacell')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(blood_end_plasmacell$celltype, blood_end_plasmacell$blood_end_plasmacell)
dfroc1[["auc"]]
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#dc7f41")+
  labs(x= '',y= '',title= 'blood_end_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.35,y=0.75,label='AUC = 0.997721',size=5)

blood_highCXCR4_plasmacell<- e
blood_highCXCR4_plasmacell$prediction[which(blood_highCXCR4_plasmacell$prediction != 'blood_highCXCR4_plasmacell')]<- 'No'
blood_highCXCR4_plasmacell$celltype[which(blood_highCXCR4_plasmacell$celltype != 'blood_highCXCR4_plasmacell')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(blood_highCXCR4_plasmacell$celltype, blood_highCXCR4_plasmacell$blood_highCXCR4_plasmacell)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
dfroc1[["auc"]]
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#478f46")+
  labs(x= '',y= '',title= 'blood_highCXCR4_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.35,y=0.75,label='AUC = 0.9952348',size=5)

blood_lowCXCR4_plasmacell<- e
blood_lowCXCR4_plasmacell$prediction[which(blood_lowCXCR4_plasmacell$prediction != 'blood_lowCXCR4_plasmacell')]<- 'No'
blood_lowCXCR4_plasmacell$celltype[which(blood_lowCXCR4_plasmacell$celltype != 'blood_lowCXCR4_plasmacell')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(blood_lowCXCR4_plasmacell$celltype, blood_lowCXCR4_plasmacell$blood_lowCXCR4_plasmacell)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#b53833")+
  labs(x= '',y= '',title= 'blood_lowCXCR4_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.35,y=0.75,label='AUC = 0.9946305',size=5)

blood_plasmablast<- e
blood_plasmablast$prediction[which(blood_plasmablast$prediction != 'blood_plasmablast')]<- 'No'
blood_plasmablast$celltype[which(blood_plasmablast$celltype != 'blood_plasmablast')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(blood_plasmablast$celltype, blood_plasmablast$blood_plasmablast)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#7f619b")+
  labs(x= '',y= '',title= 'blood_plasmablast')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.38,y=0.75,label='AUC = 0.9981277',size=5)

tonsil_cycling_plasmacell<- e
tonsil_cycling_plasmacell$prediction[which(tonsil_cycling_plasmacell$prediction != 'tonsil_cycling_plasmacell')]<- 'No'
tonsil_cycling_plasmacell$celltype[which(tonsil_cycling_plasmacell$celltype != 'tonsil_cycling_plasmacell')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(tonsil_cycling_plasmacell$celltype, tonsil_cycling_plasmacell$tonsil_cycling_plasmacell)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#81564d")+
  labs(x= '',y= '',title= 'tonsil_cycling_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.4,y=0.75,label='AUC = 0.9867113',size=5)

tonsil_highCD74_plasmacell<- e
tonsil_highCD74_plasmacell$prediction[which(tonsil_highCD74_plasmacell$prediction != 'tonsil_highCD74_plasmacell')]<- 'No'
tonsil_highCD74_plasmacell$celltype[which(tonsil_highCD74_plasmacell$celltype != 'tonsil_highCD74_plasmacell')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(tonsil_highCD74_plasmacell$celltype, tonsil_highCD74_plasmacell$tonsil_highCD74_plasmacell)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#bc76a1")+
  labs(x= '',y= '',title= 'tonsil_highCD74_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.35,y=0.75,label='AUC = 0.9973214',size=5)

tonsil_lowCD74_plasmacell<- e
tonsil_lowCD74_plasmacell$prediction[which(tonsil_lowCD74_plasmacell$prediction != 'tonsil_lowCD74_plasmacell')]<- 'No'
tonsil_lowCD74_plasmacell$celltype[which(tonsil_lowCD74_plasmacell$celltype != 'tonsil_lowCD74_plasmacell')]<- 'No'
# install.packages('pROC')
# library(pROC)
dfroc1<- roc(tonsil_lowCD74_plasmacell$celltype, tonsil_lowCD74_plasmacell$tonsil_lowCD74_plasmacell)
plotdata<- data.frame(x=c(1-dfroc1[["specificities"]]),y=dfroc1[["sensitivities"]])
ggplot(plotdata)+geom_path(aes(x=x,y=y),size=1,color= "#7d7d7f")+
  labs(x= '',y= '',title= 'tonsil_lowCD74_plasmacell')+
  geom_abline(intercept = 0, slope = 1)+
  theme_bw()+
  theme(text = element_text(size = 12))+Seurat::NoLegend()+
  annotate('text',x=0.35,y=0.75,label='AUC = 0.9950244',size=5)
