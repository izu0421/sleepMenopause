# ============================================================================
# 03 - Hormone Cell Atlas: data-overview figures
# ----------------------------------------------------------------------------
# Project : sleepMenopause - androgen decline and perimenopausal sleep
# Purpose : Summarise the single-cell Hormone Cell Atlas subset (ovary, adrenal,
#           breast) - donor counts, cell-type counts and tissue composition -
#           used to contextualise the steroidogenic-enzyme and androgen-receptor
#           analyses (manuscript Figure 6 / data overview).
# Input   : 3tissue_obs_ovary_adrenal_breast.csv, tissue_colors.csv
#           (NOT included - see "Data availability" in README.md)
# Output  : Overview bar plots (donor number and cell-type number per tissue).
# Depends : data.table, ggplot2, dplyr, stringr, scales, cowplot
# Note    : Paths below are absolute (Windows) from the original analysis
#           environment and must be updated to your local setup before running.
# ============================================================================

library(data.table)
options(datatable.fread.datatable=FALSE)
setwd('C:/FLJ/Yizhou_testerone/data_overview')
anno<-fread("C:/FLJ/Yizhou_testerone/data_overview/3tissue_obs_ovary_adrenal_breast.csv")
head(anno,2)
colnames(anno)
table(anno$tissue)

library(ggplot2)
library(dplyr)
library(stringr)
library(scales)
library(cowplot)

tissue_color<-read.csv("C:/FLJ/hormone/Color/tissue_colors.csv")
head(tissue_color)
tissue_color<-tissue_color[tissue_color$Tissue%in% anno$tissue,]
colnames(tissue_color)<-c('tissue','tissue_col')
tissue_color$tissue<-factor(tissue_color$tissue,levels = c("AdrenalGland","Ovary","Breast"))
## donor number
tmp<-as.data.frame.matrix(table(anno$tissue,anno$donor_id))
tmp[tmp>0]<-1
study<-data.frame(tissue=rownames(tmp),Number=rowSums(tmp))
study$tissue<-factor(study$tissue,levels = levels(tissue_color$tissue))#tissue_color$tissue
head(study)

p2<-ggplot(study,aes(x=tissue,y = log(Number,base = 10),fill=tissue))+geom_bar(stat = 'identity')+ #
  coord_flip()+theme_bw()+  geom_text(aes(label = Number), 
                                      hjust = -0.1) +
  theme(legend.position = 'None',plot.title  = element_text(hjust = 0.5),panel.grid = element_blank())+
  ggtitle('Donor Number')+scale_fill_manual(values = tissue_color$tissue_col)+
  scale_y_continuous(breaks = pretty_breaks())
p2

## celltype
tmp<-as.data.frame.matrix(table(anno$tissue,anno$celltype_level1))
tmp[tmp>0]<-1
study<-data.frame(tissue=rownames(tmp),Number=rowSums(tmp))
study$tissue<-factor(study$tissue,levels = levels(tissue_color$tissue))

p3<-ggplot(study,aes(x=tissue,y = Number,fill=tissue))+geom_bar(stat = 'identity')+
  coord_flip()+theme_bw()+geom_text(aes(label = Number), 
                                    hjust = -0.1) +
  theme(legend.position = 'None',plot.title  = element_text(hjust = 0.5),panel.grid = element_blank())+
  ggtitle('Celltype Number')+scale_fill_manual(values = tissue_color$tissue_col)+
  scale_y_continuous(breaks = pretty_breaks())
p3


## cell number
tmp<-as.data.frame(table(anno$tissue))
colnames(tmp)<-c('tissue','Number')
study<-tmp
#study$tissue<- str_to_title(study$tissue)
study$tissue<-factor(study$tissue,levels = levels(tissue_color$tissue))

p4<-ggplot(study,aes(x=tissue,y = log(Number,base = 10),fill=tissue))+geom_bar(stat = 'identity')+
  coord_flip()+theme_bw()+geom_text(aes(label = Number), 
                                    hjust = 0.5) +theme(legend.position = 'None',panel.grid = element_blank(),plot.title  = element_text(hjust = 0.5))+
  ggtitle('Cell Number')+scale_fill_manual(values = tissue_color$tissue_col)+
  scale_y_continuous(breaks = pretty_breaks())
p4



plot_grid(p2,p3,p4,ncol = 1)
ggsave('Data_condition_log.pdf',height =5,width = 6.5)




## age
anno['bar']=rep(0)
anno$menopause<-factor(anno$menopause,levels = c('>=50','<50'))
library(rcartocolor)
my_colors = carto_pal(7, "ArmyRose")
my_colors<-colorRampPalette(my_colors)(length(unique(anno$menopause)))
my_colors<-rev(my_colors)

age_list=list()
tissues=levels(tissue_color$tissue)

for (tissue in tissues) {
  anno1 <- anno[anno$tissue == tissue, ]
  
  p <- ggplot(anno1, aes(x = bar)) +
    geom_bar(aes(fill = menopause),
             position = 'fill',
             stat = 'count') +
    coord_flip() +
    theme_bw() +
    theme(
      panel.grid = element_blank(),
      plot.title = element_text(hjust = 0.5),
      text = element_text(size = 15),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()
      ) +
    ggtitle(paste('age_group -', tissue)) +
    scale_fill_manual(name = 'age_group', values = my_colors)
  
  age_list[[tissue]] <- p
}
plot_grid(plotlist = age_list,ncol = 1)

ggsave('age_plot.pdf',height =5,width = 6.5)



#### gender, age, and BMI

library(rcartocolor)
g1 <- gb + scale_color_carto_c(palette = "BurgYl")
g2 <- gb + scale_color_carto_c(palette = "Earth")



library(cowplot)
plot_grid(p1,p2,p3,ncol = 1)
ggsave('Data_gender_age_bmi_barplot.pdf',height =8,width = 8)



