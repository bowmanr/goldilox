setwd("/Volumes/LevineLab/Levine Lab/Bobby/Project Notes/Flt3/Mature_analysis/")
library(gridExtra)
library(cowplot)
library(magrittr)
library(RColorBrewer)
library(tidyverse)

addSmallLegend <- function(myPlot, pointSize = 0.75, textSize = 6, spaceLegend = 0.2) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
options(stringsAsFactors=FALSE)

lin<-read.delim("/Volumes/LevineLab/Levine Lab/Bobby/Project Notes/Flt3/Mature_analysis/Raw_data/BMTB4_compettitive_TKFL/lin_output.txt",sep="\t")

gg_myeloid_flt3_Time<-lin%>%mutate(myeloid=(Cd11b_posGr1_Hi+Cd11b_posGr1_Lo+Cd11b_posGr1_Neg)/Heme*100)%>%
            filter(Sample.!="Blood_BMT1123 7670_016.fcs")%>%
            filter(Tissue=="Blood")%>%
            mutate(Group=factor(gsub(" week","w",Group),levels=c("Control","2w","4w","6w","8w")))%>%
            ggplot(aes(x=Group,y=myeloid,fill=Group))+
            geom_boxplot()+
            ylab("%Cd11b+ of Cd45.2")+
            scale_y_continuous(limits = c(0,60))+
            scale_fill_manual(values=c("Control"=brewer.pal(5,"Greys")[2],
                                        "2w"=brewer.pal(5,"Reds")[2],
                                        "4w"=brewer.pal(5,"Reds")[3],
                                        "6w"=brewer.pal(5,"Reds")[4],
                                        "8w"=brewer.pal(5,"Reds")[5]))+
            theme_classic(base_size=8)+
            guides(fill="none")+
            theme(strip.background = element_blank())


gg_kit_flt3_Time<-lin%>%mutate(myeloid=(Cd11b_posGr1_Hi+Cd11b_posGr1_Lo+Cd11b_posGr1_Neg)/Heme*100)%>%
  filter(Sample.!="Blood_BMT1123 7670_016.fcs")%>%
  filter(Tissue=="Blood")%>%
  mutate(Group=factor(gsub(" week","w",Group),levels=c("Control","2w","4w","6w","8w")))%>%
  ggplot(aes(x=Group,y=Kit/Heme*100,fill=Group))+
  geom_boxplot()+
  ylab("%cKIT+ of Cd45.2")+
  scale_y_continuous(limits = c(0,20))+
  scale_fill_manual(values=c("Control"=brewer.pal(5,"Greys")[2],
                             "2w"=brewer.pal(5,"Reds")[2],
                             "4w"=brewer.pal(5,"Reds")[3],
                             "6w"=brewer.pal(5,"Reds")[4],
                             "8w"=brewer.pal(5,"Reds")[5]))+
  theme_classic(base_size=8)+
  guides(fill="none")+
  theme(strip.background = element_blank())

ggsave(plot_grid(gg_myeloid_flt3_Time,gg_kit_flt3_Time),file="./Figures/ExF1_C-gg_myeloid_flt3_Time.pdf",width=3.25,height=1.75)


lin%>%mutate(myeloid=(Cd11b_posGr1_Hi+Cd11b_posGr1_Lo+Cd11b_posGr1_Neg)/Heme*100)%>%
  filter(Sample.!="Blood_BMT1123 7670_016.fcs")%>%
  filter(Tissue=="Blood")%>%
  mutate(Group=factor(gsub(" week","w",Group),levels=c("Control","2w","4w","6w","8w")))%>%
  mutate(kit_percent=Kit/Heme*100)%>%
  group_by(Group)%>%
  summarize(mean=mean(kit_percent),
            sd=sd(kit_percent))

lin%>%mutate(myeloid=(Cd11b_posGr1_Hi+Cd11b_posGr1_Lo+Cd11b_posGr1_Neg)/Heme*100)%>%
  filter(Sample.!="Blood_BMT1123 7670_016.fcs")%>%
  filter(Tissue=="Blood")%>%
  mutate(Group=factor(gsub(" week","w",Group),levels=c("Control","2w","4w","6w","8w")))%>%
  filter(Group%in%c("Control","6w"))%>%
  mutate(kit_percent=Kit/Heme*100)%>%
  #group_by(Group)%>%
  mutate(across(kit_percent, ~ as.numeric(.) ))%>%     
  pivot_wider(id_cols = c(Group),names_from = Group,values_from = kit_percent)%>%
  mutate(mean_flt3=mean(unlist(`6w`)),
         sd_flt3=sd(unlist(`6w`)),
         mean_wt=mean(unlist(Control)),
         sd_wt=sd(unlist(Control)),
         p_value = t.test(unlist(`6w`),unlist(Control))$p.value,
         t_value = t.test(unlist(`6w`),unlist(Control))$statistic)%>%
  select(mean_flt3,sd_flt3,mean_wt,sd_wt,p_value,t_value)%>%
  data.frame
