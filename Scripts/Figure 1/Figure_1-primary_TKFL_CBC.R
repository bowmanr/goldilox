setwd("/Volumes/LevineLab/Levine Lab/Bobby/Project Notes/Flt3/Mature_analysis/")
library(ggplot2)
library(dplyr)
library(gridExtra)
library(cowplot)
library(magrittr)
library(tidyr)
library(RColorBrewer)
library(plotrix)
addSmallLegend <- function(myPlot, pointSize = 0.75, textSize = 6, spaceLegend = 0.2) {
  myPlot +
    guides(shape = guide_legend(override.aes = list(size = pointSize)),
           color = guide_legend(override.aes = list(size = pointSize))) +
    theme(legend.title = element_text(size = textSize), 
          legend.text  = element_text(size = textSize),
          legend.key.size = unit(spaceLegend, "lines"))
}
options(stringsAsFactors=FALSE)
features_of_interest <- c("Sample.ID.No","Analysis.Date","RBC.M.uL.","HGB.g.dL.","HCT...","RET..K.uL.",
                          "WBC.K.uL.","NEUT..K.uL.","LYMPH..K.uL.","MONO..K.uL.",
                          "EO..K.uL.","BASO..K.uL.","NEUT....","LYMPH....","MONO....",
                          "EO....","BASO....","PLT.K.uL.","RET....")

### Read in the data
sample_key<-data.frame(read.csv("./Raw_data/Primary_TKFL_CBC/TKFL_sample_key.csv"))
sample_key[,"Date.of.TAM"]<-ifelse(sample_key$Mouse%in%c("BLY8","BLS5","BMP4"),"2019-07-10",sample_key$Date.of.TAM)
data_files <- grep("csv",list.files("./Raw_data/Primary_TKFL_CBC",full.names=TRUE),value=TRUE)
data <- list()
for(i in 1:length(data_files)){
  data[[i]] <- read.csv(data_files[i])[,features_of_interest]
  colnames(data[[i]])[1] <-"Mouse" 
  data[[i]]<-data.frame(apply(data[[i]],2,function(x){ifelse(x=="----",NA,x)}))
}

## join the data and fix any small errors
set<-inner_join(do.call(rbind,data),sample_key)
set <-data.frame(apply(set,MARGIN=2,function(column){gsub(" ","",column)}))

set$Analysis.Date<-gsub("2015","2019",set$Analysis.Date)
set$Timepoint<-as.numeric(round((as.Date(set$Analysis.Date)-as.Date(set$Date.of.TAM))/7))
set$Group <- factor(set$Genotype,levels=c("WT","Npm1","Flt3","Flt3 Npm1"))





parameter <- "WBC.K.uL."

max_value <- set%>%
                  filter(Group%in%c("Flt3","WT")&
                         Timepoint%in%c(2,4,6,8))%>%
                  pull(all_of(parameter))%>%as.numeric(.)%>%max(.,na.rm = TRUE)*1.25

gg_WBC<-ggplot(set%>%filter(Group%in%c("Flt3","WT")&
                            Timepoint%in%c(2,4,6,8)),
              aes(x=factor(Timepoint),
                   y=as.numeric(.data[[parameter]]),
                   color=Group))  +
              ylab(parameter)  +
              stat_summary(fun= median, geom = "crossbar",
                           position = position_dodge(width = 0.90), width = 0.5) + 
              stat_summary(fun.data = mean_se, geom = "errorbar",
                           position = position_dodge(width = 0.90), width = 0.5)+
              geom_point( position = position_jitterdodge(),size=0.7)+
              xlab("Weeks post TAM")+
              ylim(0,max_value)+
              scale_color_manual(values=c("WT"=brewer.pal(5,"Greys")[5],
                                          "Flt3"=brewer.pal(5,"Reds")[5]))+
              theme_classic(base_size=8)+
              theme(strip.background = element_blank())

ggsave(addSmallLegend(gg_WBC),
       width=2.5,
       height=1.75,
       units="in",
       file="./Final_Figures/Figure 1/Figure-1B_primary_WBC.pdf")

data.frame(set)%>%
                select(Group,Timepoint,all_of(parameter))%>%
                filter(Group%in%c("Flt3","WT")&
                       Timepoint%in%c(2,4,6,8)&
                       !is.na(all_of(parameter)))%>%
               group_by(Timepoint)%>%
               mutate(across(all_of(parameter), ~ as.numeric(.) ))%>%     
               pivot_wider(id_cols = c(Timepoint),names_from = Group,values_from = all_of(parameter))%>%
               mutate(mean_flt3=mean(unlist(Flt3)),
                      sd_flt3=sd(unlist(Flt3)),
                      sem_flt3=std.error(unlist(Flt3)),
                      mean_wt=mean(unlist(WT)),
                      sd_wt=sd(unlist(WT)),
                      sem_wt=std.error(unlist(WT)),
                      p_value = t.test(unlist(Flt3),unlist(WT))$p.value,
                      t_value = t.test(unlist(Flt3),unlist(WT))$statistic)%>%
              select(mean_flt3,sd_flt3,sem_flt3,mean_wt,sd_wt,sem_wt,Timepoint,p_value,t_value)%>%
              write.csv(file="./Statistics/Supp Figure-primary_WBC.csv")


data.frame(set)%>%
  select(Group,Timepoint,all_of(parameter))%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8)&
           !is.na(all_of(parameter)))%>%
  group_by(Timepoint)%>%
  mutate(across(all_of(parameter), ~ as.numeric(.) ))%>%     
  pivot_wider(id_cols = c(Timepoint),names_from = Group,values_from = all_of(parameter))%>%
  mutate(p_value = t.test(unlist(Flt3),unlist(WT))$p.value,
         t_value = t.test(unlist(Flt3),unlist(WT))$statistic)%>%
  select(Timepoint,p_value,t_value)


parameter <- "HCT..."

max_value <- set%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8))%>%
  pull(all_of(parameter))%>%as.numeric(.)%>%max(.,na.rm = TRUE)*1.25

gg_HCT<-ggplot(set%>%filter(Group%in%c("Flt3","WT")&
                              Timepoint%in%c(2,4,6,8)),
               aes(x=factor(Timepoint),
                   y=as.numeric(.data[[parameter]]),
                   color=Group))  +
  ylab("HCT%")  +
  stat_summary(fun= median, geom = "crossbar",
               position = position_dodge(width = 0.90), width = 0.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.90), width = 0.5)+
  geom_point( position = position_jitterdodge(),size=0.7)+
  xlab("Weeks post TAM")+
  ylim(0,max_value)+
  scale_color_manual(values=c("WT"=brewer.pal(5,"Greys")[5],
                              "Flt3"=brewer.pal(5,"Reds")[5]))+
  theme_classic(base_size=8)+
  theme(strip.background = element_blank())

ggsave(addSmallLegend(gg_HCT),
       width=2.5,
       height=1.75,
       units="in",
       file="./Final_Figures/Extended Data 1/Extended Data 1c_primary_HCT.pdf")

data.frame(set)%>%
  select(Group,Timepoint,all_of(parameter))%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8)&
           !is.na(all_of(parameter)))%>%
  group_by(Timepoint)%>%
  mutate(across(all_of(parameter), ~ as.numeric(.) ))%>%     
  pivot_wider(id_cols = c(Timepoint),names_from = Group,values_from = all_of(parameter))%>%
  mutate(mean_flt3=mean(unlist(Flt3)),
         sd_flt3=sd(unlist(Flt3)),
         sem_flt3=std.error(unlist(Flt3)),
         mean_wt=mean(unlist(WT)),
         sd_wt=sd(unlist(WT)),
         sem_wt=std.error(unlist(WT)),
         p_value = t.test(unlist(Flt3),unlist(WT))$p.value,
         t_value = t.test(unlist(Flt3),unlist(WT))$statistic)%>%
  select(mean_flt3,sd_flt3,sem_flt3,mean_wt,sd_wt,sem_wt,Timepoint,p_value,t_value)%>%
  write.csv(file="./Statistics/Supp Figure-primary_HCT.csv")




parameter <- "PLT.K.uL."

max_value <- set%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8))%>%
  pull(all_of(parameter))%>%as.numeric(.)%>%max(.,na.rm = TRUE)*1.25

gg_PLT<-ggplot(set%>%filter(Group%in%c("Flt3","WT")&
                              Timepoint%in%c(2,4,6,8)),
               aes(x=factor(Timepoint),
                   y=as.numeric(.data[[parameter]]),
                   color=Group))  +
  ylab("PLT.K.uL.")  +
  stat_summary(fun= median, geom = "crossbar",
               position = position_dodge(width = 0.90), width = 0.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.90), width = 0.5)+
  geom_point( position = position_jitterdodge(),size=0.7)+
  xlab("Weeks post TAM")+
  ylim(0,max_value)+
  scale_color_manual(values=c("WT"=brewer.pal(5,"Greys")[5],
                              "Flt3"=brewer.pal(5,"Reds")[5]))+
  theme_classic(base_size=8)+
  theme(strip.background = element_blank())

ggsave(addSmallLegend(gg_PLT),
       width=2.5,
       height=1.75,
       units="in",
       file="./Final_Figures/Extended Data 1/Extended Data 1c_primary_PLT.pdf")

data.frame(set)%>%
  select(Group,Timepoint,all_of(parameter))%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8)&
           !is.na(all_of(parameter)))%>%
  group_by(Timepoint)%>%
  mutate(across(all_of(parameter), ~ as.numeric(.) ))%>%     
  pivot_wider(id_cols = c(Timepoint),names_from = Group,values_from = all_of(parameter))%>%
  mutate(mean_flt3=mean(unlist(Flt3)),
         sd_flt3=sd(unlist(Flt3)),
         sem_flt3=std.error(unlist(Flt3)),
         mean_wt=mean(unlist(WT)),
         sd_wt=sd(unlist(WT)),
         sem_wt=std.error(unlist(WT)),
         p_value = t.test(unlist(Flt3),unlist(WT))$p.value,
         t_value = t.test(unlist(Flt3),unlist(WT))$statistic)%>%
  select(mean_flt3,sd_flt3,sem_flt3,mean_wt,sd_wt,sem_wt,Timepoint,p_value,t_value)%>%
  write.csv(file="./Statistics/Supp Figure-primary_PLT.csv")



parameter <- "MONO..K.uL."

max_value <- set%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8))%>%
  pull(all_of(parameter))%>%as.numeric(.)%>%max(.,na.rm = TRUE)*1.25

gg_Mono<-ggplot(set%>%filter(Group%in%c("Flt3","WT")&
                              Timepoint%in%c(2,4,6,8)),
               aes(x=factor(Timepoint),
                   y=as.numeric(.data[[parameter]]),
                   color=Group))  +
  ylab("MONO..K.uL.")  +
  stat_summary(fun= median, geom = "crossbar",
               position = position_dodge(width = 0.90), width = 0.5) + 
  stat_summary(fun.data = mean_se, geom = "errorbar",
               position = position_dodge(width = 0.90), width = 0.5)+
  geom_point( position = position_jitterdodge(),size=0.7)+
  xlab("Weeks post TAM")+
  ylim(0,max_value)+
  scale_color_manual(values=c("WT"=brewer.pal(5,"Greys")[5],
                              "Flt3"=brewer.pal(5,"Reds")[5]))+
  theme_classic(base_size=8)+
  theme(strip.background = element_blank())

ggsave(addSmallLegend(gg_Mono),
       width=2.5,
       height=1.75,
       units="in",
       file="./Final_Figures/Extended Data 1/Extended Data 1c_primary_Mono.pdf")

data.frame(set)%>%
  select(Group,Timepoint,all_of(parameter))%>%
  filter(Group%in%c("Flt3","WT")&
           Timepoint%in%c(2,4,6,8)&
           !is.na(all_of(parameter)))%>%
  group_by(Timepoint)%>%
  mutate(across(all_of(parameter), ~ as.numeric(.) ))%>%     
  pivot_wider(id_cols = c(Timepoint),names_from = Group,values_from = all_of(parameter))%>%
  mutate(mean_flt3=mean(unlist(Flt3)),
         sd_flt3=sd(unlist(Flt3)),
         sem_flt3=std.error(unlist(Flt3)),
         mean_wt=mean(unlist(WT)),
         sd_wt=sd(unlist(WT)),
         sem_wt=std.error(unlist(WT)),
         p_value = t.test(unlist(Flt3),unlist(WT))$p.value,
         t_value = t.test(unlist(Flt3),unlist(WT))$statistic)%>%
  select(mean_flt3,sd_flt3,sem_flt3,mean_wt,sd_wt,sem_wt,Timepoint,p_value,t_value)%>%
  write.csv(file="./Statistics/Supp Figure-primary_Mono.csv")
  
