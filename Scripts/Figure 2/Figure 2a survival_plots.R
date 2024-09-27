setwd("/Volumes/LevineLab/Levine Lab/Bobby/Project Notes/Flt3/")
library(GGally)
library(survival)
library(ggplotify)
library(magick)
library(survminer)


primary <-read.csv("./Mature_analysis/Raw_data/Survival/TKFL_primary.csv")[1:36,]
primary$Genotype <- factor(primary$Genotype,levels=c("WT","Flt3","Npm1","Flt3 Npm1"))

fit <- survfit(Surv(time, status) ~ Genotype, data = primary)
primary_KM<-ggsurv(fit,size.est=1,order.legend=FALSE)+
  theme_classic(base_size=8)+
  scale_color_manual(values=c("WT"="grey60",
                              "Flt3"=brewer.pal(5,"Reds")[5],
                              "Npm1"=brewer.pal(5,"Blues")[5],
                              "Flt3 Npm1"=brewer.pal(5,"Greens")[5]))+
  scale_y_continuous(lim=c(0,1))+
  theme(plot.title = element_text(hjust=0.5),
        legend.position =  "none")+
  xlab("Weeks post TAM")

ggsave(primary_KM,width=2,height=1.75,
       file="./Mature_analysis/Final_Figures/Figure 2/Figure 2a.pdf")

fit <- survfit(Surv(time, status) ~ Genotype, data = primary%>%filter(Genotype%in%c("Flt3","Flt3 Npm1")))

surv_pvalue(fit)

