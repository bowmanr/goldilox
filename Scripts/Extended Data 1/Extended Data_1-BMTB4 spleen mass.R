
BMTB4_spleen_weights<-ggplot(spleen_mass,aes(x=Time.point,y=Spleen.Weights,fill=Time.point))+
  geom_boxplot()+
  scale_fill_manual(values=c("Control"=brewer.pal(5,"Greys")[2],
                             "2w"=brewer.pal(5,"Reds")[2],
                             "4w"=brewer.pal(5,"Reds")[3],
                             "6w"=brewer.pal(5,"Reds")[4],
                             "8w"=brewer.pal(5,"Reds")[5]))+  theme_minimal(base_size=12)+xlab("Weeks post TAM")+ylab("mass (mg)")+
  theme_classic(base_size=8)+
  guides(fill="none")+
  theme(strip.background = element_blank())
ggsave(BMTB4_spleen_weights,width = 2.25,height=1.75,
       file="./Figures/ExFig1D-BMTB4_spleen_weights.pdf")




data.frame(spleen_mass)%>%
  select(Time.point,Spleen.Weights)%>%
  filter(Time.point%in%c("8w","Control"))%>%
  group_by(Time.point)%>%
  mutate(across(Spleen.Weights, ~ as.numeric(.) ))%>%     
  pivot_wider(id_cols = c(Time.point),names_from = Time.point,values_from = Spleen.Weights)%>%
    mutate(mean_flt3=mean(unlist(`8w`)),
         sd_flt3=sd(unlist(`8w`)),
         mean_wt=mean(unlist(Control)),
         sd_wt=sd(unlist(Control)),
         p_value = t.test(unlist(`8w`),unlist(Control))$p.value,
         t_value = t.test(unlist(`8w`),unlist(Control))$statistic)%>%
  select(mean_flt3,sd_flt3,mean_wt,sd_wt,p_value,t_value)%>%
  write.csv(file="./Statistics/Extended Data 1-Spleen Mass.csv")
