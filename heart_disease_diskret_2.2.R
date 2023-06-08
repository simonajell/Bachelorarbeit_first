library(MixAll)
library(glm.predict)
library(ggplot2)
library(tidybayes)
library(latex2exp)
library(ggeasy)
library(readxl)
library(fastDummies)
library(dplyr)
library(plyr)
library(ggpubr)
library(mvtnorm)

data(HeartDisease.target,HeartDisease.cont,HeartDisease.cat)

HeartDisease<-cbind(HeartDisease.cont,HeartDisease.cat,HeartDisease.target)

HeartDisease$target<-ifelse(HeartDisease$num==0,0,1)
HeartDisease$cp_cat<-as.factor(HeartDisease$cp)
HeartDisease$restecg<-as.factor(HeartDisease$restecg)
HeartDisease$slope<-as.factor(HeartDisease$slope)
HeartDisease$cp<-as.factor(HeartDisease$cp)
levels(HeartDisease$cp_cat)<-list(asymptomatic=4, typical_angina=1,
                                  atypical_angina=2, non_anginal_pain=3)

attach(HeartDisease)

# Variablen aus anderen Kategorien: Elektrokardiographie bezogene Ergebnisse und Symptome
# Einfluss von höchstem Puls, Ruhe-EKG, ST-Senkung, Steigung des ST Abschnitts,
# Art des Brust Schmerzes und Brustschmerz durch Anstrengung auf Herzerkrankung
model4.2 <- glm(target ~ cp + restecg + thalach + exang + oldpeak + slope,
                family=binomial, data = HeartDisease)
summary(model4.2)
draws_4.2 <- rmvnorm(1000, model4.2$coefficients, vcov(model4.2))

HD_df_4 <- dplyr::select(dummy_cols(HeartDisease[,c("thalach","restecg","oldpeak","slope",
                                                    "cp", "exang", "target")]),
                         -c(restecg, slope, cp))

################
#Assumption 1 auf Brustschmerz für zufälligen oldpeak und slope = downsloping, EKG mit Hypertrophie,
# zufälligem höchstem Ruhepuls und zufälligem Aktivitätsschmerz
exang <- round(runif(303, min = 0, max = 1))
oldpeak <- round(runif(303, min = 0, max = 6.2), digits = 1)
slope_2 <- rep(0, 303)
slope_3 <- rep(1, 303)
restecg_1 <- rep(0, 303)
restecg_2 <- rep(1, 303)
thalach <- round(runif(303, min = 71, max = 202))
df_4.2 <- data.frame("cp_1" = HD_df_4$cp_1, "cp_2" = HD_df_4$cp_2,
                    "cp_3" =  HD_df_4$cp_3, "cp_4" = HD_df_4$cp_4,
                     restecg_1, restecg_2, thalach, exang, oldpeak, slope_2, slope_3)

# adjusted predictions auf selbe art und weise wie Assumption 2 nur mit den selbst ausgesuchten Werten
intval_I_4.2<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$cp_2+betas[3]*regs$cp_3+betas[4]*regs$cp_4+
    betas[5]*restecg_1+betas[6]*restecg_2 + betas[7]*thalach + betas[8]*exang +
    betas[9]*oldpeak + betas[10]*slope_2 + betas[11]*slope_3
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}

mittlerer_Ewert_I_4.2 <- list()
mittlerer_Ewert_I_4.2$cp_1 <- apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,
                                                                            dplyr::select(mutate(df_4.2,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(df_4.2)})
mittlerer_Ewert_I_4.2$cp_2 <- apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,
                                                                            dplyr::select(mutate(df_4.2,cp_1=0,cp_2=1,cp_3=0, cp_4=0),-"cp_1")))/nrow(df_4.2)})
mittlerer_Ewert_I_4.2$cp_3 <- apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,
                                                                            dplyr::select(mutate(df_4.2,cp_1=0,cp_2=0,cp_3=1, cp_4=0),-"cp_1")))/nrow(df_4.2)})
mittlerer_Ewert_I_4.2$cp_4 <- apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,
                                                                           dplyr::select(mutate(df_4.2,cp_1=0,cp_2=0,cp_3=0, cp_4=1),-"cp_1")))/nrow(df_4.2)})

mittlerer_Ewert_I_4.2_cp <- ldply(mittlerer_Ewert_I_4.2, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_I_4.2_cp)<-c("Brust_Schmerz","value")
mittlerer_Ewert_I_4.2_cp<-merge(mittlerer_Ewert_I_4.2_cp,
                                mittlerer_Ewert_I_4.2_cp %>%
                                     group_by(Brust_Schmerz)%>%
                                     dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                   by="Brust_Schmerz")
levels(mittlerer_Ewert_I_4.2_cp$Brust_Schmerz)<- list(Typical_Angina  = "cp_1", Atypical_Angina  = "cp_2", Non_Anginal_Pain  = "cp_3", Asymptomatic = "cp_4")
mittlerer_Ewert_I_4.2_cp <- data.frame(mittlerer_Ewert_I_4.2_cp, Assumption = "Annahme 1")

exp_I_4.2 <- ggplot(mittlerer_Ewert_I_4.2_cp, aes(x = value, y = Brust_Schmerz)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.I)")

# generalized marginal Effect Plot
AI_4.2 <- list()
AI_4.2$cp_2<-apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,dplyr::select(mutate(df_4.2,cp_1=0,cp_2=1,cp_3=0, cp_4=0),-"cp_1"))-
                                                  intval_I_4.2(x,dplyr::select(mutate(df_4.2,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(df_4.2)})
AI_4.2$cp_3<-apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,dplyr::select(mutate(df_4.2,cp_1=0,cp_2=0,cp_3=1, cp_4=0),-"cp_1"))-
                                                  intval_I_4.2(x,dplyr::select(mutate(df_4.2,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(df_4.2)})
AI_4.2$cp_4<-apply(draws_4.2,1,function(x){sum(intval_I_4.2(x,dplyr::select(mutate(df_4.2,cp_1=0,cp_2=0,cp_3=0, cp_4=1),-"cp_1"))-
                                                 intval_I_4.2(x,dplyr::select(mutate(df_4.2,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(df_4.2)})
AI_4.2_cp<-ldply(AI_4.2, data.frame) %>% mutate(.id=as.factor(.id))
names(AI_4.2_cp)<-c("Brust_Schmerz","value")
AI_4.2_cp<-merge(AI_4.2_cp,
                 AI_4.2_cp %>%
                      group_by(Brust_Schmerz)%>%
                      dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                    by="Brust_Schmerz")
AI_4.2_cp <- data.frame(AI_4.2_cp, Assumption = "Annahme 1")
GME_I_4.2<-ggplot(AI_4.2_cp,aes(x=value,fill=Brust_Schmerz))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Brust_Schmerz,shape=Brust_Schmerz),
                     position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  scale_shape_manual(values=c(1,15,19))+
  scale_fill_manual(values=c("gold2","green4","navyblue"))+scale_color_manual(values=c("gold2","green4","navyblue"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))

#########################
# Assumption 2
intval_II_4.2<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <- betas[1]+betas[2]*regs$cp_2+betas[3]*regs$cp_3+betas[4]*regs$cp_4+
    betas[5]*regs$restecg_1+betas[6]*regs$restecg_2 + betas[7]*regs$thalach + betas[8]*regs$exang +
    betas[9]*regs$oldpeak + betas[10]*regs$slope_2 + betas[11]*regs$slope_3
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}
mittlerer_Ewert_II_4.2 <- list()
mittlerer_Ewert_II_4.2$cp_1 <- apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,
                                                                              dplyr::select(mutate(HD_df_4,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(HD_df_4)})
mittlerer_Ewert_II_4.2$cp_2 <- apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,
                                                                               dplyr::select(mutate(HD_df_4,cp_1=0,cp_2=1,cp_3=0, cp_4=0),-"cp_1")))/nrow(HD_df_4)})
mittlerer_Ewert_II_4.2$cp_3 <- apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,
                                                                               dplyr::select(mutate(HD_df_4,cp_1=0,cp_2=0,cp_3=1, cp_4=0),-"cp_1")))/nrow(HD_df_4)})
mittlerer_Ewert_II_4.2$cp_4 <- apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,
                                                                               dplyr::select(mutate(HD_df_4,cp_1=0,cp_2=0,cp_3=0, cp_4=1),-"cp_1")))/nrow(HD_df_4)})

mittlerer_Ewert_II_4.2_cp <- ldply(mittlerer_Ewert_II_4.2, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_II_4.2_cp)<-c("Brust_Schmerz","value")
mittlerer_Ewert_II_4.2_cp<-merge(mittlerer_Ewert_II_4.2_cp,
                                 mittlerer_Ewert_II_4.2_cp %>%
                                      group_by(Brust_Schmerz)%>%
                                      dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                    by="Brust_Schmerz")
levels(mittlerer_Ewert_II_4.2_cp$Brust_Schmerz)<- list(Typical_Angina  = "cp_1", Atypical_Angina  = "cp_2", Non_Anginal_Pain  = "cp_3", Asymptomatic = "cp_4")
mittlerer_Ewert_II_4.2_cp <- data.frame(mittlerer_Ewert_II_4.2_cp, Assumption = "Annahme 2")

exp_II_4.2 <- ggplot(mittlerer_Ewert_II_4.2_cp, aes(x = value, y = Brust_Schmerz)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.II)")

# generalized marginal Effect Plot
AII_4.2 <- list()
AII_4.2$cp_2<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(mutate(HD_df_4,cp_1=0,cp_2=1,cp_3=0, cp_4=0),-"cp_1"))-
                                                   intval_II_4.2(x,dplyr::select(mutate(HD_df_4,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(HD_df_4)})
AII_4.2$cp_3<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(mutate(HD_df_4,cp_1=0,cp_2=0,cp_3=1, cp_4=0),-"cp_1"))-
                                                  intval_II_4.2(x,dplyr::select(mutate(HD_df_4,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(HD_df_4)})
AII_4.2$cp_4<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(mutate(HD_df_4,cp_1=0,cp_2=0,cp_3=0, cp_4=1),-"cp_1"))-
                                                  intval_II_4.2(x,dplyr::select(mutate(HD_df_4,cp_1=1,cp_2=0,cp_3=0, cp_4=0),-"cp_1")))/nrow(HD_df_4)})
AII_4.2_cp<-ldply(AII_4.2, data.frame) %>% mutate(.id=as.factor(.id))
names(AII_4.2_cp)<-c("Brust_Schmerz","value")
AII_4.2_cp<-merge(AII_4.2_cp,
                  AII_4.2_cp %>%
                       group_by(Brust_Schmerz)%>%
                       dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                     by="Brust_Schmerz")
AII_4.2_cp <- data.frame(AII_4.2_cp, Assumption = "Annahme 2")
GME_II_4.2 <- ggplot(AII_4.2_cp,aes(x=value,fill=Brust_Schmerz))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Brust_Schmerz,shape=Brust_Schmerz),
                     position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("gold2","green4","navyblue"))+scale_color_manual(values=c("gold2","green4","navyblue"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))

#######################
# Assumption 3
mittlerer_Ewert_III_4.2 <- list()
mittlerer_Ewert_III_4.2$cp_1<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_1==1),],-"cp_1")))/length(which(HD_df_4$cp_1==1))})
mittlerer_Ewert_III_4.2$cp_2<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_2==1),],-"cp_1")))/length(which(HD_df_4$cp_2==1))})
mittlerer_Ewert_III_4.2$cp_3<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_3==1),],-"cp_1")))/length(which(HD_df_4$cp_3==1))})
mittlerer_Ewert_III_4.2$cp_4<-apply(draws_4.2,1,function(x){sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_4==1),],-"cp_1")))/length(which(HD_df_4$cp_4==1))})


mittlerer_Ewert_III_4.2_cp <- ldply(mittlerer_Ewert_III_4.2, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_III_4.2_cp)<-c("Brust_Schmerz","value")
mittlerer_Ewert_III_4.2_cp<-merge(mittlerer_Ewert_III_4.2_cp,
                                  mittlerer_Ewert_III_4.2_cp %>%
                                       group_by(Brust_Schmerz)%>%
                                       dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                     by="Brust_Schmerz")
levels(mittlerer_Ewert_III_4.2_cp$Brust_Schmerz)<- list(Typical_Angina  = "cp_1", Atypical_Angina  = "cp_2", Non_Anginal_Pain  = "cp_3", Asymptomatic = "cp_4")
mittlerer_Ewert_III_4.2_cp <- data.frame(mittlerer_Ewert_III_4.2_cp, Assumption = "Annahme 3")

exp_III_4.2 <- ggplot(mittlerer_Ewert_III_4.2_cp, aes(x = value, y = Brust_Schmerz)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.III)")


# generalized marginal Effect Plot
AIII_4.2 <- list()
AIII_4.2$cp_2<-apply(draws_4.2,1,function(x){(sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_2==1),],-"cp_1")))/length(which(HD_df_4$cp_2==1)))-
    (sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_1==1),],-"cp_1")))/length(which(HD_df_4$cp_1==1)))})
AIII_4.2$cp_3<-apply(draws_4.2,1,function(x){(sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_3==1),],-"cp_1")))/length(which(HD_df_4$cp_3==1)))-
    (sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_1==1),],-"cp_1")))/length(which(HD_df_4$cp_1==1)))})
AIII_4.2$cp_4<-apply(draws_4.2,1,function(x){(sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_4==1),],-"cp_1")))/length(which(HD_df_4$cp_4==1)))-
    (sum(intval_II_4.2(x,dplyr::select(HD_df_4[which(HD_df_4$cp_1==1),],-"cp_1")))/length(which(HD_df_4$cp_1==1)))})

AIII_4.2_cp<-ldply(AIII_4.2, data.frame) %>% mutate(.id=as.factor(.id))
names(AIII_4.2_cp)<-c("Brust_Schmerz","value")
AIII_4.2_cp<-merge(AIII_4.2_cp,
                   AIII_4.2_cp %>%
                        group_by(Brust_Schmerz)%>%
                        dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                      by="Brust_Schmerz")
AIII_4.2_cp <- data.frame(AIII_4.2_cp, Assumption = "Annahme 3")

GME_III_4.2<-ggplot(AIII_4.2_cp,aes(x=value,fill=Brust_Schmerz))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Brust_Schmerz,shape=Brust_Schmerz),
                     position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("gold2","green4","navyblue"))+
  scale_color_manual(values=c("gold2","green4","navyblue"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))


###########
# alle generalisierten marginalen Effekte zusammen
all_A_4.2 <- rbind(AI_4.2_cp, AII_4.2_cp, AIII_4.2_cp)
all_A_4.2$Brust_Schmerz <- factor(all_A_4.2$Brust_Schmerz, labels = c("Atypischer Brustschmerz; n = 23",
                                                                      "nicht-anginöser Brustschmerz; n = 50",
                                                                      "Asymptomatischer Brustschmerz; n = 86"))

all_A_plot_4.2 <- ggplot(all_A_4.2,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_bw()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+facet_grid(Brust_Schmerz~., labeller = label_wrap_gen(width=10))+
  scale_shape_manual(values=c(1, 15,19))+
  theme(strip.text.y = element_text(angle = 0)) +
  scale_fill_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  scale_color_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  xlab(TeX("Wert des GME"))
ggsave("gme_plot_4.2.jpg", width = 7, height = 4)



# alle adjusted predictions zusammen
all_A_exp_4.2 <- rbind(mittlerer_Ewert_I_4.2_cp, mittlerer_Ewert_II_4.2_cp, mittlerer_Ewert_III_4.2_cp)
all_A_exp_4.2$Brust_Schmerz <- factor(all_A_exp_4.2$Brust_Schmerz, labels = c("Typischer Brustschmerz; n = 144",
"Atypischer Brustschmerz; n = 23", "nicht-anginöser Brustschmerz; n = 50", "Asymptomatischer Brustschmerz; n = 86"))
pred_plot_all_4.2 <- ggplot(all_A_exp_4.2,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_bw()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  facet_grid(Brust_Schmerz~., labeller = label_wrap_gen(width=10))+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("firebrick2" ,"orange", "steelblue2"))+
  scale_color_manual(values=c("firebrick2" ,"orange", "steelblue2"))+
  theme(strip.text.y = element_text(angle = 0)) +
  xlab("Wert der Angepassten Vorhersage")
ggsave("exp_plot_4.2.jpg", width = 7, height = 4)

plot_4.2 <- ggarrange(pred_plot_all_4.2, all_A_plot_4.2, nrow = 2)
annotate_figure(plot_4.2, top = text_grob("Effekt von Art des Brustschmerzes auf die Wsk. eine Herzkrankheit zu haben",
                                        face = "bold", size = 14))
