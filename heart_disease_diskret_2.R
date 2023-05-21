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
model4 <- glm(target ~ restecg + thalach + exang + oldpeak + slope+
                cp, family=binomial, data = HeartDisease)
summary(model4)
draws_4 <- rmvnorm(1000, model4$coefficients, vcov(model4))

HD_df_4 <- dplyr::select(dummy_cols(HeartDisease[,c("thalach","restecg","oldpeak","slope",
                                                    "cp", "exang", "target")]),
                         -c(restecg, slope, cp))

################
#Assumption 1 auf Ruhe EKG für zufälligen oldpeak und slope = downsloping, typischem Brustschmerz,
# zufälligem höchstem Ruhepuls und zufälligem Aktivitätsschmerz
exang <- round(runif(303, min = 0, max = 1))
oldpeak <- round(runif(303, min = 0, max = 6.2), digits = 1)
slope_2 <- rep(0, 303)
slope_3 <- rep(1, 303)
cp_2 <- rep(0, 303)
cp_3 <- rep(0, 303)
cp_4 <- rep(0, 303)
thalach <- round(runif(303, min = 71, max = 202))

# Expectation Plot
intval_I_4<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$restecg_1+betas[3]*regs$restecg_2+betas[4]*thalach+
    betas[5]*exang+betas[6]*oldpeak + betas[7]*slope_2 + betas[8]*slope_3 +
    betas[9]*cp_2 + betas[10]*cp_3 + betas[11]*cp_4
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}

mittlerer_Ewert_I_4 <- list()
mittlerer_Ewert_I_4$restecg_0 <- apply(draws_4,1,function(x){sum(intval_I_4(x,
               dplyr::select(mutate(HD_df_4,restecg_0=1,restecg_1=0,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
mittlerer_Ewert_I_4$restecg_1 <- apply(draws_4,1,function(x){sum(intval_I_4(x,
               dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=1,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
mittlerer_Ewert_I_4$restecg_2 <- apply(draws_4,1,function(x){sum(intval_I_4(x,
               dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=0,restecg_2=1),-"restecg_0")))/nrow(HD_df_4)})

mittlerer_Ewert_I_4_restecg <- ldply(mittlerer_Ewert_I_4, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_I_4_restecg)<-c("Ruhe_EKG","value")
mittlerer_Ewert_I_4_restecg<-merge(mittlerer_Ewert_I_4_restecg,
                                   mittlerer_Ewert_I_4_restecg %>%
               group_by(Ruhe_EKG)%>%
               dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
             by="Ruhe_EKG")
levels(mittlerer_Ewert_I_4_restecg$Ruhe_EKG)<- list(Normal  = "restecg_0", Abnormal  = "restecg_1", Hypertrophie  = "restecg_2")
mittlerer_Ewert_I_4_restecg <- data.frame(mittlerer_Ewert_I_4_restecg, Assumption = "Assumption 1")

exp_I_4 <- ggplot(mittlerer_Ewert_I_4_restecg, aes(x = value, y = Ruhe_EKG)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.I')")

# generalized marginal Effect Plot
AI_4 <- list()

AI_4$restecg_1<-apply(draws_4,1,function(x){sum(intval_I_4(x,dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=1,restecg_2=0),-"restecg_0"))-
                                                  intval_I_4(x,dplyr::select(mutate(HD_df_4,restecg_0=1,restecg_1=0,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
AI_4$restecg_2<-apply(draws_4,1,function(x){sum(intval_I_4(x,dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=0,restecg_2=1),-"restecg_0"))-
                                                  intval_I_4(x,dplyr::select(mutate(HD_df_4,restecg_0=1,restecg_1=0,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
AI_4_restecg<-ldply(AI_4, data.frame) %>% mutate(.id=as.factor(.id))
names(AI_4_restecg)<-c("Ruhe_EKG","value")
AI_4_restecg<-merge(AI_4_restecg,
                    AI_4_restecg %>%
                 group_by(Ruhe_EKG)%>%
                 dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
               by="Ruhe_EKG")
AI_4_restecg <- data.frame(AI_4_restecg, Assumption = "Assumption 1")
GME_I_4<-ggplot(AI_4_restecg,aes(x=value,fill=Ruhe_EKG))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Ruhe_EKG,shape=Ruhe_EKG),
      position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  scale_shape_manual(values=c(15,19))+
  scale_fill_manual(values=c("green4","navyblue"))+scale_color_manual(values=c("green4","navyblue"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))

#########################
# Assumption 2
intval_II_4<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$restecg_1+betas[3]*regs$restecg_2+betas[4]*regs$thalach+
    betas[5]*regs$exang+betas[6]*regs$oldpeak + betas[7]*regs$slope_2 + betas[8]*regs$slope_3 +
    betas[9]*regs$cp_2 + betas[10]*regs$cp_3 + betas[11]*regs$cp_4
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}
mittlerer_Ewert_II_4 <- list()
mittlerer_Ewert_II_4$restecg_0 <- apply(draws_4,1,function(x){sum(intval_II_4(x,
            dplyr::select(mutate(HD_df_4,restecg_0=1,restecg_1=0,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
mittlerer_Ewert_II_4$restecg_1 <- apply(draws_4,1,function(x){sum(intval_II_4(x,
            dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=1,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
mittlerer_Ewert_II_4$restecg_2 <- apply(draws_4,1,function(x){sum(intval_II_4(x,
            dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=0,restecg_2=1),-"restecg_0")))/nrow(HD_df_4)})

mittlerer_Ewert_II_4_restecg <- ldply(mittlerer_Ewert_II_4, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_II_4_restecg)<-c("Ruhe_EKG","value")
mittlerer_Ewert_II_4_restecg<-merge(mittlerer_Ewert_II_4_restecg,
                                   mittlerer_Ewert_II_4_restecg %>%
                                     group_by(Ruhe_EKG)%>%
                                     dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                   by="Ruhe_EKG")
levels(mittlerer_Ewert_II_4_restecg$Ruhe_EKG)<- list(Normal  = "restecg_0", Abnormal  = "restecg_1", Hypertrophie  = "restecg_2")
mittlerer_Ewert_II_4_restecg <- data.frame(mittlerer_Ewert_II_4_restecg, Assumption = "Assumption 2")

exp_II_4 <- ggplot(mittlerer_Ewert_II_4_restecg, aes(x = value, y = Ruhe_EKG)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.II)")

# generalized marginal Effect Plot
AII_4 <- list()
AII_4$restecg_1<-apply(draws_4,1,function(x){sum(intval_II_4(x,dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=1,restecg_2=0),-"restecg_0"))-
                                                  intval_II_4(x,dplyr::select(mutate(HD_df_4,restecg_0=1,restecg_1=0,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
AII_4$restecg_2<-apply(draws_4,1,function(x){sum(intval_II_4(x,dplyr::select(mutate(HD_df_4,restecg_0=0,restecg_1=0,restecg_2=1),-"restecg_0"))-
                                                  intval_II_4(x,dplyr::select(mutate(HD_df_4,restecg_0=1,restecg_1=0,restecg_2=0),-"restecg_0")))/nrow(HD_df_4)})
AII_4_restecg<-ldply(AII_4, data.frame) %>% mutate(.id=as.factor(.id))
names(AII_4_restecg)<-c("Ruhe_EKG","value")
AII_4_restecg<-merge(AII_4_restecg,
                    AII_4_restecg %>%
                      group_by(Ruhe_EKG)%>%
                      dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                    by="Ruhe_EKG")
AII_4_restecg <- data.frame(AII_4_restecg, Assumption = "Assumption 2")
GME_II_4<-ggplot(AII_4_restecg,aes(x=value,fill=Ruhe_EKG))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Ruhe_EKG,shape=Ruhe_EKG),
                     position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  scale_shape_manual(values=c(15,19))+
  scale_fill_manual(values=c("green4","navyblue"))+scale_color_manual(values=c("green4","navyblue"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))

#######################
# Assumption 3
mittlerer_Ewert_III_4 <- list()


mittlerer_Ewert_III_4$restecg_0<-apply(draws_4,1,function(x){sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_0==1),],-"restecg_0")))/length(which(HD_df_4$restecg_0==1))})
mittlerer_Ewert_III_4$restecg_1<-apply(draws_4,1,function(x){sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_1==1),],-"restecg_0")))/length(which(HD_df_4$restecg_1==1))})
mittlerer_Ewert_III_4$restecg_2<-apply(draws_4,1,function(x){sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_2==1),],-"restecg_0")))/length(which(HD_df_4$restecg_2==1))})

mittlerer_Ewert_III_4_restecg <- ldply(mittlerer_Ewert_III_4, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_III_4_restecg)<-c("Ruhe_EKG","value")
mittlerer_Ewert_III_4_restecg<-merge(mittlerer_Ewert_III_4_restecg,
                                    mittlerer_Ewert_III_4_restecg %>%
                                      group_by(Ruhe_EKG)%>%
                                      dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                    by="Ruhe_EKG")
levels(mittlerer_Ewert_III_4_restecg$Ruhe_EKG)<- list(Normal  = "restecg_0", Abnormal  = "restecg_1", Hypertrophie  = "restecg_2")
mittlerer_Ewert_III_4_restecg <- data.frame(mittlerer_Ewert_III_4_restecg, Assumption = "Assumption 3")

exp_III_4 <- ggplot(mittlerer_Ewert_III_4_restecg, aes(x = value, y = Ruhe_EKG)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.III)")


# generalized marginal Effect Plot
AIII_4 <- list()
AIII_4$restecg_1<-apply(draws_4,1,function(x){(sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_1==1),],-"restecg_0")))/length(which(HD_df_4$restecg_1==1)))-
    (sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_0==1),],-"restecg_0")))/length(which(HD_df_4$restecg_0==1)))})
AIII_4$restecg_2<-apply(draws_4,1,function(x){(sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_2==1),],-"restecg_0")))/length(which(HD_df_4$restecg_2==1)))-
    (sum(intval_II_4(x,dplyr::select(HD_df_4[which(HD_df_4$restecg_0==1),],-"restecg_0")))/length(which(HD_df_4$restecg_0==1)))})

AIII_4_restecg<-ldply(AIII_4, data.frame) %>% mutate(.id=as.factor(.id))
names(AIII_4_restecg)<-c("Ruhe_EKG","value")
AIII_4_restecg<-merge(AIII_4_restecg,
                     AIII_4_restecg %>%
                       group_by(Ruhe_EKG)%>%
                       dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                     by="Ruhe_EKG")
AIII_4_restecg <- data.frame(AIII_4_restecg, Assumption = "Assumption 3")

GME_III_4<-ggplot(AIII_4_restecg,aes(x=value,fill=Ruhe_EKG))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Ruhe_EKG,shape=Ruhe_EKG),
                     position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+
  scale_shape_manual(values=c(15,19))+
  scale_fill_manual(values=c("green4","navyblue"))+scale_color_manual(values=c("green4","navyblue"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))


###########
# alle generalisierten marginalen Effekte zusammen
all_A_4 <- rbind(AI_4_restecg, AII_4_restecg, AIII_4_restecg)
all_A_plot_4 <- ggplot(all_A_4,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+facet_grid(Ruhe_EKG~.)+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  scale_color_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))




# alle Expectation Plots zusammen
all_A_exp_4 <- rbind(mittlerer_Ewert_I_4_restecg, mittlerer_Ewert_II_4_restecg, mittlerer_Ewert_III_4_restecg)
pred_plot_all_4 <- ggplot(all_A_exp_4,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+facet_grid(Ruhe_EKG~.)+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("firebrick2" ,"orange", "steelblue2"))+
  scale_color_manual(values=c("firebrick2" ,"orange", "steelblue2"))+
  theme(legend.position = "bottom")+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"))

table(HD_df_4$restecg_1)
# Es gibt unr sehr wenige Beobachtungen, die einen abnormalen Schmerz haben,
# wodurch man eine sehr große Unischerheit hat und sehr verschiedene Werte bekommt
# Außerdem unterscheidet sich Assumption 1 stark von den anderen Assumptions, was
# vielleicht daran liegt, dass die selbst gewählten Werte unrealistisch sind.

ggarrange(pred_plot_all_4, all_A_plot_4, nrow = 2)



