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
levels(HeartDisease$cp_cat)<-list(asymptomatic=4, typical_angina=1,
                                  atypical_angina=2, non_anginal_pain=3)

attach(HeartDisease)


########## Modell für diskretes Merkmal
# Variablen, die Angaben zu physischen Eigenschaften und Allgemeinzustand geben, welche in Zusammenhang mit einer Herzkrankheit stehen
# Einfluss von Cholesterin, Geschlecht, Ruhe Blutdruck, Alter und Nüchterner Blutzucker auf Herzerkrankung
model2 <- glm(target ~ sex + chol + trestbps +
              age + fbs, family=binomial, data=HeartDisease)
summary(model2)
draws2 <- rmvnorm(1000, model2$coefficients, vcov(model2))

# neuer Datensatz mit Dummy Variablen
HD_df <- HeartDisease[,c("sex","chol","trestbps","age", "fbs", "target")]


#Assumption 1 auf Cholesterin für 30 bis 50 Jährige Personen mit hohem Blutzucker
chol_I_2 <- round(runif(303, min = 126, max = 564))
trestbps_I_2 <- round(runif(303, min = 94, max = 200))
fbs_I_2 <- round(runif(303, min = 1, max = 1))
age_I_2 <- round(runif(303, min = 30, max = 50))


# Expectation Plot
age_I_1 <- round(runif(303, min = 30, max = 50))
trestbps_I_1 <- round(runif(303, min = 94, max = 200))
fbs_I_1 <- round(runif(303, min = 0, max = 1))

intval_I_2<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$sex+betas[3]*chol_I_2+betas[4]*trestbps_I_2+
    betas[5]*age_I_1+betas[6]*fbs_I_2
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}

mittlerer_Ewert_I_2_female <- apply(draws2,1,function(x){sum(intval_I_2(x,mutate(HD_df,sex= 0)))/nrow(HD_df)})
mittlerer_Ewert_I_2_male <- apply(draws2,1,function(x){sum(intval_I_2(x,mutate(HD_df,sex = 1)))/nrow(HD_df)})
mittlerer_Ewert_I_2_sex <- data.frame()
mittlerer_Ewert_I_2_sex <- data.frame("Male" = mittlerer_Ewert_I_2_male, "Female" = mittlerer_Ewert_I_2_female)
mittlerer_Ewert_I_2_sex <- ldply(mittlerer_Ewert_I_2_sex, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_I_2_sex)<-c("Sex","value")
mittlerer_Ewert_I_2_sex <- merge(mittlerer_Ewert_I_2_sex,
                                  mittlerer_Ewert_I_2_sex %>%
                                    group_by(Sex)%>%
                                    dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                  by="Sex")
mittlerer_Ewert_I_2_sex <- data.frame(mittlerer_Ewert_I_2_sex, Assumption = "Assumption 1")
cat_1 <- ggplot(mittlerer_Ewert_I_2_sex, aes(x = value, y = Sex)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.I')")

# generalized marginal Effect Plot
AI_2 <- list()
AI_2<-apply(draws2,1,function(x){sum(intval_I_2(x,mutate(HD_df,sex=1))-intval_I_2(x,mutate(HD_df,sex=0)))/nrow(HD_df)})
AI_2_mean <- mean(AI_2)
AI_2 <- data.frame(value = AI_2, Assumption = "Assumption 1")


# Assumption 2
intval_II_2<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$sex+betas[3]*regs$chol+betas[4]*regs$trestbps+
    betas[5]*regs$age+betas[6]*regs$fbs
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}
mittlerer_Ewert_II_2_female <- apply(draws2,1,function(x){sum(intval_II_2(x,mutate(HD_df,sex = 0)))/nrow(HD_df)})
mittlerer_Ewert_II_2_male <- apply(draws2,1,function(x){sum(intval_II_2(x,mutate(HD_df,sex = 1)))/nrow(HD_df)})
mittlerer_Ewert_II_2_sex <- data.frame()
mittlerer_Ewert_II_2_sex <- data.frame("Male" = mittlerer_Ewert_II_2_male, "Female" = mittlerer_Ewert_II_2_female)

mittlerer_Ewert_II_2_sex <- ldply(mittlerer_Ewert_II_2_sex, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_II_2_sex)<-c("Sex","value")
mittlerer_Ewert_II_2_sex <- merge(mittlerer_Ewert_II_2_sex,
                              mittlerer_Ewert_II_2_sex %>%
                                group_by(Sex)%>%
                                dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                              by="Sex")
mittlerer_Ewert_II_2_sex <- data.frame(mittlerer_Ewert_II_2_sex, Assumption = "Assumption 2")
cat_2 <- ggplot(mittlerer_Ewert_II_2_sex, aes(x = value, y = Sex)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.II')")

# generalized marginal Effect Plot
AII_2 <- list()
AII_2<-apply(draws2,1,function(x){sum(intval_II_2(x,mutate(HD_df,sex=1))-intval_II_2(x,mutate(HD_df,sex=0)))/nrow(HD_df)})
AII_2_mean <- mean(AII_2)
AII_2 <- data.frame(value = AII_2, Assumption = "Assumption 2")



# Assumption 3
mittlerer_Ewert_III_2_male <- apply(draws2,1,function(x){sum(intval_II_2(x,HD_df[which(HD_df$sex==1),]))/length(which(HD_df$sex==1))})
mittlerer_Ewert_III_2_female <- apply(draws2,1,function(x){sum(intval_II_2(x,HD_df[which(HD_df$sex==0),]))/length(which(HD_df$sex==0))})
mittlerer_Ewert_III_2_sex <- data.frame()
mittlerer_Ewert_III_2_sex <- data.frame("Male" = mittlerer_Ewert_III_2_male, "Female" = mittlerer_Ewert_III_2_female)
mittlerer_Ewert_III_2_sex <- ldply(mittlerer_Ewert_III_2_sex, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_III_2_sex)<-c("Sex","value")
mittlerer_Ewert_III_2_sex <- merge(mittlerer_Ewert_III_2_sex,
                                  mittlerer_Ewert_III_2_sex %>%
                                    group_by(Sex)%>%
                                    dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
                                  by="Sex")
mittlerer_Ewert_III_2_sex <- data.frame(mittlerer_Ewert_III_2_sex, Assumption = "Assumption 3")
cat_3 <- ggplot(mittlerer_Ewert_III_2_sex, aes(x = value, y = Sex)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.III')")

# generalized marginal Effect Plot
AIII_2 <- list()
AIII_2<-apply(draws2,1,function(x){(sum(intval_II_2(x,HD_df[which(HD_df$sex==1),]))/length(which(HD_df$sex==1)))-
    (sum(intval_II_2(x,HD_df[which(HD_df$sex==0),]))/length(which(HD_df$sex==0)))})
AIII_2_mean <- mean(AIII_2)
AIII_2 <- data.frame(value = AIII_2, Assumption = "Assumption 3")



###########
# alle generalisierten marginalen Effekte zusammen
all_A_2 <- rbind(AI_2, AII_2, AIII_2)
all_A_plot_2 <- ggplot(all_A_2,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  scale_color_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))

# alle Expectation Plots zusammen
all_A_exp_2 <- rbind(mittlerer_Ewert_I_2_sex, mittlerer_Ewert_II_2_sex, mittlerer_Ewert_III_2_sex)
pred_plot_all_2 <- ggplot(all_A_exp_2,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+facet_grid(Sex~.)+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("firebrick2" ,"orange", "steelblue2"))+
  scale_color_manual(values=c("firebrick2" ,"orange", "steelblue2"))+
  theme(legend.position = "bottom")+xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))


ggarrange(pred_plot_all_2, all_A_plot_2, nrow = 2)


