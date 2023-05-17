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


# Modell schätzen
# Einfluss von Cholesterin, sex, Ruhe Blutdruck, Alter und Nüchterner Blutzucker auf Herzerkrankung
model1<-glm(target ~ chol + sex + trestbps +
              age + fbs, family=binomial, data=HeartDisease)
summary(model1)
draws1 <- rmvnorm(1000, model1$coefficients, vcov(model1))


#Assumption 1 auf Cholesterin für 30 bis 50 Jährige Männer mit dem zufälligem Ruhe Blutdruck und Blutzucker
age_I_1 <- round(runif(303, min = 30, max = 50))
trestbps_I_1 <- round(runif(303, min = 94, max = 200))
fbs_I_1 <- round(runif(303, min = 0, max = 1))
AgePredI_1 <- apply(draws1, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                            x[2]*HeartDisease$chol +
                                                            x[3]* rep(1, 303) + x[4]*trestbps_I_1 +
                                                            x[5]* age_I_1 + x[6]*fbs_I_1, x[2])))
AgePredI_1_mean <- mean(AgePredI_1)

# Expectation Plot
AgePredI_1_exp <- matrix(nrow = 439, ncol = 1000)
chol <- c(min(HeartDisease$chol) : max(HeartDisease$chol))
for(i in seq_along(chol)) {
  AgePredI_1_exp[i, ] <- apply(draws1, 1, function(x) mean(inv.logit(x[1] +
                                                                        x[2]*chol[i] +
                                                                        x[3]*rep(1, 303) + x[4]*trestbps_I_1 +
                                                                        x[5]*age_I_1 + x[6]*fbs_I_1)))
}
dtI_1_exp <- data.frame(x = chol, mean = apply(AgePredI_1_exp, 1, mean),
                         ymin=predict(loess(apply(AgePredI_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ chol),data.frame(chol = chol)),
                         ymax=predict(loess(apply(AgePredI_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ chol),data.frame(chol = chol)),
                        Assumption = "Assumption1"
)
pred_plotI_1 <- ggplot(dtI_1_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.05) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Cholesterin")+
  theme_bw()


#################
# Assumption 2 auf Cholesterin
int_leng_1 <- max(HeartDisease$chol) - min(HeartDisease$chol)
AgePredII_1 <- apply(draws1, 1, function(x) mean((1/438)*(inv.logit(x[1] +
                                              x[2]*564 +
                                              x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                              x[5]*HeartDisease$age + x[6]*HeartDisease$fbs) -
                                                            inv.logit(x[1] +
                                               x[2]*126 +
                                               x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                               x[5]*HeartDisease$age + x[6]*HeartDisease$fbs))))
AgePredII_1_mean <- mean(AgePredII_1)

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredII_1), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))


# Expectation Plot
AgePredII_1_exp <- matrix(nrow = 439, ncol = 1000)
chol <- c(min(HeartDisease$chol) : max(HeartDisease$chol))
for(i in seq_along(chol)) {
  AgePredII_1_exp[i, ] <- apply(draws1, 1, function(x) mean(inv.logit(x[1] +
                                                                            x[2]*chol[i] +
                                                                            x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                                                            x[5]*HeartDisease$age + x[6]*HeartDisease$fbs)))
}
dtII_1_exp <- data.frame(x = chol, mean = apply(AgePredII_1_exp, 1, mean),
                        ymin=predict(loess(apply(AgePredII_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ chol),data.frame(chol = chol)),
                        ymax=predict(loess(apply(AgePredII_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ chol),data.frame(chol = chol)),
                        Assumption = "Assumption2"
)
pred_plotII_1 <- ggplot(dtII_1_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.05) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Cholesterin")+
  theme_bw()

##########
# Assumption 3 auf Cholesterin
AgePredIII_1 <- apply(draws1, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                    x[2]*HeartDisease$chol +
                                                                    x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                                                    x[5]*HeartDisease$age + x[6]*HeartDisease$fbs, x[2])))
AgePredIII_1_mean <- mean(AgePredIII_1)
# wenn Alter um 1 steigt, dann verändert sich der mittlere Erwartungswert um ca.0.00553

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredIII_1), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$\\Delta_s (\\theta)$"))


# Expectation Plot
AgePredIII_1_exp <- matrix(nrow = 439, ncol = 1000)
for(i in seq_along(chol)) {
  df_exp <- HeartDisease[which(HeartDisease$chol == chol[i]), ]
  AgePredIII_1_exp[i, ] <- apply(draws1, 1, function(x) mean((inv.logit(x[1] +
                                                                          x[2]*df_exp$chol +
                                                                          x[3]*df_exp$sex + x[4]*df_exp$trestbps +
                                                                          x[5]*df_exp$age + x[6]*df_exp$fbs))))
}
dtIII_1_exp <- data.frame(x = chol, mean = apply(AgePredIII_1_exp, 1, mean),
                         ymin=predict(loess(apply(AgePredIII_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~chol),data.frame(chol=chol)),
                         ymax=predict(loess(apply(AgePredIII_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~chol),data.frame(chol=chol)),
                         Assumption = "Assumption3"
)
pred_plotIII_1 <- ggplot(dtIII_1_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.05) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Cholesterin")+
  theme_bw()

############
# durchschnittliche Steigungen vergleichen
AI_1 <- data.frame(AgePredI_1, Assumption = "Assumption 1", mean = mean(AgePredI_1))
names(AI_1)<-c("value", "Assumption", "mean")

AII_1 <- data.frame(AgePredII_1, Assumption = "Assumption 2", mean = mean(AgePredII_1))
names(AII_1)<-c("value", "Assumption", "mean")

AIII_1 <- data.frame(AgePredIII_1, Assumption = "Assumption 3", mean = mean(AgePredIII_1))
names(AIII_1)<-c("value", "Assumption", "mean")

all_A <- rbind(AI_1, AII_1, AIII_1)

all_A_plot <- ggplot(all_A,aes(x = value, fill = Assumption)) +
  geom_density(alpha=0.3) + theme_minimal() +
  stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values = c(1, 15, 19)) +
  scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
  scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
  theme(legend.position = "bottom") + xlab(TeX("$\\Delta_s (\\theta)$"))

# Expectation Plot für alle Assumptions
all_A_exp <- rbind(dtI_1_exp, dtII_1_exp, dtIII_1_exp)
pred_plot_all <- ggplot(all_A_exp) +
  geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.05) +
  scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
  scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
  labs(y=TeX("$\\Delta_s (\\theta)$"),x="Cholesterin")+
  theme_bw()

########## Modell für diskretes Merkmal
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
ggplot(all_A_2,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  scale_color_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  theme(legend.position = "bottom")+xlab(TeX("$\\Delta_s (\\theta)$"))

# alle Expectation Plots zusammen
all_A_exp_2 <- rbind(mittlerer_Ewert_I_2_sex, mittlerer_Ewert_II_2_sex, mittlerer_Ewert_III_2_sex)
ggplot(all_A_exp_2,aes(x=value,fill=Assumption))+
  geom_density(alpha=0.3)+theme_minimal()+
  stat_pointinterval(aes(color=Assumption,shape=Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  easy_remove_y_axis()+facet_grid(Sex~.)+
  scale_shape_manual(values=c(1, 15,19))+
  scale_fill_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  scale_color_manual(values=c("yellowgreen", "deepskyblue4", "darkorchid1"))+
  theme(legend.position = "bottom")+xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))


