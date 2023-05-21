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
model3 <- glm(target~  thalach + restecg + oldpeak + slope+
              cp + exang, family=binomial, data = HeartDisease)
summary(model3)
draws3 <- rmvnorm(1000, model3$coefficients, vcov(model3))

HD_df_3 <- dplyr::select(dummy_cols(HeartDisease[,c("thalach","restecg","oldpeak","slope",
                                         "cp", "exang", "target")]),
                         -c(restecg_0, slope_1, cp_1, restecg, slope, cp))
################
#Assumption 1 auf höchsten Puls für normales Ruhe-EKG, zufälligen oldpeak und slope = flat, typischem Brustschmerz und exang = 1
restecg_1 <- rep(0, 303)
restecg_2 <- rep(0, 303)
oldpeak <- round(runif(303, min = 0, max = 6.2), digits = 1)
slope_2 <- rep(1, 303)
slope_3 <- rep(0, 303)
cp_2 <- rep(0, 303)
cp_3 <- rep(0, 303)
cp_4 <- rep(0, 303)
exang <- rep(1, 303)

AgePredI_3 <- apply(draws3, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                  x[2]*HD_df_3$thalach +
                                                                  x[3]*restecg_1 + x[4]*restecg_2 +
                                                                  x[5]*oldpeak + x[6]*slope_2 +
                                                                  x[7]*slope_3 + x[8]*cp_2 +
                                                                  x[9]*cp_3 + x[10]*cp_4 +
                                                                  x[11]*exang, x[2])))
AgePredI_3_mean <- mean(AgePredI_3)

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredI_3), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))

# Expectation Plot
AgePredI_3_exp <- matrix(nrow = 132, ncol = 1000)
thalach <- c(min(HD_df_3$thalach) : max(HD_df_3$thalach))
for(i in seq_along(thalach)) {
  AgePredI_3_exp[i, ] <- apply(draws3, 1, function(x) mean(inv.logit(x[1] +
                                                                       x[2]*thalach[i] +
                                                                       x[3]*restecg_1 + x[4]*restecg_2 +
                                                                       x[5]*oldpeak + x[6]*slope_2 +
                                                                       x[7]*slope_3 + x[8]*cp_2 +
                                                                       x[9]*cp_3 + x[10]*cp_4 +
                                                                       x[11]*exang)))
}


dtI_3_exp <- data.frame(x = thalach, mean = apply(AgePredI_3_exp, 1, mean),
                         ymin=predict(loess(apply(AgePredI_3_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ thalach),data.frame(thalach = thalach)),
                         ymax=predict(loess(apply(AgePredI_3_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ thalach),data.frame(thalach = thalach)),
                         Assumption = "Assumption1"
)
pred_plotI_3 <- ggplot(dtI_3_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.15) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="höchster Pulsschlag")+
  theme_bw()

#################
# Assumption 2 auf höchstem Puls
int_leng_3 <- max(HeartDisease$thalach) - min(HeartDisease$thalach)
AgePredII_3 <- apply(draws3, 1, function(x) mean((1/131)*(inv.logit(x[1] +
                                                                      x[2]*202 +
                                                                      x[3]*HD_df_3$restecg_1 + x[4]*HD_df_3$restecg_2 +
                                                                      x[5]*HD_df_3$oldpeak + x[6]*HD_df_3$slope_2 +
                                                                      x[7]*HD_df_3$slope_3 + x[8]*HD_df_3$cp_2 +
                                                                      x[9]*HD_df_3$cp_3 + x[10]*HD_df_3$cp_4 +
                                                                      x[11]*HD_df_3$exang) -
                                                            inv.logit(x[1] +
                                                                        x[2]*71 +
                                                                        x[3]*HD_df_3$restecg_1 + x[4]*HD_df_3$restecg_2 +
                                                                        x[5]*HD_df_3$oldpeak + x[6]*HD_df_3$slope_2 +
                                                                        x[7]*HD_df_3$slope_3 + x[8]*HD_df_3$cp_2 +
                                                                        x[9]*HD_df_3$cp_3 + x[10]*HD_df_3$cp_4 +
                                                                        x[11]*HD_df_3$exang))))
AgePredII_3_mean <- mean(AgePredII_3)

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredII_3), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))


# Expectation Plot
AgePredII_3_exp <- matrix(nrow = 132, ncol = 1000)
thalach <- c(min(HD_df_3$thalach) : max(HD_df_3$thalach))
for(i in seq_along(thalach)) {
  AgePredII_3_exp[i, ] <- apply(draws3, 1, function(x) mean(inv.logit(x[1] +
                                                                        x[2]*thalach[i] +
                                                                        x[3]*HD_df_3$restecg_1 + x[4]*HD_df_3$restecg_2 +
                                                                        x[5]*HD_df_3$oldpeak + x[6]*HD_df_3$slope_2 +
                                                                        x[7]*HD_df_3$slope_3 + x[8]*HD_df_3$cp_2 +
                                                                        x[9]*HD_df_3$cp_3 + x[10]*HD_df_3$cp_4 +
                                                                        x[11]*HD_df_3$exang)))
}
dtII_3_exp <- data.frame(x = thalach, mean = apply(AgePredII_3_exp, 1, mean),
                         ymin=predict(loess(apply(AgePredII_3_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ thalach),data.frame(thalach = thalach)),
                         ymax=predict(loess(apply(AgePredII_3_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ thalach),data.frame(thalach = thalach)),
                         Assumption = "Assumption2"
)
pred_plotII_3 <- ggplot(dtII_3_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.15) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="höchster Pulsschlag")+
  theme_bw()

##########
# Assumption 3 auf Cholesterin
AgePredIII_3 <- apply(draws3, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                    x[2]*HD_df_3$thalach +
                                                                    x[3]*HD_df_3$restecg_1 + x[4]*HD_df_3$restecg_2 +
                                                                    x[5]*HD_df_3$oldpeak + x[6]*HD_df_3$slope_2 +
                                                                    x[7]*HD_df_3$slope_3 + x[8]*HD_df_3$cp_2 +
                                                                    x[9]*HD_df_3$cp_3 + x[10]*HD_df_3$cp_4 +
                                                                    x[11]*HD_df_3$exang, x[2])))
AgePredIII_3_mean <- mean(AgePredIII_3)
# wenn Alter um 1 steigt, dann verändert sich der mittlere Erwartungswert um ca.0.00553

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredIII_3), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$\\Delta_s (\\theta)$"))


# Expectation Plot
AgePredIII_3_exp <- matrix(nrow = 132, ncol = 1000)
for(i in seq_along(thalach)) {
  df_exp <- HD_df_3[which(HD_df_3$thalach == thalach[i]), ]
  AgePredIII_3_exp[i, ] <- apply(draws3, 1, function(x) mean((inv.logit(x[1] +
                                                                          x[2]*df_exp$thalach +
                                                                          x[3]*df_exp$restecg_1 + x[4]*df_exp$restecg_2 +
                                                                          x[5]*df_exp$oldpeak + x[6]*df_exp$slope_2 +
                                                                          x[7]*df_exp$slope_3 + x[8]*df_exp$cp_2 +
                                                                          x[9]*df_exp$cp_3 + x[10]*df_exp$cp_4 +
                                                                          x[11]*df_exp$exang))))
}
dtIII_3_exp <- data.frame(x = thalach, mean = apply(AgePredIII_3_exp, 1, mean),
                          ymin=predict(loess(apply(AgePredIII_3_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~thalach),data.frame(thalach=thalach)),
                          ymax=predict(loess(apply(AgePredIII_3_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~thalach),data.frame(thalach=thalach)),
                          Assumption = "Assumption3"
)
pred_plotIII_3 <- ggplot(dtIII_3_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.15) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="höchster Pulsschlag")+
  theme_bw()

############
# durchschnittliche Steigungen vergleichen
AI_3 <- data.frame(AgePredI_3, Assumption = "Assumption 1", mean = mean(AgePredI_3))
names(AI_3)<-c("value", "Assumption", "mean")

AII_3 <- data.frame(AgePredII_3, Assumption = "Assumption 2", mean = mean(AgePredII_3))
names(AII_3)<-c("value", "Assumption", "mean")

AIII_3 <- data.frame(AgePredIII_3, Assumption = "Assumption 3", mean = mean(AgePredIII_3))
names(AIII_3)<-c("value", "Assumption", "mean")

all_A_3 <- rbind(AI_3, AII_3, AIII_3)

all_A_plot_3 <- ggplot(all_A_3,aes(x = value, fill = Assumption)) +
  geom_density(alpha=0.3) + theme_minimal() +
  stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values = c(1, 15, 19)) +
  scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
  scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
  theme(legend.position = "bottom") + xlab(TeX("$\\Delta_s (\\theta)$"))

# Expectation Plot für alle Assumptions
all_A_exp_3 <- rbind(dtI_3_exp, dtII_3_exp, dtIII_3_exp)
pred_plot_all_3 <- ggplot(all_A_exp_3) +
  geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
  scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
  scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Höchster Pulsschlag")+
  theme_bw()


ggarrange(pred_plot_all_3, all_A_plot_3, nrow = 2)

