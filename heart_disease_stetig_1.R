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
# Variablen, die Angaben zu physischen Eigenschaften und Allgemeinzustand geben, welche in Zusammenhang mit einer Herzkrankheit stehen
# Einfluss von Cholesterin, Geschlecht, Ruhe Blutdruck, Alter und Nüchterner Blutzucker auf Herzerkrankung
model1<-glm(target ~ chol + sex + trestbps +
              age + fbs, family=binomial, data=HeartDisease)
summary(model1)
draws1 <- rmvnorm(1000, model1$coefficients, vcov(model1))

################
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
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.15) +
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
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.15) +
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
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.15) +
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
  scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) + xlab(TeX("$\\Delta_s (\\theta)$"))

# Expectation Plot für alle Assumptions
all_A_exp <- rbind(dtI_1_exp, dtII_1_exp, dtIII_1_exp)
pred_plot_all <- ggplot(all_A_exp) +
  geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
  scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
  scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Cholesterin")+
  theme_bw()


plot_1 <- ggarrange(pred_plot_all, all_A_plot, nrow = 2)
annotate_figure(plot_1, top = text_grob("Effekt von Cholesterin auf die Wsk. eine Herzkrankheit zu haben",
                face = "bold", size = 14))


