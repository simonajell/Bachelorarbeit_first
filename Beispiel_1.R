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

# Daten einlesen
data(HeartDisease.target,HeartDisease.cont,HeartDisease.cat)
HeartDisease<-cbind(HeartDisease.cont,HeartDisease.cat,HeartDisease.target)
# Zielvariable binär kodieren
HeartDisease$target<-ifelse(HeartDisease$num==0,0,1)
# Brustschmerz Variable zu numerischen Faktor umwandeln
HeartDisease$cp_cat<-as.factor(HeartDisease$cp)
levels(HeartDisease$cp_cat)<-list(asymptomatic=4, typical_angina=1,
                                  atypical_angina=2, non_anginal_pain=3)
attach(HeartDisease)


# Modell schätzen
# Variablen, die Angaben zu physischen Eigenschaften und Allgemeinzustand geben, welche in Zusammenhang mit einer Herzkrankheit stehen
# Einfluss von Cholesterin, Geschlecht, Ruhe-Blutdruck, Alter und nüchterner Blutzucker auf Herzerkrankung
model1<-glm(target ~ chol + sex + trestbps +
              age + fbs, family=binomial, data=HeartDisease)
summary(model1)
draws1 <- rmvnorm(1000, model1$coefficients, vcov(model1))

################
# Assumption 1 auf Cholesterin für Männer allen Alters mit zufälligem Ruhe Blutdruck und Blutzucker
# Datensatz mit allen Variablen kombiniert
age_I_1 <- round(runif(50, min = 29, max = 77))
trestbps_I_1 <- round(runif(50, min = 94, max = 200))
fbs_I_1 <- round(runif(20, min = 0, max = 1))
sex_I_1 <- 1
dt_I_1_comb <- expand.grid(age = age_I_1, trestbps = trestbps_I_1, fbs = fbs_I_1, sex = sex_I_1)
# GME schätzen
AgePredI_1 <- apply(draws1, 1, function(x) mean((1/438)*(inv.logit(x[1] +
                                                                      x[2]*564 +
                                                                      x[3]*dt_I_1_comb$sex + x[4]*dt_I_1_comb$trestbps +
                                                                      x[5]*dt_I_1_comb$age + x[6]*dt_I_1_comb$fbs) -
                                                            inv.logit(x[1] +
                                                                        x[2]*126 +
                                                                        x[3]*dt_I_1_comb$sex + x[4]*dt_I_1_comb$trestbps +
                                                                        x[5]*dt_I_1_comb$age + x[6]*dt_I_1_comb$fbs))))

AgePredI_1_mean <- mean(AgePredI_1)

# Adjusted Predictions schätzen
AgePredI_1_exp <- matrix(nrow = 439, ncol = 1000)
chol <- c(min(HeartDisease$chol) : max(HeartDisease$chol))
for(i in seq_along(chol)) {
  AgePredI_1_exp[i, ] <-apply(draws1, 1, function(x) mean(inv.logit(x[1] +
                                                                       x[2]*chol[i] +
                                                                       x[3]*dt_I_1_comb$sex + x[4]*dt_I_1_comb$trestbps +
                                                                       x[5]*dt_I_1_comb$age + x[6]*dt_I_1_comb$fbs)))
}

mean(apply(draws1, 1, function(x) mean(inv.logit(x[1] +
                                              x[2]*250 +
                                              x[3]*dt_I_1_comb$sex + x[4]*dt_I_1_comb$trestbps +
                                              x[5]*dt_I_1_comb$age + x[6]*dt_I_1_comb$fbs))))

dtI_1_exp <- data.frame(x = chol, mean = apply(AgePredI_1_exp, 1, mean),
                        ymin=predict(loess(apply(AgePredI_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ chol),data.frame(chol = chol)),
                        ymax=predict(loess(apply(AgePredI_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ chol),data.frame(chol = chol)),
                        Assumption = "Annahme 1"
)

#################
# Assumption 2 auf Cholesterin
# Intervall von Cholesterin
int_leng_1 <- max(HeartDisease$chol) - min(HeartDisease$chol)
# GME schätzen
AgePredII_1 <- apply(draws1, 1, function(x) mean((1/438)*(inv.logit(x[1] +
                                                                      x[2]*564 +
                                                                      x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                                                      x[5]*HeartDisease$age + x[6]*HeartDisease$fbs) -
                                                            inv.logit(x[1] +
                                                                        x[2]*126 +
                                                                        x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                                                        x[5]*HeartDisease$age + x[6]*HeartDisease$fbs))))
AgePredII_1_mean <- mean(AgePredII_1)

# Adjusted Predictions
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
                         Assumption = "Annahme 2"
)

##########
# Assumption 3 auf Cholesterin
# GME schätzen
AgePredIII_1 <- apply(draws1, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                    x[2]*HeartDisease$chol +
                                                                    x[3]*HeartDisease$sex + x[4]*HeartDisease$trestbps +
                                                                    x[5]*HeartDisease$age + x[6]*HeartDisease$fbs, x[2])))
AgePredIII_1_mean <- mean(AgePredIII_1)

# Adjusted Predictions schätzen
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
                          Assumption = "Annahme 3"
)

############
# GMEs vergleichen
AI_1 <- data.frame(AgePredI_1, Assumption = "Annahme 1", mean = mean(AgePredI_1))
names(AI_1)<-c("value", "Assumption", "mean")

AII_1 <- data.frame(AgePredII_1, Assumption = "Annahme 2", mean = mean(AgePredII_1))
names(AII_1)<-c("value", "Assumption", "mean")

AIII_1 <- data.frame(AgePredIII_1, Assumption = "Annahme 3", mean = mean(AgePredIII_1))
names(AIII_1)<-c("value", "Assumption", "mean")

all_A <- rbind(AI_1, AII_1, AIII_1)

all_A_plot <- ggplot(all_A,aes(x = value, fill = Assumption)) +
  geom_density(alpha=0.3) +   theme_bw() +
  stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values = c(1, 15, 19)) +
  scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
  scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) + xlab(TeX("$\\Delta_j$"))
ggsave("gme_plot_1.jpg", width = 7, height = 4)

# Adjusted Predictions unter allen Assumptions vergleichen
all_A_exp <- rbind(dtI_1_exp, dtII_1_exp, dtIII_1_exp)
pred_plot_all <- ggplot(all_A_exp) +
  geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
  geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
  scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
  scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
  labs(y="Wert der Adjusted Predictions",x="Cholesterin")+
  theme_bw()
ggsave("exp_plot_1.jpg", width = 7, height = 4)


