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
                        ymax=predict(loess(apply(AgePredII_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ chol),data.frame(chol = chol))
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
AgePredIII_1_mean <- mean(AgePredIII_1) # 0.00572 ( durchschnittliche Steigung)
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
                         ymax=predict(loess(apply(AgePredIII_1_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~chol),data.frame(chol=chol))
)
pred_plotIII_1 <- ggplot(dtIII_1_exp) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.05) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Cholesterin")+
  theme_bw()

############
# durchschnittliche Steigungen vergleichen
AII_1 <- data.frame(AgePredII_1, Assumption = "Assumption 2", mean = mean(AgePredII_1))
names(AII_1)<-c("value", "Assumption", "mean")

AIII_1 <- data.frame(AgePredIII_1, Assumption = "Assumption 3", mean = mean(AgePredIII_1))
names(AIII_1)<-c("value", "Assumption", "mean")

both_A <- rbind(AII_1, AIII_1)

both_A_plot <- ggplot(both_A,aes(x = value, fill = Assumption)) +
  geom_density(alpha=0.3) + theme_minimal() +
  stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values = c(15, 19)) +
  scale_fill_manual(values = c("yellowgreen", "mediumpurple4")) +
  scale_color_manual(values = c("yellowgreen", "mediumpurple4")) +
  theme(legend.position = "bottom") + xlab(TeX("$\\Delta_s (\\theta)$"))
