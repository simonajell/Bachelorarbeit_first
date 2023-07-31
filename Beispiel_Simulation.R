library(boot)
library(abind)
library(tibble)
library(mvtnorm)
library(ggplot2)
library(tidybayes)
library(latex2exp)
library(ggpubr)

# Analyse basiert auf der Studie von:
  # Body size and weight change over adulthood and risk of breast cancer by menopausal
  # and hormone receptor status: a pooled analysis of 20 prospective cohort studies
# Variablen: Größe, Gewicht, BMI, BMI im frühen Erwachsenenalter (18-20 J.),
  # Gewichtszunahme im Erwachsenenalter (ab 18-20 J.)
# Zielvariable: Postmenopausales Brustkrebs Risiko

set.seed(12226699)

# Helper Funktionen
inv.logit<-function(x){ #this function computes the inverse logit of a value x
  res<-numeric(length(x))
  for(i in 1:length(x)){
    res[i]<-1/(1+exp(-x[i]))
  }
  return(res)
}

inv.logit.deriv<-function(x,beta){ #this function computes the derivative of the inverse logit of a value x
  res<-numeric(length(x))
  for(i in 1:length(x)){
    res[i]<-(beta*exp(x[i]))/(1+exp(x[i]))^2
  }
  return(res)
}


# GME berechnen auf dem Datensatz und Modell für Height <= split_value
GME_Sim_smlr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
  # Mehrmals Simulieren
  sim_list <- list()
  for(i in seq_along(1:100)) {
    # Regressoren simulieren
    Height = round(rnorm(1000, 162.5, 100))
    BMI = round(rnorm(1000, 27, 28), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    Weight_Gain = (BMI - BMI_EA) * ((Height * 0.01)^2)
    X <- data.frame(Height, BMI, BMI_EA, Weight_Gain)
    # wahre Betas festlegen
    betatrue<-c(beta1, beta2, beta3, beta4, beta5)
    # Zielvariable schätzen
    Y <- apply(X,1,function(x)rbinom(1,1,inv.logit(betatrue[1]+betatrue[2]*x[1]+betatrue[3]*x[2]+betatrue[4]*x[3] + betatrue[5]*x[4])))
    # Datensatz bilden
    sim_list[[i]]<-cbind(Y,X)
  }
  sim_matrix <- abind(sim_list, along=3)
  data_S <- as.data.frame(apply(sim_matrix, c(1,2), mean))
  data_S <- data.frame(Y = round(data_S$Y),
                       Height = round(data_S$Height),
                       BMI = round(data_S$BMI, 2),
                       BMI_EA = round(data_S$BMI_EA, 2),
                       Weight_Gain = round(data_S$Weight_Gain, 2))
  # Teile Datensatz an bestimmter Stelle
  data_S_1 <- data_S[data_S$Height <= split_value, ]

  # Modell schätzen
  modelS1<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_1)
  drawsS1 <- rmvnorm(1000, modelS1$coefficients, vcov(modelS1))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI zwischen 20 und 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # GME
  Height_I_exp <- seq(min(Height_I), max(Height_I), by = 1)
  GME_I_S <<- apply(drawsS1, 1, function(x) mean((1/length(Height_I_exp))*(inv.logit(x[1] +
                                                                       x[2]*max(Height_I_exp) +
                                                                       x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                       x[5]*data_I_S$Weight_Gain) -
                                                             inv.logit(x[1] +
                                                                         x[2]*min(Height_I_exp) +
                                                                         x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                         x[5]*data_I_S$Weight_Gain))))
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  GME_II_S <<- apply(drawsS1, 1, function(x) mean((1/length(Height_S))*(inv.logit(x[1] +
                                                                        x[2]* max(Height_S)+
                                                                        x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                        x[5]*data_S_1$Weight_Gain) -
                                                              inv.logit(x[1] +
                                                                          x[2]*min(Height_S) +
                                                                          x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                          x[5]*data_S_1$Weight_Gain))))
  # Annahme 3
  GME_III_S <<- apply(drawsS1, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                      x[2]*data_S_1$Height +
                                                                      x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                      x[5]*data_S_1$Weight_Gain, x[2])))
  # GMEs vergleichen
  AI_S <- data.frame(GME_I_S, Assumption = "Annahme 1", mean = mean(GME_I_S))
  names(AI_S)<-c("value", "Assumption", "mean")

  AII_S <- data.frame(GME_II_S, Assumption = "Annahme 2", mean = mean(GME_II_S))
  names(AII_S)<-c("value", "Assumption", "mean")

  AIII_S <- data.frame(GME_III_S, Assumption = "Annahme 3", mean = mean(GME_III_S))
  names(AIII_S)<-c("value", "Assumption", "mean")

  all_A_S <- rbind(AI_S, AII_S, AIII_S)

  all_A_S_plot <<- ggplot(all_A_S,aes(x = value, fill = Assumption)) +
    geom_density(alpha=0.3) +   theme_bw() +
    stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
    scale_shape_manual(values = c(1, 15, 19)) +
    scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
    scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
    xlab(TeX("$\\Delta_j$")) +
    xlim(-0.005, 0.015)
  ggsave("exp_plot_S1.jpg", width = 7, height = 4)
  return(all_A_S_plot)
}
set.seed(12226699)
GME_Sim_smlr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

################################################################################
# GME berechnen auf dem Datensatz und Modell für Height > split_value
GME_Sim_lrgr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
  # Mehrmals Simulieren
  sim_list <- list()
  for(i in seq_along(1:100)) {
    # Regressoren simulieren
    Height = round(rnorm(1000, 162.5, 100))
    BMI = round(rnorm(1000, 27, 28), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    Weight_Gain = (BMI - BMI_EA) * ((Height * 0.01)^2)
    X <- data.frame(Height, BMI, BMI_EA, Weight_Gain)
    # wahre Betas festlegen
    betatrue<-c(beta1, beta2, beta3, beta4, beta5)
    # Zielvariable schätzen
    Y<-apply(X,1,function(x)rbinom(1,1,inv.logit(betatrue[1]+betatrue[2]*x[1]+betatrue[3]*x[2]+betatrue[4]*x[3] + betatrue[5]*x[4])))
    # Datensatz bilden
    sim_list[[i]]<-cbind(Y,X)
  }
  sim_matrix <- abind(sim_list, along=3)
  data_S <- as.data.frame(apply(sim_matrix, c(1,2), mean))
  data_S <- data.frame(Y = round(data_S$Y),
                       Height = round(data_S$Height),
                       BMI = round(data_S$BMI, 2),
                       BMI_EA = round(data_S$BMI_EA, 2),
                       Weight_Gain = round(data_S$Weight_Gain, 2))
  # Teile Datensatz an bestimmter Stelle
  # Teile Datensatz an bestimmter Stelle
  data_S_1 <- data_S[data_S$Height > split_value, ]

  # Modell schätzen
  modelS1<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_1)
  drawsS1 <- rmvnorm(1000, modelS1$coefficients, vcov(modelS1))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI größer 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # GME
  Height_I_exp <- seq(min(Height_I), max(Height_I))
  GME_I_S2 <<- apply(drawsS1, 1, function(x) mean((1/length(Height_I_exp))*(inv.logit(x[1] +
                                                                                  x[2]*max(Height_I_exp) +
                                                                                  x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                                  x[5]*data_I_S$Weight_Gain) -
                                                                        inv.logit(x[1] +
                                                                                    x[2]*min(Height_I_exp) +
                                                                                    x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                                    x[5]*data_I_S$Weight_Gain))))
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  GME_II_S2 <<- apply(drawsS1, 1, function(x) mean((1/length(Height_S))*(inv.logit(x[1] +
                                                                                  x[2]* max(Height_S)+
                                                                                  x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                                  x[5]*data_S_1$Weight_Gain) -
                                                                        inv.logit(x[1] +
                                                                                    x[2]*min(Height_S) +
                                                                                    x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                                    x[5]*data_S_1$Weight_Gain))))
  # Annahme 3
  GME_III_S2 <<- apply(drawsS1, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                    x[2]*data_S_1$Height +
                                                                    x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                    x[5]*data_S_1$Weight_Gain, x[2])))
  # GMEs vergleichen
  AI_S <- data.frame(GME_I_S2, Assumption = "Annahme 1", mean = mean(GME_I_S2))
  names(AI_S)<-c("value", "Assumption", "mean")

  AII_S <- data.frame(GME_II_S2, Assumption = "Annahme 2", mean = mean(GME_II_S2))
  names(AII_S)<-c("value", "Assumption", "mean")

  AIII_S <- data.frame(GME_III_S2, Assumption = "Annahme 3", mean = mean(GME_III_S2))
  names(AIII_S)<-c("value", "Assumption", "mean")

  all_A_S <- rbind(AI_S, AII_S, AIII_S)

  all_A_S2_plot <<- ggplot(all_A_S,aes(x = value, fill = Assumption)) +
    geom_density(alpha=0.3) +   theme_bw() +
    stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
    scale_shape_manual(values = c(1, 15, 19)) +
    scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
    scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
    xlab(TeX("$\\Delta_j$")) +
    xlim(-0.005, 0.015)
  ggsave("exp_plot_S2.jpg", width = 7, height = 4)
  return(all_A_S2_plot)
}
set.seed(12226699)
GME_Sim_lrgr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

################################################################################
# GME berechnen auf dem Datensatz  für Height <= split_value und Modell auf Height > split_value
GME_Sim_smlr_lrgr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
  # Mehrmals Simulieren
  sim_list <- list()
  for(i in seq_along(1:100)) {
    # Regressoren simulieren
    Height = round(rnorm(1000, 162.5, 100))
    BMI = round(rnorm(1000, 27, 28), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    Weight_Gain = (BMI - BMI_EA) * ((Height * 0.01)^2)
    X <- data.frame(Height, BMI, BMI_EA, Weight_Gain)
    # wahre Betas festlegen
    betatrue<-c(beta1, beta2, beta3, beta4, beta5)
    # Zielvariable schätzen
    Y <- apply(X,1,function(x)rbinom(1,1,inv.logit(betatrue[1]+betatrue[2]*x[1]+betatrue[3]*x[2]+betatrue[4]*x[3] + betatrue[5]*x[4])))
    # Datensatz bilden
    sim_list[[i]]<-cbind(Y,X)
  }
  sim_matrix <- abind(sim_list, along=3)
  data_S <- as.data.frame(apply(sim_matrix, c(1,2), mean))
  data_S <- data.frame(Y = round(data_S$Y),
                       Height = round(data_S$Height),
                       BMI = round(data_S$BMI, 2),
                       BMI_EA = round(data_S$BMI_EA, 2),
                       Weight_Gain = round(data_S$Weight_Gain, 2))
  # Teile Datensatz an bestimmter Stelle
  data_S_1 <- data_S[data_S$Height <= split_value, ]
  data_S_2 <- data_S[data_S$Height > split_value, ]

  # Modell schätzen
  modelS<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_2)
  drawsS <- rmvnorm(1000, modelS$coefficients, vcov(modelS))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI größer 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # GME
  Height_I_exp <- seq(min(Height_I), max(Height_I))
  GME_I_S3 <<- apply(drawsS, 1, function(x) mean((1/length(Height_I_exp))*(inv.logit(x[1] +
                                                                                  x[2]*max(Height_I_exp) +
                                                                                  x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                                  x[5]*data_I_S$Weight_Gain) -
                                                                        inv.logit(x[1] +
                                                                                    x[2]*min(Height_I_exp) +
                                                                                    x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                                    x[5]*data_I_S$Weight_Gain))))
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  GME_II_S3 <<- apply(drawsS, 1, function(x) mean((1/length(Height_S))*(inv.logit(x[1] +
                                                                                  x[2]* max(Height_S)+
                                                                                  x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                                  x[5]*data_S_1$Weight_Gain) -
                                                                        inv.logit(x[1] +
                                                                                    x[2]*min(Height_S) +
                                                                                    x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                                    x[5]*data_S_1$Weight_Gain))))
  # Annahme 3
  GME_III_S3 <<- apply(drawsS, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                    x[2]*data_S_1$Height +
                                                                    x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                    x[5]*data_S_1$Weight_Gain, x[2])))
  # GMEs vergleichen
  AI_S <- data.frame(GME_I_S3, Assumption = "Annahme 1", mean = mean(GME_I_S3))
  names(AI_S)<-c("value", "Assumption", "mean")

  AII_S <- data.frame(GME_II_S3, Assumption = "Annahme 2", mean = mean(GME_II_S3))
  names(AII_S)<-c("value", "Assumption", "mean")

  AIII_S <- data.frame(GME_III_S3, Assumption = "Annahme 3", mean = mean(GME_III_S3))
  names(AIII_S)<-c("value", "Assumption", "mean")

  all_A_S <- rbind(AI_S, AII_S, AIII_S)

  all_A_S3_plot <<- ggplot(all_A_S,aes(x = value, fill = Assumption)) +
    geom_density(alpha=0.3) +   theme_bw() +
    stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
    scale_shape_manual(values = c(1, 15, 19)) +
    scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
    scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
    xlab(TeX("$\\Delta_j$")) +
    xlim(-0.005, 0.015)
  ggsave("exp_plot_S3.jpg", width = 7, height = 4)
  return(all_A_S3_plot)
}
set.seed(12226699)
GME_Sim_smlr_lrgr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

################################################################################
# GME berechnen auf dem Datensatz  für Height > split_value und Modell auf Height <= split_value
GME_Sim_lrgr_smlr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
  # Mehrmals Simulieren
  sim_list <- list()
  for(i in seq_along(1:100)) {
    # Regressoren simulieren
    Height = round(rnorm(1000, 162.5, 100))
    BMI = round(rnorm(1000, 27, 28), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    BMI_EA = round(rnorm(1000, 21.5, 12), 2)
    Weight_Gain = (BMI - BMI_EA) * ((Height * 0.01)^2)
    X <- data.frame(Height, BMI, BMI_EA, Weight_Gain)
    # wahre Betas festlegen
    betatrue<-c(beta1, beta2, beta3, beta4, beta5)
    # Zielvariable schätzen
    Y <- apply(X,1,function(x)rbinom(1,1,inv.logit(betatrue[1]+betatrue[2]*x[1]+betatrue[3]*x[2]+betatrue[4]*x[3] + betatrue[5]*x[4])))
    # Datensatz bilden
    sim_list[[i]]<-cbind(Y,X)
  }
  sim_matrix <- abind(sim_list, along=3)
  data_S <- as.data.frame(apply(sim_matrix, c(1,2), mean))
  data_S <- data.frame(Y = round(data_S$Y),
                       Height = round(data_S$Height),
                       BMI = round(data_S$BMI, 2),
                       BMI_EA = round(data_S$BMI_EA, 2),
                       Weight_Gain = round(data_S$Weight_Gain, 2))
  # Teile Datensatz an bestimmter Stelle
  data_S_1 <- data_S[data_S$Height > split_value, ]
  data_S_2 <- data_S[data_S$Height <= split_value, ]

  # Modell schätzen
  modelS<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_2)
  drawsS <- rmvnorm(1000, modelS$coefficients, vcov(modelS))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI größer 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # GME
  Height_I_exp <- seq(min(Height_I), max(Height_I))
  GME_I_S4 <<- apply(drawsS, 1, function(x) mean((1/length(Height_I_exp))*(inv.logit(x[1] +
                                                                                  x[2]*max(Height_I_exp) +
                                                                                  x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                                  x[5]*data_I_S$Weight_Gain) -
                                                                        inv.logit(x[1] +
                                                                                    x[2]*min(Height_I_exp) +
                                                                                    x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                                    x[5]*data_I_S$Weight_Gain))))
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  GME_II_S4 <<- apply(drawsS, 1, function(x) mean((1/length(Height_S))*(inv.logit(x[1] +
                                                                                  x[2]* max(Height_S)+
                                                                                  x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                                  x[5]*data_S_1$Weight_Gain) -
                                                                        inv.logit(x[1] +
                                                                                    x[2]*min(Height_S) +
                                                                                    x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                                    x[5]*data_S_1$Weight_Gain))))
  # Annahme 3
  GME_III_S4 <<- apply(drawsS, 1, function(x) mean(inv.logit.deriv(x[1] +
                                                                    x[2]*data_S_1$Height +
                                                                    x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                    x[5]*data_S_1$Weight_Gain, x[2])))
  # GMEs vergleichen
  AI_S <- data.frame(GME_I_S4, Assumption = "Annahme 1", mean = mean(GME_I_S4))
  names(AI_S)<-c("value", "Assumption", "mean")

  AII_S <- data.frame(GME_II_S4, Assumption = "Annahme 2", mean = mean(GME_II_S4))
  names(AII_S)<-c("value", "Assumption", "mean")

  AIII_S <- data.frame(GME_III_S4, Assumption = "Annahme 3", mean = mean(GME_III_S4))
  names(AIII_S)<-c("value", "Assumption", "mean")

  all_A_S <- rbind(AI_S, AII_S, AIII_S)

  all_A_S4_plot <<- ggplot(all_A_S,aes(x = value, fill = Assumption)) +
    geom_density(alpha=0.3) +   theme_bw() +
    stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
    scale_shape_manual(values = c(1, 15, 19)) +
    scale_fill_manual(values = c("deepskyblue4", "yellowgreen", "darkorchid1")) +
    scale_color_manual(values = c("deepskyblue4" ,"yellowgreen", "darkorchid1")) +
    xlab(TeX("$\\Delta_j$")) +
    xlim(-0.005, 0.015)
  ggsave("exp_plot_S4.jpg", width = 7, height = 4)
  return(all_A_S4_plot)
}
set.seed(12226699)
GME_Sim_lrgr_smlr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

################################################################################
all_A_plot <- ggarrange(all_A_S_plot, all_A_S2_plot, all_A_S3_plot, all_A_S4_plot,
                        common.legend = TRUE, legend = "bottom", labels = "auto")
ggsave("GME_plot_all.jpg", width = 10.5, height = 6)


