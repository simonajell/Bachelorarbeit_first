library(boot)
library(abind)
library(tibble)
library(mvtnorm)
library(ggplot2)
library(tidybayes)
library(latex2exp)

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

AP_Sim_smlr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
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
    #betatrue<-c(beta1, beta2, beta3, beta4, beta5)
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
  # Kleinere Simulation, da der kombinierte Datensatz sonst zu groß ist.
  Height_I_exp <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I_exp, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # AP
  # Annahme 1
  AgePredI_S_exp <- matrix(nrow = length(Height_I_exp), ncol = 1000)
  for(i in seq_along(Height_I_exp)) {
    AgePredI_S_exp[i, ] <-apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                        x[2]*Height_I_exp[i] +
                                                                        x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                        x[5]*data_I_S$Weight_Gain)))
  }

  dtI_S_exp <<- data.frame(x = Height_I_exp, mean = apply(AgePredI_S_exp, 1, mean),
                          ymin=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                          ymax=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                          Assumption = "Annahme 1"
  )
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  AgePredII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                          x[2]*Height_S[i] +
                                                                          x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                          x[5]*data_S_1$Weight_Gain)))
  }
  dtII_S_exp <<- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                           ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                           ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                           Assumption = "Annahme 2"
  )
  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S_1[which(data_S_1$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean((inv.logit(x[1] +
                                                                             x[2]*df_exp$Height +
                                                                             x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                             x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S_exp <<- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                             ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                             ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                             Assumption = "Annahme 3"
  )
  # APs vergleichen
  all_A_exp <- rbind(dtI_S_exp, dtII_S_exp, dtIII_S_exp)
  pred_plot_all_S <<- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Körpergröße")+
    theme_bw()
  ggsave("ap_plot_S1.jpg", width = 7, height = 4)
  return(pred_plot_all_S)
}
set.seed(12226699)
AP_Sim_smlr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

AP_Sim_smlr_test <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
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

  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  AgePredII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                           x[2]*Height_S[i] +
                                                                           x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                           x[5]*data_S_1$Weight_Gain)))
  }
  dtII_S_exp <<- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                            ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                            ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                            Assumption = "Annahme 2"
  )
  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S_1[which(data_S_1$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean((inv.logit(x[1] +
                                                                             x[2]*df_exp$Height +
                                                                             x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                             x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S_exp <<- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                             ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                             ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                             Assumption = "Annahme 3"
  )
  # APs vergleichen
  all_A_exp <- rbind(dtII_S_exp, dtIII_S_exp)
  pred_plot_all_S <<- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Körpergröße")+
    theme_bw()
  return(pred_plot_all_S)
}
set.seed(12226699)
AP_Sim_smlr_test(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

################################################################################
AP_Sim_lrgr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
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

  # Modell schätzen
  modelS1<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_1)
  drawsS1 <- rmvnorm(1000, modelS1$coefficients, vcov(modelS1))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI größer 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I_exp <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I_exp, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # AP
  # Annahme 1
  AgePredI_S_exp <- matrix(nrow = length(Height_I_exp), ncol = 1000)
  for(i in seq_along(Height_I_exp)) {
    AgePredI_S_exp[i, ] <-apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                         x[2]*Height_I_exp[i] +
                                                                         x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                         x[5]*data_I_S$Weight_Gain)))
  }

  dtI_S2_exp <<- data.frame(x = Height_I_exp, mean = apply(AgePredI_S_exp, 1, mean),
                           ymin=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                           ymax=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                           Assumption = "Annahme 1"
  )
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  AgePredII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                           x[2]*Height_S[i] +
                                                                           x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                           x[5]*data_S_1$Weight_Gain)))
  }
  dtII_S2_exp <<- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                            ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                            ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                            Assumption = "Annahme 2"
  )
  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S_1[which(data_S_1$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean((inv.logit(x[1] +
                                                                             x[2]*df_exp$Height +
                                                                             x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                             x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S2_exp <<- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                             ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                             ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                             Assumption = "Annahme 3"
  )
  # APs vergleichen
  all_A_exp <- rbind(dtI_S2_exp, dtII_S2_exp, dtIII_S2_exp)
  pred_plot_all_S2 <<- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Körpergröße")+
    theme_bw()
  ggsave("ap_plot_S2.jpg", width = 7, height = 4)
  return(pred_plot_all_S2)
}
set.seed(12226699)
AP_Sim_lrgr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)

AP_Sim_test <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
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

  # Modell schätzen
  modelS1<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_1)
  drawsS1 <- rmvnorm(1000, modelS1$coefficients, vcov(modelS1))


  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  AgePredII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                           x[2]*Height_S[i] +
                                                                           x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                           x[5]*data_S_1$Weight_Gain)))
  }
  dtII_S2_exp <<- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                             ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                             ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                             Assumption = "Annahme 2"
  )
  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S_1[which(data_S_1$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean((inv.logit(x[1] +
                                                                             x[2]*df_exp$Height +
                                                                             x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                             x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S2_exp <<- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                              ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                              ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                              Assumption = "Annahme 3"
  )
  # APs vergleichen
  all_A_exp <- rbind(dtII_S2_exp, dtIII_S2_exp)
  pred_plot_all_S2 <<- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Körpergröße")+
    theme_bw()
  return(pred_plot_all_S2)
}
set.seed(12226699)
AP_Sim_test(-9.5, 0.009, 0.3, -0.3, 0.3, 165)

################################################################################
AP_Sim_smlr_lrgr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
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
  modelS1<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_2)
  drawsS1 <- rmvnorm(1000, modelS1$coefficients, vcov(modelS1))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI größer 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I_exp <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I_exp, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # AP
  # Annahme 1
  AgePredI_S_exp <- matrix(nrow = length(Height_I_exp), ncol = 1000)
  for(i in seq_along(Height_I_exp)) {
    AgePredI_S_exp[i, ] <-apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                         x[2]*Height_I_exp[i] +
                                                                         x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                         x[5]*data_I_S$Weight_Gain)))
  }

  dtI_S3_exp <<- data.frame(x = Height_I_exp, mean = apply(AgePredI_S_exp, 1, mean),
                            ymin=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                            ymax=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                            Assumption = "Annahme 1"
  )
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  AgePredII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                           x[2]*Height_S[i] +
                                                                           x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                           x[5]*data_S_1$Weight_Gain)))
  }
  dtII_S3_exp <<- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                             ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                             ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                             Assumption = "Annahme 2"
  )
  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S_1[which(data_S_1$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean((inv.logit(x[1] +
                                                                             x[2]*df_exp$Height +
                                                                             x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                             x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S3_exp <<- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                              ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                              ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                              Assumption = "Annahme 3"
  )
  # APs vergleichen
  all_A_exp <- rbind(dtI_S3_exp, dtII_S3_exp, dtIII_S3_exp)
  pred_plot_all_S3 <<- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Körpergröße")+
    theme_bw()
  ggsave("ap_plot_S3.jpg", width = 7, height = 4)
  return(pred_plot_all_S3)
}
set.seed(12226699)
AP_Sim_smlr_lrgr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)


################################################################################
AP_Sim_lrgr_smlr <- function(beta1, beta2, beta3, beta4, beta5, split_value) {
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
  modelS1<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S_2)
  drawsS1 <- rmvnorm(1000, modelS1$coefficients, vcov(modelS1))

  # Annahme 1
  # Subgruppe: Frauen, die einen BMI größer 25 haben, und im frühen Erwachsenenalter
  # einen BMI kleiner 21 hatten.
  # Die Gewichtszunahme ist zufällig.
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I_exp <- seq(min(data_S_1$Height), max(data_S_1$Height), by = 1)
  BMI_I <- round(runif(10, 20, 25), 2)
  BMI_EA_I = round(rnorm(10, 19, 21), 2)
  Weight_Gain_I = round(rnorm(30, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I_exp, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)
  # AP
  # Annahme 1
  AgePredI_S_exp <- matrix(nrow = length(Height_I_exp), ncol = 1000)
  for(i in seq_along(Height_I_exp)) {
    AgePredI_S_exp[i, ] <-apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                         x[2]*Height_I_exp[i] +
                                                                         x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                         x[5]*data_I_S$Weight_Gain)))
  }

  dtI_S4_exp <<- data.frame(x = Height_I_exp, mean = apply(AgePredI_S_exp, 1, mean),
                            ymin=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                            ymax=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                            Assumption = "Annahme 1"
  )
  # Annahme 2
  Height_S <- seq(min(data_S_1$Height), max(data_S_1$Height))
  AgePredII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean(inv.logit(x[1] +
                                                                           x[2]*Height_S[i] +
                                                                           x[3]*data_S_1$BMI + x[4]*data_S_1$BMI_EA +
                                                                           x[5]*data_S_1$Weight_Gain)))
  }
  dtII_S4_exp <<- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                             ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                             ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                             Assumption = "Annahme 2"
  )
  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = length(Height_S), ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S_1[which(data_S_1$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS1, 1, function(x) mean((inv.logit(x[1] +
                                                                             x[2]*df_exp$Height +
                                                                             x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                             x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S4_exp <<- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                              ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                              ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                              Assumption = "Annahme 3"
  )
  # APs vergleichen
  all_A_exp <- rbind(dtI_S4_exp, dtII_S4_exp, dtIII_S4_exp)
  pred_plot_all_S4 <<- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Körpergröße")+
    theme_bw()
  ggsave("ap_plot_S4.jpg", width = 7, height = 4)
  return(pred_plot_all_S4)
}
set.seed(12226699)
AP_Sim_lrgr_smlr(-9.5, 0.008, 0.3, -0.3, 0.3, 165)


