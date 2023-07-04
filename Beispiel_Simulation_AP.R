
AP_Sim <- function(beta1, beta2, beta3, beta4, beta5) {
  # Mehrmals Simulieren
  sim_list <- list()
  for(i in seq_along(1:100)) {
    # Regressoren simulieren
    X <- data.frame(Height = round(rnorm(1000, 1.625, 0.8), 2),
                    BMI = round(rnorm(1000, 27, 28), 2),
                    BMI_EA = round(rnorm(1000, 21.5, 12), 2),
                    Weight_Gain = round(rnorm(1000, 15, 90), 2) )
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
                       Height = round(data_S$Height, 2),
                       BMI = round(data_S$BMI, 2),
                       BMI_EA = round(data_S$BMI_EA, 2),
                       Weight_Gain = round(data_S$Weight_Gain, 2))
  # Modell schätzen
  modelS<-glm(Y ~ Height + BMI + BMI_EA + Weight_Gain, family=binomial, data=data_S)
  drawsS <- rmvnorm(1000, modelS$coefficients, vcov(modelS))

  # Annahme 1
  # Subgruppe: Frauen, die größer als 1.70m sind und einen BMI größer 25 haben.
  # Der BMI im frühen Erwachsenenalter und die Gewichtszunahme sind zufällig
  # kleinerer Datensatz, da der kombinierte Datensatz zu groß ist.
  Height_I <- round(runif(50, 1.70, 1.85), 2)
  BMI_I <- round(runif(50, 25, 35), 2)
  BMI_EA_I = round(rnorm(50, 19, 25), 2)
  Weight_Gain_I = round(rnorm(50, -15, 51), 2)
  data_I_S <- expand.grid(Height = Height_I, BMI = BMI_I, BMI_EA = BMI_EA_I, Weight_Gain = Weight_Gain_I)

  AgePredI_S_exp <- matrix(nrow = 49, ncol = 1000)
  Height_I_exp <- seq(min(Height_I), max(Height_I), by = 0.01)
  for(i in seq_along(Height_I_exp)) {
    AgePredI_S_exp[i, ] <-apply(drawsS, 1, function(x) mean(inv.logit(x[1] +
                                                                        x[2]*Height_S[i] +
                                                                        x[3]*data_I_S$BMI + x[4]*data_I_S$BMI_EA +
                                                                        x[5]*data_I_S$Weight_Change)))
  }
  dtI_S_exp <- data.frame(x = Height_S, mean = apply(AgePredI_S_exp, 1, mean),
                          ymin=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                          ymax=predict(loess(apply(AgePredI_S_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2]) ~ Height_I_exp),data.frame(Height_S = Height_I_exp)),
                          Assumption = "Annahme 1"
  )

  # Annahme 2
  AgePredII_S_exp <- matrix(nrow = 49, ncol = 1000)
  for(i in seq_along(Height_S)) {
    AgePredII_S_exp[i, ] <- apply(drawsS, 1, function(x) mean(inv.logit(x[1] +
                                                                          x[2]*Height_S[i] +
                                                                          x[3]*data_S$BMI + x[4]*data_S$BMI_EA +
                                                                          x[5]*data_S$Weight_Gain)))
  }
  dtII_S_exp <- data.frame(x = Height_S, mean = apply(AgePredII_S_exp, 1, mean),
                           ymin=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                           ymax=predict(loess(apply(AgePredII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                           Assumption = "Annahme 2"
  )

  # Annahme 3
  AgePredIII_S_exp <- matrix(nrow = 49, ncol = 1000)
  for(i in seq_along(Height_S)) {
    df_exp <- data_S[which(data_S$Height == Height_S[i]), ]
    AgePredIII_S_exp[i, ] <- apply(drawsS, 1, function(x) mean((inv.logit(x[1] +
                                                                            x[2]*df_exp$Height +
                                                                            x[3]*df_exp$BMI + x[4]*df_exp$BMI_EA +
                                                                            x[5]*df_exp$Weight_Gain))))
  }
  dtIII_S_exp <- data.frame(x = Height_S, mean = apply(AgePredIII_S_exp, 1, mean),
                            ymin=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[1]) ~ Height_S), data.frame(Height_S = Height_S)),
                            ymax=predict(loess(apply(AgePredIII_S_exp,1,function(x)HDInterval::hdi(x, credMass = 0.95)[2]) ~ Height_S), data.frame(Height_S = Height_S)),
                            Assumption = "Annahme 3"
  )
  all_A_exp <- rbind(dtI_S_exp, dtII_S_exp, dtIII_S_exp)
  pred_plot_all_S <- ggplot(all_A_exp) +
    geom_line(aes(x = x,y = mean, color = Assumption), linewidth = 1) +
    geom_point(aes(x = x,y = mean, color = Assumption), size = 0.7) +
    geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax, fill = Assumption), alpha = 0.15) +
    scale_fill_manual(values = c("firebrick2" ,"orange", "steelblue1")) +
    scale_color_manual(values = c("firebrick2" ,"orange", "steelblue2")) +
    labs(y="Wert der Adjusted Predictions",x="Cholesterin")+
    theme_bw()
  ggsave("exp_plot_1.jpg", width = 7, height = 4)
}
AP_Sim(-0.8, 3.5, -0.1, -0.15, 0.015)
