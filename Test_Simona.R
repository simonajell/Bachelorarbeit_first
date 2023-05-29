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

# Datensatz laden
utils::download.file(url="https://data.mendeley.com/public-files/datasets/25yjwbphn4/files/05601aea-ad43-4a93-8d28-adf0fa74b3c9/file_downloaded",
                     destfile = "CYPtrialData.xlsx")
data<-read_xlsx("CYPtrialData.xlsx")
# Daten filtern
data<-data[complete.cases(data),]
data<-data[which(data$LOS>24*3),] # nur Hospitalisierungen länger als 3 Tage
data$Assignment<-factor(data$Assignment,levels=c("S","G")) #Therapiegruppen in Faktor umwandeln
data$`RACE/ETHNICITY`<-factor(data$`RACE/ETHNICITY`,levels=c("W","B","L","O/U")) #Ethnizität in Faktor umwandeln
names(data)[which(names(data)=="# Psychotropic Medications")]<-"Psychotropic.Medications" # Variable umbenennen
names(data)[which(names(data)=="# Administrations")]<-"Administrations" # Variable umbenennen
names(data)[which(names(data)=="RACE/ETHNICITY")]<-"Ethnicity" # Variable umbenennen

#kleinerer Datensatz
df<-dplyr::select(dummy_cols(data[,c("AGE","Assignment","GENDER","Ethnicity")],"Ethnicity"),-Ethnicity)
names(df)<-c("age","assignmentG","genderM","ethW","ethB","ethL","ethO")
df$genderM<-ifelse(df$genderM=="M",1,0)
df$assignmentG<-ifelse(df$assignmentG=="G",1,0)

################################################################################
# Beispiel 4.1 mit frequentistischem Modell
# Modell fitten
model <- glm(LOS > median(LOS) ~ Assignment + AGE + GENDER + Ethnicity, family = "binomial", data = data)
summary(model)

# Draws:
draws <- rmvnorm(1000, model$coefficients, vcov(model))


mittlerer_Ewert_W <- apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=1,ethB=0,ethL=0,ethO=0),-"ethW")))/nrow(df)})
mittlerer_Ewert_B<-apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=1,ethL=0,ethO=0),-"ethW")))/nrow(df)})
mittlerer_Ewert_L<-apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=1,ethO=0),-"ethW")))/nrow(df)})
mittlerer_Ewert_O<-apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=0,ethO=1),-"ethW")))/nrow(df)})
mittlerer_Ewert_race <- data.frame()
mittlerer_Ewert_race <- data.frame("White" = mittlerer_Ewert_W, "Latinx" = mittlerer_Ewert_L, "Black" = mittlerer_Ewert_B,
                                   "Other" = mittlerer_Ewert_O)
mittlerer_Ewert_race<-ldply(mittlerer_Ewert_race, data.frame) %>% mutate(.id=as.factor(.id))
names(mittlerer_Ewert_race)<-c("Ethnicity","value")
mittlerer_Ewert_race <- merge(mittlerer_Ewert_race,
      mittlerer_Ewert_race %>%
        group_by(Ethnicity)%>%
        dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
      by="Ethnicity")
levels(mittlerer_Ewert_race$Ethnicity)<- list(White  = "ethW", Latinx  = "ethL", Black  = "ethB", Other  = "ethO")

cat_2 <- ggplot(mittlerer_Ewert_race, aes(x = value, y = Ethnicity)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.II')")

################################################################################



################################################################################

# im metrischen Fall: AGE ist interessierender Regressor

model_cont <- glm(LOS > median(LOS) ~ AGE + GENDER + Assignment , family = "binomial", data = data)
summary(model_cont)

# Draws:
draws_cont <- rmvnorm(1000, model_cont$coefficients, vcov(model_cont))

# mittlerer Erwartungswert für jedes Alter berechnen
intval<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$age + betas[3]*regs$genderM + betas[4]*regs$assignmentG
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}

age <- seq(18, 87, by = 1)
gender <- c(1,0)
assignment <- c(1,0)
AgePred <- apply(draws_cont,1,function(x)(inv.logit(x[1] + x[2]*age + x[3]*gender + x[4]*assignment)))
AgePred <- apply(draws_cont,1,function(x)(inv.logit(x[1] + x[2]*age + x[3]*gender + x[4]*assignment)))
AgePred_mean <- apply(AgePred,1, mean)
age_dt <- data.frame(age, "mean" = AgePred_mean)

# Plot für den erwarteten Mittelwert für jedes Alter
ggplot(age_dt, aes(x = age, y = mean)) +
  geom_line()


################################################################################
#### Assumption 2
AgePredII <- apply(draws_cont, 1, function(x) mean(1/69*((inv.logit(x[1] + x[2]*87 + x[3]*df$genderM + x[4]*df$assignmentG) -
                                                       inv.logit(x[1] + x[2]*18 + x[3]*df$genderM + x[4]*df$assignmentG)))))
AgePredII_mean <- mean(AgePredII) # 0.00549 ( durchschnittliche Steigung)

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredII), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))

# Expectation Plot
AgePredII_exp <- matrix(nrow = 70, ncol = 1000)
for(i in seq_along(age)) {
  AgePredII_exp[i, ] <- apply(draws_cont, 1, function(x) mean(inv.logit(x[1] + x[2]*age[i] + x[3]*df$genderM + x[4]*df$assignmentG)))
}



dt_exp_II <- data.frame(x = age, mean = apply(AgePredII_exp, 1, mean),
                         ymin=predict(loess(apply(AgePredII_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~age),data.frame(age=age)),
                         ymax=predict(loess(apply(AgePredII_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~age),data.frame(age=age))
)
pred_plotII <- ggplot(dt_exp_II) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.05) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Alter")+
  theme_bw()


#################################################################################
# Angenommen AGE wäre unser interessierender Regressor für die Vorhersage vom Krankenhausaufenthalt
#### Assumption 3

# Ableitung unseres Modells:
AgePredIII <- apply(draws_cont, 1, function(x) mean(inv.logit.deriv(x[1] + x[2]*df$age + x[3]*df$genderM + x[4]*df$assignmentG, x[2])))
AgePredIII_mean <- mean(AgePredIII) # 0.00572 ( durchschnittliche Steigung)
# wenn Alter um 1 steigt, dann verändert sich der mittlere Erwartungswert um ca.0.00553

# generalized marginal effect - durchschnittliche Steigung
ggplot(data.frame(rating_mean = AgePredIII), aes(x = rating_mean))+
  geom_density(alpha = 0.25,fill = 4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$\\Delta_s (\\theta)$"))


# Expectation Plot
AgePredIII_exp <- matrix(nrow = 70, ncol = 1000)
for(i in seq_along(age)) {
  df_exp <- df[which(df$age == age[i]), ]
  AgePredIII_exp[i, ] <- apply(draws_cont, 1, function(x) mean((inv.logit(x[1] + x[2]*df_exp$age + x[3]*df_exp$genderM + x[4]*df_exp$assignmentG))))
}
dt_exp_III <- data.frame(x = age, mean = apply(AgePredIII_exp, 1, mean),
                         ymin=predict(loess(apply(AgePredIII_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~age),data.frame(age=age)),
                         ymax=predict(loess(apply(AgePredIII_exp,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~age),data.frame(age=age))
                         )
pred_plotIII <- ggplot(dt_exp_III) +
  geom_line(aes(x = x,y = mean), linewidth = 1) +
  geom_ribbon(aes(x = x, ymin = ymin, ymax = ymax), alpha = 0.05) +
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Alter")+
  theme_bw()



#################################################################################
# durchschnittliche Steigungen vergleichen
As2 <- data.frame(AgePredII, Assumption = "Assumption 2", mean = mean(AgePredII))
names(As2)<-c("value", "Assumption", "mean")

As3 <- data.frame(AgePredIII, Assumption = "Assumption 3", mean = mean(AgePredIII))
names(As3)<-c("value", "Assumption", "mean")

# Assuption 2 & 3
both_As <- rbind(As2, As3)

both_cont <- ggplot(both_As,aes(x = value, fill = Assumption)) +
  geom_density(alpha=0.3) + theme_minimal() +
  stat_pointinterval(aes(color = Assumption, shape = Assumption),position = position_dodge(width = 3, preserve = "single"),point_interval = "mean_hdi",point_size=4)+
  scale_shape_manual(values = c(15, 19)) +
  scale_fill_manual(values = c("yellowgreen", "mediumpurple4")) +
  scale_color_manual(values = c("yellowgreen", "mediumpurple4")) +
  theme(legend.position = "bottom") + xlab(TeX("$\\Delta_s (\\theta)$"))


#################################################################################
# nur zwei Regressoren
model_2r <- glm(LOS > median(LOS) ~ AGE + GENDER, family = "binomial", data = data)
summary(model_2r)

# Draws:
draws_2r <- rmvnorm(1000, model_2r$coefficients, vcov(model_2r))

# Assumption 2 auf AGE
int_leng_2r <- max(data$AGE) - min(data$AGE)
AgePredII_2r <- apply(draws_2r, 1, function(x) mean((1/69)*
                 (inv.logit(x[1] + x[2]*87 +x[3]*df$genderM) -
                    inv.logit(x[1] + x[2]*18 + x[3]*df$genderM))))
AgePredII_2r_mean <- mean(AgePredII_2r)

# Assumption 3 auf AGE
AgePredIII_2r <- apply(draws_2r, 1, function(x) mean(inv.logit.deriv(x[1] +
                             x[2]*df$age + x[3]*df$genderM, x[2])))
AgePredIII_2r_mean <- mean(AgePredIII_2r)

# nur zwei Regressoren diskret
model_2rr <- glm(LOS > median(LOS) ~ GENDER + AGE, family = "binomial", data = data)
summary(model_2rr)

# Draws:
draws_2rr <- rmvnorm(1000, model_2rr$coefficients, vcov(model_2rr))

# Assumption 2 auf AGE
intval_II_2rr<-function(betas,regs,deriv=NULL){
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$genderM+betas[3]*regs$age
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}
AII_2rr<-apply(draws_2rr,1,function(x){sum(intval_II_2rr(x,mutate(df,genderM=1))-intval_II_2rr(x,mutate(df,genderM=0)))/nrow(df)})
AII_2rr_mean <- mean(AII_2rr)

# Assumption 3 auf AGE
AIII_2rr <- list()
AIII_2rr<-apply(draws_2rr,1,function(x){(sum(intval_II_2rr(x,df[which(df$genderM==1),]))/length(which(df$genderM==1)))-
    (sum(intval_II_2rr(x,df[which(df$genderM==0),]))/length(which(df$genderM==0)))})
AIII_2rr_mean <- mean(AIII_2rr)




