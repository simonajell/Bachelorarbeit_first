library(glm.predict)


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
# predicts(model, "F;all;F;F", sim.count = 1000)
draws_intercept <- rnorm(1000, -0.875915, 0.207064)
draws_assignmentG <- rnorm(1000, 0.090637, 0.127671)
draws_age <- rnorm(1000, 0.023, 0.003917)
draws_genderM <- rnorm(1000, 0.165444, 0.119941)
draws_ethnicityB <- rnorm(1000, -0.194223, 0.190011)
draws_ethnicityL <- rnorm(1000, -0.649621, 0.144719)
draws_ethnicityOU <- rnorm(1000, -0.241925, 0.289881)
draws <- data.frame(draws_intercept, draws_assignmentG, draws_age, draws_genderM,
                    draws_ethnicityB, draws_ethnicityL, draws_ethnicityOU)


mittlerer_Ewert_W <- apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=1,ethB=0,ethL=0,ethO=0),-"ethW")))/nrow(df)})
mittlerer_Ewert_B<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=1,ethL=0,ethO=0),-"ethW")))/nrow(df)})
mittlerer_Ewert_L<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=1,ethO=0),-"ethW")))/nrow(df)})
mittlerer_Ewert_O<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=0,ethO=1),-"ethW")))/nrow(df)})
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

ggplot(mittlerer_Ewert_race, aes(x = value, y = Ethnicity)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi")+
  scale_x_continuous(labels = scales::percent)+
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta,\\cdot)$"))+
  coord_flip()+theme_bw()+ggtitle("Under assumption (A.II')")+xlim(0.25,0.65)
# "White" sieht anders aus als im Beispiel? der rest ist gleich
################################################################################



################################################################################

# im metrischen Fall: AGE ist interessierender Regressor

model_cont <- glm(LOS > median(LOS) ~ AGE + GENDER , family = "binomial", data = data)
summary(model_cont)

draws_intercept_cont <- rnorm(1000, -0.971574, 0.162381)
draws_age_cont <- rnorm(1000, 0.024245, 0.003842)
draws_gender_cont <- rnorm(1000, 0.178560, 0.118659)
draws_cont <- data.frame(draws_intercept_cont, draws_age_cont, draws_gender_cont)

# mittlerer Erwartungswert für jedes Alter berechnen
age <- seq(18, 87, by = 1)
gender <- c(1,0)
AgePred <- apply(draws_cont,1,function(x)(inv.logit(x[1] + x[2]*age + x[3]*gender)))
AgePred_mean <- apply(AgePred,1, mean)
age_dt <- data.frame(age, "mean" = AgePred_mean)

# Plot für den erwarteten Mittelwert für jedes Alter
ggplot(age_dt, aes(x = age, y = mean)) +
  geom_smooth()

################################################################################
#### Assumption 2?
# mittlere Steigung zwischen 30 und 50 Jahren
dt_cont<-df[which(data$AGE > 30 & data$AGE < 50),]

AgePredII <- apply(draws_cont, 1, function(x) 1/20*((inv.logit(x[1] + x[2]*50 + x[3]*gender) - inv.logit(x[1] + x[2]*30 + x[3]*gender))))
AgePredII_mean <- mean(AgePredII)


UnifII<-apply(draws_cont,1,logistic_unif)
ggplot(data.frame(rating_mean=UnifII),aes(x=rating_mean))+
  geom_density(alpha=0.25,fill=4) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$\\Delta_s(\\theta)$"))




#################################################################################
# Angenommen AGE wäre unser interessierender Regressor für die Vorhersage vom Krankenhausaufenthalt
#### Assumption 3?

# Ableitung unseres Modells:
AgePredIII <- apply(draws_cont, 1, function(x) inv.logit.deriv(x[1] + x[2]*age + x[3]*gender, x[2]))
AgePredIII_mean <- mean(AgePredIII)
# wenn Alter um 1 steigt, dann verändert sich der mittlere Erwartungswert um ca.0.0054

# Plot dazu
UnifIII<-apply(AgePredIII,1,logistic_unif)
ggplot(data.frame(rating_mean=UnifIII),aes(x=rating_mean))+
  geom_density(alpha=0.25,fill=4) +
  xlim(0.0010, 0.002) +
  stat_pointinterval(position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4) +
  xlab(TeX("$\\Delta_s(\\theta)$"))


################
# gemeinsamer Plot zum Vergleich

data.frame(AII = UnifII, AIII = UnifIII)






