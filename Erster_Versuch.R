library(plyr)
library(dplyr)
library(tidyr)
library(tidybayes)
library(ggplot2)
library(ggpubr)
library(ggridges)
library(ggeasy)
library(readxl)
library(brms)
library(fastDummies)
library(HDInterval)
library(ggdist)
library(latex2exp)


utils::download.file(url="https://data.mendeley.com/public-files/datasets/25yjwbphn4/files/05601aea-ad43-4a93-8d28-adf0fa74b3c9/file_downloaded",
                     destfile = "CYPtrialData.xlsx")
data<-read_xlsx("CYPtrialData.xlsx")

data<-data[complete.cases(data),]
data<-data[which(data$LOS>24*3),]

data$Assignment<-factor(data$Assignment,levels=c("S","G"))
data$'RACE/ETHNICITY'<-factor(data$'RACE/ETHNICITY',levels=c("W","B","L","O/U"))

names(data)[which(names(data)=="RACE/ETHNICITY")]<-"Ethnicity"

## Section 4.1 selbst versuchen
# Modell
blogisticLOS<-brm(as.numeric(LOS>median(LOS))~Assignment+AGE+GENDER+Ethnicity,family=bernoulli(),data=data,seed=2022)
# Werte für Variablen aus Modell gewonnen
set.seed(2022)
logLOS_draws<-blogisticLOS %>%
  spread_draws(b_Intercept , b_AssignmentG ,
               b_AGE, b_GENDERM , b_EthnicityB , b_EthnicityL ,
               b_EthnicityODU ,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()

# Datensatz mit dummy Spalten (binär) erstellen
df<-dplyr::select(dummy_cols(data[,c("AGE","Assignment","GENDER","Ethnicity")],"Ethnicity"),-Ethnicity)
names(df)<-c("age","assignmentG","genderM","ethW","ethB","ethL","ethO") #Spaltennamen ändern
df$genderM<-ifelse(df$genderM=="M",1,0) #Geschlecht binär kodieren
df$assignmentG<-ifelse(df$assignmentG=="G",1,0) #assignment binär kodieren

## A2'
# arithmetisches Mittel über alle Beobachtungen wird berechnet, wenn man interessierenden Ethnizität auf 1 setzt
ethCat1<-list()
ethCat1$ethW<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=1,ethB=0,ethL=0,ethO=0),-"ethW")))/
    nrow(df)})
ethCat1$ethB<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=1,ethL=0,ethO=0),-"ethW")))/
    nrow(df)})
ethCat1$ethL<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=1,ethO=0),-"ethW")))/
    nrow(df)})
ethCat1$ethO<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=0,ethO=1),-"ethW")))/
    nrow(df)})

## A2"
# arithmetisches Mittel über alle Beobachtungen, die gewünschte Ethnizität haben -> 1000-mal
ethCat2<-list()
ethCat2$ethW<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(df[which(df$ethW==1),],-"ethW")))/length(which(df$ethW==1))})
ethCat2$ethB<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(df[which(df$ethB==1),],-"ethW")))/length(which(df$ethB==1))})
ethCat2$ethL<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(df[which(df$ethL==1),],-"ethW")))/length(which(df$ethL==1))})
ethCat2$ethO<-apply(logLOS_draws,1,function(x){sum(intval(x,dplyr::select(df[which(df$ethO==1),],-"ethW")))/length(which(df$ethO==1))})

# Dataframe der für jede Ethnizität die eben berechneten Werte angibt und den Mittelwert der Ethnizität
Pcat1<-ldply(ethCat1, data.frame) %>% mutate(.id=as.factor(.id))
names(Pcat1)<-c("Ethnicity","value") #Umbenennen der Spalten
Pcat1<-merge(Pcat1,
             Pcat1 %>%
               group_by(Ethnicity)%>%
               dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
             by="Ethnicity") #Eine Spalte für den Mittelwert der ganzen Ethnizität hinzufügen
levels(Pcat1$Ethnicity)<- list(White = "ethW", Latinx = "ethL", Black = "ethB", Other = "ethO") #Umbenennen der Ethnizitäten
# Das selbe für ethCat2
Pcat2<-ldply(ethCat2, data.frame) %>% mutate(.id=as.factor(.id))
names(Pcat2)<-c("Ethnicity","value")
Pcat2<-merge(Pcat2,
             Pcat2 %>%
               group_by(Ethnicity)%>%
               dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
             by="Ethnicity")
levels(Pcat2$Ethnicity)<- list(White = "ethW", Latinx = "ethL", Black = "ethB", Other = "ethO")

# Plots
cat_p1<-ggplot(Pcat1, aes(x = value, y = Ethnicity)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi") +
  scale_x_continuous(labels = scales::percent) +
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta ,\\cdot)$")) +
  coord_flip() +
  theme_bw() +
  ggtitle("Under assumption (A.II’)")+
  xlim(0.25,0.65)

cat_p2<-ggplot(Pcat2, aes(x = value, y = Ethnicity)) +
  stat_halfeye(alpha=0.75,point_interval = "mean_hdi") +
  scale_x_continuous(labels = scales::percent) +
  xlab(TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\theta ,\\cdot)$")) +
  coord_flip() +
  theme_bw() +
  ggtitle("Under assumption (A.II’’)") +
  xlim(0.25,0.65)

p1<-ggarrange(cat_p1,cat_p2, nrow=1,common.legend = TRUE,legend="bottom")

# Gemeinsamer Plot
# selbes Prinzip wie vorher (für die einzelnen Plots)
empEth1<-catAII1(logLOS_draws,data)
empEth2<-catAII2(logLOS_draws,data)

# Dataframe der für jede Ethnizität die eben berechneten Werte angibt und den Mittelwert der Ethnizität
Effcat1<-ldply(empEth1, data.frame) %>% mutate(.id=as.factor(.id))
names(Effcat1)<-c("Ethnicity","value")
Effcat1<-merge(Effcat1,
               Effcat1 %>%
                 group_by(Ethnicity)%>%
                 dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
               by="Ethnicity")
Effcat2<-ldply(empEth2, data.frame) %>% mutate(.id=as.factor(.id))
names(Effcat2)<-c("Ethnicity","value")
Effcat2<-merge(Effcat2,
               Effcat2 %>%
                 group_by(Ethnicity)%>%
                 dplyr::summarize_at("value",mean) %>% dplyr::rename(mean=value),
               by="Ethnicity")

# Zusammenfügen beider Effcat's
Effcat<-rbind(cbind(data.frame(assumption="A.II’"),Effcat1),
              cbind(data.frame(assumption="A.II’’"),Effcat2))
levels(Effcat$Ethnicity)<- list(White = "ethW", Latinx = "ethL", Black = "ethB", Other = "ethO")

# PLot
p2<-ggplot(Effcat,aes(x=value,fill=assumption))+
  geom_density(alpha=0.3) +
  theme_minimal()+
  stat_pointinterval(aes(color=assumption,shape=assumption),position = position_dodge(width = 3, preserve = "single"),
                     point_interval = "mean_hdi",point_size=4)+ easy_remove_y_axis() +
  facet_grid(Ethnicity~.) +
  scale_shape_manual(values=c(15,19)) +
  scale_fill_manual(values=c("green4","navyblue"))+scale_color_manual(values=c("green4","navyblue")) +
  theme(legend.position = "bottom") +
  xlab(TeX("$\\Delta_s (\\theta)$"))



