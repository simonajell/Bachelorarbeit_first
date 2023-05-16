library(MixAll)

data(HeartDisease.target,HeartDisease.cont,HeartDisease.cat)

HeartDisease<-cbind(HeartDisease.cont,HeartDisease.cat,HeartDisease.target)

HeartDisease$target<-ifelse(HeartDisease$num==0,0,1)
HeartDisease$cp_cat<-as.factor(HeartDisease$cp)
levels(HeartDisease$cp_cat)<-list(asymptomatic=4, typical_angina=1, 
                                  atypical_angina=2, non_anginal_pain=3)
HeartDisease$sex<-as.factor(ifelse(HeartDisease$sex==1,"male","female"))

attach(HeartDisease)
