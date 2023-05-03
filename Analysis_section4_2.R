library(readr)
library(ggplot2)
library(ggpubr)
library(latex2exp)
library(brms)
library(tidybayes)
library(HDInterval)
library(gridExtra)
library(grid)
library(mvtnorm)
library(ggeasy)
library(ggdist)
library(dplyr)
library(tidyr)
library(ggrepel)
library(MASS)
library(ggridges)
library(rgl)
library(viridis) 
library(zoo)
library(scales)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))#sets working directory to file location
source("helper.R")

utils::download.file(url="https://osf.io/download/fv8c3/",destfile = "SilberzahnData.csv")
SilberzahnData <- read_csv("SilberzahnData.csv")

SilberzahnData$rating_mean<-rowMeans(cbind(SilberzahnData$rater1,SilberzahnData$rater2),na.rm=TRUE)
SilberzahnData$redCard_rate<-SilberzahnData$redCards/SilberzahnData$games

lm<-brm(redCard_rate~rating_mean,data=SilberzahnData,seed=54321)
glm_bin<-brm(redCards|trials(games)~rating_mean,data=SilberzahnData,family=binomial(),seed=54321)
glm_poi<-brm(redCards~rating_mean+offset(log(games)),family = poisson(),data=SilberzahnData,seed=54321)



set.seed(600)
LinearBetas<-lm %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(LinearBetas)<-substring(names(LinearBetas), 3)

set.seed(600)
PoissonBetas<- glm_poi %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(PoissonBetas)<-substring(names(PoissonBetas), 3)


set.seed(600)
LogisticBetas<- glm_bin %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(LogisticBetas)<-substring(names(LogisticBetas), 3)


rating<-seq(0,1,by=0.01)

LinPred<-apply(LinearBetas,1,function(x)(x[1]+x[2]*rating))
BinPred<-apply(LogisticBetas,1,function(x)(inv.logit(x[1]+x[2]*rating)))
PoiPred<-apply(PoissonBetas,1,function(x)(exp(x[1]+x[2]*rating)))

Pred<-rbind(data.frame(Distribution="Normal (linear)",
                       x=rating,
                       mean=apply(LinPred,1,mean),
                       ymin=predict(loess(apply(LinPred,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~rating),data.frame(rating=rating)),
                       ymax=predict(loess(apply(LinPred,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~rating),data.frame(rating=rating))
),
data.frame(Distribution="Poisson (exponential)",
           x=rating,
           mean=apply(PoiPred,1,mean),
           ymin=predict(loess(apply(PoiPred,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~rating),data.frame(rating=rating)),
           ymax=predict(loess(apply(PoiPred,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~rating),data.frame(rating=rating))
),
data.frame(Distribution="Binomial (logistic)",
           x=rating,
           mean=apply(BinPred,1,mean),
           ymin=predict(loess(apply(BinPred,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~rating),data.frame(rating=rating)),
           ymax=predict(loess(apply(BinPred,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~rating),data.frame(rating=rating))
))


pred_plot<-ggplot(Pred)+geom_line(aes(x=x,y=mean,color=Distribution,linetype=Distribution),size=1)+
  geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax,fill=Distribution),alpha=0.05)+
  scale_colour_manual(values=c("navyblue","red3","green3"),name="") + 
  scale_linetype_manual(values=c("solid","dotted","dashed"),name="")+
  scale_fill_manual(values=c("navyblue","red3","green3"),name="")+
  labs(y=TeX("$g^{\\left[ I\\right]}_{avg}\\,(\\hat{\\theta},\\cdot)$"),x="Average skin tone rating",color="Distribution Assumption")+
  theme_bw()



LinSlope<-apply(LinearBetas,1,function(x)(rep(x[2],length(rating))))
BinSlope<-apply(LogisticBetas,1,function(x)(inv.logit.deriv(x[1]+x[2]*rating,x[2])))
PoiSlope<-apply(PoissonBetas,1,function(x)(x[2]*exp(x[1]+x[2]*rating)))

Slope<-rbind(data.frame(Distribution="Normal (linear)",
                        x=rating,
                        mean=apply(LinSlope,1,mean),
                        ymin=predict(loess(apply(LinSlope,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~rating),data.frame(rating=rating)),
                        ymax=predict(loess(apply(LinSlope,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~rating),data.frame(rating=rating))
),
data.frame(Distribution="Poisson (exponential)",
           x=rating,
           mean=apply(PoiSlope,1,mean),
           ymin=predict(loess(apply(PoiSlope,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~rating),data.frame(rating=rating)),
           ymax=predict(loess(apply(PoiSlope,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~rating),data.frame(rating=rating))
),
data.frame(Distribution="Binomial (logistic)",
           x=rating,
           mean=apply(BinSlope,1,mean),
           ymin=predict(loess(apply(BinSlope,1,function(x)HDInterval::hdi(x,credMass = 0.95)[1])~rating),data.frame(rating=rating)),
           ymax=predict(loess(apply(BinSlope,1,function(x)HDInterval::hdi(x,credMass = 0.95)[2])~rating),data.frame(rating=rating))
))


slope_plot<-ggplot(Slope)+geom_line(aes(x=x,y=mean,color=Distribution,linetype=Distribution),size=1)+
  #geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax,fill=Distribution),alpha=0.05)+
  scale_colour_manual(values=c("navyblue","red3","green3"),name="") + 
  scale_linetype_manual(values=c("solid","dotted","dashed"),name="")+
  scale_fill_manual(values=c("navyblue","red3","green3"),name="")+
  labs(y=TeX("$\\bar{s}\\,(\\hat{\\theta},\\cdot)$"),x="Average skin tone rating",color="Distribution Assumption")+
  theme_bw()

################################################################################
# figure 2
ggarrange(pred_plot,slope_plot, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
################################################################################


set.seed(55)
LinearBetas<-lm %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 4000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(LinearBetas)<-substring(names(LinearBetas), 3)

set.seed(55)
PoissonBetas<- glm_poi %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 4000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(PoissonBetas)<-substring(names(PoissonBetas), 3)


set.seed(55)
LogisticBetas<- glm_bin %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 4000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(LogisticBetas)<-substring(names(LogisticBetas), 3)





LogisticUnif<-apply(LogisticBetas,1,logistic_unif)
PoissonUnif<-apply(PoissonBetas,1,poisson_unif)



up_linear<-ggplot(LinearBetas,aes(x=rating_mean))+
  geom_density(alpha=0.25,fill=4)+
  xlim(-0.001,0.004)+theme_minimal()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_pointinterval(point_interval = "mean_hdi")


up_poisson<-ggplot(data.frame(rating_mean=PoissonUnif),aes(x=rating_mean))+
  geom_density(alpha=0.25,fill=4)+
  xlim(-0.001,0.004)+theme_minimal()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_pointinterval(point_interval = "mean_hdi")

up_logistic<-ggplot(data.frame(rating_mean=LogisticUnif),aes(x=rating_mean))+
  geom_density(alpha=0.25,fill=4)+
  xlim(-0.001,0.004)+theme_minimal()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  stat_pointinterval(point_interval = "mean_hdi")


axis1<-ggplot()+xlim(-0.001,0.004)+theme_bw()+easy_remove_y_axis()+
  xlab(TeX("$\\Delta_s(\\theta)$"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  easy_x_axis_labels_size(20)+easy_x_axis_title_size(22)


eb_poisson<-ggplot(data.frame(x=0,y=exp(summary(glm_poi)$fixed[2,1]),ymax = exp(summary(glm_poi)$fixed[2,4]),ymin = exp(summary(glm_poi)$fixed[2,3])),aes(y=y,x=x))+
  geom_point(size=4)+geom_hline(aes(yintercept=1),color="darkgray",linetype="dashed")+
  geom_errorbar(aes(ymax = ymax, ymin = ymin),size=1)+
  ylim(0.9,1.7)+xlim(-4,4)+theme_minimal()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()+
  annotate("text",y=1,x=-3.5,label="OR=1",hjust=0.5,size=6)

eb_logistic<-ggplot(data.frame(x=0,y=exp(summary(glm_bin)$fixed[2,1]),ymax = exp(summary(glm_bin)$fixed[2,4]),ymin = exp(summary(glm_bin)$fixed[2,3])),aes(y=y,x=x))+
  geom_point(size=4)+geom_hline(aes(yintercept=1),color="darkgray",linetype="dashed")+
  geom_errorbar(aes(ymax = ymax, ymin = ymin),size=1)+
  ylim(0.9,1.7)+xlim(-4,4)+theme_minimal()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()

eb_linear<-ggplot()+easy_remove_axes()+
  ylim(0,2)+xlim(0.9,1.7)+
  geom_vline(aes(xintercept=1),color="darkgray",linetype="dashed")+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  annotate("text",x=1.3,y=1,label="not\n applicable",size=10,color="darkgray")

axis2<-ggplot()+xlim(0.9,1.7)+theme_bw()+easy_remove_y_axis()+
  xlab("Odds Ratio")+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=grid::unit(c(-7,5.5,5.5,5.5), "pt"))+
  easy_x_axis_labels_size(20)+easy_x_axis_title_size(22)


beta_linear<-ggplot(data.frame(x=NA,y=NA),aes(x=x,y=y))+ylim(0,0.75)+theme_bw()+easy_remove_axes()+
  annotate("pointrange", x =  0, y =summary(lm)$fixed[2,1], ymin = summary(lm)$fixed[2,3], ymax = summary(lm)$fixed[2,4],colour = 4, size = 1.5, alpha=0.4)+
  annotate("text", x =  0, y =summary(lm)$fixed[2,1],label=expression(paste(hat(beta),"=")),vjust=-0.45,hjust=-0.001,size=6)+
  annotate("text", x =  0, y =summary(lm)$fixed[2,1],label=round(summary(lm)$fixed[2,1],5),vjust=-1.2,hjust=-0.35,size=6)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()

beta_poisson<-ggplot(data.frame(x=NA,y=NA),aes(x=x,y=y))+ylim(0,0.75)+theme_bw()+easy_remove_axes()+
  annotate("pointrange", x =  0, y =summary(glm_poi)$fixed[2,1], ymin = summary(glm_poi)$fixed[2,3], ymax = summary(glm_poi)$fixed[2,4],colour = 4, size = 1.5, alpha=0.4)+
  annotate("text", x =  0, y =summary(glm_poi)$fixed[2,1],label=expression(paste(hat(beta),"=")),vjust=-0.45,hjust=-0.001,size=6)+
  annotate("text", x =  0, y =summary(glm_poi)$fixed[2,1],label=round(summary(glm_poi)$fixed[2,1],5),vjust=-1.2,hjust=-0.35,size=6)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()

beta_logistic<-ggplot(data.frame(x=NA,y=NA),aes(x=x,y=y))+ylim(0,0.75)+theme_bw()+easy_remove_axes()+
  annotate("pointrange", x =  0, y =summary(glm_bin)$fixed[2,1], ymin = summary(glm_bin)$fixed[2,3], ymax = summary(glm_bin)$fixed[2,4],colour = 4, size = 1.5, alpha=0.4)+
  annotate("text", x =  0, y =summary(glm_bin)$fixed[2,1],label=expression(paste(hat(beta),"=")),vjust=-0.45,hjust=-0.001,size=6)+
  annotate("text", x =  0, y =summary(glm_bin)$fixed[2,1],label=round(summary(glm_bin)$fixed[2,1],5),vjust=-1.2,hjust=-0.35,size=6)+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  coord_flip()


axis3<-ggplot()+xlim(0,0.75)+theme_bw()+easy_remove_y_axis()+
  xlab(TeX("$\\beta-$coefficient"))+
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.margin=grid::unit(c(2,5.5,5.5,5.5), "pt"))+
  easy_x_axis_labels_size(20)+easy_x_axis_title_size(22)


linM<-ggplot()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  annotate("text",x=1,y=1,label="Linear\n model",size=8,color="black")
logisticM<-ggplot()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  annotate("text",x=1,y=1,label="Logistic\n model",size=8,color="black")
poiM<-ggplot()+easy_remove_axes()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())+
  annotate("text",x=1,y=1,label="Poisson\n model",size=8,color="black")


################################################################################
# figure 3
grid.arrange(linM,up_linear,eb_linear,beta_linear,
             logisticM,up_logistic,eb_logistic,beta_logistic,
             poiM,up_poisson,eb_poisson,beta_poisson,ggplot()+
               theme(panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.border = element_blank(),
                     panel.background = element_blank()),axis1,axis2,axis3,nrow=4,ncol=4,heights=c(2,2,2,1),widths=c(1,6,3,3))
################################################################################


set.seed(600)
LinearBetas_sigma<-lm %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               sigma,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(LinearBetas_sigma)[1:2]<-substring(names(LinearBetas_sigma)[1:2], 3)

set.seed(600)
PoissonBetas<- glm_poi %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(PoissonBetas)<-substring(names(PoissonBetas), 3)


set.seed(600)
LogisticBetas<- glm_bin %>%
  spread_draws(b_Intercept,
               b_rating_mean,
               ndraws = 1000) %>%
  dplyr::select(- .chain,-.iteration,-.draw) %>%
  as.data.frame()
names(LogisticBetas)<-substring(names(LogisticBetas), 3)



y<-seq(-0.2,0.2,by=0.001)

linPstart<-data.frame(group=c(rep("Individualized predictive distribution\n for skin color rating of 0",length(y)),
                              rep("Individualized predictive distribution\n for skin color rating of 1",length(y))),
                      x=rep(y,2))

linP<-data.frame()
set.seed(928)
for(i in 1:nrow(LinearBetas_sigma)){
  linP<-rbind(linP,cbind(linPstart,
                         data.frame(y=c(dnorm(y,mean=LinearBetas_sigma$Intercept[i],sd=LinearBetas_sigma$sigma[i]),
                                        dnorm(y,mean=LinearBetas_sigma$Intercept[i]+LinearBetas_sigma$rating_mean[i],sd=LinearBetas_sigma$sigma[i]))),
                         lines=paste0("beta",i)))
  
}

linPlines<-linP %>%
  dplyr::group_by(group,x) %>%
  dplyr::summarise(ymin=mean_hdi(y)$ymin,ymax=mean_hdi(y)$ymax,mean=mean(y))



fy0<-linPlines$mean[which(linPlines$group=="Individualized predictive distribution\n for skin color rating of 0")]
fy1<-linPlines$mean[which(linPlines$group=="Individualized predictive distribution\n for skin color rating of 1")]

linpp<-ggplot(linPlines,aes(x=x,color=group)) +
  geom_ribbon(aes(x=x,ymin=ymin,ymax=ymax,fill=group),alpha=0.4, colour = NA)+
  geom_line(aes(x=x,y=mean,linetype=group))+
  scale_color_manual(name='', values=c("green4","blue4"))+
  scale_fill_manual(name='', values=c("green4","blue4"))+
  theme_bw()+ theme(legend.position="bottom")+xlim(-0.2,0.2)+
  xlab("y")+ylab(TeX("$\\Pi(y,\\theta)$"))+
  geom_vline(xintercept=y[which.max(fy0)],linetype="dashed",color="grey")+
  geom_vline(xintercept=y[which.max(fy1)],linetype="dashed",color="grey")+
  annotate(x=0.0038,y=3,label="Distance between maximums:\n ca. 0.0016",vjust=1,geom="label",size=3)+ 
  labs(color  = "", linetype = "",fill="")




S<-SilberzahnData[,c("redCards","games","redCard_rate","rating_mean")]%>%drop_na()

log0<-sum(apply(LogisticBetas,1,function(x){logistic_efun(beta=x,x=0)}))/nrow(LogisticBetas)
log1<-sum(apply(LogisticBetas,1,function(x){logistic_efun(beta=x,x=1)}))/nrow(LogisticBetas)

poi0<-sum(apply(PoissonBetas,1,function(x){poisson_efun(beta=x,x=0)}))/nrow(PoissonBetas)
poi1<-sum(apply(PoissonBetas,1,function(x){poisson_efun(beta=x,x=1)}))/nrow(PoissonBetas)


nonlinP1<-data.frame(model=c(rep("logistic",6),rep("poisson",6)),sc=factor(rep(c("skin color\n rating of 0","If OR were equal to 1","skin color\n rating of 1"),4),levels=c("skin color\n rating of 0","If OR were equal to 1","skin color\n rating of 1")),
                     x=factor(rep(c(rep("P(Y>0|X)",3),rep("P(Y=0|X)",3)),2),levels=c("P(Y>0|X)",rep("P(Y=0|X)"))),
                     y=c(log0,rollmean(c(log0,log1),2),log1,
                         1-log0,rollmean(c(1-log0,1-log1),2),1-log1,
                         1-dpois(0,poi0),rollmean(c(1-dpois(0,poi0),1-dpois(0,poi1)),2),1-dpois(0,poi1),
                         dpois(0,poi0),rollmean(c(dpois(0,poi0),dpois(0,poi1)),2),dpois(0,poi1)),
                     alpha=rep(c(-1,-1,-1),4))

nonlinP2<-rbind(data.frame(model="logistic",sc="skin color\n rating of 0",x="P(Y>0|X)",y=0,alpha=-1,draws=apply(LogisticBetas,1,function(x){logistic_efun(beta=x,x=0)})),
                data.frame(model="logistic",sc="skin color\n rating of 1",x="P(Y>0|X)",y=0,alpha=-1,draws=apply(LogisticBetas,1,function(x){logistic_efun(beta=x,x=1)})),
                data.frame(model="poisson",sc="skin color\n rating of 0",x="P(Y>0|X)",y=0,alpha=-1,draws=apply(PoissonBetas,1,function(x){poisson_efun(beta=x,x=0)})),
                data.frame(model="poisson",sc="skin color\n rating of 1",x="P(Y>0|X)",y=0,alpha=-1,draws=apply(PoissonBetas,1,function(x){poisson_efun(beta=x,x=1)})),
                data.frame(model="logistic",sc="skin color\n rating of 0",x="P(Y=0|X)",y=0,alpha=-1,draws=1-apply(LogisticBetas,1,function(x){logistic_efun(beta=x,x=0)})),
                data.frame(model="logistic",sc="skin color\n rating of 1",x="P(Y=0|X)",y=0,alpha=-1,draws=1-apply(LogisticBetas,1,function(x){logistic_efun(beta=x,x=1)})),
                data.frame(model="poisson",sc="skin color\n rating of 0",x="P(Y=0|X)",y=0,alpha=-1,draws=1-apply(PoissonBetas,1,function(x){poisson_efun(beta=x,x=0)})),
                data.frame(model="poisson",sc="skin color\n rating of 1",x="P(Y=0|X)",y=0,alpha=-1,draws=1-apply(PoissonBetas,1,function(x){poisson_efun(beta=x,x=1)}))
)

nonlinP<-merge(nonlinP1 %>%
                 mutate(mscx=paste0(model,sc,x)),
               nonlinP2 %>%
                 mutate(mscx=paste0(model,sc,x))%>%
                 group_by(mscx)%>%
                 dplyr::summarise(ymax=mean_hdi(draws)$ymax,ymin=mean_hdi(draws)$ymin),
               by="mscx", all=TRUE)
nonlinP[which(nonlinP$sc=="If OR were equal to 1"),"ymin"]<-nonlinP[which(nonlinP$sc=="If OR were equal to 1"),"y"]
nonlinP[which(nonlinP$sc=="If OR were equal to 1"),"ymax"]<-nonlinP[which(nonlinP$sc=="If OR were equal to 1"),"y"]



catpp<-ggplot(data=nonlinP,aes(y=y,x=factor(model),color=sc))+
  geom_bar(aes(linetype=sc,alpha = alpha),fill="white",stat="identity",position = position_dodge(width = 0.5))+
  theme_bw()+ facet_grid(.~ factor(x))+ 
  scale_y_continuous(trans = squish_trans(0.01, 0.99, 15),
                     breaks = c(0,seq(0.05,0.95, by = 0.1),1)) + scale_alpha(guide = 'none')+
  labs(color="")+xlab("model")+ylab("")+ 
  geom_errorbar(aes(ymin=ymin,ymax=ymax),position = position_dodge(width = 0.5),width=.9,alpha=0.5)+
  scale_colour_manual(values=c("#316395","grey60","#B82E2E"),name="") + 
  scale_linetype_manual(values=c("solid","dotted","dashed"),name="")+
  theme(legend.position= "bottom",legend.justification = "left",legend.box.margin = margin(0,0,0,-0.75,"cm"))


################################################################################
# figure 4
ggarrange(linpp,catpp,nrow=1,widths = c(1.5,1))
################################################################################



