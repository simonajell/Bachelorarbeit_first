### All helper functions to produce the plots of section 4

emp_deriv<-function(x,y){ #this function determines an empirical derivative of a function f given a set of points (x,f(x))
  xres<-numeric(length(x)-1)
  yres<-numeric(length(y)-1)
  for(i in 1:length(x)-1){
    xres[i]<-x[i]+(x[i+1]-x[i])/2
    yres[i]<-(y[i+1]-y[i])/(x[i+1]-x[i])
  }
  return(data.frame(x=xres,y=yres))
}


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


logistic_unif<-function(beta){ #this function computes the proposed effect size measure for logistic regression, with X=[0,1] and mu chosen as associated with the U([0,1]) distribution
  return((exp(beta[1])*(exp(beta[2])-1))/((exp(beta[1])+1)*(exp(beta[1]+beta[2])+1)))
}

poisson_unif<-function(beta){ #this function computes the proposed effect size measure for poisson regression, with X=[0,1] and mu chosen as associated with the U([0,1]) distribution
  return(exp(beta[1])*(exp(beta[2])-1))
}

logistic_efun<-function(beta,x){ #this function computes the expectation for logistic regression
  return((exp(beta[1]+beta[2]*x))/(1+exp(beta[1]+beta[2]*x)))
}

poisson_efun<-function(beta,x){ #this function computes the expectation for poisson regression
  return(exp(beta[1]+beta[2]*x))
}

intval<-function(betas,regs,deriv=NULL){# this function computes the values of the expectation
                                        # and slope plots for the clinical trial data of Tortora (2020)
  betas <- as.numeric(betas)
  eta <-betas[1]+betas[2]*regs$assignmentG+betas[3]*regs$age+betas[4]*regs$genderM+betas[5]*regs$ethB+betas[6]*regs$ethL+betas[7]*regs$ethO
  if(is.null(deriv)){return(unlist(inv.logit(eta)))
  }else{
    return(unlist(inv.logit.deriv(eta,betas[deriv])))
  }
}

catAII1<-function(draws,data){# this function computes \Delta_s under assumption (A.I') for the clinical trial data of Tortora (2020)
  df<-dplyr::select(dummy_cols(data[,c("AGE","Assignment","GENDER","Ethnicity")],"Ethnicity"),-Ethnicity)
  names(df)<-c("age","assignmentG","genderM","ethW","ethB","ethL","ethO")
  df$genderM<-ifelse(df$genderM=="M",1,0)
  df$assignmentG<-ifelse(df$assignmentG=="G",1,0)
  out<-list()
  out$ethB<-apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=1,ethL=0,ethO=0),-"ethW"))-intval(x,dplyr::select(mutate(df,ethW=1,ethB=0,ethL=0,ethO=0),-"ethW")))/nrow(df)})
  out$ethL<-apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=1,ethO=0),-"ethW"))-intval(x,dplyr::select(mutate(df,ethW=1,ethB=0,ethL=0,ethO=0),-"ethW")))/nrow(df)})
  out$ethO<-apply(draws,1,function(x){sum(intval(x,dplyr::select(mutate(df,ethW=0,ethB=0,ethL=0,ethO=1),-"ethW"))-intval(x,dplyr::select(mutate(df,ethW=1,ethB=0,ethL=0,ethO=0),-"ethW")))/nrow(df)})
  return(out)
}

catAII2<-function(draws,data){# this function computes \Delta_s under assumption (A.I'') for the clinical trial data of Tortora (2020)
  df<-dplyr::select(dummy_cols(data[,c("AGE","Assignment","GENDER","Ethnicity")],"Ethnicity"),-Ethnicity)
  names(df)<-c("age","assignmentG","genderM","ethW","ethB","ethL","ethO")
  df$genderM<-ifelse(df$genderM=="M",1,0)
  df$assignmentG<-ifelse(df$assignmentG=="G",1,0)
  out<-list()
  out$ethB<-apply(draws,1,function(x){(sum(intval(x,dplyr::select(df[which(df$ethB==1),],-"ethW")))/length(which(df$ethB==1)))-(sum(intval(x,dplyr::select(df[which(df$ethW==1),],-"ethW")))/length(which(df$ethW==1)))})  
  out$ethL<-apply(draws,1,function(x){(sum(intval(x,dplyr::select(df[which(df$ethL==1),],-"ethW")))/length(which(df$ethL==1)))-(sum(intval(x,dplyr::select(df[which(df$ethW==1),],-"ethW")))/length(which(df$ethW==1)))}) 
  out$ethO<-apply(draws,1,function(x){(sum(intval(x,dplyr::select(df[which(df$ethO==1),],-"ethW"))/length(which(df$ethO==1))))-(sum(intval(x,dplyr::select(df[which(df$ethW==1),],-"ethW"))/length(which(df$ethW==1))))}) 
  return(out)
}



################################################################################
#### From https://stackoverflow.com/questions/35511951/r-ggplot2-collapse-or-remove-segment-of-y-axis-from-scatter-plot :

squish_trans <- function(from, to, factor) {
  trans <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < to
    ito <- x >= to
    # apply transformation
    x[isq] <- from + (x[isq] - from)/factor
    x[ito] <- from + (to - from)/factor + (x[ito] - to)
    return(x)}
  inv <- function(x) {
    if (any(is.na(x))) return(x)
    # get indices for the relevant regions
    isq <- x > from & x < from + (to - from)/factor
    ito <- x >= from + (to - from)/factor
    # apply transformation
    x[isq] <- from + (x[isq] - from) * factor
    x[ito] <- to + (x[ito] - (from + (to - from)/factor))
    return(x)}
  # return the transformation
  return(trans_new("squished", trans, inv))
}








