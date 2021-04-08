## load example data-set:
gw<-read.table("data/gibsonwu2012data.txt",
               header=TRUE)
## sum-contrast coding of predictor:
gw$so <- ifelse(
  gw$type%in%c("subj-ext"),-1,1)
## subset critical region
gw1<-subset(gw,region=="headnoun")

## load second data-set:
gw2<-read.table("data/gibsonwu2012datarepeat.txt",
                header=TRUE)
gw2$so <- ifelse(
  gw2$condition%in%c("subj-ext"),-1,1)

## frequentist analysis:
library(lme4)
m_lmer<-lmer(log(rt)~so + (1|subj),gw1)
summary(m_lmer)

library(brms)
priors <- c(set_prior("normal(6, 1.5)", class = "Intercept"),
            set_prior("normal(0, .01)", class = "b", 
                      coef = "so"),
            set_prior("normal(0, 1)", class = "sd"),
            set_prior("normal(0, 1)", class = "sigma"),
            set_prior("lkj(2)", class = "cor"))

m_gw<-brm(log(rt)~so + (1+so|subj),gw1,family=lognormal(),
          prior=priors)
summary(m_gw)

pp_check(m_gw)

## graphical visualization:
library(bayesplot)
postgw<-posterior_samples(m_gw)
## extract variances:
alpha<-postgw$b_Intercept
beta<-postgw$b_so
cor<-posterior_samples(m_gw,"^cor")
sd<-posterior_samples(m_gw,"^sd")
sigma<-posterior_samples(m_gw,"sigma")

## subject level effects:
subj_re<-posterior_samples(m_gw,"^r_subj")
meandiff<- exp(alpha + beta) - exp(alpha - beta)
mean(meandiff)
round(quantile(meandiff,prob=c(0.025,0.975)),0)

## mean effect:
hist(meandiff,freq=FALSE,
     main="Mean OR vs SR processing cost",
     xlab=expression(exp(alpha + beta)- exp(alpha - beta)))

## individual level estimates:
nsubj<-37
subjdiff<-matrix(rep(NA,nsubj*4000),nrow=nsubj)
for(i in 1:nsubj){
  subjdiff[i,]<-exp(alpha + subj_re[,i]  + (beta+subj_re[,i+nsubj])) - 
    exp(alpha + subj_re[,i] - 
          (beta+subj_re[,i+nsubj]))
}

subjdiff<-t(subjdiff)

subjdiff<-as.data.frame(subjdiff)
colnames(subjdiff)<-paste("s",c(1:nsubj),sep="")
mns <- colMeans(subjdiff)
subjdiff<-subjdiff[,order(mns)]
mcmc_areas(subjdiff)


## Bayes factors
## Bayes factor analysis:
m_gw<-brm(rt~so + (1+so|subj),gw1,family=lognormal(),
          prior=priors,warmup=5000,iter=20000,
          save_all_pars = TRUE)
summary(m_gw)

priors0 <- c(set_prior("normal(6, 1.5)", class = "Intercept"),
             set_prior("normal(0, 1)", class = "sd"),
             set_prior("normal(0, 1)", class = "sigma"),
             set_prior("lkj(2)", class = "cor"))

m_gw0<-brm(rt~1 + (1+so|subj),gw1,family=lognormal(),
           prior=priors0,,warmup=5000,iter=20000,
           save_all_pars = TRUE)

bayes_factor(m_gw,m_gw0)
