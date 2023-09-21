###############################################
## Croon's Estimation of mediation (1-1-1) effects in 2/1 partially nested designs

##
#Packages
library(MplusAutomation)
library(dplyr)
library(devtools)
library(lavaan)
library(lme4)
library(MASS)
library(psych)
library(mirt)
library(parallel)
library(foreach)
library(doParallel)
library(gdata)
library(Rfast)
library(chron)
library(powerjoin)
library(stringr)

#clean environment
rm(list=ls(all=TRUE))
# Set wd
setwd("")
#################################################################################################################################
###------------------------------------------------- Data Generation ###----------------------------------------------------
#################################################################################################################################
## Function Start
croon_partial_med <- function(datasets, n2, n1, weight1, weight2, weight3, ind_rhox,ind_rhom, ind_rho,  
                              xi00, beta0, beta0_c,
                              a, pi0, a_c,
                              covar_coef,
                              b1, B,  b1_c ,
                              rhox,rhom,rho){


# Codebook
  # datasets ## number of simulation replications or datasets created and analyzed
  
## Sample Size
  # n2 ## group level sample size
  # n1 ## individuals per group level sample size
  # nc ## total sample size in control group (n1*n2)
  
## Indicator Weights
  # weight1
  # weight2
  # weight3
## Indicator residuals (for multilevel factor models)
  # ind_rhox
  # ind_rhom
  # ind_rho
## Rho
  # rho   ## proportion of outcome indicator and factor variance at level-2
  # rhox  ## proportion of covariate indicator and factor variance at level-2
  # rhom  ## proportion of mediator indicator and factor variance at level-2
  
### Path Coefficients-all values for x and w covariates set to covar_coef
##Treat
  #OUTCOME
  # beta1    ## X  latent predictor Within
  # b1       ## M  latent predictor Within
  # B        ## M  latent predictor Between (treatment only)
  # xi1      ## X  latent predictor Between (treatment only)
  # xi2      ## W  latent predictor Between (treatment only)
#MEDIATOR
  # pi1       ## X latent predictor  Within
  # zeta1     ## X latent predictor  Between (treatment only)
  # zeta2     ## W  latent predictor Between (treatment only)
##Control
  #OUTCOME
  # beta1_c     ## X latent predictor
  # b1_c        ## M  latent predictor
#MEDIATOR
  # pi1_c       ## X latent predictor

## Intercept (group means for main effect)
  #  beta0        ##  intercept outcome Within 
  #  beta0_c      ##  intercept outcome Within 
  #  xi00         ##  intercept outcome Between 
  #  pi0          ##  intercept mediator Within
  #  a            ##  intercept mediator Between
  #  a_c          ##  intercept mediator

## Effects
# (a)*B         ## Total Mediation
# (a)*b1        ## Within Mediation
# (a)*(B-b1)    ## Between Mediation...set to 0? b1=0?
# xi00      ## Main effect
  
## Default/Test
     #rm(list=ls(all=TRUE)) 
      #datasets<-25; n2<-15; n1<-30; weight1 <-1; weight2 <- 1.5; weight3 <- 0.666; ind_rhox<-0.2; ind_rhom<-0.2; ind_rho<-0.2 
      #xi00<-0.8; beta0<-0.0; beta0_c<-0.0;
      #a<-0.5; pi0<-0.0; a_c<-0.0;
      #covar_coef<-0.3;
      #b1<-0.0; B<-0.5;  b1_c <-0.5;
      #rhox<-0.2; rhom<-0.2; rho<-0.2
      #loop<-1
  
### Data Generation Set-up ---------------------------------------
  allruns<-NULL
  runs <- NULL

  true_coefs_t <- NULL
  ml_coefs_t <- NULL
  bayes_coefs_t <- NULL
  bayesIN_coefs_t <- NULL
  fs_coefs_t <- NULL
  croon_coefs_t <- NULL
  true_coefsL1_t <- NULL
  ml_coefsL1_t <- NULL
  bayes_coefsL1_t <- NULL
  bayesIN_coefsL1_t <- NULL
  fs_coefsL1_t <- NULL
  croon_coefsL1_t <- NULL
  true_coefsL1_c <- NULL
  ml_coefsL1_c <- NULL
  bayes_coefsL1_c <- NULL
  bayesIN_coefsL1_c <- NULL
  fs_coefsL1_c <- NULL
  croon_coefsL1_c <- NULL
  
  true_conv_full <- NULL
  ml_conv_full <- NULL
  bayes_conv_full <- NULL
  bayes_psr_full <- NULL
  bayesIN_conv_full <- NULL
  bayesIN_psr_full <-NULL
  fs_conv_t <- NULL
  fs_convL1_t<- NULL
  croon_conv_t <- NULL
  croon_convL1_t <- NULL
  fs_convL1_c <- NULL
  croon_convL1_c <- NULL
  
  true_error_full <- NULL
  ml_error_full <- NULL
  bayes_error_full <- NULL
  bayesIN_error_full <- NULL
  fs_error_t <- NULL
  croon_error_t <- NULL
  fs_errorL1_t <- NULL
  croon_errorL1_t <- NULL
  fs_errorL1_c <- NULL
  croon_errorL1_c <- NULL
  
  ml_time_full<- NULL
  bayes_time_full<- NULL
  bayesIN_time_full<- NULL
  
#################################################################################################################################
###------------------------------------------------- Data Generation ---------------------------------------------------------###
#################################################################################################################################

##Loop Generate data Sets
  #numCores <- detectCores()
  #registerDoParallel(numCores)
  #registerDoParallel(cores = 8)
  #allruns<-foreach(loop=1:datasets,.packages = c('lavaan','MplusAutomation','MASS'),.combine = 'rbind') %dopar% { 

   for (loop in 1:datasets) {

runs <- NULL
d_t <-NULL
d_c <-NULL
d_full <-NULL

##Fixed Parameters and True Model
    # group ID
    id2<-rep(1:n2,each=n1) 
#treatment indicator
 tx<-1

# Control sample size 
    nc<- n1*n2
    id2_c <-rep((n2+1):(nc+n2))
#treatment indicator
 tx_c<-0 
    
### Covariates
    x_L2_t<-rnorm(n2,0,sqrt(rhox))
    x_L1_t<-rnorm(n1*n2,0,sqrt(1-rhox))
        
    w_L2_t<-rnorm(n2,0,sqrt(1))
    
    x_L1_c<-rnorm(n1*n2,0,sqrt(1))

### Mediator
  #Treat
    tau_m_L2_t<-rnorm(n2,0,sqrt(rhom))
    sig_m_L1_t<-rnorm(n1*n2,0,sqrt(1-rhom))
    m_L2_t<-a+(covar_coef*x_L2_t)+(covar_coef*w_L2_t)+tau_m_L2_t
    m_L1_t<-pi0+(covar_coef*x_L1_t)+sig_m_L1_t
  #Control
    sig_m_L1_c<-rnorm(n1*n2,0,sqrt(1))
    m_L1_c<-a_c+covar_coef*x_L1_c+sig_m_L1_c

### Outcome
  #Treat
    tau_y_L2_t<-rnorm(n2,0,sqrt(rho))
    sig_y_L1_t<-rnorm(n1*n2,0,sqrt(1-rho))
    y_L2_t<-xi00+(B*m_L2_t)+(covar_coef*x_L2_t)+(covar_coef*w_L2_t)+tau_y_L2_t
    y_L1_t<-beta0+(b1*m_L1_t)+(covar_coef*x_L1_t)+sig_y_L1_t
  #Control
    sig_y_L1_c<-rnorm(n1*n2,0,sqrt(1))
    y_L1_c<-beta0_c+(b1_c*m_L1_c)+(covar_coef*x_L1_c)+sig_y_L1_c 

##expand data from groups to individuals
    y_L2_t<-rep(y_L2_t,each=n1)
    m_L2_t<-rep(m_L2_t,each=n1)
    x_L2_t<-rep(x_L2_t,each=n1)
    w_L2_t<-rep(w_L2_t,each=n1)
    
##final variable values
    y_t<-y_L2_t+y_L1_t
    y_c<-y_L1_c
    m_t<-m_L2_t+m_L1_t
    m_c<-m_L1_c
    x_t<-x_L2_t+x_L1_t
    x_c<-x_L1_c
    w_t<-w_L2_t


###Indicator Generation
##Treat
    #creates indictors of true variables that are observable
    ##Individual level covariate indicators x
    x1_t<-weight1*x_L2_t+rep(rnorm(n2,0,sqrt(ind_rhox)),each=n1)+weight1*x_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rhox))
    x2_t<-weight2*x_L2_t+rep(rnorm(n2,0,sqrt(ind_rhox)),each=n1)+weight2*x_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rhox))
    x3_t<-weight3*x_L2_t+rep(rnorm(n2,0,sqrt(ind_rhox)),each=n1)+weight3*x_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rhox))

    ##Individual level mediator indicators m
    m1_t<-weight1*m_L2_t+rep(rnorm(n2,0,sqrt(ind_rhom)),each=n1)+weight1*m_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rhom))
    m2_t<-weight2*m_L2_t+rep(rnorm(n2,0,sqrt(ind_rhom)),each=n1)+weight2*m_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rhom))
    m3_t<-weight3*m_L2_t+rep(rnorm(n2,0,sqrt(ind_rhom)),each=n1)+weight3*m_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rhom))

    ##Cluster level covariate indicators w
    w1_t<-weight1*w_L2_t+rep(rnorm(n2,0,sqrt(1)),each=n1)             
    w2_t<-weight2*w_L2_t+rep(rnorm(n2,0,sqrt(1)),each=n1)
    w3_t<-weight3*w_L2_t+rep(rnorm(n2,0,sqrt(1)),each=n1)

    ##Individual level outcome indicators y
    y1_t<-weight1*y_L2_t+rep(rnorm(n2,0,sqrt(ind_rho)),each=n1)+weight1*y_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rho))
    y2_t<-weight2*y_L2_t+rep(rnorm(n2,0,sqrt(ind_rho)),each=n1)+weight2*y_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rho))
    y3_t<-weight3*y_L2_t+rep(rnorm(n2,0,sqrt(ind_rho)),each=n1)+weight3*y_L1_t+rnorm(n1*n2,0,sqrt(1-ind_rho))

    
##Control
    ##Individual level covariate indicators x
    x1_c<-weight1*x_L1_c+rnorm(n1*n2,0,sqrt(1))
    x2_c<-weight2*x_L1_c+rnorm(n1*n2,0,sqrt(1))
    x3_c<-weight3*x_L1_c+rnorm(n1*n2,0,sqrt(1))

    ##Individual level mediator indicators m
    m1_c<-weight1*m_L1_c+rnorm(n1*n2,0,sqrt(1))
    m2_c<-weight2*m_L1_c+rnorm(n1*n2,0,sqrt(1))
    m3_c<-weight3*m_L1_c+rnorm(n1*n2,0,sqrt(1))

    ##Individual level outcome indicators y
    y1_c<-weight1*y_L1_c+rnorm(n1*n2,0,sqrt(1))
    y2_c<-weight2*y_L1_c+rnorm(n1*n2,0,sqrt(1))
    y3_c<-weight3*y_L1_c+rnorm(n1*n2,0,sqrt(1))

### Save data frame
  ##create data frame with variables----Treatment
    d_t<-data.frame(id2, n2, n1,nc,tx, rho, rhox, rhom, weight1, weight2, weight3,ind_rhox,ind_rho,ind_rhom,  
                  beta0,xi00,pi0,a,
                  covar_coef,b1, B,
                  x1_t,x2_t,x3_t,m1_t,m2_t,m3_t,w1_t,w2_t,w3_t,y1_t,y2_t,y3_t,
                  w_t, x_t,m_t,y_t,  
                  w_L2_t,x_L2_t,x_L1_t,m_L2_t,m_L1_t,y_L2_t,y_L1_t)
    ##aggregated version of data
    dagg<-aggregate(d_t,by=list(d_t$id2),FUN=mean)
    
 ##create data frame with variables----Control
    # match ids for merge
    tx<-tx_c
    id2<-id2_c
    d_c<-data.frame(id2, n2, n1,nc,tx, rho, rhox, rhom, weight1, weight2, weight3,ind_rhox,ind_rho,ind_rhom,  
                  beta0_c, a_c,
                  covar_coef, b1_c,
                  x_c,m_c,y_c,
                  x1_c,x2_c,x3_c,m1_c,m2_c,m3_c,y1_c,y2_c,y3_c,
                  x_L1_c,m_L1_c,y_L1_c)
    

d_full<-power_full_join(d_t, d_c, by = "id2", conflict = coalesce_xy)
d_full[is.na(d_full)] <- -99

### Save Data ***Separate***------------------------------------------------------------------------------
# Save data frame as a file----------------------------for data creation and analysis together
filename_t = paste0("d_t", ".dat", sep="")
filename_dagg = paste0("dagg","1", ".dat", sep="")
filename_c = paste0("d_c", ".dat", sep="")
filename_full = paste0("d_full", ".dat", sep="")
# format for Mplus: row.names=F,col.names=F 
#write.table(d,file =filename_full,row.names=F,col.names=F)
write.table(d_t,file =filename_t,row.names=F,col.names=F)
write.table(dagg,file =filename_dagg,row.names=F,col.names=F)
write.table(d_c,file =filename_c,row.names=F,col.names=F)
write.table(d_full,file =filename_full,row.names=F,col.names=F)


#################################################################################################################################
###------------------------------------------------- Data Analysis ###----------------------------------------------------
################################################################################################################################# 
    
###############################################################################
### Estimator Comparisons ###--------------------------------------------------
###############################################################################

#############################################################
### True model FULL:Mplus
#############################################################
#treat
  true_model_full<-'
TITLE: True Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
    NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
    USEVARIABLES ARE id2 w_t x_t m_t y_t x_c  m_c  y_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;                  !identify clustering variable
    BETWEEN IS w_t;    
    WITHIN IS ;
ANALYSIS: TYPE IS TWOLEVEL;
    !use defaults
    !type = general;
    !estimator = ML;
    PROCESSORS = 16;
MODEL:               !model for tx (nested) group
%WITHIN%             !within-cluster (level-1) model for tx group
y_t ON m_t (b1); 
y_t ON x_t (beta1);
m_t ON x_t (pi1); 

y_c@0; m_c@0; x_c@0;   
y_c ON m_c@0; 
y_c ON x_c@0;
m_c ON x_c@0;
y_t with y_c@0; x_t with y_c@0; m_t with y_c@0; 
y_t with m_c@0; x_t with m_c@0; m_t with m_c@0; 
y_t with x_c@0; x_t with x_c@0; m_t with x_c@0;

%BETWEEN%                !between-cluster (level-2) model for tx group
y_t ON m_t (B); 
y_t ON x_t (xi1);
y_t ON w_t (xi2);
m_t ON x_t (zeta1);
m_t ON w_t (zeta2);

! y_c ON m_c (b1_c);  ! Mediation VS Moderated mediation 
y_c ON m_c (B);
y_c ON x_c (beta1_c);
m_c ON x_c (pi1_c);

y_t with y_c@0; x_t with y_c@0; m_t with y_c@0; w_t with y_c@0;
y_t with m_c@0; x_t with m_c@0; m_t with m_c@0; w_t with m_c@0;
y_t with x_c@0; x_t with x_c@0; m_t with x_c@0; w_t with x_c@0;

[y_c@0 m_c@0 y_t m_t];                 !estimate means

OUTPUT:
  NOCHISQUARE TECH1 TECH3;
'
### Save for Mplus- .inp
#treat
write.table(true_model_full,file="true_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="true_model_full.mplus.inp")
### Save results from Mplus
mpfit_true_full<-readModels(target="true_model_full.mplus.out")
### Check results
#mpfit_true_full$parameters$unstandardized



#############################################################
### Mplus ML                            
#############################################################
#treat
  ml_model_full<-'
TITLE: ML Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE id2 x1_t x2_t x3_t m1_t m2_t m3_t w1_t w2_t w3_t y1_t y2_t y3_t
x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c y3_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS w1_t w2_t w3_t;   !identify between-cluster (level-2) variables
ANALYSIS:
    !use defaults?
    !type = general;
    PROCESSORS = 16;
    ESTIMATOR = ML;
    TYPE = TWOLEVEL;
	  !ALGORITHM = INTEGRATION;
    !INTEGRATION = MONTECARLO(5000);
MODEL:
%WITHIN%
!Measurement model
xt1  BY x1_t  
x2_t  x3_t (Bx2 Bx3); 
mt1  BY m1_t  
m2_t  m3_t (Bm2 Bm3);
yt1  BY y1_t  
y2_t  y3_t (By2 By3);

xc1  BY x1_c@0  x2_c@0  x3_c@0; !(Wx_c1-Wx_c3); 
mc1  BY m1_c@0  m2_c@0  m3_c@0; !(Wm_c1-Wm_c3);
yc1  BY y1_c@0  y2_c@0  y3_c@0; !(Wy_c1-Wy_c3);


!Path model
yt1 ON mt1 (b1); 
yt1 ON xt1 (beta1);
mt1 ON xt1 (pi1); 

x1_c@0;  x2_c@0;  x3_c@0;       ! set indicator variance to 0
m1_c@0;  m2_c@0;  m3_c@0; 
y1_c@0;  y2_c@0;  y3_c@0;
yc1@0; mc1@0; xc1@0;            
yc1 ON mc1@0; 
yc1 ON xc1@0;
mc1 ON xc1@0;
yt1 with yc1@0; xt1 with yc1@0; mt1 with yc1@0; 
yt1 with mc1@0; xt1 with mc1@0; mt1 with mc1@0; 
yt1 with xc1@0; xt1 with xc1@0; mt1 with xc1@0;

%BETWEEN%
!Measurement model- hold loadings equal across T and C ((B.1-B.3))
xt2  BY x1_t  
x2_t  x3_t (Bx2 Bx3);
wt2  BY w1_t  
w2_t  w3_t (Bw2 Bw3);
mt2  BY m1_t  
m2_t  m3_t (Bm2 Bm3);
yt2  BY y1_t  
y2_t  y3_t (By2 By3);

xc2  BY x1_c  
x2_c  x3_c (Bx2 Bx3); 
mc2  BY m1_c  
m2_c  m3_c (Bm2 Bm3);
yc2  BY y1_c  
y2_c  y3_c (By2 By3);

! hold intercepts equal across T and C
[x1_t x1_c] (1);[m1_t m1_c] (2);[y1_t y1_c] (3);  
[x2_t x2_c] (4);[m2_t m2_c] (5);[y2_t y2_c] (6);
[x3_t x3_c] (7);[m3_t m3_c] (8);[y3_t y3_c] (9);

!Path model
yt2 ON mt2 (B); 
yt2 ON xt2 (xi1);
yt2 ON wt2 (xi2);
mt2 ON xt2 (zeta1);
mt2 ON wt2 (zeta2);

! yc2 ON mc2 (b1_c);
yc2 ON mc2 (B);
yc2 ON xc2 (beta1_c);
mc2 ON xc2 (pi1_c);

yt2 with yc2@0; xt2 with yc2@0; mt2 with yc2@0; wt2 with yc2@0;
yt2 with mc2@0; xt2 with mc2@0; mt2 with mc2@0; wt2 with mc2@0;
yt2 with xc2@0; xt2 with xc2@0; mt2 with xc2@0; wt2 with xc2@0;


[yc2@0 mc2@0 yt2 mt2];                 !estimate means	

OUTPUT:
  NOCHISQUARE tech1 tech8;
'
### Save for Mplus- .inp
write.table(ml_model_full,file="ml_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="ml_model_full.mplus.inp")
### Save results from Mplus
mpfit_ml_full<-readModels(target="ml_model_full.mplus.out")
### Check results
#mpfit_ml_full$parameters$unstandardized


#############################################################
### Mplus BAYES (uninformative prior)                            
#############################################################
#treat
  bayes_model_full<-'
TITLE: bayes Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE id2 x1_t x2_t x3_t m1_t m2_t m3_t w1_t w2_t w3_t y1_t y2_t y3_t
x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c y3_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS w1_t w2_t w3_t;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = BAYES;
    TYPE = TWOLEVEL;
    VARIANCE=0.001;
    !ALGORITHM=GIBBS(PX1);
    !PROCESS = 2;    ! define the number of chains;
    !POINT=MEAN;     ! define the estimator;
    BITERATIONS = 20000 (0); ! max number of iterations using psr;
	  PROCESSORS = 16;
MODEL:
%WITHIN%
!Measurement model
xt1  BY x1_t  
x2_t  x3_t (Bx2 Bx3); 
mt1  BY m1_t  
m2_t  m3_t (Bm2 Bm3);
yt1  BY y1_t  
y2_t  y3_t (By2 By3);

xc1  BY x1_c@0  x2_c@0  x3_c@0; !(Wx_c1-Wx_c3); 
mc1  BY m1_c@0  m2_c@0  m3_c@0; !(Wm_c1-Wm_c3);
yc1  BY y1_c@0  y2_c@0  y3_c@0; !(Wy_c1-Wy_c3);


!Path model
yt1 ON mt1 (b1); 
yt1 ON xt1 (beta1);
mt1 ON xt1 (pi1); 

x1_c@0;  x2_c@0;  x3_c@0;       ! set indicator variance to 0
m1_c@0;  m2_c@0;  m3_c@0; 
y1_c@0;  y2_c@0;  y3_c@0;
yc1@0; mc1@0; xc1@0;            
yc1 ON mc1@0; 
yc1 ON xc1@0;
mc1 ON xc1@0;
yt1 with yc1@0; xt1 with yc1@0; mt1 with yc1@0; 
yt1 with mc1@0; xt1 with mc1@0; mt1 with mc1@0; 
yt1 with xc1@0; xt1 with xc1@0; mt1 with xc1@0;

%BETWEEN%
!Measurement model- hold loadings equal across T and C
xt2  BY x1_t  
x2_t  x3_t (Bx2 Bx3);
wt2  BY w1_t  
w2_t  w3_t (Bw2 Bw3);
mt2  BY m1_t  
m2_t  m3_t (Bm2 Bm3);
yt2  BY y1_t  
y2_t  y3_t (By2 By3);

xc2  BY x1_c  
x2_c  x3_c (Bx2 Bx3); 
mc2  BY m1_c  
m2_c  m3_c (Bm2 Bm3);
yc2  BY y1_c  
y2_c  y3_c (By2 By3);

! hold intercepts equal across T and C
[x1_t x1_c] (1);[m1_t m1_c] (2);[y1_t y1_c] (3);  
[x2_t x2_c] (4);[m2_t m2_c] (5);[y2_t y2_c] (6);
[x3_t x3_c] (7);[m3_t m3_c] (8);[y3_t y3_c] (9);

!Path model
yt2 ON mt2 (B); 
yt2 ON xt2 (xi1);
yt2 ON wt2 (xi2);
mt2 ON xt2 (zeta1);
mt2 ON wt2 (zeta2);

! yc2 ON mc2 (b1_c);
yc2 ON mc2 (B);
yc2 ON xc2 (beta1_c);
mc2 ON xc2 (pi1_c);

yt2 with yc2@0; xt2 with yc2@0; mt2 with yc2@0; wt2 with yc2@0;
yt2 with mc2@0; xt2 with mc2@0; mt2 with mc2@0; wt2 with mc2@0;
yt2 with xc2@0; xt2 with xc2@0; mt2 with xc2@0; wt2 with xc2@0;


[yc2@0 mc2@0 yt2 mt2];                 !estimate means	

OUTPUT:
  NOCHISQUARE tech1 tech8;
'

### Save for Mplus- .inp
write.table(bayes_model_full,file="bayes_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="bayes_model_full.mplus.inp")
### Save results from Mplus
mpfit_bayes_full<-readModels(target="bayes_model_full.mplus.out")
### Check results
#mpfit_bayes_full$parameters$unstandardized
#mpfit_bayes_full$tech8$psr



#############################################################
### Mplus BAYES (informative prior)                            
#############################################################
#treat
  bayesIN_model_full<-'
TITLE: bayesIN Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE id2 x1_t x2_t x3_t m1_t m2_t m3_t w1_t w2_t w3_t y1_t y2_t y3_t
x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c y3_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS w1_t w2_t w3_t;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = BAYES;
    TYPE = TWOLEVEL;
    VARIANCE=0.001;
    !ALGORITHM=GIBBS(PX1);
    !PROCESS = 2;    ! define the number of chains;
    !POINT=MEAN;     ! define the estimator;
    BITERATIONS = 20000 (0); ! max number of iterations using psr;
	  PROCESSORS = 16;
MODEL:
%WITHIN%
!Measurement model
xt1  BY x1_t  
x2_t  x3_t (Bx2 Bx3); 
mt1  BY m1_t  
m2_t  m3_t (Bm2 Bm3);
yt1  BY y1_t  
y2_t  y3_t (By2 By3);

xc1  BY x1_c@0  x2_c@0  x3_c@0; 
mc1  BY m1_c@0  m2_c@0  m3_c@0; 
yc1  BY y1_c@0  y2_c@0  y3_c@0;


!Path model
yt1 ON mt1 (b1); 
yt1 ON xt1 (beta1);
mt1 ON xt1 (pi1); 

x1_c@0;  x2_c@0;  x3_c@0;       ! set indicator variance to 0
m1_c@0;  m2_c@0;  m3_c@0; 
y1_c@0;  y2_c@0;  y3_c@0;
yc1@0; mc1@0; xc1@0;            ! set variance 0
yc1 ON mc1@0; 
yc1 ON xc1@0;
mc1 ON xc1@0;
yt1 with yc1@0; xt1 with yc1@0; mt1 with yc1@0; 
yt1 with mc1@0; xt1 with mc1@0; mt1 with mc1@0; 
yt1 with xc1@0; xt1 with xc1@0; mt1 with xc1@0;

!Name for priors
yt1  (yt1lat);
mt1  (mt1lat);
xt1 (xt1lat);

%BETWEEN%
!Measurement model- hold loadings equal across T and C
xt2  BY x1_t  
x2_t  x3_t (Bx2 Bx3);
wt2  BY w1_t  
w2_t  w3_t (Bw2 Bw3);
mt2  BY m1_t  
m2_t  m3_t (Bm2 Bm3);
yt2  BY y1_t  
y2_t  y3_t (By2 By3);

xc2  BY x1_c  
x2_c  x3_c (Bx2 Bx3); 
mc2  BY m1_c  
m2_c  m3_c (Bm2 Bm3);
yc2  BY y1_c  
y2_c  y3_c (By2 By3);

! hold intercepts equal across T and C
[x1_t x1_c] (1);[m1_t m1_c] (2);[y1_t y1_c] (3);  
[x2_t x2_c] (4);[m2_t m2_c] (5);[y2_t y2_c] (6);
[x3_t x3_c] (7);[m3_t m3_c] (8);[y3_t y3_c] (9);

!Path model
yt2 ON mt2 (B); 
yt2 ON xt2 (xi1);
yt2 ON wt2 (xi2);
mt2 ON xt2 (zeta1);
mt2 ON wt2 (zeta2);

! yc2 ON mc2 (b1_c);
yc2 ON mc2 (B);
yc2 ON xc2 (beta1_c);
mc2 ON xc2 (pi1_c);

yt2 with yc2@0; xt2 with yc2@0; mt2 with yc2@0; wt2 with yc2@0;
yt2 with mc2@0; xt2 with mc2@0; mt2 with mc2@0; wt2 with mc2@0;
yt2 with xc2@0; xt2 with xc2@0; mt2 with xc2@0; wt2 with xc2@0;


[yc2@0 mc2@0 yt2 mt2]; !(xi0_c a_c  xi0 a_t);                 !estimate means

!Name for priors
yt2  (xi0); mt2  (a_t);     !xt2  (xt2lat); wt2  (xt2lat);
yc2  (xi0_c); mc2  (a_c);   !xc2  (xc2lat);

MODEL PRIORS: 
Bx2 ~ N(1.5,1);
Bx3 ~ N(.666,1);
Bm2 ~ N(1.5,1);
Bm3 ~ N(.666,1);
By2 ~ N(1.5,1);
By3 ~ N(.666,1);
Bw2 ~ N(1.5,1);
Bw3 ~ N(.666,1);

xi0_c~N(0,1); 
a_c~N(0,1); 
xi0~N(0.8,1);
a_t~N(0.5,1) ;

b1~ N(0,1); 
beta1~ N(0.3,1);
pi1~ N(0.3,1); 
B~ N(0.5,1); 
xi1~ N(0.3,1);
xi2~ N(0.3,1);
zeta1~ N(0.3,1);
zeta2~ N(0.3,1);

beta1_c~ N(0.3,1);
pi1_c~ N(0.3,1);

OUTPUT:
  NOCHISQUARE tech1 tech8;
'

### Save for Mplus- .inp
write.table(bayesIN_model_full,file="bayesIN_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="bayesIN_model_full.mplus.inp")
### Save results from Mplus
mpfit_bayesIN_full<-readModels(target="bayesIN_model_full.mplus.out")
### Check results
#mpfit_bayesIN_full$parameters$unstandardized
#mpfit_bayesIN_full$tech8$psr


#############################################################
### SEM-Sequential Approach: FS and Croon                             
#############################################################
#######################################
### Define individual measurement models ###--------------------------------
#######################################
#######################################
##factor model for X
#######################################
croon_x_model_full<-'
TITLE: Croon X Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE id2 x1_t x2_t x3_t x1_c x2_c x3_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
    TYPE = TWOLEVEL;
	  PROCESSORS = 16;

MODEL:
%WITHIN%
!Measurement model
xt1  BY x1_t  
x2_t  x3_t (Bx2 Bx3);

xc1  BY x1_c@0  x2_c@0  x3_c@0; 
x1_c@0;  x2_c@0;  x3_c@0;       ! set indicator variance to 0
xc1@0;                          ! set variance 0
xt1 with xc1@0;

%BETWEEN%
!Measurement model
xt2  BY x1_t  
x2_t  x3_t (Bx2 Bx3);

xc2  BY x1_c  
x2_c  x3_c (Bx2 Bx3);
! hold intercepts equal across T and C
[x1_t x1_c] (1);
[x2_t x2_c] (4);
[x3_t x3_c] (7);
xt2 with xc2@0;


OUTPUT:
  tech1 tech8 NOCHISQUARE; 
 Savedata: 
 file is f_x_full.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_x_model_full,file="croon_x_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_x_model_full.mplus.inp")
### Save model
x_model_full <-readModels(target="croon_x_model_full.mplus.out")
#######################################
##factor model for Xind
#######################################
croon_x_model_ind<-'
TITLE: Croon X Model_ind
DATA: FILE IS d_t.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0  xi00  pi0 a covar_coef b1 B 
        x1_t  x2_t  x3_t  m1_t  m2_t  m3_t  w1_t  w2_t  w3_t      
        y1_t  y2_t  y3_t  w_t x_t m_t y_t w_L2_t  x_L2_t  x_L1_t
        m_L2_t  m_L1_t  y_L2_t y_L1_t;
USEVARIABLES ARE id2 x1_t x2_t x3_t;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
    TYPE = TWOLEVEL;
	  PROCESSORS = 16;
MODEL:
%WITHIN%
!Measurement model
xt by x1_t x2_t  x3_t;  

%BETWEEN%
!Measurement model
xt2 by x1_t x2_t  x3_t;
OUTPUT:
  tech1 tech8 sampstat;
   Savedata:
    file is f_x_t.dat;
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_x_model_ind,file="croon_x_model_ind.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_x_model_ind.mplus.inp")
### Save model
x_model_ind <-readModels(target="croon_x_model_ind.mplus.out")
#######################################
##factor model for xind_c
#######################################
croon_x_model_ind_c<-'
TITLE: Croon x Model_ind_c
DATA: FILE IS d_c.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0_c a_c  covar_coef b1_c x_c m_c y_c
        x1_c  x2_c  x3_c  m1_c  m2_c  m3_c  y1_c  y2_c  y3_c      
        x_L1_c m_L1_c y_L1_c;
USEVARIABLES ARE x1_c  x2_c  x3_c;
    MISSING ARE ALL (-99);
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
	  PROCESSORS = 16;
MODEL:
!Measurement model
xc by x1_c x2_c  x3_c;  
OUTPUT:
  tech1 tech8 sampstat;
     Savedata:
    file is f_m_t.dat;
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_x_model_ind_c,file="croon_x_model_ind_c.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_x_model_ind_c.mplus.inp")
### Save model
x_model_ind_c <-readModels(target="croon_x_model_ind_c.mplus.out")

### Save results from Mplus------------------------------------------------------------------------
d_t<-NULL
d_c<-NULL
dagg<-NULL
### Save results from Mplus------------------------------------------------------------------------
#L1
#factor scores
d_c$fxL1_c<-try(x_model_ind_c$savedata[,"XC"],silent = T)
d_t$fxL1<-try(x_model_ind$savedata[,"XT"],silent = T)
#lambda- factor loadings
LxL1_c<-x_model_ind_c$parameters$unstandardized[1:3,"est"]   
LxL1<-x_model_ind$parameters$unstandardized[1:3,"est"]   
#variance
vxL1_c<-x_model_ind_c$parameters$unstandardized[7,"est"]     
vxL1<-x_model_ind$parameters$unstandardized[4,"est"]     

#---------------------------------------------------------------------------------------------------------------------------
# xindcovL1_t<-try(cov(x_model_full$savedata[(nc+1):(nc*2),1:3]),silent = T)
xindcovL1_t<-try(x_model_ind$sampstat$WITHIN$covariances,silent = T)
ifelse(class(xindcovL1_t)=="try-error"| is.null(xindcovL1_t)==T,
       xindcovL1_t<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),xindcovL1_t<-xindcovL1_t)
xindcovL1_t[upper.tri(xindcovL1_t)]<-t(xindcovL1_t)[upper.tri(xindcovL1_t)]
AxL1<-try(solve(xindcovL1_t) %*% LxL1  *vxL1,silent = T)                       

#xindcovL1_c<-try(cov(x_model_full$savedata[1:nc,4:6]),silent = T)
xindcovL1_c<-try(x_model_ind_c$sampstat$covariances,silent = T)
ifelse(class(xindcovL1_c)=="try-error"| is.null(xindcovL1_c)==T,
       xindcovL1_c<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),xindcovL1_c<-xindcovL1_c)
xindcovL1_c[upper.tri(xindcovL1_c)]<-t(xindcovL1_c)[upper.tri(xindcovL1_c)]
AxL1_c<-try(solve(xindcovL1_c) %*% LxL1_c  *vxL1_c,silent = T)
#---------------------------------------------------------------------------------------------------------------------------
#L2 and reliability 
#factor scores
fxL2<-try(x_model_ind$savedata[,"XT2"],silent = T)
dagg$fxL2<-try(fxL2[seq(1, length(fxL2), n1)],silent = T)                          
Lx<-x_model_ind$parameters$unstandardized[8:10,"est"]                 
vxL2<-x_model_ind$parameters$unstandardized[14,"est"]

xindcovL2<-try(x_model_ind$sampstat$BETWEEN$covariances,silent = T)
ifelse(class(xindcovL2)=="try-error"| is.null(xindcovL2)==T,
       xindcovL2<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),xindcovL2<-xindcovL2)
xindcovL2[upper.tri(xindcovL2)]<-t(xindcovL2)[upper.tri(xindcovL2)]
Ax<-try(solve(xindcovL2) %*% Lx *vxL2,silent = T)                                     
Rx<-try(xindcovL2%*%solve(xindcovL2+xindcovL1_t/n1),silent = T)                     
# check
# x_model_t$parameters$unstandardized


#######################################
##factor model for M
#######################################
croon_m_model_full<-'
TITLE: Croon M Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE id2 m1_t m2_t m3_t m1_c m2_c m3_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
    TYPE = TWOLEVEL;
	  PROCESSORS = 16;

MODEL:
%WITHIN%
!Measurement model
mt1  BY m1_t  
m2_t  m3_t (Bm2 Bm3);

mc1  BY m1_c@0  m2_c@0  m3_c@0; !(Wm_c1-Wm_c3);
m1_c@0;  m2_c@0;  m3_c@0;       ! set indicator variance to 0
mc1@0;                          ! set variance 0
mt1 with mc1@0;

%BETWEEN%
!Measurement model
mt2  BY m1_t  
m2_t  m3_t (Bm2 Bm3);

mc2  BY m1_c  
m2_c  m3_c (Bm2 Bm3);
! hold intercepts equal across T and C
[m1_t m1_c] (1);
[m2_t m2_c] (4);
[m3_t m3_c] (7);
mt2 with mc2@0;

[mc2@0 mt2];                 !estimate means


OUTPUT:
  tech1 tech8 NOCHISQUARE; 
 Savedata: 
 file is f_m_full.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_m_model_full,file="croon_m_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_m_model_full.mplus.inp")
### Save model
m_model_full <-readModels(target="croon_m_model_full.mplus.out")
# check
# m_model_full$parameters$unstandardized
#######################################
##factor model for mind
#######################################
croon_m_model_ind<-'
TITLE: Croon m Model_ind
DATA: FILE IS d_t.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0  xi00  pi0 a covar_coef b1 B 
        x1_t  x2_t  x3_t  m1_t  m2_t  m3_t  w1_t  w2_t  w3_t      
        y1_t  y2_t  y3_t  w_t x_t m_t y_t w_L2_t  x_L2_t  x_L1_t
        m_L2_t  m_L1_t  y_L2_t y_L1_t;
USEVARIABLES ARE id2 m1_t m2_t m3_t;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
    TYPE = TWOLEVEL;
	  PROCESSORS = 16;
MODEL:
%WITHIN%
!Measurement model
mt by m1_t m2_t  m3_t;  

%BETWEEN%
!Measurement model
mt2 by m1_t m2_t  m3_t;
OUTPUT:
  tech1 tech8 sampstat;
   Savedata:
    file is f_m_t.dat;
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_m_model_ind,file="croon_m_model_ind.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_m_model_ind.mplus.inp")
### Save model
m_model_ind <-readModels(target="croon_m_model_ind.mplus.out")
#######################################
##factor model for mind_c
#######################################
croon_m_model_ind_c<-'
TITLE: Croon m Model_ind_c
DATA: FILE IS d_c.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0_c a_c  covar_coef b1_c x_c m_c y_c
        x1_c  x2_c  x3_c  m1_c  m2_c  m3_c  y1_c  y2_c  y3_c      
        x_L1_c m_L1_c y_L1_c;
USEVARIABLES ARE m1_c  m2_c  m3_c;
    MISSING ARE ALL (-99);
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
	  PROCESSORS = 16;
MODEL:
!Measurement model
mc by m1_c m2_c  m3_c;  
OUTPUT:
  tech1 tech8 sampstat;
     Savedata:
    file is f_m_t.dat;
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_m_model_ind_c,file="croon_m_model_ind_c.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_m_model_ind_c.mplus.inp")
### Save model
m_model_ind_c <-readModels(target="croon_m_model_ind_c.mplus.out")

### Save results from Mplus------------------------------------------------------------------------
#L1
#factor scores
d_c$fmL1_c<-try(m_model_ind_c$savedata[,"MC"],silent = T)
d_t$fmL1<-try(m_model_ind$savedata[,"MT"],silent = T)
#lambda- factor loadings
LmL1_c<-m_model_ind_c$parameters$unstandardized[1:3,"est"]   
LmL1<-m_model_ind$parameters$unstandardized[1:3,"est"]   
#variance
vmL1_c<-m_model_ind_c$parameters$unstandardized[7,"est"]     
vmL1<-m_model_ind$parameters$unstandardized[4,"est"]     
#---------------------------------------------------------------------------------------------------

# mindcovL1_t<-try(cov(m_model_full$savedata[(nc+1):(nc*2),1:3]),silent = T)
mindcovL1_t<-try(m_model_ind$sampstat$WITHIN$covariances,silent = T)
ifelse(class(mindcovL1_t)=="try-error"| is.null(mindcovL1_t)==T,
       mindcovL1_t<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),mindcovL1_t<-mindcovL1_t)
mindcovL1_t[upper.tri(mindcovL1_t)]<-t(mindcovL1_t)[upper.tri(mindcovL1_t)]
AmL1<-try(solve(mindcovL1_t) %*% LmL1  *vmL1,silent = T)                       

#mindcovL1_c<-try(cov(m_model_full$savedata[1:nc,4:6]),silent = T)
mindcovL1_c<-try(m_model_ind_c$sampstat$covariances,silent = T)
ifelse(class(mindcovL1_c)=="try-error"| is.null(mindcovL1_c)==T,
       mindcovL1_c<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),mindcovL1_c<-mindcovL1_c)
mindcovL1_c[upper.tri(mindcovL1_c)]<-t(mindcovL1_c)[upper.tri(mindcovL1_c)]
AmL1_c<-try(solve(mindcovL1_c) %*% LmL1_c  *vmL1_c,silent = T)
#---------------------------------------------------------------------------------------------------------------------------
#L2 and reliability 
#factor scores
fmL2<-try(m_model_ind$savedata[,"MT2"],silent = T)
dagg$fmL2<-try(fmL2[seq(1, length(fmL2), n1)],silent = T)                          
Lm<-m_model_ind$parameters$unstandardized[8:10,"est"]                 
vmL2<-m_model_ind$parameters$unstandardized[14,"est"]

mindcovL2<-try(m_model_ind$sampstat$BETWEEN$covariances,silent = T)
ifelse(class(mindcovL2)=="try-error"| is.null(mindcovL2)==T,
       mindcovL2<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),mindcovL2<-mindcovL2)
mindcovL2[upper.tri(mindcovL2)]<-t(mindcovL2)[upper.tri(mindcovL2)]
Am<-try(solve(mindcovL2) %*% Lm *vmL2,silent = T)                                     
Rm<-try(mindcovL2%*%solve(mindcovL2+mindcovL1_t/n1),silent = T)                     
# check
# m_model_full$parameters$unstandardized


#######################################
##factor model for Y
#######################################
croon_y_model_full<-'
TITLE: Croon Y Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE id2 y1_t y2_t y3_t y1_c y2_c y3_c;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
    TYPE = TWOLEVEL;
	  PROCESSORS = 16;

MODEL:
%WITHIN%
!Measurement model
yt1  BY y1_t  
y2_t  y3_t (By2 By3);

yc1  BY y1_c@0  y2_c@0  y3_c@0; 
y1_c@0;  y2_c@0;  y3_c@0;       ! set indicator variance to 0
yc1@0;                          ! set variance 0
yt1 with yc1@0;

%BETWEEN%
!Measurement model
yt2  BY y1_t  
y2_t  y3_t (By2 By3);

yc2  BY y1_c  
y2_c  y3_c (By2 By3);
! hold intercepts equal across T and C
[y1_t y1_c] (1);
[y2_t y2_c] (4);
[y3_t y3_c] (7);
yt2 with yc2@0;

[yc2@0 yt2];                 !estimate means

OUTPUT:
  tech1 tech8 NOCHISQUARE; 
 Savedata: 
 file is f_y_full.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_y_model_full,file="croon_y_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_y_model_full.mplus.inp")
### Save model
y_model_full <-readModels(target="croon_y_model_full.mplus.out")
# check
# y_model_full$parameters$unstandardized
#######################################
##factor model for yind
#######################################
croon_y_model_ind<-'
TITLE: Croon y Model_ind
DATA: FILE IS d_t.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0  xi00  pi0 a covar_coef b1 B 
        x1_t  x2_t  x3_t  m1_t  m2_t  m3_t  w1_t  w2_t  w3_t      
        y1_t  y2_t  y3_t  w_t x_t m_t y_t w_L2_t  x_L2_t  x_L1_t
        m_L2_t  m_L1_t  y_L2_t y_L1_t;
USEVARIABLES ARE id2 y1_t y2_t y3_t;
    MISSING ARE ALL (-99);
    CLUSTER IS id2;              !identify clustering variable
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
    TYPE = TWOLEVEL;
	  PROCESSORS = 16;
MODEL:
%WITHIN%
!Measurement model
yt by y1_t y2_t  y3_t;  
%BETWEEN%
!Measurement model
yt2 by y1_t y2_t  y3_t;  

OUTPUT:
  tech1 tech8 sampstat; 
 Savedata: 
 file is f_y_t.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_y_model_ind,file="croon_y_model_ind.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_y_model_ind.mplus.inp")
### Save model
y_model_ind <-readModels(target="croon_y_model_ind.mplus.out")
#######################################
##factor model for yind_c
#######################################
croon_y_model_ind_c<-'
TITLE: Croon y Model_ind_c
DATA: FILE IS d_c.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0_c a_c  covar_coef b1_c x_c m_c y_c
        x1_c  x2_c  x3_c  m1_c  m2_c  m3_c  y1_c  y2_c  y3_c      
        x_L1_c m_L1_c y_L1_c;
USEVARIABLES ARE y1_c  y2_c  y3_c;
    MISSING ARE ALL (-99);
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
	  PROCESSORS = 16;
MODEL:
!Measurement model
yc by y1_c y2_c  y3_c;
OUTPUT:
  tech1 tech8 sampstat; 
 Savedata: 
 file is f_y_t.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_y_model_ind_c,file="croon_y_model_ind_c.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_y_model_ind_c.mplus.inp")
### Save model
y_model_ind_c <-readModels(target="croon_y_model_ind_c.mplus.out")

### Save results from Mplus------------------------------------------------------------------------
#L1
#factor scores
d_c$fyL1_c<-try(y_model_ind_c$savedata[,"YC"],silent = T)
d_t$fyL1<-try(y_model_ind$savedata[,"YT"],silent = T)
#lambda- factor loadings
LyL1_c<-y_model_ind_c$parameters$unstandardized[1:3,"est"]   
LyL1<-y_model_full$parameters$unstandardized[1:3,"est"]   
#variance
vyL1_c<-y_model_ind_c$parameters$unstandardized[7,"est"]     
vyL1<-y_model_ind$parameters$unstandardized[4,"est"]     
#---------------------------------------------------------------------------------------------------------------------------

# yindcovL1_t<-try(cov(y_model_full$savedata[(nc+1):(nc*2),1:3]),silent = T)
 yindcovL1_t<-try(y_model_ind$sampstat$WITHIN$covariances,silent = T)
ifelse(class(yindcovL1_t)=="try-error"| is.null(yindcovL1_t)==T,
       yindcovL1_t<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),yindcovL1_t<-yindcovL1_t)
yindcovL1_t[upper.tri(yindcovL1_t)]<-t(yindcovL1_t)[upper.tri(yindcovL1_t)]
AyL1<-try(solve(yindcovL1_t) %*% LyL1  *vyL1,silent = T)                       

#yindcovL1_c<-try(cov(y_model_full$savedata[1:nc,4:6]),silent = T)
yindcovL1_c<-try(y_model_ind_c$sampstat$covariances,silent = T)
ifelse(class(yindcovL1_c)=="try-error"| is.null(yindcovL1_c)==T,
       yindcovL1_c<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),yindcovL1_c<-yindcovL1_c)
yindcovL1_c[upper.tri(yindcovL1_c)]<-t(yindcovL1_c)[upper.tri(yindcovL1_c)]
AyL1_c<-try(solve(yindcovL1_c) %*% LyL1_c  *vyL1_c,silent = T)
#---------------------------------------------------------------------------------------------------------------------------
#L2 and reliability 
#factor scores
fyL2<-try(y_model_ind$savedata[,"YT2"],silent = T)
dagg$fyL2<-try(fyL2[seq(1, length(fyL2), n1)],silent = T)                          
Ly<-y_model_ind$parameters$unstandardized[8:10,"est"]                 
vyL2<-y_model_ind$parameters$unstandardized[14,"est"]

yindcovL2<-try(y_model_ind$sampstat$BETWEEN$covariances,silent = T)
ifelse(class(yindcovL2)=="try-error"| is.null(yindcovL2)==T,
       yindcovL2<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),yindcovL2<-yindcovL2)
yindcovL2[upper.tri(yindcovL2)]<-t(yindcovL2)[upper.tri(yindcovL2)]
Ay<-try(solve(yindcovL2) %*% Ly *vyL2,silent = T)                                     
Ry<-try(yindcovL2%*%solve(yindcovL2+yindcovL1_t/n1),silent = T)                     
# check
# y_model_full$parameters$unstandardized

#######################################
##factor model for W
#######################################
croon_w_model_full<-'
TITLE: Croon W Model_full
DATA: FILE IS d_full.dat;
VARIANCES=NOCHECK;
VARIABLE:
        NAMES ARE id2 beta0 xi00 pi0  a b1 B  x1_t x2_t  x3_t  m1_t      
m2_t m3_t w1_t w2_t  w3_t  y1_t y2_t y3_t  w_t x_t m_t       
y_t   w_L2_t x_L2_t x_L1_t m_L2_t m_L1_t  y_L2_t  y_L1_t  beta0_c a_c b1_c      
x_c  m_c  y_c  x1_c  x2_c x3_c  m1_c  m2_c  m3_c  y1_c  y2_c      
y3_c x_L1_c  m_L1_c  y_L1_c n2  n1  nc  tx  rho rhox  rhom      
weight1  weight2  weight3 ind_rhox  ind_rho ind_rhom  covar_coef;
USEVARIABLES ARE w1_t w2_t w3_t;
    MISSING ARE ALL (-99);
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
	  PROCESSORS = 16;

MODEL:
!Measurement model
wt2  BY w1_t  
w2_t  w3_t (Bw2 Bw3);

OUTPUT:
  tech1 tech8 NOCHISQUARE; 
 Savedata: 
 file is f_w_full.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_w_model_full,file="croon_w_model_full.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_w_model_full.mplus.inp")
### Save model
w_model_full <-readModels(target="croon_w_model_full.mplus.out")
# check
# w_model_full$parameters$unstandardized
#######################################
##factor model for wind
#######################################
croon_w_model_ind<-'
TITLE: Croon w Model_ind
DATA: FILE IS d_t.dat;
VARIABLE:
        NAMES ARE id2 n2  n1  nc  tx  rho rhox  rhom  weight1 weight2 weight3 
        ind_rhox  ind_rho ind_rhom beta0  xi00  pi0 a covar_coef b1 B 
        x1_t  x2_t  x3_t  m1_t  m2_t  m3_t  w1_t  w2_t  w3_t      
        y1_t  y2_t  y3_t  w_t x_t m_t y_t w_L2_t  x_L2_t  x_L1_t
        m_L2_t  m_L1_t  y_L2_t y_L1_t;
USEVARIABLES ARE w1_t w2_t w3_t;
    MISSING ARE ALL (-99);
    BETWEEN IS ;   !identify between-cluster (level-2) variables
ANALYSIS:
    ESTIMATOR = ;
	  PROCESSORS = 16;
MODEL:
!Measurement model
wt by w1_t w2_t  w3_t;  
OUTPUT:
  tech1 tech8 sampstat; 
 Savedata: 
 file is f_w_t.dat; 
 save = fscores; 
'
### Save for Mplus- .inp
write.table(croon_w_model_ind,file="croon_w_model_ind.mplus.inp",row.names=F,col.names=F,quote=F)
### Run analyis in Mplus- RUN .inp
runModels(target="croon_w_model_ind.mplus.inp")
### Save model
w_model_ind <-readModels(target="croon_w_model_ind.mplus.out")

### Save results from Mplus------------------------------------------------------------------------
#---------------------------------------------------------------------------------------------------------------------------
#L2 and reliability 
#factor scores
fwL2<-try(w_model_ind$savedata[,"WT"],silent = T)
dagg$fwL2<-try(fwL2[seq(1, length(fwL2), n1)],silent = T)                          
Lw<-w_model_ind$parameters$unstandardized[1:3,"est"]                 
vwL2<-w_model_ind$parameters$unstandardized[7,"est"]

windcovL2<-try(w_model_ind$sampstat$covariances,silent = T)
ifelse(class(windcovL2)=="try-error"| is.null(windcovL2)==T,
       windcovL2<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),windcovL2<-windcovL2)
windcovL2[upper.tri(windcovL2)]<-t(windcovL2)[upper.tri(windcovL2)]
Aw<-try(solve(windcovL2) %*% Lw *vwL2,silent = T)                                     
# check
# w_model_t$parameters$unstandardized

##############################################################################
### Level 1 Corrections-TREATMENT
##############################################################################
#unadjusted 
#construct cov matrix of relevant L1 variables and interactions
fs.cov.L1<-try((((n2*n1)-1)/(n2*n1))*cov(cbind(d_t$fyL1,d_t$fxL1,d_t$fmL1)),silent = T)
ifelse(class(fs.cov.L1)=="try-error"| is.null(fs.cov.L1)==T,
       fs.cov.L1<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),fs.cov.L1<-fs.cov.L1)

fs_cov_un.L1<-fs.cov.L1

fs_cov_un.L1[try(upper.tri(fs_cov_un.L1))]<-""
fs_cov_un.L1<-try(as.data.frame(fs_cov_un.L1),silent = T)

# Save uncorrected covariance as .dat
try(write.table(fs_cov_un.L1,file="cov1unL1.dat",row.names=F,col.names=F,quote=F),silent = T)

# croon corrected cov matrix
fsr.cov.L1<-fs.cov.L1
fsr.cov.L1[1,]<-try(fsr.cov.L1[1,]/c(t(AyL1) %*% LyL1),silent = T)
fsr.cov.L1[,1]<-try(fsr.cov.L1[,1]/c(t(AyL1) %*% LyL1),silent = T)
fsr.cov.L1[2,]<-try(fsr.cov.L1[2,]/c(t(AxL1) %*% LxL1),silent = T)
fsr.cov.L1[,2]<-try(fsr.cov.L1[,2]/c(t(AxL1) %*% LxL1),silent = T)
fsr.cov.L1[3,]<-try(fsr.cov.L1[3,]/c(t(AmL1) %*% LmL1),silent = T)
fsr.cov.L1[,3]<-try(fsr.cov.L1[,3]/c(t(AmL1) %*% LmL1),silent = T)

#corrected covs
diag(fsr.cov.L1)<-try(c(vyL1,vxL1,vmL1),silent = T)


ifelse(class(fsr.cov.L1)=="try-error"| is.null(fsr.cov.L1)==T,
       fsr.cov.L1<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),fsr.cov.L1<-fsr.cov.L1)

fsr.cov.L1[try(upper.tri(fsr.cov.L1))]<-""
fsr.cov.L1<-try(as.data.frame(fsr.cov.L1),silent = T)

write.table(fsr.cov.L1,file="cov1L1.dat",row.names=F,col.names=F,quote=F)


##############################################################################
### Level 2 Corrections-TREATMENT
##############################################################################
#unadjusted 
fs.cov.L2<-try((((n2)-1)/(n2))*cov(cbind(dagg$fyL2,dagg$fxL2,dagg$fmL2,dagg$fwL2)),silent = T)
ifelse(class(fs.cov.L2)=="try-error"| is.null(fs.cov.L2)==T,
       fs.cov.L2<-matrix(data=NA,nrow = 4, ncol = 4, byrow = FALSE,dimnames = NULL),fs.cov.L2<-fs.cov.L2)

fs_cov_un.L2<-fs.cov.L2

fs_cov_un.L2[try(upper.tri(fs_cov_un.L2))]<-""
fs_cov_un.L2<-try(as.data.frame(fs_cov_un.L2))

# Save uncorrected covariance as .dat
try(write.table(fs_cov_un.L2,file="cov1unL2.dat",row.names=F,col.names=F,quote=F))

# corrected
fsr.cov.L2<-fs.cov.L2

fsr.cov.L2[1,]<-try(fsr.cov.L2[1,]/c(t(Ay) %*% (Ry%*%Ly)),silent = T)
fsr.cov.L2[,1]<-try(fsr.cov.L2[,1]/c(t(Ay) %*% (Ry%*%Ly)),silent = T)
fsr.cov.L2[2,]<-try(fsr.cov.L2[2,]/c(t(Ax) %*% (Rx %*% Lx)),silent = T)
fsr.cov.L2[,2]<-try(fsr.cov.L2[,2]/c(t(Ax) %*% (Rx %*% Lx)),silent = T)
fsr.cov.L2[3,]<-try(fsr.cov.L2[3,]/c(t(Am) %*% (Rm %*% Lm)),silent = T)
fsr.cov.L2[,3]<-try(fsr.cov.L2[,3]/c(t(Am) %*% (Rm %*% Lm)),silent = T)
fsr.cov.L2[4,]<-try(fsr.cov.L2[4,]/c(t(Aw) %*% (Lw)),silent = T)
fsr.cov.L2[,4]<-try(fsr.cov.L2[,4]/c(t(Aw) %*% (Lw)),silent = T)

  diag(fsr.cov.L2)<-try(c(vyL2,vxL2,vmL2,vwL2),silent = T)

ifelse(class(fsr.cov.L2)=="try-error"| is.null(fsr.cov.L2)==T,
       fsr.cov.L2<-matrix(data=NA,nrow = 4, ncol = 4, byrow = FALSE,dimnames = NULL),fsr.cov.L2<-fsr.cov.L2)

fsr.cov.L2[try(upper.tri(fsr.cov.L2))]<-""
fsr.cov.L2<-try(as.data.frame(fsr.cov.L2),silent = T)

write.table(fsr.cov.L2,file="cov1L2.dat",row.names=F,col.names=F,quote=F)

##############################################################################
### Level 1 (individual cluster) Corrections-CONTROL
##############################################################################
#unadjusted 
#construct cov matrix of relevant L1 variables and interactions
fs.cov.L1_c<-try((((n2*n1)-1)/(n2*n1))*cov(cbind(d_c$fyL1_c,d_c$fxL1_c,d_c$fmL1_c)),silent = T)
ifelse(class(fs.cov.L1_c)=="try-error"| is.null(fs.cov.L1_c)==T,
       fs.cov.L1_c<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),fs.cov.L1_c<-fs.cov.L1_c)

fs_cov_un.L1_c<-fs.cov.L1_c

fs_cov_un.L1_c[try(upper.tri(fs_cov_un.L1_c))]<-""
fs_cov_un.L1_c<-try(as.data.frame(fs_cov_un.L1_c),silent = T)

# Save uncorrected covariance as .dat
try(write.table(fs_cov_un.L1_c,file="cov1unL1_c.dat",row.names=F,col.names=F,quote=F),silent = T)

# croon corrected cov matrix
fsr.cov.L1_c<-fs.cov.L1_c
fsr.cov.L1_c[1,]<-try(fsr.cov.L1_c[1,]/c(t(AyL1_c) %*% LyL1_c),silent = T)
fsr.cov.L1_c[,1]<-try(fsr.cov.L1_c[,1]/c(t(AyL1_c) %*% LyL1_c),silent = T)
fsr.cov.L1_c[2,]<-try(fsr.cov.L1_c[2,]/c(t(AxL1_c) %*% LxL1_c),silent = T)
fsr.cov.L1_c[,2]<-try(fsr.cov.L1_c[,2]/c(t(AxL1_c) %*% LxL1_c),silent = T)
fsr.cov.L1_c[3,]<-try(fsr.cov.L1_c[3,]/c(t(AmL1_c) %*% LmL1_c),silent = T)
fsr.cov.L1_c[,3]<-try(fsr.cov.L1_c[,3]/c(t(AmL1_c) %*% LmL1_c),silent = T)

#corrected covs
diag(fsr.cov.L1_c)<-try(c(vyL1_c,vxL1_c,vmL1_c),silent = T)

ifelse(class(fsr.cov.L1_c)=="try-error"| is.null(fsr.cov.L1_c)==T,
       fsr.cov.L1_c<-matrix(data=NA,nrow = 3, ncol = 3, byrow = FALSE,dimnames = NULL),fsr.cov.L1_c<-fsr.cov.L1_c)

fsr.cov.L1_c[try(upper.tri(fsr.cov.L1_c))]<-""
fsr.cov.L1_c<-try(as.data.frame(fsr.cov.L1_c),silent = T)

write.table(fsr.cov.L1_c,file="cov1L1_c.dat",row.names=F,col.names=F,quote=F)


##############################################################################
### Level 1 Path Models (uncorrected and corrected)-TREATMENT
##############################################################################
#############################################################
### Path analysis with uncorrected covariance-FS-TREATMENT                             
#############################################################

##path model specification for mplus 
fs_path_L1_t<-paste('
TITLE: fs path model_t
DATA:
    file is cov1unL1.dat;
    type is covariance;
    nobservations =', (n1*n2),';
VARIABLE:
    names are fyL1 fxL1 fmL1;
ANALYSIS:
    type = general;
    estimator = ML;
    PROCESSORS = 16;
MODEL:
fmL1 on fxL1;
fyL1 on fmL1 fxL1;

OUTPUT:
  sampstat tech1 tech3;
', sep="")

### Save for Mplus- .inp
write.table(fs_path_L1_t,file="fs_path_L1_t.mplus.inp",row.names=F,col.names=F,quote=F)

### Run uncorrected path model in Mplus
runModels(target="fs_path_L1_t.mplus.inp")
### Save corrected path model from Mplus
fs_path_L1fit_t<-readModels(target="fs_path_L1_t.mplus.out")

### Check results
# fs_path_L1fit_t$parameters$unstandardized

#############################################################
### Path analysis with corrected covariance- Croon-TREATMENT                             
#############################################################
#path model specification for mplus
croon_path_L1_t<-paste('
TITLE: croon path model_t
DATA:
    file is cov1L1.dat;
    type is covariance;
    nobservations =', (n1*n2),';
VARIABLE:
    names are fyL1 fxL1 fmL1;
ANALYSIS:
    type = general;
    estimator = ML;
    PROCESSORS = 16;
MODEL:
fmL1 on fxL1;
fyL1 on fmL1 fxL1;

OUTPUT:
 sampstat tech1 tech3;
', sep="")
### Save for Mplus- .inp
write.table(croon_path_L1_t,file="croon_path_L1_t.mplus.inp",row.names=F,col.names=F,quote=F)

### Run corrected path model in Mplus
runModels(target="croon_path_L1_t.mplus.inp")
### Save corrected path model from Mplus
croon_path_L1fit_t<-readModels(target="croon_path_L1_t.mplus.out")

### Check results
# croon_path_L1fit_t$parameters$unstandardized

##############################################################################
### Level 2 Path Models (uncorrected and corrected)-TREATMENT
##############################################################################
#############################################################
### Path analysis with uncorrected covariance-FS-TREATMENT                             
#############################################################

##path model specification for mplus 
fs_path_L2_t<-paste('
TITLE: fs path model_t
DATA:
    file is cov1unL2.dat;
    type is covariance;
    nobservations =', (n2),';
VARIABLE:
    names are fyL2 fxL2 fmL2 fwL2;
ANALYSIS:
    type = general;
    estimator = ML;
    PROCESSORS = 16;
MODEL:
fmL2 on fxL2 fwL2;
fyL2 on fmL2 fxL2 fwL2;

OUTPUT:
 sampstat tech1 tech3;
 
', sep="")

### Save for Mplus- .inp
write.table(fs_path_L2_t,file="fs_path_L2_t.mplus.inp",row.names=F,col.names=F,quote=F)

### Run uncorrected path model in Mplus
runModels(target="fs_path_L2_t.mplus.inp")
### Save corrected path model from Mplus
fs_path_L2fit_t<-readModels(target="fs_path_L2_t.mplus.out")

### Check results
# fs_path_L2fit_t$parameters$unstandardized

#############################################################
### Path analysis with corrected covariance- Croon-TREATMENT                             
#############################################################
#path model specification for mplus
croon_path_L2_t<-paste('
TITLE: croon path model_t
DATA:
    file is cov1L2.dat;
    type is covariance;
    nobservations =', (n2),';
VARIABLE:
    names are fyL2 fxL2 fmL2 fwL2;
ANALYSIS:
    type = general;
    estimator = ML;
    PROCESSORS = 16;
MODEL:
fmL2 on fxL2 fwL2;
fyL2 on fmL2 fxL2 fwL2;

OUTPUT:
  sampstat tech1 tech3;
', sep="")
### Save for Mplus- .inp
write.table(croon_path_L2_t,file="croon_path_L2_t.mplus.inp",row.names=F,col.names=F,quote=F)

### Run corrected path model in Mplus
runModels(target="croon_path_L2_t.mplus.inp")
### Save corrected path model from Mplus
croon_path_L2fit_t<-readModels(target="croon_path_L2_t.mplus.out")

### Check results
# croon_path_L2fit_t$parameters$unstandardized


##############################################################################
### Level 1 Path Models (uncorrected and corrected)-CONTROL
##############################################################################
#############################################################
### Path analysis with uncorrected covariance-FS-CONTROL                             
#############################################################
##path model specification for mplus 
fs_path_L1_c<-paste('
TITLE: fs path model_c
DATA:
    file is cov1unL1_c.dat;
    type is covariance;
    nobservations =', (nc),';
VARIABLE:
    names are fyL1_c fxL1_c fmL1_c;
ANALYSIS:
    type = general;
    estimator = ML;
    PROCESSORS = 16;
MODEL:
fmL1_c on fxL1_c;
fyL1_c on fmL1_c fxL1_c;

OUTPUT:
 sampstat tech1 tech3;
', sep="")

### Save for Mplus- .inp
write.table(fs_path_L1_c,file="fs_path_L1_c.mplus.inp",row.names=F,col.names=F,quote=F)

### Run uncorrected path model in Mplus
runModels(target="fs_path_L1_c.mplus.inp")
### Save corrected path model from Mplus
fs_path_L1fit_c<-readModels(target="fs_path_L1_c.mplus.out")

### Check results
# fs_path_L1fit_c$parameters$unstandardized

#############################################################
### Path analysis with corrected covariance- Croon-CONTROL                             
#############################################################
#path model specification for mplus
croon_path_L1_c<-paste('
TITLE: croon path model_c
DATA:
    file is cov1L1_c.dat;
    type is covariance;
    nobservations =', (nc),';
VARIABLE:
    names are fyL1_c fxL1_c fmL1_c;
ANALYSIS:
    type = general;
    estimator = ML;
    PROCESSORS = 16;
MODEL:
fmL1_c on fxL1_c;
fyL1_c on fmL1_c fxL1_c;

OUTPUT:
  sampstat tech1 tech3;
', sep="")
### Save for Mplus- .inp
write.table(croon_path_L1_c,file="croon_path_L1_c.mplus.inp",row.names=F,col.names=F,quote=F)

### Run corrected path model in Mplus
runModels(target="croon_path_L1_c.mplus.inp")
### Save corrected path model from Mplus
croon_path_L1fit_c<-readModels(target="croon_path_L1_c.mplus.out")
### Check results
# croon_path_L1fit_c$parameters$unstandardized

################################################################################################
### Save Simulation Run Results ###-------------------------------------------------------------
################################################################################################
#############################################################
### Convergence                                               
#############################################################
# True --------------------------------------------------------------------------------------------------------------
ifelse(length(mpfit_true_full$errors[]) == 0,true_conv_full_run <-0, true_conv_full_run<-1)
ifelse(true_conv_full_run==1, true_error_full_run<-paste(mpfit_true_full$errors, sep = '', collapse = ''),true_error_full_run<-NA)
true_conv_full<-rbind(true_conv_full,true_conv_full_run)  
true_error_full<-rbind(true_error_full,true_error_full_run) 

# ML --------------------------------------------------------------------------------------------------------------
ml_prob <- try(grep("PROBLEM"  , mpfit_ml_full$output), silent = T)

ifelse(length(mpfit_ml_full$errors[]) == 0 & length(ml_prob)==0,ml_conv_full_run <-0, ml_conv_full_run<-1)
ifelse(ml_conv_full_run==1, ml_error_full_run<-paste(mpfit_ml_full$errors, sep = '', collapse = ''),ml_error_full_run<-NA)
ml_conv_full<-rbind(ml_conv_full,ml_conv_full_run)  
ml_error_full<-rbind(ml_error_full,ml_error_full_run)

#chron objects store the values internally as a fraction of seconds per day. 
#Thus 1 second is equivalent to 1/(60*60*24), or 1/86400, i.e. 1.157407e-05.
et <- try(grep("Elapsed Time"  , mpfit_ml_full$output), silent = T)
ml_time_run <- try(mpfit_ml_full$output[et], silent = T)
ml_time_run<-try(chron(times= str_sub(ml_time_run,start= -8)), silent = T)
ml_time_full<-rbind(ml_time_full,ml_time_run)
## Bayes --------------------------------------------------------------------------------------------------------------
ifelse(length(mpfit_bayes_full$errors[]) == 0,bayes_conv_full_run <-0, bayes_conv_full_run<-1)
ifelse(bayes_conv_full_run==1, bayes_error_full_run<-paste(mpfit_bayes_full$errors, sep = '', collapse = ''),bayes_error_full_run<-NA)
bayes_conv_full<-rbind(bayes_conv_full,bayes_conv_full_run)  
bayes_error_full<-rbind(bayes_error_full,bayes_error_full_run)

bayes_psr_full_run <-try(mpfit_bayes_full$tech8$psr[length(mpfit_bayes_full$tech8$psr[,2]),2],silent = T)
ifelse(class(bayes_psr_full_run)=="try-error",bayes_psr_full_run<-NA, bayes_psr_full_run<-bayes_psr_full_run)
bayes_psr_full <- rbind(bayes_psr_full,bayes_psr_full_run)

et <- try(grep("Elapsed Time"  , mpfit_bayes_full$output), silent = T)
bayes_time_run <- try(mpfit_bayes_full$output[et], silent = T)
bayes_time_run<-try(chron(times= str_sub(bayes_time_run,start= -8)), silent = T)
bayes_time_full<-rbind(bayes_time_full,bayes_time_run)
## BayesIN --------------------------------------------------------------------------------------------------------------
ifelse(length(mpfit_bayesIN_full$errors[]) == 0,bayesIN_conv_full_run <-0, bayesIN_conv_full_run<-1)
ifelse(bayesIN_conv_full_run==1, bayesIN_error_full_run<-paste(mpfit_bayesIN_full$errors, sep = '', collapse = ''),bayesIN_error_full_run<-NA)
bayesIN_conv_full<-rbind(bayesIN_conv_full,bayesIN_conv_full_run)  
bayesIN_error_full<-rbind(bayesIN_error_full,bayesIN_error_full_run)

bayesIN_psr_full_run <-try(mpfit_bayesIN_full$tech8$psr[length(mpfit_bayesIN_full$tech8$psr[,2]),2],silent = T)
ifelse(class(bayesIN_psr_full_run)=="try-error",bayesIN_psr_full_run<-NA, bayesIN_psr_full_run<-bayesIN_psr_full_run)
bayesIN_psr_full <- rbind(bayesIN_psr_full,bayesIN_psr_full_run)

et <- try(grep("Elapsed Time"  , mpfit_bayesIN_full$output), silent = T)
bayesIN_time_run <- try(mpfit_bayesIN_full$output[et], silent = T)
bayesIN_time_run<-try(chron(times= str_sub(bayesIN_time_run,start= -8)), silent = T)
bayesIN_time_full<-rbind(bayesIN_time_full,bayesIN_time_run)

# FS TREATMENT--------------------------------------------------------------------------------------------------------------
ifelse(length(fs_path_L1fit_t$errors[]) == 0,fs_convL1_t_run <-0, fs_convL1_t_run<-1)
ifelse(fs_convL1_t_run==1, fs_errorL1_t_run<-paste(fs_path_L1fit_t$errors, sep = '', collapse = ''),fs_errorL1_t_run<-NA)
fs_convL1_t<-rbind(fs_convL1_t,fs_convL1_t_run)  
fs_errorL1_t<-rbind(fs_errorL1_t,fs_errorL1_t_run)

ifelse(length(fs_path_L2fit_t$errors[]) == 0,fs_conv_t_run <-0, fs_conv_t_run<-1)
ifelse(fs_conv_t_run==1, fs_error_t_run<-paste(fs_path_L2fit_t$errors, sep = '', collapse = ''),fs_error_t_run<-NA)
fs_conv_t<-rbind(fs_conv_t, fs_conv_t_run)  
fs_error_t<-rbind(fs_error_t,fs_error_t_run) 

# FS CONTROL--------------------------------------------------------------------------------------------------------------
ifelse(length(fs_path_L1fit_c$errors[]) == 0,fs_convL1_c_run <-0, fs_convL1_c_run<-1)
ifelse(fs_convL1_c_run==1, fs_convL1_c_run<-paste(fs_path_L1fit_c$errors, sep = '', collapse = ''),fs_errorL1_c_run<-NA)
fs_convL1_c<-rbind(fs_convL1_c,fs_convL1_c_run)  
fs_errorL1_c<-rbind(fs_errorL1_c,fs_errorL1_c_run)


# CROON TREATMENT--------------------------------------------------------------------------------------------------------------
ifelse(length(croon_path_L1fit_t$errors[]) == 0,croon_convL1_t_run <-0, croon_convL1_t_run<-1)
ifelse(croon_convL1_t_run==1, croon_errorL1_t_run<-paste(croon_path_L1fit_t$errors, sep = '', collapse = ''),croon_errorL1_t_run<-NA)
croon_convL1_t<-rbind(croon_convL1_t,croon_convL1_t_run)  
croon_errorL1_t<-rbind(croon_errorL1_t,croon_errorL1_t_run)

ifelse(length(croon_path_L2fit_t$errors[]) == 0,croon_conv_t_run <-0, croon_conv_t_run<-1)
ifelse(croon_conv_t_run==1, croon_error_t_run<-paste(croon_path_L2fit_t$errors, sep = '', collapse = ''),croon_error_t_run<-NA)
croon_conv_t<-rbind(croon_conv_t,croon_conv_t_run)  
croon_error_t<-rbind(croon_error_t,croon_error_t_run)

# CROON CONTROL--------------------------------------------------------------------------------------------------------------
ifelse(length(croon_path_L1fit_c$errors[]) == 0,croon_convL1_c_run <-0, croon_convL1_c_run<-1)
ifelse(croon_convL1_c_run==1, croon_convL1_c_run<-paste(croon_path_L1fit_c$errors, sep = '', collapse = ''),croon_errorL1_c_run<-NA)
croon_convL1_c<-rbind(croon_convL1_c,croon_convL1_c_run)  
croon_errorL1_c<-rbind(croon_errorL1_c,croon_errorL1_c_run)


#############################################################
### Coefficient Estimates                          
#############################################################
### For models that converged, save structural coefs

### TRUE ###
##L2
# True:Mplus
true_coefs_t_run <-try(c(mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"& 
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="W_T"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="W_T"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="Y_T"]), silent = T)
ifelse(class(true_coefs_t_run)=="try-error"| is.null(true_coefs_t_run)==T,true_coefs_t_run<-c(NA,NA,NA,NA,NA,NA,NA),
       true_coefs_t_run<-c(mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"& 
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="W_T"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="W_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="Y_T"]))


##L1
# True:Mplus
true_coefsL1_t_run <-try(c(mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"& 
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"]), silent = T)
ifelse(class(true_coefsL1_t_run)=="try-error"| is.null(true_coefsL1_t_run)==T,true_coefsL1_t_run<-c(NA,NA,NA),
       true_coefsL1_t_run<-c(mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"& 
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_T"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_T.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_T"]))

##L1
# True:control
true_coefsL1_c_run <-try(c(mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_C.ON"& 
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_C"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_C.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_C"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_C.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_C"]), silent = T)
ifelse(class(true_coefsL1_c_run)=="try-error"| is.null(true_coefsL1_c_run)==T,true_coefsL1_c_run<-c(NA,NA,NA),
       true_coefsL1_c_run<-c(mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_C.ON"& 
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="M_C"],
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="Y_C.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_C"],
  
                        mpfit_true_full$parameters$unstandardized$est[mpfit_true_full$parameters$unstandardized$paramHeader=="M_C.ON"&
                                                                    mpfit_true_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_true_full$parameters$unstandardized$param=="X_C"]))



#main and mediation effects
true_effects_full_run <-try(c(true_coefs_t_run[6]*true_coefs_t_run[1],true_coefs_t_run[7]), silent = T)
ifelse(class(true_effects_full_run)=="try-error"| is.null(true_effects_full_run)==T,true_effects_full_run<-c(NA,NA),
true_effects_full_run<-c(true_coefs_t_run[6]*true_coefs_t_run[1],true_coefs_t_run[7]))
  
true_coefs_t_run<-c(true_coefs_t_run,true_effects_full_run) 
true_coefs_t<-rbind(true_coefs_t,true_coefs_t_run)

true_coefsL1_t<-rbind(true_coefsL1_t,true_coefsL1_t_run)

true_coefsL1_c<-rbind(true_coefsL1_c,true_coefsL1_c_run)

### ml ###
##L2
# ml:Mplus
ml_coefs_t_run <-try(c(mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT2.ON"& 
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="YT2"]), silent = T)
ifelse(class(ml_coefs_t_run)=="try-error"| is.null(ml_coefs_t_run)==T,ml_coefs_t_run<-c(NA,NA,NA,NA,NA,NA,NA),
       ml_coefs_t_run<-c(mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT2.ON"& 
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MT2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="YT2"]))


##L1
# ml:Mplus
ml_coefsL1_t_run <-try(c(mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT1.ON"& 
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MT1"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT1.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT1"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MT1.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT1"]), silent = T)
ifelse(class(ml_coefsL1_t_run)=="try-error"| is.null(ml_coefsL1_t_run)==T,ml_coefsL1_t_run<-c(NA,NA,NA),
       ml_coefsL1_t_run<-c(mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT1.ON"& 
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MT1"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YT1.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT1"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MT1.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XT1"]))

##L1
# ml:control
ml_coefsL1_c_run <-try(c(mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YC2.ON"& 
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MC2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YC2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XC2"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MC2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XC2"]), silent = T)
ifelse(class(ml_coefsL1_c_run)=="try-error"| is.null(ml_coefsL1_c_run)==T,ml_coefsL1_c_run<-c(NA,NA,NA),
       ml_coefsL1_c_run<-c(mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YC2.ON"& 
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="MC2"],
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="YC2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XC2"],
  
                        mpfit_ml_full$parameters$unstandardized$est[mpfit_ml_full$parameters$unstandardized$paramHeader=="MC2.ON"&
                                                                    mpfit_ml_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_ml_full$parameters$unstandardized$param=="XC2"]))



#main and mediation effects
ml_effects_full_run <-try(c(ml_coefs_t_run[6]*ml_coefs_t_run[1],ml_coefs_t_run[7]), silent = T)
ifelse(class(ml_effects_full_run)=="try-error"| is.null(ml_effects_full_run)==T,ml_effects_full_run<-c(NA,NA),
ml_effects_full_run<-c(ml_coefs_t_run[6]*ml_coefs_t_run[1],ml_coefs_t_run[7]))
  
ml_coefs_t_run<-c(ml_coefs_t_run,ml_effects_full_run) 
ml_coefs_t<-rbind(ml_coefs_t,ml_coefs_t_run)

ml_coefsL1_t<-rbind(ml_coefsL1_t,ml_coefsL1_t_run)

ml_coefsL1_c<-rbind(ml_coefsL1_c,ml_coefsL1_c_run)

### bayes ###
##L2
# bayes:Mplus
bayes_coefs_t_run <-try(c(mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT2.ON"& 
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="YT2"]), silent = T)
ifelse(class(bayes_coefs_t_run)=="try-error"| is.null(bayes_coefs_t_run)==T,bayes_coefs_t_run<-c(NA,NA,NA,NA,NA,NA,NA),
       bayes_coefs_t_run<-c(mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT2.ON"& 
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="YT2"]))


##L1
# bayes:Mplus
bayes_coefsL1_t_run <-try(c(mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT1.ON"& 
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MT1"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT1.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT1"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MT1.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT1"]), silent = T)
ifelse(class(bayes_coefsL1_t_run)=="try-error"| is.null(bayes_coefsL1_t_run)==T,bayes_coefsL1_t_run<-c(NA,NA,NA),
       bayes_coefsL1_t_run<-c(mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT1.ON"& 
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MT1"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YT1.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT1"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MT1.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XT1"]))

##L1
# bayes:control
bayes_coefsL1_c_run <-try(c(mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YC2.ON"& 
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MC2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YC2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XC2"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MC2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XC2"]), silent = T)
ifelse(class(bayes_coefsL1_c_run)=="try-error"| is.null(bayes_coefsL1_c_run)==T,bayes_coefsL1_c_run<-c(NA,NA,NA),
       bayes_coefsL1_c_run<-c(mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YC2.ON"& 
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="MC2"],
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="YC2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XC2"],
  
                        mpfit_bayes_full$parameters$unstandardized$est[mpfit_bayes_full$parameters$unstandardized$paramHeader=="MC2.ON"&
                                                                    mpfit_bayes_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayes_full$parameters$unstandardized$param=="XC2"]))



#main and mediation effects
bayes_effects_full_run <-try(c(bayes_coefs_t_run[6]*bayes_coefs_t_run[1],bayes_coefs_t_run[7]), silent = T)
ifelse(class(bayes_effects_full_run)=="try-error"| is.null(bayes_effects_full_run)==T,bayes_effects_full_run<-c(NA,NA),
bayes_effects_full_run<-c(bayes_coefs_t_run[6]*bayes_coefs_t_run[1],bayes_coefs_t_run[7]))
  
bayes_coefs_t_run<-c(bayes_coefs_t_run,bayes_effects_full_run) 
bayes_coefs_t<-rbind(bayes_coefs_t,bayes_coefs_t_run)

bayes_coefsL1_t<-rbind(bayes_coefsL1_t,bayes_coefsL1_t_run)

bayes_coefsL1_c<-rbind(bayes_coefsL1_c,bayes_coefsL1_c_run)


### bayesIN ###
##L2
# bayesIN:Mplus
bayesIN_coefs_t_run <-try(c(mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT2.ON"& 
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="YT2"]), silent = T)
ifelse(class(bayesIN_coefs_t_run)=="try-error"| is.null(bayesIN_coefs_t_run)==T,bayesIN_coefs_t_run<-c(NA,NA,NA,NA,NA,NA,NA),
       bayesIN_coefs_t_run<-c(mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT2.ON"& 
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MT2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="WT2"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MT2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="Intercepts"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="YT2"]))


##L1
# bayesIN:Mplus
bayesIN_coefsL1_t_run <-try(c(mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT1.ON"& 
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MT1"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT1.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT1"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MT1.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT1"]), silent = T)
ifelse(class(bayesIN_coefsL1_t_run)=="try-error"| is.null(bayesIN_coefsL1_t_run)==T,bayesIN_coefsL1_t_run<-c(NA,NA,NA),
       bayesIN_coefsL1_t_run<-c(mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT1.ON"& 
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MT1"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YT1.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT1"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MT1.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Within"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XT1"]))

##L1
# bayesIN:control
bayesIN_coefsL1_c_run <-try(c(mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YC2.ON"& 
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MC2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YC2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XC2"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MC2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XC2"]), silent = T)
ifelse(class(bayesIN_coefsL1_c_run)=="try-error"| is.null(bayesIN_coefsL1_c_run)==T,bayesIN_coefsL1_c_run<-c(NA,NA,NA),
       bayesIN_coefsL1_c_run<-c(mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YC2.ON"& 
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="MC2"],
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="YC2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XC2"],
  
                        mpfit_bayesIN_full$parameters$unstandardized$est[mpfit_bayesIN_full$parameters$unstandardized$paramHeader=="MC2.ON"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    mpfit_bayesIN_full$parameters$unstandardized$param=="XC2"]))



#main and mediation effects
bayesIN_effects_full_run <-try(c(bayesIN_coefs_t_run[6]*bayesIN_coefs_t_run[1],bayesIN_coefs_t_run[7]), silent = T)
ifelse(class(bayesIN_effects_full_run)=="try-error"| is.null(bayesIN_effects_full_run)==T,bayesIN_effects_full_run<-c(NA,NA),
bayesIN_effects_full_run<-c(bayesIN_coefs_t_run[6]*bayesIN_coefs_t_run[1],bayesIN_coefs_t_run[7]))
  
bayesIN_coefs_t_run<-c(bayesIN_coefs_t_run,bayesIN_effects_full_run) 
bayesIN_coefs_t<-rbind(bayesIN_coefs_t,bayesIN_coefs_t_run)

bayesIN_coefsL1_t<-rbind(bayesIN_coefsL1_t,bayesIN_coefsL1_t_run)

bayesIN_coefsL1_c<-rbind(bayesIN_coefsL1_c,bayesIN_coefsL1_c_run)


### fs in Mplus
# treat between
fs_coefs_t_run <-try(c(fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FML2"],
                     fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
                    
                    fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
  
                    m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        m_model_full$parameters$unstandardized$BetweenWithin=="Between"& m_model_full$parameters$unstandardized$param=="MT2"],
                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        y_model_full$parameters$unstandardized$BetweenWithin=="Between"& y_model_full$parameters$unstandardized$param=="YT2"]),silent = T)
ifelse(class(fs_coefs_t_run)=="try-error"| is.null(fs_coefs_t_run)==T|class(fs_coefs_t_run)=="character"|length(fs_coefs_t_run)!=7,fs_coefs_t_run<-c(NA,NA,NA,NA,NA,NA,NA),
       fs_coefs_t_run<-c(fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FML2"],
                     fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
                    
                    fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    fs_path_L2fit_t$parameters$unstandardized$est[fs_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& fs_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
  
                    m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        m_model_full$parameters$unstandardized$BetweenWithin=="Between"& m_model_full$parameters$unstandardized$param=="MT2"],
                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        y_model_full$parameters$unstandardized$BetweenWithin=="Between"& y_model_full$parameters$unstandardized$param=="YT2"]))

# treat within
fs_coefsL1_t_run <-try(c(fs_path_L1fit_t$parameters$unstandardized$est[fs_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& fs_path_L1fit_t$parameters$unstandardized$param=="FML1"],
                    fs_path_L1fit_t$parameters$unstandardized$est[fs_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& fs_path_L1fit_t$parameters$unstandardized$param=="FXL1"],
                    
                    fs_path_L1fit_t$parameters$unstandardized$est[fs_path_L1fit_t$parameters$unstandardized$paramHeader=="FML1.ON"& fs_path_L1fit_t$parameters$unstandardized$param=="FXL1"]), silent = T)
ifelse(class(fs_coefsL1_t_run)=="try-error"| is.null(fs_coefsL1_t_run)==T|class(fs_coefsL1_t_run)=="character",fs_coefsL1_t_run<-c(NA,NA,NA),
       fs_coefsL1_t_run<-c(fs_path_L1fit_t$parameters$unstandardized$est[fs_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& fs_path_L1fit_t$parameters$unstandardized$param=="FML1"],
                    fs_path_L1fit_t$parameters$unstandardized$est[fs_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& fs_path_L1fit_t$parameters$unstandardized$param=="FXL1"],
                    
                    fs_path_L1fit_t$parameters$unstandardized$est[fs_path_L1fit_t$parameters$unstandardized$paramHeader=="FML1.ON"& fs_path_L1fit_t$parameters$unstandardized$param=="FXL1"]))

# control within
fs_coefsL1_c_run <-try(c(fs_path_L1fit_c$parameters$unstandardized$est[fs_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& fs_path_L1fit_c$parameters$unstandardized$param=="FML1_C"],
                    fs_path_L1fit_c$parameters$unstandardized$est[fs_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& fs_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"],
                    
                    fs_path_L1fit_c$parameters$unstandardized$est[fs_path_L1fit_c$parameters$unstandardized$paramHeader=="FML1_C.ON"& fs_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"]), silent = T)
ifelse(class(fs_coefsL1_c_run)=="try-error"| is.null(fs_coefsL1_c_run)==T|class(fs_coefsL1_c_run)=="character",fs_coefsL1_c_run<-c(NA,NA,NA),
       fs_coefsL1_c_run<-c(fs_path_L1fit_c$parameters$unstandardized$est[fs_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& fs_path_L1fit_c$parameters$unstandardized$param=="FML1_C"],
                    fs_path_L1fit_c$parameters$unstandardized$est[fs_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& fs_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"],
                    
                    fs_path_L1fit_c$parameters$unstandardized$est[fs_path_L1fit_c$parameters$unstandardized$paramHeader=="FML1_C.ON"& fs_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"]))

#main and mediation effects
fs_effects_full_run <-try(c(m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    m_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    m_model_full$parameters$unstandardized$param=="MT2"]*fs_coefs_t_run[1],
                                                                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    y_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    y_model_full$parameters$unstandardized$param=="YT2"]), silent = T)
ifelse(class(fs_effects_full_run)=="try-error"| is.null(fs_effects_full_run)==T,fs_effects_full_run<-c(NA,NA),
fs_effects_full_run<-c(m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    m_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    m_model_full$parameters$unstandardized$param=="MT2"]*fs_coefs_t_run[1],
                                                                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    y_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    y_model_full$parameters$unstandardized$param=="YT2"]))
  
fs_coefs_t_run<-c(fs_coefs_t_run,fs_effects_full_run) 
fs_coefs_t<-rbind(fs_coefs_t,fs_coefs_t_run)

fs_coefsL1_t<-rbind(fs_coefsL1_t,fs_coefsL1_t_run)

fs_coefsL1_c<-rbind(fs_coefsL1_c,fs_coefsL1_c_run)

### croon in Mplus
# treat between
croon_coefs_t_run <-try(c(croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FML2"],
                     croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
                    
                    croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
  
                    m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        m_model_full$parameters$unstandardized$BetweenWithin=="Between"& m_model_full$parameters$unstandardized$param=="MT2"],
                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        y_model_full$parameters$unstandardized$BetweenWithin=="Between"& y_model_full$parameters$unstandardized$param=="YT2"]), silent = T)
ifelse(class(croon_coefs_t_run)=="try-error"| is.null(croon_coefs_t_run)==T|class(croon_coefs_t_run)=="character"|length(croon_coefs_t_run)!=7,croon_coefs_t_run<-c(NA,NA,NA,NA,NA,NA,NA),
       croon_coefs_t_run<-c(croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FML2"],
                     croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FYL2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
                    
                    croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FXL2"],
                    croon_path_L2fit_t$parameters$unstandardized$est[croon_path_L2fit_t$parameters$unstandardized$paramHeader=="FML2.ON"& croon_path_L2fit_t$parameters$unstandardized$param=="FWL2"],
  
                    m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        m_model_full$parameters$unstandardized$BetweenWithin=="Between"& m_model_full$parameters$unstandardized$param=="MT2"],
                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"& 
                        y_model_full$parameters$unstandardized$BetweenWithin=="Between"& y_model_full$parameters$unstandardized$param=="YT2"]))

# treat within
croon_coefsL1_t_run <-try(c(croon_path_L1fit_t$parameters$unstandardized$est[croon_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& croon_path_L1fit_t$parameters$unstandardized$param=="FML1"],
                    croon_path_L1fit_t$parameters$unstandardized$est[croon_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& croon_path_L1fit_t$parameters$unstandardized$param=="FXL1"],
                    
                    croon_path_L1fit_t$parameters$unstandardized$est[croon_path_L1fit_t$parameters$unstandardized$paramHeader=="FML1.ON"& croon_path_L1fit_t$parameters$unstandardized$param=="FXL1"]), silent = T)
ifelse(class(croon_coefsL1_t_run)=="try-error"| is.null(croon_coefsL1_t_run)==T|class(croon_coefsL1_t_run)=="character",croon_coefsL1_t_run<-c(NA,NA,NA),
       croon_coefsL1_t_run<-c(croon_path_L1fit_t$parameters$unstandardized$est[croon_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& croon_path_L1fit_t$parameters$unstandardized$param=="FML1"],
                    croon_path_L1fit_t$parameters$unstandardized$est[croon_path_L1fit_t$parameters$unstandardized$paramHeader=="FYL1.ON"& croon_path_L1fit_t$parameters$unstandardized$param=="FXL1"],
                    
                    croon_path_L1fit_t$parameters$unstandardized$est[croon_path_L1fit_t$parameters$unstandardized$paramHeader=="FML1.ON"& croon_path_L1fit_t$parameters$unstandardized$param=="FXL1"]))

# control within
croon_coefsL1_c_run <-try(c(croon_path_L1fit_c$parameters$unstandardized$est[croon_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& croon_path_L1fit_c$parameters$unstandardized$param=="FML1_C"],
                    croon_path_L1fit_c$parameters$unstandardized$est[croon_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& croon_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"],
                    
                    croon_path_L1fit_c$parameters$unstandardized$est[croon_path_L1fit_c$parameters$unstandardized$paramHeader=="FML1_C.ON"& croon_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"]), silent = T)
ifelse(class(croon_coefsL1_c_run)=="try-error"| is.null(croon_coefsL1_c_run)==T|class(croon_coefsL1_c_run)=="character",croon_coefsL1_c_run<-c(NA,NA,NA),
       croon_coefsL1_c_run<-c(croon_path_L1fit_c$parameters$unstandardized$est[croon_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& croon_path_L1fit_c$parameters$unstandardized$param=="FML1_C"],
                    croon_path_L1fit_c$parameters$unstandardized$est[croon_path_L1fit_c$parameters$unstandardized$paramHeader=="FYL1_C.ON"& croon_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"],
                    
                    croon_path_L1fit_c$parameters$unstandardized$est[croon_path_L1fit_c$parameters$unstandardized$paramHeader=="FML1_C.ON"& croon_path_L1fit_c$parameters$unstandardized$param=="FXL1_C"]))

#main and mediation effects
croon_effects_full_run <-try(c(m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    m_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    m_model_full$parameters$unstandardized$param=="MT2"]*croon_coefs_t_run[1],
                                                                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    y_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    y_model_full$parameters$unstandardized$param=="YT2"]), silent = T)
ifelse(class(croon_effects_full_run)=="try-error"| is.null(croon_effects_full_run)==T,croon_effects_full_run<-c(NA,NA),
croon_effects_full_run<-c(m_model_full$parameters$unstandardized$est[m_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    m_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    m_model_full$parameters$unstandardized$param=="MT2"]*croon_coefs_t_run[1],
                                                                    y_model_full$parameters$unstandardized$est[y_model_full$parameters$unstandardized$paramHeader=="Means"&
                                                                    y_model_full$parameters$unstandardized$BetweenWithin=="Between"&
                                                                    y_model_full$parameters$unstandardized$param=="YT2"]))
  
croon_coefs_t_run<-c(croon_coefs_t_run,croon_effects_full_run) 

croon_coefs_t<-rbind(croon_coefs_t,croon_coefs_t_run)

croon_coefsL1_t<-rbind(croon_coefsL1_t,croon_coefsL1_t_run)

croon_coefsL1_c<-rbind(croon_coefsL1_c,croon_coefsL1_c_run)

plot(1,main=loop)
  } #end loop
#stopImplicitCluster()


#################################################################################################################################
###------------------------------------------------- Compile Results Across simulation ---------------------------------------###
#################################################################################################################################
#############################################################
### ## Converge                              
#############################################################
#Errors
error_names <-c("true_error_full",
                "ml_error_full",
                "bayes_error_full", 
                "bayesIN_error_full", 
                "fs_error_t", 
                "croon_error_t",
                "fs_errorL1_t", 
                "croon_errorL1_t", 
                "fs_errorL1_c", 
                "croon_errorL1_c")
error_list <- cbind(true_error_full[,1], 
                    ml_error_full[,1],
                    bayes_error_full[,1], 
                    bayesIN_error_full[,1], 
                    fs_error_t[,1], 
                    croon_error_t[,1], 
                    fs_errorL1_t[,1], 
                    croon_errorL1_t[,1], 
                    fs_errorL1_c[,1], 
                    croon_errorL1_c[,1])
error_list <-as.data.frame(error_list)
colnames(error_list)<-error_names
#Bayes psr
bayes_names <-c("bayes_psr_full", 
                "bayesIN_psr_full")
psr_list <- cbind(bayes_psr_full[,1], 
                    bayesIN_psr_full[,1])
psr_list <-as.data.frame(psr_list)
colnames(psr_list)<-bayes_names

# Percent convergence failure
true_conv_full<-try(mean(true_conv_full), silent = T)
ml_conv_full <-try(mean(ml_conv_full), silent = T)
bayes_conv_full <-try(mean(bayes_conv_full), silent = T)
bayesIN_conv_full <-try(mean(bayesIN_conv_full), silent = T)
fs_conv_t <-try(mean(fs_conv_t), silent = T)
fs_convL1_t<-try(mean(fs_convL1_t), silent = T)
croon_conv_t <-try(mean(croon_conv_t), silent = T)
croon_convL1_t<-try(mean(croon_convL1_t), silent = T)

fs_convL1_c <-try(mean(fs_convL1_c), silent = T)
croon_convL1_c <-try(mean(croon_convL1_c), silent = T)

# Time to converge
#chron objects store the values internally as a fraction of seconds per day. 
#Thus 1 second is equivalent to 1/(60*60*24), or 1/86400, i.e. 1.157407e-05.
ml_time_full_avg <-try(mean(ml_time_full*86400,na.rm=T), silent = T)
bayes_time_full_avg<-try(mean(bayes_time_full*86400,na.rm=T), silent = T)
bayesIN_time_full_avg<-try(mean(bayesIN_time_full*86400,na.rm=T), silent = T)

#############################################################
### ## Clean Results                              
#############################################################
#true_coefs_t_org <-true_coefs_t
#true_coefsL1_c_org <-true_coefsL1_c
true_coefs_t<-as.data.frame(true_coefs_t)
true_coefsL1_t<-as.data.frame(true_coefsL1_t)
true_coefsL1_c<-as.data.frame(true_coefsL1_c)
ml_coefs_t<-as.data.frame(ml_coefs_t)
ml_coefsL1_t<-as.data.frame(ml_coefsL1_t)
ml_coefsL1_c<-as.data.frame(ml_coefsL1_c)
bayes_coefs_t<-as.data.frame(bayes_coefs_t)
bayes_coefsL1_t<-as.data.frame(bayes_coefsL1_t)
bayes_coefsL1_c<-as.data.frame(bayes_coefsL1_c)
bayesIN_coefs_t<-as.data.frame(bayesIN_coefs_t)
bayesIN_coefsL1_t<-as.data.frame(bayesIN_coefsL1_t)
bayesIN_coefsL1_c<-as.data.frame(bayesIN_coefsL1_c)
fs_coefs_t<-as.data.frame(fs_coefs_t)
fs_coefsL1_t<-as.data.frame(fs_coefsL1_t)
fs_coefsL1_c<-as.data.frame(fs_coefsL1_c)
croon_coefs_t<-as.data.frame(croon_coefs_t)
croon_coefsL1_t<-as.data.frame(croon_coefsL1_t)
croon_coefsL1_c<-as.data.frame(croon_coefsL1_c)

true_coefs_t <- true_coefs_t %>% mutate_at(1:9, as.numeric)
true_coefsL1_t<-true_coefsL1_t %>% mutate_at(1:3, as.numeric)
true_coefsL1_c<-true_coefsL1_c %>% mutate_at(1:3, as.numeric)
ml_coefs_t <- ml_coefs_t %>% mutate_at(1:9, as.numeric)
ml_coefsL1_t<-ml_coefsL1_t %>% mutate_at(1:3, as.numeric)
ml_coefsL1_c<-ml_coefsL1_c %>% mutate_at(1:3, as.numeric)
bayes_coefs_t <- bayes_coefs_t %>% mutate_at(1:9, as.numeric)
bayes_coefsL1_t<-bayes_coefsL1_t %>% mutate_at(1:3, as.numeric)
bayes_coefsL1_c<-bayes_coefsL1_c %>% mutate_at(1:3, as.numeric)
bayesIN_coefs_t <- bayesIN_coefs_t %>% mutate_at(1:9, as.numeric)
bayesIN_coefsL1_t<-bayesIN_coefsL1_t %>% mutate_at(1:3, as.numeric)
bayesIN_coefsL1_c<-bayesIN_coefsL1_c %>% mutate_at(1:3, as.numeric)
fs_coefs_t <- fs_coefs_t %>% mutate_at(1:9, as.numeric)
fs_coefsL1_t<-fs_coefsL1_t %>% mutate_at(1:3, as.numeric)
fs_coefsL1_c<-fs_coefsL1_c %>% mutate_at(1:3, as.numeric)
croon_coefs_t <- croon_coefs_t %>% mutate_at(1:9, as.numeric)
croon_coefsL1_t<-croon_coefsL1_t %>% mutate_at(1:3, as.numeric)
croon_coefsL1_c<-croon_coefsL1_c %>% mutate_at(1:3, as.numeric)

#############################################################
### ## FS and Croon Main Effect                              
#############################################################
fs_coefs_t[,9]<-fs_coefs_t[,7]-fs_coefs_t[,1]*fs_coefs_t[,6]
croon_coefs_t[,9]<-croon_coefs_t[,7]-croon_coefs_t[,1]*croon_coefs_t[,6]

fs_coefs_t[,7]<-fs_coefs_t[,9]
croon_coefs_t[,7]<-croon_coefs_t[,9]

#############################################################
### ## Bias-                               
#############################################################
### Treatment
#L2 average value by parameter and estimator
true_coefs_t_avg<-round(colMeans( true_coefs_t,na.rm=T),3)
ml_coefs_t_avg<-round(colMeans( ml_coefs_t,na.rm=T),3)
bayes_coefs_t_avg<-round(colMeans( bayes_coefs_t,na.rm=T),3)
bayesIN_coefs_t_avg<-round(colMeans( bayesIN_coefs_t,na.rm=T),3)
fs_coefs_t_avg<-round(colMeans( fs_coefs_t,na.rm=T),3)
croon_coefs_t_avg<-round(colMeans( croon_coefs_t,na.rm=T),3)

#L1 average value by parameter and estimator
true_coefsL1_t_avg<-round(colMeans( true_coefsL1_t,na.rm=T),3)
ml_coefsL1_t_avg<-round(colMeans( ml_coefsL1_t,na.rm=T),3)
bayes_coefsL1_t_avg<-round(colMeans( bayes_coefsL1_t,na.rm=T),3)
bayesIN_coefsL1_t_avg<-round(colMeans( bayesIN_coefsL1_t,na.rm=T),3)
fs_coefsL1_t_avg<-round(colMeans( fs_coefsL1_t,na.rm=T),3)
croon_coefsL1_t_avg<-round(colMeans( croon_coefsL1_t,na.rm=T),3)

### Control
#L1 average value by parameter and estimator
true_coefsL1_c_avg<-round(colMeans( true_coefsL1_c,na.rm=T),3)
ml_coefsL1_c_avg<-round(colMeans( ml_coefsL1_c,na.rm=T),3)
bayes_coefsL1_c_avg<-round(colMeans( bayes_coefsL1_c,na.rm=T),3)
bayesIN_coefsL1_c_avg<-round(colMeans( bayesIN_coefsL1_c,na.rm=T),3)
fs_coefsL1_c_avg<-round(colMeans( fs_coefsL1_c,na.rm=T),3)
croon_coefsL1_c_avg<-round(colMeans( croon_coefsL1_c,na.rm=T),3)

###Treat
##average absolute bias across L2 structural coefs
ml_coefs_t_allavgbias<-round(mean(abs(colMeans(ml_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T))),3)
bayes_coefs_t_allavgbias<-round(mean(abs(colMeans(bayes_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T))),3)
bayesIN_coefs_t_allavgbias<-round(mean(abs(colMeans(bayesIN_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T))),3)
fs_coefs_t_allavgbias<-round(mean(abs(colMeans(fs_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T))),3)
croon_coefs_t_allavgbias<-round(mean(abs(colMeans(croon_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T))),3)
##average absolute bias across L1 structural coefs
ml_coefsL1_t_allavgbias<-round(mean(abs(colMeans(ml_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T))),3)
bayes_coefsL1_t_allavgbias<-round(mean(abs(colMeans(bayes_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T))),3)
bayesIN_coefsL1_t_allavgbias<-round(mean(abs(colMeans(bayesIN_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T))),3)
fs_coefsL1_t_allavgbias<-round(mean(abs(colMeans(fs_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T))),3)
croon_coefsL1_t_allavgbias<-round(mean(abs(colMeans(croon_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T))),3)
###Control
ml_coefsL1_c_allavgbias<-round(mean(abs(colMeans(ml_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T))),3)
bayes_coefsL1_c_allavgbias<-round(mean(abs(colMeans(bayes_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T))),3)
bayesIN_coefsL1_c_allavgbias<-round(mean(abs(colMeans(bayesIN_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T))),3)
fs_coefsL1_c_allavgbias<-round(mean(abs(colMeans(fs_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T))),3)
croon_coefsL1_c_allavgbias<-round(mean(abs(colMeans(croon_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T))),3)


#############################################################
### ## Precision/Efficiency                               
#############################################################
###Treat
#sd by param at L2
true_coefs_t_sd<-round(apply(true_coefs_t,2,sd,na.rm=T),3)
ml_coefs_t_sd<-round(apply(ml_coefs_t,2,sd,na.rm=T),3)
bayes_coefs_t_sd<-round(apply(bayes_coefs_t,2,sd,na.rm=T),3)
bayesIN_coefs_t_sd<-round(apply(bayesIN_coefs_t,2,sd,na.rm=T),3)
fs_coefs_t_sd<-round(apply(fs_coefs_t,2,sd,na.rm=T),3)
croon_coefs_t_sd<-round(apply(croon_coefs_t,2,sd,na.rm=T),3)

#sd by param at L1
true_coefsL1_t_sd<-round(apply(true_coefsL1_t,2,sd,na.rm=T),3)
ml_coefsL1_t_sd<-round(apply(ml_coefsL1_t,2,sd,na.rm=T),3)
bayes_coefsL1_t_sd<-round(apply(bayes_coefsL1_t,2,sd,na.rm=T),3)
bayesIN_coefsL1_t_sd<-round(apply(bayesIN_coefsL1_t,2,sd,na.rm=T),3)
fs_coefsL1_t_sd<-round(apply(fs_coefsL1_t,2,sd,na.rm=T),3)
croon_coefsL1_t_sd<-round(apply(croon_coefsL1_t,2,sd,na.rm=T),3)

###Control
#sd by param at L1
true_coefsL1_c_sd<-round(apply(true_coefsL1_c,2,sd,na.rm=T),3)
ml_coefsL1_c_sd<-round(apply(ml_coefsL1_c,2,sd,na.rm=T),3)
bayes_coefsL1_c_sd<-round(apply(bayes_coefsL1_c,2,sd,na.rm=T),3)
bayesIN_coefsL1_c_sd<-round(apply(bayesIN_coefsL1_c,2,sd,na.rm=T),3)
fs_coefsL1_c_sd<-round(apply(fs_coefsL1_c,2,sd,na.rm=T),3)
croon_coefsL1_c_sd<-round(apply(croon_coefsL1_c,2,sd,na.rm=T),3)


###Treat
#average sd across L2 coefs
true_coefs_t_allsd<-round(mean(apply(true_coefs_t,2,sd,na.rm=T)),3)
ml_coefs_t_allsd<-round(mean(apply(ml_coefs_t,2,sd,na.rm=T)),3)
bayes_coefs_t_allsd<-round(mean(apply(bayes_coefs_t,2,sd,na.rm=T)),3)
bayesIN_coefs_t_allsd<-round(mean(apply(bayesIN_coefs_t,2,sd,na.rm=T)),3)
fs_coefs_t_allsd<-round(mean(apply(fs_coefs_t,2,sd,na.rm=T)),3)
croon_coefs_t_allsd<-round(mean(apply(croon_coefs_t,2,sd,na.rm=T)),3)
#average sd across L1 coefs
true_coefsL1_t_allsd<-round(mean(apply(true_coefsL1_t,2,sd,na.rm=T)),3)
ml_coefsL1_t_allsd<-round(mean(apply(ml_coefsL1_t,2,sd,na.rm=T)),3)
bayes_coefsL1_t_allsd<-round(mean(apply(bayes_coefsL1_t,2,sd,na.rm=T)),3)
bayesIN_coefsL1_t_allsd<-round(mean(apply(bayesIN_coefsL1_t,2,sd,na.rm=T)),3)
fs_coefsL1_t_allsd<-round(mean(apply(fs_coefsL1_t,2,sd,na.rm=T)),3)
croon_coefsL1_t_allsd<-round(mean(apply(croon_coefsL1_t,2,sd,na.rm=T)),3)
###Control
#average sd across L1 coefs
true_coefsL1_c_allsd<-round(mean(apply(true_coefsL1_c,2,sd,na.rm=T)),3)
ml_coefsL1_c_allsd<-round(mean(apply(ml_coefsL1_c,2,sd,na.rm=T)),3)
bayes_coefsL1_c_allsd<-round(mean(apply(bayes_coefsL1_c,2,sd,na.rm=T)),3)
bayesIN_coefsL1_c_allsd<-round(mean(apply(bayesIN_coefsL1_c,2,sd,na.rm=T)),3)
fs_coefsL1_c_allsd<-round(mean(apply(fs_coefsL1_c,2,sd,na.rm=T)),3)
croon_coefsL1_c_allsd<-round(mean(apply(croon_coefsL1_c,2,sd,na.rm=T)),3)

### RMSE
###Treat
#average rmse across L2 coefs
ml_coefs_t_rmse<-round(mean(sqrt( (colMeans(ml_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T) )^2 +apply(ml_coefs_t,2,var,na.rm=T) )),3)
bayes_coefs_t_rmse<-round(mean(sqrt( (colMeans(bayes_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T) )^2 +apply(bayes_coefs_t,2,var,na.rm=T) )),3)
bayesIN_coefs_t_rmse<-round(mean(sqrt( (colMeans(bayesIN_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T) )^2 +apply(bayesIN_coefs_t,2,var,na.rm=T) )),3)
fs_coefs_t_rmse<-round(mean(sqrt( (colMeans(fs_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T) )^2 +apply(fs_coefs_t,2,var,na.rm=T) )),3)
croon_coefs_t_rmse<-round(mean(sqrt( (colMeans(croon_coefs_t,na.rm=T)-colMeans(true_coefs_t,na.rm=T) )^2 +apply(croon_coefs_t,2,var,na.rm=T) )),3)

#average rmse across L1 coefs
ml_coefsL1_t_rmse<-round(mean(sqrt( (colMeans(ml_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T) )^2 +apply(ml_coefsL1_t,2,var,na.rm=T) )),3)
bayes_coefsL1_t_rmse<-round(mean(sqrt( (colMeans(bayes_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T) )^2 +apply(bayes_coefsL1_t,2,var,na.rm=T) )),3)
bayesIN_coefsL1_t_rmse<-round(mean(sqrt( (colMeans(bayesIN_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T) )^2 +apply(bayesIN_coefsL1_t,2,var,na.rm=T) )),3)
fs_coefsL1_t_rmse<-round(mean(sqrt( (colMeans(fs_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T) )^2 +apply(fs_coefsL1_t,2,var,na.rm=T) )),3)
croon_coefsL1_t_rmse<-round(mean(sqrt( (colMeans(croon_coefsL1_t,na.rm=T)-colMeans(true_coefsL1_t,na.rm=T) )^2 +apply(croon_coefsL1_t,2,var,na.rm=T) )),3)
###Control
#average rmse across L1 coefs
ml_coefsL1_c_rmse<-round(mean(sqrt( (colMeans(ml_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T) )^2 +apply(ml_coefsL1_c,2,var,na.rm=T) )),3)
bayes_coefsL1_c_rmse<-round(mean(sqrt( (colMeans(bayes_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T) )^2 +apply(bayes_coefsL1_c,2,var,na.rm=T) )),3)
bayesIN_coefsL1_c_rmse<-round(mean(sqrt( (colMeans(bayesIN_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T) )^2 +apply(bayesIN_coefsL1_c,2,var,na.rm=T) )),3)
fs_coefsL1_c_rmse<-round(mean(sqrt( (colMeans(fs_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T) )^2 +apply(fs_coefsL1_c,2,var,na.rm=T) )),3)
croon_coefsL1_c_rmse<-round(mean(sqrt( (colMeans(croon_coefsL1_c,na.rm=T)-colMeans(true_coefsL1_c,na.rm=T) )^2 +apply(croon_coefsL1_c,2,var,na.rm=T) )),3)


#################################################################################################################################
###------------------------------------------------- End Function Save Results File ------------------------------------------###
#################################################################################################################################
# Save Errors
errorfilename = paste0("Errors_croon_partial_latent_med1_",n2,n1,rho,".csv", sep="")
#write.csv(d,file =filename )
write.csv(error_list,errorfilename)
# Save psr
psrfilename = paste0("psr_croon_partial_latent_med1_",n2,n1,rho,".csv", sep="")
#write.csv(d,file =filename )
write.csv(psr_list,psrfilename)

file.results<-c(n2,n1,nc,weight1, weight2, weight3, ind_rhox, ind_rhom, ind_rho,
                xi00, beta0, beta0_c,
                a, pi0, a_c,
                covar_coef,
                b1, B,  b1_c ,
                rhox, rhom, rho,
  
#single
  true_conv_full ,ml_conv_full, bayes_conv_full, bayesIN_conv_full ,fs_conv_t ,fs_convL1_t,croon_conv_t, croon_convL1_t, fs_convL1_c ,croon_convL1_c,
  ml_time_full_avg, bayes_time_full_avg, bayesIN_time_full_avg,
                
#multiple  
  true_coefs_t_avg,ml_coefs_t_avg,bayes_coefs_t_avg,bayesIN_coefs_t_avg,fs_coefs_t_avg,croon_coefs_t_avg,
  true_coefsL1_t_avg,ml_coefsL1_t_avg,bayes_coefsL1_t_avg,bayesIN_coefsL1_t_avg,fs_coefsL1_t_avg,croon_coefsL1_t_avg,
  true_coefsL1_c_avg,ml_coefsL1_c_avg,bayes_coefsL1_c_avg,bayesIN_coefsL1_c_avg,fs_coefsL1_c_avg,croon_coefsL1_c_avg,
  
#single  
  ml_coefs_t_allavgbias,bayes_coefs_t_allavgbias,bayesIN_coefs_t_allavgbias,fs_coefs_t_allavgbias,croon_coefs_t_allavgbias,
  ml_coefsL1_t_allavgbias,bayes_coefsL1_t_allavgbias,bayesIN_coefsL1_t_allavgbias,fs_coefsL1_t_allavgbias,croon_coefsL1_t_allavgbias,
  ml_coefsL1_c_allavgbias,bayes_coefsL1_c_allavgbias,bayesIN_coefsL1_c_allavgbias,fs_coefsL1_c_allavgbias,croon_coefsL1_c_allavgbias,

#multiple
true_coefs_t_sd,ml_coefs_t_sd,bayes_coefs_t_sd,bayesIN_coefs_t_sd,fs_coefs_t_sd,croon_coefs_t_sd,
true_coefsL1_t_sd,ml_coefsL1_t_sd,bayes_coefsL1_t_sd,bayesIN_coefsL1_t_sd,fs_coefsL1_t_sd,croon_coefsL1_t_sd,
true_coefsL1_c_sd,ml_coefsL1_c_sd,bayes_coefsL1_c_sd,bayesIN_coefsL1_c_sd,fs_coefsL1_c_sd,croon_coefsL1_c_sd,
#single
true_coefs_t_allsd,ml_coefs_t_allsd,bayes_coefs_t_allsd,bayesIN_coefs_t_allsd,fs_coefs_t_allsd,croon_coefs_t_allsd,
true_coefsL1_t_allsd,ml_coefsL1_t_allsd,bayes_coefsL1_t_allsd,bayesIN_coefsL1_t_allsd,fs_coefsL1_t_allsd,croon_coefsL1_t_allsd,
true_coefsL1_c_allsd,ml_coefsL1_c_allsd,bayes_coefsL1_c_allsd,bayesIN_coefsL1_c_allsd,fs_coefsL1_c_allsd,croon_coefsL1_c_allsd,
                
#single
ml_coefs_t_rmse,bayes_coefs_t_rmse,bayesIN_coefs_t_rmse,fs_coefs_t_rmse,croon_coefs_t_rmse,
ml_coefsL1_t_rmse,bayes_coefsL1_t_rmse,bayesIN_coefsL1_t_rmse,fs_coefsL1_t_rmse,croon_coefsL1_t_rmse,
ml_coefsL1_c_rmse,bayes_coefsL1_c_rmse,bayesIN_coefsL1_c_rmse,fs_coefsL1_c_rmse,croon_coefsL1_c_rmse
)

##variable names of results
names<-c("n2","n1","nc","weight1", "weight2", "weight3", "ind_rhox", "ind_rhom", "ind_rho",
                "xi00", "beta0", "beta0_c",
                "a", "pi0", "a_c",
                "covar_coef",
                "b1", "B",  "b1_c" ,
                "rhox", "rhom", "rho",
  
  "true_conv_full" ,"ml_conv_full", "bayes_conv_full", "bayesIN_conv_full" ,
  "fs_conv_t" ,"fs_convL1_t","croon_conv_t", "croon_convL1_t", "fs_convL1_c" ,"croon_convL1_c",
  "ml_time_full_avg", "bayes_time_full_avg", "bayesIN_time_full_avg",
  
  "true_coefs_ymL2",  "true_coefs_yxL2","true_coefs_ywL2", "true_coefs_mxL2","true_coefs_mwL2","true_coefs_ML2","true_coefs_YL2","true_coefs_MedL2","true_coefs_MainL2",
  "ml_coefs_ymL2",  "ml_coefs_yxL2","ml_coefs_ywL2", "ml_coefs_mxL2","ml_coefs_mwL2","ml_coefs_ML2","ml_coefs_YL2","ml_coefs_MedL2","ml_coefs_MainL2",
  "bayes_coefs_ymL2",  "bayes_coefs_yxL2","bayes_coefs_ywL2", "bayes_coefs_mxL2","bayes_coefs_mwL2","bayes_coefs_ML2","bayes_coefs_YL2","bayes_coefs_MedL2","bayes_coefs_MainL2",
  "bayesIN_coefs_ymL2","bayesIN_coefs_yxL2","bayesIN_coefs_ywL2", "bayesIN_coefs_mxL2","bayesIN_coefs_mwL2","bayesIN_coefs_ML2","bayesIN_coefs_YL2","bayesIN_coefs_MedL2","bayesIN_coefs_MainL2",
  "fs_coefs_ymL2",     "fs_coefs_yxL2","fs_coefs_ywL2", "fs_coefs_mxL2","fs_coefs_mwL2","fs_coefs_ML2","fs_coefs_YL2","fs_coefs_MedL2","fs_coefs_MainL2",
  "croon_coefs_ymL2",  "croon_coefs_yxL2","croon_coefs_ywL2", "croon_coefs_mxL2","croon_coefs_mwL2","croon_coefs_ML2","croon_coefs_YL2","croon_coefs_MedL2","croon_coefs_MainL2",
  
  "true_coefs_ymL1",  "true_coefs_yxL1", "true_coefs_mxL1",
  "ml_coefs_ymL1",  "ml_coefs_yxL1", "ml_coefs_mxL1",
  "bayes_coefs_ymL1",  "bayes_coefs_yxL1", "bayes_coefs_mxL1",
  "bayesIN_coefs_ymL1","bayesIN_coefs_yxL1", "bayesIN_coefs_mxL1",
  "fs_coefs_ymL1",     "fs_coefs_yxL1", "fs_coefs_mxL1",
  "croon_coefs_ymL1",  "croon_coefs_yxL1", "croon_coefs_mxL1",
  
  "true_coefs_ymL1_c",  "true_coefs_yxL1_c", "true_coefs_mxL1_c",
  "ml_coefs_ymL1_c",  "ml_coefs_yxL1_c", "ml_coefs_mxL1_c",
  "bayes_coefs_ymL1_c",  "bayes_coefs_yxL1_c", "bayes_coefs_mxL1_c",
  "bayesIN_coefs_ymL1_c","bayesIN_coefs_yxL1_c", "bayesIN_coefs_mxL1_c",
  "fs_coefs_ymL1_c",     "fs_coefs_yxL1_c", "fs_coefs_mxL1_c",
  "croon_coefs_ymL1_c",  "croon_coefs_yxL1_c", "croon_coefs_mxL1_c",
  
  "ml_coefs_t_allavgbias","bayes_coefs_t_allavgbias","bayesIN_coefs_t_allavgbias","fs_coefs_t_allavgbias","croon_coefs_t_allavgbias",
  "ml_coefsL1_t_allavgbias","bayes_coefsL1_t_allavgbias","bayesIN_coefsL1_t_allavgbias","fs_coefsL1_t_allavgbias","croon_coefsL1_t_allavgbias",
  "ml_coefsL1_c_allavgbias","bayes_coefsL1_c_allavgbias","bayesIN_coefsL1_c_allavgbias","fs_coefsL1_c_allavgbias","croon_coefsL1_c_allavgbias",

  "true_coefs_ymL2_sd",  "true_coefs_yxL2_sd","true_coefs_ywL2_sd", "true_coefs_mxL2_sd","true_coefs_mwL2_sd","true_coefs_ML2_sd","true_coefs_YL2_sd","true_coefs_MedL2_sd","true_coefs_MainL2_sd",
  "ml_coefs_ymL2_sd",  "ml_coefs_yxL2_sd","ml_coefs_ywL2_sd", "ml_coefs_mxL2_sd","ml_coefs_mwL2_sd","ml_coefs_ML2_sd","ml_coefs_YL2_sd","ml_coefs_MedL2_sd","ml_coefs_MainL2_sd",
  "bayes_coefs_ymL2_sd",  "bayes_coefs_yxL2_sd","bayes_coefs_ywL2_sd", "bayes_coefs_mxL2_sd","bayes_coefs_mwL2_sd","bayes_coefs_ML2_sd","bayes_coefs_YL2_sd","bayes_coefs_MedL2_sd","bayes_coefs_MainL2_sd",
  "bayesIN_coefs_ymL2_sd","bayesIN_coefs_yxL2_sd","bayesIN_coefs_ywL2_sd", "bayesIN_coefs_mxL2_sd","bayesIN_coefs_mwL2_sd","bayesIN_coefs_ML2_sd","bayesIN_coefs_YL2_sd","bayesIN_coefs_MedL2_sd","bayesIN_coefs_MainL2_sd",
  "fs_coefs_ymL2_sd",     "fs_coefs_yxL2_sd","fs_coefs_ywL2_sd", "fs_coefs_mxL2_sd","fs_coefs_mwL2_sd","fs_coefs_ML2_sd","fs_coefs_YL2_sd","fs_coefs_MedL2_sd","fs_coefs_MainL2_sd",
  "croon_coefs_ymL2_sd",  "croon_coefs_yxL2_sd","croon_coefs_ywL2_sd", "croon_coefs_mxL2_sd","croon_coefs_mwL2_sd","croon_coefs_ML2_sd","croon_coefs_YL2_sd","croon_coefs_MedL2_sd","croon_coefs_MainL2_sd",
  
  "true_coefs_ymL1_sd",  "true_coefs_yxL1_sd", "true_coefs_mxL1_sd",
  "ml_coefs_ymL1_sd",  "ml_coefs_yxL1_sd", "ml_coefs_mxL1_sd",
  "bayes_coefs_ymL1_sd",  "bayes_coefs_yxL1_sd", "bayes_coefs_mxL1_sd",
  "bayesIN_coefs_ymL1_sd","bayesIN_coefs_yxL1_sd", "bayesIN_coefs_mxL1_sd",
  "fs_coefs_ymL1_sd",     "fs_coefs_yxL1_sd", "fs_coefs_mxL1_sd",
  "croon_coefs_ymL1_sd",  "croon_coefs_yxL1_sd", "croon_coefs_mxL1_sd",
  
  "true_coefs_ymL1_c_sd",  "true_coefs_yxL1_c_sd", "true_coefs_mxL1_c_sd",
  "ml_coefs_ymL1_c_sd",  "ml_coefs_yxL1_c_sd", "ml_coefs_mxL1_c_sd",
  "bayes_coefs_ymL1_c_sd",  "bayes_coefs_yxL1_c_sd", "bayes_coefs_mxL1_c_sd",
  "bayesIN_coefs_ymL1_c_sd","bayesIN_coefs_yxL1_c_sd", "bayesIN_coefs_mxL1_c_sd",
  "fs_coefs_ymL1_c_sd",     "fs_coefs_yxL1_c_sd", "fs_coefs_mxL1_c_sd",
  "croon_coefs_ymL1_c_sd",  "croon_coefs_yxL1_c_sd", "croon_coefs_mxL1_c_sd",
  

"true_coefs_t_allsd","ml_coefs_t_allsd","bayes_coefs_t_allsd","bayesIN_coefs_t_allsd","fs_coefs_t_allsd","croon_coefs_t_allsd",
"true_coefsL1_t_allsd","ml_coefsL1_t_allsd","bayes_coefsL1_t_allsd","bayesIN_coefsL1_t_allsd","fs_coefsL1_t_allsd","croon_coefsL1_t_allsd",
"true_coefsL1_c_allsd","ml_coefsL1_c_allsd","bayes_coefsL1_c_allsd","bayesIN_coefsL1_c_allsd","fs_coefsL1_c_allsd","croon_coefsL1_c_allsd",
                
"ml_coefs_t_rmse","bayes_coefs_t_rmse","bayesIN_coefs_t_rmse","fs_coefs_t_rmse","croon_coefs_t_rmse",
"ml_coefsL1_t_rmse","bayes_coefsL1_t_rmse","bayesIN_coefsL1_t_rmse","fs_coefsL1_t_rmse","croon_coefsL1_t_rmse",
"ml_coefsL1_c_rmse","bayes_coefsL1_c_rmse","bayesIN_coefsL1_c_rmse","fs_coefsL1_c_rmse","croon_coefsL1_c_rmse"
)
### FULL Results
##name the results
file.results <-round(file.results,3)
names(file.results) <- names
file.results <- t(as.data.frame(file.results))


# Save individual conditions just in case
filename = paste0("croon_partial_med1",n2,n1,rho,".csv", sep="")
#write.csv(d,file =filename )
write.csv(file.results,filename)


return(file.results)

## end function 
}

#################################################################################################################################
###------------------------------------------------- Function Conditions -----------------------------------------------------###
#################################################################################################################################

file.results <- NULL #initial results
newdata<-NULL
save.results<-NULL

###Simulation Condition Function
# n1=15
newdata<-croon_partial_med(datasets=100, n2=15, n1=15, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)

newdata<-croon_partial_med(datasets=100, n2=30, n1=15, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)

newdata<-croon_partial_med(datasets=100, n2=60, n1=15, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)


newdata<-croon_partial_med(datasets=50, n2=90, n1=15, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)


# n1=30
newdata<-croon_partial_med(datasets=100, n2=15, n1=30, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)

newdata<-croon_partial_med(datasets=100, n2=30, n1=30, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)

newdata<-croon_partial_med(datasets=100, n2=60, n1=30, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)


newdata<-croon_partial_med(datasets=100, n2=90, n1=30, weight1 =1, weight2 = 1.5, weight3 = 0.666,ind_rhox=0.2,ind_rhom=0.2, ind_rho=0.2,
                                xi00<-0.8, beta0<-0.0, beta0_c<-0.0,
                                a<-0.5, pi0<-0.0, a_c<-0.0,
                                covar_coef<-0.3,
                                b1<-0.0, B<-0.5,  b1_c <-0.5,
                                rhox=0.2,rhom=0.2, rho=0.2)
save.results<- newdata; file.results<-rbind(file.results,save.results)


# Save by condition
write.csv(file.results,"croon_partial_med1.csv") 
#---------------------------------
# END