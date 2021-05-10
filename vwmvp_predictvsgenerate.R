# Test equivalence of vwmvp::generate_data and vwmvp::predict_data

library(vwmvp)
library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("grid")
library("cowplot")

# Note: in vwmvp::generate_data for VP-U models, U is defined as ranging from 0 - 2K
# for all other functions in vwmvp, U is defined as ranging from 0 - K
# for ~equivalent results any input K has to be multiplied by 2 to give the correct 
# max of 2K (and mean K). Any resultant parameter estimate has to be divided by 2
# to provide the estimate of mean K, not the max of 2K.

modelname <- "MK_U_RNminus"
paramvaluesGENDAT <- c(60,1.6,15,4) 
paramvaluesPREDDAT <- c(60,1.6,15,4*2) 

samplesize <- 6e4
setsizes <- c(1:8)

rads <- c(-180:180) * pi / 180
rads <- ifelse(rads < -pi, rads + 2 * pi, ifelse(rads > pi, rads - 2 * pi, rads))


sample_bygenerate <- vwmvp::generate_data(model = modelname,
                                          par = paramvaluesGENDAT,
                                          trials = rep(samplesize,length(setsizes)),
                                          set_size = setsizes)

probprediction <- vwmvp::predict_data(pars = paramvaluesPREDDAT,
                                      model = modelname,
                                      data = data.frame(id = "Test",
                                                      set_size = rep(setsizes,each = length(rads)),
                                                      error_0 = rep(rads,length(setsizes))))
probprediction <- probprediction %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

sample_bypredict <- NULL
for(i in seq_along(setsizes)){
  
SSi <- sample(size = samplesize,
                           probprediction %>% 
                             filter(setsize == setsizes[[i]]) %>% .$data,
                           replace=TRUE, 
                           prob =  probprediction %>% 
                             filter(setsize == setsizes[[i]]) %>% .$prediction)

SSidf <- as_tibble(data.frame(set_size = setsizes[[i]],
                     error_0 = SSi))

sample_bypredict <- sample_bypredict %>% bind_rows(SSidf)

}

sampleby <- bind_rows(sample_bypredict %>% mutate(type = "bypredict"),
          sample_bygenerate %>% mutate(type = "bygenerate")) %>% 
  mutate(model = modelname)


ggplot(sampleby,aes(x = error_0,color=type)) +
  geom_density(size=1) +
  coord_cartesian(xlim = c(-pi/2,pi/2))+
  facet_wrap(.~set_size)

sampleby %>% group_by(type,set_size) %>% summarize(mae = mean(abs(error_0)),
                                          vars =  circular::var.circular(error_0),
                                          kurt = vwmvp::circkurtosis(error_0),
                                          numtrials = length(error_0))

# Percent change in best-fit parameters in LOSSO ---------------------------------------

# Note: parameter estimates are a bit unstable in general so grain of salt
# generalizes to VP models with limited capacity (K is much less predictable)

Fits <- readRDS("Fits/FitFull.rds")
FitsLOSSO <- readRDS("Fits/CVLOSzO_Trainingfits.rds")

setsizeeffects <- bind_rows(Fits,FitsLOSSO) %>% 
  group_by(id,model,leftout) %>% 
  arrange(objective) %>% 
  slice(1) %>% 
  filter(model %in% c("MK_RNplus","MK_RNminus"))

bps <- setsizeeffects %>% group_by(exp,model,id,leftout) %>% 
  select(exp,id,model,leftout,alpha,mkappa1,tau,kappa_r) %>% 
  pivot_longer(!c(exp,id,model,leftout),names_to = "measure",values_to="value") %>% 
  arrange(exp,id,model,measure,leftout) %>% group_by(exp,id,model,measure) %>% 
  mutate(percentvalue = value/value[[1]]*100 - 100) %>% 
  filter(leftout!=0) %>% 
  mutate(leftout = factor(leftout,levels=c(1:8)))

# massive outliers / calculating mean + 95CI probably misleading

ggplot(bps, aes(x=leftout,y=percentvalue,group=leftout)) +
  geom_hline(yintercept = 0) +
  geom_boxplot() +
  coord_cartesian(ylim= c(-100,200))+
  facet_grid(model~measure) +
  scale_y_continuous(name="% change from full data set fitted")
