library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library("ggforce")
library("purrr")
library("vwmvp")
library("stringr")
library("tmvtnorm")
library("grid")
library("gridExtra")


options(pillar.sigfig = 5)
# models -----------------------------------------------------------------------

models <- c("MK_RNminus", "MK_FM_RNplus","MK_RNplus","MK_U_RNplus","MK_P_RNplus",
            "MK_P_RNminus","MK_U_RNminus","MK_FM_RNminus","J_RNminus","J_RNplus")

# Load data --------------------------------------------------------------------

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],
                            data_vdb14,data_pratte17)


FittedFull <- readRDS("Fits/FitFull.rds") %>% filter(model %in% models)

# Identify best fit

BestFull <- FittedFull %>%
  group_by(model,cvid) %>%
  filter(Dev == min(Dev, na.rm = TRUE)) %>% #finds best fit of 20 runs
  distinct() %>%
  group_by(cvid) %>%
  mutate(delLL = LL - min(LL),
         delDev = Dev - min(Dev),
         delAIC = AIC - min(AIC))


# Parameters for simulation ----------------------------------------------------
# Median

Parameters <- BestFull %>% select(c("model","mkappa1", "J1bar", "tau","alpha", "kappa_r","K")) %>%
  mutate(phi = ifelse(is.na(mkappa1),J1bar,mkappa1)) %>%
  mutate(kappa_r = ifelse(kappa_r > 4000,NA,kappa_r)) %>%
  gather(key = "Parameter",value = "value","phi", "tau","alpha", "kappa_r","K") %>%
  mutate(Parameter <- factor(Parameter,
                               levels=c("phi","alpha","tau","kappa_r","K"),
                               labels=c("phi","alpha","tau","kappa[r]","K")))


p <- c(.5)

p_names <- map_chr(p, ~paste0("qu",.x*100))

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>%
  set_names(nm = p_names)

GeneratingParameters <- Parameters %>%  group_by(model,Parameter) %>%
  summarize_at(vars(value), funs(!!!p_funs))

# Turn it into a grid.
# Long work-around away to keep it flexible for multiples
# If I knew more tidyverse, this would be less convoluted.

grid <- list()

for (i in seq_along(unique(GeneratingParameters$model))){

  modeln <- unique(GeneratingParameters$model)[[i]]

  grid[[i]] <- as_tibble(expand.grid(as.numeric(GeneratingParameters %>%
                                                  filter(model == modeln) %>%
                                                  ungroup() %>%
                                                  select(3) %>%
                                                  slice(1)),
                                     as.numeric(GeneratingParameters %>%
                                                  filter(model == modeln) %>%
                                                  ungroup() %>%
                                                  select(3) %>%
                                                  slice(2)),
                                     as.numeric(GeneratingParameters %>%
                                                  filter(model == modeln) %>%
                                                  ungroup() %>%
                                                  select(3) %>%
                                                  slice(3)),
                                     as.numeric(GeneratingParameters %>%
                                                  filter(model == modeln) %>%
                                                  ungroup() %>%
                                                  select(3) %>%
                                                  slice(4)),
                                     as.numeric(GeneratingParameters %>%
                                                  filter(model == modeln) %>%
                                                  ungroup() %>%
                                                  select(3) %>%
                                                  slice(5)))) %>%
    distinct() %>%
    set_names("kappa","alpha","tau","kappa_r","K") %>%
    tibble::rownames_to_column() %>%
    mutate(model = modeln,
           parid = paste0(modeln,"_p",rowname)) %>%
    select(-rowname)
}

grid <- bind_rows(grid) #<- generating parameters

# MVN-generated parameters -----------------------------------------------------

maxsetsizes <- experimentfile %>%
  group_by(id) %>%
  summarize(maxSs = max(set_size))


grid <- NULL
for (modeln in models){

# identify relevant parameters
selecty <- BestFull %>%
  select(c("model","mkappa1", "J1bar", "tau","alpha", "kappa_r","K")) %>%
  ungroup %>%
  gather(key = "Parameter",
         value = "value","mkappa1","J1bar", "tau","alpha", "kappa_r","K") %>%
  group_by(Parameter) %>%
  filter(model == modeln) %>%
  filter(Parameter %in% names(vwmvp::get_start_vp(modeln)))


datas <- selecty %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(Parameter,value) %>%
  drop_na() %>%
  select(-model,-cvid) %>%
  select(names(vwmvp::get_start_vp(modeln))) %>%
  rename(phi = 1)

# for each model: identify mean/sd of paramaters and correlations between them

correlations <- cor(datas)
means <- as.numeric(t(datas %>% mutate_all(mean) %>% .[1,])[,1])
sds <- as.numeric(t(datas %>% mutate_all(sd) %>% .[1,])[,1])
sigmas <- sds %*% t(sds) * correlations

# Use multivariate to generate parameters on basis of multivariate distribution
# Sample 1000 possible sets of generating parameters per model

parameters <- tmvtnorm::rtmvnorm(n=1000, mean = means, sigma=sigmas,
                                 lower=rep(0,length(means)),
                                 upper=rep(Inf,length(means)),
                                 algorithm="gibbs",
                                 burn.in.samples=1000,
                                 thinning = 100)


# Identify proportion of datasets where the best fit of the model estimates a
# K that exceeds the maximum set size in that experiment

if (modeln %in% c("MK_FM_RNminus","MK_FM_RNplus",
                  "MK_P_RNminus","MK_P_RNplus",
                  "MK_U_RNminus","MK_U_RNplus")){
# adjust for K being the mean (rather than the end point of the uniform)
critdiv <- ifelse(modeln %in% c("MK_U_RNminus","MK_U_RNplus"),2,1)

proportionlimited <- selecty %>%
  filter(value < quantile(value,.975) & value > quantile(value,.025)) %>%
  spread(Parameter,value) %>%
  bind_cols(maxSs = maxsetsizes) %>%
  mutate(Kbelow = ifelse(K/critdiv < maxSs,1,2)) %>%
  drop_na() %>%
  group_by(Kbelow) %>%
  count() %>%
  ungroup() %>%
  mutate(newn = n/sum(n))

# proportiontot <- selecty %>%
#   spread(Parameter,value) %>%
#   bind_cols(maxSs = maxsetsizes) %>%
#   mutate(Kbelow = ifelse(K/critdiv < maxSs,1,2)) %>%
#   group_by(Kbelow) %>%
#   count() %>%
#   ungroup() %>%
#   mutate(newn = n/sum(n))

# Rejection sampling for sets of generated parameters where K does not exceed 8
# 8 = max set size in the simulated data sets

parameters <- as_tibble(parameters) %>%
  set_names(names(vwmvp::get_start_vp(modeln))) %>%
  mutate(Kbelow = ifelse(K/critdiv < 8,1,2))

# aim for a total of 40 sets of generating parameters per data set
# Identify the proportion of data sets where K is below max set size
# to generate the number of sets of parameters required

if(any(proportionlimited %>% .$Kbelow == 1)){

  propbelow <- proportionlimited %>% filter(Kbelow == 1) %>% .$newn %>% .[[1]]
  Kbelowtrials <- round(propbelow * 40)
  parbelow <- parameters %>% filter(Kbelow == 1) %>% slice(1:Kbelowtrials)

} else {

# if for the model K is not below max set size for any data sets
  propbelow <- 0
  Kbelowtrials <- round(propbelow * 40)
  parbelow <- NULL
}

# Get the sets of parameters where K > max set size

parabove <-  parameters %>% filter(Kbelow == 2) %>% slice(1:(40 - Kbelowtrials))
par <- bind_rows(parbelow,parabove) %>%
  mutate(model = modeln,
         rn = 1:n(),
         parid = paste0(modeln,"_p",rn+1)) %>%
  select(-rn)

} else {


  parameters <- as_tibble(parameters) %>% set_names(names(vwmvp::get_start_vp(modeln)))
  par <- parameters %>% slice(1:40) %>%
    mutate(model = modeln,
           rn = 1:n(),
           parid = paste0(modeln,"_p",rn+1)) %>%
    select(-rn)

}
grid <- bind_rows(grid,par)
}

grid <- grid %>%
  mutate(kappa = ifelse(model %in% c("J_RNminus","J_RNplus"),J1bar,mkappa1)) %>%
  select(-mkappa1,-J1bar, -Kbelow) # <- generating parameters for all models


# Generate for each parameter combination 20 data sets -----------------------

# generating parameters are identical for different trials per set size
# number of set sizes, trials per set size need to be adapted below to match
# generating data for median parameter estimates


for (trialsSs in c(60,120,220,320,420,520,620,720,820,920,1020,2020)){

generateddata <- NULL
for (parameterid in unique(grid %>% .$parid)){

    par <- as.numeric(grid %>%
                        filter(parid == parameterid) %>%
                        select(kappa,alpha,tau,kappa_r,K))
    gendata <- vwmvp::generate_data(trials = rep(trialsSs,6),
                                    set_sizes = c(1,2,3,4,6,8),
                           model = grid %>%
                             filter(parid == parameterid) %>%
                             .$model,
                           par = par[!is.na(par)])
    gendata <- gendata %>%
      mutate(model = grid %>% filter(parid == parameterid) %>% .$model,
           id = parameterid,
           trialsSs = trialsSs) %>%
      group_nest(model,id,keep=T) %>%
      mutate(Parameters = list(grid %>% filter(parid == parameterid)))

  generateddata <- bind_rows(generateddata,gendata)

}

#saveRDS(generateddata,file = paste0("RecoveryData/ParameterRecovery_data_",trialsSs,"_rmv.rds"))


}
