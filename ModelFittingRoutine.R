library("doParallel")
library("vwmvp")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("foreach")

load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)

# test models ------------------------------------------------------------------


models <- c("SA_F_RNplus","SA_F_RNminus",
            "SA_P_RNplus", "SA_U_RNplus","EP_P_RNplus","EP_U_RNplus",
            "EP_RNplus","SA_RNplus",
            "EP_RNminus", "EP_FM_RNplus","EP_FM_RNminus",
            "EP_P_RNminus","EP_U_RNminus","SA_P_RNminus","SA_U_RNminus")

# prepare data: nested by id for foreach ---------------------------------------

# cvid, leftout to make it compatible with fitting cv data parallelized across training sets
# dat_prep <- experimentfile %>%
#   mutate(cvid = id,
#          leftout = 0) %>%
#   group_nest(exp,cvid,keep=T) %>%
#   filter(cvid == "ol17_e1_6")

dat_prep <- experimentfile %>%
  mutate(cvid = exp,
         id = exp,
         leftout = 0) %>%
  group_nest(exp,keep=T)


# Fit data ---------------------------------------------------------------------
doParallel::registerDoParallel(cores=24)
mcoptions <- list(preschedule=FALSE, set.seed=FALSE)
# tots <- NULL
# res <- foreach(Ss = unique(dat_prep %>% .$data %>% .[[1]] %>% .$set_size)
#                ,.combine = 'rbind'
#                ,.options.multicore=mcoptions) %do% {
#
#                  fitted <-vwmvp::FitVP(data = dat_prep %>%
#                                          .$data %>% .[[1]] %>% filter(set_size == Ss),
#                                        model = modelname,
#                                        method = "numint", # numerical integration
#                                        rep = 20,
#                                        startpar = NULL) # random starting parameters
#
#                  fitted <- fitted %>% mutate(setsize = Ss)
#                  tots <- tots %>% bind_rows(fitted)
#
#                }
# saveRDS(tots, file = paste0("FitVPnosetsize/", "ol17_e1_6_", modelname, "_Full20Random.rds"))


for (exps in unique(dat_prep$exp)){

  tots <- NULL
res <- foreach(modelname = models
               ,.combine = 'rbind'
               ,.options.multicore=mcoptions) %do% {

                 fitted <-vwmvp::FitVP(data = dat_prep %>% filter(exp == exps) %>% .$data %>% .[[1]],
                                       model = modelname,
                                       method = "numint", # numerical integration
                                       rep = 20,
                                       startpar = NULL) # random starting parameters

                 fitted <- fitted %>% mutate(model = modelname)
                 tots <- tots %>% bind_rows(fitted)
                 print(exps)
                 print(modelname)

               }
#saveRDS(tots, file = paste0("FullFitAgg/Full/", exps, "_nonVPmodels.rds"))
}
