#install.packages("vwmvp_1.0.tar.gz", repos = NULL, type="source",dependencies = TRUE)
library("vwmvp")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("RColorBrewer")
library("ggplot2")
library("gridExtra")
library("grid")
library("cowplot")


# Error distributions ----------------------------------------------------------


# The core prediction function is
# vwmvp::predict_data(pars = best-fit parameters of id and model,
#                     model = modelstring,
#                     data = data.frame(id = id,
#                                       set_size = unique set sizes),
#                                       error_0 = possible error range)

# In this file this function is used in various loops that depending on the type of
# data being predicted, loops across: set sizes, hold-out sets, individuals, experiments
# and produces inidvidual-level predictions, averaged parameter-estimates predictions,
# predictions based on LOSSO-CV and experiment-level predictions.

# These predictions are saved in the Prediction/ folder.

# errorspace in radian
rads <- c(-180:180) * pi / 180
rads <- ifelse(rads < -pi, rads + 2 * pi, ifelse(rads >
                                                   pi, rads - 2 * pi, rads))

# For individual-level

models <- c("J_RNminus","J_RNplus")



for (experiment in c("ol17_e1","WM4", "RTT12", "ZL8", "BCH9", "AVA11",
                     "AA12a", "AA12b", "VSCGM21a", "VSCGM21b", "VSCGM21c","pratte17")){

  LOTest <- readRDS("Fits/CVLOSzO_TestLL.rds") %>%
    filter(exp %in% experiment) %>%
    filter(model %in% models)

  FailedSzPrediction <-c() #unique(LOTest %>% filter(CVLL > 1e4) %>% .$id)

  FullFit <- readRDS("Fits/FitFull.rds") %>%
    filter(exp %in% experiment) %>%
    filter(model %in% models)

  LOTraining <- readRDS("Fits/CVLOSzO_Trainingfits.rds")  %>%
    filter(exp %in% experiment) %>%
    filter(model %in% models)

  LOTest <- LOTest %>% filter(!id %in% FailedSzPrediction)


  load("Data/ol17_prepared.rda")
  load("Data/vdb14_prepared.rda")
  load("Data/pratte17_prepared.rda")

  experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
    filter(exp == experiment) %>%
    rename(setsize = "set_size")


  expsetsize <- sort(unique(experimentfile %>% .$setsize))

  participants <- unique(FullFit$id)


  # Average


  AvPrediction <- NULL
  for (modeln in models){

    BestFull <- FullFit %>%
      filter(model == modeln) %>%
      group_by(id) %>%
      arrange(objective) %>% slice(1) %>%
      select(names(vwmvp::get_start_vp(modeln))) %>%
      ungroup() %>%
      summarise_all(mean) %>% select(-id)



    BestLOSso <- LOTraining %>%
      filter(model == modeln) %>%
      arrange(id,leftout,objective) %>%
      group_by(id,leftout) %>% slice(1) %>%
      select(names(vwmvp::get_start_vp(modeln))) %>%
      group_by(leftout) %>%
      summarise_all(mean) %>% select(-id)

    predict <- vwmvp::predict_data(pars = as.numeric(BestFull %>%
                                                       select(names(vwmvp::get_start_vp(modeln))) %>%
                                                       ungroup() ),
                                   model = modeln,
                                   data = data.frame(id = "Av",
                                                     set_size = rep(expsetsize,each=length(rads)),
                                                     error_0 = rep(rads,length(expsetsize)))
    )


    predict$prediction <- as.numeric(predict$prediction)
    predict <- predict %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

    predict$leftout <- 0


    LOSSO <- NULL
    for (leftouts in unique(BestLOSso$leftout)){

      predict2 <- vwmvp::predict_data(pars = as.numeric(BestLOSso %>% ungroup() %>%
                                                          filter(leftout == leftouts) %>%
                                                          select(-leftout)),
                                      model = modeln,
                                      data = data.frame(id = "Av",
                                                        set_size = rep(expsetsize,each = length(rads)),
                                                        error_0 = rep(rads,length(expsetsize)))
      )

      predict2$prediction <- as.numeric(predict2$prediction)
      predict2 <- predict2 %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

      predict2$leftout <- leftouts
      LOSSO <- LOSSO %>% bind_rows(predict2)
    }



    AvPrediction <- AvPrediction %>% bind_rows(predict,LOSSO)


  }

  #saveRDS(AvPrediction, file=paste0("Prediction/prediction_idav_",experiment,".rds"))



  IndividualPreds <- NULL
  for (Participant in participants) {

    print(Participant)
    IndPrediction <- NULL
    for (modeln in models){

      BestFull <- FullFit %>%
        filter(id == Participant) %>%
        filter(model == modeln) %>%
        group_by(setsize) %>%
        arrange(objective) %>% slice(1) %>%
        select(names(vwmvp::get_start_vp(modeln))) %>% ungroup()

      BestLOSso <- LOTraining %>%
        filter(id == Participant) %>%
        filter(model == modeln) %>%
        arrange(id,leftout,objective) %>%
        group_by(id,leftout) %>% slice(1) %>%
        select(names(vwmvp::get_start_vp(modeln))) %>% ungroup()

      predict <- vwmvp::predict_data(pars = as.numeric(BestFull),
                                     model = modeln,
                                     data = data.frame(id = Participant,
                                                       set_size = rep(expsetsize,each = length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))
      )



      predict <- vwmvp::predict_data(pars = as.numeric(BestFull %>% filter(setsize == indsetsize) %>%
                                                         select(names(vwmvp::get_start_vp(modeln))) %>%
                                                         ungroup() ),
                                     model = modeln,
                                     data = data.frame(id = Participant,
                                                       set_size = rep(expsetsize,each = length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))
      )


      predict$prediction <- as.numeric(predict$prediction)
      predict <- predict %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

      predict$leftout <- 0


      LOSSO <- NULL
      for (leftouts in unique(BestLOSso$leftout)){

        predict2 <- vwmvp::predict_data(pars = as.numeric(BestLOSso %>% ungroup() %>%
                                                            filter(leftout == leftouts) %>%
                                                            select(-id,-leftout)),
                                        model = modeln,
                                        data = data.frame(id = Participant,
                                                          set_size = rep(expsetsize,each = length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))
        )

        predict2$prediction <- as.numeric(predict2$prediction)
        predict2 <- predict2 %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

        predict2$leftout <- leftouts
        LOSSO <- LOSSO %>% bind_rows(predict2)
      }



      IndPrediction <- IndPrediction %>% bind_rows(predict,LOSSO)

    }

    IndividualPreds <- IndividualPreds %>% bind_rows(IndPrediction)

  }

  Prediction <- bind_rows(IndividualPreds,AvPrediction) %>% mutate(exp = experiment)
  #saveRDS(Prediction, file=paste0("Prediction/prediction_",experiment,".rds"))
}

# For experiment-level

models <- c("MK_RNminus","MK_RNplus")


for (experiment in c("ol17_e1","WM4", "RTT12", "ZL8", "BCH9", "AVA11",
                     "AA12a", "AA12b", "VSCGM21a", "VSCGM21b", "VSCGM21c","pratte17")){

  LOTest <- readRDS("Fits/FitFull_AggLOSSO_Test.rds") %>%
    filter(exp %in% experiment) %>%
    filter(model %in% models)

  FailedSzPrediction <-c() #unique(LOTest %>% filter(CVLL > 1e4) %>% .$id)

  FullFit <- readRDS("Fits/FitFull_agg.rds") %>%
    filter(exp %in% experiment) %>%
    filter(model %in% models)

  LOTraining <- readRDS("Fits/FitFull_AggLOSSO_TestLL.rds")  %>%
    filter(exp %in% experiment) %>%
    filter(model %in% models)

  LOTest <- LOTest %>% filter(!id %in% FailedSzPrediction)


  load("Data/ol17_prepared.rda")
  load("Data/vdb14_prepared.rda")
  load("Data/pratte17_prepared.rda")

  experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)  %>%
    filter(exp == experiment) %>%
    rename(setsize = "set_size")


  expsetsize <- sort(unique(experimentfile %>% .$setsize))

  participants <- unique(FullFit$id) # in the Aggregate fits exp = id



  IndividualPreds <- NULL
  for (Participant in participants) {

    print(Participant)
    IndPrediction <- NULL
    for (modeln in models){

      BestFull <- FullFit %>%
        filter(id == Participant) %>%
        filter(model == modeln) %>%
        group_by(setsize) %>%
        arrange(objective) %>% slice(1) %>%
        select(names(vwmvp::get_start_vp(modeln))) %>% ungroup()

      BestLOSso <- LOTraining %>%
        filter(id == Participant) %>%
        filter(model == modeln) %>%
        arrange(id,leftout,objective) %>%
        group_by(id,leftout) %>% slice(1) %>%
        select(names(vwmvp::get_start_vp(modeln))) %>% ungroup()

      predict <- vwmvp::predict_data(pars = as.numeric(BestFull),
                                     model = modeln,
                                     data = data.frame(id = Participant,
                                                       set_size = rep(expsetsize,each = length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))
      )



      predict <- vwmvp::predict_data(pars = as.numeric(BestFull %>% filter(setsize == indsetsize) %>%
                                                         select(names(vwmvp::get_start_vp(modeln))) %>%
                                                         ungroup() ),
                                     model = modeln,
                                     data = data.frame(id = Participant,
                                                       set_size = rep(expsetsize,each = length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))
      )


      predict$prediction <- as.numeric(predict$prediction)
      predict <- predict %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

      predict$leftout <- 0


      LOSSO <- NULL
      for (leftouts in unique(BestLOSso$leftout)){

        predict2 <- vwmvp::predict_data(pars = as.numeric(BestLOSso %>% ungroup() %>%
                                                            filter(leftout == leftouts) %>%
                                                            select(-id,-leftout)),
                                        model = modeln,
                                        data = data.frame(id = Participant,
                                                          set_size = rep(expsetsize,each = length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))
        )

        predict2$prediction <- as.numeric(predict2$prediction)
        predict2 <- predict2 %>% group_by(setsize) %>% mutate(normprediction = prediction/sum(prediction))

        predict2$leftout <- leftouts
        LOSSO <- LOSSO %>% bind_rows(predict2)
      }



      IndPrediction <- IndPrediction %>% bind_rows(predict,LOSSO)

    }

    IndividualPreds <- IndividualPreds %>% bind_rows(IndPrediction)

  }

  Prediction <- bind_rows(IndividualPreds) %>% mutate(exp = experiment)
  #saveRDS(Prediction, file=paste0("Prediction/prediction_agg_",experiment,".rds"))
}

# Summary statistics -----------------------------------------------------

# For the graphs in the manuscript, this routine is replicated in the individual files
# for the compound figures

# This routine here produces files that are used in the the experiment-/individual-level
# .html files. Note: for RQ3 we present averaged predictions. We excluded participants that showed
# extreme LOSSO out-of-sample deviance for at least one of the models in the comparison
# accordingly, we also removed that individual's data before aggregating. If the below routine is
# run across all models, the data removed will differ.

# MakeSumStat function is wrapped in various loops to pull the correct set of
# predictions for the specific case.
# Results are put into SummaryStat/ folder

MakeSumStat <- function(Prediction, modelname){


  CV <- NULL
  for(SS in unique(Prediction %>% .$setsize)){


    CVid <- NULL
    for (subj in unique(Prediction$id)) {

      CVleft <- NULL

      for(leftouts in unique(Prediction$leftout)) {


        print(modelname)
        print(subj)
        print(leftouts)

        testNA <- any(is.na(Prediction %>%
                              filter(id == subj) %>%
                              filter(setsize == SS) %>%
                              filter(leftout == leftouts) %>%
                              filter(model == modelname) %>%
                              .$normprediction) == TRUE)
        if (testNA){

          CircVar <- NA
          MAE <- NA
          kurtosis <- NA

        } else {

          samples <- sample(size=6e4,Prediction %>%
                              filter(id == subj) %>%
                              filter(leftout == leftouts) %>%
                              filter(setsize == SS) %>%
                              filter(model == modelname) %>%
                              .$data,replace=TRUE,
                            prob =  Prediction %>%
                              filter(id == subj) %>%
                              filter(leftout == leftouts) %>%
                              filter(setsize == SS) %>%
                              filter(model == modelname) %>%
                              .$normprediction)
          CircVar <- circular::var.circular(samples)

          MAE <- mean(abs(samples))
          kurtosis <- vwmvp::circkurtosis(samples)

        }

        out <- as_tibble(data.frame(id = subj,
                                    leftout = leftouts,
                                    setsize = SS,
                                    model = modelname,
                                    Var = CircVar,

                                    meanErr = MAE,
                                    kurt = kurtosis))
        CVleft <- CVleft %>%  bind_rows(out)

      }

      CVid <- CVid %>%  bind_rows(CVleft)
    }

    CV <- CV %>%  bind_rows(CVid)

  }
  return(CV)
}



# Load observed data


load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")

experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],
                            data_vdb14,data_pratte17)  %>%
  rename(setsize = "set_size")

expsetsize <- sort(unique(experimentfile %>% .$setsize))

models <- c("MK_RNminus","MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus",
            "MK_RNplus","MK_FM_RNplus","MK_P_RNplus","MK_U_RNplus")

for (experiment in c("ol17_e1","WM4", "RTT12", "ZL8", "BCH9", "AVA11",
                     "AA12a", "AA12b", "VSCGM21a", "VSCGM21b", "VSCGM21c","pratte17")){


  LOTest <- readRDS("Fits/CVLOSzO_TestLL.rds") %>%
    filter(exp == experiment) %>%

    filter(model %in% models)

  FailedSzPrediction <- unique(LOTest %>% filter(CVLL > 1e4) %>% .$id) # extreme LOSSO ids

  PredictionInd <- bind_rows(readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_ind.rds")),
                             readRDS(paste0("Prediction/prediction_nonVP_full_",experiment,"_ind.rds")))

  PredictionIDAV <- readRDS(paste0("Prediction/prediction_idavExclFailed_",experiment,".rds"))

  PredictionAgg <- bind_rows(readRDS(paste0("Prediction/prediction_VP_fullLOSSO_",experiment,"_agg.rds")),
                             readRDS(paste0("Prediction/prediction_nonVP_full_",experiment,"_agg.rds")))


  FRPredictionAvP <- PredictionInd %>%
    filter(model %in% models) %>%
    #mutate(model = factor(model,levels=models)) %>%
    filter(id != "Av") %>%
    filter(!id %in% FailedSzPrediction) %>%
    group_by(model,leftout,setsize,data) %>%
    summarize(mprediction = mean(prediction)) %>%
    rename("prediction" = mprediction) %>%
    mutate(normprediction = prediction/sum(prediction)) %>%
    mutate(id = "Averagepred")

  CV_aggdata <- NULL

  for(SS in unique(experimentfile %>%
                   filter(exp == experiment) %>%  .$setsize)){

    expdata <- experimentfile %>%
      filter(exp == experiment) %>%
      #  filter(!id %in% FailedSzPrediction) %>%
      filter(setsize == SS) %>%
      .$error_0

    cv <- circular::var.circular(expdata)
    mae <- mean(abs(expdata))
    kurtosis_data <- vwmvp::circkurtosis(expdata)

    out2 <- as_tibble(data.frame(id = "Agg",
                                 leftout = 9,
                                 setsize = SS,
                                 model = "Data",
                                 Var = cv,
                                 meanErr = mae,
                                 kurt = kurtosis_data))

    CV_aggdata <- CV_aggdata %>% bind_rows(out2)
  }

  CV_avdata <- NULL

  for(SS in unique(experimentfile %>%
                   filter(exp == experiment) %>%  .$setsize)){

    expdata <- experimentfile %>%
      filter(exp == experiment) %>%
      filter(!id %in% FailedSzPrediction) %>%
      filter(setsize == SS) %>%
      .$error_0

    cv <- circular::var.circular(expdata)
    mae <- mean(abs(expdata))
    kurtosis_data <- vwmvp::circkurtosis(expdata)

    out2 <- as_tibble(data.frame(id = "Averagepred",
                                 leftout = 9,
                                 setsize = SS,
                                 model = "Data",
                                 Var = cv,
                                 meanErr = mae,
                                 kurt = kurtosis_data))

    CV_avdata <- CV_avdata %>% bind_rows(out2)
  }

  CV_inddata <- NULL


  for(SS in unique(experimentfile %>%
                   filter(exp == experiment) %>%  .$setsize)){

    CVid <- NULL
    for (subj in unique(experimentfile %>%
                        filter(exp == experiment) %>%  .$id)){


      expdata <- experimentfile %>%
        filter(exp == experiment) %>%
        filter(id == subj) %>%
        filter(setsize == SS) %>%
        .$error_0

      cv <- circular::var.circular(expdata)
      mae <- mean(abs(expdata))
      kurtosis_data <- vwmvp::circkurtosis(expdata)

      out2 <- as_tibble(data.frame(id = subj,
                                   leftout=9,
                                   setsize = SS,
                                   model = "Data",
                                   Var = cv,
                                   meanErr = mae,
                                   kurt = kurtosis_data))

      CVid <- CVid %>% bind_rows(out2)
    }

    CV_inddata <- CV_inddata %>% bind_rows(CVid)
  }

  CV_ind <- NULL
  CV_agg <- NULL
  CV_averagepred <- NULL
  CV_avidexcl <- NULL

  for(modelname in models){

    CV_ind <- CV_ind %>% bind_rows(MakeSumStat(PredictionInd,modelname))
    CV_averagepred <- CV_averagepred %>% bind_rows(MakeSumStat(FRPredictionAvP,modelname))
    CV_agg <- CV_agg %>% bind_rows(MakeSumStat(PredictionAgg,modelname))
    CV_avidexcl <- CV_avidexcl %>% bind_rows(MakeSumStat(PredictionIDAV,modelname))
  }


  # CVind <- bind_rows(CV_ind,CV_inddata) %>% mutate(exp = experiment)
  # saveRDS(CVind,file=paste0("SummaryStat/SumStat_VP_",experiment,"_ind.rds"))
  # CVagg <- bind_rows(CV_agg,CV_aggdata) %>% mutate(exp = experiment)
  # saveRDS(CVagg,file=paste0("SummaryStat/SumStat_VP_",experiment,"_agg.rds"))
  # CVaveragepred <- bind_rows(CV_averagepred,CV_avdata) %>% mutate(exp = experiment)
  # saveRDS(CVaveragepred,file=paste0("SummaryStat/SumStat_VP_",experiment,"_averagepred.rds"))
  #
  # CVavidexcl <- bind_rows(CV_avidexcl,CV_avdata %>% mutate(id = "AvExcl"))
  # saveRDS(CVaver,file=paste0("SummaryStat/SumStatLOSSO_",experiment,"_idAvExcl.rds"))

}