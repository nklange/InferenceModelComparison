library("tidyverse")
library("xtable")
library("circular")
library("vwmvp")

models <- c("MK_RNminus","MK_RNplus","J_RNminus","J_RNplus","MK_F_RNplus",
            "MK_P_RNplus","MK_U_RNplus","EP_F_RNminus","SA_F_RNminus")
models2 <- c("MK_RNminus","MK_RNplus",
                        "MK_FM_RNplus","MK_U_RNplus","MK_FM_RNminus",
                        "MK_U_RNminus")

# Prediction wrapper ----------------------------------------

PredictHoldOutSz <- function(TrainData, modelN, HoldOutData){
  
  pars <- TrainData %>% filter(model == modelN) %>% 
    select(names(vwmvp::get_start_vp(modelN))) %>% unlist()
  
  dp <- vwmvp::prep_data(HoldOutData)
  set_sizes <- dp$set_sizes
  error_list <- dp$datalist
  
  out <- vwmvp::ll_vp_numint(pars = pars, model = modelN, error_list = error_list,set_sizes = set_sizes)
  return(out)
  
}


# Predict Hold-out Sets: 10F CV ------------------------------------------------------------
  
# Load all training fits
CVTraining <- readRDS("Fits/CV10Fold_TrainingfitsAlt.rds")

CVTraining <- CVTraining %>% 
  group_by(model,cvid) %>% 
  filter(Dev == min(Dev, na.rm = TRUE)) %>% 
  distinct() 
CVTraining$CVLL <- NA



# Load data to be fitted 

CVTest <- readRDS("CVData/testforCV10Fold.rds") %>% unnest() %>% filter(set_size != 1) %>% 
  nest_by(exp,id,keep=T)

# loops 
PredictCV <- NULL


for (subj in unique(CVTraining$id)){
 
  SubjTrain <- CVTraining %>% ungroup() %>% filter(id == subj)
  SubjTest <- CVTest %>% unnest() %>% filter(id == subj)
  
  #unique(CVTraining$model)
  for (modelN in c("MK_RNminus","MK_RNplus")){ # add models to be predicted here

    holdouts <- unique(SubjTrain$leftout)
    for (i in seq_along(holdouts)){


      TrainData <- SubjTrain %>% filter(model ==  modelN) %>%  
        filter(leftout == holdouts[[i]])
      if (nrow(TrainData) > 0){
        HoldOutData <- SubjTest %>% filter(Fold == i)
      
        TrainData$CVLL <- PredictHoldOutSz(TrainData,modelN,HoldOutData)
        print(TrainData$CVLL)
        PredictCV <- bind_rows(PredictCV,TrainData)
      }
    }
    
  }
  
}


saveRDS(PredictCV,file = "Fits/CV10FoldAlt_TestLL.rds")

# Predict Hold-out Sets: 5x2 CV ------------------------------------------------------------

# Load best fit parameters

CVTraining <- readRDS("Fits/CV5x2Fold_Trainingfits.rds")

CVTraining <- CVTraining %>% 
  group_by(model,cvid) %>% 
  filter(Dev == min(Dev, na.rm = TRUE)) %>% 
  distinct() 
CVTraining$CVLL <- NA

# Load data to be fitted 

CVTest <- readRDS("CVData/testforCV5times2.rds")

# loops 
PredictCV <- NULL

for (subj in unique(CVTraining %>% .$id)){
  
  SubjTrain <- CVTraining %>% ungroup() %>% filter(id == subj)
  SubjTest <- CVTest %>% unnest() %>% filter(id == subj)
  
  for (modelN in unique(CVTraining$model)){
    
    holdouts <- unique(SubjTrain$leftout)
    for (i in seq_along(holdouts)){
      
      
      TrainData <- SubjTrain %>% 
        filter(model ==  modelN) %>% 
        filter(leftout == holdouts[[i]])
      HoldOutData <- SubjTest %>%
        filter(.,grepl(substring(holdouts[[i]],1,1),CVRep)) %>% 
        filter(.,grepl(substring(holdouts[[i]],2,2),Fold))
      
      TrainData$CVLL <- PredictHoldOutSz(TrainData,modelN,HoldOutData)
      PredictCV <- bind_rows(PredictCV,TrainData)
    }
    
  }
  
}


#saveRDS(PredictCV,file = "Fits/CV5x2_TestLL.rds")

# Predict Holdout sets: leave-one-setsize-out ----------------------------

# Load best fit parameters

CVTraining <- readRDS("Fits/CVLOSzO_Trainingfits.rds")
CVTraining <- CVTraining %>% 
  group_by(model,cvid) %>% 
  filter(Dev == min(Dev, na.rm = TRUE)) %>% 
  distinct() 
CVTraining$CVLL <- NA


# Load data to be fitted 

CVTest <- readRDS("CVData/testforCVLOSZO.rds")

# loops 

PredictCV <- NULL

for (subj in unique(CVTraining$id)){
  
  SubjTrain <- CVTraining %>% ungroup() %>% filter(id == subj)
  SubjTest <- CVTest %>% unnest() %>% filter(id == subj)
  
  for (modelN in unique(CVTraining$model)){
    
    set_sizes <- unique(SubjTrain$leftout)
    for (i in seq_along(set_sizes)){
      
      TrainData <- SubjTrain %>% 
        filter(model ==  modelN) %>%  
        filter(leftout == set_sizes[[i]])
      HoldOutData <- SubjTest %>% 
        filter(set_size == set_sizes[[i]])
      
      TrainData$CVLL <- PredictHoldOutSz(TrainData,modelN,HoldOutData)
      TrainData$set_size <- set_sizes[[i]]
      PredictCV <- bind_rows(PredictCV,TrainData)
    }
    
  }
  
  
  
}



#saveRDS(PredictCV,file = "Fits/CVLOSzO_TestLL.rds")


# Predict Holdout sets: leave-one-setsize-out AGG ----------------------------

# Load best fit parameters

CVTraining <- readRDS("Fits/FitFull_AggLOSSOTraining_secondbatch.rds")
CVTraining <- CVTraining %>% 
  group_by(model,cvid) %>% 
  filter(Dev == min(Dev, na.rm = TRUE)) %>% 
  distinct() 
CVTraining$CVLL <- NA


# Load data to be fitted 

CVTest <- readRDS("CVData/testforCVLOSZO.rds")

# loops 

PredictCV <- NULL

for (subj in unique(CVTraining$exp)){

  
  SubjTrain <- CVTraining %>% ungroup() %>% filter(exp == subj)
  SubjTest <- CVTest %>% unnest() %>% filter(exp == subj)
  
  for (modelN in unique(CVTraining$model)){
    
    set_sizes <- unique(SubjTrain$leftout)
    for (i in seq_along(set_sizes)){
      
      TrainData <- SubjTrain %>% 
        filter(model ==  modelN) %>%  
        filter(leftout == set_sizes[[i]])
      HoldOutData <- SubjTest %>% 
        filter(set_size == set_sizes[[i]])
      
      TrainData$CVLL <- PredictHoldOutSz(TrainData,modelN,HoldOutData)
      TrainData$set_size <- set_sizes[[i]]
      PredictCV <- bind_rows(PredictCV,TrainData)
    }
    
  }
  
  
  
}



saveRDS(PredictCV,file = "Fits/FitFull_AggLOSSO_TestLL_secondbatch.rds")


# PredictCV for inidvidual set-sizes -------------------------------------------
# if hold-out set prediction is split by set sizes to compare directly to LOSsO-CV

# # Load best fit parameters
# 
# CVTraining <- readRDS("CVSz_TrainingBestFit.rds")
# CVTraining$CVLL <- NA
# 
# # Load data to be fitted 
# 
# CVTest <- readRDS("testforCVLOSZO.rds")
# 
# PredictCV <- NULL
# 
# 
# for (subj in unique(CVTraining$id)){
# 
#   
#   SubjTrain <- CVTraining %>% ungroup() %>% filter(id == subj)
#   SubjTest <- CVTest %>% unnest() %>% filter(id == subj)
#   
#   for (modelN in models){
#     
# 
#     holdouts <- unique(SubjTrain$leftout)
#     for (i in seq_along(holdouts)){
#       
#       for (sz in unique(SubjTest$set_size)){
#       
#       TrainData <- SubjTrain %>% filter(model ==  modelN) %>%  filter(leftout == holdouts[[i]])
#       HoldOutData <- SubjTest %>% filter(Fold == i) %>% filter(set_size == sz)
#       
#       TrainData$CVLL <- PredictHoldOutSz(TrainData,modelN,HoldOutData)
#       TrainData$set_size <- sz
#       PredictCV <- bind_rows(PredictCV,TrainData)
#       }
#     }
#     
#   }
#   
# }
# 
# 
# 
# 
# 
# # Load best fit parameters
# 
# CVTraining <- readRDS("CV5x2_TrainingBestFit.rds")
# CVTraining$CVLL <- NA
# 
# # Load data to be fitted 
# 
# CVTest <- readRDS("testforCV5times2.rds")
# 
# # loops 
# PredictCV <- NULL
# 
# 
# for (subj in unique(CVTraining$id)){
#   
#   SubjTrain <- CVTraining %>% ungroup() %>% filter(id == subj)
#   SubjTest <- CVTest %>% unnest() %>% filter(id == subj)
#   
#   for (modelN in models){
#     
#     holdouts <- unique(SubjTrain$leftout)
#     for (i in seq_along(holdouts)){
#       for (sz in unique(SubjTest$set_size)){
#       
#       TrainData <- SubjTrain %>% filter(model ==  modelN) %>%  filter(leftout == holdouts[[i]])
#       HoldOutData <- SubjTest %>% 
#         filter(.,grepl(substring(holdouts[[i]],1,1),CVRep)) %>% 
#         filter(.,grepl(substring(holdouts[[i]],2,2),Fold)) %>% 
#         filter(set_size == sz)
#       
#       TrainData$set_size <- sz
#       TrainData$CVLL <- PredictHoldOutSz(TrainData,modelN,HoldOutData)
#       PredictCV <- bind_rows(PredictCV,TrainData)
#       }
#     }
#     
#   }
#   
# }
# 
# 
# 
# 
