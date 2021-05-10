library("tidyverse")
library("circular")


load("Data/ol17_prepared.rda")
load("Data/vdb14_prepared.rda")
load("Data/pratte17_prepared.rda")
experimentfile <- bind_rows(data_ol17[data_ol17$exp == "ol17_e1",],data_vdb14,data_pratte17)

# prep data -----

datprep <- experimentfile %>% 
  mutate(cvid = id,
         leftout = 0) %>% 
  group_nest(exp,cvid,keep=T)


## prep CV -----------------
# Leave one id out

trainforCV <- NULL

for (experiment in unique(experimentfile$exp)){

  dat_prep <- experimentfile %>%
    filter(exp == experiment)
  
  subj <- sort(unique(dat_prep$id))
  
  dtrain <-NULL
  

  for (i in seq_along(subj)) {
    
    dtrain1 <- dat_prep %>%  
      filter(id != subj[[i]]) %>%
      mutate(leftout = subj[[i]]) %>% 
      mutate(cvid = paste0(exp,"_LeftOut",leftout))
    
    dtrain <- bind_rows(dtrain,dtrain1) 
    
    
  }

  trainforCV <- bind_rows(trainforCV,dtrain) 
  
}

trainforCV <- trainforCV %>% 
  group_nest(exp,cvid,keep=T)
#saveRDS(trainforCV, file = "trainforCVLOIDO.rds")



# Leave one setsize out


# create dataframe of id_leftout combination to fit separately
trainforCV <- NULL

for (subj in unique(experimentfile$id)){
  
  dat_prep <- experimentfile %>% 
    filter(id == subj)
  
  expsetsizes <- sort(unique(dat_prep$set_size))
  
  dtrain <-NULL

  
  for (i in seq_along(expsetsizes)) {
    
    dtrain1 <- dat_prep %>% 
      filter(set_size != expsetsizes[[i]]) %>%
      mutate(leftout = expsetsizes[[i]]) %>% 
      mutate(cvid = paste0(id,"_LeftOut",leftout))
    
    dtrain <- bind_rows(dtrain,dtrain1) 

    
  }
  
  trainforCV <- bind_rows(trainforCV,dtrain) 

}

trainforCV <- trainforCV %>% 
  group_nest(exp,cvid,keep=T)
#saveRDS(trainforCV, file = "trainforCVLOSZO.rds")


# 5 x 2 CV

trainforCV <- NULL

for (subj in unique(experimentfile$id)){
  
  dat_prep <- experimentfile %>% 
    filter(id == subj)
  
  expsetsizes <- sort(unique(dat_prep$set_size))
  
  train <- NULL
  for (j in c(1:5)){
    dtrain <- NULL
    for (i in seq_along(expsetsizes)) {
      
      dtrain1 <- dat_prep %>% 
        filter(set_size == expsetsizes[[i]])
      
      dtrain2 <- dtrain1 %>% 
        mutate(Fold = sample(factor(rep(c("A","B"), length.out=nrow(dtrain1)))),
               CVRep = j) %>% 
        mutate(leftout = paste0(j,ifelse(Fold == "A","B","A"))) %>% 
        mutate(cvid = paste0(id,"_Leftout",leftout))
      dtrain <- bind_rows(dtrain,dtrain2) 
    }
    
    train <- bind_rows(train,dtrain)
    
    
  }
  trainforCV <- bind_rows(trainforCV,train)
}

trainforCV5times2 <- trainforCV %>% 
  group_nest(exp,cvid,keep=T)
#saveRDS(trainforCV5times2, file = "trainforCV5times2.rds")


# 10Fold CV

trainforCV <- NULL

for (subj in unique(experimentfile$id)){
  
  dat_prep <- experimentfile %>% 
    filter(id == subj)
  
  expsetsizes <- sort(unique(dat_prep$set_size))
  
  train <- NULL
  
  for (i in seq_along(expsetsizes)) {
    
    dtrain1 <- dat_prep %>% 
      filter(set_size == expsetsizes[[i]])
    
    dtrain2 <- dtrain1 %>% 
      mutate(Fold = sample(factor(rep(c(1:10), length.out=nrow(dtrain1)))))
    
    dtrain <- NULL
    
    for (j in c(1:10)) {
      
      dtrain1 <- dtrain2 %>% 
        filter(Fold != j) %>%
        mutate(leftout = j) %>% 
        mutate(cvid = paste0(id,"_LeftOutFold",leftout))
      
      dtrain <- bind_rows(dtrain,dtrain1) 
      
    }
    
    train <- bind_rows(train,dtrain) 
  }
  trainforCV <- bind_rows(trainforCV,train)
}

trainforCV10Fold <- trainforCV %>% 
  group_nest(exp,cvid,keep=T)
#saveRDS(trainforCV10Fold, file = "trainforCV10Fold.rds")

# Remove set sizes from CV10Fold for altCV10F

train10Fold <- readRDS("CVData/trainforCV10Fold.rds")

train10Fold_noSs1 <- train10Fold %>% unnest() %>% filter(set_size > 1)
train10Fold_noSs2 <- train10Fold %>% unnest() %>% filter(set_size != 2)
#saveRDS(train10Fold_noSs1, file = "trainforCV10Fold_noSs1.rds")
#saveRDS(train10Fold_noSs2, file = "CVData/trainforCV10Fold_noSs2.rds")


# Make appropriate test set files ---------------------------------------
#Train id


#TrainSz

testforCVSz <- experimentfile %>%
  mutate(cvid = id) %>% 
  mutate(leftout = NA) %>% 
  group_nest(exp,id,keep=T)
#saveRDS(testforCVSz, file = "testforCVLOSZO.rds")

# CV 5 x2

testforCV5times2 <- readRDS("trainforCV5times2.rds")
#saveRDS(testforCV5times2, file = "testforCV5times2.rds")

# CV 10 Fold

train10Fold <- readRDS("trainforCV10Fold.rds")

trainmainfold <- train10Fold %>% filter(.,grepl('LeftOutFold1',cvid)) %>% filter(.,!grepl('LeftOutFold10',cvid)) %>% 
  unnest()
trainfinalfold <- train10Fold %>% filter(.,grepl('LeftOutFold10',cvid)) %>% unnest() %>% filter(Fold == 1) 

test10Fold <- bind_rows(trainmainfold,trainfinalfold) %>% select(-exp1,-cvid1,-leftout,-cvid)

testforCV10Fold <- test10Fold %>% 
  group_nest(exp,id,keep=T)
#saveRDS(testforCV10Fold, file = "testforCV10Fold.rds")
   

# Prep by exp LOSsO -----------------------------

LOSSO<- readRDS("CVData/trainforCVLOSZO.rds")
LOSSOExp <- NULL
for (expe in unique(LOSSO$exp)){

  Exps <- LOSSO %>% filter(exp == expe) %>% unnest()
  expsetsize <- sort(unique(Exps %>% .$set_size))
  selectleftouts <- expsetsize[c(which(expsetsize == 1),which(expsetsize == min(expsetsize[expsetsize > 1])))]
  
  ExpC <- Exps %>% filter(leftout %in% selectleftouts) %>% 
    mutate(cvid = paste0(exp,"_Leftout",leftout),
           id = exp)
  
  LOSSOExp <- LOSSOExp %>%  bind_rows(ExpC)
}
models <- c("MK_RNplus","MK_RNminus","EP_RNplus","EP_RNminus")


dat_prep <- LOSSOExp %>% group_nest(exp,cvid,keep=T) %>% 
  slice(rep(1:n(),4)) %>% 
  mutate(modeln = rep(models,each=length(unique(LOSSOExp$cvid)))) %>% 
  mutate(uniqueid = paste0(cvid,"_",modeln))
