#RNTheoPred graph

library("doParallel")
library("vwmvp")
library("magrittr")
library("dplyr")
library("tibble")
library("tidyr")
library("ggplot2")
library(magick)
library(cowplot)
library(ungeviz)


# Sumstat function --------------------------------

MakeSumStat <- function(Prediction){
  CV <- NULL
  for(SS in unique(Prediction %>% .$setsize)){

    CircVars <- NULL

    for(modeln in unique(Prediction$model)){


      testNA <- any(is.na(Prediction %>%
                            filter(setsize == SS) %>%
                            filter(model == modeln) %>%
                            .$normprediction) == TRUE)
      if (testNA){

        CircVar <- NA
        CircSD <- NA
        MAE <- NA

      } else {

        samples <- sample(size=6e4,Prediction %>%
                            filter(setsize == SS) %>%
                            filter(model == modeln) %>%
                            .$data,replace=TRUE,
                          prob =  Prediction %>%
                            filter(setsize == SS) %>%
                            filter(model == modeln) %>%
                            .$normprediction)
        CircVar <- circular::var.circular(samples)
        CircSD <- circular::sd.circular(samples)
        MAE <- mean(abs(samples))
        kurtosis <- vwmvp::circkurtosis(samples)

      }

      out <- as_tibble(data.frame(setsize = SS,
                                  model = modeln,
                                  Var = CircVar,
                                  SD = CircSD,
                                  meanErr = MAE,
                                  kurt = kurtosis))

      CircVars <- CircVars %>%  bind_rows(out)
    }


    CV <- CV %>% bind_rows(CircVars)
  }


  sumstat <- CV %>% select(-SD) %>%
    pivot_longer(!c(model,setsize),names_to = "measure",values_to="value") %>%
    mutate(measure = factor(measure, levels = c("Var","meanErr","kurt"),
                            labels = c(expression(paste("circular ",sigma^2)),
                                       expression(paste("mean abs. deviation (radian)")),
                                       expression(paste("kurtosis")))))

  return(sumstat)
}

# Generate theoretical predictions ---------------------------------------------

## Note: for the uniform capacity limit, the limit is defined formally as sampling from values from 0 to 2K
## with the uniform distribution being defined by a mean of K (shown in the supplementary material).
## In the code in the vwmvp package it is practically implemented as ranging from 0 to K, hence a mean of K/2.
## To simulate data from F, P and U models with the equivalent mean-K, we therefore multiplied the input K parameter
## by 2 in the model simulations below.

rads <- c(-180:180) * pi / 180
rads <- ifelse(rads < -pi, rads + 2 * pi, ifelse(rads >
                                                   pi, rads - 2 * pi, rads))

## Fixed limit

dat_MKRNminus <- vwmvp::predict_data(model = "MK_RNminus",
                                     pars = c(60,1.6,15),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNminus")

dat_MKFMRNminus1 <- vwmvp::predict_data(model = "MK_FM_RNminus",
                                     pars = c(60,1.6,15,4),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_FM_RNminus1")

dat_MKFMRNminus2 <- vwmvp::predict_data(model = "MK_FM_RNminus",
                                        pars = c(60,1.6,15,2),
                                        data = data.frame(id = "Test",
                                                          set_size = rep(expsetsize,each=length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_FM_RNminus2")

dat_MKFMRNminus3 <- vwmvp::predict_data(model = "MK_FM_RNminus",
                                        pars = c(60,1.9,15,4),
                                        data = data.frame(id = "Test",
                                                          set_size = rep(expsetsize,each=length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_FM_RNminus3")


gendatF <- bind_rows(dat_MKRNminus,dat_MKFMRNminus1,dat_MKFMRNminus2,dat_MKFMRNminus3) %>%
  mutate(model = factor(model,levels = c("MK_RNminus","MK_FM_RNminus1","MK_FM_RNminus2","MK_FM_RNminus3"))) %>%
  mutate(normprediction = prediction/sum(prediction)) %>%
  mutate(setsize = paste0("Set size ", setsize))

# Variable limit (uniform)

dat_MKRNminus <- vwmvp::predict_data(model = "MK_RNminus",
                                     pars = c(60,1.6,15),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNminus")

dat_MKURNminus1 <- vwmvp::predict_data(model = "MK_U_RNminus",
                                        pars = c(60,1.6,15,4*2),
                                        data = data.frame(id = "Test",
                                                          set_size = rep(expsetsize,each=length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_U_RNminus1")

dat_MKURNminus2 <- vwmvp::predict_data(model = "MK_U_RNminus",
                                        pars = c(60,1.6,15,2*2),
                                        data = data.frame(id = "Test",
                                                          set_size = rep(expsetsize,each=length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_U_RNminus2")

dat_MKURNminus3 <- vwmvp::predict_data(model = "MK_U_RNminus",
                                        pars = c(60,1.9,15,4*2),
                                        data = data.frame(id = "Test",
                                                          set_size = rep(expsetsize,each=length(rads)),
                                                          error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_U_RNminus3")

gendatU <- bind_rows(dat_MKRNminus,dat_MKURNminus1,dat_MKURNminus2,dat_MKURNminus3) %>%
  mutate(model = factor(model,levels = c("MK_RNminus","MK_U_RNminus1","MK_U_RNminus2","MK_U_RNminus3"))) %>%
  mutate(normprediction = prediction/sum(prediction)) %>%
  mutate(setsize = paste0("Set size ", setsize))

## Variable limit (poisson)

dat_MKRNminus <- vwmvp::predict_data(model = "MK_RNminus",
                                     pars = c(60,1.6,15),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNminus")

dat_MKPRNminus1 <- vwmvp::predict_data(model = "MK_P_RNminus",
                                       pars = c(60,1.6,15,4),
                                       data = data.frame(id = "Test",
                                                         set_size = rep(expsetsize,each=length(rads)),
                                                         error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_P_RNminus1")

dat_MKPRNminus2 <- vwmvp::predict_data(model = "MK_P_RNminus",
                                       pars = c(60,1.6,15,2),
                                       data = data.frame(id = "Test",
                                                         set_size = rep(expsetsize,each=length(rads)),
                                                         error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_P_RNminus2")

dat_MKPRNminus3 <- vwmvp::predict_data(model = "MK_P_RNminus",
                                       pars = c(60,1.9,15,4),
                                       data = data.frame(id = "Test",
                                                         set_size = rep(expsetsize,each=length(rads)),
                                                         error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_P_RNminus3")

gendatP <- bind_rows(dat_MKRNminus,dat_MKPRNminus1,dat_MKPRNminus2,dat_MKPRNminus3) %>%
  mutate(model = factor(model,levels = c("MK_RNminus","MK_P_RNminus1","MK_P_RNminus2","MK_P_RNminus3"))) %>%
  mutate(normprediction = prediction/sum(prediction)) %>%
  mutate(setsize = paste0("Set size ", setsize))

## One of each vs unlimited

dat_MKRNminus <- vwmvp::predict_data(model = "MK_RNminus",
                                     pars = c(60,1.6,15),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNminus")

dat_MKFMRNminus <- vwmvp::predict_data(model = "MK_FM_RNminus",
                                       pars = c(60,1.6,15,4),
                                       data = data.frame(id = "Test",
                                                         set_size = rep(expsetsize,each=length(rads)),
                                                         error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_FM_RNminus")

dat_MKPRNminus <- vwmvp::predict_data(model = "MK_P_RNminus",
                                       pars = c(60,1.6,15,4),
                                       data = data.frame(id = "Test",
                                                         set_size = rep(expsetsize,each=length(rads)),
                                                         error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_P_RNminus")

dat_MKURNminus <- vwmvp::predict_data(model = "MK_U_RNminus",
                                       pars = c(60,1.6,15,4 * 2),
                                       data = data.frame(id = "Test",
                                                         set_size = rep(expsetsize,each=length(rads)),
                                                         error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_U_RNminus")


gendatAll <- bind_rows(dat_MKRNminus,dat_MKFMRNminus,dat_MKPRNminus,dat_MKURNminus) %>%
  mutate(model = factor(model,levels = c("MK_RNminus","MK_FM_RNminus","MK_P_RNminus","MK_U_RNminus"))) %>%
  mutate(normprediction = prediction/sum(prediction)) %>%
  mutate(setsize = paste0("Set size ", setsize))


# Pminus --------

sumstatP <- MakeSumStat(gendatP)

errordistributionP <-ggplot(data = gendatP,
                            aes(x = data,y=prediction,color=model,linetype=model)) +
  geom_line(size = 1) +

  facet_wrap(.~setsize,ncol=4) +
  coord_cartesian(xlim = c(-pi/2,pi/2),ylim=c(0,NA)) +

  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_point(aes(shape=model,color=model,x=0,y=-1),size=3) +
  scale_color_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  scale_shape_manual(values = c(16,15,15,15),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +
  scale_linetype_manual(values = c("solid",rep(c("dashed"),3)),
                        labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = ",infinity,")")),
                                   expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 4)")),
                                   expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 2)")),
                                   expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.9,",
                                                    tau," = 15,",
                                                    K," = 4)")))) +
  theme(axis.text.x = element_text(size=12),
        plot.title = element_text(size = 12,hjust=0.5),
        axis.text.y = element_text(size=12,color="white"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill="transparent"),
        legend.position = c(0.6,0.27),
        legend.box.background = element_rect(fill="white",color="lightgrey"),
        legend.text = element_text(size=12),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))


sumstat_plotP <- ggplot(sumstatP, aes(x = setsize, y = value)) +
  geom_point(data = sumstatP,
             aes(x = setsize, y = value,shape = model, fill = model,alpha=model),
             size=2.5,color="black",position=position_dodge(0.6)) +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_wrap(.~measure,ncol=3,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25),limits=c(0,NA))+
  scale_alpha_manual(values=c(1,0.7,0.7,0.7))+
  scale_fill_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                    labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = ",infinity,")")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 4)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 2)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.9,",
                                                tau," = 15,",
                                                K," = 4)")))) +
  scale_shape_manual(values = c(21,22,22,22),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  theme(axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "white",color="black"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")



CapP <- plot_grid(errordistributionP,sumstat_plotP,labels = c('(1)','(2)'), label_size = 12,ncol=1,
                  rel_heights=c(2,1), label_x = 0, label_y = 0.95,hjust = 0.2)

#Fminus -------------

sumstatF <- MakeSumStat(gendatF)

errordistributionF <- ggplot(data = gendatF,
                             aes(x = data,y=prediction,color=model,linetype=model)) +
  geom_line(size = 1) +

  facet_wrap(.~setsize,ncol=4) +
  coord_cartesian(xlim = c(-pi/2,pi/2),ylim=c(0,NA)) +

  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_point(aes(shape=model,color=model,x=0,y=-1),size=3) +
  scale_color_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  scale_shape_manual(values = c(16,15,15,15),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +
  scale_linetype_manual(#values = rep(c("solid"),4),
    values = c("solid",rep("dashed",3)),
    labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                alpha," = 1.6,",
                                tau," = 15,",
                                K," = ",infinity,")")),
               expression(paste("VP(",kappa,")F-, fixed limit (",
                                bar(kappa)[1]," = 60,",
                                alpha," = 1.6,",
                                tau," = 15,",
                                K," = 4)")),
               expression(paste("VP(",kappa,")F-, fixed limit (",
                                bar(kappa)[1]," = 60,",
                                alpha," = 1.6,",
                                tau," = 15,",
                                K," = 2)")),
               expression(paste("VP(",kappa,")F-, fixed limit (",
                                bar(kappa)[1]," = 60,",
                                alpha," = 1.9,",
                                tau," = 15,",
                                K," = 4)")))) +
  theme(axis.text.x = element_text(size=12),
        plot.title = element_text(size = 12,hjust=0.5),
        axis.text.y = element_text(size=12,color="white"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill="transparent",color="lightgrey"),
        legend.position = c(0.65,0.27),
        legend.box.background = element_rect(fill="white"),
        legend.text = element_text(size=12),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))


sumstat_plotF <- ggplot(sumstatF, aes(x = setsize, y = value)) +
  geom_point(data = sumstatF,
             aes(x = setsize, y = value,shape = model, fill = model,alpha=model),
             size=2.5,color="black",position=position_dodge(0.6)) +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_wrap(.~measure,ncol=3,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25),limits=c(0,NA))+
  scale_alpha_manual(values=c(1,0.7,0.7,0.7))+
  scale_fill_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                    labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = ",infinity,")")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 4)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 2)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.9,",
                                                tau," = 15,",
                                                K," = 4)")))) +
  scale_shape_manual(values = c(21,22,22,22),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  theme(axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "white",color="black"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")



CapF <- plot_grid(errordistributionF,sumstat_plotF,labels = c('(1)','(2)'), label_size = 12,ncol=1,
                  rel_heights=c(2,1), label_x = 0, label_y = 0.95,hjust = 0.2)

#Uminus -------------


sumstatU <- MakeSumStat(gendatU)

errordistributionU <- ggplot(data = gendatU,
                             aes(x = data,y=prediction,color=model,linetype=model)) +
  geom_line(size = 1) +

  facet_wrap(.~setsize,ncol=4) +
  coord_cartesian(xlim = c(-pi/2,pi/2),ylim=c(0,NA)) +

  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_point(aes(shape=model,color=model,x=0,y=-1),size=3) +
  scale_color_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  scale_shape_manual(values = c(16,15,15,15),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +
  scale_linetype_manual(values = c("solid",rep("dashed",3)),
                        labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = ",infinity,")")),
                                   expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 4)")),
                                   expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 2)")),
                                   expression(paste("VP(",kappa,")U-, unif.-dist. limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.9,",
                                                    tau," = 15,",
                                                    K," = 4)")))) +
  theme(axis.text.x = element_text(size=12),
        plot.title = element_text(size = 12,hjust=0.5),
        axis.text.y = element_text(size=12,color="white"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill="transparent",color="lightgrey"),
        legend.position = c(0.6,0.27),
        legend.box.background = element_rect(fill="white"),
        legend.text = element_text(size=12),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))


sumstat_plotU <- ggplot(sumstatU, aes(x = setsize, y = value)) +
  geom_point(data = sumstatU,
             aes(x = setsize, y = value,shape = model, fill = model,alpha=model),
             size=2.5,color="black",position=position_dodge(0.6)) +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_wrap(.~measure,ncol=3,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25),limits=c(0,NA))+
  scale_alpha_manual(values=c(1,0.7,0.7,0.7))+
  scale_fill_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                    labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = ",infinity,")")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 4)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 2)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.9,",
                                                tau," = 15,",
                                                K," = 4)")))) +
  scale_shape_manual(values = c(21,22,22,22),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  theme(axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "white",color="lightgrey"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")



CapU <- plot_grid(errordistributionU,sumstat_plotU,labels = c('(1)','(2)'), label_size = 12,ncol=1,
                  rel_heights=c(2,1), label_x = 0, label_y = 0.95,hjust = 0.2)


#Allminus -------------

sumstatA <- MakeSumStat(gendatAll)

errordistributionA <- ggplot(data = gendatAll,
                             aes(x = data,y=prediction,color=model,linetype=model)) +
  geom_line(size = 1) +

  facet_wrap(.~setsize,ncol=4) +
  coord_cartesian(xlim = c(-pi/2,pi/2),ylim=c(0,NA)) +

  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +
  geom_point(aes(shape=model,color=model,x=0,y=-1),size=3) +

  scale_color_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")U-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  scale_shape_manual(values = c(16,15,15,15),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")F-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")U-, fixed limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +
  scale_linetype_manual(values = rep(c("solid"),4),
                        # values = c("solid",rep("dashed",3)),
                        labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = ",infinity,")")),
                                   expression(paste("VP(",kappa,")F-, fixed limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 4)")),
                                   expression(paste("VP(",kappa,")P-, fixed limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 4)")),
                                   expression(paste("VP(",kappa,")U-, fixed limit (",
                                                    bar(kappa)[1]," = 60,",
                                                    alpha," = 1.6,",
                                                    tau," = 15,",
                                                    K," = 4)")))) +
  theme(axis.text.x = element_text(size=12),
        plot.title = element_text(size = 12,hjust=0.5),
        axis.text.y = element_text(size=12,color="white"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=12),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill="transparent",color="lightgrey"),
        legend.position = c(0.65,0.27),
        legend.box.background = element_rect(fill="white"),
        legend.text = element_text(size=12),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))


sumstat_plotA <- ggplot(sumstatA, aes(x = setsize, y = value)) +
  geom_point(data = sumstatA,
             aes(x = setsize, y = value,shape = model, fill = model,alpha=model),
             size=2.5,color="black",position=position_dodge(0.6)) +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_wrap(.~measure,ncol=3,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25),limits=c(0,NA))+
  scale_alpha_manual(values=c(1,0.7,0.7,0.7))+
  scale_fill_manual(values = c("#332288","#E69F00","#44AA99","#CC79A7"),
                    labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = ",infinity,")")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 4)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.6,",
                                                tau," = 15,",
                                                K," = 2)")),
                               expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                bar(kappa)[1]," = 60,",
                                                alpha," = 1.9,",
                                                tau," = 15,",
                                                K," = 4)")))) +
  scale_shape_manual(values = c(21,22,22,22),
                     labels = c(expression(paste("VP(",kappa,")A-, unlimited (",bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = ",infinity,")")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 4)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.6,",
                                                 tau," = 15,",
                                                 K," = 2)")),
                                expression(paste("VP(",kappa,")P-, poisson-dist. limit (",
                                                 bar(kappa)[1]," = 60,",
                                                 alpha," = 1.9,",
                                                 tau," = 15,",
                                                 K," = 4)")))) +

  theme(axis.text = element_text(size=12),
        axis.title.y = element_blank(),
        axis.title.x = element_text(size=12),
        plot.title = element_text(size = 18,hjust=0.5),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 12),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill = "white",color="black"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = "none",
        legend.text.align = 0,
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")



CapA <- plot_grid(errordistributionA,sumstat_plotA,labels = c('(1)','(2)'), label_size = 11,ncol=1,
                  rel_heights=c(2,1), label_x = 0, label_y = 0.95,hjust = 0.2)



# CP fin ---------------------

plot_grid(CapA,CapF,CapP,CapU,labels="AUTO",ncol=2,scale=0.95)

#ggsave("RQ2_TheoPred.png", units="cm", width=40, height=40, dpi=600)
