#install.packages("vwmvp_1.0.tar.gz", repos = NULL, type="source",dependencies = TRUE)


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
library(ggforce)


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

# RQ1 Theoretical Predictions  ---------------------


# # Plot Fisher information vs kappa in introduction

FisherKappa <- tibble(J = vwmvp::J_map,kappa = vwmvp::kappa_map)  %>% mutate(Jpercent = (J )/kappa)

FisherKappa <- tibble(J = c(seq(0.1,1,0.1),seq(1.5,5,0.5),c(1:100)),kappa = KappafromJ((J))) %>%
  mutate(Jpercent = (J )/kappa)


rawFisher <- ggplot(FisherKappa %>% filter(J<=50), aes(x=kappa,y =J)) +
  geom_point() +
  geom_abline(slope=1)+
  scale_x_continuous(name = as.expression(bquote(kappa~"="~Phi(J)))
                    ) +
  scale_y_continuous(name = expression(paste(kappa~"="~f(phi)))) +
  facet_zoom(xlim = c(0,5), ylim = c(0,5), zoom.size = 1) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = c(0.9,0.75))
tail(FisherKappa)


percentFisher <- ggplot(FisherKappa , aes(x=kappa,y =Jpercent)) +
  geom_point() +
  #geom_abline(slope=1)+
  scale_y_continuous(name = expression(paste(frac(kappa~"="~f(phi), kappa~"="~Phi(J)))),
                     limits = c(0,1),breaks=seq(0,1,.25),labels=seq(0,1,.25)) +
  scale_x_continuous(name =  as.expression(bquote(kappa~"="~Phi(J))),
                     limits = c(0,50),breaks=seq(0,50,10)) +
  #facet_zoom(xlim = c(0,5), ylim = c(0,5), zoom.size = 1) +
  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        panel.background = element_rect(fill = "white"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        legend.position = c(0.9,0.75))

FisherKappa <- plot_grid(rawFisher,percentFisher,ncol=2, labels = c('(1)','(2)'),
                         label_size = 12,rel_widths=c(1.9,1),scale=.95,label_x = 0, label_y = 0.95)

# Generate predictions

rads <- c(-180:180) * pi / 180
rads <- ifelse(rads < -pi, rads + 2 * pi, ifelse(rads >
                                                   pi, rads - 2 * pi, rads))

dat_MKRNminus <- vwmvp::predict_data(model = "MK_RNminus",
                                    pars = c(30,1.5,15),
                                    data = data.frame(id = "Test",
                                                      set_size = rep(expsetsize,each=length(rads)),
                                                      error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNminus")

dat_MKRNplus0 <- vwmvp::predict_data(model = "MK_RNplus",
                                     pars = c(30,1.5,15,15),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNplus0")

dat_MKRNplus1 <- vwmvp::predict_data(model = "MK_RNplus",
                                     pars = c(20,1.5,15,35),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNplus1")

dat_MKRNplus2 <- vwmvp::predict_data(model = "MK_RNplus",
                                     pars = c(90,1.8,15,35),
                                     data = data.frame(id = "Test",
                                                       set_size = rep(expsetsize,each=length(rads)),
                                                       error_0 = rep(rads,length(expsetsize)))) %>%
  mutate(model = "MK_RNplus2")


gendat <- bind_rows(dat_MKRNminus,dat_MKRNplus0,dat_MKRNplus1,dat_MKRNplus2) %>%
  mutate(setsize = paste0("Set size ", setsize)) %>% mutate(normprediction = prediction/sum(prediction))

errordistribution <- ggplot(data = gendat,
                            aes(x = data,y=prediction,color=model,linetype=model)) +
  geom_line(size = 1) +

  facet_wrap(.~setsize,ncol=4) +
  #geom_point(aes(shape=model,color=model,x=0,y=-1),size=3) +
  coord_cartesian(xlim = c(-pi/2,pi/2),ylim=c(0,NA)) +

  scale_x_continuous(limits=c(-pi,pi),name = "Error (radian)",breaks = c(pi,-pi/2,0,pi/2,pi),
                     labels = c(expression(-pi),
                                expression(-frac(pi,2)),
                                0,
                                expression(frac(pi,2)),
                                expression(pi))) +

  scale_color_manual(values =  c("#000000","#E69F00", "#56B4E9", "#009E73"),
                     labels = c(expression(paste("VP(",kappa,")A- (",kappa," = 30,",
                                                 tau," = 15,",
                                                 alpha," = 1.5)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 30,",
                                                 tau," = 15,",
                                                 alpha," = 1.5,",
                                                 kappa[r]," = 15)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 20,",
                                                 tau," = 15,",
                                                 alpha," = 1.5,",
                                                 kappa[r]," = 35)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 90,",
                                                 tau," = 15,",
                                                 alpha," = 1.8,",
                                                 kappa[r]," = 35)")))
  ) +
  scale_shape_manual(values = c(16,15,15,15),
                     labels = c(expression(paste("VP(",kappa,")A- (",kappa," = 30,",
                                                 tau," = 15,",
                                                 alpha," = 1.5)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 30,",
                                                 tau," = 15,",
                                                 alpha," = 1.5,",
                                                 kappa[r]," = 15)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 20,",
                                                 tau," = 15,",
                                                 alpha," = 1.5,",
                                                 kappa[r]," = 35)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 90,",
                                                 tau," = 15,",
                                                 alpha," = 1.8,",
                                                 kappa[r]," = 35)")))
                     ) +

  scale_linetype_manual(values = c("solid","dashed","dashed","dashed"),
                        labels = c(expression(paste("VP(",kappa,")A- (",kappa," = 30,",
                                                    tau," = 15,",
                                                    alpha," = 1.5)")),
                                   expression(paste("VP(",kappa,")A+ (",kappa," = 30,",
                                                    tau," = 15,",
                                                    alpha," = 1.5,",
                                                    kappa[r]," = 15)")),
                                   expression(paste("VP(",kappa,")A+ (",kappa," = 20,",
                                                    tau," = 15,",
                                                    alpha," = 1.5,",
                                                    kappa[r]," = 35)")),
                                   expression(paste("VP(",kappa,")A+ (",kappa," = 90,",
                                                    tau," = 15,",
                                                    alpha," = 1.8,",
                                                    kappa[r]," = 35)")))
  ) +

  theme(axis.text.x = element_text(size=12),
        plot.title = element_text(size = 12,hjust=0.5),
        axis.text.y = element_text(size=12,color="white"),
        axis.ticks.y = element_blank(),
        axis.title.x = element_text(size=14),
        axis.title.y = element_blank(),
        panel.background = element_rect(fill = "transparent"),
        panel.border = element_rect(colour = "black", fill=NA, size=0.5),
        strip.background = element_rect(fill = "transparent"),
        plot.background = element_rect(fill = "transparent"),
        strip.text = element_text(size = 14),
        legend.key = element_rect(colour = "transparent", fill = "white"),
        legend.background = element_rect(fill="white",color="lightgrey"),
        legend.position = c(0.7,0.27),
        legend.box.background = element_rect(fill="white"),
        legend.text = element_text(size=14),
        legend.text.align = 0,
        legend.title = element_blank(),
        legend.key.width = unit(3,"line"))

errordistribution


sumstat <- MakeSumStat(gendat)


sumstat_plot <- ggplot(sumstat, aes(x = set_size, y = value)) +
  geom_point(data = sumstat,
             aes(x = setsize, y = value,shape = model, fill = model),
             size=3.5,alpha=0.7,color="black") +
  scale_x_discrete(labels=c(c(1:8)),name="Set size") +
  facet_grid(.~measure,scales="free_y",labeller=label_parsed) +
  scale_y_continuous(name = "Summary Statistic",breaks = c(0,0.25,0.5,0.75,1,1.25,1.5),limits=c(-0.01,1.5))+


  scale_fill_manual(values = c("#000000","#E69F00", "#56B4E9", "#009E73"),
                    labels = c(expression(paste("VP(",kappa,")A- (",kappa," = 30,",
                                                tau," = 15,",
                                                alpha," = 1.5)")),
                               expression(paste("VP(",kappa,")A+ (",kappa," = 30,",
                                                tau," = 15,",
                                                alpha," = 1.5,",
                                                kappa[r]," = 15)")),
                               expression(paste("VP(",kappa,")A+ (",kappa," = 20,",
                                                tau," = 15,",
                                                alpha," = 1.5,",
                                                kappa[r]," = 35)")),
                               expression(paste("VP(",kappa,")A+ (",kappa," = 90,",
                                                tau," = 15,",
                                                alpha," = 1.8,",
                                                kappa[r]," = 35)")))
                    ) +
  scale_shape_manual(values = c(21,22,22,22),
                     labels = c(expression(paste("VP(",kappa,")A- (",kappa," = 30,",
                                                 tau," = 15,",
                                                 alpha," = 1.5)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 30,",
                                                 tau," = 15,",
                                                 alpha," = 1.5,",
                                                 kappa[r]," = 15)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 20,",
                                                 tau," = 15,",
                                                 alpha," = 1.5,",
                                                 kappa[r]," = 35)")),
                                expression(paste("VP(",kappa,")A+ (",kappa," = 90,",
                                                 tau," = 15,",
                                                 alpha," = 1.8,",
                                                 kappa[r]," = 35)")))
                     ) +
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
        legend.background = element_rect(fill = "transparent"),
        legend.text = element_text(size = 12,face="plain"),
        legend.position = c(0.15,0.8),
        legend.text.align = 0,
        legend.box.background = element_rect(fill="white",color="lightgrey"),
        legend.key.height = unit(0.5,"line"),
        legend.title = element_blank(),
        legend.box = "vertical")
sumstat_plot

plot_grid(FisherKappa,errordistribution,sumstat_plot,labels = c('A','B','C'), ncol=1,
          rel_heights=c(1,1.5,1))

#ggsave("RQ1_TheoPred.png", units="cm", width=35, height=35, dpi=600)


