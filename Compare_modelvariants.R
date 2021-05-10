# Compare model implementation in manuscript with alternate model implementations

library("plyr")
library("magrittr")
library("tidyr")
library("xtable")
library("ggplot2")
library("dplyr")
library("RColorBrewer")
library(cowplot)

# Models to compare ------------------------------------------------------------

models <- c("MK_F_RNminus","MK_P_RNminus","MK_U_RNminus","MK_FM_RNminus",
            "MK_F_RNplus","MK_P_RNplus","MK_U_RNplus","MK_FM_RNplus",
            "MK_P2_RNminus","MK_U2_RNminus","MK_FM2_RNminus",
            "MK_P2_RNplus","MK_U2_RNplus","MK_FM2_RNplus")

# Load Fits --------------------------------------------------------------------

FittedFull <- readRDS("Fits/FitFull.rds") %>% filter(model %in% models)


# Identify best fit ------------------------------------------------------------


BestFull <- FittedFull %>%
  group_by(model,cvid) %>%
  filter(Dev == min(Dev, na.rm = TRUE)) %>% #finds best fit of 20 runs
  arrange(cvid) %>%
  distinct() %>%
  group_by(cvid) %>%
  mutate(delLL = LL - min(LL),
         delDev = Dev - min(Dev),
         delAIC = AIC - min(AIC)) %>%
  group_by(model)

## globally determines order of bars in graphs
modelsinorder <- BestFull %>%
  group_by(model) %>%
  summarize(means = mean(AIC)) %>%
  arrange(-means) %>%
  .$model

BestFull$model <- factor(BestFull$model,levels = modelsinorder)


# Manuscript-style graph -------------------------------------------------------

## Fitted Full

AICmeans <- BestFull %>%
  group_by(model) %>%
  summarize(means = mean(AIC)) %>%
  arrange(means)

#sum all aics
AICtot <- BestFull %>%
  group_by(cvid,model) %>%
  mutate(deltaAIC = AIC - AICmeans %>% slice(1) %>% .$means) %>%
  group_by(model) %>%
  summarize(mdeltaAIC = mean(deltaAIC)) %>%
  mutate(mdeltaAICadj = ifelse(mdeltaAIC > 20,20,mdeltaAIC)) %>%
  arrange(model)


allAIC <- ggplot(AICtot,aes(y = mdeltaAICadj,x = model,fill = mdeltaAIC)) +
  geom_rect(aes(xmin=0, xmax=14.5, ymin=0, ymax=1), fill = "grey") +
  geom_bar(stat = "identity") +
  coord_flip(ylim = c(0,20)) +
  scale_fill_gradientn(name = expression(paste(Delta," AIC")),
                       colors = rev(brewer.pal(n = 9, name = "YlOrRd")),
                       limits=c(0,20))+

  scale_y_continuous(name = expression(paste(Delta," AIC")),
                     breaks = seq(0,20,2),
                     labels = c(seq(0,18,2),"...")) +
  scale_x_discrete(name = "Model", labels = c(
    expression(paste("VP(",kappa,")U2-")),
    expression(paste("VP(",kappa,")U-")),
    expression(paste("VP(",kappa,")P-")),
    expression(paste("VP(",kappa,")P2-")),
    expression(paste("VP(",kappa,")F3+")),
    expression(paste("VP(",kappa,")F2+")),
    expression(paste("VP(",kappa,")F+")),
    expression(paste("VP(",kappa,")F2-")),
    expression(paste("VP(",kappa,")F-")),
    expression(paste("VP(",kappa,")F3-")),
    expression(paste("VP(",kappa,")U2+")),
    expression(paste("VP(",kappa,")P+")),
    expression(paste("VP(",kappa,")U+")),
    expression(paste("VP(",kappa,")P2+"))
  )) +

  theme(axis.text.x = element_text(size=12),
        axis.text.y = element_text(size=12),
        axis.title = element_text(size=12),
        legend.text = element_text(size = 12),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.position = c(0.7,0.7),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size=12),
        strip.placement = "outside",
        strip.text = element_text(size=12))


# Split by individual differences ----------------------------------------------

Refcompare <- BestFull %>%
  mutate(modeltype =
           case_when(model %in% c("MK_FM_RNminus", "MK_F_RNminus") ~ "F-",
                     model %in% c("MK_FM_RNplus","MK_F_RNplus") ~ "F+",
                     model %in% c("MK_FM2_RNplus") ~ "FM+",
                     model %in% c("MK_FM2_RNminus") ~ "FM-",
                     model %in% c("MK_P_RNplus","MK_P2_RNplus") ~ "P+",
                     model %in% c("MK_U_RNplus","MK_U2_RNplus") ~ "U+",
                     model %in% c("MK_P_RNminus","MK_P2_RNminus") ~ "P-",
                     model %in% c("MK_U_RNminus","MK_U2_RNminus") ~ "U-")) %>%
  mutate(reference =
           case_when(model %in% c("MK_FM_RNminus","MK_P_RNminus",
                                  "MK_U_RNminus","MK_FM_RNplus",
                                  "MK_P_RNplus","MK_U_RNplus") ~ "aRef",
                     TRUE ~ "bAlt"))

RefFp <- Refcompare %>% filter(model %in% c("MK_FM_RNplus")) %>%
  mutate(modeltype = "FM+")
RefFm <- Refcompare %>% filter(model %in% c("MK_FM_RNminus")) %>%
  mutate(modeltype = "FM-")


Refcompare <- bind_rows(Refcompare,RefFp,RefFm) %>%
  mutate(modeltype = factor(modeltype, levels = c("F-","F+","FM-","FM+",
                                                  "P-","P+","U-","U+"),
                            labels = c("F2- - F-",
                                       "F2+ - F+",
                                       "F3- - F-",
                                       "F3+ - F+",
                                       "P2- - P-",
                                       "P2+ - P+",
                                       "U2- - U-",
                                       "U2+ - U+")))
Allmodels <- Refcompare %>% group_by(modeltype,exp,cvid) %>%
  arrange(modeltype,reference) %>% summarize(diffs = diff(AIC))

Refcomp <- ggplot(data = Allmodels,
       aes(x = cvid, y = diffs, colour = exp)) +
  scale_y_continuous(limits=c(-50,50),
                     name=expression(paste(Delta, "AIC (alternative model - reference model)"))) +
  geom_point() +
  facet_wrap(.~modeltype, nrow=4) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size=12),
        axis.title.y = element_text(size=12),
        axis.title.x = element_blank(),
        legend.text = element_text(size = 12),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = "white"),
        legend.key = element_rect(fill="transparent"),
        legend.title = element_blank(),
        panel.border = element_rect(colour = "black", fill = "transparent"),
        strip.background = element_rect(fill = "white"),
        strip.text.x = element_text(size=12),
        strip.placement = "outside",
        strip.text = element_text(size=12))

plot_grid(allAIC,Refcomp,nrow=1,rel_widths=c(1,2))
#ggsave("VP_modelvariants.png", units="cm", width=35, height=20, dpi=600)

