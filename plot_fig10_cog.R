library(tidyverse)
library(ggplot2)
library("RColorBrewer")
library(wesanderson)
library("PerformanceAnalytics")
library(GGally)
library(ggcorrplot)
library(reshape2)
library(gridExtra)
library(cowplot)

get_legend<-function(myggplot){
  tmp <- ggplot_gtable(ggplot_build(myggplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)
}


setwd('~/Box Sync/Connectome/')

cog <- read.csv("reference/CognitiveScore_2019-08-06_16-54-17.csv", fileEncoding = "UTF-8-BOM")
cog %>% summary()

cog$Visit[cog$Visit.Id == 0] = 'BSL'
cog$Visit[cog$Visit.Id == 1] = 'Wk12'
cog$Visit[cog$Visit.Id == 2] = 'Yr1'
cog$Visit[cog$Visit.Id == 3] = 'Yr2'
cog$cohort <- as.factor(cog$cohort)
cog$Visit <- as.factor(cog$Visit)
cog %>% summary()


# read age_matched subjects, who is in image analysis
age_match_list <- read.csv('../Connectome/reference/age_matched_list.csv', fileEncoding = "UTF-8-BOM")
cog_age_matched <- inner_join(cog,age_match_list, by = c("Participant.Id" = "subject.id"))
cog_age_matched$ParticipantId <- as.factor(cog_age_matched$Participant.Id)
cog_age_matched %>% 
  select (ParticipantId, site, cohort.x, EXE.ZS, SPE.ZS, ATT.ZS, LEA.ZS, MEM.ZS, MOT.ZS, OVERALL.ZS, age, Visit) %>% 
  summary()


# compare hiv bsl wk12, and hc  
cog_hiv_bsl_wk12_hc <- cog %>% 
  filter( Visit == 'BSL' | Visit == 'Wk12') %>%
  filter(cohort != 'LTNP') %>% 
  filter(site == 'ROC' ) %>%
  group_by(cohort)


cog_hiv_bsl_wk12_hc$Visit
cog_hiv_bsl_wk12_hc$cohort

cog_hiv_bsl_wk12_hc <- unite(cog_hiv_bsl_wk12_hc, "group", c("cohort","Visit"))
cog_hiv_bsl_wk12_hc$group <- as.factor(cog_hiv_bsl_wk12_hc$group)



p1 <- ggplot(data = cog_hiv_bsl_wk12_hc, aes(x = group, y = OVERALL.ZS, fill = group))
p1 <- p1 + scale_color_manual(values = brewer.pal(n = 8, name = "Set1")) + 
  scale_fill_brewer(palette = "Set1") + 
  geom_jitter(width = .2)+
  geom_boxplot(alpha = .8) +
  ylab("Overall Z-Score") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 15)) +
  scale_x_discrete(labels=c("HC_BSL" = "HC", "HIV+_BSL" = "HIV_BSL",
                            "HIV+_Wk12" = "HIV_12wk"))

p1

p3 <- ggplot(data = cog_hiv_bsl_wk12_hc, aes(x = group, y = MOT.ZS, fill = group))
p3 <- p3 + scale_color_manual(values = brewer.pal(n = 8, name = "Set1")) + 
  scale_fill_brewer(palette = "Set1") + 
  geom_jitter(width = .2)+
  geom_boxplot(alpha = .8) +
  ylab("Motor Z-Score") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        axis.title.x=element_blank(),
        legend.text=element_text(size=24),
        axis.title.y = element_text(size = 24),
        axis.text = element_text(size = 15)) +
  scale_x_discrete(labels=c("HC_BSL" = "HC", "HIV+_BSL" = "HIV_BSL",
                            "HIV+_Wk12" = "HIV_12wk"))


p3


ggsave("Overall_Z-score.pdf",p1, width = 8, height = 6)
ggsave("Motor_Z_score.pdf",p3, width = 8, height = 6)

