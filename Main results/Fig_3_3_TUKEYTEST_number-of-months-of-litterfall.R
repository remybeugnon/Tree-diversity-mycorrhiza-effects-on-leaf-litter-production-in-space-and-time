#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-03-03
# Tukey test: Fig. 3 Number of months of litterfall 
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)
# anova Type  I

#============================ Packages ===============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)

#============================ Dataset ===============================

df.all.wide.info <- read.csv

# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition) %>% # removed trap and month
  dplyr::summarise(Ac_mean = mean(Ac_m2, na.rm = TRUE),	
                   Ae_mean = mean(Ae_m2, na.rm = TRUE),
                   Be_mean = mean(Be_m2, na.rm = TRUE),	
                   Ca_mean = mean(Ca_m2, na.rm = TRUE),
                   Fa_mean = mean(Fa_m2, na.rm = TRUE),	
                   Fr_mean = mean(Fr_m2, na.rm = TRUE),	
                   Pr_mean = mean(Pr_m2, na.rm = TRUE),	
                   Qu_mean = mean(Qu_m2, na.rm = TRUE), 
                   So_mean = mean(So_m2, na.rm = TRUE),	
                   Ti_mean = mean(Ti_m2, na.rm = TRUE), 
                   cont_mean = mean(cont_m2, na.rm = TRUE))

df.all.wide.mean$div<-as.factor(df.all.wide.mean$tree_species_richness)
df.all.wide.mean$blk<-as.factor(df.all.wide.mean$block)
df.all.wide.mean$myc<-as.factor(df.all.wide.mean$mycorrhizal_type)
df.all.wide.mean$sr<-df.all.wide.mean$tree_species_richness
df.all.wide.mean$sr_myc<-paste(df.all.wide.mean$sr,df.all.wide.mean$myc,sep="_")#interaction terms
df.all.wide.mean$myc <- recode_factor(df.all.wide.mean$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")


# 2) long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c(13:22),
               names_to="species",
               values_to="dryweight")


# 3)sum across species per plot
df.litter.sum = df.all.long |> 
  group_by(block, plotID, div, sr, myc, composition, month1) |> 
  summarise(litterfall_sum = sum(dryweight, na.rm = T),
            litterfall_mean = mean(dryweight, na.rm =T),
            litterfall_sd = sd(dryweight, na.rm = T))

# 4) monthly coverage - sum per plot, if litter dryweight larger then 0 ####
df.litter.cover = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall_sum >0) |>
  summarise(number_month_litterfall = n()) 

df_avg <- df.litter.cover %>%
  group_by(sr, myc) %>%
  summarise(avg_litterfall_months = mean(number_month_litterfall, na.rm = TRUE))


ggplot(df.litter.cover, aes(x=number_month_litterfall, y=div, color=myc, fill=myc))+
  geom_bar(position = "dodge", stat = "summary", aes(fill=myc, color=myc), alpha=0.5)+
  labs(y=bquote("Tree species richness"), 
       x = "Number of months of litterfall")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  scale_x_continuous(breaks= c(0:12))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))

temp.cov <- ggplot(df.litter.cover, aes(y=number_month_litterfall, x=sr, color=myc, fill=myc))+
  geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Number of months of litterfall"), 
       x = "Tree species richness")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  scale_y_continuous(breaks= c(0:12))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))

temp.cov

# 5) Model ####
#library(nlme)

df.litter.cover = df.litter.sum |>
  group_by(block, plotID, div, sr, myc, composition) |>
  filter(litterfall_sum >0) |>
  summarise(number_month_litterfall = n())

hist(df.litter.cover$number_month_litterfall)

# new model (Feb 2025)
mod.coverage =
  lme(number_month_litterfall ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.litter.cover)


# 6) Check the model quality ####
#library(performance)
performance::check_model(mod.coverage)

# 7) Summary ####
summary(mod.coverage)

# 8) Anova (Type III SS) ####
anova(object = mod.coverage, type = 'sequential')


# 9) Tukey HSD to test what mycorrhizal types differ 
library(emmeans)
emmeans(mod.coverage, list(pairwise ~ "myc"), adjust = "tukey")


# Get average CV per richness level and mycorrhizal type
d.mean_nr <- df.litter.cover %>%
  group_by(sr, myc) %>%
  summarise(mean_nr = mean(number_month_litterfall, na.rm = TRUE), .groups = "drop")

d.mean_nr_wide <- d.mean_nr %>%
  pivot_wider(names_from = myc, values_from = mean_nr)

d.mean_nr_wide <- d.mean_nr_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per mycorrhizal type
d.mean_nr_myc <- df.litter.cover %>%
  group_by(myc) %>%
  summarise(mean_nr = mean(number_month_litterfall, na.rm = TRUE), .groups = "drop")

d.mean_nr_myc_wide <- d.mean_nr_myc %>%
  pivot_wider(names_from = myc, values_from = mean_nr)

d.mean_nr_myc_wide <- d.mean_nr_myc_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per tree species richness
d.mean_nr_sr <- df.litter.cover %>%
  group_by(sr) %>%
  summarise(mean_nr = mean(number_month_litterfall, na.rm = TRUE), .groups = "drop")

d.mean_nr_sr_wide <- d.mean_nr_sr %>%
  pivot_wider(names_from = sr, values_from = mean_nr)

d.mean_nr_sr_wide <- d.mean_nr_sr_wide %>%
  mutate(
    `1_vs_2` = (`1` - `2`) / `2` * 100,
    `1_vs_4` = (`1` - `4`) / `4` * 100,
    `2_vs_4` = (`2` - `4`) / `4` * 100)


# whats a good threshhold - 0g, 1g, 5g, 10g
ggplot(df.litter.sum, aes(x=div, y=litterfall))+
  geom_boxplot(alpha=0.5)+
  facet_grid(div~myc)+
  stat_summary(fun=mean, geom="point", shape=20, size=3, color="red", fill="red") +
  labs(y=bquote("Tree species richness"), 
       x = "Number of months of litterfall")+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", size = 12),
        axis.title.y = element_text(size = 12),
        axis.title.x = element_text(size=12),
        axis.ticks.x = element_line(),
        strip.text.x = element_text(12),
        panel.border = element_rect(colour="black", fill=NA),
        plot.background = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.title = element_text(size =12),
        plot.subtitle = element_text(size=12),
        legend.position = "right",
        legend.direction = "vertical",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))


### end ###
