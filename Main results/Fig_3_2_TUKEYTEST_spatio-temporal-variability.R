#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-02-27
# Tukey test: Fig. 3 spatio-temporal variability of litter fall
# Coefficient of variation
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)
# anova Type I

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

# 1) NO! average across traps ####
df.all.wide.trap <- df.all.wide.info 

df.all.wide.trap$div<-as.factor(df.all.wide.trap$tree_species_richness)
df.all.wide.trap$blk<-as.factor(df.all.wide.trap$block)
df.all.wide.trap$myc<-as.factor(df.all.wide.trap$mycorrhizal_type)
df.all.wide.trap$sr<-df.all.wide.trap$tree_species_richness
df.all.wide.trap$sr_myc<-paste(df.all.wide.trap$sr,df.all.wide.trap$myc,sep="_")#interaction terms
df.all.wide.trap$myc <- recode_factor(df.all.wide.trap$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")
df.all.wide.trap$block <- as.factor(df.all.wide.trap$block)


# 2) long format ####
df.all.long <- df.all.wide.trap %>% 
  pivot_longer(cols=c(20:29),
               names_to="species",
               values_to="litterfall")


#### 3.1) Yearly spatial variability (Coefficient of variation CV) ####
df.year.spat.variability = df.all.long |>
  group_by(block, sr, div, myc, plotID, composition, month1) |> 
  summarise(litterfall_CV = sd(litterfall, na.rm=T)/mean(litterfall, na.rm=T))|>  # mean between traps and sd between traps
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, plotID, composition) |> 
  summarise(litterfall_CV = mean(litterfall_CV, na.rm=T))

# distribution of data
hist(log(df.year.spat.variability$litterfall_CV))

# 3.2) Plot: yearly spatial variability  ####
spa.var.yearly <- ggplot()+
  geom_point(data = df.year.spat.variability , 
             aes(x=sr, y=litterfall_CV, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.spat.variability, 
              aes(x=sr, y=litterfall_CV, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Annual spatial variability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  #scale_y_log10() +
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
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
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(c)")
show(spa.var.yearly)

#3.3) Model ####
#library(nlme)
mod.year.spa.variability =
  lme(litterfall_CV ~ log2(sr) * myc,
      random = ~1|composition,
      data = df.year.spat.variability)

# 3.4) Check the model quality ####
#library(performance)
#performance::check_model(mod.year.spa.variability)

# 3.5) Summary ####
summary(mod.year.spa.variability)

# 3.6) Anova (Type III SS) ####
anova(object = mod.year.spa.variability, type = 'sequential')

library(emmeans)
emmeans(mod.year.spa.variability, list(pairwise ~ "myc"), adjust = "tukey")


# Get average CV per richness level and mycorrhizal type
mean_CV <- df.year.spat.variability %>%
  group_by(sr, myc) %>%
  summarise(meanCV = mean(litterfall_CV, na.rm = TRUE), .groups = "drop")

spat_wide <- mean_CV %>%
  pivot_wider(names_from = myc, values_from = meanCV)

spat_wide <- spat_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per mycorrhizal type
mean_CV_myc <- df.year.spat.variability %>%
  group_by(myc) %>%
  summarise(meanCV = mean(litterfall_CV, na.rm = TRUE), .groups = "drop")

mean_CV_myc_wide <- mean_CV_myc %>%
  pivot_wider(names_from = myc, values_from = meanCV)

mean_CV_myc_wide <- mean_CV_myc_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per tree species richness
mean_CV_sr <- df.year.spat.variability %>%
  group_by(sr) %>%
  summarise(meanCV = mean(litterfall_CV, na.rm = TRUE), .groups = "drop")

mean_CV_sr_wide <- mean_CV_sr %>%
  pivot_wider(names_from = sr, values_from = meanCV)

mean_CV_sr_wide <- mean_CV_sr_wide %>%
  mutate(
    `1_vs_2` = (`1` - `2`) / `2` * 100,
    `1_vs_4` = (`1` - `4`) / `4` * 100,
    `2_vs_4` = (`2` - `4`) / `4` * 100)








#### 4.1) Monthly spatial variabilty (CV coefficient of variation) ####
df.month.spat.variability = df.all.long |>
  group_by(block, sr, div, myc, plotID, composition, month1) |> 
  summarise(litterfall_CV = sd(litterfall, na.rm=T)/mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]),labels=month.abb))

df.month.spat.variability$month2 <- dplyr::recode_factor(df.month.spat.variability$month1, 
                                             "March"="Mar", 
                                             "April"="Apr",
                                             "May"="May",
                                             "June"="Jun",
                                             "July"="Jul",
                                             "August"="Aug",
                                             "September"="Sep",
                                             "October"="Oct",
                                             "November"="Nov",
                                             "December"="Dec",
                                             "January"="Jan",
                                             "February"="Feb")

df.month.spat.variability$month1 = factor(df.month.spat.variability$month2, 
                              levels = c(month.abb[3:12], month.abb[1:2]))

# 4.2) Check months individually ####
M = map_df( .x = unique(df.month.spat.variability$month1),
            .f = ~ {
              mod = 
                lme(log(litterfall_CV) ~ log2(sr) * myc, 
                    random= ~1|composition,
                    data= df.month.spat.variability |>
                      filter(month1 == .x), 
                    na.action = na.omit) |> 
                anova(type = 'sequential') |> 
                data.frame() |> 
                mutate(month = .x)
              mod$exp = rownames(mod)
              rownames(mod) <- NULL
              mod |> 
                select(month, explanatory = exp, 
                       numDF, denDF, F.value, p.value)
            }) |> 
  mutate(sign = if_else(p.value < 0.001, '***', 
                        if_else(p.value < 0.01, "**", 
                                if_else(p.value < 0.05, '*',
                                        if_else(p.value<0.1, '.', 'n.s.')))))
M

M$sign[M$explanatory == 'log2(sr)' & M$month == 'March'] = 
  paste0("log2(sr) = ",M$sign[M$explanatory == 'log2(sr)' & M$month == 'March'])
M$sign[M$explanatory == 'myc' & M$month == 'March'] = 
  paste0("myc = ",M$sign[M$explanatory == 'myc' & M$month == 'March'])
M$sign[M$explanatory == 'log2(sr):myc' & M$month == 'March'] = 
  paste0("log2(sr):myc = ",M$sign[M$explanatory == 'log2(sr):myc' & M$month == 'March'])


# 4.3) Plot: Monthly spatial variability ####
spa.var.monthly <- ggplot()+
  geom_point(data = df.month.spat.variability, 
             aes(x=sr, y=litterfall_CV, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.month.spat.variability, 
              aes(x=sr, y=litterfall_CV, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  facet_grid(.~month1)+
  labs(y=bquote("Monthly spatial variability"), 
       x = "Tree species richness")+
  scale_y_log10(breaks=c(0.0001,0.001,0.01,0.1,1,10,100),labels=c(0.0001,0.001,0.01,0.1,1,10,100))+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type")+
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr)'),
            aes(x=1, y=20,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'myc'),
            aes(x=1, y=15,
                label = sign,
                hjust = 0),
            color = 'black') +
  geom_text(data = M |>
              mutate(month1 = month) |>
              filter(explanatory == 'log2(sr):myc'),
            aes(x=1, y=10,
                label = sign,
                hjust = 0),
            color = 'black') +
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
        legend.position = "bottom",
        legend.direction = "horizontal",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(d)")

spa.var.monthly

# 4.4) Model ####
# with correlation structure
#library(nlme)

hist(df.month.spat.variability$litterfall_CV)

mod.monthly.spat.variability =
  lme(log(litterfall_CV) ~ month1 * log2(sr) * myc,
      random = ~1|composition,
      data = df.month.spat.variability,
      correlation=corCAR1(), na.action = na.omit)

# 4.5) Check the model quality ####
# library(performance)
performance::check_model(mod.monthly.spat.variability)

# 4.6) Summary ####
summary(mod.monthly.spat.variability)

# 4.7) Anova (Type I SS) ####
anova(mod.monthly.spat.variability, type="sequential")


#### TUKEY TEST MONTHLY SPATIAL VARIABILITY LITTERFALL ####

library(nlme)
library(emmeans)
library(dplyr)
library(purrr)
library(writexl)

# Calculate raw group means per month and myc
group_means <- df.month.spat.variability %>%
  group_by(month1, myc) %>%
  summarise(mean_litter = mean(litterfall_CV, na.rm = TRUE), .groups = "drop") %>%
  rename(month = month1)

library(stringr)

clean_contrast <- function(x) {
  x %>%
    str_replace_all("[()]", "") %>%   # remove parentheses
    str_squish()                      # remove extra whitespace
}

# Get all pairwise combinations of myc groups per month, calculate relative % difference
pairwise_raw_diff <- group_means %>%
  group_by(month) %>%
  summarise(
    comb = list(combn(myc, 2, simplify = FALSE))
  ) %>%
  unnest(cols = c(comb)) %>%
  rowwise() %>%
  mutate(
    group1 = comb[[1]],
    group2 = comb[[2]],
    mean1 = group_means$mean_litter[group_means$month == month & group_means$myc == group1],
    mean2 = group_means$mean_litter[group_means$month == month & group_means$myc == group2],
    rel_change_percent = round(100 * (mean1 - mean2) / mean2, 1),
    contrast = paste(group1, "-", group2)
  ) %>%
  select(month, contrast, rel_change_percent) %>%
  ungroup()

myc_sig_months <- M %>%
  filter(explanatory == "myc", p.value < 0.05) %>%
  pull(month)

Tukey_results_monthly_lit <- map_df(myc_sig_months, function(mo) {
  dat <- df.month.spat.variability %>% filter(month1 == mo)
  
  mod <- lme(litterfall_CV ~ log2(sr) * myc, random = ~1|composition, na.action = na.omit, data = dat)
  
  out <- tryCatch({
    emmeans(mod, pairwise ~ myc, adjust = "tukey")$contrasts %>%
      as.data.frame() %>%
      mutate(month = mo) %>%
      select(month, contrast, df, t.ratio, p.value)
    
  }, error = function(e) {
    message("Skipping month: ", mo, " due to error: ", e$message)
    NULL
  })
  
  out
})

Tukey_results_monthly_lit <- Tukey_results_monthly_lit %>%
  mutate(contrast = clean_contrast(contrast))

pairwise_raw_diff <- pairwise_raw_diff %>%
  mutate(contrast = clean_contrast(contrast))


final_results3 <- left_join(pairwise_raw_diff, Tukey_results_monthly_lit, by = c("month", "contrast"))



#### 5.1) Yearly temporal variability  ####
df.year.temp.variability = df.all.long |>
  group_by(block, sr, div, myc, plotID, composition, month1) |> 
  summarise(litterfall_mean = mean(litterfall, na.rm=T))|>
  mutate(month1 = factor(month1, levels = c(month.name[3:12], month.name[1:2]))) |>
  group_by(block, sr, div, myc, composition, plotID) |> 
  summarise(litterfall_CV = sd(litterfall_mean, na.rm=T)/mean(litterfall_mean, na.rm=T))

hist(log(df.year.temp.variability$litterfall_CV))

# 5.2) Plot: Yearly temporal variability ####
temp.var.yearly <- ggplot()+
  geom_point(data = df.year.temp.variability , 
             aes(x=sr, y=litterfall_CV, 
                 color = myc, fill = myc), shape =21, size = 1, alpha=0.5)+
  geom_smooth(data = df.year.temp.variability, 
              aes(x=sr, y=litterfall_CV, 
                  color = myc, fill = myc),
              method="lm", alpha=0.3)+
  labs(y=bquote("Annual temporal variability"), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
  #scale_y_log10()+
  scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                    name = "Mycorrhizal type",
                    guide="none")+ 
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                     name = "Mycorrhizal type",
                     guide="none")+
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
        legend.box.background = element_rect(color="transparent"))+
  labs(tag = "(a)")
temp.var.yearly

# 5.3) Model ####
#library(nlme)
mod.yearly.temp.var =
  lme(litterfall_CV ~ log2(sr) * myc ,
      random = ~1|composition,
      data = df.year.temp.variability)

# 5.4) Check the model quality ####
#library(performance)
#performance::check_model(mod.yearly.temp.var)

# 5.5) Summary ####
summary(mod.yearly.temp.var)

# 5.6) Anova (Type I SS) ####
anova(mod.yearly.temp.var, type="sequential")



emmeans(mod.yearly.temp.var, list(pairwise ~ "myc"), adjust = "tukey")

# Get average CV per richness level and mycorrhizal type
mean_CV <- df.year.temp.variability %>%
  group_by(sr, myc) %>%
  summarise(meanCV = mean(litterfall_CV, na.rm = TRUE), .groups = "drop")

spat_wide <- mean_CV %>%
  pivot_wider(names_from = myc, values_from = meanCV)

spat_wide <- spat_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per mycorrhizal type
mean_CV_myc <- df.year.temp.variability %>%
  group_by(myc) %>%
  summarise(meanCV = mean(litterfall_CV, na.rm = TRUE), .groups = "drop")

mean_CV_myc_wide <- mean_CV_myc %>%
  pivot_wider(names_from = myc, values_from = meanCV)

mean_CV_myc_wide <- mean_CV_myc_wide %>%
  mutate(
    EM_vs_AM      = (EM - AM) / AM * 100,
    AMEM_vs_AM    = (`AM + EM` - AM) / AM * 100,
    AMEM_vs_EM    = (`AM + EM` - EM) / EM * 100
  )

# Get average CV per tree species richness
mean_CV_sr <- df.year.temp.variability %>%
  group_by(sr) %>%
  summarise(meanCV = mean(litterfall_CV, na.rm = TRUE), .groups = "drop")

mean_CV_sr_wide <- mean_CV_sr %>%
  pivot_wider(names_from = sr, values_from = meanCV)

mean_CV_sr_wide <- mean_CV_sr_wide %>%
  mutate(
    `1_vs_2` = (`1` - `2`) / `2` * 100,
    `1_vs_4` = (`1` - `4`) / `4` * 100,
    `2_vs_4` = (`2` - `4`) / `4` * 100)


### end ###
