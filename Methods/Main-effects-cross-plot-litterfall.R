#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-09-04
# cross-plot litterfall (i.e., leaf litter originating from neighboring plots and non-target species)  in g per m2
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)
# anova Type  I

#============================ Packages ===============================

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

#============================ Dataset ===============================

# df.all.wide.info <- read_csv

# 1) average across traps ####
df.all.wide.mean <- df.all.wide.info %>%
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, block, month, composition) %>% # removed trap and month
  dplyr::summarise(Ac_mean = mean(Ac, na.rm = TRUE),	
                   Ae_mean = mean(Ae, na.rm = TRUE),
                   Be_mean = mean(Be, na.rm = TRUE),	
                   Ca_mean = mean(Ca, na.rm = TRUE),
                   Fa_mean = mean(Fa, na.rm = TRUE),	
                   Fr_mean = mean(Fr, na.rm = TRUE),	
                   Pr_mean = mean(Pr, na.rm = TRUE),	
                   Qu_mean = mean(Qu, na.rm = TRUE), 
                   So_mean = mean(So, na.rm = TRUE),	
                   Ti_mean = mean(Ti, na.rm = TRUE), 
                   cont_mean = mean(cont, na.rm = TRUE))

df.all.wide.mean$div<-as.factor(df.all.wide.mean$tree_species_richness)
df.all.wide.mean$blk<-as.factor(df.all.wide.mean$block)
df.all.wide.mean$myc<-as.factor(df.all.wide.mean$mycorrhizal_type)
df.all.wide.mean$sr<-df.all.wide.mean$tree_species_richness
df.all.wide.mean$sr_myc<-paste(df.all.wide.mean$sr,df.all.wide.mean$myc,sep="_")#interaction terms
df.all.wide.mean$myc <- recode_factor(df.all.wide.mean$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")


# 2) long format ####
df.all.long <- df.all.wide.mean %>% 
  pivot_longer(cols=c("Ac_mean":"cont_mean"),
               names_to="species",
               values_to="dryweight")


# 3) plot main interaction effects ####
df.annual.litter = 
  df.all.long |> 
  group_by(plotID, myc, tree_species_richness, block, composition, species) |>
  summarise(litter.prod = sum(dryweight, na.rm = T), 
            sd.litter.prod = sd(dryweight, na.rm = T))

species_names <- c(
  "Ac_mean" = "Acer pseudoplatanus",	
  "Ae_mean" = "Aesculus hippocastanum",
  "Be_mean" = "Betula pendula",	
  "Ca_mean" = "Carpinus betulus",
  "Fa_mean" = "Fagus sylvatica",	
  "Fr_mean" = "Fraxinus excelsior",	
  "Pr_mean" = "Prunus avium",	
  "Qu_mean" = "Quercus petraea", 
  "So_mean" = "Sorbus aucuparia",	
  "Ti_mean" = "Tilia platyphyllos", 
  "cont_mean" = "contaminants"
)

df.annual.litter$species2 <- recode(df.annual.litter$species, !!!species_names)

lit <- 
  ggplot(df.annual.litter, aes(x = species2, y = litter.prod, fill = species)) +
  geom_col() +
  labs(y=bquote("Leaf litter dryweight"~(g~m^-2)), 
       x = "Tree species")+
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size=12),
        axis.line = element_line(color='black'),
        axis.text.y = element_text(color="black", size = 12),
        axis.text.x = element_text(color="black", angle = 45, hjust = 1, face="italic"),
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
        legend.position = "none",
        legend.key = element_rect(color="transparent"),   
        legend.title = element_text("Biodiversity effects", size = 12),
        legend.text = element_text(size=12),
        legend.background = element_rect(colour=NA),
        legend.box= NULL,
        legend.box.background = element_rect(color="transparent"))
lit


df.crossplot <- df.all.wide.mean %>%
  rowwise() %>%
  mutate(total_litter = sum(c_across(ends_with("_mean")), na.rm = TRUE),
         frac_crossplot = cont_mean / total_litter) %>%
  ungroup()

df.crossplot_summary <- df.crossplot %>%
  group_by(sr, myc) %>%
  summarise(mean_frac = mean(frac_crossplot, na.rm = TRUE),
            se_frac = sd(frac_crossplot, na.rm = TRUE) / sqrt(n()))

cross <- 
  ggplot(df.crossplot_summary, aes(x = sr, y = mean_frac, color = myc)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = mean_frac - se_frac, ymax = mean_frac + se_frac),
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(y=bquote("Cross-plot litterfall (contaminants)"~~"% of total litterfall"), ### in g per m2!!!
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
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
cross

cross <- 
  ggplot(df.crossplot_summary, aes(x = sr, y = mean_frac, color = myc)) +
  geom_point(position = position_dodge(width = 0.3), size = 3) +
  geom_errorbar(aes(ymin = mean_frac - se_frac, ymax = mean_frac + se_frac),
                width = 0.2, position = position_dodge(width = 0.3)) +
  labs(
    y = "Cross-plot litterfall (contaminants) \n % of total litterfall",
    x = "Tree species richness") +
  scale_x_continuous(
    trans = 'log2',
    breaks = c(1, 2, 4)
  ) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  scale_fill_manual(
    values = c("#71b540", "#4c8ecb", "#febf00"),
    name = "Mycorrhizal type",
    guide = "none"
  ) +
  scale_color_manual(
    values = c("#71b540", "#4c8ecb", "#febf00"), 
    name = "Mycorrhizal type"
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "sans"),  # default safe font (fixes warnings)
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.line = element_line(color = 'black'),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.ticks.x = element_line(),
    strip.text.x = element_text(size = 12),
    panel.border = element_rect(colour = "black", fill = NA),
    plot.background = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    plot.title = element_text(size = 12),
    plot.subtitle = element_text(size = 12),
    legend.position = "right",
    legend.direction = "vertical",
    legend.key = element_rect(color = "transparent"),   
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_rect(colour = NA),
    legend.box = NULL,
    legend.box.background = element_rect(color = "transparent")
  )

cross



df.crossplot_plotlevel <- df.crossplot %>%
  group_by(plotID, plotName, sr, myc, block, composition) %>%
  summarise(frac_crossplot_mean = mean(frac_crossplot, na.rm = TRUE),
            .groups = "drop")

df.crossplot_plotlevel_month <- df.crossplot %>%
  group_by(plotID, plotName, sr, myc, block, composition, month) %>%
  summarise(frac_crossplot_mean = mean(frac_crossplot, na.rm = TRUE),
            .groups = "drop")


main.eff <- ggplot(df.annual.litter, 
       aes(x=tree_species_richness, y=litter.prod, 
           color = myc, fill = myc))+
  geom_point(shape =21, size = 1, alpha=0.5)+
  geom_smooth(method="lm", alpha=0.3)+
  labs(y=bquote("Leaf litter dryweight of contaminants"~(g~m^-2)), 
       x = "Tree species richness")+
  scale_x_continuous(trans='log2',
                     breaks=c(1,2,4))+
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
main.eff


# 4) Model ####


# new model with lme, composition as random factor and Type 1 anova
mod.total.litterfall2 =
  lme(litter.prod ~ log2(tree_species_richness) * myc,
      random = ~1|composition,
      data = df.annual.litter)

# 5) Check the model quality ####
#library(performance)
performance::check_model(mod.total.litterfall2)
#dev.off()

# 6) Summary ####
summary(mod.total.litterfall2)

# 7) Anova (Type I SS) ####
anova(object = mod.total.litterfall2, type = 'sequential')

# 8) plot boxplots
main.eff.box <- ggplot(df.annual.litter, 
                   aes(x=as.factor(tree_species_richness), y=litter.prod, 
                       color = myc, fill = myc))+
  geom_boxplot(shape =21, size = 1, alpha=0.5)+
  labs(y=bquote("Leaf litter dryweight"~(g~m^-2)), 
       x = "Tree species richness")+
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
main.eff.box


### end ###