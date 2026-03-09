#---------------------------------------------------------------------
# MyDiv experiment; Litterfall project
# 2025-11-12
# directly use NIRS output
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)

# Propagate error from the NIRS-inferred nutrient concentrations:

#============================ Packages ==============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)

#============================ Dataset ===============================

# CNP_predicted <- read_csv
# df.mass.info <- read_csv

#============================ Data sorting ============================

# get kg/kg or percentage out of the mg/g concentrations
CNP_predicted.p <- CNP_predicted %>%
  dplyr::mutate(Cp = C * 0.001,
                Np = N * 0.001,
                Pp = P * 0.001)

CNP_predicted.p <-  CNP_predicted.p %>%
  rename(
    C_pred = C, # mg/g
    N_pred = N,
    P_pred = P
  )

# g/m2 in kg/ha
df.mass.info.kgha <- df.mass.info %>%
  dplyr::mutate(Ac_kgha = Ac_m2 * 10,	
                Ae_kgha = Ae_m2 * 10,
                Be_kgha = Be_m2 * 10,
                Ca_kgha = Ca_m2 * 10,
                Fa_kgha = Fa_m2 * 10,
                Fr_kgha = Fr_m2 * 10,
                Pr_kgha = Pr_m2 * 10,	
                Qu_kgha = Qu_m2 * 10, 
                So_kgha = So_m2 * 10,
                Ti_kgha = Ti_m2 * 10,
                cont_kgha = cont_m2 * 10)

df.mass.info.kgha <- df.mass.info.kgha %>%
  dplyr::select(cols=-c(8:19))


# average across traps ####
df.mass.info.mean <- df.mass.info.kgha %>%
  dplyr::group_by(plotID, plotName, tree_species_richness, mycorrhizal_type, myc, sr, div, block, blk, month1, month, composition) %>% # removed trap and month
  dplyr::summarise(Ac = mean(Ac_kgha, na.rm = TRUE),	
                   Ae = mean(Ae_kgha, na.rm = TRUE),
                   Be = mean(Be_kgha, na.rm = TRUE),	
                   Ca = mean(Ca_kgha, na.rm = TRUE),
                   Fa = mean(Fa_kgha, na.rm = TRUE),	
                   Fr = mean(Fr_kgha, na.rm = TRUE),	
                   Pr = mean(Pr_kgha, na.rm = TRUE),	
                   Qu = mean(Qu_kgha, na.rm = TRUE), 
                   So = mean(So_kgha, na.rm = TRUE),	
                   Ti = mean(Ti_kgha, na.rm = TRUE), 
                   cont = mean(cont_kgha, na.rm = TRUE))%>%
  dplyr::select(col=-c(cont))


# long format ####
df.mass.long <- df.mass.info.mean %>% 
  pivot_longer(cols=c("Ac":"Ti"),
               names_to="species",
               values_to="mass")

df.mass.long <- na.omit(df.mass.long)

df.mass.long <- df.mass.long %>% 
  dplyr::mutate(month2 = recode(month1,
                                January = "Jan",
                                February = "Feb",
                                March = "Mar",
                                April = "Apr",
                                May = "May",
                                June = "Jun",
                                July = "Jul",
                                August = "Aug",
                                September = "Sep",
                                October = "Oct",
                                November = "Nov",
                                December = "Dec")) %>%
  ungroup() %>%
  dplyr::select(-c("month","month1","myc","sr","div","blk"))

df.mass.long$myc <- recode_factor(df.mass.long$mycorrhizal_type, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")

df.mass.long.2 <- df.mass.long %>%
  dplyr::select(-c("mycorrhizal_type"))

CNP_predicted.p$myc <- recode_factor(CNP_predicted.p$mycorrhizal_type, "A" ="AM", "E" ="EM", "AE" = "AM + EM")

CNP_predicted.p.2 <- CNP_predicted.p %>%
  dplyr::select(-c("mycorrhizal_type"))

# join both datasets
df.mass.CNP <- full_join(df.mass.long.2, CNP_predicted.p.2, by=c("plotID","plotName","tree_species_richness","myc","month2"="month","species","composition","block"))


df.mass.CNP$div<-as.factor(df.mass.CNP$tree_species_richness)
df.mass.CNP$blk<-as.factor(df.mass.CNP$block)
df.mass.CNP$sr<-df.mass.CNP$tree_species_richness
df.mass.CNP$block <- as.factor(df.mass.CNP$block)


# sum(concentration_species * biomass_species)/total biomass

# we need to consider the proportion of litter per species on each plot first
# otherwise we don't know the true backflow of nutrients to the ground

# Cflux and Nflux and Pflux per species
df.mass.CNP$Cflux <- df.mass.CNP$C_pred/1000 * df.mass.CNP$mass # C pred in mg/m-2
df.mass.CNP$Nflux <- df.mass.CNP$N_pred/1000 * df.mass.CNP$mass
df.mass.CNP$Pflux <- df.mass.CNP$P_pred/1000 * df.mass.CNP$mass

### annual effects - sum per plot (no month) ####

# Define residual means and SDs (from NIRS model output)
m.c  <- -0.0082
sd.c <- 0.5179

m.n  <- -0.0072
sd.n <- 0.1207

m.p  <- 0
sd.p <- 0.68


run_one_sim_lme <- function(data, iter_id, m.c, sd.c, m.n, sd.n, m.p, sd.p) {
  n <- nrow(data)
  
  # Simulate new observed values using residual means & SDs
  sim_df <- data %>%
    rowwise() %>%
    mutate(
      C_obs = C_pred + rnorm(1, mean = m.c, sd = sd.c), #Cpred= C concentration based on predicted NIRS data
      N_obs = N_pred + rnorm(1, mean = m.n, sd = sd.n),
      P_obs = P_pred + rnorm(1, mean = m.p, sd = sd.p)
    )

  sim_df <-
    sim_df %>%
    group_by(plotID, blk, sr, myc, composition) |>
    summarise(Cflux.tot = sum((C_obs/1000) * mass, na.rm = T),
              Nflux.tot = sum((N_obs/1000) * mass, na.rm = T),
              Pflux.tot = sum((P_obs/1000) * mass, na.rm = T),
              mass.tot = sum(mass, na.rm = T)) %>%
    mutate(Cc=(Cflux.tot/mass.tot)*1000
           ) %>%   # C concentration *1000 so the unit is mg/g
    mutate(Nc=(Nflux.tot/mass.tot)*1000
           ) %>%   # N concentration
    mutate(Pc=(Pflux.tot/mass.tot)*1000
           ) %>%   # P concentration
    mutate(CN=Cc/Nc) %>% # C:N ratio
    mutate(CP=Cc/Pc) %>% # C:P ratio
    mutate(NP=Nc/Pc) %>% # N:P ratio
    mutate(across(c(Cc, Nc, Pc), ~na_if(., 0))) #turn zeros in C, N, and P concentration to NA
  
  
  # Refit paper models
  modC <- lme(Cc ~ log2(sr) * myc,  # change Cc etc
              random = ~1 | composition,
              data = sim_df,
              method = "REML",
              na.action = na.exclude)
  
  modN <- lme(Nc ~ log2(sr) * myc,
              random = ~1 | composition,
              data = sim_df,
              method = "REML",
              na.action = na.exclude)
  
  modP <- lme(Pc ~ log2(sr) * myc,
              random = ~1 | composition,
              data = sim_df,
              method = "REML",
              na.action = na.exclude)
  
  # Extract fixed effects
  extract_fixed <- function(mod, response) {
    coefs <- summary(mod)$tTable
    varcorr <- VarCorr(mod)
    df <- data.frame(
      term = rownames(coefs),
      estimate = coefs[, "Value"],
      std.error = coefs[, "Std.Error"],
      t.value = coefs[, "t-value"],
      p.value = coefs[, "p-value"],
      re_sd = as.numeric(varcorr[1, "StdDev"]),   # random effect SD
      resid_sd = as.numeric(varcorr[nrow(varcorr), "StdDev"]), # residual SD
      response = response,
      iter = iter_id,
      stringsAsFactors = FALSE
    )
    rownames(df) <- NULL
    return(df)
  }
  
  bind_rows(
    extract_fixed(modC, "C"),
    extract_fixed(modN, "N"),
    extract_fixed(modP, "P")
  )
}

n_iter <- 100  

sim_results <- do.call(
  rbind,
  lapply(1:n_iter, function(i) run_one_sim_lme(
    df.mass.CNP, i, 
    m.c = m.c, sd.c = sd.c, 
    m.n = m.n, sd.n = sd.n, 
    m.p = m.p, sd.p = sd.p
  ))
)

run_one_sim <- function(data, iter_id, m.c, sd.c, m.n, sd.n, m.p, sd.p) {
  n <- nrow(data)
  
  # Simulate new observed values using residual means & SDs
  sim_df <- data %>%
    rowwise() %>%
    mutate(
      C_obs = C_pred + rnorm(1, mean = m.c, sd = sd.c), #Cc= C concentration based on predicted NIRS data
      N_obs = N_pred + rnorm(1, mean = m.n, sd = sd.n),
      P_obs = P_pred + rnorm(1, mean = m.p, sd = sd.p)
    )
  
  sim_df <-
    sim_df %>%
    group_by(plotID, blk, sr, myc, composition) |>
    summarise(Cflux.tot = sum((C_obs/1000) * mass, na.rm = T),
              Nflux.tot = sum((N_obs/1000) * mass, na.rm = T),
              Pflux.tot = sum((P_obs/1000) * mass, na.rm = T),
              mass.tot = sum(mass, na.rm = T)) %>%
    mutate(Cc=(Cflux.tot/mass.tot)*1000
           ) %>%   # C concentration *1000 so the unit is mg/g
    mutate(Nc=(Nflux.tot/mass.tot)*1000
           ) %>%   # N concentration
    mutate(Pc=(Pflux.tot/mass.tot)*1000
           ) %>%   # P concentration
    mutate(CN=Cc/Nc) %>% # C:N ratio
    mutate(CP=Cc/Pc) %>% # C:P ratio
    mutate(NP=Nc/Pc) %>% # N:P ratio
    mutate(across(c(Cc, Nc, Pc), ~na_if(., 0))) #turn zeros in C, N, and P concentration to NA
  
  sim_df

}


n_iter <- 100 

sim_data <- map_df(.x=1:n_iter, 
                   .f= ~{
                     run_one_sim(
                       data = df.mass.CNP, 
                       .x, 
                       m.c = m.c, sd.c = sd.c, 
                       m.n = m.n, sd.n = sd.n,
                       m.p = m.p, sd.p = sd.p) |> 
                       mutate(n_iter = .x)})

pl <- ggpubr::ggarrange(
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = Cc)) + 
    geom_violin()+
    labs(y= expression(atop("C concentration", (mg ~ g^-1))),
         x = "plot ID")+
    labs(tag = "(a)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = Nc)) + 
    geom_violin()+
    labs(y= expression(atop("N concentration", (mg ~ g^-1))),
         x = "plot ID")+
    labs(tag = "(b)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = Pc)) + 
    geom_violin()+
    labs(y= expression(atop("P concentration", (mg ~ g^-1))),
         x = "plot ID")+
    labs(tag = "(c)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = Cflux.tot)) + 
    geom_violin()+
    labs(y=expression(atop("C flux", (kg ~ ha^-1 ~ yr^-1))),
         x = "plot ID")+
    labs(tag = "(d)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = Nflux.tot)) + 
    geom_violin()+
    labs(y=expression(atop("N flux", (kg ~ ha^-1 ~ yr^-1))),
         x = "plot ID")+
    labs(tag = "(e)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = Pflux.tot)) + 
    geom_violin()+
    labs(y=expression(atop("P flux", (kg ~ ha^-1 ~ yr^-1))),
         x = "plot ID")+
    labs(tag = "(f)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = CN)) + 
    geom_violin()+
    labs(y= "C : N",
         x = "plot ID")+
    labs(tag = "(g)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = CP)) + 
    geom_violin()+
    labs(y= "C : P",
         x = "plot ID")+
    labs(tag = "(h)"),
  ggplot(data = sim_data, 
         aes(x = factor(plotID), y = NP)) + 
    geom_violin()+
    labs(y= "N : P",
         x = "plot ID")+
    labs(tag = "(i)"),
  
  ncol = 3, nrow = 3,
  align ="hv", common.legend = T, legend="bottom"
)

pl


p.C = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=Cc, color=myc))
p.N = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=Nc, color=myc))
p.P = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=Pc, color=myc))
p.Cf = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=Cflux.tot, color=myc))
p.Nf = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=Nflux.tot, color=myc))
p.Pf = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=Pflux.tot, color=myc))
p.CN = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=CN, color=myc))
p.CP = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=CP, color=myc))
p.NP = 
  ggplot(data=sim_data, 
         aes(x=sr, 
             y=NP, color=myc))

# # remove P concentration outlier
# sim_data <- sim_data[-1062, ]

for(i in 1:100){
  p.C = 
    p.C +
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=Cc, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=Cc, 
                    color=myc),
                se=F, size = 0.5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y= expression(atop("C concentration", (mg ~ g^-1))),
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
  
  p.N = 
    p.N +
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=Nc, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=Nc, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y= expression(atop("N concentration", (mg ~ g^-1))),
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(b)")
  
  p.P = 
    p.P +
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=Pc, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=Pc, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y= expression(atop("P concentration", (mg ~ g^-1))),
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
  
  p.Cf = 
    p.Cf +
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=Cflux.tot, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=Cflux.tot, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y=expression(atop("C flux", (kg ~ ha^-1 ~ yr^-1))),
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(d)")
  
  p.Nf = 
    p.Nf +
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=Nflux.tot, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=Nflux.tot, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y=expression(atop("N flux", (kg ~ ha^-1 ~ yr^-1))),
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(e)")
  
  p.Pf = 
    p.Pf +
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=Pflux.tot, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=Pflux.tot, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y=expression(atop("P flux", (kg ~ ha^-1 ~ yr^-1))),
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(f)")
  
  p.CN = 
    p.CN+
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=CN, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=CN, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y = "C : N",
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(g)")
  
  p.CP = 
    p.CP+
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=CP, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=CP, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y = "C : P",
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(h)")
  
  p.NP = 
    p.NP+
    geom_jitter(data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, 
                    y=NP, color=myc),
                shape =21, size = 1, alpha=0.5, width= 0.1)+
    geom_smooth(method="lm",
                data=sim_data %>% 
                  filter(n_iter==i), 
                aes(x=sr, y=NP, 
                    color=myc),
                se=F, size = .5)+
    scale_x_continuous(trans='log2',
                       breaks=c(1,2,4))+
    labs(y = "N : P",
         x = "Tree species richness")+
    scale_fill_manual(values= c("#71b540","#4c8ecb","#febf00"),
                      name = "Mycorrhizal type",
                      guide="none")+ 
    scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"), 
                       name = "Mycorrhizal type",
                       guide="none")+
    theme(strip.background = element_blank(),
          strip.text = element_text(size=12),
          axis.line = element_line(color='black'),
          axis.text.y = element_text(color="black", size = 12),
          axis.text.x = element_text(color="black", size = 12),
          axis.title.y = element_text(color="black", size = 12),
          axis.title.x = element_blank(),
          axis.ticks = element_line(color="black"),
          strip.text.x = element_text(12),
          panel.border = element_rect(colour="black", fill=NA),
          panel.background = element_rect(fill="white", colour=NA),
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
    labs(tag = "(i)")
  
}


library(ggpubr)

annual.plot <- ggarrange(p.C, p.N, p.P,
                         p.Cf, p.Nf, p.Pf,
                         p.CN, p.CP, p.NP,
                         ncol = 3, nrow = 3, align ="hv", common.legend = T, legend="bottom")

annual.plot


### end ###



