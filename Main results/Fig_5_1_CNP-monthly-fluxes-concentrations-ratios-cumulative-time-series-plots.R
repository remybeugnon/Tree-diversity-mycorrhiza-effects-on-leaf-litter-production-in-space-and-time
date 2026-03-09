#---------------------------------------------------------------------
# MyDiv experiment; Litterfall data March 2023 - February 2024
# 2025-06-25
# annual and monthly effects (tree diversity and mycorrhizal type) on leaf litter C, N, P concentration, flux and CN, CP, NP
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)
# anova Type  I

#============================ Packages ===============================

install.packages(c(
  "readr",
  "readxl",
  "tidyverse",
  "nlme",
  "lme4",
  "lmerTest",
  "ggpubr",
  "patchwork"
))

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(ggpubr)
library(patchwork)  # for combining plots

#============================ Datasets ===============================

# df.CNP.info <- read_csv
# df.mass.info <- read_csv
#============================ Data sorting ============================

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


# (1) average across traps ####
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

# exclude rows with NAs/missing data
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

df.CNP.info <- df.CNP.info |>
    dplyr::select(-c("block","plotName","tree_species_richness","mycorrhizal_type")) %>%
    dplyr::rename(month2=month)

# get kg/kg or percentage out of the mg/g concentrations

df.CNP.info.percent <- df.CNP.info %>%
  dplyr::mutate(Cp = C * 0.001,	
                Np = N * 0.001,
                Pp = P * 0.001)

# join both datasets
df.mass.CNP <- full_join(df.mass.long, df.CNP.info.percent, by=c("plotID","month2","species","composition"))

df.mass.CNP$div<-as.factor(df.mass.CNP$tree_species_richness)
df.mass.CNP$blk<-as.factor(df.mass.CNP$block)
df.mass.CNP$myc<-as.factor(df.mass.CNP$mycorrhizal_type)
df.mass.CNP$sr<-df.mass.CNP$tree_species_richness
df.mass.CNP$myc <- recode_factor(df.mass.CNP$myc, "AMF" ="AM", "EMF" ="EM", "AMF+EMF" = "AM + EM")



# sum(concentration_species * biomass_species)/total biomass

# we need to consider the proportion of litter per species on each plot first
# otherwise we don't know the true backflow of nutrients to the ground

# Cflux and Nflux and Pflux per species
df.mass.CNP$Cflux <- df.mass.CNP$Cp * df.mass.CNP$mass
df.mass.CNP$Nflux <- df.mass.CNP$Np * df.mass.CNP$mass
df.mass.CNP$Pflux <- df.mass.CNP$Pp * df.mass.CNP$mass




# (2) sum per plot and month (across species) & cumulative effects ####


# exclude rows with NAs/missing data
df.mass.CNP.nona <- df.mass.CNP %>%
  filter(!(is.na(C) & is.na(N) & is.na(P)))

df.flux.plotmonth <-
  #df.mass.CNP.nona %>%
  df.mass.CNP %>%
  group_by(plotID, blk, month2, sr, myc, composition) |>
  summarise(Cflux.tot = sum(Cflux, na.rm = T),
            Nflux.tot = sum(Nflux, na.rm = T),
            Pflux.tot = sum(Pflux, na.rm = T),
            mass.tot = sum(mass, na.rm = T)) %>%
  mutate(Cc=(Cflux.tot/mass.tot)*1000)%>%
  mutate(Nc=(Nflux.tot/mass.tot)*1000)%>%
  mutate(Pc=(Pflux.tot/mass.tot)*1000)%>%
  mutate(CN=Cc/Nc)%>%
  mutate(CP=Cc/Pc)%>%
  mutate(NP=Nc/Pc)%>%
  mutate(across(c(Cc, Nc, Pc), ~na_if(., 0))) #turn zeros in C, N, and P concentration to NA


df.flux.plotmonth$month2 = factor(df.flux.plotmonth$month2, levels = 
                                    c( month.abb[3:12],  month.abb[1:2]))

df.flux.plotmonth$myc <- factor(df.flux.plotmonth$myc, levels = c("AM", "EM", "AM + EM"))

df.flux.plotmonth2 <-
  df.flux.plotmonth |>
  group_by(plotID, blk, sr, myc, composition) |>
  arrange(month2) |>
  mutate(Cflux.cs = cumsum(Cflux.tot),
         Nflux.cs = cumsum(Nflux.tot),
         Pflux.cs = cumsum(Pflux.tot)) |>
  ungroup()


### PLOTTING with stats for the asterisks in figure ####


df.flux.plotmonth2$months <- as.numeric(df.flux.plotmonth2$month2)
months <- sort(unique(df.flux.plotmonth2$months))

df.flux.plotmonth2$div<-as.factor(df.flux.plotmonth2$sr)

# Vector of month labels starting from March
month_labels <- c("Mar", "Apr", "May", "Jun", "Jul", "Aug",
                  "Sep", "Oct", "Nov", "Dec", "Jan", "Feb")

# Calculate monthly means per mycorrhizal type ####
df.monthly.mean.myc<- df.flux.plotmonth2 %>%
  group_by(months, myc) %>%
  summarise(Cc_mean = mean(Cc, na.rm = TRUE),
            Nc_mean = mean(Nc, na.rm = TRUE),
            Pc_mean = mean(Pc, na.rm = TRUE),
            Cflux.tot_mean = mean(Cflux.tot, na.rm = TRUE),
            Nflux.tot_mean = mean(Nflux.tot, na.rm = TRUE),
            Pflux.tot_mean = mean(Pflux.tot, na.rm = TRUE),
            Cflux.cs_mean = mean(Cflux.cs, na.rm = TRUE),
            Nflux.cs_mean = mean(Nflux.cs, na.rm = TRUE),
            Pflux.cs_mean = mean(Pflux.cs, na.rm = TRUE),
            .groups = "drop")

# Calculate monthly means per tree species richness ####
df.monthly.mean.div <- df.flux.plotmonth2 %>%
  group_by(months, div) %>%
  summarise(Cc_mean = mean(Cc, na.rm = TRUE),
            Nc_mean = mean(Nc, na.rm = TRUE),
            Pc_mean = mean(Pc, na.rm = TRUE),
            Cflux.tot_mean = mean(Cflux.tot, na.rm = TRUE),
            Nflux.tot_mean = mean(Nflux.tot, na.rm = TRUE),
            Pflux.tot_mean = mean(Pflux.tot, na.rm = TRUE),
            Cflux.cs_mean = mean(Cflux.cs, na.rm = TRUE),
            Nflux.cs_mean = mean(Nflux.cs, na.rm = TRUE),
            Pflux.cs_mean = mean(Pflux.cs, na.rm = TRUE),
            .groups = "drop")


library(ggplot2)
library(dplyr)
library(patchwork)  # for combining plots

# Define month labels and ensure numeric months
month_labels <- c("Mar", "Apr", "May", "Jun", "Jul", "Aug", 
                  "Sep", "Oct", "Nov", "Dec", "Jan", "Feb")

response_vars <- c("Cc", "Nc", "Pc", "Cflux.tot", "Nflux.tot", "Pflux.tot", "Cflux.cs", "Nflux.cs", "Pflux.cs") 
#"CN", "CP", "NP")

# Initialize an empty data frame to store results
p_values <- data.frame(
  Month = character(),
  Month_Name = character(),
  Response = character(),
  Variable = character(),
  DF = numeric(),
  DDF = numeric(),
  F_Value = numeric(),
  P_Value = numeric(),
  Stars = character(),
  stringsAsFactors = FALSE
)

# Loop over each response variable
for (response in response_vars) {
  for (month in unique(df.flux.plotmonth2$months)) {
    
    dd_month <- df.flux.plotmonth2 %>%
      filter(months == month & !is.na(.data[[response]]))
    
    if (nrow(dd_month) > 0) {
      model <- tryCatch({
        lme(as.formula(paste(response, "~ I(log2(sr)) * myc")), 
            random = ~1 | composition, 
            data = dd_month)
      }, error = function(e) {
        warning(paste("Error fitting model for month:", month, 
                      "and response:", response, ":", e$message))
        return(NULL)
      })
      
      if (is.null(model)) next
      
      anova_res <- anova(model)
      terms <- c("I(log2(sr))", "myc", "I(log2(sr)):myc")
      
      for (term in terms) {
        if (term %in% rownames(anova_res)) {
          p_val <- anova_res[term, "p-value"]
          F_val <- anova_res[term, "F-value"]
          df_val <- anova_res[term, "numDF"]
          ddf_val <- anova_res[term, "denDF"]
          
          if (!is.na(p_val)) {
            stars <- ifelse(p_val < 0.001, "***",
                            ifelse(p_val < 0.01, "**",
                                   ifelse(p_val < 0.05, "*",
                                          ifelse(p_val < 0.1, "·", "n.s."))))
            month_name <- month_labels[as.integer(month)]
            
            new_row <- data.frame(
              Month = month,
              Month_Name = month_name,
              Response = response,
              Variable = term,
              DF = df_val,
              DDF = ddf_val,
              F_Value = F_val,
              P_Value = p_val,
              Stars = stars,
              stringsAsFactors = FALSE
            )
            
            p_values <- rbind(p_values, new_row)
          }
        }
      }
    }
  }
}

# Optional: order by month and response
p_values <- p_values %>% arrange(Response, Month)

# Check results
head(p_values)

# Print the output data frame to see all p-values and significance
print(p_values)

p_values$Variable <- recode(p_values$Variable,
                            "I(log2(sr))" = "div",
                            "I(log2(sr)):myc" = "div:myc")



library(ggplot2)
library(rlang)
install.packages("scales") 
library(scales)

plot_biodiv_effect <- function(df_raw, df_mean, yvar, predictor, predictor_label, colors, ylab_text, tag_label, p_values_df) {
  
  # Filter p_values for this plot's response variable and predictor
  pvals_subset <- p_values_df %>%
    filter(Response == yvar, Variable == predictor)
  
  y_min <- min(df_raw[[yvar]], na.rm = TRUE)
  y_max <- max(df_raw[[yvar]], na.rm = TRUE)
  y_buffer <- (y_max - y_min) * 0.5
  y_star_pos <- y_max + y_buffer
  
  # Create data frame for star annotations per month
  stars_df <- data.frame(
    months = pvals_subset$Month,
    stars = pvals_subset$Stars,
    y = y_star_pos
  )
  
  # Add size category: larger for *, **, ***; small for n.s.
  stars_df$size <- ifelse(stars_df$stars == "n.s.", 4, 6)  # Adjust sizes as needed
  
  ggplot() +
    geom_point(data = df_raw, aes(x = months, y = .data[[yvar]], color = .data[[predictor]]), alpha = 0.5) +
    geom_jitter(data = df_raw, aes(x = months, y = .data[[yvar]], color = .data[[predictor]]),
      width = 0.3,   # small horizontal jitter
      height = 0,    # no vertical jitter (optional)
      alpha = 0.3)+
    geom_line(data = df_mean, aes(x = months, y = .data[[paste0(yvar, "_mean")]], color = .data[[predictor]]), size = 1.2) +
    geom_text(data = stars_df, aes(x = months, y = y, label = stars, size = factor(size)), color = "black") +  # add stars
    scale_color_manual(values = colors, name = predictor_label) +
    scale_y_continuous(limits = c(y_min, y_star_pos * 1.05),
                       labels = scales::label_number(big.mark = "")) +  
    scale_x_discrete(limits = month_labels)+
    scale_size_manual(values = c("4" = 4, "6" = 6), guide = "none") +
    ylab(ylab_text) +
    xlab(NULL) +
    theme_minimal() +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(size = 12),
      axis.line = element_line(color = 'black'),
      axis.text.y = element_text(color = "black", size = 12),
      axis.text.x = element_text(color = "black", size = 12),
      axis.title.y = element_text(size = 12),
      axis.title.x = element_blank(),
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
      legend.title = element_text("Biodiversity effects", size = 12),
      legend.text = element_text(size = 12),
      legend.background = element_rect(colour = NA),
      legend.box = NULL,
      plot.tag = element_text()
    ) +
    labs(tag = tag_label)
}


# Variable groups
concentration_vars <- c("Cc", "Nc", "Pc")
flux_vars <- c("Cflux.tot", "Nflux.tot", "Pflux.tot")
cs_vars <- c("Cflux.cs", "Nflux.cs", "Pflux.cs")

# Colors
div_colors <- c("#dda1ba", "#a93c69", "#5e223a")
myc_colors <- c("#71b540","#4c8ecb","#febf00")

var_labels <- c(
  Cc = expression(atop("C concentration", (mg ~ g^-1))),
  Nc = expression(atop("N concentration", (mg ~ g^-1))),
  Pc = expression(atop("P concentration", (mg ~ g^-1))),
  Cflux.tot = expression(atop("C flux", (kg ~ ha^-1 ~ yr^-1))),
  Nflux.tot = expression(atop("N flux", (kg ~ ha^-1 ~ yr^-1))),
  Pflux.tot = expression(atop("P flux", (kg ~ ha^-1 ~ yr^-1))),
  Cflux.cs = expression(atop("C flux (cumulative)", (kg ~ ha^-1 ~ yr^-1))),
  Nflux.cs = expression(atop("N flux (cumulative)", (kg ~ ha^-1 ~ yr^-1))),
  Pflux.cs = expression(atop("P flux (cumulative)", (kg ~ ha^-1 ~ yr^-1)))
)

# Assume you have monthly means calculated:
# df.monthly.mean.div: for species richness
# df.monthly.mean.myc: for mycorrhizal type

# Helper function to generate plots
make_plots <- function(var_list, y_suffix, start_letter = 1) {
  plot_list <- list()
  
  for (i in seq_along(var_list)) {
    var <- var_list[i]
    
    # Richness effect
    p1 <- plot_biodiv_effect(
      df_raw = df.flux.plotmonth2,
      df_mean = df.monthly.mean.div,
      yvar = var,
      predictor = "div",
      predictor_label = "Tree species richness",
      colors = div_colors,
      ylab_text = var_labels[[var]],
      tag_label = paste0("(", letters[i], ")"),
      p_values_df = p_values
    )
    
    # Mycorrhizal effect
    p2 <- plot_biodiv_effect(
      df_raw = df.flux.plotmonth2,
      df_mean = df.monthly.mean.myc,
      yvar = var,
      predictor = "myc",
      predictor_label = "Mycorrhizal type",
      colors = myc_colors,
      ylab_text = var_labels[[var]],
      tag_label = paste0("(", letters[start_letter + i + length(var_list) - 1], ")"),
      p_values_df = p_values
    )
    
    plot_list[[length(plot_list) + 1]] <- p1
    plot_list[[length(plot_list) + 1]] <- p2
  }
  
  return(plot_list)
}


# Generate plots
plots_conc <- make_plots(concentration_vars, "concentration", start_letter = 1)
plots_flux <- make_plots(flux_vars, "flux", start_letter = 1)
plots_cs   <- make_plots(cs_vars, "cumulative flux", start_letter = 1)

# Combine using patchwork
library(patchwork)

concentration_figure <- wrap_plots(plots_conc, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

flux_figure <- wrap_plots(plots_flux, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

cs_figure <- wrap_plots(plots_cs, ncol = 2) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom", legend.direction = "horizontal")

# Show one or all as needed
concentration_figure
flux_figure
cs_figure


### Boxplots Mycorrhizal type - sum sum ####
#### (1) C flux cumulative sum ####
Ccs.myc <-ggplot() +
  geom_point(data = df.flux.plotmonth2, aes(x = months, y = Cflux.cs, color = myc), alpha = 0.5) +
  geom_line(data = df.monthly.mean.myc, aes(x = months, y = Cflux.cs_mean, color = myc), size = 1.2) +
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"),
                     name = "Mycorrhizal type") +
  scale_x_continuous(
    name = "Month",
    breaks = 1:12,
    labels = month_labels
  ) +
  ylab("C flux cumulative sum") +
  labs(fill = "Tree species richness") +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.line = element_line(color = 'black'),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
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
    legend.title = element_text("Biodiversity effects", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_rect(colour = NA),
    legend.box = NULL,
    legend.box.background = element_rect(color = "transparent")
  )+ labs(
    tag = "(d)")
Ccs.myc

#### (2) N flux cumulative sum ####
Ncs.myc <-ggplot() +
  geom_point(data = df.flux.plotmonth2, aes(x = months, y = Nflux.cs, color = myc), alpha = 0.5) +
  geom_line(data = df.monthly.mean.myc, aes(x = months, y = Nflux.cs_mean, color = myc), size = 1.2) +
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"),
                     name = "Mycorrhizal type") +
  scale_x_continuous(
    name = "Month",
    breaks = 1:12,
    labels = month_labels
  ) +
  ylab("N flux cumulative sum") +
  labs(fill = "Tree species richness") +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.line = element_line(color = 'black'),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
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
    legend.title = element_text("Biodiversity effects", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_rect(colour = NA),
    legend.box = NULL,
    legend.box.background = element_rect(color = "transparent")
  )+ labs(
    tag = "(e)")
Ncs.myc

#### (3) flux cumulative sum ####
Pcs.myc 
ggplot() +
  #geom_boxplot(data = df.flux.plotmonth2, 
               # aes(x = months, y = Pflux.cs, color = myc, group = interaction(months, myc)), 
               # alpha = 0.5)+
  geom_point(data = df.flux.plotmonth2, aes(x = months, y = Pflux.cs, color = myc), alpha = 0.5) +
  geom_line(data = df.monthly.mean.myc, aes(x = months, y = Pflux.cs_mean, color = myc), size = 1.2) +
  scale_color_manual(values = c("#71b540","#4c8ecb","#febf00"),
                     name = "Mycorrhizal type") +
  scale_x_continuous(
    name = "Month",
    breaks = 1:12,
    labels = month_labels
  ) +
  ylab("P flux cumulative sum") +
  labs(fill = "Tree species richness") +
  theme_minimal() +
  theme(
    strip.background = element_blank(),
    strip.text = element_text(size = 12),
    axis.line = element_line(color = 'black'),
    axis.text.y = element_text(color = "black", size = 12),
    axis.text.x = element_text(color = "black", size = 12),
    axis.title.y = element_text(size = 12),
    axis.title.x = element_blank(),
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
    legend.title = element_text("Biodiversity effects", size = 12),
    legend.text = element_text(size = 12),
    legend.background = element_rect(colour = NA),
    legend.box = NULL,
    legend.box.background = element_rect(color = "transparent")
  )+ labs(
    tag = "(f)")
Pcs.myc


