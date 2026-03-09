#---------------------------------------------------------------------
# MyDiv experiment; Litterfall project
# 2025-10-28
# Precipitation 2018-2023
# by Elisabeth Boenisch (elisabeth.boenisch@idiv.de)

#============================ Packages ==============================

rm(list=ls())

library(readr)
library(readxl)
library(tidyverse)
library(nlme)
library(lme4)
library(lmerTest)
library(ggplot2)

#============================ Dataset ===============================

precipitation <- read.delim2 #data from DWD (Deutscher Wetterdienst) for Bad Lauchstädt (Wetterstation 02878; 2018-2024, monthly data)
# https://opendata.dwd.de/climate_environment/CDC/observations_germany/climate/


#================= Definition of Abbreviations ======================

# MO_NSH  Monatssumme der Neuschneehoehe, cm 
# MO_RR   Monatssumme der Niederschlagshoehe, mm
# MO_SH_S Monatssumme der Schneehoehe, cm
# MX_RS   Maximale Niderschlagshoehe des Monats, mm

head(precipitation)

precipitation_recent <- precipitation %>%
  mutate(MESS_DATUM_BEGINN = as.Date(as.character(MESS_DATUM_BEGINN), format = "%Y%m%d")) %>%
  filter(MESS_DATUM_BEGINN >= as.Date("2018-01-01"))

precipitation_recent$MO_RR <- as.numeric(precipitation_recent$MO_RR)

ggplot(precipitation_recent, aes(x = MESS_DATUM_BEGINN, y = MO_RR)) +
  geom_line(color = "blue") +
  labs(x = "Date", y = "Monthly precipitation (mm)", 
       title = "Precipitation over time") +
  theme_minimal()

line.plot <- ggplot(precipitation_recent, aes(x = MESS_DATUM_BEGINN, y = MO_RR)) +
  geom_line(color = "steelblue", linewidth = 0.8) +
  geom_point(size = 1.2, color = "steelblue") +
  geom_smooth(method = "loess", se = FALSE, color = "darkred", linetype = "dashed") +
  scale_x_date(date_breaks = "1 year", date_labels = "%Y") +
  labs(x = "Date", y = "Monthly precipitation (mm)",
       title = "Monthly Precipitation Over Time") +
  theme_minimal(base_size = 13)+
  labs(tag = "(a)")


###########

library(dplyr)
library(tidyr)
library(lubridate)

# Generate full sequence of months
full_months <- data.frame(
  MESS_DATUM_BEGINN = seq(min(precipitation_recent$MESS_DATUM_BEGINN, na.rm = TRUE),
                          max(precipitation_recent$MESS_DATUM_BEGINN, na.rm = TRUE),
                          by = "month")
)

# Join with your data (preserves all months, fills missing with NA)
precipitation_full <- full_months %>%
  left_join(precipitation_recent, by = "MESS_DATUM_BEGINN")

############

# Make sure date and MO_RR are correctly formatted
precipitation_full <- precipitation_full %>%
  mutate(MESS_DATUM_BEGINN = as.Date(MESS_DATUM_BEGINN),
         year = format(MESS_DATUM_BEGINN, "%Y"),
         MO_RR = as.numeric(MO_RR))

# Sum precipitation per year
precip_annual <- precipitation_full %>%
  group_by(year) %>%
  summarise(annual_sum = sum(MO_RR, na.rm = TRUE))

bar.plot <- ggplot(precip_annual, aes(x = as.integer(year), y = annual_sum)) +
  geom_col(fill = "steelblue", color = "gray40") +
  geom_text(aes(label = round(annual_sum, 1)), 
            vjust = -0.3, size = 3, color = "black") +
  labs(x = NULL, y = "Annual precipitation (mm)", 
       title = "Sum of Monthly Precipitation per Year") +
  theme_minimal(base_size = 13) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank()) +
  scale_x_continuous(breaks = seq(min(as.integer(precip_annual$year)),
                                  max(as.integer(precip_annual$year)), 1))+
  labs(tag = "(b)")

library(patchwork)


plot <- line.plot / bar.plot

### end ###
