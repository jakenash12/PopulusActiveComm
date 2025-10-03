#This script calculates the mean absolute change of soil moisture and temperature.
#This is a metric of seasonal variability in soil moisture/temp
#It is calculated by summing the absolute value of the difference in soil temp or moisture
#between consecutive days and then dividing this sum by the total number of days in the series
#First, daily means of soil moisture and temp are calculated and then mean absolute change
#is calculated from these daily means

#Need to run the script SoilSensorDataProcessing.R prior to this script

library(tidyverse)
library(car)

###############################Soil Temp#########################################

#calculates daily means of soil temp
interpolated_df_soiltemp_daily=
  interpolated_df_soiltemp %>%
  mutate(Site=str_extract(Plot, "^.."),
         Date=as.Date((Time))) %>% 
  group_by(Plot, Date) %>%
  dplyr::summarise(MeanDailySoilTemp=mean(NewVals)) %>%
  filter(Date>"2021-05-24", Date<"2021-10-08")

#initializes empty data frame where values will be deposited
Summary_soiltemp=data.frame(
  Plot = character(0),
  MAC_soiltemp = numeric(0),
  min_soiltemp = numeric(0),
  max_soiltemp = numeric(0),
  mean_soiltemp = numeric(0),
  diel_soiltemp = numeric(0),
  seasonal_range_soiltemp=numeric(0)
)

#uses loop to calculate summary stats for each plot
for (i in unique(interpolated_df_soiltemp_daily$Plot)){
  temp_df=
    interpolated_df_soiltemp_daily %>%
    filter(Plot==i, !(is.na(MeanDailySoilTemp))) %>%
    arrange(Date) %>%  # Ensure the data is ordered by Date
    mutate(
      AbsDifference = abs(MeanDailySoilTemp - lag(MeanDailySoilTemp))
    )
  
  diel_var_df =
    interpolated_df_soiltemp %>%
    mutate(Site=str_extract(Plot, "^.."),
           Date=as.Date((Time))) %>% 
    group_by(Plot, Date) %>%
    filter(Date>"2021-05-24", Date<"2021-10-08") %>%
    filter(Plot==i, !(is.na(NewVals))) %>%
    group_by(Date) %>%
    summarise(Min=min(NewVals),
              Max=max(NewVals)) %>%
    mutate(Range=abs(Max-Min))
  
  diel_range=mean(diel_var_df$Range)
  min_soiltemp=min(temp_df$MeanDailySoilTemp)
  max_soiltemp=max(temp_df$MeanDailySoilTemp)
  mean_soiltemp=mean(temp_df$MeanDailySoilTemp)
  mean_abs_ch=(sum(temp_df$AbsDifference, na.rm=TRUE))/(as.numeric(as.Date("2021-10-08") - as.Date("2021-05-24")))
  new_df=data.frame(Plot=i, 
                    MAC_soiltemp=mean_abs_ch, 
                    min_soiltemp=min_soiltemp,
                    max_soiltemp=max_soiltemp,
                    mean_soiltemp=mean_soiltemp,
                    diel_soiltemp=diel_range,
                    seasonal_range_soiltemp=max_soiltemp-min_soiltemp)
  Summary_soiltemp=rbind(Summary_soiltemp, new_df)
}

###############################Soil Moisture#########################################
#calculates daily means of soil moisture
interpolated_df_moist_daily=
  interpolated_df_moist %>%
  mutate(Site=str_extract(Plot, "^.."),
         Date=as.Date((Time))) %>% 
  group_by(Plot, Date) %>%
  dplyr::summarise(MeanDailyMoist=mean(NewVals)) %>%
  filter(Date>"2021-05-24", Date<"2021-10-08")

#initializes empty dataframe where summary stats will be deposited
Summary_moist=data.frame(
  Plot = character(0),
  MAC_moist = numeric(0),
  min_moist = numeric(0),
  max_moist = numeric(0),
  mean_moist = numeric(0),
  diel_moist = numeric(0),
  seasonal_range_moist=numeric(0)
)

#uses loop to calculate summary stats for each plot
for (i in unique(interpolated_df_moist_daily$Plot)){
  temp_df=
    interpolated_df_moist_daily %>%
    filter(Plot==i, !(is.na(MeanDailyMoist))) %>%
    arrange(Date) %>%  # Ensure the data is ordered by Date
    mutate(
      AbsDifference = abs(MeanDailyMoist - lag(MeanDailyMoist))
    )
  
  diel_var_df =
    interpolated_df_moist %>%
    mutate(Site=str_extract(Plot, "^.."),
           Date=as.Date((Time))) %>% 
    group_by(Plot, Date) %>%
    filter(Date>"2021-05-24", Date<"2021-10-08") %>%
    filter(Plot==i, !(is.na(NewVals))) %>%
    group_by(Date) %>%
    summarise(Min=min(NewVals),
              Max=max(NewVals)) %>%
    mutate(Range=abs(Max-Min))
  
  diel_range=mean(diel_var_df$Range)
  min_soilmoist=min(temp_df$MeanDailyMoist)
  max_soilmoist=max(temp_df$MeanDailyMoist)
  mean_soilmoist=mean(temp_df$MeanDailyMoist)
  mean_abs_ch=(sum(temp_df$AbsDifference, na.rm=TRUE))/(as.numeric(as.Date("2021-10-08") - as.Date("2021-05-24")))
  new_df=data.frame(Plot=i, 
                    MAC_moist=mean_abs_ch, 
                    min_moist=min_soilmoist,
                    max_moist=max_soilmoist,
                    mean_moist=mean_soilmoist,
                    diel_moist=diel_range,
                    seasonal_range_moist=max_soilmoist-min_soilmoist)
  Summary_moist=rbind(Summary_moist, new_df)
}

################################################################################

