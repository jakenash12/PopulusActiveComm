library(stringr)
library(readxl)
library(dplyr)
library(anytime)
library(tidyverse)
library(lubridate)
library(tidyr)
library(simtimer)
library(tibbletime)
library(ggplot2)
library(cowplot)
library(car)
library(magrittr)

path="./InputFiles/Teros11DataDownload"
filelist=list.files(path)

#custom function to rename columns as port # + measurement type
renamecols <- function(df){
  for (g in ncol(df):2){
    colnames(df)[g] = 
      if ((g%%2)==0){
        paste(colnames(df)[g],df[1,g])
      } else {
        paste(colnames(df)[g-1],df[1,g])
      }
  }
  colnames(df)[1]="Time"
  return(df)
}

filtercols <- function(df){
  vect = !grepl("Sensor Output", df[1,])
  return(df[vect])
}

Fulldf=data.frame()
origin_date <- as.POSIXct("1900-01-01 00:00:00", tz = "UTC")
for (i in filelist) {
  Site=str_extract(i, "^..")
  sheets=
    paste(path, i, sep="/") %>%
    excel_sheets
  sheets=sheets[sheets!="Metadata"]
  sheets=sheets[!str_detect(sheets, regex("Raw", ignore_case = TRUE))]
  
  for (z in sheets){
    print(i)
    WorkingDF1=
      read_excel((paste(path, i, sep="/")), sheet=z,guess_max = 1048576) %>%
      slice(-1) %>% #deletes the first row
      filtercols %>%
      renamecols  %>%
      slice(-1, -2) %>%
      mutate(Time = round_date(as.datetime(as.numeric(Time)*86400, origin_date = origin_date)), unit="minute") %>%
      select(-unit)
    print(colnames(WorkingDF1))
    WorkingDF=
      WorkingDF1 %>%
      gather(measurement,value,2:`Port8 °C Logger Temperature`) %>%
      mutate(Site = Site, 
             Port = str_extract(.$measurement, "Port."),
             DataType = str_extract(.$measurement, "(?<=Port.)[^+]+")) %>%
      mutate(Plot=paste(Site, str_extract(.$Port, "(?<=Port)[^+]+"),sep=""),
             SiteType=str_extract(.$Site, "^."))
    Fulldf=rbind(Fulldf, WorkingDF)
  }
}

#filters df by time window from 5/19/21 to 10/11/21
Fulldf %<>%
  filter(Time<"2021-10-11", Time >"2021-05-19") %>% #subsets data by a start and end date 
  mutate(value=as.numeric(value))

#subsets dataframe to only include data collected by sensors
#which were plugged into ports 2, 3, or 4
Sensordf=
  Fulldf %>%
  filter((Port=="Port2"|Port=="Port3"|Port=="Port4"))


###################Interpolating Data################################
#Here I use a data interpolation method to fill in missing data based on 
#linear models
#The basic idea is that missing data points at a give plot and time point
#are predicted based on values at other plots at the same point
#
#Ideally, data from a missing timepoint are predicted based on the data values
#at the other two plots within that same site. However, in some cases
#data are missing from all three plots within a site. 
#


##################Data interpolation for soil moisture#################
#creates wide format soil moisture data frame
Moistdf_wide =
  Sensordf %>%
  filter(DataType==" m³/m³ Water Content") %>% #selects only data that is of type soil moisture
  filter(!(Plot=="HB4"&Time>"2021-09-14" & Time<"2021-10-12")) %>% #filters out data from plot HB4 during a period when sensor was dug up and sitting on surface of soil
  select(Time, Plot, value) %>% #selects time, plot, and value columns
  distinct(Plot, Time, .keep_all = TRUE) %>% #removes duplicate rows
  filter(value > 0.01) %>% #removes data that is less than 1% soil moisture because it is likely incorrect (mostly just removing zeroes)
  spread(Plot, value) #converts from long to wide format, with time as column and each of the plots as 27 more columns containing soil moisture data

Moistdf_wide=Moistdf_wide[rowSums(is.na(Moistdf_wide))<27,] #removes timepoints for which all plots had NA values - these cannot be interpolated

#makes list of plot names
SiteNames=
  Moistdf_wide %>% 
  select(-Time) %>%
  colnames

#creates an empty dataframe where data for each plot will be appended
colnames_soilmoist=c("Time","OldVals","SameSite1","SameSite2",
                      "PredicType", "Plot", "SameSitesPredic", "SameSite1Predic",
                      "SameSite2Predic", "NewVals", "DataOrigin")
interpolated_df_moist=
  data.frame(matrix(ncol=length(colnames_soilmoist),nrow=0)) %>%
  set_colnames(colnames_soilmoist)

#uses loop to run data interpolation script for each plot
for (i in SiteNames) {
  ActiveSite=str_extract(i, "^..") #identifies which site the plot of interest belongs to by extracting first two characters of plot ID
  SameSites= #identifies names of two other plots from same site as plot of interest 
    SiteNames %>%
    str_subset(pattern = i, negate = TRUE) %>%
    str_subset(pattern=ActiveSite)
  SameSite1=SameSites[1]
  SameSite2=SameSites[2]
  
  tempdf_moist=data.frame(Time=Moistdf_wide$Time,
                    OldVals=Moistdf_wide[i],
                    SameSite1=Moistdf_wide[SameSite1],
                    SameSite2=Moistdf_wide[SameSite2]) %>%
    set_colnames(c("Time","OldVals","SameSite1","SameSite2")) %>%
    mutate(PredicType=ifelse(!is.na(OldVals),"None", 
                             ifelse((!is.na(SameSite1)&!is.na(SameSite2)),"SameSites",
                                    ifelse(!is.na(SameSite1),"SameSite1",
                                           ifelse(!is.na(SameSite2),"SameSite2","Other")))),
           Plot=i)
  
  ############Prediction when 2 plots from same site are available
  predic_mod_SameeSites=lm(OldVals~SameSite1+SameSite2, tempdf_moist)
  tempdf_moist$SameSitesPredic= #predicts data to be interpolated based on linear model
    predict(predic_mod_SameeSites, newdata = tempdf_moist) 
  
  ###########Prediction when only first plot from site is available
  predic_mod_SameeSite1=lm(OldVals~SameSite1, tempdf_moist)
  tempdf_moist$SameSite1Predic= #predicts data to be interpolated based on linear model
    predict(predic_mod_SameeSite1, newdata = tempdf_moist) 
  
  ###########Prediction when only second plot from site is available
  predic_mod_SameeSite2=lm(OldVals~SameSite2, tempdf_moist)
  tempdf_moist$SameSite2Predic= #predicts data to be interpolated based on linear model
    predict(predic_mod_SameeSite2, newdata = tempdf_moist) 
  
  tempdf_moist %<>%
    mutate(NewVals=ifelse(!is.na(OldVals),OldVals,
                          ifelse(!is.na(SameSitesPredic),SameSitesPredic, 
                                 ifelse(!is.na(SameSite1Predic),SameSite1Predic,
                                        ifelse(!is.na(SameSite2Predic),SameSite2Predic,NA)))),
           DataOrigin=ifelse(!is.na(OldVals),"Measured",ifelse((!is.na(SameSitesPredic)|!is.na(SameSite1Predic)|!is.na(SameSite2Predic)),"Interpolated","Missing")))
  
  #Appends temp df from current plot to full df
  interpolated_df_moist=rbind(interpolated_df_moist,tempdf_moist)
}

#creates plot of soil moisture at each plot over time
#with averages by site and lines colored by site type
SoilMoistPlot=
  interpolated_df_moist %>%
  mutate(Site=str_extract(Plot, "^.."),
         Date=as.Date((Time))) %>% 
  group_by(Site, Date) %>%
  dplyr::summarise(MeanDailyMoist=mean(NewVals)) %>% 
  mutate(SiteType=str_extract(Site,"^.")) %>%
  filter(Date>"2021-05-24", Date<"2021-10-08") %>%
  ggplot(aes(x=Date, y=MeanDailyMoist, color=SiteType, group=Site)) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(size=2) +
  geom_vline(xintercept = as.Date("2021-06-04"), linetype=2, size=2) +
  geom_vline(xintercept = as.Date("2021-08-17"), linetype=2, size=2) +
  geom_vline(xintercept = as.Date("2021-10-07"), linetype=2, size=2) +
  theme_test() +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text=element_blank(),
        plot.title = element_text(hjust = 0.5, size=16),
        panel.border=element_rect(size=2),
        axis.ticks=element_line(size=2)) +
  ggtitle("")

##################Data interpolation for soil temperature#################
#converts soil temperature data to wide format
SoilTempdf_wide =
  Sensordf %>%
  filter(DataType==" °C Soil Temperature") %>% #selects only data that is of type soil temperature
  filter(!(Plot=="HB4"&Time>"2021-09-14" & Time<"2021-10-12")) %>% #filters out data from plot HB4 during a period when sensor was dug up and sitting on surface of soil
  select(Time, Plot, value) %>% #selects time, plot, and value columns
  distinct(Plot, Time, .keep_all = TRUE) %>% #removes duplicate rows
  filter(value > 0.01) %>% #removes data that is less than 0.01 C because it is likely incorrect (mostly just removing zeroes)
  spread(Plot, value) #converts from long to wide format, with time as column and each of the plots as 27 more columns containing soil temperature data

SoilTempdf_wide=SoilTempdf_wide[rowSums(is.na(SoilTempdf_wide))<27,] #removes timepoints for which all plots had NA values - these cannot be interpolated

#creates a dataframe of data to exclude. These data looked suspect
#upon visual examination and are filtered prior to interpolation

#creates list of plot names to loop through
SiteNames=
  SoilTempdf_wide %>% 
  select(-Time) %>%
  colnames

#creates an empty dataframe where data for each plot will be appended
colnames_soiltemp=c("Time","OldVals","SameSite1","SameSite2",
                     "PredicType", "Plot", "SameSitesPredic", "SameSite1Predic",
                     "SameSite2Predic", "NewVals", "DataOrigin")
interpolated_df_soiltemp=
  data.frame(matrix(ncol=length(colnames_soiltemp),nrow=0)) %>%
  set_colnames(colnames_soiltemp)

#uses loop to run data interpolation for each plot
for (i in SiteNames) {
  ActiveSite=str_extract(i, "^..") #identifies which site the plot of interest belongs to by extracting first two characters of plot ID
  SameSites= #identifies names of two other plots from same site as plot of interest 
    SiteNames %>%
    str_subset(pattern = i, negate = TRUE) %>%
    str_subset(pattern=ActiveSite)
  SameSite1=SameSites[1]
  SameSite2=SameSites[2]
  
  tempdf_temp=data.frame(Time=SoilTempdf_wide$Time,
                          OldVals=SoilTempdf_wide[i],
                          SameSite1=SoilTempdf_wide[SameSite1],
                          SameSite2=SoilTempdf_wide[SameSite2]) %>%
    set_colnames(c("Time","OldVals","SameSite1","SameSite2")) %>%
    mutate(PredicType=ifelse(!is.na(OldVals),"None", 
                             ifelse((!is.na(SameSite1)&!is.na(SameSite2)),"SameSites",
                                    ifelse(!is.na(SameSite1),"SameSite1",
                                           ifelse(!is.na(SameSite2),"SameSite2","Other")))),
           Plot=i)
  
  ############Prediction when 2 plots from same site are available
  predic_mod_SameeSites=lm(OldVals~SameSite1+SameSite2, tempdf_temp)
  tempdf_temp$SameSitesPredic= #predicts data to be interpolated based on linear model
    predict(predic_mod_SameeSites, newdata = tempdf_temp) 
  
  ###########Prediction when only first plot from site is available
  predic_mod_SameeSite1=lm(OldVals~SameSite1, tempdf_temp)
  tempdf_temp$SameSite1Predic= #predicts data to be interpolated based on linear model
    predict(predic_mod_SameeSite1, newdata = tempdf_temp) 
  
  ###########Prediction when only second plot from site is available
  predic_mod_SameeSite2=lm(OldVals~SameSite2, tempdf_temp)
  tempdf_temp$SameSite2Predic= #predicts data to be interpolated based on linear model
    predict(predic_mod_SameeSite2, newdata = tempdf_temp) 
  
  tempdf_temp %<>%
    mutate(NewVals=ifelse(!is.na(OldVals),OldVals,
                          ifelse(!is.na(SameSitesPredic),SameSitesPredic, 
                                 ifelse(!is.na(SameSite1Predic),SameSite1Predic,
                                        ifelse(!is.na(SameSite2Predic),SameSite2Predic,NA)))),
           DataOrigin=ifelse(!is.na(OldVals),"Measured",ifelse((!is.na(SameSitesPredic)|!is.na(SameSite1Predic)|!is.na(SameSite2Predic)),"Interpolated","Missing")))
  
  #Appends temp df from current plot to full df
  interpolated_df_soiltemp=rbind(interpolated_df_soiltemp,tempdf_temp)
}

#makes a plot of average soil temperature over time 
#average by site, colored by site type
SoilTempPlot=
  interpolated_df_soiltemp %>%
  mutate(Site=str_extract(Plot, "^.."),
         Date=as.Date((Time))) %>% 
  group_by(Site, Date) %>%
  dplyr::summarise(MeanDailySoilTemp=mean(NewVals)) %>% 
  mutate(SiteType=str_extract(Site,"^.")) %>%
  filter(Date>"2021-05-24", Date<"2021-10-08") %>%
  ggplot(aes(x=Date, y=MeanDailySoilTemp, color=SiteType, group=Site)) +
  scale_color_brewer(palette = "Dark2") +
  geom_line(size=2) +
  geom_vline(xintercept = as.Date("2021-06-04"), linetype=2, size=2) +
  geom_vline(xintercept = as.Date("2021-08-17"), linetype=2, size=2) +
  geom_vline(xintercept = as.Date("2021-10-07"), linetype=2, size=2) +
  theme_test()  +
  theme(legend.position = "none",
        axis.title=element_blank(),
        axis.text=element_blank(),
        plot.title = element_text(hjust = 0.5, size=16),
        panel.border=element_rect(size=2),
        axis.ticks=element_line(size=2)) +
  ggtitle("")

#makes a 2-panel plot showing soil temperature and moisture
#over the season
plot_grid(SoilTempPlot, SoilMoistPlot, nrow=2)

#calculates the % of data points that are measured,
#interpolated, and missing (i.e. could not be interpolated)
table(interpolated_df_soiltemp$DataOrigin) %>%
  as.data.frame %>%
  mutate(Proportion=Freq/sum(Freq))
#########################################################

##############Creating summary dataframes#############################
#Here, interpolated data is formatted into summary dataframes that include
#seasonal averages, as well as averages in time windows before each sampling date

#Start and end date cutoffs for data filtering
#Note that these dates are actually a discrete cutoff at 12 am
#This means that start date is inclusive, and end date is non-inclusive
StartDate="2021-05-24"
EndDate="2021-10-08"

#creates a dataframe with seasonal averages of 
#soil temperature by plot
PlotAvgSoilTemp_temp=
  interpolated_df_soiltemp %>%
  filter(Time>StartDate,Time<EndDate) %>%
  group_by(Plot) %>%
  summarise(plotavg_soiltemp=mean(NewVals, na.rm=TRUE),
            sd_soiltemp=sd(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."))

PlotAvgSoilTemp=
  interpolated_df_soiltemp %>%
  filter(Time>StartDate,Time<EndDate) %>%
  mutate(Day=as.Date(Time)) %>%
  group_by(Plot, Day) %>%
  summarise(SoilTempDailyAvg=mean(NewVals)) %>%
  group_by(Plot) %>%
  summarise(sd_dailyavgs_soiltemp=sd(SoilTempDailyAvg, na.rm=TRUE)) %>%
  left_join(PlotAvgSoilTemp_temp)

#creates dataframe of seasonal averages of soil
#moisture grouped by plot
PlotAvgSoilMoist_temp=
  interpolated_df_moist %>%
  filter(Time>StartDate,Time<EndDate) %>%
  group_by(Plot) %>%
  summarise(plotavg_moist=mean(NewVals, na.rm=TRUE),
            sd_moist=sd(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."))

PlotAvgSoilMoist=
  interpolated_df_moist %>%
  filter(Time>StartDate,Time<EndDate) %>%
  mutate(Day=as.Date(Time)) %>%
  group_by(Plot, Day) %>%
  summarise(MoistDailyAvg=mean(NewVals)) %>%
  group_by(Plot) %>%
  summarise(sd_dailyavgs_moist=sd(MoistDailyAvg, na.rm=TRUE)) %>%
  left_join(PlotAvgSoilMoist_temp)
  
#creates a dataframe with the average within-day
#fluctuation in soil moisture
Diurnal_Moisture_Summary =
  Fulldf %>%
  filter(DataType==" m³/m³ Water Content", Time<"2021-10-07", Time>"2021-09-07") %>%
  mutate(value=as.numeric(value),
         Date = as.Date(Time)) %>%
  group_by(Plot, Date) %>%
  summarise(avg_moist=mean(value),
            min_moist=min(value),
            max_moist=max(value),
            DiurnalVar=max(value)-min(value)) %>%
  filter(avg_moist!=0) %>%
  ungroup %>%
  group_by(Plot) %>%
  summarise(avg_moist=mean(avg_moist),
            min_moist=mean(min_moist),
            max_moist=mean(max_moist),
            DiurnalVar=mean(DiurnalVar)) %>%
  mutate(SiteType=str_extract(.$Plot, "^."), 
         PlotNum=str_extract(.$Plot, ".$"),
         Site=str_extract(.$Plot, "^..")) %>%
  filter(PlotNum<5) %>%
  mutate(SiteType=dplyr::recode(.$SiteType, H="High Elevation", 
                                M="Mid Elevation",
                                R="Riparian"))

#creates a dataframe with the average within-day
#fluctuation in soil temperature
Diurnal_Temp_Summary =
  Fulldf %>%
  filter(DataType==" °C Soil Temperature", Time < "2021-08-04") %>%
  mutate(value=as.numeric(value),
         Date = as.Date(Time)) %>%
  group_by(Plot, Date) %>%
  summarise(avg_temp=mean(value),
            min_temp=min(value),
            max_temp=max(value),
            DiurnalVar=max(value)-min(value)) %>%
  filter(avg_temp!=0) %>%
  ungroup %>%
  group_by(Plot) %>%
  summarise(avg_temp=mean(avg_temp),
            min_temp=mean(min_temp),
            max_temp=mean(max_temp),
            DiurnalVar=mean(DiurnalVar)) %>%
  mutate(SiteType=str_extract(.$Plot, "^."), 
         PlotNum=str_extract(.$Plot, ".$"),
         Site=str_extract(.$Plot, "^..")) %>%
  filter(PlotNum<5) %>%
  mutate(SiteType=dplyr::recode(.$SiteType, H="High Elevation", 
                               M="Mid Elevation",
                               R="Riparian"))

#plots average soil temp by site type
AvgTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=SiteType, y = avg_temp, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("Average Soil Temperature (C)") +
  ggtitle("Average Temp") +
  theme_test() +
  theme(legend.position="none")

#plots average maximum daily soil temp by site type
DiurnalTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=SiteType, y = DiurnalVar, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Maximum Daily Temp") +
  theme_test() +
  theme(legend.position="none")

#plots average maximum daily soil temp by site type
MaxTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=SiteType, y = max_temp, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Maximum Daily Temp") +
  theme_test() +
  theme(legend.position="none")

#plots average minimum daily soil temp by site type
MinTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=SiteType, y = min_temp, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Minimum Daily Temp") +
  theme_test() +
  theme(legend.position="none") 
######################################Moisture Plots#####################
DiurnalMoistPlot=
  ggplot(Diurnal_Moisture_Summary, aes(x=SiteType, y = DiurnalVar, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Diurnal Variability in Moisture") +
  theme_test() +
  theme(legend.position="none")

AvgMoistPlot=
  Diurnal_Moisture_Summary %>%
  mutate(Site=factor(Site, level=c("MC", "MB","RB","RC","HA","HB","HC","RA","MA"))) %>%
  ggplot(aes(x=Site, y = avg_moist, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("Average Volumetric Moisture") +
  ggtitle("Average Moist") +
  theme_test() +
  theme(legend.position="none")


MaxMoistPlot=
  ggplot(Diurnal_Moisture_Summary, aes(x=SiteType, y = max_moist, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Maximum Daily Moisture") +
  theme_test() +
  theme(legend.position="none")


MinMoistPlot=
  ggplot(Diurnal_Moisture_Summary, aes(x=SiteType, y = min_moist, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Minimum Daily Moisture") +
  theme_test() +
  theme(legend.position="none") 

#makes multipanel plots showing moist and temp parameters
plot_grid(AvgTempPlot, MaxTempPlot, MinTempPlot, DiurnalTempPlot)
plot_grid(AvgMoistPlot, MaxMoistPlot, MinMoistPlot, DiurnalMoistPlot)
