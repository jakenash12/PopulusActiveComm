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

path="C:/Users/akeja/OneDrive - Duke University/Documents/Duke PhD/Projects/PMI/AUE_2021/Teros11DataDownload"
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
Fulldf %<>%
  filter(Time<"2021-10-11", Time >"2021-05-19") %>% #subsets data by a start and end date 
  mutate(value=as.numeric(value))
#subsets dataframe to only include air temperature data
AirTemp_df=
  Fulldf %>%
  filter(DataType==" °C Logger Temperature") 

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
Moistdf_wide =
  Sensordf %>%
  filter(DataType==" m³/m³ Water Content") %>% #selects only data that is of type soil moisture
  filter(!(Plot=="HB4"&Time>"2021-09-14" & Time<"2021-10-12")) %>% #filters out data from plot HB4 during a period when sensor was dug up and sitting on surface of soil
  select(Time, Plot, value) %>% #selects time, plot, and value columns
  distinct(Plot, Time, .keep_all = TRUE) %>% #removes duplicate rows
  filter(value > 0.01) %>% #removes data that is less than 1% soil moisture because it is likely incorrect (mostly just removing zeroes)
  spread(Plot, value) #converts from long to wide format, with time as column and each of the plots as 27 more columns containing soil moisture data

Moistdf_wide=Moistdf_wide[rowSums(is.na(Moistdf_wide))<27,] #removes timepoints for which all plots had NA values - these cannot be interpolated

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

ggplot(tempdf_moist, aes(x=Time)) +
  geom_point(aes(y=OldVals),color="black") +
  geom_point(aes(y=SameSite1),color="green") +
  geom_point(aes(y=SameSite2),color="blue") 

ggplot(interpolated_df_moist, aes(x=Time, y=NewVals,color=DataOrigin)) +
  geom_point() +
  facet_wrap(.~Plot) +
  theme(strip.text.x = element_blank())

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



ggplot(filter(interpolated_df_moist,Plot=="MB2", Time>"2021-10-01"), aes(x=Time, y=NewVals,color=DataOrigin)) +
  geom_point()

interpolated_df_moist %>%
  filter(Plot=="", DataOrigin=="Interpolated") %>% View

##################Data interpolation for soil temperature#################
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

ggplot(interpolated_df_soiltemp, aes(x=Time, y=NewVals,color=DataOrigin)) +
  geom_point() +
  facet_wrap(.~Plot)

SoilTempPlot=
  interpolated_df_soiltemp %>%
  mutate(Site=str_extract(Plot, "^.."),
         Date=as.Date((Time))) %>% 
  group_by(Site, Date) %>%
  dplyr::summarise(MeanDailySoilTemp=mean(NewVals)) %>% 
  mutate(SiteType=str_extract(Site,"^.")) %>%
  filter(Date>"2021-05-24", Date<"2021-10-08") %>%
  ggplot(aes(x=Date, y=MeanDailySoiLTemp, color=SiteType, group=Site)) +
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

plot_grid(SoilTempPlot, SoilMoistPlot, nrow=2)

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

#creates a dataframe with seasonal averages of soil temperature
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

#Dataframes with averages of soil temp over 10 days before each sampling date
Start_T1="2021-05-24"
End_T1="2021-06-02"

SoilTempT1=
  interpolated_df_soiltemp %>%
  filter(Time>Start_T1,Time<End_T1) %>%
  group_by(Plot) %>%
  summarise(seasonal_soiltemp=mean(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."),
         Time="June")

SoilMoistT1=
  interpolated_df_moist %>%
  filter(Time>Start_T1,Time<End_T1) %>%
  group_by(Plot) %>%
  summarise(seasonal_moist=mean(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."),
         Time="June")

Start_T2="2021-08-07"
End_T2="2021-08-17"

SoilTempT2=
  interpolated_df_soiltemp %>%
  filter(Time>Start_T2,Time<End_T2) %>%
  group_by(Plot) %>%
  summarise(seasonal_soiltemp=mean(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."),
         Time="August")

SoilMoistT2=
  interpolated_df_moist %>%
  filter(Time>Start_T2,Time<End_T2) %>%
  group_by(Plot) %>%
  summarise(seasonal_moist=mean(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."),
         Time="August")

Start_T3="2021-09-27"
End_T3="2021-10-07"

SoilTempT3=
  interpolated_df_soiltemp %>%
  filter(Time>Start_T3,Time<End_T3) %>%
  group_by(Plot) %>%
  summarise(seasonal_soiltemp=mean(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."),
         Time="October")

SoilMoistT3=
  interpolated_df_moist %>%
  filter(Time>Start_T3,Time<End_T3) %>%
  group_by(Plot) %>%
  summarise(seasonal_moist=mean(NewVals, na.rm=TRUE)) %>%
  mutate(Site=str_extract(Plot, "^.."),
         SiteType=str_extract(Plot, "^."),
         Time="October")

SeasonalSoilTempAvg=rbind(SoilTempT1, SoilTempT2, SoilTempT3)
SeasonalSoilMoistAvg=rbind(SoilMoistT1, SoilMoistT2, SoilMoistT3)

######################Processing of air temp data#####################

#filters out air temperature data that is within growing season start/end dates
AirTemp_df=
  Fulldf %>%
  filter(DataType==" °C Logger Temperature",
         Time>StartDate,Time<EndDate)

#averages over the whole season are calculated, as well as standard deviations
#First, a sd on raw data (which includes variation induced by diel cycles) is calculated
#Then, an sd is calculated on a dataframe of daily averages, to focus
#on inter-day variation (i.e. seasonality)
AirTemp_avg_df_temp=
  AirTemp_df %>%
  group_by(Site) %>% 
  summarise(avg_airtemp=mean(value),
            sd_airtemp=sd(value)) 

AirTemp_avg_df =
  AirTemp_df %>%
  mutate(Day=as.Date(Time)) %>%
  group_by(Site, Day) %>%
  summarise(dailyavg_airtemp=mean(value)) %>%
  group_by(Site) %>%
  summarise(sd_airtemp_daily_avg=sd(dailyavg_airtemp)) %>%
  left_join(AirTemp_avg_df_temp)

AirTempT1=
  AirTemp_df %>%
  filter(Time>Start_T1,Time<End_T1) %>%
  group_by(Site) %>%
  summarise(seasonal_airtemp=mean(value, na.rm=TRUE)) %>%
  mutate(Time="June")

AirTempT2=
  AirTemp_df %>%
  filter(Time>Start_T2,Time<End_T2) %>%
  group_by(Site) %>%
  summarise(seasonal_airtemp=mean(value, na.rm=TRUE)) %>%
  mutate(Time="August")

AirTempT3=
  AirTemp_df %>%
  filter(Time>Start_T3,Time<End_T3) %>%
  group_by(Site) %>%
  summarise(seasonal_airtemp=mean(value, na.rm=TRUE)) %>%
  mutate(Time="October")

SeasonalAirTempAvg=rbind(AirTempT1, AirTempT2, AirTempT3)

test=
  SeasonalSoilMoistAvg %>%
  spread(Time,seasonal_moist)

ggplot(test, aes(x=August, y =October)) +
  geom_point()
  
Matric=
  Fulldf %>%
  filter(DataType==" kPa Matric Potential") %>%
  mutate(value=as.numeric(value))

WaterDF =
  Fulldf %>%
  filter((DataType==" m³/m³ Water Content" | 
           DataType==" kPa Matric Potential")) %>%
  mutate(value=as.numeric(value)) %>%
  filter((Site=="RA" |
           Site=="MA" |
           Site=="HA" |
           Site=="HC"))

WaterPotentialDF =
  WaterDF %>%
  filter(DataType==" kPa Matric Potential") %>%
  rename(MatricPotential=value)

VolWaterDf=
  WaterDF %>%
  filter(DataType==" m³/m³ Water Content") %>%
  rename(VolWater=value) %>%
  filter((Plot=="RA3" |
            Plot=="MA3" |
            Plot=="HA3" |
            Plot=="HC3"))

FullwaterDF =
  WaterPotentialDF %>%
  left_join(VolWaterDf, by=c("Time", "Site"))

ggplot(FullwaterDF, aes(x=VolWater, y = MatricPotential)) +
  geom_point(aes(color=Site))

Fulldf %>%
  filter(DataType==" m³/m³ Water Content") %>%
  mutate(value=as.numeric(value)) %>%
  mutate(Drought=value<0.15) %>%
  filter(value!=0, value<0.5) %>%
  ggplot(aes(x=Time, y=value, color=Drought, group=Plot)) +
  geom_line(size=1) +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-06-04")), col = "black") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-08-18")), col = "black") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-10-08")), col = "black") +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Month") +
  facet_wrap(.~Plot,ncol=3) +
  scale_color_manual(values=c("black", "red")) +
  theme_bw() +
  theme(legend.position="none") +
  theme(text=element_text(size=16)) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()) 

Fulldf_wide=
  Fulldf %>%
  filter(DataType==" m³/m³ Water Content") %>%
  filter(value!=0, value<0.5) %>%
  select(Time, Plot, value) %>% 
  distinct(Time, Plot, .keep_all = TRUE) %>%
  pivot_wider(id_cols=Time,names_from=Plot, values_from=value)
  
Filldf_gathered=
  Fulldf_wide %>%
  gather(Plot, value, -Time) %>%
  mutate(Site=str_extract(Plot, "^.."),
         Missing=is.na(value),
         value=as.numeric(value),
         PlotID=str_extract(Plot,".$")) %>% 
  filter(PlotID=="2" |
           PlotID=="3" |
           PlotID=="4")

ggplot(Filldf_gathered) +
  geom_point(aes(x=Time, y=0.1, color=Missing)) +
  facet_grid(Site~PlotID) 

Fulldf_w_missing=
  data.frame(Time=unique(Fulldf$Time)) %>%
  

Fulldf %>%
  filter(DataType==" m³/m³ Water Content") %>%
  mutate(value=as.numeric(value)) %>%
  mutate(Drought=value<0.15) %>%
  filter(value!=0, value<0.5) %>%
  ggplot(aes(x=Time, y=value, group=Plot)) +
  geom_line(size=1) +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Day") +
  facet_wrap(.~Plot, ncol=3) 

Fulldf %>%
  filter(DataType==" m³/m³ Water Content", Plot=="HA3") %>%
  mutate(value=as.numeric(value)) %>%
  mutate(Drought=value<0.15) %>%
  filter(value!=0, value<0.5) %>%
  ggplot(aes(x=Time, y=value, group=Plot)) +
  geom_line(size=1) +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Day") +
  facet_wrap(.~Plot, ncol=3) +
  scale_color_manual(values=c("black", "red")) +
  theme(legend.position="none",
        text = element_text(size=18))

Fulldf %>%
  filter(DataType==" m³/m³ Water Content", Plot=="HA3") %>%
  mutate(value=as.numeric(value)) %>%
  mutate(Drought=value<0.15) %>%
  filter(value!=0, value<0.5) %>%
  ggplot(aes(x=Time, y=value, group=Plot)) +
  geom_line(size=1) +
  #geom_vline(xintercept = as.POSIXct(as.Date("2021-08-18")), 
             #color = "blue", 
             #lwd = 1,
             #linetype=2) +
  #geom_vline(xintercept = as.POSIXct(as.Date("2021-06-03")), 
             #color = "blue", 
             #lwd = 1,
             #linetype=2) +
  #geom_hline(yintercept=0.15, color="red", size=1) +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Month") +
  scale_color_manual(values=c("black", "red")) +
  theme(legend.position="none",
        text = element_text(size=18))

Fulldf %>%
  filter(DataType==" °C Soil Temperature") %>%
  mutate(value=as.numeric(value),
         Day=as.Date(Time)) %>%
  filter(value>0, Day>"2021-05-20") %>% 
  group_by(SiteType, Day) %>% 
  summarise(value=mean(value)) %>%
  ggplot(aes(x=Day, y=value, color=SiteType)) +
  theme_test() +
  geom_line(size=1) +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-08-18")), 
             color = "blue", 
             lwd = 1,
             linetype=2) +
  geom_vline(xintercept = as.POSIXct(as.Date("2021-06-03")), 
             color = "blue", 
             lwd = 1,
             linetype=2) +
  ylab("Soil Temperature (C)") +
  xlab("Month") +
  theme(text = element_text(size=18))
  
Fulldf %>%
  filter(DataType==" m³/m³ Water Content") %>%
  mutate(value=as.numeric(value),
         Day=as.Date(Time)) %>%
  filter(Day>"2021-05-20"&Day<"2021-10-06") %>% 
  group_by(SiteType, Day) %>% 
  summarise(value=mean(value)) %>%
  ggplot(aes(x=Day, y=value, color=SiteType)) +
  theme_test() +
  geom_line(size=1) +
  ylab("Soil Moisture") +
  xlab("Month") +
  theme(text = element_text(size=18))

Fulldf %>%
  filter(DataType==" m³/m³ Water Content") %>%
  mutate(value=as.numeric(value)) %>%
  filter(value!=0, value<0.4, Time >=as.Date("2021-05-27"), Time<=("2021-10-09")) %>%
  group_by(SiteType, Time) %>%
  summarise(AvgTemp=mean(value)) %>%
  ungroup %>%
  ggplot(aes(x=Time, y=AvgTemp, color=SiteType)) +
  geom_line(size=1) +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-06-04")), col = "red", size=1,linetype="dashed") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-08-19")), col = "red", size=1,linetype="dashed") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-10-08")), col = "red", size=1,linetype="dashed") +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Day") +
  theme_test()

Fulldf %>%
  filter(DataType==" m³/m³ Water Content") %>%
  mutate(value=as.numeric(value)) %>%
  filter(value!=0, value<0.4, Time >=as.Date("2021-05-27"), Time<=("2021-10-09")) %>%
  group_by(Site, Time) %>%
  summarise(AvgMoist=mean(value)) %>%
  ungroup %>%
  mutate(SiteType=str_extract(.$Site, "^.")) %>%
  ggplot(aes(x=Time, y=AvgMoist, color=SiteType, group=Site)) +
  geom_line(size=1) +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-06-04")), col = "red", size=1, linetype = "longdash") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-08-19")), col = "red", size=1, linetype = "longdash") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-10-08")), col = "red", size=1, linetype = "longdash") +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Date") +
  theme_classic() 
  theme(legend.position="none")

Fulldf %>%
  filter(DataType==" °C Soil Temperature") %>%
  mutate(value=as.numeric(value),
         Date=as.Date(Time)) %>%
  filter(Time >=as.Date("2021-05-27"), Time<=("2021-10-09")) %>%
  group_by(Site, Date) %>%
  summarise(AvgTemp=mean(value)) %>%
  ungroup %>%
  mutate(SiteType=str_extract(.$Site, "^.")) %>%
  ggplot(aes(x=Date, y=AvgTemp, color=SiteType, group=Site)) +
  geom_line(size=1) +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-06-04")), col = "black") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-08-19")), col = "black") +
  geom_vline(xintercept = as.integer(as.POSIXct("2021-10-08")), col = "black") +
  ylab("Volumetric Soil Moisture (m3/m3)") +
  xlab("Day") +
  theme_classic()

Fulldf %>%
  filter(DataType==" °C Soil Temperature") %>%
  mutate(value=as.numeric(value)) %>%
  filter(value!=0) %>%
  group_by(Site, Time) %>%
  summarise(AvgTemp=mean(value)) %>%
  ungroup %>%
  mutate(SiteType=str_extract(.$Site, "^.")) %>%
  ggplot(aes(x=Time, y=AvgTemp)) +
  geom_line(size=1) +
  ylab("Temp (C)") +
  xlab("Month") +
  theme_classic() 
  facet_wrap(.~Plot, ncol=3)

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

AirTempSummary =
  AirTemp_df %>%
  mutate(value=as.numeric(value),
         Date = as.Date(Time)) %>%
  group_by(Site, Date) %>%
  summarise(avg_temp=mean(value),
            min_temp=min(value),
            max_temp=max(value),
            DiurnalVar=max(value)-min(value)) %>%
  filter(avg_temp!=0) %>%
  ungroup %>%
  group_by(Site) %>%
  summarise(avg_airtemp=mean(avg_temp),
            min_airtemp=mean(min_temp),
            max_airtemp=mean(max_temp),
            DiurnalVar_airtemp=mean(DiurnalVar)) %>%
  mutate(SiteType=str_extract(.$Site, "^.")) %>%
  mutate(SiteType=dplyr::recode(.$SiteType, H="High Elevation", 
                                M="Mid Elevation",
                                R="Riparian"))

mod=lm(avg_moist~Site, Diurnal_Moisture_Summary)
Anova(mod)  


pairwise.t.test(Diurnal_Moisture_Summary$avg_moist,Diurnal_Moisture_Summary$SiteType)

AvgTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=Site, y = avg_temp, colour=SiteType)) +
  geom_boxplot(size=1) +
  ggtitle("Average Soil Temp") +
  ylab("Temperature (C)") +
  theme_test() +
  theme(legend.position="none")

ggplot(AirTempSummary, aes(x=Site, y=avg_temp, color=SiteType)) +
  geom_point(aes(size=3)) +  
  theme_test() +
  theme(legend.position="none")

AvgTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=SiteType, y = avg_temp, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("Average Soil Temperature (C)") +
  ggtitle("Average Temp") +
  theme_test() +
  theme(legend.position="none")


MaxTempPlot=
  ggplot(Diurnal_Temp_Summary, aes(x=SiteType, y = max_temp, colour=SiteType)) +
  geom_boxplot(size=1) +
  ylab("") +
  ggtitle("Maximum Daily Temp") +
  theme_test() +
  theme(legend.position="none")


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
  ggtitle("") +
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

plot_grid(AvgTempPlot, MaxTempPlot, MinTempPlot, DiurnalTempPlot)
plot_grid(AvgMoistPlot, MaxMoistPlot, MinMoistPlot, DiurnalMoistPlot)

ggplot(WorkingDF, aes(x=(as.numeric(`Port6 kPa Matric Potential`)), y=as.numeric(`Port3 m³/m³ Water Content`))) +
  geom_point()

mod=lm(avg_temp~SiteType, Diurnal_Temp_Summary)
Anova(mod)
summary(mod)
pairwise.t.test(Diurnal_Temp_Summary$avg_temp, Diurnal_Temp_Summary$SiteType)
