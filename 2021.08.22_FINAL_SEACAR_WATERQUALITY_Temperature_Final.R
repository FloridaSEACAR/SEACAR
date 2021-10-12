#SEACAR Water Column_Quality: Temperature
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton

#Seasonal Kendall Tau for long term trends in Temperature at the MA level

#About KENDALL TAU:
#Is there a monotonic trend, and 
##what is the significance?
#tests are non-parametric (no distribution assumptions, based on signs)
#tests require a user-defined interval
#output a p-value AND 
##the direction of the trend (as tau, or T)
##T with a hat=direction and magnitude of trend
##B1 with a hat is the LINEAR rate of change
#give a rate of change as the slope
#in Seasonal Kendall Tau, different tests by month across
##years have pooled results

#####---------------------------------------------------------------
#clear the environment
rm(list = ls())
# Load required libraries
library(tidyverse)
library(dplyr)
library(tidyr)
library(gt)
library(tidymodels)

#read in data
Temperature <- read.csv("~/Downloads/Water Temperature and pH-2021-Jul-26.csv")
#correct the MA names
Temperature$ManagedAreaName<-ifelse(Temperature$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                                    ifelse(Temperature$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                           ifelse(Temperature$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                                  ifelse(Temperature$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                         ifelse(Temperature$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                                paste(Temperature$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.	
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(Temperature)

#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
colnames(Temperature)[colnames(Temperature) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(Temperature)[colnames(Temperature) == "X.WaterTemperature_DegreesC."] <- "Temperature"

unique(Temperature$Temperature)
class(Temperature$Temperature)
nrow(Temperature)

#remove NAs in Temperature
Temperature<-Temperature[!is.na(Temperature$Temperature),]
#check the minimum and maximum values
min(Temperature$Temperature)
max(Temperature$Temperature)
#there are negative values. They are wrong. Make them NAs.
Temperature$Temperature[Temperature$Temperature<0] <- NA
#remove NAs again in Temperature
Temperature<-Temperature[!is.na(Temperature$Temperature),]
#check the minimum and maximum values again
min(Temperature$Temperature)
max(Temperature$Temperature)
########################
#Prepare date to be read in properly
unique(Temperature$Month)
unique(Temperature$Year)

#We need date to be read, by R, as a date. See what class it is.
class(Temperature$sampledate)
unique(Temperature$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
Temperature$sampledate2<-as.character(Temperature$sampledate)
class(Temperature$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
Temperature$sampledate <- as.Date(Temperature$sampledate2, format="%m/%d/%y")
###works to here with old file
class(Temperature$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(Temperature)
#remove NAs
#make a new df (Temperature2), and remove any remaining Temperature NAs. This way we can always look at the 
#original without reloading files.
Temperature2<-Temperature[!is.na(Temperature$Temperature),]
colnames(Temperature2)
nrow(Temperature2)

#Remove date NAs
Temperature2<-Temperature2[!is.na(Temperature2$Year),]
Temperature2<-Temperature2[!is.na(Temperature2$Month),]
Temperature2<-Temperature2[!is.na(Temperature2$sampledate),]
noNAN<-nrow(Temperature2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(Temperature2$Value_Qualifier_STORET_WIN)
nrow(Temperature2)
Temperature2$Qinv<-grepl("Q",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Hinv<-grepl("H",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Yinv<-grepl("Y",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Jinv<-grepl("J",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Ninv<-grepl("N",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Dinv<-grepl("D",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Ginv<-grepl("G",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$Tinv<-grepl("T",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#Temperature2$Uinv<-grepl("U",Temperature2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Temperature2$STORETWIN_remove<-
  ifelse(Temperature2$Qinv=="TRUE","remove",
         ifelse(Temperature2$Hinv=="TRUE","remove",
                ifelse(Temperature2$Yinv=="TRUE","remove",
                       ifelse(Temperature2$Jinv=="TRUE","remove",
                              ifelse(Temperature2$Ninv=="TRUE","remove",
                                     ifelse(Temperature2$Dinv=="TRUE","remove",
                                            ifelse(Temperature2$Ginv=="TRUE","remove",
                                                   ifelse(Temperature2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(Temperature2)
Temperature2<-Temperature2%>%filter(STORETWIN_remove=="KEEP")
nrow(Temperature2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(Temperature2$Value_Qualifier_NERR)
Temperature2$Ainv2<-grepl("A",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Dinv2<-grepl("D",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Hinv2<-grepl("H",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Kinv2<-grepl("K",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Minv2<-grepl("M",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Ninv2<-grepl("N",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Uinv2<-grepl("U",Temperature2$Value_Qualifier_NERR,fixed=TRUE)
Temperature2$Sinv2<-grepl("S",Temperature2$Value_Qualifier_NERR,fixed=TRUE)

Temperature2$NERR_remove<-
  ifelse(Temperature2$Ainv2=="TRUE","remove",
         ifelse(Temperature2$Dinv2=="TRUE","remove",
                ifelse(Temperature2$Hinv2=="TRUE","remove",
                       ifelse(Temperature2$Kinv2=="TRUE","remove",
                              ifelse(Temperature2$Minv2=="TRUE","remove",
                                     ifelse(Temperature2$Ninv2=="TRUE","remove",
                                            ifelse(Temperature2$Uinv2=="TRUE","remove",
                                                   ifelse(Temperature2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(Temperature2)
Temperature2<-Temperature2%>%filter(NERR_remove=="KEEP")
nrow(Temperature2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(Temperature2$QAQCFlag)
unique(Temperature2$QAQCFlag)
Temperature2$flag5<-grepl("<-5>",Temperature2$QAQCFlag,fixed=TRUE)
Temperature2$flag4<-grepl("<-4>",Temperature2$QAQCFlag,fixed=TRUE)
Temperature2$flag3<-grepl("<-3>",Temperature2$QAQCFlag,fixed=TRUE)
Temperature2$flag2<-grepl("<-2>",Temperature2$QAQCFlag,fixed=TRUE)
Temperature2$flag1<-grepl("<1>",Temperature2$QAQCFlag,fixed=TRUE)


Temperature2$FLAG_remove<-
  ifelse(Temperature2$flag5=="TRUE","remove",
         ifelse(Temperature2$flag4=="TRUE","remove",
                ifelse(Temperature2$flag3=="TRUE","remove",
                       ifelse(Temperature2$flag2=="TRUE","remove",
                              ifelse(Temperature2$flag1=="TRUE","remove","KEEP")))))

nrow(Temperature2)
Temperature2<-Temperature2%>%filter(FLAG_remove=="KEEP")
nrow(Temperature2)

#remove units outside of range
min(Temperature2$Temperature)
max(Temperature2$Temperature)
a<-nrow(Temperature2)
Temperature2<-Temperature2%>%filter(Temperature<=50)
Temperature2<-Temperature2%>%filter(Temperature>=-2)
b<-nrow(Temperature2)
a-b
#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
Temperature3surf<-Temperature2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
Temperature3surf$yearmonth<-ifelse(Temperature3surf$Month==10|Temperature3surf$Month==11|Temperature3surf$Month==12, (paste(Temperature3surf$Year, Temperature3surf$Month, sep=".")),
                                   (paste(Temperature3surf$Year, Temperature3surf$Month, sep=".0")))
colnames(Temperature3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(Temperature ~ yearmonth, data = Temperature3surf, xlab = "yearmonth",
        ylab = "Temperature, degrees C", main = "Temperature boxplot", notch=FALSE)
#then at year level
boxplot(Temperature ~ Year, data = Temperature3surf, xlab = "Year",
        ylab = "Temperature, degrees C", main = "Temperature boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(Temperature3surf$ManagedAreaName), function (i) {
  dat <- filter(Temperature3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = Temperature)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, y="Temperature (degrees C)") 
})
#make a pdf with all of these
pdf("Quality_Temperature_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(Temperature3surf$ManagedAreaName), function (i) {
  dat <- filter(Temperature3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = Temperature, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, y="Temperature (degrees C)") 
})
#make a pdf with all of these
pdf("Quality_Temperature_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average Temperature by month, year, depth, and area. We will use this in analyses.
colnames(Temperature2)
Temperatureavg<-
  Temperature2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(Temperature))
#rename columns
colnames(Temperatureavg) 
colnames(Temperatureavg) [5]<-"Temperature"

#Make a new df that will let us count years
Temperatureavgyrct<-
  Temperature2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(Temperature))
#rename columns
colnames(Temperatureavgyrct) 
colnames(Temperatureavgyrct) [5]<-"Temperature"
#################
#################
#Run analyses for Surface Temperature:
#Make a new df with only Surface Temperatures so we can run analyses for those first 
#we have already removed NAs and qualifier codes from Temperatureavg and Temperatureavgyrct
colnames(Temperatureavg)
SurfaceTemperatureavg<-
  Temperatureavg%>%
  filter(Relative_Depth=="Surface")
colnames(SurfaceTemperatureavg)

SurfaceTemperatureavgyrct<-
  Temperatureavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(SurfaceTemperatureavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(SurfaceTemperatureavg$Month)
class(SurfaceTemperatureavg$Year)
SurfaceTemperatureavg$mth<-as.numeric(SurfaceTemperatureavg$Month)
SurfaceTemperatureavg$yr<-as.numeric(SurfaceTemperatureavg$Year)
class(SurfaceTemperatureavg$mth)
class(SurfaceTemperatureavg$Temperature)

class(SurfaceTemperatureavgyrct$Month)
class(SurfaceTemperatureavgyrct$Year)
SurfaceTemperatureavgyrct$mth<-as.numeric(SurfaceTemperatureavgyrct$Month)
class(SurfaceTemperatureavgyrct$mth)
class(SurfaceTemperatureavgyrct$Year)
SurfaceTemperatureavgyrct$yr<-as.numeric(SurfaceTemperatureavgyrct$Year)
class(SurfaceTemperatureavgyrct$yr)
class(SurfaceTemperatureavgyrct$Temperature)
colnames(SurfaceTemperatureavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(SurfaceTemperatureavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  SurfaceTemperatureavg%>%
  group_by(ManagedAreaName)%>%
  summarise(length(unique(Year)))
colnames(nyearsS) [2]<-"numyrs"
eliminateS<-
  nyearsS%>%
  filter(numyrs<10)
eliminateS
colnames(eliminateS)
class(eliminateS)
#use anti_join to include only programs>=10 yrs of data
SurfaceTemperatureavg<-anti_join(SurfaceTemperatureavg,eliminateS,by=c("ManagedAreaName"))
unique(SurfaceTemperatureavg$ManagedAreaName)
SurfaceTemperatureavgyrct<-SurfaceTemperatureavg

colnames(SurfaceTemperatureavg)
colnames(SurfaceTemperatureavgyrct)

# PLOTTING
colnames(SurfaceTemperatureavgyrct)
SurfaceTemperatureavgyrct$yearmonth<-ifelse(SurfaceTemperatureavgyrct$mth==10|SurfaceTemperatureavgyrct$mth==11|SurfaceTemperatureavgyrct$mth==12, (paste(SurfaceTemperatureavgyrct$yr, SurfaceTemperatureavgyrct$mth, sep=".")),
                                            (paste(SurfaceTemperatureavgyrct$yr, SurfaceTemperatureavgyrct$mth, sep=".0")))

unique(SurfaceTemperatureavgyrct$yearmonth)
###
plots_list <- lapply(unique(SurfaceTemperatureavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(SurfaceTemperatureavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = Temperature)) +
    geom_point() +
    labs(title = i, y="Temperature (degrees C)") +
    scale_x_discrete(breaks = SurfaceTemperatureavgyrct$yearmonth[seq(1, length(SurfaceTemperatureavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Quality_Temperature_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(SurfaceTemperatureavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(Temperature~mth+yr, data = subset(SurfaceTemperatureavgyrct, ManagedAreaName == i))

#First, view the vector
unique(SurfaceTemperatureavgyrct$ManagedAreaName)
colnames(SurfaceTemperatureavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(SurfaceTemperatureavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(Temperature~mth + yr, data = subset(SurfaceTemperatureavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(SurfaceTemperatureavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Quality_Temperature_MA_surf_KToutput.csv")
kendall_ManagedAreaNameal_list_Surface
sink()
#

#********************************************************
#GENERATE RESULTS TABLES

#first, results table (after seasonal test)
kendall_map_surface <- map(.x = kendall_ManagedAreaNameal_list_Surface, .f = broom::tidy)

### excel export version
formatted_kendall_surface <- enframe(kendall_map_surface) %>% 
  unnest() %>% 
  rename(`tau statistic` = estimate1) %>% 
  rename(`senn slope` = estimate2) %>% 
  rename(intercept = estimate3) %>% #put the z-trend (Kendal tau) p-value then chi-square stat and p-value 
  group_by(name) %>% 
  mutate(p_value_stat = case_when(
    row_number() == 1 ~ "Chi_sq",
    row_number() == 2 ~ "z_(trend)"
  )) %>% 
  mutate(Trend_direction = case_when(
    `tau statistic` > 0  ~ "Increase",
    `tau statistic` == 0  ~ "No trend",
    `tau statistic` < 0  ~ "Decrease",
  )) %>% 
  ungroup() %>% 
  pivot_wider(names_from = p_value_stat,
              values_from = c(p.value, statistic)) %>% 
  rename(`Managed Area` = name) %>%
  rename(`z statistic p-value` = `p.value_z_(trend)`) %>% 
  rename(`chi-sq` = `statistic_Chi_sq`) %>% 
  rename(`chi-sq p-value` = `p.value_Chi_sq`) %>% 
  rename(`z statistic` = `statistic_z_(trend)`) %>% 
  rename(`Trend Direction` = Trend_direction) %>% 
  select(`Managed Area`, `tau statistic`, `senn slope`,
         intercept, `z statistic`, `z statistic p-value`,
         `chi-sq`, `chi-sq p-value`, `Trend Direction`)


kendall_with_regions_surface <- left_join(formatted_kendall_surface, regions, by = c("Managed Area" = "ManagedAreaName"))

kendall_with_regions_final_surface <- kendall_with_regions_surface %>% 
  select(-X) %>% 
  relocate(Region, .before = `Managed Area`)

kendall_with_regions_final_surface %>%
  gt() 

colnames(kendall_with_regions_final_surface)
colnames(kendall_with_regions_final_surface) [7] <- "zpvalue" 
kendall_with_regions_final_surface$significance<-ifelse(kendall_with_regions_final_surface$zpvalue<=0.05,"significant","not significant")
colnames(kendall_with_regions_final_surface) [7] <- "z statistic p-value" 
colnames(kendall_with_regions_final_surface) [10] <- "tau statistic trend direction" 
kendall_with_regions_final_surface<-kendall_with_regions_final_surface %>% select(-11) 
write_csv(kendall_with_regions_final_surface, "Quality_Temperature_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(Temperature3surf)
SUMMARYTemperaturesurfmeanmedstdev<-
  Temperature3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(Temperature),max=max(Temperature),median=median(Temperature),mean=mean(Temperature),sd=sd(Temperature),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYTemperaturesurfmeanmedstdev)

colnames(SUMMARYTemperaturesurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [4] <- "Minimum Value (Temperature, degrees C)" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [5] <- "Maximum Value (Temperature, degrees C)" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [6] <- "Median (Temperature, degrees C)" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [7] <- "Mean (Temperature, degrees C)" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [8] <- "Standard Deviation (Temperature, degrees C)" 
colnames(SUMMARYTemperaturesurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYTemperaturesurfmeanmedstdev,"Quality_Temperature_MA_surf_resultssummarystats2.csv")
