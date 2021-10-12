#SEACAR Water Column_Quality: DO
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton'

#Seasonal Kendall Tau for long term trends in DO at the MA level

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
DO <- read.csv("~/Downloads/Dissolved Oxygen-2021-Jul-26.csv")
#correct the MA names
DO$ManagedAreaName<-ifelse(DO$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                           ifelse(DO$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                  ifelse(DO$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                         ifelse(DO$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                ifelse(DO$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                       paste(DO$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.	
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(DO)

#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
colnames(DO)[colnames(DO) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(DO)[colnames(DO) == "X.DissolvedOxygen_mg.L."] <- "DO"

unique(DO$DO)
class(DO$DO)
nrow(DO)

#remove NAs in DO
DO<-DO[!is.na(DO$DO),]
#check the minimum and maximum values
min(DO$DO)
max(DO$DO)
#there are negative values. They are wrong. Make them NAs.
DO$DO[DO$DO<0] <- NA
#remove NAs again in DO
DO<-DO[!is.na(DO$DO),]
#check the minimum and maximum values again
min(DO$DO)
max(DO$DO)
########################
#Prepare date to be read in properly
unique(DO$Month)
unique(DO$Year)

#We need date to be read, by R, as a date. See what class it is.
class(DO$sampledate)
unique(DO$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
DO$sampledate2<-as.character(DO$sampledate)
class(DO$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
DO$sampledate <- as.Date(DO$sampledate2, format="%m/%d/%y")
###works to here with old file
class(DO$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(DO)
#remove NAs
#make a new df (DO2), and remove any remaining DO NAs. This way we can always look at the 
#original without reloading files.
DO2<-DO[!is.na(DO$DO),]
colnames(DO2)
nrow(DO2)

#Remove date NAs
DO2<-DO2[!is.na(DO2$Year),]
DO2<-DO2[!is.na(DO2$Month),]
DO2<-DO2[!is.na(DO2$sampledate),]
noNAN<-nrow(DO2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(DO2$Value_Qualifier_STORET_WIN)
nrow(DO2)
DO2$Qinv<-grepl("Q",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Hinv<-grepl("H",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Yinv<-grepl("Y",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Jinv<-grepl("J",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Ninv<-grepl("N",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Dinv<-grepl("D",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Ginv<-grepl("G",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$Tinv<-grepl("T",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#DO2$Uinv<-grepl("U",DO2$Value_Qualifier_STORET_WIN,fixed=TRUE)
DO2$STORETWIN_remove<-
  ifelse(DO2$Qinv=="TRUE","remove",
         ifelse(DO2$Hinv=="TRUE","remove",
                ifelse(DO2$Yinv=="TRUE","remove",
                       ifelse(DO2$Jinv=="TRUE","remove",
                              ifelse(DO2$Ninv=="TRUE","remove",
                                     ifelse(DO2$Dinv=="TRUE","remove",
                                            ifelse(DO2$Ginv=="TRUE","remove",
                                                   ifelse(DO2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(DO2)
DO2<-DO2%>%filter(STORETWIN_remove=="KEEP")
nrow(DO2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(DO2$Value_Qualifier_NERR)
DO2$Ainv2<-grepl("A",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Dinv2<-grepl("D",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Hinv2<-grepl("H",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Kinv2<-grepl("K",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Minv2<-grepl("M",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Ninv2<-grepl("N",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Uinv2<-grepl("U",DO2$Value_Qualifier_NERR,fixed=TRUE)
DO2$Sinv2<-grepl("S",DO2$Value_Qualifier_NERR,fixed=TRUE)

DO2$NERR_remove<-
  ifelse(DO2$Ainv2=="TRUE","remove",
         ifelse(DO2$Dinv2=="TRUE","remove",
                ifelse(DO2$Hinv2=="TRUE","remove",
                       ifelse(DO2$Kinv2=="TRUE","remove",
                              ifelse(DO2$Minv2=="TRUE","remove",
                                     ifelse(DO2$Ninv2=="TRUE","remove",
                                            ifelse(DO2$Uinv2=="TRUE","remove",
                                                   ifelse(DO2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(DO2)
DO2<-DO2%>%filter(NERR_remove=="KEEP")
nrow(DO2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(DO2$QAQCFlag)
unique(DO2$QAQCFlag)
DO2$flag5<-grepl("<-5>",DO2$QAQCFlag,fixed=TRUE)
DO2$flag4<-grepl("<-4>",DO2$QAQCFlag,fixed=TRUE)
DO2$flag3<-grepl("<-3>",DO2$QAQCFlag,fixed=TRUE)
DO2$flag2<-grepl("<-2>",DO2$QAQCFlag,fixed=TRUE)
DO2$flag1<-grepl("<1>",DO2$QAQCFlag,fixed=TRUE)


DO2$FLAG_remove<-
  ifelse(DO2$flag5=="TRUE","remove",
         ifelse(DO2$flag4=="TRUE","remove",
                ifelse(DO2$flag3=="TRUE","remove",
                       ifelse(DO2$flag2=="TRUE","remove",
                              ifelse(DO2$flag1=="TRUE","remove","KEEP")))))

nrow(DO2)
DO2<-DO2%>%filter(FLAG_remove=="KEEP")
nrow(DO2)

#remove units outside of range
min(DO2$DO)
max(DO2$DO)
a<-nrow(DO2)
DO2<-DO2%>%filter(DO<=22)
#DO2<-DO2%>%filter(DO>=0)
b<-nrow(DO2)
a-b
#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
DO3surf<-DO2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
DO3surf$yearmonth<-ifelse(DO3surf$Month==10|DO3surf$Month==11|DO3surf$Month==12, (paste(DO3surf$Year, DO3surf$Month, sep=".")),
                          (paste(DO3surf$Year, DO3surf$Month, sep=".0")))
colnames(DO3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(DO ~ yearmonth, data = DO3surf, xlab = "yearmonth",
        ylab = "DO, mg/L", main = "DO boxplot", notch=FALSE)
#then at year level
boxplot(DO ~ Year, data = DO3surf, xlab = "Year",
        ylab = "DO, mg/L", main = "DO boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(DO3surf$ManagedAreaName), function (i) {
  dat <- filter(DO3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = DO)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, y="DO (mg/L)") 
})
#make a pdf with all of these
pdf("Quality_DOmgL_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(DO3surf$ManagedAreaName), function (i) {
  dat <- filter(DO3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = DO, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, y="DO (mg/L)") 
})
#make a pdf with all of these
pdf("Quality_DOmgL_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average DO by month, year, depth, and area. We will use this in analyses.
colnames(DO2)
DOavg<-
  DO2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(DO))
#rename columns
colnames(DOavg) 
colnames(DOavg) [5]<-"DO"

#Make a new df that will let us count years
DOavgyrct<-
  DO2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(DO))
#rename columns
colnames(DOavgyrct) 
colnames(DOavgyrct) [5]<-"DO"
#################
#################
#Run analyses for Surface DO:
#Make a new df with only Surface DOs so we can run analyses for those first 
#we have already removed NAs and qualifier codes from DOavg and DOavgyrct
colnames(DOavg)
SurfaceDOavg<-
  DOavg%>%
  filter(Relative_Depth=="Surface")
colnames(SurfaceDOavg)

SurfaceDOavgyrct<-
  DOavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(SurfaceDOavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(SurfaceDOavg$Month)
class(SurfaceDOavg$Year)
SurfaceDOavg$mth<-as.numeric(SurfaceDOavg$Month)
SurfaceDOavg$yr<-as.numeric(SurfaceDOavg$Year)
class(SurfaceDOavg$mth)
class(SurfaceDOavg$DO)

class(SurfaceDOavgyrct$Month)
class(SurfaceDOavgyrct$Year)
SurfaceDOavgyrct$mth<-as.numeric(SurfaceDOavgyrct$Month)
class(SurfaceDOavgyrct$mth)
class(SurfaceDOavgyrct$Year)
SurfaceDOavgyrct$yr<-as.numeric(SurfaceDOavgyrct$Year)
class(SurfaceDOavgyrct$yr)
class(SurfaceDOavgyrct$DO)
colnames(SurfaceDOavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(SurfaceDOavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  SurfaceDOavg%>%
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
SurfaceDOavg<-anti_join(SurfaceDOavg,eliminateS,by=c("ManagedAreaName"))
unique(SurfaceDOavg$ManagedAreaName)
SurfaceDOavgyrct<-SurfaceDOavg

colnames(SurfaceDOavg)
colnames(SurfaceDOavgyrct)

# PLOTTING
colnames(SurfaceDOavgyrct)
SurfaceDOavgyrct$yearmonth<-ifelse(SurfaceDOavgyrct$mth==10|SurfaceDOavgyrct$mth==11|SurfaceDOavgyrct$mth==12, (paste(SurfaceDOavgyrct$yr, SurfaceDOavgyrct$mth, sep=".")),
                                   (paste(SurfaceDOavgyrct$yr, SurfaceDOavgyrct$mth, sep=".0")))

unique(SurfaceDOavgyrct$yearmonth)
###
plots_list <- lapply(unique(SurfaceDOavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(SurfaceDOavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = DO)) +
    geom_point() +
    labs(title = i, y="DO (mg/L)") +
    scale_x_discrete(breaks = SurfaceDOavgyrct$yearmonth[seq(1, length(SurfaceDOavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Quality_DOmgL_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(SurfaceDOavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(DO~mth+yr, data = subset(SurfaceDOavgyrct, ManagedAreaName == i))

#First, view the vector
unique(SurfaceDOavgyrct$ManagedAreaName)
colnames(SurfaceDOavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(SurfaceDOavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(DO~mth + yr, data = subset(SurfaceDOavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(SurfaceDOavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Quality_DOmgL_MA_surf_KToutput.csv")
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
write_csv(kendall_with_regions_final_surface, "Quality_DOmgL_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(DO3surf)
SUMMARYDOsurfmeanmedstdev<-
  DO3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(DO),max=max(DO),median=median(DO),mean=mean(DO),sd=sd(DO),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYDOsurfmeanmedstdev)

colnames(SUMMARYDOsurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYDOsurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYDOsurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYDOsurfmeanmedstdev) [4] <- "Minimum Value (DO, mg/L)" 
colnames(SUMMARYDOsurfmeanmedstdev) [5] <- "Maximum Value (DO, mg/L)" 
colnames(SUMMARYDOsurfmeanmedstdev) [6] <- "Median (DO, mg/L)" 
colnames(SUMMARYDOsurfmeanmedstdev) [7] <- "Mean (DO, mg/L)" 
colnames(SUMMARYDOsurfmeanmedstdev) [8] <- "Standard Deviation (DO, mg/L)" 
colnames(SUMMARYDOsurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYDOsurfmeanmedstdev,"Quality_DOmgL_MA_surf_resultssummarystats2.csv")


