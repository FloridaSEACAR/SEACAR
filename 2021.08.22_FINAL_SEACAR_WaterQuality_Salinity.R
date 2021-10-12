#SEACAR Water Column_Quality: Salinity
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton

#Seasonal Kendall Tau for long term trends in Salinity at the MA level

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
Salinity <- read.csv("~/Desktop/Salinity and Specific Conductivity-2021-Jul-26.csv")
#correct the MA names
Salinity$ManagedAreaName<-ifelse(Salinity$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                                 ifelse(Salinity$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                        ifelse(Salinity$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                               ifelse(Salinity$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                      ifelse(Salinity$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                             paste(Salinity$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.	
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(Salinity)

#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
colnames(Salinity)[colnames(Salinity) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(Salinity)[colnames(Salinity) == "Salinity_ppt"] <- "Salinity"

unique(Salinity$Salinity)
class(Salinity$Salinity)
nrow(Salinity)

#remove NAs in Salinity
Salinity<-Salinity[!is.na(Salinity$Salinity),]
#check the minimum and maximum values
min(Salinity$Salinity)
max(Salinity$Salinity)
#there are negative values. They are wrong. Make them NAs.
Salinity$Salinity[Salinity$Salinity<0] <- NA
#remove NAs again in Salinity
Salinity<-Salinity[!is.na(Salinity$Salinity),]
#check the minimum and maximum values again
min(Salinity$Salinity)
max(Salinity$Salinity)
########################
#Prepare date to be read in properly
unique(Salinity$Month)
unique(Salinity$Year)

#We need date to be read, by R, as a date. See what class it is.
class(Salinity$sampledate)
unique(Salinity$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
Salinity$sampledate2<-as.character(Salinity$sampledate)
class(Salinity$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
Salinity$sampledate <- as.Date(Salinity$sampledate2, format="%m/%d/%y")
###works to here with old file
class(Salinity$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(Salinity)
#remove NAs
#make a new df (Salinity2), and remove any remaining Salinity NAs. This way we can always look at the 
#original without reloading files.
Salinity2<-Salinity[!is.na(Salinity$Salinity),]
colnames(Salinity2)
nrow(Salinity2)

#Remove date NAs
Salinity2<-Salinity2[!is.na(Salinity2$Year),]
Salinity2<-Salinity2[!is.na(Salinity2$Month),]
Salinity2<-Salinity2[!is.na(Salinity2$sampledate),]
noNAN<-nrow(Salinity2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(Salinity2$Value_Qualifier_STORET_WIN)
nrow(Salinity2)
Salinity2$Qinv<-grepl("Q",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Hinv<-grepl("H",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Yinv<-grepl("Y",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Jinv<-grepl("J",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Ninv<-grepl("N",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Dinv<-grepl("D",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Ginv<-grepl("G",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$Tinv<-grepl("T",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#Salinity2$Uinv<-grepl("U",Salinity2$Value_Qualifier_STORET_WIN,fixed=TRUE)
Salinity2$STORETWIN_remove<-
  ifelse(Salinity2$Qinv=="TRUE","remove",
         ifelse(Salinity2$Hinv=="TRUE","remove",
                ifelse(Salinity2$Yinv=="TRUE","remove",
                       ifelse(Salinity2$Jinv=="TRUE","remove",
                              ifelse(Salinity2$Ninv=="TRUE","remove",
                                     ifelse(Salinity2$Dinv=="TRUE","remove",
                                            ifelse(Salinity2$Ginv=="TRUE","remove",
                                                   ifelse(Salinity2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(Salinity2)
Salinity2<-Salinity2%>%filter(STORETWIN_remove=="KEEP")
nrow(Salinity2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(Salinity2$Value_Qualifier_NERR)
Salinity2$Ainv2<-grepl("A",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Dinv2<-grepl("D",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Hinv2<-grepl("H",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Kinv2<-grepl("K",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Minv2<-grepl("M",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Ninv2<-grepl("N",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Uinv2<-grepl("U",Salinity2$Value_Qualifier_NERR,fixed=TRUE)
Salinity2$Sinv2<-grepl("S",Salinity2$Value_Qualifier_NERR,fixed=TRUE)

Salinity2$NERR_remove<-
  ifelse(Salinity2$Ainv2=="TRUE","remove",
         ifelse(Salinity2$Dinv2=="TRUE","remove",
                ifelse(Salinity2$Hinv2=="TRUE","remove",
                       ifelse(Salinity2$Kinv2=="TRUE","remove",
                              ifelse(Salinity2$Minv2=="TRUE","remove",
                                     ifelse(Salinity2$Ninv2=="TRUE","remove",
                                            ifelse(Salinity2$Uinv2=="TRUE","remove",
                                                   ifelse(Salinity2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(Salinity2)
Salinity2<-Salinity2%>%filter(NERR_remove=="KEEP")
nrow(Salinity2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(Salinity2$QAQCFlag)
unique(Salinity2$QAQCFlag)
Salinity2$flag5<-grepl("<-5>",Salinity2$QAQCFlag,fixed=TRUE)
Salinity2$flag4<-grepl("<-4>",Salinity2$QAQCFlag,fixed=TRUE)
Salinity2$flag3<-grepl("<-3>",Salinity2$QAQCFlag,fixed=TRUE)
Salinity2$flag2<-grepl("<-2>",Salinity2$QAQCFlag,fixed=TRUE)
Salinity2$flag1<-grepl("<1>",Salinity2$QAQCFlag,fixed=TRUE)


Salinity2$FLAG_remove<-
  ifelse(Salinity2$flag5=="TRUE","remove",
         ifelse(Salinity2$flag4=="TRUE","remove",
                ifelse(Salinity2$flag3=="TRUE","remove",
                       ifelse(Salinity2$flag2=="TRUE","remove",
                              ifelse(Salinity2$flag1=="TRUE","remove","KEEP")))))

nrow(Salinity2)
Salinity2<-Salinity2%>%filter(FLAG_remove=="KEEP")
nrow(Salinity2)

#remove units outside of range
min(Salinity2$Salinity)
max(Salinity2$Salinity)
a<-nrow(Salinity2)
Salinity2<-Salinity2%>%filter(Salinity<=70)
Salinity2<-Salinity2%>%filter(Salinity>=0)
b<-nrow(Salinity2)
a-b
#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
Salinity3surf<-Salinity2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
Salinity3surf$yearmonth<-ifelse(Salinity3surf$Month==10|Salinity3surf$Month==11|Salinity3surf$Month==12, (paste(Salinity3surf$Year, Salinity3surf$Month, sep=".")),
                                (paste(Salinity3surf$Year, Salinity3surf$Month, sep=".0")))
colnames(Salinity3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(Salinity ~ yearmonth, data = Salinity3surf, xlab = "yearmonth",
        ylab = "Salinity, PSU", main = "Salinity boxplot", notch=FALSE)
#then at year level
boxplot(Salinity ~ Year, data = Salinity3surf, xlab = "Year",
        ylab = "Salinity, PSU", main = "Salinity boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(Salinity3surf$ManagedAreaName), function (i) {
  dat <- filter(Salinity3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = Salinity)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, y="Salinity (PSU)") 
})
#make a pdf with all of these
pdf("Quality_Salinity_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(Salinity3surf$ManagedAreaName), function (i) {
  dat <- filter(Salinity3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = Salinity, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, y="Salinity (PSU)") 
})
#make a pdf with all of these
pdf("Quality_Salinity_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average Salinity by month, year, depth, and area. We will use this in analyses.
colnames(Salinity2)
Salinityavg<-
  Salinity2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(Salinity))
#rename columns
colnames(Salinityavg) 
colnames(Salinityavg) [5]<-"Salinity"

#Make a new df that will let us count years
Salinityavgyrct<-
  Salinity2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(Salinity))
#rename columns
colnames(Salinityavgyrct) 
colnames(Salinityavgyrct) [5]<-"Salinity"
#################
#################
#Run analyses for Surface Salinity:
#Make a new df with only Surface Salinitys so we can run analyses for those first 
#we have already removed NAs and qualifier codes from Salinityavg and Salinityavgyrct
colnames(Salinityavg)
SurfaceSalinityavg<-
  Salinityavg%>%
  filter(Relative_Depth=="Surface")
colnames(SurfaceSalinityavg)

SurfaceSalinityavgyrct<-
  Salinityavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(SurfaceSalinityavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(SurfaceSalinityavg$Month)
class(SurfaceSalinityavg$Year)
SurfaceSalinityavg$mth<-as.numeric(SurfaceSalinityavg$Month)
SurfaceSalinityavg$yr<-as.numeric(SurfaceSalinityavg$Year)
class(SurfaceSalinityavg$mth)
class(SurfaceSalinityavg$Salinity)

class(SurfaceSalinityavgyrct$Month)
class(SurfaceSalinityavgyrct$Year)
SurfaceSalinityavgyrct$mth<-as.numeric(SurfaceSalinityavgyrct$Month)
class(SurfaceSalinityavgyrct$mth)
class(SurfaceSalinityavgyrct$Year)
SurfaceSalinityavgyrct$yr<-as.numeric(SurfaceSalinityavgyrct$Year)
class(SurfaceSalinityavgyrct$yr)
class(SurfaceSalinityavgyrct$Salinity)
colnames(SurfaceSalinityavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(SurfaceSalinityavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  SurfaceSalinityavg%>%
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
SurfaceSalinityavg<-anti_join(SurfaceSalinityavg,eliminateS,by=c("ManagedAreaName"))
unique(SurfaceSalinityavg$ManagedAreaName)
SurfaceSalinityavgyrct<-SurfaceSalinityavg

colnames(SurfaceSalinityavg)
colnames(SurfaceSalinityavgyrct)

# PLOTTING
colnames(SurfaceSalinityavgyrct)
SurfaceSalinityavgyrct$yearmonth<-ifelse(SurfaceSalinityavgyrct$mth==10|SurfaceSalinityavgyrct$mth==11|SurfaceSalinityavgyrct$mth==12, (paste(SurfaceSalinityavgyrct$yr, SurfaceSalinityavgyrct$mth, sep=".")),
                                         (paste(SurfaceSalinityavgyrct$yr, SurfaceSalinityavgyrct$mth, sep=".0")))

unique(SurfaceSalinityavgyrct$yearmonth)
###
plots_list <- lapply(unique(SurfaceSalinityavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(SurfaceSalinityavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = Salinity)) +
    geom_point() +
    labs(title = i, y="Salinity (PSU)") +
    scale_x_discrete(breaks = SurfaceSalinityavgyrct$yearmonth[seq(1, length(SurfaceSalinityavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Quality_Salinity_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(SurfaceSalinityavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(Salinity~mth+yr, data = subset(SurfaceSalinityavgyrct, ManagedAreaName == i))

#First, view the vector
unique(SurfaceSalinityavgyrct$ManagedAreaName)
colnames(SurfaceSalinityavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(SurfaceSalinityavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(Salinity~mth + yr, data = subset(SurfaceSalinityavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(SurfaceSalinityavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Quality_Salinity_MA_surf_KToutput.csv")
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
write_csv(kendall_with_regions_final_surface, "Quality_Salinity_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(Salinity3surf)
SUMMARYSalinitysurfmeanmedstdev<-
  Salinity3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(Salinity),max=max(Salinity),median=median(Salinity),mean=mean(Salinity),sd=sd(Salinity),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYSalinitysurfmeanmedstdev)

colnames(SUMMARYSalinitysurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [4] <- "Minimum Value (Salinity PSU)" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [5] <- "Maximum Value (Salinity PSU)" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [6] <- "Median (Salinity PSU)" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [7] <- "Mean (Salinity PSU)" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [8] <- "Standard Deviation (Salinity PSU)" 
colnames(SUMMARYSalinitysurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYSalinitysurfmeanmedstdev,"Quality_Salinity_MA_surf_resultssummarystats2.csv")
