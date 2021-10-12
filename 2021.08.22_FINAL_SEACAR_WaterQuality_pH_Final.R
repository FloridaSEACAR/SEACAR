#SEACAR Water Column_Quality: pH
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton'

#Seasonal Kendall Tau for long term trends in pH at the MA level

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
pH <- read.csv("~/Downloads/Water Temperature and pH-2021-Jul-26.csv")
#correct the MA names
pH$ManagedAreaName<-ifelse(pH$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                           ifelse(pH$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                  ifelse(pH$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                         ifelse(pH$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                ifelse(pH$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                       paste(pH$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.	
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(pH)

#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
colnames(pH)[colnames(pH) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(pH)[colnames(pH) == "pH"] <- "pH"

unique(pH$pH)
class(pH$pH)
nrow(pH)

#remove NAs in pH
pH<-pH[!is.na(pH$pH),]
#check the minimum and maximum values
min(pH$pH)
max(pH$pH)
#there are negative values. They are wrong. Make them NAs.
pH$pH[pH$pH<0] <- NA
#remove NAs again in pH
pH<-pH[!is.na(pH$pH),]
#check the minimum and maximum values again
min(pH$pH)
max(pH$pH)
########################
#Prepare date to be read in properly
unique(pH$Month)
unique(pH$Year)

#We need date to be read, by R, as a date. See what class it is.
class(pH$sampledate)
unique(pH$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
pH$sampledate2<-as.character(pH$sampledate)
class(pH$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
pH$sampledate <- as.Date(pH$sampledate2, format="%m/%d/%y")
###works to here with old file
class(pH$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(pH)
#remove NAs
#make a new df (pH2), and remove any remaining pH NAs. This way we can always look at the 
#original without reloading files.
pH2<-pH[!is.na(pH$pH),]
colnames(pH2)
nrow(pH2)

#Remove date NAs
pH2<-pH2[!is.na(pH2$Year),]
pH2<-pH2[!is.na(pH2$Month),]
pH2<-pH2[!is.na(pH2$sampledate),]
noNAN<-nrow(pH2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(pH2$Value_Qualifier_STORET_WIN)
nrow(pH2)
pH2$Qinv<-grepl("Q",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Hinv<-grepl("H",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Yinv<-grepl("Y",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Jinv<-grepl("J",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Ninv<-grepl("N",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Dinv<-grepl("D",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Ginv<-grepl("G",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$Tinv<-grepl("T",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#pH2$Uinv<-grepl("U",pH2$Value_Qualifier_STORET_WIN,fixed=TRUE)
pH2$STORETWIN_remove<-
  ifelse(pH2$Qinv=="TRUE","remove",
         ifelse(pH2$Hinv=="TRUE","remove",
                ifelse(pH2$Yinv=="TRUE","remove",
                       ifelse(pH2$Jinv=="TRUE","remove",
                              ifelse(pH2$Ninv=="TRUE","remove",
                                     ifelse(pH2$Dinv=="TRUE","remove",
                                            ifelse(pH2$Ginv=="TRUE","remove",
                                                   ifelse(pH2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(pH2)
pH2<-pH2%>%filter(STORETWIN_remove=="KEEP")
nrow(pH2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(pH2$Value_Qualifier_NERR)
pH2$Ainv2<-grepl("A",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Dinv2<-grepl("D",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Hinv2<-grepl("H",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Kinv2<-grepl("K",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Minv2<-grepl("M",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Ninv2<-grepl("N",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Uinv2<-grepl("U",pH2$Value_Qualifier_NERR,fixed=TRUE)
pH2$Sinv2<-grepl("S",pH2$Value_Qualifier_NERR,fixed=TRUE)

pH2$NERR_remove<-
  ifelse(pH2$Ainv2=="TRUE","remove",
         ifelse(pH2$Dinv2=="TRUE","remove",
                ifelse(pH2$Hinv2=="TRUE","remove",
                       ifelse(pH2$Kinv2=="TRUE","remove",
                              ifelse(pH2$Minv2=="TRUE","remove",
                                     ifelse(pH2$Ninv2=="TRUE","remove",
                                            ifelse(pH2$Uinv2=="TRUE","remove",
                                                   ifelse(pH2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(pH2)
pH2<-pH2%>%filter(NERR_remove=="KEEP")
nrow(pH2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(pH2$QAQCFlag)
unique(pH2$QAQCFlag)
pH2$flag5<-grepl("<-5>",pH2$QAQCFlag,fixed=TRUE)
pH2$flag4<-grepl("<-4>",pH2$QAQCFlag,fixed=TRUE)
pH2$flag3<-grepl("<-3>",pH2$QAQCFlag,fixed=TRUE)
pH2$flag2<-grepl("<-2>",pH2$QAQCFlag,fixed=TRUE)
pH2$flag1<-grepl("<1>",pH2$QAQCFlag,fixed=TRUE)


pH2$FLAG_remove<-
  ifelse(pH2$flag5=="TRUE","remove",
         ifelse(pH2$flag4=="TRUE","remove",
                ifelse(pH2$flag3=="TRUE","remove",
                       ifelse(pH2$flag2=="TRUE","remove",
                              ifelse(pH2$flag1=="TRUE","remove","KEEP")))))

nrow(pH2)
pH2<-pH2%>%filter(FLAG_remove=="KEEP")
nrow(pH2)

#remove units outside of range
min(pH2$pH)
max(pH2$pH)
a<-nrow(pH2)
pH2<-pH2%>%filter(pH<=14)
pH2<-pH2%>%filter(pH>=0)
b<-nrow(pH2)
a-b
#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
pH3surf<-pH2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
pH3surf$yearmonth<-ifelse(pH3surf$Month==10|pH3surf$Month==11|pH3surf$Month==12, (paste(pH3surf$Year, pH3surf$Month, sep=".")),
                          (paste(pH3surf$Year, pH3surf$Month, sep=".0")))
colnames(pH3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(pH ~ yearmonth, data = pH3surf, xlab = "yearmonth",
        ylab = "pH, SU", main = "pH boxplot", notch=FALSE)
#then at year level
boxplot(pH ~ Year, data = pH3surf, xlab = "Year",
        ylab = "pH, SU", main = "pH boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(pH3surf$ManagedAreaName), function (i) {
  dat <- filter(pH3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = pH)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, y="pH (SU)") 
})
#make a pdf with all of these
pdf("Quality_pH_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(pH3surf$ManagedAreaName), function (i) {
  dat <- filter(pH3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = pH, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, y="pH (SU)") 
})
#make a pdf with all of these
pdf("Quality_pH_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average pH by month, year, depth, and area. We will use this in analyses.
colnames(pH2)
pHavg<-
  pH2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(pH))
#rename columns
colnames(pHavg) 
colnames(pHavg) [5]<-"pH"

#Make a new df that will let us count years
pHavgyrct<-
  pH2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(pH))
#rename columns
colnames(pHavgyrct) 
colnames(pHavgyrct) [5]<-"pH"
#################
#################
#Run analyses for Surface pH:
#Make a new df with only Surface pHs so we can run analyses for those first 
#we have already removed NAs and qualifier codes from pHavg and pHavgyrct
colnames(pHavg)
SurfacepHavg<-
  pHavg%>%
  filter(Relative_Depth=="Surface")
colnames(SurfacepHavg)

SurfacepHavgyrct<-
  pHavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(SurfacepHavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(SurfacepHavg$Month)
class(SurfacepHavg$Year)
SurfacepHavg$mth<-as.numeric(SurfacepHavg$Month)
SurfacepHavg$yr<-as.numeric(SurfacepHavg$Year)
class(SurfacepHavg$mth)
class(SurfacepHavg$pH)

class(SurfacepHavgyrct$Month)
class(SurfacepHavgyrct$Year)
SurfacepHavgyrct$mth<-as.numeric(SurfacepHavgyrct$Month)
class(SurfacepHavgyrct$mth)
class(SurfacepHavgyrct$Year)
SurfacepHavgyrct$yr<-as.numeric(SurfacepHavgyrct$Year)
class(SurfacepHavgyrct$yr)
class(SurfacepHavgyrct$pH)
colnames(SurfacepHavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(SurfacepHavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  SurfacepHavg%>%
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
SurfacepHavg<-anti_join(SurfacepHavg,eliminateS,by=c("ManagedAreaName"))
unique(SurfacepHavg$ManagedAreaName)
SurfacepHavgyrct<-SurfacepHavg

colnames(SurfacepHavg)
colnames(SurfacepHavgyrct)

# PLOTTING
colnames(SurfacepHavgyrct)
SurfacepHavgyrct$yearmonth<-ifelse(SurfacepHavgyrct$mth==10|SurfacepHavgyrct$mth==11|SurfacepHavgyrct$mth==12, (paste(SurfacepHavgyrct$yr, SurfacepHavgyrct$mth, sep=".")),
                                   (paste(SurfacepHavgyrct$yr, SurfacepHavgyrct$mth, sep=".0")))

unique(SurfacepHavgyrct$yearmonth)
###
plots_list <- lapply(unique(SurfacepHavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(SurfacepHavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = pH)) +
    geom_point() +
    labs(title = i, y="pH (SU)") +
    scale_x_discrete(breaks = SurfacepHavgyrct$yearmonth[seq(1, length(SurfacepHavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Quality_pH_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(SurfacepHavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(pH~mth+yr, data = subset(SurfacepHavgyrct, ManagedAreaName == i))

#First, view the vector
unique(SurfacepHavgyrct$ManagedAreaName)
colnames(SurfacepHavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(SurfacepHavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(pH~mth + yr, data = subset(SurfacepHavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(SurfacepHavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Quality_pH_MA_surf_KToutput.csv")
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
write_csv(kendall_with_regions_final_surface, "Quality_pH_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(pH3surf)
SUMMARYpHsurfmeanmedstdev<-
  pH3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(pH),max=max(pH),median=median(pH),mean=mean(pH),sd=sd(pH),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYpHsurfmeanmedstdev)

colnames(SUMMARYpHsurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYpHsurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYpHsurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYpHsurfmeanmedstdev) [4] <- "Minimum Value (pH, SU)" 
colnames(SUMMARYpHsurfmeanmedstdev) [5] <- "Maximum Value (pH, SU)" 
colnames(SUMMARYpHsurfmeanmedstdev) [6] <- "Median (pH, SU)" 
colnames(SUMMARYpHsurfmeanmedstdev) [7] <- "Mean (pH, SU)" 
colnames(SUMMARYpHsurfmeanmedstdev) [8] <- "Standard Deviation (pH, SU)" 
colnames(SUMMARYpHsurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYpHsurfmeanmedstdev,"Quality_pH_MA_surf_resultssummarystats2.csv")
