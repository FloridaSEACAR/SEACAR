#SEACAR Water Column_Nutrient Scoring: Total Phosphorus
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton'

#Seasonal Kendall Tau for long term trends in Phoshorus at the MA level

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
phos <- read.csv("~/Desktop/Nutrients - Nitrogen and Phosphorus-2021-Jul-26.csv")
#correct the MA names
phos$ManagedAreaName<-ifelse(phos$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                            ifelse(phos$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                   ifelse(phos$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                          ifelse(phos$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                 ifelse(phos$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                        paste(phos$ManagedAreaName, "Aquatic Preserve", sep=" "))))))

unique(regions$ManagedAreaName)

#rename columns
colnames(phos) 
#read in a "regions" df to help in table-building, later.
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(phos)

#********************************************************
#1. PREPARE DATA FOR ANALYSES
#rename columns
colnames(phos)[colnames(phos) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(phos)[colnames(phos) == "X.TotalPhosphorus_mg.L."] <- "TP"

unique(phos$TP)
class(phos$TP)
nrow(phos)

#remove NAs in TP
phos<-phos[!is.na(phos$TP),]
#check the minimum and maximum values
min(phos$TP)
max(phos$TP)
#there are negative values. They are wrong. Make them NAs.
phos$TP[phos$TP<0] <- NA
#remove NAs again in TP
phos<-phos[!is.na(phos$TP),]
#check the minimum and maximum values again
min(phos$TP)
max(phos$TP)
########################
#Prepare date to be read in properly
unique(phos$Month)
unique(phos$Year)

#We need date to be read, by R, as a date. See what class it is.
class(phos$sampledate)
unique(phos$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
phos$sampledate2<-as.character(phos$sampledate)
class(phos$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
phos$sampledate <- as.Date(phos$sampledate2, format="%m/%d/%y")
###works to here with old file
class(phos$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(phos)
#remove NAs
#make a new df (phos2), and remove any remaining phos NAs. This way we can always look at the 
#original without reloading files.
phos2<-phos[!is.na(phos$TP),]
colnames(phos2)
nrow(phos2)

#Remove date NAs
phos2<-phos2[!is.na(phos2$Year),]
phos2<-phos2[!is.na(phos2$Month),]
phos2<-phos2[!is.na(phos2$sampledate),]
noNAN<-nrow(phos2)
#how many rows have we removed so far?
originalN-noNAN
#
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(phos2$Value_Qualifier_STORET_WIN)
nrow(phos2)
phos2$Qinv<-grepl("Q",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Hinv<-grepl("H",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Yinv<-grepl("Y",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Jinv<-grepl("J",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Ninv<-grepl("N",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Dinv<-grepl("D",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Ginv<-grepl("G",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$Tinv<-grepl("T",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#phos2$Uinv<-grepl("U",phos2$Value_Qualifier_STORET_WIN,fixed=TRUE)
phos2$STORETWIN_remove<-
  ifelse(phos2$Qinv=="TRUE","remove",
         ifelse(phos2$Hinv=="TRUE","remove",
                ifelse(phos2$Yinv=="TRUE","remove",
                       ifelse(phos2$Jinv=="TRUE","remove",
                              ifelse(phos2$Ninv=="TRUE","remove",
                                     ifelse(phos2$Dinv=="TRUE","remove",
                                            ifelse(phos2$Ginv=="TRUE","remove",
                                                   ifelse(phos2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(phos2)
phos2<-phos2%>%filter(STORETWIN_remove=="KEEP")
nrow(phos2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(phos2$Value_Qualifier_NERR)
phos2$Ainv2<-grepl("A",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Dinv2<-grepl("D",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Hinv2<-grepl("H",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Kinv2<-grepl("K",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Minv2<-grepl("M",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Ninv2<-grepl("N",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Uinv2<-grepl("U",phos2$Value_Qualifier_NERR,fixed=TRUE)
phos2$Sinv2<-grepl("S",phos2$Value_Qualifier_NERR,fixed=TRUE)

phos2$NERR_remove<-
  ifelse(phos2$Ainv2=="TRUE","remove",
         ifelse(phos2$Dinv2=="TRUE","remove",
                ifelse(phos2$Hinv2=="TRUE","remove",
                       ifelse(phos2$Kinv2=="TRUE","remove",
                              ifelse(phos2$Minv2=="TRUE","remove",
                                     ifelse(phos2$Ninv2=="TRUE","remove",
                                            ifelse(phos2$Uinv2=="TRUE","remove",
                                                   ifelse(phos2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(phos2)
phos2<-phos2%>%filter(NERR_remove=="KEEP")
nrow(phos2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(phos2$QAQCFlag)
unique(phos2$QAQCFlag)
phos2$flag5<-grepl("<-5>",phos2$QAQCFlag,fixed=TRUE)
phos2$flag4<-grepl("<-4>",phos2$QAQCFlag,fixed=TRUE)
phos2$flag3<-grepl("<-3>",phos2$QAQCFlag,fixed=TRUE)
phos2$flag2<-grepl("<-2>",phos2$QAQCFlag,fixed=TRUE)
phos2$flag1<-grepl("<1>",phos2$QAQCFlag,fixed=TRUE)


phos2$FLAG_remove<-
  ifelse(phos2$flag5=="TRUE","remove",
         ifelse(phos2$flag4=="TRUE","remove",
                ifelse(phos2$flag3=="TRUE","remove",
                       ifelse(phos2$flag2=="TRUE","remove",
                              ifelse(phos2$flag1=="TRUE","remove","KEEP")))))

nrow(phos2)
phos2<-phos2%>%filter(FLAG_remove=="KEEP")
nrow(phos2)

#********************************************************
#2. GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
phos3surf<-phos2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
phos3surf$yearmonth<-ifelse(phos3surf$Month==10|phos3surf$Month==11|phos3surf$Month==12, (paste(phos3surf$Year, phos3surf$Month, sep=".")),
                            (paste(phos3surf$Year, phos3surf$Month, sep=".0")))
colnames(phos3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(TP ~ yearmonth, data = phos3surf, xlab = "year.month",
        ylab = "Total Phosphorus mg/L", main = "Total Phosphorus boxplot", notch=FALSE)
#then at year level
boxplot(TP ~ Year, data = phos3surf, xlab = "Year",
        ylab = "Total Phosphorus mg/L", main = "Total Phosphorus boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(phos3surf$ManagedAreaName), function (i) {
  dat <- filter(phos3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = TP)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, x="year.month",y="Total Phosphorus mg/L") 
   
})
#make a pdf with all of these
pdf("Nutrients_P_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(phos3surf$ManagedAreaName), function (i) {
  dat <- filter(phos3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = TP, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, x="year.month",y="Total Phosphorus mg/L") 
})
#make a pdf with all of these
pdf("Nutrients_P_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#4a. BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average phos by month, year, depth, and area. We will use this in analyses.
colnames(phos2)
phosavg<-
  phos2%>%
   
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(TP))
#rename columns
colnames(phosavg) 
colnames(phosavg) [5]<-"TP"

#Make a new df that will let us count years
phosavgyrct<-
  phos2%>%
   
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(TP))
#rename columns
colnames(phosavgyrct) 
colnames(phosavgyrct) [5]<-"TP"
#################
#################
#Run analyses for Surface phos:
#Make a new df with only Surface phoss so we can run analyses for those first 
#we have already removed NAs and qualifier codes from phosavg and phosavgyrct
colnames(phosavg)
Surfacephosavg<-
  phosavg%>%
  filter(Relative_Depth=="Surface")
colnames(Surfacephosavg)

Surfacephosavgyrct<-
  phosavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(Surfacephosavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(Surfacephosavg$Month)
class(Surfacephosavg$Year)
Surfacephosavg$mth<-as.numeric(Surfacephosavg$Month)
Surfacephosavg$yr<-as.numeric(Surfacephosavg$Year)
class(Surfacephosavg$mth)
class(Surfacephosavg$TP)

class(Surfacephosavgyrct$Month)
class(Surfacephosavgyrct$Year)
Surfacephosavgyrct$mth<-as.numeric(Surfacephosavgyrct$Month)
class(Surfacephosavgyrct$mth)
class(Surfacephosavgyrct$Year)
Surfacephosavgyrct$yr<-as.numeric(Surfacephosavgyrct$Year)
class(Surfacephosavgyrct$yr)
class(Surfacephosavgyrct$TP)
colnames(Surfacephosavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(Surfacephosavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  Surfacephosavg%>%
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
Surfacephosavg<-anti_join(Surfacephosavg,eliminateS,by=c("ManagedAreaName"))
unique(Surfacephosavg$ManagedAreaName)
Surfacephosavgyrct<-Surfacephosavg

colnames(Surfacephosavg)
colnames(Surfacephosavgyrct)

# PLOTTING
colnames(Surfacephosavgyrct)
Surfacephosavgyrct$yearmonth<-ifelse(Surfacephosavgyrct$mth==10|Surfacephosavgyrct$mth==11|Surfacephosavgyrct$mth==12, (paste(Surfacephosavgyrct$yr, Surfacephosavgyrct$mth, sep=".")),
                                     (paste(Surfacephosavgyrct$yr, Surfacephosavgyrct$mth, sep=".0")))

unique(Surfacephosavgyrct$yearmonth)
###
plots_list <- lapply(unique(Surfacephosavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(Surfacephosavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = TP)) +
    geom_point() +
    labs(title = i, x="year.month",y="Total Phosphorus mg/L") +
    scale_x_discrete(breaks = Surfacephosavgyrct$yearmonth[seq(1, length(Surfacephosavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Nutrients_P_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(Surfacephosavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(TP~mth+yr, data = subset(Surfacephosavgyrct, ManagedAreaName == i))

#First, view the vector
unique(Surfacephosavgyrct$ManagedAreaName)
colnames(Surfacephosavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(Surfacephosavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(TP~mth + yr, data = subset(Surfacephosavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(Surfacephosavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Nutrients_P_MA_surf_KToutput.csv")
kendall_ManagedAreaNameal_list_Surface
sink()
#

#********************************************************
#4b. GENERATE RESULTS TABLES

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
write_csv(kendall_with_regions_final_surface, "Nutrients_P_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(phos3surf)
SUMMARYphossurfmeanmedstdev<-
  phos3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(TP),max=max(TP),median=median(TP),mean=mean(TP),sd=sd(TP),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYphossurfmeanmedstdev)

colnames(SUMMARYphossurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYphossurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYphossurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYphossurfmeanmedstdev) [4] <- "Minimum Value (Total Phosphorus mg/L)" 
colnames(SUMMARYphossurfmeanmedstdev) [5] <- "Maximum Value (Total Phosphorus mg/L)" 
colnames(SUMMARYphossurfmeanmedstdev) [6] <- "Median (Total Phosphorus mg/L)" 
colnames(SUMMARYphossurfmeanmedstdev) [7] <- "Mean (Total Phosphorus mg/L)" 
colnames(SUMMARYphossurfmeanmedstdev) [8] <- "Standard Deviation (Total Phosphorus mg/L)" 
colnames(SUMMARYphossurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYphossurfmeanmedstdev,"Nutrients_P_MA_surf_resultssummarystats2.csv")
