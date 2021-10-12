#SEACAR Water Column_Clarity: chla
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton

#Seasonal Kendall Tau for long term trends in chla at the MA level

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
chlA <- read.csv("~/Desktop/Chlorophyll-2021-Jul-26.csv")
#correct the MA names
chlA$ManagedAreaName<-ifelse(chlA$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                             ifelse(chlA$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                    ifelse(chlA$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                           ifelse(chlA$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                  ifelse(chlA$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                         paste(chlA$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(chlA)

#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
colnames(chlA)[colnames(chlA) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(chlA)[colnames(chlA) == "X.Chlorophyllauncorrectedforpheophytin_ug.L."] <- "chla"

unique(chlA$chla)
class(chlA$chla)
nrow(chlA)

#remove NAs in chla
chlA<-chlA[!is.na(chlA$chla),]
#check the minimum and maximum values
min(chlA$chla)
max(chlA$chla)
#there are negative values. They are wrong. Make them NAs.
chlA$chla[chlA$chla<0] <- NA
#remove NAs again in chla
chlA<-chlA[!is.na(chlA$chla),]
#check the minimum and maximum values again
min(chlA$chla)
max(chlA$chla)
########################
#Prepare date to be read in properly
unique(chlA$Month)
unique(chlA$Year)

#We need date to be read, by R, as a date. See what class it is.
class(chlA$sampledate)
unique(chlA$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
chlA$sampledate2<-as.character(chlA$sampledate)
class(chlA$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
chlA$sampledate <- as.Date(chlA$sampledate2, format="%m/%d/%y")
###works to here with old file
class(chlA$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(chlA)
#remove NAs
#make a new df (chlA2), and remove any remaining chlA NAs. This way we can always look at the 
#original without reloading files.
chlA2<-chlA[!is.na(chlA$chla),]
colnames(chlA2)
nrow(chlA2)

#Remove date NAs
chlA2<-chlA2[!is.na(chlA2$Year),]
chlA2<-chlA2[!is.na(chlA2$Month),]
chlA2<-chlA2[!is.na(chlA2$sampledate),]
noNAN<-nrow(chlA2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(chlA2$Value_Qualifier_STORET_WIN)
nrow(chlA2)
chlA2$Qinv<-grepl("Q",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Hinv<-grepl("H",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Yinv<-grepl("Y",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Jinv<-grepl("J",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Ninv<-grepl("N",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Dinv<-grepl("D",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Ginv<-grepl("G",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$Tinv<-grepl("T",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#chlA2$Uinv<-grepl("U",chlA2$Value_Qualifier_STORET_WIN,fixed=TRUE)
chlA2$STORETWIN_remove<-
  ifelse(chlA2$Qinv=="TRUE","remove",
         ifelse(chlA2$Hinv=="TRUE","remove",
                ifelse(chlA2$Yinv=="TRUE","remove",
                       ifelse(chlA2$Jinv=="TRUE","remove",
                              ifelse(chlA2$Ninv=="TRUE","remove",
                                     ifelse(chlA2$Dinv=="TRUE","remove",
                                            ifelse(chlA2$Ginv=="TRUE","remove",
                                                   ifelse(chlA2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(chlA2)
chlA2<-chlA2%>%filter(STORETWIN_remove=="KEEP")
nrow(chlA2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(chlA2$Value_Qualifier_NERR)
chlA2$Ainv2<-grepl("A",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Dinv2<-grepl("D",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Hinv2<-grepl("H",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Kinv2<-grepl("K",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Minv2<-grepl("M",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Ninv2<-grepl("N",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Uinv2<-grepl("U",chlA2$Value_Qualifier_NERR,fixed=TRUE)
chlA2$Sinv2<-grepl("S",chlA2$Value_Qualifier_NERR,fixed=TRUE)

chlA2$NERR_remove<-
  ifelse(chlA2$Ainv2=="TRUE","remove",
         ifelse(chlA2$Dinv2=="TRUE","remove",
                ifelse(chlA2$Hinv2=="TRUE","remove",
                       ifelse(chlA2$Kinv2=="TRUE","remove",
                              ifelse(chlA2$Minv2=="TRUE","remove",
                                     ifelse(chlA2$Ninv2=="TRUE","remove",
                                            ifelse(chlA2$Uinv2=="TRUE","remove",
                                                   ifelse(chlA2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(chlA2)
chlA2<-chlA2%>%filter(NERR_remove=="KEEP")
nrow(chlA2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(chlA2$QAQCFlag)
unique(chlA2$QAQCFlag)
chlA2$flag5<-grepl("<-5>",chlA2$QAQCFlag,fixed=TRUE)
chlA2$flag4<-grepl("<-4>",chlA2$QAQCFlag,fixed=TRUE)
chlA2$flag3<-grepl("<-3>",chlA2$QAQCFlag,fixed=TRUE)
chlA2$flag2<-grepl("<-2>",chlA2$QAQCFlag,fixed=TRUE)
chlA2$flag1<-grepl("<1>",chlA2$QAQCFlag,fixed=TRUE)


chlA2$FLAG_remove<-
  ifelse(chlA2$flag5=="TRUE","remove",
         ifelse(chlA2$flag4=="TRUE","remove",
                ifelse(chlA2$flag3=="TRUE","remove",
                       ifelse(chlA2$flag2=="TRUE","remove",
                              ifelse(chlA2$flag1=="TRUE","remove","KEEP")))))

nrow(chlA2)
chlA2<-chlA2%>%filter(FLAG_remove=="KEEP")
nrow(chlA2)

#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
chlA3surf<-chlA2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
chlA3surf$yearmonth<-ifelse(chlA3surf$Month==10|chlA3surf$Month==11|chlA3surf$Month==12, (paste(chlA3surf$Year, chlA3surf$Month, sep=".")),
                            (paste(chlA3surf$Year, chlA3surf$Month, sep=".0")))
colnames(chlA3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(chla ~ yearmonth, data = chlA3surf, xlab = "yearmonth",
        ylab = "chl a, ug/L", main = "chla boxplot", notch=FALSE)
#then at year level
boxplot(chla ~ Year, data = chlA3surf, xlab = "Year",
        ylab = "chl a, ug/L", main = "chla boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(chlA3surf$ManagedAreaName), function (i) {
  dat <- filter(chlA3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = chla)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, y="chla (ug/L)") 
})
#make a pdf with all of these
pdf("Clarity_chla_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(chlA3surf$ManagedAreaName), function (i) {
  dat <- filter(chlA3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = chla, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, y="chla (ug/L)") 
})
#make a pdf with all of these
pdf("Clarity_chla_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average chlA by month, year, depth, and area. We will use this in analyses.
colnames(chlA2)
chlAavg<-
  chlA2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(chla))
#rename columns
colnames(chlAavg) 
colnames(chlAavg) [5]<-"chla"

#Make a new df that will let us count years
chlAavgyrct<-
  chlA2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(chla))
#rename columns
colnames(chlAavgyrct) 
colnames(chlAavgyrct) [5]<-"chla"
#################
#################
#Run analyses for Surface chlA:
#Make a new df with only Surface chlAs so we can run analyses for those first 
#we have already removed NAs and qualifier codes from chlAavg and chlAavgyrct
colnames(chlAavg)
SurfacechlAavg<-
  chlAavg%>%
  filter(Relative_Depth=="Surface")
colnames(SurfacechlAavg)

SurfacechlAavgyrct<-
  chlAavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(SurfacechlAavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(SurfacechlAavg$Month)
class(SurfacechlAavg$Year)
SurfacechlAavg$mth<-as.numeric(SurfacechlAavg$Month)
SurfacechlAavg$yr<-as.numeric(SurfacechlAavg$Year)
class(SurfacechlAavg$mth)
class(SurfacechlAavg$chla)

class(SurfacechlAavgyrct$Month)
class(SurfacechlAavgyrct$Year)
SurfacechlAavgyrct$mth<-as.numeric(SurfacechlAavgyrct$Month)
class(SurfacechlAavgyrct$mth)
class(SurfacechlAavgyrct$Year)
SurfacechlAavgyrct$yr<-as.numeric(SurfacechlAavgyrct$Year)
class(SurfacechlAavgyrct$yr)
class(SurfacechlAavgyrct$chla)
colnames(SurfacechlAavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(SurfacechlAavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  SurfacechlAavg%>%
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
SurfacechlAavg<-anti_join(SurfacechlAavg,eliminateS,by=c("ManagedAreaName"))
unique(SurfacechlAavg$ManagedAreaName)
SurfacechlAavgyrct<-SurfacechlAavg

colnames(SurfacechlAavg)
colnames(SurfacechlAavgyrct)

# PLOTTING
colnames(SurfacechlAavgyrct)
SurfacechlAavgyrct$yearmonth<-ifelse(SurfacechlAavgyrct$mth==10|SurfacechlAavgyrct$mth==11|SurfacechlAavgyrct$mth==12, (paste(SurfacechlAavgyrct$yr, SurfacechlAavgyrct$mth, sep=".")),
                                     (paste(SurfacechlAavgyrct$yr, SurfacechlAavgyrct$mth, sep=".0")))

unique(SurfacechlAavgyrct$yearmonth)
###
plots_list <- lapply(unique(SurfacechlAavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(SurfacechlAavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = chla)) +
    geom_point() +
    labs(title = i, y="chla (ug/L)") +
    scale_x_discrete(breaks = SurfacechlAavgyrct$yearmonth[seq(1, length(SurfacechlAavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Clarity_chla_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(SurfacechlAavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(chla~mth+yr, data = subset(SurfacechlAavgyrct, ManagedAreaName == i))

#First, view the vector
unique(SurfacechlAavgyrct$ManagedAreaName)
colnames(SurfacechlAavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(SurfacechlAavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(chla~mth + yr, data = subset(SurfacechlAavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(SurfacechlAavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Clarity_chla_MA_surf_KToutput.csv")
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
write_csv(kendall_with_regions_final_surface, "Clarity_chla_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(chlA3surf)
SUMMARYchlAsurfmeanmedstdev<-
  chlA3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(chla),max=max(chla),median=median(chla),mean=mean(chla),sd=sd(chla),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYchlAsurfmeanmedstdev)

colnames(SUMMARYchlAsurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYchlAsurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYchlAsurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYchlAsurfmeanmedstdev) [4] <- "Minimum Value (chla ug/L)" 
colnames(SUMMARYchlAsurfmeanmedstdev) [5] <- "Maximum Value (chla ug/L)" 
colnames(SUMMARYchlAsurfmeanmedstdev) [6] <- "Median (chla ug/L)" 
colnames(SUMMARYchlAsurfmeanmedstdev) [7] <- "Mean (chla ug/L)" 
colnames(SUMMARYchlAsurfmeanmedstdev) [8] <- "Standard Deviation (chla ug/L)" 
colnames(SUMMARYchlAsurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYchlAsurfmeanmedstdev,"Clarity_chla_MA_surf_resultssummarystats2.csv")
