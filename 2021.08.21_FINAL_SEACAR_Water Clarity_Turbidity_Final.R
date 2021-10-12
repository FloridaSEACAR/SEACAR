#SEACAR Water Column_Clarity: Turbidity
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton

#Seasonal Kendall Tau for long term trends in Turbidity at the MA level

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
turbid <- read.csv("~/Desktop/Turbidity and TSS-2021-Jul-26.csv")
#correct the MA names
turbid$ManagedAreaName<-ifelse(turbid$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                               ifelse(turbid$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                      ifelse(turbid$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                             ifelse(turbid$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                    ifelse(turbid$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                           paste(turbid$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(turbid)


#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
colnames(turbid)[colnames(turbid) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(turbid)[colnames(turbid) == "Turbidity_NTU"] <- "Turbid"

unique(turbid$Turbid)
class(turbid$Turbid)
nrow(turbid)

#remove NAs in Turbid
turbid<-turbid[!is.na(turbid$Turbid),]
#check the minimum and maximum values
min(turbid$Turbid)
max(turbid$Turbid)
#there are negative values. They are wrong. Make them NAs.
turbid$Turbid[turbid$Turbid<0] <- NA
#remove NAs again in Turbid
turbid<-turbid[!is.na(turbid$Turbid),]
#check the minimum and maximum values again
min(turbid$Turbid)
max(turbid$Turbid)
########################
#Prepare date to be read in properly
unique(turbid$Month)
unique(turbid$Year)

#We need date to be read, by R, as a date. See what class it is.
class(turbid$sampledate)
unique(turbid$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
turbid$sampledate2<-as.character(turbid$sampledate)
class(turbid$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
turbid$sampledate <- as.Date(turbid$sampledate2, format="%m/%d/%y")
###works to here with old file
class(turbid$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(turbid)
#remove NAs
#make a new df (turbid2), and remove any remaining turbid NAs. This way we can always look at the 
#original without reloading files.
turbid2<-turbid[!is.na(turbid$Turbid),]
colnames(turbid2)
nrow(turbid2)

#Remove date NAs
turbid2<-turbid2[!is.na(turbid2$Year),]
turbid2<-turbid2[!is.na(turbid2$Month),]
turbid2<-turbid2[!is.na(turbid2$sampledate),]
noNAN<-nrow(turbid2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(turbid2$Value_Qualifier_STORET_WIN)
nrow(turbid2)
turbid2$Qinv<-grepl("Q",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Hinv<-grepl("H",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Yinv<-grepl("Y",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Jinv<-grepl("J",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Ninv<-grepl("N",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Dinv<-grepl("D",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Ginv<-grepl("G",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$Tinv<-grepl("T",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#turbid2$Uinv<-grepl("U",turbid2$Value_Qualifier_STORET_WIN,fixed=TRUE)
turbid2$STORETWIN_remove<-
  ifelse(turbid2$Qinv=="TRUE","remove",
         ifelse(turbid2$Hinv=="TRUE","remove",
                ifelse(turbid2$Yinv=="TRUE","remove",
                       ifelse(turbid2$Jinv=="TRUE","remove",
                              ifelse(turbid2$Ninv=="TRUE","remove",
                                     ifelse(turbid2$Dinv=="TRUE","remove",
                                            ifelse(turbid2$Ginv=="TRUE","remove",
                                                   ifelse(turbid2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(turbid2)
turbid2<-turbid2%>%filter(STORETWIN_remove=="KEEP")
nrow(turbid2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(turbid2$Value_Qualifier_NERR)
turbid2$Ainv2<-grepl("A",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Dinv2<-grepl("D",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Hinv2<-grepl("H",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Kinv2<-grepl("K",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Minv2<-grepl("M",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Ninv2<-grepl("N",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Uinv2<-grepl("U",turbid2$Value_Qualifier_NERR,fixed=TRUE)
turbid2$Sinv2<-grepl("S",turbid2$Value_Qualifier_NERR,fixed=TRUE)

turbid2$NERR_remove<-
  ifelse(turbid2$Ainv2=="TRUE","remove",
         ifelse(turbid2$Dinv2=="TRUE","remove",
                ifelse(turbid2$Hinv2=="TRUE","remove",
                       ifelse(turbid2$Kinv2=="TRUE","remove",
                              ifelse(turbid2$Minv2=="TRUE","remove",
                                     ifelse(turbid2$Ninv2=="TRUE","remove",
                                            ifelse(turbid2$Uinv2=="TRUE","remove",
                                                   ifelse(turbid2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(turbid2)
turbid2<-turbid2%>%filter(NERR_remove=="KEEP")
nrow(turbid2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(turbid2$QAQCFlag)
unique(turbid2$QAQCFlag)
turbid2$flag5<-grepl("<-5>",turbid2$QAQCFlag,fixed=TRUE)
turbid2$flag4<-grepl("<-4>",turbid2$QAQCFlag,fixed=TRUE)
turbid2$flag3<-grepl("<-3>",turbid2$QAQCFlag,fixed=TRUE)
turbid2$flag2<-grepl("<-2>",turbid2$QAQCFlag,fixed=TRUE)
turbid2$flag1<-grepl("<1>",turbid2$QAQCFlag,fixed=TRUE)


turbid2$FLAG_remove<-
  ifelse(turbid2$flag5=="TRUE","remove",
         ifelse(turbid2$flag4=="TRUE","remove",
                ifelse(turbid2$flag3=="TRUE","remove",
                       ifelse(turbid2$flag2=="TRUE","remove",
                              ifelse(turbid2$flag1=="TRUE","remove","KEEP")))))

nrow(turbid2)
turbid2<-turbid2%>%filter(FLAG_remove=="KEEP")
nrow(turbid2)

#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
turbid3surf<-turbid2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
turbid3surf$yearmonth<-ifelse(turbid3surf$Month==10|turbid3surf$Month==11|turbid3surf$Month==12, (paste(turbid3surf$Year, turbid3surf$Month, sep=".")),
                              (paste(turbid3surf$Year, turbid3surf$Month, sep=".0")))
colnames(turbid3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(Turbid ~ yearmonth, data = turbid3surf, xlab = "yearmonth",
        ylab = "Turbidity (NTU)", main = "Turbidity boxplot", notch=FALSE)
#then at year level
boxplot(Turbid ~ Year, data = turbid3surf, xlab = "Year",
        ylab = "Turbidity (NTU)", main = "Turbidity boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(turbid3surf$ManagedAreaName), function (i) {
  dat <- filter(turbid3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = Turbid)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, x="year.month", y="Turbidity (NTU)")
})
#make a pdf with all of these
pdf("Clarity_Turbidity_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(turbid3surf$ManagedAreaName), function (i) {
  dat <- filter(turbid3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = Turbid, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, y="Turbidity (NTU)") 
})
#make a pdf with all of these
pdf("Clarity_Turbidity_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average turbid by month, year, depth, and area. We will use this in analyses.
colnames(turbid2)
turbidavg<-
  turbid2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(Turbid))
#rename columns
colnames(turbidavg) 
colnames(turbidavg) [5]<-"Turbid"

#Make a new df that will let us count years
turbidavgyrct<-
  turbid2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(Turbid))
#rename columns
colnames(turbidavgyrct) 
colnames(turbidavgyrct) [5]<-"Turbid"
#################
#################
#Run analyses for Surface turbid:
#Make a new df with only Surface turbids so we can run analyses for those first 
#we have already removed NAs and qualifier codes from turbidavg and turbidavgyrct
colnames(turbidavg)
Surfaceturbidavg<-
  turbidavg%>%
  filter(Relative_Depth=="Surface")
colnames(Surfaceturbidavg)

Surfaceturbidavgyrct<-
  turbidavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(Surfaceturbidavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(Surfaceturbidavg$Month)
class(Surfaceturbidavg$Year)
Surfaceturbidavg$mth<-as.numeric(Surfaceturbidavg$Month)
Surfaceturbidavg$yr<-as.numeric(Surfaceturbidavg$Year)
class(Surfaceturbidavg$mth)
class(Surfaceturbidavg$Turbid)

class(Surfaceturbidavgyrct$Month)
class(Surfaceturbidavgyrct$Year)
Surfaceturbidavgyrct$mth<-as.numeric(Surfaceturbidavgyrct$Month)
class(Surfaceturbidavgyrct$mth)
class(Surfaceturbidavgyrct$Year)
Surfaceturbidavgyrct$yr<-as.numeric(Surfaceturbidavgyrct$Year)
class(Surfaceturbidavgyrct$yr)
class(Surfaceturbidavgyrct$Turbid)
colnames(Surfaceturbidavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(Surfaceturbidavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  Surfaceturbidavg%>%
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
Surfaceturbidavg<-anti_join(Surfaceturbidavg,eliminateS,by=c("ManagedAreaName"))
unique(Surfaceturbidavg$ManagedAreaName)
Surfaceturbidavgyrct<-Surfaceturbidavg

colnames(Surfaceturbidavg)
colnames(Surfaceturbidavgyrct)

# PLOTTING
colnames(Surfaceturbidavgyrct)
Surfaceturbidavgyrct$yearmonth<-ifelse(Surfaceturbidavgyrct$mth==10|Surfaceturbidavgyrct$mth==11|Surfaceturbidavgyrct$mth==12, (paste(Surfaceturbidavgyrct$yr, Surfaceturbidavgyrct$mth, sep=".")),
                                       (paste(Surfaceturbidavgyrct$yr, Surfaceturbidavgyrct$mth, sep=".0")))

unique(Surfaceturbidavgyrct$yearmonth)
###
plots_list <- lapply(unique(Surfaceturbidavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(Surfaceturbidavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = Turbid)) +
    geom_point() +
    labs(title = i) +
    scale_x_discrete(breaks = Surfaceturbidavgyrct$yearmonth[seq(1, length(Surfaceturbidavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Clarity_Turbidity_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(Surfaceturbidavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(Turbid~mth+yr, data = subset(Surfaceturbidavgyrct, ManagedAreaName == i))

#First, view the vector
unique(Surfaceturbidavgyrct$ManagedAreaName)
colnames(Surfaceturbidavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(Surfaceturbidavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(Turbid~mth + yr, data = subset(Surfaceturbidavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(Surfaceturbidavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Clarity_Turbidity_MA_surf_KToutput.csv")
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
write_csv(kendall_with_regions_final_surface, "Clarity_Turbidity_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(turbid3surf)
SUMMARYturbidsurfmeanmedstdev<-
  turbid3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(Turbid),max=max(Turbid),median=median(Turbid),mean=mean(Turbid),sd=sd(Turbid),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYturbidsurfmeanmedstdev)

colnames(SUMMARYturbidsurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYturbidsurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYturbidsurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYturbidsurfmeanmedstdev) [4] <- "Minimum Value (Turbidity, NTU)" 
colnames(SUMMARYturbidsurfmeanmedstdev) [5] <- "Maximum Value (Turbidity, NTU)" 
colnames(SUMMARYturbidsurfmeanmedstdev) [6] <- "Median (Turbidity, NTU)" 
colnames(SUMMARYturbidsurfmeanmedstdev) [7] <- "Mean (Turbidity, NTU)" 
colnames(SUMMARYturbidsurfmeanmedstdev) [8] <- "Standard Deviation (Turbidity, NTU)" 
colnames(SUMMARYturbidsurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYturbidsurfmeanmedstdev,"Clarity_Turbidity_MA_surf_resultssummarystats2.csv")
