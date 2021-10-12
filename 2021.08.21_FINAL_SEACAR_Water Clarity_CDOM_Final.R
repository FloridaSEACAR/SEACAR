#SEACAR Water Column_Clarity: CDOM
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.08.21 
#Script written by Katie May Laumann, klaumann@umces.edu, on 2020.09.04
#QAQC'd by J. Edgerton

#Seasonal Kendall Tau for long term trends in CDOM at the MA level

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
cdom <- read.csv("~/Desktop/CDOM and FDOM-2021-Jul-26.csv")
#correct the MA names
cdom$ManagedAreaName<-ifelse(cdom$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                             ifelse(cdom$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                    ifelse(cdom$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                           ifelse(cdom$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                  ifelse(cdom$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                         paste(cdom$ManagedAreaName, "Aquatic Preserve", sep=" "))))))


#read in a "regions" df to help in table-building, later.
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(cdom)

#********************************************************
#PREPARE DATA FOR ANALYSES
#rename columns
unique(cdom$X.Coloreddissolvedorganicmatter)
colnames(cdom)[colnames(cdom) == "Activity_Start_Date_Time"] <- "sampledate"
colnames(cdom)[colnames(cdom) == "X.Coloreddissolvedorganicmatter"] <- "CDOM"

unique(cdom$CDOM)
class(cdom$CDOM)
nrow(cdom)

#remove NAs in CDOM
cdom<-cdom[!is.na(cdom$CDOM),]
#check the minimum and maximum values
min(cdom$CDOM)
max(cdom$CDOM)
#there are negative values. They are wrong. Make them NAs.
cdom$CDOM[cdom$CDOM<0] <- NA
#remove NAs again in CDOM
cdom<-cdom[!is.na(cdom$CDOM),]
#check the minimum and maximum values again
min(cdom$CDOM)
max(cdom$CDOM)
########################
#Prepare date to be read in properly
unique(cdom$Month)
unique(cdom$Year)

#We need date to be read, by R, as a date. See what class it is.
class(cdom$sampledate)
unique(cdom$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
cdom$sampledate2<-as.character(cdom$sampledate)
class(cdom$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
cdom$sampledate <- as.Date(cdom$sampledate2, format="%m/%d/%y")
###works to here with old file
class(cdom$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(cdom)
#remove NAs
#make a new df (cdom2), and remove any remaining cdom NAs. This way we can always look at the 
#original without reloading files.
cdom2<-cdom[!is.na(cdom$CDOM),]
colnames(cdom2)
nrow(cdom2)

#Remove date NAs
cdom2<-cdom2[!is.na(cdom2$Year),]
cdom2<-cdom2[!is.na(cdom2$Month),]
cdom2<-cdom2[!is.na(cdom2$sampledate),]
noNAN<-nrow(cdom2)
#how many rows have we removed so far?
originalN-noNAN
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(cdom2$Value_Qualifier_STORET_WIN)
nrow(cdom2)
cdom2$Qinv<-grepl("Q",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Hinv<-grepl("H",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Yinv<-grepl("Y",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Jinv<-grepl("J",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Ninv<-grepl("N",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Dinv<-grepl("D",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Ginv<-grepl("G",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$Tinv<-grepl("T",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#cdom2$Uinv<-grepl("U",cdom2$Value_Qualifier_STORET_WIN,fixed=TRUE)
cdom2$STORETWIN_remove<-
  ifelse(cdom2$Qinv=="TRUE","remove",
         ifelse(cdom2$Hinv=="TRUE","remove",
                ifelse(cdom2$Yinv=="TRUE","remove",
                       ifelse(cdom2$Jinv=="TRUE","remove",
                              ifelse(cdom2$Ninv=="TRUE","remove",
                                     ifelse(cdom2$Dinv=="TRUE","remove",
                                            ifelse(cdom2$Ginv=="TRUE","remove",
                                                   ifelse(cdom2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(cdom2)
cdom2<-cdom2%>%filter(STORETWIN_remove=="KEEP")
nrow(cdom2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(cdom2$Value_Qualifier_NERR)
cdom2$Ainv2<-grepl("A",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Dinv2<-grepl("D",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Hinv2<-grepl("H",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Kinv2<-grepl("K",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Minv2<-grepl("M",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Ninv2<-grepl("N",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Uinv2<-grepl("U",cdom2$Value_Qualifier_NERR,fixed=TRUE)
cdom2$Sinv2<-grepl("S",cdom2$Value_Qualifier_NERR,fixed=TRUE)

cdom2$NERR_remove<-
  ifelse(cdom2$Ainv2=="TRUE","remove",
         ifelse(cdom2$Dinv2=="TRUE","remove",
                ifelse(cdom2$Hinv2=="TRUE","remove",
                       ifelse(cdom2$Kinv2=="TRUE","remove",
                              ifelse(cdom2$Minv2=="TRUE","remove",
                                     ifelse(cdom2$Ninv2=="TRUE","remove",
                                            ifelse(cdom2$Uinv2=="TRUE","remove",
                                                   ifelse(cdom2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(cdom2)
cdom2<-cdom2%>%filter(NERR_remove=="KEEP")
nrow(cdom2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(cdom2$QAQCFlag)
unique(cdom2$QAQCFlag)
cdom2$flag5<-grepl("<-5>",cdom2$QAQCFlag,fixed=TRUE)
cdom2$flag4<-grepl("<-4>",cdom2$QAQCFlag,fixed=TRUE)
cdom2$flag3<-grepl("<-3>",cdom2$QAQCFlag,fixed=TRUE)
cdom2$flag2<-grepl("<-2>",cdom2$QAQCFlag,fixed=TRUE)
cdom2$flag1<-grepl("<1>",cdom2$QAQCFlag,fixed=TRUE)


cdom2$FLAG_remove<-
  ifelse(cdom2$flag5=="TRUE","remove",
         ifelse(cdom2$flag4=="TRUE","remove",
                ifelse(cdom2$flag3=="TRUE","remove",
                       ifelse(cdom2$flag2=="TRUE","remove",
                              ifelse(cdom2$flag1=="TRUE","remove","KEEP")))))

nrow(cdom2)
cdom2<-cdom2%>%filter(FLAG_remove=="KEEP")
nrow(cdom2)

#********************************************************
#GENERATE SUMMARY PLOTS
#All-depth DATA boxplots
#make a new df with just the data you want. 
cdom3surf<-cdom2#%>%
  #filter(Relative_Depth=="All-depth")

#create a yearmonth column for fine-scale plotting by year and month
cdom3surf$yearmonth<-ifelse(cdom3surf$Month==10|cdom3surf$Month==11|cdom3surf$Month==12, (paste(cdom3surf$Year, cdom3surf$Month, sep=".")),
                            (paste(cdom3surf$Year, cdom3surf$Month, sep=".0")))
colnames(cdom3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(CDOM ~ yearmonth, data = cdom3surf, xlab = "yearmonth",
        ylab = "CDOM", main = "CDOM boxplot", notch=FALSE)
#then at year level
boxplot(CDOM ~ Year, data = cdom3surf, xlab = "Year",
        ylab = "CDOM", main = "CDOM boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(cdom3surf$ManagedAreaName), function (i) {
  dat <- filter(cdom3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = CDOM)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, x="year.month",y="CDOM") 
})
#make a pdf with all of these
pdf("Clarity_CDOM_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(cdom3surf$ManagedAreaName), function (i) {
  dat <- filter(cdom3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = CDOM, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, x="year.month",y="CDOM") 
})
#make a pdf with all of these
pdf("Clarity_CDOM_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average cdom by month, year, depth, and area. We will use this in analyses.
colnames(cdom2)
cdomavg<-
  cdom2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(CDOM))
#rename columns
colnames(cdomavg) 
colnames(cdomavg) [5]<-"CDOM"

#Make a new df that will let us count years
cdomavgyrct<-
  cdom2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(CDOM))
#rename columns
colnames(cdomavgyrct) 
colnames(cdomavgyrct) [5]<-"CDOM"
#################
#################
#Run analyses for All-depth cdom:
#Make a new df with only All-depth cdoms so we can run analyses for those first 
#we have already removed NAs and qualifier codes from cdomavg and cdomavgyrct
colnames(cdomavg)
All-depthcdomavg<-
  cdomavg#%>%
  #filter(Relative_Depth=="All-depth")
colnames(All-depthcdomavg)

All-depthcdomavgyrct<-
  cdomavgyrct#%>%
  #filter(Relative_Depth=="All-depth")
colnames(All-depthcdomavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(All-depthcdomavg$Month)
class(All-depthcdomavg$Year)
All-depthcdomavg$mth<-as.numeric(All-depthcdomavg$Month)
All-depthcdomavg$yr<-as.numeric(All-depthcdomavg$Year)
class(All-depthcdomavg$mth)
class(All-depthcdomavg$CDOM)

class(All-depthcdomavgyrct$Month)
class(All-depthcdomavgyrct$Year)
All-depthcdomavgyrct$mth<-as.numeric(All-depthcdomavgyrct$Month)
class(All-depthcdomavgyrct$mth)
class(All-depthcdomavgyrct$Year)
All-depthcdomavgyrct$yr<-as.numeric(All-depthcdomavgyrct$Year)
class(All-depthcdomavgyrct$yr)
class(All-depthcdomavgyrct$CDOM)
colnames(All-depthcdomavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(All-depthcdomavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  All-depthcdomavg%>%
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
All-depthcdomavg<-anti_join(All-depthcdomavg,eliminateS,by=c("ManagedAreaName"))
unique(All-depthcdomavg$ManagedAreaName)
All-depthcdomavgyrct<-All-depthcdomavg

colnames(All-depthcdomavg)
colnames(All-depthcdomavgyrct)

# PLOTTING
colnames(All-depthcdomavgyrct)
All-depthcdomavgyrct$yearmonth<-ifelse(All-depthcdomavgyrct$mth==10|All-depthcdomavgyrct$mth==11|All-depthcdomavgyrct$mth==12, (paste(All-depthcdomavgyrct$yr, All-depthcdomavgyrct$mth, sep=".")),
                                     (paste(All-depthcdomavgyrct$yr, All-depthcdomavgyrct$mth, sep=".0")))

unique(All-depthcdomavgyrct$yearmonth)
###
plots_list <- lapply(unique(All-depthcdomavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(All-depthcdomavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = CDOM)) +
    geom_point() +
    labs(title = i, x="year.month",y="CDOM") +
    scale_x_discrete(breaks = All-depthcdomavgyrct$yearmonth[seq(1, length(All-depthcdomavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Clarity_CDOM_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(All-depthcdomavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(CDOM~mth+yr, data = subset(All-depthcdomavgyrct, ManagedAreaName == i))

#First, view the vector
unique(All-depthcdomavgyrct$ManagedAreaName)
colnames(All-depthcdomavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_All-depth<-lapply(unique(All-depthcdomavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(CDOM~mth + yr, data = subset(All-depthcdomavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_All-depth) <- unique(All-depthcdomavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Clarity_CDOM_MA_surf_KToutput.csv")
kendall_ManagedAreaNameal_list_All-depth
sink()
#

#********************************************************
#GENERATE RESULTS TABLES

#first, results table (after seasonal test)
kendall_map_All-depth <- map(.x = kendall_ManagedAreaNameal_list_All-depth, .f = broom::tidy)

### excel export version
formatted_kendall_All-depth <- enframe(kendall_map_All-depth) %>% 
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


kendall_with_regions_All-depth <- left_join(formatted_kendall_All-depth, regions, by = c("Managed Area" = "ManagedAreaName"))

kendall_with_regions_final_All-depth <- kendall_with_regions_All-depth %>% 
  select(-X) %>% 
  relocate(Region, .before = `Managed Area`)

kendall_with_regions_final_All-depth %>%
  gt() 

colnames(kendall_with_regions_final_All-depth)
colnames(kendall_with_regions_final_All-depth) [7] <- "zpvalue" 
kendall_with_regions_final_All-depth$significance<-ifelse(kendall_with_regions_final_All-depth$zpvalue<=0.05,"significant","not significant")
colnames(kendall_with_regions_final_All-depth) [7] <- "z statistic p-value" 
colnames(kendall_with_regions_final_All-depth) [10] <- "tau statistic trend direction" 
kendall_with_regions_final_All-depth<-kendall_with_regions_final_All-depth %>% select(-11) 
write_csv(kendall_with_regions_final_All-depth, "Clarity_CDOM_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(cdom3surf)
SUMMARYcdomsurfmeanmedstdev<-
  cdom3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(CDOM),max=max(CDOM),median=median(CDOM),mean=mean(CDOM),sd=sd(CDOM),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYcdomsurfmeanmedstdev)

colnames(SUMMARYcdomsurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYcdomsurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYcdomsurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYcdomsurfmeanmedstdev) [4] <- "Minimum Value (CDOM)" 
colnames(SUMMARYcdomsurfmeanmedstdev) [5] <- "Maximum Value (CDOM)" 
colnames(SUMMARYcdomsurfmeanmedstdev) [6] <- "Median (CDOM)" 
colnames(SUMMARYcdomsurfmeanmedstdev) [7] <- "Mean (CDOM)" 
colnames(SUMMARYcdomsurfmeanmedstdev) [8] <- "Standard Deviation (CDOM)" 
colnames(SUMMARYcdomsurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYcdomsurfmeanmedstdev,"Clarity_CDOM_MA_surf_resultssummarystats2.csv")

