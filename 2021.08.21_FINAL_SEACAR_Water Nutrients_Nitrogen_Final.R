#SEACAR Water Column_Nutrient Scoring: Total Nitrogen- #STOPPED AT LINE 629
#Last edited by  Katie May Laumann, klaumann@umces.edu, on 2021.05.17 to quantify and remove applicable qualifier codes
#edited by edited Katie May Laumann, klaumann@umces.edu, on 2021.08.25 to add in table formatting and remove MAs with incorrect data
#QAQC'd by J. Edgerton


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
nitro <- read.csv("~/Desktop/Nitrogen-2021-Jul-26.csv")
#correct the MA names
nitro$ManagedAreaName<-ifelse(nitro$ManagedAreaName=="Apalachicola NERR","Apalachicola National Estuarine Research Reserve",
                             ifelse(nitro$ManagedAreaName=="Florida Keys NMS","Florida Keys National Marine Sanctuary",
                                    ifelse(nitro$ManagedAreaName=="Coral ECA","Southeast Florida Coral Reef Ecosystem Conservation Area",
                                           ifelse(nitro$ManagedAreaName=="Guana Tolomato Matanzas NERR","Guana Tolomato Matanzas National Estuarine Research Reserve",
                                                  ifelse(nitro$ManagedAreaName=="Rookery Bay NERR","Rookery Bay National Estuarine Research Reserve",
                                                         paste(nitro$ManagedAreaName, "Aquatic Preserve", sep=" "))))))

unique(nitro$ManagedAreaName)

#read in a "regions" df to help in table-building, later.
regions <- read.csv("~/Desktop/SEACARMARegionList.csv")
unique(regions$ManagedAreaName)
glimpse(regions)

#view the column names
colnames(nitro)

#********************************************************
#1. PREPARE DATA FOR ANALYSES
colnames(nitro)[colnames(nitro) == "Activity_Start_Date_Time"] <- "sampledate"
unique(nitro$X.TotalNitrogen_mg.L)
colnames(nitro)[colnames(nitro) == "X.TotalNitrogen_mg.L."] <- "Total.Nitrogen_mg.L"
unique(nitro$Total.Nitrogen_mg.L)
colnames(nitro)[colnames(nitro) == "X.TotalKjeldahlNitrogenTKN_mg.L."] <- "Total.Kjeldahl.Nitrogen.TKN_mg.L"
colnames(nitro)[colnames(nitro) == "X.NO2.3Filtered_mg.L."] <- "NO2.3.Filtered_mg.L"
colnames(nitro)[colnames(nitro) == "X.NH4Filtered_mg.L."] <- "NH4.Filtered_mg.L"
nrow(nitro)
#first, check what the values are. We can't just remove NAs or negative values, as we ultimately
#need to add together 3 columns. So:
  #get rid of negative values by replacing them with "999999999" 
  #Then, convert NAs to 0.
  #Do this for the 3 different measures we will be using to sum up to Total Nitrogen 
  #(Total.Nitrogen_mg.L, Total.Kjeldahl.Nitrogen.TKN_mg.L, and NO2.3.Filtered_mg.L)
  #Then remove "999999999" values (previously negative values)
colnames(nitro)
unique(nitro$Total.Nitrogen_mg.L)
nitro$Total.Nitrogen_mg.L<-ifelse(nitro$Total.Nitrogen_mg.L>0,nitro$Total.Nitrogen_mg.L,999999999)
nitro$Total.Nitrogen_mg.L[is.na(nitro$Total.Nitrogen_mg.L)] <- 0
unique(nitro$Total.Nitrogen_mg.L)

unique(nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L)
nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L<-ifelse(nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L>0,nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L,999999999)
nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L[is.na(nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L)] <- 0
unique(nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L)

nitro$NO2.3.Filtered_mg.L<-ifelse(nitro$NO2.3.Filtered_mg.L>0,nitro$NO2.3.Filtered_mg.L,999999999)
nitro$NO2.3.Filtered_mg.L[is.na(nitro$NO2.3.Filtered_mg.L)] <- 0
unique(nitro$NO2.3.Filtered_mg.L)

#filter out formerly negative values
nitro<-
  nitro%>%
  filter(nitro$Total.Nitrogen_mg.L!=999999999)

nitro<-
  nitro%>%
  filter(nitro$NO2.3.Filtered_mg.L!=999999999)

nitro<-
  nitro%>%
  filter(Total.Kjeldahl.Nitrogen.TKN_mg.L!=999999999)

#check values again
unique(nitro$Total.Nitrogen_mg.L)
unique(nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L)
unique(nitro$NO2.3.Filtered_mg.L)

#sum the 3 columns with nitrogen measures to calculate TN, or Total Nitrogen
nitro$TN<-nitro$Total.Nitrogen_mg.L+nitro$Total.Kjeldahl.Nitrogen.TKN_mg.L+nitro$NO2.3.Filtered_mg.L
min(nitro$TN)
max(nitro$TN)
########################
#Prepare date to be read in properly
unique(nitro$Month)
unique(nitro$Year)

#We need date to be read, by R, as a date. See what class it is.
class(nitro$sampledate)
#need to change it from a character to a factor, then a date.
##first, make new column "sampledate2", where it is a character
nitro$sampledate2<-as.character(nitro$sampledate)
class(nitro$sampledate2)
##then, in the original column, paste the "sampledate2" characters as a DATE object in yyy-mm-dd format.
nitro$sampledate <- as.Date(nitro$sampledate2, format="%m/%d/%y")
###works to here with old file
class(nitro$sampledate)

#we'll remove NAs (and invalid measurements). When we're done, we want to see how much of the data we
##have omitted. Store the current number of rows in an object "originalN"
originalN<-nrow(nitro)
#remove NAs
#make a new df (nitro2), and remove any remaining nitro NAs. This way we can always look at the 
  #original without reloading files.
nitro2<-nitro[!is.na(nitro$TN),]
colnames(nitro2)
nrow(nitro2)
#remove Relative_Depth NAs
#nitro2<-nitro2[!is.na(nitro2$Relative_Depth),]  #try also doing calculations for 'NAs'
#nrow(nitro2)
#Remove date NAs
nitro2<-nitro2[!is.na(nitro2$Year),]
nitro2<-nitro2[!is.na(nitro2$Month),]
nitro2<-nitro2[!is.na(nitro2$sampledate),]
noNAN<-nrow(nitro2)
#how many rows have we removed so far?
originalN-noNAN
colnames(nitro2)
#
#REMOVE CERTAIN QUALIFIER CODES AND FLAGS
#Identify in a new column, then remove, qualifier codes and flagged codes as decided with SMEs
##Starting with Value_Qualifier_STORET_WIN
#remove: Q, H, Y, J, N, D, G, T
#for U, the indicator value is already replaced with MDL value, so we don't need to do anything with it here.
class(nitro2$Value_Qualifier_STORET_WIN)
nrow(nitro2)
nitro2$Qinv<-grepl("Q",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Hinv<-grepl("H",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Yinv<-grepl("Y",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Jinv<-grepl("J",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Ninv<-grepl("N",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Dinv<-grepl("D",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Ginv<-grepl("G",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$Tinv<-grepl("T",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
#nitro2$Uinv<-grepl("U",nitro2$Value_Qualifier_STORET_WIN,fixed=TRUE)
nitro2$STORETWIN_remove<-
  ifelse(nitro2$Qinv=="TRUE","remove",
         ifelse(nitro2$Hinv=="TRUE","remove",
                ifelse(nitro2$Yinv=="TRUE","remove",
                       ifelse(nitro2$Jinv=="TRUE","remove",
                              ifelse(nitro2$Ninv=="TRUE","remove",
                                     ifelse(nitro2$Dinv=="TRUE","remove",
                                            ifelse(nitro2$Ginv=="TRUE","remove",
                                                   ifelse(nitro2$Tinv=="TRUE","remove","KEEP"))))))))
nrow(nitro2)
nitro2<-nitro2%>%filter(STORETWIN_remove=="KEEP")
nrow(nitro2)

##Next, Value_Qualifier_NERR
#remove: A, D, H, K, M, N, U, S
class(nitro2$Value_Qualifier_NERR)
nitro2$Ainv2<-grepl("A",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Dinv2<-grepl("D",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Hinv2<-grepl("H",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Kinv2<-grepl("K",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Minv2<-grepl("M",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Ninv2<-grepl("N",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Uinv2<-grepl("U",nitro2$Value_Qualifier_NERR,fixed=TRUE)
nitro2$Sinv2<-grepl("S",nitro2$Value_Qualifier_NERR,fixed=TRUE)

nitro2$NERR_remove<-
  ifelse(nitro2$Ainv2=="TRUE","remove",
         ifelse(nitro2$Dinv2=="TRUE","remove",
                ifelse(nitro2$Hinv2=="TRUE","remove",
                       ifelse(nitro2$Kinv2=="TRUE","remove",
                              ifelse(nitro2$Minv2=="TRUE","remove",
                                     ifelse(nitro2$Ninv2=="TRUE","remove",
                                            ifelse(nitro2$Uinv2=="TRUE","remove",
                                                   ifelse(nitro2$Sinv2=="TRUE","remove","KEEP"))))))))

nrow(nitro2)
nitro2<-nitro2%>%filter(NERR_remove=="KEEP")
nrow(nitro2)


##Next, QAQCFlag
#remove: <-5> Outside high sensor range; <-4> Outside low sensor range; <-3> Data rejected due to QA/QC; <-2> Missing Data; <1>  Suspect data

class(nitro2$QAQCFlag)
unique(nitro2$QAQCFlag)
nitro2$flag5<-grepl("<-5>",nitro2$QAQCFlag,fixed=TRUE)
nitro2$flag4<-grepl("<-4>",nitro2$QAQCFlag,fixed=TRUE)
nitro2$flag3<-grepl("<-3>",nitro2$QAQCFlag,fixed=TRUE)
nitro2$flag2<-grepl("<-2>",nitro2$QAQCFlag,fixed=TRUE)
nitro2$flag1<-grepl("<1>",nitro2$QAQCFlag,fixed=TRUE)


nitro2$FLAG_remove<-
  ifelse(nitro2$flag5=="TRUE","remove",
         ifelse(nitro2$flag4=="TRUE","remove",
                ifelse(nitro2$flag3=="TRUE","remove",
                       ifelse(nitro2$flag2=="TRUE","remove",
                              ifelse(nitro2$flag1=="TRUE","remove","KEEP")))))

nrow(nitro2)
nitro2<-nitro2%>%filter(FLAG_remove=="KEEP")
nrow(nitro2)


#********************************************************
#********************************************************
#GENERATE SUMMARY PLOTS
#SURFACE DATA boxplots
#make a new df with just the data you want. 
nitro3surf<-nitro2%>%
  filter(Relative_Depth=="Surface")

#create a yearmonth column for fine-scale plotting by year and month
nitro3surf$yearmonth<-ifelse(nitro3surf$Month==10|nitro3surf$Month==11|nitro3surf$Month==12, (paste(nitro3surf$Year, nitro3surf$Month, sep=".")),
                             (paste(nitro3surf$Year, nitro3surf$Month, sep=".0")))
colnames(nitro3surf)

#Test by making boxplots of ALL data from all ManagedAreaNames/Managed Areas together
#first at yearmonth level
boxplot(TN ~ yearmonth, data = nitro3surf, xlab = "year.month",
        ylab = "Total Nitrogen mg/L", main = "Total Nitrogen boxplot", notch=FALSE)
#then at year level
boxplot(TN ~ Year, data = nitro3surf, xlab = "Year",
        ylab = "Total Nitrogen mg/L", main = "Total Nitrogen boxplot", notch=FALSE)
#Now make a list of boxplots for each ManagedAreaName individually
#at yearmonth level
plots_listyrmth <- lapply(unique(nitro3surf$ManagedAreaName), function (i) {
  dat <- filter(nitro3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = yearmonth, y = TN)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )+
    labs(title = i, x="year.month",y="Total Nitrogen mg/L") 
})
#make a pdf with all of these
pdf("Nutrients_N_MA_surf_MonthlyBoxplots.pdf", width = 16, height = 4)
print(plots_listyrmth)
dev.off()

#and at Year level
plots_listyr <- lapply(unique(nitro3surf$ManagedAreaName), function (i) {
  dat <- filter(nitro3surf, ManagedAreaName == i)
  ggplot(data = dat, aes(x = Year, y = TN, group=Year)) +
    geom_boxplot(notch=FALSE) +
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=12)
    )+
    labs(title = i, x="year.month",y="Total Nitrogen mg/L") 
})
#make a pdf with all of these
pdf("Nutrients_N_MA_surf_AnnualBoxplots.pdf", width = 16, height = 4)
print(plots_listyr)
dev.off()



#********************************************************
#BEGIN ANALYSES

#We will run seasonal Kendall Tau, using monthly averages and being sure to account for
#Relative Depth. Make a new df that gives us average nitro by month, year, depth, and area. We will use this in analyses.
colnames(nitro2)
nitroavg<-
  nitro2%>%
  
  group_by(ManagedAreaName,Month, Year, Relative_Depth)%>%
  summarise(mean(TN))
#rename columns
colnames(nitroavg) 
colnames(nitroavg) [5]<-"TN"

#Make a new df that will let us count years
nitroavgyrct<-
  nitro2%>%
  
  group_by(ManagedAreaName,Month, Year,Relative_Depth)%>%
  summarise(mean(TN))
#rename columns
colnames(nitroavgyrct) 
colnames(nitroavgyrct) [5]<-"TN"
#################
#################
#Run analyses for Surface nitro:
#Make a new df with only Surface nitros so we can run analyses for those first 
#we have already removed NAs and qualifier codes from nitroavg and nitroavgyrct
colnames(nitroavg)
Surfacenitroavg<-
  nitroavg%>%
  filter(Relative_Depth=="Surface")
colnames(Surfacenitroavg)

Surfacenitroavgyrct<-
  nitroavgyrct%>%
  filter(Relative_Depth=="Surface")
colnames(Surfacenitroavgyrct)

#We need numerics to run kendallSeasonalTrendTest. 
#check to see if these are numerics, and if not, convert them.
class(Surfacenitroavg$Month)
class(Surfacenitroavg$Year)
Surfacenitroavg$mth<-as.numeric(Surfacenitroavg$Month)
Surfacenitroavg$yr<-as.numeric(Surfacenitroavg$Year)
class(Surfacenitroavg$mth)
class(Surfacenitroavg$TN)

class(Surfacenitroavgyrct$Month)
class(Surfacenitroavgyrct$Year)
Surfacenitroavgyrct$mth<-as.numeric(Surfacenitroavgyrct$Month)
class(Surfacenitroavgyrct$mth)
class(Surfacenitroavgyrct$Year)
Surfacenitroavgyrct$yr<-as.numeric(Surfacenitroavgyrct$Year)
class(Surfacenitroavgyrct$yr)
class(Surfacenitroavgyrct$TN)
colnames(Surfacenitroavgyrct)
#Make sure we only use Managed Areas with enough data
#what are the managed areas?
unique(Surfacenitroavg$ManagedAreaName)
#remove Managed Areas that don't have at least 10 years of data. 
##First, ID which ones have fewer than 10 yrs of data
nyearsS<-
  Surfacenitroavg%>%
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
Surfacenitroavg<-anti_join(Surfacenitroavg,eliminateS,by=c("ManagedAreaName"))
unique(Surfacenitroavg$ManagedAreaName)
Surfacenitroavgyrct<-Surfacenitroavg

colnames(Surfacenitroavg)
colnames(Surfacenitroavgyrct)

# PLOTTING
colnames(Surfacenitroavgyrct)
Surfacenitroavgyrct$yearmonth<-ifelse(Surfacenitroavgyrct$mth==10|Surfacenitroavgyrct$mth==11|Surfacenitroavgyrct$mth==12, (paste(Surfacenitroavgyrct$yr, Surfacenitroavgyrct$mth, sep=".")),
                                      (paste(Surfacenitroavgyrct$yr, Surfacenitroavgyrct$mth, sep=".0")))

unique(Surfacenitroavgyrct$yearmonth)
###
plots_list <- lapply(unique(Surfacenitroavgyrct$ManagedAreaName), function (i) {
  
  dat <- filter(Surfacenitroavgyrct, ManagedAreaName == i)
  
  ggplot(data = dat, aes(x = yearmonth, y = TN)) +
    geom_point() +
    labs(title = i, x="year.month",y="Total Nitrogen mg/L") +
    scale_x_discrete(breaks = Surfacenitroavgyrct$yearmonth[seq(1, length(Surfacenitroavgyrct$yearmonth), by = 3)])+
    theme(
      axis.text.x = element_text(angle = 90, hjust = 1,size=6)
    )
  
})
##
pdf("2Nutrients_N_MA_surf_plots.pdf", width = 16, height = 4)
print(plots_list)
dev.off()



#TESTING WITH ONE Managed Area before we run Kendall Tau for all at once
#*******************
library(EnvStats)
i=unique(Surfacenitroavgyrct$ManagedAreaName)[1]
kendallSeasonalTrendTest(TN~mth+yr, data = subset(Surfacenitroavgyrct, ManagedAreaName == i))

#First, view the vector
unique(Surfacenitroavgyrct$ManagedAreaName)
colnames(Surfacenitroavgyrct)

library(EnvStats)
kendall_ManagedAreaNameal_list_Surface<-lapply(unique(Surfacenitroavgyrct$ManagedAreaName), function(i)kendallSeasonalTrendTest(TN~mth + yr, data = subset(Surfacenitroavgyrct, ManagedAreaName == i)) )
names(kendall_ManagedAreaNameal_list_Surface) <- unique(Surfacenitroavgyrct$ManagedAreaName)
detach(package:EnvStats)
sink("Nutrients_N_MA_surf_KToutput.csv")
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
write_csv(kendall_with_regions_final_surface, "Nutrients_N_MA_surf_resultstables.csv")

#then write a summary stats table, including: mean, med, and standard deviation table
colnames(nitro3surf)
SUMMARYnitrosurfmeanmedstdev<-
  nitro3surf%>%
  group_by(ManagedAreaName,Year)%>%
  summarize(n=(length(unique(RowID))),min=min(TN),max=max(TN),median=median(TN),mean=mean(TN),sd=sd(TN),ProgramIDs=paste(unique(ProgramID), collapse = ";"))
colnames(SUMMARYnitrosurfmeanmedstdev)

colnames(SUMMARYnitrosurfmeanmedstdev) [1] <- "Managed Area Name" 
colnames(SUMMARYnitrosurfmeanmedstdev) [2] <- "Year" 
colnames(SUMMARYnitrosurfmeanmedstdev) [3] <- "Number of Samples" 
colnames(SUMMARYnitrosurfmeanmedstdev) [4] <- "Minimum Value (Total Nitrogen mg/L)" 
colnames(SUMMARYnitrosurfmeanmedstdev) [5] <- "Maximum Value (Total Nitrogen mg/L)" 
colnames(SUMMARYnitrosurfmeanmedstdev) [6] <- "Median (Total Nitrogen mg/L)" 
colnames(SUMMARYnitrosurfmeanmedstdev) [7] <- "Mean (Total Nitrogen mg/L)" 
colnames(SUMMARYnitrosurfmeanmedstdev) [8] <- "Standard Deviation (Total Nitrogen mg/L)" 
colnames(SUMMARYnitrosurfmeanmedstdev) [9] <- "Programs with Data" 

write.csv(SUMMARYnitrosurfmeanmedstdev,"Nutrients_N_MA_surf_resultssummarystats2.csv")
