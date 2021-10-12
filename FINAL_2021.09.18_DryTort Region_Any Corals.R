#FLORIDA SEACAR
#Coral Scripts, Percent Cover Analyses
#Script written by Katie May Laumann, klaumann@umces.edu, on 25 June 2020
#Script last modified: 19 Sept 2021 by Katie May Laumann, klaumann@umces.edu

#This script uses linear mixed effects models to find the effect of time (Year) on Percent Cover 
#in each Coral Region. Notes on modifying this script for Managed Area level analyses are included.
#Notes for changing/modifying model parameters, which should be done based on the data used for each 
#analysis, are provided.

#clear the environment
rm(list = ls())
#load packages
library(tidyverse)
library(dplyr)

#read in data
CorDRYTORTUGAS<-read.csv("~/Downloads/Percent Cover  - DRY TORT-2021-Jul-26.csv")
colnames(CorDRYTORTUGAS)

#run all analyses on fixed and random samples independently. 
#At the time of code production, no column is available in data to do so.
#Once this is corrected, simply run the two following lines of code, and run the remaining script 
#once for each (Corfixed and Corrand).
#CorDryfixed<-CorDRY%>%filter(samplingmethod=="fixed")
#CorDryfixed<-CorDRY%>%filter(samplingmethod=="random")
#if a certain program (eg 62) should be left out, remove it from analyses using the following:
#CorDRY<-CorDRY%>%filter(ProgramID!="62")

#Explore the data and prepare them for analysis.
#check unique values in some columns
unique(CorDRYTORTUGAS$SpeciesGroup1)
unique(CorDRYTORTUGAS$SpeciesGroup2)

#We have more than corals here: see Groups 1 and 2.
#We only want corals, so filter by Group 2.
CorDRYTORTUGAS<-CorDRYTORTUGAS%>%filter(SpeciesGroup2=="Octocoral"|SpeciesGroup2=="Unspecifiedcoral")
#if you want to exclude Octocorals, run the next line instead:
#CorDRYTORTUGAS<-CorDRYTORTUGAS%>%filter(SpeciesGroup2=="Unspecifiedcoral")

#Current datasets do not include sampling method, which was requested. If sampling method is added to the dataset,
#seperate analyses by sampling method by first building a df for each method:
#CorDRYTORTUGASMethod1<-CorDRYTORTUGAS%>%MethodRowName=="Method1"
#CorDRYTORTUGASMethod2<-CorDRYTORTUGAS%>%MethodRowName=="Method2"
#then by running all of the following analyses for each method individually. 
#This could be simplified by using lapply, if desired.

#rename columns to more easily work with them
colnames(CorDRYTORTUGAS)[colnames(CorDRYTORTUGAS) == "X.PercentCover.SpeciesComposition_.."] <- "Perccov"
nrow(CorDRYTORTUGAS)
#remove % cover NAs
CorDRYTORTUGAS <- CorDRYTORTUGAS[!is.na(CorDRYTORTUGAS$Perccov), ]  
nrow(CorDRYTORTUGAS)

#!!! don't need the next 3 lines for Region-level analysis. You will need them for Managed Area Level analysis.
#Create a Region List that we will write to CorDRYTORTUGASPerccovRegionList.csv to help in building our final results table, later.
#CorDRYTORTUGASRegionList<-unique(CorDRYTORTUGAS[c("ManagedAreaName","Region")])
#write.csv(CorDRYTORTUGASRegionList,"CorDRYTORTUGASPerccovRegionList.csv")
#*******************************
#check unique values
unique(CorDRYTORTUGAS$Perccov)
#convert to a numeric to be able to do math
CorDRYTORTUGAS$Perccov<-as.numeric(CorDRYTORTUGAS$Perccov)
CorDRYTORTUGAS <- CorDRYTORTUGAS[!is.na(CorDRYTORTUGAS$Perccov), ]  
unique(CorDRYTORTUGAS$SpeciesName)
unique(CorDRYTORTUGAS$GenusName)
nrow(CorDRYTORTUGAS)

#combine genus and spp into one column in case you want to refine analyses to the species level.
CorDRYTORTUGAS$gensp<-paste(CorDRYTORTUGAS$GenusName, CorDRYTORTUGAS$SpeciesName, sep=" ")
nrow(CorDRYTORTUGAS)
unique(CorDRYTORTUGAS$gensp)
CorDRYTORTUGAS<-CorDRYTORTUGAS%>%filter(gensp!="NA NA")

#how many rows?
allsamples<-nrow(CorDRYTORTUGAS)
class(CorDRYTORTUGAS$Year)
unique(CorDRYTORTUGAS$Year)

#GENERATE SUMMARY STATS TABLE
#BUILD A SUMMARY STATS TABLES with: # samples, min and max values, median, avg, and stdev of values, and program IDs
colnames(CorDRYTORTUGAS)
summaryDRYTORTUGASAnyCoral<-
  CorDRYTORTUGAS%>%
  group_by(Coral_Region, Year)%>%
  summarize(n=(length(unique(RowID))),min=min(Perccov),max=max(Perccov),median=median(Perccov),mean=mean(Perccov),sd=sd(Perccov),ProgramIDs=paste(unique(ProgramID), collapse = ","))
colnames(summaryDRYTORTUGASAnyCoral)
#rename columns in table we just created
colnames(summaryDRYTORTUGASAnyCoral) [1] <- "Coral Region" 
colnames(summaryDRYTORTUGASAnyCoral) [2] <- "Year" 
colnames(summaryDRYTORTUGASAnyCoral) [3] <- "Number of Samples" 
colnames(summaryDRYTORTUGASAnyCoral) [4] <- "Minimum Value (Percent Cover)" 
colnames(summaryDRYTORTUGASAnyCoral) [5] <- "Maximum Value (Percent Cover)" 
colnames(summaryDRYTORTUGASAnyCoral) [6] <- "Median (Percent Cover)" 
colnames(summaryDRYTORTUGASAnyCoral) [7] <- "Mean (Percent Cover)" 
colnames(summaryDRYTORTUGASAnyCoral) [8] <- "Standard Deviation (Percent Cover)" 
colnames(summaryDRYTORTUGASAnyCoral) [9] <- "Programs with Data" 
#write a summary stats .csv
write.csv(summaryDRYTORTUGASAnyCoral,"CoralSummaryStatisticsTable_AnyCoral_DRYTORTUGASRegion.csv")


##############################
#Analyses for ANY CORAL
#1. Summary statistics plots
par(mfrow=c(2,2))
plot(lm(CorDRYTORTUGAS$Year ~ CorDRYTORTUGAS$Perccov))
#print to file
pdf("CoralSummaryStatisticsPlots_AnyCoral_DRYTORTUGASRegion.pdf", width = 16, height = 4)
plot(lm(CorDRYTORTUGAS$Year ~ CorDRYTORTUGAS$Perccov))
dev.off()


#2. Boxplots
#BOXPLOTS
CorDRYTORTUGASboxplots<-CorDRYTORTUGAS# if you are going by MA, you will specifiy the desired MA here: %>% filter(ManagedAreaName=="")
boxplotsCorDRYTORTUGAS<- lapply(unique(CorDRYTORTUGAS$Coral_Region), function (i) {
  dat <- filter(CorDRYTORTUGASboxplots, Coral_Region == i)
  ggplot(data = dat, aes(group=Year, x = Year, y = Perccov)) +
    geom_boxplot() +
    labs(title = "Percent cover of any coral, Dry Tortugas Coral Region",subtitle=i, y="Percent Cover", x="Year") })
#view
boxplotsCorDRYTORTUGAS
#print to file
pdf("CoralBoxplots_AnyCoral_DRYTORTUGASRegion.pdf", width = 16, height = 4) 
print(boxplotsCorDRYTORTUGAS)
dev.off()
write.csv(CorDRYTORTUGAS,"CorDryTortugas.csv")

#3. Linear Mixed Effects Models
#two options to run lmes. Use lme4 OR nlme. Here, lme4 is hashtagged out, and we are running nlme.
#library(lme4)
#test1=lmer(Perccov ~ Year + (1|ProgramLocationID),data=CorDRYTORTUGAS)
#summary(test1)

#to run just one lme, do the following
library(nlme)
CorDRYTORTUGAS_AnyCoral<-lme(Perccov ~ Year,
                             random =~1|ProgramLocationID,
                             na.action = na.omit,
                             data = CorDRYTORTUGAS)
summary(CorDRYTORTUGAS_AnyCoral)

unique(CorDRYTORTUGAS$ProgramID)
unique(CorDRYTORTUGAS$ProgramLocationID)

#Here is where you examine your options for setting random variables and modifying your model.
#For example, there are only 2 programs sampling, but many program location IDs. For example, 
#you may want to allow for different starting amounts of coral in each program location ID. 
#To run multiple lmes at different levels, for example by SpeciesGroup2, do the following:
#Use the same type of model as above, but run it for each level (specified as i in the script)
#and print them to a list using lapply.

#In the following, we will allow for variation by program location ID. 
CorDRYTORTUGAS_AnyCoral_lmelist2<-lapply(unique(CorDRYTORTUGAS$Coral_Region), 
                                         function(i)summary(lme(Perccov ~ Year,
                                                                random =~1|ProgramLocationID,
                                                                na.action = na.omit, data = subset(CorDRYTORTUGAS, Coral_Region == i))))

#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
names(CorDRYTORTUGAS_AnyCoral_lmelist2) <- unique(CorDRYTORTUGAS$Coral_Region)
CorDRYTORTUGAS_AnyCoral_lmelist2

#If you want to run the same analysis at a different level, for example the SpeciesGroup2 level, run the following
#8 lines instead.
#CorDRYTORTUGAS_AnyCoral_lmelist2<-lapply(unique(CorDRYTORTUGAS$SpeciesGroup2), 
#function(i)summary(lme(Perccov ~ Year,
#random =~1|ProgramLocationID,
#na.action = na.omit, data = subset(CorDRYTORTUGAS, SpeciesGroup2 == i))))
#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
#names(CorDRYTORTUGAS_AnyCoral_lmelist2) <- unique(CorDRYTORTUGAS$SpeciesGroup2)
#CorDRYTORTUGAS_AnyCoral_lmelist2

#print to a csv file
sink("Corallme_AnyCoral_DRYTORTUGASRegion.csv")
CorDRYTORTUGAS_AnyCoral_lmelist2
sink()

#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(CorDRYTORTUGAS_AnyCoral_lmelist2), function(i) {
  mod <- CorDRYTORTUGAS_AnyCoral_lmelist2[[i]]
  pred <- data.frame(Coral_Region = i, Year = mod$data$Year, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_CorDRYTORTUGAS_AnyCoral<-ggplot() +
  geom_jitter(data = CorDRYTORTUGAS, aes(x = Year, y = Perccov, col = Coral_Region),size=0.1) +
  geom_line(data = preds, aes(x = Year, y = pred, group = Coral_Region, col = Coral_Region))
#label the plot and y axis (x was already labeled above)
Plot_CorDRYTORTUGAS_AnyCoral<-Plot_CorDRYTORTUGAS_AnyCoral + labs(y = "Percent Cover")
Plot_CorDRYTORTUGAS_AnyCoral<-Plot_CorDRYTORTUGAS_AnyCoral + labs(title = "Percent cover of any coral, Dry Tortugas Region")
#color code your points and lines by taxon
Plot_CorDRYTORTUGAS_AnyCoral<-Plot_CorDRYTORTUGAS_AnyCoral + labs(colour = "Coral Region Code")
Plot_CorDRYTORTUGAS_AnyCoral
#print to a pdf
pdf("CorallmePlot_AnyCoral_DRYTORTUGASRegion.pdf", width = 16, height = 4)
print(Plot_CorDRYTORTUGAS_AnyCoral)
dev.off()

#TABLES
#non-summarized lme for broom::tidy()
CorDRYTORTUGAS_lmelist2no_summary<-lapply(unique(CorDRYTORTUGAS$Coral_Region), 
                                          function(i)(lme(Perccov ~ Year,
                                                          random =~1|ProgramLocationID,
                                                          na.action = na.omit, data = subset(CorDRYTORTUGAS, Coral_Region== i)) ))

names(CorDRYTORTUGAS_lmelist2no_summary) <- unique(CorDRYTORTUGAS$Coral_Region)
#

CorDRYTORTUGAS_lmelist2summary <-  CorDRYTORTUGAS_AnyCoral_lmelist2 %>% 
  enframe() %>% 
  group_by(name) %>% 
  mutate(`Log-restricted-likelihood` =  map(.x = value, .f = ~.x$logLik)) %>% #For each species in list, pull info from each lmeObject summary
  unnest(`Log-restricted-likelihood`) %>% 
  mutate(`Number of Observations` = map (.x = value, .f = ~.x$dims$N)) %>% 
  unnest(`Number of Observations`) %>% 
  mutate(`Number of Groups` = map (.x = value, .f = ~rownames_to_column(as.data.frame(.x$dims$ngrps)))) %>% 
  unnest(`Number of Groups`) %>% 
  filter(rowname == "ProgramLocationID") %>% #specify this since there are also X and Y columns to choose from, but not needed here
  select(-rowname) %>% 
  rename(`Number of Groups` = `.x$dims$ngrps`) %>% 
  mutate(`stats` = map(.x = value, .f = ~rownames_to_column(as.data.frame(.x$tTable)))) %>%
  unnest(`stats`) %>% 
  rename(RegionID = name) %>% #more cleaning and formatting at this point
  rename(`Effect of Year` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)
#
###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject

CorDRYTORTUGAS_lmelist_tidy <- map(.x = CorDRYTORTUGAS_lmelist2no_summary, .f = broom.mixed::tidy) %>%  
  enframe() %>% 
  unnest(cols = value) %>%
  pivot_wider(names_from = term, values_from =  estimate) %>%  
  select(-c(std.error,
            df,
            statistic,
            p.value,
            effect,
            group)) %>% 
  group_by(name) %>% 
  summarise(across(.cols = c(`(Intercept)`,
                             Year,
                             `sd_(Intercept)`,
                             # `cor_Year.(Intercept)`, ****
                             # sd_Year, ****
                             sd_Observation),
                   .fns = ~sum(.x, na.rm = T))) %>% 
  rename(RegionID = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

#
###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(CorDRYTORTUGAS_lmelist_tidy, CorDRYTORTUGAS_lmelist2summary, by = "RegionID") %>% 
  select(-value) %>% 
  filter(!Type == "(Intercept)") %>% #Need to select just the year rows for Fixed Effects output stats (e.g. p-value, etc.)
  select(-Type) %>% #Past this point is just cleanup and rearragning columns
  rename(`Random Effects Intercept (SD)`= `sd_(Intercept)`) %>% 
  # rename(`Random Effects Year (SD)` = sd_Year) %>% ****
  rename(`Random Effects Residual (SD)` = sd_Observation) %>% 
  # rename(`Random Effects Year (corr)` = `cor_Year.(Intercept)`) %>% ****
  rename(FE_year = Year) %>% 
  rename(FE_intercept = `(Intercept)`) %>% 
  rename(`Fixed Effects Intercept` = FE_intercept) %>% 
  relocate(c(`Random Effects Intercept (SD)`,
             # `Random Effects Year (SD)`, ****
             `Random Effects Residual (SD)`
             # `Random Effects Year (corr)` ****
  ),
  .after = `p-value`) %>% 
  relocate(c(`Log-restricted-likelihood`,
             `Number of Observations`,
             `Number of Groups`),
           .after = `Random Effects Residual (SD)`) %>% #**** `Random Effects Year (corr)` to `Random Effects Residual (SD)`
  select(-FE_year) %>% 
  rename(`Fixed Effects Std.Error` = Std.Error) %>% 
  rename(`Fixed Effects df` = df) %>% 
  rename(`Fixed Effects t-value` = `t-value`) %>% 
  rename(`Fixed Effects p-value` = `p-value`)
library(gt)
joined_output_table %>% #check the table
  gt()

#in the following, if you are using Managed Area, change "Region and Level" to Managed Area"
nrow(joined_output_table)
joined_output_tablevec<- c("Dry Tortugas Region, All Corals")            

joined_output_table["Region and Level"] <- joined_output_tablevec   
joined_output_table<-joined_output_table%>%
  relocate("Region and Level",.before = "RegionID")
joined_output_table
table1<-joined_output_table

write.csv(joined_output_table,"Coral_lmeResults_AnyCoral_DRYTORTUGASRegion.csv")
