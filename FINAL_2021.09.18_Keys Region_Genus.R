#FLORIDA SEACAR
#Coral Scripts, Percent Cover Analyses
#Script written by Katie May Laumann, klaumann@umces.edu, on 25 June 2020
#Script last modified: 18 Sept 2021 by Katie May Laumann, klaumann@umces.edu


#This script uses mixed effects models to find the effect of time (Year) on Percent Cover 
#in each Region. 

#clear the environment
rm(list = ls())
#load packages
library(tidyverse)
library(dplyr)
#read in data
CorKEYS<-read.csv("~/Downloads/Percent Cover  - FLA KEYS-2021-Jul-26.csv")
colnames(CorKEYS)
#run all analyses on fixed and random samples independently. 
#At the time of code production, no column is available in data to do so.
#Once this is corrected, simply run the two following lines of code, and run the remaining script 
#once for each (Corfixed and Corrand).
#CorFLKEYSfixed<-CorFLKEYS%>%filter(samplingmethod=="fixed")
#CorFLKEYSrandom<-CorFLKEYS%>%filter(samplingmethod=="random")
#if a certain program (eg 62) should be left out, remove it from analyses using the following:
#CorKEYS<-CorKEYS%>%filter(ProgramID!="62")

unique(CorKEYS$ManagedAreaName)
CorKEYS<-CorKEYS%>%filter(ManagedAreaName!="NA")
#check unique values in some columns
unique(CorKEYS$SpeciesGroup1)
unique(CorKEYS$SpeciesGroup2)

#We have more than corals here: see Group 1
#We only want corals, so filter by Group 2
CorKEYS<-CorKEYS%>%filter(SpeciesGroup2=="Octocoral"|SpeciesGroup2=="Unspecifiedcoral")
#if you want to exclude Octocorals, run the next line instead:
#CorKEYS<-CorKEYS%>%filter(SpeciesGroup2=="Unspecifiedcoral")

#Current datasets do not include sampling method, which was requested. If sampling method is added to the dataset,
#seperate analyses by sampling method by first building a df for each method:
#CorFLKEYSMethod1<-CorFLKEYS%>%MethodRowName=="Method1"
#CorFLKEYSMethod2<-CorFLKEYS%>%MethodRowName=="Method2"
#then by running all of the following analyses for each method individually. 
#This could be simplified by using lapply, if desired.

#rename columns to more easily work with them
colnames(CorKEYS)[colnames(CorKEYS) == "X.PercentCover.SpeciesComposition_.."] <- "Perccov"
nrow(CorKEYS)
#remove % cover NAs
CorKEYS <- CorKEYS[!is.na(CorKEYS$Perccov), ]  
nrow(CorKEYS)

#!!! don't need the next 3 lines for Region-level analysis. You will need them for Managed Area Level analysis.
#Create a Region List that we will write to CorKEYSPerccovRegionList.csv to help in building our final results table, later
#CorKEYSRegionList<-unique(CorKEYS[c("ManagedAreaName","Region")])
#write.csv(CorKEYSRegionList,"CorKEYSPerccovRegionList.csv")
#*******************************
#check unique values
unique(CorKEYS$Perccov)
#convert to a numeric to be able to do math
CorKEYS$Perccov<-as.numeric(CorKEYS$Perccov)
CorKEYS <- CorKEYS[!is.na(CorKEYS$Perccov), ]  
unique(CorKEYS$SpeciesName)
unique(CorKEYS$GenusName)
nrow(CorKEYS)

#combine genus and spp into one column
CorKEYS$gensp<-paste(CorKEYS$GenusName, CorKEYS$SpeciesName, sep=" ")
nrow(CorKEYS)
unique(CorKEYS$gensp)
CorKEYS<-CorKEYS%>%filter(gensp!="NA NA")

#how many rows?
allsamples<-nrow(CorKEYS)
##STOPPED
unique(CorKEYS$SpeciesGroup1)
unique(CorKEYS$SpeciesGroup2)
colnames(CorKEYS)

class(CorKEYS$Year)
unique(CorKEYS$Year)

#BUILD A SUMMARY STATS TABLES with: # samples, min and max values, median, avg, and stdev of values, and program IDs
summaryKEYS<-
  CorKEYS%>%
  group_by(Coral_Region, GenusName, Year)%>%
  summarize(n=(length(unique(RowID))),min=min(Perccov),max=max(Perccov),median=median(Perccov),mean=mean(Perccov),sd=sd(Perccov),ProgramIDs=paste(unique(ProgramID), collapse = ","))
colnames(summaryKEYS)
#rename columns in table we just created
colnames(summaryKEYS) [1] <- "Coral Region" 
colnames(summaryKEYS) [2] <- "Genus" 
colnames(summaryKEYS) [4] <- "Number of Samples" 
colnames(summaryKEYS) [5] <- "Minimum Value (Percent Cover)" 
colnames(summaryKEYS) [6] <- "Maximum Value (Percent Cover)" 
colnames(summaryKEYS) [7] <- "Median (Percent Cover)" 
colnames(summaryKEYS) [8] <- "Mean (Percent Cover)" 
colnames(summaryKEYS) [9] <- "Standard Deviation (Percent Cover)" 
colnames(summaryKEYS) [10] <- "Programs with Data" 
#write a summary stats .csv
write.csv(summaryKEYS,"CoralSummaryStatisticsTable_Genus_KEYSRegion.csv")


##############################
#Analyses for GENUS LEVEL
#1. Summary statistics plots 
par(mfrow=c(2,2))
plot(lm(CorKEYS$Year ~ CorKEYS$Perccov))


#2. Boxplots
#BOXPLOTS
CorKEYSboxplots<-CorKEYS# if you are going by MA, you will specifiy the desired MA here: %>% filter(ManagedAreaName=="")
boxplotsCorKEYS<- lapply(unique(CorKEYS$GenusName), function (i) {
  dat <- filter(CorKEYSboxplots, GenusName == i)
  ggplot(data = dat, aes(group=Year, x = Year, y = Perccov)) +
    geom_boxplot() +
    labs(title = "Percent cover coral by genus, Florida Keys Region",subtitle=i, y="Percent Cover", x="Year") })
#view
boxplotsCorKEYS
#print to file
pdf("CoralBoxplots_ByGenus_KEYSRegion.pdf", width = 16, height = 4) 
print(boxplotsCorKEYS)
dev.off()

#3. Linear Mixed Effects Models
#two options to run lmes. Use lme4 OR nlme. Here, lme4 is hashtagged out, and we are running nlme.
#library(lme4)
#test1=lmer(Perccov ~ Year + (1|ProgramLocationID),data=CorKEYS)
#summary(test1)

#Here is where you examine your options for setting random variables and modifying your model.
#This may, for example, depend on the # of programs sampling, program location IDs, etc.  
#To run multiple lmes at different levels, for example by SpeciesGroup1, do the following:

#Here, we will allow for variation by program location ID. 
#if there are not enough data for a specific group to run analysis, remove that group
##First, ID which ones have fewer than 10 yrs of data
unique(CorKEYS$GenusName)
nyears<-
  CorKEYS%>%
  group_by(GenusName)%>%
  summarise(length(unique(Year)))
colnames(nyears) [2]<-"nyrs"
eliminate<-
  nyears%>%
  filter(nyrs<10)
eliminate
colnames(CorKEYS)
#use anti_join to include only programs>=10 yrs of data
CorKEYS<-anti_join(CorKEYS,eliminate,by=c("GenusName"))
unique(CorKEYS$GenusName)
colnames(CorKEYS)
#to run just one lme, do the following
library(nlme)
CorKEYS_ByGenus<-lme(Perccov ~ Year,
                     random =~1|ProgramLocationID,
                     na.action = na.omit,
                     data = CorKEYS)
summary(CorKEYS_ByGenus)

#to run multiple lmes by group, for example my GenusName, do the following
#Use the same type of model as above, but run it for each taxon (specified as i in the script)
#and print them to a list using lapply.
#Program location ID will be allowed to vary
#2 genera don't have enough data; remove them before analyses
CorKEYS<-CorKEYS%>%filter(GenusName!="Eusmilia")
CorKEYS<-CorKEYS%>%filter(GenusName!="Meandrina")

CorKEYS_ByGenus_lmelist2<-lapply(unique(CorKEYS$GenusName), 
                                 function(i)summary(lme(Perccov ~ Year,
                                                        random =~1|ProgramLocationID,
                                                        na.action = na.omit, data = subset(CorKEYS, GenusName == i))))

#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
names(CorKEYS_ByGenus_lmelist2) <- unique(CorKEYS$GenusName)
CorKEYS_ByGenus_lmelist2

#print to a csv file
sink("Corallme_ByGenus_KEYSRegion.csv")
CorKEYS_ByGenus_lmelist2
sink()

#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(CorKEYS_ByGenus_lmelist2), function(i) {
  mod <- CorKEYS_ByGenus_lmelist2[[i]]
  pred <- data.frame(GenusName = i, Year = mod$data$Year, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_CorKEYS_ByGenus<-ggplot() +
  geom_jitter(data = CorKEYS, aes(x = Year, y = Perccov, col = GenusName),size=0.1) +
  geom_line(data = preds, aes(x = Year, y = pred, group = GenusName, col = GenusName))
#label the plot and y axis (x was already labeled above)
Plot_CorKEYS_ByGenus<-Plot_CorKEYS_ByGenus + labs(y = "Percent Cover")
Plot_CorKEYS_ByGenus<-Plot_CorKEYS_ByGenus + labs(title = "Percent cover coral by genus, Florida Keys Region")
#color code your points and lines by taxon
Plot_CorKEYS_ByGenus<-Plot_CorKEYS_ByGenus + labs(colour = "GenusName")
Plot_CorKEYS_ByGenus
#print to a pdf
pdf("CorallmePlot_ByGenus_KEYSRegion.pdf", width = 16, height = 4)
print(Plot_CorKEYS_ByGenus)
dev.off()

##
##
##
#TABLES
#non-summarized lme for broom::tidy()
CorKEYS_lmelist2no_summary<-lapply(unique(CorKEYS$GenusName), 
                                   function(i)(lme(Perccov ~ Year,
                                                   random =~1|ProgramLocationID,
                                                   na.action = na.omit, data = subset(CorKEYS, GenusName== i)) ))

names(CorKEYS_lmelist2no_summary) <- unique(CorKEYS$GenusName)
#

CorKEYS_lmelist2summary <-  CorKEYS_ByGenus_lmelist2 %>% 
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
  rename(Genus = name) %>% #more cleaning and formatting at this point
  rename(`Effect of Year` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)

###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
CorKEYS_lmelist_tidy <- map(.x = CorKEYS_lmelist2no_summary, .f = broom.mixed::tidy) %>%  
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
  rename(Genus = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(CorKEYS_lmelist_tidy, CorKEYS_lmelist2summary, by = "Genus") %>% 
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


nrow(joined_output_table)
joined_output_tablevec<- c("Florida Keys","Florida Keys","Florida Keys","Florida Keys",
                           "Florida Keys","Florida Keys","Florida Keys","Florida Keys",
                           "Florida Keys","Florida Keys","Florida Keys","Florida Keys")            

joined_output_table["Region"] <- joined_output_tablevec   
joined_output_table<-joined_output_table%>%
  relocate("Region",.before = "Genus")
joined_output_table
table2<-joined_output_table

write.csv(joined_output_table,"Coral_lmeResults_ByGenus_KEYSRegion.csv")
