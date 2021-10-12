#FLORIDA SEACAR
#Wetlands Percent Cover by group- all, then by MA
#Script written by Katie May Laumann, klaumann@umces.edu, on 19 November 2020
#Script last modified: 19 Sept 2021
#Script last modified by:Katie May Laumann

#clear the environment
rm(list = ls())

library(hrbrthemes)
#first, load packages
library(tidyverse)
library(dplyr)

#read in subsetted file that the above code produced and saved
wetlands <- read.csv("~/Downloads/All Parameters but Hecatres-2021-Jul-26.csv")
colnames(wetlands)

#colnames(wetlands) [23]<-"perccov"
colnames(wetlands)[colnames(wetlands) == "X.PercentCover.SpeciesComposition_.."] <- "perccov"

unique(wetlands$ProgramID)
unique(wetlands$Region)

unique(wetlands$SpeciesGroup1)
wetlands<-
  wetlands%>%
  filter(SpeciesGroup1=="Mangroves and associate"|SpeciesGroup1=="Marsh succulents"|SpeciesGroup1=="Marsh"|SpeciesGroup1=="Invasive")

wetlands<-
  wetlands%>%
  filter(Month!="NA")

wetlands<-
  wetlands%>%
  filter(Year!="NA")

#You may decide to remove certain programs, for example here we remove 
#program 5015, which only looks at a specific site, and program 651, at request of SMEs
wetlands<-
  wetlands%>%
  filter(ProgramID!="5015")
unique(wetlands$ProgramID)

wetlands<-
  wetlands%>%
  filter(ProgramID!="651")


#remove duplicates and percent cover 0s
class(wetlands$MADup)
wetlands<-
  wetlands%>%
  filter(MADup==1)

wetlands<-
  wetlands%>%
  filter(perccov!=0)

wetlands$mth<-as.numeric(wetlands$Month)
class(wetlands$Year)
wetlands$yr<-as.numeric(wetlands$Year)

#BUILD A SUMMARY STATS TABLES with: # samples, min and max values, median, avg, and stdev of values, and program IDs
summarywetlands<-
  wetlands%>%
  group_by(ManagedArea, SpeciesGroup1, Year)%>%
  summarize(n=(length(unique(RowID))),min=min(perccov),max=max(perccov),median=median(perccov),mean=mean(perccov),sd=sd(perccov),ProgramIDs=paste(unique(ProgramID), collapse = ","))
colnames(summarywetlands)
#rename columns in table we just created
colnames(summarywetlands) [1] <- "Wetlands Managed Area" 
colnames(summarywetlands) [2] <- "Wetlands Group" 
colnames(summarywetlands) [4] <- "Number of Samples" 
colnames(summarywetlands) [5] <- "Minimum Value (Percent Cover)" 
colnames(summarywetlands) [6] <- "Maximum Value (Percent Cover)" 
colnames(summarywetlands) [7] <- "Median (Percent Cover)" 
colnames(summarywetlands) [8] <- "Mean (Percent Cover)" 
colnames(summarywetlands) [9] <- "Standard Deviation (Percent Cover)" 
colnames(summarywetlands) [10] <- "Programs with Data" 
#write a summary stats .csv
write.csv(summarywetlands,"WetlandsSummaryStatisticsTable_PercentCover_AllWetlands.csv")

##############################
#Analyses for GROUP 1 LEVEL
#1. Summary statistics plots for full data set
par(mfrow=c(2,2))
plot(lm(wetlands$Year ~ wetlands$perccov))

#2. Boxplots
#BOXPLOTS
wetlandsboxplots<-wetlands# if you are going by MA, you will specifiy the desired MA here: %>% filter(ManagedAreaName=="")
boxplotswetlands<- lapply(unique(wetlands$SpeciesGroup1), function (i) {
  dat <- filter(wetlandsboxplots, SpeciesGroup1 == i)
  ggplot(data = dat, aes(group=Year, x = Year, y = perccov)) +
    geom_boxplot() +
    labs(title = "Percent cover by group, all wetlands",subtitle=i, y="Percent Cover", x="Year") })
#view
boxplotswetlands

#view
boxplotswetlands
#print to file
pdf("PercentcoverBoxplots_AllWetlands.pdf", width = 16, height = 4) 
print(boxplotswetlands)
dev.off()

#3. Linear Mixed Effects Models
#two options to run lmes. Use lme4 OR nlme. Here, lme4 is hashtagged out, and we are running nlme.
#library(lme4)
#test1=lmer(perccov ~ Year + (1|ProgramLocationID),data=wetlands)
#summary(test1)

#to run just one lme, do the following
library(nlme)
#wetlands_ByGroup1<-lme(perccov ~ Year,
#                     random =~1|ProgramLocationID,
#                     na.action = na.omit,
#                     data = wetlands)
#summary(wetlands_ByGroup1)

#Here is where you examine your options for setting random variables and modifying your model.
#THis may, for example, depend on the # of programs sampling, program location IDs, etc.  
#To run multiple lmes at different levels, for example by SpeciesGroup1, do the following:

#Here, we will allow for variation by program location ID. 
#if there are not enough data for a specific group to run analysis, remove that group
#
wetlands_ByGroup1_lmelist2<-lapply(unique(wetlands$SpeciesGroup1), 
                                   function(i)summary(lme(perccov ~ Year,
                                                          random =~1|ProgramLocationID,
                                                          na.action = na.omit, data = subset(wetlands, SpeciesGroup1 == i))))

#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
names(wetlands_ByGroup1_lmelist2) <- unique(wetlands$SpeciesGroup1)
wetlands_ByGroup1_lmelist2

#print to a csv file
sink("Wetlandslme_ByGroup1_AllWetlands.csv")
wetlands_ByGroup1_lmelist2
sink()

#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(wetlands_ByGroup1_lmelist2), function(i) {
  mod <- wetlands_ByGroup1_lmelist2[[i]]
  pred <- data.frame(SpeciesGroup1 = i, Year = mod$data$Year, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_wetlands_ByGroup1<-ggplot() +
  geom_jitter(data = wetlands, aes(x = Year, y = perccov, col = SpeciesGroup1),size=0.1) +
  geom_line(data = preds, aes(x = Year, y = pred, group = SpeciesGroup1, col = SpeciesGroup1))
#label the plot and y axis (x was already labeled above)
Plot_wetlands_ByGroup1<-Plot_wetlands_ByGroup1 + labs(y = "Percent Cover")
Plot_wetlands_ByGroup1<-Plot_wetlands_ByGroup1 + labs(title = "Percent cover by group, All Wetlands Data")
#color code your points and lines by taxon
Plot_wetlands_ByGroup1<-Plot_wetlands_ByGroup1 + labs(colour = "Species Group")
Plot_wetlands_ByGroup1
#print to a pdf
pdf("WetlandslmePlot_ByGroup1_AllWetlands.pdf", width = 16, height = 4)
print(Plot_wetlands_ByGroup1)
dev.off()

##
#TABLES
#non-summarized lme for broom::tidy()
wetlands_lmelist2no_summary<-lapply(unique(wetlands$SpeciesGroup1), 
                                    function(i)(lme(perccov ~ Year,
                                                    random =~1|ProgramLocationID,
                                                    na.action = na.omit, data = subset(wetlands, SpeciesGroup1== i)) ))

names(wetlands_lmelist2no_summary) <- unique(wetlands$SpeciesGroup1)
#

wetlands_lmelist2summary <-  wetlands_ByGroup1_lmelist2 %>% 
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
  rename(Group1 = name) %>% #more cleaning and formatting at this point
  rename(`Effect of Year` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)
#
###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
wetlands_lmelist_tidy <- map(.x = wetlands_lmelist2no_summary, .f = broom.mixed::tidy) %>%  
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
  rename(Group1 = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

#
###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(wetlands_lmelist_tidy, wetlands_lmelist2summary, by = "Group1") %>% 
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

write.csv(joined_output_table,"Wetlands_lmeResults_ByGroup1_AllWetlands.csv")

##NOW RUN ANALYSES BY MA. for those with >=5 years of data
unique(wetlands$ManagedArea)
nyears<-
  wetlands%>%
  group_by(ManagedArea)%>%
  summarise(length(unique(Year)))
colnames(nyears) [2]<-"nyrs"
nyears

#only run analyses for Guana River Marsh AND Guana Tolomato Matanzas NERR
#Make a new df for each
wetlands1<-wetlands%>%filter(ManagedArea=="Guana River Marsh")
wetlands2<-wetlands%>%filter(ManagedArea=="Guana Tolomato Matanzas NERR")

#BUILD A SUMMARY STATS TABLES with: # samples, min and max values, median, avg, and stdev of values, and program IDs
summarywetlands1<-
  wetlands1%>%
  group_by(ManagedArea, SpeciesGroup1, Year)%>%
  summarize(n=(length(unique(RowID))),min=min(perccov),max=max(perccov),median=median(perccov),mean=mean(perccov),sd=sd(perccov),ProgramIDs=paste(unique(ProgramID), collapse = ","))
colnames(summarywetlands1)
#rename columns in table we just created
colnames(summarywetlands1) [1] <- "Wetlands Managed Area" 
colnames(summarywetlands1) [2] <- "Wetlands Group" 
colnames(summarywetlands1) [4] <- "Number of Samples" 
colnames(summarywetlands1) [5] <- "Minimum Value (Percent Cover)" 
colnames(summarywetlands1) [6] <- "Maximum Value (Percent Cover)" 
colnames(summarywetlands1) [7] <- "Median (Percent Cover)" 
colnames(summarywetlands1) [8] <- "Mean (Percent Cover)" 
colnames(summarywetlands1) [9] <- "Standard Deviation (Percent Cover)" 
colnames(summarywetlands1) [10] <- "Programs with Data" 
#write a summary stats .csv
write.csv(summarywetlands1,"WetlandsSummaryStatisticsTable_PercentCover_GuanaRiverMarsh.csv")

##############################
#Analyses for GROUP 1 LEVEL
#1. Summary statistics plots for full data set
par(mfrow=c(2,2))
plot(lm(wetlands1$Year ~ wetlands1$perccov))

#2. Boxplots
#BOXPLOTS
wetlands1boxplots<-wetlands1# if you are going by MA, you will specifiy the desired MA here: %>% filter(ManagedAreaName=="")
boxplotswetlands1<- lapply(unique(wetlands1$SpeciesGroup1), function (i) {
  dat <- filter(wetlands1boxplots, SpeciesGroup1 == i)
  ggplot(data = dat, aes(group=Year, x = Year, y = perccov)) +
    geom_boxplot() +
    labs(title = "Percent cover by group, all wetlands1",subtitle=i, y="Percent Cover", x="Year") })
#view
boxplotswetlands1

#view
boxplotswetlands1
#print to file
pdf("PercentcoverBoxplots_GuanaRiverMarsh.pdf", width = 16, height = 4) 
print(boxplotswetlands1)
dev.off()

#3. Linear Mixed Effects Models
#two options to run lmes. Use lme4 OR nlme. Here, lme4 is hashtagged out, and we are running nlme.
#library(lme4)
#test1=lmer(perccov ~ Year + (1|ProgramLocationID),data=wetlands1)
#summary(test1)

#to run just one lme, do the following
library(nlme)
#wetlands1_ByGroup1<-lme(perccov ~ Year,
#                     random =~1|ProgramLocationID,
#                     na.action = na.omit,
#                     data = wetlands1)
#summary(wetlands1_ByGroup1)

#Here is where you examine your options for setting random variables and modifying your model.
#THis may, for example, depend on the # of programs sampling, program location IDs, etc.  
#To run multiple lmes at different levels, for example by SpeciesGroup1, do the following:

#Here, we will allow for variation by program location ID. 
#if there are not enough data for a specific group to run analysis, remove that group
#
wetlands1_ByGroup1_lmelist2<-lapply(unique(wetlands1$SpeciesGroup1), 
                                    function(i)summary(lme(perccov ~ Year,
                                                           random =~1|ProgramLocationID,
                                                           na.action = na.omit, data = subset(wetlands1, SpeciesGroup1 == i))))

#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
names(wetlands1_ByGroup1_lmelist2) <- unique(wetlands1$SpeciesGroup1)
wetlands1_ByGroup1_lmelist2

#print to a csv file
sink("Wetlandslme_ByGroup1_GuanaRiverMarsh.csv")
wetlands1_ByGroup1_lmelist2
sink()

#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(wetlands1_ByGroup1_lmelist2), function(i) {
  mod <- wetlands1_ByGroup1_lmelist2[[i]]
  pred <- data.frame(SpeciesGroup1 = i, Year = mod$data$Year, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_wetlands1_ByGroup1<-ggplot() +
  geom_jitter(data = wetlands1, aes(x = Year, y = perccov, col = SpeciesGroup1),size=0.1) +
  geom_line(data = preds, aes(x = Year, y = pred, group = SpeciesGroup1, col = SpeciesGroup1))
#label the plot and y axis (x was already labeled above)
Plot_wetlands1_ByGroup1<-Plot_wetlands1_ByGroup1 + labs(y = "Percent Cover")
Plot_wetlands1_ByGroup1<-Plot_wetlands1_ByGroup1 + labs(title = "Percent cover by group, Guana River Marsh Data")
#color code your points and lines by taxon
Plot_wetlands1_ByGroup1<-Plot_wetlands1_ByGroup1 + labs(colour = "Species Group")
Plot_wetlands1_ByGroup1
#print to a pdf
pdf("WetlandslmePlot_ByGroup1_GuanaRiverMarsh.pdf", width = 16, height = 4)
print(Plot_wetlands1_ByGroup1)
dev.off()

##
#TABLES
#non-summarized lme for broom::tidy()
wetlands1_lmelist2no_summary<-lapply(unique(wetlands1$SpeciesGroup1), 
                                     function(i)(lme(perccov ~ Year,
                                                     random =~1|ProgramLocationID,
                                                     na.action = na.omit, data = subset(wetlands1, SpeciesGroup1== i)) ))

names(wetlands1_lmelist2no_summary) <- unique(wetlands1$SpeciesGroup1)
#

wetlands1_lmelist2summary <-  wetlands1_ByGroup1_lmelist2 %>% 
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
  rename(Group1 = name) %>% #more cleaning and formatting at this point
  rename(`Effect of Year` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)
#
###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
wetlands1_lmelist_tidy <- map(.x = wetlands1_lmelist2no_summary, .f = broom.mixed::tidy) %>%  
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
  rename(Group1 = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

#
###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(wetlands1_lmelist_tidy, wetlands1_lmelist2summary, by = "Group1") %>% 
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

write.csv(joined_output_table,"Wetlands_lmeResults_ByGroup1_GuanaRiverMarsh.csv")

#BUILD A SUMMARY STATS TABLES with: # samples, min and max values, median, avg, and stdev of values, and program IDs
summarywetlands2<-
  wetlands2%>%
  group_by(ManagedArea, SpeciesGroup1, Year)%>%
  summarize(n=(length(unique(RowID))),min=min(perccov),max=max(perccov),median=median(perccov),mean=mean(perccov),sd=sd(perccov),ProgramIDs=paste(unique(ProgramID), collapse = ","))
colnames(summarywetlands2)
#rename columns in table we just created
colnames(summarywetlands2) [1] <- "Wetlands Managed Area" 
colnames(summarywetlands2) [2] <- "Wetlands Group" 
colnames(summarywetlands2) [4] <- "Number of Samples" 
colnames(summarywetlands2) [5] <- "Minimum Value (Percent Cover)" 
colnames(summarywetlands2) [6] <- "Maximum Value (Percent Cover)" 
colnames(summarywetlands2) [7] <- "Median (Percent Cover)" 
colnames(summarywetlands2) [8] <- "Mean (Percent Cover)" 
colnames(summarywetlands2) [9] <- "Standard Deviation (Percent Cover)" 
colnames(summarywetlands2) [10] <- "Programs with Data" 
#write a summary stats .csv
write.csv(summarywetlands2,"WetlandsSummaryStatisticsTable_PercentCover_GuanaTolomatoMatanzasNERR.csv")

##############################
#Analyses for GROUP 1 LEVEL
#1. Summary statistics plots for full data set
par(mfrow=c(2,2))
plot(lm(wetlands2$Year ~ wetlands2$perccov))

#2. Boxplots
#BOXPLOTS
wetlands2boxplots<-wetlands2# if you are going by MA, you will specifiy the desired MA here: %>% filter(ManagedAreaName=="")
boxplotswetlands2<- lapply(unique(wetlands2$SpeciesGroup1), function (i) {
  dat <- filter(wetlands2boxplots, SpeciesGroup1 == i)
  ggplot(data = dat, aes(group=Year, x = Year, y = perccov)) +
    geom_boxplot() +
    labs(title = "Percent cover by group, all wetlands2",subtitle=i, y="Percent Cover", x="Year") })
#view
boxplotswetlands2

#view
boxplotswetlands2
#print to file
pdf("PercentcoverBoxplots_GuanaTolomatoMatanzasNERR.pdf", width = 16, height = 4) 
print(boxplotswetlands2)
dev.off()

#3. Linear Mixed Effects Models
#two options to run lmes. Use lme4 OR nlme. Here, lme4 is hashtagged out, and we are running nlme.
#library(lme4)
#test1=lmer(perccov ~ Year + (1|ProgramLocationID),data=wetlands2)
#summary(test1)

#to run just one lme, do the following
library(nlme)
#wetlands2_ByGroup1<-lme(perccov ~ Year,
#                     random =~1|ProgramLocationID,
#                     na.action = na.omit,
#                     data = wetlands2)
#summary(wetlands2_ByGroup1)

#Here is where you examine your options for setting random variables and modifying your model.
#THis may, for example, depend on the # of programs sampling, program location IDs, etc.  
#To run multiple lmes at different levels, for example by SpeciesGroup1, do the following:

#Here, we will allow for variation by program location ID. 
#if there are not enough data for a specific group to run analysis, remove that group
#
wetlands2_ByGroup1_lmelist2<-lapply(unique(wetlands2$SpeciesGroup1), 
                                    function(i)summary(lme(perccov ~ Year,
                                                           random =~1|ProgramLocationID,
                                                           na.action = na.omit, data = subset(wetlands2, SpeciesGroup1 == i))))

#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
names(wetlands2_ByGroup1_lmelist2) <- unique(wetlands2$SpeciesGroup1)
wetlands2_ByGroup1_lmelist2

#print to a csv file
sink("Wetlandslme_ByGroup1_GuanaTolomatoMatanzasNERR.csv")
wetlands2_ByGroup1_lmelist2
sink()

#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(wetlands2_ByGroup1_lmelist2), function(i) {
  mod <- wetlands2_ByGroup1_lmelist2[[i]]
  pred <- data.frame(SpeciesGroup1 = i, Year = mod$data$Year, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_wetlands2_ByGroup1<-ggplot() +
  geom_jitter(data = wetlands2, aes(x = Year, y = perccov, col = SpeciesGroup1),size=0.1) +
  geom_line(data = preds, aes(x = Year, y = pred, group = SpeciesGroup1, col = SpeciesGroup1))
#label the plot and y axis (x was already labeled above)
Plot_wetlands2_ByGroup1<-Plot_wetlands2_ByGroup1 + labs(y = "Percent Cover")
Plot_wetlands2_ByGroup1<-Plot_wetlands2_ByGroup1 + labs(title = "Percent cover by group, Guana Tolomato Matanzas NERR Data")
#color code your points and lines by taxon
Plot_wetlands2_ByGroup1<-Plot_wetlands2_ByGroup1 + labs(colour = "Species Group")
Plot_wetlands2_ByGroup1
#print to a pdf
pdf("WetlandslmePlot_ByGroup1_GuanaTolomatoMatanzasNERR.pdf", width = 16, height = 4)
print(Plot_wetlands2_ByGroup1)
dev.off()

##
#TABLES
#non-summarized lme for broom::tidy()
wetlands2_lmelist2no_summary<-lapply(unique(wetlands2$SpeciesGroup1), 
                                     function(i)(lme(perccov ~ Year,
                                                     random =~1|ProgramLocationID,
                                                     na.action = na.omit, data = subset(wetlands2, SpeciesGroup1== i)) ))

names(wetlands2_lmelist2no_summary) <- unique(wetlands2$SpeciesGroup1)
#

wetlands2_lmelist2summary <-  wetlands2_ByGroup1_lmelist2 %>% 
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
  rename(Group1 = name) %>% #more cleaning and formatting at this point
  rename(`Effect of Year` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)
#
###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
wetlands2_lmelist_tidy <- map(.x = wetlands2_lmelist2no_summary, .f = broom.mixed::tidy) %>%  
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
  rename(Group1 = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

#
###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(wetlands2_lmelist_tidy, wetlands2_lmelist2summary, by = "Group1") %>% 
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

write.csv(joined_output_table,"Wetlands_lmeResults_ByGroup1_GuanaTolomatoMatanzasNERR.csv")

