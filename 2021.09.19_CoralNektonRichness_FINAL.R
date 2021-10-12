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
CorSE<-read.csv("~/Downloads/CorSE.csv")
colnames(CorSE)
unique(CorSE$ManagedAreaName)
CorSE<-CorSE%>%filter(ManagedAreaName=="Coral ECA")

CorKeys<-read.csv("~/Downloads/CorFLKeys.csv")
colnames(CorKeys)
unique(CorKeys$ManagedAreaName)
CorKeys<-CorKeys%>%filter(ManagedAreaName=="Coral ECA"|ManagedAreaName=="Florida Keys NMS"|ManagedAreaName=="Coupon Bight")

CorDry<-read.csv("~/Downloads/CorDryTortugas.csv")
colnames(CorDry)
unique(CorDry$ManagedAreaName)
CorDry<-CorDry%>%filter(ManagedAreaName=="Florida Keys NMS")

Cor<-rbind(CorSE,CorKeys)
Cor<-rbind(Cor,CorDry)

unique(Cor$MADup)
#remove duplicates
Cor<-Cor%>%filter(MADup!=2)
colnames(Cor)
unique(Cor$ProgramName)

nrow(Cor)
unique(Cor$ManagedAreaName)
#so we will do analyses for three managed areas.
#there are multiple variables, and the datasets are not built to work together well. 
#the best we can do, aside from re-building the datasets to work together, is to
#join columns that are shared between the coral and the nekton datasets we will be using.
#
#Now that we have the coral dataset, we need to put together the nekton dataset
#The file is large, so I have preprocessed it using the following code, and now start from loading Nekton2.csv
#Nekton1<-read.csv("~/Downloads/Nekton_SE-2021-Jul-26.csv")
#colnames(Nekton1)
#unique(Nekton1$ProgramName)
#unique(Nekton1$SpeciesGroup1)
#Nekton1<-Nekton1%>%filter(SpeciesGroup1=="Grazers and reef dependent species")

#unique(Nekton1$GenusName)
#Nekton2 <- Nekton1[!is.na(Nekton1$GenusName), ]  
#unique(Nekton2$GenusName)
#write.csv(Nekton2,"Nekton2.csv")
Nekton<-read.csv("~/Downloads/Nekton2.csv")
colnames(Nekton)
colnames(Cor)
unique(Nekton$Year)

#What do both datasets share?
unique(Cor$ProgramName)
unique(Nekton$ProgramName)
#Programs are different, so location IDs will also differ.
#Looking at column names, it seems the finest analysis level we can get to is to compare the monthly spp richness and percent covers.

#first, set up the summary datasets to merge
#
#We have a database problem here. The nekton dataset only has "NULL" for month, so in this example we will not include it. The final dataset should include it,
#in which case you should replace "group_by(ManagedAreaName,Year)" with "group_by(ManagedAreaName,Year,Month)"
#We will allow this to vary by Managed Area.
spprichness<-
  Nekton%>%
  group_by(ManagedAreaName,Year)%>%
  summarise(length(unique(GenusName)))
colnames(spprichness) [3] <- "GenericRichness" 

coralperccov<-
  Cor%>%
  group_by(ManagedAreaName,Year)%>%
  summarise(mean(Perccov))
colnames(coralperccov) [3] <- "MeanPercCov" 

df <- inner_join(spprichness, coralperccov, by = c("ManagedAreaName", "Year"))

#To run just one lme at the Region level (we are in the SE region currently), 
#looking to see if there is a relationship between richness and percent cover, run the following.
#we will allow for variation around Year
library(nlme)
dfmodel<-lme(GenericRichness ~ MeanPercCov,
                    random =~1|Year,
                    na.action = na.omit,
                    data = df)
summary(dfmodel)

sink("CoralNektonlmeSERegion.csv")
dfmodel
sink()
#Breaking this down, run it at the MA level. Many more considerations will need to be taken into account once the database
#is corrected to include month and Coral Region in the Nekton dataset.

#In the following, we will allow for variation by program location ID. 
unique(df$ManagedAreaName)
df<-df%>%filter(ManagedAreaName=="Florida Keys NMS"|ManagedAreaName=="Coral ECA")
dfmodel_lmelist2<-lapply(unique(df$ManagedAreaName), 
                                function(i)summary(lme(GenericRichness ~ MeanPercCov,
                                                       random =~1|Year,
                                                       na.action = na.omit, data = subset(df, ManagedAreaName == i))))

#create and add a names list to the output, with each name corresponding to one
#unique taxon for the MA/in the model
names(dfmodel_lmelist2) <- unique(df$ManagedAreaName)
dfmodel_lmelist2

#print to a csv file
sink("CoralNektonlmebyMA.csv")
dfmodel_lmelist2
sink()

#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(dfmodel_lmelist2), function(i) {
  mod <- dfmodel_lmelist2[[i]]
  pred <- data.frame(ManagedAreaName = i, MeanPercCov = mod$data$MeanPercCov, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_SpeciesRichness<-ggplot() +
  geom_jitter(data = df, aes(x = MeanPercCov, y = GenericRichness, col = ManagedAreaName),size=0.1) +
  geom_line(data = preds, aes(x = MeanPercCov, y = pred, group = ManagedAreaName, col = ManagedAreaName))
#label the plot and y axis
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(y = "Species Richness")
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(x = "Mean Percent Cover")
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(title = "Species Richness by Managed Area")
#color code your points and lines by taxon
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(colour = "ManagedAreaName")
Plot_SpeciesRichness

#print to a pdf
pdf("NektonCoral_SppRichnessbyMAmodelPlot.pdf", width = 16, height = 4)
print(Plot_SpeciesRichness)
dev.off()

#TABLES
#non-summarized lme for broom::tidy()
dfmodel_lmelist2no_summary<-lapply(unique(df$ManagedAreaName), 
                                   function(i)(lme(GenericRichness~ MeanPercCov,
                                                   random =~1|Year,
                                                   na.action = na.omit, data = subset(df, ManagedAreaName== i)) ))

names(dfmodel_lmelist2no_summary) <- unique(df$ManagedAreaName)
#
df_lmelist2summary <-  dfmodel_lmelist2%>% 
  enframe() %>% 
  group_by(name) %>% 
  mutate(`Log-restricted-likelihood` =  map(.x = value, .f = ~.x$logLik)) %>% #For each species in list, pull info from each lmeObject summary
  unnest(`Log-restricted-likelihood`) %>% 
  mutate(`Number of Observations` = map (.x = value, .f = ~.x$dims$N)) %>% 
  unnest(`Number of Observations`) %>% 
  mutate(`Number of Groups` = map (.x = value, .f = ~rownames_to_column(as.data.frame(.x$dims$ngrps)))) %>% 
  unnest(`Number of Groups`) %>% 
  filter(rowname == "Year") %>% #specify this since there are also X and Y columns to choose from, but not needed here
  select(-rowname) %>% 
  rename(`Number of Groups` = `.x$dims$ngrps`) %>% 
  mutate(`stats` = map(.x = value, .f = ~rownames_to_column(as.data.frame(.x$tTable)))) %>%
  unnest(`stats`) %>% 
  rename(ManagedArea = name) %>% #more cleaning and formatting at this point
  rename(`Effect of MeanPercCov` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)

###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
df_lmelist_tidy <- map(.x = dfmodel_lmelist2no_summary, .f = broom.mixed::tidy) %>%  
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
                             MeanPercCov,
                             `sd_(Intercept)`,
                             # `cor_MeanPercCov.(Intercept)`, ****
                             # sd_MeanPercCov, ****
                             sd_Observation),
                   .fns = ~sum(.x, na.rm = T))) %>% 
  rename(ManagedArea = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(df_lmelist_tidy, df_lmelist2summary, by = "ManagedArea") %>% 
  select(-value) %>% 
  filter(!Type == "(Intercept)") %>% #Need to select just the year rows for Fixed Effects output stats (e.g. p-value, etc.)
  select(-Type) %>% #Past this point is just cleanup and rearragning columns
  rename(`Random Effects Intercept (SD)`= `sd_(Intercept)`) %>% 
  # rename(`Random Effects MeanPercCov (SD)` = sd_MeanPercCov) %>% ****
  rename(`Random Effects Residual (SD)` = sd_Observation) %>% 
  # rename(`Random Effects MeanPercCov (corr)` = `cor_MeanPercCov.(Intercept)`) %>% ****
  rename(FE_year = MeanPercCov) %>% 
  rename(FE_intercept = `(Intercept)`) %>% 
  rename(`Fixed Effects Intercept` = FE_intercept) %>% 
  relocate(c(`Random Effects Intercept (SD)`,
             # `Random Effects MeanPercCov (SD)`, ****
             `Random Effects Residual (SD)`
             # `Random Effects MeanPercCov (corr)` ****
  ),
  .after = `p-value`) %>% 
  relocate(c(`Log-restricted-likelihood`,
             `Number of Observations`,
             `Number of Groups`),
           .after = `Random Effects Residual (SD)`) %>% #**** `Random Effects MeanPercCov (corr)` to `Random Effects Residual (SD)`
  select(-FE_year) %>% 
  rename(`Fixed Effects Std.Error` = Std.Error) %>% 
  rename(`Fixed Effects df` = df) %>% 
  rename(`Fixed Effects t-value` = `t-value`) %>% 
  rename(`Fixed Effects p-value` = `p-value`)
library(gt)
joined_output_table %>% #check the table
  gt()


write.csv(joined_output_table,"NektonCoral_SppRichnessbyMAmodelTable.csv")

