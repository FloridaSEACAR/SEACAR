#FLORIDA SEACAR
#Water Quality Scripts, Nekton species reichness
#Script originally written by Katie May Laumann, klaumann@umces.edu, in May 2020
#Script last modified on 19 Sept 2021 by: Katie May Laumann, klaumann@umces.edu
#QAQCd by J. Edgerton

#This script calculates and plots annual taxonomic richness of nekton

#PREPARE DATA
#clear the environment
rm(list = ls())
#load packages
library(tidyverse)
library(dplyr)
library(RColorBrewer)

#Load your data. In order to look at the data in different regions, simply load the data for the region you want to look at.
#Here, we use NE region data as an example.
Nekton<-read.csv("~/Downloads/Nekton_NE-2021-Jul-26.csv")
colnames(Nekton)
unique(Nekton$SamplingMethod1)
unique(Nekton$SamplingMethod2)
#make one column where Genus and species are combined to correctly display species names
Nekton$gensp<-paste(Nekton$GenusName, Nekton$SpeciesName, sep=" ")
nrow(Nekton)
Nekton <- Nekton[!is.na(Nekton$GenusName), ] 
Nekton <- Nekton[!is.na(Nekton$ManagedAreaName), ]
unique(Nekton$gensp)
unique(Nekton$Genus)

#we are going to calculate species diversity, aka # of species. Ideally you would have no NAs in species name.
##since we do, you have 3 options : 1) move forward and assume "NA" means some species other than the individuals we
##are actually IDing to spp, 2) remove all spp NAs, or 3) calculate generic diversity. Here we are going with option 1,
##but the code can be easily adapted for any of the other options.

#We need to calculate diversity by MA. First, we want to calculate the number of unique species per managed area per year.
#make a new column so, out of curiosity, we can see how many times each spp was seen in each MA and each Year.
Nekton$NumberPresent<-ifelse(Nekton$gensp!="NA",1,0)
unique(Nekton$NumberPresent)

#out of curiosity, calculate the number of rows of data in each region and year
Nekton2<-
  Nekton%>%
  group_by(ManagedAreaName,gensp,Year,Region)%>%
  summarise(sum(Nekton$NumberPresent))

#Now, for our analysis, we need to know the number of unique spp seen ea year in ea MA
spprichness<-
  Nekton2%>%
  group_by(ManagedAreaName,Year,Region)%>%
  summarise(length(unique(gensp)))

#Moving forward with analyses, we only want MAs with >=5 yrs of data. Let's eliminate those with insufficient data.
nyears<-
  spprichness%>%
  group_by(ManagedAreaName)%>%
  summarise(length(unique(Year)))
colnames(nyears) [2]<-"nyrs"
eliminate<-
  nyears%>%
  filter(nyears$nyrs<5)
eliminate
#use anti_join to include only programs>=10 yrs of data
spprichness<-anti_join(spprichness,eliminate,by=c("ManagedAreaName"))
unique(spprichness$ManagedAreaName)

colnames(spprichness)
colnames(spprichness) [4]<-"numspp"

#we're going to look at linear trends, running an lme as we have done in many analyses previously. 
#We will allow for variation on Region; this will apply when you look at MAs that occur in multiple regions
#(with the full dataset, rather than the split ones we are using here)
library(nlme)
spprichnessmodel<-lapply(unique(spprichness$ManagedAreaName), 
                         function(i)summary(lme(numspp ~ Year,
                                                random =~1|Region,
                                                na.action = na.omit, data = subset(spprichness, ManagedAreaName == i))))

names(spprichnessmodel) <- unique(spprichness$ManagedAreaName)
spprichnessmodel

#print to a csv file
sink("Nekton_SppRichnessbyMAmodel.csv")
spprichnessmodel
sink()

#####
#testing graphs
#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(spprichnessmodel), function(i) {
  mod <- spprichnessmodel[[i]]
  pred <- data.frame(ManagedAreaName = i, Year = mod$data$Year, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_SpeciesRichness<-ggplot() +
  geom_jitter(data = spprichness, aes(x = Year, y = numspp, col = ManagedAreaName),size=0.1) +
  geom_line(data = preds, aes(x = Year, y = pred, group = ManagedAreaName, col = ManagedAreaName))
#label the plot and y axis (x was already labeled above)
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(y = "Species Richness")
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(title = "Species Richness by Managed Area")
#color code your points and lines by taxon
Plot_SpeciesRichness<-Plot_SpeciesRichness + labs(colour = "ManagedAreaName")
Plot_SpeciesRichness
#print to a pdf
pdf("Nekton_SppRichnessbyMAmodelPlot.pdf", width = 16, height = 4)
print(Plot_SpeciesRichness)
dev.off()
##
##
##
#TABLES
#non-summarized lme for broom::tidy()
spprichnessmodelno_summary<-lapply(unique(spprichness$ManagedAreaName), 
                                   function(i)(lme(numspp~ Year,
                                                   random =~1|Region,
                                                   na.action = na.omit, data = subset(spprichness, ManagedAreaName== i)) ))

names(spprichnessmodelno_summary) <- unique(spprichness$ManagedAreaName)
#
spprichness_lmelist2summary <-  spprichnessmodel%>% 
  enframe() %>% 
  group_by(name) %>% 
  mutate(`Log-restricted-likelihood` =  map(.x = value, .f = ~.x$logLik)) %>% #For each species in list, pull info from each lmeObject summary
  unnest(`Log-restricted-likelihood`) %>% 
  mutate(`Number of Observations` = map (.x = value, .f = ~.x$dims$N)) %>% 
  unnest(`Number of Observations`) %>% 
  mutate(`Number of Groups` = map (.x = value, .f = ~rownames_to_column(as.data.frame(.x$dims$ngrps)))) %>% 
  unnest(`Number of Groups`) %>% 
  filter(rowname == "Region") %>% #specify this since there are also X and Y columns to choose from, but not needed here
  select(-rowname) %>% 
  rename(`Number of Groups` = `.x$dims$ngrps`) %>% 
  mutate(`stats` = map(.x = value, .f = ~rownames_to_column(as.data.frame(.x$tTable)))) %>%
  unnest(`stats`) %>% 
  rename(ManagedArea = name) %>% #more cleaning and formatting at this point
  rename(`Effect of Year` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)

###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
spprichness_lmelist_tidy <- map(.x = spprichnessmodelno_summary, .f = broom.mixed::tidy) %>%  
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
  rename(ManagedArea = name) %>%
  rownames_to_column() %>% 
  select(-rowname)

###PART 3: join the stats from the summary() and tidy() methods into one table to be formatted

joined_output_table <- left_join(spprichness_lmelist_tidy, spprichness_lmelist2summary, by = "ManagedArea") %>% 
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


write.csv(joined_output_table,"Nekton_SppRichnessbyMAmodelTable.pdf.csv")

