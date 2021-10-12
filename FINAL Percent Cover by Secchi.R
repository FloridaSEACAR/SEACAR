#FLORIDA SEACAR
#SAV Scripts, % cover by secchi
#Script written by Katie May Laumann, klaumann@umces.edu, on 25 June 2020
#Script last modified: 19 Aug 2021
#Script last modified: 20 Sept 2021 by Katie May Laumann, klaumann@umces.edu

#This script uses mixed effects models to find the effect of secchi on Percent Cover 
#in each Managed Area. 

#clear the environment
rm(list = ls())
#load packages
library(tidyverse)
library(dplyr)

#read in data
SAV<-read.csv("~/Desktop/SAVforlme.csv")
colnames(SAV)

secchi<-read.csv("~/Desktop/secchiforlm.csv")
colnames(secchi)

df <- inner_join(SAV, secchi, by = c("Region", "ManagedAreaName","Year","Month","yearmonth"))  # Applying inner_join() function
df       
colnames(df)
colnames(df) [7]<-"PercentCover"
colnames(df) [9]<-"secchi"

#MODEL
#we're going to look at linear trends, running an lme as we have done in many analyses previously. 
#We will allow for variation on Region; this will apply when you look at MAs that occur in multiple regions
#(with the full dataset, rather than the split ones we are using here)
unique(df$ManagedAreaName)
##First, ID which ones have fewer than 10 yrs of data
dfyrs<-
  df%>%
  group_by(ManagedAreaName)%>%
  summarise(length(unique(Year)))
colnames(dfyrs) [2]<-"numyrs"
eliminate<-
  dfyrs%>%
  filter(numyrs<10)
eliminate
colnames(eliminate)
#use anti_join to include only programs>=10 yrs of data
df<-anti_join(df,eliminate,by=c("ManagedAreaName"))
unique(df$ManagedAreaName)
library(nlme)
colnames(df)


model<-lapply(unique(df$ManagedAreaName), 
               function(i)summary(lme(PercentCover ~ secchi,
                                      random =~1|Year,
                                      na.action = na.omit, data = subset(df, ManagedAreaName == i))))

names(model) <- unique(df$ManagedAreaName)
model
#NOTE: models are overfitted. We need more data.

#we're going to look at linear trends, running an lme as we have done in many analyses previously. 
#We will allow for variation on Year; this will apply when you look at MAs that occur in multiple regions
#(with the full dataset, rather than the split ones we are using here)
library(nlme)
dfmodel<-lapply(unique(df$ManagedAreaName), 
                function(i)summary(lme(PercentCover~ secchi,
                                       random =~1|Year,
                                       na.action = na.omit, data = subset(df, ManagedAreaName == i))))

names(dfmodel) <- unique(df$ManagedAreaName)
dfmodel

#print to a csv file
sink("secchiandSAVmodel.csv")
dfmodel
sink()

#####
#testing graphs
#4. Graphs
#first, convert your list of analyses to a dataframe (named preds here)
preds <- do.call(rbind, lapply(names(dfmodel), function(i) {
  mod <- dfmodel[[i]]
  pred <- data.frame(ManagedAreaName = i, secchi = mod$data$secchi, pred = predict(mod, level = 0))
}))
#using that dataframe, plot, for each taxon, the observed values in each programlocationID
#for each year. Plot these as 'jittered' points. 
#create a line that is the fitted model for each taxon, as well.
Plot_secchiandSAV<-ggplot() +
  geom_jitter(data = df, aes(x = secchi, y = PercentCover, col = ManagedAreaName),size=0.1) +
  geom_line(data = preds, aes(x = secchi, y = pred, group = ManagedAreaName, col = ManagedAreaName))
#label the plot and y axis (x was already labeled above)
Plot_secchiandSAV<-Plot_secchiandSAV + labs(y = "Percent Cover")
Plot_secchiandSAV<-Plot_secchiandSAV + labs(title = "Percent Cover x secchi by Managed Area")
#color code your points and lines by taxon
Plot_secchiandSAV<-Plot_secchiandSAV + labs(colour = "ManagedAreaName")
Plot_secchiandSAV
#print to a pdf
pdf("secchiandSAVmodelPlot.pdf", width = 16, height = 4)
print(Plot_secchiandSAV)
dev.off()
##
##
##
#TABLES
#non-summarized lme for broom::tidy()
dfmodelno_summary<-lapply(unique(df$ManagedAreaName), 
                          function(i)(lme(PercentCover~ secchi,
                                          random =~1|Year,
                                          na.action = na.omit, data = subset(df, ManagedAreaName== i)) ))

names(dfmodelno_summary) <- unique(df$ManagedAreaName)
#
df_lmelist2summary <-  dfmodel%>% 
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
  rename(`Effect of secchi` = Value) %>% 
  rename(df = DF) %>% 
  rename(Type = rowname)

###PART 2 of tables: Used the non_summary model to get other info

#not using summary() since it's easier to get some info with tidy() on a non-summarized lmeObject
df_lmelist_tidy <- map(.x = dfmodelno_summary, .f = broom.mixed::tidy) %>%  
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
                             secchi,
                             `sd_(Intercept)`,
                             # `cor_secchi.(Intercept)`, ****
                             # sd_secchi, ****
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
  # rename(`Random Effects secchi (SD)` = sd_secchi) %>% ****
  rename(`Random Effects Residual (SD)` = sd_Observation) %>% 
  # rename(`Random Effects secchi (corr)` = `cor_secchi.(Intercept)`) %>% ****
  rename(FE_year = secchi) %>% 
  rename(FE_intercept = `(Intercept)`) %>% 
  rename(`Fixed Effects Intercept` = FE_intercept) %>% 
  relocate(c(`Random Effects Intercept (SD)`,
             # `Random Effects secchi (SD)`, ****
             `Random Effects Residual (SD)`
             # `Random Effects secchi (corr)` ****
  ),
  .after = `p-value`) %>% 
  relocate(c(`Log-restricted-likelihood`,
             `Number of Observations`,
             `Number of Groups`),
           .after = `Random Effects Residual (SD)`) %>% #**** `Random Effects secchi (corr)` to `Random Effects Residual (SD)`
  select(-FE_year) %>% 
  rename(`Fixed Effects Std.Error` = Std.Error) %>% 
  rename(`Fixed Effects df` = df) %>% 
  rename(`Fixed Effects t-value` = `t-value`) %>% 
  rename(`Fixed Effects p-value` = `p-value`)
library(gt)
joined_output_table %>% #check the table
  gt()


write.csv(joined_output_table,"secchiandSAVmodelTable.pdf.csv")



