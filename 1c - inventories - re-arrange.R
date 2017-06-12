---
  title: "FunDivSynthesis - prepare Inventory"
output: html_notebook
---
  
rm(list=ls())
gc()
memory.size(max=T)
memory.limit(size = 10000)

if(!require(data.table)){
  install.packages("data.table")
  library(data.table)}
if(!require(vegan)){
  install.packages("vegan")
  library(vegan)}
if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)}
if(!require(reshape2)){
  install.packages("reshape2")
  library(reshape2)}
if(!require(raster)){
  install.packages("raster")
  library(raster)}
if(!require(cluster)){
  install.packages("cluster")
  library(cluster)}
if(!require(optmatch)){
  install.packages("optmatch")
  library(optmatch)}
if(!require(rgdal)){
  install.packages("rgdal")
  library(rgdal)}
if(!require(ggplot2)){
  install.packages("ggplot2")
  library(ggplot2)}
if(!require(rworldmap)){
  install.packages("rworldmap")
  library(rworldmap)}
if(!require(scales)){
  install.packages("scales")
  library(scales)}



setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Inventories\\NFI with French data")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Inventories\\NFI with French data")

species_name<-data.table(read.csv("FunDivEUROPE_species.csv",na.strings=c("","NA")))
plotdata<-data.table(read.csv("FunDivEUROPE_plot_data.csv",na.strings=c("","NA")))
treedata<-data.table(read.csv("FunDivEUROPE_all_trees_10cm.csv",na.strings=c("","NA")))

#Converting plotcode to a character and removing the trailing space:
trim.trailing <- function (x) sub("\\s+$", "", x)
plotdata$plotcode <- trim.trailing (as.character(plotdata$plotcode))
treedata$plotcode <- trim.trailing (as.character(treedata$plotcode))

# remove plots that do not match between plot and tree data
plotdata = subset(plotdata, plotcode %in% treedata$plotcode)
treedata = subset(treedata, plotcode %in% plotdata$plotcode)
#-> 132281 plots remain

#remove all plots that have harvested trees (status == 3) or are managed ( management != 0)
remove.plotcode = treedata$plotcode[which(treedata$treestatus_th == 3)]
remove.plotcode = c(remove.plotcode,
                    plotdata$plotcode[which(!(plotdata$management2 %in% "0" | plotdata$management1 %in% "0"))])
stable.plots = subset(plotdata, !(plotcode %in% remove.plotcode))

# remove all trees that are dead (status = 4,5 or 6?)
stable.trees = subset(treedata, !(treestatus_th %in% c(4,5,6)))
stable.plots = subset(stable.plots, (plotcode %in% stable.trees$plotcode))

#change zeros in dhb and ba2 to NA
stable.trees$dbh1[which(stable.trees$dbh1 == 0)] = NA
stable.trees$dbh2[which(stable.trees$dbh2 == 0)] = NA
stable.trees$ba1[which(stable.trees$ba1 == 0)] = NA
stable.trees$ba2[which(stable.trees$ba2 == 0)] = NA

#filter out plots for which no growth data is available
stable.trees.with.growth = stable.trees[(!(is.na(ba_ha1 )) & !(is.na(ba_ha2))) | !(is.na(bachange_ha)) | !(is.na(bachange_ha_yr))]
stable.plots = stable.plots[plotcode %in% unique(as.character(stable.trees.with.growth$plotcode))]
stable.trees = stable.trees[plotcode %in% stable.plots$plotcode]

# check
table(stable.plots$country)

# rename forest types
forest.type.rename = read.table("european_forest_types.txt",header=T,sep="\t")
stable.plots$europeanforesttype = as.character(stable.plots$europeanforesttype)

for(i in 1:nrow(forest.type.rename)){
  stable.plots$europeanforesttype[which(stable.plots$europeanforesttype %in% as.character(forest.type.rename$forest.code[i]))] =
    as.character(forest.type.rename$forest.type_short[i])
}

# check
table(stable.plots$europeanforesttype)
table(stable.plots$country)

###{r homogenize species names}
species_rename = data.table(read.table("species_names_correction.txt",header=T,sep="\t"))
names(species_rename) = c("code_old","code_new","acceptedname")
stable.trees = merge(stable.trees, species_name[,.(id,code)], by.x="speciesid",by.y="id",all.x=T)
stable.trees = merge(stable.trees, species_rename[,.(code_old,code_new)], by.x="code",by.y="code_old",all.x=T)

stable.trees$code = NULL

#remove plots with unkown species ids
if(length(which(is.na(stable.trees$code_new))) >= 1){
  stable.trees = stable.trees[-(which(stable.trees$plotcode %in% stable.trees$plotcode[which(is.na(stable.trees$code_new))])),]
}

# check
table(stable.plots$europeanforesttype)
table(stable.plots$country)

###r calculate basal area und dbh heterogeneity for each plot}

# add years between surveys to stable.trees
stable.trees = merge(stable.trees,stable.plots[,.(plotcode,yearsbetweensurveys)],by.x="plotcode",by.y="plotcode")

#fill NA's in ba1 with 0
#stable.trees$ba1[which(is.na(stable.trees$ba1))] = 0

#not neccessary
#stable.trees$ba2[which(is.na(stable.trees$ba2))] = 0 #fill NA's in ba2 with 0
#stable.trees = stable.trees[which(stable.trees$ba2 > 0),] #remove dead trees

#calculate ba_change_per_year increment per year and tree
stable.trees$ba_change_per_year = (stable.trees$ba2 - stable.trees$ba1) / stable.trees$yearsbetweensurveys
stable.trees$ba_change_per_year[is.na(stable.trees$ba_change_per_year)] = stable.trees$bachange_ha_yr[is.na(stable.trees$ba_change_per_year)]

#per species: sum of ba per, number of trees and mean ba increment
species.ba.n.growth = stable.trees[,.(species_sum_ba = sum(ba_ha2),
                                      species_n_trees=length(ba_ha2),
                                      species_growth_per_tree = mean(ba_change_per_year)),
                                   by=.(plotcode,code_new)]

#per plot: sum of bam ba heterogeneity
plot.statistics = stable.trees[,.(plot_sum_ba = sum(ba2),
                                  plot_ba_heterogeneity = sd(ba2) / mean(ba2),
                                  plot_n_trees = length(ba2)),
                               by=.(plotcode)]
#####r classify species richness}
#build a wide-format table to classify
plot.species.ba = data.frame(acast(species.ba.n.growth, plotcode ~ code_new, value.var="species_sum_ba",fill=0))
plot.species.ba$plotcode = rownames(plot.species.ba)

classify.with.uneven.proportions = function(x){
  x = x / sum(x)
  threshold = 0.1
  x.sort = sort(x,decreasing = T)
  x.sum = sum(x)
  
  #mono
  if(x.sort[1] >= (x.sum*0.9) & x.sort[2] <= (x.sum*threshold)){return(1)}
  #two
  if((x.sort[1] + x.sort[2]) >= (x.sum*0.9) & x.sort[2] > 0.1 & x.sort[3] <= (x.sum*threshold)){return(2)}
  #three
  if((x.sort[1] + x.sort[2] + x.sort[3]) >= (x.sum*0.9) & x.sort[3] > 0.1 & x.sort[4] <= (x.sum*threshold)){return(3)}
  #four
  if((x.sort[1] + x.sort[2] + x.sort[3] + x.sort[4]) >= (x.sum*0.9) & x.sort[4] > 0.1  & x.sort[5] <= (x.sum*threshold)){return(4)}
  #five
  if((x.sort[1] + x.sort[2] + x.sort[3] + x.sort[4] + x.sort[5])  >= (x.sum*0.9) & x.sort[5] > 0.1 & x.sort[6] <= (x.sum*threshold)){return(5)}
  #six
  if((x.sort[1] + x.sort[2] + x.sort[3] + x.sort[4] + x.sort[5] + x.sort[6]) >= (x.sum*0.9) & x.sort[6] > 0.1 & x.sort[7] <= (x.sum*threshold)){return(6)}
  else{return(NA)}
}

classified.uneven.proportions = apply(plot.species.ba[,c(1:(ncol(plot.species.ba)-1))],1,classify.with.uneven.proportions)

classified.uneven.proportions = data.frame("plotcode" = names(classified.uneven.proportions),"SR"=classified.uneven.proportions)

species.ba.n.growth = merge(species.ba.n.growth,classified.uneven.proportions,by.x="plotcode",by.y="plotcode")
plot.statistics = merge(plot.statistics,classified.uneven.proportions,by.x="plotcode",by.y="plotcode")

# remove unclassified
species.ba.n.growth = species.ba.n.growth[!(is.na(SR))]
plot.statistics = plot.statistics[!(is.na(SR))]

#####r clean up: remove plots with < 5 trees and all neglegible species per plot}
plot.statistics = subset(plot.statistics, plot_n_trees >= 5)
species.ba.n.growth = subset(species.ba.n.growth, plotcode %in% plot.statistics$plotcode)

#this might take 3-4 hours
species.ba.n.growth.clean = species.ba.n.growth
species.ba.n.growth.clean$plot_paste_code =   paste(species.ba.n.growth.clean$plotcode, species.ba.n.growth.clean$code_new,sep="-")
i
pb <- txtProgressBar(min=0,max=length(unique(plot.statistics$plotcode)),initial=0,style=3)
for(i in 1:length(unique(plot.statistics$plotcode))){
  subset.temp = species.ba.n.growth.clean[plotcode %in% unique(plot.statistics$plotcode)[i]]
  subset.temp = subset.temp[order(subset.temp$species_sum_ba,decreasing=T),]
  SR.temp = subset.temp$SR[1]
  if(SR.temp < nrow(subset.temp)){
    plotcode.to.remove = subset.temp$plot_paste_code[(SR.temp+1):nrow(subset.temp)]
    species.ba.n.growth.clean = 
      species.ba.n.growth.clean[-which(species.ba.n.growth.clean$plot_paste_code %in% plotcode.to.remove),]}
  setTxtProgressBar(pb, i)
}

stable.plots.clean = subset(stable.plots, plotcode %in% plot.statistics$plotcode)
stable.trees.clean = subset(stable.trees, plotcode %in% species.ba.n.growth.clean$plotcode)

#check
length(unique(stable.plots.clean$plotcode))
length(unique(stable.trees.clean$plotcode))
length(unique(plot.statistics$plotcode))
length(unique(species.ba.n.growth.clean$plotcode))

#####{r calculate diversity indices}

plot.statistics2 = species.ba.n.growth.clean[,.(
  focal_species = paste(code_new, collapse = "-"),
  focal_n_trees = sum(species_n_trees),
  focal_growth_per_tree = sum(species_n_trees*species_growth_per_tree)/sum(species_n_trees),
  focal_shannon_div = diversity(species_sum_ba, index="shannon")),
  by = .(plotcode)]

plot.statistics3 = merge(plot.statistics,plot.statistics2,by.x="plotcode",by.y="plotcode") 
plot.statistics3 = merge(plot.statistics3,plotdata[,.(plotcode,longitude, latitude)],by.x="plotcode",by.y="plotcode") 


write.table(plot.statistics3,"intermediate_output\\plot_statistics_clean_without_environment.txt",sep="\t",dec=".",row.names=F,quote=F)
write.table(stable.trees.clean,"intermediate_output\\stable_trees_clean.txt",sep="\t",dec=".",row.names=F,quote=F)
write.table(stable.plots.clean,"intermediate_output\\stable_plots_clean.txt",sep="\t",dec=".",row.names=F,quote=F)
write.table(species.ba.n.growth.clean,"intermediate_output\\species_statistics_clean.txt",sep="\t",dec=".",row.names=F,quote=F)


#####r add climate variables and forest type}

plot.statistics3 = data.table(read.table("intermediate_output\\plot_statistics_clean_without_environment.txt",sep="\t",dec=".",header=T))
species.ba.n.growth.clean = data.table(read.table("intermediate_output\\species_statistics_clean.txt",sep="\t",dec=".",header=T))
stable.plots.clean = data.table(read.table("intermediate_output\\stable_plots_clean.txt",sep="\t",dec=".",header=T))

dem = read.table("FunDivEUROPE_inventories_dem.csv",sep=",",header=T)
clim = read.table("FunDivEUROPE_inventories_climate.csv",sep=",",header=T)
clim$country =NULL

# climate = getData('worldclim', var='bio', res=2.5)
# bioclim = raster::extract(climate, data.frame(plot.statistics3[,.(longitude,latitude)]), method="bilinear")
# bioclim = as.data.frame(bioclim)
# plot.statistics3 = cbind(plot.statistics3,bioclim)

plot.statistics3 = merge(plot.statistics3,stable.plots.clean[,.(plotcode,europeanforesttype)],by.x="plotcode",by.y="plotcode")
plot.statistics3 = merge(plot.statistics3,dem,by.x="plotcode",by.y="plotcode")
plot.statistics4 = merge(plot.statistics3,clim,by.x="plotcode",by.y="plotcode")

species.ba.n.growth.clean = merge(species.ba.n.growth.clean,stable.plots.clean[,.(plotcode,europeanforesttype)],by.x="plotcode",by.y="plotcode")
species.ba.n.growth.clean = merge(species.ba.n.growth.clean,dem,by.x="plotcode",by.y="plotcode")
species.ba.n.growth.clean = merge(species.ba.n.growth.clean,clim,by.x="plotcode",by.y="plotcode")

write.table(plot.statistics4,"intermediate_output\\plot_statistics_clean_with_environment.txt",sep="\t",dec=".",row.names=F,quote=F)
write.table(species.ba.n.growth.clean,"intermediate_output\\species_statistics_clean_with_environment.txt",sep="\t",dec=".",row.names=F,quote=F)

#######r match monocultures and mixtures for meta-analysis}
species.ba.n.growth.clean = data.table(read.table("intermediate_output\\species_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))
plot.statistics4 = data.table(read.table("intermediate_output\\plot_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))

# remove plots with negative growth
plots.to.remove = species.ba.n.growth.clean$plotcode[which(species.ba.n.growth.clean$species_growth_per_tree <= 0 )]
species.ba.n.growth.clean = species.ba.n.growth.clean[-which(species.ba.n.growth.clean$plotcode %in% plots.to.remove),]

do.match <-function(data,species, forest_type){
  
  #construct a df to store results
  df.results.temp = data.frame("code_new" = NA, "forest_type" = NA,"monoculture" = NA, "mixture" = NA, 
                               "SR_mixture" = NA,"distance" = NA)[0,]
  
  #take the subset of data
  df.temp = as.data.frame(data[code_new == species & europeanforesttype == forest_type])
  
  # across each richness level of mixtures
  for(richness in 2:length(table(df.temp$SR))){
    
    #subset for richness level
    df.temp.richness = subset(df.temp, SR == 1 | SR == richness)
    
    #check for inf and missing values
    df.temp.richness = df.temp.richness[complete.cases(df.temp.richness),]
    
    # only continue if at least two plots in each SR class
    if(length(which(table(df.temp.richness$SR) > 1)) == 2){
      row.names(df.temp.richness) = df.temp.richness$plotcode
      
      # distance matrix
      dist.matrix = as.matrix(daisy(df.temp.richness[,c("bio1","bio4","bio12","bio15","slope","plot_ba_heterogeneity","plot_ba_heterogeneity")],metric="euclidean"))
      
      #assign a treatment variable
      df.temp.richness$treat = NA
      if(table(df.temp.richness$SR)[1] > table(df.temp.richness$SR)[2]){
        df.temp.richness$treat[which(df.temp.richness$SR == 1)] = 0
        df.temp.richness$treat[which(df.temp.richness$SR == richness)] = 1
      }else{
        df.temp.richness$treat[which(df.temp.richness$SR == 1)] = 1
        df.temp.richness$treat[which(df.temp.richness$SR == richness)] = 0}
      
      # Match algorithm
      temp.formula = as.formula(paste("treat ~ ",
                                paste(c("bio1","bio4","bio12","bio15","slope","plot_ba_heterogeneity","plot_ba_heterogeneity"),collapse = "+")))
      
      dist.temp = match_on(temp.formula, data=df.temp.richness, method="euclidean")
      match.temp = pairmatch(dist.temp, data = df.temp.richness, remove.unmatchables = TRUE)
      
      #create match list
      for(match.levels in as.character(levels(match.temp))){
        one.match = which(match.temp %in% match.levels)
        one.match.df = df.temp.richness[c(one.match),]
        one.match.df = one.match.df[order(one.match.df$SR),]
        
        match.df.temp = data.frame("code_new" = species, "forest_type" = forest_type,
                                   "monoculture" = one.match.df$plotcode[1],
                                   "mixture" = one.match.df$plotcode[2], 
                                   "SR_mixture" = richness,
                                   "distance"= dist.matrix[which(rownames(dist.matrix) %in% one.match.df$plotcode[1]),which(colnames(dist.matrix) %in% one.match.df$plotcode[2])])
        
        #append results
        df.results.temp = rbind(df.results.temp,match.df.temp)
      }
    }
  }
  return(df.results.temp)
}


# environment to determine distance = bio1, bio4, bio12, bio15, slope,sum_ba_plot, ba_heterogeneity
data = species.ba.n.growth.clean
data = merge(data, plot.statistics4[,.(plotcode,plot_sum_ba,plot_ba_heterogeneity)],
             by.x="plotcode",by.y="plotcode")
data = data[,.(plotcode,code_new,SR,europeanforesttype,
               bio1 = scale(bio1), bio4=scale(bio4), bio12=scale(bio12),bio15 = scale(bio15),slope=scale(slope),
               plot_sum_ba = scale(plot_sum_ba),plot_ba_heterogeneity=scale(plot_ba_heterogeneity))]

# get all posible species forest combinations
species_forest_type_combinations = unique(melt(species.ba.n.growth.clean[,.(code_new,europeanforesttype)], id.vars=c("code_new", "europeanforesttype")))

matched.results = data.frame("code_new" = NA, "forest_type" = NA,"monoculture" = NA, "mixture" = NA,
                             "SR_mixture" = NA,"distance" = NA)[0,]

# run the matchit algorithm - 30 min
pb = txtProgressBar(min = 0, max = nrow(species_forest_type_combinations),style = 3)
for(i in 1:nrow(species_forest_type_combinations)){
  matched.results.temp = do.match(data = data,
                                  species = species_forest_type_combinations$code_new[i],
                                  forest_type = species_forest_type_combinations$europeanforesttype[i])
  matched.results = rbind(matched.results,matched.results.temp)
  setTxtProgressBar(pb, i)
}

write.table(matched.results,"intermediate_output\\matched_for_meta-analysis.txt",sep="\t",dec=".",row.names=F,quote=F)

######r append growth, remove 10% outliers, calculated effect sizes and plot the location of the monos/mixtures}
matched.results = data.table(read.table("intermediate_output\\matched_for_meta-analysis.txt",sep="\t",dec=".",head=T))
species.ba.n.growth.clean = data.table(read.table("intermediate_output\\species_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))

# append growth to monocultures and mixtures
matched.results.merged = merge(matched.results,species.ba.n.growth.clean[,.(plotcode,code_new,country, species_growth_per_tree)],by.x=c("monoculture","code_new"),by.y=c("plotcode","code_new"))
names(matched.results.merged)[c(7,8)] =c("country_mono", "growth_mono")

matched.results.merged = merge(matched.results.merged,species.ba.n.growth.clean[,.(plotcode,code_new,country,species_growth_per_tree)],by.x=c("mixture","code_new"),by.y=c("plotcode","code_new"))
names(matched.results.merged)[c(9,10)] =c("country_mix","growth_mix")

# rename forest type
# plot distances histogram
matched.results.merged$forest_type = as.character(matched.results.merged$forest_type)
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "acidophilous")] = "Acidophilous"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "alpine")] = "Alpine"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "beech")] = "Beech"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "boreal")] = "Boreal"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "broadleaved_evergreen")] = "Broadleaved evergreen"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "floodplain")] = "Floodplain"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "hemiboreal")] = "Hemiboreal"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "introduced")] = "Introduced"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "mediterranean_coniferous")] = "Mediterranean coniferous"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "mesophytic_deciduous")] = "Mesophytic deciduous"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "mountain_beech")] = "Mountain beech"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "non_riverine")] = "Non-riverine"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "swamp")] = "Swamp"
matched.results.merged$forest_type[which(matched.results.merged$forest_type %in% "thermophilous_deciduous")] = "Thermophilous deciduous"
matched.results.merged$forest_type = factor(matched.results.merged$forest_type, levels = c(
  "Mediterranean coniferous","Thermophilous deciduous","Broadleaved evergreen","Beech","Introduced","Non-riverine",
  "Mesophytic deciduous","Floodplain","Swamp","Acidophilous","Mountain beech","Hemiboreal","Boreal","Alpine"))

# cut distances > 10%
quantile(matched.results$distance ,c(0.9))
matched.results.cut_off = subset(matched.results.merged, distance < quantile(matched.results.merged$distance ,0.9))


# plot histogram of distances ########
png(filename = "C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\manuscript 1_meta-analysis\\output\\Histogram_distances.png",
    height = 1500, width = 3000, res = 300)
ggplot() + geom_histogram(aes(x = matched.results.merged$distance, fill = matched.results.merged$forest_type), bins = 100) +
  geom_vline(aes(xintercept = 1.273256),linetype="dotted", size = 2) +
  xlab("Euclidean distance between\nmonospecific and mixed plots") +
  ylab("Count") +
  theme_bw() + theme(panel.grid = element_blank()) +
  scale_fill_manual(name = "Forest type", values = c("#D11A10","#9F5504","#D67521","#D69224","#D6AE27","#D6C92A",
                    "#7FD635","#0B7D08","#4FD5CC","#55B6D5","#5A92D5","#6071D5","#7866D5","#9D6CD5")) +
  scale_x_sqrt(breaks = c(0,0.5,1,2,4,6,8)) 
graphics.off() 

# get all potential species combinations
all.species.combinations = as.character(matched.results.cut_off$code_new)
pb <- txtProgressBar(min = 0, max = nrow(matched.results.cut_off), style = 3)

for(i in 1:nrow(matched.results.cut_off)){
  species.temp = species.ba.n.growth.clean[plotcode %in% as.character(matched.results.cut_off$mixture[i])]
  all.species.combinations = c(all.species.combinations,paste(sort(species.temp$code_new), collapse="_"))  
  setTxtProgressBar(pb, i)
}

all.species.combinations.with.frequency = as.data.frame(table(all.species.combinations))
names(all.species.combinations.with.frequency) = c("species_combinations","n_replicates")

write.table(all.species.combinations.with.frequency, "all_species_compositions.csv",sep="\t",dec=".",quote=F,row.names=F)

# calculate effect sizes - per country comparison needed
matched.results.cut_off$counties_compared = paste(matched.results.cut_off$country_mono,matched.results.cut_off$country_mix, sep="_")
inventories.effect.sizes = matched.results.cut_off[,.(
  mono_mean = mean(growth_mono), mono_sd = sd(growth_mono),mono_n = length(growth_mono),
  mix_mean = mean(growth_mix), mix_sd = sd(growth_mix), mix_n = length(growth_mix)),
  by= .(code_new,forest_type,SR_mixture,counties_compared)]

write.table(inventories.effect.sizes,"intermediate_output\\inventories_effect_sizes.txt",sep="\t",dec=".",row.names=F,quote=F)

# plot the number of comparisons within and between countries
png(filename = "C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\manuscript 1_meta-analysis\\output\\Forest_types_per_country.png",
    height = 2000, width = 3500, res = 300)
ggplot(data =matched.results.cut_off ) + 
  geom_bar(aes(x=counties_compared, fill = forest_type)) +
  scale_fill_manual(name = "Forest type", values = c("#D11A10","#9F5504","#D67521","#D69224","#D6AE27","#D6C92A",
                                                     "#7FD635","#0B7D08","#4FD5CC","#55B6D5","#5A92D5","#6071D5","#7866D5","#9D6CD5")) +
  ylab("Count") + xlab("Origin of the monospecifec and\norigin of the mixed plot that are compared")+ 
  theme_bw() + theme(panel.grid = element_blank()) 
graphics.off()

# plot the number of species comparisons within and between countries
inventories.effect.sizes.per.species = matched.results.cut_off[,.(
  mono_mean = mean(growth_mono), mono_sd = sd(growth_mono),mono_n = length(growth_mono),
  mix_mean = mean(growth_mix), mix_sd = sd(growth_mix), mix_n = length(growth_mix)),
  by= .(code_new,forest_type,counties_compared)]

png(filename = "C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\manuscript 1_meta-analysis\\output\\Forest_types_per_country_species.png",
    height = 2000, width = 3500, res = 300)
ggplot(data =inventories.effect.sizes.per.species ) + 
  geom_bar(aes(x=counties_compared, fill = forest_type)) +
  scale_fill_manual(name = "Forest type", values = c("#D11A10","#9F5504","#D67521","#D69224","#D6AE27","#D6C92A",
                                                     "#7FD635","#0B7D08","#4FD5CC","#55B6D5","#5A92D5","#6071D5","#7866D5","#9D6CD5")) +
  ylab("Count") + xlab("Origin of the monospecifec and\norigin of the mixed plot that are compared")+ 
  theme_bw() + theme(panel.grid = element_blank()) 
graphics.off()


# plotting the position of monos and mixtures ###########
data(countriesLow)
countries = countriesLow

plot.data = matched.results.cut_off[,.(monoculture,mixture,SR_mixture, growth_mono,growth_mix)]
plot.statistics4 = data.table(read.table("intermediate_output\\plot_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))
stable.plots.clean = read.table("intermediate_output\\stable_plots_clean.txt",sep="\t",dec=".",header=T)
  
plot.data = merge(plot.data,stable.plots.clean[,c("plotcode","longitude","latitude")],
                  by.x = "monoculture", by.y="plotcode")
names(plot.data)[which(names(plot.data) %in% c("longitude","latitude"))] = c("long_mono","lat_mono")
plot.data = merge(plot.data,stable.plots.clean[,c("plotcode","longitude","latitude")],
                  by.x = "mixture", by.y="plotcode")
names(plot.data)[which(names(plot.data) %in% c("longitude","latitude"))] = c("long_mix","lat_mix")

setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\manuscript 1_meta-analysis\\output")
png(filename= "Appendix_map_monos2.png",height = 1500, width = 1500, res= 300)
ggplot(data=countries) + geom_polygon(aes(x=long,y=lat,group = group),fill="white",color="grey",size=0.3) +
  coord_cartesian(xlim=c(-20,30),ylim=c(25,70)) +
  geom_point(data=plot.data[!duplicated(plot.data$monoculture),],
             aes(x=long_mono,y=lat_mono),alpha=0.15,size=0.6) +
  xlab("Longitude") + ylab("Latitude") + ggtitle("Monoculture plots")
graphics.off()

png(filename= "Appendix_map_2_mix2.png",height = 1500, width = 1500, res= 300)
ggplot(data=countries) + geom_polygon(aes(x=long,y=lat,group = group),fill="white",color="grey",size=0.3) +
  coord_cartesian(xlim=c(-20,30),ylim=c(25,70)) +
  geom_point(data=plot.data[plot.data$SR_mixture == 2,],
             aes(x=long_mix,y=lat_mix),alpha=0.15,size=0.6) +
  xlab("Longitude") + ylab("Latitude") + ggtitle("2-species mixture plots")
graphics.off()

png(filename= "Appendix_map_3_mix2.png",height = 1500, width = 1500, res= 300)
ggplot(data=countries) + geom_polygon(aes(x=long,y=lat,group = group),fill="white",color="grey",size=0.3) +
  coord_cartesian(xlim=c(-20,30),ylim=c(25,70)) +
  geom_point(data=plot.data[plot.data$SR_mixture == 3,],
             aes(x=long_mix,y=lat_mix),alpha=0.15,size=0.6) +
  xlab("Longitude") + ylab("Latitude") + ggtitle("3-species mixture plots")
graphics.off()

png(filename= "Appendix_map_4_mix2.png",height = 1500, width = 1500, res= 300)
ggplot(data=countries) + geom_polygon(aes(x=long,y=lat,group = group),fill="white",color="grey",size=0.3) +
  coord_cartesian(xlim=c(-20,30),ylim=c(25,70)) +
  geom_point(data=plot.data[plot.data$SR_mixture == 4,],
             aes(x=long_mix,y=lat_mix),alpha=0.15,size=0.6) +
  xlab("Longitude") + ylab("Latitude") + ggtitle(">= 4 species mixture plots")
graphics.off()

# save stand conditions for later plotting with the other approaches
for.export.histograms = rbind(plot.data,plot.data)
for.export.histograms$SR = for.export.histograms$SR_mixture
for.export.histograms$SR[1:nrow(for.export.histograms)/2] = 1
for.export.histograms$plotcode = as.character(for.export.histograms$mixture)
for.export.histograms$plotcode[1:nrow(for.export.histograms)/2] = as.character(for.export.histograms$monoculture[1:nrow(for.export.histograms)/2])
for.export.histograms = for.export.histograms[,.(plotcode,SR)]
for.export.histograms = merge(for.export.histograms,
                              plot.statistics4[,.(plotcode,plot_sum_ba,plot_ba_heterogeneity,bio1,bio4,bio12,bio15)],
                              by.x="plotcode",by.y = "plotcode")

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts")

write.table(for.export.histograms,"data\\Inventories\\NFI with French data\\intermediate_output\\export_for_histogram_plotting.txt", quote=F,row.names=F, sep="\t")




```{r match monocultures and mixtures for over-yielding}

species.ba.n.growth.clean = data.table(read.table("intermediate_output\\2017\\species_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))
plot.statistics4 = data.table(read.table("intermediate_output\\2017\\plot_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))

# environment to determine distance = bio1, bio4, bio12, bio15, slope,sum_ba_plot, ba_heterogeneity
data = species.ba.n.growth.clean
data = merge(data, plot.statistics4[,.(plotcode,plot_sum_ba,plot_ba_heterogeneity)],
             by.x="plotcode",by.y="plotcode")
data = data[,.(plotcode,code_new,SR,europeanforesttype,species_growth_per_tree,species_sum_ba,species_n_trees,
               bio1 = scale(bio1), bio4=scale(bio4), bio12=scale(bio12),bio15 = scale(bio15),slope=scale(slope),
               plot_sum_ba = scale(plot_sum_ba),plot_ba_heterogeneity=scale(plot_ba_heterogeneity))]

get.best.monos.and.mixtures = function(data,plotcode_mix){
  df.mixture = data[plotcode %in% plotcode_mix]
  
  #check if monocultures of all species in this forest_type available
  all.monos = subset(data, europeanforesttype %in% df.mixture$europeanforesttype &
                       code_new %in%  df.mixture$code_new & SR ==1)
  
  # only continue if for each mixed species there is a monoculture at least
  if(length(unique(all.monos$code_new)) == nrow(df.mixture)){
    mixture.and.best.monos = df.mixture
    mixture.and.best.monos$distance = 0
    
    # calculate euclidean distance
    df.for.distance = rbind(df.mixture,all.monos)
    df.for.distance$distance = 0
    dist.matrix = as.matrix(daisy(df.for.distance[,.(bio1,bio4,bio12,bio15,slope,plot_ba_heterogeneity,plot_ba_heterogeneity)],metric="euclidean"))
    
    #extract monos with lowest distance
    for(mixture.species in df.mixture$code_new){
      
      monos.rownumbers = which(df.for.distance$code_new %in% mixture.species &
                                 df.for.distance$SR == 1)
      mix.rownumber = which(df.for.distance$code_new %in% mixture.species &
                              df.for.distance$SR > 1)
      if(length(monos.rownumbers) > 1){
        best.mono = which(dist.matrix[mix.rownumber,monos.rownumbers] ==
                            min(dist.matrix[mix.rownumber,monos.rownumbers]))[1]
        mixture.and.best.monos = rbind(mixture.and.best.monos,
                                       df.for.distance[as.integer(names(best.mono)),])
        mixture.and.best.monos$distance[nrow(mixture.and.best.monos)] =
          dist.matrix[mix.rownumber,as.integer(names(best.mono)[1])]
      }else{
        mixture.and.best.monos = rbind(mixture.and.best.monos,
                                       df.for.distance[which(df.for.distance$code_new %in% mixture.species &
                                                               df.for.distance$SR == 1),])
        mixture.and.best.monos$distance[nrow(mixture.and.best.monos)] =
          dist.matrix[mix.rownumber,monos.rownumbers]}
    }
    mixture.and.best.monos$plotcode_mix = mixture.and.best.monos$plotcode[1]
    mixture.and.best.monos$plotcodes_monos = paste(c(as.character(mixture.and.best.monos$plotcode)[which(mixture.and.best.monos$SR == 1)]),collapse="-")
    return(mixture.and.best.monos)
  }
}

calculate.overyield = function(mixture.and.best.monos){
  growth.mix.estimated.from.mono = 0
  growth.mix = 0
  species.temp = ""
  
  #here overyielding is calculated but how, if different numbers of trees, probably shit
  for(species.growth in unique(as.character(mixture.and.best.monos$code_new))){
    growth.mix.estimated.from.mono = growth.mix.estimated.from.mono +
      mixture.and.best.monos[code_new %in% species.growth & SR == 1]$species_growth_per_tree*
      mixture.and.best.monos[code_new %in% species.growth & SR > 1]$species_n_trees
    growth.mix = growth.mix +
      mixture.and.best.monos[code_new %in% species.growth & SR > 1]$species_growth_per_tree*
      mixture.and.best.monos[code_new %in% species.growth & SR > 1]$species_n_trees
  }
  
  #fill results table
  results.overyield.temp = 
    data.frame("species"= paste(unique(as.character(mixture.and.best.monos$code_new)),collapse="-"),
               "forest_type" = mixture.and.best.monos$europeanforesttype[1],
               "plotcode_mix" = mixture.and.best.monos$plotcode[which(mixture.and.best.monos$SR > 1)][1],
               "growth_mix_estimated" = growth.mix.estimated.from.mono,
               "growth_mix" = growth.mix,
               "trees_n_mono" = sum(mixture.and.best.monos$species_n_trees[
                 which(mixture.and.best.monos$SR == 1)]),
               "trees_n_mix" = sum(mixture.and.best.monos$species_n_trees[
                 which(mixture.and.best.monos$SR > 1)]),
               "mean_distance" = mean(mixture.and.best.monos$distance[
                 which(mixture.and.best.monos$SR == 1)]))
  
  return(results.overyield.temp)
}


# get all mixture plotcodes
mixture.plotcodes = unique(subset(plot.statistics4, SR > 1)$plotcode)

# build results tables
all.mixture.and.best.monos=data.frame(plotcode=NA,code_new=NA,SR=NA,europeanforesttype=NA, species_growth_per_tree=NA,species_sum_ba=NA,species_n_trees=NA,bio1=NA,bio4=NA,bio12=NA,bio15=NA,slope=NA,plot_sum_ba=NA,plot_ba_heterogeneity=NA,distance=NA,plotcode_mix=NA,plotcodes_monos=NA)[0,]
all.results.overyield = data.frame(species=NA,forest_type=NA,plotcode_mix=NA,growth_mix_estimated=NA,growth_mix=NA,trees_n_mono=NA,trees_mix=NA,mean_distance=NA)[0,]

# run algorithm
pb <- txtProgressBar(min=0,max=length(mixture.plotcodes),initial=0,style=3)
for(i in 1:length(mixture.plotcodes)){
  mixture.and.best.monos = get.best.monos.and.mixtures(data,mixture.plotcodes[i])
  if(is.null(mixture.and.best.monos) == FALSE){
    all.mixture.and.best.monos = rbind(all.mixture.and.best.monos,mixture.and.best.monos)
    all.results.overyield = rbind(all.results.overyield,calculate.overyield(mixture.and.best.monos))
  }
  setTxtProgressBar(pb, i)
}

# append SR, sum_ba and ba_heterogeneity bioclim
all.results.overyield.merged = merge(all.results.overyield,plot.statistics4[,.(plotcode,SR,plot_sum_ba,plot_ba_heterogeneity,bio1,bio4,bio12,bio15,pet)],by.x="plotcode_mix",by.y="plotcode")

write.table(all.results.overyield.merged,"intermediate_output\\2017\\overyield_results_for_Helge.txt",sep="\t",dec=".",row.names=F,quote=F)
write.table(all.mixture.and.best.monos,"intermediate_output\\2017\\Mono_Mix_matched_for_helge.txt",sep="\t",dec=".",row.names=F,quote=F)

```
```{r check all the output}
#setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Inventories")
setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Inventories")

results.overyield = data.table(read.table("intermediate_output\\2017\\overyield_results_for_Helge.txt",sep="\t",dec=".",header=T))
mixture.and.best.monos = data.table(read.table("intermediate_output\\2017\\Mono_Mix_matched_for_helge.txt",sep="\t",dec=".",header=T))

```



