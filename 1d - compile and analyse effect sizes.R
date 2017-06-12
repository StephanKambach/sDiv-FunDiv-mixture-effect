### FunDiv Synthesis:
### Diversity ~ Growth across platforms
### By: Stephan Kambach
### Date: 1.03.2016

rm(list=ls())
gc()

#libraries
if(!require(metafor)){install.packages("metafor")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(scales)){install.packages("scales")}
if(!require(lmerTest)){install.packages("lmerTest")}
if(!require(corrplot)){install.packages("corrplot")}
if(!require(extrafont)){install.packages("extrafont")}
if(!require(raster)){install.packages("raster")}
if(!require(data.table)){install.packages("data.table")}

library(metafor)
library(ggplot2)
library(scales)
library(extrafont)
library(raster)
library(data.table)
library(lmerTest)
font_import(pattern ="[T/t]imes")
y
loadfonts(device="win")
fonts()

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts")

experiments = data.table(read.table("data\\Experiments\\all_experiments_re-arranged.csv",sep="\t",header=T))
str(experiments)
exploratories = data.table(read.table("data\\Exploratories\\all_exploratories.csv",sep="\t",header=T))
str(exploratories)
inventories = data.table(read.table("data\\Inventories\\NFI with French data\\inventories_effect_sizes.txt",sep="\t",header=T))
str(inventories)

# rename and column names
inventories$platform  = "inventories"
names(inventories) = c("species","forest_type","SR_mixture","countries_compared","mono_mean","mono_sd","mono_n",
                       "mix_mean","mix_sd","mix_n","platform") 
experiments$counties_compared = paste(experiments$platform,experiments$exp, sep = "_")
exploratories$counties_compared = paste(exploratories$platform, exploratories$exp, sep = "_")

names(experiments) = c("species","forest_type","mono_mean","mono_sd","mono_n","SR_mixture",
                       "mix_mean","mix_sd","mix_n","platform","countries_compared")

names(exploratories) = c("species","forest_type","SR_mixture","mono_mean","mono_sd","mono_n",
                         "mix_mean","mix_sd","mix_n","platform","countries_compared")

# rename exploratory forest types
exploratories$forest_type = as.character(exploratories$forest_type)
exploratories$forest_type[which(exploratories$forest_type %in% "acidophilous")] = "Acidophilous"
exploratories$forest_type[which(exploratories$forest_type %in% "alpine")] = "Alpine"
exploratories$forest_type[which(exploratories$forest_type %in% "beech")] = "Beech"
exploratories$forest_type[which(exploratories$forest_type %in% "boreal")] = "Boreal"
exploratories$forest_type[which(exploratories$forest_type %in% "broadleaved_evergreen")] = "Broadleaved evergreen"
exploratories$forest_type[which(exploratories$forest_type %in% "floodplain")] = "Floodplain"
exploratories$forest_type[which(exploratories$forest_type %in% "hemiboreal")] = "Hemiboreal"
exploratories$forest_type[which(exploratories$forest_type %in% "introduced")] = "Introduced"
exploratories$forest_type[which(exploratories$forest_type %in% "mediterranean_coniferous")] = "Mediterranean coniferous"
exploratories$forest_type[which(exploratories$forest_type %in% "mesophytic_deciduous")] = "Mesophytic deciduous"
exploratories$forest_type[which(exploratories$forest_type %in% "mountain_beech")] = "Mountain beech"
exploratories$forest_type[which(exploratories$forest_type %in% "non_riverine")] = "Non-riverine"
exploratories$forest_type[which(exploratories$forest_type %in% "swamp")] = "Swamp"
exploratories$forest_type[which(exploratories$forest_type %in% "thermophilous_deciduous")] = "Thermophilous deciduous"
exploratories$forest_type = factor(exploratories$forest_type, levels = c(
    "Mediterranean coniferous","Thermophilous deciduous","Broadleaved evergreen","Beech","Introduced","Non-riverine",
    "Mesophytic deciduous","Floodplain","Swamp","Acidophilous","Mountain beech","Hemiboreal","Boreal","Alpine"))

# assign row weightening due to multiple comparisons
inventories.weights = inventories[,.(weight = 1 / length(countries_compared)),
                                  by = .(species, forest_type)]
inventories.merged = merge(inventories, inventories.weights, by.x = c("species","forest_type"),
                           by.y = c("species","forest_type"))

experiments.weights = experiments[,.(weight = 1 / length(countries_compared)),
                                  by = .(species, forest_type)]
experiments.merged = merge(experiments, experiments.weights, by.x = c("species","forest_type"),
                           by.y = c("species","forest_type"))

exploratories.weights = exploratories[,.(weight = 1 / length(countries_compared)),
                                  by = .(species, forest_type)]
exploratories.merged = merge(exploratories, exploratories.weights, by.x = c("species","forest_type"),
                           by.y = c("species","forest_type"))

# put data together
all.platforms = rbind(inventories.merged,experiments.merged,exploratories.merged)
str(all.platforms)
table(all.platforms$forest_type)

# calculate Log-RR and remove NAS
all.platforms$RR = log(all.platforms$mix_mean / all.platforms$mono_mean)
all.platforms.complete = all.platforms[!(is.na(all.platforms$RR)),]

# check
str(all.platforms.complete)

# has countries_compared an effect?
countries_compared.test = lm(RR ~ countries_compared, data = subset(all.platforms.complete, platform %in% "inventories"),
                             weights = weight)
anova(countries_compared.test)
summary(countries_compared.test)

# grand mean
lmer.grand.mean = lmer(RR ~ (1|countries_compared), data = all.platforms.complete,
                       weights = weight)

summary(lmer.grand.mean)
# approximated CI = 0.08665 +/- 1.96*0.05770
# CI = -0.026442 to 0.199742

# SR as predictor variable
lmer.SR = lmer(RR ~ SR_mixture + (1|countries_compared), data = all.platforms.complete,
                       weights = weight)

summary(lmer.SR)
anova(lmer.grand.mean,lmer.SR)

# mean effect size per platform
lmer.platform = lmer(RR ~ platform + (1|countries_compared) - 1, data = all.platforms.complete,
                     weights = weight)
summary(lmer.platform)

grand.mean.inventories = lmer(RR ~ (1|countries_compared), data = subset(all.platforms.complete, platform %in% "inventories"),
                        weights = weight)
summary(grand.mean.inventories)

grand.mean.experiments = lm(RR ~ 1, data = subset(all.platforms.complete, platform %in% "experiments"),
                            weights = weight)
summary(grand.mean.experiments)

grand.mean.exploratories = lm(RR ~ 1, data = subset(all.platforms.complete, platform %in% "exploratories"),
                        weights = weight)
summary(grand.mean.exploratories)

# mean effect size per forest_type
lmer.forest.type = lmer(RR ~ forest_type + (1|countries_compared) - 1, data = all.platforms.complete,
                        weights = weight)
summary(lmer.forest.type)

# mean effect size per exp per platform
lmer.inventories.forest.type = lmer(RR ~ forest_type + (1|countries_compared) - 1, 
                                    weights = weight, data = all.platforms.complete, subset= platform %in% "inventories")
summary(lmer.inventories.forest.type)
lmer.experiments.forest.type = lm(RR ~ forest_type - 1, 
                                    weights = weight, data = all.platforms.complete, subset= platform %in% "experiments")
summary(lmer.experiments.forest.type)
lmer.exploratories.forest.type = lm(RR ~ forest_type - 1, 
                                      weight = weight, data = all.platforms.complete, subset= platform %in% "exploratories")
summary(lmer.exploratories.forest.type)


# get the number of mono-mix-comparisons per platform-exp - combination
to.count.replicates = all.platforms.complete[,.(n_species = length(unique(species)),n_comparisons = sum(mix_n)),
                                          by= .(platform,forest_type)]
# all species
length(unique(all.platforms.complete$species[all.platforms.complete$platform %in% "inventories"]))
sum(all.platforms.complete$mix_n[all.platforms.complete$platform %in% "inventories"])

length(unique(all.platforms.complete$species[all.platforms.complete$platform %in% "experiments"]))
sum(all.platforms.complete$mix_n[all.platforms.complete$platform %in% "experiments"])

length(unique(all.platforms.complete$species[all.platforms.complete$platform %in% "exploratories"]))
sum(all.platforms.complete$mix_n[all.platforms.complete$platform %in% "exploratories"])

# plotting ##########################
plot.data1 = all.platforms.complete[,.(platform,species,forest_type,SR_mixture,RR)] 

plot.data.inventories = subset(plot.data1, platform %in% "inventories")
plot.data.inventories = plot.data.inventories[order(plot.data.inventories$forest_type,plot.data.inventories$RR,decreasing=F),]
plot.data.inventories$exp = factor(plot.data.inventories$forest_type)

plot.data.experiments = subset(plot.data1, platform %in% "experiments")
plot.data.experiments = plot.data.experiments[order(plot.data.experiments$forest_type,plot.data.experiments$RR,decreasing=F),]
plot.data.experiments$exp = factor(plot.data.experiments$forest_type)

plot.data.exploratories = subset(plot.data1, platform %in% "exploratories")
plot.data.exploratories = plot.data.exploratories[order(plot.data.exploratories$forest_type,plot.data.exploratories$RR,decreasing=F),]
plot.data.exploratories$exp = factor(plot.data.exploratories$forest_type)


#insert empty row between different levels of exp
insert.empty.rows = function(data.to.split,factor.to.split){
    empty.row = data.to.split[NA,][1,]
    factor.table = table(data.to.split[,which(names(data.to.split) %in% factor.to.split)])
    
    new.table = data.to.split[0,]
    
    start.row = 1
    end.row = 0
    
    for(i in factor.table){
        
        end.row = end.row + i
        new.table= rbind(new.table,
                         data.to.split[start.row : end.row,],
                         empty.row)
        start.row=start.row + i
    }
    return(new.table)
}

plot.data.experiments = insert.empty.rows(plot.data.experiments,"forest_type")
plot.data.exploratories = insert.empty.rows(plot.data.exploratories,"forest_type")
plot.data.inventories = insert.empty.rows(plot.data.inventories,"forest_type")

plot.data.experiments$uniqueID = seq(1:nrow(plot.data.experiments))
plot.data.exploratories$uniqueID = seq(1:nrow(plot.data.exploratories))
plot.data.inventories$uniqueID = seq(1:nrow(plot.data.inventories))

inventories.estimates = coefficients(summary(lmer.inventories.forest.type))[,"Estimate"]
experiments.estimates = coefficients(summary(lmer.experiments.forest.type))[,"Estimate"]
exploratories.estimates = coefficients(summary(lmer.exploratories.forest.type))[,"Estimate"]

inventories.sds = coefficients(summary(lmer.inventories.forest.type))[,"Std. Error"]
experiments.sds = coefficients(summary(lmer.experiments.forest.type))[,"Std. Error"]
exploratories.sds = coefficients(summary(lmer.exploratories.forest.type))[,"Std. Error"]

inventories.replicates = data.frame(table(as.character(subset(all.platforms, platform %in% "inventories")$forest_type)))
inventories.replicates = rbind(data.frame("Var1" = "grand mean",
                                          "Freq" = nrow(subset(all.platforms, platform %in% "inventories"))),
                               inventories.replicates)
experiments.replicates = data.frame(table(as.character(subset(all.platforms, platform %in% "experiments")$forest_type)))
experiments.replicates = rbind(data.frame("Var1" = "grand mean",
                                          "Freq" = nrow(subset(all.platforms, platform %in% "experiments"))),
                               experiments.replicates)
exploratories.replicates = data.frame(table(as.character(subset(all.platforms, platform %in% "exploratories")$forest_type)))
exploratories.replicates = rbind(data.frame("Var1" = "grand mean",
                                          "Freq" = nrow(subset(all.platforms, platform %in% "exploratories"))),
                               exploratories.replicates)

# create dataframe for plotting grand means
plot.inventories.grand.means = data.frame(
    "forest_type"= c("grand mean",NA,names(inventories.estimates)),
    "RR" = c(coefficients(summary(grand.mean.inventories))[1,"Estimate"],NA,inventories.estimates),
    "RR_SD" = c(coefficients(summary(grand.mean.inventories))[1,"Std. Error"],NA,inventories.sds))
plot.experiments.grand.means = data.frame(
    "forest_type"= c("grand mean",NA,names(experiments.estimates)),
    "RR" = c(coefficients(summary(grand.mean.experiments))[1,"Estimate"],NA,experiments.estimates),
    "RR_SD" = c(coefficients(summary(grand.mean.experiments))[1,"Std. Error"],NA,experiments.sds))
plot.exploratories.grand.means = data.frame(
    "forest_type"= c("grand mean",NA,names(exploratories.estimates)),
    "RR" = c(coefficients(summary(grand.mean.exploratories))[1,"Estimate"],NA,exploratories.estimates),
    "RR_SD" = c(coefficients(summary(grand.mean.exploratories))[1,"Std. Error"],NA,exploratories.sds))


plot.inventories.grand.means$forest_type = gsub("forest_type","",plot.inventories.grand.means$forest_type)
plot.experiments.grand.means$forest_type = gsub("forest_type","",plot.experiments.grand.means$forest_type)
plot.exploratories.grand.means$forest_type = gsub("forest_type","",plot.exploratories.grand.means$forest_type)

plot.inventories.grand.means$forest_type = factor(plot.inventories.grand.means$forest_type, levels = c(plot.inventories.grand.means$forest_type))
plot.inventories.grand.means = merge(plot.inventories.grand.means,inventories.replicates,
                                     by.x = "forest_type", by.y = "Var1", all.x=T)
plot.experiments.grand.means$forest_type = factor(plot.experiments.grand.means$forest_type, levels = c(plot.experiments.grand.means$forest_type))
plot.experiments.grand.means = merge(plot.experiments.grand.means,experiments.replicates,
                                     by.x = "forest_type", by.y = "Var1", all.x=T)
plot.exploratories.grand.means$forest_type = factor(plot.exploratories.grand.means$forest_type, levels = c(plot.exploratories.grand.means$forest_type))
plot.exploratories.grand.means = merge(plot.exploratories.grand.means,exploratories.replicates,
                                     by.x = "forest_type", by.y = "Var1", all.x=T)

plot.inventories.grand.means$forest_type = as.character(plot.inventories.grand.means$forest_type)
plot.experiments.grand.means$forest_type = as.character(plot.experiments.grand.means$forest_type)
plot.exploratories.grand.means$forest_type = as.character(plot.exploratories.grand.means$forest_type)

# sort order
#sort.order  = c("grand mean",NA,"thermophilous_deciduous","mediterranean_coniferous","introduced","broadleaved_evergreen",
#               "beech","acidophilous","mesophytic_deciduous","non_riverine","floodplain","swamp",
#               "mountain_beech","boreal","hemiboreal","alpine")
plot.inventories.grand.means = plot.inventories.grand.means[c(1,nrow(plot.inventories.grand.means),
                                                              order(plot.inventories.grand.means$RR)),]
plot.inventories.grand.means = plot.inventories.grand.means[!duplicated(plot.inventories.grand.means),]
plot.experiments.grand.means = plot.experiments.grand.means[c(1,nrow(plot.experiments.grand.means),
                                                              order(plot.experiments.grand.means$RR)),]
plot.experiments.grand.means = plot.experiments.grand.means[!duplicated(plot.experiments.grand.means),]
plot.exploratories.grand.means = plot.exploratories.grand.means[c(1,nrow(plot.exploratories.grand.means),
                                                              order(plot.exploratories.grand.means$RR)),]
plot.exploratories.grand.means = plot.exploratories.grand.means[!duplicated(plot.exploratories.grand.means),]

plot.inventories.grand.means$uniqueID = seq(1:nrow(plot.inventories.grand.means))
plot.experiments.grand.means$uniqueID = seq(1:nrow(plot.experiments.grand.means))
plot.exploratories.grand.means$uniqueID = seq(1:nrow(plot.exploratories.grand.means))

# factor levels
plot.inventories.grand.means$forest_type = factor(plot.inventories.grand.means$forest_type, levels = plot.inventories.grand.means$forest_type)
plot.experiments.grand.means$forest_type = factor(plot.experiments.grand.means$forest_type, levels = plot.experiments.grand.means$forest_type)
plot.exploratories.grand.means$forest_type = factor(plot.exploratories.grand.means$forest_type, levels = plot.exploratories.grand.means$forest_type)

#############
# inventories
svg(filename="manuscript 1_meta-analysis\\output\\inventories_grand_means.svg",height=5,width=5)  
ggplot(data=plot.inventories.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_SD), ymax = RR + (1.96* RR_SD),color=forest_type),size = 1) +
    geom_point(aes(x=factor(uniqueID),y=RR, fill = forest_type), shape = 21, color = "black",size = 3) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() +

    theme_bw() + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left", text = element_text(family = "Times New Roman"),
          panel.grid = element_blank()) +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.inventories.grand.means$forest_type))),
                       labels = c("Non-riverine","Alpine","Acidophilous","Broadleaved-evergreen","Hemiboreal",
                                  "Boreal","Introduced","Thermophilous-deciduous","Mediterranean-coniferous",
                                  "Beech","Mountain beech","Mesophytic-deciduous","Swamp","Floodplain","Grand mean"),
                       values= c("red",rep("#efc406",14))) +
    scale_fill_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.inventories.grand.means$forest_type))),
                       labels = c("Non-riverine","Alpine","Acidophilous","Broadleaved-evergreen","Hemiboreal",
                                  "Boreal","Introduced","Thermophilous-deciduous","Mediterranean-coniferous",
                                  "Beech","Mountain beech","Mesophytic-deciduous","Swamp","Floodplain","Grand mean"),
                      values= c("red",rep("#efc406",14))) +
    
    ggtitle("a)") + ylab("Mixture effect (log RR) - inventory approach")
graphics.off()


# experiments
svg(filename="manuscript 1_meta-analysis\\output\\experiments_grand_means.svg",height=2.5,width=5)  
ggplot(data=plot.experiments.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_SD), ymax = RR + (1.96* RR_SD),color=forest_type),size = 1) +
    geom_point(aes(x=factor(uniqueID),y=RR, fill = forest_type), shape = 21, color = "black",size = 3) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() +
    
    theme_bw() + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),legend.position="right", text = element_text(family = "Times New Roman"),
                       panel.grid = element_blank()) +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.experiments.grand.means$forest_type))),
                       labels = c("BIOTREE","Kreinitz","FORBIO-Zedelgem","Satakunta","FORBIO-Gedinne",
                                  "ORPHEE","Grand mean"),
                       values= c("red",rep("#82B5CE",7))) +
    scale_fill_manual(name ="Forest type",
                      breaks = c(rev(levels(plot.experiments.grand.means$forest_type))),
                      labels = c("BIOTREE","Kreinitz","FORBIO-Zedelgem","Satakunta","FORBIO-Gedinne",
                                 "ORPHEE","Grand mean"),
                      values= c("red",rep("#82B5CE",7))) +
    
    ggtitle("b)") + ylab("Mixture effect (log RR) - experimental approach")
graphics.off()

# exploratories
svg(filename="manuscript 1_meta-analysis\\output\\exploratories_grand_means.svg",height=2.5,width=5)  
ggplot(data=plot.exploratories.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_SD), ymax = RR + (1.96* RR_SD),color=forest_type),size = 1) +
    geom_point(aes(x=factor(uniqueID),y=RR, fill = forest_type), shape = 21, color = "black",size = 3) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() +
    
    theme_bw() + theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
                       axis.title.y = element_blank(),legend.position="right", text = element_text(family = "Times New Roman"),
                       panel.grid = element_blank()) +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.exploratories.grand.means$forest_type))),
                       labels = c("Mediterranean-coniferous","Thermophilous-deciduous","Hemiboreal","Beech",
                                  "Mountain beech","Boreal","Grand mean" ),
                       values= c("red",rep("#96B94D",7))) +
    scale_fill_manual(name ="Forest type",
                      breaks = c(rev(levels(plot.exploratories.grand.means$forest_type))),
                      labels = c("Mediterranean-coniferous","Thermophilous-deciduous","Hemiboreal","Beech",
                                 "Mountain beech","Boreal","Grand mean" ),
                      values= c("red",rep("#96B94D",7))) +
    
    ggtitle("c)") + ylab("Mixture effect (log RR) - exploratory approach")
graphics.off()



###### get the stand conditions to plot histograms alltogether
stand.cond.inv = read.table("data\\Inventories\\NFI with French data\\intermediate_output\\export_for_histogram_plotting.txt",
                            header=T, sep="\t")
locations.exper.explo = read.csv("intermediate output\\locations of experiments and exploratories.csv",sep=";",dec=",")
str(locations.exper.explo)
clim = getData("worldclim",var="bio",res=2.5 )
bioclim.data = extract(clim,locations.exper.explo[,c("Longitude","Latitude")])
locations.exper.explo = cbind(locations.exper.explo,bioclim.data)

locations.exper.explo$bio1 = locations.exper.explo$bio1 / 10
locations.exper.explo$bio4 = locations.exper.explo$bio4 / 10
stand.cond.inv$bio1 = stand.cond.inv$bio1 / 10
stand.cond.inv$bio4 = stand.cond.inv$bio4 / 10

# merge data
stand.cond.inv$platform = "inventories"
locations.exper.explo = locations.exper.explo[,c("Forest_Type","Type","bio1","bio4","bio12","bio15")]
names(locations.exper.explo) = c("plotcode","platform","bio1","bio4","bio12","bio15" )
locations.exper.explo$plot_sum_ba = NA
locations.exper.explo$plot_ba_heterogeneity = NA
locations.exper.explo$SR = NA

stand.cond.all = rbind(stand.cond.inv,locations.exper.explo)

svg(filename = "manuscript 1_meta-analysis\\output\\Appendix_hist_bio1.svg", height = 4, width = 12)
ggplot() + geom_histogram(data=subset(stand.cond.all, platform %in% "inventories"), aes(x = bio1), bins = 100, color= NA, fill = "#DFB738", alpha = 0.5) +
    geom_point(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio1, y = rep(1,6)), fill = "#82B5CE", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio1, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    geom_point(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio1, y = rep(1,6)), size =4, fill = "#96B94D", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio1, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    ggtitle("Mean annual temperature") +
    theme_bw() + theme(panel.grid = element_blank(), text = element_text(family = "Times New Roman"))
graphics.off()

svg(filename = "manuscript 1_meta-analysis\\output\\Appendix_hist_bio4.svg", height = 4, width = 12)
ggplot() + geom_histogram(data=subset(stand.cond.all, platform %in% "inventories"), aes(x = bio4), bins = 100, color= NA, fill = "#DFB738", alpha = 0.5) +
    geom_point(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio4, y = rep(1,6)), fill = "#82B5CE", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio4, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    geom_point(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio4, y = rep(1,6)), size =4, fill = "#96B94D", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio4, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    ggtitle("Temperature Seasonality") +
    theme_bw() + theme(panel.grid = element_blank(), text = element_text(family = "Times New Roman"))
graphics.off()

svg(filename = "manuscript 1_meta-analysis\\output\\Appendix_hist_bio12.svg", height = 4, width = 12)
ggplot() + geom_histogram(data=subset(stand.cond.all, platform %in% "inventories"), aes(x = bio12), bins = 100, color= NA, fill = "#DFB738", alpha = 0.5) +
    geom_point(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio12, y = rep(1,6)), fill = "#82B5CE", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio12, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    geom_point(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio12, y = rep(1,6)), size =4, fill = "#96B94D", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio12, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    ggtitle("Annual precipitation") +
    theme_bw() + theme(panel.grid = element_blank(), text = element_text(family = "Times New Roman"))
graphics.off()

svg(filename = "manuscript 1_meta-analysis\\output\\Appendix_hist_bio15.svg", height = 4, width = 12)
ggplot() + geom_histogram(data=subset(stand.cond.all, platform %in% "inventories"), aes(x = bio15), bins = 100, color= NA, fill = "#DFB738", alpha = 0.5) +
    geom_point(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio15, y = rep(1,6)), fill = "#82B5CE", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Experiment"), aes(x = bio15, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman") +
    geom_point(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio15, y = rep(1,6)), size =4, fill = "#96B94D", pch =21, color ="black") +
    geom_text(data = subset(stand.cond.all, platform %in% "Forest"), aes(x = bio15, y= rep(1,6), label = plotcode), angle = 90, hjust = -0.1, family = "Times New Roman")+
    ggtitle("Precipitation seasonality") +
    theme_bw() + theme(panel.grid = element_blank(), text = element_text(family = "Times New Roman"))
graphics.off()



##################################################################
# correlate model coefficient estimates between research approaches

# calculate coefficient estimates
lmer.species.inv = lmer(RR ~ species + (1|countries_compared) -1 , data = all.platforms.complete,
                       weights = weight, subset = platform %in% "inventories")
lmer.species.exper = lmer(RR ~ species + (1|countries_compared) - 1, data = all.platforms.complete,
                        weights = weight, subset = platform %in% "experiments")
lmer.species.explo = lmer(RR ~ species + (1|countries_compared) - 1, data = all.platforms.complete,
                        weights = weight, subset = platform %in% "exploratories")

# extract coefficient estimates
coef.inv = as.data.frame(summary(lmer.species.inv)$coefficients)
coef.exper = as.data.frame(summary(lmer.species.exper)$coefficients)
coef.explo = as.data.frame(summary(lmer.species.explo)$coefficients)

names(coef.inv)[1] = "Estimate_inv"
names(coef.exper)[1] = "Estimate_exper"
names(coef.explo)[1] = "Estimate_explo"

coef.inv$species = gsub(patter = "species",replacement="",rownames(coef.inv))
coef.exper$species = gsub(patter = "species",replacement="",rownames(coef.exper))
coef.explo$species = gsub(patter = "species",replacement="",rownames(coef.explo))

# create dataframe to contrast shared coefficient estimates
coef.inv.exper = merge(coef.inv,coef.exper, by.x = "species",by.y="species")
coef.inv.explo = merge(coef.inv,coef.explo, by.x = "species",by.y="species")
coef.exper.explo = merge(coef.exper,coef.explo, by.x = "species",by.y="species")

# correlations
cor.test(coef.inv.exper[,"Estimate_inv"],coef.inv.exper[,"Estimate_exper"], method="kendall")
cor.test(coef.inv.explo[,"Estimate_inv"],coef.inv.explo[,"Estimate_explo"], method="kendall")
cor.test(coef.exper.explo[,"Estimate_exper"],coef.exper.explo[,"Estimate_explo"], method="kendall")

# plotting
plot_inv_exper = ggplot(data=coef.inv.exper)
svg(filename = "manuscript 1_meta-analysis\\output\\comparison_inv_exper.svg",height=5, width = 5)
plot_inv_exper + geom_point(aes(x = Estimate_inv, y = Estimate_exper), size = 2) +
    geom_text(aes(x = Estimate_inv, y = Estimate_exper-0.02, family = "Times New Roman",
                  label = species), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"),
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Mixture effect (log RR) - inventory approach") +
    ylab("Mixture effect (log RR) - experimental approach") + ggtitle("a)")
graphics.off()

plot_inv_explo = ggplot(data=coef.inv.explo)
svg(filename = "manuscript 1_meta-analysis\\output\\comparison_inv_explo.svg",height=5, width = 5)
plot_inv_explo + geom_point(aes(x = Estimate_inv, y = Estimate_explo), size = 2) +
    geom_text(aes(x = Estimate_inv, y = Estimate_explo-0.02, family = "Times New Roman",
                  label = species), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"),
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Mixture effect (log RR) - inventory approach") +
    ylab("Mixture effect (log RR) - exploratory approach") + ggtitle("a)")
graphics.off()

plot_exper_explo = ggplot(data=coef.exper.explo)
svg(filename = "manuscript 1_meta-analysis\\output\\comparison_exper_explo.svg",height=5, width = 5)
plot_exper_explo + geom_point(aes(x = Estimate_explo, y = Estimate_exper), size = 2) +
    geom_text(aes(x = Estimate_explo, y = Estimate_exper-0.02, family = "Times New Roman",
                  label = species), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"),
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Mixture effect (log RR) - exploratory approach") +
    ylab("Mixture effect (log RR) - experimental approach") + ggtitle("a)")
graphics.off()


#######################################################################################
# correlate model coefficient estimates between inv and explo / separated by fores type

# get mean effect sizes per species and forest type
inv.species.forest.types = unique(all.platforms.complete[platform %in% "inventories",.(species,forest_type)])
explo.species.forest.types = unique(all.platforms.complete[platform %in% "exploratories",.(species,forest_type)])
inv.species.forest.types$lmer_intercept_inv = NA
explo.species.forest.types$lmer_intercept_explo = NA

# run models for all subsets of species and forest types
# inventories
for(i in 1:nrow(inv.species.forest.types)){
    data.subset = subset(all.platforms.complete,platform %in% "inventories" & 
                             species %in% as.character(inv.species.forest.types$species)[i] &
                             forest_type %in% as.character(inv.species.forest.types$forest_type)[i])
    if(length(unique(data.subset$countries_compared)) >=2){
        lmer.temp = lmer(RR ~ (1|countries_compared), data = data.subset, weights = weight)
        inv.species.forest.types$lmer_intercept_inv[i] = summary(lmer.temp)$coefficients[,"Estimate"]
    }else{
        lmer.temp = lm(RR ~ 1, data = data.subset, weights = weight)
        inv.species.forest.types$lmer_intercept_inv[i] = summary(lmer.temp)$coefficients[,"Estimate"]}}

# exploratories
for(i in 1:nrow(explo.species.forest.types)){
    data.subset = subset(all.platforms.complete,platform %in% "exploratories" & 
                             species %in% as.character(explo.species.forest.types$species)[i] &
                             forest_type %in% as.character(explo.species.forest.types$forest_type)[i])
    if(length(unique(data.subset$countries_compared)) >=2){
        lmer.temp = lmer(RR ~ (1|countries_compared), data = data.subset, weights = weight)
        explo.species.forest.types$lmer_intercept_explo[i] = summary(lmer.temp)$coefficients[,"Estimate"]
    }else{
        lmer.temp = lm(RR ~ 1, data = data.subset, weights = weight)
        explo.species.forest.types$lmer_intercept_explo[i] = summary(lmer.temp)$coefficients[,"Estimate"]}}

# test correlation in coefficients
inv.explo.species.forest.types = merge(inv.species.forest.types, explo.species.forest.types,
                                       by.x = c("species","forest_type"), by.y = c("species","forest_type"))
cor.test(inv.explo.species.forest.types$lmer_intercept_inv,inv.explo.species.forest.types$lmer_intercept_explo,
         method="kendall")

plot_inv_explo_forest_type = ggplot(data=inv.explo.species.forest.types)
svg(filename = "manuscript 1_meta-analysis\\output\\comparison_inv_explo_per_forest_type.svg",height=5, width = 5)
plot_inv_explo_forest_type + geom_point(aes(x = lmer_intercept_inv, y = lmer_intercept_explo, color = forest_type), size = 2) +
    geom_text(aes(x = lmer_intercept_inv, y = lmer_intercept_explo-0.02, color = forest_type, label = species),
              family = "Times New Roman", size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"), legend.position = "bottom",
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    scale_color_manual(name="Forest type",
                       values = c("grey","black","#2c7bb6","#fdae61","#abd9e9","#d7191c"),
                       breaks = c("Thermophilous deciduous","Mediterranean coniferous","Beech",
                                  "Mountain beech","Hemiboreal","Boreal"),
                       labels = c("Thermophilous deciduous","Mediterranean coniferous","Beech",
                                  "Mountain beech","Hemiboreal","Boreal")) +
    xlab("Mixture effect (log RR) - inventory approach") +
    ylab("Mixture effect (log RR) - exploratory approach") + ggtitle("d)")
graphics.off()














lmer.species.inv.forest.type = lmer(RR ~ species*forest_type + (1|countries_compared) -1 , data = all.platforms.complete,
                        weights = weight, subset = platform %in% "inventories")
lmer.species.explo.forest.type = lmer(RR ~ species*forest_type + (1|countries_compared) -1 , data = all.platforms.complete,
                                    weights = weight, subset = platform %in% "exploratories")

# extract coefficient estimates
coef.inv.forest.type = as.data.frame(summary(lmer.species.inv.forest.type)$coefficients)
coef.explo.forest.type = as.data.frame(summary(lmer.species.inv.forest.type)$coefficients)



svg(filename="figures\\experiments_exp_as_mod_RR.svg",height=10,width=5)  
ggplot(data=plot.data.experiments) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * sd), ymax = RR + 1.96* (RR_Var),colour=exp)) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp)) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2))) +
    
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name="Tree diversity experiment",
                       breaks = c(rev(levels(plot.data.experiments$exp))),
                       values= c("#5971C9","#8ACDFF","#5F98D8","#85E0E7","#007BB5")) +
    ggtitle("b) Experimental approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


# grand means of all platforms
svg(filename="figures\\experiments_exp_as_mod_grand_means_RR.svg",height=4,width=4)  
ggplot(data=plot.experiments.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp),size=6) +
    annotate("text",x = factor(plot.experiments.grand.means$uniqueID), y = plot.experiments.grand.means$RR +1.96* (plot.experiments.grand.means$RR_Var) + 0.035, label = plot.experiments.grand.means$replicates,size =3.5) + 
    
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2)))+
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Tree diversity experiment",
                       breaks = c(rev(row.names(rma.experiments.exp$b)),"grand mean"),
                       values= c("#5971C9","#8ACDFF","#5F98D8","#85E0E7","#007BB5","#AA0000")) +
    ggtitle("b) Experimental approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()

#############
# exploratories

svg(filename="figures\\exploratories_exp_as_mod_RR.svg",height=10,width=5)  
ggplot(data=plot.data.exploratories) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp)) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp)) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2))) +
    
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.data.exploratories$exp))),
                       values= c("#10580A","#21b014","#98A95D","#71D45A","#186F22","#91FB5C")) +
    ggtitle("c) Exploratory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


# grand means of all platforms
svg(filename="figures\\exploratories_exp_as_mod_grand_means_RR.svg",height=4,width=4)  
ggplot(data=plot.exploratories.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp),size=6) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    annotate("text",x = factor(plot.exploratories.grand.means$uniqueID), y = plot.exploratories.grand.means$RR +1.96* (plot.exploratories.grand.means$RR_Var) + 0.05, label = plot.exploratories.grand.means$replicates,size =3.5) + 
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2)))+
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name = "Forest type",
                       breaks = c(rev(row.names(rma.exploratories.exp$b)),"grand mean"),
                       values= c("#10580A","#21b014","#98A95D","#71D45A","#186F22","#91FB5C","#AA0000")) +
    ggtitle("c) Exploratory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


#############
# inventories

svg(filename="figures\\inventories_exp_as_mod_RR.svg",height=10,width=5)  
ggplot(data=plot.data.inventories) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp)) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp)) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2))) +
    
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.data.inventories$exp))),
                       values= c(rep(c("#634C00", "#957201","#C79801","#F9BE02","#FAD14D","#E4EAC7"),4),"#634C00", "#957201")) +
    ggtitle("a) Inventory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


# grand means of all platforms
svg(filename="figures\\inventories_exp_as_mod_grand_means_RR.svg",height=9,width=6)  
ggplot(data=plot.inventories.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp),size=6) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    annotate("text",x = factor(plot.inventories.grand.means$uniqueID), y = plot.inventories.grand.means$RR +1.96* (plot.inventories.grand.means$RR_Var) + 0.2, label = plot.inventories.grand.means$replicates,size =3.5) + 
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2)))+
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(row.names(rma.inventories.exp$b)),"grand mean"),
                       values= c(rep(c("#634C00", "#957201","#C79801","#F9BE02","#FAD14D","#E4EAC7"),4),"#634C00", "#957201","#AA0000")) +
    ggtitle("a) Inventory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()



# platform*forest type
lmer.inventories.forest.type = lmer(RR ~ exp + (1|species) - 1, data = subset(all.platforms.complete, platform %in% "inventories"))
lmer.experiments.forest.type = lmer(RR ~ exp + (1|species) - 1, data = subset(all.platforms.complete, platform %in% "experiments"))
lmer.exploratories.forest.type = lmer(RR ~ exp + (1|species) - 1, data = subset(all.platforms.complete, platform %in% "exploratories"))
summary(lmer.inventories.forest.type)
summary(lmer.experiments.forest.type)
summary(lmer.exploratories.forest.type)

# species
lmer.species = lmer(RR ~ species + (1|species) - 1, data = all.platforms.complete)
summary(lmer.species)

# platform*species
lmer.inventories.species = lmer(RR ~ species + (1|exp) - 1, data = subset(all.platforms.complete, platform %in% "inventories"))
lmer.experiments.species = lmer(RR ~ species + (1|exp) - 1, data = subset(all.platforms.complete, platform %in% "experiments"))
lmer.exploratories.species = lmer(RR ~ species + (1|exp) - 1, data = subset(all.platforms.complete, platform %in% "exploratories"))
summary(lmer.inventories.species)
summary(lmer.experiments.species)
summary(lmer.exploratories.species)

###################
# rank correlations

# extract estimates
coef.inventories.forest.type = data.frame("forest_type" = names(coef(summary(lmer.inventories.forest.type))[,"Estimate"]),
                                          "estimate_inventories" = coef(summary(lmer.inventories.forest.type))[,"Estimate"],
                                          "error_inventories" = coef(summary(lmer.inventories.forest.type))[,"Std. Error"])

coef.experiments.forest.type = data.frame("forest_type" = names(coef(summary(lmer.experiments.forest.type))[,"Estimate"]),
                                          "estimate_experiments" = coef(summary(lmer.experiments.forest.type))[,"Estimate"],
                                          "error_experiments" = coef(summary(lmer.experiments.forest.type))[,"Std. Error"])

coef.exploratories.forest.type = data.frame("forest_type" = names(coef(summary(lmer.inventories.forest.type))[,"Estimate"]),
                                            "estimate_exploratories" = coef(summary(lmer.inventories.forest.type))[,"Estimate"],
                                            "error_exploratories" = coef(summary(lmer.inventories.forest.type))[,"Std. Error"])

coef.inventories.species = data.frame("species" = names(coef(summary(lmer.inventories.species))[,"Estimate"]),
                                      "estimate_inventories" = summary(lmer.inventories.species)$coefficients[,"Estimate"],
                                      "error_inventories" = summary(lmer.inventories.species)$coefficients[,"Std. Error"])

coef.experiments.species = data.frame("species" = names(coef(summary(lmer.experiments.species))[,"Estimate"]),
                                      "estimate_experiments" = summary(lmer.experiments.species)$coefficients[,"Estimate"],
                                      "error_experiments" = summary(lmer.experiments.species)$coefficients[,"Std. Error"])

coef.exploratories.species = data.frame("species" = names(coef(summary(lmer.exploratories.species))[,"Estimate"]),
                                      "estimate_exploratories" = summary(lmer.exploratories.species)$coefficients[,"Estimate"],
                                      "error_exploratories" = summary(lmer.exploratories.species)$coefficients[,"Std. Error"])
# combine across platforms
forest.type.estimates = merge(coef.inventories.forest.type,coef.experiments.forest.type,
                              by.x = "forest_type",by.y= "forest_type",all.x=T,all.y=T)
forest.type.estimates = merge(forest.type.estimates,coef.exploratories.forest.type,
                              by.x = "forest_type",by.y= "forest_type",all.x=T,all.y=T)


species.estimates = merge(coef.inventories.species,coef.experiments.species,
                          by.x = "species",by.y= "species",all.x=T,all.y=T)
species.estimates = merge(species.estimates,coef.exploratories.species,
                          by.x = "species",by.y= "species",all.x=T,all.y=T)

# scale - not neccessary
#forest.type.estimates.scaled = cbind(forest.type.estimates$forest_type,scale(forest.type.estimates[,2:4]))
#species.estimates.scaled = cbind(species.estimates$species, scale(species.estimates[,2:4]))

# correlations
cor.test(forest.type.estimates$estimate_inventories, forest.type.estimates$estimate_experiments,method = "spearman")
cor.test(forest.type.estimates$estimate_inventories, forest.type.estimates$estimate_exploratories,method = "spearman")
cor.test(forest.type.estimates$estimate_experiments, forest.type.estimates$estimate_exploratories,method = "spearman")

cor.test(species.estimates$estimate_inventories, species.estimates$estimate_experiments, method = "spearman",use= )
cor.test(species.estimates$estimate_inventories, species.estimates$estimate_exploratories, method = "spearman")
cor.test(species.estimates$estimate_experiments, species.estimates$estimate_exploratories, method = "spearman")


# figures - species
plot.data = species.estimates
plot.data$species = gsub("species","",plot.data$species)

svg("manuscript 1_meta-analysis\\output\\Scatterplot_experiments-inventories.svg",height=5,width=6)
ggplot(data = plot.data) +
    geom_point(aes(x = estimate_inventories, y = estimate_experiments)) +
    geom_text(aes(x = estimate_inventories, y = estimate_experiments - 0.03, label = species)) +
    coord_cartesian(xlim = c(-0.25,0.9),ylim = c(-0.5,0.9)) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Inventory approach - mixture effect (log RR)") +
    ylab("Experimental approach - mixture effect (log RR)")
graphics.off()

svg("manuscript 1_meta-analysis\\output\\Scatterplot_exploratories-inventories.svg",height=5,width=6)
ggplot(data = plot.data) +
    geom_point(aes(x = estimate_inventories, y = estimate_exploratories)) +
    geom_text(aes(x = estimate_inventories, y = estimate_exploratories-0.03, label = species)) +
    coord_cartesian(xlim = c(-0.25,0.9),ylim = c(-0.5,0.9)) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Inventory approach - mixture effect (log RR)") +
    ylab("Exploratory approach - mixture effect (log RR)") 
graphics.off()

svg("manuscript 1_meta-analysis\\output\\Scatterplot_experiments-exploratories.svg",height=5,width=6)
ggplot(data = plot.data) +
    geom_point(aes(x = estimate_exploratories, y = estimate_experiments)) +
    geom_text(aes(x = estimate_exploratories, y = estimate_experiments-0.03, label = species)) +
    coord_cartesian(xlim = c(-0.25,0.9),ylim = c(-0.5,0.9)) +
    geom_hline(yintercept=0, linetype="dotted") +
    geom_vline(xintercept=0, linetype="dotted") +
    theme_bw() + theme(panel.grid = element_blank()) +
    xlab("Experimental approach - mixture effect (log RR)") +
    ylab("Exploratory approach - mixture effect (log RR)") 
graphics.off()

# correlation
cor.species = cor(cbind(species.estimates$estimate_inventories,species.estimates$estimate_experiments,species.estimates$estimate_exploratories),
                  use = "pairwise.complete.obs", method="spearman")
rownames(cor.species) = c("Inventories","Experiments","Exploratories")
colnames(cor.species) = c("Inventories","Experiments","Exploratories")

svg("manuscript 1_meta-analysis\\output\\corrplot_species.svg",height=6,width=6)
corrplot(cor.species,type= "lower",method="color",addCoef.col="black",tl.col="black")
graphics.off()

# plot in scatterplot - species
plot(estimate_inventories ~ estimate_experiments, data = species.estimates,xlim = c(-0.5,0.5),ylim=c(-0.5,0.5))
text(estimate_inventories ~ estimate_experiments, data = species.estimates, pos = 1,
     labels = gsub("species","",species.estimates$species),xlim = c(-0.5,0.5),ylim=c(-0.5,0.5))

plot(estimate_inventories ~ estimate_exploratories, data = species.estimates,xlim = c(-0.5,0.5),ylim=c(-0.5,0.5))
text(estimate_inventories ~ estimate_exploratories, data = species.estimates, pos = 1,
     labels = gsub("species","",species.estimates$species), xlim = c(-0.5,0.5),ylim=c(-0.5,0.5))


############################################################################
# MISSING = Comparison of mixtures that are in all three research approaches





#calculate effect sizes
ROM_effect_sizes = escalc(m1i=mix_mean,sd1i=mix_sd,n1i=mix_n,
              m2i=mono_mean,sd2i=mono_sd,n2i=mono_n,
              data=all.platforms,measure="ROM",append=F)
names(ROM_effect_sizes) = c("RR","RR_Var")

SMDH_effect_sizes = escalc(m1i=mix_mean,sd1i=mix_sd,n1i=mix_n,
                          m2i=mono_mean,sd2i=mono_sd,n2i=mono_n,
                          data=all.platforms,measure="SMDH",append=F)
names(SMDH_effect_sizes) = c("SMDH","SMDH_Var")

all.platforms = cbind(all.platforms,as.data.frame(ROM_effect_sizes),as.data.frame(SMDH_effect_sizes))
all.platforms$unique_id = factor(seq(1:nrow(all.platforms)))

all.platforms.complete = all.platforms[which(complete.cases(all.platforms)),]
########################################################
# Meta-Analysis ----------------------------------------


####################
# with lmer---------
# grand mean
lmer.grand.mean = lmer(all.platforms.complete$RR ~ (1|all.platforms.complete$species),weights=(1/all.platforms.complete$RR_Var))
summary(lmer.grand.mean)

#lmer with moderators
lmer.mods.full = lmer(all.platforms.complete$RR ~ all.platforms.complete$platform + all.platforms.complete$exp + all.platforms.complete$SR -1 +
                     (1|all.platforms.complete$species),weights=(1/all.platforms.complete$RR_Var))
lmer.mods.full
summary(lmer.mods.full)
anova(lmer.mods.full)

lmer.mods.no.SR = lmer(all.platforms.complete$RR ~ all.platforms.complete$platform + all.platforms.complete$exp  -1 +
                      (1|all.platforms.complete$species),weights=(1/all.platforms.complete$RR_Var))
summary(lmer.mods.no.SR)
anova(lmer.mods.no.SR)

rma.grand.mean.no.STlmer = lmer(all.platforms.complete$RR ~ all.platforms.complete$platform + all.platforms.complete$exp + all.platforms.complete$SR -1 +
                               (1|all.platforms.complete$species),weights=(1/all.platforms.complete$RR_Var))

##########################
# with metafor-----------

################
#grand mean-----
rma.grand.mean = rma(yi=RR,vi=RR_Var,data=all.platforms)
rma.grand.mean

# platform effect
rma.platforms = rma(yi=RR,vi=RR_Var,data=all.platforms,
                    mods=~ platform - 1)
rma.platforms

#exp effect within each platform
rma.inventories.exp = rma(yi=RR,vi=RR_Var,data=all.platforms,
                    mods=~ exp -1,subset= platform == "inventories")
rma.experiments.exp = rma(yi=RR,vi=RR_Var,data=all.platforms,
                          mods=~ exp-1,subset= platform == "experiments")
rma.exploratories.exp = rma(yi=RR,vi=RR_Var,data=all.platforms,
                          mods=~ exp-1,subset= platform == "exploratories")
rma.inventories.exp
rma.experiments.exp
rma.exploratories.exp

#SR effect within exp
rma.inventories.exp.species = rma(yi=RR,vi=RR_Var,data=all.platforms,
                          mods=~ factor(SR),subset= platform == "inventories")
rma.experiments.exp.species = rma(yi=RR,vi=RR_Var,data=all.platforms,
                          mods=~ factor(SR),subset= platform == "experiments")
rma.exploratories.exp.species = rma(yi=RR,vi=RR_Var,data=all.platforms,
                            mods=~ factor(SR),subset= platform == "exploratories")
rma.inventories.exp.species
rma.experiments.exp.species
rma.exploratories.exp.species


##############################################
# Figures ------------------------------------

# one figure for each research approach - otherwise its too messy
plot.data.inven = all.platforms[,c("platform","species","exp","SR","RR","RR_Var")][which(all.platforms$platform %in% "inventories"),]
plot.data.exper = all.platforms[,c("platform","species","exp","SR","RR","RR_Var")][which(all.platforms$platform %in% "experiments"),]
plot.data.explo = all.platforms[,c("platform","species","exp","SR","RR","RR_Var")][which(all.platforms$platform %in% "exploratories"),]

# inventories
plot.data.inven = plot.data.inven[complete.cases(plot.data.inven),]
plot.data.inven = plot.data.inven[order(plot.data.inven$RR),]
plot.data.inven$uniqueID = seq(1:nrow(plot.data.inven))

svg(filename="figures\\forest_plot_inventories.svg",height=10,width=6)  
ggplot(data=plot.data.inven)  +
    #geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var)),color="#FFC838") +
    geom_point(aes(x=factor(uniqueID),y=RR,size=1/RR_Var),fill="#FFC838",alpha=0.5,shape=21) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + theme_cowplot() + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),legend.position = "none") +
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1)) +
    geom_hline(aes(yintercept=-1), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=-0.5), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.5), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=1), lty=2,size=0.2,colour="grey") +
    ylab("\nResponse ratio of species growth\nin mixtures vs. monocultures") +
    ggtitle("a) Inventory approach")
graphics.off()

# experiments
plot.data.exper = plot.data.exper[complete.cases(plot.data.exper),]
plot.data.exper = plot.data.exper[order(plot.data.exper$RR),]
plot.data.exper$uniqueID = seq(1:nrow(plot.data.exper))

svg(filename="figures\\forest_plot_experiments.svg",height=4.6,width=6)  
ggplot(data=plot.data.exper)  +
    #geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var)),color="#FFC838") +
    geom_point(aes(x=factor(uniqueID),y=RR,size=1/RR_Var),fill="#4090DB",alpha=0.5,shape=21) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + theme_cowplot() + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),legend.position = "none") +
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1)) +
    geom_hline(aes(yintercept=-1), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=-0.5), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.5), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=1), lty=2,size=0.2,colour="grey") +
    ylab("\nResponse ratio of species growth\nin mixtures vs. monocultures") +
    ggtitle("b) Experimental approach")
graphics.off()

# exploratories
plot.data.explo = plot.data.explo[complete.cases(plot.data.explo),]
plot.data.explo = plot.data.explo[order(plot.data.explo$RR),]
plot.data.explo$uniqueID = seq(1:nrow(plot.data.explo))

svg(filename="figures\\forest_plot_exploratories.svg",height=4.6,width=6)  
ggplot(data=plot.data.explo)  +
    #geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var)),color="#FFC838") +
    geom_point(aes(x=factor(uniqueID),y=RR,size=1/RR_Var),fill="#6DC993",alpha=0.5,shape=21) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + theme_cowplot() + 
    theme(axis.text.y = element_blank(),
          axis.title.y = element_blank(),legend.position = "none") +
    scale_y_continuous(breaks = c(-1,-0.5,0,0.5,1)) +
    geom_hline(aes(yintercept=-1), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=-0.5), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.5), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=1), lty=2,size=0.2,colour="grey") +
    ylab("\nResponse ratio of species growth\nin mixtures vs. monocultures") +
    ggtitle("c) Exploratory approach")
graphics.off()

#grand means
plot.grand.means = rbind( data.frame("platform"="grand mean","species"="grand mean","exp"="grand mean", "SR"="grand mean",
                                     "RR"=rma.grand.mean$b, "RR_Var" = rma.grand.mean$se),
                          empty.row,
                          data.frame("platform"="exploratories","species"="exploratories","exp"="exploratories", "SR"="grand mean",
                                     "RR"=rma.platforms$b[which(row.names(rma.platforms$b) %in% "platformexploratories")], "RR_Var" = rma.platforms$se[which(row.names(rma.platforms$b) %in% "platformexploratories")] ),
                          data.frame("platform"="experiments","species"="experiments","exp"="experiments", "SR"="grand mean",
                                     "RR"=rma.platforms$b[which(row.names(rma.platforms$b) %in% "platformexperiments")], "RR_Var" = rma.platforms$se[which(row.names(rma.platforms$b) %in% "platformexperiments")] ),
                          data.frame("platform"="inventories","species"="inventories","exp"="inventories", "SR"="grand mean",
                                     "RR"=rma.platforms$b[which(row.names(rma.platforms$b) %in% "platforminventories")], "RR_Var" = rma.platforms$se[which(row.names(rma.platforms$b) %in% "platforminventories")] ))

rownames(plot.grand.means) = seq(1:nrow(plot.grand.means))                   
plot.grand.means$uniqueID = factor(rownames(plot.grand.means),levels = rownames(plot.grand.means))                   
plot.grand.means$replicates = c(rma.platforms$k,
                                NA,
                                sum(rma.platforms$X[,which(colnames(rma.platforms$X) %in% "platformexperiments")]),
                                sum(rma.platforms$X[,which(colnames(rma.platforms$X) %in% "platformexploratories")]),
                                sum(rma.platforms$X[,which(colnames(rma.platforms$X) %in% "platforminventories")]))


svg(filename="figures\\forest_plot_grand_means.svg",height=3,width=6)  
ggplot(data=plot.grand.means) + 
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=platform),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=platform),size=6) +
    annotate("text",x = factor(plot.grand.means$uniqueID), y = plot.grand.means$RR +1.96* (plot.grand.means$RR_Var) + 0.035, label = plot.grand.means$replicates,size =3.5) + 
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + theme_cowplot()+
    theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Research Approach",
                       breaks = c("inventories", "experiments", "exploratories","grand mean"),
                       values=c("#FB3958","#4090DB","#6DC993","#FFC838")) +
    ylab("Response ratio of species growth\nin mixtures vs. monocultures") +
    ylim(c(-0.1,0.4)) +
    geom_hline(aes(yintercept=-0.1), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.1), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.2), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.3), lty=2,size=0.2,colour="grey") +
    geom_hline(aes(yintercept=0.4), lty=2,size=0.2,colour="grey") 
graphics.off()







# B: separate for each platform with exp as moderator ------------------

plot.data1 = all.platforms[,c("platform","species","exp","SR","RR","RR_Var")] 

plot.data1 = plot.data1[-which(is.na(plot.data1$RR) | is.na(plot.data1$RR_Var)),]

plot.data.experiments = subset(plot.data1, platform %in% "experiments")
plot.data.experiments = plot.data.experiments[order(plot.data.experiments$exp,plot.data.experiments$RR,decreasing=F),]
plot.data.experiments$exp = factor(plot.data.experiments$exp )

plot.data.exploratories = subset(plot.data1, platform %in% "exploratories")
plot.data.exploratories = plot.data.exploratories[order(plot.data.exploratories$exp,plot.data.exploratories$RR,decreasing=F),]
plot.data.exploratories$exp = factor(plot.data.exploratories$exp )

plot.data.inventories = subset(plot.data1, platform %in% "inventories")
plot.data.inventories = plot.data.inventories[order(plot.data.inventories$exp,plot.data.inventories$RR,decreasing=F),]
plot.data.inventories$exp = factor(plot.data.inventories$exp )

#insert empty row between different levels of exp
insert.empty.rows = function(data.to.split,factor.to.split){
    empty.row = data.to.split[NA,][1,]
    factor.table = table(data.to.split[,which(names(data.to.split) %in% factor.to.split)])
    
    new.table = data.to.split[0,]
    
    start.row = 1
    end.row = 0
    
    for(i in factor.table){
        
        end.row = end.row + i
        new.table= rbind(new.table,
                         data.to.split[start.row : end.row,],
                         empty.row)
        start.row=start.row + i
    }
    return(new.table)
}


plot.data.experiments = insert.empty.rows(plot.data.experiments,"exp")
plot.data.exploratories = insert.empty.rows(plot.data.exploratories,"exp")
plot.data.inventories = insert.empty.rows(plot.data.inventories,"exp")

plot.data.experiments$uniqueID = seq(1:nrow(plot.data.experiments))
plot.data.exploratories$uniqueID = seq(1:nrow(plot.data.exploratories))
plot.data.inventories$uniqueID = seq(1:nrow(plot.data.inventories))

# create dataframe for plotting grand means

rma.experiments = rma(yi=RR,vi=RR_Var,data=all.platforms,subset=(platform == "experiments"))
rma.exploratories = rma(yi=RR,vi=RR_Var,data=all.platforms,subset=(platform == "exploratories"))
rma.inventories = rma(yi=RR,vi=RR_Var,data=all.platforms,subset=(platform == "inventories"))


plot.experiments.grand.means = data.frame(
    "exp"=c("grand mean",NA,row.names(rma.experiments.exp$b)),
    "RR" =c(rma.experiments$b,NA,rma.experiments.exp$b),
    "RR_Var" =c(rma.experiments$se,NA,rma.experiments.exp$se))
plot.exploratories.grand.means = data.frame(
    "exp"=c("grand mean",NA,row.names(rma.exploratories.exp$b)),
    "RR" =c(rma.exploratories$b,NA,rma.exploratories.exp$b),
    "RR_Var" =c(rma.exploratories$se,NA,rma.exploratories.exp$se))
plot.inventories.grand.means = data.frame(
    "exp"=c("grand mean",NA,row.names(rma.inventories.exp$b)),
    "RR" =c(rma.inventories$b,NA,rma.inventories.exp$b),
    "RR_Var" =c(rma.inventories$se,NA,rma.inventories.exp$se))
plot.experiments.grand.means$uniqueID = seq(1:nrow(plot.experiments.grand.means))
plot.exploratories.grand.means$uniqueID = seq(1:nrow(plot.exploratories.grand.means))
plot.inventories.grand.means$uniqueID = seq(1:nrow(plot.inventories.grand.means))

plot.experiments.grand.means$replicates = c(rma.experiments$k,NA,colSums(as.data.frame(rma.experiments.exp$X))) 
plot.exploratories.grand.means$replicates = c(rma.exploratories$k,NA,colSums(as.data.frame(rma.exploratories.exp$X))) 
plot.inventories.grand.means$replicates = c(rma.inventories$k,NA,colSums(as.data.frame(rma.inventories.exp$X))) 

#############
# experiments

svg(filename="figures\\experiments_exp_as_mod_RR.svg",height=10,width=5)  
ggplot(data=plot.data.experiments) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp)) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp)) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2))) +
    
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name="Tree diversity experiment",
                       breaks = c(rev(levels(plot.data.experiments$exp))),
                       values= c("#5971C9","#8ACDFF","#5F98D8","#85E0E7","#007BB5")) +
    ggtitle("b) Experimental approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


# grand means of all platforms
svg(filename="figures\\experiments_exp_as_mod_grand_means_RR.svg",height=4,width=4)  
ggplot(data=plot.experiments.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp),size=6) +
    annotate("text",x = factor(plot.experiments.grand.means$uniqueID), y = plot.experiments.grand.means$RR +1.96* (plot.experiments.grand.means$RR_Var) + 0.035, label = plot.experiments.grand.means$replicates,size =3.5) + 
    
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2)))+
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Tree diversity experiment",
                       breaks = c(rev(row.names(rma.experiments.exp$b)),"grand mean"),
                       values= c("#5971C9","#8ACDFF","#5F98D8","#85E0E7","#007BB5","#AA0000")) +
    ggtitle("b) Experimental approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()

#############
# exploratories

svg(filename="figures\\exploratories_exp_as_mod_RR.svg",height=10,width=5)  
ggplot(data=plot.data.exploratories) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp)) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp)) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2))) +
    
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.data.exploratories$exp))),
                       values= c("#10580A","#21b014","#98A95D","#71D45A","#186F22","#91FB5C")) +
    ggtitle("c) Exploratory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


# grand means of all platforms
svg(filename="figures\\exploratories_exp_as_mod_grand_means_RR.svg",height=4,width=4)  
ggplot(data=plot.exploratories.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp),size=6) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    annotate("text",x = factor(plot.exploratories.grand.means$uniqueID), y = plot.exploratories.grand.means$RR +1.96* (plot.exploratories.grand.means$RR_Var) + 0.05, label = plot.exploratories.grand.means$replicates,size =3.5) + 
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2)))+
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name = "Forest type",
                       breaks = c(rev(row.names(rma.exploratories.exp$b)),"grand mean"),
                       values= c("#10580A","#21b014","#98A95D","#71D45A","#186F22","#91FB5C","#AA0000")) +
    ggtitle("c) Exploratory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


#############
# inventories

svg(filename="figures\\inventories_exp_as_mod_RR.svg",height=10,width=5)  
ggplot(data=plot.data.inventories) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp)) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp)) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2))) +
    
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(levels(plot.data.inventories$exp))),
                       values= c(rep(c("#634C00", "#957201","#C79801","#F9BE02","#FAD14D","#E4EAC7"),4),"#634C00", "#957201")) +
    ggtitle("a) Inventory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()


# grand means of all platforms
svg(filename="figures\\inventories_exp_as_mod_grand_means_RR.svg",height=9,width=6)  
ggplot(data=plot.inventories.grand.means) + 
    
    geom_linerange(aes(x=factor(uniqueID),ymin=RR - (1.96 * RR_Var), ymax = RR + 1.96* (RR_Var),colour=exp),size=3) +
    geom_point(aes(x=factor(uniqueID),y=RR,colour=exp),size=6) +
    geom_hline(aes(yintercept=0), lty=2,size=0.5,colour="black") +
    annotate("text",x = factor(plot.inventories.grand.means$uniqueID), y = plot.inventories.grand.means$RR +1.96* (plot.inventories.grand.means$RR_Var) + 0.2, label = plot.inventories.grand.means$replicates,size =3.5) + 
    
    coord_flip() + scale_y_continuous(labels=trans_format("exp",comma_format(digits=2)))+
    theme(panel.background = element_blank(), axis.text.y = element_blank(),axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),legend.position="left") +
    scale_color_manual(name ="Forest type",
                       breaks = c(rev(row.names(rma.inventories.exp$b)),"grand mean"),
                       values= c(rep(c("#634C00", "#957201","#C79801","#F9BE02","#FAD14D","#E4EAC7"),4),"#634C00", "#957201","#AA0000")) +
    ggtitle("a) Inventory approach") +
    ylab("\nResponse ratio of species growth\n in mixtures vs. monocultures")
graphics.off()

