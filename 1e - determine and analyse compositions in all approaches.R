### FunDiv Synthesis:
### Diversity ~ Growth across platforms
### By: Stephan Kambach
### Date: 26.04.2017
### check with mixtures are comparable between all the platforms

if(!require(data.table)){install.packages("data.table")}
if(!require(cluster)){install.packages("cluster")}
if(!require(optmatch)){install.packages("optmatch")}
if(!require(metafor)){install.packages("metafor")}
if(!require(ggplot2)){install.packages("ggplot2")}
if(!require(scales)){install.packages("scales")}
if(!require(lmerTest)){install.packages("lmerTest")}
if(!require(corrplot)){install.packages("corrplot")}
if(!require(extrafont)){install.packages("extrafont")}

library(data.table)
library(cluster)
library(optmatch)
library(metafor)
library(ggplot2)
library(lmerTest)
library(corrplot)
library(extrafont)
font_import(pattern ="[T/t]imes")
y
loadfonts(device="win")
fonts()

#########################################################
# PART 1 - determine comparable species compositions----
rm(list=ls())
gc()

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")

inventories = read.table("Inventories\\NFI with French data\\all_species_compositions.csv",header=T)
experiments = read.table("Experiments\\all_species_compositions.csv",header=T)
exploratories = read.table("Exploratories\\all_species_compositions.csv",header=T)
experiments$in_experiments = "yes"
exploratories$in_exploratories = "yes"

all.platforms = merge(experiments,exploratories,by.x="x",by.y="x", all.x=T,all.y=T)
all.platforms = merge(all.platforms,inventories,by.x = "x", by.y = "species_combinations", all.x = T)

in.inventories.experiments = all.platforms[which(!is.na(all.platforms$n_replicates) & !(is.na(all.platforms$in_experiments))),][,c("x","n_replicates","in_experiments")]
in.inventories.experiments = in.inventories.experiments[in.inventories.experiments$n_replicates >= 1,]
in.inventories.exploratories = all.platforms[which(!is.na(all.platforms$n_replicates) & !(is.na(all.platforms$in_exploratories))),][,c("x","n_replicates","in_exploratories")]
in.inventories.exploratories = in.inventories.exploratories[in.inventories.exploratories$n_replicates >= 1,]
in.experiments.exploratories = all.platforms[which(!is.na(all.platforms$in_experiments) & !(is.na(all.platforms$in_exploratories))),][,c("x","in_experiments","in_exploratories")]

# check to keep only monocultures that are also in mixtures and only mixtures for which monocultures are available

monos = unique(gsub(pattern ="_.*",replace="",in.inventories.experiments$x))
for(i in 1:length(monos)){
    if(length(grep(pattern = monos[i], x = in.inventories.experiments$x)) == 1){
        in.inventories.experiments = in.inventories.experiments[-grep(pattern = monos[i], x = in.inventories.experiments$x),]}
}

monos = unique(gsub(pattern ="_.*",replace="",in.inventories.exploratories$x))
for(i in 1:length(monos)){
    if(length(grep(pattern = monos[i], x = in.inventories.exploratories$x)) == 1){
        in.inventories.exploratories = in.inventories.exploratories[-grep(pattern = monos[i], x = in.inventories.exploratories$x),]}
}

# create table with all the mono-mix for comparison
composition.summary = data.frame("composition" = unique(c(as.character(in.inventories.experiments$x),
                                                        as.character(in.inventories.exploratories$x),
                                                        as.character(in.experiments.exploratories$x))),
                                 "inv_exper" = "no",
                                 "inv_explo" = "no",
                                 "exper_explo" = "no")
composition.summary$inv_exper = as.character(composition.summary$inv_exper)
composition.summary$inv_explo = as.character(composition.summary$inv_explo)
composition.summary$exper_explo = as.character(composition.summary$exper_explo)

composition.summary$inv_exper[which(composition.summary$composition %in% as.character(in.inventories.experiments$x))] = "yes"
composition.summary$inv_explo[which(composition.summary$composition %in% as.character(in.inventories.exploratories$x))] = "yes"
composition.summary$exper_explo[which(composition.summary$composition %in% as.character(in.experiments.exploratories$x))] = "yes"

write.table(composition.summary,"compositions_comparable_across_approaches.txt",sep="\t",quote=F,row.names=F)

##########################################
# PART 2: Calculate ES for experiments----

rm(list=ls())
gc()

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")

comparable.compositions = read.table("compositions_comparable_across_approaches.txt",sep = "\t",header=T)
str(comparable.compositions)

correct.names = read.table("species_names_correction.txt", header = T, sep="\t")
names(correct.names)[1] = "code_old"
correct.names$code_old = as.character(correct.names$code_old)
correct.names$code_new = as.character(correct.names$code_new)

#functions names(data.to.order) = exp,species,SR,growth_mean,growth_sd,growth_sd
re.order.function = function(data.to.order){
    mix.data = subset(data.to.order, SR >= 2)
    mono.data = subset(data.to.order, SR == 1)
    mono.data[,SR:=NULL]
    mix.data[,exp:=NULL]
    merge.data = merge(mono.data,mix.data,by.x="species",by.y="species",all.y=T)
    names(merge.data) = c("species","exp","mono_mean","mono_sd","mono_n","SR","mix_mean","mix_sd","mix_n")
    return(merge.data)
}

# input should be as character vector
rename.species.compositions = function(species.compositions){
    for(i in 1:length(species.compositions)){
        species.comp.temp = species.compositions[i]
        species.temp = strsplit(species.comp.temp, split = "_|\\+")[[1]]
        for(j in 1:length(species.temp)){
            if(species.temp[j] %in% correct.names$code_old){
            species.temp[j] = correct.names$code_new[which(correct.names$code_old %in% species.temp[j])]}
        }
        species.compositions[i] = paste(sort(species.temp), collapse = "_")
    }
    return(species.compositions)
}

rename.species = function(species){
    for(i in 1:length(species)){
        if(species[i] %in% correct.names$code_old){
            species[i] = correct.names$code_new[which(correct.names$code_old %in% species[i])]}

    }
    return(species)
}

####### experiments 1 ######## growth as basal diameter after 5 years, also height available
experiment1 = data.table(read.table("Experiments\\update_from_Simon\\forbio.zed.plot.csv",header=T,sep=","))
experiment1$SR[which(experiment1$SR > 3)] = 4

# rename
experiment1$sp.comp = rename.species.compositions(as.character(experiment1$sp.comp))
experiment1$sp.code = rename.species(as.character(experiment1$sp.code))

# only include plots that are matched in inventory and exploratory
experiment1.for.inv.comp = experiment1[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]
experiment1.for.explo.comp = experiment1[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$exper_explo %in% "yes")])]

# for inv
experiment1.summarized.inv = experiment1.for.inv.comp[,.("exp"="forbio.zed",
                                        "growth_mean"=mean(basaldiam.mm,na.rm=T),"growth_sd"=sd(basaldiam.mm,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment1.summarized.inv)[1] = "species"                                       
experiment1.summarized.inv= re.order.function(data.to.order = experiment1.summarized.inv)

# for explo
experiment1.summarized.explo = experiment1.for.explo.comp[,.("exp"="forbio.zed",
                                                         "growth_mean"=mean(basaldiam.mm,na.rm=T),"growth_sd"=sd(basaldiam.mm,na.rm=T),"growth_n"=length(unique(X))),
                                                      by=.(sp.code,SR)]
names(experiment1.summarized.explo)[1] = "species"                                       
experiment1.summarized.explo= re.order.function(data.to.order = experiment1.summarized.explo)

####### experiments 2 ######## growth as basal area.cm2 after 8 years, also height available
experiment2 = data.table(read.table("Experiments\\update_from_Simon\\kreinitz.plot.csv",header=T,sep=","))
experiment2$SR[which(experiment2$SR > 3)] = 4

# rename
experiment2$sp.comp = rename.species.compositions(as.character(experiment2$sp.comp))
experiment2$sp.code = rename.species(as.character(experiment2$sp.code))

# only include plots that are matched in inventory and exploratory
experiment2.for.inv.comp = experiment2[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]
experiment2.for.explo.comp = experiment2[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$exper_explo %in% "yes")])]

# for inv
experiment2.summarized.inv = experiment2.for.inv.comp[,.("exp"="kreinitz",
                                        "growth_mean"=mean(basalarea.cm2,na.rm=T),"growth_sd"=sd(basalarea.cm2,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment2.summarized.inv)[1] = "species"                                       
experiment2.summarized.inv= re.order.function(data.to.order = experiment2.summarized.inv)
experiment2.summarized.inv = experiment2.summarized.inv[species %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]

# for explo
experiment2.summarized.explo = experiment2.for.explo.comp[,.("exp"="kreinitz",
                                                         "growth_mean"=mean(basalarea.cm2,na.rm=T),"growth_sd"=sd(basalarea.cm2,na.rm=T),"growth_n"=length(unique(X))),
                                                      by=.(sp.code,SR)]
names(experiment2.summarized.explo)[1] = "species"                                       
experiment2.summarized.explo= re.order.function(data.to.order = experiment2.summarized.explo)
experiment2.summarized.explo = experiment2.summarized.explo[species %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_explo %in% "yes")])]

####### experiments 3 ######## growth as dbh.mmm after 12 years, also basal area.cm2 and height available
experiment3 = data.table(read.table("Experiments\\update_from_Simon\\satakunta.plot.csv",header=T,sep=","))
experiment3$SR[which(experiment3$SR > 3)] = 4

# rename
experiment3$sp.comp = rename.species.compositions(as.character(experiment3$sp.comp))
experiment3$sp.code = rename.species(as.character(experiment3$sp.code))

# only include plots that are matched in inventory and exploratory
experiment3.for.inv.comp = experiment3[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]
experiment3.for.explo.comp = experiment3[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$exper_explo %in% "yes")])]

# for inv
experiment3.summarized.inv = experiment3.for.inv.comp[,.("exp"="satakunta",
                                        "growth_mean"=mean(dbh.mm,na.rm=T),"growth_sd"=sd(dbh.mm,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment3.summarized.inv)[1] = "species"                                       
experiment3.summarized.inv= re.order.function(data.to.order = experiment3.summarized.inv)

#for explo
experiment3.summarized.explo = experiment3.for.explo.comp[,.("exp"="satakunta",
                                                         "growth_mean"=mean(dbh.mm,na.rm=T),"growth_sd"=sd(dbh.mm,na.rm=T),"growth_n"=length(unique(X))),
                                                      by=.(sp.code,SR)]
names(experiment3.summarized.explo)[1] = "species"                                       
experiment3.summarized.explo= re.order.function(data.to.order = experiment3.summarized.explo)

####### experiments 4 ######## growth as basaldiameter.mean after 5 years, also height.mean available
experiment4 = data.table(read.table("Experiments\\update_from_Simon\\ged.plot.csv",header=T,sep=","))
experiment4$SR[which(experiment4$SR > 3)] = 4

# rename
experiment4$sp.comp = rename.species.compositions(as.character(experiment4$sp.comp))
experiment4$sp.code = rename.species(as.character(experiment4$sp.code))

experiment4.for.inv.comp = experiment4[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]
experiment4.for.explo.comp = experiment4[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$exper_explo %in% "yes")])]

# for inv
experiment4.summarized.inv = experiment4.for.inv.comp[,.("exp"="forbio.ged",
                                        "growth_mean"=mean(basaldiameter.mean,na.rm=T),"growth_sd"=sd(basaldiameter.mean,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment4.summarized.inv)[1] = "species"                                       
experiment4.summarized.inv= re.order.function(data.to.order = experiment4.summarized.inv)

# for explo
experiment4.summarized.explo = experiment4.for.explo.comp[,.("exp"="forbio.ged",
                                                         "growth_mean"=mean(basaldiameter.mean,na.rm=T),"growth_sd"=sd(basaldiameter.mean,na.rm=T),"growth_n"=length(unique(X))),
                                                      by=.(sp.code,SR)]
names(experiment4.summarized.explo)[1] = "species"                                       
experiment4.summarized.explo= re.order.function(data.to.order = experiment4.summarized.explo)

####### experiments 5 ######## growth as basalarea.mean after 11 years
experiment5 = data.table(read.table("Experiments\\update_from_Simon\\biotree\\biotree.subplot.csv",header=T,sep=","))
experiment5$SR[which(experiment5$SR > 3)] = 4

# rename
experiment5$sp.comp = rename.species.compositions(as.character(experiment5$sp.comp))
experiment5$sp.code = rename.species(as.character(experiment5$sp.code))

experiment5.for.inv.comp = experiment5[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]
experiment5.for.explo.comp = experiment5[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$exper_explo %in% "yes")])]

# for inv
experiment5.summarized.inv = experiment5.for.inv.comp[,.("exp"="biotree",
                                        "growth_mean"=mean(basalarea.cm2,na.rm=T),"growth_sd"=sd(basalarea.cm2,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment5.summarized.inv)[1] = "species"                                       
experiment5.summarized.inv= re.order.function(data.to.order = experiment5.summarized.inv)

# for explo
experiment5.summarized.explo = experiment5.for.explo.comp[,.("exp"="biotree",
                                                         "growth_mean"=mean(basalarea.cm2,na.rm=T),"growth_sd"=sd(basalarea.cm2,na.rm=T),"growth_n"=length(unique(X))),
                                                      by=.(sp.code,SR)]
names(experiment5.summarized.explo)[1] = "species"                                       
experiment5.summarized.explo= re.order.function(data.to.order = experiment5.summarized.explo)

####### experiments 6 ######## growth as height in 2014 after 6 years
experiment6 = data.table(read.table("Experiments\\update_from_Simon\\gro.sp.orphee_new_2015.csv",sep=",",header=T,row.names=1))
experiment6$richness[which(experiment6$richness > 3)] = 4
experiment6$X = seq(1:nrow(experiment6))
names(experiment6)[which(names(experiment6) == "treatment")] = "sp.comp"
names(experiment6)[which(names(experiment6) == "SP")] = "sp.code"
names(experiment6)[which(names(experiment6) == "richness")] = "SR"

# rename
experiment6$sp.comp = rename.species.compositions(as.character(experiment6$sp.comp))
experiment6$sp.code = rename.species(as.character(experiment6$sp.code))

experiment6.for.inv.comp = experiment6[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$inv_exper %in% "yes")])]
experiment6.for.explo.comp = experiment6[sp.comp %in% as.character(comparable.compositions$composition[which(comparable.compositions$exper_explo %in% "yes")])]

# for inv
experiment6.summarized.inv = experiment6.for.inv.comp[,.("exp"="biotree",
                                                         "growth_mean"=mean(mean_plot,na.rm=T),"growth_sd"=sd(mean_plot,na.rm=T),"growth_n"=length(unique(X))),
                                                      by=.(sp.code,SR)]
names(experiment6.summarized.inv)[1] = "species"                                       
experiment6.summarized.inv= re.order.function(data.to.order = experiment6.summarized.inv)

# for explo
experiment6.summarized.explo = experiment6.for.explo.comp[,.("exp"="biotree",
                                                             "growth_mean"=mean(mean_plot,na.rm=T),"growth_sd"=sd(mean_plot,na.rm=T),"growth_n"=length(unique(X))),
                                                          by=.(sp.code,SR)]
names(experiment6.summarized.explo)[1] = "species"                                       
experiment6.summarized.explo= re.order.function(data.to.order = experiment6.summarized.explo)

####### merge all experiment tables 

all.experiments.ordered.inv = rbind(experiment1.summarized.inv,experiment2.summarized.inv,experiment3.summarized.inv,
                                    experiment4.summarized.inv,experiment5.summarized.inv,experiment6.summarized.inv)

all.experiments.ordered.explo = rbind(experiment1.summarized.explo,experiment2.summarized.explo,experiment3.summarized.explo,
                                      experiment4.summarized.explo,experiment5.summarized.explo,experiment6.summarized.explo)

all.experiments.ordered.inv$platform = "experiments"
all.experiments.ordered.explo$platform = "experiments"

write.table(all.experiments.ordered.inv, "Experiments\\all_experiments_re-arranged_for_comparison_with_inventories.csv",sep="\t",dec=".",quote=F,row.names=F)
write.table(all.experiments.ordered.explo, "Experiments\\all_experiments_re-arranged_for_comparison_with_exploratories.csv",sep="\t",dec=".",quote=F,row.names=F)

############################################
# PART 3: Calculate ES for exploratories----

rm(list=ls())
gc()

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")

files.plots = list.files("Exploratories\\Plot_descriptors\\")
files.growth = list.files("Exploratories\\Growth_data\\")
plot.alldata = read.table("Exploratories\\alldata_sDiv_FunDiv_exploratories20150611.txt",sep="\t",header=T)
rename_species = read.table("species_names_correction.txt",header=T, sep="\t")
names(rename_species)[1]  = "code_old"
comparable.compositions = read.table("compositions_comparable_across_approaches.txt",sep = "\t",header=T)

growth.es.final.inv = data.frame("SPCode"=NA,"exp"=NA,"mono_mean"=NA,"mono_sd"=NA,"mono_n"=NA,
                             "SR"=NA,"mix_mean"=NA,"mix_sd"=NA,"mix_n"=NA)[0,]
growth.es.final.exper = data.frame("SPCode"=NA,"exp"=NA,"mono_mean"=NA,"mono_sd"=NA,"mono_n"=NA,
                             "SR"=NA,"mix_mean"=NA,"mix_sd"=NA,"mix_n"=NA)[0,]


# a table to extract the growth per species and plot for the supplement
growth.table.for.supplement.inv = data.frame(Increase_species_ba = NA, PlotID = NA, SpeciesCode= NA)[0,]
growth.table.for.supplement.exper = data.frame(Increase_species_ba = NA, PlotID = NA, SpeciesCode= NA)[0,]

for(i in 1:6){
    plot.raw = read.csv(paste("Exploratories\\Plot_descriptors\\",files.plots[i],sep=""))
    growth.raw = read.csv(paste("Exploratories\\Growth_data\\",files.growth[i],sep=""))
    
    #rename species codes
    plot.raw$SpeciesCode = as.character(plot.raw$SpeciesCode)
    
    for(j in 1: nrow(plot.raw)){
        if(plot.raw$SpeciesCode[j] %in% rename_species$code_old){
            plot.raw$SpeciesCode[j] = as.character(rename_species$code_new[which(rename_species$code_old %in%  plot.raw$SpeciesCode[j])])
        }}
    
    plot.raw = subset(plot.raw, Target == 1)
    plot.raw$UNIQUE_ID = paste(plot.raw$PlotID, plot.raw$TreeNumber,sep="")
    
    #calculate the mean growth rate over 1990-till end
    growth.mean = data.frame("PlotID" = growth.raw$PlotID,
                             "TreeNumber" = growth.raw$TreeNumber, 
                             "UNIQUE_ID" = paste(growth.raw$PlotID,growth.raw$TreeNumber,sep=""),
                             "Increase_ba" =  ((rowSums(growth.raw[,c(3 : which(names(growth.raw) %in% "RG_2010"))],na.rm=TRUE)^2)*pi) - ((rowSums(growth.raw[,c(3 : which(names(growth.raw) %in% "RG_2000"))],na.rm=TRUE)^2)*pi))
    
    growth.mean = subset(growth.mean, UNIQUE_ID %in% plot.raw$UNIQUE_ID)
    
    #Attach species code to each species
    growth.mean = merge(growth.mean, plot.raw[,c("PlotID","TreeNumber","SpeciesCode")],
                        by.x=c("PlotID","TreeNumber"),by.y=c("PlotID","TreeNumber"),
                        all.x=T, all.y=F)
    growth.mean = na.omit(growth.mean)
    
    #calculate the mean growth rate per species over all trees per plot
    growth.mean = data.table(growth.mean)
    
    growth.summary = growth.mean[,.(Increase_species_ba = mean(Increase_ba)),by=.(PlotID,SpeciesCode)]
    growth.table.for.supplement.inv = rbind(growth.table.for.supplement.inv,growth.summary)
    growth.table.for.supplement.exper = rbind(growth.table.for.supplement.exper,growth.summary)
    
    growth.summary = merge(growth.summary,plot.alldata[,c("PlotID","Restricted_Species_Richness")],
                           by.x="PlotID",by.y="PlotID",
                           all.x=TRUE,all.y=FALSE)
    
    #rename SR and cut to 4
    growth.summary$SR = growth.summary$Restricted_Species_Richness
    growth.summary$SR[which(growth.summary$SR > 3)] = 4
    
    #get the location of the site
    location = sub("_growth_data.csv*","",files.growth[i])
    
    # subset only combinations that also occur in the experiments or inventories
    growth.summary.inv = growth.summary
    growth.summary.exper = growth.summary
    
    plots.to.check = growth.summary[,.(species.combination = paste(sort(SpeciesCode),collapse ="_")),
                                                              by=.(PlotID)]
    for(k in 1:length(plots.to.check$PlotID)){
        species.comp.temp = plots.to.check$species.combination[k]
    
        if(species.comp.temp %in% comparable.compositions$composition){
          if(comparable.compositions$inv_explo[which(comparable.compositions$composition %in% species.comp.temp)] == "yes"){
            
          }else{growth.summary.inv = 
            growth.summary.inv[!(growth.summary.inv$PlotID %in%  plots.to.check$PlotID[k]),]}
        }else{growth.summary.inv = 
          growth.summary.inv[!(growth.summary.inv$PlotID %in%  plots.to.check$PlotID[k]),]}

        
        if(species.comp.temp %in% comparable.compositions$composition){
          if(comparable.compositions$exper_explo[which(comparable.compositions$composition %in% species.comp.temp)] == "yes"){
            
          }else{growth.summary.exper = 
            growth.summary.exper[!(growth.summary.exper$PlotID %in%  plots.to.check$PlotID[k]),]}
        }else{growth.summary.exper = 
          growth.summary.exper[!(growth.summary.exper$PlotID %in%  plots.to.check$PlotID[k]),]}
    }
    
    #calculate growth per Species and Richness level
    growth.summary.inv = growth.summary.inv[,.(exp=location, 
                                       SPECIES_MEAN = mean(Increase_species_ba),
                                       SPECIES_SD = sd(Increase_species_ba),
                                       SPECIES_N = length(unique(Increase_species_ba))),
                                    by=.(SpeciesCode,SR)]
    
    growth.summary.exper = growth.summary.exper[,.(exp=location, 
                                               SPECIES_MEAN = mean(Increase_species_ba),
                                               SPECIES_SD = sd(Increase_species_ba),
                                               SPECIES_N = length(unique(Increase_species_ba))),
                                            by=.(SpeciesCode,SR)]
    
    # remove incomplete rows (with one one replicate)
    #growth.summary.inv = growth.summary.inv[complete.cases(growth.summary.inv),]
    #growth.summary.exper = growth.summary.exper[complete.cases(growth.summary.exper),]
    
    # delete monocultures that have no mixtures to compare
    growth.summary.inv = growth.summary.inv[SpeciesCode %in% names(which(table(growth.summary.inv$SpeciesCode) >= 1))]
    growth.summary.exper = growth.summary.exper[SpeciesCode %in% names(which(table(growth.summary.exper$SpeciesCode) >= 1))]
    
    # order table so that each row contains monoculture and a mixture
    growth.es.inv = merge(growth.summary.inv[growth.summary.inv$SR %in% 1,],growth.summary.inv[!(growth.summary.inv$SR %in% 1),],
                      by.x="SpeciesCode",by.y="SpeciesCode",all.y=T)
    growth.es.inv[,SR.x := NULL]
    growth.es.inv[,exp.y := NULL]
    names(growth.es.inv) = c("SPCode","exp","mono_mean","mono_sd","mono_n","SR","mix_mean","mix_sd","mix_n")
    
    growth.es.exper = merge(growth.summary.exper[growth.summary.exper$SR %in% 1,],growth.summary.exper[!(growth.summary.exper$SR %in% 1),],
                          by.x="SpeciesCode",by.y="SpeciesCode",all.y=T)
    growth.es.exper[,SR.x := NULL]
    growth.es.exper[,exp.y := NULL]
    names(growth.es.exper) = c("SPCode","exp","mono_mean","mono_sd","mono_n","SR","mix_mean","mix_sd","mix_n")
    
    # ones again delete mixtures that have no monoculture to compare
    growth.es.inv = growth.es.inv[!is.na(growth.es.inv$mono_mean) & !is.na(growth.es.inv$mix_mean),]
    growth.es.exper = growth.es.exper[!is.na(growth.es.exper$mono_mean) & !is.na(growth.es.exper$mix_mean),]
    
    #put all together in result dataframe with all sites, species, mixtures
    growth.es.final.inv = rbind(growth.es.final.inv,growth.es.inv)
    growth.es.final.exper = rbind(growth.es.final.exper,growth.es.exper)
}


exploratories.inv = data.frame("species"=growth.es.final.inv$SPCode,"exp"=growth.es.final.inv$exp,"SR"=growth.es.final.inv$SR,
                           "mono_mean"=growth.es.final.inv$mono_mean,"mono_sd"=growth.es.final.inv$mono_sd,"mono_n"=growth.es.final.inv$mono_n,
                           "mix_mean"=growth.es.final.inv$mix_mean,"mix_sd"=growth.es.final.inv$mix_sd,"mix_n"=growth.es.final.inv$mix_n,
                           "platform"= "exploratories")
exploratories.exper = data.frame("species"=growth.es.final.exper$SPCode,"exp"=growth.es.final.exper$exp,"SR"=growth.es.final.exper$SR,
                               "mono_mean"=growth.es.final.exper$mono_mean,"mono_sd"=growth.es.final.exper$mono_sd,"mono_n"=growth.es.final.exper$mono_n,
                               "mix_mean"=growth.es.final.exper$mix_mean,"mix_sd"=growth.es.final.exper$mix_sd,"mix_n"=growth.es.final.exper$mix_n,
                               "platform"= "exploratories")


exploratories.inv = exploratories.inv[!(is.na(exploratories.inv$mono_mean)),]
exploratories.exper = exploratories.exper[!(is.na(exploratories.exper$mono_mean)),]

#rename forest types
exploratories.inv$exp = as.character(exploratories.inv$exp)
exploratories.inv$exp[which(exploratories.inv$exp %in% "Finland")] = "boreal"
exploratories.inv$exp[which(exploratories.inv$exp %in% "Poland")] = "hemiboreal"
exploratories.inv$exp[which(exploratories.inv$exp %in% "Germany")] = "beech"
exploratories.inv$exp[which(exploratories.inv$exp %in% "Romania")] = "mountain_beech"
exploratories.inv$exp[which(exploratories.inv$exp %in% "Italy")] = "thermophilous_deciduous"
exploratories.inv$exp[which(exploratories.inv$exp %in% "Spain")] = "mediterranean_coniferous"

exploratories.exper$exp = as.character(exploratories.exper$exp)
exploratories.exper$exp[which(exploratories.exper$exp %in% "Finland")] = "boreal"
exploratories.exper$exp[which(exploratories.exper$exp %in% "Poland")] = "hemiboreal"
exploratories.exper$exp[which(exploratories.exper$exp %in% "Germany")] = "beech"
exploratories.exper$exp[which(exploratories.exper$exp %in% "Romania")] = "mountain_beech"
exploratories.exper$exp[which(exploratories.exper$exp %in% "Italy")] = "thermophilous_deciduous"
exploratories.exper$exp[which(exploratories.exper$exp %in% "Spain")] = "mediterranean_coniferous"

table(exploratories.inv$exp)
table(exploratories.exper$exp)

write.table(exploratories.inv,"Exploratories\\all_exploratories_re-arranged_for_comparison_with_inventories.csv", quote=F,sep="\t",row.names=F,dec=".")
write.table(exploratories.exper,"Exploratories\\all_exploratories_re-arranged_for_comparison_with_experiments.csv", quote=F,sep="\t",row.names=F,dec=".")

##########################################
# PART 4: Calculate ES for inventories----
rm(list=ls())
gc()

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data")

matched.results = data.table(read.table("Inventories\\NFI with French data\\intermediate_output\\matched_for_meta-analysis.txt",sep="\t",dec=".",head=T))
species.ba.n.growth.clean = data.table(read.table("Inventories\\NFI with French data\\intermediate_output\\species_statistics_clean_with_environment.txt",sep="\t",dec=".",header=T))
comparable.compositions = read.table("compositions_comparable_across_approaches.txt",sep = "\t",header=T)
str(comparable.compositions)

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


# only keep that comparisons that are in experiments or exploratories
matched.results.cut_off.exper = matched.results.cut_off
matched.results.cut_off.explo = matched.results.cut_off

pb <- txtProgressBar(min = 0, max = nrow(matched.results.cut_off), style = 3)

for(i in 1:nrow(matched.results.cut_off)){
    species.comp.temp = paste(sort(species.ba.n.growth.clean$code_new[
        which(species.ba.n.growth.clean$plotcode %in% matched.results.cut_off$mixture[i])]), collapse = "_")
    if(species.comp.temp %in% comparable.compositions$composition){
        if(comparable.compositions$inv_exper[which(comparable.compositions$composition %in% species.comp.temp)] == "yes"){
            
        }else{matched.results.cut_off.exper = 
            matched.results.cut_off.exper[!(matched.results.cut_off.exper$mixture %in% matched.results.cut_off$mixture[i]),]}
    }else{matched.results.cut_off.exper = 
        matched.results.cut_off.exper[!(matched.results.cut_off.exper$mixture %in% matched.results.cut_off$mixture[i]),]}
    
    if(species.comp.temp %in% comparable.compositions$composition){
        if(comparable.compositions$inv_explo[which(comparable.compositions$composition %in% species.comp.temp)] == "yes"){
            
        }else{matched.results.cut_off.explo = 
            matched.results.cut_off.explo[!(matched.results.cut_off.explo$mixture %in% matched.results.cut_off$mixture[i]),]}
    }else{matched.results.cut_off.explo = 
        matched.results.cut_off.explo[!(matched.results.cut_off.explo$mixture %in% matched.results.cut_off$mixture[i]),]}
    
    setTxtProgressBar(pb, i)
}


matched.results.cut_off.exper$counties_compared = paste(matched.results.cut_off.exper$country_mono,matched.results.cut_off.exper$country_mix, sep="_")
matched.results.cut_off.explo$counties_compared = paste(matched.results.cut_off.explo$country_mono,matched.results.cut_off.explo$country_mix, sep="_")

# calculate effect sizes

inventories.effect.sizes.exper = matched.results.cut_off.exper[,.(
    mono_mean = mean(growth_mono), mono_sd = sd(growth_mono),mono_n = length(growth_mono),
    mix_mean = mean(growth_mix), mix_sd = sd(growth_mix), mix_n = length(growth_mix)),
    by= .(code_new,forest_type,SR_mixture,counties_compared)]

inventories.effect.sizes.explo = matched.results.cut_off.explo[,.(
    mono_mean = mean(growth_mono), mono_sd = sd(growth_mono),mono_n = length(growth_mono),
    mix_mean = mean(growth_mix), mix_sd = sd(growth_mix), mix_n = length(growth_mix)),
    by= .(code_new,forest_type,SR_mixture,counties_compared)]

write.table(inventories.effect.sizes.exper, "Inventories\\all_inventories_re-arranged_for_comparison_with_experiments.txt",
            sep="\t",dec=".",quote=F, row.names=F)
write.table(inventories.effect.sizes.explo, "Inventories\\all_inventories_re-arranged_for_comparison_with_exploratories.txt",
            sep="\t",dec=".",quote=F, row.names=F)

#####################################
# PART 5: Analyse output together----

rm(list=ls())
gc()

setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts")

inventories.exper = read.table("data\\Inventories\\all_inventories_re-arranged_for_comparison_with_experiments.txt",sep="\t",header=T)
inventories.explo = read.table("data\\Inventories\\all_inventories_re-arranged_for_comparison_with_exploratories.txt",sep="\t",header=T)
experiments.inv = read.table("data\\Experiments\\all_experiments_re-arranged_for_comparison_with_inventories.csv",sep="\t",header=T)
experiments.explo = read.table("data\\Experiments\\all_experiments_re-arranged_for_comparison_with_exploratories.csv",sep="\t",header=T)
exploratories.inv = read.table("data\\Exploratories\\all_exploratories_re-arranged_for_comparison_with_inventories.csv",sep="\t",header=T)
exploratories.exper = read.table("data\\Exploratories\\all_exploratories_re-arranged_for_comparison_with_experiments.csv",sep="\t",header=T)

comparable.compositions = read.table("data\\compositions_comparable_across_approaches.txt",sep = "\t",header=T)

exploratories.inv$exp = as.character(exploratories.inv$exp)
exploratories.inv$exp[which(exploratories.inv$exp %in% "acidophilous")] = "Acidophilous"
exploratories.inv$exp[which(exploratories.inv$exp %in% "alpine")] = "Alpine"
exploratories.inv$exp[which(exploratories.inv$exp %in% "beech")] = "Beech"
exploratories.inv$exp[which(exploratories.inv$exp %in% "boreal")] = "Boreal"
exploratories.inv$exp[which(exploratories.inv$exp %in% "broadleaved_evergreen")] = "Broadleaved evergreen"
exploratories.inv$exp[which(exploratories.inv$exp %in% "floodplain")] = "Floodplain"
exploratories.inv$exp[which(exploratories.inv$exp %in% "hemiboreal")] = "Hemiboreal"
exploratories.inv$exp[which(exploratories.inv$exp %in% "introduced")] = "Introduced"
exploratories.inv$exp[which(exploratories.inv$exp %in% "mediterranean_coniferous")] = "Mediterranean coniferous"
exploratories.inv$exp[which(exploratories.inv$exp %in% "mesophytic_deciduous")] = "Mesophytic deciduous"
exploratories.inv$exp[which(exploratories.inv$exp %in% "mountain_beech")] = "Mountain beech"
exploratories.inv$exp[which(exploratories.inv$exp %in% "non_riverine")] = "Non-riverine"
exploratories.inv$exp[which(exploratories.inv$exp %in% "swamp")] = "Swamp"
exploratories.inv$exp[which(exploratories.inv$exp %in% "thermophilous_deciduous")] = "Thermophilous deciduous"
exploratories.inv$exp = factor(exploratories.inv$exp, levels = c(
    "Mediterranean coniferous","Thermophilous deciduous","Broadleaved evergreen","Beech","Introduced","Non-riverine",
    "Mesophytic deciduous","Floodplain","Swamp","Acidophilous","Mountain beech","Hemiboreal","Boreal","Alpine"))

inventories.exper$platform = "inventories"
inventories.explo$platform = "inventories"

names(inventories.exper) = c("species","forest_type","SR_mixture","countries_compared","mono_mean","mono_sd","mono_n","mix_mean","mix_sd","mix_n","platform")
names(inventories.explo) = c("species","forest_type","SR_mixture","countries_compared","mono_mean","mono_sd","mono_n","mix_mean","mix_sd","mix_n","platform")

names(experiments.inv) = c("species","forest_type","mono_mean","mono_sd","mono_n","SR_mixture","mix_mean","mix_sd","mix_n","platform")
names(experiments.explo) = c("species","forest_type","mono_mean","mono_sd","mono_n","SR_mixture","mix_mean","mix_sd","mix_n","platform")

names(exploratories.inv) = c("species","forest_type","SR_mixture","mono_mean","mono_sd","mono_n","mix_mean","mix_sd","mix_n","platform")
names(exploratories.exper) = c("species","forest_type","SR_mixture","mono_mean","mono_sd","mono_n","mix_mean","mix_sd","mix_n","platform")

# add random effects
experiments.inv$countries_compared = paste(experiments.inv$platform, experiments.inv$forest_type, sep="_")
experiments.explo$countries_compared = paste(experiments.explo$platform, experiments.explo$forest_type, sep="_")
exploratories.inv$countries_compared = paste(exploratories.inv$platform, exploratories.inv$forest_type, sep="_")
exploratories.exper$countries_compared = paste(exploratories.exper$platform, exploratories.exper$forest_type, sep="_")

# no model with countries compared in the experiments and exploratories here
str(inventories.exper)
str(experiments.inv)
str(exploratories.inv)

# ckeck if species overlap -> works
sort(unique(inventories.exper$species))
sort(unique(experiments.inv$species))

sort(unique(inventories.explo$species))
sort(unique(exploratories.inv$species))

sort(unique(experiments.explo$species))
sort(unique(exploratories.exper$species))
#-> QURO is not measured in monoculture and in mixture partly only QURO and not the second species was measured
experiments.explo = experiments.explo[which(!(experiments.explo$species %in% c("QURO"))),]

# calculate log RR
inventories.exper$RR = log(inventories.exper$mix_mean / inventories.exper$mono_mean)
inventories.explo$RR = log(inventories.explo$mix_mean / inventories.explo$mono_mean)
experiments.inv$RR = log(experiments.inv$mix_mean / experiments.inv$mono_mean)
experiments.explo$RR = log(experiments.explo$mix_mean / experiments.explo$mono_mean)
exploratories.inv$RR = log(exploratories.inv$mix_mean / exploratories.inv$mono_mean)
exploratories.exper$RR = log(exploratories.exper$mix_mean / exploratories.exper$mono_mean)

#assign weights
inventories.exper.weights = data.table(inventories.exper)
inventories.exper.weights = inventories.exper.weights[,.(weight = length(RR)),by = .(species)]
inventories.exper = merge(inventories.exper, inventories.exper.weights, by.x = "species",by.y = "species")

inventories.explo.weights = data.table(inventories.explo)
inventories.explo.weights = inventories.explo.weights[,.(weight = length(RR)),by = .(species)]
inventories.explo = merge(inventories.explo, inventories.explo.weights, by.x = "species",by.y = "species")

exploratories.inv.weights = data.table(exploratories.inv)
exploratories.inv.weights = exploratories.inv.weights[,.(weight = length(RR)),by = .(species)]
exploratories.inv = merge(exploratories.inv, exploratories.inv.weights, by.x = "species",by.y = "species")

exploratories.exper.weights = data.table(exploratories.exper)
exploratories.exper.weights = exploratories.exper.weights[,.(weight = length(RR)),by = .(species)]
exploratories.exper = merge(exploratories.exper, exploratories.exper.weights, by.x = "species",by.y = "species")

experiments.inv.weights = data.table(experiments.inv)
experiments.inv.weights = experiments.inv.weights[,.(weight = length(RR)),by = .(species)]
experiments.inv = merge(experiments.inv, experiments.inv.weights, by.x = "species",by.y = "species")

experiments.explo.weights = data.table(experiments.explo)
experiments.explo.weights = experiments.explo.weights[,.(weight = length(RR)),by = .(species)]
experiments.explo = merge(experiments.explo, experiments.explo.weights, by.x = "species",by.y = "species")

########### plot and analysis

# inventories-experiments, compare species RR
lmer.inventories.exper = lmer(RR ~ species + (1|countries_compared) - 1, data = inventories.exper, weights = weight)
lmer.experiments.inv = lmer(RR ~ species + (1|countries_compared) - 1, data = experiments.inv, weights = weight)
coef.inv.exper = summary(lmer.inventories.exper)$coefficients[,"Estimate"]
coef.exper.inv = summary(lmer.experiments.inv)$coefficients[,"Estimate"]

cor.test(coef.inv.exper,coef.exper.inv, method = "kendall")

plot_inv_exper = ggplot()
svg(filename = "manuscript 1_meta-analysis\\output\\restricted_comparison_inv_exper.svg",height=5, width = 5)
plot_inv_exper + geom_point(aes(x = coef.inv.exper, y = coef.exper.inv), size = 2) +
    geom_text(aes(x = coef.inv.exper, y = coef.exper.inv-0.02, family = "Times New Roman",
                  label = gsub(patter = "species", replacement = "", x = names(coef.inv.exper))), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"),
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Mixture effect (log RR) - inventory approach") +
    ylab("Mixture effect (log RR) - experimental approach") + ggtitle("a)")
graphics.off()

# inventories-experiments, compare species RR
lmer.inventories.explo = lmer(RR ~ species + (1|countries_compared) - 1, data = inventories.explo, weights = weight)
lmer.exploratories.inv = lmer(RR ~ species + (1|countries_compared) - 1, data = exploratories.inv, weights = weight)
coef.inv.explo = summary(lmer.inventories.explo)$coefficients[,"Estimate"]
coef.explo.inv = summary(lmer.exploratories.inv)$coefficients[,"Estimate"]

cor.test(coef.inv.explo,coef.explo.inv, method = "kendall")

plot_inv_explo = ggplot()
svg(filename = "manuscript 1_meta-analysis\\output\\restricted_comparison_inv_explo.svg",height=5, width = 5)
plot_inv_explo + geom_point(aes(x = coef.inv.explo, y = coef.explo.inv), size = 2) +
    geom_text(aes(x = coef.inv.explo, y = coef.explo.inv-0.05, family = "Times New Roman", 
                  label = gsub(patter = "species", replacement = "", x = names(coef.explo.inv))), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"),
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Mixture effect (log RR) - inventory approach") +
    ylab("Mixture effect (log RR) - exploratory approach") +ggtitle("b)")
graphics.off()

# experiments-exploratories, compare species RR
lmer.experiments.explo = lmer(RR ~ species + (1|countries_compared) - 1, data = experiments.explo, weights = weight)
lmer.exploratories.exper = lmer(RR ~ species + (1|countries_compared) - 1, data = exploratories.exper, weights= weight)
coef.exper.explo = summary(lmer.experiments.explo)$coefficients[,"Estimate"]
coef.explo.exper = summary(lmer.exploratories.exper)$coefficients[,"Estimate"]

cor.test(coef.exper.explo,coef.explo.exper, method = "kendall")

plot_exper_explo = ggplot()
svg(filename = "manuscript 1_meta-analysis\\output\\restricted_comparison_exper_explo.svg",height=5, width = 5)
plot_exper_explo + geom_point(aes(x = coef.exper.explo, y = coef.explo.exper), size = 2) +
    geom_text(aes(x = coef.exper.explo, y = coef.explo.exper-0.03, family = "Times New Roman",
                  label = gsub(patter = "species", replacement = "", x = names(coef.exper.explo))), size = 4) +
    geom_hline(yintercept = 0, linetype = "dotted") +
    geom_vline(xintercept = 0, linetype = "dotted") +
    theme_bw() + theme(text = element_text(family = "Times New Roman"),
                       panel.grid.minor = element_blank(),panel.grid.major = element_blank()) +
    xlab("Mixture effect (log RR) - exploratory approach") +
    ylab("Mixture effect (log RR) - experimental approach") + ggtitle("c)")
graphics.off()


###############
# inventories-exploratories but now with the overlapping species:exp combinations
# correlate model coefficient estimates between inv and explo / separated by fores type

# get mean effect sizes per species and forest type
inventories.explo.forest.types = unique(inventories.explo[,c("species","forest_type")])
exploratories.inv.forest.types = unique(exploratories.inv[,c("species","forest_type")])
inventories.explo.forest.types$lmer_intercept_inv = NA
exploratories.inv.forest.types$lmer_intercept_explo = NA

# run models for all subsets of species and forest types
# inventories
for(i in 1:nrow(inventories.explo.forest.types)){
    data.subset = subset(inventories.explo, species %in% as.character(inventories.explo.forest.types$species)[i] &
                                            forest_type %in% as.character(inventories.explo.forest.types$forest_type)[i])
    if(length(unique(data.subset$countries_compared)) >=2 & length(unique(data.subset$countries_compared)) < nrow(data.subset)){
        lmer.temp = lmer(RR ~ (1|countries_compared), data = data.subset, weights = weight)
        inventories.explo.forest.types$lmer_intercept_inv[i] = summary(lmer.temp)$coefficients[,"Estimate"]
    }else{
        lmer.temp = lm(RR ~ 1, data = data.subset, weights = weight)
        inventories.explo.forest.types$lmer_intercept_inv[i] = summary(lmer.temp)$coefficients[,"Estimate"]}}

# exploratories
for(i in 1:nrow(exploratories.inv.forest.types)){
    data.subset = subset(exploratories.inv, species %in% as.character(exploratories.inv.forest.types$species)[i] &
                                            forest_type %in% as.character(exploratories.inv.forest.types$forest_type)[i])
    if(length(unique(data.subset$countries_compared)) >=2){
        lmer.temp = lmer(RR ~ (1|countries_compared), data = data.subset, weights = weight)
        exploratories.inv.forest.types$lmer_intercept_explo[i] = summary(lmer.temp)$coefficients[,"Estimate"]
    }else{
        lmer.temp = lm(RR ~ 1, data = data.subset, weights = weight)
        exploratories.inv.forest.types$lmer_intercept_explo[i] = summary(lmer.temp)$coefficients[,"Estimate"]}}

# test correlation in coefficients
restricted.inv.explo.species.forest.types = merge(inventories.explo.forest.types, exploratories.inv.forest.types,
                                       by.x = c("species","forest_type"), by.y = c("species","forest_type"))
cor.test(restricted.inv.explo.species.forest.types$lmer_intercept_inv,restricted.inv.explo.species.forest.types$lmer_intercept_explo,
         method="kendall")

plot_restricted_inv_explo_forest_Type = ggplot(data=restricted.inv.explo.species.forest.types)
svg(filename = "manuscript 1_meta-analysis\\output\\restricted_comparison_inv_explo_per_forest_type.svg",height=5, width = 5)
plot_restricted_inv_explo_forest_Type + geom_point(aes(x = lmer_intercept_inv, y = lmer_intercept_explo, color = forest_type), size = 2) +
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





















plot(coef.inv.explo.overlap,coef.explo.inv.overlap, col = as.factor(exp.vector), pch=16); abline(h=0); abline(v=0)


# create datasets with same species_exp combinations for comparison of logRR
#-> not for comparison with experiments as these are no forest types
RR_inv_exper = rbind(inventories.exper,experiments.inv)
RR_inv_explo = rbind(inventories.exper,experiments.inv)

RR_inv_exper_replicates = RR_inv_exper[duplicated(RR_inv_exp[,c("species","exp")]),c("species","exp")]




# analysis with exp as random effect
lmer.inventories.exper = lmer(RR ~ species + (1|exp) - 1, data = inventories.exper)
lmer.inventories.explo = lmer(RR ~ species + (1|exp) - 1, data = inventories.explo)
lmer.experiments.inv = lmer(RR ~ species + (1|exp) - 1, data = experiments.inv)
lmer.experiments.explo = lmer(RR ~ species + (1|exp) - 1, data = experiments.explo)
lmer.exploratories.inv = lmer(RR ~ species + (1|exp) - 1 , data = exploratories.inv)
lmer.exploratories.exper = lmer(RR ~ species + (1|exp) - 1, data = exploratories.exper)

# extract coefficient estimates ################################
coef.inv.exper = summary(lmer.inventories.exper)$coefficients[,"Estimate"]
coef.exper.inv = summary(lmer.experiments.inv)$coefficients[,"Estimate"]
coef.inv.explo = summary(lmer.inventories.explo)$coefficients[,"Estimate"]
coef.explo.inv = summary(lmer.exploratories.inv)$coefficients[,"Estimate"]
coef.exper.explo = summary(lmer.exploratories.exper)$coefficients[,"Estimate"]
coef.explo.exper = summary(lmer.experiments.explo)$coefficients[,"Estimate"]

# check if species overlap
names(coef.inv.exper) == names(coef.exper.inv)
names(coef.inv.explo) == names(coef.explo.inv)
names(coef.exper.explo) == names(coef.explo.exper)

# correlations
cor.test(coef.inv.exper,coef.exper.inv, method="spearman")
cor.test(coef.inv.explo,coef.explo.inv, method="spearman")
cor.test(coef.explo.exper,coef.exper.explo, method="spearman")

# plot
plot(coef.inv.exper,coef.exper.inv); abline(h=0); abline(v=0)
plot(coef.inv.explo,coef.explo.inv); abline(h=0); abline(v=0)
plot(coef.explo.exper,coef.exper.explo); abline(h=0); abline(v=0)

# analysis as exp:species interaction between inventories and exploratories ##########################
lm.inventories.explo = lm(RR ~ species:exp - 1, data = inventories.explo)
lm.exploratories.inv = lm(RR ~ species:exp - 1, data = exploratories.inv)

# extract coefficient estimates
coef.inv.explo.overlap = summary(lm.inventories.explo)$coefficients[,"Estimate"]
coef.explo.inv.overlap = summary(lm.exploratories.inv)$coefficients[,"Estimate"]

# only include coefficients that are in both models
coef.inv.explo.overlap = coef.inv.explo.overlap[names(coef.inv.explo.overlap) %in% names(coef.explo.inv.overlap)]
coef.explo.inv.overlap = coef.explo.inv.overlap[names(coef.explo.inv.overlap) %in% names(coef.inv.explo.overlap)]
coef.inv.explo.overlap = coef.inv.explo.overlap[order(names(coef.inv.overlap))]
coef.explo.inv.overlap = coef.explo.inv.overlap[order(names(coef.explo.inv.overlap))]

length(coef.inv.explo.overlap)
length(coef.explo.inv.overlap)

# make a vector for the experiment
exp.vector = gsub(pattern=".*:exp", replacement="",x = names(coef.explo.inv.overlap))

# correlations
cor.test(coef.inv.explo.overlap,coef.explo.inv.overlap, method="kendall")

# plot
plot(coef.inv.explo.overlap,coef.explo.inv.overlap, col = as.factor(exp.vector), pch=16); abline(h=0); abline(v=0)

############### Plotting ######################
library(ggplot2)

plot1 = ggplot()
plot1 + geom_point(aes(x = coef.explo.inv.overlap, y = coef.inv.explo.overlap, 
                       col = exp.vector))
