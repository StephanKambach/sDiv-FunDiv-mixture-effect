### FunDiv Synthesis:
### Diversity ~ Growth across platforms
### By: Stephan Kambach
### Date: 29.02.2016

#!there was some error data in the spain growth data and I deleted it
#spain = read.csv("C:\\Users\\kambach\\Desktop\\aktuelle Arbeiten\\FunDiv WS-GRoup1\\R_scripts\\Exploratories\\Raw data\\Growth_data\\Spain_growth_data.csv")
#spain = spain[c(1:488),]
#spain = spain[,-1]
#write.csv(spain,"C:\\Users\\kambach\\Desktop\\aktuelle Arbeiten\\FunDiv WS-GRoup1\\R_scripts\\Exploratories\\Raw data\\Growth_data\\Spain_growth_data.csv",row.names=F)

rm(list=ls())
# libraries
if(!require(data.table)){
    install.packages("data.table")
    library(data.table)}
if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)}

#data
setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Exploratories")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Exploratories")

files.plots = list.files("Plot_descriptors\\")
files.growth = list.files("Growth_data\\")
plot.alldata = read.table("alldata_sDiv_FunDiv_exploratories20150611.txt",sep="\t",header=T)
str(plot.alldata)
#dataframe to compile all effect sizes in one
growth.es.final = data.frame("SPCode"=NA,"exp"=NA,"mono_mean"=NA,"mono_sd"=NA,"mono_n"=NA,
                             "SR"=NA,"mix_mean"=NA,"mix_sd"=NA,"mix_n"=NA)[0,]

rename_species = read.table("species_names_correction.txt",header=T, sep="\t")
#str(rename_species)
names(rename_species)[1]  = "code_old"

# a table to extract the growth per species and plot for the supplement
growth.table.for.supplement = data.frame(Increase_species_ba = NA, PlotID = NA, SpeciesCode= NA)[0,]

# collect all species compositions
all.species.combinations = c()
#this loop takes the growth and plot description for each site 
#and calculates effect sizes for each species in each mixture
for(i in 1:6){
    plot.raw = read.csv(paste("Plot_descriptors\\",files.plots[i],sep=""))
    growth.raw = read.csv(paste("Growth_data\\",files.growth[i],sep=""))

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

    # get all species combinations
    for(k in unique(growth.mean$PlotID)){
        growth.mean.temp = subset(growth.mean, PlotID %in% k)
        species.temp = sort(unique(growth.mean.temp$SpeciesCode))
        all.species.combinations = c(all.species.combinations,paste(species.temp, collapse="_"))
    }
    #calculate the mean growth rate per species over all trees per plot
    growth.mean = data.table(growth.mean)
    
    growth.summary = growth.mean[,.(Increase_species_ba = mean(Increase_ba)),by=.(PlotID,SpeciesCode)]
    growth.table.for.supplement = rbind(growth.table.for.supplement,growth.summary)
    
    growth.summary = merge(growth.summary,plot.alldata[,c("PlotID","Restricted_Species_Richness")],
                    by.x="PlotID",by.y="PlotID",
                    all.x=TRUE,all.y=FALSE)
    
    #rename SR and cut to 4
    growth.summary$SR = growth.summary$Restricted_Species_Richness
    growth.summary$SR[which(growth.summary$SR > 3)] = 4
    
    #get the location of the site
    location = sub("_growth_data.csv*","",files.growth[i])

    #calculate growth per Species and Richness level
    growth.summary = growth.summary[,.(exp=location, 
                             SPECIES_MEAN = mean(Increase_species_ba),
                             SPECIES_SD = sd(Increase_species_ba),
                             SPECIES_N = length(unique(Increase_species_ba))),
                          by=.(SpeciesCode,SR)]
                            
    #order table so that each row contains monoculture and a mixture
    growth.es = merge(growth.summary[growth.summary$SR %in% 1,],growth.summary[!(growth.summary$SR %in% 1),],
             by.x="SpeciesCode",by.y="SpeciesCode",all.y=T)

    
    growth.es[,SR.x := NULL]
    growth.es[,exp.y := NULL]
    
    names(growth.es) = c("SPCode","exp","mono_mean","mono_sd","mono_n","SR","mix_mean","mix_sd","mix_n")

    #put all together in result dataframe with all sites, species, mixtures
    growth.es.final = rbind(growth.es.final,growth.es)

}

exploratories = data.frame("species"=growth.es.final$SPCode,"exp"=growth.es.final$exp,"SR"=growth.es.final$SR,
                           "mono_mean"=growth.es.final$mono_mean,"mono_sd"=growth.es.final$mono_sd,"mono_n"=growth.es.final$mono_n,
                           "mix_mean"=growth.es.final$mix_mean,"mix_sd"=growth.es.final$mix_sd,"mix_n"=growth.es.final$mix_n,
                           "platform"= "exploratories")

exploratories = exploratories[- which(is.na(exploratories$mono_mean)),]

#rename forest types
exploratories$exp = as.character(exploratories$exp)
exploratories$exp[which(exploratories$exp %in% "Finland")] = "boreal"
exploratories$exp[which(exploratories$exp %in% "Poland")] = "hemiboreal"
exploratories$exp[which(exploratories$exp %in% "Germany")] = "beech"
exploratories$exp[which(exploratories$exp %in% "Romania")] = "mountain_beech"
exploratories$exp[which(exploratories$exp %in% "Italy")] = "thermophilous_deciduous"
exploratories$exp[which(exploratories$exp %in% "Spain")] = "mediterranean_coniferous"

table(exploratories$exp)

write.table(exploratories,"all_exploratories.csv",
            quote=F,sep="\t",row.names=F,dec=".")

write.table(growth.table.for.supplement, "intermediate_output\\growth_per_species_and_plot_for_supplement.csv",sep="\t",quote=F)

all.species.combinations = as.data.frame(table(all.species.combinations))
all.species.combinations = all.species.combinations$all.species.combinations[which(all.species.combinations$Freq >= 1)]
write.table(all.species.combinations, "all_species_compositions.csv",sep="\t",dec=".",quote=F,row.names=F)


########### count species pool 
species.nr = c()
species = c()
for(i in 1:6){
    plot.raw = read.csv(paste("Plot_descriptors\\",files.plots[i],sep=""))
    species.nr = c(species.nr,length(unique(plot.raw$Species)))
    species = c(species,as.character(plot.raw$Species))
    }
length(unique(species))
