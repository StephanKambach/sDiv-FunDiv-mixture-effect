# ### FunDiv Synthesis:
# ### Diversity ~ Growth across platforms
# ### By: Stephan Kambach
# ### Date: 29.02.2016
# 


#################################################
# use update from Simon -----------------------------

# libraries
if(!require(data.table)){
    install.packages("data.table")
    library(data.table)}
if(!require(ggplot2)){
    install.packages("ggplot2")
    library(ggplot2)}
if(!require(cowplot)){
    install.packages("cowplot")
    library(cowplot)}
if(!require(raster)){
    install.packages("raster")
    library(raster)}


rm(list=ls())
setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Experiments")
#setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\data\\Experiments")

correct.names = read.table("species_names_correction.txt", header = T, sep="\t")
names(correct.names)[1] = "code_old"
correct.names$code_old = as.character(correct.names$code_old)
correct.names$code_new = as.character(correct.names$code_new)

all.species.compositions = c()

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
experiment1 = data.table(read.table("update_from_Simon\\forbio.zed.plot_update03.05.2017.csv",header=T,sep=","))
experiment1$SR[which(experiment1$SR > 3)] = 4
str(experiment1)

# rename
experiment1$sp.comp = rename.species.compositions(as.character(experiment1$sp.comp))
experiment1$sp.code = rename.species(as.character(experiment1$sp.code))

# all.species.compositions
all.species.compositions = c(all.species.compositions,as.character(experiment1$sp.code),as.character(experiment1$sp.comp))

experiment1.summarized = experiment1[,.("exp"="forbio.zed",
                                        "growth_mean"=mean(basaldiam.mm,na.rm=T),"growth_sd"=sd(basaldiam.mm,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment1.summarized)[1] = "species"                                       

experiment1.ordered= re.order.function(data.to.order = experiment1.summarized)

####### experiments 2 ######## growth as basal area.cm2 after 8 years, also height available
experiment2 = data.table(read.table("update_from_Simon\\kreinitz.plot.csv",header=T,sep=","))
experiment2$SR[which(experiment2$SR > 3)] = 4
str(experiment2)

# rename
experiment2$sp.comp = rename.species.compositions(as.character(experiment2$sp.comp))
experiment2$sp.code = rename.species(as.character(experiment2$sp.code))

# all.species.compositions
all.species.compositions = c(all.species.compositions,as.character(experiment2$sp.code),as.character(experiment2$sp.comp))

experiment2.summarized = experiment2[,.("exp"="kreinitz",
                                        "growth_mean"=mean(basalarea.cm2,na.rm=T),"growth_sd"=sd(basalarea.cm2,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment2.summarized)[1] = "species"                                       

experiment2.ordered= re.order.function(data.to.order = experiment2.summarized)

####### experiments 3 ######## growth as dbh.mmm after 12 years, also basal area.cm2 and height available
experiment3 = data.table(read.table("update_from_Simon\\satakunta.plot.csv",header=T,sep=","))
experiment3$SR[which(experiment3$SR > 3)] = 4
str(experiment3)

# rename
experiment3$sp.comp = rename.species.compositions(as.character(experiment3$sp.comp))
experiment3$sp.code = rename.species(as.character(experiment3$sp.code))

# all.species.compositions
all.species.compositions = c(all.species.compositions,as.character(experiment3$sp.code),as.character(experiment3$sp.comp))

experiment3.summarized = experiment3[,.("exp"="satakunta",
                                        "growth_mean"=mean(dbh.mm,na.rm=T),"growth_sd"=sd(dbh.mm,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment3.summarized)[1] = "species"                                       

experiment3.ordered= re.order.function(data.to.order = experiment3.summarized)

####### experiments 4 ######## growth as basaldiameter.mean after 5 years, also height.mean available
experiment4 = data.table(read.table("update_from_Simon\\ged.plot.csv",header=T,sep=","))
experiment4$SR[which(experiment4$SR > 3)] = 4
str(experiment4)

# rename
experiment4$sp.comp = rename.species.compositions(as.character(experiment4$sp.comp))
experiment4$sp.code = rename.species(as.character(experiment4$sp.code))

# all.species.compositions
all.species.compositions = c(all.species.compositions,as.character(experiment4$sp.code),as.character(experiment4$sp.comp))

experiment4.summarized = experiment4[,.("exp"="forbio.ged",
                                        "growth_mean"=mean(basaldiameter.mean,na.rm=T),"growth_sd"=sd(basaldiameter.mean,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment4.summarized)[1] = "species"                                       

experiment4.ordered= re.order.function(data.to.order = experiment4.summarized)
str(experiment4.ordered)

####### experiments 5 ######## growth as basalarea.mean after 11 years
experiment5 = data.table(read.table("update_from_Simon\\biotree\\biotree.subplot.csv",header=T,sep=","))
experiment5$SR[which(experiment5$SR > 3)] = 4
str(experiment5)

# rename
experiment5$sp.comp = rename.species.compositions(as.character(experiment5$sp.comp))
experiment5$sp.code = rename.species(as.character(experiment5$sp.code))

# all.species.compositions
all.species.compositions = c(all.species.compositions,as.character(experiment5$sp.code),as.character(experiment5$sp.comp))

experiment5.summarized = experiment5[,.("exp"="biotree",
                                        "growth_mean"=mean(basalarea.cm2,na.rm=T),"growth_sd"=sd(basalarea.cm2,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]
names(experiment5.summarized)[1] = "species"                                       

experiment5.ordered= re.order.function(data.to.order = experiment5.summarized)
str(experiment5.ordered)

####### experiments 6 ######## growth as height in 2014 after 6 years
experiment6 = data.table(read.table("update_from_Simon\\gro.sp.orphee_new_2015.csv",sep=",",header=T,row.names=1))
experiment6$richness[which(experiment6$richness > 3)] = 4
experiment6$X = seq(1:nrow(experiment6))
str(experiment6)
names(experiment6)[which(names(experiment6) == "treatment")] = "sp.comp"
names(experiment6)[which(names(experiment6) == "SP")] = "sp.code"
names(experiment6)[which(names(experiment6) == "richness")] = "SR"
experiment6$sp.comp =  gsub(pattern="\\+",replacement="_",as.character(experiment6$sp.comp))

# rename
experiment6$sp.comp = rename.species.compositions(as.character(experiment6$sp.comp))
experiment6$sp.code = rename.species(as.character(experiment6$sp.code))

# all.species.compositions
all.species.compositions = c(all.species.compositions,as.character(experiment6$sp.code),as.character(experiment6$sp.comp))

experiment6.summarized = experiment6[,.("exp"="orphee",
                                        "growth_mean"=mean(mean_plot,na.rm=T),"growth_sd"=sd(mean_plot,na.rm=T),"growth_n"=length(unique(X))),
                                     by=.(sp.code,SR)]

names(experiment6.summarized) = c("species","SR","exp","growth_mean","growth_sd","growth_n")

experiment6.ordered= re.order.function(data.to.order = experiment6.summarized)

################## combine all experiments ###########################
all.experiments.ordered = rbind(experiment1.ordered,experiment2.ordered,experiment3.ordered,
                                experiment4.ordered,experiment5.ordered,experiment6.ordered)
all.experiments.ordered$platform = "experiments"

write.table(all.experiments.ordered, "all_experiments_re-arranged.csv",sep="\t",dec=".",quote=F,row.names=F)

all.species.compositions = unique(all.species.compositions)

write.table(all.species.compositions, "all_species_compositions.csv",sep="\t",dec=".",quote=F,row.names=F)

############
#figures----
setwd("C:\\Users\\kambach\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\")
setwd("C:\\Users\\Agando\\Dropbox\\current_tasks\\FunDiv WS-GRoup1\\R_scripts\\")

locations = read.csv("intermediate output\\locations of experiments and exploratories.csv",sep=";",dec=",")
str(locations)
locations = subset(locations, Type %in% "Experiment")
clim = getData("worldclim",var="bio",res=2.5 )
bioclim.data = extract(clim,locations[,c("Longitude","Latitude")])
locations = cbind(locations,bioclim.data)

locations$bio1 = locations$bio1 / 10
locations$bio7 = locations$bio7 / 10

png("data\\Experiments\\intermediate_output\\data_overview\\bio1.png",height=600,width=800)
ggplot()+
    geom_histogram(aes(x=locations$bio1)) +
    theme_cowplot() +  xlab("\nMean annual temperature °C") + ylab("Count\n") +
    theme(text = element_text(size=20),axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) +
    xlim(c(-3,20))
graphics.off()

png("data\\Experiments\\intermediate_output\\data_overview\\bio12.png",height=600,width=800)
ggplot()+
    geom_histogram(aes(x=locations$bio12)) +
    theme_cowplot() +  xlab("\nAnnual precipitation mm") + ylab("Count\n") +
    theme(text = element_text(size=20),axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) +
    xlim(c(200,1700))
graphics.off()

png("data\\Experiments\\intermediate_output\\data_overview\\bio7.png",height=600,width=800)
ggplot()+
    geom_histogram(aes(x=locations$bio7)) +
    theme_cowplot() +  xlab("\nTemperature annual range °C") + ylab("Count\n") +
    theme(text = element_text(size=20),axis.text.x = element_text(size=20),axis.text.y = element_text(size=20)) +
    xlim(c(14,40))
graphics.off()
