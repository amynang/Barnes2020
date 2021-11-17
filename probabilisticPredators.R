rm(list=ls())
library(tidyverse)
library(reshape2)
library(fluxweb)
library(vegan) 
library(ggplot2)
library(patchwork)
library(parallel)
library(rlist)
library(brms)
library(glmmTMB)
library(DHARMa)
n.cores<-detectCores()

# https://figshare.com/articles/software/R_Code_for_Barnes_et_al_Biodiversity_enhances_the_multi-trophic_control_of_arthropod_herbivory_/12909962/1
# The code in the link above is efficient but (to me) it is difficult to see 
# what is going on
# So this is a "rephrasing" that works for me. 
# Probably just a matter of personal taste.

####################### the data ######################## 
# Jena
download.file("https://figshare.com/ndownloader/files/24558608", 
              destfile = "JenaExp_arthro.csv")
je.com.data = read.csv("JenaExp_arthro.csv")

download.file("https://figshare.com/ndownloader/files/24558617", 
              destfile = "JenaExp_matrix.csv")
je.foodwebmat = as.matrix(read.csv("JenaExp_matrix.csv", row.names = 1))

# Cedar Creek
download.file("https://figshare.com/ndownloader/files/24558605", 
              destfile = "CedCrExp_arthro.csv")
cc.com.data = read.csv("CedCrExp_arthro.csv")

download.file("https://figshare.com/ndownloader/files/24558614", 
              destfile = "CedCrExp_matrix.csv")
cc.foodwebmat = as.matrix(read.csv("CedCrExp_matrix.csv", row.names = 1))

# individual attribute dataframes for each year-plot combination in Jena
je.att = split(je.com.data, with(je.com.data, sample))
# individual matrices for each year-plot combination in Jena
je.web = vector(mode = "list", length=length(je.att))
for (i in 1:length(je.att)) {
  je.web[[i]] = je.foodwebmat[je.att[[i]]$trophic.group,
                              je.att[[i]]$trophic.group]
}

# individual attribute dataframes for each year-plot combination in Cedar Creek
cc.att = split(cc.com.data, with(cc.com.data, sample))
# individual matrices for each year-plot combination in Cedar Creek
cc.web = vector(mode = "list", length=length(cc.att))
for (i in 1:length(cc.att)) {
  cc.web[[i]] = cc.foodwebmat[cc.att[[i]]$trophic.group,
                              cc.att[[i]]$trophic.group]
}


thousand = vector(mode = "list",length=1000)

for (k in 1:1000) { 

# a list for feeding preference matrices
mat.prefs = vector(mode = "list",length=length(je.att))
# a list for flux matrices
fluxes = vector(mode = "list",length=length(je.att))

allmetrics1 = data.frame(experiment = factor("Jena"),
                         year = integer(length(je.att)),
                         plot = character(length(je.att)),
                         plant.rich = integer(length(je.att)), # richness gradient
                         tot.flux = numeric(length(je.att)),   # total energy flux
                         pred.flux = numeric(length(je.att)),  # predation flux
                         herb.flux = numeric(length(je.att)),  # herbivory flux
                         top.down = numeric(length(je.att)),   # from herbivores per unit herbivore biomass
                         bot.up = numeric(length(je.att)),     # to herbivores per unit herbivore biomass
                         herb.press = numeric(length(je.att))) # to herbivores per unit plant biomass   

# This loop works through all the year-plot foodwebs of the respective experiment.
# In each one first it creates a feeding preference matrix such that omnivores are 
# feeding equally from their resource channels (animals, plants, detritus) and within 
# each channel they eat from individual nodes based on the relative availability 
# of each node. Predators are eating from animal nodes also based on relative 
# availability (biomass).
# Then it calculates fluxes for each foodweb and places the resulting flux matrix
# in the respective slot in the list "fluxes".
# Then it calculates flux aggregates of interest like fluxes from herbivores and
# standardizes them by the appropriate biomass. These end up in the "allmetrics"
# dataframe with one row per foodweb (plot-year combination)

for (i in 1:length(je.att)) {
  
  ####################   Omnivores' Balanced Diet Plan   #########################
  # OK, so this bit is not easy to read but it works(?)
  # also, it borrows heavily from hand-me-down code
  # that was probably originally written by Barnes
  # (Not sure where this happens in the code of the
  # paper, but it is mentioned in the methods)
  
  # add biomass values in the matrix to 'manually' define the preferences
  # first create a matrix with species biomasses
  mat.bioms = replicate(length(je.att[[i]]$biomass), je.att[[i]]$biomass)
  # mat.prefs contains preference of predators based on their prey biomasses
  mat.prefs[[i]] = je.web[[i]] * mat.bioms
  # then identify omnivorous species
  # these are species feeding at least on:
  # one basal species and at least one non basal species
  # we want them to feed equally on plant animals and detritus
  basals = which(je.att[[i]]$trophic.level == "basal")
  animals = which(je.att[[i]]$trophic.level != "basal") ; #S_animals[i]=length(animals)
  plants = which(je.att[[i]]$trophic.group == "plants")        #; S_plants[i]=length(plants)
  detritus = which(je.att[[i]]$trophic.group == "detritus")
  #not relevant here but we will need to keep track of herbivores and predators
  herbivores = which(je.att[[i]]$trophic.level == "herbivore")
  predators = which(je.att[[i]]$trophic.level == "carnivore")
  
  #omnivores that feed on detritus plants and animals
  omnivores.c.h.d = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])>0 &
                            colSums(mat.prefs[[i]][plants,,drop = FALSE])>0 &
                            colSums(mat.prefs[[i]][animals,])>0)
  
  #omnivores that feed on plants and animals
  omnivores.c.h = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])==0 &
                          colSums(mat.prefs[[i]][plants,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][animals,])>0)
  
  #omnivores that feed on detritus and animals
  omnivores.c.d = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][plants,,drop = FALSE])==0 &
                          colSums(mat.prefs[[i]][animals,])>0)
  
  #omnivores that feed on detritus and plants
  omnivores.h.d = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][plants,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][animals,])==0)
  
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][animals, omnivores.c.h.d] = mat.prefs[[i]][animals, omnivores.c.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores.c.h.d,drop=FALSE])),
         length(omnivores.c.h.d),length(omnivores.c.h.d)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][plants, omnivores.c.h.d] = mat.prefs[[i]][plants, omnivores.c.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][plants, omnivores.c.h.d,drop=FALSE])),
         length(omnivores.c.h.d),length(omnivores.c.h.d))
  
  # normalize preferences of omnivores over detritus to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][detritus, omnivores.c.h.d] = mat.prefs[[i]][detritus, omnivores.c.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][detritus, omnivores.c.h.d,drop=FALSE])),
         length(omnivores.c.h.d),length(omnivores.c.h.d))
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][animals, omnivores.c.h] = mat.prefs[[i]][animals, omnivores.c.h,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores.c.h,drop=FALSE])),
         length(omnivores.c.h),length(omnivores.c.h)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][plants, omnivores.c.h] = mat.prefs[[i]][plants, omnivores.c.h,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][plants, omnivores.c.h,drop=FALSE])),
         length(omnivores.c.h),length(omnivores.c.h))
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][animals, omnivores.c.d] = mat.prefs[[i]][animals, omnivores.c.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores.c.d,drop=FALSE])),
         length(omnivores.c.d),length(omnivores.c.d)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][detritus, omnivores.c.d] = mat.prefs[[i]][detritus, omnivores.c.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][detritus, omnivores.c.d,drop=FALSE])),
         length(omnivores.c.d),length(omnivores.c.d)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][plants, omnivores.h.d] = mat.prefs[[i]][plants, omnivores.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][plants, omnivores.h.d,drop=FALSE])),
         length(omnivores.h.d),length(omnivores.h.d))
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][detritus, omnivores.h.d] = mat.prefs[[i]][detritus, omnivores.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][detritus, omnivores.h.d,drop=FALSE])),
         length(omnivores.h.d),length(omnivores.h.d))
  
  mat.prefs[[i]] = decostand(mat.prefs[[i]], "total", MARGIN = 2)
  
  # then make total biomass of animal prey species equal to the biomass of detritus:
  #mat.prefs[[i]][!basals, omnivores] =  mat.prefs[[i]][!basals, omnivores] * biomass.plot[1]
  ################################################################################
  
  ################################# Uncertainty ##################################
  # Here we take each consumer in the foodweb and replace its fixed preferences 
  # with a random sample from a dirichlet distribution whose component probabilities 
  # are given by the vector of the original preferences. The vector is multiplied 
  # by a scalar that modifies the shape of the distribution (larger = less uncertainty)
  # Across several iterations our expectations regarding what consumers feed on 
  # are met, on average. But in each iteration consumer preferences deviate somewhat
  # from those expected based on intrinsic preference and/or relative availability.
  for (j in 1: nrow(mat.prefs[[i]])) { 
    mat.prefs[[i]][,j] = rdirichlet(1, mat.prefs[[i]][,j]*100)
  }
  mat.prefs[[i]][is.nan(mat.prefs[[i]])] = 0 #removes NaNs from basal node "preferences"
  
  ################################################################################
  
  fluxes[[i]] <- fluxing(mat.prefs[[i]],
                         je.att[[i]]$biomass, 
                         je.att[[i]]$com.metabolism,
                         je.att[[i]]$efficiency,
                         bioms.prefs = F,
                         bioms.losses = F,
                         ef.level = "prey")
  
  fluxes[[i]] = fluxes[[i]] * 744 #Joules per month
  
  
  allmetrics1[i,]$year = unique(je.att[[i]]$year)
  allmetrics1[i,]$plot = unique(je.att[[i]]$plot)
  allmetrics1[i,]$plant.rich = unique(je.att[[i]]$plant.div)          # richness gradient
  allmetrics1[i,]$tot.flux = sum(fluxes[[i]])                         # total energy flux              
  allmetrics1[i,]$pred.flux = sum(fluxes[[i]][animals, ])             # predation flux
  allmetrics1[i,]$herb.flux = sum(fluxes[[i]][plants, ])              # herbivory flux
  allmetrics1[i,]$top.down = sum(fluxes[[i]][herbivores, predators])/ # from herbivores per unit herbivore biomass
    sum(je.att[[i]][herbivores,"biomass"])
  allmetrics1[i,]$bot.up = sum(fluxes[[i]][plants, herbivores])/      # to herbivores per unit herbivore biomass
    sum(je.att[[i]][herbivores,"biomass"])
  allmetrics1[i,]$herb.press = sum(fluxes[[i]][plants, herbivores])/  # to herbivores per unit plant biomass 
    sum(unique(je.att[[i]]$plant.biomass))     
  
}




######################### Same deal for Cedar Creek ############################

# a list for feeding preference matrices
mat.prefs = vector(mode = "list",length=length(cc.att))
# a list for flux matrices
fluxes = vector(mode = "list",length=length(cc.att))

allmetrics2 = data.frame(experiment = factor("Cedar"),
                         year = integer(length(cc.att)),
                         plot = character(length(cc.att)),
                         plant.rich = integer(length(cc.att)), # richness gradient
                         tot.flux = numeric(length(cc.att)),   # total energy flux
                         pred.flux = numeric(length(cc.att)),  # predation flux
                         herb.flux = numeric(length(cc.att)),  # herbivory flux
                         top.down = numeric(length(cc.att)),   # from herbivores per unit herbivore biomass
                         bot.up = numeric(length(cc.att)),     # to herbivores per unit herbivore biomass
                         herb.press = numeric(length(cc.att))) # to herbivores per unit plant biomass   

# This loop works through all the year-plot foodwebs of the respective experiment.
# In each one, first it creates a feeding preference matrix such that omnivores are 
# feeding equally from their resource channels (animals, plants, detritus) and within 
# each channel they eat from individual nodes based on the relative availability 
# of each node. Predators are eating from animal nodes also based on relative 
# availability (biomass).
# Then it calculates fluxes for each foodweb and places the resulting flux matrix
# in the respective slot in the list "fluxes".
# Then it calculates flux aggregates of interest like fluxes from herbivores and
# standardizes them by the appropriate biomass. These end up in the "allmetrics"
# dataframe with one row per foodweb (plot-year combination)

for (i in 1:length(cc.att)) {
  
  ####################   Omnivores' Balanced Diet Plan   #######################
  # OK, so this bit is not easy to read but it works(?)
  # also, it borrows heavily from hand-me-down code
  # that was probably originally written by Barnes
  # (Not sure where this happens in the code of the
  # paper, but it is mentioned in the methods)
  
  # add biomass values in the matrix to 'manually' define the preferences
  # first create a matrix with species biomasses
  mat.bioms = replicate(length(cc.att[[i]]$biomass), cc.att[[i]]$biomass)
  # mat.prefs contains preference of predators based on their prey biomasses
  mat.prefs[[i]] = cc.web[[i]] * mat.bioms
  # then identify omnivorous species
  # these are species feeding at least on:
  # one basal species and at least one non basal species
  # we want them to feed equally on plant animals and detritus
  basals = which(cc.att[[i]]$trophic.level == "basal")
  animals = which(cc.att[[i]]$trophic.level != "basal") ; #S_animals[i]=length(animals)
  plants = which(cc.att[[i]]$trophic.group == "plants")        #; S_plants[i]=length(plants)
  detritus = which(cc.att[[i]]$trophic.group == "detritus")
  #not relevant here but we will need to keep track of herbivores and predators
  herbivores = which(cc.att[[i]]$trophic.level == "herbivore")
  predators = which(cc.att[[i]]$trophic.level == "carnivore")
  
  #omnivores that feed on detritus plants and animals
  omnivores.c.h.d = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])>0 &
                            colSums(mat.prefs[[i]][plants,,drop = FALSE])>0 &
                            colSums(mat.prefs[[i]][animals,])>0)
  
  #omnivores that feed on plants and animals
  omnivores.c.h = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])==0 &
                          colSums(mat.prefs[[i]][plants,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][animals,])>0)
  
  #omnivores that feed on detritus and animals
  omnivores.c.d = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][plants,,drop = FALSE])==0 &
                          colSums(mat.prefs[[i]][animals,])>0)
  
  #omnivores that feed on detritus and plants
  omnivores.h.d = which(colSums(mat.prefs[[i]][detritus,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][plants,,drop = FALSE])>0 &
                          colSums(mat.prefs[[i]][animals,])==0)
  
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][animals, omnivores.c.h.d] = mat.prefs[[i]][animals, omnivores.c.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores.c.h.d,drop=FALSE])),
         length(omnivores.c.h.d),length(omnivores.c.h.d)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][plants, omnivores.c.h.d] = mat.prefs[[i]][plants, omnivores.c.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][plants, omnivores.c.h.d,drop=FALSE])),
         length(omnivores.c.h.d),length(omnivores.c.h.d))
  
  # normalize preferences of omnivores over detritus to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][detritus, omnivores.c.h.d] = mat.prefs[[i]][detritus, omnivores.c.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][detritus, omnivores.c.h.d,drop=FALSE])),
         length(omnivores.c.h.d),length(omnivores.c.h.d))
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][animals, omnivores.c.h] = mat.prefs[[i]][animals, omnivores.c.h,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores.c.h,drop=FALSE])),
         length(omnivores.c.h),length(omnivores.c.h)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][plants, omnivores.c.h] = mat.prefs[[i]][plants, omnivores.c.h,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][plants, omnivores.c.h,drop=FALSE])),
         length(omnivores.c.h),length(omnivores.c.h))
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][animals, omnivores.c.d] = mat.prefs[[i]][animals, omnivores.c.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][animals, omnivores.c.d,drop=FALSE])),
         length(omnivores.c.d),length(omnivores.c.d)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over animals to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][detritus, omnivores.c.d] = mat.prefs[[i]][detritus, omnivores.c.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][detritus, omnivores.c.d,drop=FALSE])),
         length(omnivores.c.d),length(omnivores.c.d)) #diag(4)!=diag(4,1,1) important if single omnivore
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][plants, omnivores.h.d] = mat.prefs[[i]][plants, omnivores.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][plants, omnivores.h.d,drop=FALSE])),
         length(omnivores.h.d),length(omnivores.h.d))
  
  # normalize preferences of omnivores over plants to 1: (sum of prey prefs for omn is 1)
  mat.prefs[[i]][detritus, omnivores.h.d] = mat.prefs[[i]][detritus, omnivores.h.d,drop=FALSE] %*%
    diag(1/colSums(as.matrix(mat.prefs[[i]][detritus, omnivores.h.d,drop=FALSE])),
         length(omnivores.h.d),length(omnivores.h.d))
  
  mat.prefs[[i]] = decostand(mat.prefs[[i]], "total", MARGIN = 2)
  
  # then make total biomass of animal prey species equal to the biomass of detritus:
  #mat.prefs[[i]][!basals, omnivores] =  mat.prefs[[i]][!basals, omnivores] * biomass.plot[1]
  ##############################################################################
  
  ################################# Uncertainty ##################################
  # Here we take each consumer in the foodweb and replaced its fixed preferences 
  # with a random sample from a dirichlet distribution whose component probabilities 
  # are given by the vector of the original preferences. The vector is multiplied 
  # by a scalar that modifies the shape of the distribution (larger = less uncertainty)
  for (j in 1: nrow(mat.prefs[[i]])) { 
    mat.prefs[[i]][,j] = rdirichlet(1, mat.prefs[[i]][,j]*100)
  }
  mat.prefs[[i]][is.nan(mat.prefs[[i]])] = 0 #removes NaNs from basal node "preferences"
  
  ################################################################################
  
  fluxes[[i]] <- fluxing(mat.prefs[[i]],
                         cc.att[[i]]$biomass, 
                         cc.att[[i]]$com.metabolism,
                         cc.att[[i]]$efficiency,
                         bioms.prefs = F,
                         bioms.losses = F,
                         ef.level = "prey")
  
  fluxes[[i]] = fluxes[[i]] * 744 #Joules per month
  
  
  allmetrics2[i,]$year = unique(cc.att[[i]]$year)
  allmetrics2[i,]$plot = unique(cc.att[[i]]$plot)
  allmetrics2[i,]$plant.rich = unique(cc.att[[i]]$plant.div)          # richness gradient
  allmetrics2[i,]$tot.flux = sum(fluxes[[i]])                         # total energy flux              
  allmetrics2[i,]$pred.flux = sum(fluxes[[i]][animals, ])             # predation flux
  allmetrics2[i,]$herb.flux = sum(fluxes[[i]][plants, ])              # herbivory flux
  allmetrics2[i,]$top.down = sum(fluxes[[i]][herbivores, predators])/ # from herbivores per unit herbivore biomass
    sum(cc.att[[i]][herbivores,"biomass"])
  allmetrics2[i,]$bot.up = sum(fluxes[[i]][plants, herbivores])/      # to herbivores per unit herbivore biomass
    sum(cc.att[[i]][herbivores,"biomass"])
  allmetrics2[i,]$herb.press = sum(fluxes[[i]][plants, herbivores])/  # to herbivores per unit plant biomass 
    sum(unique(cc.att[[i]]$plant.biomass))     
  
}

#combine Jena and Cedar Creek
allmetrics = rbind(allmetrics1,allmetrics2)

################################ Multiverse ####################################
# So we had 487 food webs (plot-year-experiment) and for each one we created 1000
# versions. In each one of them, consumers' realised preferences deviate somewhat
# from our expectations (See "Uncertainty"). As a result of this variation, 
# fluxes also vary across the thousand versions of our foodwebs.
# We can now run 1000 regression models to see wether the effect of our predictor,
# plant richness, is consistent given this variation. 
# The simplest way to do this would be to examine the distribution of our point 
# estimates as well as the lower bound of our 95% credible interval in relation 
# to zero. (This may or may not be what Benjamin had in mind) 

thousand[[k]] = allmetrics

############################## Show loop progress ##############################
cat(paste0(round((k/1000)*100), '% completed'))
Sys.sleep(.05)
if (k == 1000) cat(': Done')
else cat('\014')
################################################################################
}

#list.save(thousand, "thousand.rds") #saves the list to wd
#thousand <- readRDS("thousand.rds") #loads the list from wd

for (i in 1:1000) {
  thousand[[i]]$year = as.factor(thousand[[i]]$year)
  thousand[[i]]$plot = as.factor(thousand[[i]]$plot)
  thousand[[i]]$plant.rich = as.integer(thousand[[i]]$plant.rich)
}


# g1 = ggplot(allmetrics, 
#             aes(plant.rich, log(tot.flux), color = experiment)) +     
#   geom_point(position = position_jitterdodge(jitter.width = 0.30), alpha = 0.3)+
#   geom_smooth(method=lm) +
#   theme_classic()
# g2 = ggplot(allmetrics, 
#             aes(plant.rich, log(pred.flux), color = experiment)) +     
#   geom_point(position = position_jitterdodge(jitter.width = 0.30), alpha = 0.3)+
#   geom_smooth(method=lm) +
#   theme_classic()
# g3 = ggplot(allmetrics, 
#             aes(plant.rich, log(herb.flux), color = experiment)) +     
#   geom_point(position = position_jitterdodge(jitter.width = 0.30), alpha = 0.3)+
#   geom_smooth(method=lm) +
#   theme_classic()
# 
# (g1+g2+g3)

confInt = vector(mode = "list",length=1000)
for (i in 1:1000) { 
m = glmmTMB(log(tot.flux) ~ plant.rich + (1|year),
             data = thousand[[i]])
confInt[[i]] = confint(m)
print(confInt[[i]])
}
results = data.frame(intercept = numeric(1000),
                     intercept_2.5 = numeric(1000),
                     intercept_97.5 = numeric(1000),
                     slope = numeric(1000),
                     slope_2.5 = numeric(1000),
                     slope_97.5 = numeric(1000))
for (i in 1:1000) { 
 results[i,1] = confInt[[i]][1,3]
 results[i,2] = confInt[[i]][1,1]
 results[i,3] = confInt[[i]][1,2]
 results[i,4] = confInt[[i]][2,3]
 results[i,5] = confInt[[i]][2,1]
 results[i,6] = confInt[[i]][2,2]
}

long1 = melt(results[,1:3])
long2 = melt(results[,4:6])

p1 = ggplot(long1, aes(x=value, fill = variable)) + 
     geom_density() +
  facet_wrap(~variable, scales="free")

p2 = ggplot(long2, aes(x=value, fill = variable)) + 
  geom_density() +
  facet_wrap(~variable, scales="free")

(p1/p2)

# Because this script was not running forever enough I will now do the 1000 
# models in brms


m.1 <- brm(bf(log(tot.flux) ~ log(plant.rich) + (1|year)),
           chains = 4,
           iter = 2000,
           control = list(adapt_delta = .99),
           cores = n.cores,
           data = allmetrics)

pp = brms::pp_check(m.1)
pp + theme_bw()
plot(m.1)
#plot_model(m.1)
brms::conditional_effects(m.1)
summary(m.1)
