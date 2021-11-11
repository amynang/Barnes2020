library(fluxweb)

#https://figshare.com/articles/software/R_Code_for_Barnes_et_al_Biodiversity_enhances_the_multi-trophic_control_of_arthropod_herbivory_/12909962/1
#The code in the link above is efficient but (to me) it is difficult to see 
#what is going on
#So this is a reframing that works for me. 
#Probably just a matter of personal taste.

#the data
#Jena
download.file("https://figshare.com/ndownloader/files/24558608", 
              destfile = "JenaExp_arthro.csv")
je.com.data = read.csv("JenaExp_arthro.csv")

download.file("https://figshare.com/ndownloader/files/24558617", 
              destfile = "JenaExp_matrix.csv")
je.foodwebmat = as.matrix(read.csv("JenaExp_matrix.csv", row.names = 1))

#Cedar Creek
download.file("https://figshare.com/ndownloader/files/24558605", 
              destfile = "CedCrExp_arthro.csv")
cc.com.data = read.csv("CedCrExp_arthro.csv")

download.file("https://figshare.com/ndownloader/files/24558614", 
              destfile = "CedCrExp_matrix.csv")
cc.foodwebmat = as.matrix(read.csv("CedCrExp_matrix.csv", row.names = 1))

#individual attribute dataframes for each year-plot combination in Jena
je.att = split(je.com.data, with(je.com.data, sample))
#individual matrices for each year-plot combination in Jena
je.web = vector(mode = "list", length=length(je.att))
for (i in 1:length(je.att)) {
  je.web[[i]] = je.foodwebmat[je.att[[i]]$trophic.group,
                              je.att[[i]]$trophic.group]
}

#individual attribute dataframes for each year-plot combination in Cedar Creek
cc.att = split(cc.com.data, with(cc.com.data, sample))
#individual matrices for each year-plot combination in Cedar Creek
cc.web = vector(mode = "list", length=length(cc.att))
for (i in 1:length(cc.att)) {
  cc.web[[i]] = cc.foodwebmat[cc.att[[i]]$trophic.group,
                              cc.att[[i]]$trophic.group]
}

