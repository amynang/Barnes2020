library(fluxweb)


download.file("https://figshare.com/ndownloader/files/24558608", destfile = "JenaExp_arthro.csv")
je.com.data = read.csv("JenaExp_arthro.csv")

download.file("https://figshare.com/ndownloader/files/24558617", destfile = "JenaExp_matrix.csv")
je.foodwebmat = as.matrix(read.csv("JenaExp_matrix.csv", row.names = 1))

download.file("https://figshare.com/ndownloader/files/24558605", destfile = "CedCrExp_arthro.csv")
cc.com.data = read.csv("CedCrExp_arthro.csv")

download.file("https://figshare.com/ndownloader/files/24558614", destfile = "CedCrExp_matrix.csv")
cc.foodwebmat = as.matrix(read.csv("CedCrExp_matrix.csv", row.names = 1))