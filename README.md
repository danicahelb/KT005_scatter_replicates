# KT005_scatter_replicates

```r
rm(list=ls())
setwd("~/Desktop/Microarray Analyses/Kevin's array/151118 data/") #set working directory (folder where Apac array files are stored)
```

##########################################################################################
##Import data
```r
print(load("./Data/Processed Data/0004.KTarray.Apac.transformed_normalized.RData"))
```

##########################################################################################
##clean up annAnti antigen and Unique IDs
```r
annAnti$Unique.ID[annAnti$antigen=="1534"] <- gsub("1534", "15-34", annAnti$Unique.ID[annAnti$antigen=="1534"])
annAnti$Unique.ID[annAnti$antigen=="1538"] <- gsub("1538", "15-38", annAnti$Unique.ID[annAnti$antigen=="1538"])
annAnti$Unique.ID[annAnti$antigen=="1734"] <- gsub("1734", "17-34", annAnti$Unique.ID[annAnti$antigen=="1734"])
annAnti$Unique.ID[annAnti$antigen=="5242"] <- gsub("5242", "52-42", annAnti$Unique.ID[annAnti$antigen=="5242"])

annAnti$antigen[annAnti$antigen=="1534"] <- "15-34"
annAnti$antigen[annAnti$antigen=="1538"] <- "15-38"
annAnti$antigen[annAnti$antigen=="1734"] <- "17-34"
annAnti$antigen[annAnti$antigen=="5242"] <- "52-42"
```

##########################################################################################
##define the spot type
```r
annAnti$type <- "Plasmodium"
annAnti$type[c(grep("blank", annAnti$antigen), grep("PBS", annAnti$antigen))] <- "Neg Ctrl"
annAnti$type[c(grep("std", annAnti$antigen))] <- "Pos Ctrl"
annAnti$type[c(grep("GST", annAnti$antigen))] <- "GST"
```

##########################################################################################
##determine which replicate
```r
annAnti$replicate <- 1
annAnti$replicate[duplicated(annAnti$antigen)] <- 2
annAnti$replicate[annAnti$type %in% c("Neg Ctrl", "GST")] <- NA
```

##########################################################################################
##determine how many times that antigen was repeated
```r
for (i in unique(annAnti$antigen)){
  annAnti[annAnti$antigen==i, "n_replicates"] <- length(annAnti$antigen[annAnti$antigen==i])
}
```


##########################################################################################
##scatter plot of repeats
```r
pdf(file="./Figures/Exploratory/0005.KTarray.Apac.scatter_replicates.pdf",width=12,height=9)
par(mfrow=c(3,4))
for(i in unique(annAnti$antigen[annAnti$n_replicates==2])){
  x <- annAnti$Unique.ID[annAnti$antigen==i & annAnti$replicate==1]
  y <- annAnti$Unique.ID[annAnti$antigen==i & annAnti$replicate==2]
  for(j in c("X1", "X2", "X3")){
    plot(log(MNA.raw[rownames(MNA.raw)==x, grep(j,colnames(MNA.raw))]+0.1),log(MNA.raw[rownames(MNA.raw)==y, grep(j,colnames(MNA.raw))]+0.1),xlab="Replicate 1", ylab="Replicate 2", main=paste("Mean-raw,",i, paste("Apac", j, sep=""), sep=" "), ylim=c(log(0.1),log(80000)), xlim=c(log(0.1),log(80000)), col="red")
    plot(intensity.MNA[rownames(intensity.MNA)==x, grep(j,colnames(intensity.MNA))], intensity.MNA[rownames(intensity.MNA)==y, grep(j,colnames(intensity.MNA))], xlab=NA, ylab=NA, main=paste("Mean-Normalized,",i, sep=" "), ylim=c(0,20), xlim=c(0,20))
    plot(log(MedA.raw[rownames(MedA.raw)==x, grep(j,colnames(MedA.raw))]+0.1),log(MedA.raw[rownames(MedA.raw)==y, grep(j,colnames(MedA.raw))]+0.1),xlab="Replicate 1", ylab="Replicate 2", main=paste("Median-raw,",i, paste("Apac", j, sep=""), sep=" "), ylim=c(log(0.1),log(80000)), xlim=c(log(0.1),log(80000)), col="red")
    plot(intensity.MedA[rownames(intensity.MedA)==x, grep(j,colnames(intensity.MedA))], intensity.MedA[rownames(intensity.MedA)==y, grep(j,colnames(intensity.MedA))] ,xlab=NA, ylab=NA, main=paste("Median-Normalized,", i, sep=" "), ylim=c(0,18), xlim=c(0,18))
  }
}
dev.off()
```

##########################################################################################
##std curves (normalized data)
```r
std <- intensity.MNA[grep("std", rownames(intensity.MNA)),]
std <- std[,order(as.numeric(gsub("b", "", gsub("^X._", "", colnames(std)))))]

pdf(file="./Figures/Exploratory/0005.KTarray.Apac.std_curves.by_slide.normalized.pdf",width=20,height=10)
par(mfrow=c(3,5))
for(i in unique(sort(pheSera$Slide))){
  for (j in c(2:6)){
    temp <- std[,colnames(std) %in% pheSera$SampleID[pheSera$Slide==i]]
    temp <- merge(temp, annAnti[,c("Unique.ID", "description")], all.x=T, by.x=0, by.y="Unique.ID")
    temp$dilution <- as.numeric(gsub(" [(]PBS[])]", "", temp$description))
    plot(x=log(temp$dilution + 0.1), 
         y=temp[,j],
         ylim=c(6,18),
         cex=1.5,
         #col="red", 
         #pch=19,
         xaxt='n',
         #yaxt="n",
         xlab="IgG Concentration (pg/nL)",
         ylab="Intensity",
         main=paste("Slide = ", i, ", Sample = ", names(temp[j]), sep=""))
    axis(side=1,
         at=sort(unique(log(temp$dilution + 0.1))), 
         las=0,
         labels=c("0","0.4", "0.8","1.56","3.13","6.25","12.5", "25", "50", "100"),
         cex.axis=1)
  }
}
dev.off()    
```

##########################################################################################
##std curves (raw data)
```r
std <- MNA.raw[grep("std", rownames(MNA.raw)),]
std <- std[,order(as.numeric(gsub("b", "", gsub("^X._", "", colnames(std)))))]

pdf(file="./Figures/Exploratory/0005.KTarray.Apac.std_curves.by_slide.raw.pdf",width=20,height=10)
par(mfrow=c(3,5))
for(i in unique(sort(pheSera$Slide))){
  for (j in c(2:6)){
    temp <- std[,colnames(std) %in% pheSera$SampleID[pheSera$Slide==i]]
    temp <- merge(temp, annAnti[,c("Unique.ID", "description")], all.x=T, by.x=0, by.y="Unique.ID")
    temp$dilution <- as.numeric(gsub(" [(]PBS[])]", "", temp$description))
    plot(x=log(temp$dilution + 0.1), 
         y=log(temp[,j]+0.1),
         ylim=log(c(500,70000)),
         cex=1.5,
         #col="red", 
         #pch=19,
         xaxt='n',
         #yaxt="n",
         xlab="IgG Concentration (pg/nL)",
         ylab="Log Intensity",
         main=paste("Slide = ", i, ", Sample = ", names(temp[j]), sep=""))
    axis(side=1,
         at=sort(unique(log(temp$dilution + 0.1))), 
         las=0,
         labels=c("0","0.4", "0.8","1.56","3.13","6.25","12.5", "25", "50", "100"),
         cex.axis=1)
  }
}
dev.off()     
```
