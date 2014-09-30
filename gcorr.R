setwd("genetic_corr/")
library(gtools)
read.table("TN_N34.txt",header=TRUE,na.string=".")->TN

TN[TN$Environment == "N3",]->TN3
TN3[,c("Environment","Genotype","Block","SN","SP")]->SNP3
SNP3$T <- SNP3[,length(SNP3)-1] + SNP3[,length(SNP3)]

aov(SP~Block + Genotype, data = SNP3)->SP3
aov(SN~Block + Genotype, data = SNP3)->SN3
aov(T~Block + Genotype, data = SNP3)->T3
summary(SP3)->SP3
summary(SN3)->SN3
summary(T3)->T3

TN[,c("Environment","Genotype","Block","SN","SP")]->SNP
SNP[!(is.na(SNP$SN) | is.na(SNP$SP)),]->SNP
SNP$T <- SNP[,length(SNP)-1] + SNP[,length(SNP)]
aov(SP~Block + Environment * Genotype, data=SNP)->SP
aov(SN~Block + Environment * Genotype, data=SNP)->SN
aov(T~Block + Environment * Genotype, data=SNP)->T
summary(T)[[1]]->T
summary(SP)[[1]]->SP
summary(SN)[[1]]->SN

TN[,c("Environment","Genotype","Block","SN","SN")]->SNN
SNN[!(is.na(SNN$SN) | is.na(SNN$SN)),]->SNN
SNN$T <- SNN[,length(SNN)-1] + SNN[,length(SNN)]
aov(SN~Block + Environment * Genotype, data=SNN)->SPN
#aov(SN~Block + Environment * Genotype, data=SNN)->SNK
aov(T~Block + Environment * Genotype, data=SNN)->TNs
summary(TNs)[[1]]->TNs
summary(SPN)[[1]]->SPN


gcor <- function(data, factors,values,fla,corrf){
  if(is.numeric(values)){
    values<-names(data)[values]
  }
  if(is.numeric(factors)){
    factors<-names(data)[factors]
  }
  geneticCorr <- matrix(nrow=length(values),ncol=length(values))
  rownames(geneticCorr)<-values
  colnames(geneticCorr)<-values
  for(i in 1:(length(values))){
    for(j in 1:length(values)){
  #for(i in 1:(length(values)-1)){
    #for(j in (i+1):length(values)){
      tmp <- data[,c(factors,values[c(i,j)])]
      tmp <- tmp[!(is.na(tmp[,length(tmp)-1]) | is.na(tmp[,length(tmp)])),]
      tmp$T <- tmp[,length(tmp)-1] + tmp[,length(tmp)]
      f <- as.formula(paste(names(tmp)[length(tmp)-2],"~",fla))
      tmpAov<-summary(aov(f,data=tmp))[[1]]
      print(tmpAov)
      r1 <- getGvar(tmpAov)
      f <- as.formula(paste(names(tmp)[length(tmp)-1],"~",fla))
      tmpAov<-summary(aov(f,data=tmp))[[1]]
      #print(tmpAov)
      r2 <- getGvar(tmpAov)
      #r2 <- getGvar(summary(aov(f,data=tmp))[[1]])
      f <- as.formula(paste(names(tmp)[length(tmp)],"~",fla))
      tmpAov<-summary(aov(f,data=tmp))[[1]]
      #print(tmpAov)
      cb <- getGvar(tmpAov)
      #cb <- getGvar(summary(aov(f,data=tmp))[[1]])
      gcov <- (cb -r1 -r2)/2
      gcorr <- gcov/sqrt(r1)/sqrt(r2)
      geneticCorr[values[i],values[j]]<-gcorr[corrf]
      #aovFram = summary(aov(f,data=tmp))[[1]]
      #getGvar(aovFram,"Genotype",c("Block"))
      #print(values[c(i,j)])
      #print(r1)
      #print(r2)
      #print(cb)
      #print(gcov)
      #print(gcorr)
    }
  }
  return(geneticCorr)
}

getGvar <- function(aovFram, Geno,rm){
  effname <- gsub("\\s+","", rownames(aovFram))
  for(i in 1:length(effname)){
    effname[i] <- paste(sort(unlist(strsplit(effname[i],":"))),collapse=":")
  }
  fs <- effname[-grep(":",effname)]
  fs <- fs[fs!="Residuals"]
  ta <- t(aovFram)
  colnames(ta) <- effname
  residuals <- ta[3,c("Residuals")]
  gvar <- vector("numeric",2^length(fs)-1)
  l <- 1;
  for( i in length(fs):1){
    combs <- combinations(length(fs),i)
    for(j in 1:length(combs[,1])){
      combfs <- fs[combs[j,]]
      remains <- fs[-combs[j,]]
      calfs <-paste(sort(combfs),collapse=":")
      if(any(effname==calfs)){
        MSe <- ta[3,calfs] - residuals
        for(k in length(remains):1){
          recombfs <- combinations(length(remains),k)
          for(m in 1 : length(recombfs[,1])){
            subfs <- paste(sort(c(combfs,remains[recombfs[m,]])),collapse=":")
            if(any(names(gvar)==subfs,na.rm=TRUE)){
              MSe <- MSe - gvar[subfs]
            }
          }
        }
        gvar[l] <- MSe
        names(gvar)[l] <- calfs
        l <- l +1
      }
    }
  }
  gvar <- gvar[1:(l-1)]
  #print(gvar)
  for(name in names(gvar)){
    remains <- fs
    for(invar in unlist(split(name,":"))){
      remains <- fs[fs != invar]
    }
    for(r in remains){
      gvar[name] <- gvar[name]/(round(ta[1,r]+1))
    }
  }
 return(gvar)
}


