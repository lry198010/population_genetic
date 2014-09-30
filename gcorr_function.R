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
      #print(tmpAov)
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
