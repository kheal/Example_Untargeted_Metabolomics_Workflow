#get libraries----
library(plyr) 
library(ggplot2) 
library(gridExtra) 
library(xcms) 
#library(XLConnect) 
library(dplyr) 
library(seqinr) 
library(lubridate)
library(reshape2)

#getTIC----
getTIC <- function(file,rtcor=NULL) {
  object <- xcmsRaw(file)
  cbind(if (is.null(rtcor)) object@scantime else rtcor, rawEIC(object,mzrange=range(object@env$mz))$intensity) 
}

#Calculate CV----
CV <- function(mean, sd){
  (sd/mean)*100
}

#overlay TIC from all files in current folder or from xcmsSet, create pdf----
getTICs <- function(xcmsSet=NULL,files=NULL, pdfname="TICs.pdf",rt=c("raw","corrected")) {
  if (is.null(xcmsSet)) {
    filepattern <- c("[Cc][Dd][Ff]", "[Nn][Cc]", "([Mm][Zz])?[Xx][Mm][Ll]",
                     "[Mm][Zz][Dd][Aa][Tt][Aa]", "[Mm][Zz][Mm][Ll]")
    filepattern <- paste(paste("\\.", filepattern, "$", sep = ""), collapse = "|")
    if (is.null(files))
      files <- getwd()
    info <- file.info(files)
    listed <- list.files(files[info$isdir], pattern = filepattern,
                         recursive = TRUE, full.names = TRUE)
    files <- c(files[!info$isdir], listed)
  } else {
    files <- filepaths(xcmsSet)
  }
  
  N <- length(files)
  TIC <- vector("list",N)
  
  for (i in 1:N) {
    cat(files[i],"n")
    if (!is.null(xcmsSet) && rt == "corrected")
      rtcor <- xcmsSet@rt$corrected[[i]] else 
        rtcor <- NULL
      TIC[[i]] <- getTIC(files[i],rtcor=rtcor)
  }
  
  pdf(pdfname,w=16,h=10)
  cols <- rainbow(N)
  lty = 1:N
  pch = 1:N
  xlim = range(sapply(TIC, function(x) range(x[,1])))
  ylim = range(sapply(TIC, function(x) range(x[,2])))
  plot(0, 0, type="n", xlim = xlim, ylim = ylim, main = "Total Ion Chromatograms", xlab = "Retention Time", ylab = "TIC")
  for (i in 1:N) {
    tic <- TIC[[i]]
    points(tic[,1], tic[,2], col = cols[i], pch = pch[i], type="l")
  }
  legend("topright",paste(basename(files)), col = cols, lty = lty, pch = pch)
  dev.off()
  
  invisible(TIC)
}



#doRetcorPG----
doRetcorPG<- function(xset1){
  retcor(xset1, 
         method="peakgroups", 
         plottype = NULL, 
         missing=1, 
         extra=16,
         smooth="loess",
         family = "symmetric")
}


#doPeakPicking----
doPeakPick <- function(DatFiles){
  xcmsSet(DatFiles,
          method='centWave',
          ppm=Params["PPM", Fraction], 
          peakwidth=c(Params["PEAKWIDTHlow", Fraction], Params["PEAKWIDTHhigh", Fraction]),
          snthresh=Params["SNthresh", Fraction] , 
          mzdiff=Params["MZDIFF", Fraction],
          prefilter=c(Params["PREFILTERlow", Fraction], Params["PREFILTERhigh", Fraction]))
}

#doGrouping----
dogrouping<- function(xset1){
  group(xset1,
        method = "density", 
        bw = Params["BW", Fraction], 
        minfrac = Params["MINFRAC", Fraction],
        minsamp = Params["MINSAMP", Fraction], 
        mzwid = Params["MZWID", Fraction], 
        max = Params["MAX", Fraction])
}

#doRetcor----
doRetcor<- function(xset1){
  retcor(xset1, 
         method="obiwarp", 
         plottype = c("none", "deviation"), 
         profStep=Params["profStep", Fraction], 
         gapInit=Params["GAPInit", Fraction], 
         gapExtend=Params["GAPExtend", Fraction])
}

#checkit ----
checkit <- function(x, ranges){
  which(x>=ranges[,1] & x<=ranges[,2])}


###########Function to find isotope matches-------
isotopehunter8<-function(mzxcms,isotopefile,scantime){
  
  starttime<-Sys.time()
  
  if(scantime[1]=="all"){
    startscan<-1
    endscan<-length(mzxcms@scantime)
  }else{
    startscan<-which(mzxcms@scantime>scantime[1])[1]
    endscan<-tail(which(mzxcms@scantime<scantime[2]),1)
  }
  
  readpattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  
  pattern<-readpattern[(readpattern[,6]=='Y'),]
  uppermass<-pattern[-1,2]+pattern[-1,4]-pattern[1,2]
  lowermass<-pattern[-1,2]-pattern[-1,4]-pattern[1,2]
  upperratio<-pattern[-1,3]*pattern[-1,5]/pattern[1,3]
  lowerratio<-pattern[-1,3]/pattern[-1,5]/pattern[1,3]
  nisotope<-length(uppermass)
  
  columnnames<-c(sapply(pattern[,1],FUN=function(x) c(paste(x,'mass'),paste(x,'intensity'))))
  
  results<-matrix(0,ncol=(4+nisotope*2),nrow=1E6)
  r<-1
  colnames(results) = c('scan','time', columnnames)
  
  
  mzint<-20
  mzbuffer<-uppermass[length(uppermass)]
  
  for(i in startscan:endscan) {
    tscan<-getScan(mzxcms,i)
    tscan<-tscan[which(tscan[,2]>0),]
    
    #print(paste(i,'of',endscan))
    
    first<-tscan[1,1]
    last<-first+mzint
    
    while(last<tscan[nrow(tscan),1]){
      
      scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
      
      if(length(scan)>2){
        final<-length(which(scan[,1]<(last-mzbuffer)))
        for(j in 1:final){
          x<-scan[j,1]
          y<-scan[j,2]
          isotopes<-matrix(NA,ncol=nisotope)
          k<-1
          while(k<nisotope+1){
            isotopes[k]<-(which(
              scan[,1]>(x+lowermass[k]) & 
                scan[,1]<(x+uppermass[k]) & 
                scan[,2]>(y*lowerratio[k]) & 
                scan[,2]<(y*upperratio[k])
            )[1])
            if(is.na(isotopes[k])){k<-(nisotope+2)
            } else {k<-k+1
            }
          }
          if(k==(nisotope+1)){
            results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
            r<-r+1
          }
          
        }
      }
      first<-last-mzbuffer
      last<-last+mzint
    }
    scan<-tscan[which(tscan[,1]>first & tscan[,1]<last),]
    if(length(scan)>2){
      for(j in 1:nrow(scan)){
        x<-scan[j,1]
        y<-scan[j,2]
        isotopes<-matrix(NA,ncol=nisotope)
        k<-1
        while(k<nisotope+1){
          isotopes[k]<-(which(
            scan[,1]>(x+lowermass[k]) & 
              scan[,1]<(x+uppermass[k]) & 
              scan[,2]>(y*lowerratio[k]) & 
              scan[,2]<(y*upperratio[k])
          )[1])
          if(is.na(isotopes[k])){k<-(nisotope+2)
          } else {k<-k+1
          }
        }
        if(k==(nisotope+1)){
          results[r,]<-c(i,mzxcms@scantime[i],x,y,c(t(scan[isotopes,])))
          r<-r+1
        }
      }
    }
  }
  if('O' %in% readpattern[,6]){
    optionalpattern<-readpattern[(readpattern[,6]=='O')]   
    
    for(i in 1:nrow(results)){
      
    }
    
  }
  
  
  print(Sys.time()-starttime)
  
  return(results[1:(r-1),])
}

###########Function to remove noise from isotope matches-------
Removenoise4<-function(result,background,twidth){
  
  masslist2<-c()
  
  if(length(which(result[,4]>background))>1){
    
    #remove elements below x counts (x=background)
    results2<-result[(result[,4]>background),]
    
    #round data for binning
    results3<-round(results2, digits=2)
    
    #remove elements that don't appear twice or more in x chromatographic seconds (x=twidth) and which don't have a modest peak if there are more than 10 points.
    masslist<-unique(results3[which(duplicated(results3[,3])),3])
    
    
    #remove elements that appear in x or more 30 second intervals (x = noise).
    for(i in 1:length(masslist)){
      scantime<-0
      intense<-NULL
      intense2<-NULL
      times<-results3[which(results3[,3]==masslist[i]),2]
      timediff<-times[2:length(times)]-times[1:length(times)-1]
      
      if(min(timediff)<twidth){
        while(scantime<max(result[,2])){
          resultint<-results3[which(results3[,3]==masslist[i] & results3[,2]>scantime & results3[,2]<(scantime+30)),4]
          resultint2<-results3[which(results3[,3]==masslist[i] & results3[,2]>scantime & results3[,2]<(scantime+30)),6]
          if(length(resultint)>0){
            intense<-c(intense,max(resultint))
            intense2<-c(intense2,max(resultint2))
          }
          scantime<-(scantime+30)
        }
        intense<-sort(intense)
        intense2<sort(intense2)
        if(length(intense)<7.5){masslist2<-c(masslist2,masslist[i])}
        if(length(intense)>7.5){
          quartile<-as.integer(length(intense)*3/4)
          threshold1<-mean(intense[1:quartile]+4*sd(intense[1:quartile]))
          threshold2<-mean(intense2[1:quartile]+4*sd(intense[1:quartile]))
          if(max(intense)>threshold1 & max(intense2)>threshold2){masslist2<-c(masslist2,masslist[i])}
        }
      }
    }
  }
  
  ##Now return only the result with the maximum intensity for each element:
  
  result4<-matrix(data=0,nrow=length(masslist2),ncol=ncol(result))
  colnames(result4) = colnames(result)
  
  if(length(masslist2)>0){
    
    for(i in 1:length(masslist2)){
      resultint<-results2[(results3[,3]==masslist2[i]),]
      result4[i,]<-resultint[which.max(resultint[,4]),]
    }
  }
  
  return(result4)
}


#Generate PDF report for isotope search----
XCMSreport<-function(savefile,mzxcms,results,isotopefile,timerange){
  
  masslist<-results[,3]
  
  pdf(paste(savefile,'.pdf', sep = ""))
  for(i in 1:length(masslist)){
    finalplotter_MSonly2(mzxcms,masslist[i],results,isotopefile,timerange,i)
  }
  dev.off()
  write.table(results,paste(savefile,'.txt', sep = ""),sep="\t")
  return()
}

#Final plot function for generating report for any pattern. Includes MS spectra as last panel----
finalplotter_MSonly<-function(mzxcms,mass,results,isotopefile,timerange,i){
  par(mfrow=c(2,1))
  par(mar=c(4,4,3,3))
  mzxcms@scantime<-mzxcms@scantime/60
  timerange<-timerange/60
  isotopeplotter3(mzxcms,isotopefile,mass,timerange,(results[i,2]/60))
  MSplotter2(mzxcms,results,i,isotopefile)
}

#Isotope plotter for final report -----
#Makes one plot with EIC's of every isotope found, scaled to the 'same' intensity
isotopeplotter3<-function(mzxcms,isotopefile,lowmass,timerange,maxscan){
  pattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  
  plotrange<-pattern[,2]-pattern[1,2]+lowmass
  isotoperatio<-pattern[,3]/pattern[1,3]
  namerange<-as.vector(pattern[,1])
  
  colors<-c('blue4','darkorange2','burlywood4','black','red')
  usedcolors<-colors[1:ncol(pattern)]
  
  maxy<-0
  times<-mzxcms@scantime
  
  for(i in 1:nrow(pattern)){
    mzrange<-c(-0.005,0.005)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    EIC1<-unlist(EIC[2])/isotoperatio[i]
    newmax<-max(EIC1[which(times>timerange[1]&times<timerange[2])])
    maxy<-max(c(maxy,newmax))
  }
  
  
  mzrange<-c(-0.005,0.005)+lowmass
  EIC<-rawEIC(mzxcms,mzrange)
  plot(times,EIC[[2]],
       type='l',
       lwd=1,
       xlim=timerange,
       ylim=c(0,maxy*1.1),
       ylab='Scaled Intensity',
       xlab='Retention Time (min)',
       col=usedcolors[1],
       bty='n',
       cex.axis=1,
       cex.lab=1
  )
  
  for(i in 2:ncol(pattern)){
    mzrange<-c(-0.005,0.005)+plotrange[i]
    EIC<-rawEIC(mzxcms,mzrange)
    lines(times,EIC[[2]]/isotoperatio[i],col=usedcolors[i],lwd=1)
  }
  
  abline(v=maxscan, lty=3,col='gray48')
  title(paste('  Isotope peaks LC-ESIMS.  EIC = ',toString(round(plotrange,digits=3))))
  legend('topright',bty="n",namerange,lwd=2,col=usedcolors,cex=0.7,horiz=TRUE)
  
}

#Mass spectra plotter-----
MSplotter2<-function(mzxcms,results,i,isotopefile){
  
  #determine correct isotope pattern
  readpattern<-read.csv(isotopefile, header = TRUE, sep = ",")
  pattern<-readpattern[(readpattern[,6]=='Y'),]
  theorymass<-pattern[,2]-pattern[1,2]+results[i,3]
  theoryratio<-pattern[,3]/pattern[1,3]*results[i,4]
  
  isotopes<-matrix(results[i,-1:-2],byrow=TRUE,ncol=2)
  massrange<-c(mean(isotopes[,1])-10,mean(isotopes[,1])+10)
  scanrange<-getScan(mzxcms,results[i,1],mzrange=massrange)
  #scanrange<-scan[which(scan[,1]>mzrange[1] & scan[,1]<mzrange[2]),]
  plot(scanrange,type='h',ylim=c(0,max(scanrange[,2])))
  lines(theorymass,theoryratio,type='h',col='cornsilk3',lwd=3)
  lines(isotopes,type='h',col='red')
  title(paste('MS: ',toString(round(results[i,2]/60,digits=2)),' min'))
  
}

#KRH Peak Peeker Function----------
KRHS1picker <- function(DFFiltered){
  print(paste("Get ready to look at", length(DFFiltered$scan), "peaks!"))
  DFFiltered <- split(DFFiltered, DFFiltered$SampleID)
  for (j in 1:length(DFFiltered)){ #Loop through # of mzxml files to open and shut
    print(paste("Opening next sample, please wait, this is number", j, "of", length(DFFiltered)))
    cleanresult <- DFFiltered[[j]]
    mzdatafiles <- DFFiltered[[j]]$mzFile
    mzxcms <- xcmsRaw(mzdatafiles[1],profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
    masslist32S<-cleanresult$X32S.mass
    masslist34S<-cleanresult$X34S.mass
    timelist <- cleanresult$scan
    cleanlist <- cleanresult$GoodPeak
    intensitylist <- cleanresult$X32S.intensity
    namelist <- cleanresult$SampleID
    
    for (i in 1:length(timelist)) {   
      EICS32<-rawEIC(mzxcms,mzrange=c(masslist32S[i]-0.005,masslist32S[i]+0.005))#, timerange = c(timelist[i]-100, timelist[i]+100))
      EICS34<-rawEIC(mzxcms,mzrange=c(masslist34S[i]-0.005,masslist34S[i]+0.005))#, timerange = c(timelist[i]-100, timelist[i]+100)
      plot(EICS32[[2]],
           xlim = c(timelist[i]-300, timelist[i]+300),
           ylim = c(0, (2*intensitylist[i])),
           type='l',
           lty=3,
           lwd=2,
           ylab='Scaled Intensity',
           xlab='scan number',
           col="black",
           bty='n',
           cex.axis=1,
           cex.lab=1, 
           main = namelist[i])
      lines(22.12821*EICS34[[2]],
            xlim = c(timelist[i]-300, timelist[i]+300),
            ylim = c(0, (2*intensitylist[i])),
            type='l',
            lwd=2,
            col="orange",
            bty='n')
      abline(v=timelist[i], col='cyan',lty=3, lwd=6)
      ask<-readline(prompt="Enter 'y' if this is a good hit: ")
      if(ask=='y'){cleanlist[i]="GOOD"}else{cleanlist[i]="BAD"}
    }
    
    DFFiltered[[j]]$GoodPeak <- cleanlist
  }
  DFFiltered2 <- do.call(rbind, DFFiltered)
  write.csv(DFFiltered2,  paste(as.character(Dirs[as.character(Fraction) , "SResults"]), 
                                "/Culled_S1_List_wtihTargets_Checked.csv", sep = "", collapse = NULL))
}

KRHS1pickerOrg <- function(DFFiltered){
  print(paste("Get ready to look at", length(DFFiltered$scan), "peaks!"))
  DFFiltered <- split(DFFiltered, DFFiltered$SampleID)
  for (j in 1:length(DFFiltered)){ #Loop through # of mzxml files to open and shut
    print(paste("Opening next sample, please wait, this is number", j, "of", length(DFFiltered)))
    cleanresult <- DFFiltered[[j]]
    mzdatafiles <- DFFiltered[[j]]$mzFile
    mzxcms <- xcmsRaw(mzdatafiles[1],profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
    masslist32S<-cleanresult$X32S.mass
    masslist34S<-cleanresult$X34S.mass
    timelist <- cleanresult$scan
    cleanlist <- cleanresult$GoodPeak
    intensitylist <- cleanresult$X32S.intensity
    namelist <- cleanresult$SampleID
    
    for (i in 1:length(timelist)) {   
      EICS32<-rawEIC(mzxcms,mzrange=c(masslist32S[i]-0.005,masslist32S[i]+0.005))#, timerange = c(timelist[i]-100, timelist[i]+100))
      EICS34<-rawEIC(mzxcms,mzrange=c(masslist34S[i]-0.005,masslist34S[i]+0.005))#, timerange = c(timelist[i]-100, timelist[i]+100)
      plot(EICS32[[2]],
           xlim = c(timelist[i]-300, timelist[i]+300),
           ylim = c(0, (2*intensitylist[i])),
           type='l',
           lty=3,
           lwd=2,
           ylab='Scaled Intensity',
           xlab='scan number',
           col="black",
           bty='n',
           cex.axis=1,
           cex.lab=1, 
           main = namelist[i])
      lines(22.12821*EICS34[[2]],
            xlim = c(timelist[i]-300, timelist[i]+300),
            ylim = c(0, (2*intensitylist[i])),
            type='l',
            lwd=2,
            col="orange",
            bty='n')
      abline(v=timelist[i], col='cyan',lty=3, lwd=6)
      ask<-readline(prompt="Enter 'y' if this is a good hit: ")
      if(ask=='y'){cleanlist[i]="GOOD"}else{cleanlist[i]="BAD"}
    }
    
    DFFiltered[[j]]$GoodPeak <- cleanlist
  }
  DFFiltered2 <- do.call(rbind, DFFiltered)
  write.csv(DFFiltered2,  paste(as.character(Dirs[as.character(Fraction) , "SResults"]), "/", Org ,
                                "_Culled_S1_List_wtihTargets_Checked.csv", sep = "", collapse = NULL))
}

#Final plot function for generating report for any pattern. Includes MS spectra as last panel
finalplotter_MSonly<-function(mzxcms,mass,results,isotopefile,timerange,i){
  par(mfrow=c(2,1))
  par(mar=c(4,4,3,3))
  mzxcms@scantime<-mzxcms@scantime/60
  timerange<-timerange/60
  isotopeplotter3(mzxcms,isotopefile,mass,timerange,(results[i,2]/60))
  MSplotter2(mzxcms,results,i,isotopefile)
}

#Final plot function for generating report for any pattern. Includes Apo form MS spectra as last panel
finalplotter_MSonly2<-function(mzxcms,mass,results,isotopefile,timerange,i){
  par(mfrow=c(3,1))
  par(mar=c(4,4,3,3))
  mzxcms@scantime<-mzxcms@scantime/60
  timerange<-timerange/60
  isotopeplotter3(mzxcms,isotopefile,mass,timerange,(results[i,2]/60))
  if(isotopefile=='./Feisotope.csv'){EICplot1(mzxcms,mass-50.9161,'Apo form',timerange,(results[i,2]/60))}
  if(isotopefile=='./Cuisotope.csv'){EICplot1(mzxcms,mass-60.9139,'Apo form',timerange,(results[i,2]/60))}
  MSplotter2(mzxcms,results,i,isotopefile)
}




#Check known compounds
checkKnownCompounds <- function(MFs, ppmtol, rttolRP, rttolHILIC) {
  if(missing(ppmtol)) {
    ppmtol <- 15
  }
  if(missing(rttolRP)) {
    rttolRP <- 0.5
  }
  if(missing(rttolHILIC)) {
    rttolHILIC <- 1
  }
  matchedKnownCompounds <- list()
  knownShortCompounds <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/Ingalls_Lab_Standards.csv"), header = T)  %>%
    mutate(mz = m.z, RT = RT..min.) %>% select(Compound.Name, mz, RT, Fraction1, Fraction2) %>% 
    mutate(Fraction1 = as.character(Fraction1),
           Fraction2 = as.character(Fraction2))
  MFstry <- MFs %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_") %>% select(-MF)  
  matchedKnownCompounds[[1]] <- difference_inner_join(x= MFstry, y = knownShortCompounds, 
                                                      by = "mz", max_dist = .01,  distance_col = NULL) %>%
    mutate(MZdiff = abs(mz.x-mz.y ))%>%
    mutate(RTdiff = abs(RT.x-RT.y),
           ppm = (MZdiff/mz.x *10^6)) %>%
    filter(ppm < ppmtol) %>%
    filter(Frac == Fraction1 | Frac == Fraction2) %>%
    mutate(Frac = as.factor(Frac), mz = mz.x) %>%
    filter(!((Frac == "CyanoAq" | Frac == "CyanoDCM") &  RTdiff > rttolRP)) %>%
    filter(!((Frac == "HILICNeg" | Frac == "HILICPos") &  RTdiff > rttolHILIC))   %>%
    select(MF_Frac, mz, rt, Compound.Name, RTdiff, ppm)
  matchedKnownCompounds[[2]] <- matchedKnownCompounds[[1]] %>% 
    arrange(RTdiff) %>%
    group_by(MF_Frac) %>%
    summarise(TargetMatches = as.character(paste(Compound.Name,  collapse="; ")),
              ppmMatches = as.character(paste(ppm,  collapse="; ")),
              RTdiffMatches = as.character(paste(RTdiff,  collapse="; ")))
  return(matchedKnownCompounds)
  
}

#flag for known Contaminants -------
checkContaminants <- function(MFs, ppmtol) {
  if(missing(ppmtol)) {
    ppmtol <- 15
  }
  matchedKnownContaminants <- list()
  knownContaminants <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/CommonContams.csv"), comment.char = "#", header = T)%>%
    mutate(mz = m.z) %>% select(Fraction1:Fraction3, mz, Compound) %>% mutate(Flag = "PossibleContamin")
  MFstry <- MFs %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_") %>% select(-MF)
  matchedKnownContaminants[[1]] <- difference_inner_join(x= knownContaminants, y = MFstry, 
                                                         by = "mz", max_dist = .01,  distance_col = NULL) %>%
    mutate(ppm = (abs(mz.x-mz.y )/mz.x *10^6)) %>%
    filter(ppm < ppmtol, Frac == Fraction1 | Frac == Fraction2 | Frac == Fraction3) %>%
    mutate(Frac = as.factor(Frac), mz = mz.x) %>% 
    select(MF_Frac, Compound, ppm, Flag ) %>% 
    rename(ContamMatches = Compound, Contamppm = ppm, ContamFlag = Flag)
  
  matchedKnownContaminants[[2]] <- matchedKnownContaminants[[1]] %>% 
    group_by(MF_Frac) %>%
    summarise(ContamMatches = as.character(paste(ContamMatches,  collapse="; ")),
              Contamppm = as.character(paste(Contamppm,  collapse="; ")))
  return(matchedKnownContaminants)
}

#Check against the KEGG list of compounds -----
checkKEGG <- function(MFs, ppmtol) {
  if(missing(ppmtol)) {
    ppmtol <- 5
  }
  matchedKEGGs <- list()
  keggCompounds <- read.csv(text = getURL("https://raw.githubusercontent.com/kheal/Example_Untargeted_Metabolomics_Workflow/master/KEGGCompounds_withMasses.csv"), header = T) %>% rename(Compound= OtherCmpds)
  keggPos <- keggCompounds %>% select(Compound, PosMZ) %>% 
    mutate(Fraction1 = "HILICPos", Fraction2="CyanoAq", Fraction3= "CyanoDCM") %>% unique() %>% 
    rename(mz = PosMZ)
  keggCompoundsShort <- keggCompounds %>% select(Compound, NegMZ) %>% 
    mutate(Fraction1 = "HILICNeg", Fraction2=NA, Fraction3= NA) %>% unique() %>%
    rename(mz = NegMZ) %>% rbind(keggPos)
  MFstry <- MFs %>% mutate(MF_Frac2 = MF_Frac) %>% separate(MF_Frac2, c("MF", "Frac"), sep =  "_") %>% select(-MF)
  matched <- difference_inner_join(x= keggCompoundsShort, y = MFstry, 
                                   by = "mz", max_dist = .01,  distance_col = NULL) %>%
    mutate(ppm = (abs(mz.x-mz.y )/mz.x *10^6)) %>%
    filter(ppm < ppmtol, Frac == Fraction1 | Frac == Fraction2 | Frac == Fraction3) %>%
    mutate(Frac = as.factor(Frac),mz = mz.x) %>% 
    select(MF_Frac, Compound, ppm ) %>% rename(Keggppm = ppm)
  
  matchedNames <- left_join(matched, keggCompounds) %>% select(Compound, Name) %>%
    group_by(Compound) %>%
    summarise(KEGGMatchesNames = as.character(paste(Name,  collapse=" ")))
  matchedKEGGs[[1]] <- left_join(matched, matchedNames) %>% rename(KeggMatches = Compound)
  matchedKEGGs[[2]] <- matchedKEGGs[[1]] %>%
    arrange(Keggppm) %>%
    group_by(MF_Frac) %>%
    summarise(KeggMatches = as.character(paste(KeggMatches,  collapse="; ")),
              Keggppm = as.character(paste(Keggppm,  collapse="; ")),
              KeggNames = as.character(paste(KEGGMatchesNames,  collapse="; ")))
  return(matchedKEGGs)}



#Get max RT function----- 
maxRT <- function(mz, rt, mzxcms, compound, xset3, XcmsIndex){
  EICinfo<-rawEIC(mzxcms, mzrange=c(mz-0.005,mz+0.005))
  SubEICdat <- data.frame("intensity"=EICinfo[[2]],"time"=mzxcms@scantime, "correcttime" = xset3@rt[[2]][[XcmsIndex]]) %>%
    filter(time > (rt-100) & time < (rt+100))
  if(sum(SubEICdat$intensity)==0){return(NA)}else{
    maxtime <- SubEICdat %>%
      filter(intensity == max(intensity)) %>%
      select(time) %>%
      as.numeric()
    maxcorrecttime <- SubEICdat %>%
      filter(intensity == max(intensity)) %>%
      select(correcttime) %>%
      as.numeric()
    plot(x = SubEICdat$time, y = SubEICdat$intensity,
         type='l', 
         xlim = c(rt-300, rt+300), 
         ylab='Intensity',
         xlab='time',
         main = compound)
    abline(v=rt, col='cyan',lty=3, lwd=6)
    abline(v= maxtime, col = 'red', lty = 3, lwd =6)
    ask<-readline(prompt="Enter 'y' if this is a good match: ")
    if(ask=='y'){return(maxcorrecttime)}else{return(NA)}}
}


                        
                        
#Pull an MS2 spectra from a DDA file according to the a targeted mass and time------
                        #Sums intensities and filters out fragments like in Tabb2003
getms2spectra <-
  function(xs=xs, mass=mass, time=time, timetol=15, masstolLR=0.4, masstolHR=0.02){
  #xs is msn2xcmsRaw(xcmsRaw(DDAFILE, includeMSn=TRUE))
  masslow<-mass-masstolLR
  masshigh<-mass+masstolLR
  masslowhr <- mass - masstolHR
  masshighr <- mass + masstolHR
  timelow<-time-timetol
  timehigh<-time+timetol
  scanlist<-which(xs@msnPrecursorMz>masslow & xs@msnPrecursorMz<masshigh & xs@scantime>timelow & xs@scantime<timehigh)
  allMS2list <- list()
  if(length(scanlist)>0){
    for (l in 1:length(scanlist)){
      ms2dat <- as.data.frame(getScan(xs,scanlist[l])) %>% mutate(scannumber = scanlist[l])
      allMS2list[[l]] <- ms2dat} #This loop gets all matching scans with LR
      allMS2_matched <- do.call(rbind, allMS2list) %>%
        filter(mz > masslowhr & mz < masshighr) 
      matchedscanlist <- allMS2_matched$scannumber} else {matchedscanlist <- c()} #Do we find any HR matches?  
      if(length(matchedscanlist > 0)){
        MS2Summed <- do.call(rbind, allMS2list) %>%
          filter(scannumber %in% allMS2_matched$scannumber) %>%
          mutate(mzRND = as.factor(round(mz, digits = 2))) %>%
          group_by(mzRND) %>%
          summarise(mz = mean(mz),
                    intensity = sum(intensity),
                    sqrtintensity = sqrt(intensity),
                    scancount = n())
        bestscan <- MS2Summed %>%
          arrange(desc(intensity)) %>%
          mutate(intensity = round(intensity/max(intensity)*100, digits = 1))%>%
          filter(intensity > 0.5) %>%
          mutate(mz = round (mz, digits = 5),
                 mash = paste(mz, intensity, sep = ", " ))
        sortedscanrange <- paste(bestscan$mash, collapse = "; ")
        return(sortedscanrange)}else{return(NA)} #Close if/then for scanlist >1
  
  }
#Close Funtion

                        
#Get a scan table from a concatenated scanlist---- (Scan is output from getms2spectra function)
 scantable <- function(Scan){
  Try <- read.table(text=as.character(Scan),
                  col.names=c('mz','intensity')) %>% 
          mutate(mz = as.numeric(mz %>% str_replace(",", "")),
                 intensity = as.numeric(intensity %>% str_replace(";", "")))
return(Try)
  }                       

                        
############## Cosine scoring w/ spectral peaks-----
                        #From Rene Boiteau
MSMScosine_1<-function(scan1,scan2,mass1,mass2){

  mztolerance<-0.02
  
  w1<-(scan1[,1]^2)*sqrt(scan1[,2])
  w2<-(scan2[,1]^2)*sqrt(scan2[,2])
  
  diffmatrix<-sapply(scan1[,1], function(x) scan2[,1]-x)
  sameindex<-which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity)
}



############## Cosine scoring w/ neutral losses-----
                     #From Rene Boiteau
MSMScosine_2<-function(scan1,scan2,mass1,mass2){

  mztolerance<-0.02
  
  loss1<-scan1
  loss1[,1]<-mass1-loss1[,1]
  loss2<-scan2
  loss2[,1]<-mass2-loss2[,1]
  w1<-(loss1[,1]^2)*sqrt(loss1[,2])
  w2<-(loss2[,1]^2)*sqrt(loss2[,2])
  
  diffmatrix<-sapply(loss1[,1], function(x) loss2[,1]-x)
  sameindex<-which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity2<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity2)
}
                     
                     
#  QUICKLY LOOK AT MANY PEAKS TO SEE IF THEY ARE HIGH QUALITY ------
### xs = xcmsRaw(DATAFILE,profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
### XCMSIndex = index # in xset3 that matches your DATAFILE
### xset3 = RT corrected xset that includes the DATAFILE
### rt, rtmin, rtmax = all in seconds
                     
peakPeeker <- function(xs, mz, rt, XcmsIndex, xset3){
    EICinfo<-rawEIC(xs, mzrange=c(mz-0.005,mz+0.005))
    SubEICdat <- data.frame("intensity"=EICinfo[[2]],
                            "time"=xs@scantime, 
                            "correcttime" = xset3@rt[[2]][[XcmsIndex]]) %>%
      filter(time > (rt-100) & time < (rt+100))
    plot(x = SubEICdat$correcttime, y = SubEICdat$intensity,
         type='l', 
         ylab='Intensity',
         xlab='time')
    abline(v=rt, col='cyan',lty=3, lwd=2)
    ask<-readline(prompt="Enter 'y' if this is a good peak: ")
    if(ask=='y'){return("YES")}else{return("NO")}
}
                     
#Check against MassBank Spectrum (m/z and spectra match for putative ID)
#Need to have the Spectra in your working set, which is parsed from MONA, KRH has copy, but can't upload -- too big
massbankMS2MatchPOSITIVE <- function(ShortestDat){
  mz <- as.numeric(ShortestDat["mz"])
  MS2 <- as.character(ShortestDat["MS2"])
  MF_Frac <- as.character(ShortestDat["MF_Frac"])
  
  Candidates <- Spectra %>%
    mutate(MH_mass = M_mass + 1.0072766) %>%
    filter(near(MH_mass, mz, tol= 0.02)) %>%
    mutate(scan1 = spectrum_KRHform_filtered, 
           scan2 = MS2, mass1 = MH_mass, mass2 = mz)%>%
    filter(!is.na(scan1))
  
   NoMatchReturn <- Spectra %>% 
        mutate(MassBankMatch = NA,
               MassBankppm = NA,
               MassBankCosine1 = NA,
               MF_Frac = MF_Frac) %>%
        select(MF_Frac, MassBankMatch:MassBankCosine1) %>% head(1)
    
  if (length(Candidates$ID > 1)){
      Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMSconsine1_df(x))
      Candidates$Cosine2 <- apply(Candidates, 1, FUN=function(x) MSMSconsine2_df(x))
      Candidates <- Candidates %>% filter(Cosine1 > 0.8) %>% arrange(desc(Cosine1)) %>% 
        mutate(MF_Frac = MF_Frac) %>% head(1)
      Candidates <- Candidates %>% 
        mutate(MassBankMatch = paste(Names, ID, sep= " ID:"),
               MassBankppm = abs(mass1-mass2)/mass2 *10^6,
               MassBankCosine1 = Cosine1) %>%
        filter(MassBankppm < 5) %>%
        select(MF_Frac, MassBankMatch:MassBankCosine1)
      }
      if (length(Candidates$MF_Frac) == 0){
      Candidates <-NoMatchReturn
      }
      return(Candidates)
  }

massbankMS2MatchNEGATIVE <- function(ShortestDat){
  mz <- as.numeric(ShortestDat["mz"])
  MS2 <- as.character(ShortestDat["MS2"])
  MF_Frac <- as.character(ShortestDat["MF_Frac"])
  
  Candidates <- Spectra %>%
    mutate(MH_mass = M_mass - 1.0072766) %>%
    filter(near(MH_mass, mz, tol= 0.02)) %>%
    mutate(scan1 = spectrum_KRHform_filtered, 
           scan2 = MS2, mass1 = MH_mass, mass2 = mz)%>%
    filter(!is.na(scan1))
  
  NoMatchReturn <- Spectra %>% 
        mutate(MassBankMatch = NA,
               MassBankppm = NA,
               MassBankCosine1 = NA,
               MF_Frac = MF_Frac) %>%
        select(MF_Frac, MassBankMatch:MassBankCosine1) %>% head(1)
    
  if (length(Candidates$ID > 1)){
      Candidates$Cosine1 <- apply(Candidates, 1, FUN=function(x) MSMSconsine1_df(x))
      Candidates$Cosine2 <- apply(Candidates, 1, FUN=function(x) MSMSconsine2_df(x))
      Candidates <- Candidates %>% filter(Cosine1 > 0.8) %>% arrange(desc(Cosine1)) %>% 
        mutate(MF_Frac = MF_Frac) %>% head(1)
      Candidates <- Candidates %>% 
        mutate(MassBankMatch = paste(Names, ID, sep= " ID:"),
               MassBankppm = abs(mass1-mass2)/mass2 *10^6,
               MassBankCosine1 = Cosine1) %>%
        filter(MassBankppm < 5) %>%
        select(MF_Frac, MassBankMatch:MassBankCosine1)
      }
      if (length(Candidates$MF_Frac) == 0){
      Candidates <-NoMatchReturn
      }
      return(Candidates)
  }

#These functions can be used with aply functions, otherwise the same as MSMScosine1 and 2
MSMSconsine1_df <- function(df){
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- df["mass1"]
  mass2 <- df["mass2"]

  mztolerance<-0.02
  
  w1<-(scan1[,1]^2)*sqrt(scan1[,2])
  w2<-(scan2[,1]^2)*sqrt(scan2[,2])
  
  diffmatrix<-sapply(scan1[,1], function(x) scan2[,1]-x)
  sameindex<-which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity)
}

MSMSconsine2_df <- function(df){
  scan1 <- scantable(df["scan1"])
  scan2 <- scantable(df["scan2"])
  mass1 <- as.numeric(df["mass1"])
  mass2 <- as.numeric(df["mass2"])

  mztolerance<-0.02
  
  loss1<-scan1
  loss1[,1]<-mass1-loss1[,1]
  loss2<-scan2
  loss2[,1]<-mass2-loss2[,1]
  w1<-(loss1[,1]^2)*sqrt(loss1[,2])
  w2<-(loss2[,1]^2)*sqrt(loss2[,2])
  
  diffmatrix<-sapply(loss1[,1], function(x) loss2[,1]-x)
  sameindex<-which(abs(diffmatrix)<mztolerance,arr.ind=T)
  
  similarity2<-sum(w1[sameindex[,2]]*w2[sameindex[,1]])/(sqrt(sum(w2^2))*sqrt(sum(w1^2)))
  
  return(similarity2)
}
                     
####For 2nd search for S-containing compounds            
KRHS1pickerOrg_SecondSearch  <- function(DFFiltered){
  print(paste("Get ready to look at", length(DFFiltered$scan), "peaks!"))
  DFFiltered <- split(DFFiltered, DFFiltered$SampleID)
  for (j in 1:length(DFFiltered)){ #Loop through # of mzxml files to open and shut
    print(paste("Opening next sample, please wait, this is number", j, "of", length(DFFiltered)))
    cleanresult <- DFFiltered[[j]]
    cleanresult$GoodPeak <- NA
    mzdatafiles <- DFFiltered[[j]]$mzFile
    mzxcms <- xcmsRaw(mzdatafiles[1],profstep=0.01,profmethod="bin",profparam=list(),includeMSn=FALSE,mslevel=NULL, scanrange=NULL)
    masslist32S<-cleanresult$X32S.mass
    masslist34S<-cleanresult$X34S.mass
    timelist <- cleanresult$scan
    cleanlist <- cleanresult$GoodPeak
    intensitylist <- cleanresult$X32S.intensity
    namelist <- cleanresult$SampleID
    
    for (i in 1:length(timelist)) {   
      EICS32<-rawEIC(mzxcms,mzrange=c(masslist32S[i]-0.005,masslist32S[i]+0.005))#, timerange = c(timelist[i]-100, timelist[i]+100))
      EICS34<-rawEIC(mzxcms,mzrange=c(masslist34S[i]-0.005,masslist34S[i]+0.005))#, timerange = c(timelist[i]-100, timelist[i]+100)
      plot(EICS32[[2]],
           xlim = c(timelist[i]-200, timelist[i]+200),
           ylim = c(0, (2*intensitylist[i])),
           type='l',
           lwd=2,
           ylab='Scaled Intensity',
           xlab='scan number',
           col="black",
           bty='n',
           cex.axis=1,
           cex.lab=1, 
           main = namelist[i])
      lines(22.12821*EICS34[[2]],
            xlim = c(timelist[i]-300, timelist[i]+300),
            ylim = c(0, (2*intensitylist[i])),
            type='l',
            lwd=1,
            col="orange",
            bty='n')
      abline(v=timelist[i], col='cyan',lty=3, lwd=6)
      ask<-readline(prompt="Enter 'y' if this is a good hit: ")
      if(ask=='y'){cleanlist[i]="GOOD"}else{cleanlist[i]="BAD"}
    }
    
    DFFiltered[[j]]$GoodPeak <- cleanlist
  }
  DFFiltered2 <- do.call(rbind, DFFiltered)
  return(DFFiltered2)
}
                                  
                                  
                                  
                                  
