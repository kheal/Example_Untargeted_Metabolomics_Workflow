# CAMERA post-processing function --------------------------
# This function processes a CAMERA generated peak table (*.annot) and organizes neutral mass and isotope information.  
# Input is:
#       1. a CAMERA object (xsa; most likely xsaFA)
#       2. a data.frame of the output from CAMERA (.annot). There will be one row for
#       every ion in the original data 
#       Columns:
#             a. the mass feature name (MassFeature), 
#             b. the m/z (mz), 
#             c. the retention time in minutes (RT), 
#             d. any isotopes CAMERA thinks it has found (isotopes), 
#             e. any adducts it thinks it has found (adduct), and 
#             f. the pseudospectrum group (pcgroup). "Pseudospectra" are the 
#             chromatograms that CAMERA hypothesizes could be related, i.e. they
#             follow the same pattern in intensity. Pseudospectra do not necessarily
#             have m/z that would actually make them the same compound, though, so
#             think of grouping the ions into pseudospectra as simply the first 
#             step in making a hypothesis about which ions might be related. 
# Output is a named list containing:
#       1. the same CAMERA object used as input (root part of the name is whatever the input 
#       xcmsSet object was called, suffix is .xsa)
#       2. a data.frame of the output from CAMERA (.annot). 
#       3. a data.frame (suffix is .otherions) where each row is an ion that had 
#       at least one other potentially related ion, so >= 1 other isotope or 
#       adduct. The columns are:
#             a. MassFeature, 
#             b. mz, 
#             c. RT, 
#             d. pseudospectrum group (pcgroup), 
#             e. the type of ion CAMERA thinks this ion might be (IonType), 
#             f. the possible charge of the ion (Charge), 
#             g. the m/z of the other, potentially related ion (MassOfM), 
#             h. and a group number where other ions with the same number
#             are potentially isotopes of this one (IsoGroup). 


camerapostprocess <- function(xset, xset.annot) {
  # Making a column with the mass feature name.
  xset.annot$MassFeature <- paste0("I", round((xset.annot$mz),
                                              digits = 4), "R", 
                                   round(xset.annot$rt/60, 
                                         digits = 2))
  
  # Making a column with the RT in minutes. Note that this is different
  # from the column "rt", which is the RT in seconds. 
  xset.annot$RT <- round(xset.annot$rt/60, digits = 2)
  
  # Retaining only columns of interest in a logical order
  xset.annot <- xset.annot[, c("MassFeature", "mz", "RT", "isotopes",
                               "adduct", "pcgroup")]
  
  # Using regex to extract meaningful info about isotopes and adducts.
  
  # Isotopes
  # Make the column "isotopes" be character data instead of the default, factor.
  # Replace empty cells with "NA".
  xset.annot$isotopes <- as.character(xset.annot$isotopes)
  xset.annot$isotopes[xset.annot$isotopes == ""] <- NA
  
  
  # Take a subset of the data that is only the mass features that appear
  # to have isotopes. 
  Iso <- xset.annot[complete.cases(xset.annot$isotopes), 
                    c("MassFeature", "mz", "RT", "isotopes", 
                      "pcgroup")]
  
  
  # Split the isotopes column into three pieces:
  # 1. xset.annot$IsoGroup (numeric isotope group)
  # 2. xset.annot$IonType (the isotope of that particular ion s/a "M" or "M+1")
  # 3. xset.annot$Charge (the charge of the isotope)
  IsoString <- str_split(Iso$isotopes, "\\]")
  IsoGroup <- sapply(IsoString, function(x) x[[1]])
  IsoGroup <- gsub("\\[", "", IsoGroup)
  
  IonType <- sapply(IsoString, function(x) x[[2]])
  IonType <- gsub("\\[", "", IonType)
  
  Charge <- sapply(IsoString, function(x) x[[3]])
  
  IsoList <- data.frame(MassFeature = Iso$MassFeature,
                        mz = Iso$mz,
                        RT = Iso$RT, 
                        pcgroup = as.numeric(as.character(Iso$pcgroup)),
                        IsoGroup = as.numeric(as.character(IsoGroup)), 
                        IonType = as.character(IonType), 
                        Charge = as.character(Charge),
                        stringsAsFactors = FALSE)
  
  for (i in 1:nrow(IsoList)){
    
    if (IsoList$IonType[i] == "M") {
      IsoList$mzOfI[i] <- IsoList$mz[i]
    } else {
      n <- as.numeric(str_sub(IsoList$IonType[i], 3, 3))
      
      IsoList$mzOfI[i] <- IsoList$mz[i] - n * 1.00866
      
    }
    
  }
  
  IsoList$IonType <- sub("M", "I", IsoList$IonType)
  
  
  # Adducts
  # Make the column "adduct" be character data instead of the default, factor.
  # Replace empty cells with "NA".
  xset.annot$adduct <- as.character(xset.annot$adduct)
  xset.annot$adduct[xset.annot$adduct == ""] <- NA
  
  # Take a subset of the data that is only the mass features that appear
  # to have adduct. 
  Adduct <- xset.annot[complete.cases(xset.annot$adduct), 
                       c("MassFeature", "mz", "RT", "adduct", 
                         "pcgroup")]
  
  # Split the adduct column into one list for every possible adduct with 
  # each list having 4 pieces:
  # 1. IonType = type of adduct, eg. M+Cl
  # 2. MassOfM = the neutral mass of that particular adduct
  # 3. MassFeature 
  # 4. pcgroup
  
  AdSplit <- str_split(Adduct$adduct, "\\[")
  AdList <- list()
  
  for (i in 1:nrow(Adduct)){
    
    IonType <- c()
    Charge <- c()
    MassOfM <- c()
    
    for (m in 2:length(AdSplit[[i]])){
      
      IonType[m-1] <- unlist(str_extract(AdSplit[[i]][m], ".*\\]"))
      IonType[m-1] <- gsub("\\]", "", IonType[m-1])
      
      Charge[m-1] <- unlist(str_extract(AdSplit[[i]][m], "\\].{1,2}"))
      Charge[m-1] <- str_trim(gsub("\\]", "", Charge[m-1]))
      
      MassOfM[m-1] <- str_trim(str_extract(AdSplit[[i]][m], " .*"))
    }
    
    AdList[[i]] <- data.frame(MassFeature = Adduct$MassFeature[i],
                              mz = Adduct$mz[i],
                              RT = Adduct$RT[i],
                              pcgroup = as.numeric(as.character(
                                Adduct$pcgroup[i])),
                              IonType = IonType,
                              Charge = Charge,
                              MassOfM = as.numeric(MassOfM),
                              stringsAsFactors = FALSE)
  }
  
  AdList <- rbind.fill(AdList)
  
  
  OtherIons <- rbind.fill(AdList, IsoList)
  
  IonList <- list(xset, xset.annot, OtherIons)
  #names(IonList) <- c(
  #  paste0(as.character(deparse(substitute(xset))),".xsa"),
  #   paste0(as.character(deparse(substitute(xset))),".annot"),
  #  paste0(as.character(deparse(substitute(xset))),".otherions"))
  
  return(invisible(IonList))
}