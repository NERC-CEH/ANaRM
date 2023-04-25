
##### Adapted Nature-based-solutions Rational Method (ANaRM) - Main program file #####
# 
# Available for use under the terms of the included MIT licence. 
# Please consider citing:  
# J.D Miller, G. Vesuviano, J.R. Wallbank, D.H. Fletcher, L. Jones (2023). 
# "Hydrological assessment of urban Nature-Based Solutions for urban planning using Ecosystem Service toolkit applications"
# Landscape and Urban Planning 234  104737


##### Instructions #####  
#1. Check the options file is in the directory specified by getwd() or else edit the variable "optionsFN" below
#2. Check the options file (options.inp) and runnoff.csv are correct.
#3. Run this script.
#4. See Miller 2023 and README.md for further detail, or contact johwal@ceh.ac.uk

#### Be aware that:    #####
# 1.  World cover and out-flow may disagree - e.g. lake catchment can be missed due to how outflow directions work for lake (see BARTLEY RSV), OR the lake catchment may need enlarging especially if option "FARL" (original) is used (min ccar = 1.2 * lake_area is used in that case). TO DO: A function to modify outf based on lake locations could be produced. 
# 1b. This could be especially problematic for a lake change scenario - the outf and ccar grids won't be updated.
# 2.  world cover can count broad rivers as "water". Not a problem for Rea at	Calthorpe Park though.
# 3.  DTMs may or may not be particularly "hydrological" .
# 4.  There is a small modification to the definition of FARL so that all lakes can be considered "online".

#### Changes : ####
# Apr23 : Added the option "removeEdgeCats". This removes all cells downstream from the edge of either the outflow grid or cells with no outflow direction.
# Apr23 : Fixed bug occuring when multiple lakeId's specified in options file.

##### John R Wallbank, UK Centre for Ecology & Hydrology (UKCEH), johwal@ceh.ac.uk, April 2023.


rm(list=ls()) 
Sys.time()
###################################################
## EDIT LINE BELOW TO NAME OPTIONS FILE IF NEEDED:
optionsFN = paste0(getwd(),"/options.inp")
###################################################

options(max.print=900)
options(stringsAsFactors = FALSE)

library(raster)#install.packages("raster")
library(matlab)#install.packages("matlab")

#reads an options file and evaluates all options. 
#If "FN" (File Name) in a variable (and not empty) then this is made relative to baseDir.
readOptions<-function(optionsFN,verbose=F){
  optLines =trimws(sapply(readLines(optionsFN),function(x) gsub("#.*", "", x),USE.NAMES = F))
  optLines=optLines[optLines!=""]
  optLines=optLines[grepl("=",optLines)]
  for(ln in optLines){#ln = optLines[5]
    tmp = trimws(strsplit(ln,"=")[[1]])
    if( grepl("FN",tmp[1]) ){#ie a file
      if(length(tmp)==1){
        assign( tmp[1],"" , envir = .GlobalEnv)
      }else{
        assign( tmp[1], paste0(baseDir,"/",tmp[2]), envir = .GlobalEnv)
      }
    }else{
      numericalNow = suppressWarnings(as.numeric(tmp[2]))
      logicalNow = suppressWarnings(as.logical(tmp[2]))
      if(!is.na(numericalNow)){
        assign( tmp[1], numericalNow, envir = .GlobalEnv)
      }else if(!is.na(logicalNow)){
        assign( tmp[1], logicalNow, envir = .GlobalEnv)
      }else{
        if(tmp[1]=="baseDir" & is.na(tmp[2]) ) tmp[2] = getwd()
        assign( tmp[1], tmp[2], envir = .GlobalEnv)
        if(tmp[1]=="lakeId"  )assign( tmp[1], as.numeric( strsplit( gsub( ")", "",gsub("c(", "",tmp[2],fixed=T),fixed=T),"," )[[1]] ), envir = .GlobalEnv) 
      }
    }
  }
  cat("Read",length(optLines),"options from",optionsFN,"\n")
  if(verbose)print(optLines)
}

################### read options #####################
optDF = readOptions(optionsFN,verbose=T)
cat("\n**** Running with options file:",optionsFN,"****\n\n")
source(paste0(baseDir,"/ANaRM_functions.R"))

#convert rain intensity to m^3 s^-1 (= mm/hr * 1/60^2 hr/s * 1/1000 m/mm )
iRain = iRain*( (1/60^2) * (1/1000)) 

#load runoff coef (percentage) for each wrldcover class
cvDF = read.csv(runoffFN)
if(any(cvDF[lakeId,"cv"]!=100))print("WARNING: cv for lake is not 100")


############# Part 1 ###############
## Load OUTF, and base run-off grid. Crop and pad them

#Outflow directions are used as reference grid
refRast = raster(outfFN)
refCrs  = toString(crs(refRast))
refBnds = getBnds(refRast) 

#load/make outf, lcmMP,  and roff rasters 
outfMP = loadRaster(refRast,refBnds=refBnds,res=res,pad=T,refCrs=refCrs) #outflow directions
lcmMP  = loadRaster(lcmFN,refBnds=refBnds,res=res,pad=T,padWith=NA,refCrs=refCrs) #Land cover map (for lakes and roff)
roffMP = lcmToRoff(lcmMP,cvDF) #run-off coefs 

#Some roffMP cells maybe NA - turn them into "sea" (class 0).
outfMP[is.na(outfMP)]=0 
outfMP[is.na(roffMP)]=0 
roffMP[is.na(roffMP)]=0 

#outfMP follows esri convention for flows: 1,2,4...,128 where 1 is East, 2 is SE, 4 is S etc
#outfIndx uses a number that can be added to the 1D matrix index. This goes down the columns and then across the rows.
outfIndx = getOutfIndxMatrix(outfMP)
rm(outfMP)

#Remove all cells downstream from the edge of either the outflow grid or  no outflow direction
if(removeEdgeCats){
  incompMP = zeros(dim(outfIndx))
  incompMP[2,]=1;   incompMP[,2]=1
  incompMP[nrow(incompMP)-1,]=1;   incompMP[,ncol(incompMP)-1]=1
  incompMP[1:10,1:10]
  incompCats = sumGridOverUpstreamArea( incompMP ,outfIndx,maxCnt = Inf,verbose=T) #no-zero downstream of edge of outf
  outfIndx[incompCats>0] = 0
  excludeMat = (outfIndx==0)
  rm(incompCats)
}else{
  excludeMat=NULL
}


## load/calculate ccar
#calculate ccar .....
if(ccarFN==""){
  ccarMP = sumGridOverUpstreamArea( ones(dim(outfIndx)) ,outfIndx,maxCnt = Inf,verbose=T)
  #....and save it.
  if(ccarOutFN!="")exportTIFF( ccarMP ,refRast,"DEM_CCAR",ccarOutFN,dePad=T,excludeMat = excludeMat)
}
#or load it ...
if( ccarFN!=""){
  ccarMP = loadRaster(ccarFN,refBnds=refBnds,res=res,pad=T,refCrs=refCrs)
  ccarMP[is.na(ccarMP)] = 0
}




############# Part 2a ###############
## Create/save/load roffSum for base run-off grid
#Calculate roffSum ...
if( roffSumFN==""){
  print("CALCULATING accumulated base run-off grid. Consider specifying loadRoff0SumFN to load in future.")
  roffSum = sumGridOverUpstreamArea(roffMP,outfIndx,maxCnt = Inf,verbose=T) # < 5 seconds
  # ...and save roffSum ...
  if(outRoffSumFN!="") exportTIFF( roffSum ,refRast,"DEM_ROFF_SUM",outRoffSumFN,dePad=T,excludeMat = excludeMat)
}
# ... or load roffSum.
if( roffSumFN!=""){
  roffSum = loadRaster(roffSumFN,refBnds=refBnds,res=res,pad=T,refCrs=refCrs)
  roffSum[is.na(roffSum)] = 0
}


############# Part 2b ###############
## Create roffScnSum for scenario run-off grid roffScnFN. 
## Uses difference with base run-off grid so that this can be done quickly "on the fly"
if(lcmScnFN!=""){
  ########## Load grid :
  lcmScnMP = loadRaster(lcmScnFN,refBnds=refBnds,res=res,pad=T,padWith=NA,refCrs=refCrs)
  roffScnMP = lcmToRoff(lcmScnMP,cvDF)
  
  #Some roffScnMP cells are NA - turn them into "sea".
  roffScnMP[is.na(roffScnMP)]=0 
  
  ######## calculate roffScnSum :
  diff = roffScnMP - roffMP
  T0 = Sys.time()
  diffSum = sumGridOverUpstreamArea(diff,outfIndx,maxCnt = Inf,verbose=F)
  printp("Total time for roffScnSum is ", round(Sys.time() -T0, 4), "seconds" ) # < 0.5 seconds
  roffScnSum = diffSum + roffSum 
  
  ########  save roffScnSum (not neccessary)
  if(outRoffScnSumFN!="" )exportTIFF( roffScnSum ,refRast,"DEM_ROFF_SCN_SUM",outRoffScnSumFN,dePad=T,excludeMat = excludeMat)
}else{
  print("*no scenario run-off*")
}


############# Part 3a ###############
## Create/Load base FARL grid.
if(farlFN!=""){#load FARL.....
  farlMP = loadRaster(farlFN,refBnds=refBnds,res=res,pad=T,refCrs=refCrs)
  farlMP[is.na(farlMP)]=1
}else if(calcBaseFARL){#... or calculate FARL ...
  print("*calculate FARL*")
  lakeMP = matrix(F, nrow = nrow(lcmMP), ncol = ncol(lcmMP)) 
  for(ii in 1:length(lakeId)){
    lakeMP = lakeMP | (lcmMP == lakeId[ii])
  }
  lakeMP[is.na(lakeMP)| outfIndx==0] = FALSE
  lakesDF = mkGriddedLakeDFv2(lakeMP,ccarMP,outfIndx)
  farlMP = calcFarl(lakesDF,outfIndx,ccarMP,opt =farlOpt)
  if(farlOutFN!="" ) exportTIFF( farlMP ,refRast,farlOpt,farlOutFN,dePad=T,excludeMat = excludeMat)  #... and save FARL...
}else{#.. or just say no to FARL.
  print("*no FARL*")
  farlMP = ones(  dim(outfIndx) )
}


############# Part 3b ###############
## Create scenario FARL grid
# for useLcmScnForFarlScn=T all lakes (including those in base) are specified (by lakeId)
# for useLcmScnForFarlScn=F only new lakes are specified (by any value greater than zero)
if(  lakesScnFN!=""  | useLcmScnForFarlScn   ){
  
  ##create lakeScnMP - a matrix of lake cells
  #if using lcmScnMP for lake location (all lakes) :
  if(useLcmScnForFarlScn & lakesScnFN==""){
    ##lakeScnMP = !is.na(lcmScnMP) & (lcmScnMP == lakeId)
    lakeScnMP = matrix(F, nrow = nrow(lcmScnMP), ncol = ncol(lcmScnMP)) 
    for(ii in 1:length(lakeId)){
      lakeScnMP = lakeScnMP | (lcmScnMP == lakeId[ii])
    }
    lakeScnMP[is.na(lakeScnMP)| outfIndx==0] = FALSE
  }
  
  #if using lakesScnFN for lake location (just the new lakes) : 
  if( lakesScnFN!="" ){ 
    lakeScnMP = loadRaster(lakesScnFN,refBnds=refBnds,res=res,pad=T,refCrs=refCrs)
    lakeScnMP =  !is.na(lakeScnMP) & (lakeScnMP>0)
  }
  
  ##create lakesScnDF (a DF of lake outlet location, area, and ccar) and use for calc FARL (farlScnMP)
  T0 = Sys.time()
  lakesScnDF = mkGriddedLakeDFv2(lakeScnMP,ccarMP,outfIndx)
  farlScnMP = calcFarl(lakesScnDF,outfIndx,ccarMP,opt =farlOpt)
  if(lakesScnFN!="" )farlScnMP = farlScnMP*farlMP #only new lakes where included
  printp("Total time for farlScnMP is ", round(difftime(Sys.time(),T0,units="secs"), 2) ,"seconds") # 
 
  #and export
  if(farlScnOutFN!="" ) exportTIFF( farlScnMP ,refRast,farlOpt,farlScnOutFN,dePad=T,excludeMat = excludeMat)
  
}else{print("*no scenario FARL*")
  farlScnMP = farlMP
}


############# Part 4  ###############
##calculate flows and output

#Export base scenario flow. 
if( flowBaseFN  != "" ){
  flowBase = Cr* iRain * (res^2 * roffSum/100) * (farlMP^farlPower)#runoff was in percent (hence divide by 100)
  exportTIFF( flowBase ,refRast,"FLOW_BASE_M3S",flowBaseFN,dePad=T,excludeMat = excludeMat)
  
  ### print outlet flow ###
  iOut = which.max(ccarMP)
  printp("outlet flow (FARL ignored)  = ", round(Cr * res^2 * roffSum[iOut]/100 * iRain ,4 ) , " m^3/s")
  printp("outlet flow (FARL included) = ", round(Cr * res^2 * roffSum[iOut]/100 * iRain * (farlMP[iOut]^farlPower) ,4 ), " m^3/s")
}



#Export new scenario flow. 
if( flowScnFN  != "" ){
  flowScn = Cr*iRain * (res^2 * roffScnSum/100) * (farlScnMP^farlPower)
  exportTIFF( flowScn ,refRast,"FLOW_SCN_M3S",flowScnFN,dePad=T,excludeMat = excludeMat)
  
  ### print outlet flow ###
  iOut = which.max(ccarMP)
  printp("outlet flow (FARL ignored)  = ", round( Cr * res^2 * roffScnSum[iOut]/100 * iRain ,4 ) , " m^3/s")
  printp("outlet flow (FARL included) = ", round( Cr * res^2 * roffScnSum[iOut]/100 * iRain * (farlScnMP[iOut]^farlPower) ,4 ), " m^3/s")
}


#Export the ratio of the new scenario to the base flow. 
if( flowRatioFN  != "" ){
  flowRatio = ((farlScnMP/farlMP)^farlPower) * roffScnSum/roffSum - 1
  exportTIFF( flowRatio ,refRast,"FLOW_RATIO_SCN_BASE",flowRatioFN,dePad=T,excludeMat = excludeMat)
}

#take a csv with X and Y of bng coords and output various data
if(outletsInFN!="" & outletsOutFN!="" ) mkOutletDF(outletsInFN,outletsOutFN)

Sys.time()

