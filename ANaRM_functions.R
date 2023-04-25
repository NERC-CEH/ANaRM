##### Adapted Nature-based-solutions Rational Method’ (ANaRM) - supporting functions #####
# 
# Available for use under the terms of the included MIT licence. 
# Please consider citing:  
# J.D Miller, G. Vesuviano, J.R. Wallbank, D.H. Fletcher, L. Jones (2023). 
# "Hydrological assessment of urban Nature-Based Solutions for urban planning using Ecosystem Service toolkit applications"
# Landscape and Urban Planning 234  104737
#
##### John R Wallbank, UK Centre for Ecology & Hydrology (UKCEH), johwal@ceh.ac.uk, April 2023.

printp<-function(...)print(paste0(...))


#use a DF of classes and runoff values (percentages) to create a run-off grid from a Land Cover Map
lcmToRoff <-function(lcmMP,cvDF){
  roffMP = lcmMP
  if( any(! unique(lcmMP[!is.na(lcmMP)]) %in% cvDF[,"wrdclass"])) stop("Unrecognised classes in LCM. cvDF needs updating? ")
  for( i  in 1:nrow(cvDF)){
    roffMP[lcmMP == cvDF[i,"wrdclass"]] = cvDF[i,"cv"]  
  }
  return(roffMP)
}


#return vector of xmin,xmax,ymin,ymax from a raster file.
getBnds<-function(rasterFle){
  return(c(xmin(extent(rasterFle) ),xmax(extent(rasterFle)),ymin(extent(rasterFle) ),ymax(extent(rasterFle))))
}


#supply either file name or the raster itself to load it and return a matrix padded on the outer edge with "padWith"
loadRaster<-function(FN,refBnds=NA,res=10,pad=T,padWith=0,refCrs=NA){
  printp("LOADING grid :",FN)
  
  if(class(FN) == "character"){
    matrixPadded = raster(FN)
    if(!is.na(refCrs) & toString(crs(matrixPadded))!=refCrs ) cat("WARNING:\n CRS of",FN," is:\n",toString(crs(matrixPadded)),"\n refCrs is:\n",refCrs,"\n")
  } else if(class(FN) == "RasterLayer"){
    matrixPadded = FN
  }else{
    stop("Problem in loadRaster input") 
  }
  
  bnds    = getBnds(matrixPadded)
  matrixPadded = as.matrix(matrixPadded)
  
  #cropping to refBnds..
  if(all(!is.na(refBnds)))matrixPadded = crop2RefGrid(bnds,matrixPadded,refBnds,res=res,padWith=padWith)
  
  #pad with "sea" (zero value)  to prevent routing off the grid
  if(pad)matrixPadded = padarray(matrixPadded, c(1, 1), direction="both",padval=padWith)
  
  return(matrixPadded)
}


#strip the padding (from loadRaster) and export
exportTIFF <- function(datMat,rasterTemplate,nameStr,FN,dePad,excludeMat = NULL,excludeVal = NA,dType = 'FLT4S'){#
  printp("EXPORTING : ",FN)
  
  if(!is.null(excludeMat)) datMat[excludeMat]=excludeVal

  
  if(dePad){
    values(rasterTemplate) = datMat[2:(nrow(datMat)-1),2:(ncol(datMat)-1)]
  }else{
    values(rasterTemplate) = datMat
  }
  names(rasterTemplate)  = nameStr
  writeRaster(rasterTemplate,FN,overwrite=TRUE,datatype = dType)#default datatype = 'FLT4', v. minor inacuracies
  #Could use INT4S (but max of 2147483647 is only 40 times higher than max in current roffSum)
}


#crop dat (a Matrix) with bounds (bnds, [L,R,B,T]) to refBnds
crop2RefGrid<-function(bnds,dat,refBnds,res=10,padWith=0){
  dimDiff=c(1,-1,1,-1)*(refBnds - bnds)/res
  
  if( sum(abs( round(dimDiff)-dimDiff ))>res/10000 ){
    print("WARNING: grids off by non-integer multiple of the resolution!")
  }
  if( sum(abs( round(dimDiff)-dimDiff ) -c(0.5,0.5,0.5,0.5))>res/10000 ){
    print("WARNING: dat grid will be shifted up and right by half a grid.")
    dimDiff = dimDiff + 0.5*c(1,-1,1,-1)
  }
  dimDiff = round(dimDiff)
  
  cropBy = sapply(dimDiff, function(x)max(0,x))#[L,R,B,T]
  dat = dat[(cropBy[4]+1) : (nrow(dat)-cropBy[3]),(cropBy[1]+1):(ncol(dat)-cropBy[2])]
  padBy = sapply(dimDiff, function(x)max(0,-x))#[L,R,B,T]
  dat = padarray(dat, c(padBy[4], padBy[1]),direction="pre", padval=padWith) 
  dat = padarray(dat, c(padBy[3], padBy[2]),direction="post", padval=padWith) 
  return(dat)
}


#outfMP follows esri convention for flows: 1,2,4...,128 where 1 is East, 2 is SE, 4 is S etc
#outfIndx uses a number that can be added to the 1D matrix index. 
#This goes down the columns and then across the rows (hence add nrow(outfMP) to add 1 to row).
getOutfIndxMatrix<-function(outfMP){
  outfIndx = matrix(data = 0, nrow = nrow(outfMP), ncol = ncol(outfMP))
  outfIndx[outfMP==64] = -1
  outfIndx[outfMP==0 ] =  0 
  outfIndx[outfMP==4 ] =  1
  outfIndx[outfMP==128]=  nrow(outfMP) - 1
  outfIndx[outfMP==1 ] =  nrow(outfMP)
  outfIndx[outfMP==2 ] =  nrow(outfMP) + 1
  outfIndx[outfMP==32] = -nrow(outfMP) - 1
  outfIndx[outfMP==16] = -nrow(outfMP)
  outfIndx[outfMP==8 ] = -nrow(outfMP) + 1
  return(outfIndx)
}


#Takes a matrix of data (datGrid, e.g. roff coefs), and the OUTF directions (outfIndx, in index format) and returns matrix of datGrid summed over all upstream grid cells. Options: maxCnt = Inf (limit number of times each cell is routed), verbose=F (display timing/debuging info).
#Method: (1) Take a grid of ones (tmp0). (2) Make a list of all non-zero cells in that grid (nzIndx). (3) Use outf direction to route them all down stream, add them to the CCAR, and update nzIndx. (4) nzIndx will need some cleaning (could use nzIndx = which((tmp0>0)) here but would be slower). (5) loop until everything is routed.
#The benefit of this method is that the time needed for an itteration of the while loop is small when there are few cells to route. In Method 2 the time for each itteration becomes dominated by steps using the entire matrix.
sumGridOverUpstreamArea<-function(datGrid,outfIndx,maxCnt = Inf,verbose=F){
  LandCells = (outfIndx!=0)
  tmp0 = datGrid
  tmp0[!LandCells] = 0 #(S1)
  accum = zeros( dim(outfIndx) )
  nzIndx = which((tmp0!=0))#(S2) single integer indexing down columns then across rows.
  if(verbose)printp("initial non-zero cells ",length(nzIndx))
  #repeatedly propogate tmp0 downstream ...
  cnt=0 ; T0 = Sys.time();T1 = Sys.time()
  while( (cnt < maxCnt) & (length(nzIndx)>0) ){
    T1 = Sys.time()
    nzIndxNxt = nzIndx #nzIndx for the next loop
    for( ii in 1:length(nzIndx) ){#(S3)
      i  = nzIndx[ii]
      iN = i + outfIndx[i]#cell into which cell i flows
      tmp0[iN] = tmp0[iN] + tmp0[i]
      accum[i] = accum[i] + tmp0[i]
      tmp0[i] = 0
      nzIndxNxt[ii] = iN
    }
    T2 = Sys.time()
    #update indices to route (S4)
    nzIndxNxt = nzIndxNxt[tmp0[nzIndxNxt]!=0] #remove from cells to route if nothing to route
    nzIndxNxt = nzIndxNxt[LandCells[nzIndxNxt]]#only keep land cells
    nzIndxNxt = unique(nzIndxNxt)#remove from cells to route not unique
    nzIndx = nzIndxNxt
    T3 = Sys.time()
    cnt=cnt+1
    if(verbose & (cnt<10 | cnt%%100==0)){
      printp("* loop ",cnt," non-zero cells ",length(nzIndx)," *")
      printp("T2-T1 = ",T2-T1," T3-T2 = ",T3-T2)
    }
  }
  if(verbose)printp("T_total = ",as.numeric(T2-T0,  units= "secs"))# 8 seconds
  return(accum)
}


#For a lake point find all connected lakes points, setting lakesDatTmp=False to avoid double counting. 
#Was getConnectedOnlineLakeV2 in catchment descriptors python code. It checked if the river was online also (assumed online here).
getConnectedOnlineLake <- function(i,lakesDatTmp){
  if( !lakesDatTmp[i] ){
    printp("getConnectedOnlineLake : no lake at ",i)
    return(c())
  }
  nr=nrow(lakesDatTmp)
  nhd = c( -nr-1,  -nr, -nr+1,  -1, 1,   nr-1, nr , nr+1 )#indices of 8 neighbours to i in lakesDatTmp
  newi  = c(i) # lake cells found in last itteration of while loop
  lakei = c(i) #  all lake cells found so far
  #onLine=False #could change to T if river found in vicinity.
  while(length(newi) != 0){
    lakesDatTmp[newi]=F
    newiTmp = c()
    for(i in newi ){
      tmp = i+nhd[ lakesDatTmp[ i+nhd ] ]#indices of lake cells neighbouring cell i
      newiTmp = c(newiTmp,tmp)#add neighbourng lake to newiTmp
      lakesDatTmp[tmp]=F  #and set false in lakesDatTmp to avoid finding it again
      #if(!onLine) ..... then can check for river cells here. 
    }
    newi = newiTmp  
    lakei = c(lakei,newi)#call to unique removed (CHECK)
  }
  return(lakei)
}


#Input : lakei - the indices of lake cells (or a closed bdry) within ccarMP and outfMP. 
#Output: if there is a index with the highest ccar return it as the outlet. Otherwise use downstream ccar to tie-break. 
getOutletCoords<-function(lakei,ccarMP,outfIndx){ #(oldX,oldY,ccarDat,outfDat,flowDirectionDic):
  cLakeNow = ccarMP[lakei]
  oLakeNow = outfIndx[lakei]
  #check if the lake is all on bad data (e.g.) in sea
  if( all( cLakeNow<1 ) ){
    print("Lake is fully on non-land cells")
    return(-1) # ... then return an err flag 
  }
  maxcLakeNow= max(cLakeNow) #which.max can't deal with ties
  lMaxNow=lakei[ cLakeNow==maxcLakeNow ]
  if( length(lMaxNow) ==1 ){
    return(lMaxNow)
  }else{#tie break using downstream ccar
    ccarDownStearm = ccarMP[ lMaxNow+outfIndx[lMaxNow] ] 
    return(lMaxNow[which.max(ccarDownStearm)])
  }
}


#given a T/F grid for lakes this makes a DataFrame (lakesDF) of outlet, area, adn ccar [area units = number of cells]
#slightly faster than mkGriddedLakeDF
mkGriddedLakeDFv2 <- function(lakeMP,ccarMP,outfIndx){
  lakesOutletI= c()  #list of outlet-indices of lakes
  lakesCCAR   = c()  #list of catchment areas of lakes [units of res^2]
  lakesArea   = c()  #list of areas of lakes [units of res^2]
  allLakeCells = which(lakeMP)
  while( length(allLakeCells)!=0 ){
    printp("Lake cells remaining: ",length(allLakeCells))
    i=allLakeCells[1]
    #For a given lake cell (i) find all lake cells connected to it, ...
    lakei = getConnectedOnlineLake(i,lakeMP)
    #... remove them from allLakeCells
    allLakeCells = allLakeCells[!allLakeCells %in% lakei]
    # ....and record their outlet index, area, and catchment area.
    lakeout = getOutletCoords(lakei,ccarMP,outfIndx)
    lakesOutletI = c(lakesOutletI , lakeout  )
    lakesArea = c( lakesArea ,  length(lakei) )
    lakesCCAR = c( lakesCCAR , ifelse(lakeout>0, ccarMP[lakeout],1))
  }
  return( data.frame(outI = lakesOutletI, area = lakesArea, ccar = lakesCCAR ) )
}

#superceeded by mkGriddedLakeDFv2
mkGriddedLakeDF<-function(lakeMP,ccarMP,outfIndx){
  lakeMPTmp=lakeMP#Working array recording Lakes no yet processed (gets set to F)
  lakesOutletI=c() #list of outlet-indices of lakes
  lakesCCAR =c()   #list of catchment areas of lakes [units of res^2]
  lakesArea=c()    #list of areas of lakes [units of res^2]
  for( i in 1:length(lakeMP) ){
    if( i %% 100000 == 0 )printp("element ", i/1000,"k out of ", round(length(lakeMP)/1000),"k")
    if( lakeMPTmp[i] ){#if a lake cell is found ...
      #... then find all lake cells connected to it, ...
      lakei = getConnectedOnlineLake(i,lakeMPTmp)
      #...set them all to F to avoid finding again ...
      lakeMPTmp[lakei]=F #(would be beter to do this in getConnectedOnlineLake)
      # and record their outlet index, area, and catchemtn area.
      lakeout = getOutletCoords(lakei,ccarMP,outfIndx)
      lakesOutletI = c(lakesOutletI , lakeout  )
      lakesArea = c( lakesArea ,  length(lakei) )
      lakesCCAR = c( lakesCCAR , ifelse(lakeout>0, ccarMP[lakeout],1))
    }
  }
  return( data.frame(outI = lakesOutletI, area = lakesArea, ccar = lakesCCAR ) )
}


#for every lake identified in lakesDF move progressivley downstream calculating FARL.
#note : lakes' catchment area may need expanding as don't re-draw outf directions when lake is added
#hence the lake may be bigger than its catchment area! For the default option opt=="FARL" we insist
#the catchment is at least maxFF = 1.2 times larger than the lake otherwise the FARL effect 
#will not decay appropriately as we move downstream.
#For opt = "FARL" the standard definition of Flood Attenuation by Reservoits and Lakes is used from the
#UKCEH Flood Estimation handbook. For "FARL2-SL" or "FARL2-GL" a new proposed definition is used
#that should work better for small lakes. It is based on a simple model of how long a lake may delay 
#water for, but is not yet written up or published anywhere (USE AT OWN RISK) April 2022. The "FARL2-GL" 
#version should work for any size of lake (lake area compared to its catchment). "FARL2-SL" is a simpler 
#approximation to "FARL2-GL" for smaller lakes (relative to their catchments).
calcFarl<-function(lakesDF,outfIndx,ccarMP,opt = "FARL", maxFF = 1.2){
  ### OPTIONS ###
  maxCatFrac = ifelse(opt=="FARL",maxFF,1)  # for opt = "FARL" : there is a problem if area = AreaSubCat
  kappa = 30       # for opt = "FARL2-SL" or "FARL2-GL"
  mLL   = 1.5      # for opt = "FARL2-SL" or "FARL2-GL"
  ###############
  
  if( ! opt%in%c("FARL","FARL2-SL","FARL2-GL"))stop("bad option in calcFarl")
  
  #check if any lakes are bigger than their catchment (or similar in size for opt = "FARL")
  cond = maxCatFrac*lakesDF[,"area"] > lakesDF[,"ccar"]
  if(any(cond)){
    print("Increasing catchment area (ccar) of following lakes:")
    print(lakesDF[cond,])
    lakesDF[cond,"ccar"] = ceiling(maxCatFrac*lakesDF[cond,"area"])
  }
  
  LandCells = (outfIndx!=0)
  farlMP = ones(  dim(outfIndx) )
  
  for( n in 1:nrow(lakesDF) ){
    AreaSubCat = lakesDF[n,"ccar"]
    indx       = lakesDF[n,"outI"]
    area       = lakesDF[n,"area"]
    
    if(opt=="FARL"){ #ORIGINAL DEFINITION
      farlPt = 1.0 - sqrt( area/AreaSubCat )
      #if(n==36)printp("AreaSubCat=",AreaSubCat,"  indx=",indx,"  area=",area," farlPt=",farlPt," ccarMP",ccarMP[indx])
      while( LandCells[indx] ){
        farlMP[indx] = farlMP[indx]*(farlPt^(AreaSubCat/ccarMP[indx]) ) #could have FARL !=1 from previous lake
        indx = indx + outfIndx[indx]
      }
    }
    
    if(opt%in%c("FARL2-SL","FARL2-GL") ){ #PROPOSED DEFINITIONS - should work better for small lakes as may be added in urban environments  
      Ar = kappa*(area/AreaSubCat)
      if(opt == "FARL2-SL") TLdivT = Ar                     #FARL2 in "Small" Lake approximation (Alake<< Acatchment as most are)
      if(opt == "FARL2-GL") TLdivT = ( Ar^4 + (0.5^4)*Ar^(4*mLL) )^(0.25)      #FARL2 for General (area) Lake)
      farlPt = 1.0/( 1.0 + 1.0/TLdivT) #This is (1-f(i))
      while( LandCells[indx] ){
        farlMP[indx] = farlMP[indx]* (1.0 - farlPt*(AreaSubCat/ccarMP[indx]) ) 
        indx = indx + outfIndx[indx]
      }
    }
  }
  return(farlMP)
}


#A function to take a csv with X and Y of bng coords and output various data
mkOutletDF <- function(outletsInFN,outletsOutFN=""){
  hdr ="#iOutlet        = index of outlet on the grid
#ccar           = ccar (units of gridcells)
#flowBase       = flow_base_scenario  (includes effect of FARL)
#C_base = mean runoff ratio; base scenario (percent)
#C_Scn  = mean runoff ratio; new  scenario (percent)
#FARL_base = farl; base scenario
#FARL_Scn  = farl; new  scenario
#flowBaseNoFARL = flow_base_scenario (does not include effect of FARL)
#flowScn        = flow_new_scenario
#flowRatio      = flow_ratio_scenario_to_base (percentage)
"
  outletDF = read.csv(outletsInFN)
  outletDF[,c("iOutlet","ccar")] = -9999
  if(exists("roffSum"))   outletDF[,"C_base"]          = -9999
  if(exists("roffScnSum"))outletDF[,"C_Scn"]           = -9999
  
  if(exists("farlMP"))    outletDF[,"FARL_base"]       = -9999
  if(exists("farlScnMP")) outletDF[,"FARL_Scn"]        = -9999
  
  if(exists("flowBase"))  outletDF[,"flowBase"]        = -9999
  if(exists("flowBase"))  outletDF[,"flowBaseNoFARL"]  = -9999
  
  if(exists("flowScn"))   outletDF[,"flowScn"]         = -9999
  
  if(exists("flowRatio")) outletDF[,"flowRatio"]       = -9999
  
  nr = nrow(flowBase)
  
  for(n in 1:nrow(outletDF)){
    x=outletDF[n,"X"];y=outletDF[n,"Y"]
    iOutlet = floor((x - refBnds[1])/res + 1 ) *nr + ceiling((refBnds[4] - y)/res+1)
    
    outletDF[n,"iOutlet"]    = iOutlet
    if(ccarMP[iOutlet]<=0)next 
      
    outletDF[n,"ccar"]       = ccarMP[iOutlet]
    
    if(exists("roffSum"))   outletDF[n,"C_base"]         = roffSum[iOutlet]/ccarMP[iOutlet]
    if(exists("roffScnSum"))outletDF[n,"C_Scn"]          = roffScnSum[iOutlet]/ccarMP[iOutlet]
    
    if(exists("farlMP"))    outletDF[n,"FARL_base"]      = farlMP[iOutlet]
    if(exists("farlScnMP")) outletDF[n,"FARL_Scn"]       = farlScnMP[iOutlet]
    
    if(exists("flowBase"))  outletDF[n,"flowBase"]       = flowBase[iOutlet]
    if(exists("flowBase"))  outletDF[n,"flowBaseNoFARL"] = (flowBase/(farlMP^farlPower))[iOutlet]
    if(exists("flowScn"))   outletDF[n,"flowScn"]        = flowScn[iOutlet]
    if(exists("flowRatio")) outletDF[n,"flowRatio"]      = 100*flowRatio[iOutlet]
  }
  if(outletsOutFN!=""){
    write(hdr, file=outletsOutFN)
    suppressWarnings(write.table(outletDF, outletsOutFN,quote = F, sep = ",",row.names = F,append =T) )
  }
  return(outletDF)
}


