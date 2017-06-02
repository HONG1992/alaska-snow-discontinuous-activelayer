#Snow Seasonality and Discontinuous Active Layer Depth Processing
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(ggplot2)
library(psych)

#open folder where all files (shape and raster and data) are stored
    setwd("~/Desktop/modis_data/MODISSnowSeasonality")
    dataBaseDir <- c("~/Desktop/")
    lstDir <- paste(dataBaseDir,c("modis_data/MODISSnowSeasonality/"),sep="")

#Function for cropping raster using shape file
#source: https://stat.ethz.ch/pipermail/r-sig-geo/2013-July/018912.html
    homemadeClip<-function(raster,shape)
    {
      cropped<-crop(raster,shape)
      cells<-rasterize(shape,cropped)
      cropped*cells
    }

#in folder, fileNames are only the ones that end in .tif
    fileNames <- list.files(pattern = ".tif")

#create data frames
#first snow day full snow season
    fsSnow = data.frame("2001"=rep(0, 5), "2002"=rep(0,5), "2003"=rep(0,5),
                      "2004"=rep(0,5), "2005"=rep(0,5),"2006"=rep(0,5),
                      "2007"=rep(0,5), "2008"=rep(0,5), "2009"=rep(0,5),
                      "2010"=rep(0,5), "2011"=rep(0,5), "2012"=rep(0,5),
                      "2013"=rep(0,5),"2014"=rep(0,5), "2015"=rep(0,5), 
                    "2016"=rep(0,5));
#last snow day full snow season
    lsSnow = data.frame("2001"=rep(0, 5), "2002"=rep(0,5), "2003"=rep(0,5),
                      "2004"=rep(0,5), "2005"=rep(0,5), "2006"=rep(0,5),
                      "2007"=rep(0,5), "2008"=rep(0,5), "2009"=rep(0,5),
                      "2010"=rep(0,5), "2011"=rep(0,5), "2012"=rep(0,5),
                      "2013"=rep(0,5),"2014"=rep(0,5), "2015"=rep(0,5),
                    "2016"=rep(0,5));
#first snow day continuous snow season
    fsSnowCont = data.frame("2001"=rep(0, 5), "2002"=rep(0,5), "2003"=rep(0,5),
                          "2004"=rep(0,5), "2005"=rep(0,5), "2006"=rep(0,5),
                          "2007"=rep(0,5), "2008"=rep(0,5), "2009"=rep(0,5),
                          "2010"=rep(0,5), "2011"=rep(0,5), "2012"=rep(0,5),
                          "2013"=rep(0,5),"2014"=rep(0,5), "2015"=rep(0,5), 
                        "2016"=rep(0,5));
#last snow day continuous snow season
    lsSnowCont = data.frame("2001"=rep(0, 5), "2002"=rep(0,5), "2003"=rep(0,5),
                          "2004"=rep(0,5), "2005"=rep(0,5), "2006"=rep(0,5),
                          "2007"=rep(0,5), "2008"=rep(0,5), "2009"=rep(0,5),
                          "2010"=rep(0,5), "2011"=rep(0,5), "2012"=rep(0,5),
                          "2013"=rep(0,5), "2014"=rep(0,5), "2015"=rep(0,5),
                        "2016"=rep(0,5));

#change row names of each data frame (MEDIAN MEAN AND SD don't populate)
    row.names(fsSnow) <- c("min","max","median","mean","SD")
    row.names(lsSnow) <- c("min","max","median","mean","SD")
    row.names(fsSnowCont) <- c("min","max","median","mean","SD")
    row.names(lsSnowCont) <- c("min","max","median","mean","SD")

#change column name of each data frame
    colnames(fsSnow)<-c("2001","2002","2003","2004","2005","2006","2007","2008",
                    "2009","2010","2011","2012","2013","2014","2015","2016")
    colnames(lsSnow)<-c("2001","2002","2003","2004","2005","2006","2007","2008",
                    "2009","2010","2011","2012","2013","2014","2015","2016")
    colnames(fsSnowCont)<-c("2001","2002","2003","2004","2005","2006","2007","2008",
                        "2009","2010","2011","2012","2013","2014","2015","2016")
    colnames(lsSnowCont)<-c("2001","2002","2003","2004","2005","2006","2007","2008",
                        "2009","2010","2011","2012","2013","2014","2015","2016")


#iterate through all of the Snow seasonality files
    for (i in 1:length(fileNames))
    {
      savedFile = file.path('~/Desktop/modis_data/MODISSnowSeasonality/AKNew', fileNames[i]);
      show(savedFile)
      
      #if the file already exists, don't clip it
      if (file.exists(savedFile))
      {
        #Need these for later regardless
        AK_raster_stack <- stack(paste(lstDir, fileNames[i],sep=""))
        AKTest <- readOGR(dsn="./Alaska/Alaska.shp",layer ="Alaska")
        AK <-AKTest[AKTest$STATE_NAME == "Alaska",]
        AK_proj <- spTransform(AK, projection(AK_raster_stack))
        
        #cropped file
        AK_raster_stack_clip <- stack(paste(savedFile,sep=""))
      }
      #else, crop it and continue
      else
      {
        
        #read stack
        AK_raster_stack <- stack(paste(lstDir, fileNames[i],sep=""))
        
        #read Alaska shp file
        AKTest <- readOGR(dsn="./Alaska/Alaska.shp",layer ="Alaska")
        AK <-AKTest[AKTest$STATE_NAME == "Alaska",]
        
        #change coordinate reference system
        AK_proj <- spTransform(AK, projection(AK_raster_stack))
        
        #clip raster using shape file
        clip_snow_depth <- homemadeClip(raster = AK_raster_stack, shape = AK_proj)
        AK_raster_stack_clip <- clip_snow_depth
        
        #Write cropped raster to new file
        writeRaster(AK_raster_stack_clip, savedFile, options="INTERLEAVE=BAND", overwrite=TRUE)
      } 
      
      #populate data frame
      #first stack layer - first snow day full snow season
      extremaData <- AK_raster_stack_clip[[1]]
      extremaData[extremaData == 0] <- NA
      fsSnow[,i] = c("Min" = extremaData@data@min, "Max" = extremaData@data@max,  
                     median(extremaData@data@values,na.rm=TRUE), 
                     round(mean(extremaData@data@values,na.rm=TRUE)),
                     sd(extremaData@data@values,na.rm=TRUE))
      #second stack layer - last snow day full snow season
      extremaData <- AK_raster_stack_clip[[2]]
      extremaData[extremaData == 0] <- NA
      lsSnow[,i] = c("Min" = extremaData@data@min, "Max" = extremaData@data@max,  
                     median(extremaData@data@values,na.rm=TRUE), 
                     round(mean(extremaData@data@values,na.rm=TRUE)),
                     sd(extremaData@data@values,na.rm=TRUE))
      #third stack layer - first snow day continuous snow season
      extremaData <- AK_raster_stack_clip[[3]]
      extremaData[extremaData == 0] <- NA
      fsSnowCont[,i] = c("Min" = extremaData@data@min, "Max" = extremaData@data@max, 
                         median(extremaData@data@values,na.rm=TRUE), 
                         round(mean(extremaData@data@values,na.rm=TRUE)),
                         sd(extremaData@data@values,na.rm=TRUE))
      #fourth stack layer - last snow day continuous snow season
      extremaData <- AK_raster_stack_clip[[4]]
      extremaData[extremaData == 0] <- NA
      lsSnowCont[,i] = c("Min" = extremaData@data@min, "Max" = extremaData@data@max, 
                         median(extremaData@data@values,na.rm=TRUE), 
                         round(mean(extremaData@data@values,na.rm=TRUE)),
                         sd(extremaData@data@values,na.rm=TRUE))
      
      #check if it did anything
      show(i)
    }


#-------------------------------------------------------------------------------------------------------------
#Analysis of Permafrost Data
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

#open folder where permafrost data is
    setwd("~/Desktop/modis_data/calm_sites/CALM_SitesData_Cleaned")
    dataBaseDir <- c("~/Desktop/")
    listDir <- paste(dataBaseDir,c("modis_data/calm_sites/CALM_SitesData_Cleaned/"),sep="")
    boundDir <- paste(dataBaseDir,c("modis_data/calm_sites/CALM_SitesData_Cleaned/LCC_Arctic_Alaska_Yukon_copy/"),sep="")

#Load North Slope shape file (to exclude from AK later)
    NS <- readOGR(dsn=path.expand(paste(boundDir, c("LCC_Arctic_Alaska_Yukon.shp"),sep="")),layer ="LCC_Arctic_Alaska_Yukon")

#Load active layer file
    activeIn <- shapefile(paste(listDir,c("CALM_SitesData_Cleaned.shp"),sep=""),stringsAsFactors=FALSE)

#make sure they are in the same coordinate system
    ntSlopeProj <-projection(NS)
    activeAK <- spTransform(activeIn,ntSlopeProj)
    activeAKInd <- is.na(over(activeAK, as(NS, "SpatialPolygons"))) & !is.na(over(activeAK, as(AK_proj, "SpatialPolygons")))


#Create dataframe that has permafrost data for Alaska excluding NS
    var <- 25
    ALDataFrame = data.frame("SiteCode"=rep(0, var),"SiteName"=rep(0, var),
                           "Latitude"=rep(0, var),"Longitude"=rep(0, var),
                           "1990"=rep(0, var),"1991"=rep(0, var),"1992"=rep(0, var),
                           "1993"=rep(0, var),"1994"=rep(0, var), "1995"=rep(0, var),
                           "1996"=rep(0, var),"1997"=rep(0, var),"1998"=rep(0, var),
                           "1999"=rep(0, var),"2000"=rep(0, var),"2001"=rep(0, var),
                           "2002"=rep(0,var), "2003"=rep(0,var), "2004"=rep(0,var),
                           "2005"=rep(0,var), "2006"=rep(0,var), "2007"=rep(0,var),
                           "2008"=rep(0,var), "2009"=rep(0,var),"2010"=rep(0,var),
                           "2011"=rep(0,var), "2012"=rep(0,var), "2013"=rep(0,var),
                           "2014"=rep(0,var), "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(ALDataFrame)<-c("SiteCode","SiteName","Latitude","Longitude",
                           "1990","1991","1992","1993","1994","1995","1996",
                           "1997","1998","1999","2000","2001","2002","2003",
                           "2004","2005","2006","2007","2008","2009","2010",
                           "2011","2012","2013","2014","2015","2016")

#populate dataframe
    ALDataFrame[,1:4] <- activeAK@data[activeAKInd,1:4]
    ALDataFrame[,5:31] <- activeAK@data[activeAKInd,7:33]

#write out
    setwd("~/Desktop/modis_data/MODISSnowSeasonality/AKNew/dataFrames/")
    write.csv(ALDataFrame, file = "ALData.csv")


#-------------------------------------------------------------------------------------------------------------
#Finds snow extent where active layer point data is - extract and populate data frames
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

#Create data frames for first snow day (full) median snow depth at active layer depth locations
    fsSnowAtLoc = data.frame("2001"=rep(0, var), "2002"=rep(0,var), "2003"=rep(0,var),
                           "2004"=rep(0,var),"2005"=rep(0,var), "2006"=rep(0,var),
                           "2007"=rep(0,var), "2008"=rep(0,var),"2009"=rep(0,var),
                           "2010"=rep(0,var), "2011"=rep(0,var), "2012"=rep(0,var),
                         "2013"=rep(0,var),"2014"=rep(0,var), "2015"=rep(0,var),
                         "2016"=rep(0,var));
    colnames(fsSnowAtLoc) <- c("2001", "2002", "2003", "2004",
                           "2005", "2006", "2007", "2008",
                           "2009", "2010", "2011", "2012",
                           "2013","2014", "2015", "2016")

    row.names(fsSnowAtLoc) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])
#Create data frames for last snow day (full) median snow depth at active layer depth locations
    lsSnowAtLoc = data.frame("2001"=rep(0, var), "2002"=rep(0,var), "2003"=rep(0,var),
                           "2004"=rep(0,var),"2005"=rep(0,var), "2006"=rep(0,var),
                           "2007"=rep(0,var), "2008"=rep(0,var),"2009"=rep(0,var),
                           "2010"=rep(0,var), "2011"=rep(0,var), "2012"=rep(0,var),
                         "2013"=rep(0,var),"2014"=rep(0,var), "2015"=rep(0,var),
                         "2016"=rep(0,var));
    colnames(lsSnowAtLoc) <- c("2001", "2002", "2003", "2004",
                           "2005", "2006", "2007", "2008",
                           "2009", "2010", "2011", "2012",
                           "2013","2014", "2015", "2016")
    row.names(lsSnowAtLoc) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Create data frames for first snow day (continuous) median snow depth at active layer depth locations
    fsSnowContAtLoc = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                               "2003"=rep(0,var), "2004"=rep(0,var),
                             "2005"=rep(0,var), "2006"=rep(0,var),
                             "2007"=rep(0,var), "2008"=rep(0,var),
                             "2009"=rep(0,var), "2010"=rep(0,var),
                             "2011"=rep(0,var), "2012"=rep(0,var),
                             "2013"=rep(0,var),"2014"=rep(0,var),
                             "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(fsSnowContAtLoc) <- c("2001", "2002", "2003", "2004",
                               "2005", "2006", "2007", "2008",
                               "2009", "2010", "2011", "2012",
                               "2013","2014", "2015", "2016")
    row.names(fsSnowContAtLoc) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])
#Create data frames for last snow day (continuous) median snow depth at active layer depth locations
    lsSnowContAtLoc = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                               "2003"=rep(0,var), "2004"=rep(0,var),
                             "2005"=rep(0,var), "2006"=rep(0,var),
                             "2007"=rep(0,var), "2008"=rep(0,var),
                             "2009"=rep(0,var), "2010"=rep(0,var),
                             "2011"=rep(0,var), "2012"=rep(0,var),
                             "2013"=rep(0,var),"2014"=rep(0,var),
                             "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(lsSnowContAtLoc) <- c("2001", "2002", "2003", "2004",
                               "2005", "2006", "2007", "2008",
                               "2009", "2010", "2011", "2012",
                               "2013","2014", "2015", "2016")
    row.names(lsSnowContAtLoc) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Create data frames for duration of snow period (continuous) median snow depth at active layer depth locations
    durContSnowPer = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                              "2003"=rep(0,var), "2004"=rep(0,var),
                            "2005"=rep(0,var), "2006"=rep(0,var),
                            "2007"=rep(0,var), "2008"=rep(0,var),
                            "2009"=rep(0,var), "2010"=rep(0,var),
                            "2011"=rep(0,var), "2012"=rep(0,var),
                            "2013"=rep(0,var),"2014"=rep(0,var),
                            "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(durContSnowPer) <- c("2001", "2002", "2003", "2004",
                              "2005", "2006", "2007", "2008",
                              "2009", "2010", "2011", "2012",
                              "2013","2014", "2015", "2016")
    row.names(durContSnowPer) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Create data frames for the number of snow days at active layer depth locations
    numSDAtLoc = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                          "2003"=rep(0,var), "2004"=rep(0,var),
                        "2005"=rep(0,var), "2006"=rep(0,var),
                        "2007"=rep(0,var), "2008"=rep(0,var),
                        "2009"=rep(0,var), "2010"=rep(0,var),
                        "2011"=rep(0,var), "2012"=rep(0,var),
                        "2013"=rep(0,var),"2014"=rep(0,var),
                        "2015"=rep(0,var), "2016"=rep(0,var),
                        "SiteName"=rep(0, var));
    colnames(numSDAtLoc) <- c("2001", "2002", "2003", "2004",
                          "2005", "2006", "2007", "2008",
                          "2009", "2010", "2011", "2012",
                          "2013","2014", "2015", "2016","SiteName")
    row.names(numSDAtLoc) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Create data frames for total number of all days within CSS segments at active layer depth locations
    totNumDaysCSSAtLoc = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                                  "2003"=rep(0,var), "2004"=rep(0,var),
                                "2005"=rep(0,var), "2006"=rep(0,var),
                                "2007"=rep(0,var), "2008"=rep(0,var),
                                "2009"=rep(0,var), "2010"=rep(0,var),
                                "2011"=rep(0,var), "2012"=rep(0,var),
                                "2013"=rep(0,var),"2014"=rep(0,var),
                                "2015"=rep(0,var), "2016"=rep(0,var),
                                "SiteName"=rep(0, var));
    colnames(totNumDaysCSSAtLoc) <- c("2001", "2002", "2003", "2004",
                                  "2005", "2006", "2007", "2008",
                                  "2009", "2010", "2011", "2012",
                                  "2013","2014", "2015", "2016", "SiteName")
    row.names(totNumDaysCSSAtLoc) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])


    buf <- c(3)

#extracts snow data at active layer point locations
    for (i in 1:length(fileNames))
    {
      
      savedFile = file.path('~/Desktop/modis_data/MODISSnowSeasonality/AKNew', fileNames[i]);
      AK_raster_stack_clip <- stack(paste(savedFile,sep=""))
      
      #raster layers
      layer1 <- AK_raster_stack_clip[[1]] # first day of the full snow season (FSS start day)
      layer2 <- AK_raster_stack_clip[[2]] #last day of the full snow season (FSS end day)
      layer4 <- AK_raster_stack_clip[[4]] #first day of the longest CSS segment (CSS start day)
      layer5 <- AK_raster_stack_clip[[5]] #last day of the longest CSS segment (CSS end day)
      layer6 <- AK_raster_stack_clip[[6]] #longest_css_last_day-longest_css_first_day +1
      layer7 <- AK_raster_stack_clip[[7]] #the number of snow days
      layer12<- AK_raster_stack_clip[[12]] # total number of all days within CSS segments
      
      #if there are 0's set them to NA
      layer1[layer1 == 0] <- NA 
      layer2[layer2 == 0] <- NA
      layer4[layer4 == 0] <- NA
      layer5[layer5 == 0] <- NA
      layer6[layer6 == 0] <- NA
      layer7[layer7 == 0] <- NA
      layer12[layer12 == 0] <- NA
      
      #extract snow at active layer locations
      snowAtLocFFS <-extract(layer1,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      snowAtLocLFS <-extract(layer2,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      snowAtLocFCS <-extract(layer4,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      snowAtLocLCS <-extract(layer5,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      contSnowPerAtLoc <-extract(layer6,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      snowDaysAtLoc <-extract(layer7,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      totCSSAtLoc <- extract(layer12,activeAK@coords[activeAKInd,],buffer=buf,fun=median,na.rm=TRUE)
      
      #populate dataframes
      fsSnowAtLoc[,i] = c(snowAtLocFFS)
      lsSnowAtLoc[,i] = c(snowAtLocLFS)
      fsSnowContAtLoc[,i] = c(snowAtLocFCS)
      lsSnowContAtLoc[,i] = c(snowAtLocLCS)
      durContSnowPer[,i] = c(contSnowPerAtLoc)
      numSDAtLoc[,i] = c(snowDaysAtLoc)
      totNumDaysCSSAtLoc[,i] = c(totCSSAtLoc)
    }

    numSDAtLoc[17]<-ALDataFrame$SiteName
    totNumDaysCSSAtLoc[17]<-ALDataFrame$SiteName

#write out
    setwd("~/Desktop/modis_data/MODISSnowSeasonality/AKNew/dataFrames/")
    write.csv(fsSnowAtLoc, file = "fsSnowAtLoc.csv")
    write.csv(lsSnowAtLoc, file = "lsSnowAtLoc.csv")
    write.csv(fsSnowContAtLoc, file = "fsSnowContAtLoc.csv")
    write.csv(lsSnowContAtLoc, file = "lsSnowContAtLoc.csv")


#-------------------------------------------------------------------------------------------------------------
#Find dates of snow free period 2 ways
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

#Dataframe for snow free data (full season) -> Snow Year
    FULLSnowFreeSY = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                              "2003"=rep(0,var), "2004"=rep(0,var),
                            "2005"=rep(0,var), "2006"=rep(0,var),
                            "2007"=rep(0,var), "2008"=rep(0,var),
                            "2009"=rep(0,var), "2010"=rep(0,var),
                            "2011"=rep(0,var), "2012"=rep(0,var),
                            "2013"=rep(0,var),"2014"=rep(0,var),
                            "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(FULLSnowFreeSY) <- c("2001", "2002", "2003", "2004",
                              "2005", "2006", "2007", "2008",
                              "2009", "2010", "2011", "2012",
                              "2013","2014", "2015", "2016")
    row.names(FULLSnowFreeSY) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Dataframe for snow free data (continuous season) -> Snow Year
    CONTSnowFreeSY = data.frame("2001"=rep(0, var), "2002"=rep(0,var),
                              "2003"=rep(0,var), "2004"=rep(0,var),
                            "2005"=rep(0,var), "2006"=rep(0,var),
                            "2007"=rep(0,var), "2008"=rep(0,var),
                            "2009"=rep(0,var), "2010"=rep(0,var),
                            "2011"=rep(0,var), "2012"=rep(0,var),
                            "2013"=rep(0,var),"2014"=rep(0,var),
                            "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(CONTSnowFreeSY) <- c("2001", "2002", "2003", "2004",
                              "2005", "2006", "2007", "2008",
                              "2009", "2010", "2011", "2012",
                              "2013","2014", "2015", "2016")
    row.names(CONTSnowFreeSY) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Dataframe for snow free data (full season) -> Calendar Year
    FULLSnowFreeCY = data.frame("2002"=rep(0,var), "2003"=rep(0,var),
                              "2004"=rep(0,var), "2005"=rep(0,var),
                              "2006"=rep(0,var), "2007"=rep(0,var),
                              "2008"=rep(0,var), "2009"=rep(0,var),
                              "2010"=rep(0,var), "2011"=rep(0,var),
                              "2012"=rep(0,var), "2013"=rep(0,var),
                              "2014"=rep(0,var), "2015"=rep(0,var),
                              "2016"=rep(0,var));
    colnames(FULLSnowFreeCY) <- c("2002", "2003", "2004",
                              "2005", "2006", "2007", "2008",
                              "2009", "2010", "2011", "2012",
                              "2013","2014", "2015", "2016")
    row.names(FULLSnowFreeCY) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Dataframe for snow free data (continuous season) -> Calendar Year
    CONTSnowFreeCY = data.frame("2002"=rep(0,var), "2003"=rep(0,var),
                              "2004"=rep(0,var), "2005"=rep(0,var),
                              "2006"=rep(0,var), "2007"=rep(0,var),
                              "2008"=rep(0,var), "2009"=rep(0,var),
                              "2010"=rep(0,var), "2011"=rep(0,var),
                              "2012"=rep(0,var), "2013"=rep(0,var),
                              "2014"=rep(0,var), "2015"=rep(0,var),
                              "2016"=rep(0,var));
    colnames(CONTSnowFreeCY) <- c("2002", "2003", "2004",
                              "2005", "2006", "2007", "2008",
                              "2009", "2010", "2011", "2012",
                              "2013","2014", "2015", "2016")
    row.names(CONTSnowFreeCY) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])


#------------------------------------------------------------------------------------
#Merging the dataframes
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

#new dataframe with all of the data in one
    n <- 375
    totalsData <- data.frame("SiteName"=rep(0,n), "Year"=rep(0,n),
                           "ActiveLayer"=rep(0,n), "CFullSnow"=rep(0,n),
                           "CContSnow"=rep(0,n), "SFullSnow"=rep(0,n),
                           "SContSnow"=rep(0,n));
    colnames(totalsData) <- c("SiteName","Year","ActiveLayer","FULLSnowFreeCY",
                            "CONTSnowFreeCY", "FULLSnowFreeSY", "CONTSnowFreeSY")

#populate respectively
    for (i in 1:15)
    {
      a <- 1+(var*(i-1))
      b <- var+(var*(i-1))
      #Site Names
      totalsData[a:b,1] <- ALDataFrame[1:var,2]
      #Years
      totalsData[a:b,2] <- i+2001
      #Active Layer Depth
      totalsData[a:b,3] <- ALDataFrame[1:var, i+16]
      #C Full Snow
      totalsData[a:b,4] <- FULLSnowFreeCY[1:var,i]
      #C Cont Snow
      totalsData[a:b,5] <- CONTSnowFreeCY[1:var,i]
      #S Full Snow
      totalsData[a:b,6] <- FULLSnowFreeSY[1:var,i+1]
      #S Cont Snow
      totalsData[a:b,7] <- CONTSnowFreeSY[1:var,i+1]
    }
  
    totalsData[,3:7][totalsData[,3:7] == 0] <- NA

#------------------------------------------------------------------------------------
#First attempt at plotting
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

#This section doesn't really have a purpose since there isn't exactly any complete record

#Plot of largest complete set of active layer points (1996-2016, Barrow-Franklin Buff (1-10))
    #r=210
    #mostCompleteAL <- data.frame("SiteName"=rep(0,r), "Year"=rep(0,r),
    #                           "ActiveLayer"=rep(0,r));
    #colnames(mostCompleteAL) <- c("SiteName","Year","ActiveLayer")
    #for (i in 1:21)
    #{
    #  a <- 1+(10*(i-1))
    #  b <- 10+(10*(i-1))
    #  #Site Names
    #  mostCompleteAL[a:b,1] <- ALDataFrame[1:10,2]
    #  #Years
    #  mostCompleteAL[a:b,2] <- i+1995
    #  #Active Layer Depth
    #  mostCompleteAL[a:b,3] <- ALDataFrame[1:10, i+10]
    #}

#------------------------------------------------------------------------------------
#Freeze and Melt periods
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

#Freeze period - freeze period = start snow - 213
    CONTFreezeSY = data.frame("2002"=rep(0,var), "2003"=rep(0,var),
                            "2004"=rep(0,var), "2005"=rep(0,var),
                            "2006"=rep(0,var), "2007"=rep(0,var),
                            "2008"=rep(0,var), "2009"=rep(0,var),
                            "2010"=rep(0,var), "2011"=rep(0,var),
                            "2012"=rep(0,var), "2013"=rep(0,var),
                            "2014"=rep(0,var), "2015"=rep(0,var),
                            "2016"=rep(0,var));
    colnames(CONTFreezeSY) <- c("2002", "2003", "2004",
                            "2005", "2006", "2007", "2008",
                            "2009", "2010", "2011", "2012",
                            "2013","2014", "2015", "2016")
    row.names(CONTFreezeSY) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])

#Melt period - freeze period = 578 - last snow
    CONTMeltSY = data.frame("2002"=rep(0,var), "2003"=rep(0,var), "2004"=rep(0,var),
                        "2005"=rep(0,var), "2006"=rep(0,var), "2007"=rep(0,var),
                        "2008"=rep(0,var), "2009"=rep(0,var), "2010"=rep(0,var),
                        "2011"=rep(0,var), "2012"=rep(0,var), "2013"=rep(0,var),
                        "2014"=rep(0,var), "2015"=rep(0,var), "2016"=rep(0,var));
    colnames(CONTMeltSY) <- c("2002", "2003", "2004",
                          "2005", "2006", "2007", "2008",
                          "2009", "2010", "2011", "2012",
                          "2013","2014", "2015", "2016")
    row.names(CONTMeltSY) <-c(activeAK@data$Site_Name[activeAKInd == TRUE])


#Combine Active layer depth, freeze period, and melt period into one dataframe
    m <- 375
    meltFreezeAL <- data.frame("SiteName"=rep(0,m), "Year"=rep(0,m),
                             "ActiveLayer"=rep(0,m), "d_Melt"=rep(0,m),
                             "d_Freeze"=rep(0,m),"durationContSnowPer"=rep(0,m));
    colnames(meltFreezeAL) <- c("SiteName","Year","ActiveLayer","Melt", "Freeze", "durationContSnowPer") 

#for 2002-2016 populate the meltFreezeAL dataframe
    for (i in 1:15)
    {
      a <- 1+(var*(i-1))
      b <- var+(var*(i-1))
      #Site Names
      meltFreezeAL[a:b,1] <- ALDataFrame[1:var,2]
      #Years
      meltFreezeAL[a:b,2] <- i+2001
      #Active Layer Depth
      meltFreezeAL[a:b,3] <- ALDataFrame[1:var, i+16]
      #Melt
      meltFreezeAL[a:b,4] <- CONTMeltSY[1:var,i]
      #Freeze
      meltFreezeAL[a:b,5] <- CONTFreezeSY[1:var,i]
      #durationContSnowPer
      meltFreezeAL[a:b,6] <- durContSnowPer[1:var,i+1]
    }
    meltFreezeAL[,3:6][meltFreezeAL[,3:6] == 0] <- NA

#Combine RED (look at original data online) --> All of them are red hence doesn't apply
#Active layer depth, freeze period, and melt period into one dataframe
    q <- 375
    temp <- data.frame("SiteName"=rep(0,q), "Year"=rep(0,q), "ActiveLayer"=rep(0,q),
                   "d_Melt"=rep(0,q), "d_Freeze"=rep(0,q),
                   "durationContSnowPer"=rep(0,q));
    colnames(temp) <- c("SiteName","Year","ActiveLayer","Melt", "Freeze",
                        "durationContSnowPer") 

#for 2002-2016 populate Dataframe
    for (i in 1:15)
    {
      #look only at RED (see original data) data (also see temp for the dataframe)
      a <- 1+(25*(i-1))
      b <- 25+(25*(i-1))
      #Site Names
      temp[a:b,1] <- ALDataFrame[1:25,2]
      #Years
      temp[a:b,2] <- i+2001
      #Active Layer Depth
      temp[a:b,3] <- ALDataFrame[1:25, i+16]
      #Melt
      temp[a:b,4] <- CONTMeltSY[1:25,i]
      #Freeze
      temp[a:b,5] <- CONTFreezeSY[1:25,i]
      #durationContSnowPer
      temp[a:b,6] <- durContSnowPer[1:25,i+1]
    }

#new dataframe looking at the data collected using the 1000 method (only 2 points)
    q <- 30
    thousandDF <- data.frame("SiteName"=rep(0,q), "Year"=rep(0,q),
                           "ActiveLayer"=rep(0,q), "durationContSnowPer"=rep(0,q))
#1000 method (The red ones in the data using the 1000 method -- only 2 points)
    tempvec <- c(24,25) #the position of the 1000 method ones in the temp DF

    for (i in 1:15)
    {
      a <- 1+(2*(i-1))
      b <- 2+(2*(i-1))
      c<- tempvec+(25*i-25)
      thousandDF[a:b,1] <- temp[c,1]
      thousandDF[a:b,2] <- i+2001
      thousandDF[a:b,3] <- temp[c, 3]
      thousandDF[a:b,4] <- temp[c,6] 
    }
