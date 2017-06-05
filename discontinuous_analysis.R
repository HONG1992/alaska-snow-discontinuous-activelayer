#This script contains the analysis of snow seasonality descriptors and permafrost
#active layer depth (discontinuous zone) data using the dataframes created in the
#file snow_seasonality_proc.R. See snow_seasonality_proc.R and/or comments for
#more info.


#Snow Duration and Active Layer Depth (Discontinuous Zone) Analysis
#Earth Lab - Project Permafrost
#By: Ksenia Lepikhina and Jeffery Thompson (mentor)
#Copy Right

library(sp)
library(rgdal)
library(raster)
library(rgeos)
library(ggplot2)
library(psych)
library(dplyr)
library(tidyr)
library(boot)
#-------------------------------------------------------------------------------------------------------------
#Plot north slope data

#label plots and plot
#would need to create a for loop in order to see each year. Else it'll just plot
#whatever file was the last one used.
    labels <- c("first_snow_day", "last_snow_day", "fss_range",
                "longest_css_first_day", "longest_css_last_day",
                "longest_css_day_range", "snow_days", "no_snow_days",
                "css_segment_num", "mflag", "cloud_days", "tot_css_days")
    names(AK_raster_stack_clip) <- labels[1:12]
    plot(AK_raster_stack_clip)
#-------------------------------------------------------------------------------------------------------------
#Analysis of Permafrost Data

#Plot the North Slope -**WARNING** takes forever
    plot(coordinates(activeAK),type="n")
    plot(NS,border="blue",add=TRUE)
    plot(AK_proj,border="green",add=TRUE)

#Plot active layer sites inside and outside shape file of North Slope 
    points(activeAK[!activeAKInd, ], pch=1, col="gray") #outside
    points(activeAK[activeAKInd, ], pch=16, col="red") #inside
#-------------------------------------------------------------------------------------------------------------
#Find dates of snow free period (full and continuous) 2 ways

#Snow Year Analysis (FULL snow season)
    dSFP <- vector("numeric", 25L) #Vector to populate
#data goes from 2001-2016
    for (i in 1:16)
    {
      #for 42 locations
      for (j in 1:25)
      {
        #dSFP = 365-(End –Start)
        k <- 365 -(lsSnowAtLoc[j,i]-fsSnowAtLoc[j,i])
        dSFP[j] <- k
      }
      #populate data frame
      FULLSnowFreeSY[,i] = c(dSFP)
    }

#Snow Year Analysis (CONTINUOUS snow season) (same as above w/ continuous instead of full)
    dSFP2 <- vector("numeric", 25L)
    for (i in 1:16)
    {
      for (j in 1:25)
      {
        k <- 365 -(lsSnowContAtLoc[j,i]-fsSnowContAtLoc[j,i])
        dSFP2[j] <- k
      }
      CONTSnowFreeSY[,i] = c(dSFP2)
    }

#Calendar Year Analysis FULL (lose a year in the math so data goes from 2002-2016)
    dSFP3 <- vector("numeric", 25L)
    for (i in 1:15)
    {
      for (j in 1:25)
      {
        #check if NA
        if (is.na(fsSnowAtLoc[j,i+1]) || is.na(lsSnowAtLoc[j,i]))
        {
          dSFP3[j] <- NA
        }
        # start(y+1) < 365 & end(y) > 365
        else if (fsSnowAtLoc[j,i+1] < 365 && lsSnowAtLoc[j,i] > 365)
        {
          #dSFP = start(y+1) – (end(y) -365)
          k <- fsSnowAtLoc[j,i+1] - (lsSnowAtLoc[j,i] -365)
          dSFP3[j] <- k
        }
        #start(y+1) < 365 & end(y) < 365
        else if (fsSnowAtLoc[j,i+1] < 365 && lsSnowAtLoc[j,i] < 365)
        {
          #dSFP = start(y+1)
          dSFP3[j] <- fsSnowAtLoc[j,i+1] 
        }
        #start(y+1) > 365 & end(y) > 365
        else if (fsSnowAtLoc[j,i+1] > 365 && lsSnowAtLoc[j,i] > 365)
        {
          #dSFP = 365 – (end(y) -365)
          k <- 365 - (lsSnowAtLoc[j,i] -365)
          dSFP3[j] <- k
        }
        #start(y+1) > 365 & end(y) < 365
        else if (fsSnowAtLoc[j,i+1] > 365 && lsSnowAtLoc[j,i] < 365)
        {
          #dSFP = 365
          dSFP3[j] <- 365
        }
      }
      #populate dataframe
      FULLSnowFreeCY[,i] = c(dSFP3)
    }

#Calendar Year Analysis CONTINUOUS (same format as above except continuous instead of full)
    dSFP4 <- vector("numeric", 25L)
    for (i in 1:15)
    {
      for (j in 1:25)
      {
        if (is.na(fsSnowContAtLoc[j,i+1]) || is.na(lsSnowContAtLoc[j,i]))
        {
          dSFP4[j] <- NA
        }
        else if (fsSnowContAtLoc[j,i+1] < 365 && lsSnowContAtLoc[j,i] > 365)
        {
          k <- fsSnowContAtLoc[j,i+1] - (lsSnowContAtLoc[j,i] -365)
          dSFP4[j] <- k
        }
        else if (fsSnowContAtLoc[j,i+1] < 365 && lsSnowContAtLoc[j,i] < 365)
        {
          dSFP4[j] <- fsSnowContAtLoc[j,i+1] 
        }
        else if (fsSnowContAtLoc[j,i+1] > 365 && lsSnowContAtLoc[j,i] > 365)
        {
          k <- 365 - (lsSnowContAtLoc[j,i] -365)
          dSFP4[j] <- k
        }
        else if (fsSnowAtLoc[j,i+1] > 365 && lsSnowContAtLoc[j,i] < 365)
        {
          dSFP4[j] <- 365
        }
      }
      CONTSnowFreeCY[,i] = c(dSFP4)
    }

#------------------------------------------------------------------------------------
#Merging the dataframes

#------------------------------------------------------------------------------------
#First attempt at plotting

#test box plots
#AL depth vs Years
    boxplot(ActiveLayer ~ Year, data=totalsData,main="Active Layer Depths by Year")
#Continuous snow free period (snow year) vs year
    boxplot(CONTSnowFreeSY ~ Year, data=totalsData,main="Continous Snow Free Period - Snow Year")
#Continuous snow free period (calendar year) vs year
    boxplot(CONTSnowFreeCY ~ Year, data=totalsData,main="Continous Snow Free Period - Cal. Year")
#Full snow free period (snow year) vs year
    boxplot(FULLSnowFreeSY ~ Year, data=totalsData,main="Full Snow Free Period - Snow Year")
#Full snow free period (calendar year) vs year
    boxplot(FULLSnowFreeCY ~ Year, data=totalsData,main="Full Snow Free Period - Cal. Year")

#Active layer vs Continuous Snow Free Period (snow year)
    contSnowFreeSnowYearInd <- totalsData$CONTSnowFreeSY >0 #& totalsData$CONTSnowFreeSY <250
    plot(totalsData[contSnowFreeSnowYearInd,7],totalsData[contSnowFreeSnowYearInd,3],main="Continuous Snow Free Period (Snow Year) vs Active Layer Depth", xlab = "Continuous Snow Free Period in days (Snow Year)", ylab = "Active Layer Depth (cm)")
    contSnowFreeSYReg <- lm(totalsData[contSnowFreeSnowYearInd,3] ~totalsData[contSnowFreeSnowYearInd,7])
    abline(contSnowFreeSYReg)
    summary(contSnowFreeSYReg)

#same as above but with more stats
    contSnowFreeSYReg <- lm(ActiveLayer ~ CONTSnowFreeSY, data=totalsData)
    summary(contSnowFreeSYReg)
    plot(contSnowFreeSYReg)

#same as above but with better labels
    #plot(totalsData$CONTSnowFreeSY,totalsData$ActiveLayer)
    #abline(contSnowFreeSYReg)


#checking if there is any correlation between years - FULL Snow Year
    fullSnowFreeSnowYearInd <- totalsData$FULLSnowFreeSY<300
    par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE) #sets border size
    plot(totalsData[fullSnowFreeSnowYearInd,6],totalsData[fullSnowFreeSnowYearInd,3],
         pch = c(0,1,2,5,6,15,16,17,18,19,20,21,22,23,24),bg = c(
           "red","red","red","red"), main = "Full Snow Year, AL Depth", xlab = "Full Snow Free Period in days (Snow Year)", ylab = "Active Layer Depth (cm)")
    legend('topright', inset=c(-0.2,0),names(FULLSnowFreeSY[2:16]), 
           pch=c(0,1,2,5,6,15,16,17,18,19,20,21,22,23,24), pt.bg=c("red","red","red","red"), bty='n', cex=.75)
    fullSnowFreeSYReg <- lm(totalsData[fullSnowFreeSnowYearInd,3] ~ totalsData[fullSnowFreeSnowYearInd,6])
    abline(fullSnowFreeSYReg)
    summary(fullSnowFreeSYReg)

#Looking at individual years -FULL Snow Year
    fullSnowFreeSnowYearInd <- totalsData$FULLSnowFreeSY<300 #& 100<totalsData$FULLSnowFreeSY[1:42,1]
    for (i in 2:16)
    {
      #set each row to go from a:b each iteration
      a <- 1+(25*(i-2))
      b <- 25+(25*(i-2))
      par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(totalsData[a:b,6],totalsData[a:b,3])
      legend('topright', inset=c(-0.2,0),names(FULLSnowFreeSY[i]), bty='n', cex=.75)
      fullSnowFreeSYReg <- lm(totalsData[a:b,3] ~ totalsData[a:b,6])
      abline(fullSnowFreeSYReg)
      show(summary(fullSnowFreeSYReg))
    } #not great R^2 for any

#checking if there is any correlation between locations - FULL Snow Year
    fullSnowFreeSnowYearInd <- totalsData$FULLSnowFreeSY<300 #& 67<totalsData$FULLSnowFreeSY
    par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
    myshapes <- c(21,22)#,21,22,23,24,21,22,23,24,
    #21,22,23,24,21,22,23,24,21,22,23,24,
    #21,22,23,24,21,22,23,24,21,22,23,24,
    #21,22,23,24,21,22)
    mycolors <- c("red","red",
                  "green","green",
                  "orange","orange",
                  "blue","blue",
                  "cyan","cyan",
                  "black","black", "black",
                  "grey","grey",
                  "yellow","yellow",
                  "pink","pink",
                  "darkgreen","darkgreen",
                  "white","white")
    plot(totalsData[fullSnowFreeSnowYearInd,6], totalsData[fullSnowFreeSnowYearInd,3],
         pch=myshapes,bg=mycolors, main = "Full Snow Year, Locations")
    legend(275,225, inset=c(-0.1,0),totalsData$SiteName[1:25], 
           pch= myshapes, pt.bg = mycolors,bty='n', cex=.55, ncol=1)

#Looking at individual locations - FULL Snow Year
    fullSnowFreeSnowYearInd <- totalsData$FULLSnowFreeSY<300 #& 100<totalsData$FULLSnowFreeSY[1:42,1]
#plot each location on a separate plot
    for (i in 1:25)
    {
      locData <- c(i,i+25,i+(2*25), i+(3*25),i+(4*25),i+(5*25),i+(6*25),i+(7*25),
                   i+(8*25),i+(9*25),i+(10*25),i+(11*25),i+(12*25),i+(13*25),i+(14*25))
      #a <- 1+(42*(i-1))
      #b <- 42+(42*(i-1))
      par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
      if (all(is.na(totalsData[locData,3])))
      {
        next
      }
      plot(totalsData[locData,6],totalsData[locData,3], main = "Individual Locations", xlab = "Full Snow Free Period (Snow Year)", ylab = "Active Layer Depth" )
      legend('topright', inset=c(-0.2,0),ALDataFrame$SiteName[i], bty='n', cex=.75)
      fullSnowFreeSYReg <- lm(totalsData[locData,3] ~ totalsData[locData,6])
      abline(fullSnowFreeSYReg)
      show(summary(fullSnowFreeSYReg))
    }

#Plot of largest complete set of active layer points (1996-2016, Barrow-Franklin Buff (1-10))
#(Not quite as applicable for the discontinuous dataset)    
    #pred <- pre.frame$mostCompleteAL...2.
    
    completeALReg <- lm(mostCompleteAL[,3] ~ mostCompleteAL[,2])
    #pre.frame <- data.frame(mostCompleteAL[,2])
    #pp<- predict(completeALReg, int = "p", newdata = pre.frame)
    #pc<- predict(completeALReg, int = "c", newdata = pre.frame)
    #plot(mostCompleteAL[,2], mostCompleteAL[,3], main = "Active Layer Depth
    #     vs Years", xlab = "Years(1996-2016)", ylab = "Active Layer Depth (cm)", 
    #     ylim=range(mostCompleteAL[,3], pp, na.rm = T))
    #matlines(pred, pc, lty = c(1,2,2), col = "blue") #confidence interval
    #matlines(pred,pp, lty = c(1,3,3), col = "black") #prediction interval
    #show(summary(completeALReg))
    #intcp <- coef(completeALReg)[1] 
    #slp <-  coef(completeALReg)[2] 
    
#------------------------------------------------------------------------------------
#Freeze and Melt periods
    dSFP5 <- vector("numeric", 25L) #new vector
    #2002-2016
    for (i in 1:15)
    {
      for (j in 1:25)
      {
        if (is.na(fsSnowContAtLoc[j,i]) || is.na(lsSnowContAtLoc[j,i]))
        {
          dSFP5[j] <- NA
        }
        else
        {
          #check
          show(j)
          show(i)
          #calculate the freeze period
          d_freeze <- fsSnowContAtLoc[j,i+1]-213
          dSFP5[j] <- d_freeze
          #check
          show(d_freeze)
        }
      }
      #place in Dataframe
      CONTFreezeSY[,i] = c(dSFP5) 
    }
    
    #almost the same as above
    dSFP6 <- vector("numeric", 25L)
    for (i in 1:15)
    {
      for (j in 1:25)
      {
        if (is.na(fsSnowContAtLoc[j,i]) || is.na(lsSnowContAtLoc[j,i]))
        {
          dSFP6[j] <- NA
        }
        else
        {
          #calculate the melt period 
          d_melt <- 578 - lsSnowContAtLoc[j,i+1]
          dSFP6[j] <- d_melt
        }
      }
      CONTMeltSY[,i] = c(dSFP6) 
    }

#checking if there is any correlation between years - AL vs Freeze
    freeze <- meltFreezeAL$Freeze <150 #& 100<totalsData$FULLSnowFreeSY[1:42,1]
    par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
    plot(meltFreezeAL[freeze,5],meltFreezeAL[freeze,3],
         pch = c(0,1,2,5,6,15,16,17,18,19,20,21,22,23,24),bg = c(
           "red","red","red","red"), main = "AL vs Freeze", xlab = "Duration of 'Freeze' Period", ylab = "Active Layer Depth (cm)")
    legend('topright', inset=c(-0.2,0),names(FULLSnowFreeSY[2:16]),
           pch=c(0,1,2,5,6,15,16,17,18,19,20,21,22,23,24), pt.bg=c("red","red","red","red"), bty='n', cex=.75)
    freezeReg <- lm(meltFreezeAL[freeze,3] ~ (meltFreezeAL[freeze,4]+meltFreezeAL[freeze,5] +meltFreezeAL[freeze,6]))
    abline(freezeReg)
    summary(freezeReg) 

#Looking at individual years -AL and duration of continuous snow period <-(meltFreezeLineIndiv)
    for (i in 2:16)
    {
      a <- 1+(25*(i-2))
      b <- 25+(25*(i-2))
      par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
      plot(meltFreezeAL[a:b,6],meltFreezeAL[a:b,3], xlab = "Number of Snow Days", ylab = "Active Layer Depth")
      legend('topright', inset=c(-0.2,0),names(FULLSnowFreeSY[i]), bty='n', cex=.75)
      meltFreezeLineIndiv <- lm(meltFreezeAL[a:b,3] ~ meltFreezeAL[a:b,6])#meltFreezeAL[a:b,4]+ meltFreezeAL[a:b,5]+meltFreezeAL[a:b,6])
      abline(meltFreezeLineIndiv)
      show(summary(meltFreezeLineIndiv))
    } 
    
    
#Looking at individual locations -AL and duration of continuous snow period <-(meltFreezeLineIndiv)
    for (i in 1:25)
    {
      locData2 <- c(i,i+25,i+(2*25), i+(3*25),i+(4*25),i+(5*25),i+(6*25),i+(7*25),
                   i+(8*25),i+(9*25),i+(10*25),i+(11*25),i+(12*25),i+(13*25),i+(14*25))
      par(mar=c(4.1, 4.1, 4.1, 8.1), xpd=TRUE)
      if (all(is.na(meltFreezeAL[locData2,3])))
      {
        next
      }
      if (all(is.na(meltFreezeAL[locData2,6])))
      {
        next
      }
      if (all(!is.na(meltFreezeAL[locData2,6])))
      {
        if (all(!is.na(meltFreezeAL[locData2,3])))
        {
        plot(meltFreezeAL[locData2,6],meltFreezeAL[locData2,3])
        legend('topright', inset=c(-0.2,0),meltFreezeAL$SiteName[i], bty='n', cex=.75)
        meltFreezeLineIndiv <- lm(meltFreezeAL[locData2,3] ~ meltFreezeAL[locData2,6])#meltFreezeAL[a:b,4]+ meltFreezeAL[a:b,5]+meltFreezeAL[a:b,6])
        abline(meltFreezeLineIndiv)
        #show(summary(meltFreezeLineIndiv))
        }
      }
    }
    
   

#AL depth vs Years
    boxplot(ActiveLayer ~ Year, data=meltFreezeAL,main="Active Layer Depths by Year")
#Melt vs Years
    boxplot(Melt ~ Year, data=meltFreezeAL,main="Melt - Year")
#Freeze vs Years
    boxplot(Freeze ~ Year, data=meltFreezeAL,main="Freeze - Year")
#duration of continuous snow period vs Years
    boxplot(durationContSnowPer ~ Year, data=meltFreezeAL,main="Duration of Cont Snow Per - Year")

#look only at RED (see original data) (Doesn't matter they are all red)
    plot(temp[,6], temp[,3], main = "ActiveLayer vs Duration of Continuous
         Snow Period", xlab = "Duration of Continuous Snow Period", ylab = "Active Layer Depth")
    tempReg <- lm(temp[,3] ~ (temp[,6]))
    abline(tempReg, method="spearman")
    summary(tempReg)
    corr.test(temp[3],temp[6], method="spearman")

#1000 method (The red ones in the data using the 1000 method)
    plot(thousandDF[,4], thousandDF[,3], main = "1000x1000m grid:
          Active Layer Depth vs Duration of Continuous Snow Period", xlab = 
           "Duration of Continuous Snow Period", ylab = "Active Layer Depth")
    thousandReg <- lm(thousandDF[,3] ~ (thousandDF[,4]))
    abline(thousandReg, method= "spearman")
    summary(thousandReg)
    corr.test(thousandDF[3],thousandDF[4],method= "spearman") #NOT a very good correlation AT ALL using the 1000 method

 
    
    vars <- c("1990", "1991","1992","1993","1994","1995","1996","1997","1998","1999","2000","2001",
              "2002","2003","2004","2005","2006","2007","2008","2009","2010","2011","2012","2013",
              "2014","2015","2016")   
#------------------------------------------------------------------------------------
# make a pairs plot for numeric data

#active layer temporal
    ALDataFrame %>%
      select(starts_with("1"), starts_with("2"))  %>% #one_of(vars)) %>% #
      pairs
    
      
#active layer spatial    
    ALDataFrame %>%
      gather(year, value, -SiteCode, -SiteName, -Latitude, -Longitude) %>%
      ggplot(aes(as.numeric(year), value, color = SiteName)) + 
      geom_line() +
      labs(title = "Active Layer Spatial")
    
#longest continuous snow season day range (layer 7) temporal
    numSDAtLoc %>%
      select(starts_with("2")) %>%
      pairs
    
#longest continuous snow season day range (layer 7) spatial    
    numSDAtLoc %>%
      gather(year, value, -SiteName) %>%
      ggplot(aes(as.numeric(year), value, color = SiteName)) + 
      geom_line()+
      labs(title = "Number of Snow Days Temporal")
    
#total continuous snow season days (layer 12) temporal
    totNumDaysCSSAtLoc %>%
      select(starts_with("2")) %>%
      pairs

#total continuous snow season days (layer 12) spatial    
    totNumDaysCSSAtLoc %>%
      gather(year, value, -SiteName) %>%
      ggplot(aes(as.numeric(year), value, color = SiteName)) + 
      geom_line() +
      labs(title = "Number of days with CSS Segments Temporal")
    
#------------------------------------------------------------------------------------
#test bootstrap

#continuous snow free snow year and year
    bootContSYvsYear <- boot(data=totalsData, statistic = rsq, 
                    R=1000, formula=CONTSnowFreeSY ~ Year)
    summary(bootContSYvsYear)
    plot(bootContSYvsYear, rm.na=TRUE)
    boot.ci(bootContSYvsYear, type="bca")
    
#AL and year
    bootALvsYear <- boot(data=totalsData, statistic=rsq, 
                    R=1000, formula= ActiveLayer ~ Year)
    summary(bootALvsYear)
    plot(bootALvsYear, rm.na=TRUE)
    boot.ci(bootALvsYear, type="bca")

    
    