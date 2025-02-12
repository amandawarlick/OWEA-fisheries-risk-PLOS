### Functions

# =======================================
# = Make dt into SpatialPointsDataFrame =
# =======================================
### If want to do all sppcode then make sp=="all"

library(sf)

make_spdf <- function(dt, sp = sp, com_por = com_por, yr = yr, save = T, output = F, proj.type, lon.var = "declon", lat.var = "declat"){
  print("step1")
  dt2 <- as.data.frame(dt)
  print("step2")
  
  # convert to sf object, the CRS is 4326 (WGS84 Lat/Lon)
  locs <- dt2 %>% 
    st_as_sf(coords = c('Best_Long','Best_Lat'), crs = 4326) %>% # may want to change 'Best_Long','Best_Lat' to lon.var, lat.var
    st_transform(crs = proj.type) # change to proj.type

  print("step5")
  ### Save results as shape file
  if(save == T){
    name <- paste("spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", yr, sep="")
	
    # Write to shape file
    wd <- getwd()
    st_write(locs, paste0(wd,'/Confidential/footprint_files/', sp, '/PVC/', name,'.shp'), append=FALSE)
    
    name2 <- paste("csv_","spp_", sp, "_com_port_", unique(com_por), "_year_", yr, ".csv", sep="")
    write.csv(locs, name2)
  }
  if(output==T){
    return(locs)
  }
  
}

# ========================================
# = Calculate kde and output into raster =
# ========================================
calc_kde <- function(locs, sp, com_por, yr, save = T, output = F, proj.type, sigma.fixed = F, var = "qty", var.name = "qtykept") 	{
  ### Get rid of NA in gear
  gear.types <- unique(locs$geargroup)
  gear.types <- gear.types[gear.types != "NA"]
  
  if(sigma.fixed == F){
	  sig_gear <- 
	    bandwidth_df$bandwidth[
	      which(bandwidth_df$gear%in%gear.types)]
  
	  sig_gear2 <- max(sig_gear) #if multiple gear types, take the maximum bandwidth for any one gear type
  
  } else{sig_gear2 <- sigma.fixed}
  
  
  if(var == "qty"){
	  print("step1")
	  # Convert x into spatial points data frame 
    
	  mypattern.qty <- ppp(x = jitter(st_coordinates(locs)[,1]), 
	                                 y = jitter(st_coordinates(locs)[,2]),
	                                 # round extent to nearest thousandth so resolution will work
	                                 xrange=round(st_bbox(locs)[c(1,3)], -3) + c(-200000, 200000), #+ c(-200000, 200000)
	                                 yrange=round(st_bbox(locs)[c(2,4)],-3) + c(-400000, 200000), #+ c(-400000, 200000)
	                                 marks = locs[[eval(var.name)]] )
  
  #
	  print("step2")
	  ### Calculate kde for qty_kept
	  #Don't convert to km first, just use eps=1000 for pixel size
	  kde.test1 <- density.ppp(mypattern.qty, kernel = "gaussian", sigma = sig_gear2, eps = 1000,
	                           weights = mypattern.qty$marks, positive = TRUE)
  
	  ###Scale values by 1000*1000 so that sum to the sum of the original qtykept
	  ### Output= qtykept per km2
	  kde.test1$v <- 1000*1000*kde.test1$v
	  
	  rast1 <- raster::raster(kde.test1) #based on reading here. https://maczokni.github.io/crime_mapping_textbook/studying-spatial-point-patterns.html
    crs(rast1) <- projalt
	  
	  # First expand bounding box to consistent regional scale
	  # extent (xmin, xmax, ymin, ymax)
	  regional <- extent(70000, 900000, 3500000, 5500000) #Note need to change this so as to not hardcode; extent(70000, 900000, 3500000, 5500000) - but having some errors of grid too small
	  rast1 <- extend(rast1, regional, value = 0)
	  
	  if(save == T){
		  print("step3")

	    ### Write raster to tif file
	    name1 <- paste("Confidential/footprint_files/", sp, "/KDE/", var.name, "_spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", unique(yr),"_logbook", ".tif", sep="")
	    writeRaster(x=rast1, filename=name1, format="GTiff", overwrite=T)
	    
	    ### Write raster to .rst file
	    name1b <- paste("Confidential/footprint_files/", sp, "/KDE/",var.name,"_spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", unique(yr),"_logbook",".rst", sep="")
	    #crs(rast1) <- "+proj=longlat" # JS added 03-30-2022. https://stackoverflow.com/questions/69014085/saving-raster-not-updated-for-proj-6
	    writeRaster(x=rast1, filename=name1b, format="IDRISI", overwrite=T)
	  }
	  print("step4")
  	#rast.out <- vector("list", length=1)  
  	  if(output==T){
  	  	#rast.out[[1]] <- rast1
  	  	rast.out <- rast1
  	  	
  	  }
    	return(rast.out)
  
  }#end of if statement for var==qty

  
  if(var == "days"){
	  print("step1")
	  mypattern.days <- spatstat::ppp(x=jitter(st_coordinates(locs)[,1]),
	                                  y=jitter(st_coordinates(locs)[,2]),
	                                  # round extent to nearest thousandth so resolution will work
	                                  xrange=round(st_bbox(locs)[c(1,3)], -3)+ c(-50000,50000),
	                                  yrange=round(st_bbox(locs)[c(2,4)],-3) + c(-50000,50000),
	                                  marks=locs[[eval(var.name)]])
	  
	  print("step2")
	  ### Calculate kde for adj_fisher_days
	  kde.test2 <- spatstat::density.ppp(mypattern.days, kernel = "gaussian", sigma= sig_gear2, eps = 1000,
	                           weights = mypattern.days$marks, positive = TRUE)
  
	  ###Scale values by 1000*1000 so that sum to the sum of the original adj_fisher_days
	  ### Output= adj_fisher_days per km2
	  kde.test2$v <- 1000*1000*kde.test2$v
	  rast2 <- raster::raster(kde.test2)#, crs = CRS(proj.type))
	  crs(rast2) <- sp::CRS(projalt) # this was the last line jameal added based on reading here. https://maczokni.github.io/crime_mapping_textbook/studying-spatial-point-patterns.html
	  
	  # First expand bounding box to consistent regional scale
	  # extent (xmin, xmax, ymin, ymax)
	  regional <- extent(70000, 640000, 3500000, 5500000) #Note need to change this so as to not hardcode
	  
	  rast2 <- extend(rast2, regional, value=0)
	  
	  print("step3")
	  if(save==T){
	    
	    # name1 <- paste("KDE/","qty_","spp_", sp, "_com_port_", unique(com_por), "_year_", unique(yr),".tif", sep="")
	    name2 <- paste("Confidential/footprint_files/", sp, "/KDE/","days_","spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", unique(yr), "_logbook", ".tif", sep="")
	    
	    ### Write raster to tif file
	    # writeRaster(rast1, name1, overwrite=T)
	    #crs(rast2) <- "+proj=longlat" # JS added 03-30-2022. https://stackoverflow.com/questions/69014085/saving-raster-not-updated-for-proj-6
	    writeRaster(rast2, name2, overwrite=T)
	    
	    ## Write raster to .rst file
	    name2b <- paste("Confidential/footprint_files/", sp, "/KDE/","days_","spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", unique(yr), "_logbook", ".rst", sep="")
	    writeRaster(x=rast2, filename=name2b, format="IDRISI", overwrite=T)
	  }
	  
	  print("step4")
    
	#rast.out <- vector("list", length=1)  
	  if(output==T){
	  	#rast.out[[1]] <- rast2
	  	rast.out <- rast2
	  }
  	return(rast.out)

  } # end of if var=days
}


# =========================================================
# = Use Percent Volume Contour Function from adehabitatHR =
# =========================================================
library(adehabitatHR)

calc_pvc <- function(kde.out, cont, sp, com_por, yr, save=T, output=F, var="qty", var.name="qtykept"){
  
  if(var == "qty"){
	  #qty.r <- kde.out[[1]]
	  #qty.r <- rast1
    qty.r <- kde.out
	  
	  kde.sp <- as(qty.r, "SpatialPixelsDataFrame")
	  qty.ud <- new("estUD", kde.sp)
	  qty.ud@vol = FALSE
	  qty.ud@h$meth="Plug-in Bandwidth"
	  qty.ud.vol <- getvolumeUD(qty.ud, standardize=TRUE)
	  qty.ud.vol.raster <- raster(qty.ud.vol)

	  qty.percvol <- getverticeshr(qty.ud, percent = cont, ida = NULL, unin = "m",
	  unout = "km2", standardize = TRUE) %>%
	    st_as_sf
  
	  if(save==T){

	 	name1 <- paste(var.name, "_", cont, "_spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", unique(yr), "_logbook", sep="")
	 	
	 	wd <- getwd()
	 	st_write(qty.percvol, paste0(wd,'/Confidential/footprint_files/', sp, '/PVC/', name1,'.shp'), append=FALSE)
	 	
	  }
	  pvc.out <- NA
	  
	  if(output==T){
	    # stk <- stack(pvc.qty, pvc.days)
	    # names(stk) <- c("qty", "days")
	    # return(stk)
		#pvc.out[[1]] <- qty.percvol
		pvc.out <- qty.percvol
		
	  }
	return(pvc.out)
	  
  }

  if(var=="days"){
	 #days.r <- kde.out[[1]]
	  days.r <- kde.out
	  
	  kde.sp <- as(days.r, "SpatialPixelsDataFrame")
	  days.ud <- new("estUD", kde.sp)
	  days.ud@vol = FALSE
	  days.ud@h$meth="Plug-in Bandwidth"
	  days.ud.vol <- getvolumeUD(days.ud, standardize=TRUE)
	  days.ud.vol.raster <- raster(days.ud.vol)

	  days.percvol <- getverticeshr(days.ud, percent = cont,ida = NULL, unin = "m",
	  unout = "km2", standardize=TRUE)
  
	  if(save==T){
	    name1 <- paste(var.name, "_", cont, "_spp_", sp, "_", yrSrt, "_", yrEnd, "_com_port_", unique(com_por), "_year_", unique(yr), sep="")
	    
	    wd <- getwd()
	    st_write(days.percvol, paste0(wd,'/Confidential/footprint_files/', sp, '/PVC/', name1,'.shp'), append=FALSE)

	  }
  	
	#pvc.out <- vector("list", length=1)  
	  pvc.out <- NA
	  
	  if(output==T){
	    # stk <- stack(pvc.qty, pvc.days)
	    # names(stk) <- c("qty", "days")
	    # return(stk)
		#pvc.out[[1]] <- days.percvol
		pvc.out <- days.percvol
	  }
	  return(pvc.out)
  }
}


