### Calculate fishing footprints in each year, port community, and species target group

library(data.table)
library(here)
library(maps)
library(mapdata)
library(maptools)
library(raster)
library(spatstat)
library(tidyverse)
library(janitor)
library(sf)

# ===========================
# Load Port Information 
# ===========================
ports <- read.csv(here::here('Non_Confidential', 'iopac_logbook_conversion_table.csv'))

# ===========================
# Identify years of interest 
# ===========================
yrSrt <- 1994 # options are '1987' or '1994' or '2011'
yrEnd <- 2020 # options are 2019 (exclude the new Amendment 28 changes) or 2020

# ===========================
# Load and Clean Fishticket/Logbook Data 
# ===========================
##Load raw joined fishticket/logbook data from Samhouri et al. 2024
load(here('Confidential/data/LB.ShortForm.with.Hake.Strat 26 Apr 2022.RData')) #data from PSMFC PacFIN database

new_df <- LB.ShortForm.with.Hake.Strat %>%
  filter(Strategy != 'HAKE', #remove fishing events targeting Pacific hake
         RYEAR >= yrSrt, # >=1994 so we have consistent info from all 3 states and this consists of vessels that separated for hake vs. other bottom trawling
         Total_All_FMP.kg != 0 | is.na(Total_All_FMP.kg) == TRUE) # drop rows with zero or NA catch of FMP species

# join port information
new_df2 <- new_df %>%
  left_join(ports, by = c('RPCID' = 'PACFIN_PORT_CODE')) %>%
  filter(is.na(IOPAC)!= TRUE) %>%
  filter(IOPAC != 'Unknown Ports')

#Retain necessary attributes
dat_subset <- new_df2 %>%
  dplyr::select(!c(DRVID, TOWNUM, AGID, RPCID, RDAY, DDATE, RDATE, TOWDATE, SET_LAT, SET_LONG, UP_LAT, UP_LONG, LATLONG_TYPE, SET_TIME, DURATION, ADJ_TOWTIME, ARID_PSMFC, BLOCK, BLOCK_OR, TICKET_DATE, DEPTH1, PERMID_1, LEN_ENDOR_1, bimo, InsideEEZ, InsideAllPoly, BlkCentrd_Lat, BlkCentrd_Long, Best_Lat, Best_Long, DepthGIS.m, Month, BestBtmDepth.m, GoodTow, Pacfin.Port.Description, Pacfin.Port.Name, PACFIN_LAT, PACFIN_LONG, Map_order, in.LB.ShortForm.with.Hake.Strat..24.Mar.2022, PRIMARY_LATITUDE, PRIMARY_LONGITUDE, STATE, UP_TIME, State.Waters.PSMFC, Agency.Name, Strategy, TRIP_ID))

#melt to group by haul and sum per species
dat_long <- dat_subset %>%
  reshape2::melt(id.vars = c('RYEAR', 'GRID', 'PACFIN_TARGET', 'Key', 'FTID', 'IOPAC', 'Total_All_FMP.kg')) %>%
  transform(FTID = as.character(FTID)) %>%
  filter(value > 0)

####Target species to use for defining each haul into subfisheries (e.g. DTS bottom trawl, non-DTS bottom trawl, and Midwater rockfish trawl)
DTS_names <- c('DOVR.kg',  'LSP1.kg',  'THDS.kg',  'SSP1.kg',  'SABL.kg',  'SBL1.kg', 'THD1.kg')
nDTS_names <- c('ART1.kg',  'REX1.kg',  'BSOL.kg', 'EGL1.kg',  'SLNS.kg',  'ARRA.kg',  'CLPR.kg',  'DBRK.kg', 'STL1.kg',  'YEYE.kg',  'LCD1.kg', 'CSL1.kg',  'STRK.kg',  'CNRY.kg',  'SBLY.kg',  'PDB1.kg',  'PTR1.kg', 'SRK1.kg',  'SRKR.kg',  'SNOS.kg',  'BCAC.kg', 'ARTH.kg',  'CSOL.kg',  'PDAB.kg',  'PTRL.kg',  'REX.kg', 'EGLS.kg',  'ARR1.kg',  'CLP1.kg',  'DBR1.kg', 'CNR1.kg',  'SBL1.kg',  'SNS1.kg',  'BCC1.kg',  'POP.kg',  'YEY1.kg', 'LCOD.kg', 'BGL1.kg', 'SSOL.kg', 'DSRK.kg', 'PCOD.kg', 'USKT.kg', 'STRY.kg', 'OSKT.kg', 'NUSP.kg', 'RSOL.kg', 'BNK1.kg', 'SSRK.kg', 'LSKT.kg', 'VRM1.kg', 'URCK.kg', 'BLK1.kg', 'NUSF.kg', 'RCK1.kg', 'LSRK.kg', 'BRW1.kg', 'KLPG.kg', 'GSP1.kg', 'BLU1.kg', 'USLF.kg', 'GSR1.kg', 'RATF.kg', 'CSKT.kg', 'OSKT.kg', 'RST1.kg', 'GRS1.kg', 'RCK7.kg', 'ROS1.kg', 'CBZN.kg', 'CHN1.kg', 'BSKT.kg', 'STR1.kg', 'ORCK.kg', 'CWC1.kg', 'SCR1.kg', 'USLP.kg', 'GRDR.kg', 'COP1.kg')
MDT_names <- c('WDOW.kg', 'YTRK.kg', 'WDW1.kg', 'YTR1.kg', 'PWHT.kg') 


####Aggregate by fishing fleet - sum species first - and assign targeted species based on which species has the highest proportion of landings
dat_long_agg <- dat_long %>%
  transform(subfishery = ifelse(GRID == 'MDT', 'MDT',
                                ifelse(variable %in% c(DTS_names), 'DTS', 'nDTS'))) %>%
  distinct() %>%
  transform(subfishery = ifelse(!is.na(GRID), subfishery, 
                                ifelse(variable %in% c(DTS_names), 'DTS',
                                       ifelse(variable %in% c(MDT_names), 'MDT', 'nDTS')))) %>%
  distinct() %>%
  #sum of landings by subfishery
  group_by(FTID, Key, subfishery) %>% 
  dplyr::summarize(subfishery_kgs = sum(value), .groups = 'keep') %>%
  ungroup() %>%
  
  #sum of landings by haul
  group_by(FTID, Key) %>%
  mutate(tot_kgs_haul = sum(subfishery_kgs)) %>%
  ungroup() %>%
  
  #calculate proportions of total catch per subfishery
  group_by(FTID, Key) %>%
  mutate(prop_lbs_subfish = subfishery_kgs/tot_kgs_haul) %>%
  
  #assign a TARGET to the haul, defined as the species with the largest proportion of landings for that haul
  mutate(first = dplyr::first(prop_lbs_subfish, order_by = desc(prop_lbs_subfish)),
         second = dplyr::nth(prop_lbs_subfish, n = 2, order_by = desc(prop_lbs_subfish)),
         target_sub = dplyr::first(subfishery, order_by = desc(prop_lbs_subfish)),
         second_sub = dplyr::nth(subfishery, n = 2, order_by = desc(prop_lbs_subfish)))

ftid_attributes <- new_df2 %>%
  dplyr::select(FTID, Key, RYEAR, DRVID, TRIP_ID, Best_Lat, Best_Long, IOPAC) %>% distinct() 

dat_w_targets <- dat_long_agg %>%
  left_join(ftid_attributes, by=c('FTID', 'Key')) %>%
  filter(subfishery == target_sub) %>%
  group_by(FTID, Key, RYEAR, IOPAC, DRVID, TRIP_ID, target_sub, Best_Lat, Best_Long) %>%
  dplyr::summarize(qtykept = tot_kgs_haul,
                   target_species_landings = subfishery_kgs, #let's keep subfishery landings weight also for the Sensitivity calculation
                   .groups = 'keep')

# use all IOPAC port groups in logbook data from years of interest
focal <- new_df2 %>% 
  filter(RYEAR %in% seq(yrSrt, yrEnd, 1)) %>% 
  dplyr::select(IOPAC) %>% 
  distinct() %>% 
  pull() 

# ==========
# Parameters for kernel density and percent volume contour functions
# ==========
source('scripts/logbook_fishshed_functions.R') #functions for creating KDE and PVC's

bandwidth_df <- data.frame(gear = c('Trawl'),
                           bandwidth = c(10000))

# Use UTM 10
proj1 <- '+proj=utm +zone=10 +datum=WGS84 +units=m' # originally +init=EPSG:26910. could use 26910 with sf()
projalt <- 26910




##############################################
####Calculate annual footprints using for loop across all fisheries, but do each Percent Volume Contour-of-interest separately####
##############################################


##############################################
####Percent Volume Contours of interest #1####
##############################################
cont <- 50 #could be any value; we compared 50, 75 and 95

#run across all three fisheries
sp <- c('DTS', 'nDTS', 'MDT')
for (i in 1:length(sp)) {
  new_df3 <- dat_w_targets %>%
    dplyr::select(Best_Lat, Best_Long, TRIP_ID, DRVID, FTID, target_sub, qtykept, target_species_landings) %>%
    mutate(geargroup = 'Trawl') %>%
    filter(target_sub == sp[i]) %>%
    rename(fish_group = target_sub)
  
  agg_comm_focal4 <- data.table(copy(new_df3))
  
  # 3. save haul-level data for summary statistics for within vs. outside WEA portion of footprint
  haul_dat <- new_df3 %>% filter(qtykept > 0)
  write_rds(haul_dat, file = here('Confidential', 'processed logbook', paste0(sp[i], '_', yrSrt, '_', yrEnd, '_logbook_hauls_ports_IOPACs_depths.rds')))
  
  # look at num vessels by IOPAC and era
  vessels_iopac_era <- agg_comm_focal4 %>% 
    dplyr::select(IOPAC, RYEAR, DRVID) %>%
    group_by(IOPAC, RYEAR) %>%
    dplyr::summarise(num_unique_vessels = n_distinct(DRVID), .groups = 'keep')
  
  # 4. confidentiality: drop IOPAC-era/year combos with <3 vessels.... 
  agg_comm_focal5 <- agg_comm_focal4 %>% left_join(vessels_iopac_era) %>%
    #filter(num_unique_vessels >= 3) %>% #let's not do this until, or if, we need to do it for showing stats/maps
    filter(qtykept > 0)  
  
  
  # ==============================================
  # = Calculate KDE and PVC for fishing location =
  # ==============================================
  # Do all years at once BY fishing fleet and year
  
  # By year
  ### footprint for qty kept
  setorder(agg_comm_focal5, IOPAC, RYEAR)
  
  ports.all <- unique(agg_comm_focal5$IOPAC)
  
  qty.kde.era <- agg_comm_focal5[IOPAC %in% ports.all, j = {
    print(paste0(IOPAC, '_', RYEAR, '_', sp[i]))
    t.dt <- .SD
    print('Convert to spdf')
    
    out <- make_spdf(t.dt, sp = sp[i], com_por = IOPAC, lon.var = 'Best_Long', lat.var = 'Best_Lat', 
                     yr = RYEAR,
                     save = F, output = T, proj.type = projalt)
    print('Kernel Density Estimation')
    kde.qty <- calc_kde(out, sp = sp[i], com_por = IOPAC, yr = RYEAR,
                        save = T, output = T, 
                        proj.type = projalt, sigma.fixed = F, var = 'qty', var.name = 'qtykept')
    
    print('Percent Volume Contours')
    pvc.qty <- calc_pvc(kde.qty, cont = cont, sp = sp[i], com_por = IOPAC, yr = RYEAR,
                        save = T, output = T, var = 'qty')
    
    list('done')
    
  }, by = c('IOPAC', 'geargroup', 'RYEAR'), .SDcols = c('geargroup', 'qtykept', 'Best_Lat', 'Best_Long')
  ]
  #warnings are okay based on whether a year x port x fishery occurred
  
  # 5. combine all of the shape files into one
  
  pvcs <- list.files(here('Confidential', 'footprint_files', sp[i], 'PVC'), full.names = T) %>% 
    str_subset('qtykept') %>% str_subset(paste0(as.character(cont), '_spp')) %>% str_subset(paste0(yrSrt, '_', yrEnd)) %>% str_subset('.shp')
  
  
  # import all shapes in the given directory
  pvcs_sf <- pvcs %>% purrr::map_df(read_sf)
  
  #fix names for ports, years, contours
  port_names <- pvcs %>% map_chr(~str_match(.x, 'port_\\s*(.*?)\\s*_year')[,2])
  year_names <- pvcs %>% map_chr(~str_match(.x, 'year_\\s*(.*?)\\s*_logbook')[,2]) #years
  cont_names <- pvcs %>% map_chr(~str_match(.x, 'qtykept_\\s*(.*?)\\s*_spp')[,2]) #contours
  
  pvcs_sf <- pvcs_sf %>% mutate(port_name=port_names,year=year_names, cont=cont_names) %>%
    filter(year >= yrSrt & year <= yrEnd)#years
  
  # 6. save combined shape files
  write_sf(pvcs_sf, dsn = paste0('Confidential/', 'data/', 'PVC/', 'logbook_', sp[i], '_', yrSrt, '_', yrEnd, '_', cont, '_pvcs_yrs_target.shp')) #target species designation
}

##############################################
####Percent Volume Contours of interest #2####
##############################################
cont <- 75 #could be any value; we compared 50, 75 and 95

#run across all three fisheries
sp <- c('DTS', 'nDTS', 'MDT')
for (i in 1:length(sp)) {
  new_df3 <- dat_w_targets %>%
    dplyr::select(Best_Lat, Best_Long, TRIP_ID, DRVID, FTID, target_sub, qtykept, target_species_landings) %>%
    mutate(geargroup = 'Trawl') %>%
    filter(target_sub == sp[i]) %>%
    rename(fish_group = target_sub)
  
  agg_comm_focal4 <- data.table(copy(new_df3))
  
  # 3. save haul-level data for summary statistics for within vs. outside WEA portion of footprint
  haul_dat <- new_df3 %>% filter(qtykept > 0)
  write_rds(haul_dat, file = here('Confidential', 'processed logbook', paste0(sp[i], '_', yrSrt, '_', yrEnd, '_logbook_hauls_ports_IOPACs_depths.rds')))
  
  # look at num vessels by IOPAC and era
  vessels_iopac_era <- agg_comm_focal4 %>% 
    dplyr::select(IOPAC, RYEAR, DRVID) %>%
    group_by(IOPAC, RYEAR) %>%
    dplyr::summarise(num_unique_vessels = n_distinct(DRVID), .groups = 'keep')
  
  # 4. confidentiality: drop IOPAC-era/year combos with <3 vessels.... 
  agg_comm_focal5 <- agg_comm_focal4 %>% left_join(vessels_iopac_era) %>%
    #filter(num_unique_vessels >= 3) %>% #let's not do this until, or if, we need to do it for showing stats/maps
    filter(qtykept > 0)  
  
  
  # ==============================================
  # = Calculate KDE and PVC for fishing location =
  # ==============================================
  # Do all years at once BY fishing fleet and year
  
  # By year
  ### footprint for qty kept
  setorder(agg_comm_focal5, IOPAC, RYEAR)
  
  ports.all <- unique(agg_comm_focal5$IOPAC)
  
  qty.kde.era <- agg_comm_focal5[IOPAC %in% ports.all, j = {
    print(paste0(IOPAC, '_', RYEAR, '_', sp[i]))
    t.dt <- .SD
    print('Convert to spdf')
    
    out <- make_spdf(t.dt, sp = sp[i], com_por = IOPAC, lon.var = 'Best_Long', lat.var = 'Best_Lat', 
                     yr = RYEAR,
                     save = F, output = T, proj.type = projalt)
    print('Kernel Density Estimation')
    kde.qty <- calc_kde(out, sp = sp[i], com_por = IOPAC, yr = RYEAR,
                        save = T, output = T, 
                        proj.type = projalt, sigma.fixed = F, var = 'qty', var.name = 'qtykept')
    
    print('Percent Volume Contours')
    pvc.qty <- calc_pvc(kde.qty, cont = cont, sp = sp[i], com_por = IOPAC, yr = RYEAR,
                        save = T, output = T, var = 'qty')
    
    list('done')
    #}) 
  }, by = c('IOPAC', 'geargroup', 'RYEAR'), .SDcols = c('geargroup', 'qtykept', 'Best_Lat', 'Best_Long')
  ]
  #warnings are okay
  
  # 5. combine all of the shape files into one
  
  pvcs <- list.files(here('Confidential', 'footprint_files', sp[i], 'PVC'), full.names = T) %>% 
    str_subset('qtykept') %>% str_subset(paste0(as.character(cont), '_spp')) %>% str_subset(paste0(yrSrt, '_', yrEnd)) %>% str_subset('.shp')
  
  
  # import all shapes in the given directory
  pvcs_sf <- pvcs %>% purrr::map_df(read_sf)
  
  #fix names for ports, years, contours
  port_names <- pvcs %>% map_chr(~str_match(.x, 'port_\\s*(.*?)\\s*_year')[,2])
  year_names <- pvcs %>% map_chr(~str_match(.x, 'year_\\s*(.*?)\\s*_logbook')[,2]) #years
  cont_names <- pvcs %>% map_chr(~str_match(.x, 'qtykept_\\s*(.*?)\\s*_spp')[,2]) #contours
  
  
  # pvcs_sf <- pvcs_sf %>% mutate(port_name=port_names,period=period_names) #eras
  pvcs_sf <- pvcs_sf %>% mutate(port_name=port_names,year=year_names, cont=cont_names) %>%
    filter(year >= yrSrt & year <= yrEnd)#years
  
  # 6. save combined shape files
  write_sf(pvcs_sf, dsn = paste0('Confidential/', 'data/', 'PVC/', 'logbook_', sp[i], '_', yrSrt, '_', yrEnd, '_', cont, '_pvcs_yrs_target.shp')) #target species designation
}

##############################################
####Percent Volume Contours of interest #3####
##############################################
cont <- 95 #could be any value; we compared 50, 75 and 95

sp <- c('DTS', 'nDTS', 'MDT')
for (i in 1:length(sp)) {
  new_df3 <- dat_w_targets %>%
    dplyr::select(Best_Lat, Best_Long, TRIP_ID, DRVID, FTID, target_sub, qtykept, target_species_landings) %>%
    mutate(geargroup = 'Trawl') %>%
    filter(target_sub == sp[i]) %>%
    rename(fish_group = target_sub)
  
  agg_comm_focal4 <- data.table(copy(new_df3))
  
  # 3. save haul-level data for summary statistics for within vs. outside WEA portion of footprint
  haul_dat <- new_df3 %>% filter(qtykept > 0)
  write_rds(haul_dat, file = here('Confidential', 'processed logbook', paste0(sp[i], '_', yrSrt, '_', yrEnd, '_logbook_hauls_ports_IOPACs_depths.rds')))
  
  # look at num vessels by IOPAC and era
  vessels_iopac_era <- agg_comm_focal4 %>% 
    dplyr::select(IOPAC, RYEAR, DRVID) %>%
    group_by(IOPAC, RYEAR) %>%
    dplyr::summarise(num_unique_vessels = n_distinct(DRVID), .groups = 'keep')
  
  # 4. confidentiality: drop IOPAC-era/year combos with <3 vessels.... 
  agg_comm_focal5 <- agg_comm_focal4 %>% left_join(vessels_iopac_era) %>%
    #filter(num_unique_vessels >= 3) %>% #let's not do this until, or if, we need to do it for showing stats/maps
    filter(qtykept > 0)  
  
  
  # ==============================================
  # = Calculate KDE and PVC for fishing location =
  # ==============================================
  # Do all years at once BY fishing fleet and year
  
  # By year
  ### footprint for qty kept
  setorder(agg_comm_focal5, IOPAC, RYEAR)
  
  ports.all <- unique(agg_comm_focal5$IOPAC)
  
  qty.kde.era <- agg_comm_focal5[IOPAC %in% ports.all, j = {
    print(paste0(IOPAC, '_', RYEAR, '_', sp[i]))
    t.dt <- .SD
    print('Convert to spdf')
    
    out <- make_spdf(t.dt, sp = sp[i], com_por = IOPAC, lon.var = 'Best_Long', lat.var = 'Best_Lat', 
                     yr = RYEAR,
                     save = F, output = T, proj.type = projalt)
    print('Kernel Density Estimation')
    kde.qty <- calc_kde(out, sp = sp[i], com_por = IOPAC, yr = RYEAR,
                        save = T, output = T, 
                        proj.type = projalt, sigma.fixed = F, var = 'qty', var.name = 'qtykept')
    
    print('Percent Volume Contours')
    pvc.qty <- calc_pvc(kde.qty, cont = cont, sp = sp[i], com_por = IOPAC, yr = RYEAR,
                        save = T, output = T, var = 'qty')
    
    list('done')
    #}) 
  }, by = c('IOPAC', 'geargroup', 'RYEAR'), .SDcols = c('geargroup', 'qtykept', 'Best_Lat', 'Best_Long')
  ]
  #warnings are okay
  
  # 5. combine all of the shape files into one
  
  pvcs <- list.files(here('Confidential', 'footprint_files', sp[i], 'PVC'), full.names = T) %>% 
    str_subset('qtykept') %>% str_subset(paste0(as.character(cont), '_spp')) %>% str_subset(paste0(yrSrt, '_', yrEnd)) %>% str_subset('.shp')
  
  
  # import all shapes in the given directory
  pvcs_sf <- pvcs %>% purrr::map_df(read_sf)
  
  #fix names for ports, years, contours
  port_names <- pvcs %>% map_chr(~str_match(.x, 'port_\\s*(.*?)\\s*_year')[,2])
  year_names <- pvcs %>% map_chr(~str_match(.x, 'year_\\s*(.*?)\\s*_logbook')[,2]) #years
  cont_names <- pvcs %>% map_chr(~str_match(.x, 'qtykept_\\s*(.*?)\\s*_spp')[,2]) #contours
  
  
  # pvcs_sf <- pvcs_sf %>% mutate(port_name=port_names,period=period_names) #eras
  pvcs_sf <- pvcs_sf %>% mutate(port_name=port_names,year=year_names, cont=cont_names) %>%
    filter(year >= yrSrt & year <= yrEnd)#years
  
  # 6. save combined shape files
  write_sf(pvcs_sf, dsn = paste0('Confidential/', 'data/', 'PVC/', 'logbook_', sp[i], '_', yrSrt, '_', yrEnd, '_', cont, '_pvcs_yrs_target.shp')) #target species designation
}



############################################
#Create footprint across all years of data; changing 'yr = RYEAR' to 'yr = fish_group' in the code below
# 1. instead of by species, across all years
##############################################
sp = 'all'

new_df3 <- dat_w_targets %>%
  dplyr::select(Best_Lat, Best_Long, TRIP_ID, DRVID, FTID, target_sub, qtykept) %>%
  mutate(geargroup = 'Trawl') %>%
  rename(fish_group = target_sub) %>%
  filter(RYEAR >= yrSrt)

agg_comm_focal4 <- data.table(copy(new_df3))

# 3. save haul-level data for summary statistics for within vs. outside WEA portion of footprint
haul_dat <- new_df3 %>% filter(qtykept > 0)
write_rds(haul_dat, file = here('Confidential', 'processed logbook', paste0(sp, '_', yrSrt, '_', yrEnd, '_logbook_hauls_ports_IOPACs_depths.rds')))

# look at num vessels by IOPAC and era
vessels_iopac_era <- agg_comm_focal4 %>% 
  dplyr::select(IOPAC, RYEAR, DRVID, fish_group) %>%
  group_by(IOPAC, RYEAR, fish_group) %>%
  dplyr::summarise(num_unique_vessels = n_distinct(DRVID), .groups = 'keep')

# 4. confidentiality: drop IOPAC-era/year combos with <3 vessels.... 
agg_comm_focal5 <- agg_comm_focal4 %>% left_join(vessels_iopac_era) %>%
  #filter(num_unique_vessels >= 3) %>% #let's not do this until, or if, we need to do it for showing stats/maps
  filter(qtykept > 0)  

# By Target
### footprint for qty kept
setorder(agg_comm_focal5, IOPAC, fish_group)

ports.all <- unique(agg_comm_focal5$IOPAC)

qty.kde.era <- agg_comm_focal5[IOPAC %in% ports.all, j = {
  print(paste0(IOPAC, '_', fish_group))
  t.dt <- .SD
  print('Convert to spdf')
  
  out <- make_spdf(t.dt, sp=sp, com_por=IOPAC, lon.var='Best_Long', lat.var='Best_Lat', 
                   yr = fish_group,
                   save = F, output = T, proj.type = projalt)
  print('Kernel Density Estimation')
  kde.qty <- calc_kde(out, sp = sp, com_por = IOPAC, yr = fish_group, #yr = period,
                      save = T, output = T, 
                      proj.type = projalt, sigma.fixed = F, var = 'qty', var.name = 'qtykept')
  
  print('Percent Volume Contours')
  pvc.qty <- calc_pvc(kde.qty, cont=cont, sp=sp, com_por=IOPAC, yr = fish_group, #yr = period, 
                      save=T, output=T, var='qty')
  
  list('done')
  
}, by=c('IOPAC', 'geargroup', 'fish_group'
), .SDcols=c('geargroup', 'qtykept', 'Best_Lat', 'Best_Long')
]
#warnings are okay

pvcs <- list.files(here('Confidential', 'footprint_files', sp, 'PVC'), full.names = T) %>% 
  str_subset('qtykept') %>% str_subset(paste0(as.character(cont), '_spp')) %>% str_subset(paste0(yrSrt, '_', yrEnd)) %>% str_subset('.shp')


# import all shapes in the given directory
pvcs_sf <- pvcs %>% purrr::map_df(read_sf)

#fix names for ports, years, contours
port_names <- pvcs %>% map_chr(~str_match(.x, 'port_\\s*(.*?)\\s*_year')[,2])
year_names <- pvcs %>% map_chr(~str_match(.x, 'year_\\s*(.*?)\\s*_logbook')[,2]) #fishery
cont_names <- pvcs %>% map_chr(~str_match(.x, 'qtykept_\\s*(.*?)\\s*_spp')[,2]) #contour

pvcs_sf <- pvcs_sf %>% 
  mutate(port_name = port_names, 
         year = year_names, 
         cont = cont_names) %>%
  rename('subfishery' = 'year')

# 6. save combined shape files
write_sf(pvcs_sf, dsn = paste0('Confidential/', 'data/', 'PVC/', 'logbook_', sp, '_', yrSrt, '_', yrEnd, '_', cont, '_pvcs_yrs_target.shp')) #target species designation




