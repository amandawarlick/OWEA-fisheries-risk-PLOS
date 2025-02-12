####Script uses footprint files generated in 'create_footprints_target_sp.R' script to calculate and map commercial fishing footprints, calculate risk (exposure and vulnerability) of displacement from Offshore Wind Energy Areas, and identify potential changes in this risk when the effects of climate change (via any shifts in targeted species' spatial distributions) are accounted for.

####Load packages####
library(tidyr)
library(reshape2)
library(data.table)
library(ggplot2)
library(cowplot)
library(knitr)
library(stringr)
library(here)
library(Hmisc)
library(lubridate)
library(coda)
library(readr)
library(captioner)
library(latex2exp) #TeX
library(dplyr)
library(visreg) 
library(tictoc) 
library(tidyverse)
library(magrittr)
library(rnaturalearth)
library(sf)
library(RANN)
library(paletteer)
library(ggsci)
library(patchwork)
library(glue)
library(ggh4x)
library(ggnewscale)
library(ggrepel)
library(zoo)
library(slider)

######################################
####User inputs####
######################################
#Period of interest for footprints
yrSrt <- 1994 #earliest year of data to be used to calculate footprints; options considered were "1987" or "1994" or "2011"
yrEnd <- 2020

#PVC contour of interest
cont = "75" #options used for main manuscript ("75") and for sensitivity comparisons in appendix were ("50" and "95")

#OWE areas used
oweas <- "ORCAS" #either Oregon Call Areas: "ORCAS"; or Oregon WEAs: "ORWEAS"

#Ports of interest and factor level order
ports.all <- c('Puget Sound', 'North WA Coast', 'WA Coast', 'Astoria', 'Tillamook', 'Newport', 'Coos Bay', 'Brookings', 'Crescent City', 'Eureka', 'Fort Bragg','Bodega Bay', 'San Francisco', 'Monterey', 'Morro Bay', 'Santa Barbara', 'Los Angeles', 'San Diego')

######################################
####map objects and WEA boundaries####
######################################
projection_extent <- readRDS(here('Confidential', 'data', 'cropped_domain_for_projected_sdms.rds'))
bbox = st_bbox(projection_extent)

coast_outline <- read_sf(here('GIS layers', 'NWFSC_LOADER_TMER', 'NWFSC_LOADER_TMER.shp')) %>% st_as_sf() %>%
  filter(NAME %in% c('California','Oregon','Washington','Nevada')) %>%
  st_transform(crs = "+proj=utm +zone=10 +datum=WGS84 +units=km")

gr <- projection_extent %>% st_make_grid(cellsize = 10, what = 'centers') %>% st_as_sf() %>%
  st_intersection(projection_extent)
gr_xy <- st_coordinates(gr)

#Oregon wind area shapefiles
if(oweas == "ORCAS"){
  OR_owe_shp <- read_sf(here('GIS layers','ORwindenergyareas','CallAreas_OR_Outline_2022_04_22.shp')) %>%
    st_transform(st_crs(projection_extent)) %>%
    dplyr::rename(Call_Area = Call_Area_)
  OR_owe_shp <- st_zm(OR_owe_shp, drop = TRUE, what = "ZM")
  
} else{
  OR_owe_shp <- read_sf(here('GIS layers','Oregon_Lease_Areas_0','Oregon_Lease_Areas_2024_08_21.shp')) %>% 
    st_transform(st_crs(projection_extent)) %>%
    dplyr::rename(Call_Area = Name,
                  Acres = ACRES,
                  Sq_Miles = SQ_MI)
}

#California wind energy shapefiles
hum_owe_shp <- read_sf(here('GIS layers','CAwindenergyareas','LeaseAreas_CA_Outlines_2022_09_08.shp')) %>% 
  st_transform(st_crs(projection_extent)) %>%
  dplyr::select(-c(Lease_Numb)) %>%
  transform(Call_Area = ifelse(Sq_Miles > 115, 'Morro Bay', 'Humboldt')) %>%
  filter(Call_Area == 'Humboldt')
hum_owe_shp <- st_union(hum_owe_shp[1,], hum_owe_shp[2,])
hum_owe_shp <- st_zm(hum_owe_shp, drop = TRUE, what = "ZM")

mb_owe_shp <- read_sf(here('GIS layers','CAwindenergyareas','LeaseAreas_CA_Outlines_2022_09_08.shp')) %>% 
  st_transform(st_crs(projection_extent)) %>%
  dplyr::select(-c(Lease_Numb)) %>%
  transform(Call_Area = ifelse(Sq_Miles > 115, 'Morro Bay', 'Humboldt')) %>%
  filter(Call_Area == 'Morro Bay')
mb_owe_shp <- st_union(mb_owe_shp[1,], st_union(mb_owe_shp[2,], mb_owe_shp[3,]))
mb_owe_shp <- st_zm(mb_owe_shp, drop = TRUE, what = "ZM")

#use this if want to filter or summarize by WEA
owe_shps <- bind_rows(OR_owe_shp, hum_owe_shp, mb_owe_shp)

#use this if differencing (or intersecting), otherwise treats geometries separately and thinks there's no overlap
owe_shps_union <- st_union(owe_shps)



######################################
####stats for confidentiality########
######################################

log_dat <- read_rds(here('Confidential', 'processed logbook', glue('MDT_{yrSrt}_{yrEnd}_logbook_hauls_ports_IOPACs_depths.rds'))) %>%
  bind_rows(read_rds(here('Confidential', 'processed logbook', glue('nDTS_{yrSrt}_{yrEnd}_logbook_hauls_ports_IOPACs_depths.rds')))) %>%
  bind_rows(read_rds(here('Confidential', 'processed logbook', glue('DTS_{yrSrt}_{yrEnd}_logbook_hauls_ports_IOPACs_depths.rds'))))

#convert to spatial object
log_dat_sf <- log_dat %>% 
  st_as_sf(coords = c("Best_Long", "Best_Lat"),crs = 4326, remove = F) %>% 
  st_transform(st_crs(owe_shps))

#spatial join
log_dat_wea_sf <- log_dat_sf %>% st_join(owe_shps)

#remove geometry 
log_dat_wea <- log_dat_wea_sf %>% st_set_geometry(NULL)

#calculate subfishery/port group unique vessel, trips characteristics
log_dat_wea <- log_dat_wea %>%
  transform(IOPAC = factor(IOPAC, levels = ports.all),
            fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl')))
log_dat_wea <- droplevels(log_dat_wea)


#for statistics
#all hauls
conf_all_combs <- log_dat_wea %>%
  expand(RYEAR, IOPAC, fish_group)
conf_port_fish_group_combs <- log_dat_wea %>%
  expand(IOPAC, fish_group)
conf_fish_group_combs <- log_dat_wea %>%
  expand(fish_group)

conf_dat_year_port_fish_group <- log_dat_wea %>% 
  group_by(RYEAR, IOPAC, fish_group) %>% 
  # total catch
  summarise(all_areas_total_mt = sum(qtykept)/1000,
            all_areas_subfishery_mt = sum(target_species_landings)/1000,
            total_hauls = n(),
            total_trips = n_distinct(TRIP_ID),
            total_unique_vessels = n_distinct(DRVID), .groups = 'keep') %>%
  right_join(conf_all_combs, by = c('RYEAR', 'IOPAC', 'fish_group')) %>%
  arrange(RYEAR, IOPAC, fish_group)

#catch within proposed OWEs
conf_dat_port_fish_group_wea <- log_dat_wea %>%
  filter(!is.na(Call_Area)) %>%
  group_by(IOPAC, fish_group) %>%
  summarise(hauls_in_OWEAs = n(),
            inside_no_years_fished_logbook = n_distinct(RYEAR),
            inside_unique_vessels = n_distinct(DRVID),
            inside_total_mt = sum(qtykept)/1000,
            inside_subfishery_mt = sum(target_species_landings)/1000) %>%
  mutate(across(where(is.numeric), round, 1))

conf_dat_port_fish_group <- log_dat_wea %>% 
  group_by(IOPAC, fish_group) %>% 
  # total catch
  summarise(all_areas_total_mt = sum(qtykept)/1000,
            all_areas_subfishery_mt = sum(target_species_landings)/1000,
            no_of_years_fished = n_distinct(RYEAR),
            total_hauls = n(),
            total_trips = n_distinct(TRIP_ID),
            total_unique_vessels = n_distinct(DRVID), .groups = 'keep') %>%
  right_join(conf_port_fish_group_combs, by = c('IOPAC', 'fish_group')) %>%
  arrange(IOPAC, fish_group) %>%
  left_join(conf_dat_port_fish_group_wea) %>%
  mutate(across(where(is.numeric), round, 1))



######################################
####Figure 1: footprint map across all years####
######################################
footprints <- read_sf(here('Confidential', 'data', 'PVC', glue('logbook_all_{yrSrt}_{yrEnd}_{cont}_pvcs_yrs_target.shp'))) %>% 
  st_transform(st_crs(projection_extent)) %>% 
  transform(port_name = factor(port_name, levels = ports.all),
            subfishery = factor(subfishery, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  group_by(port_name, subfishery) %>%
  left_join(conf_dat_port_fish_group, by = c("port_name" = "IOPAC", "subfishery" = "fish_group"))
footprints$port_name <- droplevels(footprints$port_name)

bb <- st_bbox(footprints)

#colors for footprints and port locations - needs to correspond to the number of ports that have footprints for the time period of interest
cols <- c('deeppink', 'cyan','gray70', 'violetred4', 'orchid3',  'dodgerblue3', 'deepskyblue1', '#4aaaa5', '#a3d39c', 
          '#f6b61c', 'chocolate2', 'orangered2', 'red3', 'red4', 'gray30', 'black', 'mintcream', 'white')

#coordinates for ports for maps
port_locs <- read.csv(here::here('Non_Confidential', 'iopac_logbook_conversion_table.csv')) %>%
  dplyr::select(IOPAC, Pacfin.Port.Name, PACFIN_LAT, PACFIN_LONG, PRIMARY_LATITUDE, PRIMARY_LONGITUDE) %>%
  filter(IOPAC %in% unique(footprints$port_name)) %>%
  filter(Pacfin.Port.Name != "N PUGET S") %>%
  mutate(PACFIN_LONG = ifelse(IOPAC == "San Francisco" | Pacfin.Port.Name == "SF", -122.4, PACFIN_LONG + 0.1),
         PACFIN_LONG = ifelse(IOPAC == "San Diego", PRIMARY_LONGITUDE, PACFIN_LONG),
         PACFIN_LAT = ifelse(IOPAC == "San Diego", PRIMARY_LATITUDE, PACFIN_LAT)) %>%
  filter(!is.na(PACFIN_LONG)) %>%
  dplyr::select(IOPAC, Pacfin.Port.Name, PACFIN_LAT, PACFIN_LONG) %>%
  rowwise() %>%
  mutate(Pacfin.Port.Name = str_to_title(Pacfin.Port.Name),
         use_this_one = if_else(Pacfin.Port.Name %in% c("Westport", "Sf", "Bellingham", "Neah Bay", "La"), "yes", if_else(any(grep(Pacfin.Port.Name, IOPAC)), "yes", "no"))) %>%
  filter(use_this_one == "yes",
         !Pacfin.Port.Name %in% c("Bragg", "Morro", "Bodega")) %>%
  distinct() %>%
  st_as_sf(coords = c("PACFIN_LONG", "PACFIN_LAT"),
           crs = 4326) %>%
  transform(IOPAC = factor(IOPAC, levels = ports.all))
port_locs$IOPAC <- droplevels(port_locs$IOPAC)


foot_map_all <- ggplot() +
  geom_sf(data = footprints %>% filter(total_unique_vessels > 2,
                                       !port_name %in% c("Los Angeles", "San Diego", "Santa Barbara")), aes(fill = port_name), alpha = 0.8) + ##no overlap and makes the map really difficult to manage
  scale_fill_manual(name = "Port footprints", values = cols) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf_text(data = coast_outline, aes(label = NAME),
               nudge_x = c(-100, -120, -70, -100),
               nudge_y = c(145,  300,  0,   -150),
               size = 2) +
  geom_sf(data = owe_shps, linewidth = 0.5, color = "black", fill = 'transparent') +
  new_scale_fill() +
  geom_sf(data = port_locs, aes(fill = IOPAC), size = 1.5, shape = 22, color = "gray40", show.legend = FALSE) +
  scale_fill_manual(values = cols) +
  facet_grid( ~ subfishery) +
  xlim(280, 750) +
  ylim(3800, 5350) +
  coord_sf() +
  theme(panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.7),
        panel.grid = element_blank(),
        panel.background = element_blank(),
        axis.title = element_blank(),
        axis.text.y = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 45, vjust = 0.5))
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_1_whole_coast.tiff"), width = 6.5, height = 5)
#foot_map_all


######################################
####Figure 2: Overlap of footprints annually (Exposure axis)####
######################################

##load footprints
footprints_yrs <- read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_{cont}_pvcs_yrs_target.shp'))) %>% 
  transform(fish_group = 'DTS') %>% 
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_ndts_{yrSrt}_{yrEnd}_{cont}_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'nDTS')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_mdt_{yrSrt}_{yrEnd}_{cont}_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'MDT')) %>%
  st_transform(st_crs(projection_extent)) %>%
  filter(cont == cont) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT')))

#overlap shapes
intersect_shp <- st_intersection(footprints_yrs %>% filter(area < 100000), owe_shps)

#how many distinct year x port_name x fish_group combinations?
combos <- intersect_shp %>%
  distinct(year, port_name, fish_group)

#port-year combinations where there was no overlap; return all year-port-group combos not in the overlap but in the footprints
no_over <- anti_join(footprints_yrs %>%
                       filter(area < 100000) %>%
                       distinct(year, port_name, fish_group),
                     intersect_shp %>% distinct(year, port_name, fish_group)) %>%
  transform(year = as.integer(year))


intersect_grp <- intersect_shp %>% 
  mutate(intersect_area = st_area(.)) %>%   # calculate area of intersections
  transform(perc = as.numeric(intersect_area/area)) %>%
  dplyr::select(area, port_name, year, fish_group, Call_Area, intersect_area, perc) %>%   
  st_drop_geometry() %>%
  transform(Call_Area = factor(Call_Area, levels = c('Coos Bay', 'Brookings', 'Humboldt', 'Morro Bay'),
                               labels = c('Coos Bay', 'Brookings', 'Humboldt', 'Morro Bay')),
            year = as.integer(year)) %>%
  
  #add in all rows, not just intersection rows
  bind_rows(no_over) %>%
  transform(perc = ifelse(is.na(perc), 0, perc)) #want to identify whether the fishery occurred but that there was no overlap, so change NAs to 0's


#Now want time series plot to have disconnected lines when there is data (years) missing, so need to have every year for each port and fish_group in the dataframe
displ_grp_all_combs <- intersect_grp %>%
  expand("year" = yrSrt:yrEnd, port_name, fish_group)

#sum across WEAs for the total percent of each footprint that is taken up wind
displ_grp_exposure <- intersect_grp %>%
  group_by(year, port_name, fish_group) %>%
  dplyr::summarize(perc = sum(perc, na.rm = T),
                   summed_area = sum(intersect_area, na.rm = T)) %>%
  group_by(port_name, fish_group) %>%
  mutate(years_of_overlap_footprints = sum(perc != 0, na.rm = T)) %>%
  ungroup() %>%
  full_join(displ_grp_all_combs, by = c("year", "port_name", "fish_group")) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  left_join(conf_dat_year_port_fish_group, by = c("year" = "RYEAR", "port_name" = "IOPAC", "fish_group")) %>%
  group_by(port_name, fish_group) %>%
  arrange(port_name, fish_group, year)
displ_grp_exposure <- droplevels(displ_grp_exposure)

#for mean dotted line
displ_grp_avg <- displ_grp_exposure %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(perc = round(mean(perc, na.rm = T), 2)) %>%
  arrange(perc) 

displ_grp_sd <- displ_grp_exposure %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(sd = round(sd(perc, na.rm = T), 2)) %>%
  reshape2::dcast(port_name ~ fish_group, value.var = 'sd')

displ_tab <- displ_grp_avg %>%
  reshape2::dcast(port_name ~ fish_group, value.var = 'perc') %>%
  merge(displ_grp_sd, by ='port_name', suffixes = c('.mean', '.sd'))

#Figure 2: Exposure - %>% filter(total_unique_vessels < 3), #this removes year-subfishery-port group combinations with <3 vessels for confidentiality
if(yrSrt > 2010){
  exp_breaks <- seq(yrSrt, yrEnd, by = 1)
} else{
  exp_breaks <- seq(yrSrt, yrEnd, by = 5)
}

overlap_plot_dat <- displ_grp_exposure %>% 
  mutate(perc = ifelse(total_unique_vessels < 3, NA, perc)) %>%
  filter(years_of_overlap_footprints > 4, #only use fleets in which 5 or more years of overlap with OWEAs occurred
         !all(is.na(years_of_overlap_footprints))) %>%
  mutate(port_name = droplevels(port_name))

overlap_plot <- ggplot(overlap_plot_dat, aes(x = year, y = perc, color = fish_group)) +
  geom_point(size = 1.5) +
  geom_line() +
  geom_hline(data = displ_grp_avg %>%
               filter(port_name %in% unique(overlap_plot_dat$port_name)), aes(yintercept = perc, color = fish_group), linetype = 'dashed') +
  geom_text(data = displ_grp_exposure %>% 
              filter(port_name %in% unique(overlap_plot_dat$port_name),
                     total_unique_vessels < 3), aes(x = year, y = -Inf, color = fish_group, label = "#"), size = 2, vjust = -0.4, show.legend = FALSE) +
  xlab('') + 
  ylab('Footprint overlap with OWE areas') +
  facet_grid(fish_group ~ port_name) +
  theme_bw() +
  theme(legend.position = "none",
        legend.margin = margin(5, 2, 0, 2, unit = "pt"),
        legend.box.background = element_rect(color = "black"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.7),
        panel.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(name = "Fishery", values = cols[c(6,7,9)]) +
  scale_x_continuous(breaks = exp_breaks) +
  scale_y_continuous(limits = c(-0.02, max(displ_grp_exposure$perc))) #change to c(0, 0.6) if using OR Call Areas; 0.2 if using OR WEAs
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_2_exposure.tiff"), width = 6.5, height = 4.5)


######################################
####Table 1: statistics summary####
######################################

summary_table <- conf_dat_port_fish_group %>%
  left_join(displ_grp_exposure %>% dplyr::select("IOPAC" = "port_name", fish_group, years_of_overlap_footprints) %>% filter(!is.na(years_of_overlap_footprints)) %>% distinct()) %>%
  mutate(cont = cont)
write.csv(summary_table, here('Confidential', 'output', glue('{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_summary_stats.csv')))


######################################
####Figure 3: Fishing fidelity (Adaptive Capacity component of Vulnerability axis)####
######################################

#Uses footprint_yrs from above...this takes a while to run
#sum across WEAs for the total percent of each footprint that overlaps with WEAs
displ_grp <- intersect_grp %>%
  group_by(year, port_name, fish_group) %>%
  dplyr::summarize(perc = sum(perc, na.rm = T),
                   summed_area = sum(intersect_area, na.rm = T)) %>%
  group_by(port_name, fish_group) %>%
  mutate(years_of_overlap_footprints = sum(perc != 0, na.rm = T)) %>%
  ungroup() %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  left_join(conf_dat_year_port_fish_group, by = c("year" = "RYEAR", "port_name" = "IOPAC", "fish_group"))
displ_grp <- droplevels(displ_grp) 


fid_all <- data.frame() #storage

sp_grp <- c('DTS', 'nDTS', 'MDT')
for (s in 1:length(sp_grp)) {
  
  #ports with >4 years of footprints overlapping with OWE areas in a fishery
  ports.fid <- displ_grp %>% distinct(port_name)
  
  for (p in 1:dim(ports.fid)[1]) {
    
    #subset footprints for port-species combinations
    foot_simp <- footprints_yrs %>%
      filter(port_name == ports.fid[p,] & fish_group == sp_grp[s])
    
    #find the unique area of each footprint not included in the others combined
    foot_yrs <- unique(foot_simp$year)
    
    yr_grid <- expand_grid(year1 = c(foot_yrs), year2 = c(foot_yrs)) %>%
      filter(year1 != year2) %>%
      mutate(row_id = row_number()) %>%
      pivot_longer(year1:year2) %>%
      arrange(row_id, value) %>%
      group_by(row_id) %>%
      mutate(name = paste0("year", row_number())) %>%
      ungroup() %>%
      pivot_wider(id_cols = row_id, names_from = name, values_from = value) %>%
      dplyr::select(-row_id) %>%
      distinct()
    
    fid_yrs <- data.frame()
    for (i in 1:dim(yr_grid)[1]) {
      try({
        numer <- st_intersection(foot_simp %>% filter(year == as.character(yr_grid[i,1])),
                                 foot_simp %>% filter(year == as.character(yr_grid[i,2])))
        
        denom <- st_union(foot_simp %>% filter(year == as.character(yr_grid[i,1])),
                          foot_simp %>% filter(year == as.character(yr_grid[i,2])))
        
        #look at individual footprints
        ggplot() +
          geom_sf(data = foot_simp %>% filter(year == as.character(yr_grid[i,1])), fill = 'dodgerblue', alpha = 0.5) +
          geom_sf(data = foot_simp %>% filter(year == as.character(yr_grid[i,2])), fill = 'lightblue', alpha = 0.5)
        
        #look at numerator/denominator
        ggplot() +
          geom_sf(data = numer, fill = 'dodgerblue', alpha = 0.5) +
          geom_sf(data = denom, fill = 'lightblue', alpha = 0.5)
        
        
        
        if (dim(numer)[1] == 0) #if footprints don't overlap, fidelity measure == 0
        {fid <- data.frame(fid = 0,
                           port_name = ports.fid[p,],
                           fish_group = sp_grp[s])}
        else
          
        {fid <- data.frame(fid = as.numeric(st_area(numer)/st_area(denom)),
                           port_name = ports.fid[p,],
                           fish_group = sp_grp[s])}
        
        fid_yrs <- bind_rows(fid, fid_yrs)
      })
      
    } # i years
    
    fid_all <- bind_rows(fid_yrs, fid_all)
    
  } #port
} #species group

write.csv(fid_all, file = here('Confidential', 'output', glue('fidelity_{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_target.csv')), row.names = F)


#Can comment out the above fid_all for-loop if these data are already saved in the write.csv step directly above this
fid_all <- read.csv(file = here('Confidential', 'output', glue('fidelity_{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_target.csv')), header = T, stringsAsFactors = F)

fid_dat <- fid_all %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(mean = mean(fid, na.rm = TRUE), 
                   sd = sd(fid, na.rm = TRUE)) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  mutate(port_name = droplevels(port_name))

#let's add the number of years fished to the x-axis of the figure and number of years footprint overlaps with OWEAs
years.fished.labels <- conf_dat_port_fish_group %>%
  st_drop_geometry() %>%
  dplyr::select("port_name" = "IOPAC", fish_group, no_of_years_fished) %>%
  left_join(displ_grp_exposure %>% 
              dplyr::select(port_name, fish_group, years_of_overlap_footprints, total_unique_vessels), by = c("port_name", "fish_group"))

fid_plot_dat <- fid_dat %>%
  left_join(years.fished.labels, by = c("port_name", "fish_group")) %>%
  filter(years_of_overlap_footprints > 4, #only keep fleets that overlapped in at least 5 years of the period
         total_unique_vessels > 2)

fid_plot <- ggplot(fid_plot_dat, aes(x = port_name, y = mean, color = fish_group, shape = fish_group, group = fish_group)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ifelse(mean-sd < 0, 0, mean-sd), ymax = ifelse(mean+sd > 1, 1, mean+sd)), width = 0.3, position = position_dodge(width = 0.5)) +
  xlab('') + 
  ylab('Fishing site fidelity') +
  geom_text(aes(x = port_name, -Inf, label = no_of_years_fished), position = position_dodge(width = 0.5), size = 2, vjust = -0.4, show.legend = FALSE) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  scale_color_manual(values = cols[c(6,7,9)], name = 'Fishery') +
  scale_shape_manual(values = c(19, 17, 15), name = 'Fishery')
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_3_adaptive_capacity.tiff"), width = 6.5, height = 4)


######################################
####Figure 4: logbook derived stats and proportion of targeted species caught within WEAS (Sensitivity component of Vulnerability axis)####
######################################

#for statistics
#all hauls
log_dat_wea_all_combs <- log_dat_wea %>%
  expand(RYEAR, IOPAC, fish_group) 
log_dat_wea_all_avg_combs <- log_dat_wea %>%
  expand(IOPAC, fish_group)

dat_summ_all <- log_dat_wea %>% 
  group_by(RYEAR, IOPAC, fish_group) %>% 
  # total catch
  summarise(all_areas_total_mt = sum(qtykept)/1000,
            all_areas_subfishery_mt = sum(target_species_landings)/1000,
            total_hauls = n(),
            total_trips = n_distinct(TRIP_ID),
            total_unique_vessels = n_distinct(DRVID), .groups = 'keep') %>%
  right_join(log_dat_wea_all_combs, by = c('RYEAR', 'IOPAC', 'fish_group')) %>%
  replace(is.na(.), 0) %>%
  arrange(RYEAR, IOPAC, fish_group)

dat_summ_all_avg <- log_dat_wea %>% 
  group_by(IOPAC, fish_group) %>% 
  # total catch
  summarise(all_areas_total_mt = sum(qtykept)/1000,
            all_areas_subfishery_mt = sum(target_species_landings)/1000,
            total_hauls = n(),
            total_trips = n_distinct(TRIP_ID),
            total_unique_vessels = n_distinct(DRVID), .groups = 'keep') %>%
  right_join(log_dat_wea_all_avg_combs, by = c('IOPAC', 'fish_group')) %>%
  replace(is.na(.), 0) %>%
  arrange(IOPAC, fish_group)

#inside/outside specific WEA proportions
dat_weas_summ <- log_dat_wea %>% 
  # designate NAs as "Outside" wind energy areas
  mutate(wind_area = replace_na(Call_Area, "Outside WEA")) %>% 
  group_by(RYEAR, IOPAC, fish_group, wind_area) %>% 
  summarise(within_areas_total_mt = sum(qtykept)/1000,
            within_areas_subfishery_mt = sum(target_species_landings)/1000,
            trips = n_distinct(TRIP_ID),
            hauls = n(),
            vessels = n_distinct(DRVID), .groups = 'keep') %>% 
  #join the overall data
  left_join(dat_summ_all, by = c("RYEAR", "IOPAC", 'fish_group')) %>% 
  mutate(prop_total_mt = within_areas_total_mt/all_areas_total_mt,
         prop_subfishery_mt = within_areas_subfishery_mt/all_areas_subfishery_mt,
         prop_hauls = hauls/total_hauls,
         prop_trips = trips/total_trips)


#inside outside any WEA portions -- annual; this captures and calculates proportions for all years in which there was no overlap with OWE areas - i.e. retains 0's for mean calculations
dat_summ <- log_dat_wea %>% 
  # designate NAs as "Outside" wind energy areas
  mutate(wind_area = replace_na(Call_Area, "Outside WEA")) %>% 
  transform(out_in = factor(ifelse(wind_area == 'Outside WEA', 'Outside', 'Inside'))) %>%
  group_by(RYEAR, IOPAC, fish_group, out_in) %>% 
  summarise(within_areas_total_mt = sum(qtykept)/1000,
            within_areas_subfishery_mt = sum(target_species_landings)/1000,
            trips = n_distinct(TRIP_ID),
            hauls = n(),
            vessels = n_distinct(DRVID)) %>%
  ungroup()

dat_summ_combs <- dat_summ %>%
  expand(RYEAR, IOPAC, fish_group, out_in)

dat_summ_all_combs <- dat_summ_combs %>%
  left_join(dat_summ, by = c('RYEAR', 'IOPAC', 'fish_group', 'out_in')) %>%
  #join the overall data
  left_join(dat_summ_all, by = c("RYEAR","IOPAC", 'fish_group')) %>%
  replace(is.na(.), 0) %>%
  mutate(prop_total_mt = within_areas_total_mt/all_areas_total_mt,
         prop_subfishery_mt = within_areas_subfishery_mt/all_areas_subfishery_mt,
         prop_hauls = hauls/total_hauls,
         prop_trips = trips/total_trips) %>%
  mutate(across(where(is.numeric), \(x) ifelse(is.nan(x), NA, x))) %>%
  arrange(RYEAR, IOPAC, fish_group) %>%
  filter(out_in == 'Inside',
         all_areas_total_mt != 0) %>%
  group_by(IOPAC, fish_group) %>%
  mutate(no_years_overlap_logbook = sum(prop_total_mt > 0)) %>%
  #filter(no_years_overlap_logbook > 2) %>%
  group_by(IOPAC, fish_group) %>%
  #filter(!all(total_unique_vessels < 3)) %>%
  mutate(IOPAC = droplevels(IOPAC))


dat_summ_mean_sd <- dat_summ_all_combs %>%
  group_by(IOPAC, fish_group, no_years_overlap_logbook) %>%
  summarise(n_years_fished = n(),
            mean_prop_mt = mean(prop_total_mt, na.rm = TRUE),
            sd_prop_mt = sd(prop_total_mt, na.rm = TRUE),
            subfishery_mean_prop_mt = mean(prop_subfishery_mt, na.rm = TRUE),
            subfishery_sd_prop_mt = sd(prop_subfishery_mt, na.rm = TRUE)) %>%
  left_join(displ_grp %>% dplyr::select(port_name, fish_group, years_of_overlap_footprints), by = c("IOPAC" = "port_name", "fish_group")) %>%
  distinct() %>%
  left_join(dat_summ_all_avg %>% dplyr::select(IOPAC, fish_group, "total_unique_vessels_logbook" = "total_unique_vessels"), by = c("IOPAC", "fish_group")) %>%
  filter(!is.na(years_of_overlap_footprints))


#Sensitivity plot
subfishery_prop_mt_plot <- ggplot(dat_summ_mean_sd %>% filter(total_unique_vessels_logbook > 2,
                                                              years_of_overlap_footprints > 4), aes(x = IOPAC, y = subfishery_mean_prop_mt, col = fish_group, shape = fish_group)) + #only keep fleets that overlapped in at least 5 years of the period
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ifelse(subfishery_mean_prop_mt - subfishery_sd_prop_mt < 0, 0, subfishery_mean_prop_mt - subfishery_sd_prop_mt), ymax = subfishery_mean_prop_mt + subfishery_sd_prop_mt), width = 0.3, position = position_dodge(width = 0.5)) +
  geom_text(aes(x = IOPAC, -Inf, label = n_years_fished), position = position_dodge(width = 0.5), size = 2, vjust = -0.4, show.legend = FALSE) +
  ylab('Mean proportion of landings within OWE areas') + xlab('') +
  scale_y_continuous(breaks = c(seq(0,0.6,0.1))) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  scale_color_manual(values = cols[c(6,7,9)], name = 'Fishery') +
  scale_shape_manual(values = c(19, 17, 15), name = 'Fishery')
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_4_sensitivity.tiff"), width = 6.5, height = 4)



######################################
####Figure 5: Risk calculations ####
######################################

#Join Exposure, Sensitivity and Adaptive Capacity values for plotting
exposure.mean <- displ_grp_avg %>%
  rename("exposure_mean" = "perc")
exposure.sd <- displ_grp_sd %>%
  pivot_longer(cols = -port_name) %>%
  rename("fish_group" = "name",
         "exposure_sd" = "value")
exposure <- exposure.mean %>%
  left_join(exposure.sd, by = c("port_name", "fish_group"))

vulnerability <- dat_summ_all_combs %>% #sensitivity
  dplyr::select(RYEAR, IOPAC, fish_group, prop_subfishery_mt) %>%
  rename("port_name" = "IOPAC") %>%
  merge(fid_dat %>% rename(fidelity.mean = mean) %>% dplyr::select(-sd), 
        by = c('port_name', 'fish_group')) %>%
  mutate(sensitivity = prop_subfishery_mt,
         adapt_cap = fidelity.mean,
         vulnerability_euclid = sqrt(sensitivity^2 + adapt_cap^2)) %>%
  group_by(port_name, fish_group) %>%
  summarise(vul_mean_euclid = mean(vulnerability_euclid, na.rm = TRUE),
            vul_sd_euclid = sd(vulnerability_euclid, na.rm = TRUE)) %>%
  left_join(dat_summ_mean_sd, by = c('port_name' = 'IOPAC', 'fish_group')) %>%
  rename('sensitivity' = 'subfishery_mean_prop_mt',
         'sensitivity_sd' = 'subfishery_sd_prop_mt') %>%
  left_join(fid_dat %>% rename("adapt_cap" = "mean", "adapt_cap_sd" = "sd"), by = c('port_name', 'fish_group')) %>%
  left_join(exposure) %>%
  mutate(fish_group = factor(fish_group, levels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl')),
         cont = cont,
         relative_risk = round(sqrt(vul_mean_euclid^2 + exposure_mean^2), 2)) %>%
  dplyr::select(port_name, fish_group, n_years_fished, "n_years_overlap_logbook" = "no_years_overlap_logbook", "n_years_overlap_footprints" = "years_of_overlap_footprints", total_unique_vessels_logbook, exposure_mean, exposure_sd, adapt_cap, adapt_cap_sd, sensitivity, sensitivity_sd, vul_mean_euclid, vul_sd_euclid, mean_prop_mt, sd_prop_mt, cont, relative_risk)
write.csv(vulnerability, glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_risk_values_and_summary_stats.csv"))

risk_plot_vul_x_exp <- ggplot(vulnerability %>% filter(n_years_overlap_footprints > 4,
                                                       n_years_overlap_logbook > 4,
                                                       total_unique_vessels_logbook > 2), aes(exposure_mean, vul_mean_euclid, color = port_name, fill = port_name, shape = port_name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ifelse(vul_mean_euclid - vul_sd_euclid < 0, 0, vul_mean_euclid - vul_sd_euclid), ymax = vul_mean_euclid + vul_sd_euclid)) +
  geom_errorbar(aes(xmin = ifelse(exposure_mean - exposure_sd < 0, 0, exposure_mean - exposure_sd), xmax = exposure_mean + exposure_sd)) +
  xlab('Exposure') + ylab('Vulnerability') +
  facet_grid(~fish_group, drop = T) +
  theme_bw() +
  theme(#panel.grid = element_blank(),
    legend.position = c(0.91, 0.77),
    legend.box.background = element_rect(color = "black"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.spacing = unit(0.45, "cm"),
    legend.key.size = unit(0.45, "cm")) +
  scale_y_continuous(limits = c(0, 1.0)) +
  # scale_x_continuous(limits = c(0, 0.55)) +
  scale_color_manual(values = cols[c(6,7,8,9,10,15)], name = 'Port group') +
  scale_fill_manual(values = cols[c(6,7,8,9,10,15)], name = 'Port group') +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 8), name = 'Port group')
risk_plot_vul_x_exp
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_5_risk_plot.tiff"), width = 6.5, height = 4)



######################################
####SDM Projections and climate change effects ####
######################################

######################################
####Gather data from projections####
######################################

#load all available groundfish CPUE projections - data from Liu et al. 2023 & 2024)#
spp_othr <- c('sablefish', 'lingcod', 'longspine thornyhead', 'pacific grenadier', 'shortspine thornyhead')
othr_dat <- data.frame()
for (s in spp_othr) {
  print(s)
  dat_temp <- readRDS(file = here::here('Non_Confidential', 'Liu_projections', 'Other', s, 'projection_3ESMs.rds')) %>%
    transform(sp = s)
  
  othr_dat <- bind_rows(othr_dat, dat_temp)
}

spp_rock <- c('aurora rockfish', 'chilipepper', 'darkblotched rockfish', 'greenstriped rockfish', 'sharpchin rockfish', 'shortbelly rockfish', 'shortraker rockfish', 'splitnose rockfish', 'stripetail rockfish', 'canary rockfish', 'bocaccio', 'pacific ocean perch', 'widow rockfish', 'yelloweye rockfish', 'yellowtail rockfish')

rock_dat <- data.frame()
for (s in spp_rock) {
  print(s)
  dat_temp <- readRDS(file = here::here('Non_Confidential', 'Liu_projections', 'Rockfish', s, 'projection_3ESMs.rds')) %>%
    transform(sp = s)
  
  rock_dat <- bind_rows(rock_dat, dat_temp)
}

spp_flat <- c('Arrowtooth Flounder', 'Dover Sole', 'Curlfin Sole', 'Pacific Sanddab',
              'Petrale Sole', 'Rex Sole', 'English Sole', 'Slender Sole')
flat_dat <- data.frame()
for (s in spp_flat) {
  print(s)
  dat_temp <- readRDS(file = here::here('Non_Confidential', 'Liu_projections', 'Flatfish', s, 'projection_3ESMs.rds')) %>%
    transform(sp = s)
  
  flat_dat <- bind_rows(flat_dat, dat_temp)
}

all_spp <- bind_rows(flat_dat, rock_dat, othr_dat)

#identify target species for each fishery
dts_trawl <- c('Dover Sole', 'shortspine thornyhead', 'longspine thornyhead', 'sablefish')

ndts_trawl <- c('Arrowtooth Flounder', 'Curlfin Sole', 'Pacific Sanddab', 'Petrale Sole', 'Rex Sole', 'English Sole', 'Slender Sole', 'aurora rockfish', 'chilipepper', 'darkblotched rockfish', 'greenstriped rockfish', 'sharpchin rockfish', 'shortbelly rockfish', 'shortraker rockfish', 'splitnose rockfish', 'stripetail rockfish', 'canary rockfish', 'bocaccio', 'pacific ocean perch', 'yelloweye rockfish', 'lingcod', 'pacific grenadier')

mid_trawl <- c('widow rockfish', 'yellowtail rockfish')

all_spp <- all_spp %>% 
  transform(fishery = ifelse(sp %in% dts_trawl, 'dts_trawl', 
                             ifelse(sp %in% mid_trawl, 'mid_trawl', 'ndts_trawl')))

#Mapping functions#
map_year <- function(spp_preds, yr_vec = yr_vec, return_pred_df = F, plot_leg = T){
  # scale for the legend, common within species
  df <- spp_preds %>%
    filter(year%in%yr_vec) %>%
    group_by(year, sp, longitude,latitude) %>%
    summarise(est = mean(est, na.rm=T, .groups = 'keep') %>% exp()) %>% #averaging across the 3 Earth System Models
    group_by(year, longitude, latitude) %>%
    summarise(est = sum(est, na.rm = T)) %>% #summing over all species in the fishery (i.e. total biomass of targeted species)
    group_by(longitude, latitude) %>%
    summarise(est = mean(est, na.rm = T, .groups = 'keep')) #mean sum across years (i.e. average total biomass of targeted species)
  
  scl <- c(0, quantile(df$est, 0.99)) # rescale super large positive outliers for mapping purposes
  
  df <- df %>%
    mutate(est = ifelse(est > scl[2], scl[2], est))
  
  # match nearest neighbors from predictions to grid
  pred_points <- df %>% 
    dplyr::select(longitude, latitude) %>% 
    as.matrix()
  
  nns <- nn2(pred_points, gr_xy, k=1)$nn.idx
  
  gr_pred <- gr_xy %>% 
    as_tibble() %>% 
    mutate(est=df$est[nns])
  
  bbox = st_bbox(projection_extent)
  
  if(return_pred_df) {out <- gr_pred %>% as_tibble()}
  
  else{
    out<-ggplot(coast_outline) +
      geom_sf() +
      geom_raster(data = gr_pred, aes(x = X, y = Y, fill = est), interpolate = F) +
      scale_fill_viridis_c(limits = scl) +
      coord_sf(datum = NA) +
      xlim(bbox[1], bbox[3]) + ylim(bbox[2], bbox[4])  +
      labs(x = "", y = "", fill = "CPUE", title = '')
    if(!plot_leg) out <- out + theme(legend.position = 'None')
  }
  out
}


######################################
####Figure 6: CPUE maps (example from each fishery)#####
######################################

##DTS trawl predictions 
spp_preds <- all_spp %>% filter(fishery == 'dts_trawl')

basepred <- map_year(spp_preds, yr_vec = 2020, return_pred_df = T)
newpred <- map_year(spp_preds, yr_vec = 2050:2100, return_pred_df = T) %>%  #check times/years
  mutate(est_comp = est) %>% dplyr::select(-est)
both <- basepred %>% 
  left_join(newpred, by=c('X','Y')) %>% 
  mutate(diffpred = est_comp-est) %>%
  mutate(perc = (est_comp-est)/est*100)

bbox = st_bbox(projection_extent)

dts_cpue_map <- ggplot(both) +
  geom_raster(aes(x = X, y = Y, fill = perc), interpolate = F) + 
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red3", limits = c(-100,100)) +
  geom_sf(data = owe_shps, fill = NA, linewidth = 0.7) +
  geom_sf(data = footprints %>% filter(port_name == ifelse(oweas == "ORCAS", 'Eureka', 'Coos Bay') & cont == cont & subfishery == 'DTS trawl'), fill = NA, color = "tomato3", linewidth = 0.7) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf(data = port_locs %>% filter(IOPAC == ifelse(oweas == "ORCAS", 'Eureka', 'Coos Bay')), size = 3, shape = 22, fill = "red3", color = "gray40", show.legend = FALSE) +
  geom_sf_text(data = coast_outline, aes(label = NAME),
               nudge_x = c(-300, -220, -230, -200),
               nudge_y = c(450,  300,  -200,   -150),
               size = 3) +
  xlim(bbox[1],550) +
  ylim(4410,5100) +
  labs(title = ifelse(oweas == "ORCAS", "a) Eureka DTS", "a) Coos Bay DTS"), x = "", y = "", fill = "%\u0394\nbiomass\nindex") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.position = "none",
        #legend.position = c(0.75, 0.4),
        legend.title.align = 0.5, 
        legend.background = element_rect(fill = "white")) +
  coord_sf(expand = FALSE)
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_6_DTS_cpue.tiff"), width = 3, height = 6.5)
dts_cpue_map

##Non-DTS trawl predictions 
spp_preds <- all_spp %>% filter(fishery == 'ndts_trawl')

basepred <- map_year(spp_preds, yr_vec = 2020, return_pred_df = T)
newpred <- map_year(spp_preds, yr_vec = 2050:2100, return_pred_df = T) %>%  #check times/years
  mutate(est_comp = est) %>% dplyr::select(-est)
both <- basepred %>% 
  left_join(newpred, by=c('X','Y')) %>% 
  mutate(diffpred = est_comp-est) %>%
  mutate(perc = (est_comp-est)/est*100)

bbox = st_bbox(projection_extent)
non_dts_cpue_map <- ggplot(both) +
  geom_raster(aes(x = X, y = Y, fill = perc), interpolate = F) + 
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red3", limits = c(-100, 100)) +
  geom_sf(data = owe_shps, fill = NA, linewidth = 0.7) +
  geom_sf(data = footprints %>% filter(port_name == ifelse(oweas == "ORCAS", 'Coos Bay', 'Eureka') & cont == cont & subfishery == 'non-DTS trawl'), fill = NA, color = "tomato3", linewidth = 0.7) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf(data = port_locs %>% filter(IOPAC == ifelse(oweas == "ORCAS", 'Coos Bay', 'Eureka')), size = 3, shape = 22, fill = "red3", color = "gray40", show.legend = FALSE) +
  xlim(bbox[1],550) +
  ylim(4410,5100) +
  labs(title = ifelse(oweas == "ORCAS", "b) Coos Bay non-DTS", "b) Eureka non-DTS"), x = "", y = "", fill = "%\u0394\nbiomass\nindex") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(size = 8, angle = 90, vjust = 0.4),
        #legend.position = c(0.75, 0.4),
        legend.position = "none",
        legend.title.align = 0.5,
        legend.background = element_rect(fill = "white")) +
  coord_sf(expand = FALSE)
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_6_non_dts_cpue.tiff"), width = 3, height = 6.5)
non_dts_cpue_map


##Midwater trawl predictions 
spp_preds <- all_spp %>% filter(fishery == 'mid_trawl')

basepred <- map_year(spp_preds, yr_vec = 2020, return_pred_df = T)
newpred <- map_year(spp_preds, yr_vec = 2050:2100, return_pred_df = T) %>%  #check times/years
  mutate(est_comp = est) %>% dplyr::select(-est)
both <- basepred %>% 
  left_join(newpred, by=c('X','Y')) %>% 
  mutate(diffpred = est_comp-est) %>%
  mutate(perc = (est_comp-est)/est*100)

bbox = st_bbox(projection_extent)
midwater_cpue_map <- ggplot(both) +
  geom_raster(aes(x = X, y = Y, fill = perc), interpolate = F) + 
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red3", limits = c(-100,100)) +
  geom_sf(data = owe_shps, fill = NA, linewidth = 0.7) +
  geom_sf(data = footprints %>% filter(port_name == ifelse(oweas == "ORCAS", 'Brookings', 'Crescent City') & cont == cont & subfishery == 'Midwater trawl'), fill = NA, color = "tomato3", linewidth = 0.7) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf(data = port_locs %>% filter(IOPAC == ifelse(oweas == "ORCAS", 'Brookings', 'Crescent City')), size = 3, shape = 22, fill = "red3", color = "gray40", show.legend = FALSE) +
  xlim(bbox[1],550) +
  ylim(4410,5100) +
  labs(title = ifelse(oweas == 'ORCAS', 'c) Brookings midwater', 'c) Crescent City midwater'), x = "", y = "", fill = "%\u0394\nbiomass\nindex") +
  theme_bw() +
  theme(plot.title = element_text(size = 10),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        axis.text.y = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.position = c(0.75, 0.4),
        legend.title.align = 0.5, 
        legend.background = element_rect(fill = "white")) +
  coord_sf(expand = FALSE)
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_6_midwater_cpue.tiff"), width = 3, height = 6.5)
midwater_cpue_map


dts_cpue_map + non_dts_cpue_map + midwater_cpue_map
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_6_all_three_fisheries_cpue.tiff"), width = 6.5, height = 5.5)



######################################
####Figure 7: CPUE inside and outside of WEAs within historical footprints#####
######################################
# CPUE inside/outside WEAs
# Logic- do a sequential join of footprints and WEAs, keeping track of:
# inside vs. outside footprints
# inside vs. outside WEAs
# then just do algebra

# first, summed CPUE per fishery (across species) per year and grid cell
cpue_gridcell_fishery_year <- all_spp %>%
  mutate(est = exp(est)) %>%
  group_by(year, fishery, sp, latitude, longitude) %>%
  dplyr::summarize(est = mean(est, na.rm = T), .groups = 'keep') %>% # this should summarize across the 3 ESMs
  group_by(year, fishery, longitude, latitude) %>%
  summarise(est = sum(est, na.rm = T)) %>% # summed CPUE across species per fishery
  ungroup()

sdm_id <- cpue_gridcell_fishery_year %>%
  dplyr::select(longitude, latitude) %>%
  distinct() %>%
  mutate(sdm_grid_id = row_number())


# now, join footprints from period of interest
cpue_inside_outside_fps <- cpue_gridcell_fishery_year %>% 
  left_join(sdm_id) %>% 
  st_as_sf(coords=c("longitude","latitude"), crs=st_crs(footprints)) %>% 
  st_join(footprints_yrs %>%
            mutate(fish_group = ifelse(fish_group == "DTS", "dts_trawl", ifelse(fish_group == "nDTS", "ndts_trawl", "mid_trawl")),
                   year = as.numeric(year)) %>% #get the two fishery categories into the same values
            filter(port_name %in% ports.all)) %>% 
  mutate(inside_fp = ifelse(!is.na(port_name),"Yes","No")) %>% # indicates whether or not a grid cell is within a footprint
  filter(fishery == fish_group) #need to make sure only footprints for the specific fishery species is kept

# for this analysis, we only care about values inside footprints (and inside vs. outside WEAs) and fleets that overlapped > 4 years
cpue_inside_fps <- cpue_inside_outside_fps %>% 
  filter(inside_fp == "Yes") 

# now, we join these to the unioned WEAs (we don't care about the identity/name of the WEA here)
owe_shps_union <- owe_shps_union %>% st_as_sf() %>% mutate(is_wea = "Yes")

cpue_fps_inside_outside_weas <- cpue_inside_fps %>% 
  st_join(owe_shps_union %>% st_as_sf) %>% 
  mutate(is_wea = ifelse(is.na(is_wea), "No", is_wea)) # just like the footprints above, we indicate whether a point is inside or outside WEAs

#count number of SDM points inside and outside WEAs
cpue_fps_inside_outside_weas_counts <- cpue_fps_inside_outside_weas %>%
  group_by(fishery, port_name, is_wea) %>%
  summarise(SDM_points = length(unique(sdm_grid_id)))

sdm_points_plot <- ggplot(cpue_fps_inside_outside_weas_counts) +
  geom_bar(aes(SDM_points)) +
  facet_grid(fishery ~ port_name) +
  ylim(0, 25)

cpue_fps_inside_outside_weas_edit <- cpue_fps_inside_outside_weas %>%
  left_join(cpue_fps_inside_outside_weas_counts %>% 
              st_drop_geometry() %>%
              filter(is_wea == "Yes") %>%
              dplyr::select(fishery, port_name, SDM_points), by = c("fishery", "port_name")) %>%
  filter(SDM_points > 1) %>%
  left_join(displ_grp %>% dplyr::select(port_name, fish_group, years_of_overlap_footprints) %>% mutate(fish_group = ifelse(fish_group == "DTS trawl", "dts_trawl", ifelse(fish_group == "non-DTS trawl", "ndts_trawl", "mid_trawl"))) %>% distinct(), by = c("port_name", "fish_group")) %>%
  filter(!is.na(years_of_overlap_footprints))

# make time series plot of biomass index (Fig. 7 - Biomass index is CPUE (in kg/km2))
#new names for subfishery facet text
subfishery.labs <- c("DTS trawl", "non-DTS trawl", "Midwater trawl")
names(subfishery.labs) <- c("dts_trawl", "ndts_trawl", "mid_trawl")

cpue_time_series <- cpue_fps_inside_outside_weas_edit %>% 
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(cpue_fps_inside_outside_weas_edit) %>%
  st_drop_geometry() %>% 
  dplyr::select(-c(year.y, area)) %>% # get single CPUE values for each set of coordinates by year/fishery/port/inside vs. outside wea
  distinct() %>% 
  st_drop_geometry() %>% 
  group_by(year.x, fishery, port_name, is_wea) %>% # take a mean CPUE by year/fishery/port/inside vs. outside wea
  summarise(mean_cpue = mean(est, na.rm = T),
            median_cpue = median(est, na.rm = T),
            sd_cpue = sd(est, na.rm = T),
            se_cpue = sd_cpue/sqrt(n())) %>%
  group_by(fishery, port_name, is_wea) %>% 
  mutate(mean_5year = slide_dbl(mean_cpue, mean, .before = 2, .after = 2),
         se_5year = slide_dbl(se_cpue, mean, .before = 2, .after = 2),
         sd_5year = slide_dbl(sd_cpue, mean, .before = 2, .after = 2)) %>%
  ungroup() %>%
  mutate(port_name = factor(port_name, levels = ports.all),
         fishery = factor(fishery, levels = c('dts_trawl', 'ndts_trawl', 'mid_trawl')),
         port_name = fct_drop(port_name))

cpue_time_series_fig <- cpue_time_series %>%
  filter(year.x >= 2020) %>%
  ggplot(aes(x = year.x, y = mean_5year, color = is_wea)) +
  geom_point(size = 0.2) +
  geom_line(linewidth = 0.25) +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = F) +
  geom_ribbon(aes(ymin = ifelse(mean_5year-se_5year < 0, 0, mean_5year-se_5year), ymax = mean_5year+se_5year), alpha = 0.2, linewidth = 0.1) +
  scale_color_discrete(labels = c("Outside OWEAs", "Inside OWEAs")) +
  facet_grid(fishery ~ port_name, scales = 'free',
             labeller = labeller(fishery = subfishery.labs)) +
  #facetted_pos_scales(y = list(NULL, NULL, scale_y_continuous(limits = c(-200, 200)))) +
  theme_bw() +
  theme(legend.position = c(0.85, 0.2),
        legend.box.background = element_rect(colour = "black"),
        #legend.key = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  labs(x = 'Year', y = 'Biomass index', col = "Portion of footprints")
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_7_time_series_mean.tiff"), width = 6.5, height = 6)



######################################
####Figure 8: CPUE ratio of inside vs outside of WEAs within historical footprints#####
######################################
cpue_ratio_fig <- cpue_time_series %>%
  pivot_wider(names_from = is_wea, values_from = c(median_cpue, mean_cpue, sd_cpue, se_cpue, mean_5year, sd_5year, se_5year)) %>%
  mutate(ratio = mean_5year_Yes/mean_5year_No,
         sd_ratio = ratio * sqrt((sd_5year_Yes/mean_5year_Yes)^2 + (sd_5year_No/mean_5year_No)^2),
         se_ratio = ratio * sqrt((se_5year_Yes/mean_5year_Yes)^2 + (se_5year_No/mean_5year_No)^2)) %>%
  filter(year.x >= 2020) %>%
  ggplot(aes(x = year.x, y = ratio)) +
  geom_line() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = F) +
  geom_ribbon(aes(ymin = ifelse(ratio-se_ratio < 0, 0, ratio-se_ratio), ymax = ratio+se_ratio), alpha = 0.2, linewidth = 0.1) +
  geom_hline(yintercept = 1, linetype=2) +
  facet_grid(fishery ~ port_name,
             labeller = labeller(fishery = subfishery.labs)) +
  theme_bw() +
  theme(#legend.position = 'none',
    #legend.key = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab('Year') + 
  ylab('Ratio of biomass index')
ggsave(glue("Confidential/output/{oweas}_{cont}PVC_{yrSrt}_{yrEnd}_Figure_8_time_series_ratio.tiff"), width = 6.5, height = 6)





######################################
#### Appendix 1: Sensitivity to choice of kernel density percent volume contour####
######################################

######################################
#### Appendix 1: Figure A1 - Example footprints map across PVC contours####
######################################

footprints_all_pvcs <- read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_95_pvcs_yrs_target.shp'))) %>% 
  transform(cont = factor(cont)) %>% 
  bind_rows(read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_75_pvcs_yrs_target.shp')))) %>% 
  transform(cont = factor(cont)) %>% 
  bind_rows(read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_50_pvcs_yrs_target.shp')))) %>% 
  transform(cont = factor(cont)) %>% 
  st_transform(st_crs(projection_extent))

bb <- st_bbox(footprints_all_pvcs)

#faceted
map_contours <- ggplot() +
  geom_sf(data = coast_outline, fill = 'grey90')+
  geom_sf(data = owe_shps, aes(linetype = 'linetype'), size = 0.5, color = 'gray40', fill = 'transparent') +
  
  geom_sf(data = footprints_all_pvcs %>% filter(port_name == 'Coos Bay' & year %in% c(seq(1994, 2020, 2))), 
          aes(fill = cont, col = cont), alpha = 0.2) +
  geom_sf(data = port_locs %>% filter(IOPAC == 'Coos Bay'), size = 1.5, shape = 22, fill = 'red3', color = "gray40", show.legend = FALSE) +
  facet_wrap(.~year) +
  ggtitle('Fishing footprints for DTS landings to Coos Bay, OR') +
  #whole coast
  # xlim(bb[1],bb[3]) +
  # ylim(bb[2],bb[4]) +
  
  #OR/northern CA
  xlim(100,750) +
  ylim(4500, 5050) +
  
  #central CA
  # xlim(bbox[1],750) +
  # ylim(3600,4200) +
  coord_sf(datum = NA) +
  
  theme(legend.position = 'top',
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.7),
        panel.grid = element_blank(),
        panel.background = element_blank()) +
  scale_linetype_manual(values = c('linetype' = 'solid'), labels = 'OWEAs', name = '') +
  scale_color_manual(values = c('red', 'black', '#f6b61c'), name = 'PVC (%)') +
  scale_fill_manual(values = c('red', 'black', '#f6b61c'), name = 'PVC (%)')
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A1_footprints.tiff"), width = 6.5, height = 6.5)




######################################
#### Appendix 1: Figures A2 and A3 - Example of Exposure across PVC contours (area and percentages)####
######################################

##load footprints
footprints_yrs_app <- read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_50_pvcs_yrs_target.shp'))) %>% 
  transform(fish_group = 'DTS') %>% 
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_ndts_{yrSrt}_{yrEnd}_50_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'nDTS')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_mdt_{yrSrt}_{yrEnd}_50_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'MDT')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_75_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'DTS')) %>% 
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_ndts_{yrSrt}_{yrEnd}_75_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'nDTS')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_mdt_{yrSrt}_{yrEnd}_75_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'MDT')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_dts_{yrSrt}_{yrEnd}_95_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'DTS')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_ndts_{yrSrt}_{yrEnd}_95_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'nDTS')) %>%
  bind_rows(
    read_sf(here('Confidential', 'data', 'PVC', glue('logbook_mdt_{yrSrt}_{yrEnd}_95_pvcs_yrs_target.shp'))) %>% 
      transform(fish_group = 'MDT')) %>%
  st_transform(st_crs(projection_extent))

#overlap shapes
intersect_shp_app <- st_intersection(footprints_yrs_app %>% filter(area < 100000), owe_shps)

#port-year combinations where there was no overlap; return all year-port-group combos not in the overlap but in the footprints
no_over_app <- anti_join(footprints_yrs_app %>%
                           filter(area < 100000) %>%
                           distinct(year, port_name, fish_group, cont),
                         intersect_shp %>% distinct(year, port_name, fish_group, cont)) %>%
  transform(year = as.integer(year))


intersect_grp_app <- intersect_shp_app %>% 
  mutate(intersect_area = st_area(.)) %>%   # calculate area of intersections
  transform(perc = as.numeric(intersect_area/area)) %>%
  dplyr::select(area, port_name, year, cont, fish_group, Call_Area, intersect_area, perc) %>%   
  st_drop_geometry() %>%
  transform(Call_Area = factor(Call_Area, levels = c('Coos Bay', 'Brookings', 'Humboldt', 'Morro Bay'),
                               labels = c('Coos Bay', 'Brookings', 'Humboldt', 'Morro Bay')),
            year = as.integer(year)) %>%
  
  #add in all rows, not just intersection rows
  bind_rows(no_over_app) %>%
  transform(perc = ifelse(is.na(perc), 0, perc)) #want to identify whether the fishery occurred but that there was no overlap, so change NAs to 0's


#sum across WEAs for the total percent of each footprint that is taken up wind
over_area_app <- intersect_grp_app %>%
  group_by(year, port_name, fish_group, cont) %>%
  dplyr::summarize(perc = sum(perc, na.rm = T),
                   summed_area = as.vector(sum(intersect_area, na.rm = T))) %>%
  group_by(port_name, fish_group, cont) %>%
  mutate(years_of_overlap_footprints = sum(perc != 0, na.rm = T)) %>%
  filter(port_name %in% c("Newport", "Coos Bay", "Brookings", "Crescent City", "Eureka", "Morro Bay")) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  left_join(conf_dat_year_port_fish_group, by = c("year" = "RYEAR", "port_name" = "IOPAC", "fish_group")) %>%
  # group_by(port_name, fish_group) %>%
  filter(years_of_overlap_footprints > 4) %>% #only keep fleets that overlapped in at least 5 years of the period
  arrange(port_name, fish_group, year)
over_area_app <- droplevels(over_area_app)

#for mean dotted line
over_area_avg_app <- over_area_app %>%
  group_by(port_name, fish_group, cont) %>%
  dplyr::summarize(perc_mean = round(mean(perc, na.rm = T), 2),
                   perc_sd = round(sd(perc, na.rm = T), 2),
                   summed_area_mean = mean(summed_area, na.rm = T),
                   summed_area_sd = sd(summed_area, na.rm = T),
                   .groups = 'keep')


#average displacement over the study period
#overlap based on square kilometers overlap
overlap_plot_sqkm_app <- ggplot(over_area_app %>% mutate(summed_area = ifelse(total_unique_vessels < 3, NA, summed_area)), aes(x = year, y = summed_area, color = cont)) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(data = over_area_avg_app, aes(yintercept = summed_area_mean, color = cont), linewidth = 0.3, linetype = 'dashed') +
  geom_text(data = over_area_app %>% filter(total_unique_vessels < 3), aes(x = year, y = -Inf, color = cont, label = "#"), size = 2, vjust = -0.4, show.legend = FALSE) +
  xlab('') + 
  ylab('Footprint overlap with OWE areas (km^2)') +
  facet_grid(fish_group ~ port_name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top",
        legend.margin = margin(5, 2, 0, 2, unit = "pt"),
        legend.box.background = element_rect(color = "black"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.7),
        panel.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(name = "PVC (%)", values = c('red', 'black', '#f6b61c')) +
  scale_x_continuous(breaks = exp_breaks) +
  scale_y_continuous(limits = c(-0.05, max(over_area_app$summed_area)))
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A2_overlap_sqkm.tiff"), width = 6.5, height = 4.5)


#overlap based on percent overlap
overlap_plot_percent_app <- ggplot(over_area_app %>% mutate(perc = ifelse(total_unique_vessels < 3, NA, perc)), aes(x = year, y = perc, color = cont)) +
  geom_point(size = 0.5) +
  geom_line() +
  geom_hline(data = over_area_avg_app, aes(yintercept = perc_mean, color = cont), linewidth = 0.3, linetype = 'dashed') +
  geom_text(data = over_area_app %>% filter(total_unique_vessels < 3), aes(x = year, y = -Inf, color = cont, label = "#"), size = 2, vjust = -0.4, show.legend = FALSE) +
  xlab('') + 
  ylab('Footprint overlap with OWE areas') +
  facet_grid(fish_group ~ port_name, scales = "free_y") +
  theme_bw() +
  theme(legend.position = "top",
        legend.margin = margin(5, 2, 0, 2, unit = "pt"),
        legend.box.background = element_rect(color = "black"),
        legend.title = element_text(size = 7),
        legend.text = element_text(size = 6),
        axis.text.x = element_text(size = 6, angle = 90, vjust = 0.5),
        panel.border = element_rect(color = 'black', fill = NA, linewidth = 0.7),
        panel.background = element_blank(),
        panel.grid.minor.x = element_blank(),
        legend.key = element_blank()) +
  scale_color_manual(name = "PVC (%)", values = c('red', 'black', '#f6b61c')) +
  scale_x_continuous(breaks = exp_breaks) +
  scale_y_continuous(limits = c(-0.05, 1))
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A3_overlap_percent.tiff"), width = 6.5, height = 4.5)




######################################
#### Appendix 1: Figure A4 - Example of Fishing Fidelity across PVC contours (Adaptive Capacity comparisons)####
######################################
contours <- c('50', '95')
for (c in 1:length(contours)) {
  fid_all <- data.frame() #storage
  
  sp_grp <- c('DTS', 'nDTS', 'MDT')
  for (s in 1:length(sp_grp)) {
    
    #ports with >4 years of footprints overlapping with OWE areas in a fishery
    ports.fid <- displ_grp %>% distinct(port_name)
    
    for (p in 1:dim(ports.fid)[1]) {
      
      #subset footprints for port-species combinations
      foot_simp <- footprints_yrs_app %>%
        filter(port_name == ports.fid[p,] & cont == contours[c] & fish_group == sp_grp[s])
      
      #find the unique area of each footprint not included in the others combined
      foot_yrs <- unique(foot_simp$year)
      
      yr_grid <- expand_grid(year1 = c(foot_yrs), year2 = c(foot_yrs)) %>%
        filter(year1 != year2) %>%
        mutate(row_id = row_number()) %>%
        pivot_longer(year1:year2) %>%
        arrange(row_id, value) %>%
        group_by(row_id) %>%
        mutate(name = paste0("year", row_number())) %>%
        ungroup() %>%
        pivot_wider(id_cols = row_id, names_from = name, values_from = value) %>%
        dplyr::select(-row_id) %>%
        distinct()
      
      fid_yrs <- data.frame()
      for (i in 1:dim(yr_grid)[1]) {
        try({
          numer <- st_intersection(foot_simp %>% filter(year == as.character(yr_grid[i,1])),
                                   foot_simp %>% filter(year == as.character(yr_grid[i,2])))
          
          denom <- st_union(foot_simp %>% filter(year == as.character(yr_grid[i,1])),
                            foot_simp %>% filter(year == as.character(yr_grid[i,2])))
          
          #look at individual footprints
          ggplot() +
            geom_sf(data = foot_simp %>% filter(year == as.character(yr_grid[i,1])), fill = 'dodgerblue', alpha = 0.5) +
            geom_sf(data = foot_simp %>% filter(year == as.character(yr_grid[i,2])), fill = 'lightblue', alpha = 0.5)
          
          #look at numerator/denominator
          ggplot() +
            geom_sf(data = numer, fill = 'dodgerblue', alpha = 0.5) +
            geom_sf(data = denom, fill = 'lightblue', alpha = 0.5)
          
          
          
          if (dim(numer)[1] == 0) #if footprints don't overlap, fidelity measure == 0
          {fid <- data.frame(fid = 0,
                             port_name = ports.fid[p,],
                             fish_group = sp_grp[s])}
          else
            
          {fid <- data.frame(fid = as.numeric(st_area(numer)/st_area(denom)),
                             port_name = ports.fid[p,],
                             fish_group = sp_grp[s])}
          
          fid_yrs <- bind_rows(fid, fid_yrs)
        })
        
      } # i years
      
      fid_all <- bind_rows(fid_yrs, fid_all)
      
    } #port
  } #species group
  
  write.csv(fid_all, file = here('Confidential', 'output', glue('fidelity_{oweas}_{contours[c]}PVC_{yrSrt}_{yrEnd}_target.csv')), row.names = F)
  
} #contours

fid_cont_all <- read.csv(file = here('Confidential', 'output', glue('fidelity_{oweas}_50PVC_{yrSrt}_{yrEnd}_target.csv')), header = T, stringsAsFactors = F) %>%
  mutate(cont = factor("50")) %>%
  filter(port_name %in% c("Newport", "Coos Bay", "Brookings", "Crescent City", "Eureka", "Morro Bay")) %>%
  bind_rows(read.csv(file = here('Confidential', 'output', glue('fidelity_{oweas}_75PVC_{yrSrt}_{yrEnd}_target.csv')), header = T, stringsAsFactors = F) %>%
              mutate(cont = factor("75"))) %>%
  bind_rows(read.csv(file = here('Confidential', 'output', glue('fidelity_{oweas}_95PVC_{yrSrt}_{yrEnd}_target.csv')), header = T, stringsAsFactors = F) %>%
              mutate(cont = factor("95")))

fid_cont_dat <- fid_cont_all %>%
  group_by(port_name, fish_group, cont) %>%
  dplyr::summarize(mean = mean(fid, na.rm = TRUE), 
                   sd = sd(fid, na.rm = TRUE)) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl')))
fid_cont_dat <- droplevels(fid_cont_dat)

#let's add the number of years fished to the x-axis of the figure
years.fished.labels <- conf_dat_port_fish_group %>%
  st_drop_geometry() %>%
  dplyr::select("port_name" = "IOPAC", fish_group, no_of_years_fished)

fid_cont_dat <- fid_cont_dat %>%
  left_join(years.fished.labels, by = c("port_name", "fish_group"))

fid_cont_plot <- ggplot(fid_cont_dat %>% filter(port_name %in% c("Newport", "Coos Bay", "Brookings", "Crescent City", "Eureka", "Morro Bay")), aes(x = port_name, y = mean, color = cont, shape = cont, group = cont)) +
  geom_point(size = 2, position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = ifelse(mean-sd < 0, 0, mean-sd), ymax = ifelse(mean+sd > 1, 1, mean+sd)), width = 0.3, position = position_dodge(width = 0.5)) +
  xlab('') + 
  ylab('Fishing site fidelity') +
  geom_text(aes(x = port_name, -Inf, label = no_of_years_fished), color = 'black', size = 2, vjust = -0.4, show.legend = FALSE) +
  facet_grid(fish_group ~ .) +
  theme_bw() +
  theme(legend.position = 'top',
        legend.key = element_blank(),
        axis.text.x = element_text(angle = 0, vjust = 0.5)) +
  scale_color_manual(values = c('red', 'black', '#f6b61c'), name = 'PVC (%)') +
  scale_shape_manual(values = c(19, 17, 15), name = 'PVC (%)') +
  scale_y_continuous(limits = c(-0.05, max(fid_cont_dat$mean + fid_cont_dat$sd)))
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A4_adaptive_capacity.tiff"), width = 6.5, height = 4)




######################################
#### Appendix 1: Figure A5 - Risk values across PVC contours####
######################################
exposure_all <- exposure %>% #grab the main text exposure values for the 75 PVC
  mutate(contour = "75") #will add the other two PVC to this file

contours <- c('50', '95')
for (c in 1:length(contours)) {
  exposure_all <- exposure_all
  
#overlap shapes
intersect_shp_all <- st_intersection(footprints_yrs_app %>% filter(cont == contours[c] & area < 100000), owe_shps)

#how many distinct year x port_name x fish_group combinations?
combos_all <- intersect_shp_all %>%
  distinct(year, port_name, fish_group)

#port-year combinations where there was no overlap; return all year-port-group combos not in the overlap but in the footprints
no_over_all <- anti_join(footprints_yrs_app %>%
                       filter(area < 100000) %>%
                       distinct(year, port_name, fish_group),
                     intersect_shp_all %>% distinct(year, port_name, fish_group)) %>%
  transform(year = as.integer(year))


intersect_grp_all <- intersect_shp_all %>% 
  mutate(intersect_area = st_area(.)) %>%   # calculate area of intersections
  transform(perc = as.numeric(intersect_area/area)) %>%
  dplyr::select(area, port_name, year, fish_group, Call_Area, intersect_area, perc) %>%   
  st_drop_geometry() %>%
  transform(Call_Area = factor(Call_Area, levels = c('Coos Bay', 'Brookings', 'Humboldt', 'Morro Bay'),
                               labels = c('Coos Bay', 'Brookings', 'Humboldt', 'Morro Bay')),
            year = as.integer(year)) %>%
  
  #add in all rows, not just intersection rows
  bind_rows(no_over_all) %>%
  transform(perc = ifelse(is.na(perc), 0, perc)) #want to identify whether the fishery occurred but that there was no overlap, so change NAs to 0's


#Now want time series plot to have disconnected lines when there is data (years) missing, so need to have every year for each port and fish_group in the dataframe
displ_grp_all_combs_all <- intersect_grp_all %>%
  expand("year" = yrSrt:yrEnd, port_name, fish_group)

#sum across WEAs for the total percent of each footprint that is taken up wind
displ_grp_exposure_all <- intersect_grp_all %>%
  group_by(year, port_name, fish_group) %>%
  dplyr::summarize(perc = sum(perc, na.rm = T),
                   summed_area = sum(intersect_area, na.rm = T)) %>%
  group_by(port_name, fish_group) %>%
  mutate(years_of_overlap_footprints = sum(perc != 0, na.rm = T)) %>%
  ungroup() %>%
  full_join(displ_grp_all_combs_all, by = c("year", "port_name", "fish_group")) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  left_join(conf_dat_year_port_fish_group, by = c("year" = "RYEAR", "port_name" = "IOPAC", "fish_group")) %>%
  group_by(port_name, fish_group) %>%
  arrange(port_name, fish_group, year)
displ_grp_exposure <- droplevels(displ_grp_exposure)

#for mean dotted line
displ_grp_avg_all <- displ_grp_exposure_all %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(perc = round(mean(perc, na.rm = T), 2)) %>%
  arrange(perc) 

displ_grp_sd_all <- displ_grp_exposure_all %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(sd = round(sd(perc, na.rm = T), 2)) %>%
  reshape2::dcast(port_name ~ fish_group, value.var = 'sd')

#Join Exposure, Sensitivity and Adaptive Capacity values for plotting
exposure.mean_all <- displ_grp_avg_all %>%
  rename("exposure_mean" = "perc")
exposure.sd_all <- displ_grp_sd %>%
  pivot_longer(cols = -port_name) %>%
  rename("fish_group" = "name",
         "exposure_sd" = "value")

exposure_contours <- exposure.mean_all %>%
  left_join(exposure.sd_all, by = c("port_name", "fish_group")) %>%
  mutate(contour = contours[c])

exposure_all <- bind_rows(exposure_contours, exposure_all)

}


fid_dat_50 <- read.csv(file = here('Confidential', 'output', glue('fidelity_{oweas}_50PVC_{yrSrt}_{yrEnd}_target.csv')), header = T, stringsAsFactors = F) %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(mean = mean(fid, na.rm = TRUE), 
                   sd = sd(fid, na.rm = TRUE)) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  mutate(port_name = droplevels(port_name))

fid_dat_95 <- read.csv(file = here('Confidential', 'output', glue('fidelity_{oweas}_95PVC_{yrSrt}_{yrEnd}_target.csv')), header = T, stringsAsFactors = F) %>%
  group_by(port_name, fish_group) %>%
  dplyr::summarize(mean = mean(fid, na.rm = TRUE), 
                   sd = sd(fid, na.rm = TRUE)) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS', 'nDTS', 'MDT'),
                                labels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl'))) %>%
  mutate(port_name = droplevels(port_name))

vulnerability_50 <- dat_summ_all_combs %>% #sensitivity
  dplyr::select(RYEAR, IOPAC, fish_group, prop_subfishery_mt) %>%
  rename("port_name" = "IOPAC") %>%
  merge(fid_dat_50 %>% rename(fidelity.mean = mean) %>% dplyr::select(-sd), 
        by = c('port_name', 'fish_group')) %>%
  mutate(sensitivity = prop_subfishery_mt,
         adapt_cap = fidelity.mean,
         vulnerability_euclid = sqrt(sensitivity^2 + adapt_cap^2)) %>%
  group_by(port_name, fish_group) %>%
  summarise(vul_mean_euclid = mean(vulnerability_euclid, na.rm = TRUE),
            vul_sd_euclid = sd(vulnerability_euclid, na.rm = TRUE)) %>%
  left_join(dat_summ_mean_sd, by = c('port_name' = 'IOPAC', 'fish_group')) %>%
  rename('sensitivity' = 'subfishery_mean_prop_mt',
         'sensitivity_sd' = 'subfishery_sd_prop_mt') %>%
  left_join(fid_dat_50 %>% rename("adapt_cap" = "mean", "adapt_cap_sd" = "sd"), by = c('port_name', 'fish_group')) %>%
  left_join(exposure_all %>% filter(contour == "50")) %>%
  mutate(fish_group = factor(fish_group, levels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl')),
         cont = "50",
         relative_risk = round(sqrt(vul_mean_euclid^2 + exposure_mean^2), 2)) %>%
  dplyr::select(port_name, fish_group, n_years_fished, "n_years_overlap_logbook" = "no_years_overlap_logbook", "n_years_overlap_footprints" = "years_of_overlap_footprints", total_unique_vessels_logbook, exposure_mean, exposure_sd, adapt_cap, adapt_cap_sd, sensitivity, sensitivity_sd, vul_mean_euclid, vul_sd_euclid, mean_prop_mt, sd_prop_mt, cont, relative_risk)
write.csv(vulnerability_50, glue("Confidential/output/{oweas}_50PVC_{yrSrt}_{yrEnd}_risk_values_and_summary_stats.csv"))

vulnerability_95 <- dat_summ_all_combs %>% #sensitivity
  dplyr::select(RYEAR, IOPAC, fish_group, prop_subfishery_mt) %>%
  rename("port_name" = "IOPAC") %>%
  merge(fid_dat_95 %>% rename(fidelity.mean = mean) %>% dplyr::select(-sd), 
        by = c('port_name', 'fish_group')) %>%
  mutate(sensitivity = prop_subfishery_mt,
         adapt_cap = fidelity.mean,
         vulnerability_euclid = sqrt(sensitivity^2 + adapt_cap^2)) %>%
  group_by(port_name, fish_group) %>%
  summarise(vul_mean_euclid = mean(vulnerability_euclid, na.rm = TRUE),
            vul_sd_euclid = sd(vulnerability_euclid, na.rm = TRUE)) %>%
  left_join(dat_summ_mean_sd, by = c('port_name' = 'IOPAC', 'fish_group')) %>%
  rename('sensitivity' = 'subfishery_mean_prop_mt',
         'sensitivity_sd' = 'subfishery_sd_prop_mt') %>%
  left_join(fid_dat_95 %>% rename("adapt_cap" = "mean", "adapt_cap_sd" = "sd"), by = c('port_name', 'fish_group')) %>%
  left_join(exposure_all %>% filter(contour == "95")) %>%
  mutate(fish_group = factor(fish_group, levels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl')),
         cont = "95",
         relative_risk = round(sqrt(vul_mean_euclid^2 + exposure_mean^2), 2)) %>%
  dplyr::select(port_name, fish_group, n_years_fished, "n_years_overlap_logbook" = "no_years_overlap_logbook", "n_years_overlap_footprints" = "years_of_overlap_footprints", total_unique_vessels_logbook, exposure_mean, exposure_sd, adapt_cap, adapt_cap_sd, sensitivity, sensitivity_sd, vul_mean_euclid, vul_sd_euclid, mean_prop_mt, sd_prop_mt, cont, relative_risk)
write.csv(vulnerability_95, glue("Confidential/output/{oweas}_95PVC_{yrSrt}_{yrEnd}_risk_values_and_summary_stats.csv"))

risk_values <- read.csv(here('Confidential', 'output', 'ORCAS_95PVC_1994_2020_risk_values_and_summary_stats.csv')) %>%
  bind_rows(read.csv(here('Confidential', 'output', 'ORCAS_75PVC_1994_2020_risk_values_and_summary_stats.csv'))) %>%
  bind_rows(read.csv(here('Confidential', 'output', 'ORCAS_50PVC_1994_2020_risk_values_and_summary_stats.csv'))) %>%
  transform(port_name = factor(port_name, levels = ports.all)) %>%
  transform(fish_group = factor(fish_group, levels = c('DTS trawl', 'non-DTS trawl', 'Midwater trawl')))

risk_plot_contours <- ggplot(risk_values %>% filter(n_years_overlap_footprints > 4,
                                                    n_years_overlap_logbook > 4,
                                                    total_unique_vessels_logbook > 2), aes(exposure_mean, vul_mean_euclid, color = port_name, fill = port_name, shape = port_name)) +
  geom_point(size = 3) +
  geom_errorbar(aes(ymin = ifelse(vul_mean_euclid - vul_sd_euclid < 0, 0, vul_mean_euclid - vul_sd_euclid), ymax = vul_mean_euclid + vul_sd_euclid)) +
  geom_errorbar(aes(xmin = ifelse(exposure_mean - exposure_sd < 0, 0, exposure_mean - exposure_sd), xmax = exposure_mean + exposure_sd)) +
  xlab('Exposure') + ylab('Vulnerability') +
  facet_grid(cont ~ fish_group, drop = T) +
  theme_bw() +
  theme(#panel.grid = element_blank(),
    legend.position = c(0.91, 0.75),
    legend.box.background = element_rect(color = "black"),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 6),
    legend.spacing = unit(0.5, "cm"),
    legend.key.size = unit(0.5, "cm")) +
  scale_y_continuous(limits = c(0, 1.0)) +
  #scale_x_continuous(limits = c(0, 0.75)) +
  scale_color_manual(values = cols[c(6,7,8,9,10,15)], name = 'Port group') +
  scale_fill_manual(values = cols[c(6,7,8,9,10,15)], name = 'Port group') +
  scale_shape_manual(values = c(21, 22, 23, 24, 25, 8), name = 'Port group')
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A5_risk_plots.tiff"), width = 6.5, height = 4)


######################################
#### Appendix 1: Figure A6 - Example of projected CPUE in WEAs across PVC contours####
######################################

#with all PVC contours
footprints_app <- read_sf(here('Confidential', 'data', 'PVC', glue('logbook_all_{yrSrt}_{yrEnd}_50_pvcs_yrs_target.shp'))) %>%
  bind_rows(read_sf(here('Confidential', 'data', 'PVC', glue('logbook_all_{yrSrt}_{yrEnd}_75_pvcs_yrs_target.shp')))) %>%
  bind_rows(read_sf(here('Confidential', 'data', 'PVC', glue('logbook_all_{yrSrt}_{yrEnd}_95_pvcs_yrs_target.shp'))))

spp_preds <- all_spp %>% filter(fishery == 'dts_trawl')

basepred <- map_year(spp_preds, yr_vec = 2020, return_pred_df = T)
newpred <- map_year(spp_preds, yr_vec = 2050:2100, return_pred_df = T) %>%
  mutate(est_comp = est) %>% dplyr::select(-est)
both <- basepred %>% 
  left_join(newpred, by=c('X','Y')) %>% 
  mutate(diffpred = est_comp-est) %>%
  mutate(perc = (est_comp-est)/est*100)

bbox = st_bbox(projection_extent)

dts_cpue_map_app <- ggplot(both) +
  geom_raster(aes(x = X, y = Y, fill = perc), interpolate = F) + 
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red3", limits = c(-100,100)) +
  geom_sf(data = owe_shps, fill = NA, color = 'gray50', linewidth = 0.7) +
  geom_sf(data = footprints_app %>% filter(port_name == ifelse(oweas == "ORCAS", 'Eureka', 'Coos Bay') & subfishery == 'DTS'), aes(col = cont), fill = NA, linewidth = 0.7) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf(data = port_locs %>% filter(IOPAC == ifelse(oweas == "ORCAS", 'Eureka', 'Coos Bay')), size = 3, shape = 22, fill = "red3", color = "gray40", show.legend = FALSE) +
  geom_sf_text(data = coast_outline, aes(label = NAME),
               nudge_x = c(-300, -220, -230, -200),
               nudge_y = c(450,  300,  -200,   -150),
               size = 3) +
  xlim(bbox[1],550) +
  ylim(4410,5100) +
  scale_color_manual(values = c('red', 'black', '#f6b61c'), name = 'Fishing\nfootprint\nPVC%') +
  labs(title = ifelse(oweas == "ORCAS", "a) Eureka DTS", "a) Coos Bay DTS"), x = "", y = "", fill = "%\u0394\nbiomass\nindex") +
  theme_bw() +
  theme(plot.title = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.position = "none",
        legend.title.align = 0.5, 
        legend.background = element_rect(fill = "white")) +
  coord_sf(expand = FALSE)

spp_preds <- all_spp %>% filter(fishery == 'ndts_trawl')

basepred <- map_year(spp_preds, yr_vec = 2020, return_pred_df = T)
newpred <- map_year(spp_preds, yr_vec = 2050:2100, return_pred_df = T) %>%
  mutate(est_comp = est) %>% dplyr::select(-est)
both <- basepred %>% 
  left_join(newpred, by=c('X','Y')) %>% 
  mutate(diffpred = est_comp-est) %>%
  mutate(perc = (est_comp-est)/est*100)

ndts_cpue_map_app <- ggplot(both) +
  geom_raster(aes(x = X, y = Y, fill = perc), interpolate = F) + 
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red3", limits = c(-100,100)) +
  geom_sf(data = owe_shps, fill = NA, color = 'gray50', linewidth = 0.7) +
  geom_sf(data = footprints_app %>% filter(port_name == ifelse(oweas == "ORCAS", 'Coos Bay', 'Eureka') & subfishery == 'nDTS'), aes(col = cont), fill = NA, linewidth = 0.7) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf(data = port_locs %>% filter(IOPAC == ifelse(oweas == "ORCAS", 'Coos Bay', 'Eureka')), size = 3, shape = 22, fill = "red3", color = "gray40", show.legend = FALSE) +
  geom_sf_text(data = coast_outline, aes(label = NAME),
               nudge_x = c(-300, -220, -230, -200),
               nudge_y = c(450,  300,  -200,   -150),
               size = 3) +
  xlim(bbox[1],550) +
  ylim(4410,5100) +
  scale_color_manual(values = c('red', 'black', '#f6b61c'), name = 'Fishing\nfootprint\nPVC%') +
  labs(title = ifelse(oweas == "ORCAS", "b) Coos Bay non-DTS", "b) Eureka non-DTS"), x = "", y = "", fill = "%\u0394\nbiomass\nindex") +
  theme_bw() +
  theme(plot.title = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.position = "none",
        #legend.position = c(0.75, 0.4),
        legend.title.align = 0.5, 
        legend.background = element_rect(fill = "white")) +
  coord_sf(expand = FALSE)

spp_preds <- all_spp %>% filter(fishery == 'mid_trawl')

basepred <- map_year(spp_preds, yr_vec = 2020, return_pred_df = T)
newpred <- map_year(spp_preds, yr_vec = 2050:2100, return_pred_df = T) %>%
  mutate(est_comp = est) %>% dplyr::select(-est)
both <- basepred %>% 
  left_join(newpred, by=c('X','Y')) %>% 
  mutate(diffpred = est_comp-est) %>%
  mutate(perc = (est_comp-est)/est*100)

mid_cpue_map_app <- ggplot(both) +
  geom_raster(aes(x = X, y = Y, fill = perc), interpolate = F) + 
  scale_fill_gradient2(low = "dodgerblue2", mid = "white", high = "red3", limits = c(-100,100)) +
  geom_sf(data = owe_shps, fill = NA, color = 'gray50', linewidth = 0.7) +
  geom_sf(data = footprints_app %>% filter(port_name == ifelse(oweas == "ORCAS", 'Brookings', 'Crescent City') & subfishery == 'MDT'), aes(col = cont), fill = NA, linewidth = 0.7) +
  geom_sf(data = coast_outline, fill = 'grey90') +
  geom_sf(data = port_locs %>% filter(IOPAC == ifelse(oweas == "ORCAS", 'Brookings', 'Crescent City')), size = 3, shape = 22, fill = "red3", color = "gray40", show.legend = FALSE) +
  geom_sf_text(data = coast_outline, aes(label = NAME),
               nudge_x = c(-300, -220, -230, -200),
               nudge_y = c(450,  300,  -200,   -150),
               size = 3) +
  xlim(bbox[1],550) +
  ylim(4410,5100) +
  scale_color_manual(values = c('red', 'black', '#f6b61c'), name = 'Fishing\nfootprint\nPVC%') +
  labs(title = ifelse(oweas == "ORCAS", "c) Brookings midwater", "c) Crescent City midwater"), x = "", y = "", fill = "%\u0394\nbiomass\nindex") +
  theme_bw() +
  theme(plot.title = element_text(size = 8),
        plot.margin = margin(0, 0, 0, 0, "pt"),
        panel.grid = element_blank(),
        axis.text = element_text(size = 8),
        axis.text.x = element_text(angle = 90, vjust = 0.4),
        legend.title.align = 0.5, 
        legend.background = element_rect(fill = "white")) +
  coord_sf(expand = FALSE)

dts_cpue_map_app + ndts_cpue_map_app + mid_cpue_map_app
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A6_all_three_fisheries_cpue.tiff"), width = 6.5, height = 4.5)



######################################
#### Appendix 1: Figure A7 - Figure 7 from main text across all PVC contours####
######################################

cpue_inside_outside_fps_app <- cpue_gridcell_fishery_year %>% 
  left_join(sdm_id) %>% 
  st_as_sf(coords=c("longitude","latitude"), crs=st_crs(footprints_yrs_app)) %>% 
  st_join(footprints_yrs_app %>%
            mutate(fish_group = ifelse(fish_group == "DTS", "dts_trawl", ifelse(fish_group == "nDTS", "ndts_trawl", "mid_trawl")),
                   year = as.numeric(year)) %>% 
            filter(port_name %in% ports.all)) %>% 
  mutate(inside_fp = ifelse(!is.na(port_name),"Yes","No")) %>% # indicates whether or not a grid cell is within a footprint
  filter(fishery == fish_group) #need to make sure only footprints for the specific fishery species is kept

# for this analysis, we only care about values inside footprints (and inside vs. outside WEAs) and fleets that overlapped > 4 years
cpue_inside_fps_app <- cpue_inside_outside_fps_app %>% 
  filter(inside_fp == "Yes") 

# now, we join these to the unioned WEAs (we don't care about the identity/name of the WEA here)
owe_shps_union <- owe_shps_union %>% st_as_sf() %>% mutate(is_wea = "Yes")

cpue_fps_inside_outside_weas_app <- cpue_inside_fps_app %>% 
  st_join(owe_shps_union %>% st_as_sf) %>% 
  mutate(is_wea = ifelse(is.na(is_wea), "No", is_wea)) # just like the footprints above, we indicate whether a point is inside or outside WEAs

#how many SDM points are included in the estimate inside and outside WEAs?
cpue_fps_inside_outside_weas_counts_app <- cpue_fps_inside_outside_weas_app %>%
  group_by(fishery, cont, port_name, is_wea) %>%
  summarise(SDM_points = length(unique(sdm_grid_id)))

sdm_points_plot_app <- ggplot(cpue_fps_inside_outside_weas_counts_app) +
  geom_bar(aes(SDM_points, color = cont)) +
  facet_grid(fishery ~ port_name) +
  ylim(0, 25)

cpue_fps_inside_outside_weas_edit_app <- cpue_fps_inside_outside_weas_app %>%
  left_join(cpue_fps_inside_outside_weas_counts_app %>% 
              st_drop_geometry() %>%
              filter(is_wea == "Yes") %>%
              dplyr::select(cont, fishery, port_name, SDM_points), by = c("cont", "fishery", "port_name")) %>%
  filter(SDM_points > 1) %>%
  left_join(over_area_app %>% dplyr::select(port_name, fish_group, cont, years_of_overlap_footprints) %>% mutate(fish_group = ifelse(fish_group == "DTS trawl", "dts_trawl", ifelse(fish_group == "non-DTS trawl", "ndts_trawl", "mid_trawl"))) %>% distinct(), by = c("port_name", "fish_group", "cont")) %>%
  filter(!is.na(years_of_overlap_footprints))


cpue_time_series_app <- cpue_fps_inside_outside_weas_edit_app %>% 
  st_coordinates() %>%
  as.data.frame() %>%
  bind_cols(cpue_fps_inside_outside_weas_edit_app) %>%
  st_drop_geometry() %>% 
  dplyr::select(-c(year.y, area)) %>% # get single CPUE values for each set of coordinates by year/fishery/port/inside vs. outside wea
  distinct() %>% 
  st_drop_geometry() %>% 
  group_by(year.x, fishery, cont, port_name, is_wea) %>% # take a mean CPUE by year/fishery/port/inside vs. outside wea
  summarise(mean_cpue = mean(est, na.rm = T),
            median_cpue = median(est, na.rm = T),
            sd_cpue = sd(est, na.rm = T),
            se_cpue = sd_cpue/sqrt(n())) %>%
  group_by(fishery, cont, port_name, is_wea) %>% 
  mutate(mean_5year = slide_dbl(mean_cpue, mean, .before = 2, .after = 2),
         se_5year = slide_dbl(se_cpue, mean, .before = 2, .after = 2),
         sd_5year = slide_dbl(sd_cpue, mean, .before = 2, .after = 2)) %>%
  ungroup() %>%
  mutate(port_name = factor(port_name, levels = ports.all),
         fishery = factor(fishery, levels = c('dts_trawl', 'ndts_trawl', 'mid_trawl')),
         port_name = fct_drop(port_name))


cpue_ratio_fig_app <- cpue_time_series_app %>%
  pivot_wider(names_from = is_wea, values_from = c(median_cpue, mean_cpue, sd_cpue, se_cpue, mean_5year, sd_5year, se_5year)) %>%
  mutate(ratio = mean_5year_Yes/mean_5year_No,
         sd_ratio = ratio * sqrt((sd_5year_Yes/mean_5year_Yes)^2 + (sd_5year_No/mean_5year_No)^2),
         se_ratio = ratio * sqrt((se_5year_Yes/mean_5year_Yes)^2 + (se_5year_No/mean_5year_No)^2)) %>%
  filter(year.x >= 2020) %>%
  ggplot(aes(x = year.x, y = ratio, color = cont)) +
  geom_line() +
  geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"), se = F) +
  geom_ribbon(aes(ymin = ifelse(ratio-se_ratio < 0, 0, ratio-se_ratio), ymax = ratio+se_ratio), alpha = 0.2, linewidth = 0.1) +
  geom_hline(yintercept = 1, linetype=2) +
  facet_grid(fishery ~ port_name,
             labeller = labeller(fishery = subfishery.labs)) +
  scale_color_manual(values = c('red', 'black', '#f6b61c'), name = 'PVC (%)') +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) +
  xlab('Year') + 
  ylab('Ratio of biomass index')
ggsave(glue("Confidential/output/{oweas}_all_PVC_{yrSrt}_{yrEnd}_Figure_A7_time_series_ratio.tiff"), width = 6.5, height = 6)
