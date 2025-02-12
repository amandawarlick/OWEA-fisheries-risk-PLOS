
## *A framework to evaluate dynamic social and ecological interactions between offshore wind energy development and commercial fisheries: a US West Coast perspective* 

#### Amanda J. Warlick, Owen Liu, Janelle Layton, Chris J. Harvey, Jameal Samhouri, Elliott Hazen, Kelly Andrews

##### Please contact Kelly Andrews (kelly.andrews@noaa.gov) or Amanda Warlick (amanda.warlick@noaa.gov) for questions about the code, analysis, or underlying data.

DOI: https:

________________________________________________________________________________

## Abstract

Offshore wind energy (OWE) planning is occurring rapidly alongside efforts to understand the potential effects of long-term environmental variability and climate change on social-ecological systems. To minimize potential conflicts between current and new ocean-use sectors, there is a need to identify tradeoffs between OWE development and other ocean users under dynamic environmental conditions. Here, we present a framework for evaluating the risk of groundfish fisheries being displaced from traditional fishing grounds by the designation of proposed OWE areas (OWEAs) and how risk may be affected by climate change impacts on targeted species. Specifically, we use fishery-dependent catch data from three groundfish fisheries to derive annual fishing “footprints” for port groups along the U.S. West Coast (1994-2020). We calculate the historical risk of these fleets being displaced from fishing grounds that have been proposed as sites for OWE development using an exposure-vulnerability framework. Risk varies across fishing fleets, but generally corresponds to a fleet’s target species and distance to proposed OWEAs. We then use existing climate-driven projections to map the spatial distribution of targeted species biomass for each of the three fisheries from 2020 to 2100. In some cases, future target species biomass indices have higher predicted values inside proposed OWEAs compared to outside OWEAs, indicating that incorporating climate change impacts may increase the perceived risk of displacement for these fleets. These results indicate that tradeoffs between commercial fishing and OWE development will not be fully understood unless the effects of climate change are incorporated into marine spatial planning efforts and efforts to develop appropriately-scaled mitigation measures. 

### Table of Contents 

#### [Scripts](./scripts)

Contains the script to process raw data and run all analyses. 

1. create_footprint_target_sp.R  
This file pulls in cleaned logbook data ('data/LB.ShortForm.with.Hake.Strat 26 Apr 2022.RData') and generates and saves (1) summarized data for within/outside WEA-overlap area statistics ('data/logbook/sp_logbook_hauls_ports_IOPAC_depths.rds') and (2) annual fishing footprints for each landings port and subfishery ('data/PVC/logbook_sp_pvcs_yrs.shp'). 

Users specify which sp (subfishery) is of interest and which PVC (cont) in the kernel density calculation. The script must be run separately for each subfishery and/or choice of PVC. The script relies on functions that can be found in logbook_fishshed_functions.R  

- This script relies on 'LB.ShortForm.with.Hake.Strat.RData' file and the IOPAC look-up table 'iopac_logbook_conversion_table.csv', both in the 'data' folder. 

# 2. fishticket_rev.Rmd: no longer using
# This file contains code adapted from Owen Liu to clean fish ticket data and assign a target species according to fish ticket.
# - This script relies on two large .csv files containing raw fish ticket data.

3. Liu et al. groundfish species distribution projections are combined with footprint summaries to generate manuscript summary statistics and figures (script available upon request). 
 
#### [Data](./Data) 

Contains raw and processed data necessary for developing summary statistics and figures. 


### Details of Article 

Warlick AJ, O Liu, J Layton, CJ Harvey, J Samhouri, EL Hazen, K Andrews. 2025. A framework to evaluate dynamic social and ecological interactions between offshore wind energy development and commercial fisheries: a US West Coast perspective. PLOS Climate. 


### Disclaimer

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project content is provided on an "as is" basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.

Created by Kelly Andrews (February 2024) via transfer from amandawarlick/WEA_mapping

