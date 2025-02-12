####This script examines confidentiality of the "North WA Coast" IOPAC port group to ensure confidentiality is maintained

north_wa_coast <- read.csv("Confidential/data/North WA Coast trawl landings.csv")

nwc_summ <- north_wa_coast %>%
  group_by(PARTICIPATION_GROUP_NAME) %>%
  summarise(n_years = n_distinct(PACFIN_YEAR),
            n_dealer_id = n_distinct(DEALER_ID),
            n_dealer_num = n_distinct(DEALER_NUM),
            n_vessels = n_distinct(VESSEL_ID),
            n_ftid = n_distinct(FTID))

