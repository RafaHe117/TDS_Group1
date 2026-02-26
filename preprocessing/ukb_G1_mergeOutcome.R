library(data.table)
rm(list = ls())

setwd("/rds/general/user/rh725/projects/hda_25-26/live/TDS/TDS_Group1")

ukb <- readRDS("ukb_G1_final.rds")
cvd <- readRDS("cvd_events.rds")

ukb$eid <- as.integer(ukb$eid)
cvd$eid <- as.integer(cvd$eid)
cvd$date <- as.Date(cvd$date)

first_cvd <- aggregate(date ~ eid, data = cvd, FUN = function(x) min(x, na.rm = TRUE))
names(first_cvd)[2] <- "cvd_first_date"

ukb$cvd_first_date <- first_cvd$cvd_first_date[match(ukb$eid, first_cvd$eid)]

ukb$date_recr <- as.Date(ukb$date_recr)

ukb$cvd_event     <- as.integer(!is.na(ukb$cvd_first_date))
ukb$cvd_prevalent <- as.integer(ukb$cvd_event == 1 & ukb$cvd_first_date <  ukb$date_recr)
ukb$cvd_incident  <- as.integer(ukb$cvd_event == 1 & ukb$cvd_first_date >= ukb$date_recr)

stopifnot(all(ukb$cvd_prevalent <= ukb$cvd_event))
stopifnot(all(ukb$cvd_incident  <= ukb$cvd_event))
stopifnot(all(ukb$cvd_prevalent + ukb$cvd_incident == ukb$cvd_event))

ukb <- ukb[, c(
  "eid",
  "cvd_event",
  "cvd_prevalent",
  "cvd_incident",
  "cvd_first_date",
  setdiff(names(ukb), c(
    "eid",
    "cvd_event",
    "cvd_prevalent",
    "cvd_incident",
    "cvd_first_date"
  ))
), drop = FALSE]

saveRDS(ukb, "ukb_G1_final.rds")