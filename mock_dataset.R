##################
## Below are codes to generate mock dataset for the synergistic effect analysis
## Note each variable within the dataset was simulated randomly thus no  
## meaningful association is expected in the dataset
##################

## To note, I shared a few other processed datasets without privacy concern
## please download them into the data folder before running this code
## https://drive.google.com/file/d/1UCctnt9r4C7FM1QuVtHQYpYtQ2Dd8HGD/view?usp=sharing

## "eh85_wf15_binary_1772zcta_0619.rds" for exposure time-series data
## "zcta_list_eh85_wf15_binary_1772zcta.csv" for a list of California zctas with exposure and health data
## "census_2010_pop.csv" zcta centroids and popualtion (2010 census)
## "zip_selected_hpi3_hpi_wt_041823_clean.csv" zcta-specific Healthy Places Index data after cleaning

## setting the stage
library(data.table)
outdir1 <- "" ## working directory for the project

### code to create necessary folders
if (!dir.exists(file.path(outdir1, "data"))) dir.create(file.path(outdir1, "data"))
if (!dir.exists(file.path(outdir1, "results"))) dir.create(file.path(outdir1, "results"))
if (!dir.exists(file.path(outdir1, "figures"))) dir.create(file.path(outdir1, "figures"))

### code to create a mock health dataset
#################
## read in a list of zctas with health data
zctas <- fread(file.path(outdir1, "data", "zcta_list_eh85_wf15_binary_1772zcta.csv"))

set.seed(1234)
ha <- data.frame(
  # zcta = rep(c(90006, 90007, 90008, 90010, 90011, 90012, 90013, 90014, 90015, 90016), each=5113),
  zcta = rep(zctas$zcta, each=5113),
  date = rep(seq(from = as.Date("2006-01-01"), to = as.Date("2019-12-31"), by = 1), time = nrow(zctas)),
  circulatory = rpois(nrow(zctas)*5113, 1), 
  respiratory = rpois(nrow(zctas)*5113, 1.2))
ha <- ha[ha$respiratory!=0 | ha$circulatory != 0, ]
write.csv(ha, file.path(outdir1, "data", "combo_99_19_patzip_cvd_resp_EM.csv"), row.names = FALSE)
#################
