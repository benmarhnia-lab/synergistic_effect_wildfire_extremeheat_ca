################################
### These are codes used to estimate synergistic effects of wildfire smoke
### and extreme heat on cardiorespiratory hospitalizations in California
### at state-level (case-crossover design) and ZIP code-level (within community 
### matched design + spatial Bayesian hierarchical model). 
### Effect modification by community characteristics were evaluated as well.
### Codes written and organized by Chen Chen on 3/13/2023
################################

outdir1 <- "" ## working directory for the project

## run case-crossover analysis at state-level
library(survival)
library(msm)
##################
if (!file.exists(file.path(outdir1, "results", "state_models"))) dir.create(file.path(outdir1, "results", "state_models"))
dataset <- "eh85_wf15_binary_1772zcta"

## calculation of reri--based on SAS codes from VanderWeele and Knol 2014 Appendix
## variance for joint effect based on deltamethod (car package)
## Note "se" in print-outs are se of transformed coefficients (delta method for joint effect, VanderWeele & Knol 2014 method for RERI). 
## The rest are se for coefficients of the model.
reri_interaction <- function(b, v, exposure1, exposure2) {
  reri <- exp(b[1] + b[2] + b[3]) - exp(b[1]) - exp(b[2]) + 1
  k1 <- exp(b[1] + b[2] + b[3]) - exp(b[1])
  k2 <- exp(b[1] + b[2] + b[3]) - exp(b[2])
  k3 <- exp(b[1] + b[2] + b[3])
  vreri <- v[1, 1] * k1 * k1 + 
    v[2, 2] * k2 * k2 + 
    v[3, 3] * k3 * k3 + 
    2 * v[1, 2] * k1 * k2 +
    2 * v[1, 3] * k1 * k3 +
    2 * v[2, 3] * k2 * k3
  se_reri <- sqrt(vreri)
  reri_ci95_l <- reri - 1.96 * se_reri
  reri_ci95_u <- reri + 1.96 * se_reri
  est <- exp(b)
  se <- sqrt(c(v[1, 1], v[2, 2], v[3, 3]))
  est_l <- exp(b - 1.96 * se)
  est_u <- exp(b + 1.96 * se)
  k3_se <- deltamethod(~exp(x1+x2+x3), b, v)
  k3_l <- k3 - 1.96 * k3_se
  k3_u <- k3 + 1.96 * k3_se
  temp <- data.frame( 
    est=c(est, reri, k3),
    ll=c(est_l, reri_ci95_l, k3_l),
    ul=c(est_u, reri_ci95_u, k3_u),
    se=c(se, se_reri, k3_se))
  rownames(temp) <- c(exposure1, exposure2, "interaction_multiplicative", "interaction_reri", "joint_effet")
  return(temp)
}

## read in exposure dataset
exposure <- readRDS(file.path(outdir1, "data", paste0(dataset, "_0619.rds")))

## clean health data to create dataset for case-crossover
ha <- fread(file.path(outdir1, "data", "combo_99_19_patzip_cvd_resp_EM.csv"))
names(ha)[1:2] <- c("zcta", "case_date")
## Focus on overall population
ha <- ha[ha$case_date > as.Date("2005-12-31") & ha$case_date < as.Date("2020-01-01"), .(zcta, case_date, circulatory, respiratory)]
length(unique(ha$zcta)) ## 1778 zipcodes
unique(ha$zcta[!(ha$zcta %in% exposure$zcta)]) ## no exposure for 93042, no health data for 91719
ha <- ha[ha$zcta %in% exposure$zcta, ] ## keep of data from zipcodes with exposure data
length(unique(ha$zcta)) ## 1772 zipcodes
ha$csdresp <- ha$circulatory + ha$respiratory

## create variables for merging
ha$wday <- wday(ha$case_date)
ha$month <- month(ha$case_date)
ha$year <- year(ha$case_date)

## create empty variable 
out <- numeric()

lags <- c("same_day", "lag1")
outcomes <- c("csdresp", "circulatory", "respiratory")
for (lag in lags) {
  if (lag!="same_day") { ## create lags in exposure
    nn <- as.numeric(gsub("lag", "", lag))
    new <- copy(exposure)
    new <- new[order(new$date), ]
    new[, eh:=shift(eh, n = nn, fill = NA, type = "lag"), by=.(zcta)]
    new[, wf:=shift(wf, n = nn, fill = NA, type = "lag"), by=.(zcta)]
  } else {
    new <- copy(exposure)
    new <- new[order(new$date), ]
  }
  new$wday <- wday(new$date)
  for (outcome in outcomes) {
    if (outcome=="csdresp") {
      ha$id <- 1:nrow(ha)
      dt <- merge(ha, new, by = c("zcta", "wday", "month", "year"), all.x = TRUE, allow.cartesian=TRUE)
    } else if (outcome=="circulatory") {
      csd <- ha[circulatory>0 & !is.na(circulatory), ]
      csd$respiratory <- NULL
      csd$id <- 1:nrow(csd)
      dt <- merge(csd, new, by = c("zcta", "wday", "month", "year"), all.x = TRUE, allow.cartesian=TRUE)
    } else if(outcome=="respiratory") {
      resp <- ha[respiratory>0 & !is.na(respiratory), ]
      resp$circulatory <- NULL
      resp$id <- 1:nrow(resp)
      dt <- merge(resp, new, by = c("zcta", "wday", "month", "year"), all.x = TRUE, allow.cartesian=TRUE)
    }
    dt$case <- ifelse(dt$case_date==dt$date, 1, 0)
    loc <- grep(outcome, names(dt))
    names(dt)[loc] <- "outcome"
    
    m <- clogit(case ~ eh*wf + strata(id), data=dt, weights=outcome, method="approximate")
    saveRDS(m, file.path(outdir1, "results", "state_models", paste0(dataset, "_", lag, "_", outcome, ".rds")))

    dt <- dt[!is.na(dt$wf & !is.na(dt$eh)), ] ## remove missing values
    m.eh <- clogit(case ~ eh + strata(id), data=dt, weights=outcome, method="approximate")
    saveRDS(m.eh, file.path(outdir1, "results", "state_models", paste0(dataset, "_", lag, "_", outcome, "_eh.rds")))

    m.wf <- clogit(case ~ wf + strata(id), data=dt, weights=outcome, method="approximate")
    saveRDS(m.wf, file.path(outdir1, "results", "state_models", paste0(dataset, "_", lag, "_", outcome, "_wf.rds")))

    ## create results for plotting
    temp <- data.frame(rbind(summary(m.eh)$conf.int[, c(1, 3, 4)], 
                             summary(m.wf)$conf.int[, c(1, 3, 4)]))
    row.names(temp) <- c("extreme_heat", "wildfire")
    names(temp) <- c("est", "ll", "ul")
    temp$se <- c(sqrt(m.eh$var), sqrt(m.wf$var))
    temp <- rbind(reri_interaction(b=m$coefficients, v=m$var, "extreme_heat only", "wildfire only"),
                  temp)
    temp$estimate <- row.names(temp)
    out <- rbind(out, cbind(outcome = outcome, lag = lag, temp))
    
    m <- m.wf <- m.eh <- dt <- temp <- NULL
  }
}

write.csv(out, file.path(outdir1, "results", paste0(dataset, "_state_model_summary.csv")), row.names = FALSE)

##################

## run within-community matched design (Schwarz et al. 2021): controls identified as in the same summer and not an event day, use inverse distance weighting of all control days
## matched time-series (Liu et al. 2017): controls identified as 1) within the window of 7 calendar days before or after the exposed day in another year and 2) separated from any other exposed day for more than 2 days
## Bobb et al. 2014: controls identified as 1) within the window of 3 days before and after the exposed day in another year and 2) separated from any other exposed day for more than 2 days
# library(lme4)
library(survival)
dataset <- "eh85_wf15_binary_1772zcta"
##################
## identify potential control days for each exposed day based on Liu et al. method with controls identified as: 
## 1) within the window of nbuffer calendar days before or after the exposed day in another year and 
## 2) separated from any other exposed day for more than 2 days
year.control <- function(exposed, control, years, nbuffer) { 
  if (length(exposed) > 0) {
    out <- lapply(exposed, function(e_day) {
      e_buffer <- yday(e_day) + (-nbuffer:nbuffer) ## potential days of year for control (7 day window from exposure)
      buffer_year <- years[years!=year(e_day)] ## potential years for control (not the year of event)
      c_day <- sapply(buffer_year, function(yr) { ## potential control days based on exposed date
        
        ndays <- yday(as.Date(paste0(yr, "-12-31"))) ## days in the year of exposure
        ndays_bf <- yday(as.Date(paste0(yr - 1, "-12-31"))) ## days in the year before exposure
        
        sapply(e_buffer, function(dd) {
          if(dd < 0) {
            as.IDate(paste(yr-1, dd + ndays_bf), format = "%Y %j")
          } else if (dd > ndays) {
            as.IDate(paste(yr+1, dd - ndays), format = "%Y %j")
          } else {
            as.IDate(paste(yr, dd), format = "%Y %j")
          }
        })
      })
      c_pool <- control[control %in% c(c_day)]  ## potential control days after excluding days close to other exposure
      return(c(e_day, c_pool))
    })
    names(out) <- exposed
  } else {
    out <- numeric()
  }
  return(out)
}

## identify potential control days for each exposed day baed on schwartz et al. and Liu et al. method: controls identified as 
## 1) within the window of buffer calendar days (nbuffer) before or after the exposed day and 
## 2) separated from any other exposed day for more than 2 days
month.control <- function(exposed, control, nbuffer) { 
  if (length(exposed) > 0) {
    out <- lapply(exposed, function(e_day) {
      e_buffer <- e_day + (-nbuffer:nbuffer) ## potential days for control (buffer day window from exposure)
      c_pool <- control[control %in% e_buffer]  ## potential control days after excluding days close to other exposure
      return(c(e_day, c_pool))
    })
    names(out) <- exposed
  } else {
    out <- numeric()
  }
  return(out)
}

## randomly selected n.bobb.control controls for each exposed day, 
## use glm for each zcta
glmer.analysis.zcta <- function(outcome.dt, event) { 
  if (nrow(outcome.dt) > 0) {
    f <- reformulate("exposed", response = event)
    m <- tryCatch({
      glm(f , data = outcome.dt, family = "poisson")
    }, condition = function(cond) {
      # cat("\t", "glmer", zcta, exposure, event, as.character(cond))
      cond$call <- NULL
      cond
    })
    # }
  } else {
    m <- simpleError(paste("No available exposed day"))
    m$call <- NULL
  }
  
  if (!inherits(m, what = "condition")) {
    est <- summary(m)$coefficients["exposed",1]
    se <- summary(m)$coefficients["exposed",2]
    temp <- c(exp(est), exp(est - 1.96*se),  exp(est + 1.96*se), length(unique(outcome.dt$id)))
  } else {
    temp <- m
  }
  return(temp)
}

## weight the outcome of all potential controls using the inverse distance in year 
## and directly calculate the incidence ratio for each match exposed and control group, then calculate the average
year.wt.analysis <- function(exposure, event, c.list, outcome.dt) { 
  if (length(c.list) > 0) {
    out <- numeric()
    for (i in 1:length(c.list)) {
      if (length(c.list[[i]]) > 1) {
        e_day <- c.list[[i]][1]
        c_pool <- data.table(date=c.list[[i]], exposed = c(1, rep(0, length(c.list[[i]])-1)))
        c_pool$wt <- 1/round(abs(as.numeric(difftime(c_pool$date, e_day, units = "days"))/365), digits = 0)
        c_pool$wt[1] <- 0 ## assign 0 weight to the day with exposure
        c_pool <- merge(c_pool, outcome.dt, by="date", all.x=TRUE)
        
        rr <- c_pool[exposed==1, eval(as.name(event))]/c_pool[, sum(eval(as.name(event))*wt)/sum(wt)]
        rr <- ifelse(is.infinite(rr), NA, rr)
        out <- c(out, rr)
      }
    }
    
    if (length(out)==0) {
      temp <- simpleError(paste("No available control day"))
    } else {
      temp <- c(mean(out, na.rm=TRUE), sum(!is.na(out)))
      if (is.na(temp[1])) {
        temp <- simpleError(paste("No non-infinite RR"))
        temp$call <- NULL
      }
    }
    
  } else {
    temp <- simpleError(paste("No available exposed day"))
    temp$call <- NULL
  }
  return(temp)
}

## weight the outcome of all potential controls using the distance in day and 
## directly calculate the incidence ratio for each match exposed and control group, then calculate the average
month.wt.analysis <- function(exposure, event, c.list, outcome.dt) { 
  if (length(c.list) > 0) {
    out <- numeric()
    for (i in 1:length(c.list)) {
      if (length(c.list[[i]]) > 1) {
        e_day <- c.list[[i]][1]
        c_pool <- data.table(date=c.list[[i]], exposed = c(1, rep(0, length(c.list[[i]])-1)))
        c_pool$wt <- 1/round(abs(as.numeric(difftime(c_pool$date, e_day, units = "days"))), digits = 0)
        c_pool$wt[1] <- 0 ## assign 0 weight to the day with exposure
        c_pool <- merge(c_pool, outcome.dt, by="date", all.x=TRUE)
        
        rr <- c_pool[exposed==1, eval(as.name(event))]/c_pool[, sum(eval(as.name(event))*wt)/sum(wt)]
        rr <- ifelse(is.infinite(rr), NA, rr)
        out <- c(out, rr)
      }
    }
    
    if (length(out)==0) {
      temp <- simpleError(paste("No available control day"))
    } else {
      temp <- c(mean(out, na.rm=TRUE), sum(!is.na(out)))
      if (is.na(temp[1])) {
        temp <- simpleError(paste("No non-infinite RR"))
        temp$call <- NULL
      }
    }
    
    
  } else {
    temp <- simpleError(paste("No available exposed day"))
    temp$call <- NULL
  }
  return(temp)
}

dt <- readRDS(file.path(outdir1, "data", paste0(dataset, "_0619.rds")))
dt$yday <- yday(dt$date)

## Focus on overall population
ha <- fread(file.path(outdir1, "data", "combo_99_19_patzip_cvd_resp_EM.csv"))
names(ha)[1:2] <- c("zcta", "date")
ha <- ha[ha$date > as.Date("2005-12-31") & ha$date < as.Date("2020-01-01"), .(zcta, date, circulatory, respiratory)]

bobb.control.n <- tarik.control.n <- 4
events <- c("csdresp")
nms <- c(paste(rep(c("year_glmer", "month_glmer"), each=4), rep(c("rr", "rr_ll", "rr_ul", "ngrp"), times=2), sep="_"), 
         paste(rep(c("year_wt", "month_wt"), each=2), rep(c("rr", "ngrp"), times=2), sep="_"))
nms <- paste(rep(nms, times=3), rep(c("eh1wf1", "eh1wf0", "eh0wf1"), each=length(nms)), sep="_")
bar <- data.frame(zcta=rep(unique(dt$zcta), each = length(events)), event = events, 
                  eh1wf1 = NA, eh1wf0 = NA, eh0wf1 = NA, 
                  bobb_control_pool = NA, 
                  setNames(replicate(length(nms), NA, simplify = F), nms)
)
fail <- numeric()
set.seed(824)


for (i in 1:nrow(bar)) {
  event_ <- bar$event[i]
  ha_ <- ha[ha$zcta==bar$zcta[i], ]
  foo <- dt[dt$zcta==bar$zcta[i], ]
  foo <- merge(foo, ha_[, .(date, circulatory, respiratory)], by="date", all.x = TRUE)
  setnafill(foo, fill=0, cols = c("circulatory", "respiratory")) ## fill in zeros for days without either csd or resp
  foo$csdresp <- foo$circulatory + foo$respiratory
  
  ## potential control identification method from Bobb et al--remove those close to exposure
  buffer_days <- foo$date[foo$eh==1 | foo$wf==1]
  buffer_days <- unique(c(buffer_days-3, buffer_days-2, buffer_days-1, buffer_days, buffer_days + 1, buffer_days + 2, buffer_days + 3))
  day00_bobb <- foo$date[!(foo$date %in% buffer_days)]
  bar$bobb_control_pool[i] <- length(day00_bobb) ## total number of potential control days excluding those close to exposure
  
  for (exposure_ in c("eh1wf1", "eh1wf0", "eh0wf1")) {
    ## identify exposed days using Bobb method
    days <- foo$date[foo$eh==as.numeric(substring(exposure_, 3, 3)) & foo$wf==as.numeric(substring(exposure_, 6, 6))]
    bar[i, exposure_] <- length(days)
    
    ## identify potential control days for each exposed day based on Bobb et al. method
    baz <- year.control(days, day00_bobb, 2006:2019, 7)
    
    ## weighted analysis method adopted from Schwarz et al.
    temp2 <- year.wt.analysis(exposure_, event_, baz, foo)
    if (!inherits(temp2, what = "condition")) {
      bar[i, paste0(c("year_wt_rr_", "year_wt_ngrp_"), exposure_)] <- temp2
    } else {
      fail <- rbind(fail, data.frame(zcta=bar$zcta[i], event = event_, exposure = exposure_, ngrps = length(baz), method = "year_wt", error=as.character(temp2)))
      bar[i, paste0("year_wt_ngrp_", exposure_)] <- length(baz)
    }
    
    ## control identification method modified from Schwarz et al.: same year within 30 days before and after the exposure and exclude these within three days of an event (second is the same as bobb et al. method)
    baz2 <- month.control(days, day00_bobb, 30)
    
    ## weighted analysis method adopted from Schwarz et al.
    temp4 <- month.wt.analysis(exposure_, event_, baz2, foo)
    if (!inherits(temp4, what = "condition")) {
      bar[i, paste0(c("month_wt_rr_", "month_wt_ngrp_"), exposure_)] <- temp4
    } else {
      fail <- rbind(fail, data.frame(zcta=bar$zcta[i], event = event_, exposure = exposure_, ngrps = length(baz2), method = "month_wt", error=as.character(temp4)))
      bar[i, paste0("month_wt_ngrp_", exposure_)] <- length(baz2)
    }
  }
  
  if (i%%500==0) cat("\n", "first three methods finished zcta #", i, "\n")
}

pop <- fread(file.path(outdir1, "data", "census_2010_pop.csv"))[, .(pop, zcta)]

for (exposure_ in c("eh1wf1", "eh1wf0", "eh0wf1")) {
  ## use dataset created in state-wise sensitivity analysis--it involves sampling so better use the same dataset
  out.year <- readRDS(file.path(outdir1, "data", paste0(dataset, "_", exposure_, "_yearcontrol.rds")))
  out.month <- readRDS(file.path(outdir1, "data", paste0(dataset, "_", exposure_, "_monthcontrol.rds")))
  
  ## create necessary varibles
  out.year <- merge(out.year, pop, by = "zcta", all.x = TRUE)
  out.year$dow <- wday(out.year$date)
  out.month <- merge(out.month, pop, by = "zcta", all.x = TRUE)
  out.month$dow <- wday(out.month$date)
  
  for (i in 1:nrow(bar)) {
    event_ <- bar$event[i]
    foo.year <- out.year[out.year$zcta==bar$zcta[i], ]
    foo.month <- out.month[out.month$zcta==bar$zcta[i], ]
    
    ## analysis method from Bobb et al using yearly controls
    temp <- glmer.analysis.zcta(foo.year, event_)
    if (!inherits(temp, what = "condition")) {
      bar[i, paste0(c("year_glmer_rr_", "year_glmer_rr_ll_", "year_glmer_rr_ul_", "year_glmer_ngrp_"), exposure_)] <- temp
    } else {
      fail <- rbind(fail, data.frame(zcta=bar$zcta[i], event = event_, exposure = exposure_, ngrps = length(unique(foo.year$id)), method = "year_glmer", error=as.character(temp)))
      bar[i, paste0("year_glmer_ngrp_", exposure_)] <- length(unique(foo.year$id))
    }
    
    ## analysis method from Bobb et al using monthly controls
    temp3 <- glmer.analysis.zcta(foo.month, event_)
    if (!inherits(temp3, what = "condition")) {
      bar[i, paste0(c("month_glmer_rr_", "month_glmer_rr_ll_", "month_glmer_rr_ul_", "month_glmer_ngrp_"), exposure_)] <- temp3
    } else {
      fail <- rbind(fail, data.frame(zcta=bar$zcta[i], event = event_, exposure = exposure_, ngrps = length(unique(foo.month$id)), method = "month_glmer", error=as.character(temp3)))
      bar[i, paste0("month_glmer_ngrp_", exposure_)] <- length(unique(foo.month$id))
    }
    
    if (i%%500==0) cat("\n", exposure_, "glm methods finished zcta #", i, "\n")
  }
  
}

write.csv(bar, file.path(outdir1, "results", paste0(dataset, "_specific.csv")), row.names = FALSE)
write.csv(fail, file.path(outdir1, "results", paste0(dataset, "_specific_fail.csv")), row.names = FALSE)
#################

## Apply Bayesian spatial analysis to incorporate spatial heterogeneity into the estimates for reri
library(sp)
library(gstat) 
library("spBayes")
library(coda)
library(MBA)

library(rgdal)
library(sf)
library(ggplot2)
library(FRK)
library(fields)
library(RColorBrewer)

outdir2 <- file.path(outdir1, "figures", "spatial")
if (!dir.exists(outdir2)) dir.create(outdir2)
##################
dataset <- "eh85_wf15_binary_1772zcta"
if (!dir.exists(file.path(outdir2, dataset))) dir.create(file.path(outdir2, dataset))
sink(file.path(outdir2, dataset, "summary of running spatial bayesian.txt"))

methods <- c("month_wt", "year_wt", "month_glm", "year_glm")

## read in data for analysis
bar <- fread(file.path(outdir1, "results", paste0(dataset, "_specific.csv")))
fail <- fread(file.path(outdir1, "results", paste0(dataset, "_specific_fail.csv")))
coords <- read.csv(file.path(outdir1, "results", paste0("zcta_list_", dataset, ".csv")), as.is = TRUE) ## same as above but with population info
bar <- merge(bar, coords, by = "zcta", all.x=TRUE)

## change glmer to glm
names(bar) <- gsub("glmer", "glm", names(bar))
fail$method <- gsub("glmer", "glm", fail$method)

for (m in methods) {
  set.seed(824)
  cat("\n\n", m, "\n")
  fail_m <- fail[fail$method==m, ]
  fail.zcta <- unique(fail_m$zcta[fail_m$ngrps==0])
  baz <- bar[!(zcta %in% fail.zcta), ]
  cat("# of zctas removed due to no exposed day is", length(fail.zcta), "\n")
  nn <- nrow(baz)
  baz <- baz[baz$pop>1000 &!is.na(baz$pop), ]
  cat("# of zctas removed due to no pop or pop<=1000 is", nn-nrow(baz), "\n")
  cat("# of zctas left after 1k pop and exposure day criteria is", nrow(baz), "\n")
  ex.fail <- fail_m[fail_m$ngrps>0]
  ex.fail <- ex.fail[!ex.fail$zcta %in% fail.zcta, ]
  ex.fail <- ex.fail[ex.fail$zcta %in% baz$zcta, ]
  cat("# of zctas removed due to failures other than no exposed day or pop is", length(unique(ex.fail$zcta)), "\n")
  cat("zctas removed due to failures other than no exposed day or pop", "\n")
  print(merge(ex.fail, coords[,c("zcta", "pop")], by="zcta", all.x=TRUE))
  baz <- baz[!baz$zcta %in% ex.fail$zcta, ]
  
  ## calculate reri
  baz[, paste0(m, "_reri"):= eval(as.name(paste0(m, "_rr_eh1wf1"))) - eval(as.name(paste0(m, "_rr_eh1wf0"))) -
        eval(as.name(paste0(m, "_rr_eh0wf1"))) + 1]
  
  ## drop extreme values in glmer results
  if (sum(abs(baz[, paste0(m, "_reri"), with=FALSE]) > 50) >0) {
    fail.zcta2 <- unique(baz$zcta[abs(baz[, paste0(m, "_reri"), with=FALSE]) > 50])
    cat(length(fail.zcta2), "zcta removed due to large reri in", m, "\n")
    print(baz[baz$zcta %in% fail.zcta2, .(zcta, eval(as.name(paste0(m, "_reri"))), pop, eh1wf1, eh1wf0, eh0wf1)])
    baz <- baz[!(baz$zcta %in% fail.zcta2), ]
  }
  cat("# of zctas with estimates for all three exposures and", m, "is", nrow(baz), "\n")
  
  bayesDF <- baz[, c("zcta", paste0(m, "_reri"), "FinalLon", "FinalLat"), with=FALSE]
  names(bayesDF)[2] <- "est"
  
  # US census bureau datum
  spdf <- SpatialPointsDataFrame(coords = bayesDF[,.(FinalLon, FinalLat)],
                                 data = bayesDF)
  
  v1 <- variogram(est ~ 1, data = spdf)
  png(file.path(outdir2, dataset, paste0(m, nrow(bayesDF), "zcta_reri_variogram.png")))
  print(plot(v1, main=paste0(dataset, "_", m, nrow(bayesDF), "zcta"), cex=1.5))
  #Assumes Isotropy
  dev.off()
  
  ## conduct Bayesian model with flat priors
  # tau sq is nugget
  # sigma sq sill
  # phi range
  n.samples = 10000
  if (m=="year_wt"){
    bef.sp <- spLM(est ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,.(FinalLon, FinalLat)]),
                   starting = list("phi" = 2, "sigma.sq" = 0.8, "tau.sq" = .5),
                   tuning = list("phi" = 0.1, "sigma.sq" = 0.04, "tau.sq" = 0.025),
                   priors = list("phi.Unif" = c(0.001, 6), "sigma.sq.IG" = c(2, 1/0.8)
                                 , "tau.sq.IG" = c(2, 1/.5)
                   ), ## actually used--(2, 1/mean) priors for IG
                   cov.model = "spherical", n.samples = n.samples, verbose = TRUE, n.report=2000)
  } else if (m=="month_wt") {
    bef.sp <- spLM(est ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,.(FinalLon, FinalLat)]),
                   starting = list("phi" = 2, "sigma.sq" = 2, "tau.sq" = .75),
                   tuning = list("phi" = 0.1, "sigma.sq" = 0.1, "tau.sq" = 0.0375),
                   priors = list("phi.Unif" = c(0.001, 6), "sigma.sq.IG" = c(2, 1/2)
                                 , "tau.sq.IG" = c(2, 1/.75)
                   ), ## actually used--(2, 1/mean) priors for IG
                   cov.model = "spherical", n.samples = n.samples, verbose = TRUE, n.report=2000)
  } else if (m=="year_glm") {
    bef.sp <- spLM(est ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,.(FinalLon, FinalLat)]),
                   starting = list("phi" = 3, "sigma.sq" = 2, "tau.sq" = 0.7), ## glm results
                   tuning = list("phi" = 0.15, "sigma.sq" = 0.1, "tau.sq" = 0.035),
                   priors = list("phi.Unif" = c(0.001, 6), "sigma.sq.IG" = c(2, 1/2)
                                 , "tau.sq.IG" = c(2, 1/0.7)
                   ), ## actually used--(2, 1/mean) priors for IG
                   cov.model = "spherical", n.samples = n.samples, verbose = TRUE, n.report=2000)
  } else if (m=="month_glm") {
    bef.sp <- spLM(est ~ 1, data = bayesDF, coords = as.matrix(bayesDF[,.(FinalLon, FinalLat)]),
                   starting = list("phi" = 2, "sigma.sq" = 1, "tau.sq" = 1), ## glm results
                   tuning = list("phi" = 0.1, "sigma.sq" = 0.05, "tau.sq" = 0.05),
                   priors = list("phi.Unif" = c(0.001, 6), "sigma.sq.IG" = c(2, 1/1)
                                 , "tau.sq.IG" = c(2, 1/1)
                   ), ## actually used--(2, 1/mean) priors for IG
                   cov.model = "spherical", n.samples = n.samples, verbose = TRUE, n.report=2000)
  }
  
  
  cat("summary of thetas before burn-in", "\n")
  print(round(summary(mcmc(bef.sp$p.theta.samples))$quantiles, 3))
  png(file.path(outdir2, dataset, paste0(m, nrow(bayesDF), "zcta_reri_mcmctrace.png")))
  plot(bef.sp$p.theta.samples)
  dev.off()
  
  ## exclude 75% samples as burn-in
  burn.in <- floor(0.75*n.samples)
  bef.sp <- spRecover(bef.sp, start = burn.in, n.report=1000)
  
  cat("summary of thetas after burn-in", "\n")
  print(round(summary(mcmc(bef.sp$p.theta.recover.samples))$quantiles, 3))
  png(file.path(outdir2, dataset, paste0(m, nrow(bayesDF), "zcta_reri_mcmctrace_afterburnin.png")))
  plot(bef.sp$p.theta.recover.samples)
  dev.off()
  
  beta.samples <- bef.sp$p.beta.recover.samples
  w.samples <- bef.sp$p.w.recover.samples
  
  bayesDF$overall_mu <- mean(beta.samples)
  bayesDF$overall_sd <- sd(beta.samples)
  cat("Statewise mean is", mean(beta.samples), ", sd is", sd(beta.samples), "\n") ## state-level estimates after spatial pooling
  bayesDF$w_hat_mu <- apply(w.samples, 1, mean)
  bayesDF$w_hat_sd <- apply(w.samples, 1, sd)
  bayesDF$SNR <- bayesDF$w_hat_mu/bayesDF$w_hat_sd
  bayesDF$truncSNR <- ifelse(bayesDF$SNR < 2 & bayesDF$SNR > -2, NA, bayesDF$SNR )
  print(summary(bayesDF$est))
  print(summary(bayesDF$w_hat_mu))
  print(summary(bayesDF$w_hat_sd))
  print(summary(bayesDF$SNR))
  print(summary(bayesDF$truncSNR))
  write.csv(bayesDF, file.path(outdir1, "results", paste0(dataset, "_reri_", m, nrow(bayesDF), "zcta.csv")))
}
sink()
##################

## Meta regression
## I also did linear regression for comparison
## All indicators were coded so that higher value represents better socioeconomical status (SES) based on previous believes, except for the ethnicity groups.  
library(meta)
##################
dataset <- "eh85_wf15_binary_1772zcta"

methods <- c("month_wt", "year_wt", "month_glm", "year_glm")
nc <- c(991, 994, 981, 990) #zipcodes included for each method
names(nc) <- methods

## read in variables for meta regression
hpi <- fread(file.path(outdir1, "zip_selected_hpi3_hpi_wt_041823.csv"))

## population density
zcta.sp <- readOGR(file.path(outdir1, "tl_2010_06_zcta510",
                             "tl_2010_06_zcta510.shp"), stringsAsFactors = FALSE)
area <- zcta.sp@data[, c("ZCTA5CE10", "ALAND10")]
area$ZCTA5CE10 <- as.numeric(area$ZCTA5CE10)
hpi <- merge(hpi, area, by.x="ZIP", by.y="ZCTA5CE10", all.x=TRUE) 
hpi$ALAND10 <- as.numeric(hpi$ALAND10)
hpi$popden <- hpi$pop/hpi$ALAND10 * 10000 ## so that it's per 10k
## after transformation, all variables are true to its meaning and higher better
hpi_vbs <- c("employed", "abovepoverty", "bachelorsed", 
             "treecanopy", "AC", "insured", "percapitaincome", "automobile",
             "white", "black", "asian", "latino", "NativeAm","PacificIsl",
             "popden")
hpi <- hpi[!is.na(employed) & !is.na(treecanopy), c("ZIP", hpi_vbs), with=FALSE]
### using more EM variables
out <- numeric()
sink(file.path(outdir1, "summary of running meta regression reri_022723.txt"))
for (m in methods) {
  ## reri after pooling
  bayesDF <- fread(file.path(outdir1, "results", paste0(dataset, "_reri_", m, nc[m], "zcta.csv")))
  ntemp <- nrow(bayesDF)
  cat("\n\n", "# of zctas with", m, "based reri is", ntemp, "\n")
  bayesDF <- merge(bayesDF, hpi, by.x="zcta", by.y="ZIP", all.x=TRUE)
  bayesDF <- bayesDF[!is.na(employed), ]
  cat("# of zctas included in meta regression of", m, "based reri is", nrow(bayesDF), "\n")
  cat("# of zctas removed due to missing hpi in meta regression of", m, "based reri is", ntemp - nrow(bayesDF), "\n")
  
  ## reri before pooling--there is no CI here...
  bar[, paste0(m, "_reri"):= eval(as.name(paste0(m, "_rr_eh1wf1"))) - eval(as.name(paste0(m, "_rr_eh1wf0"))) -
        eval(as.name(paste0(m, "_rr_eh0wf1"))) + 1]
  baz <- bar[, c("zcta", paste0(m, "_reri")), with=FALSE]
  names(baz)[2] <- "est_raw"
  bayesDF <- merge(bayesDF, baz, by="zcta", all.x=TRUE) ## keep zctas included in Bayesian spatial pooling--removed extreme values
  
  nms <- paste(rep(c("lm", "meta", "lm_raw"), each=2), rep(c("coef", "se"), time=3), sep="_")
  estimate <- data.frame(vb=hpi_vbs,
                         min=NA, q1=NA, median=NA, q3=NA, max=NA, zcta5.used=NA, 
                         setNames(replicate(length(nms), NA, simplify = FALSE), nms))
  for (i in 1:nrow(estimate)) {
    if (estimate$vb[i]=="AC") {
      cache <- copy(bayesDF)
      bayesDF <- bayesDF[!is.na(bayesDF$AC), ]
      cat("# of zctas removed due to missing AC in meta regression of", m, "based reri is", nrow(cache) - nrow(bayesDF), "\n")
    }
    estimate[i, "zcta5.used"] <- nrow(bayesDF)
    
    ## linear regression of associations after bayesian spatial pooling
    g <- lm(as.formula(paste0("w_hat_mu ~ ", estimate$vb[i])), data=bayesDF)
    estimate[i, 2:6] <- quantile(bayesDF[, eval(as.name(estimate$vb[i]))], probs = c(0, 0.25, 0.5, 0.75, 1))
    estimate[i, nms[1:2]] <- summary(g)$coefficients[2, 1:2]
    
    ## meta-regression of reri after bayesian spatial pooling
    m2 <- metagen(TE = w_hat_mu, seTE = w_hat_sd, studlab = zcta, data = bayesDF)
    g2 <- metareg(m2, as.formula(paste0(" ~ ", estimate$vb[i])))
    estimate[i, nms[3:4]] <- c(g2$beta[2], g2$se[2])
    
    ## linear regression of associations before bayesian spatial pooling
    g3 <- lm(as.formula(paste0("est_raw ~ ", estimate$vb[i])), data=bayesDF)
    estimate[i, nms[5:6]] <- summary(g3)$coefficients[2, 1:2]
    
    if (estimate$vb[i]=="AC") {
      bayesDF <- copy(cache)
    }
  }
  estimate$iqr <- estimate$q3 - estimate$q1
  temp <- data.frame(method=m)
  temp <- cbind(rep(temp), estimate)
  out <- rbind(out, temp)
  estimate <- nms <- g <- g2 <- m2 <- temp <- NULL
}
write.csv(out, file.path(outdir1, "results",  paste0(dataset, "_reri_metaregression_022723.csv")), row.names = FALSE)
sink()
##################