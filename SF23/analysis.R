# install swfdr
install.packages("BiocManager")
BiocManager::install("swfdr")

# load the scraped p-values
abstracts <- rbind(readRDS("re-analysis/abstracts.RDS"), readRDS("re-analysis/abstracts_PLOS_Medicine.RDS"))
abstracts <- abstracts[abstracts$journal != "Epidemiologic reviews",]
abstracts <- abstracts[abstracts$journal != "American journal of epidemiology",]

p_values  <- rbind(readRDS("re-analysis/p_values.RDS"),  readRDS("re-analysis/p_values_PLOS_Medicine.RDS"),
                   readRDS("re-analysis/confidence_intervals.RDS"), readRDS("re-analysis/confidence_intervals_PLOS_Medicine.RDS"))

p_values  <- merge(abstracts[,c("pmid", "journal", "year", "RCT", "CT")], p_values, by = "pmid", all.x = FALSE, all.y = TRUE)

p_values$year <- as.numeric(p_values$year)
p_values      <- p_values[!is.na(p_values$year), ]
p_values      <- p_values[p_values$RCT | p_values$CT, ]

xtable::xtable(cbind(
  "Abstracts" = table(abstracts$journal[!duplicated(abstracts$pmid)]),
  "Abstracts" = table(abstracts$journal[!duplicated(abstracts$pmid) & (abstracts$CT | abstracts$RCT)]),
  "Elligible" = table(p_values$journal[!duplicated(p_values$pmid)]),
  "p-values"  = table(p_values$journal)
))

t(t(table(p_values$journal[!duplicated(p_values$pmid)])))
cbind(t(t(table(p_values$journal))), t(t(table(p_values$journal))))
t(t(table(p_values$journal)))
t(t(table(ifelse(p_values$year <= 2010, "<= 2010", "> 2010"))))
table(p_values$journal, ifelse(p_values$year <= 2010, "<= 2010", "> 2010"))
table(p_values$journal, p_values$year)

# some helper functions for z-curve
get_lb <- function(p){
  decimals       <- nchar(p) - 2
  lb             <- ifelse(decimals %in% c(2,3), p - 0.5 * 10^(-decimals), NA)
  lb[!is.na(lb)] <- ifelse(lb[!is.na(lb)] < 0, 0, lb[!is.na(lb)])
  return(lb)
}
get_ub <- function(p){
  decimals       <- nchar(p) - 2
  ub             <- ifelse(decimals %in% c(2,3), p + 0.5 * 10^(-decimals), NA)
  ub[!is.na(ub)] <- ifelse(ub[!is.na(ub)] > 0.05, 0.05, ub[!is.na(ub)])
  return(ub)
}

### collect the data ----
# estimate FDR (combining everything)
job::job({

  data <- p_values

  # select data
  tmp <- data[data$p_value <= 0.05,]

  # create rounding / censoring indicators for swfdr (according to the original script)
  tmp$tt              <- as.numeric(tmp$truncated)
  tmp$rr              <- rep(0,length(tmp$tt))
  tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
  tmp$pp              <- tmp$p_value

  cl   <- parallel::makeCluster(20)
  parallel::clusterExport(cl, c("tmp"), envir = environment())
  fit_swfdr_FDR <- parallel::parSapplyLB(cl, 1:1000, function(i) {
    set.seed(i)
    temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
      this_pmid <- tmp[tmp$pmid == pmid, ]
      return(this_pmid[sample(nrow(this_pmid), 1), ])
    }))
    return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
  })
  parallel::stopCluster(cl)


  # create rounding / censoring bounds for z-curve
  tmp$lb <- get_lb(tmp$p_value)
  tmp$ub <- get_ub(tmp$p_value)

  tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                       ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
  tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]


  cl   <- parallel::makeCluster(20)
  parallel::clusterExport(cl, c("tmp"), envir = environment())
  fit_zcurve <- parallel::parLapplyLB(cl, 1:1000, function(i) {
    set.seed(i)
    temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
      this_pmid <- tmp[tmp$pmid == pmid, ]
      return(this_pmid[sample(nrow(this_pmid), 1), ])
    }))
    temp_fit <- zcurve::zcurve(
      p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
      p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
      p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
      bootstrap = FALSE)
    return(data.frame(t(c(
      ERR = zcurve::ERR(temp_fit)$Estimate,
      EDR = zcurve::EDR(temp_fit)$Estimate,
      FDR = zcurve::Soric(temp_fit)$Estimate,
      weights   = temp_fit$fit$weights,
      prop_high = temp_fit$fit$prop_high
    ))))
  })
  fit_zcurve <- do.call(rbind, fit_zcurve)
  parallel::stopCluster(cl)
  saveRDS(fit_zcurve, file = "re-analysis/fit_zcurve_all.RDS")

  journal_FDR <- data.frame(
    N              = length(unique(tmp$pmid)),
    zcurve_ERR     = median(fit_zcurve$ERR),
    zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
    zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
    zcurve_EDR     = median(fit_zcurve$EDR),
    zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
    zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
    zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
    zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
    zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
    swfdr_FDR      = median(fit_swfdr_FDR),
    swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
    swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
  )

  saveRDS(journal_FDR, file = "re-analysis/FDR.RDS")
},import = "all", title = "all")

# estimate FDR by journal
job::job({

  data <- p_values

  cl   <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR.RDS")
},import = "all", title = "CT + RCT")

job::job({

  data <- p_values

  cl   <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05 & data$CT,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_RCT.RDS")
},import = "all", title = "CT")

job::job({

  data <- p_values

  cl   <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05 & data$RCT,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_RCT.RDS")
},import = "all", title = "RCT")

# estimate FDR by journal past 2010
job::job({

  data <- p_values[p_values$year > 2010, ]

  cl   <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))



    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_2020.RDS")
},import = "all", title = "CT + RCT > 2010")

job::job({
  data <- p_values[p_values$year > 2010 & p_values$CT, ]

  cl <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_CT_2020.RDS")},import = "all", title = "CT > 2010")

job::job({

  data <- p_values[p_values$year > 2010  & p_values$RCT, ]

  cl <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_RCT_2020.RDS")
},import = "all", title = "RCT > 2010")

# estimate FDR by journal before 2010
job::job({

  data <- p_values[p_values$year <= 2010, ]

  cl   <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_2010.RDS")
},import = "all", title = "CT + RCT <= 2010")

job::job({
  data <- p_values[p_values$year <= 2010 & p_values$CT, ]

  cl <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal)[unique(data$journal) != "PLoS medicine"], function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_CT_2010.RDS")},import = "all", title = "CT <= 2010")
#---
job::job({

  data <- p_values[p_values$year <= 2010 & p_values$RCT, ]

  cl <- parallel::makeCluster(5)
  parallel::clusterExport(cl, c("data", "get_lb", "get_ub"), envir = environment())
  journal_FDR <- parallel::parLapplyLB(cl, unique(data$journal), function(journal) {

    # select data
    tmp <- data[data$journal == journal & data$p_value <= 0.05,]

    # create rounding / censoring indicators for swfdr (according to the original script)
    tmp$tt              <- as.numeric(tmp$truncated)
    tmp$rr              <- rep(0,length(tmp$tt))
    tmp$rr[tmp$tt == 0] <- (tmp$p_value[tmp$tt==0] == round(tmp$p_value[tmp$tt==0],2))
    tmp$pp              <- tmp$p_value

    fit_swfdr_FDR       <- sapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      return(swfdr::calculateSwfdr(pValues = temp_sample$pp, truncated = temp_sample$tt, rounded = temp_sample$rr, numEmIterations = 100)$pi0)
    })

    # create rounding / censoring bounds for z-curve
    tmp$lb <- get_lb(tmp$p_value)
    tmp$ub <- get_ub(tmp$p_value)

    tmp$lb[tmp$truncated == 1] <- ifelse(tmp$p_value[tmp$truncated == 1] <= 0.001, 0,
                                         ifelse(tmp$p_value[tmp$truncated == 1] <= 0.01, 0.001, 0.01))
    tmp$ub[tmp$truncated == 1] <- tmp$p_value[tmp$truncated == 1]

    fit_zcurve                 <- do.call(rbind, lapply(1:1000, function(i){
      set.seed(i)
      temp_sample <- do.call(rbind, lapply(unique(tmp$pmid), function(pmid){
        this_pmid <- tmp[tmp$pmid == pmid, ]
        return(this_pmid[sample(nrow(this_pmid), 1), ])
      }))
      temp_fit <- zcurve::zcurve(
        p    = as.vector(temp_sample$p_value)[is.na(temp_sample$lb)],
        p.lb = temp_sample$lb[!is.na(temp_sample$lb)],
        p.ub = temp_sample$ub[!is.na(temp_sample$ub)],
        bootstrap = FALSE)
      return(data.frame(
        ERR = zcurve::ERR(temp_fit)$Estimate,
        EDR = zcurve::EDR(temp_fit)$Estimate,
        FDR = zcurve::Soric(temp_fit)$Estimate
      ))
    }))


    return(data.frame(
      journal        = journal,
      N              = length(unique(tmp$pmid)),
      zcurve_ERR     = median(fit_zcurve$ERR),
      zcurve_ERR_lCI = max(c(quantile(fit_zcurve$ERR, 0.025) - 0.03, 0.025)),
      zcurve_ERR_uCI = min(c(quantile(fit_zcurve$ERR, 0.975) + 0.03, 1)),
      zcurve_EDR     = median(fit_zcurve$EDR),
      zcurve_EDR_lCI = max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)),
      zcurve_EDR_uCI = min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),
      zcurve_FDR     = zcurve:::.get_Soric_FDR(median(fit_zcurve$EDR), 0.05),
      zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(fit_zcurve$EDR, 0.975) + 0.05, 1)),    0.05),
      zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(fit_zcurve$EDR, 0.025) - 0.05, 0.05)), 0.05),
      swfdr_FDR      = median(fit_swfdr_FDR),
      swfdr_FDR_lCI  = quantile(fit_swfdr_FDR, 0.025),
      swfdr_FDR_uCI  = quantile(fit_swfdr_FDR, 0.975)
    ))
  })
  parallel::stopCluster(cl)
  journal_FDR <- do.call(rbind, journal_FDR)
  saveRDS(journal_FDR, file = "re-analysis/journal_FDR_RCT_2010.RDS")
},import = "all", title = "RCT <= 2010")




### analyze the results ----
FDR                  <- readRDS(file = "re-analysis/FDR.RDS")
journal_FDR          <- readRDS(file = "re-analysis/journal_FDR.RDS")
journal_FDR_2010     <- readRDS(file = "re-analysis/journal_FDR_2010.RDS")
journal_FDR_2020     <- readRDS(file = "re-analysis/journal_FDR_2020.RDS")
journal_FDR_CT_2010  <- readRDS(file = "re-analysis/journal_FDR_CT_2010.RDS")
journal_FDR_CT_2020  <- readRDS(file = "re-analysis/journal_FDR_CT_2020.RDS")
journal_FDR_RCT_2010 <- readRDS(file = "re-analysis/journal_FDR_RCT_2010.RDS")
journal_FDR_RCT_2020 <- readRDS(file = "re-analysis/journal_FDR_RCT_2020.RDS")


journal_name <- function(journal){
  return(unname(sapply(journal, function(j) switch(
    j,
    "Lancet (London, England)"             = "Lancet",
    "BMJ (Clinical research ed.)"          = "BMJ",
    "PLoS medicine"                        = "PLoS",
    "The New England journal of medicine"  = "NEJM",
    "JAMA"                                 = "JAMA"
  ))))
}

journal_FDR_2020$journal       <- journal_name(journal_FDR_2020$journal)
journal_FDR_CT_2020$journal    <- journal_name(journal_FDR_CT_2020$journal)
journal_FDR_RCT_2020$journal   <- journal_name(journal_FDR_RCT_2020$journal)
journal_FDR_2010$journal       <- journal_name(journal_FDR_2010$journal)
journal_FDR_CT_2010$journal    <- journal_name(journal_FDR_CT_2010$journal)
journal_FDR_RCT_2010$journal   <- journal_name(journal_FDR_RCT_2010$journal)

### Figures for journal / decade
pdf("re-analysis/fdr_by_journal_and_decade.pdf", width = 12, height = 4)
{
  journals <- c("Lancet", "BMJ","NEJM", "JAMA", "PLoS")
  x_sec    <- seq(0, 1, length.out = 5 + 2)

  par(mfcol = c(1, 2), mar = c(3, 4.5, 2, 1))
  plot(NA, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), main = "Before 2010")
  x_labels <- c(sapply(journals[journals != "PLoS"], function(j) paste0(j, "\n", "(", journal_FDR_2010$N[journal_FDR_2010$journal == j], ")") ), " ")
  axis(1, x_sec[-c(1, length(x_sec))], x_labels, col = NA, col.ticks = NA)
  axis(2, seq(0, 1, 0.2), las = 1)
  mtext("FDR", 2, 3)
  i <- 1
  x_adj <- 0.01
  for(j in journals){

    if(j == "PLoS")
      next

    i <- i + 1
    points(x_sec[i] - x_adj, journal_FDR_2010$zcurve_FDR[journal_FDR_2010$journal == j], pch = 15)
    arrows(
      x0 = x_sec[i] - x_adj, x1 = x_sec[i] - x_adj,
      y0 = journal_FDR_2010$zcurve_FDR_lCI[journal_FDR_2010$journal == j],
      y1 = journal_FDR_2010$zcurve_FDR_uCI[journal_FDR_2010$journal == j],
      code = 3, angle = 90, length = 0.02
    )

    points(x_sec[i] + x_adj, journal_FDR_2010$swfdr_FDR[journal_FDR_2010$journal == j], pch = 1)
    arrows(
      x0 = x_sec[i] + x_adj, x1 = x_sec[i] + x_adj,
      y0 = journal_FDR_2010$swfdr_FDR_lCI[journal_FDR_2010$journal == j],
      y1 = journal_FDR_2010$swfdr_FDR_uCI[journal_FDR_2010$journal == j],
      code = 3, angle = 90, length = 0.02
    )
  }


  plot(NA, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), main = "After 2010")
  x_labels <- sapply(journals, function(j) paste0(j, "\n", "(", journal_FDR_2020$N[journal_FDR_2020$journal == j], ")") )
  axis(1, x_sec[-c(1, length(x_sec))], x_labels, col = NA, col.ticks = NA)
  axis(2, seq(0, 1, 0.2), las = 1)
  mtext("FDR", 2, 3)
  i <- 1
  x_adj <- 0.01
  for(j in journals){

    i <- i + 1
    points(x_sec[i] - x_adj, journal_FDR_2020$zcurve_FDR[journal_FDR_2020$journal == j], pch = 15)
    arrows(
      x0 = x_sec[i] - x_adj, x1 = x_sec[i] - x_adj,
      y0 = journal_FDR_2020$zcurve_FDR_lCI[journal_FDR_2020$journal == j],
      y1 = journal_FDR_2020$zcurve_FDR_uCI[journal_FDR_2020$journal == j],
      code = 3, angle = 90, length = 0.02
    )

    points(x_sec[i] + x_adj, journal_FDR_2020$swfdr_FDR[journal_FDR_2020$journal == j], pch = 1)
    arrows(
      x0 = x_sec[i] + x_adj, x1 = x_sec[i] + x_adj,
      y0 = journal_FDR_2020$swfdr_FDR_lCI[journal_FDR_2020$journal == j],
      y1 = journal_FDR_2020$swfdr_FDR_uCI[journal_FDR_2020$journal == j],
      code = 3, angle = 90, length = 0.02
    )
  }

  legend("topright", legend = c("z-curve", "swfdr"), pch = c(15, 1), bty = "n", cex = 1.20)
  dev.off()
}

pdf("re-analysis/fdr_by_journal_decade_and_type.pdf", width = 12, height = 4)
{
  journals <- c("Lancet", "BMJ","NEJM", "JAMA", "PLoS")
  x_sec    <- seq(0, 1, length.out = 5 + 2)

  par(mfcol = c(1, 2), mar = c(3, 4.5, 2, 1))
  plot(NA, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), main = "Before 2010")
  x_labels <- c(sapply(journals[journals != "PLoS"], function(j) paste0(j, "\n", "(", journal_FDR_RCT_2010$N[journal_FDR_RCT_2010$journal == j], ")\n(", journal_FDR_CT_2010$N[journal_FDR_CT_2010$journal == j], ")") ), " ")
  axis(1, x_sec[-c(1, length(x_sec))], x_labels, col = NA, col.ticks = NA, line = 0.5)
  axis(2, seq(0, 1, 0.2), las = 1)
  mtext("FDR", 2, 3)
  i <- 1
  x_adj <- 0.01
  for(j in journals){

    if(j == "PLoS")
      next

    i <- i + 1
    points(x_sec[i] - x_adj, journal_FDR_RCT_2010$zcurve_FDR[journal_FDR_RCT_2010$journal == j], pch = 15)
    arrows(
      x0 = x_sec[i] - x_adj, x1 = x_sec[i] - x_adj,
      y0 = journal_FDR_RCT_2010$zcurve_FDR_lCI[journal_FDR_RCT_2010$journal == j],
      y1 = journal_FDR_RCT_2010$zcurve_FDR_uCI[journal_FDR_RCT_2010$journal == j],
      code = 3, angle = 90, length = 0.02
    )

    points(x_sec[i] + x_adj, journal_FDR_CT_2010$zcurve_FDR[journal_FDR_CT_2010$journal == j], pch = 0)
    arrows(
      x0 = x_sec[i] + x_adj, x1 = x_sec[i] + x_adj,
      y0 = journal_FDR_CT_2010$zcurve_FDR_lCI[journal_FDR_CT_2010$journal == j],
      y1 = journal_FDR_CT_2010$zcurve_FDR_uCI[journal_FDR_CT_2010$journal == j],
      code = 3, angle = 90, length = 0.02
    )
  }


  plot(NA, type = "n", xlab = "", ylab = "", axes = FALSE, xlim = c(0, 1), ylim = c(0, 1), main = "After 2010")
  x_labels <- sapply(journals, function(j) paste0(j, "\n", "(", journal_FDR_RCT_2020$N[journal_FDR_RCT_2020$journal == j], ")\n(", journal_FDR_CT_2020$N[journal_FDR_CT_2020$journal == j], ")") )
  axis(1, x_sec[-c(1, length(x_sec))], x_labels, col = NA, col.ticks = NA, line = 0.5)

  axis(2, seq(0, 1, 0.2), las = 1)
  mtext("FDR", 2, 3)
  i <- 1
  x_adj <- 0.01
  for(j in journals){

    i <- i + 1
    points(x_sec[i] - x_adj, journal_FDR_RCT_2020$zcurve_FDR[journal_FDR_RCT_2020$journal == j], pch = 15)
    arrows(
      x0 = x_sec[i] - x_adj, x1 = x_sec[i] - x_adj,
      y0 = journal_FDR_RCT_2020$zcurve_FDR_lCI[journal_FDR_RCT_2020$journal == j],
      y1 = journal_FDR_RCT_2020$zcurve_FDR_uCI[journal_FDR_RCT_2020$journal == j],
      code = 3, angle = 90, length = 0.02
    )

    points(x_sec[i] + x_adj, journal_FDR_CT_2020$zcurve_FDR[journal_FDR_CT_2020$journal == j], pch = 0)
    arrows(
      x0 = x_sec[i] + x_adj, x1 = x_sec[i] + x_adj,
      y0 = journal_FDR_CT_2020$zcurve_FDR_lCI[journal_FDR_CT_2020$journal == j],
      y1 = journal_FDR_CT_2020$zcurve_FDR_uCI[journal_FDR_CT_2020$journal == j],
      code = 3, angle = 90, length = 0.02
    )
  }

  legend("topleft", legend = c("RCT", "CT"), pch = c(15, 0), bty = "n", cex = 1.20)
  dev.off()
}

### Table 3
T3 <- rbind(
  cbind(journal = journal_FDR$journal, round(journal_FDR[,c("zcurve_FDR", "zcurve_FDR_lCI", "zcurve_FDR_uCI", "swfdr_FDR", "swfdr_FDR_lCI", "swfdr_FDR_uCI")], 2)),
  c(journal = "Combined", round(FDR[,c("zcurve_FDR", "zcurve_FDR_lCI", "zcurve_FDR_uCI", "swfdr_FDR", "swfdr_FDR_lCI", "swfdr_FDR_uCI")], 2)))
T3

### maximum FDR as a function of alpha
fit_zcurve <- readRDS(file = "re-analysis/fit_zcurve_all.RDS")
alpha      <- c(0.05, 0.01, 0.005, 0.001, 0.0001)

est_alpha <- do.call(rbind, lapply(alpha, function(temp_alpha) {

  new_estimates   <- data.frame(do.call(rbind, lapply(1:nrow(fit_zcurve), function(i){
    zcurve:::.get_estimates(0:6, unlist(fit_zcurve[i, c("weights1", "weights2", "weights3", "weights4", "weights5", "weights6", "weights7")]), fit_zcurve$prop_high[i], temp_alpha, qnorm(0.975))[1:2]
  })))

  new_FDR <- data.frame(
    alpha          = temp_alpha,
    zcurve_ERR     = median(new_estimates$ERR),
    zcurve_ERR_lCI = max(c(quantile(new_estimates$ERR, 0.025) - 0.03, temp_alpha / 2)),
    zcurve_ERR_uCI = min(c(quantile(new_estimates$ERR, 0.975) + 0.03, 1)),
    zcurve_EDR     = median(new_estimates$EDR),
    zcurve_EDR_lCI = max(c(quantile(new_estimates$EDR, 0.025) - 0.05, temp_alpha)),
    zcurve_EDR_uCI = min(c(quantile(new_estimates$EDR, 0.975) + 0.05, 1)),
    zcurve_FDR     = zcurve:::.get_Soric_FDR(median(new_estimates$EDR), temp_alpha),
    zcurve_FDR_lCI = zcurve:::.get_Soric_FDR(min(c(quantile(new_estimates$EDR, 0.975) + 0.05, 1)), temp_alpha),
    zcurve_FDR_uCI = zcurve:::.get_Soric_FDR(max(c(quantile(new_estimates$EDR, 0.025) - 0.05, temp_alpha)), temp_alpha)
  )

  return(new_FDR)
}))
cbind(alpha = as.character(est_alpha$alpha), round(est_alpha[,-1], 2))
