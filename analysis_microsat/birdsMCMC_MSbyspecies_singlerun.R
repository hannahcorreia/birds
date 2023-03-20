## For large nbroods, will likely run into "vector memory exhausted" error
## On Mac systems, 
## Open terminal,and type in the following three lines:
## cd ~
## touch .Renviron
## open .Renviron
## Save the following as the first line of .Renviron:
## R_MAX_VSIZE=100Gb 
## Note: This limit includes both physical and virtual memory.


# clear the workspace
rm(list = ls(all = T))

library(rjags)
library(R2OpenBUGS)
library(coda)
library(extraDistr)
library(ggplot2)

##### specify Ben's post.summ function
post.summ = function(post, var) {
  
  # coerce to matrix for easy subsetting
  post.samp = as.matrix(post)
  
  # if parameter is indexed
  if(substr(var, nchar(var), nchar(var)) == "[") {
    # extract columns with headers equal to the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate desired quantities
    summ = apply(post.sub, 2, function(x) c(mean = mean(x), sd = sd(x), quantile(x, c(0.5, 0.025, 0.975))))
    return(summ)
  }
  
  # if parameter is not indexed
  if(substr(var, nchar(var), nchar(var)) != "[") {
    # extract the column with the same header as the desired variable
    post.sub = post.samp[,substr(colnames(post.samp), 1, nchar(var)) == var]
    # calculate the desired quantities
    summ = c(mean = mean(post.sub), sd = sd(post.sub), quantile(post.sub, c(0.5, 0.025, 0.975)))
    return(summ)
  }
}

#########################################################
#### ESTIMATION OF Q (RATE OF SUCCESS OF EACH SIRE) 
#### AND CALCULATION OF M (# MATES) FROM R (# SIRES)
#########################################################

##### read and prepare data #####
dat1 <- read.csv("paternity_birds_ms.csv")

## remove birds with nbrood = NA and avgbrood = NA (can't run without avgbrood)
dat <- dat1[!is.na(dat1$nbrood) & !is.na(dat1$avgbrood),]

# expand data to create multiple broods within each species (sample size * 10)
for(i in 1:nrow(dat)){
  temp.sp <- expand.grid(species = rep(dat[i,]$species, dat[i,]$nbrood*10))
  if(i == 1){
    temp <- temp.sp
  } else {
    temp <- rbind(temp, temp.sp)
  } 
}
mates.df <- merge(temp, dat, all.x = TRUE, all.y = FALSE)


##### specify model code #####
model_byspecies_q.file <- "model_byspecies_q.txt"
jagsscript.byspecies_q <- cat("
  model{
  # define prior: probability of q (prob of success for each sire)
  q ~ dbeta(1,1)
  
  # likelihood: stochastic relationship
  # sires ~ (truncated binomial)
  for(i in 1:N){
    nsires[i] ~ dbinom(q, littersize[i]) T(0,)
  }
  
  # DERIVED QUANTITIES
  # Calculate the number of mates (deterministic relation?)
  # q is probability of success here
  # for(i in 1:N){
  #   nmates[i] <- (q * nsires[i]) / (1-q)
  # }
  
  # mates ~ (negative binomial); 
  ##
  # NOTE: in R, representation of NB dist is same as Casella-Berger stat theory book
  # whereas, Ash uses the alternative NB pmf seen in the note for NB dist in Casella-Berger
  # Therefore, P(X=x|r,p) used in R is P(X=M-r|r,q) in terms of number of failures M-r
  # mean of NB is μ = n(1-p)/p and variance n(1-p)/p^2, where p is the prob of SUCCESS
  # R is calculating the number of failures M-r, and we generated the number of successes (sires, R)
  # therefore the number of mates M = failures + successes = failed.m + nsires
  ##
  for(i in 1:N){
    failed.m[i] ~ dnbinom(q, nsires[i]) T(nsires[i],)
    nmates[i] <- failed.m[i] + nsires[i]
  }
  for(i in 1:N){
    est.M[i] <- nsires[i]/q
  }
  
  # Calculated quantities from original data
  for(j in 1:J){
    avg.R[j] <- avgbrood[j]*q/(1-(1-q)^avgbrood[j])     # avg sires, given avgbrood and estimated q
  }
  for(j in 1:J){
    avg.M[j] <- avgsire[j]/q     # avg mates, given avgsire and estimated q
  }
  for(j in 1:J){
    est.pmult.kq[j] <- 1 - avgbrood[j] * q * (1 - q)^(avgbrood[j] - 1) / (1 - (1 - q)^(avgbrood[j]))   
    # estimated prob. of mult. paternity, given avgbrood and estimated q
  }
}
",file = model_byspecies_q.file)



#### BEGINNING OF REGENERATION OF DATA LOOP ####
# Remove any datapoints with NA for avgsire, since Poisson data gen will fail
mates.dat <- mates.df[!is.na(mates.df$avgsire),]

set.seed(36849)
  
#### (1) generate broods and sires from Poisson for each species ####
mates <- cbind(mates.dat, littersize = NA, nsires = NA)
for(i in 1:nrow(mates)){
  # generate broods
  nsample <- mates[i,]$nbrood   # number of broods
  avg.bs <- mates[i,]$avgbrood  # average brood size
  gen.k <- rtpois(1, avg.bs, a = 0, b = Inf)
  mates[i,]$littersize <- gen.k
  # generate sires
  avg.sire <- mates[i,]$avgsire  # average number of sires
  sire.max <- mates[i,]$littersize # max number of sires per litter is litter size
  gen.r <- rtpois(1, avg.sire, a = 0, b = sire.max)
  mates[i,]$nsires <- gen.r
}
## remove multiple paternity proportion of 1 or 0
#mates <- mates[mates$pmult != 0 & mates$pmult != 1,]
## remove an average sire value of NA
mates <- mates[!is.na(mates$avgsire), ]

## save data generation (b/c it takes too long to redo each time)
# save(mates, file = "mates.rda")
# load("mates.rda")

# ##### (2) initial values for Q #####
# # Calculate an approximate initial q based on p (multiple paternity)
# ## Negative log-likelihood of zero-trunc binomial
# truncbin.lik <- function(q, s, k){
#   logl <- sum(s)*log(q) + sum(k-s)*(log(1-q)) - sum(log(1 - (1-q)^k))
#   return(-logl)
# }
# fq <- function(q, k, p) 1 - k*q*(1-q)^(k-1)/(1-(1-q)^k) - p
# mates.p <- mates
# mates.p$q0 <- mates.p$qs <- NA
# for(i in sort(unique(mates.p$species))){
#   for(j in 1:nrow(mates.p)){
#     mates.p[j,]$q0 <- uniroot(fq, c(0.01,.99), k=mates.p[j,]$littersize, p=mates.p[j,]$pmult)$root    # estimated prob success per sire
#   }
#   mates.p[mates.p$species==i,]$qs <- optim(.5, truncbin.lik, s=mates[mates.p$species==i,]$nsires, k=mates[mates.p$species==i,]$littersize, method = "Brent", hessian = TRUE, lower=0, upper=1)$par
# }
# q.init1 <- mean(mates.p$qs)
#   
# # obtain different q.init for 2nd chain
# q.init2 <- mean(sample(mates.p$qs, size = 25, replace = TRUE))
# #q.init2 <- mean(sample(mates.p$qs, size = 10, replace = TRUE))
# #q.init2 <- mean(sample(mates.p$qs, size = 10, replace = FALSE))
# #q.init2 <- 0.5084281
# #q.init2 <- 0.6295649
#   
# inits <- list(list(q=q.init1), list(q=q.init2))
  
##### (3) MCMC dimensions #####
ni = 10000
nb = 1000
nc = 2 # needs to match number of initial values proposed
nt = 2
n.iter = ni + nb
  
##### (4) parameters to monitor #####
params = c("q", "nsires", "failed.m", "nmates", "est.M", "avg.R", "avg.M", "est.pmult.kq")
  
  
##### CREATE LIST TO HOLD ALL RESULTS FROM JAGS
# # max number of iterations through all combinations of sires and litters
# maxlitter.for.minsire <- length(min(mates.dat$nsires):max(mates.dat$littersize))
# maxlitter.for.maxsire <- length(max(mates.dat$nsires):max(mates.dat$littersize))
# nlist <- sum(maxlitter.for.minsire:maxlitter.for.maxsire) 
# # slightly more slots in list than necessary, e.g. sires 5 & 6 only have max of 13 littersize
# # could change with random generation of litter sizes though
jags.results <- list() 
  
##### (5) run the model in JAGS #####
# Initial conditions for loop
# Likelihood needs to change for each nsire value, so q calculated for each # nsire
mates.p <- mates
species <- sort(unique(mates.p$species))
b <- 1
# Start JAGS
#starttime <- Sys.time()
for(r in species){
  # data containing only relevant values
  m.dat <- mates.p[mates.p$species==r,]
  orig.dat <- dat[dat$species==r,]
  jags.dat <- list(nsires = m.dat$nsires, 
                   littersize = m.dat$littersize, 
                   N = nrow(m.dat), 
                   avgbrood = orig.dat$avgbrood,
                   avgsire = orig.dat$avgsire,
                   J = nrow(orig.dat))
  
  # run JAGS (not using initial values for q)
  jmod <- jags.model(file = model_byspecies_q.file, data = jags.dat, n.chains = nc, n.adapt = 1000)  
  update(jmod, n.iter = nb, by = 1, progress.bar = 'text')
  post <- coda.samples(jmod, params, n.iter = ni, thin = nt) 
      
  # save current JAGS output
  name <- paste0("species=",r)
  jags.results[[name]] <- post
      
  # Save quantities from data
  MCMC_temp <- cbind(
    as.character(orig.dat$species),
    as.numeric(orig.dat$avgbrood),
    as.numeric(orig.dat$avgsire),
    as.numeric(orig.dat$pmult),
    post.summ(post, "q")[[1]],  # mean of the est. param. q
    post.summ(post, "q")[[2]],  # sd of the est. param. q
    post.summ(post, "q")[[4]],  # 2.5% est. param. q
    post.summ(post, "q")[[5]],  # 97.5% est. param. q
    post.summ(post, "nsires")[[1]],  # mean of nsires
    post.summ(post, "nsires")[[2]],  # sd of nsires
    post.summ(post, "nsires")[[4]],  # 2.5% nsires
    post.summ(post, "nsires")[[5]],  # 97.5% nsires
    post.summ(post, "nmates")[[1]],  # mean of nmates
    post.summ(post, "nmates")[[2]],  # sd of nmates
    post.summ(post, "nmates")[[4]],  # 2.5% nmates
    post.summ(post, "nmates")[[5]],  # 97.5% nmates
    post.summ(post, "est.M")[[1]],  # mean of est.M
    post.summ(post, "est.M")[[2]],  # sd of est.M
    post.summ(post, "est.M")[[4]],  # 2.5% est.M
    post.summ(post, "est.M")[[5]],  # 97.5% est.M
    post.summ(post, "est.pmult.kq")[[1]],  # mean of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[2]],  # sd of the est.pmult.kq
    post.summ(post, "est.pmult.kq")[[4]],  # 2.5% est.pmult.kq
    post.summ(post, "est.pmult.kq")[[5]],  # 97.5% est.pmult.kq
    post.summ(post, "avg.R")[[1]],  # mean of avg.R
    post.summ(post, "avg.R")[[2]],  # sd of avg.R
    post.summ(post, "avg.R")[[4]],  # 2.5% avg.R
    post.summ(post, "avg.R")[[5]],  # 97.5% avg.R
    post.summ(post, "avg.M")[[1]],  # mean of avg.M
    post.summ(post, "avg.M")[[2]],  # sd of avg.M
    post.summ(post, "avg.M")[[4]],  # 2.5% avg.M
    post.summ(post, "avg.M")[[5]]  # 97.5% avg.M
    )
   
if(b==1) MCMC_sumtemp <- MCMC_temp else MCMC_sumtemp <- rbind(MCMC_sumtemp, MCMC_temp)
  
# move to next step
b <- b + 1
}
  
MCMC_summary <- data.frame(species = as.character(MCMC_sumtemp[,1]),
                             avgbrood = as.numeric(MCMC_sumtemp[,2]),
                             avgsire = as.numeric(MCMC_sumtemp[,3]),
                             pmult = as.numeric(MCMC_sumtemp[,4]),
                             mean.q = as.numeric(MCMC_sumtemp[,5]),
                             sd.q = as.numeric(MCMC_sumtemp[,6]),
                             q2.5 = as.numeric(MCMC_sumtemp[,7]),
                             q97.5 = as.numeric(MCMC_sumtemp[,8]),
                             mean.nsires = as.numeric(MCMC_sumtemp[,9]),
                             sd.nsires = as.numeric(MCMC_sumtemp[,10]),
                             nsires2.5 = as.numeric(MCMC_sumtemp[,11]),
                             nsires97.5 = as.numeric(MCMC_sumtemp[,12]),
                             mean.nmates = as.numeric(MCMC_sumtemp[,13]),
                             sd.nmates = as.numeric(MCMC_sumtemp[,14]),
                             nmates2.5 = as.numeric(MCMC_sumtemp[,15]),
                             nmates97.5 = as.numeric(MCMC_sumtemp[,16]),
                             mean.estM = as.numeric(MCMC_sumtemp[,17]),
                             sd.estM = as.numeric(MCMC_sumtemp[,18]),
                             estM2.5 = as.numeric(MCMC_sumtemp[,19]),
                             estM97.5 = as.numeric(MCMC_sumtemp[,20]),
                             mean.est.pmult = as.numeric(MCMC_sumtemp[,21]),
                             sd.est.pmult = as.numeric(MCMC_sumtemp[,22]),
                             est.pmult2.5 = as.numeric(MCMC_sumtemp[,23]),
                             est.pmult97.5 = as.numeric(MCMC_sumtemp[,24]),
                             mean.avgR = as.numeric(MCMC_sumtemp[,25]),
                             sd.avgR = as.numeric(MCMC_sumtemp[,26]),
                             avgR2.5 = as.numeric(MCMC_sumtemp[,27]),
                             avgR97.5 = as.numeric(MCMC_sumtemp[,28]),
                             mean.avgM = as.numeric(MCMC_sumtemp[,29]),
                             sd.avgM = as.numeric(MCMC_sumtemp[,30]),
                             avgM2.5 = as.numeric(MCMC_sumtemp[,31]),
                             avgM97.5 = as.numeric(MCMC_sumtemp[,32]))
  
save(MCMC_summary, file = paste0("MCMC_summary_singlerun.rda"))
  

# plot single run with original pmult data (from Avis)
g1a <- ggplot() +
  stat_smooth(data = MCMC_summary, aes(avgbrood, mean.est.pmult), color = "red", method = "loess") +
  geom_point(data = MCMC_summary, aes(avgbrood, mean.est.pmult), color = "red", alpha = 0.2, size = 2.5) +
  geom_point(data = dat, aes(avgbrood, pmult)) +
  #stat_smooth(data = dat, aes(avgbrood, pmult), color = "red") +
  labs(x="Brood size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  theme_bw(base_size = 16)
ggsave("p_singlerun.eps", plot = g1a, device = cairo_ps, width = 11, height = 7.5, dpi = 320)



# ################ DIAGNOSTICS BELOW ###################
# current <- 1 #change for different species (61 total)
# nam.current <- names(jags.results[current])
# jags.post <- jags.results[[nam.current]]
# gelman.diag(jags.post, multivariate = F) # values below 1.1 should be okay
# gelman.plot(jags.post, multivariate = F)
# n.eff <- effectiveSize(jags.post)
# 
# # visualize trace and posterior plots
# par(mar=c(2,2,1,1))
# plot(jags.post)
# 
# ##### make inference #####
# post.summ(jags.post, "q") # probability of success for each sire for a given data-gen run
# post.summ(jags.post, "nsires")
# post.summ(jags.post, "nmates")
# post.summ(jags.post, "avg.R") # avg number of sires given the est prob of success for each sire, using avg litter size
# post.summ(jags.post, "avg.M")
# post.summ(jags.post, "est.pmult.kq")

# #### Compare the est.sires (calculated with est. q) to true.sires
# sires.results <- sires.results[order(sires.results$species),]
# head(sires.results)
# for(i in unique(sires.results$true.sires)){
#   print(apply(sires.results[sires.results$true.sires==i, 2:5], 2, mean))
# }


################ CALCULATE RESIDUALS ###################
MCMC_resids <- MCMC_summary
loess_pmult <- loess(mean.est.pmult ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult <- predict(loess_pmult, newdata = MCMC_resids$avgbrood)
# lower limit
loess_pmult_lcl <- loess(est.pmult2.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult_lcl <- predict(loess_pmult_lcl, newdata = MCMC_resids$avgbrood)
# upper limit
loess_pmult_ucl <- loess(est.pmult97.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_pmult_ucl <- predict(loess_pmult_ucl, newdata = MCMC_resids$avgbrood)
# residual calculations
MCMC_resids$resid_pmult <- MCMC_resids$pred_pmult - MCMC_resids$pmult
MCMC_resids$resid_pmult_lcl <- MCMC_resids$pred_pmult_lcl - MCMC_resids$pmult
MCMC_resids$resid_pmult_ucl <- MCMC_resids$pred_pmult_ucl - MCMC_resids$pmult

# save dataframe with all residual calculations
save(MCMC_resids, file = "MCMC_resids.rda")
write.csv(MCMC_resids, file = "MCMC_resids.csv")

# plot pred_pmult with lower and upper limits
library(wesanderson)

g1b <- ggplot() +
  geom_point(data = dat, aes(avgbrood, pmult), alpha = 0.5, size = 2.5) +
  #### below is for regular "p_CL.eps" 
  # stat_smooth(data = dat, aes(avgbrood, pmult), se = FALSE, color = "black", method = "lm") +
  # stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult), color = "red", se = FALSE, method = "loess") +
  # geom_point(data = MCMC_resids, aes(avgbrood, mean.est.pmult), color = "red", alpha = 0.5, size = 2.5) +
  # stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  # stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  #### below is for "p_CL_alt.eps" plot for publication
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult), color = wes_palette("Darjeeling1")[1], se = FALSE, method = "loess") +
  geom_point(data = MCMC_resids, aes(avgbrood, mean.est.pmult), color = wes_palette("Darjeeling1")[1], alpha = 0.5, size = 2.5) +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_ucl), color = wes_palette("Darjeeling1")[2], se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_pmult_lcl), color = wes_palette("Darjeeling1")[2], se = FALSE, size = 0.5, method = "loess") +
  labs(x="Brood size", y="Probability of multiple paternity") +
  lims(y = c(0,1)) +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  theme_bw(base_size = 16)
ggsave("p_CL_alt.eps", plot = g1b, device = cairo_ps, width = 11, height = 7.5, dpi = 320)


################ CALCULATE EFFECT SIZES ###################
birds_1a <- read.csv("paternity_birds_ms.csv", header = TRUE)
## remove birds with nbrood = NA
birds_1b <- birds_1a[!is.na(birds_1a$nbrood),]
birds_1 <- birds_1b[!is.na(birds_1b$avgsire),]

# yi = Bayes pmult - true pmult (i.e. MCMC_total_resids$resid_pmult)
# vi = pmult*(1-pmult)/nbrood

birds_singlerun1 <- merge(birds_1[,c(1:3)], MCMC_resids, by = c("species", "avgbrood"))
# order by brood size
birds_singlerun <- birds_singlerun1[order(birds_singlerun1$avgbrood, birds_singlerun1$pmult, birds_singlerun1$nbrood),]

### variance p(1-p)/n for fixed p assumption
birds_singlerun$vi <- birds_singlerun$pmult*(1-birds_singlerun$pmult)/birds_singlerun$nbrood

## labels
bird_species <- birds_singlerun[order(birds_singlerun$avgbrood),]$species

library(metafor)
paternity.meta.bayes <- rma(resid_pmult, vi, data = birds_singlerun, slab = bird_species)
paternity.meta.bayes

pdf(file = paste0("forest_bayes_singlerun-w-labels.pdf"), width = 10, height = 25)
print(
  forest(paternity.meta.bayes,  
         order = order(birds_singlerun$avgbrood), cex = 1,
         xlim = c(-2.5,2), 
         xlab = expression(paste(p[B], " - p"))),
  quote = FALSE,
  text(-2.5, 139, "Species",  pos=4),
  text(2, 139, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()

#### Need to plot forest plot over three pages (way too long otherwise)
## Part 1
res1 <- paternity.meta.bayes
res1$vi.f <- res1$vi.f[1:69]
res1$yi.f <- res1$yi.f[1:69]
res1$slab <- bird_species[1:69]

# part1.species <- birds_singlerun[order(birds_singlerun$avgbrood),]$species[1:68]

## Species should already be ordered by avgbrood and labelled correctly - see above
pdf(file = paste0("forest_bayes_singlerun1-w-labels.pdf"), width = 9, height = 16)
print(
  forest(res1, #slab = part1.species, 
         #order = order(birds_singlerun[birds_singlerun$species %in% part1.species,]$avgbrood), 
         cex = 1,
         xlim = c(-2.5,2), alim = c(-0.6, 0.8), 
         xlab = expression(paste(p[B], " - p")),
         mlab = ""),  ## mlab controls the "RE Model" text at the bottom left
  quote = FALSE,
  text(-2.5, 71, "Species",  pos=4),
  text(2, 71, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()

## Part 2
res2 <- paternity.meta.bayes
res2$vi.f <- res2$vi.f[70:138]
res2$yi.f <- res2$yi.f[70:138]
res2$slab <- bird_species[70:138]

# part2.species <- birds_singlerun[order(birds_singlerun$avgbrood),]$species[70:138]

pdf(file = paste0("forest_bayes_singlerun2-w-labels.pdf"), width = 9, height = 16)
print(
  forest(res2, #slab = part2.species, 
         #order = order(birds_singlerun[birds_singlerun$species %in% part2.species,]$avgbrood), 
         cex = 1,
         xlim = c(-2.5,2), alim = c(-0.6, 0.8), 
         xlab = expression(paste(p[B], " - p")),
         mlab = ""),  ## mlab controls the "RE Model" text at the bottom left
  quote = FALSE,
  text(-2.5, 71, "Species",  pos=4),
  text(2, 71, expression(paste(p[B], " - p [95% CI]")), pos=2)
)
dev.off()

# ## Part 3
# res3 <- paternity.meta.bayes
# res3$vi.f <- res3$vi.f[92:136]
# res3$yi.f <- res3$yi.f[92:136]
# res3$slab <- bird_species[92:136]
# 
# # part3.species <- birds_singlerun[order(birds_singlerun$avgbrood),]$species[163:243]
# 
# pdf(file = paste0("forest_bayes_singlerun3-w-labels.pdf"), width = 15, height = 10)
# print(
#   forest(res3, #slab = part3.species, 
#          #order = order(birds_singlerun[birds_singlerun$species %in% part3.species,]$avgbrood), 
#          cex = 1,
#          xlim = c(-2.5,2),
#          xlab = expression(paste(p[B], " - p"))),
#   quote = FALSE,
#   text(-2.5, 47, "Species",  pos=4),
#   text(2, 47, expression(paste(p[B], " - p [95% CI]")), pos=2)
# )
# dev.off()

# k as moderator 
paternity.meta.k <- rma(resid_pmult, vi, mods = ~ avgbrood, data = birds_singlerun)
paternity.meta.k
pdf(file = paste0("forest_bayes_singlerun_k-mod.pdf"), width = 14, height = 25)
print(
  forest(paternity.meta.k, slab = birds_singlerun$species, 
         order = order(birds_singlerun$avgbrood), cex = 1)
)
dev.off()

# k as moderator (data without pmult = 0 or 1)
paternity.meta.k <- rma(resid_pmult, vi, mods = ~ avgbrood, data = birds_singlerun[!birds_singlerun$pmult==0 & !birds_singlerun$pmult==1,])
paternity.meta.k
forest(paternity.meta.k, slab = birds_singlerun[!birds_singlerun$pmult==0 & !birds_singlerun$pmult==1,]$species, 
       order = order(birds_singlerun[!birds_singlerun$pmult==0 & !birds_singlerun$pmult==1,]$avgbrood), cex = 1)


################ E(M) AND E(S) VS K PLOT ###################
## Plot estimated M from data using Bayesian p overlaid on observed avg M vs obs. littersize
loess_R <- loess(mean.avgR ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_R <- predict(loess_R, newdata = MCMC_resids$avgbrood)
# lower limit
loess_R_lcl <- loess(avgR2.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_R_lcl <- predict(loess_R_lcl, newdata = MCMC_resids$avgbrood)
# upper limit
loess_R_ucl <- loess(avgR97.5 ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_R_ucl <- predict(loess_R_ucl, newdata = MCMC_resids$avgbrood)

g2 <- ggplot() +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_R), color = "red", se = FALSE, method = "loess") +
  geom_point(data = MCMC_resids, aes(avgbrood, mean.avgR), color = "red", alpha = 0.5, size = 2.5) +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_R_ucl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_R_lcl), color = "blue", se = FALSE, size = 0.5, method = "loess") +
  stat_smooth(data = birds_1, aes(avgbrood, avgsire), color = "black", method = "lm", se = FALSE) +
  geom_point(data = birds_1, aes(avgbrood, avgsire), alpha = 0.5, size = 2.5) +
  labs(x="Brood size", y="Number of sires") +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  theme_bw(base_size = 16)
ggsave("R_CL.eps", plot = g2, device = cairo_ps, width = 11, height = 7.5, dpi = 320)


## Replicate Ash's plot of #mates/sires vs brood size with Bayesian estimates

# Combinatorics avgmates calculated using pred_pmult and true k
birds_1$avgmate <- exp(log(1-birds_1$pmult)/(1-birds_1$avgbrood))

# Used the Avise data avgbrood and avgsire to obtain avgR = k*q/(1-(1-q)^k) and avgM = r/q
# calculate predicted M
loess_M <- loess(mean.avgM ~ avgbrood, data = MCMC_resids)
MCMC_resids$pred_M <- predict(loess_M, newdata = MCMC_resids$avgbrood)
# plot
g3 <- ggplot() +
  # avgR in red is calculated from mean of ZTB dist, since S ~ ZTB
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_R, linetype = "Estimated # sires using\ndata K and Bayesian q"), color = "red", method = "loess", se = TRUE) +
  # points from which the above line was generated:
  #geom_point(data = MCMC_resids, aes(avgbrood, mean.avgR), color = "red", alpha = 0.2) +
  # True number of average sires from the Avise data
  geom_point(data = birds_1, aes(avgbrood, avgsire, color = "Observed number of sires"), alpha = 0.5, size = 2.5) +
  # Avise used avgsires as estimated mates (i.e. he does not distinguish b/w mates and sires?)
  stat_smooth(data = birds_1, aes(avgbrood, avgsire, linetype = "Observed number of sires"), color = "black", method = "lm", se = FALSE) +
  # avgM in blue is calculated from mean of NB dist knowing S and q, since M ~ NB
  stat_smooth(data = MCMC_resids, aes(avgbrood, pred_M, linetype = "Estimated # mates using\ndata S and Bayesian q"), color = "blue", method = "lm", se = TRUE) +
  # points from which the above line was generated:
  #geom_point(data = MCMC_resids, aes(avgbrood, mean.avgM), color = "blue", alpha = 0.2) +
  # could instead plot avgmates calculated using true pmult and true k from combinatorics formula?
  stat_smooth(data = birds_1[!is.infinite(birds_1$avgmate),], aes(avgbrood, avgmate, linetype = "Estimated # mates using\ncombinatorial formula"), 
              color = "green3", method = "lm", se = TRUE) +
  geom_point(data = birds_1[!is.infinite(birds_1$avgmate),], aes(avgbrood, avgmate, color = "Estimated # mates using\ncombinatorial formula"), alpha = 0.5, size = 2.5) +
  labs(x="Brood size", y="Number of mates/sires") +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_y_continuous(limits = c(0,10)) +
  scale_colour_manual(name = "", values = c("Observed number of sires" = "black",
                                            "Estimated # mates using\ncombinatorial formula" = "green3")) +
  scale_linetype_manual(values = c("Estimated # mates using\ndata S and Bayesian q" = 1, # blue line
                                   "Estimated # mates using\ncombinatorial formula" = 1, # dark green line
                                   "Estimated # sires using\ndata K and Bayesian q" = 1, # red line
                                   "Observed number of sires" = 1), # black line
                        name = "",
                        guide = guide_legend(override.aes = list(color = c("blue", "green3", "red", "black"), 
                                                                 fill = c(NA,NA,NA,NA)))) +
  theme_bw(base_size = 16) + theme(legend.key.height=unit(15, "mm"))
ggsave("Avise_MS-k.eps", plot = g3, device = cairo_ps, width = 12, height = 7.5, dpi = 320)


# Used the generated data to create nsires and calculate estM = nsires/q
g4 <- ggplot() +
  # NB estimated mates using different values of q
  stat_function(fun = function(x) x/(1-0.2^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
  stat_function(fun = function(x) x/(1-0.3^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
  stat_function(fun = function(x) x/(1-0.4^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
  stat_function(fun = function(x) x/(1-0.5^x), aes(x = x, linetype = "Expected # mates from NB distribution"), data = data.frame(x = c(2,10)), color = "magenta", size = 1.2) +
  # # nsires in red is generated from Poission dist with lambda = avgsire in Avise data
  # stat_smooth(data = MCMC_resids, aes(avgbrood, mean.nsires, linetype = "Generated # sires using Poisson(avgsire)"), color = "red", method = "lm", se = TRUE) +
  # # points from which the above line was generated:
  # geom_point(data = MCMC_resids, aes(avgbrood, mean.nsires, color = "Generated # sires using Poisson(avgsire)"), alpha = 0.5) +
  # True number of average sires from the Avise data
  geom_point(data = birds_1, aes(avgbrood, avgsire, color = "Observed number of sires"), alpha = 0.5, size = 2.5) +
  # Avise used avgsires as estimated mates (i.e. he does not distinguish b/w mates and sires?)
  stat_smooth(data = birds_1, aes(avgbrood, avgsire, linetype = "Observed number of sires"), color = "black", method = "lm", se = FALSE) +
  # estM in blue is calculated from mean of NB dist knowing generated S and q, since M ~ NB
  stat_smooth(data = MCMC_resids, aes(avgbrood, mean.estM, linetype = "Bayesian MCMC estimated # mates"), color = "blue", method = "lm", se = FALSE) +
  # points from which the above line was generated:
  geom_point(data = MCMC_resids, aes(avgbrood, mean.estM, color = "Bayesian MCMC estimated # mates"), alpha = 0.5, size = 2.5) +
  # could instead plot avgmates calculated using true pmult and true k from combinatorics formula?
  stat_smooth(data = birds_1[!is.infinite(birds_1$avgmate),], aes(avgbrood, avgmate, linetype = "Estimated # mates using\ncombinatorial formula"), 
              color = "green3", method = "lm", se = FALSE) +
  geom_point(data = birds_1[!is.infinite(birds_1$avgmate),], aes(avgbrood, avgmate, color = "Estimated # mates using\ncombinatorial formula"), alpha = 0.5, size = 2.5) +
  labs(x="Brood size", y="Number of mates/sires") +
  scale_x_continuous(breaks = c(seq(1,12, by = 1)) , limits = c(1,12)) +
  scale_y_continuous(breaks = c(seq(0,13, by = 2)) ,limits = c(0,13)) +
  scale_colour_manual(values = c("Bayesian MCMC estimated # mates" = "blue",
                                 "Estimated # mates using\ncombinatorial formula" = "green3",
                                 #"Generated # sires using Poisson(avgsire)" = "red",
                                 "Observed number of sires" = "black"),
                      name = "") +
  scale_linetype_manual(values = c("Bayesian MCMC estimated # mates" = 1, # blue line,
                                   "Estimated # mates using\ncombinatorial formula" = 1, # green3 line
                                   "Expected # mates from NB distribution" = 1, # magenta line
                                   #"Generated # sires using Poisson(avgsire)" = 1, # red line
                                   "Observed number of sires" = 1), # black line
                        name = "",
                        guide = guide_legend(override.aes = list(color = c("blue","green3", 
                                                                           "magenta", 
                                                                           #"red", 
                                                                           "black"), 
                                                                 fill = c(NA,NA,
                                                                          NA,
                                                                          #NA,
                                                                          NA)))) +
  #theme_bw(base_size = 16) + theme(legend.key.height=unit(10, "mm")) + # For legend on right
  theme_bw(base_size = 16) + theme(legend.key.height=unit(10, "mm"), 
                                   legend.direction = "vertical",
                                   legend.position = "bottom")
#ggsave("Bayes_gen_MS-k.eps", plot = g4, device = cairo_ps, width = 12, height = 7.5, dpi = 320) #Legend on right
# Edited plot with legend on bottom for
ggsave("g4plot.eps", plot = g4, device = cairo_ps, width = 8, height = 10, dpi = 320) #Legend on bottom
