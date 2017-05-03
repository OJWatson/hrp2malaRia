
#' Individual Based Malaria Simulator
#'
#' \code{hrp2_Simulation} simulates an individual based malaria system considering
#' multiplicity of infection and infection status of individuals.
#'
#' The simulator considers individuals to progress through the Griffin 2016 model.
#'
#' The model also introduces a further acknowledgement of the hrp2 strains present in hosts, and
#' the effects of this on detection.
#'
#' The model first allocates a fixed length for all storage and proceeds in a ring like fashion allowing for delays to
#' be incorporated.
#'
#' ## Simulation parameters
#' @param years Length of Simulation (years). Default = 0.5
#' @param time.step Simulation time interval considered in days. Has to be 1.
#' @param states Vector of the compartmental states. Default = c(1, 2, 3, 4, 5, 6), which are c("S","D","A","U","T","P")
#' @param storage Number of last time steps to store. If storage = 0 (default), then all results are stored
#' ## Demographic paramaters
#' @param N Population Size. Default = 1000
#' @param max.age Maximum age (years). Default = 100
#' @param prev.0 Assumed starting prevalence. Default = 0.3
#' @param strains.0 Maximum starting MOI possible. Default = 5
#' @param desired.freq Desired hrp2 deletion frequency to be enforced when a previous simulation is reloaded
#' @param d.death Average life expectancy (years). Deafult = 21
#' ## Epidemiological Parameters
#' @param EIR Maximum EIR / day. Default = 200/365,
#' @param a.0 Age-dependent biting parameter (days). Default = 8*365
#' @param rho Age-dependent biting parameter (days). Default = 0.85
#' @param zeta.sd sd of lognormal distribution used to generate indiviudal specific biting. Default = 1.67
#' @param omega BM age scaling parameter. Default = 0.75
#' ## Entomological parameters
#' @param mu.0 Baseline mosquito mortality. Default = 0.132
#' @param beta Birthrate in form of y = beta[1]*EIR + beta[2]. Default = c(0.0572289,0.0696121)
#' @param theta Seasonal effects. Default = 1
#' @param a.k Daily biting rate. Default = 0.30667
#' ## Diagnostics and nmf
#' @param ft Probability of receiving assesment. Default = 0.4
#' @param rdt.det Probability of detecting malaria given mono infection of hrp2' strains. Default = 1
#' @param rdt.nonadherence Probability of non adherence to test result. Default = 0, i.e. perfect adherence
#' @param microscopy.use Probability of microscopy use. Default = 0, i.e. no microscopy
#' @param include.nmf Boolean detailing whether to include nmf section. Default = FALSE
#' @param fever.age.breaks Age breaks for fevers
#' @param annual.nmf.per.age.bracket Annual number of non malarial fevers for each age bracket in fever.age.breaks.
#' Default is mean across 5 representative surveys from Burundi 2012, Liberia 2009/2011 and Nigeria 2010/2015
#' @param nmf.multiplier Multiplication for annual.nmf.per.age.bracket to introduce sensitivity. Default = 1
#' ## delays and durations
#' @param delay.mos Mosquito extrinsic incubation period (days). Default = 10
#' @param delay.gam Delay from emergence of blood-stage parasites to onward infectivity (days). Default = 12.5
#' @param d.E Latent period delay in pre-erythrocytic stage (days). Default = 12
#' @param d.I Patent infection duration with no immunity (including disease) (days). Default = 200
#' @param d.T Clinical disease duration with treatment (days). Default  = 5
#' @param d.D Clinical disease duration without treatment (days). Default  = 5
#' @param d.U Sub-patent infection duration (days). Default = 100
#' @param d.P Prophylaxis with SP following treatment duration (days). Default = 25
#' @param d.A Asymptomatic infection duration (days). Default = d.I - d.D
#' ## Pre-erythrocytic immunity parameters
#' @param d.1 = 0.161 Probability with maximum immunity
#' @param d.ID = 10*365 Inverse of decay rate
#' @param I.D0 Scale parameter. Default = 1.58
#' @param k.D Shape parameter. Default = 0.477
#' @param u.D Duration in which immunity is not boosted (days). Default = 9.45
#' @param a.D Scale parameter relating age to immunity (days). Default = 21.923*365
#' @param f.D0 Parameter relating age to immunity. Default = 0.0071
#' @param g.D Shape parameter relating age to immunity. Default = 4.81
#' ## Blood stage immunity parameters
#' @param b.0 Probability with no immunity. Default = 0.590
#' @param b.1 Maximum relative reduction. Default = 0.5
#' @param d.B Inverse of decay rate (days). Default = 10*365
#' @param I.B0 Scale parameter. Default = 43.9
#' @param k.B  Shape parameter. Default = 2.16
#' @param u.B Duration in which immunity is not boosted (days). Default = 7.2
#' ## Acquired and maternal immunity parameters
#' @param phi.0 Probability with no immunity. Default = 0.7916
#' @param phi.1 Maximum relative reduction. Default = 0.00074
#' @param d.CA Inverse of decay rate (days). Default = 30*365
#' @param I.C0 Scale parameter. Default = 18.0
#' @param k.C Shape parameter. Default = 2.37
#' @param u.C Duration in which immunity is not boosted (days). Default = 6.06
#' @param P.M New-born immunity relative to mother’s. Default = 0.774
#' @param d.CM Inverse of decay rate of maternal immunity (days). Default = 67.7
#' ## contributions to infectious reservoir by state and age
#' @param c.S Contribution from S. Default = 0
#' @param c.D Contribution from D. Default = 0.068
#' @param c.T Contribution from T. Default = 0.32 * c.D
#' @param c.U Contribution from U. Default = 0.0062
#' @param c.A Contribution from A. Default = c.U + (c.D - c.U)*q^g.1, where q is
#' the probability of detection of parasites as calculated using pre-erythrcytic immunity and age
#' @param c.P Contribution from P. Default = 0
#' @param g.1 Relates infectiousness probability of detection. Default = 1.82
#'
#' @return Result a list of 59 elements including the following
#' Time - Vector of time steps considered
#' Status - Matrix of individuals' infection status by Time
#' Last_Infection_Time - Vector of individual's last infection
#' Last_Bite_Time - Vector of individual's last bite
#' Age - Vector of individual's ages
#' N.Sens - Matrix of individuals' MOI of RDT sensitive strains by Time
#' N.Dels - Matrix of individuals' MOI of RDT resistant strains by Time
#' I.B - Vector of individual's blood stage immunity
#' I.CA - Vector of individual's acquired immunity
#' I.CM - Vector of individual's maternal immunity
#' I.D - Vector of individual's pre-erythrocytic immunity
#' Counter - Simulation step counter
#' d.EIR - daily EIR in individuals > 18 years
#' I.Reservoir - Vector of each strains weighted contribution to the last step's I.Reservoir
#' Sv - Susceptible mosquito population
#' Ev - Exposed mosquito population
#' Iv - Infected mosquito population
#' Prev - Daily prevalence
#' N.Sens/N.Dels - sum of sensitive to diagnostic and deletion strains
#' Incidence - Clinical cases / population / per day
#' Mono - monoinfected individual proportions
#' FOIv - Force of infection to vectors
#' S. ... - Series storage of various of the above


hrp2_Simulation <- function(
  ## simulation parameters
  ID=NULL,
  root=NULL,
  direct.input = NULL,
  years = 0.5,
  time.step = 1,
  states = c(1, 2, 3, 4, 5, 6),
  storage = 50,
  #states = c("S","D","A","U","T","P"),
  ## demographic parameters
  N = 1000,
  max.age = 60,
  prev.0 = 0.3,
  strains.0 = 9,
  desired.freq = NULL,
  d.death = 21,
  ## epidemiological parameters
  EIR = 200/365,
  a.0 = 8*365,
  rho = 0.85,
  zeta.sd = 1.67,
  omega = 0.8405912,
  ## country / diagnostic parameters
  theta = 1,
  ft = 0.4,
  rdt.det = 1,
  rdt.nonadherence = 0,
  microscopy.use = 0,
  include.nmf = FALSE,
  fever.age.breaks = c(-0.1,1,2,3,4,5,7,9,11,13,15,20,25,30,101),
  annual.nmf.per.age.bracket = c(2.235,1.841,1.503,1.175,0.928,0.67,0.486,0.512,0.475,0.429,0.652,0.747,0.804,0.789),
  nmf.multiplier=1,
  region = NULL,
  country = NULL,
  fitness = 1,
  ## entomological parameters
  mu.0 = 0.132,
  beta = c(0.05879106,0.07108894),
  a.k = 0.30677,
  ## delays and durations
  delay.mos = 10,
  delay.gam = 12.5,
  d.E = 12,
  d.I = 200,
  d.T = 5,
  d.D = 5,
  d.U = 110,
  d.P = 25,
  d.A = d.I - d.D,
  ## Pre-erythrocytic immunity parameters
  d.1 = 0.160527,
  d.ID = 10*365,
  I.D0 = 1.577533,
  k.D = 0.476614,
  u.D = 9.44512,
  a.D = 8001.99,
  f.D0 = 0.007055,
  g.D = 4.8183,
  ## Blood stage immunity parameters
  b.0 = 0.590076,
  b.1 = 0.5,
  d.B = 10*365,
  I.B0 = 43.8787,
  k.B = 2.15506,
  u.B = 7.19919,
  ## Acquired and maternal immunity parameters
  phi.0 = 0.791666,
  phi.1 = 0.000737,
  d.CA = 10950,
  I.C0 = 18.02366,
  k.C = 2.36949,
  u.C = 6.06349,
  P.M = 0.774368,
  d.CM = 67.6952,
  ## contributions to infectious reservoir by state and age
  g.1 = 1.82425,
  c.S = 0,
  c.D = 0.0676909,
  c.T = 0.322 * c.D,
  c.U = 0.0062,
  c.A = 0, # calculated internally
  c.P = 0)
{
  ##############################################################################################################################################################################
  ## EXTRA FUNCTIONS ##
  ##############################################################################################################################################################################

  ## function to createa annual season patterns
  seasonality <- function(country.values,region){

    # define vector of times spanning one year
    tvec = seq(0,1,l=365/time.step)

    # calculate Fourier series
    seasonality <- 1
    seasonality <- seasonality + country.values["seasonal_a1",region]*cos(1 * 2*pi*tvec) + country.values["seasonal_b1",region]*sin(1 * 2*pi*tvec)
    seasonality <- seasonality + country.values["seasonal_a2",region]*cos(2 * 2*pi*tvec) + country.values["seasonal_b2",region]*sin(2 * 2*pi*tvec)
    seasonality <- seasonality + country.values["seasonal_a3",region]*cos(3 * 2*pi*tvec) + country.values["seasonal_b3",region]*sin(3 * 2*pi*tvec)

    # ensure that scaling factor never goes below zero (this can happen in practice because we are only using the first few terms in an infinite series)
    seasonality[seasonality<0] <- 0

    return(seasonality)

  }

  ## function to find positions of infected individuals and randomly assign a strain type
  random_strains <- function(res) {
    infected <- which(res$Status[,res$Buffer[res$Counter]] %in% c(2,3,4,5))
    ## 2 in following exponential taken from mean MOI from Omer 2011 being roughly 2
    res$N.Sens[infected,res$Buffer[res$Counter]] <- sample(0:strains.0,size = length(infected),replace = TRUE,prob = rep(1/strains.0,strains.0+1))
    res$N.Dels[infected,res$Buffer[res$Counter]] <- sample(0:1,size = length(infected),replace = TRUE,prob = rep(1,2))
    check <- res$N.Dels[infected,res$Buffer[res$Counter]] + res$N.Sens[infected,1]
    res$N.Sens[infected[which(check==0)],res$Buffer[res$Counter]] <- 1
    return(res)
  }

  ## functon to determine new infection states of confirmed infections
  new_states <- function(prob.matrix){
    return (sample(x = c(2,5,3),size = 1,prob = prob.matrix))
  }

  ## function to update the result matrix in anticipation of next time step
  update <- function(res,i){

    ## initially assume new time-step has same states and strains as before

    ## find buffer position
    im1 <- which(res$Buffer==i-1)

    res$Status[,res$Counter] <- res$Status[,im1]
    res$N.Dels[,res$Counter] <- res$N.Dels[,im1]
    res$N.Sens[,res$Counter] <- res$N.Sens[,im1]

    ## Adjust immunity exponential decline
    res$I.B <- res$I.B*exp(-r.B)
    res$I.CA <- res$I.CA*exp(-r.CA)
    res$I.CM <- res$I.CM*exp(-r.CM)
    res$I.D <- res$I.D*exp(-r.ID)

    ## Calculate new individual contributions

    # individual infectivity to mosquitos contribution

    psi <- (1 - rho*exp(-res$Age/a.0))

    #psi <- psi / (mean(psi) / omega)
    psi <- psi/mean(psi)

    pi <- psi * zeta

    ## find buffer position
    img <- which(res$Buffer==i-grouped.delay)

    # calculate strains present in humans proportional to their individual bite rates and weighted MOI
    tots <- res$N.Dels[,img]+res$N.Sens[,img]
    tots[tots==0] <- NaN

    h.del.reservoir <- pi*res$N.Dels[,img]
    h.sen.reservoir <- pi*res$N.Sens[,img]

    # calculate contribution based on their infection status
    ## first non age dependent, i.e. non asymptomatic
    cont <- vector(mode = "numeric", length = N)
    non.A <- which(res$Status[,img] %in% c(2,4,5))
    cont[non.A] <- contribution[res$Status[non.A,img]]

    ## then asymptomatics
    A <- which(res$Status[,img] %in% 3)
    # calculate the f.D function
    f.D <- 1 - ((1 - f.D0) / (1 + ((res$Age[A] / a.D)^g.D)))

    # calculate q, which examines asymptomatics in terms of their parasetemia and thus contribution
    q <- d.1 + ( (1 - d.1) / (1 + f.D*((res$I.D[A] / I.D0)^k.D)))
    cont[A] <- c.U + (c.D - c.U)*(q^g.1)


    # multiply contribution by human reservoir
    h.del.reservoir <- (h.del.reservoir * cont)/tots
    h.sen.reservoir <- (h.sen.reservoir * cont)/tots


    # calculate human reservoir
    res$I.Reservoir[,res$Counter] <- c(sum(h.del.reservoir,na.rm=TRUE)/N,sum(h.sen.reservoir,na.rm=TRUE)/N)

    # mosquito dynamics
    surv.0 <- exp(-mu.0 * delay.mos)
    ince <- a.k*(sum(res$I.Reservoir[,res$Counter]))*res$Sv[im1]
    incv <- ince * surv.0

    res$FOIv[res$Counter] <- a.k*sum(res$I.Reservoir[,res$Counter])

    mv <- res$Sv[im1]+res$Ev[im1]+res$Iv[im1]
    beta <- mv * mu.0 * theta[i]

    res$Sv[res$Counter] <- max(0,res$Sv[im1] + (-ince - (mu.0*res$Sv[im1]) + beta))
    res$Ev[res$Counter] <- max(0,res$Ev[im1] + (ince - incv - (mu.0*theta[i]*res$Ev[im1])))
    res$Iv[res$Counter] <- max(0,res$Iv[im1] + (incv - (mu.0*theta[i]*res$Iv[im1])))

    return(res)
  }

  ## function to delete the last column from results as the last time step is not updated
  end_curation <- function(res){
    res$Time <- head(res$Time,-1)
    res$N.Dels <- res$N.Dels[,-res$Counter]
    res$N.Sens <- res$N.Sens[,-res$Counter]
    res$Status <- res$Status[,-res$Counter]
    res$Infections <- head(res$Infections, -1)
    res$Kidz.Infections <- head(res$Kidz.Infections, -1)
    res$dEIR <- head(res$dEIR, -1)
    res$Sv <- head(res$Sv, -1)
    res$Ev <- head(res$Ev, -1)
    res$Iv <- head(res$Iv, -1)
    res$FOIv <- head(res$FOIv,-1)
    res$Counter <- res$Counter - 1
    return(res)
  }

  # function to reload specified result
  fetch <- function(res,ID,direct.input){

    if(!is.null(direct.input)){
      reload.res <- direct.input
    } else if (length(unlist(strsplit(ID,"/")))>1){
      reload.res <- readRDS(ID)
    } else {
      reload.res <- context::task_result(id = ID,db = root)
    }

    res$Time <- c(reload.res$Time,res$Time+max(reload.res$Time)+1)
    res$Buffer <- reload.res$Buffer
    res$Status <- reload.res$Status
    res$Age <- reload.res$Age
    res$Zeta <- reload.res$Zeta
    ## Immunities
    res$Last_Infection_Time.CA <- reload.res$Last_Infection_Time.CA
    res$Last_Infection_Time.D <- reload.res$Last_Infection_Time.D
    res$Last_Bite_Time <- reload.res$Last_Bite_Time
    res$I.CA <- reload.res$I.CA
    res$I.CM <- reload.res$I.CM
    res$I.D <- reload.res$I.D
    res$I.B <- reload.res$I.B
    ## Mosquito / infection
    res$Sv <- reload.res$Sv
    res$Ev <- reload.res$Ev
    res$Iv <- reload.res$Iv
    res$I.Reservoir <- reload.res$I.Reservoir
    ## Prev / incidence outputs
    res$N.Sens <- reload.res$N.Sens
    res$N.Dels <- reload.res$N.Dels

    ## series
    res$S.Times <- reload.res$S.Times
    res$S.Status <- cbind(reload.res$S.Status,matrix(0,nrow=6,ncol=length(save.times)))
    res$S.dEIR <- c(reload.res$S.dEIR,rep(0,length(save.times)))
    res$S.I.B <- c(reload.res$S.I.B,rep(0,length(save.times)))
    res$S.I.CA <- c(reload.res$S.I.CA,rep(0,length(save.times)))
    res$S.I.CM <- c(reload.res$S.I.CM,rep(0,length(save.times)))
    res$S.I.D <- c(reload.res$S.I.D,rep(0,length(save.times)))
    res$S.Prev.All <- c(reload.res$S.Prev.All,rep(0,length(save.times)))
    res$S.Prev.05 <- c(reload.res$S.Prev.05,rep(0,length(save.times)))
    res$S.N.Sens <- c(reload.res$S.N.Sens,rep(0,length(save.times)))
    res$S.N.Dels <- c(reload.res$S.N.Dels,rep(0,length(save.times)))
    res$S.N.Sens.05 <- c(reload.res$S.N.Sens.05,rep(0,length(save.times)))
    res$S.N.Dels.05 <- c(reload.res$S.N.Dels.05,rep(0,length(save.times)))
    res$S.Incidence <- c(reload.res$S.Incidence,rep(0,length(save.times)))
    res$S.Incidence.05 <- c(reload.res$S.Incidence.05,rep(0,length(save.times)))
    res$S.Prev.Mono.D <- c(reload.res$S.Prev.Mono.D,rep(0,length(save.times)))
    res$S.Prev.Mono.S <- c(reload.res$S.Prev.Mono.S,rep(0,length(save.times)))
    res$S.Prev.Mono.D.05 <- c(reload.res$S.Prev.Mono.D.05,rep(0,length(save.times)))
    res$S.Prev.Mono.S.05 <- c(reload.res$S.Prev.Mono.S.05,rep(0,length(save.times)))
    res$S.Clin.Mono.D <- c(reload.res$S.Prev.Mono.D,rep(0,length(save.times)))
    res$S.Clin.Mono.S <- c(reload.res$S.Prev.Mono.S,rep(0,length(save.times)))
    res$S.Clin.Mono.D.05 <- c(reload.res$S.Clin.Mono.D.05,rep(0,length(save.times)))
    res$S.Clin.Mono.S.05 <- c(reload.res$S.Clin.Mono.S.05,rep(0,length(save.times)))

    return(res)
  }

  # function to change which strains infected individuals posses given a desired frequency of deleted strains
  freq.adjust <- function(res,target){

    former.hrp2 <- colSums(res$N.Sens)
    former.hrp2d <- colSums(res$N.Dels)
    former.ratio <- former.hrp2d/colSums(rbind(former.hrp2d,former.hrp2))

    hrp2d.tot <- colSums(res$N.Dels)
    hrp2.tot <- colSums(res$N.Sens)
    tot.tot <- colSums(rbind(hrp2d.tot,hrp2.tot))

    desired.d.tot <- round(tot.tot*target)
    n.change <- hrp2d.tot -  desired.d.tot

    if(sum(n.change) > 0){

      del.pos <- which(res$N.Dels!=0)
      dels <- res$N.Dels[del.pos]


      to.be.changed <- table(sample(rep(del.pos,dels),size = sum(n.change),replace=F))

      res$N.Dels[as.numeric(names(to.be.changed))] <- res$N.Dels[as.numeric(names(to.be.changed))] - to.be.changed
      res$N.Sens[as.numeric(names(to.be.changed))] <- res$N.Sens[as.numeric(names(to.be.changed))] + to.be.changed

    } else if (sum(n.change) < 0){

      sen.pos <- which(res$N.Sens!=0)
      sens <- res$N.Sens[sen.pos]


      to.be.changed <- table(sample(rep(sen.pos,sens),size = abs(sum(n.change)),replace=F))

      res$N.Dels[as.numeric(names(to.be.changed))] <- res$N.Dels[as.numeric(names(to.be.changed))] + to.be.changed
      res$N.Sens[as.numeric(names(to.be.changed))] <- res$N.Sens[as.numeric(names(to.be.changed))] - to.be.changed

    }

    hrp2 <- colSums(res$N.Sens)
    hrp2d <- colSums(res$N.Dels)
    ratio <- hrp2d/colSums(rbind(hrp2d,hrp2))

    hrp2.multipliers <- (hrp2/former.hrp2)
    hrp2d.multipliers <- (hrp2d/former.hrp2d)

    res$I.Reservoir <- t(cbind(res$I.Reservoir[1,]*hrp2d.multipliers,res$I.Reservoir[2,]*hrp2.multipliers))

    return(res)

  }


  ########################################################################################################################################################################################################################################
  ## HANDLE VARIABLES ##
  ########################################################################################################################################################################################################################################

  ## Handle arguments and create simulation variables
  #################################

  # variables to determine size preallocation
  max.time <- round(365 * years)
  its <- (max.time/time.step)+1

  ### res.cols <- (max.time/time.step)+1

  ## buffer set up
  res.cols <- storage
  buffer <- 1:res.cols
  ## round its up to the nearest res.cols so that result represents the last res.cols steps nicely
  its <- ceiling(its/res.cols)*res.cols
  max.time <- its*time.step
  logtenth <- floor(seq(2,its,its/10))


  # Sample of ages from given age distribution
  age <- sample(0:(365*max.age),size = N,replace = TRUE,prob = (365*max.age)*exp((0:(365*max.age))/-(d.death*365)))

  # lognormal individual biting heterogeneity
  zeta <- rlnorm(n = N,meanlog = -zeta.sd/2, sdlog = sqrt(zeta.sd))

  # beta given the roughly desired EIR
  beta <- beta[1]*(EIR*365) + beta[2]

  # group infectious state contributions for later
  contribution <- as.numeric(rbind(c.S,c.D,0,c.U,c.T,c.P))*time.step

  # seasonal effects
  theta <- rep(theta,length.out=max.time+1)

  # if alternative country parameterisation is set load from data frame
  if (!is.null(region)){
    load(paste("data/",country,".rda",sep=""))
    country.values <- get(country)
    ft <- country.values["fT",region]
    theta <- rep(seasonality(country.values,region),length.out=max.time+1)
    theta <- rep(1,length.out=max.time+1)
  }

  # individual ids
  ids <- 1:N


  # Rough equilibrium starting points relationship
  load("data/eq.0.rda")
  dat <- eq.0

  ## Create fixed variables
  #################################

  # total delay
  grouped.delay <- round((delay.gam)/time.step)
  rounded.dE <- round((d.E)/time.step)
  # per time step rate for treatment
  r.T = time.step /d.T
  # per time step rate for clinincal disease clearance
  r.D = time.step /d.D
  # per time step rate for becoming newly susceptible
  r.U = time.step /d.U
  # per time step rate for prophylaxis clearance
  r.P = time.step /d.P
  # per time step rate for sub-patent infection development
  r.A = time.step / d.A
  # per time step deccrease in blood stage immunity
  r.ID = time.step / d.ID
  # per time step decrease in pre-erythrocytic immunity
  r.B = time.step / d.B
  # per time step deccrease in acquired immunity
  r.CA = time.step / d.CA
  # per time step deccrease in maternal immunity
  r.CM = time.step / d.CM
  # per time step rate of death
  r.Death = time.step/(d.death*365)
  # per time step fever rate
  annual.nmf.per.age.bracket = (nmf.multiplier*annual.nmf.per.age.bracket) / (365/time.step)

  srs <- years*12
  save.times <- floor(seq(30,its,its/srs))


  ## Initiate Results Matrix ##
  ##########################################################
  res <- list("Buffer" = buffer,
              "Time" = seq(from=0,to=max.time,by=time.step),
              "Status" = matrix(nrow = N,ncol = res.cols,data = 0),
              "Age" = age,
              "Zeta" = zeta,
              "Counter" = 1,

              ## Immunity outputs

              "Last_Infection_Time.CA" = runif(N,0,1),
              "Last_Infection_Time.D" = runif(N,0,1),
              "Last_Bite_Time" = runif(N,0,1),
              "I.B" = vector(mode = "numeric", length = N),
              "I.CA" = vector(mode = "numeric", length = N),
              "I.CM" = vector(mode = "numeric", length = N),
              "I.D" = vector(mode = "numeric", length = N),

              ## Mosquito / infection outputs

              "I.Reservoir" = matrix(0,nrow = 2,ncol=res.cols),
              "Sv" = vector(mode = "numeric", length = res.cols),
              "Ev" = vector(mode = "numeric", length = res.cols),
              "Iv" = vector(mode = "numeric", length = res.cols),

              ## Prevalence / Incidence outputs

              "Prev.All" = vector(mode = "numeric", length = res.cols),
              "Prev.05" = vector(mode = "numeric", length = res.cols),
              "N.Sens" = matrix(nrow = N,ncol = res.cols,data = 0),
              "N.Dels" = matrix(nrow = N,ncol = res.cols,data = 0),
              "N.Sens.05" = vector(mode = "numeric", length = res.cols),
              "N.Dels.05" = vector(mode = "numeric", length = res.cols),
              "Incidence" = vector(mode = "numeric", length = res.cols),
              "Incidence.05" = vector(mode = "numeric", length = res.cols),
              "Prev.Mono.D" = vector(mode = "numeric", length = res.cols),
              "Prev.Mono.S" = vector(mode = "numeric", length = res.cols),
              "Prev.Mono.D.05" = vector(mode = "numeric", length = res.cols),
              "Prev.Mono.S.05" = vector(mode = "numeric", length = res.cols),
              "Clin.Mono.D" = vector(mode = "numeric", length = res.cols),
              "Clin.Mono.S" = vector(mode = "numeric", length = res.cols),
              "Clin.Mono.D.05" = vector(mode = "numeric", length = res.cols),
              "Clin.Mono.S.05" = vector(mode = "numeric", length = res.cols),
              "dEIR" = vector(mode = "numeric", length = res.cols),
              "FOIv" = vector(mode = "numeric", length = res.cols),

              ## Series Collection also

              "S.Times" = save.times,
              "S.Status" = matrix(nrow = 6,ncol = srs,data = 0),
              "S.I.B" = vector(mode = "numeric", length = srs),
              "S.I.CA" = vector(mode = "numeric", length = srs),
              "S.I.CM" = vector(mode = "numeric", length = srs),
              "S.I.D" = vector(mode = "numeric", length = srs),
              "S.Prev.All" = vector(mode = "numeric", length = srs),
              "S.Prev.05" = vector(mode = "numeric", length = srs),
              "S.N.Sens" = vector(mode = "numeric", length = srs),
              "S.N.Dels" = vector(mode = "numeric", length = srs),
              "S.N.Sens.05" = vector(mode = "numeric", length = srs),
              "S.N.Dels.05" = vector(mode = "numeric", length = srs),
              "S.Incidence" = vector(mode = "numeric", length = srs),
              "S.Incidence.05" = vector(mode = "numeric", length = srs),
              "S.Prev.Mono.D" = vector(mode = "numeric", length = srs),
              "S.Prev.Mono.S" = vector(mode = "numeric", length = srs),
              "S.Prev.Mono.D.05" = vector(mode = "numeric", length = srs),
              "S.Prev.Mono.S.05" = vector(mode = "numeric", length = srs),
              "S.Clin.Mono.D" = vector(mode = "numeric", length = srs),
              "S.Clin.Mono.S" = vector(mode = "numeric", length = srs),
              "S.Clin.Mono.D.05" = vector(mode = "numeric", length = srs),
              "S.Clin.Mono.S.05" = vector(mode = "numeric", length = srs),
              "S.dEIR" = vector(mode = "numeric", length = srs)

  )

  if(!is.null(ID)){

    # reload previous simulations
    res <- fetch(res,ID,direct.input)

    # adjust allele frequencies if required
    if(!is.null(desired.freq)){
      res <- freq.adjust(res,desired.freq)
    }
    # collect non stored varables and prepare results for new simulation
    zeta <- res$Zeta
    max.time <- max(res$Time)
    start <- res$Buffer[1]
    save.times <- floor(seq(start+28,max.time-1,(its-1)/srs))
    logtenth <- floor(seq(start+1,max.time-1,its/10))
    res$S.Times <- c(res$S.Times,save.times)
    its <- max.time-1
    theta <- rep(theta,length.out=its+1)
    true_grouped.delay <- grouped.delay
    true_rounded.dE <- rounded.dE
    ptm <- proc.time()

  } else {

    res$Sv[1] <- 0.4212434*(EIR*365) + 0.5143586
    res$Ev[1] <- 0.018553*(EIR*365)
    res$Iv[1] <- 0.006762714*(EIR*365)

    dat.x <- c(1,2,3,4,5,6,7,8,9,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90,95,100,105,
               110,115,120,125,130,135,140,145,150,155,160,165,170,175,180,185,190,195,200)

    eq <- which.min(abs(dat.x-(EIR*365)))

    ## initialise day 1 from previous equilibriums
    res$Status[,1]<- sample(states,N,replace = TRUE,
                            prob = c(dat$Eq.States[eq,])
    )

    ## randomly assign strains to those infected
    res <- random_strains(res)
    ##########################################
    ## INITIALISATION ##
    ####################

    # individual biting heterogeneity scaled by new ages
    psi <- (1 - rho*exp(-res$Age/a.0))
    psi <- psi/mean(psi)

    ##  immunities from previous data equilibriums
    res$I.D <- 1.4*zeta * dat$Eq.I.D[eq,1]*age + dat$Eq.I.D[eq,2]
    res$I.CA <- 1.3*zeta * dat$Eq.I.CA[eq,1]*age + dat$Eq.I.CA[eq,2]
    res$I.B <- 0.9*zeta * dat$Eq.I.B[eq,1]*age + dat$Eq.I.B[eq,2]

    mum.pos <- res$Age %in% seq(from=20*365, to=21*365, by=0.5) # "mums"
    res$I.CM <- (P.M * mean(res$I.CA[mum.pos])) * exp(-res$Age*r.CM) # maternal immunity

    # biting success rate dependent on age and pre-erythocytic immunity
    b <- b.0 * ( b.1 + ( (1 - b.1) / (1 + (res$I.B / I.B0)^k.B) ) )

    # symptom success rate dependent on age and pre-erythocytic immunity
    phi <- phi.0 * ( phi.1 + ( (1 - phi.1) / (1 + ((res$I.CA +res$I.CM) / I.C0)^k.C) ) )

    ## set true delay for initialisation hack
    true_grouped.delay <- grouped.delay
    true_rounded.dE <- rounded.dE
    grouped.delay <- 1
    rounded.dE <- 1


    ########################
    ## END INITIALISATION ##
    ##########################################

    ## prepare next step and start simulation from there
    ## buffering

    i <- 1
    res$Counter <- res$Counter + 1
    if (res$Counter > res.cols) res$Counter <- 1
    res$Buffer[res$Counter] <- i + 1

    res <- update(res,i+1)
    res$I.Reservoir[,1] <- res$I.Reservoir[,2]

    start <- res$Counter
    ## start timer
    ptm <- proc.time()

  }

  ########################################################################################################################################################################################################################################
  ## INITIATE LOOP ##
  ########################################################################################################################################################################################################################################

  for (i in start:(its)){

    t <- res$Time[i]

    ## initial hack to allow the first delay time steps to iterate
    if ((res$Counter - true_grouped.delay)>0){
      grouped.delay <- true_grouped.delay
      rounded.dE <- true_rounded.dE
    }

    ## prepare next time step ##
    ############################

    ## parameters that are needed outside the update loop ##

    # individual biting heterogeneity scaled by new ages
    psi <- (1 - rho*exp(-res$Age/a.0))

    #psi <- psi / (mean(psi) / omega)
    psi <- psi/mean(psi)

    # first calculate which individuals have been bitten
    pi <- psi * zeta

    ## find buffer position
    ime <- which(res$Buffer==i-rounded.dE)

    num.bites <- rbinom(n = N,1,prob = 1-exp(-(a.k * time.step * pi * res$Iv[ime])))
    pos.bites <- which(num.bites>0)

    res$dEIR[res$Counter] <- sum(num.bites[res$Age>(18*365)])/(sum(res$Age>(18*365))*time.step)
    test <- sum(num.bites[res$Age>(58*365)])/(sum(res$Age>(58*365))*time.step)
    ## logging
    if(is.element(i,logtenth)) {
      print(c(as.integer(i),res$Sv[res$Counter],res$Ev[res$Counter],res$Iv[res$Counter],res$dEIR[res$Counter]*365,test*365))
    }

    # update immunity due to bites
    bite.boostables.pos <- pos.bites[which(res$Last_Bite_Time[pos.bites] < t - u.B)]
    res$I.B[bite.boostables.pos] <- res$I.B[bite.boostables.pos] + 1

    # update time of last bite
    res$Last_Bite_Time[bite.boostables.pos] <- t - (res$Last_Bite_Time[bite.boostables.pos] -
                                                      floor(res$Last_Bite_Time[bite.boostables.pos] ))

    # biting success rate dependent on age and pre-erythocytic immunity
    b <- b.0 * ( b.1 + ( (1 - b.1) / (1 + ((res$I.B / I.B0)^k.B) ) ) )

    # find the positions of individuals who were bitten and can be infected
    I_pos <- pos.bites[res$Status[pos.bites,ime] %in% c(1,2,3,4)]

    # apply force of infection to all individuals who have been bitten and can be infected
    foni <-  (1-b[I_pos])^num.bites[I_pos]
    new.infections <- rbinom(n = length(foni),1,prob = foni)

    pos.new.infections <- I_pos[which(new.infections==0)]
    num.new.infections <- length(pos.new.infections)


    # catch if no infections
    if(num.new.infections==0) {
      no.infections <- TRUE
      pos.no.infections <- ids
    } else {
      no.infections <- FALSE
      pos.no.infections <- ids[-pos.new.infections]

      # update immunity due to inoculating bites
      I.D.boostables.pos <- pos.new.infections[which(res$Last_Infection_Time.D[pos.new.infections] < t - u.D)]
      I.CA.boostables.pos <- pos.new.infections[which(res$Last_Infection_Time.CA[pos.new.infections] < t - u.C)]
      res$I.D[I.D.boostables.pos] <- res$I.D[I.D.boostables.pos] + 1
      res$I.CA[I.CA.boostables.pos] <- res$I.CA[I.CA.boostables.pos] + 1

      # symptom success rate dependent on age and acquired immunity
      phi <- phi.0 * ( phi.1 + ( (1 - phi.1) / (1 + ( ( (res$I.CA +res$I.CM) / I.C0 )^k.C) ) ) )

      # update time of last infection
      res$Last_Infection_Time.CA[I.CA.boostables.pos] <- t - (res$Last_Infection_Time.CA[I.CA.boostables.pos] -
                                                                floor(res$Last_Infection_Time.CA[I.CA.boostables.pos]))
      res$Last_Infection_Time.D[I.D.boostables.pos] <- t - (res$Last_Infection_Time.D[I.D.boostables.pos] -
                                                              floor(res$Last_Infection_Time.D[I.D.boostables.pos]))
    }

    if(!no.infections){

      # U are subpatenet so when considered for reinfection they are essentially gaining a monoinfection in terms of when
      # we consider them for treatment by RDT
      pos.new.infections.U <- pos.new.infections[res$Status[pos.new.infections,res$Counter] == 4]
      pos.new.infections.non.U <- setdiff(pos.new.infections,pos.new.infections.U)

      # If there are any assumed fitness effects we assume that it causes an overall effect on the deletion strain contribution to
      # the infectious reservoir
      if(fitness != 1){

        fitness.effect <- res$I.Reservoir[,ime] * c(fitness,1)

        # draw which strains are passed on
        new.strains.U <- sample(x = c(FALSE,TRUE), size = length(pos.new.infections.U),replace = TRUE, prob = fitness.effect)
        new.strains.non.U <- sample(x = c(FALSE,TRUE), size = length(pos.new.infections.non.U),replace = TRUE, fitness.effect)
      } else {

        # draw which strains are passed on
        new.strains.U <- sample(x = c(FALSE,TRUE), size = length(pos.new.infections.U),replace = TRUE, prob = res$I.Reservoir[,ime])
        new.strains.non.U <- sample(x = c(FALSE,TRUE), size = length(pos.new.infections.non.U),replace = TRUE, prob = res$I.Reservoir[,ime])

      }
      # Add strain to those who are non U previously
      res$N.Sens[pos.new.infections.non.U[new.strains.non.U],res$Counter] <- res$N.Sens[pos.new.infections.non.U[new.strains.non.U],res$Counter] + 1
      res$N.Dels[pos.new.infections.non.U[!new.strains.non.U],res$Counter] <- res$N.Dels[pos.new.infections.non.U[!new.strains.non.U],res$Counter] + 1

      # Clear all strains from those who are U before assigning strain
      res$N.Sens[pos.new.infections.U,res$Counter] <- 0
      res$N.Dels[pos.new.infections.U,res$Counter] <- 0
      res$N.Sens[pos.new.infections.U[new.strains.U],res$Counter] <-  1
      res$N.Dels[pos.new.infections.U[!new.strains.U],res$Counter] <- 1

      total.strains <- sum(res$N.Sens[,res$Counter]) + sum(res$N.Dels[,res$Counter])

      # identify any hrp2' only individuals
      rdt.det.new.infections <- rep(1,num.new.infections)

      # probability that hrp2' only individual will still be treated, 1 - chance of not being hrp3, microscopied or nonadherence
      rdt.det.new.infections[which(res$N.Sens[pos.new.infections,res$Counter] == 0)] <- 1 - ((1-rdt.det)*(1-rdt.nonadherence)*(1-microscopy.use))

      # infection outcome

      # N -> D
      p.D <- phi[pos.new.infections] * (1 - (rdt.det.new.infections * ft))
      # N -> T
      p.T <- phi[pos.new.infections] * rdt.det.new.infections * ft
      # N -> A
      p.A <- 1 - phi[pos.new.infections]

      ## Need to find those who are to be infected again who are D as they should not then go to being
      ## asymptomatic trough the additonal infection. Going to T however seems plausible.
      pos.in.pos.new.infections.D <- which(res$Status[pos.new.infections,res$Counter] == 2)
      p.A[pos.in.pos.new.infections.D] <- 0

      # create infection outcome probability vector
      p.tot <- rbind(p.D,p.T,p.A)
      p.tot.col.cumsum <- apply(p.tot, 2, cumsum)
      normalised.p.tot.sum <- apply(p.tot.col.cumsum,2,function(x){return(x/x[3])})

      ## sample from cumulative distribution
      res$Status[pos.new.infections,res$Counter] <- c(2,5,3)[colSums(t(matrix(rep(runif(num.new.infections,0,1),3),ncol=3)) >(normalised.p.tot.sum))+1]



      if(sum(is.na(res$Status[,res$Counter]))>0){
        catach = 1
        fev <- catach * 5
        browser()
        im1 <- which(res$Buffer==i-1)
      }

      #####################################################################
      ## Record Prevalence / Incidence / Sensitive strain related values ##
      #####################################################################

      ## find buffer position
      im1 <- which(res$Buffer==i-1)

      ### Prevalence ########################################################################

      abs.prev <- sum(res$Status[,res$Counter] %in% c(2,3,4,5))
      res$Prev.All[res$Counter] <- abs.prev / N

      ## find location of under fives  and assign under five prevalence
      pos.05 <- which(res$Age<(365*5))
      abs.prev.05 <- sum(res$Status[pos.05,res$Counter] %in% c(2,3,4,5))
      res$Prev.05[res$Counter] <- abs.prev.05 / (length(pos.05))

      ### Incidence ########################################################################

      ## Incidence measurements, thus require locations of individuas who could be infected, i.e. S, D, A, U
      previous.y.pos <- pos.new.infections[res$Status[pos.new.infections,im1] %in% c(1,2,3,4)]

      ## Total New Infections
      res$Incidence[res$Counter] <- sum((res$Status[previous.y.pos,res$Counter] %in% c(2,5))) / N

      ## Total under 5s infection
      kids.y.pos <- previous.y.pos[res$Age[previous.y.pos] <= (5*365)]
      n.kids <- length(pos.05)
      res$Incidence.05[res$Counter] <- sum((res$Status[kids.y.pos,res$Counter] %in% c(2,5)))/n.kids

      ### Strain Values ########################################################################

      ## Total strains in under 5s
      total.strains.05 <- sum(res$N.Sens[pos.05,res$Counter]) + sum(res$N.Dels[pos.05,res$Counter])
      res$N.Sens.05[res$Counter] <- sum(res$N.Sens[pos.05,res$Counter]) / total.strains.05
      res$N.Dels.05[res$Counter] <- sum(res$N.Dels[pos.05,res$Counter]) / total.strains.05


      ## Monoinfection Prevalence ##
      pos.no.sens <- which(res$N.Sens[,res$Counter] == 0)
      res$Prev.Mono.D[res$Counter] <- sum((res$N.Dels[pos.no.sens,res$Counter] > 0)) / abs.prev
      pos.no.dels <- which(res$N.Dels[,res$Counter] == 0)
      res$Prev.Mono.S[res$Counter] <- sum((res$N.Sens[pos.no.dels,res$Counter] > 0)) / abs.prev

      ## Monoinfection Under 5s Prevalence ##
      pos.no.sens.05 <- pos.05[(res$N.Sens[pos.05,res$Counter] == 0)]
      res$Prev.Mono.D.05[res$Counter] <- sum((res$N.Dels[pos.no.sens.05,res$Counter] > 0))/abs.prev.05
      pos.no.dels.05 <- pos.05[(res$N.Dels[pos.05,res$Counter] == 0)]
      res$Prev.Mono.S.05[res$Counter] <- sum((res$N.Sens[pos.no.dels.05,res$Counter] > 0))/abs.prev.05

      ## Monoinfection Incidence ##
      res$Clin.Mono.D[res$Counter] <- sum((res$N.Dels[intersect(pos.no.sens,previous.y.pos)] > 0))/(N*365)
      res$Clin.Mono.S[res$Counter] <- sum((res$N.Sens[intersect(pos.no.dels,previous.y.pos)] > 0))/(N*365)

      ## Monoinfection Incidence Under 5s ##
      res$Clin.Mono.D.05[res$Counter] <- sum((res$N.Dels[intersect(pos.no.sens.05,previous.y.pos)] > 0))/n.kids
      res$Clin.Mono.S.05[res$Counter] <- sum((res$N.Sens[intersect(pos.no.dels.05,previous.y.pos)] > 0))/n.kids

      #####################################################################

    }

    ################################
    ## Non infection changes ##
    ################################

    # T -> P, i.e. 5 -> 6

    Tr.pos <- intersect(which(res$Status[,res$Counter] %in% 5),pos.no.infections)
    Tr.changes <- Tr.pos[rbinom(n = length(Tr.pos),1, prob = 1-exp(-r.T)) > 0]

    # update locations of prophylaxis
    res$Status[Tr.changes,res$Counter] <- 6
    res$N.Dels[Tr.changes,res$Counter] <- 0
    res$N.Sens[Tr.changes,res$Counter] <- 0

    # P -> S, i.e. 6 -> 1

    P.pos <- intersect(which(res$Status[,res$Counter] %in% 6),pos.no.infections)
    P.changes <- P.pos[rbinom(n = length(P.pos),1, prob = 1-exp(-r.P)) > 0]
    res$Status[P.changes,res$Counter] <- 1

    # D -> A, i.e. 2 -> 3

    D.pos <- intersect(which(res$Status[,res$Counter] %in% 2),pos.no.infections)
    D.changes <- D.pos[rbinom(n = length(D.pos),1, prob = 1-exp(-r.D)) > 0]
    res$Status[D.changes,res$Counter] <- 3

    # A -> U, i.e. 3 -> 4

    A.pos <- intersect(which(res$Status[,res$Counter] %in% 3),pos.no.infections)
    A.changes <- A.pos[rbinom(n = length(A.pos),1, prob = 1-exp(-r.A)) > 0]
    res$Status[A.changes,res$Counter] <- 4

    # U -> S, i.e. 4 -> 1

    U.pos <- intersect(which(res$Status[,res$Counter] %in% 4),pos.no.infections)
    U.changes <- U.pos[rbinom(n = length(U.pos),1, prob = 1-exp(-r.U)) > 0]

    # update locations of susceptibles

    res$Status[U.changes,res$Counter] <- 1
    res$N.Dels[U.changes,res$Counter] <- 0
    res$N.Sens[U.changes,res$Counter] <- 0

    ################################
    ## Non compartmental changes ##
    ################################

    ################################
    ## NMF changes ##
    ################################

    if(include.nmf==TRUE){

      band_ages <- cut(res$Age,fever.age.breaks*365)
      nmf.pos <- which(sapply(band_ages,function(x){return(rbinom(1,1,1-exp(-annual.nmf.per.age.bracket[x])))})==1)

      if(length(nmf.pos)>0){

        # first find those in S as the only chance that they will be treated upon seeking it for the nmf is due to nonadherence and we move them to prophylaxis
        # as they don't have any strains and should thus not be in treated
        S.pos.nmf <- nmf.pos[(res$Status[nmf.pos,res$Counter] %in% c(1))]
        S.treated.pos.nmf <- S.pos.nmf[which(rbinom(length(S.pos.nmf),size = 1,rdt.nonadherence*ft)==1)]
        res$Status[S.treated.pos.nmf,res$Counter] <- 6

        # next find those in U as the only chance that they will be treated upon seeking it for the nmf is due to nonadherence and we move them to treated
        U.pos.nmf <- nmf.pos[(res$Status[nmf.pos,res$Counter] %in% c(4))]
        U.treated.pos.nmf <- U.pos.nmf[which(rbinom(length(U.pos.nmf),size = 1,rdt.nonadherence*ft)==1)]
        res$Status[U.treated.pos.nmf,res$Counter] <- 5

        # next find those in D or A as they will move to treated unless they are hrp2 etc
        D.or.A.pos.nmf <- nmf.pos[(res$Status[nmf.pos,res$Counter] %in% c(2,3))]

        # identify any hrp2' only individuals
        rdt.det.D.or.A.nmf <- rep(1,length(D.or.A.pos.nmf))
        # probability that hrp2' only individual will still be treated, 1 - chance of not being hrp3, microscopied or nonadherence
        rdt.det.D.or.A.nmf[which(res$N.Sens[D.or.A.pos.nmf,res$Counter] == 0)] <- ft * (1 - ((1-rdt.det)*(1-rdt.nonadherence)*(1-microscopy.use)))

        D.or.A.treated.pos.nmf <- D.or.A.pos.nmf[which(rbinom(length(D.or.A.pos.nmf),size = 1,rdt.det.D.or.A.nmf)==1)]
        res$Status[D.or.A.treated.pos.nmf,res$Counter] <- 5

      }

    }

    ## loss of individual strains
    #########################

    Strain_pos <- which(res$Status[,res$Counter] %in% c(2,3,4,5))

    # calculate individual probabilities of losing a strain for those with strains
    probs <- 1-(exp(-time.step *(res$N.Dels[Strain_pos,res$Counter] + res$N.Sens[Strain_pos,res$Counter]) / (d.A+d.U)))
    # Identify individuals who have been selected for strain loss
    queue_drop_pos <- Strain_pos[which(rbinom(length(probs),1,probs)==1)]
    # calculate ratio of selective strains to deletion strains for those selected
    p <- res$N.Sens[queue_drop_pos,res$Counter] / (res$N.Sens[queue_drop_pos,res$Counter]+res$N.Dels[queue_drop_pos,res$Counter])
    # Use ratios to sample strains to be lost
    sens.drop.pos <- which(rbinom(length(p),1,p)==1)
    # Apply loss of selected strain type
    res$N.Sens[queue_drop_pos[sens.drop.pos],res$Counter] <- res$N.Sens[queue_drop_pos[sens.drop.pos],res$Counter] - 1
    if(length(sens.drop.pos)==0){
      res$N.Dels[queue_drop_pos,res$Counter] <- res$N.Dels[queue_drop_pos,res$Counter] - 1
    } else {
      res$N.Dels[queue_drop_pos[-sens.drop.pos],res$Counter] <- res$N.Dels[queue_drop_pos[-sens.drop.pos],res$Counter] - 1
    }


    ## death/age update
    #########################

    # first update everyone's age
    res$Age <- res$Age + time.step

    # identify those who die
    death.pos <- unique(c(which(res$Age>max.age*365),which(rbinom(n=N,size=1,prob=1-exp(-r.Death))==1)))

    # update those individuals to appear as susceptible newborns and assign immunities
    res$Status[death.pos,res$Counter] <- 1 # susceptible status
    res$N.Sens[death.pos,res$Counter] <- 0 # no strains
    res$N.Dels[death.pos,res$Counter] <- 0
    res$I.D[death.pos] <- 0 # no blood stage immunity
    res$I.CA[death.pos] <- 0 # no acquired immunity
    mum.pos <- res$Age %in% seq(from=20*365, to=21*365, by=0.5) # "mums"
    res$I.CM[death.pos] <- P.M * mean(res$I.CA[mum.pos]) # maternal immunity
    res$I.B[death.pos] <- 0 # no preerythrocytic immunity
    # assume time to be current time and no boosting period is assumed within maternal immunity
    res$Last_Infection_Time.CA[death.pos] <- t - runif(n = length(death.pos),min = 0,max = 1)
    res$Last_Infection_Time.D[death.pos] <- t - runif(n = length(death.pos),min = 0,max = 1)
    res$Last_Bite_Time[death.pos] <- t - runif(n = length(death.pos),min = 0,max = 1)
    res$Age[death.pos] <- 0
    res$zeta[death.pos] <- rlnorm(n = length(death.pos),meanlog = -zeta.sd/2, sdlog = sqrt(zeta.sd))
    zeta[death.pos] <- res$zeta[death.pos]


    #######################################
    ## END DEMOGRAPHIC/GENETIC CHANGES ##
    #######################################

    # store series

    if(is.element(t,res$S.Times)) {

      ## use mean over last week for storage to help with capturing low incidence settings
      last.7 <- match((i-6) : i,res$Buffer)
      # which series column are the results going into
      pos <- match(t,res$S.Times)

      ### immunity storage ###

      res$S.I.B[pos] <- mean(res$I.B)
      res$S.I.CA[pos] <- mean(res$I.CA)
      res$S.I.CM[pos] <- mean(res$I.CM)
      res$S.I.D[pos] <- mean(res$I.D)

      ### status and EIR storage ###
      res$S.Status[,pos] <- table(factor(res$Status[,last.7],levels=1:6)) / (7*N)
      res$S.dEIR[pos] <- mean(res$dEIR[last.7])

      ### Prevalence / Incidence Storage ###
      res$S.Prev.All[pos] <- mean(res$Prev.All[last.7])
      res$S.Prev.05[pos] <- mean(res$Prev.05[last.7])
      res$S.Incidence[pos] <- mean(res$Incidence[last.7])
      res$S.Incidence.05[pos] <- mean(res$Incidence.05[last.7])

      ### Strain Storage ###
      total.last.7.strains <- sum(res$N.Sens[,last.7]) + sum(res$N.Dels[,last.7])
      res$S.N.Sens[pos] <- sum(res$N.Sens[,last.7]) / total.last.7.strains
      res$S.N.Dels[pos] <- sum(res$N.Dels[,last.7]) / total.last.7.strains
      res$S.N.Sens.05[pos] <- mean(res$N.Sens.05[last.7])
      res$S.N.Dels.05[pos] <- mean(res$N.Dels.05[last.7])

      res$S.Prev.Mono.D[pos] <- mean(res$Prev.Mono.D[last.7])
      res$S.Prev.Mono.S[pos] <- mean(res$Prev.Mono.S[last.7])
      res$S.Prev.Mono.D.05[pos] <- mean(res$Prev.Mono.D.05[last.7])
      res$S.Prev.Mono.S.05[pos] <- mean(res$Prev.Mono.S.05[last.7])
      res$S.Clin.Mono.D[pos] <- mean(res$Clin.Mono.D[last.7])
      res$S.Clin.Mono.S[pos] <- mean(res$Clin.Mono.S[last.7])
      res$S.Clin.Mono.D.05[pos] <- mean(res$Clin.Mono.D.05[last.7])
      res$S.Clin.Mono.S.05[pos] <- mean(res$Clin.Mono.S.05[last.7])

    }

    ## buffering
    res$Counter <- res$Counter + 1
    if (res$Counter > res.cols) res$Counter <- 1
    res$Buffer[res$Counter] <- i + 1



    # update results
    res <- update(res,i + 1)

  }

  ## print times
  print(paste(N,"individuals analysed, over a period of ",max.time, " days at ",time.step,
              " day intervals within ", ((proc.time() - ptm)[3]),"secs"))

  res$Status <- matrix(as.integer(res$Status),nrow = N)
  res$N.Sens <- matrix(as.integer(res$N.Sens),nrow = N)
  res$N.Dels <- matrix(as.integer(res$N.Dels),nrow = N)



  return(res)
}
