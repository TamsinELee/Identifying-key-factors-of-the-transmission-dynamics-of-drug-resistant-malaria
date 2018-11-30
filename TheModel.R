# Identifying key factors of the transmission dynamics of drug-resistant malaria
# Lee and Penny 2018 Journal of Theoretical Biology
# The model: Equations (1) to (8) 

dev.off()
rm(list = ls()) 

library(ggplot2)
library(reshape2)

# Time step is one day
NYr          <- 3
NDays        <- 365 * NYr
BurnIn       <- 0
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# "Known" parameters (from Mandal et al. Malaria Journal 2011)
mvec         <- 20.25+(39.5/2)*sin(seq(0,2*pi*NYr,2*pi*NYr/NDays)) # Number of female mosquitoes to humans (seasonal)
mvec         <- mvec/mvec * 28                          # Switched off seasonality
a            <- 0.25 #runif(1, 0.01, 0.5)               # Biting rate (changed after review)
b0           <- 0.3 #runif(1, 0.2, 0.5)                 # Probability of transmitting from mosq to human
b1           <- 0.3#0.28#0.3 #runif(1, 0.2, 0.5)        # Probability of transmitting from mosq to human
b2           <- 0.3#0.2 #runif(1, 0.2, 0.5)             # Probability of transmitting from mosq to human
c0           <- 0.5                                     # Probability of transmitting from human to mosq
c1           <- 0.5#0.4                                 # Probability of transmitting from human to mosq
c2           <- 0.5#0.3                                 # Probability of transmitting from human to mosq

mu           <- 0.017/365                               # Background death rate of humans
muhat        <- 0.2 #runif(1, 0.05, 0.5)                # Background death rate of mosq
tau          <- 11 #floor(runif(1, 5, 15))              # Latency period of mosquitoes


# Model parameters
alphaI       <- mu * 2                                  # Death rate of infecteds
alphaT       <- mu                                      # Death rate of treated (mu means perfect treatment)
rI           <- 0.02 #0.03 #runif(1, 0.005, 0.05)       # Recovery rate of infecteds
rT0          <- rI * 3                                  # Recovery rate of treated, sensitive infection
rT1          <- rI * 2                                  # Recovery rate of treated, partially resistant infections
rx           <- 0.03                                    # Treatment rate per day (0.036 is 40% chance in 2 week period)

phi          <- 1/110                                   # Replacement rate


#-------------------------------------------------------------------------
# Intial total population sizes
#-------------------------------------------------------------------------

HH           <- matrix(nrow = NDays, ncol = 1)
H            <- matrix(nrow = NDays, ncol = 1)

Nhuman       <- 100
HH[1:tau]    <- Nhuman                                 # Number of humans
H[1:tau]     <- Nhuman                                 # Number of mosquitoes

#-------------------------------------------------------------------------
# Initial population spread and setting up storage.
# Double letters is human, single is mosquito. 
#-------------------------------------------------------------------------

# Storage for number of infected and treated humans
II0          <- matrix(nrow = NDays, ncol = 1)          # Sensitive
II1          <- matrix(nrow = NDays, ncol = 1)          # Partially resistant
II2          <- matrix(nrow = NDays, ncol = 1)          # Resistant
TT0          <- matrix(nrow = NDays, ncol = 1)          # Sensitive
TT1          <- matrix(nrow = NDays, ncol = 1)          # Partially resistant 
TT2          <- matrix(nrow = NDays, ncol = 1)          # Resistant

# Storage of number of exposed and infected mosquitoes
E0           <- matrix(nrow = NDays, ncol = 1)          # Sensitive
E1           <- matrix(nrow = NDays, ncol = 1)          # Partially resistant
E2           <- matrix(nrow = NDays, ncol = 1)          # Resistant
I0           <- matrix(nrow = NDays, ncol = 1)          # Sensitive
I1           <- matrix(nrow = NDays, ncol = 1)          # Partially resistant
I2           <- matrix(nrow = NDays, ncol = 1)          # Resistant

# Storage of rate of change for number of infected and treated humans
dII0         <- matrix(nrow = NDays, ncol = 1)          # Sensitive
dII1         <- matrix(nrow = NDays, ncol = 1)          # Partially resistant
dII2         <- matrix(nrow = NDays, ncol = 1)          # Resistant
dTT0         <- matrix(nrow = NDays, ncol = 1)          # Sensitive
dTT1         <- matrix(nrow = NDays, ncol = 1)          # Partially resistant
dTT2         <- matrix(nrow = NDays, ncol = 1)          # Resistant

II0          <- 0.005 * HH                              # 
II1          <- 0 * HH                                  # No initial partially resistant infections
II2          <- 0 * HH                                  # No initial resistant infections

TT0          <- 0.001 * HH * 0.99   # proportion of infected with Sensitive
TT1          <- 0.001 * HH * 0.008  # proportion of infected with Partially resistant - initally v.small
TT2          <- 0.001 * HH * 0.002  # proportion of infected with Resistant - initally v.small

SS           <- HH - II0 - II1 - II2 - TT0 - TT1 - TT2 # susceptible human

E0           <- 0.1 * H      # exposed mosq.Res 0
E1           <- 0 * H        # exposed mosq.Res 1
E2           <- 0 * H        # exposed mosq.Res 2

I0           <- 0.1 * H      # infected mosq. with Res0 
I1           <- 0 * H        # infected mosq. with Res1 
I2           <- 0 * H        # infected mosq. with Res1 

S            <- H - I0 - I1 - I2 - E0 - E1 - E2   # susceptible mosq.

phiAdded     <- matrix(NA, nrow = NDays, ncol = 1)     

#-------------------------------------------------------------------------
# Time stepping
#-------------------------------------------------------------------------
  for (day in (tau+1):NDays){
    
    #if (day > 365){phi <- 0}  # switch off replacement rate - like doing a cycle of drugs after a while
    
    m        <- mvec[day] # m changes from low to high (to replicate seasonality)
    
    #-------------------------------------------------------------------------
    # The rates of change in mosquitoes
    #-------------------------------------------------------------------------
    
    dE0dt    <- a * c0 * (II0[day-1] + TT0[day-1])/HH[day-1] * S[day-1] - # susceptible become exposed
      a * c0 * (II0[day-tau] + TT0[day-tau])/HH[day-tau] * S[day-tau] * exp(-muhat * tau) - # exposed become infected
      muhat * E0[day-1]  # background death
    
    dE1dt    <- a * c1 * (II1[day-1] + TT1[day-1])/HH[day-1] * S[day-1] - # susceptible become exposed
      a * c1 * (II1[day-tau] + TT1[day-tau])/HH[day-tau] * S[day-tau] * exp(-muhat * tau) - # exposed become infected
      muhat * E1[day-1]  # background death
    
    dE2dt    <- a * c2 * (II2[day-1] + TT2[day-1])/HH[day-1] * S[day-1] - # susceptible become exposed
      a * c2 * (II2[day-tau] + TT2[day-tau])/HH[day-tau] * S[day-tau] * exp(-muhat * tau) - # exposed become infected
      muhat * E2[day-1]  # background death
    
    dI0dt    <- a * c0 * (II0[day-tau] + TT0[day-tau])/HH[day-tau] * S[day-tau] * exp(-muhat * tau) - # exposed become infected
      muhat * I0[day-1]  # background death
    
    dI1dt    <- a * c1 * (II1[day-tau] + TT1[day-tau])/HH[day-tau] * S[day-tau] * exp(-muhat * tau) - # exposed become infected
      muhat * I1[day-1]  # background death
    
    dI2dt    <- a * c2 * (II2[day-tau] + TT2[day-tau])/HH[day-tau] * S[day-tau] * exp(-muhat * tau) - # exposed become infected
      muhat * I2[day-1]  # background death
    
    #-------------------------------------------------------------------------
    # The rates of change in humans
    #-------------------------------------------------------------------------
    
    # To keep the population constant, set lambda to lambdaConst
    lambdaConst     <- alphaI * (II0[day-1] + II1[day-1] + II2[day-1] + TT2[day-1]) + 
                       alphaT * (TT0[day-1] + TT1[day-1]) 
    lambda          <- lambdaConst
    
    totalRecover    <- rI * (II0[day-1] + II1[day-1] + II2[day-1] + TT2[day-1]) + 
                      rT0 * TT0[day-1] +  rT1 * TT1[day-1]
    
    dSSdt           <- lambda + totalRecover - 
      (a * b0 * m * SS[day-1]/HH[day-1] * I0[day-1]) -
      (a * b1 * m * SS[day-1]/HH[day-1] * I1[day-1]) - 
      (a * b2 * m * SS[day-1]/HH[day-1] * I2[day-1])
    
    dII0dt          <- a * b0 * m * SS[day-1]/HH[day-1] * I0[day-1] - # become infected 
      rx * II0[day-1] - # treated
      rI * II0[day-1] - # recover
      alphaI * II0[day-1] # virulence
    
    dII1dt <- a * b1 * m * SS[day-1]/HH[day-1] * I1[day-1] - # become infected       
      rx * II1[day-1] - # treated
      rI * II1[day-1] - # recover
      alphaI * II1[day-1] # virulence
    
    dII2dt          <- a * b2 * m * SS[day-1]/HH[day-1] * I2[day-1] - # become infected 
      rx * II2[day-1] - # treated
      rI * II2[day-1] - # recover
      alphaI * II2[day-1] # virulence
    
    dTT0dt          <- rx * II0[day-1] -
      rT0 * TT0[day-1] - # recover
      alphaT * TT0[day-1]  # virulence
    
    dTT1dt          <-  rx * II1[day-1] -
      rT1 * TT1[day-1] - # recover
      alphaT * TT1[day-1]  # virulence
    
    dTT2dt          <- rx * II2[day-1] -
      rI * TT2[day-1] - # recover - note same rI as for infected group
      alphaI * TT2[day-1] # virulence - note same alpha as for infected group
    
    #if (day > (BurnIn * 365)){
    #if (dTT2dt < 0.00001){
      
      phiAdded[day] <- day
      
      dTT0dt        <- dTT0dt - phi * TT0[day-1] # replaced by Partially resistant
      
      dTT1dt        <- dTT1dt + phi * TT0[day-1] - # from Sensitive
        phi * TT1[day-1] # replaced by Resistant
      
      dTT2dt        <- dTT2dt + phi * TT1[day-1] # from Resistant

    #} # end if loop 
  
    #-------------------------------------------------------------------------
    # Updating 
    #-------------------------------------------------------------------------
    
    H[day]        <- H[day-1]  # constant population 
    
    E0[day]       <- E0[day-1] + dE0dt
    E1[day]       <- E1[day-1] + dE1dt
    E2[day]       <- E2[day-1] + dE2dt
    I0[day]       <- I0[day-1] + dI0dt
    I1[day]       <- I1[day-1] + dI1dt
    I2[day]       <- I1[day-1] + dI2dt
    S             <- H - E0 - E1 - E2 - I0 - I1 - I2   # susceptible mosq - this eqn only works for constant popln
    
    II0[day]      <- II0[day-1] + dII0dt
    II1[day]      <- II1[day-1] + dII1dt
    II2[day]      <- II2[day-1] + dII2dt
    
    TT0[day]      <- TT0[day-1] + dTT0dt
    TT1[day]      <- TT1[day-1] + dTT1dt
    TT2[day]      <- TT2[day-1] + dTT2dt
    
    SS[day]       <- SS[day-1] + dSSdt
    HH[day]       <- SS[day] + II0[day] + II1[day] + II2[day] + TT0[day] + TT1[day] + TT2[day] # constant population 

    dII0[day]     <- dII0dt
    dII1[day]     <- dII1dt
    dII2[day]     <- dII2dt
    dTT0[day]     <- dII0dt
    dTT1[day]     <- dTT1dt
    dTT2[day]     <- dTT2dt
    
} # End day loop


#-------------------------------------------------------------------------
# Mosquitoes store and plot
#-------------------------------------------------------------------------

MosqDF            <- data.frame(matrix(ncol=8,nrow=NDays))
colnames(MosqDF)  <- c("Day","S","E0","E1","E2","I0","I1","I2")
MosqDF$Day        <- seq(1,NDays,1)
MosqDF$S          <- as.numeric(S)
MosqDF$E0         <- as.numeric(E0)
MosqDF$E1         <- as.numeric(E1)
MosqDF$E2         <- as.numeric(E2)
MosqDF$I0         <- as.numeric(I0)
MosqDF$I1         <- as.numeric(I1)
MosqDF$I2         <- as.numeric(I2)

#-------------------------------------------------------------------------

MosqDFmelt        <- melt(MosqDF, id=c("Day"))
pm                <- ggplot(MosqDFmelt, aes(x=Day, y=value, fill=variable)) + geom_area()
pm  + scale_fill_manual(values=alpha(c("grey83", "tan1","tan3", "tan4",
                                "violet", "violetred3","violetred4"), 0.64)) +
 theme(axis.text.x = element_text(size=14),
       axis.text.y = element_text(size=14)) +
 scale_x_continuous(name="Years", limits=c((BurnIn * 365 + 1), NDays), breaks=seq((BurnIn * 365 + 1),NDays,365), labels=seq(BurnIn + 1,NYr,1)) +
 scale_y_continuous(name="Percentage of mosquitoes") +
 guides(fill=guide_legend(title="")) +
 theme(legend.text=element_text(size=16)) +
 theme(axis.title=element_text(size=18))

#-------------------------------------------------------------------------
# Human store and plot
#-------------------------------------------------------------------------

HumanDF           <- data.frame(matrix(ncol=8,nrow=NDays))
colnames(HumanDF) <- c("Day","S","I0","I1","I2","T0","T1","T2")
HumanDF$Day       <- seq(1,NDays,1)
HumanDF$S         <- as.numeric(SS)
HumanDF$I0        <- as.numeric(II0)
HumanDF$I1        <- as.numeric(II1)
HumanDF$I2        <- as.numeric(II2)
HumanDF$T0        <- as.numeric(TT0)
HumanDF$T1        <- as.numeric(TT1)
HumanDF$T2        <- as.numeric(TT2)

#-------------------------------------------------------------------------

HumanDFmelt       <- melt(HumanDF, id=c("Day")) 
ph                <- ggplot(HumanDFmelt, aes(x=Day, y=value, fill=variable)) + geom_area()
ph  + scale_fill_manual(values = alpha(c("grey83", "violet", "violetred1","violetred3", 
                                 "cadetblue1", "cadetblue3", "violetred4"), 0.78)) + 
  theme_bw() +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  scale_x_continuous(name="Years", limits=c((BurnIn * 365), NDays), breaks=seq((BurnIn * 365),NDays,365), labels=seq(BurnIn,NYr,1)) + 
  scale_y_continuous(name="Percentage of hosts") + 
  guides(fill=guide_legend(title="")) + 
  guides(fill=FALSE) + 
  theme(legend.text=element_text(size=30)) + 
  theme(axis.title=element_text(size=30))

#-------------------------------------------------------------------------
# Treated population store and plot
#-------------------------------------------------------------------------
 
ResDF             <- data.frame(matrix(ncol=4,nrow=NDays))
colnames(ResDF)   <- c("Day","T0","T1","T2")
ResDF$Day         <- seq(1,NDays,1)
ResDF$T0          <- as.numeric(TT0)
ResDF$T1          <- as.numeric(TT1)
ResDF$T2          <- as.numeric(TT2)
 
ResDFmelt         <- melt(ResDF, id=c("Day")) 
pt                <- ggplot(ResDFmelt, aes(x=Day, y=value, fill=variable)) + geom_area(alpha=1)
pt  + scale_fill_manual(values=c("cadetblue1", "cadetblue3", "cadetblue4")) + 
 theme(axis.text.x = element_text(size=14),
       axis.text.y = element_text(size=14)) + 
 scale_x_continuous(name="Years", limits=c((BurnIn * 365 + 1), NDays), breaks=seq((BurnIn * 365 + 1),NDays,365), labels=seq(BurnIn + 1,NYr,1)) + 
 scale_y_continuous(name="Treated hosts only") + 
 guides(fill=guide_legend(title="")) + 
 theme(legend.text=element_text(size=16)) + 
 theme(axis.title=element_text(size=18))

#=====================================================================================================

# Reproductive numbers

#===================================================================================================== 

OutsideBrackets    <- a^2 * m / muhat * exp(-muhat * tau)

InsideBrackets00   <- b0 * c0 * (rT0 + rx + phi + mu)/((rI + rx + alphaI) * (rT0 + phi + mu))
InsideBrackets11   <- b1 * c1 * (rT1 + rx + phi + mu)/((rI + rx + alphaI) * (rT1 + phi + mu))
InsideBrackets22   <- b2 * c2 * 1/(rI + alphaI)
InsideBrackets01   <- b1 * c0 * (rx * phi)/((rI + rx + alphaI) * (rT0 + phi + mu) * (rT1 + phi + mu))
InsideBrackets12   <- b2 * c1 * (rx * phi)/((rI + rx + alphaI) * (rT1 + phi + mu) * (rI + alphaI))
InsideBrackets02   <- b2 * c0 * (rx * phi * phi)/((rI + rx + alphaI) * (rT0 + phi + mu) * (rT1 + phi + mu) * (rI + alphaI))

R00                <- OutsideBrackets * InsideBrackets00
R11                <- OutsideBrackets * InsideBrackets11
R22                <- OutsideBrackets * InsideBrackets22
R01                <- OutsideBrackets * InsideBrackets01
R12                <- OutsideBrackets * InsideBrackets12
R02                <- OutsideBrackets * InsideBrackets02

RR                 <- R00 + R11 + R22 + R01 + R02 + R12
paste("Total reproductive number:", RR)