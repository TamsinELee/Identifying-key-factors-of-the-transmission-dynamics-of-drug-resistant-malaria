
dev.off()
rm(list = ls()) 

library(ggplot2)
library(reshape2)

a                   <- 0.25
b                   <- 0.2
c                   <- 0.5
m                   <- 28
muhat               <- 0.2
rI                  <- 0.2
alpha               <- 0.0004657534
rT0                 <- 3 * rI
rT1                 <- 2 * rI
mu                  <- 0.017/365
lambda              <- 0.017/365 # assume lambda = mu for constant pop.in disease free equilibrium
tauhat              <- 11

rx                  <- 0.5
phi                 <- 0.5

Npts                <- 1001

R012012             <- matrix(nrow = Npts, ncol = 1)
R00                 <- matrix(nrow = Npts, ncol = 1)
R11                 <- matrix(nrow = Npts, ncol = 1)
R22                 <- matrix(nrow = Npts, ncol = 1)
R01                 <- matrix(nrow = Npts, ncol = 1)
R12                 <- matrix(nrow = Npts, ncol = 1)
R02                 <- matrix(nrow = Npts, ncol = 1)


# ====================================================================================
# rT0
# ====================================================================================

vec                  <- seq(1,10,length.out=Npts)

for (k1 in 1:Npts){
  
  rT0                <- vec[k1] * rI

  OutsideBrackets    <- a^2 * b * c * m / muhat * exp(-muhat * tauhat)
  
  InsideBrackets00   <- (rT0 + rx + phi + mu)/((rI + rx + alpha) * (rT0 + phi + mu))
  InsideBrackets11   <- (rT1 + rx + phi + mu)/((rI + rx + alpha) * (rT1 + phi + mu))
  InsideBrackets22   <- 1/(rI + alpha)
  InsideBrackets01   <- (rx * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu))
  InsideBrackets12   <- (rx * phi)/((rI + rx + alpha) * (rT1 + phi + mu) * (rI + alpha))
  InsideBrackets02   <- (rx * phi * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu) * (rI + alpha))
  
  R00[k1]            <- OutsideBrackets * InsideBrackets00
  R11[k1]            <- OutsideBrackets * InsideBrackets11
  R22[k1]            <- OutsideBrackets * InsideBrackets22
  R01[k1]            <- OutsideBrackets * InsideBrackets01
  R12[k1]            <- OutsideBrackets * InsideBrackets12
  R02[k1]            <- OutsideBrackets * InsideBrackets02
  
}

R012012              <- R00 + R11 + R22 + R01 + R12 + R02

store1               <- cbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

storea               <- rbind(R012012/R012012[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

storea1              <- rbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1])

storea2              <- rbind(R01/R01[1], R12/R12[1], R02/R02[1])

rT0                  <- 3 * rI # set back to default

# ====================================================================================
# rT1
# ====================================================================================

vec                  <- seq(1,10,length.out=Npts)

for (k1 in 1:Npts){
  
  rT1                <- vec[k1] * rI
  
  OutsideBrackets    <- a^2 * b * c * m / muhat * exp(-muhat * tauhat)
  
  InsideBrackets00   <- (rT0 + rx + phi + mu)/((rI + rx + alpha) * (rT0 + phi + mu))
  InsideBrackets11   <- (rT1 + rx + phi + mu)/((rI + rx + alpha) * (rT1 + phi + mu))
  InsideBrackets22   <- 1/(rI + alpha)
  InsideBrackets01   <- (rx * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu))
  InsideBrackets12   <- (rx * phi)/((rI + rx + alpha) * (rT1 + phi + mu) * (rI + alpha))
  InsideBrackets02   <- (rx * phi * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu) * (rI + alpha))
  
  R00[k1]            <- OutsideBrackets * InsideBrackets00
  R11[k1]            <- OutsideBrackets * InsideBrackets11
  R22[k1]            <- OutsideBrackets * InsideBrackets22
  R01[k1]            <- OutsideBrackets * InsideBrackets01
  R12[k1]            <- OutsideBrackets * InsideBrackets12
  R02[k1]            <- OutsideBrackets * InsideBrackets02
  
}

R012012              <- R00 + R11 + R22 + R01 + R12 + R02

store2               <- cbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

storeb               <- rbind(R012012/R012012[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

storeb1              <- rbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1])

storeb2              <- rbind(R01/R01[1], R12/R12[1], R02/R02[1])

rT1                  <- 2 * rI # set back to default

# ====================================================================================
# rx 
# ====================================================================================

vec                  <- seq(0.01,1,length.out=Npts) # rx and phi

for (k1 in 1:Npts){
  
  rx                 <- vec[k1]
  
  OutsideBrackets    <- a^2 * b * c * m / muhat * exp(-muhat * tauhat)
  
  InsideBrackets00   <- (rT0 + rx + phi + mu)/((rI + rx + alpha) * (rT0 + phi + mu))
  InsideBrackets11   <- (rT1 + rx + phi + mu)/((rI + rx + alpha) * (rT1 + phi + mu))
  InsideBrackets22   <- 1/(rI + alpha)
  InsideBrackets01   <- (rx * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu))
  InsideBrackets12   <- (rx * phi)/((rI + rx + alpha) * (rT1 + phi + mu) * (rI + alpha))
  InsideBrackets02   <- (rx * phi * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu) * (rI + alpha))
  
  R00[k1]            <- OutsideBrackets * InsideBrackets00
  R11[k1]            <- OutsideBrackets * InsideBrackets11
  R22[k1]            <- OutsideBrackets * InsideBrackets22
  R01[k1]            <- OutsideBrackets * InsideBrackets01
  R12[k1]            <- OutsideBrackets * InsideBrackets12
  R02[k1]            <- OutsideBrackets * InsideBrackets02  
  
}

R012012              <- R00 + R11 + R22 + R01 + R12 + R02

store3               <- cbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

storec               <- rbind(R012012/R012012[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

storec1              <- rbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1])

storec2              <- rbind(R01/R01[1], R12/R12[1], R02/R02[1])

rx                   <- 0.5 # set back to default

# ====================================================================================
# phi
# ====================================================================================

vec                  <- seq(0.01,1,length.out=Npts) # rx and phi

for (k1 in 1:Npts){
  
  phi                <- vec[k1]
  
  OutsideBrackets    <- a^2 * b * c * m / muhat * exp(-muhat * tauhat)
  
  InsideBrackets00   <- (rT0 + rx + phi + mu)/((rI + rx + alpha) * (rT0 + phi + mu))
  InsideBrackets11   <- (rT1 + rx + phi + mu)/((rI + rx + alpha) * (rT1 + phi + mu))
  InsideBrackets22   <- 1/(rI + alpha)
  InsideBrackets01   <- (rx * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu))
  InsideBrackets12   <- (rx * phi)/((rI + rx + alpha) * (rT1 + phi + mu) * (rI + alpha))
  InsideBrackets02   <- (rx * phi * phi)/((rI + rx + alpha) * (rT0 + phi + mu) * (rT1 + phi + mu) * (rI + alpha))
  
  R00[k1]            <- OutsideBrackets * InsideBrackets00
  R11[k1]            <- OutsideBrackets * InsideBrackets11
  R22[k1]            <- OutsideBrackets * InsideBrackets22
  R01[k1]            <- OutsideBrackets * InsideBrackets01
  R12[k1]            <- OutsideBrackets * InsideBrackets12
  R02[k1]            <- OutsideBrackets * InsideBrackets02
  
}

R012012              <- R00 + R11 + R22 + R01 + R12 + R02

store4               <- cbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

stored               <- rbind(R012012/R012012[1], 
                        R01/R01[1], R12/R12[1], R02/R02[1])

stored1              <- rbind(R012012/R012012[1], 
                        R00/R00[1], R11/R11[1], R22/R22[1])

stored2              <- rbind(R01/R01[1], R12/R12[1], R02/R02[1])

# ====================================================================================
# put together for x axis to be reprod. no.
# ====================================================================================

Results              <- as.data.frame(rbind(store1, store2, store3, store4))
Variable             <- c(rep("rT_0",Npts), rep("rT_0.5",Npts), rep("rx",Npts), rep("phi",Npts))
X                    <- rep(seq(0,1,length.out=Npts),4)
Results              <- cbind(Variable, X, Results)
colnames(Results)    <- c("Variable", "X", "R012012","R00","R11","R22","R01","R12","R02")

# ====================================================================================

presults             <- ggplot(Results, aes(x=X, y=R012012, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_continuous(name="Change in reprod. no.") + 
  scale_x_continuous(name="Change in variable") + 
  scale_color_discrete(labels = unname(TeX(c("$\\phi", "$r_{T_P}", "$r_{T_S}", "$r_{x}")))) + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.title=element_blank()) +
  theme(legend.text.align = 0) +
  theme(legend.text = element_text(size=30))

# ====================================================================================

presults             <- ggplot(Results, aes(x=X, y=R00, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_continuous(name="Change in R0_0") + 
  scale_x_continuous(name="Change in variable") + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.position="none")

# ====================================================================================

presults             <- ggplot(Results, aes(x=X, y=R11, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_continuous(name="Change in R1_1") + 
  scale_x_continuous(name="Change in variable") + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.position="none")

# ====================================================================================

presults             <- ggplot(Results, aes(x=X, y=R22, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_continuous(name="Change in R2_2") + 
  scale_x_continuous(name="Change in variable") + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.position="none")

# ====================================================================================

presults             <- ggplot(Results, aes(x=X, y=R01, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_continuous(name="Change in R0_1") + 
  scale_x_continuous(name="Change in variable") + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.position="none")

# ====================================================================================

presults             <- ggplot(Results, aes(x=X, y=R12, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_continuous(name="Change in R1_2") + 
  scale_x_continuous(name="Change in variable") + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.position="none")

# ====================================================================================

Results              <- as.data.frame(rbind(store1, store2, store3, store4))
Variable             <- c(rep("rT_S*",Npts), rep("rT_P*",Npts), rep("rx",Npts), rep("phi",Npts))
X                    <- rep(seq(0,1,length.out=Npts),4)
Results              <- cbind(Variable, X, Results)
colnames(Results)    <- c("Variable", "X", "R012012","R00","R11","R22","R01","R12","R02")

presults             <- ggplot(Results, aes(x=X, y=R02, group=Variable)) + geom_line(aes(color=Variable), size=2)
presults + scale_y_log10(name= unname(TeX("Change in $R_{S \\rightarrow R}") )) + 
  scale_x_log10(name="Change in variable") + 
  scale_color_discrete(labels = unname(TeX(c("$\\phi", "$r_{T_P}*", "$r_{T_S}*", "$r_{x}")))) + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  theme(legend.text.align = 0) +
guides(color=guide_legend(title="")) + 
theme(legend.text=element_text(size=30)) 

# ====================================================================================
# put together for x axis to be reprod. no.
# ====================================================================================

Results2              <- as.data.frame(cbind(storea, storeb, storec, stored))
R                     <- c(rep("Overall",Npts), rep("SS_SR",Npts), rep("SR_RR",Npts), rep("SS_RR",Npts))
X                     <- rep(seq(0.01,1,length.out=Npts), 7)
Results2              <- cbind(R, X, Results2)
colnames(Results2)    <- c("R", "X", "rT_SS", "rT_SR","rx","phi")

# ====================================================================================

Results2              <- as.data.frame(cbind(storea, storeb, storec, stored))
R                     <- c(rep("Overall",Npts), rep("S_P*",Npts), rep("P_R*",Npts), rep("S_R*",Npts))
X                     <- rep(seq(0.01,1,length.out=Npts), 4)
Results2              <- cbind(R, X, Results2)
colnames(Results2)    <- c("R", "X", "rT_S", "rT_P","rx","phi")
# So that legend has stars to indicate which lines are indentical
temp                  <- levels(Results2$R)
Results2$R            <- factor(Results2$R, levels = temp[c(4,2,3,1)])
presults              <- ggplot(Results2, aes(x=X, y=rx, group=R, order = as.numeric(Results2$R))) + geom_line(aes(color=R), size = 1.2)
presults + scale_y_continuous(name="Change in reprod. no.") + 
  scale_x_continuous(name= unname(TeX("$r_x"))) +
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  scale_color_brewer(palette = "Set1", labels = unname(TeX(c("$R_{S \\rightarrow R}*", "$R_{P \\rightarrow R}*", "$R_{S \\rightarrow P}*", "$Overall")))) + 
  theme(legend.title=element_blank()) + 
  theme(legend.text.align = 0) +
  theme(legend.text=element_text(size=30))

# ====================================================================================

Results2              <- as.data.frame(cbind(storea, storeb, storec, stored))
R                     <- c(rep("Overall",Npts), rep("S_P",Npts), rep("P_R",Npts), rep("S_R",Npts))
X                     <- rep(seq(0.01,1,length.out=Npts), 4)
Results2              <- cbind(R, X, Results2)
colnames(Results2)    <- c("R", "X", "rT_S", "rT_P","rx","phi")
# So that legend has stars to indicate which lines are indentical - but nothing identical here
temp                  <- levels(Results2$R)
Results2$R            <- factor(Results2$R, levels = temp[c(4,2,3,1)])

presults              <- ggplot(Results2, aes(x=X, y=phi, group=R)) + geom_line(aes(color=R), size=2)
presults + scale_y_continuous(name="Change in reprod. no.") + 
  scale_x_continuous(name= unname(TeX("$\\phi"))) + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  scale_color_brewer(palette = "Set1", labels = unname(TeX(c("$R_{S \\rightarrow R}", "$R_{P \\rightarrow R}", "$R_{S \\rightarrow P}", "$Overall")))) + 
  theme(legend.title=element_blank()) + 
  theme(legend.text.align = 0) +
  theme(legend.text=element_text(size=30))

# ====================================================================================
Results2               <- as.data.frame(cbind(storea, storeb, storec, stored))
R                      <- c(rep("Overall",Npts), rep("S_P*",Npts), rep("P_R",Npts), rep("S_R*",Npts))
X                      <- rep(seq(0.01,1,length.out=Npts), 4)
Results2               <- cbind(R, X, Results2)
colnames(Results2)     <- c("R", "X", "rT_S", "rT_P","rx","phi")
# So that legend has stars to indicate which lines are indentical - but nothing identical here
temp                   <- levels(Results2$R)
Results2$R             <- factor(Results2$R, levels = temp[c(4,2,3,1)])

presults               <- ggplot(Results2, aes(x=X, y=rT_S, group=R)) + geom_line(aes(color=R), size=2)
presults + scale_y_continuous(name="Change in reprod. no.") + 
  scale_x_continuous(name = unname(TeX("$r_{T_S}"))) + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  scale_color_brewer(palette = "Set1", labels = unname(TeX(c("$R_{S \\rightarrow R}*", "$R_{P \\rightarrow R} ", "$R_{S \\rightarrow P}*", "$Overall")))) + 
  theme(legend.title=element_blank()) + 
  theme(legend.text.align = 0) +
  theme(legend.text=element_text(size=30))

# ====================================================================================

Results2                <- as.data.frame(cbind(storea, storeb, storec, stored))
R                       <- c(rep("Overall",Npts), rep("S_P*",Npts), rep("P_R*",Npts), rep("S_R*",Npts))
X                       <- rep(seq(0.01,1,length.out=Npts), 4)
Results2                <- cbind(R, X, Results2)
colnames(Results2)      <- c("R", "X", "rT_S", "rT_P","rx","phi")
# So that legend has stars to indicate which lines are indentical - but nothing identical here
temp                    <- levels(Results2$R)
Results2$R              <- factor(Results2$R, levels = temp[c(4,2,3,1)])

colnames(Results2)      <- c("R", "X", "rT_S", "rT_P","rx","phi")
# So that legend has stars to indicate which lines are indentical
presults                <- ggplot(Results2, aes(x=X, y=rT_P, group=R)) + geom_line(aes(color=R), size=2)
presults + scale_y_continuous(name="Change in reprod. no.") + 
  scale_x_continuous(name = unname(TeX("$r_{T_P}"))) + 
  theme(axis.title.x = element_text(margin = margin(t = 20), size = 30)) +
  theme(axis.title.y = element_text(margin = margin(r = 20), size = 30)) +
  theme(axis.text.x = element_text(size=30), 
        axis.text.y = element_text(size=30)) + 
  scale_color_brewer(palette = "Set1", labels = unname(TeX(c("$R_{S \\rightarrow R}*", "$R_{P \\rightarrow R}*", "$R_{S \\rightarrow P}*", "$Overall")))) + 
  theme(legend.title=element_blank()) + 
  theme(legend.text.align = 0) +
  theme(legend.text=element_text(size=30))

# ====================================================================================
# put together for x axis to be reprod. no. - ONLY STATIC ONES
# ====================================================================================

Results4                <- as.data.frame(cbind(storea1, storeb1, storec1, stored1))
R                       <- c(rep("R012012",Npts), rep("R00",Npts), rep("R11",Npts), rep("R02",Npts))
X                       <- rep(seq(0,1,length.out=Npts), 4)
Results4                <- cbind(R, X, Results4)
colnames(Results4)      <- c("R", "X", "rT_SS", "rT_SR", "rx","phi")

# ====================================================================================
