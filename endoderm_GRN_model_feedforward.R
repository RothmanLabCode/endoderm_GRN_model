# C. elegans endoderm development GRN model - Feed-Forward model
# From "Recursive feedforward loops govern robust cell fate lockdown during C. elegans embryogenesis"
# by Chee Kiang Ewe, Erica M. Sommermann, Josh Kenchel, Sagen E. Flowers, Morris F. Maduro, and Joel H. Rothman
# UC Santa Barbara, UC Riverside, and UC Los Angeles
# Updated 2021 August 23 by Josh Kenchel
# jkenchel@ucsb.edu

# Program should be ready to run as-is.
# Flip logical vectors to make knockout mutants or adjust parameters as needed.
# Program will create a plot of gene expression levels and print elt-2 expression metrics to the console.

#Switches to make knockout mutants. Make TRUE to knock out. wt is all FALSE.
skn1ko <- FALSE
pop1ko <- FALSE
med2ko <- FALSE
med1ko <- FALSE
end3ko <- FALSE
end1ko <- FALSE
elt7ko <- FALSE
elt2ko <- FALSE

# Model parameters described in supplementary table
delta <- 9.213993e-02 # Native degradation rate
a1 <- 4.867228e-02 #SKN-1 activation of med-2
a2 <- 1.771412e-02 #SKN-1 activation of med-1
a3 <- 8.393703e-18 #SKN-1 activation of end-3
a4 <- 4.830789e-17 #SKN-1 activation of end-1
a5 <- 2.733470e-02 #MED-2 activation of med-1
a6 <- 2.333903e-02 #MED-2 activation of end-3 adjusted from fitted value of 8.333903e-03
a7 <- 2.272943e-01 #MED-1 activation of end-3
a8 <- 1.950453e-01 #MED-1 activation of end-1 adjusted from fitted value of 9.550453e-18
a9 <- 1.084313e-01 #POP-1 activation of end-3
a10 <- 9.011876e-16 #POP-1 activation of end-1
a11 <- 2.666866e-01 #END-3 activation of end-1
a12 <- 2.380570e-01 #END-3 activation of elt-7 adjusted from fitted value of 8.280570e-02
a13 <- 5.893009e-02 #END-1 activation of elt-7
a14 <- 2.862568e-16 #END-3 activation of elt-2
a15 <- 4.412201e-02 #END-1 activation of elt-2 adjusted from fitted value of 1.912201e-16
a16 <- 5.574331e-02 #ELT-7 activation of elt-2
f1 <- 7.026551e-02 #ELT-7 feedback
f2 <- 2.699666e-01 #ELT-2 feedback
theta1 <- 4.403810e+01 #Threshold required for ELT-7 feedback (ppm)
theta2 <- 4.539943e+00 #Threshold required for ELT-2 feedback (ppm)
tau <- 3.567803e+01 #Transcription-activation delay (min)

#Initial conditions (ppm)
skn1_0 <- 23
pop1_0 <- 3
med2_0 <- 13.8
med1_0 <- 5.6
end3_0 <- 5.3
end1_0 <- 0
elt7_0 <- 0
elt2_0 <- 0

#SKN-1 and POP-1 conditions for square waves
skn1on <- TRUE #at t=0
skn1offTime <- 23 #min
pop1on <- FALSE #at t=0
pop1onTime <- 23 #min
pop1offTime <- 41 #min

# Time settings
tEnd <- 1000 #total time of run, in min
dt <- 0.01 #Time step for Euler approximation
nTimepoints <- tEnd/dt + 1 #number of timepoints
times <- seq(0, tEnd, by = dt) #vector of all timepoints
txnDelayTimepoints <- round(tau/dt)

#Creating vectors to store expression level time series for all genes
med2 <- numeric(nTimepoints)
med1 <- numeric(nTimepoints)
end3 <- numeric(nTimepoints)
end1 <- numeric(nTimepoints)
elt7 <- numeric(nTimepoints)
elt2 <- numeric(nTimepoints)
med2[1] <- med2_0*!med2ko
med1[1] <- med1_0*!med1ko
end3[1] <- end3_0*!end3ko
end1[1] <- end1_0*!end1ko
elt7[1] <- elt7_0*!elt7ko
elt2[1] <- elt2_0*!elt2ko

# Loop through the number of timepoints given by the run time and time step dt
# The main part of the model
for(i in 1:(nTimepoints-1)){
  
  # Logical gates to create square waves for SKN-1 and POP-1
  if(times[i] > skn1offTime){
    skn1on <- FALSE
  }
  if(times[i] > pop1onTime){
    pop1on <- TRUE
  }
  if(times[i] > pop1offTime){
    pop1on <- FALSE
  }
  
  # conditions at t_i
  skn1 <- skn1_0*skn1on*!skn1ko
  pop1 <- pop1_0*pop1on*!pop1ko
  med2_i <- med2[i]
  med1_i <- med1[i]
  end3_i <- end3[i]
  end1_i <- end1[i]
  elt7_i <- elt7[i]
  elt2_i <- elt2[i]
  
  # Calculation of timepoint used for activation due to transcription-activation delay
  txnEffectTimepoint <- max(c(1, i - txnDelayTimepoints))
  med2_prev <- med2[txnEffectTimepoint]
  med1_prev <- med1[txnEffectTimepoint]
  end3_prev <- end3[txnEffectTimepoint]
  end1_prev <- end1[txnEffectTimepoint]
  elt7_prev <- elt7[txnEffectTimepoint]
  elt2_prev <- elt2[txnEffectTimepoint]
  
  # System of differential equations described in supplemental text
  # Conditions at t_i used to calculate rates of change
  dmed2 <- (skn1*a1 - delta*med2_i)*!med2ko
  dmed1 <- (skn1*a2 + med2_prev*a5 - delta*med1_i)*!med1ko
  dend3 <- (skn1*a3 + med1_prev*a7 + med2_prev*a6 + pop1*a9 - delta*end3_i)*!end3ko
  dend1 <- (skn1*a4 + end3_prev*a11 + med1_prev*a8 + pop1*a10 - delta*end1_i)*!end1ko
  delt7 <- (end1_prev*a13 + end3_prev*a12 + f1*(elt7_prev>theta1)*elt7_prev - delta*elt7_i)*!elt7ko
  delt2 <- (elt7_prev*a16 + end1_prev*a15 + end3_prev*a14 + f2*(elt2_prev>theta2)*elt2_prev - delta*elt2_i)*!elt2ko
  
  # Calculation of conditions at t_i+1 given the derivatives and time step dt
  med2[i+1] <- med2_i + dt*dmed2
  med1[i+1] <- med1_i + dt*dmed1
  end3[i+1] <- end3_i + dt*dend3
  end1[i+1] <- end1_i + dt*dend1
  elt7[i+1] <- elt7_i + dt*delt7
  elt2[i+1] <- elt2_i + dt*delt2
}

#Make graph of gene expression levels
plot(y=med2, x=times, type = "l",
     ylab = "[Transcript] (ppm)", xlab = "Time post 4-cell stage (min)", 
     col = "red", ylim = c(0,200), xlim = c(0,200),
     main = "Gene expression levels")
lines(y=med1, x=times, type = "l", col = "orange")
lines(y=end3, x=times, type = "l", col = "green")
lines(y=end1, x=times, type = "l", col = "skyblue")
lines(y=elt7, x=times, type = "l", col = "blue")
lines(y=elt2, x=times, type = "l", col = "purple")
legend("topleft", legend = c("med-2", "med-1", "end-3", "end-1", "elt-7", "elt-2"),
       col = c("red", "orange", "green", "skyblue", "blue", "purple"), lwd = 1)

#Print metrics of elt-2 timing
#Arbitrary concentration threshold used for comparative metrics, based on transcriptomics data
elt2timeToArbitraryConc <- min(which(floor(elt2) == 134))*dt
#Arbitrary timepoint chosen to quantify elt-2 timing, again based on transcriptomics data
elt2atArbitraryTimepoint <- elt2[min(which(floor(times) == 186))]
elt2at260min <- elt2[min(which(floor(times) == 260))]
cat("elt-2 time to 134 ppm transcript: ", elt2timeToArbitraryConc, "\n",
    "[elt-2] at 186 min: ", elt2atArbitraryTimepoint, "\n",
    "[elt-2] at 260 min: ", elt2at260min, "\n")