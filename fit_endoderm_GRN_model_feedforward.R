# "Modified Gauss-Newton algorithm" to fit endoderm GRN Feed-Forward model to transcriptomics data
# From "Recursive feedforward loops govern robust cell fate lockdown during C. elegans embryogenesis"
# by Chee Kiang Ewe, Erica M. Sommermann, Josh Kenchel, Sagen E. Flowers, Morris F. Maduro, and Joel H. Rothman
# UC Santa Barbara, UC Riverside, and UC Los Angeles
# Updated 13 March 2022 by Josh Kenchel
# jkenchel@ucsb.edu

# *** Recommended to run in R STUDIO to keep all plots in one window. ***
# Program takes parameter guesses, changes some parameter values randomly, runs model,
# calculates sum of squared differences, and keeps new values if less than previous sum of squares.
# Program iterates nGen times, changing parameter values with probability mu.
# New parameter values are drawn from normal distribution with mean of current value and SD of mutationSD multiplied by current value.

# IF RUNNING PROGRAM FOR FIRST TIME:
#   It should be okay to run as-is.
# IF ITERATING ON PREVIOUSLY EVOLVED VALUES STORED IN bestParamGuesses IN THE WORKSPACE:
# set this value to TRUE
iterating <- TRUE

#Sets least squares value to Infinite if this is the first time running the program in this workspace
if(!iterating){
  leastSquares <- Inf
}

# Settings for evolution (fitting) of model parameters
nGen <- 1000 #number of generations (cycles)
mu <- 0.05 #probability that a parameter gets mutated to a new value

# Given mutation, the standard deviation (as a percent of current value) of the normal distribution
# with mean of the current value from which to draw a new value
mutationSD <- 0.1

# Parameter names described in the supplementary text and files
paramNames <- c("delta", #Native degradation rate
                "a1", #SKN-1 activation of med-2
                "a2", #SKN-1 activation of med-1
                "a3", #SKN-1 activation of end-3
                "a4", #SKN-1 activation of end-1
                "a5", #MED-2 activation of med-1
                "a6", #MED-2 activation of end-3
                "a7", #MED-1 activation of end-3
                "a8", #MED-1 activation of end-1
                "a9", #POP-1 activation of end-3
                "a10", #POP-1 activation of end-1
                "a11", #END-3 activation of end-1
                "a12", #END-3 activation of elt-7
                "a13", #END-1 activation of elt-7
                "a14", #END-3 activation of elt-2
                "a15", #END-1 activation of elt-2
                "a16", #ELT-7 activation of elt-2
                "f1", #ELT-7 feedback
                "f2", #ELT-2 feedback
                "theta1", #Threshold required for ELT-7 feedback (ppm)
                "theta2", #Threshold required for ELT-2 feedback (ppm)
                "tau") #Transcription-activation delay (min)

# Parameter values used in the paper and listed in the supplementary file
bestParams <- c(
  9.213993e-02, #delta, Native degradation rate
  4.867228e-02, #a1, SKN-1 activation of med-2
  1.771412e-02, #a2, SKN-1 activation of med-1
  8.393703e-18, #a3, SKN-1 activation of end-3
  4.830789e-17, #a4, SKN-1 activation of end-1
  2.733470e-02, #a5, MED-2 activation of med-1
  8.333903e-03, #a6, MED-2 activation of end-3
  2.272943e-01, #a7, MED-1 activation of end-3
  9.550453e-18, #a8, MED-1 activation of end-1
  1.084313e-01, #a9, POP-1 activation of end-3
  9.011876e-16, #a10, POP-1 activation of end-1
  2.666866e-01, #a11, END-3 activation of end-1
  8.280570e-02, #a12, END-3 activation of elt-7
  5.893009e-02, #a13, END-1 activation of elt-7
  2.862568e-16, #a14, END-3 activation of elt-2
  1.912201e-16, #a15, END-1 activation of elt-2
  5.574331e-02, #a16, ELT-7 activation of elt-2
  7.026551e-02, #f1, ELT-7 feedback
  2.699666e-01, #f2, ELT-2 feedback
  4.403810e+01, #theta1, Threshold required for ELT-7 feedback (ppm)
  4.539943e+00, #theta2, Threshold required for ELT-2 feedback (ppm)
  3.567803e+01 #tau, Transcription-activation delay (min)
)

#If iterating on previously fitted values, parameter guesses will be drawn from bestParamGuesses data frame in the workspace.
#If running for the first time, bestParamGuesses will be reset to the published values.
if(!iterating){
  bestParamGuesses <- data.frame(parameter = paramNames, guess = bestParams)
}

# This line draws the parameter values back out of the bestParamGuesses data frame (to allow for iteration)
bestParams <- bestParamGuesses[, "guess"]
nParams <- length(bestParams)

#Switches to make knockout mutants
#Currently not used -- leave as FALSE
skn1ko <- FALSE
pop1ko <- FALSE
med2ko <- FALSE
med1ko <- FALSE
end3ko <- FALSE
end1ko <- FALSE
elt7ko <- FALSE
elt2ko <- FALSE

# Time settings
tEnd <- 200 #min
dt <- 0.01 #time step for Euler approximation
nTimepoints <- tEnd/dt + 1
times <- seq(0, tEnd, by = dt)

#Initial conditions (ppm)
skn1_0 <- 23
pop1_0 <- 3
med2_0 <- 13.8
med1_0 <- 5.6
end3_0 <- 5.3
end1_0 <- 0
elt7_0 <- 0
elt2_0 <- 0

#SKN-1 and POP-1 settings for square wave
skn1offTime <- 23 #min
pop1onTime <- 23 #min
pop1offTime <- 41 #min

####################################################################################################

# Load transcriptomics data used to fit model (units of time and ppm)
time <- c(0,23,41,53,66,83,101,122,143,186)
med2 <- c(13.8,13.8,0,0,0,0,0,0,0,0)
med1 <- c(5.6,5.3,0.4,0,0,0,0,0,0,0)
end3 <- c(5.3,11.0,18.0,17.0,17.0,23.0,0,0,0,0)
end1 <- c(0,7,17,26,35,49,62,36,45,0)
elt7 <- c(0,0.5,9,12,27,24,38,47,36,27)
elt2 <- c(0,0,0,0,0,0,0,0,53,134)
FFdataFit <- data.frame(time,med2,med1,end3,end1,elt7,elt2)

####################################################################################################

#Begin evolution (fitting) of model with nGen generations (cycles) of iteration
for(gen in 1:nGen){
  
  #Print progress to the console every 100 generations
  if(gen/100 == round(gen/100)){
    cat("Simulation ", gen, " of ", nGen, ", sum of squares = ", leastSquares, "\n")
  }
  
  #Transfer parameter values to another vector for mutation and testing
  testParams <- bestParams
  
  # Introduce mutations (change random parameter values randomly)
  mutation <- runif(nParams) < mu #Creates logical vector of length nParams dictating whether each parameter will be changed (TRUE) or not (FALSE)
  newParams <- abs(rnorm(nParams, mean = 1, sd = mutationSD)*testParams) # new values to be assigned to parameters for which mutation is TRUE
  testParams <- testParams*(!mutation) + newParams*mutation #combine current values for unchanged parameters with new values for mutated parameters
  
  # Store parameter values to be tested in their respective named vectors
  delta <- testParams[1]
  a1 <- testParams[2]
  a2 <- testParams[3]
  a3 <- testParams[4]
  a4 <- testParams[5]
  a5 <- testParams[6]
  a6 <- testParams[7]
  a7 <- testParams[8]
  a8 <- testParams[9]
  a9 <- testParams[10]
  a10 <- testParams[11]
  a11 <- testParams[12]
  a12 <- testParams[13]
  a13 <- testParams[14]
  a14 <- testParams[15]
  a15 <- testParams[16]
  a16 <- testParams[17]
  f1 <- testParams[18]
  f2 <- testParams[19]
  theta1 <- testParams[20]
  theta2 <- testParams[21]
  tau <- testParams[22]
  
  #Calculate number of timepoints represented by transcription-activation delay tau
  txnDelayTimepoints <- round(tau/dt)
  
  #SKN-1 and POP-1 initial on or off settings
  skn1on <- TRUE
  pop1on <- FALSE

  #Creating vectors to store expression level time series for all genes
  med2 <- numeric(nTimepoints)
  med1 <- numeric(nTimepoints)
  end3 <- numeric(nTimepoints)
  end1 <- numeric(nTimepoints)
  elt7 <- numeric(nTimepoints)
  elt2 <- numeric(nTimepoints)
  med2[1] <- med2_0*(-1*(med2ko - 1))
  med1[1] <- med1_0*(-1*(med1ko - 1))
  end3[1] <- end3_0*(-1*(end3ko - 1))
  end1[1] <- end1_0*(-1*(end1ko - 1))
  elt7[1] <- elt7_0*(-1*(elt7ko - 1))
  elt2[1] <- elt2_0*(-1*(elt2ko - 1))
  
  # This is the model code, mostly the same as in endoderm_GRN_model.R
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
    
    #Conditions at t_i
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
  
  # Calculate sum of squared differences between model predictions and transcriptomics data
  simData <- data.frame(med2, med1, end3, end1, elt7, elt2)
  sumOfSquares <- 0
  for(i in 1:nrow(FFdataFit)){
    fitTime <- FFdataFit[i,"time"]
    sumOfSquares <- sumOfSquares + sum((FFdataFit[i, 2:7] - simData[which(times == fitTime),])^2)
  }
  
  # If sum of squares is less than before, replace parameter values with testParams and plot the model run
  # Otherwise, bestParams remains intact
  if(sumOfSquares < leastSquares){
    leastSquares <- sumOfSquares
    bestParams <- testParams
    plot(y=med2, x=times, type = "l", ylab = "[Transcript] (ppm)", xlab = "Time post 4-cell stage (min)", 
         col = "red", ylim = c(0,200), xlim = c(0,200),
         main = graphTitle)
    lines(y=med1, x=times, type = "l", col = "orange")
    lines(y=end3, x=times, type = "l", col = "green")
    lines(y=end1, x=times, type = "l", col = "skyblue")
    lines(y=elt7, x=times, type = "l", col = "blue")
    lines(y=elt2, x=times, type = "l", col = "purple")
    legend("topleft", legend = c("med-2", "med-1", "end-3", "end-1", "elt-7", "elt-2"),
           col = c("red", "orange", "green", "skyblue", "blue", "purple"), lwd = 1)
  }
}

#Make graph of gene expression levels for final parameter values
plot(y=med2, x=times, type = "l",
     ylab = "[Transcript] (ppm)", xlab = "Time post 4-cell stage (min)", 
     col = "red", ylim = c(0,200), xlim = c(0,200),
     main = graphTitle)
lines(y=med1, x=times, type = "l", col = "orange")
lines(y=end3, x=times, type = "l", col = "green")
lines(y=end1, x=times, type = "l", col = "skyblue")
lines(y=elt7, x=times, type = "l", col = "blue")
lines(y=elt2, x=times, type = "l", col = "purple")
legend("topleft", legend = c("med-2", "med-1", "end-3", "end-1", "elt-7", "elt-2"),
       col = c("red", "orange", "green", "skyblue", "blue", "purple"), lwd = 1)

# Store evolved parameters and print to console
bestParamGuesses <- data.frame(parameter = paramNames, guess = bestParams)
print(bestParamGuesses)