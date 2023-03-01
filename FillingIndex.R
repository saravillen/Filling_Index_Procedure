# Information ------------------------------------------------------------------
#     Authors: Fernanda Alves-Martins, Sara Villén-Pérez, Ignacio Morales-Castilla
#     Title: Proof-of-concept Filling index procedure
#     Date: June 2022
#     Last update: November 11 2022
#     
#     Workflow: 
#       1. fi: Function to calculate plot filling index
#       2. null_model_fi: Function that simulate line-shaped plots, generate a 
#          null distribution of filling indexes, and compare the observed 
#          filling index to the simulated plots  
#       3. Running "null_model_fi" function



# Libraries --------------------------------------------------------------------
libs = c("dplyr","raster","quantreg","LaplacesDemon","extras","beepr")
lapply(libs, require, character.only = TRUE)
rm(libs)

# Filling index function--------------------------------------------------------
fi = function(df, Xvariable, Yvariable){
  options(warn = -1)
  xmin = min(df[, Xvariable], na.rm = TRUE)
  xmax = max(df[, Xvariable], na.rm = TRUE)
  ymin = 0
  ymax = max(df[, Yvariable], na.rm = TRUE)
  
  # This function is to create the raster comprising the X and Y
  env_space = raster(xmn = xmin, xmx = xmax,
                     ymn = ymin, ymx = ymax,
                     res = c((xmax-xmin)/50, (ymax-ymin)/50)) # Set the raster resolution
  values(env_space) = 0
  
  coords_IB = raster::extract(env_space, y = na.omit(df[,c(Xvariable,Yvariable)]),cellnumbers = TRUE)
  n_IB = table(coords_IB)
  env_space_IB = env_space
  values(env_space_IB)[as.numeric(names(n_IB))] = n_IB
  IB_values = values(env_space_IB)
  IB_values[IB_values == 0] = NA
  values(env_space_IB) = IB_values
  
  # Convert the raster to a dataframe
  d = data.frame(coordinates(env_space_IB), count=env_space_IB[])
  d$X_class = as.integer(cut(d$x, breaks=50,labels=c(1:50)))
  d$Y_class = as.integer(cut(d$y, breaks=50, labels=c(1:50)))
  
  df_pixel = data.frame()
  
  for (j in 1:50){
    
    df.j = d[which(d$X_class == j),]
    id=max(sort(df.j$Y_class[which(!is.na(df.j$count))], decreasing = T))
    df.j = subset(df.j, df.j$Y_class <= id)
    pix_color = sum(!is.na(df.j$count))
    pix_no_color = sum(is.na(df.j$count))
    df_pixel2 = as.data.frame(cbind(pix_color,pix_no_color))
    df_pixel = rbind(df_pixel,df_pixel2)
  }
  
  prop_pixels_V1 =  sum(df_pixel$pix_color)/sum(df_pixel)
  return(prop_pixels_V1)
}

# Null model Filling index function --------------------------------------------

null_model_fi = function(data,Xvariable,Yvariable,iterations){
  
  fi_sim_all = data.frame()
  filling_index = data.frame()
  
  spdata = datos[ ,c(Yvariable,Xvariable)]
  
  # Observed filling index
  fi_obs = fi(spdata,Xvariable,Yvariable)
  
  n_spsi = nrow(spdata) # sample size of the species (number of observations)
  
  # Transform abundance data to avoid negative predicted abundances
  z_spdata = spdata
  z_spdata[,Yvariable] = z_spdata[,Yvariable]/1000000
  small = min(z_spdata[z_spdata[,Yvariable] > 0,Yvariable])/2 # small value to add before transformation
  z_spdata[,Yvariable] =  z_spdata[,Yvariable] + small
  z_spdata[,Yvariable] = logit( z_spdata[,Yvariable])
  
  # Create quadratic and cubic terms of predictors
  z_spdata = within(z_spdata,{
    X2 = z_spdata[,Xvariable]^2})
  colnames (z_spdata) = c(Yvariable, Xvariable, paste(Xvariable,"_2", sep=""))
  
  # Formula
  passFormula = function(df, Xvariable, Yvariable){
    fit = rq(df[,Yvariable] ~ df[,Xvariable] + df[,paste(Xvariable,"_2", sep="")], tau = 0.5)
  }
  
  qr50 = passFormula(z_spdata,Xvariable, Yvariable)
  
  ## 2 ## Simulation of line-shaped pattern
  ## 2.1. Predictors for partial effect (sample size equal to original sample size)
  predictors = data.frame(spdata[,Xvariable],
                          spdata[,Xvariable]^2)
  colnames (predictors) = c(paste(Xvariable),
                            paste(Xvariable,"_2", sep=""))
  
  
  ## 2.2. Simulation of average response (simulated line-shaped pattern without error)
  pred_varj_50 = predict(qr50, predictors)
  pred_varj_50 = (invlogit(pred_varj_50)-small)*1000000 #re-transform abundance data to the original scale
  
  
  for(k in 1:iterations){
    
    ## 2.3. Simulation of error around simulated linear pattern
    ## 2.3.1. Estimation of the maximum predicted abundance at median (tau 0.5)
    max_ab_pred = sqrt(max(pred_varj_50))
    
    ## 2.3.2. Error around simulated linear pattern
    error = rnorm(n=n_spsi, mean=0, sd=max_ab_pred)
    
    ## 2.4. Average response + error (x iterations possible errors)
    pred_varj_50_err = pred_varj_50 + error
    pred_varj_50_err[pred_varj_50_err<0] = 0 # convert possible negative abundances to zero
    pred_varj_50_err = as.data.frame(cbind(z_spdata[,Xvariable], pred_varj_50_err))
    colnames (pred_varj_50_err)  = c(Xvariable,Yvariable)
    
    ## 3 ## Simulated filling indexes
    fi_sim = fi(pred_varj_50_err,Xvariable,Yvariable)
    fi_sim_all = rbind(fi_sim_all,fi_sim)
    colnames (fi_sim_all) = "fi_sim"
    
  }
  
  ## 4 ## Probability that observed filling index is significantly higher than that of simulated linear patterns (i.e., that the observed pattern is not linear)
  pValue = (length(which(fi_sim_all>=fi_obs))+1)/(iterations+1)
  
  if (pValue < 0.05){
    print(paste("p value =", pValue))
    print("Polygonal-shaped pattern: Observed filling index is significantly higher than the simulated ones")
  }
  
  else{
    print(paste("p value =", pValue))
    print("No Polygonal-shaped pattern: Observed filling index is statistically similar to the simulated ones")
  }
  
}


# Load data --------------------------------------------------------------------
datos = read.csv("fagus_grand.csv", header= T, sep=";")

# Running the function ---------------------------------------------------------
## the function arguments to be declared are:
## - dataframe
## - predictor column name within quotes
## - response variable column name within quotes
## - iterations
null_model_fi(datos,"gdd","abundance",999) #It should take some time
beep(2) # notifies when the function return results
