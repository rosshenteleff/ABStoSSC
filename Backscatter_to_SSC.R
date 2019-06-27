##Installing packages and libraries...

install.packages("tidyverse")
library(tidyr)

install.packages("plot3D")
library(plot3D)

install.packages("plotly")
library(plotly)
library(RColorBrewer)

install.packages("animation")
library(animation)

install.packages("readr")
library(readr)

install.packages("purrr")                    
library(purrr)                                          ##Note: This seems to fail sometimes. If the console reads "Error in library(purrr) : there is no package called ‘purrr’".
                                                        ##If this happens, go to C:\Users\s2448519\Documents\R\R-3.6.0\library (change user as required) and delete the purrr folder, and rerun the program.
install.packages("dplyr")
library(dplyr)
                                                        
##Bringing the data into the program...

##Amplitude Data:
Converse_A1 <- read.table("C:/Users/s2448519/Downloads/NS044_December2018/NS044_December2018/ADCP/ADPCon01.a1", header = FALSE, nrows = -1)
Converse_A2 <- read.table("C:/Users/s2448519/Downloads/NS044_December2018/NS044_December2018/ADCP/ADPCon01.a2", header = FALSE, nrows = -1)
Converse_A3 <- read.table("C:/Users/s2448519/Downloads/NS044_December2018/NS044_December2018/ADCP/ADPCon01.a3", header = FALSE, nrows = -1)

##Date-Time and Depth Data:
Converse_SEN <- read.table("C:/Users/s2448519/Downloads/NS044_December2018/NS044_December2018/ADCP/ADPCon01.sen", header = FALSE, nrows = -1)


##Settings of instrument, change as needed.
Blanking_Distance = 0.2  ##(m)
Bin_Size = 0.2           ##(m)
Bin_Number = 20          ##(count)
Salinity = 30            ##(ppt)
Freq = 2                 ##(MHz)
Noise_Floor = 13         ##(counts)
Max_Depth = 2.879        ##(m)


##Extracting depth values from SEN file...

Depth <- matrix(Converse_SEN[,14], nrow = nrow(Converse_SEN))
Depth_Matrix <- matrix(Depth, nrow = nrow(Depth), ncol = ncol(Converse_A1), byrow = FALSE)

##Note: Although the units of Col. 14 in Converse_SEN are actually dbar, this works out to be equal to the depth of seawater above the sensor in m.
##Note: Depth never goes above 2.879m, and frequently goes below blanking distance - cannot measure SSC in air or within blanking distance.


##Calculating SNR per beam...

SNR_Matrix1 <- matrix(0, nrow = nrow(Converse_A1), ncol = ncol(Converse_A1))

for(row in 1:nrow(Converse_A1)){
  for(col in 1:ncol(Converse_A1)){
    SNR_Matrix1[row,col] = 20 * log(Converse_A1[row,col] / Noise_Floor)
  }
}

SNR_Matrix2 <- matrix(0, nrow = nrow(Converse_A2), ncol = ncol(Converse_A2))

for(row in 1:nrow(Converse_A2)){
  for(col in 1:ncol(Converse_A2)){
    SNR_Matrix2[row,col] = 20 * log(Converse_A2[row,col] / Noise_Floor)   
  }
}

SNR_Matrix3 <- matrix(0, nrow = nrow(Converse_A3), ncol = ncol(Converse_A3))

for(row in 1:nrow(Converse_A3)){
  for(col in 1:ncol(Converse_A3)){
    SNR_Matrix3[row,col] = 20 * log(Converse_A3[row,col] / Noise_Floor)   
  }
}


##Filtering out data with SNR < 3dB...

for(row in 1:nrow(Converse_A1)){
  for(col in 1:ncol(Converse_A1)){
    if(SNR_Matrix1[row,col] < 3){
      Converse_A1[row,col] = 0
    }
  }
}

for(row in 1:nrow(Converse_A2)){
  for(col in 1:ncol(Converse_A2)){
    if(SNR_Matrix2[row,col] < 3){
      Converse_A3[row,col] = 0
    }
  }
}

for(row in 1:nrow(Converse_A3)){
  for(col in 1:ncol(Converse_A3)){
    if(SNR_Matrix3[row,col] < 3){
      Converse_A3[row,col] = 0
    }
  }
}


##Extracting time data from SEN file...

Month <- matrix(Converse_SEN[,1], nrow = nrow(Converse_SEN))
Day <- matrix(Converse_SEN[,2], nrow = nrow(Converse_SEN))
Year <- matrix(Converse_SEN[,3], nrow = nrow(Converse_SEN))
Hour <- matrix(Converse_SEN[,4], nrow = nrow(Converse_SEN))
Minute <- matrix(Converse_SEN[,5], nrow = nrow(Converse_SEN))
Second <- matrix(Converse_SEN[,6], nrow = nrow(Converse_SEN))

Date <- as.data.frame(cbind(Month, Day, Year))
Time <- as.data.frame(cbind(Hour, Minute, Second))

Date_comb <- unite(Date[1:3], comb_date, sep = "-", remove = TRUE)
Time_comb <- unite(Time[1:3], comb_time, sep = ":", remove = TRUE)

Date_Time <- as.data.frame(cbind(Date_comb, Time_comb))
Date_Time_comb <- unite(Date_Time[1:2], comb_date_time, sep = " ", remove = TRUE)


##Calculating range of each bin...

Range_Vector = seq((Blanking_Distance + (Bin_Size / 2) + (Bin_Size / 4)),((Bin_Size * Bin_Number) + 0.25), by = Bin_Size)
Range_Matrix = matrix(Range_Vector, nrow = nrow(Converse_A1), ncol = length(Range_Vector), byrow = TRUE)

##Note: Range per bin = distance to centre of bin + bin size/4.
##Note: Maximum length 4.15 - seems too long, but isn't, as blanking distance taken into account.


##Converting ranges to x-y-z coordinates for each beam...

##Beam 1:
x1 = -(Range_Vector) * sin(25)
y1 = Range_Vector * 0
z1 = Range_Vector * cos(25)

##Beam 2:
x2 = Range_Vector * sin(30) * -sin(25)
y2 = -(Range_Vector) * cos(30)
z2 = Range_Vector * cos(25)

##Beam 3:
x3 = Range_Vector * sin(30) * -sin(25)
y3 = Range_Vector * cos(30)
z3 = Range_Vector * cos(25)

##Note: Using coordinate system as described on p. 23 of manual.


##Calculating general spread of acoustic beam...

Acoustic_Spread = 20*log(Range_Vector)
AS_Matrix = matrix(Acoustic_Spread, nrow = nrow(Converse_A1), ncol = ncol(Converse_A1), byrow = TRUE)


##Calculating spread of acoustic beam due to water...

WA_30sal_3MHz = ((2.9 - 2.4)/(35 - 0))*Salinity + 2.4
WA_30sal_1.5MHz = ((0.7 - 0.6)/(35 - 0))*Salinity + 0.6
Water_Absorption = ((WA_30sal_3MHz - WA_30sal_1.5MHz)/(3 - 1.5))*Freq + WA_30sal_1.5MHz
Water_Spread = 2 * Water_Absorption * Range_Vector
WS_Matrix = matrix(Water_Spread, nrow = nrow(Converse_A1), ncol = ncol(Converse_A1), byrow = TRUE)

##Note: This is for Freq between 1.5 and 3 MHz, and Salinity between 0 and 35 ppt - if otherwise consult table @ https://www.nortekgroup.com/assets/documents/Sediments.pdf and change values as needed.


##Filtering out data for which amplitude = 0...

for(row in 1:nrow(Converse_A1)){
  for(col in 1:ncol(Converse_A1)){
    if(Converse_A1[row,col] == 0){
      AS_Matrix[row,col] = 0
      WS_Matrix[row,col] = 0
    }
  }
}


##Calculating fluid-corrected backscatter matrices...

Converse_FCB1 = 0.43 * Converse_A1 + AS_Matrix + WS_Matrix
Converse_FCB2 = 0.43 * Converse_A2 + AS_Matrix + WS_Matrix
Converse_FCB3 = 0.43 * Converse_A3 + AS_Matrix + WS_Matrix


##Calculating slopes of each row to get particle attenuation factor...        [Could this be done with map() rather than for()?]

Slope_Matrix1 <- matrix(0, nrow = nrow(Converse_FCB1))

for(row in 1:nrow(Converse_FCB1)){
  y <- as.numeric(Converse_FCB1[row,])
  x <- Range_Vector
  model <- lm(y ~ x)
  Slope_Matrix1[row,] <- as.numeric(model$coef[2])
}

Slope_Matrix2 <- matrix(0, nrow = nrow(Converse_FCB2))

for(row in 1:nrow(Converse_FCB2)){
  y <- as.numeric(Converse_FCB2[row,])
  x <- Range_Vector
  model <- lm(y ~ x)
  Slope_Matrix2[row,] <- as.numeric(model$coef[2])
}

Slope_Matrix3 <- matrix(0, nrow = nrow(Converse_FCB3))

for(row in 1:nrow(Converse_FCB3)){
  y <- as.numeric(Converse_FCB3[row,])
  x <- Range_Vector
  model <- lm(y ~ x)
  Slope_Matrix3[row,] <- as.numeric(model$coef[2])
}


##Calculating particle attenuation from slopes...

alphaP_1 = mean(Slope_Matrix1) * -0.5
alphaP_2 = mean(Slope_Matrix2) * -0.5
alphaP_3 = mean(Slope_Matrix3) * -0.5

##Note: This calculation is done following Wright et al.'s 2010 paper from the 2nd Joint Federal Interagency Conference, "Disciminating silt-and-clay from suspended sand in rivers using side-facing acoustic profilers".
##This is as follows: alphaP = -0.5*(dFCB/dr).

##Calculating particle spread matrices...

PS_Matrix1 = alphaP_1 * matrix(Range_Vector, nrow = nrow(Converse_FCB1), ncol = ncol(Converse_FCB1)) * 20
PS_Matrix2 = alphaP_2 * matrix(Range_Vector, nrow = nrow(Converse_FCB2), ncol = ncol(Converse_FCB2)) * 20
PS_Matrix3 = alphaP_3 * matrix(Range_Vector, nrow = nrow(Converse_FCB3), ncol = ncol(Converse_FCB3)) * 20


##Filtering out data for which amplitude = 0...

for(row in 1:nrow(PS_Matrix1)){
  for(col in 1:ncol(PS_Matrix1)){
    if(Converse_A1[row,col] == 0){
      PS_Matrix1[row,col] = 0
      PS_Matrix2[row,col] = 0
      PS_Matrix3[row,col] = 0
    }
  }
}

##Calculating final backscatter matrices...

BS_Matrix1 = Converse_FCB1 + PS_Matrix1
BS_Matrix2 = Converse_FCB2 + PS_Matrix2
BS_Matrix3 = Converse_FCB3 + PS_Matrix3


##Calculating suspended sediment concentration (SSC) from backscatter...           [Good candidate for map() rather than for() - what about if() statments?]

SSC_Matrix1 <- matrix(0, nrow = nrow(BS_Matrix1), ncol = ncol(BS_Matrix1))

for(row in 1:nrow(SSC_Matrix1)){
  for(col in 1:ncol(SSC_Matrix1)){
    SSC_Matrix1[row,col] <- 10^((BS_Matrix1[row,col]) * 0.025408151)
  }
}

max(SSC_Matrix1)

SSC_Matrix2 <- matrix(0, nrow = nrow(BS_Matrix2), ncol = ncol(BS_Matrix2))

for(row in 1:nrow(SSC_Matrix2)){
  for(col in 1:ncol(SSC_Matrix2)){
    SSC_Matrix2[row,col] <- 10^((BS_Matrix2[row,col]) * 0.025408151)
  }
}

max(SSC_Matrix2)

SSC_Matrix3 <- matrix(0, nrow = nrow(BS_Matrix3), ncol = ncol(BS_Matrix3))

for(row in 1:nrow(SSC_Matrix3)){
  for(col in 1:ncol(SSC_Matrix3)){
    SSC_Matrix3[row,col] <- 10^((BS_Matrix3[row,col]) * 0.025408151)
  }
}

max(SSC_Matrix3)

##Note: BS to SSC calculation done according to J.W. Gartner's 2004 article in Marine Geology 211, "Estimating suspended solids concentrations from backscatter intensity measured by acoustic Doppler current profiler in San Francisco Bay, California" (p. 181)
##This is as follows: SSC = 10^(a*BS + b), where a = 0.025408151, from linear modelling assuming the estimated and measured maxima correspond, and b = 0, as actual BS values were used.


##Creating x-y-z matrices...

X_Matrix1 <- matrix(x1, nrow = nrow(SSC_Matrix1), ncol = ncol(SSC_Matrix1), byrow = TRUE)
Y_Matrix1 <- matrix(y1, nrow = nrow(SSC_Matrix1), ncol = ncol(SSC_Matrix1), byrow = TRUE)
Z_Matrix1 <- matrix(z1, nrow = nrow(SSC_Matrix1), ncol = ncol(SSC_Matrix1), byrow = TRUE)

X_Matrix2 <- matrix(x2, nrow = nrow(SSC_Matrix2), ncol = ncol(SSC_Matrix2), byrow = TRUE)
Y_Matrix2 <- matrix(y2, nrow = nrow(SSC_Matrix2), ncol = ncol(SSC_Matrix2), byrow = TRUE)
Z_Matrix2 <- matrix(z2, nrow = nrow(SSC_Matrix2), ncol = ncol(SSC_Matrix2), byrow = TRUE)

X_Matrix3 <- matrix(x3, nrow = nrow(SSC_Matrix3), ncol = ncol(SSC_Matrix3), byrow = TRUE)
Y_Matrix3 <- matrix(y3, nrow = nrow(SSC_Matrix3), ncol = ncol(SSC_Matrix3), byrow = TRUE)
Z_Matrix3 <- matrix(z3, nrow = nrow(SSC_Matrix3), ncol = ncol(SSC_Matrix3), byrow = TRUE)


##Combining matrices...

X_Matrix <- cbind(X_Matrix1, X_Matrix2, X_Matrix3)
Y_Matrix <- cbind(Y_Matrix1, Y_Matrix2, Y_Matrix3)
Z_Matrix <- cbind(Z_Matrix1, Z_Matrix2, Z_Matrix3)

SSC_Matrix <- matrix(0, nrow = nrow(SSC_Matrix1), ncol = ncol(SSC_Matrix1))

for(row in 1:nrow(SSC_Matrix)){
  for(col in 1:ncol(SSC_Matrix)){
    SSC_Matrix[row,col] = mean(SSC_Matrix1[row,col], SSC_Matrix2[row,col], SSC_Matrix3[row,col])
  }
}

SSC_df <- as.data.frame(SSC_Matrix)
Z_df <- as.data.frame(Z_Matrix1)
SSC_df_long <- gather(SSC_df, value = "SSC", key = "bin")
Z_df_long <- gather(Z_df, value = "Z", key = "Z")
SSC_df_date_long <- as.data.frame(cbind(Date_Time_comb, SSC_df_long, Depth, Z_df_long[,2]))

colnames(SSC_df_date_long) <- c("Date", "bin", "SSC", "Depth", "Z")

for(row in 1:nrow(SSC_df_date_long)){
  if(SSC_df_date_long[row,"Z"] > SSC_df_date_long[row,"Depth"]){
    SSC_df_date_long[row,"SSC"] = 0
  }
}

##Creating plots...

palette <- colorRampPalette(c("black", "blue", "lightblue1", "lightblue2", "green", "yellow", "orange", "red", "darkred"))

SSC_plot = plot_ly(data = SSC_df_date_long, type = 'heatmap', x = ~Date, y = ~Z, z = ~SSC, colors = palette(nrow(SSC_df_date_long)))

SSC_plot



