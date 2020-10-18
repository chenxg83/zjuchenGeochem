# Seawater Density
#' value calculation
#'
#' @param Temperature the temperature of seawater. default: 25 degree C
#' @param P the pressure condition of seawater. default: 1 bar
#' @param S seawater salinity. default: 35
#'
#' @return seawater density at this P-T-S condition
#' @export
#'
#' @examples
#' seawater_density(Temperature = 100, P = 20, S = 35)
#' The output is: 0.9846648 with an unit of g cm-3
#' @description calculate seawater density as a function of T, P, and salinity using the equations proposed by Safarov et al., (2012,2013)

seawater_density <- function (Temperature = 25, P = 1, S = 35) {
  TK <- Temperature+273.15            #transform temperature in degree C to Calvin
  P_mpa <- P/10                       #transform pressure in bar to MPa

  a10 <- 16.8260119036                #all the following coefficients are shown in
  a11 <- -0.0141313351238             #Table 2 of Safarov et al., 2012
  a12 <- 0.247406490575/1000
  a20 <- -0.174551209401
  a21 <- 0.336290751342/10^4
  a22 <- -0.102404029138/10^5
  a30 <- 0.590549251552/10^3
  a31 <- 0
  a32 <- 0
  a40 <- -0.843286627505/10^6
  a41 <- 0
  a42 <- 0.482122273652/10^11
  a50 <- 0.444848916045/10^9
  a51 <- -0.346870313329/10^13
  a52 <- -0.534550469142/10^14
  A1 <- TK*(a10*S^0+a11*S+a12*S^2)
  A2 <- TK^2*(a20*S^0+a21*S+a22*S^2)
  A3 <- TK^3*(a30*S^0+a31*S+a32*S^2)
  A4 <- TK^4*(a40*S^0+a41*S+a42*S^2)
  A5 <- TK^5*(a50*S^0+a51*S+a52*S^2)
  A <- A1+A2+A3+A4+A5

  b00 <- 1780.94028549
  b01 <- -8.44625953277
  b02 <- -0.0203022580799
  b10 <- -23.9619103076
  b11 <- 0.0616367739239
  b12 <- 0
  b20 <- 0.123588063860
  b21 <- -0.115367580107/10^3
  b22 <- 0
  b30 <- -0.230808294277/10^3
  b31 <- 0
  b32 <- 0
  b40 <- 0.142119761663/10^6
  b41 <- 0
  b42 <- 0.534482054293/10^12
  B1 <- (b00+b01*S+b02*S^2)
  B2 <- TK*(b10+b11*S+b12*S^2)
  B3 <- TK^2*(b20+b21*S+b22*S^2)
  B4 <- TK^3*(b30+b31*S+b32*S^2)
  B5 <- TK^4*(b40+b41*S+b42*S^2)
  B <- B1+B2+B3+B4+B5

  c00 <- -1197.00357946
  c01 <- 8.48513299409
  c02 <- 0
  c10 <- 11.7622451237
  c11 <- -0.0571254359702
  c12 <- 0
  c20 <- -0.0402338524280
  c21 <- 0.958997239703/10^4
  c22 <- 0.185255385513/10^6
  c30 <- 0.418262515254/10^4
  c31 <- 0
  c32 <- 0
  c40 <- 0
  c41 <- 0
  c42 <- -0.986850131470/10^12
  C1 <- (c00+c01*S+c02*S^2)
  C2 <- TK*(c10+c11*S+c12*S^2)
  C3 <- TK^2*(c20+c21*S+c22*S^2)
  C4 <- TK^3*(c30+c31*S+c32*S^2)
  C5 <- TK^4*(c40+c41*S+c42*S^2)
  C <- C1+C2+C3+C4+C5

  # P_mpa = A*\rho^2+B*\rho^8+C*\rho^12
  # using polyroot to find the density of seawater, d (g cm-3)

  if (Temperature<0|Temperature>373|P<0|P>1400|S<0|S>55.529)
    density <- ("Temperature/pressure/salinity out of range!")
  else {
    densityfunction <- function (x) {
      A*x^2+B*x^8+C*x^12-P_mpa
      }

    #find the seawater density in a range of 0.8 ~ 1.2 g cm-3

    density <- uniroot(densityfunction,lower = 0.8,upper = 1.2)$root

  }

  return(density)
}

# Critical Pressure
#' data frame calculation
#'
#' @param Temperature the temperature of seawater. default: 25 degree C
#' @param type should be either "saturated" or "full". default: "saturated"
#'
#' @return seawater pressure at critical curve
#' @export
#'
#' @examples
#' critical_pressure(Temperature = 100, type="saturated")
#' The output is: 1.901956 with an unit of bar
#' @description calculate the critical pressure (bar) of seawater at given temperature using the equation proposed by Driesner and Heinrich (2007)

critical_pressure <- function (Temperature = 25, type = "saturated") {
  criticalPressureH2O <- 220.54915
  criticalTempH2O <- 373.976
  c1 <- -2.36
  c2 <- 1.28534/10
  c3 <- -2.3707/100
  c4 <- 3.20089/1000
  c5 <- -1.38917/10^4
  c6 <- 1.02789/10^7
  c7 <- -4.8376/10^11
  c1a <- 1
  c2a <- 1.5
  c3a <- 2
  c4a <- 2.5
  c5a <- 3
  c6a <- 4
  c7a <- 5
  if (type=="full") {
    Critical_P <- criticalPressureH2O+
      c1*(criticalTempH2O-Temperature)^c1a+
      c2*(criticalTempH2O-Temperature)^c2a+
      c3*(criticalTempH2O-Temperature)^c3a+
      c4*(criticalTempH2O-Temperature)^c4a+
      c5*(criticalTempH2O-Temperature)^c5a+
      c6*(criticalTempH2O-Temperature)^c6a+
      c7*(criticalTempH2O-Temperature)^c7a
  } else if(type=="saturated") {
    if (Temperature>=100) {
      Critical_P <- criticalPressureH2O+
        c1*(criticalTempH2O-Temperature)^c1a+
        c2*(criticalTempH2O-Temperature)^c2a+
        c3*(criticalTempH2O-Temperature)^c3a+
        c4*(criticalTempH2O-Temperature)^c4a+
        c5*(criticalTempH2O-Temperature)^c5a+
        c6*(criticalTempH2O-Temperature)^c6a+
        c7*(criticalTempH2O-Temperature)^c7a
    } else Critical_P <- 1
  } else Critical_P <- ("Wrong input: The type should be either 'full' or 'saturated'!")
  return(Critical_P)
}

# Quartz solubility
#' data frame calculation
#'
#' @param Temperature the temperature of seawater. default: 25 degree C
#' @param S the salinity of seawater. default: 35
#' @param ref the reference that used to calculate the quartz solubility. "vonDamm1991", "Gunnarsson2000". default: "vonDamm1991"
#'
#'
#' @return dissolved SiO2 concentration in the fluid that equilibrate with quartz at giving condition
#' @export
#'
#' @examples
#' quartz_sol(Temperature = 100, S = 35, ref = "vonDamm1991")
#' The output is: 1.003644 with an unit of mM
#' @description calculate the quartz solubility


quartz_sol <- function (Temperature = 25, S = 35, ref = "vonDamm1991") {

   TK <- Temperature+273.15
   if (ref == "vonDamm1991") {

    P <- critical_pressure(Temperature=Temperature, #calculate the critical pressure
                           type = "saturated")
    d <- seawater_density(Temperature, P, S)        #seawater density

                          #transform temperature in degree C to calvin
    ln_Sol <- -2.32888 + 1.79547 * log(d, base = exp(1)) + (-2263.62 + 0.00407350 * TK * TK) / TK + 0.0398808 * P / TK
    SiO2_mM <- 1000 * exp(ln_Sol)                  #SiO2 content in mM
  } else if (ref =="Gunnarsson2000") {

    logK_ppm <- -34.188+197.47/TK-5.851*10^-6*TK^2+12.245*log10(TK)
    SiO2_mM <- (10^logK_ppm)*1000

  } else print("wrong input")

  return(SiO2_mM)

}

# Dissolution constant of minerals
#' data frame calculation
#'
#' @param Temperature temperature condition. default: 25 degree C
#' @param mineral the specified minerals: anhydrite, calcite, barite, celestite
#'
#' @return dissolution constant of the specified mineral
#' @export
#'
#' @examples
#' dissolution_constant(Temperature = 100, mineral = "anhydrite")
#' The output is: -2.497612
#' @description calculate the dissolution constant of a specified mineral at given temperature

dissolution_constant <- function (Temperature = 25, mineral) {

  TK <- Temperature+273.15     #transform temperature in degree C to calvin

  anhydrite <- function (TK) {
    3.94-677/TK-12.39/1000*TK
  }

  calcite <- function(TK) {
    -1.46-41/TK-17.41/10^6*TK^2
  }

  barite <- function(TK) {
    lnK <- 275.053-43.014*log(TK)-15806.3/TK
    return(log10(exp(lnK)))
  }

  celestite <- function(TK) {
    lnK <- 224.069-35.9422*log(TK)-10302.32/TK
    return (log10(exp(lnK)))
  }

  if (mineral == "anhydrite")
    logK <- anhydrite(TK)
  else if (mineral == "calcite")
    logK <- calcite(TK)
  else if (mineral == "barite")
    logK <- barite(TK)
  else if (mineral == "celestite")
    logK <- celestite(TK)
  else logK <- ("Wrong input: The mineral should be ")

  return (logK)
}

# Henry's constants
#' data frame calculation
#'
#' @param Temperature temperature condition. default: 25 degree C
#' @param gas the gas species. default: "N2"
#'
#' @return KH constant of a gas at given temperature
#' @export
#'
#' @examples
#' Henry_constant(Temperature = 100, gas = "Ar")
#' The output is: 64421.35
#' @description calculate the Henry's constants

Henry_constant <- function (Temperature = 25, gas = 'N2') {

  Henry_Ar <- function (Temperature) {
    Tcl <- 647
    T0 <- (Temperature+273)/Tcl
    a0 <- 4.289125
    a1 <- 19.225988
    a2 <- -21.603721
    ln_KH <- a0+a1*(1-T0)^(1/3)/T0^2+a2*(1-T0)^(2/3)/T0^2
    KH <- exp(ln_KH)
    return (KH)
  }

  Henry_N2 <- function (Temperature) {
    pB <- 9.755-1121/(Temperature+273)-0.009169*(Temperature+273)
    Hk <- 10^(-pB)/55.55
    KH <- 1/Hk
    return (KH)
  }

  N2_Ar <- function (Temperature) { #N2/Ar ratio at given temperature
    Ar_air <- 9.32
    N2_air <- 781
    ratio <- N2_air/Ar_air*Henry_Ar(Temperature)/Henry_N2(Temperature)
    return (ratio)
  }

  if (gas=="N2") return(Henry_N2(Temperature))
  else if (gas=="Ar") return(Henry_Ar(Temperature))
  else if (gas=="N2/Ar") return(N2_Ar(Temperature))
  else print("Wrong input")

}

# CO2-CH4 equilibration temperature
#' data frame calculation
#'
#' @param Temperature temperature condition. default: 25 degree C
#'
#'
#' @return fractionation of δ13C between CO2 and CH4
#' @export
#'
#' @examples
#' Horita_T (Temperature = 50)
#' The output is: 61.45526
#' @description calculate stable carbon fractionation between CO2 and CH4

Horita_T <- function (Temperature = 25) { #input temperature in degree C, return dCO2-CH4

  TK <- Temperature + 273.15

  dCO2_CH4 <- 0.16+11.754*10^6/TK^2-2.3655*10^9/TK^3+0.2054*10^12/TK^4
  #the delta (CO2-CH4) if isotopic equilibrium achieved, Horita, 2001
  dd <- data.frame(Temperature=Temperature,delta=dCO2_CH4)
  #made a dataframe for plotting
  return (dCO2_CH4)
}

# CO2-CH4 equilibration to calculate temperature
#' data frame calculation
#'
#' @param dC ΔCO2-CH4
#'
#'
#' @return fractionation of δ13C between CO2 and CH4
#' @export
#'
#' @examples
#' Horita (dC = 50)
#' The output is: 212 with  unit of degree C
#' @description calculate the equilibration temperature based on the stable carbon fractionation between CO2 and CH4

Horita <- function (dC) {                     #input dCO2-CH4, return temperature in Celcius
  df <- function (temperature) {              #input temperature in K, return dCO2-CH4
    dCO2_CH4 <- 0.16+11.754*10^6/(temperature+273.15)^2
    -2.3655*10^9/(temperature+273.15)^3
    +0.2054*10^12/(temperature+273.15)^4         #the delta (CO2-CH4) if isotopic equilibrium
    #achieved, Horita, 2001
    dd <- data.frame(Temp=temperature,delta=dCO2_CH4)     #made a dataframe for plotting
    return (dd)
  }

  temperature <- seq(1,1000,1)
  df <- df(temperature)

  if (df$delta[1]<dC) {
    return("Warning: the input is larger than the maximum theoretical value!")
  } else if (df$delta[1000]>dC) {
    return("Warning: the input is smaller than the minimum theoretical value!")
  } else {
    n <- which.min(abs(df$delta-dC))             #return the nearest delta value at temperature
    return (n)
  }

}

# carbonate veins
#' data frame calculation
#'
#' @param mineral the mineral type of the carbonate veins.
#' @param d18O the d18O value of the carbonate veins. default: 20
#' @param ref the d18O of precipitation fluid. default: 0
#'
#'
#' @return Precipitation temperature of carbonate veins
#' @export
#'
#' @examples
#' d18O_temperature(d18O = 20, mineral = "dolomite", ref = -1)
#' The output is: 88.6529 with an unit of degree C
#' @description calculate the precipitation temperature of carbonate veins

d18O_temperature <- function (d18O = 20, mineral, ref = 0) {

  calcite <- function (d18O) { #input the d18O of calcite veins, return the precipitation temperture in degree C
    temp <- (d18O+2.89-ref)/2.78
    temperature <- sqrt(10^6/temp)-273.15 #Friedman and O'Neil, 1977, Kim and O'Neil, 1997
    return (temperature)
  }

  aragonite <- function (d18O) { #Kim et al., 2007
    17.88*10^3/(d18O+31.14-ref)-273.15 #input d18O of aragonite, return temperature
  }

  dolomite <- function (d18O) { #Horita, 2014
    temp <- (d18O+3.14-ref)/3.16
    temperature <- sqrt(10^6/temp)-273.15
    return (temperature)
  }

  if (mineral == "calcite")
    calculated_temperature <- calcite(d18O)
  else if (mineral == "aragonite")
    calculated_temperature <- aragonite(d18O)
  else if (mineral =="dolomite")
    calculated_temperature <- dolomite(d18O)
  else calculated_temperature <- ("Wrong input: the mineral should be either calcite, aragonite, or dolomite!")

  return (calculated_temperature)
}

# Normalize Rare Earth Elements (REEs) data
#' data frame calculation
#'
#' @param data the raw data that need to be normalized
#' @param method the normalization reference. default: "CI"
#'
#'
#' @return a data frame with REEs data normalized by CI chondrite or PAAS
#' @export
#'
#' @examples
#' REE_normalization (data = temp, method = "PAAS")
#' @description Normalize the REEs data

REE_normalization <- function (data, method = "CI") {
  #method should be either CI or PAAS

  #create a temporary data frame exclude the REEs

  data_temp <- data %>% select(-La,-Ce,-Pr,-Nd,-Sm,-Eu,
                               -Gd,-Tb,-Dy,-Y,-Ho,-Er,
                               -Tm,-Yb,-Lu)

  if (method=="CI") { #CI: CI chondrite
    data$La=data$La/0.237
    data$Ce=data$Ce/0.613
    data$Pr=data$Pr/0.0928
    data$Nd=data$Nd/0.457
    data$Sm=data$Sm/0.148
    data$Eu=data$Eu/0.0563
    data$Gd=data$Gd/0.199
    data$Tb=data$Tb/0.0361
    data$Dy=data$Dy/0.246
    data$Y=data$Y/1.57
    data$Ho=data$Ho/0.0546
    data$Er=data$Er/0.16
    data$Tm=data$Tm/0.0247
    data$Yb=data$Yb/0.161
    data$Lu=data$Lu/0.0246
    data_temp2 <- data.frame(La=data$La,Ce=data$Ce,Pr=data$Pr,Nd=data$Nd,
                             Sm=data$Sm,Eu=data$Eu,Gd=data$Gd,Tb=data$Tb,
                             Dy=data$Dy,Y=data$Y,Ho=data$Ho,Er=data$Er,
                             Tm=data$Tm,Yb=data$Yb,Lu=data$Lu)
    data_normalized <- bind_cols(data_temp,data_temp2)
    return(data_normalized)

  } else if (method=="PAAS") { #PAAS
    data$La=data$La/38.2
    data$Ce=data$Ce/79.6
    data$Pr=data$Pr/8.83
    data$Nd=data$Nd/33.9
    data$Sm=data$Sm/5.55
    data$Eu=data$Eu/1.08
    data$Gd=data$Gd/4.66
    data$Tb=data$Tb/0.774
    data$Dy=data$Dy/4.68
    data$Y=data$Y/27
    data$Ho=data$Ho/0.991
    data$Er=data$Er/2.85
    data$Tm=data$Tm/0.405
    data$Yb=data$Yb/2.82
    data$Lu=data$Lu/0.433
    data_temp2 <- data.frame(La=data$La,Ce=data$Ce,Pr=data$Pr,Nd=data$Nd,
                             Sm=data$Sm,Eu=data$Eu,Gd=data$Gd,Tb=data$Tb,
                             Dy=data$Dy,Y=data$Y,Ho=data$Ho,Er=data$Er,
                             Tm=data$Tm,Yb=data$Yb,Lu=data$Lu)
    data_normalized <- bind_cols(data_temp,data_temp2)
    return(data_normalized)

  } else return("Please input the correct normalization method, CI or PAAS")
}

# CO2 partial dissolution
#' data frame calculation
#'
#' @param Temperature the temperature of fluid. default: 25 degree C
#' @param pH the pH value of the fluid. default: 7
#'
#'
#' @return molar ratio between CO2 dissolved in aqueous solution and in gas phase
#' @export
#'
#' @examples
#' CO2_dissolution(Temperature = 100, pH = 8)
#' The output is: 1.61824, indicating that the molar ratio between dissolved CO2 and gas phase CO2 is 1.62
#' @description calculate molar ratio between dissolved CO2 and in gas phase

CO2_dissolution <- function (Temperature = 25, pH = 7) { #input the Temperature in degree C
  ka1 <- 4.45*10^-7
  ka2 <- 4.67*10^-11
  Hi <- 10^(-pH)
  vant1 <- 2400*(1/(Temperature+273.15)-1/298.15)
  vant2 <- exp(vant1)
  CO2_pH <- 1+ka1*vant2/Hi+ka1*ka2*vant2*vant2/(Hi^2)
  CO2_aq <- 0.83*vant2
  total_CO2 <- CO2_pH*CO2_aq
  return (total_CO2) #ratios between CO~2~ in aqueous solution and in gas phase.
}
