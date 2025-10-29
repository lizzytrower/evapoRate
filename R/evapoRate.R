#' evapoRate: An Evaporation Modeling Function
#'
#' @description
#' Predicts aqueous solution change as a result of evaporation, including CO<sub>2</sub> exchange with the atmosphere and CaCO<sub>3</sub> precipitation
#'
#' @details
#' This function is an implementation of PHREEQC that can be used to predict the changes in major ion concentrations of a solution as the result of evaporation (i.e., loss of H2O). The function also includes the ability to model the effects of CO<sub>2</sub> gas exchange (invasion or degassing) on the evolving solution chemistry, and the C isotopic evolution of dissolved inorganic carbon (DIC). evapoRate includes 4 modes for CO<sub>2</sub> gas exchange:
#' (1) no gas exchange
#' (2) gas exchange based on the "stagnant film" diffusion model
#' (3) gas exchange based on the "solid wall" model
#' (4) gas exchange parameterized as a function of wind speed
#' Modes 2 and 3 are most appropriate for application to indoor evaporation experiments; mode 4 is most appropriate for application to the environment.
#'
#' This function also includes the ability to model the effects of CaCO<sub>3</sub> precipitation, allowing the user to provide a precipitation flux and set the ratios of Ca, Mg, and Alkalinity consumed per mole of CaCO<sub>3</sub> precipitated.
#'
#' Finally, this function can also model changes in \eqn{\delta}<sup>13</sup>C values of dissolved inorganic carbon (DIC) driven by gas exchange. Note that changes to \eqn{\delta}<sup>13</sup>C<sub>DIC</sub> values related to precipitation is not currently enabled, but may be included in a future version of the model.
#'
#' Similarly, precipitation of other evaporite minerals is not currently enabled, but could be added in future model versions, or could be implemented by the user via changing concentrations of relevant elements between timesteps.
#'
#' Overall, this function is designed to be run within a for-loop over a series of timesteps through which water is removed incrementally.
#'
#' @param tempC temperature in degrees Celsius; default is 25 C
#' @param P_atm atmospheric pressure in atm; default is 1 atm
#' @param xCO2 mole fraction of CO<sub>2</sub> in ppm; default is 420 ppm
#' @param water_initial_mass initial mass of water in kg; default is 1 kg
#' @param water_removed_mass amount of water removed by evaporation in kg
#' @param Na_i initial concentration of sodium (Na) in mmol/kg, default is 0
#' @param Mg_i initial concentration of magnesium (Mg) in mmol/kg, default is 0
#' @param Ca_i initial concentration of calcium (Ca) in mmol/kg, default is 0
#' @param K_i initial concentration of potassium (K) in mmol/kg, default is 0
#' @param Cl_i initial concentration of chloride (Cl) in mmol/kg, default is 0
#' @param SO4_i initial concentration of sulfate (SO4) in mmol/kg, default is 0
#' @param Si_i initial concentration of dissolved silica (Si) in mmol/kg, default is 0
#' @param Alk_i initial alkalinity in mmol/kg, default is 0
#' @param DIC_i initial concentration of dissolved inorganic carbon (DIC) in mmol/kg, default is 0
#' @param gas_exchange_mode one of 4 possible text strings selecting the option for CO<sub>2</sub> gas exchange mode; the default option is "none" - no CO<sub>2</sub> exchange; the other options are "stagnant_film", "solid_wall", and "wind"
#' @param timestep the timestep duration in hr, default is 1 hr; in theory this should be related to the amount of water removed via an observed or estimated evaporation rate
#' @param surface_area a surface area over which evaporation is occurring, in m<sup>2</sup>; required for any gas exchange mode other than "none"; default is 1 m<sup>2</sup>
#' @param z_film stagnant film thickness, in m; required for "stagnant_film" gas exchange mode but unused for other gas exchange modes
#' @param rho_air air density in kg/m<sup>3</sup>; required for "solid_wall" gas exchange mode but unused for other gas exchange modes
#' @param u_air indoor air velocity in m/s; required for "solid_wall" gas exchange mode but unused for other gas exchange modes
#' @param surface_length length of surface over which evaporation is occurring; required for "solid_wall" gas exchange mode but unused for other gas exchange modes
#' @param alpha_exchange a dimensionless factor describing the enhancement of CO<sub>2</sub> exchange due to CO<sub>2</sub> hydration or dehydration reactions, where a value of 1 means no enhancement
#' @param u_wind wind velocity in m/s; required for "wind" gas exchange mode bu unused for other gas exchange modes
#' @param d13_C_tracking one of two possible strings determining whether to track carbon isotopes; default is "no" (\eqn{\delta}<sup>13</sup>C tracking is off); the other option is "yes" (\eqn{\delta}<sup>13</sup>C tracking is on). \eqn{\delta}<sup>13</sup>C tracking is only available when gas exchange is turned on.
#' @param d13C_i_invasion initial \eqn{\delta}<sup>13</sup>C<sub>DIC</sub> value for gas invasion (value from previous timestep)
#' @param d13C_i_Rayleigh initial \eqn{\delta}<sup>13</sup>C<sub>DIC</sub> value for degassing (value when degassing starts)
#' @param eps_g_DIC carbon isotope fractionation (permil) associated with CO<sub>2</sub> invasion
#' @param alpha_DIC_g Rayleigh carbon isotope fractionation factor associated with CO<sub>2</sub degassing >
#' @param f_remaining fraction of solution remaining - value between 1 and 0, needed for predicted carbon isotope evolution during CO<sub>2</sub> degassing
#' @param CaCO3_precip one of two possible strings determining whether to include CaCO<sub>3</sub> precipitation; default is "no" (no CaCO<sub>3</sub> precipitation); "yes" turns on CaCO<sub>3</sub> precipitation
#' @param Fcarb_mol_kg CaCO<sub>3</sub> precipitation flux, in mol/kg; note that this should be scaled for the timestep duration
#' @param Alk_Ca_carb ratio of alkalinity consumed per mole of CaCO<sub>3</sub> precipitated; default is 2
#' @param Ca_Ca_carb ratio of calcium consumed per mole of CaCO<sub>3</sub> precipitated; default is 1 (note this might be appropriate to change for high-Mg calcite)
#' @param Mg_Ca_carb ratio of magnesium consumed per mole of CaCO<sub>3</sub> precipitated; default is 0 (note this might be appropriate to change for high-Mg calcite)
#' @param MgCO3_precip one of two possible strings determining whether to include MgCO<sub>3</sub> precipitation; default is "no" (no MgCO<sub>3</sub> precipitation); "yes" turns on MgCO<sub>3</sub> precipitation
#' @param FMgcarb_mol_kg MgCO<sub>3</sub> precipitation flux, in mol/kg; note that this should be scaled for the timestep duration
#' @param Alk_Mgcarb ratio of alkalinity consumed per mole of MgCO<sub>3</sub> precipitated; default is 10
#' @param DIC_Mgcarb ratio of DIC consumed per mole of MgCO<sub>3</sub> precipitated; default is 4
#' @param Mg_Mgcarb ratio of magnesium consumed per mole of MgCO<sub>3</sub> precipitated; default is 5
#' @param SO4_precip one of two possible strings determining whether to include SO<sub>4</sub> precipitation; default is "no" (no SO<sub>4</sub> precipitation); "yes" turns on SO<sub>4</sub> precipitation
#' @param FSO4_mol_kg SO<sub>4</sub> precipitation flux, in mol/kg; note that this should be scaled for the timestep duration
#' @param Ca_SO4 ratio of calcium consumed per mole of SO<sub>4</sub> precipitated; default is 1 (can set to zero to simulate sulfate diffusion into sediment or microbial sulfate reduction)
#' @param SO4_SO4 ratio of sulfate consumed per mole of SO<sub>4</sub> precipitated; default is 1
#' @param Mg_SO4 ratio of magnesium consumed per mole of SO<sub>4</sub> precipitated; default is 1 (can set to zero to simulate sulfate diffusion into sediment or microbial sulfate reduction)
#' @param MgSi_precip one of two possible strings determining whether to include MgSi precipitation; default is "no" (no MgSi precipitation); "yes" turns on MgSi precipitation
#' @param FMgSi_mol_kg MgSi precipitation flux, in mol/kg; note that this should be scaled for the timestep duration
#' @param Mg_MgSi ratio of magnesium consumed per mole of MgSi precipitated; default is 4
#' @param Si_MgSi ratio of silicon consumed per mole of MgSi precipitated; default is 6
#' @param Alk_MgSi ratio of alkalinity consumed per mole of MgSi precipitated; default is 8
#'
#' @returns output dataframe that includes major ion concentrations, saturation states for a number of relevant minerals, and the calculated gas flux
#' @examples
#' output <- evapoRate(water_removed_mass = 0.01, Na_i = 1, Cl_i = 1)
#'
#' @export

evapoRate <- function(
    tempC = 25, #temperature in degrees C
    P_atm = 1, #atmospheric pressure in atm
    xCO2 = 420, #mole fraction of CO2 in ppm
    water_initial_mass = 1, #initial water mass in kg
    water_removed_mass, #amount of water removed in kg
    Na_i = 0, #initial [Na] (mmol/kg)
    Mg_i = 0, #initial [Mg] (mmol/kg)
    Ca_i = 0, #initial [Ca] (mmol/kg)
    K_i = 0, #initial [K] (mmol/kg)
    Cl_i = 0, #initial [Cl] (mmol/kg)
    SO4_i = 0, #initial [SO4] (mmol/kg)
    Si_i = 0, #initial [dSi] (mmol/kg)
    Alk_i = 0, #initial alkalinity (mmol/kg)
    DIC_i = 0, #initial [DIC] (mmol/kg)
    gas_exchange_mode = "none",
    timestep = 1, #timestep duration (hr), if using any gas exchange mode other than "none"
    surface_area = 1, #surface area over which evaporation is occurring (m^2), for any gas exchange mode other than "none"
    z_film = 300*10^-6, #stagnant film thickness, if using the "stagnant_film" gas exchange mode (m)
    rho_air = 0.984, #air density, if using the "solid_wall" gas exchange mode (kg/m^3)
    u_air = 0.2, #indoor air velocity, if using the "solid_wall" gas exchange mode (m/s)
    surface_length = 0, #length of surface over which evaporation is occurring, if using the "solid_wall" gas exchange mode (m)
    alpha_exchange = 1, #a factor describing the enhancement of CO2 exchange due to CO2 hydration or dehydration reactions, where a value of 1 means no enhancement
    u_wind = 0, #wind speed, if using the "wind" gas exchange mode (m/s)
    d13C_tracking = "no", #text string, either "yes" or "no" that directs whether or not to track d13C_DIC change
    d13C_DIC_i_invasion = 0, #initial d13C_DIC value for gas invasion (value from previous timestep)
    d13C_DIC_i_Rayleigh = 0, #initial d13C_DIC value for degassing (value when degassing starts)
    eps_g_DIC = 0,
    alpha_DIC_g = 1,
    f_remaining = 1, #fraction of solution remaining
    CaCO3_precip = "no",
    Fcarb_mol_kg = 0, #calculated CaCO3 precipitation flux (mol), where Fcarb > 0 corresponds to net precipitation
    Alk_Ca_carb = 2, #ratio of moles of alkalinity consumed to moles of CaCO3 produced
    Mg_Ca_carb = 0, #ratio of moles of Mg consumed to moles of CaCO3 produced (helpful for exploring precipitation of high-Mg calcite or protodolomite)
    Ca_Ca_carb = 1,
    MgCO3_precip = "no",
    FMgcarb_mol_kg = 0, #calculated MgCO3 precipitation flux (mol)
    Alk_Mgcarb = 10, #ratio of moles of alkalinity consumed to moles of MgCO3 produced
    DIC_Mgcarb = 4, #ratio of moles of DIC consumed to moles of MgCO3 produced
    Mg_Mgcarb = 5, #ratio of moles of Mg consumed to moles of MgCO3 produced
    SO4_precip = "no",
    FSO4_mol_kg = 0, #calculated sulfate mineral precipitation flux (mol)
    Ca_SO4 = 1, #ratio of moles of Ca consumed to moles of sulfate mineral produced (can set to zero to simulate sulfate diffusion into sediment or microbial sulfate reduction)
    SO4_SO4 = 1, #ratio of moles of SO4 consumed per moles of sulfate mineral produced
    Mg_SO4 = 0, #ratio of moles of Mg consumed to moles of sulfate mineral produced (can set to zero to simulate sulfate diffusion into sediment or microbial sulfate reduction)
    MgSi_precip = "no",
    FMgSi_mol_kg = 0, #calculated Mg-silicate precipitation flux (mol), where Fcarb > 0 corresponds to net precipitation
    Mg_MgSi = 4, #ratio of moles of Mg consumed to moles of Mg-silicate produced
    Si_MgSi = 6, #ratio of moles of Si consumed to moles of Mg-silicate produced
    Alk_MgSi = 8 #ratio of moles of Alk consumed to moles of Mg-silicate produced
) {

  #warnings to catch incorrect inputs
  if (water_removed_mass > water_initial_mass) {
    stop("Warning: mass of water removed cannot be greater than original water mass.")
  }

  if (water_removed_mass < 0) {
    stop("Warning: mass of water removed must be a positive number.")
  }

  if (Na_i < 0 || Mg_i < 0 || Ca_i < 0 || K_i < 0 || Cl_i < 0 || SO4_i < 0 || Si_i < 0 || Alk_i < 0 || DIC_i < 0) {
    stop("Warning: solution concentration inputs cannot be negative.")
  }

  #warnings for combinations of modes that aren't available yet
  if (d13C_tracking == "yes" && CaCO3_precip == "yes") {
    print("Warning: d13C tracking currently not enabled in combination with CaCO3 precipitation. Turning off d13C tracking mode.")
    d13C_tracking <- "no"
  }

  if (CaCO3_precip == "yes" && Fcarb_mol_kg < 0) {
    print("Warning: CaCO3 dissolution kinetics not available yet. Setting Fcarb to zero.")
    Fcarb_mol_kg <- 0
  }

  if (gas_exchange_mode == "none" && CaCO3_precip == "yes") {
    print("Warning: CaCO3 precipitation not available when gas exchange mode is set to none. Turning off CaCO3 precipitation.")
    CaCO3_precip <- "no"
  }

  if (gas_exchange_mode == "none" && d13C_tracking == "yes") {
    print("Warning: d13C tracking not available when gas exchange mode is set to none. Turning off d13C tracking mode.")
    d13C_tracking <- "no"
  }

  water_removed_mol <- water_removed_mass/18.01528*1000
  pCO2_obs <- xCO2*P_atm #{uatm}

  phreeqc::phrLoadDatabaseString(pitzer.dat)
  phreeqc::phrSetOutputStringsOn(TRUE)

  input1 <- c(
    ' PHASES',
    ' Hydromagnesite',
    '      Mg5(CO3)4(OH)2:4H2O +6.0000 H+  =  + 4.0000 HCO3- + 5.0000 Mg++ + 6.0000 H2O',
    '      log_k           30.8539',
    ' -delta_H	-289.696	kJ/mol',	# Calculated enthalpy of reaction
    ' Monohydrocalcite',
    '    CaCO3:H2O +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 HCO3-',
    '      log_k           2.6824',
    ' -delta_H	-20.5648	kJ/mol',	# Calculated enthalpy of reaction

    '  SOLUTION 1',
    '  units         mmol/kgw',
    paste('  temp              ',as.character(tempC)),
    paste('  pressure          ',as.character(P_atm)),
    paste('  Alkalinity        ',as.character(Alk_i)),
    paste('  C(4)              ',as.character(DIC_i)),
    paste('  Ca                ',as.character(Ca_i)),
    paste('  Mg                ',as.character(Mg_i)),
    paste('  Na                ',as.character(Na_i)),
    paste('  K                 ',as.character(K_i)),
    paste('  Cl                ',as.character(Cl_i)),
    paste('  S(6)              ',as.character(SO4_i)),
    paste('  Si                ',as.character(Si_i)),
    paste('  water             ',as.character(water_initial_mass)),

    'SELECTED_OUTPUT       ',
    '  -high precision   TRUE',
    '  -pH               TRUE',
    '  -Alkalinity       TRUE',
    '  -molalities       CO2',
    '  -totals           C(4) Na Ca Mg K Cl S(6) Si',
    '  -si               aragonite calcite dolomite magnesite monohydrocalcite hydromagnesite gypsum sepiolite sylvite hexahydrite epsomite CO2(g)',
    '  -activities       H2O',

    ' REACTION 1',
    '    H2O     -1.0',
    paste(as.character(water_removed_mol), 'moles'),
    ' SAVE solution 2',
    ' END',

    ' EQUILIBRIUM_PHASES',
    paste( ' CO2(g) ',as.character(round(log10(pCO2_obs*10^-6),digits = 3)),' 100'),
    ' USE solution 2',
    ' END')

  phreeqc::phrRunString(input1)
  phreeqc_output1 <- phreeqc::phrGetSelectedOutput()

  output <- as.data.frame(pCO2_obs)

  #save the conservative ions and water activity
  output$Na_mmol_kg <- phreeqc_output1$n1$Na.mol.kgw.[2]*1000
  output$K_mmol_kg <- phreeqc_output1$n1$K.mol.kgw.[2]*1000
  output$Cl_mmol_kg <- phreeqc_output1$n1$Cl.mol.kgw.[2]*1000
  output$salinity <- phreeqc_output1$n1$Na.mol.kgw.[2]*22.989769 + phreeqc_output1$n1$Ca.mol.kgw.[2]*40.078 + phreeqc_output1$n1$Mg.mol.kgw.[2]*24.305 + phreeqc_output1$n1$K.mol.kgw.[2]*39.0983 + phreeqc_output1$n1$Cl.mol.kgw.[2]*35.453 + phreeqc_output1$n1$S.6..mol.kgw.[2]*96.06 + phreeqc_output1$n1$C.4..mol.kgw.[2]*61.0168
  output$H2O_activity <- 10^phreeqc_output1$n1$la_H2O[2]

  #CO2 exchange
  if (gas_exchange_mode == "none") {
    #if there is no gas exchange, we just save the outputs of the initial calculation
    output$pH <- phreeqc_output1$n1$pH[2]
    output$DIC_mmol_kg <- phreeqc_output1$n1$C.4..mol.kgw.[2]*1000
    output$Alk_mmol_kg <- phreeqc_output1$n1$Alk.eq.kgw.[2]*1000
    output$pCO2 <- 10^phreeqc_output1$n1$si_CO2.g.[2]*10^6
    output$SO4_mmol_kg <- phreeqc_output1$n1$S.6..mol.kgw.[2]*1000
    output$Ca_mmol_kg <- phreeqc_output1$n1$Ca.mol.kgw.[2]*1000
    output$Mg_mmol_kg <- phreeqc_output1$n1$Mg.mol.kgw.[2]*1000
    output$Si_mmol_kg <- phreeqc_output1$n1$Si.mol.kgw.[2]*1000
    output$Omega_ar <- 10^phreeqc_output1$n1$si_aragonite[2]
    output$Omega_cc <- 10^phreeqc_output1$n1$si_calcite[2]
    output$Omega_dol <- 10^phreeqc_output1$n1$si_dolomite[2]
    output$Omega_mag <- 10^phreeqc_output1$n1$si_magnesite[2]
    output$Omega_mhc <- 10^phreeqc_output1$n1$si_monohydrocalcite[2]
    output$Omega_hmag <- 10^phreeqc_output1$n1$si_hydromagnesite[2]
    output$Omega_gy <- 10^phreeqc_output1$n1$si_gypsum[2]
    output$Omega_sep <- 10^phreeqc_output1$n1$si_sepiolite[2]
    output$Omega_hexahydrite <- 10^phreeqc_output2$n1$si_hexahydrite[2]
    output$Omega_epsomite <- 10^phreeqc_output2$n1$si_epsomite[2]
    output$Omega_sylvite <- 10^phreeqc_output2$n1$si_sylvite[2]
    output$Fgas <- 0

  } else {
    if (gas_exchange_mode == "stagnant_film") {

      mu <- seawaterdynamicviscosity(tempC,output$salinity) #{kg/m*s}
      kB <- 1.380649*10^-23 #{kg*m^2/K*s^2}
      nSE <- 4 #Stokes-Einstein number {dimensionless}
      a298 <- 168*10^-12 #{m}
      alpha <- 2*10^-3
      tempK <- tempC + 273
      a <- a298*(1 + alpha*(tempK - 298)) #{m}
      D_CO2 <- kB*tempK/(nSE*pi*mu*a) #{m^2/s}
      D_CO2 <- D_CO2*3600 #{m^2/hr}
      k_CO2 <- D_CO2/z_film

    } else if (gas_exchange_mode == "solid_wall") {

      mu <- seawaterdynamicviscosity(tempC,output$salinity) #{kg/m*s}
      kB <- 1.380649*10^-23 #{kg*m^2/K*s^2}
      nSE <- 4 #Stokes-Einstein number {dimensionless}
      a298 <- 168*10^-12 #{m}
      alpha <- 2*10^-3
      tempK <- tempC + 273
      a <- a298*(1 + alpha*(tempK - 298)) #{m}
      D_CO2 <- kB*tempK/(nSE*pi*mu*a) #{m^2/s}
      D_CO2 <- D_CO2*3600 #{m^2/hr}

      dens <- fluiddensity_TS(tempC,output$salinity)
      R_l <- u_air*surface_length/(mu/dens)
      u_star <- 0.477*(log10(R_l))^(-1.29)*u_air
      Sc <- (mu/dens)/(D_CO2/3600) #Schmidt number {dimensionless}
      k_CO2 <- 0.082*(rho_air/dens)^0.5*Sc^(-2/3)*u_star #{m/s}
      k_CO2 <- k_CO2*3600 #{m/hr}

    } else if (gas_exchange_mode == "wind") {

      mu <- seawaterdynamicviscosity(tempC,output$salinity) #{kg/m*s}
      kB <- 1.380649*10^-23 #{kg*m^2/K*s^2}
      nSE <- 4 #Stokes-Einstein number {dimensionless}
      a298 <- 168*10^-12 #{m}
      alpha <- 2*10^-3
      tempK <- tempC + 273
      a <- a298*(1 + alpha*(tempK - 298)) #{m}
      D_CO2 <- kB*tempK/(nSE*pi*mu*a) #{m^2/s}
      D_CO2 <- D_CO2*3600 #{m^2/hr}

      dens <- fluiddensity_TS(tempC,output$salinity)
      Sc <- (mu/dens)/(D_CO2/3600) #Schmidt number {dimensionless}
      k_CO2 <- 0.251*u_wind^2*(Sc/600)^0.5 #{cm/hr}
      k_CO2 <- k_CO2/100 #{m/hr}

    }

    #gas flux calculations
    dens <- fluiddensity_TS(tempC,output$salinity)
    KH_kg <- phreeqc_output1$n1$m_CO2.mol.kgw.[3]/(pCO2_obs*10^-6)
    KH_m3 <- KH_kg*dens #{mol/m^3/atm}
    delta_pCO2 <- 10^phreeqc_output1$n1$si_CO2.g.[2]*10^6 - pCO2_obs #{uatm}
    delta_pCO2 <- delta_pCO2*10^-6 #{atm}
    Fgas <- k_CO2*KH_m3*delta_pCO2*alpha_exchange #{mol/m^2/hr}
    Fgas <- Fgas*timestep*surface_area #{mol}

    water_final_mass <- water_initial_mass - water_removed_mass

    #CaCO3 precipitation calculations (Si and SO4 are in here as part of the framing for future expansion to include precipitation of additional minerals)
    if (Fcarb_mol_kg <= 0 || CaCO3_precip == "no") {
      Alk_new <- phreeqc_output1$n1$Alk.eq.kgw.[2]*1000
      DIC_new <- phreeqc_output1$n1$C.4..mol.kgw.[2]*1000
      Ca_new <- phreeqc_output1$n1$Ca.mol.kgw.[2]*1000
      Mg_new <- phreeqc_output1$n1$Mg.mol.kgw.[2]*1000
      SO4_new <- phreeqc_output1$n1$S.6..mol.kgw.[2]*1000
      Si_new <- phreeqc_output1$n1$Si.mol.kgw.[2]*1000
      Cl_new <- phreeqc_output1$n1$Cl.mol.kgw.[2]*1000
      K_new <- phreeqc_output1$n1$K.mol.kgw.[2]*1000
    } else if (CaCO3_precip == "yes") {
      Fcarb_mmol_kg <- Fcarb_mol_kg*10^3
      Alk_new <- phreeqc_output1$n1$Alk.eq.kgw.[2]*1000 - Alk_Ca_carb*Fcarb_mmol_kg
      DIC_new <- phreeqc_output1$n1$C.4..mol.kgw.[2]*1000 - Fcarb_mmol_kg
      Ca_new <- phreeqc_output1$n1$Ca.mol.kgw.[2]*1000 - Ca_Ca_carb*Fcarb_mmol_kg
      Mg_new <- phreeqc_output1$n1$Mg.mol.kgw.[2]*1000 - Mg_Ca_carb*Fcarb_mmol_kg
      SO4_new <- phreeqc_output1$n1$S.6..mol.kgw.[2]*1000
      Si_new <- phreeqc_output1$n1$Si.mol.kgw.[2]*1000
      Cl_new <- phreeqc_output1$n1$Cl.mol.kgw.[2]*1000
      K_new <- phreeqc_output1$n1$K.mol.kgw.[2]*1000
    }

    #MgCO3 precipitation calculations
    if (FMgcarb_mol_kg <= 0 || MgCO3_precip == "no") {
      Alk_new <- Alk_new
      DIC_new <- DIC_new
      Mg_new <- Mg_new
    } else {
      FMgcarb_mmol_kg <- FMgcarb_mol_kg*10^3
      Alk_new <- Alk_new - Alk_Mgcarb*FMgcarb_mmol_kg
      DIC_new <- DIC_new - DIC_Mgcarb*FMgcarb_mmol_kg
      Mg_new <- Mg_new - Mg_Mgcarb*FMgcarb_mmol_kg
    }

    #sulfate mineral precipitation calculations
    if (FSO4_mol_kg <= 0 || SO4_precip == "no") {
      Mg_new <- Mg_new
      SO4_new <- SO4_new
      Ca_new <- Ca_new
    } else {
      FSO4_mmol_kg <- FSO4_mol_kg*10^3
      Mg_new <- Mg_new - Mg_SO4*FSO4_mmol_kg
      Ca_new <- Ca_new - Ca_SO4*FSO4_mmol_kg
      SO4_new <- SO4_new - SO4_SO4*FSO4_mmol_kg
    }

    #Mg-silicate mineral precipitation calculations
    if (FMgSi_mol_kg <= 0 || MgSi_precip == "no") {
      Mg_new <- Mg_new
      Alk_new <- Alk_new
      Si_new <- Si_new
    } else {
      FMgSi_mmol_kg <- FMgSi_mol_kg*10^3
      Mg_new <- Mg_new - Mg_MgSi*FMgSi_mmol_kg
      Alk_new <- Alk_new - Alk_MgSi*FMgSi_mmol_kg
      Si_new <- Si_new - Si_MgSi*FMgSi_mmol_kg
    }

    #chloride mineral precipitation calculations
    if (FCl_mol_kg <= 0 || Cl_precip == "no") {
      K_new <- K_new
      Cl_new <- Cl_new
    } else {
      FCl_mmol_kg <- FCl_mol_kg*10^3
      K_new <- K_new - K_Cl*FCl_mmol_kg
      Cl_new <- Cl_new - Cl_Cl*FCl_mmol_kg
    }

    #second PHREEQC calculation to recalculate carbonate speciation incorporating gas exchange and CaCO3 precipitation
    input2 <- c(
      ' PHASES',
      ' Hydromagnesite',
      '      Mg5(CO3)4(OH)2:4H2O +6.0000 H+  =  + 4.0000 HCO3- + 5.0000 Mg++ + 6.0000 H2O',
      '      log_k           30.8539',
      ' -delta_H	-289.696	kJ/mol',	# Calculated enthalpy of reaction
      ' Monohydrocalcite',
      '    CaCO3:H2O +1.0000 H+  =  + 1.0000 Ca++ + 1.0000 H2O + 1.0000 HCO3-',
      '      log_k           2.6824',
      ' -delta_H	-20.5648	kJ/mol',	# Calculated enthalpy of reaction

      '  SOLUTION 2',
      '  units         mmol/kgw',
      paste('  temp              ',as.character(tempC)),
      paste('  pressure          ',as.character(P_atm)),
      paste('  Alkalinity        ',as.character(Alk_new)),
      paste('  C(4)              ',as.character(DIC_new)),
      paste('  Ca                ',as.character(Ca_new)),
      paste('  Mg                ',as.character(Mg_new)),
      paste('  Na                ',as.character(output$Na_mmol_kg)),
      paste('  K                 ',as.character(K_new)),
      paste('  Cl                ',as.character(Cl_new)),
      paste('  S(6)              ',as.character(SO4_new)),
      paste('  Si                ',as.character(Si_new)),
      paste('  water             ',as.character(water_final_mass)),

      'SELECTED_OUTPUT       ',
      '  -high precision   TRUE',
      '  -pH               TRUE',
      '  -Alkalinity       TRUE',
      '  -molalities       CO2',
      '  -activities       H+ Mg+2 H4SiO4 H2O',
      '  -totals           C(4) Na Ca Mg K Cl S(6) Si',
      '  -si               aragonite calcite dolomite magnesite monohydrocalcite hydromagnesite gypsum sepiolite sylvite hexahydrite epsomite CO2(g)',

      ' REACTION 1',
      '    CO2     -1.0',
      paste(as.character(Fgas), 'moles'),
      ' SAVE solution 2',
      ' END')

    phreeqc::phrRunString(input2)
    phreeqc_output2 <- phreeqc::phrGetSelectedOutput()

    #save carbonate chemistry outputs
    output$Ca_mmol_kg <- phreeqc_output2$n1$Ca.mol.kgw.[2]*1000
    output$Mg_mmol_kg <- phreeqc_output2$n1$Mg.mol.kgw.[2]*1000
    output$SO4_mmol_kg <- phreeqc_output2$n1$S.6..mol.kgw.[2]*1000
    output$Si_mmol_kg <- phreeqc_output2$n1$Si.mol.kgw.[2]*1000
    output$pH <- phreeqc_output2$n1$pH[2]
    output$DIC_mmol_kg <- phreeqc_output2$n1$C.4..mol.kgw.[2]*1000
    output$Alk_mmol_kg <- phreeqc_output2$n1$Alk.eq.kgw.[2]*1000
    output$pCO2 <- 10^phreeqc_output2$n1$si_CO2.g.[2]*10^6
    output$Omega_ar <- 10^phreeqc_output2$n1$si_aragonite[2]
    output$Omega_cc <- 10^phreeqc_output2$n1$si_calcite[2]
    output$Omega_dol <- 10^phreeqc_output2$n1$si_dolomite[2]
    output$Omega_mag <- 10^phreeqc_output2$n1 $si_magnesite[2]
    output$Omega_mhc <- 10^phreeqc_output2$n1$si_monohydrocalcite[2]
    output$Omega_hmag <- 10^phreeqc_output2$n1$si_hydromagnesite[2]
    output$Omega_gy <- 10^phreeqc_output2$n1$si_gypsum[2]
    output$Omega_sep <- 10^phreeqc_output2$n1$si_sepiolite[2]
    output$Omega_hexahydrite <- 10^phreeqc_output2$n1$si_hexahydrite[2]
    output$Omega_epsomite <- 10^phreeqc_output2$n1$si_epsomite[2]
    output$Omega_sylvite <- 10^phreeqc_output2$n1$si_sylvite[2]
    output$Fgas <- Fgas
    output$a_H <- 10^phreeqc_output2$n1$la_H.[2]
    output$a_Mg <- 10^phreeqc_output2$n1$la_Mg.2[2]
    output$a_Si <- 10^phreeqc_output2$n1$la_H4SiO4[2]
    output$a_H2O <- 10^phreeqc_output2$n1$la_H2O[2]


  if (d13C_tracking == "yes") {

    if (Fgas > 0) {
      eps_DIC_g <- (alpha_DIC_g - 1)*1000*log(f_remaining)
      output$d13C_DIC <- eps_DIC_g + d13C_DIC_i_Rayleigh
    } else {
        output$d13C_DIC <- (d13C_DIC_i_invasion*phreeqc_output1$n1$C.4..mol.kgw.[2]*water_final_mass - (-8.5 + eps_g_DIC)*Fgas)/(output$DIC_mmol_kg*water_final_mass*10^-3)
      }
    }
  }

  return(output)
}
