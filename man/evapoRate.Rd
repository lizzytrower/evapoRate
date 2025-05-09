% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/evapoRate.R
\name{evapoRate}
\alias{evapoRate}
\title{evapoRate: An Evaporation Modeling Function}
\usage{
evapoRate(
  tempC = 25,
  P_atm = 1,
  xCO2 = 420,
  water_initial_mass = 1,
  water_removed_mass,
  Na_i = 0,
  Mg_i = 0,
  Ca_i = 0,
  K_i = 0,
  Cl_i = 0,
  SO4_i = 0,
  Si_i = 0,
  Alk_i = 0,
  DIC_i = 0,
  gas_exchange_mode = "none",
  timestep = 1,
  surface_area = 1,
  z_film = 300 * 10^-6,
  rho_air = 0.984,
  u_air = 0.2,
  surface_length = 0,
  alpha_exchange = 1,
  u_wind = 0,
  d13C_tracking = "no",
  d13C_DIC_i_invasion = 0,
  d13C_DIC_i_Rayleigh = 0,
  eps_g_DIC = 0,
  alpha_DIC_g = 1,
  f_remaining = 1,
  CaCO3_precip = "no",
  Fcarb_mol_kg = 0,
  Alk_Ca_carb = 2,
  Mg_Ca_carb = 0,
  Ca_Ca_carb = 1
)
}
\arguments{
\item{tempC}{temperature in degrees Celsius; default is 25 C}

\item{P_atm}{atmospheric pressure in atm; default is 1 atm}

\item{xCO2}{mole fraction of CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} in ppm; default is 420 ppm}

\item{water_initial_mass}{initial mass of water in kg; default is 1 kg}

\item{water_removed_mass}{amount of water removed by evaporation in kg}

\item{Na_i}{initial concentration of sodium (Na) in mmol/kg, default is 0}

\item{Mg_i}{initial concentration of magnesium (Mg) in mmol/kg, default is 0}

\item{Ca_i}{initial concentration of calcium (Ca) in mmol/kg, default is 0}

\item{K_i}{initial concentration of potassium (K) in mmol/kg, default is 0}

\item{Cl_i}{initial concentration of chloride (Cl) in mmol/kg, default is 0}

\item{SO4_i}{initial concentration of sulfate (SO4) in mmol/kg, default is 0}

\item{Si_i}{initial concentration of dissolved silica (Si) in mmol/kg, default is 0}

\item{Alk_i}{initial alkalinity in mmol/kg, default is 0}

\item{DIC_i}{initial concentration of dissolved inorganic carbon (DIC) in mmol/kg, default is 0}

\item{gas_exchange_mode}{one of 4 possible text strings selecting the option for CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} gas exchange mode; the default option is "none" - no CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} exchange; the other options are "stagnant_film", "solid_wall", and "wind"}

\item{timestep}{the timestep duration in hr, default is 1 hr; in theory this should be related to the amount of water removed via an observed or estimated evaporation rate}

\item{surface_area}{a surface area over which evaporation is occurring, in m\if{html}{\out{<sup>}}2\if{html}{\out{</sup>}}; required for any gas exchange mode other than "none"; default is 1 m\if{html}{\out{<sup>}}2\if{html}{\out{</sup>}}}

\item{z_film}{stagnant film thickness, in m; required for "stagnant_film" gas exchange mode but unused for other gas exchange modes}

\item{rho_air}{air density in kg/m\if{html}{\out{<sup>}}3\if{html}{\out{</sup>}}; required for "solid_wall" gas exchange mode but unused for other gas exchange modes}

\item{u_air}{indoor air velocity in m/s; required for "solid_wall" gas exchange mode but unused for other gas exchange modes}

\item{surface_length}{length of surface over which evaporation is occurring; required for "solid_wall" gas exchange mode but unused for other gas exchange modes}

\item{alpha_exchange}{a dimensionless factor describing the enhancement of CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} exchange due to CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} hydration or dehydration reactions, where a value of 1 means no enhancement}

\item{u_wind}{wind velocity in m/s; required for "wind" gas exchange mode bu unused for other gas exchange modes}

\item{eps_g_DIC}{carbon isotope fractionation (permil) associated with CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} invasion}

\item{alpha_DIC_g}{Rayleigh carbon isotope fractionation factor associated with CO\if{html}{\out{<sub>}}2</sub degassing >}

\item{f_remaining}{fraction of solution remaining - value between 1 and 0, needed for predicted carbon isotope evolution during CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} degassing}

\item{CaCO3_precip}{one of two possible strings determining whether to include CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitation; default is "no" (no CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitation); "yes" turns on CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitation}

\item{Fcarb_mol_kg}{CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitation flux, in mol/kg; note that this should be scaled for the timestep duration}

\item{Alk_Ca_carb}{ratio of alkalinity consumed per mole of CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitated; default is 2}

\item{Mg_Ca_carb}{ratio of magnesium consumed per mole of CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitated; default is 0 (note this might be appropriate to change for high-Mg calcite)}

\item{Ca_Ca_carb}{ratio of calcium consumed per mole of CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitated; default is 1 (note this might be appropriate to change for high-Mg calcite)}

\item{d13_C_tracking}{one of two possible strings determining whether to track carbon isotopes; default is "no" (\eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C tracking is off); the other option is "yes" (\eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C tracking is on). \eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C tracking is only available when gas exchange is turned on.}

\item{d13C_i_invasion}{initial \eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C\if{html}{\out{<sub>}}DIC\if{html}{\out{</sub>}} value for gas invasion (value from previous timestep)}

\item{d13C_i_Rayleigh}{initial \eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C\if{html}{\out{<sub>}}DIC\if{html}{\out{</sub>}} value for degassing (value when degassing starts)}
}
\value{
output dataframe that includes major ion concentrations, saturation states for a number of relevant minerals, and the calculated gas flux
}
\description{
Predicts aqueous solution change as a result of evaporation, including CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} exchange with the atmosphere and CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitation
}
\details{
This function is an implementation of PHREEQC that can be used to predict the changes in major ion concentrations of a solution as the result of evaporation (i.e., loss of H2O). The function also includes the ability to model the effects of CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} gas exchange (invasion or degassing) on the evolving solution chemistry, and the C isotopic evolution of dissolved inorganic carbon (DIC). evapoRate includes 4 modes for CO\if{html}{\out{<sub>}}2\if{html}{\out{</sub>}} gas exchange:
(1) no gas exchange
(2) gas exchange based on the "stagnant film" diffusion model
(3) gas exchange based on the "solid wall" model
(4) gas exchange parameterized as a function of wind speed
Modes 2 and 3 are most appropriate for application to indoor evaporation experiments; mode 4 is most appropriate for application to the environment.

This function also includes the ability to model the effects of CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitation, allowing the user to provide a precipitation flux and set the ratios of Ca, Mg, and Alkalinity consumed per mole of CaCO\if{html}{\out{<sub>}}3\if{html}{\out{</sub>}} precipitated.

Finally, this function can also model changes in \eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C values of dissolved inorganic carbon (DIC) driven by gas exchange. Note that changes to \eqn{\delta}\if{html}{\out{<sup>}}13\if{html}{\out{</sup>}}C\if{html}{\out{<sub>}}DIC\if{html}{\out{</sub>}} values related to precipitation is not currently enabled, but may be included in a future version of the model.

Similarly, precipitation of other evaporite minerals is not currently enabled, but could be added in future model versions, or could be implemented by the user via changing concentrations of relevant elements between timesteps.

Overall, this function is designed to be run within a for-loop over a series of timesteps through which water is removed incrementally.
}
\examples{
output <- evapoRate(water_removed_mass = 0.01, Na_i = 1, Cl_i = 1)

}
