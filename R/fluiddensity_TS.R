#' Seawater Density Calculator
#'
#' @description
#' Calculates fluid density as a function of temperature in Celsius and salinity in g/kg, using equations developed for seawater
#'
#' @details
#' This function contains two empirical equations, one from:
#'  F. J. Millero, and A. Poisson, International one-atmosphere equation of state of seawater. Deep-Sea Research, 28A (6), 625 – 629, 1981.
#'  and another from:
#'  Sharqawy, M. H., Lienhard, J. H., & Zubair, S. M. (2010). Thermophysical properties of seawater: a review of existing correlations and data. Desalination and Water Treatment, 16(1–3), 354–380.
#'
#' @param tempC temperature in degrees Celsius; valid for -2 to 40 degrees C ("low" range from Millero & Poisson) or 0 to 180 degrees C ("high" range from Sharqawy et al)
#' @param sal practical salinity in g/kg; valid for 0 to 42 g/kg ("low" range from Millero & Poisson) or 0 to 160 g/kg ("high" range from Sharqawy et al)
#' @param range either "low" for Millero & Poisson or "high" for Sharqawy et al
#' @returns density in kg/m^3
#' @examples
#' rho_sw1 <- seawaterdensity(25,35)
#' rho_sw2 <- seawaterdensity(c(25,30,35),35)
#'
#' @export

fluiddensity_TS <- function(tempC,sal,range = "low") {

  if (range == "low") {
    A <- 0.824496 - 4.0899*10^-3*tempC + 7.6438*10^-5*tempC^2 - 8.2467*10^-7*tempC^2 + 5.3875*10^-9*tempC^4
    B <- -5.72466*10^-3 + 1.0227*10^-4*tempC - 1.6546*10^-6*tempC^2
    C <- 4.8314*10^-4
    rho_w <- 999.842594 + 6.793952*10^-2*tempC - 9.09529*10^-3*tempC^2 + 1.001685*10^-4*tempC^3 - 1.120083*10^-6*tempC^4 + 6.536336*10^-9*tempC^5
    rho_sw <- rho_w + A*sal + B*sal^1.5 + C*sal^2
  } else if (range == "high") {
    a1 <- 9.999*10^2
    a2 <- 2.034*10^-2
    a3 <- -6.162*10^-3
    a4 <- 2.261*10^-5
    a5 <- -4.657*10^-8
    b1 <- 8.020*10^2
    b2 <- -2.001
    b3 <- 1.677*10^-2
    b4 <- -3.060*10^-5
    b5 <- -1.613*10^-5

    sal <- sal/1000

    rho_sw <- a1 + a2*tempC + a3*tempC^2 + a4*tempC^3 + a5*tempC^4 + b1*sal + b2*sal*tempC + b3*sal*tempC^2 + b4*sal*tempC^3 + b5*sal^2*tempC^2
  }

  return(rho_sw)
}
