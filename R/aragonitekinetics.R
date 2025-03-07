#' Aragonite Kinetics Interpolation Function
#'
#' @description
#' Calculates the precipitation rate parameters (k and n) for aragonite as a function of temperature.
#'
#' @details
#' This function interpolates between the experimental constrains of the rate constant (k) and reaction order (n) for aragonite precipitation as observed at 5,25, and 37 C in Burton and Walter (1987):
#'
#' Burton, E. A., & Walter, L. M. (1987). Relative precipitation rates of aragonite and Mg calcite from seawater: Temperature or carbonate ion control? Geology, 15, 111–114.
#' or as observed at 10, 25, and 40 C in Romanek et al. (2011):
#' Romanek, C. S., Morse, J. W., & Grossman, E. L. (2011). Aragonite Kinetics in Dilute Solutions. Aquatic Geochemistry, 17(4), 339.
#'
#' @param tempC temperature in degrees Celsius
#' @param reference string with the source of kinetics data - must be either "BurtonWalter" or "Romanek"
#' @returns a list with two items, k (umol/m<sup>2</sup>/hr) and n (unitless)
#' @examples
#' arag1 <- aragonitekinetics(25)
#' arag2 <- mapply(aragonitekinetics,c(25,26,27))
#'
#' @export

aragonitekinetics <- function(tempC,
                              reference = "Romanek") {

  if (reference == "BurtonWalter") {
    caln <- c(0.4, 1.7, 2.4)
    calk <- c(21.8, 40.6, 45.1)
    calT <- c(5, 25, 37)

    if (tempC > 37) {
      k <- 45.1
      n <- 2.4
    }

    if (tempC < 5) {
      k <- 21.8
      n <- 0.4
    }

    if (tempC >= 5 && tempC <= 37) {
      k <- pracma::interp1(calT, calk, tempC)
      n <- pracma::interp1(calT, caln, tempC)
    }
  } else if (reference == "Romanek") {
    caln <- c(1.74, 1.74, 1.52)
    calk <- c(10^1.27, 10^1.85, 10^2.55)
    calT <- c(10, 25, 40)

    if (tempC > 40) {
      k <- 10^2.55
      n <- 2.55
    }

    if (tempC < 10) {
      k <- 10^1.27
      n <- 1.74
    }

    if (tempC >= 10 && tempC <= 40) {
      k <- pracma::interp1(calT, calk, tempC)
      n <- pracma::interp1(calT, caln, tempC)
    }
  } else {
    k <- NA
    n <- NA
  }

  return(list(k = k, n = n))
}
