% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/main.R
\name{CO2_dissolution}
\alias{CO2_dissolution}
\title{CO2 partial dissolution}
\usage{
CO2_dissolution(Temperature = 25, pH = 7)
}
\arguments{
\item{Temperature}{the temperature of fluid. default: 25 degree C}

\item{pH}{the pH value of the fluid. default: 7}
}
\value{
molar ratio between CO2 dissolved in aqueous solution and in gas phase
}
\description{
calculate molar ratio between dissolved CO2 and in gas phase
}
\examples{
CO2_dissolution(Temperature = 100, pH = 8)
The output is: 1.61824, indicating that the molar ratio between dissolved CO2 and gas phase CO2 is 1.62
}
