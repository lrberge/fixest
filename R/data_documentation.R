

#' Trade data sample
#'
#' This data reports trade information between countries of the European Union (EU15).
#'
#' @usage
#' data(trade, package = "fixest")
#'
#' @format
#' `trade` is a data frame with 38,325 observations and 6 variables named `Destination`, `Origin`, `Product`, `Year`, `dist_km` and `Euros`.
#'
#' * `Origin`: 2-digits codes of the countries of origin of the trade flow.
#' * `Destination`: 2-digits codes of the countries of destination of the trade flow.
#' * `Products`: Number representing the product categories (from 1 to 20).
#' * `Year`: Years from 2007 to 2016
#' * `dist_km`: Geographic distance in km between the centers of the countries of origin and destination.
#' * Euros: The total amount in euros of the trade flow for the specific year/product category/origin-destination country pair.
#'
#'
#' @source
#' This data has been extrated from Eurostat on October 2017.
#'
#'
#'
"trade"


#' Sample data for difference in difference
#'
#' This data has been generated to illustrate the use of difference in difference functions in 
#' package \pkg{fixest}. This is a balanced panel of 104 individuals and 10 periods. 
#' About half the individuals are treated, the treatment having a positive effect on 
#' the dependent variable `y` after the 5th period. The effect of the treatment on `y` is gradual.
#'
#' @usage
#' data(base_did, package = "fixest")
#'
#' @format
#' `base_did` is a data frame with 1,040 observations and 6 variables named 
#' `y`, `x1`, `id`, `period`, `post` and `treat`.
#'
#' \describe{
#' \item{y}{The dependent variable affected by the treatment.}
#' \item{x1}{ An explanatory variable.}
#' \item{id}{ Identifier of the individual.}
#' \item{period}{ From 1 to 10}
#' \item{post}{ Indicator taking value 1 if the period is strictly greater than 5, 0 otherwise.}
#' \item{treat}{ Indicator taking value 1 if the individual is treated, 0 otherwise.}
#'
#' }
#'
#' @source
#' This data has been generated from \pkg{R}.
#'
#'
#'
#'
"base_did"





#' Sample data for staggered difference in difference
#'
#' This data has been generated to illustrate the Sun and Abraham (Journal of Econometrics, 2021) method for staggered difference-in-difference. This is a balanced panel of 95 individuals and 10 periods. Half the individuals are treated. For those treated, the treatment date can vary from the second to the last period. The effect of the treatment depends on the time since the treatment: it is first negative and then increasing.
#'
#' @usage
#' data(base_stagg, package = "fixest")
#'
#' @format
#' `base_stagg` is a data frame with 950 observations and 7 variables:
#'
#' * id: panel identifier.
#' * year: from 1 to 10.
#' * year_treated: the period at which the individual is treated.
#' * time_to_treatment: different between the year and the treatment year.
#' * treated: indicator taking value 1 if the individual is treated, 0 otherwise.
#' * treatment_effect_true: true effect of the treatment.
#' * x1: explanatory variable, correlated with the period.
#' * y: the dependent variable affected by the treatment.
#'
#'
#' @source
#' This data has been generated from \pkg{R}.
#'
"base_stagg"


#' Publication data sample
#'
#' This data reports the publication output (number of articles and number of citations received) 
#' for a few scientists from the start of their career to 2000. 
#' Most of the variables are processed from the Microsoft Academic Graph (MAG) data set. A few variables are randomly generated.
#' 
#'
#' @usage
#' data(base_pub, package = "fixest")
#'
#' @format
#' `base_pub` is a data frame with 4,024 observations and 10 variables. There are 200 different scientists and 51 different years (ends in 2000).
#'
#' * `author_id`: scientist identifier
#' * `year`: current year
#' * `affil_id`: affiliation ID of the scientist's current affiliation
#' * `affil_name`: affiliation name of the scientist's current affiliation (character)
#' * `field`: field name of the scientist (character), time invariant
#' * `nb_pub`: number of publications of the scientist for the current year
#' * `nb_cites`: number of citations received by the publications of the scientist in the current year. Accounts for the citations received from articles published up to 2020.
#' * `birth_year`: birth year of the scientist (this is randomly generated)
#' * `is_woman`: 1 if the scientist is a woman, 0 otherwise (this is randomly generated)
#' * `age`: current age of the scientist (formally `year - birth_year`)
#'
#'
#' @source
#' The source of this data set is the Microsoft Academic Graph data set, extracted in 2020. Now a defunct project, you can find similar data on [OpenAlex](https://docs.openalex.org/).
#' 
#' The variables `birth_year`, `is_woman` and `age` were randomly generated. All other variables have created from the raw MAG files.
#'
#'
#'
"base_pub"

