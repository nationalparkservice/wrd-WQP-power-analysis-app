#' Generate vector of yearly data
#'
#' @param df A data frame with columns for Value, Date, Year, Month, Day, Quarter
#' @return A numeric vector
#' @export
#'
#' @examples
#' yearly_subsample(data.frame(Value = rnorm(10, 2, 0.1), Year = sample(2000:2005, 10, replace = TRUE), Month = c(1:10), Day = sample(1:30, 10)) |> dplyr::mutate(Date = paste0(Year, "-", Month, "-", Day), Quarter = c(rep(c(1, 2, 3), each = 3), rep(4, 1))))

yearly_subsample <- function(df){
  ## yearly subsample
  year_vec <- c()

  for (i in 1:length(unique(df$Year))){
    year <- unique(df$Year)[i]
    tmp <- df |> dplyr::filter(Year == year) |>
      dplyr::slice(1)
    year_vec <- c(year_vec, tmp$Value)
  }

  year_vec <- as.numeric(year_vec)

  return(year_vec)
}


#' Generate vector of biannual data
#'
#' @param df A data frame with columns for Value, Date, Year, Month, Day, Quarter, Half
#' @return A numeric vector
#' @export
#'
#' @examples
#' biannually_subsample(data.frame(Value = rnorm(10, 2, 0.1), Year = sample(2000:2005, 10, replace = TRUE), Month = c(1:10), Day = sample(1:30, 10)) |> dplyr::mutate(Date = paste0(Year, "-", Month, "-", Day), Half = c(rep(1, 6), rep(2, 4)), Quarter = c(rep(c(1, 2, 3), each = 3), rep(4, 1))))

biannually_subsample <- function(df){
  ## quarterly subsample
  biannual_vec <- c()

  for (i in 1:length(unique(df$Year))){
    for (j in 1:length(unique(df$Half))){
      year <- unique(df$Year)[i]
      half <- unique(df$Half)[j]
      tmp <- df |> dplyr::filter(Year == year & Half == half) |>
        dplyr::slice(1)
      biannual_vec <- c(biannual_vec, tmp$Value)
    }
  }

  biannual_vec <- as.numeric(biannual_vec)
  return(biannual_vec)
}

#' Generate vector of quarterly data
#'
#' @param df A data frame with columns for Value, Date, Year, Month, Day, Quarter
#' @return A numeric vector
#' @export
#'
#' @examples
#' quarterly_subsample(data.frame(Value = rnorm(10, 2, 0.1), Year = sample(2000:2005, 10, replace = TRUE), Month = c(1:10), Day = sample(1:30, 10)) |> dplyr::mutate(Date = paste0(Year, "-", Month, "-", Day), Quarter = c(rep(c(1, 2, 3), each = 3), rep(4, 1))))
quarterly_subsample <- function(df){
  ## quarterly subsample
  quarter_vec <- c()

  for (i in 1:length(unique(df$Year))){
    for (j in 1:length(unique(df$Quarter))){
      year <- unique(df$Year)[i]
      quarter <- unique(df$Quarter)[j]
      tmp <- df |> dplyr::filter(Year == year & Quarter == quarter) |>
        dplyr::slice(1)
      quarter_vec <- c(quarter_vec, tmp$Value)
    }
  }

  quarter_vec <- as.numeric(quarter_vec)
  return(quarter_vec)
}

#' Generate vector of monthly data
#'
#' @param df A data frame with columns for Value, Date, Year, Month, Day, Quarter
#' @return A numeric vector
#' @export
#'
#' @examples
#' monthly_subsample(data.frame(Value = rnorm(10, 2, 0.1), Year = sample(2000:2005, 10, replace = TRUE), Month = c(1:10), Day = sample(1:30, 10)) |> dplyr::mutate(Date = paste0(Year, "-", Month, "-", Day), Quarter = c(rep(c(1, 2, 3), each = 3), rep(4, 1))))

monthly_subsample <- function(df){
  ## monthly subsample
  month_vec <- c()

  for (i in 1:length(unique(df$Year))){
    for (j in 1:length(unique(df$Month))){
      year <- unique(df$Year)[i]
      month <- unique(df$Month)[j]
      tmp <- df |> dplyr::filter(Year == year & Month == month) |>
        dplyr::slice(1)
      month_vec <- c(month_vec, tmp$Value)
    }
  }

  month_vec <- as.numeric(month_vec)

  return(month_vec)
}

#' Calculate coefficient of variation
#'
#' @param vec a vector of numbers
#' @return A number
#' @export
#'
#' @examples
#' get_cv(c(1, 2, 1, 0, 1, 3))

get_cv <- function(vec){
  round(sd(vec)/mean(vec), 3)
}


#' Generate MC simulated data
#'
#' @param x a data frame with columns nobs, slope, var_ratio, var_r, cv, nsim
#' @return data frame list
#' @export
#'
#' @examples
#' get_data(data.frame(nobs = 10, slope = 2, var_ratio = .1, var_r = .05, cv = .2, nsim = 10))
get_data <- function(x) {
  t <- seq(from = 1, to = x[[1]], by = 1)
  trend_comp <- x[[2]] * t
  y <- trend_comp + rnorm(x[[1]], mean = 100, sd = 100 * x[[3]])
  data.frame(t, y)
}

#' run tests on nsim datasets
#'
#' @param x a data frame with columns nobs, slope, var_ratio, var_r, cv, nsim
#' @return Power_MK vector
#' @export
#'
#' @examples
#' simulate(data.frame(nobs = 10, slope = 2, var_ratio = .1, var_r = .05, cv = .2, nsim = 10))
simulate <- function(x){

  # Get data
  dt <- t(replicate(1000, get_data(x)))

  # Run Mann-Kendall test
  mktest <- apply(dt, 1, function(x) Kendall::MannKendall(x[[2]])$sl)
  mktest_power <- sum(mktest < 0.05)/1000
  list(Power_MK = mktest_power)
}


#' Create hypothetical power plots
#'
#' @param df0 a data frame with columns nobs, slope, var_ratio, var_r, cv, nsim, Power_MK
#' @param frequency sampling frequency
#' @param cv coefficient of variation for respective frequency
#' @param sample_size sample size for respective frequency
#' @return ggplot object
#' @export
#'
#' @examples
#' make_plot(data.frame(nobs = c(2, 10), slope = c(2, 2), var_ratio = c(.1, .1), var_r = c(0.05, .05), cv = c(.2, .2), nsim = c(10, 10), Power_MK = c(.1, .311)), "Annually", 0.5, 12)
make_plot <- function(df0, frequency, cv, sample_size) {

  library(ggplot2)

  p <- ggplot(data = df0, aes(x = nobs, y = as.numeric(Power_MK),
                              color = as.factor(slope))) +
    geom_point() +
    geom_line() +
    geom_point(aes(x = sample_size, y = 0),
               colour = "black",
               size = 3,
               shape = 8) +
    labs(title = paste0(frequency, " Sampling: CV = ", cv),
         x = 'Sample Size',
         y = 'Hypothetical Power',
         caption = "* represents the corresponding frequency's sample size.") +
    scale_color_hue(l = 50) +
    theme_bw(base_size = 16) +
    theme(plot.title = element_text(size = 16, hjust = 0.5),
          axis.title.y = element_text(size = 16, color = 'black'),
          axis.title.x = element_text(size = 16, color = 'black'),
          axis.text.x = element_text(size = 16, color = 'black'),
          axis.text.y = element_text(size = 16, color = 'black'),
          legend.key.height = unit(0.75, 'line')) +
    guides(color = guide_legend(title = "Trend Slope"))

  return(p)
}

#' Find Mode of a vector
#'
#' @param x a vector
#' @return mode of the vector
#' @export
#'
#' @examples
#' Mode(c("cat", "dog", "cat", "fish"))
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

