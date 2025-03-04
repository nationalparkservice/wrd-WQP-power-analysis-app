library(shiny)
library(Kendall)
library(tidyverse)
library(dataRetrieval)
library(kableExtra)
library(npsMK)

set.seed(1999)

ui <- fluidPage(
  theme = shinythemes::shinytheme("paper"),
  titlePanel("Mann-Kendall Test"),
  sidebarLayout(
    sidebarPanel(
      selectInput(inputId = "frequency",
                  label = "Frequency",
                  choices = c("Annually", "Biannually", "Quarterly", "Monthly"),
                  selected = "Annually"),
      selectizeInput(inputId = "network",
                  label = "Network",
                  choices = unique(ids$Network),
                  selected = unique(ids$Network)[1]),
      selectizeInput(inputId = "sitename",
                     label = "Site",
                     choices = unique(ids$Name),
                     selected = unique(ids$Name)[1]),
      selectInput(inputId = "parameter_cd",
                  label = "Parameter",
                  choices = "Nitrate",
                  selected = "Nitrate"),
      sliderInput(inputId = "trend_slope",
                  label = "Trend Slope",
                  min = 0.5,
                  max = 10,
                  step = 0.5,
                  value = 2),
      actionButton("go", "Perform Analysis")
    ),
    mainPanel(
      htmlOutput("description"),
      htmlOutput("site_info"),
      plotOutput("test_plot")
    )
  )
)



server <- function(input, output, session) {


  ## update parameter input
  observe({
    selectedNetwork = input$network

    Sites.selected = unique((ids |> dplyr::filter(Network == selectedNetwork))$Name)

    shiny::validate(need(length(Sites.selected) > 0,
                         "No nutrient data available for selected site."))

    updateSelectizeInput(session = session,
                      inputId = "sitename",
                      choices = Sites.selected,
                      selected = Sites.selected[1])
  })

  observe({
    selectedSite = input$sitename
    SiteId = ids[ids$Name == selectedSite,]$ID
    SiteParams_all = unique(dataRetrieval::readWQPqw(siteNumbers = SiteId,
                                                     "",
                                                     "",
                                                     "")$CharacteristicName)
    toMatch <- c("nitr", "phos", "ammon")

    SiteParams <- unique(grep(paste(toMatch,collapse="|"),
                              SiteParams_all, value=TRUE,
                              ignore.case = TRUE))

    shiny::validate(need(length(SiteParams) > 0,
                         "No nutrient data available for selected site."))

    updateSelectInput(session = session,
                      inputId = "parameter_cd",
                      choices = SiteParams,
                      selected = SiteParams[1])
  })

  siteID <- reactive({

    ids[ids$Name == input$sitename,]$ID

  })

  Sample <- eventReactive(input$go, {

    df <- dataRetrieval::readWQPqw(siteNumbers = siteID(),
                                   "",
                                   "",
                                   parameterCd = input$parameter_cd) |>
      data.frame() |>
      dplyr::filter(ResultValueTypeName == "Actual") |>
      dplyr::filter(ActivityTypeCode == "Sample-Routine")

    units = Mode(df$ResultMeasure.MeasureUnitCode)
    df <- df |>
      dplyr::filter(ResultMeasure.MeasureUnitCode == units) |>
      dplyr::rename(Date = ActivityStartDate,
                    Value = ResultMeasureValue)

    return(df)

  })

  df <- eventReactive(input$go, {
    Sample() |>
      dplyr::select(Value, Date) |>
      dplyr::mutate(Month = lubridate::month(Date),
                    Year = lubridate::year(Date)
      ) |>
      dplyr::mutate(Quarter = ifelse(as.numeric(Month) %in% c(1:3), "Q1",
                                     ifelse(as.numeric(Month) %in% c(4:6), "Q2",
                                            ifelse(as.numeric(Month) %in% c(7:9), "Q3", "Q4"))),
                    Half = ifelse(as.numeric(Month) %in% c(1:6), "H1", "H2")) |>
      dplyr::arrange(Date) |>
      na.omit()
  })

  freq_list <- eventReactive(input$go, {

    # shiny::validate(shiny::need(nrow(Sample()) > 0,
    #                 "No data for selected site and parameter."))

    list(annual_vec = yearly_subsample(df()),
         quarterly_vec = quarterly_subsample(df()),
         biannual_vec = biannually_subsample(df()),
         monthly_vec = monthly_subsample(df())
    )
  })

  frequencies_table <- eventReactive(input$go, {
    data.frame("Frequency" = c("Annually", "Biannually", "Quarterly", "Monthly"),
               "N" = c(length(freq_list()$annual_vec),
                       length(freq_list()$biannual_vec),
                       length(freq_list()$quarterly_vec),
                       length(freq_list()$monthly_vec)),
               "CV" = c(get_cv(freq_list()$annual_vec),
                        get_cv(freq_list()$biannual_vec),
                        get_cv(freq_list()$quarterly_vec),
                        get_cv(freq_list()$monthly_vec)))
  })

  sample_sizes <- reactive({
    shiny::validate(need(dplyr::all_of(frequencies_table()$N) > 5,
                         "Not enough data."))
    samp_size0 <- frequencies_table()[frequencies_table()$Frequency == input$frequency,]$N

    if(samp_size0 < 100){
      samp_size = samp_size0*3.5
    } else if (samp_size0 >= 100 & samp_size0 < 200){
      samp_size = samp_size0*2
    } else if (samp_size0 >= 200 & samp_size0 < 500){
      samp_size = samp_size0*1.5
    } else{
      samp_size = samp_size0*1.25
    }

    if (samp_size > 100 & samp_size < 300){
      nobs <- seq(from = 5,
                  to = samp_size,
                  by = 10)
    } else if (samp_size >= 300 & samp_size < 500){
      nobs <- seq(from = 5,
                  to = samp_size,
                  by = 20)
    } else if (samp_size >= 500 & samp_size < 1000){
      nobs <- seq(from = 5,
                  to = samp_size,
                  by = 50)
    }else if (samp_size >= 1000){
      nobs <- seq(from = 5,
                  to = samp_size,
                  by = 100)
    }else{
      nobs <- seq(from = 5,
                  to = samp_size,
                  by = 5)
    }

    return(nobs)
  })

  trend_slopes <- eventReactive(input$go, {
    unique(c(.1, 1, 10, input$trend_slope))
  })

  test0 <- eventReactive(input$go, {
    rbind(data.frame(
      nobs = sample_sizes(),
      slope = 0.1,
      cv = frequencies_table()[frequencies_table()$Frequency == input$frequency,]$CV,
      var_ratio = NA,
      var_r = NA,
      nsim = 10000
    ),
    data.frame(
      nobs = sample_sizes(),
      slope = 1,
      cv = frequencies_table()[frequencies_table()$Frequency == input$frequency,]$CV,
      var_ratio = NA,
      var_r = NA,
      nsim = 10000
    ),
    data.frame(
      nobs = sample_sizes()    ,
      slope = 10,
      cv = frequencies_table()[frequencies_table()$Frequency == input$frequency,]$CV,
      var_ratio = NA,
      var_r = NA,
      nsim = 10000
    ),
    data.frame(
      nobs = sample_sizes(),
      slope = input$trend_slope,
      cv = frequencies_table()[frequencies_table()$Frequency == input$frequency,]$CV,
      var_ratio = NA,
      var_r = NA,
      nsim = 10000 )
    )  |>
      dplyr::mutate(
        var_r = (100 * cv)^2) |>
      dplyr::mutate(
        var_ratio = (nobs^2 - 1) * slope^2 / (12 * var_r)
      )

  })

  Power_calc <- eventReactive(input$go, {

    shiny::withProgress(message = "Analyzing...", {

      # cl <- parallel::makeCluster(ncores)
      # parallel::clusterExport(
      #   cl,
      #   list("get_data", "MannKendall")
      # )

      set.seed(1999)

      # zz <- parallel::parApply(cl, test0(), 1, simulate)

      # parallel::stopCluster(cl)

      # return(zz)
      zz <- apply(test0(), 1, npsMK::simulate) |>
        unlist() |>
        as.vector()

    })

    # parallel::stopCluster(cl)

    return(zz)

  })

  Power_calc_df <- eventReactive(input$go, {

    # do.call(rbind, lapply(Power_calc(), as.data.frame))
    data.frame(Power_MK = Power_calc())
  })

  test1 <- eventReactive(input$go, {
    cbind(test0(), Power_calc_df())
  })

  # output$site_info <- renderText({
  #   paste0(site_info()$station_nm, " - ", site_info()$param.nm)
  # })

  description_tab <- eventReactive(input$go, {
    data.frame(x = input$network, y = input$sitename, z = input$parameter_cd)
  })


  output$description <- renderText({
    kableExtra::kable(description_tab(), col.names = NULL)  |>
      kableExtra::kable_styling(bootstrap_options= c("striped","condensed"))
  })

  output$site_info <- renderText({
    kableExtra::kable(frequencies_table()) |>
      kableExtra::kable_styling(bootstrap_options= c("striped","condensed"))

  })

  freq_plot <- eventReactive(input$go, {
    make_plot(test1(), input$frequency,
              frequencies_table()[frequencies_table()$Frequency == input$frequency,]$CV,
              frequencies_table()[frequencies_table()$Frequency == input$frequency,]$N)
  })
  output$test_plot <- renderPlot({

    freq_plot()
  })

}

shinyApp(ui = ui, server = server)






