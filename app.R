rm(list = ls())

library(LBSPR)
library(shiny)
library(reshape2)
library(ReporteRs)
library(ggplot2)

source("utilities.R")

# Define the UI
ui <- tagList(
  tags$style(HTML("div#lb_indicator_plot img {width: auto; height: auto; max-width: auto; max-height: 800px;}")),
  navbarPage(title = "LBIndicator Application",
             id = "navbar",
             tabPanel(title = "Instructions",
                      p("To use the LBI Application you will need: 1) a length frequency distribution 2) weight at length data, and 3) estimates of the life history parameters. The following paragraphs outline the steps to use the LBI Application.  Each heading refers to a tab on the menu."),
                      h4("Upload Data"),
                      p("The first step is to upload two CSV (comma separated variable) files containing the length frequency distribution and mean weight at length.  The file must be in CSV format and contain only numeric values except for the header row which can contain labels.  Multiple years of data should be placed in seperate columns."),
                      p("Length frequency data must have the midpoints of the length classes (the length bins) in the first column, and numeric values for all counts (i.e., all columns are the same length). Length measurements should be raw numbers, each column representing a different year."),
                      p("Two example stocks have been included. Click the 'I want to look at sample data' radio button to make see the exact format."),
                      h4("Plot length frequency distribution"),
                      p("Slide the bar according to the bin width that you want to aggregate the data. Both length frequency distribution and weight at length data frames are modified concurrently. Only the LFD is displayed"),
                      h4("LBI plot"),
                      p("Enter the life history parameters for your stock. The application will automatically generate the LBI indicators and reference points."),
                      h4("LBI table"),
                      p("The last 3 years of LBI indicators are generated in a table relative to reference points. The shading of the cell indicates the relative status."),
                      h4("Downloads"),
                      p("All figures and tables can be downloaded as a .docx report or as individual .png images.")
             ), #tabPanel
             tabPanel(title = "Upload data",
                      sidebarLayout(
                        sidebarPanel(
                          radioButtons(inputId = "input_type", 
                                       label = "How do you want to upload data?",
                                       list("I have a length frequency distribution I want to explore." = "input_file",
                                            "I want to look at sample data" = "example_file"),
                                       selected = ""),
                          
                          conditionalPanel(
                            condition = "input.input_type == 'input_file'", 
                            fileInput(inputId = 'choose_input_data',
                                      label = "Choose length frequency distribution to upload",
                                      accept = c(
                                        'text/csv',
                                        'text/comma-separated-values',
                                        'text/tab-separated-values',
                                        'text/plain',
                                        '.csv',
                                        '.tsv'
                                      )
                            ), #fileInput
                            fileInput(inputId = 'choose_weight_data',
                                      label = "Choose weight at length data to upload",
                                      accept = c(
                                        'text/csv',
                                        'text/comma-separated-values',
                                        'text/tab-separated-values',
                                        'text/plain',
                                        '.csv',
                                        '.tsv'
                                      )
                            ), #fileInput
                            radioButtons(inputId = 'separator_button',
                                         label = 'Separator',
                                         selected = ",",
                                         choices =  c(
                                           Comma = ',',
                                           Semicolon = ';',
                                           Tab = '\t')
                            ), #radioButtons
                            actionButton(inputId = "load_input_data",
                                         label = "Load Data",
                                         icon("cloud-upload")
                            ) #fileInput
                          ), #conditionalPanel
                          
                          conditionalPanel(
                            condition = "input.input_type == 'example_file'", 
                            selectInput(inputId = 'choose_example_data',
                                        label = "Example Data",
                                        choices = c(
                                          "",
                                          "Frequency - Nephrops" = "freqNeph",
                                          "Frequency - Gadoid" = "freqGad")
                            ) #selectInput
                          ) #ConditionalPanel
                        ), #sidebarPanel
                        mainPanel(
                          h3("Length frequency distribution"),
                          h4("Steps:"),
                          tags$ol(
                            tags$li("Upload a CSV data file or select an example data set from the menu on the left."),
                            tags$li("Length frequency distributions should have the length bin in the first column and annual frequency in subsequent columns."),
                            tags$li("Double check the first 6 rows to make sure everything looks OK.")
                          ),
                          dataTableOutput("data_table")
                        ) # mainPanel
                      ) # sidebar layout
             ), #tabPanel
             tabPanel(title = "Plot length frequency data",
                      sidebarLayout(
                        sidebarPanel(
                          sliderInput(inputId = "binwidth_value",
                                      label = "Bin width for length frequencies",
                                      min = 1,
                                      max = 100,
                                      value = ""
                          ), #sliderInput
                          radioButtons(inputId = 'length_units', 
                                       label = "",
                                       selected = ",",
                                       choices =  c(
                                         "cm" = 'cm',
                                         "mm" = 'mm',
                                         "in" = 'in')
                          ) #radioButtons
                        ), #sidebarPanel
                        mainPanel(
                          h3("Plots go here"),
                          h4("Steps:"),
                          tags$ol(
                            tags$li("Select the optimal bin width for length frequencies and the proper units. Note: units are only used for plot labels."),
                            tags$li("If the length frequencies displayed are not correct, please review the raw data")
                          ),
                          # Reactive length frequencies
                          plotOutput("histogram_plot")
                        ) # mainPanel
                      ) # sidebar layout
             ), #tabPanel
             tabPanel(title = "Length based indicator - plots",
                      sidebarLayout(
                        sidebarPanel(
                          # Input Linf and Lmat parameters
                          tags$div(
                            HTML("L<sub>&#x221e</sub>")
                          ),
                          numericInput(inputId = 'linf_value',
                                       label = "",
                                       value = "",
                                       min = 1,
                                       max = 100
                          ), # textInput
                          tags$div(
                            HTML("L<sub>mat</sub>")
                          ),
                          numericInput(inputId = 'lmat_value',
                                       label = "",
                                       value = "",
                                       min = 1,
                                       max = 100
                          ) # textInput
                        ), #sidebarPanel
                        mainPanel(
                          h3("Plots go here"),
                          h4("Steps:"),
                          tags$ol(
                            tags$li("Fill in the name of your stock if you want to download summary "),
                            tags$li("If the length frequencies displayed are not correct, please review the raw data")
                          ),
                          # Reactive table 
                          plotOutput("lb_indicator_plot")   
                        ) # mainPanel
                      ) # sidebar layout
             ), #tabPanel
             tabPanel(title = "Length based indicator - table",
                      # sidebarLayout(
                      # sidebarPanel(
                      
                      # ), #sidebarPanel
                      # mainPanel(
                      h3("Stock status relative to reference points"),
                      tableOutput(outputId = "lb_indicator_table")
                      # ) # mainPanel
                      # ) # sidebar layout
             ), #tabPanel
             tabPanel("Downloads",
                      # h3(":"),
                      # tags$ol(
                      #   tags$li("Fill in the name of your stock if you want to download summary "),
                      #   tags$li("")
                      # ),
                      # h4("Binned Data"),
                      # downloadButton('download_data',
                      #                'Download binned data (.csv)'
                      # ), # downloadButton
                      # br(),
                      h4("Plots"),
                      downloadButton('download_lb_plot',
                                     'Download indicator time series plots (.png)' 
                      ), # downloadButton
                      # downloadButton('download_histogram',
                      #                'Download histogram (.png)' 
                      # ), # downloadButton
                      br(),
                      h4("Summary sheet"),
                      downloadButton('download_doc',
                                     'Download summary document (.docx)'
                      ), # downloadButton
                      br(),
                      br(),
                      textInput(inputId = 'stock_name',
                                label = "Enter stock name as filename",
                                value = "My stock"
                      ), # textInput
                      h3("Sources:"),
                      p("The application was created by Scott Large from example R code generously supplied to WKLIFE5 by T. Miethe and C. Silva. Inspiration for the application was drawn from J. Cope and A. Hordyk."),
                      p("The application should be used as a tool to explore data and a thorough audit of the length based indicators should be conducted before they can be used for decision making."),
                      p("Source code can be found at https://github.com/ices-tools-dev/LBI_shiny. Questions and issues should be directed to the github repository.")
             ) #tabPanel
  ) # NavbarPage
) #tagList

server <- function(input, output, session) {
  
  # LOAD YOUR OWN DATA
  data_from_user <- reactive({
    
    req(input$input_type,
        input$choose_input_data,
        input$load_input_data
    )
    
    datapath <- input$choose_input_data$datapath
    
    dat <- read.csv(file = datapath, 
                    header = TRUE,
                    sep = input$separator_button,
                    stringsAsFactors = FALSE,
                    check.names = FALSE)
    
    if(class(dat) == "data.frame" | class(dat) == "matrix") {
      if(ncol(dat) > 1) {
        chk_NAs <- apply(dat, 2, is.na) # check NAs
        dat <- dat[!apply(chk_NAs, 1, prod), , drop = FALSE]
        dat <- dat[, !apply(chk_NAs, 2, prod), drop = FALSE]
      }
    }
    if(class(dat) == "numeric" | class(dat) == "integer") {
      dat <- dat[!is.na(dat)]
    }
    return(dat)
    
  })
  
  weights_from_user <- reactive({
    
    req(input$input_type,
        input$choose_weight_data,
        input$load_input_data
    )
    
    # weights
    weightpath <- input$choose_weight_data$datapath
    
    wght <- read.csv(file = weightpath, 
                     header = TRUE,
                     sep = input$separator_button,
                     stringsAsFactors = FALSE,
                     check.names = FALSE)
    
    if(class(wght) == "data.frame" | class(wght) == "matrix") {
      if(ncol(wght) > 1) {
        chk_NAs <- apply(wght, 2, is.na) # check NAs
        wght <- wght[!apply(chk_NAs, 1, prod), , drop = FALSE]
        wght <- wght[, !apply(chk_NAs, 2, prod), drop = FALSE]
      }
    }
    if(class(wght) == "numeric" | class(wght) == "integer") {
      wght <- wght[!is.na(wght)]
    }
    
    return(wght)
  })
  
  
  
  data_from_example <- reactive({
    req(
      input$input_type,
      input$choose_example_data
    )
    
    weightpath <- input$choose_weight_data$datapath
    
    
    read.csv(switch(input$choose_example_data,
                    freqNeph = "data/freqNeph.csv",
                    freqGad= "data/freqGad.csv"),
             header = TRUE,
             stringsAsFactors = FALSE,
             check.names = FALSE)
  })
  
  weights_from_example <- reactive({
    req(
      input$input_type,
      input$choose_example_data
    )
    
    read.csv(switch(input$choose_example_data,
                    freqNeph = "data/walNeph.csv",
                    freqGad = "data/walGad.csv"),
             header = TRUE,
             stringsAsFactors = FALSE,
             check.names = FALSE)
  })
  
  
  output$data_table <- renderDataTable({
    
    # Print out first 6 observations
    validate(
      need(input$input_type, "")
    )
    
    if(input$input_type == "input_file"){
      validate(
        need(input$choose_input_data, "Load length frequency data to explore"),
        need(input$choose_weight_data, "Load the corresponding weight at length data, too.")
      )
      dat <- data_from_user()
      wght <- weights_from_user()
      
      validate(
        need(all(dim(dat) == dim(wght)), "Dimensions of length frequency data and weight at length data do not agree")
      )
    }
    if(input$input_type == "example_file"){
      validate(
        need(input$choose_example_data, "Select the example data to explore")
      )
      dat <- data_from_example()
      wght <- weights_from_example()
    }
    head(dat)
  }, options = list(pageLength = -1,
                    searching = FALSE,
                    paging = FALSE,
                    ordering = FALSE,
                    info = FALSE)
  )
  
  output$histogram_plot <- renderPlot({
    
    validate(
      need(input$input_type, "")
    )
    
    if(input$input_type == "input_file"){
      validate(
        need(input$choose_input_data, "Load length frequency data to explore"),
        need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
      )
      
      dat <- data_from_user()
      wght <- weights_from_user()
    }
    if(input$input_type == "example_file"){
      validate(
        need(input$choose_example_data, "Select the example data to explore"),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
      )
      
      dat <- data_from_example()
      wght <- weights_from_example()
    }
    
    bin_plot(dat, 
             binwidth = input$binwidth_value, 
             l_units = input$length_units)
    
  })
  
  
  
  output$lb_indicator_data <- renderDataTable({
    
    validate(
      need(input$input_type, "")
    )
    
    if(input$input_type == "input_file"){
      validate(
        need(input$choose_input_data, "Load length frequency data to explore"),
        need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
      ) 
      dat <- data_from_user()
      wght <- weights_from_user()
    }
    if(input$input_type == "example_file"){
      validate(
        need(input$choose_example_data, "Select the example data to explore"),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
        
      )
      dat <- data_from_example()
      wght <- weights_from_example()
    }
    
    validate(
      need(input$linf_value, "Enter Linf"),
      need(input$lmat_value, "Enter Lmat")
    )
    
    lb_dat <- lb_table(data = dat,
                       binwidth = input$binwidth_value,
                       linf = input$linf_value,
                       lmat = input$lmat_value,
                       weight = wght)
  })
  
  
  
  output$lb_indicator_plot <- renderPlot({
    
    validate(
      need(input$input_type, "")
    )
    
    if(input$input_type == "input_file"){
      validate(
        need(input$choose_input_data, "Load length frequency data to explore"),
        need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
      ) 
      dat <- data_from_user()
      wght <- weights_from_user()
    }
    if(input$input_type == "example_file"){
      validate(
        need(input$choose_example_data, "Select the example data to explore"),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
        
      )
      dat <- data_from_example()
      wght <- weights_from_example()
    }
    
    validate(
      need(input$linf_value, "Enter Linf"),
      need(input$lmat_value, "Enter Lmat")
    )
    
    lb_plot(data = dat,
            binwidth = input$binwidth_value,
            l_units = input$length_units,
            linf = input$linf_value,
            lmat = input$lmat_value,
            weight = wght)
    
  },
  width = 1000,
  height = 600)
  
  
  output$lb_indicator_table <- renderFlexTable({
    
    validate(
      need(input$input_type, "")
    )
    
    if(input$input_type == "input_file"){
      validate(
        need(input$choose_input_data, "Load length frequency data to explore"),
        need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
      ) 
      dat <- data_from_user()
      wght <- weights_from_user()
    }
    if(input$input_type == "example_file"){
      validate(
        need(input$choose_example_data, "Select the example data to explore"),
        need(input$binwidth_value, "Select the binwidth for the histograms"),
        need(input$length_units, "Select the proper units.")
      )
      dat <- data_from_example()
      wght <- weights_from_example()
    }
    
    validate(
      need(input$linf_value, "Enter Linf"),
      need(input$lmat_value, "Enter Lmat")
    )
    
    lb_table(data = dat,
             binwidth = input$binwidth_value,
             l_units = input$length_units,
             linf = input$linf_value,
             lmat = input$lmat_value,
             weight = wght)
  })
  # }
  
# 
#   output$download_data <- downloadHandler(
#     
#     filename = function() {
#       nm <- input$stock_name
#       nm <- gsub(" ", "", nm)
#       if (nchar(nm) < 1) nm <- "stock_name"
#       paste0(nm, '_binned_data.docx')
#     },
#     
#     content = function(file) {
#       validate(
#         need(input$input_type, "")
#       )
#       
#       if(input$input_type == "input_file"){
#         validate(
#           need(input$choose_input_data, "Load length frequency data to explore"),
#           need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
#           need(input$binwidth_value, "Select the binwidth for the histograms"),
#           need(input$length_units, "Select the proper units.")
#         )
#         dat <- data_from_user()
#         wght <- weights_from_user()
#       }
#       if(input$input_type == "example_file"){
#         validate(
#           need(input$choose_example_data, "Select the example data to explore"),
#           need(input$binwidth_value, "Select the binwidth for the histograms"),
#           need(input$length_units, "Select the proper units.")
#           
#         )
#         dat <- data_from_example()
#         wght <- weights_from_example()
#       }
#       
#       bin_mat(data = dat,
#              binwidth = input$binwidth_value,
#              l_units = input$length_units)
#     }
#   )
# 
#   output$download_histogram <- downloadHandler(
#     
#   )
  
  output$download_doc <- downloadHandler(
    
    filename = function() {
      nm <- input$stock_name
      nm <- gsub(" ", "", nm)
      if (nchar(nm) < 1) nm <- "stock_name"
      paste0(nm, '_LBI_summary.docx')
    },
    
    content = function(file) {
      validate(
        need(input$input_type, "")
      )
      
      if(input$input_type == "input_file"){
        validate(
          need(input$choose_input_data, "Load length frequency data to explore"),
          need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
          need(input$binwidth_value, "Select the binwidth for the histograms"),
          need(input$length_units, "Select the proper units.")
        )
        dat <- data_from_user()
        wght <- weights_from_user()
      }
      if(input$input_type == "example_file"){
        validate(
          need(input$choose_example_data, "Select the example data to explore"),
          need(input$binwidth_value, "Select the binwidth for the histograms"),
          need(input$length_units, "Select the proper units.")
          
        )
        dat <- data_from_example()
        wght <- weights_from_example()
      }
      
      validate(
        need(input$linf_value, "Enter Linf"),
        need(input$lmat_value, "Enter Lmat")
      )
      
      lb_doc(data = dat,
             binwidth = input$binwidth_value,
             l_units = input$length_units,
             linf = input$linf_value,
             lmat = input$lmat_value,
             stock = input$stock_name,
             weight = wght,
             filename = file)
      
    })
  
  output$download_lb_plot <- downloadHandler(
    
    filename = function() {
      nm <- input$stock_name
      nm <- gsub(" ", "", nm)
      if (nchar(nm) < 1) nm <- "stock_name"
      paste0(nm, '_LBI.png')
    },
    
    content = function(file) {
      validate(
        need(input$input_type, "")
      )
      
      if(input$input_type == "input_file"){
        validate(
          need(input$choose_input_data, "Load length frequency data to explore"),
          need(input$choose_weight_data, "Load the corresponding weight at length data, too."),
          need(input$binwidth_value, "Select the binwidth for the histograms"),
          need(input$length_units, "Select the proper units.")
        ) 
        dat <- data_from_user()
        wght <- weights_from_user()
      }
      if(input$input_type == "example_file"){
        validate(
          need(input$choose_example_data, "Select the example data to explore"),
          need(input$binwidth_value, "Select the binwidth for the histograms"),
          need(input$length_units, "Select the proper units.")
          
        )
        dat <- data_from_example()
        wght <- weights_from_example()
      }
      
      validate(
        need(input$linf_value, "Enter Linf"),
        need(input$lmat_value, "Enter Lmat")
      )
      
      png(file,
          bg = "white",
          pointsize = 5,
          units = "cm",
          width = 16,
          height = 15,
          res = 400)
      
      lb_plot(data = dat,
              binwidth = input$binwidth_value,
              l_units = input$length_units,
              linf = input$linf_value,
              lmat = input$lmat_value,
              weight = wght)
      
      dev.off()
    })
  
}
# 

# Return a Shiny app object
shinyApp(ui = ui, server = server)

