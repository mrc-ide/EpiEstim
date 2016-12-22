# TODO: Add comment
#
# Author: annecori
###############################################################################

library(shiny)

###############################################################################

# Define UI for miles per gallon application
shinyUI(pageWithSidebar(

				# Application title
				headerPanel("Estimation of the reproduction number for an outbreak"),

				# Sidebar with controls to select the variable to plot against mpg
				# and to specify whether outliers should be included
				sidebarPanel(

						####################################################################################
						h3("Data"),
						selectInput("dataset", "Choose dataset:",
								list("Upload my own data" = 0,"Measles, 1861, Hagelloch, Germany" = 1,
										"Influenza, 1918, Baltimore" = 2,
										"Smallpox, 1972, Kosovo" = 3,
										"Severe acute respiratory syndrome (SARS), 2003, Hong Kong" = 4,
										"Influenza, 2009, scholl in Pennsylvania" = 5
								)),

						conditionalPanel(
								condition = "input.dataset == 0",
								fileInput('DataFile', 'Choose CSV File (The file must have 2 columns called "Time" and "Incidence"):',
										accept=c('text/csv', 'text/comma-separated-values,text/plain'))
						),
						tags$hr(),
						h3("Double Censored"),
						selectInput("whichCensor", "Double Censored?",
						            list("No" = 1,
						                 "Yes" = 2)),
						conditionalPanel(
						  condition = "input.whichCensor == 2",
						  selectInput("dataset1", "Choose dataset for CDT:",
						              list("Upload my own data" = 0, "Measles, 1861, Hagelloch, Germany" = 1,
						                   "Simulation" = 2
						              ))
						),
# 						selectInput("dataset1", "Choose dataset for CDT:",
# 						            list("Upload my own data" = 0, "Measles, 1861, Hagelloch, Germany" = 1,
# 						                 "Simulation" = 2
# 						            )),

						conditionalPanel(
						  condition = "input.whichCensor == 2 & input.dataset1 == 0",
						  fileInput('DataFile1', 'Choose CSV File (The file must have 5 columns called "EL", "ER", "SL", "SR" and "type"):',
						            accept=c('text/csv', 'text/comma-separated-values,text/plain'))
						),

						tags$hr(),

						####################################################################################
						h3("Instantaneous or case reproduction number"),
						selectInput("whichR", "Choose reproduction number to estimate:",
								list("Instantaneous reproduction number (Cori et al.)" = 1,
										"Case reproduction number (Wallinga and Teunis)" = 2)),
						tags$hr(),

						####################################################################################
						h3("Serial interval distribution"),
						conditionalPanel(
								condition = "input.whichR == 1",
								selectInput("method1", "Serial interval specification:",
										list("Parametric" = 1,"Non parametric" = 2,
												"Parametric, accounting for uncertainty" = 3
										))),

						conditionalPanel(
								condition = "input.whichR == 2",
								selectInput("method2", "Serial interval specification:",
										list("Parametric" = 1,"Non parametric" = 2
										))),

						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 1) || (input.whichR == 2 && input.method2 == 1)",
								textInput("meanSI", "Mean serial interval:", "3")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 1) || (input.whichR == 2 && input.method2 == 1)",
								textInput("sdSI", "Standard deviation of the serial interval:", "2")
						),

						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 2) || (input.whichR == 2 && input.method2 == 2)",
								fileInput('SIFile', 'Choose CSV File (The file must have 2 columns called "k" and "w[k]", and start with k=0, w[k]=0):',
										accept=c('text/csv', 'text/comma-separated-values,text/plain'))
						),

						conditionalPanel(condition = "(input.whichR == 1 && input.method1 == 3)",
								h4("Distribution of the mean serial interval")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("meanmeanSI", "Mean:", "3")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("stdmeanSI", "Standard deviation:", "1")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("minmeanSI", "Lower bound:", "1")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("maxmeanSI", "Upper bound:", "5")
						),

						conditionalPanel(condition = "(input.whichR == 1 && input.method1 == 3)",
								h4("Distribution of the standard deviation of the serial interval")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("meanstdSI", "Mean:", "1")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("stdstdSI", "Standard deviation:", "0.5")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("minstdSI", "Lower bound:", "0.25")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								textInput("maxstdSI", "Upper bound:", "1.75")
						),

						conditionalPanel(condition = "(input.whichR == 1 && input.method1 == 3)",
								h4("Sample sizes for simulations")
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								sliderInput("n1SI", "Number of SI distributions sampled (large number will lead to large computing times):", min = 10,
										max = 200,
										value = 50)
						),
						conditionalPanel(
								condition = "(input.whichR == 1 && input.method1 == 3)",
								sliderInput("n2SI", "Posterior sample size for each SI distribution (large number will lead to large computing times):", min = 10,
										max = 200,
										value = 50)
						),


						tags$hr(),

						####################################################################################
						h3("Time windows"),
						#textInput("width", "Width of the time window used for estimation:", "7"),
						sliderInput("width",
								"Width of the time window used for estimation:",
								min = 1,
								max = 30,
								value = 7),

						tags$hr(),

						####################################################################################
						conditionalPanel(
								condition = "input.whichR == 1",
								h3("Prior distribution")
						),
						conditionalPanel(
								condition = "input.whichR == 1",
								textInput("meanPrior", "Mean prior for R:", "5")
						),
						conditionalPanel(
								condition = "input.whichR == 1",
								textInput("sdPrior", "Standard deviation of the prior for R:", "5")
						),
						conditionalPanel(
								condition = "input.whichR == 1",
								tags$hr()
						),

						####################################################################################
						conditionalPanel(
								condition = "input.whichR == 2",
								h3("Number of simulations")
						),
						conditionalPanel(
								condition = "input.whichR == 2",
								sliderInput("nSim", "Number of simulations used to estimate variability around central estimate of Rc:", min = 10,
										max = 100,
										value = 20)
						),
						conditionalPanel(
								condition = "input.whichR == 2",
								tags$hr()
						),

						####################################################################################
						h3("Reference"),
						helpText("A. Cori, N.M. Ferguson, C. Fraser and S. Cauchemez, A new framework and software to estimate time-varying reproduction numbers during epidemics, AJE, 2013.")



#submitButton("Update Results")
				),

# Show a plot of the generated distribution
				mainPanel(
						tabsetPanel(
								tabPanel("Data", downloadButton('downloadData', 'Download'),
										helpText("The downloaded file will be empty if no dataset has been selected."),
										tableOutput("viewData")
								),
								tabPanel("CDT Data",
								         helpText("The downloaded file will be empty if no dataset has been selected."),
								         tableOutput("viewData1")
								),
								tabPanel("Plots", downloadButton('savePlot', 'Save image'),
										helpText("The downloaded file will be empty if no dataset has been selected."),
										plotOutput("plot")
								),
								tabPanel("CDT Estimation",
								         helpText("Estimates from CDT; It might take a few minutes :)"),
								         tableOutput("viewCDT")),
								tabPanel("Estimated Reproduction number", downloadButton('downloadR', 'Download'),
										helpText("The downloaded file will be empty if no dataset has been selected."),
										tableOutput("viewR")),
								tabPanel("Serial interval distribution",
										helpText("For parametric serial intervals, this gives the corresponding discrete distribution,"),
										helpText("For non parametric serial intervals, this gives the corresponding mean and standard deviation,"),
										helpText("For uncertain serial intervals, this is empty."),
										downloadButton('downloadSI', 'Download'),
										helpText("The downloaded file will be empty if no dataset has been selected."),
										tableOutput("viewSI"))
						)
				)
		))