library(shiny)
library(shinydashboard)

# Define UI for application that draws a histogram



ui <- dashboardPage(
  dashboardHeader(title="Exploration of Coevolutionary Models"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Gavrilets & Hastings 1998", tabName = "Gav_hast", icon=icon("bar-chart-o")),
      menuItem("Nuismer, Thompson, and Gomulkiewicz 1999", tabName="Nuismer_99",icon=icon("bar-chart-o"))
    )
  ),
  dashboardBody(
    tabItems(
      tabItem(tabName="Gav_hast",
              fluidRow(
                plotOutput("Gav_Hastings_phase_plot")
              ),
              fluidRow(
                shinydashboard::box(    sliderInput("selected_a",
                                                    "Indirect within species interaction strength, spp 1: (a)",
                                                    min = -2,
                                                    max = 2,
                                                    step = 0.1,
                                                    value=0.05),
                                        sliderInput("selected_b",
                                                    "Direct, between species interaction strenght, spp 1: (b)",
                                                    min = -2,
                                                    max = 2,
                                                    step=0.1,
                                                    value = 0.05),
                                        sliderInput("selected_c",
                                                    "Direct, between species interaction strength, spp 2: (c)",
                                                    min = -2,
                                                    max = 2,
                                                    step = 0.1,
                                                    value=0.05),
                                        sliderInput("selected_d",
                                                    "Indirect, within species interaction strength, spp 2: (d)",
                                                    min = -2,
                                                    max = 2,
                                                    step = 0.1,
                                                    value=0.05),
                )
                )
              ),
      tabItem(tabName="Nuismer_99",
              fluidRow(
                plotOutput("Nuismer_99_plot")
              ),
              fluidRow(
                shinydashboard::box(title="Starting allele frequencies",
                                    sliderInput("selected_x1",
                                                "x1:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.5),
                                    sliderInput("selected_y1",
                                                "y1:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.5),
                                    sliderInput("selected_x2",
                                                "x2:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.5),
                                    sliderInput("selected_y2",
                                                "y2:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.5)),
                shinydashboard::box(title="Model parameters",
                                    sliderInput("selected_C",
                                                "Sensitivity of fitness to change in other species trait, C:",
                                                min=0,
                                                max=1,
                                                step=0.02,
                                                value=0.02),
                                    sliderInput("selected_Km",
                                                "Strength of mutualistic interaction, Km:",
                                                min=0,
                                                max=1,
                                                step=0.02,
                                                value=0.02),
                                    sliderInput("selected_Ka",
                                                "Strength of antagonistic interaction, Ka:",
                                                min=0,
                                                max=1,
                                                step=0.02,
                                                value=0.02),
                                    sliderInput("migration",
                                                "Migration rate, m:",
                                                min=0,
                                                max=1,
                                                step=0.02,
                                                value=0.02),
                                    sliderInput("generations",
                                                "Number of generations",
                                                min=0,
                                                max=4000,
                                                step=50,
                                                value=500),
                                    )
              )
              ))
    )
  )


# 
# C <- 0.02
# C2 <- 0.02
# Km <- 0.04
# Km2 <- 0.04
# Ka <- 0.02
# Ka2 <- 0.02
# m <- 0.05
# gens <- 4000
# 
# 


# ui <- shinyUI(fluidPage(
#   
#   # Application title
#   titlePanel("Gavrilets & Hastings 1998"),
#   
#   # Sidebar with a slider input for number of bins
#   sidebarLayout(
#     sidebarPanel(
#       sliderInput("selected_a",
#                   "a:",
#                   min = -2,
#                   max = 2,
#                   step = 0.1,
#                   value=0.05),
#       sliderInput("selected_b",
#                              "b:",
#                              min = -2,
#                              max = 2,
#                   step=0.1,
#                              value = 0.05),
#       sliderInput("selected_c",
#                              "c:",
#                              min = -2,
#                              max = 2,
#                              step = 0.1,
#                   value=0.05),
#     sliderInput("selected_d",
#                              "d:",
#                              min = -2,
#                              max = 2,
#                              step = 0.1,
#                 value=0.05),
#     ),
#     # Show a plot of the generated distribution
#      mainPanel(
#          plotOutput("Gav_Hastings_phase_plot")
#      )
#   )
# ))

library(phaseR)
library(deSolve)

#source("Gavrilets.Hastings.R")
chase_test <- function(t,state,parameters){
  p1 <- state[1]
  q1 <- 1-p1
  p2 <- state[2]
  q2 <- 1-p2
  a <- parameters[1]
  b <- parameters[2]
  c <- parameters[3]
  d <- parameters[4]
  dp1dt <- p1*q1*((a*(p1-q1))+(b*(p2-q2)))
  dp2dt <- p2*q2*((c*(p1-q1))+(d*(p2-q2)))
  dt <- numeric(2)
  print(dp1dt)
  
  dt[1] <- dp1dt
  dt[2] <- dp2dt
  list(dt)
}
migrate_primes <- function(xy_df, C, Km, Ka,m){
  #starting allele frequencies
  y1 <- xy_df$y1
  x1 <- xy_df$x1
  y2 <- xy_df$y2
  x2 <- xy_df$x2
  
  
  #relative fitness values
  #mutualist community
  #benefits for x from y
  
  WA1 <- 1+(C*y1)
  Wa1 <- 1+(C*(1-y1))
  #benefits for y from x
  WB1 <- 1+(Km*(x1))
  Wb1 <- 1+(Km*(1-x1))
  
  #antagonist community:
  WA2 <- 1+(C*y2)
  Wa2 <- 1+(C*(1-y2))
  WB2 <- 1-(Ka*x2)
  Wb2 <- 1-(Ka*(1-x2))
  
  #post migration allele frequencies
  
  xstar1 <- (x1*(1-m))+(m*x2)
  xstar2 <- (x2*(1-m))+(m*x1)
  ystar1 <- (y1*(1-m))+(m*y2)
  ystar2 <- (y2*(1-m))+(m*y1)
  
  #next generations allele freqs, based on fitness and post-migration freqs
  
  next_y1 <- (ystar1*WB1)/((ystar1*WB1)+((1-ystar1)*Wb1))
  next_x1 <- (xstar1*WA1)/((xstar1*WA1)+((1-xstar1)*Wa1))
  next_y2 <- (ystar2*WB2)/((ystar2*WB2)+((1-ystar2)*Wb2))
  next_x2 <- (xstar2*WA2)/((xstar2*WA2)+((1-xstar2)*Wa2))
  
  return(c(next_y1, next_x1, next_y2, next_x2))
  
}

buncha_gens <- function(gens,xy_df,C,Km,Ka,m){
  
  
  
  
  
  for(generation in 1:gens){
    xy_df[generation+1,] <- migrate_primes(xy_df[generation,],C ,Km, Ka,m)
  }
  return(xy_df)
}

test_df <- data.frame(matrix(NA,nrow=gens+1,ncol=4))
names(test_df) <- c("y1","x1","y2","x2")


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output) {
  
  output$Gav_Hastings_phase_plot <- renderPlot({
    
    # generate bins based on input$bins from ui.R
    chase_test.flow_field <- flowField(chase_test,xlim=c(0,1),ylim=c(0,1),
                                       parameters=c(input$selected_a,input$selected_b,input$selected_c,input$selected_d),
                                       points=20,add=F, col='black')
    chase_test.null <- nullclines(chase_test,xlim=c(0,1),ylim=c(0,1),
                                  parameters=c(input$selected_a,input$selected_b,input$selected_c,input$selected_d),
                                  points=20,add=T)
    
  })
  output$Nuismer_99_plot <- renderPlot({
    test_df <- data.frame(matrix(NA,nrow=input$generations+1,ncol=4))
    names(test_df) <- c("y1","x1","y2","x2")
    test_df[1,] <- c(input$selected_y1,input$selected_x1, input$selected_y2, input$selected_x2)
    
    
    results <- buncha_gens(input$generations, test_df, input$selected_C, input$selected_Km, input$selected_Ka, input$migration)
    long_results <- gather(results)
    
    long_results$community <- substr(long_results$key,2,2)

    g1 <- input$generations+1
    long_results$generation <- rep(1:g1,4)
    long_results$species <- substr(long_results$key,1,1)
    
    lab <- c("Mutualistic","Antagonistic")
    names(lab) <- c("1","2")
    
    glong <- ggplot(data=long_results,aes(x=generation,color=species))+geom_point(aes(y=value))+
      facet_wrap(~community, labeller=labeller(community=lab))+ylim(0,1)+
      ggtitle("Frequency of x1 and y1 in mutualistic and antagonistic communities")
    glong
    
  })
  
})
