library(shiny)
library(shinydashboard)
library(tidyr)
library(ggplot2)
library(igraph)
library(Matrix)
library(matlib)
library(gganimate)
library(chorddiag)
library(viridis)
library(phaseR)
library(deSolve)

# Define UI for application that draws a histogram



ui <- dashboardPage(
  dashboardHeader(title="Exploration of Coevolutionary Models"),
  dashboardSidebar(
    sidebarMenu(
      menuItem("Introduction to the Application",tabName="Intro_tab",icon=icon("fas fa-info-circle")),
      menuItem("Gavrilets & Hastings 1998", tabName = "Gav_hast", icon=icon("bar-chart-o")),
      menuItem("Nuismer 1999", tabName="Nuismer_99",icon=icon("bar-chart-o")),
      menuItem("Guimarães et al. 2017", tabName="Guimaraes", icon=icon("bar-chart-o"))
    )
  ),
  dashboardBody(
    
    tabItems(
      
      
      tabItem(tabName="Intro_tab",
      h2("Coevolutionary Modeling Exploration"),
      h3("Introduction: The Evolution of Interacting Species"),
      p(h4("Hi, hello, and welcome! This application is built to explore mathematical models of coevoluationary systems.To start
      out, I'll give a quick overview of what we mean when we say coevolution and the difficulties surrounding its study. As John
      Thompson points out in the opening chapter of his book", em('The Coevolutionary Process')," (1994), the complex web of interactions
      within and between species was a central concern of Darwin's thinking on the evolution of biodiversity. Darwin's 
      closing paragraph in",em("On the Origin of Species"), "describes the 'tangled bank,' a vision of life as a dense network
      of interacting organisms, all shaping each other's evolution. The evolution of interacting
      species is a complex problem that continues to fascinate evolutionary biologists, and coevolution is an important
      process within the broader field of community evolution. Coevolution is, in some ways, a sort of 'special case'
      of community evolution where two or more species mutually affect each other's evolutionary trajectories. So, for example, predators may
      drive evolution of toxins as a chemical defense in their prey. These defenses then drive evolution of resistance in 
      predators. As predator resistance increases, prey with higher levels of toxin survive better, thus driving toxin levels
      higher in surviving prey populations, which of course then selects for even higher resistance in predators This example, often called a 'coevolutionary arms
      race' is only one example of how interacting species can shape each other's evolutionary trajectories. In this case, we expect
      coevolutionary forces to drive each species' trait to larger and larger values until some limit is reached. This limit may be
      abiotic (by which we mean it doesn't relate to another living species and is instead defined by environmental factors) or 
      biotic (some other species limits how toxic prey can get or how resistant predators can become). Coevolution can
      occur in antagonistic relationships (like the predator prey example), in mutualistic systems (like flowering plants and 
      pollinating insects), in competitive systems (where the interaction is competition between species for a shared resource), and in systems where
      the nature of the relationship may differ from site-to-site (for example, the", em("Greya-Lithophragma")," relationship discussed in
      Nuismer et al. 1999, one of the models you'll find as part of this web app).")),
      
      h3("The Difficulties of Reciprocal Selection and the Geographic Mosaic Model"),
      p(h4("While it is possible for any of the above relationships to represent coevolution, none of these types of interactions
           are inherently coevolutionary. None are inherently reciprocal. Toxic prey may exert selective pressure on predator resistance, but predator resistance
           may not drive toxicity higher (see Benkman et al. 2019, where resistance in garter snakes correlates to prey (newt)
           toxicity, but prey toxicity is best explained by neutral population structure). While the concept of reciprocal 
           selection can be quite simple, disentangling the complex patterns of how traits are distributed in space into clear, discrete evolutionary
           processes is much more difficult. Even in situations where traits (i.e. phenotypes) are highly correlated among
           interacting species and patterns and mechanisms of reciprocal selection line up perfectly, it remains possible
           for non-reciprocal processes to drive trait distributions.")),
      
      #okay, so next part is digging into the geographic mosaic of coevolution - why that's the 
      #prevailing model now, how that leads to "meaningful" patterns in trait distributions (ie
      #patterns that suggest a strong or diffuse coevolutionary pattern) across species ranges,
      #and the ambiguity that arises when we consider alternative hypotheses - shared abiotic
      #optima for traits, "trait remixing" via demographic factors, one species adapting to
      #the other but no reciprocity involved (in conjunction with spatially structured remixing
      #this could create very similar patterns to spatially structured coevolution from a classic
      #mosaic). 
      
      p(h4("This difficulty of interpreting patterns stems from a number of factors inherent to community evolution and more specifically to coevolution.")),
      
      #Next would be why models are useful - "utility of models." These generate explicit 
      #predictions that can be tested against empirical systems. these can range from tests
      #of whether correlated phenotypes are meaningful to how complex trait distributions
      #can arise from reciprocity to how demographic models may generate specific patterns of
      #genetic differentiaton among differently evolving populations.
      
      h3("The Basics of Coevolution"),
      p(h4())
      ),
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
                                                    value=0.1),
                                        sliderInput("selected_b",
                                                    "Direct, between species interaction strength, spp 1: (b)",
                                                    min = -2,
                                                    max = 2,
                                                    step=0.1,
                                                    value = 0.1),
                                        sliderInput("selected_c",
                                                    "Direct, between species interaction strength, spp 2: (c)",
                                                    min = -2,
                                                    max = 2,
                                                    step = 0.1,
                                                    value=0.1),
                                        sliderInput("selected_d",
                                                    "Indirect, within species interaction strength, spp 2: (d)",
                                                    min = -2,
                                                    max = 2,
                                                    step = 0.1,
                                                    value=0.1),
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
                                                value=0.4),
                                    sliderInput("selected_y1",
                                                "y1:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.4),
                                    sliderInput("selected_x2",
                                                "x2:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.4),
                                    sliderInput("selected_y2",
                                                "y2:",
                                                min=0,
                                                max=1,
                                                step=0.1,
                                                value=0.4)),
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
              ),
      tabItem(tabName="Guimaraes",

                chorddiagOutput("Guim_net", height = 900),
              
              fluidRow(
                plotOutput("Guim_eq")
              ),
              
              fluidRow(
                shinydashboard::box(title="Plot Controls",
                                    
                                    
                                    uiOutput("slider")
                                    
                                    
                                    # sliderInput(inputId = "guim_time",
                                    #             label="Timestep",
                                    #             min = 1, 
                                    #             max=1, 
                                    #             value=1,
                                    #             step=1,
                                    #             )),
              ))
              
              
            )
      )
    )
  )



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

r <- read.table("Bcor.txt")

adaptive_peak <- function(matriz_row,zs,theta){
  qs_that_matter <- matriz_row[which(matriz_row!=0)]
  zs_that_matter <- zs[which(matriz_row!=0)]
  muts <- sum(qs_that_matter*zs_that_matter)
  abio <- (1-sum(qs_that_matter)*theta)
  peak <- muts+abio
  return(peak)
}

eq_z <- function(q_mat, phi_mat_diag, theta, dims){
  
  I_mat <- diag(dims)
  
  temp <- solve(I_mat-q_mat)
  
  
  z_star <- temp%*%(phi_mat_diag%*%theta)
  
  
  return(z_star)
  
}


guim_sim <- function(incidence_mat=NULL, #an incedence matrix for a bipartite graph
                       alpha=0.2,  #scaling parameter for trait matching, larger is more sensitive, lower is less sensitive.
                       timesteps=6000, #maximum number of timesteps
                       th_medio=10, #maximal value for abiotic selection, this uses a uniform distro from 0:1, so I'm not letting users set a bound for that.
                       m_medio=0.7, #mean for mutualistic strength
                       h_medio=0.2, #product of additive genetic variance and slope of selection gradient. Basically controls how quickly evolution can move
                       h_sd=1){ #standard deviation for generating random h parameters for each species. normally distributed
    
    
    if(is.null(incidence_mat)){
      r <- read.table("Bcor.txt")
    }else{
      r <- incidence_mat
    }
    
    
    lines <- dim(r)[1] #number of, we'll say, plants?
    columns <- dim(r)[2] #number of, we'll say, pollinators?
    
    dims <- lines+columns #total number of species in network
    
    F1 <- matrix(0,nrow=dims,ncol=dims)
    
    
    
    #fill q matrix (F1) with binary interactions
    F1[1:lines,lines+1:columns] <- as.matrix(r)
    F1[lines+1:columns,1:lines] <- as.matrix(r)
    F1[F1>0] <- 1
    
    #generate our random h's
    h <- h_medio+0.01*rnorm(dims,0,h_sd)
    
    #this is just a lil while loop that redoes the h vector above if any of the values are negative or above 1.
    

    while(min(h)<0||max(h)>1){
      h <- h_medio+0.01*rnorm(dims,0,h_sd)
    }
    
    #gimme them thetas
    theta <- th_medio*runif(dims,0,1)
    
    #generate abiotic selection values
    p_medio <- 1-m_medio
    pvec <- p_medio+0.01*rnorm(dims,0,h_sd)
    

    while(min(pvec)<0||max(pvec)>1){
      pvec <- p_medio+(0.01*rnorm(dims,0,1))
    }
    

    z = th_medio*runif(dims,0,1)
    
    zn <- z
    
    step <- 0
    matriz <- F1
    tau <- timesteps
    
    one_mat <- matrix(1,nrow=dims,ncol=dims)
    int_mat <- diag(1-pvec)
    
    
    peaks_z_list <- list()
    
    z_star_frame <- data.frame(matrix(0,nrow=dims, ncol=1))
    
    
    #MODEL
    for(iter in 1:timesteps){
      peaks <- numeric(dims)
      #check is simulation is running based on step value
      if(step ==0){
        #loop through rows
        for(i in 1:dims){
          #loop through columns
          for(j in 1:dims){
            
            matriz[i,j] <- F1[i,j]*exp(-alpha*((z[i]-z[j])^2))
            
            
          }
        }
        for(i in 1:dims){
          matriz[i,] <- matriz[i,]/sum(matriz[i,])
          
          peaks[i] <- adaptive_peak(matriz[i,],z,theta[i])
          
        }
        
        
        
        
        #ABOVE - the two parts above are equation 2 from the methods in paper. First part is numerator, second is denom. We use F1 as the m
        #parameter, which just indicates whether a link exists, and thus whether we should calculate selection pressure q for species i and j.
        #exp and -alpha are pretty standard, but basically we calculate the squared diff between the two traits, multiply by scaling constant
        #that shapes sensitivity of trait matching, then exponentiate. This gives us a roughly normally distributed 0:1 value that decreases towards
        #0 as diff increases and approaches 1 as traits match.
        
        #second part is just normalizing those links to sum to 1 for each species, so basically each row has to sum to 1.
        
        #honestly...not sure what this does. Int mat is a diagonal of 1-pvec, so basically a diagonal matrix describing the complement to 
        matriz <- int_mat%*%matriz
        
        
        zmatriz <- diag(z)
        
        zmatriz <- (one_mat%*%zmatriz)-(zmatriz%*%one_mat)
        
        echanges <- matriz*zmatriz
        
        echanges <- rowSums(echanges)
        
        z2 <- z
        
        for(i in 1:dims){
          #this is just the non-zi component of eq 1 in paper. Let's try and walk through. So h is the vector of phi parameters. That's additive
          #genetic variance and the slope of the selection gradient - compound var? qij from the eq is echanges here I think. Then in the eq
          #it's 1-qij (so presumably the p vec here is that??), times that phi param, times the diff between abiotic opt and the trait. So this computes
          #just the change based on additive variance (so what's genetically feasible), the changes from our interactions, and the changes from
          #abiotic. p vec is our non-mutualistic selection, here taken as enviornmental selection towards some optimum, as opposed to say, a further
          #interaction shaping that trait - like predation or a different mutualism, etc.
          
          zn[i] <- (h[i]*echanges[i])+(pvec[i]*h[i]*(theta[i]-z[i]))
          
          #this is the final part of eq 1 in paper, zi + the rest of that equation
          z[i] <- z[i]+zn[i]
        }
        
        #is change occuring?
        peaks <- as.data.frame(peaks)
        peaks$type <- c(rep("Plant",lines),rep("Pollinator",columns))
        peaks$time <- iter
        peaks$index <- 1:dims
        peaks$z <- z
        peaks$z_star <- eq_z(matriz,diag(pvec), theta, dims)
        peaks_z_list[[iter]] <- peaks
        d <- z2-z
        
        if(mean(abs(d))<0.001){
          tau <- iter
          step <- 1
        }
        
        
        
        
      }
      
    }
    
    big_frame <- do.call("rbind", peaks_z_list)  
    return(list(Q_Mat=matriz, peaks_z=big_frame, r=r,tau=tau))
  }


# Define server logic required to draw a histogram
server <- shinyServer(function(input, output, session) {
  
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
  
  output$Guim_net <- renderChorddiag({
    r <- as.matrix(r)
    
    plant_stuff <- paste("Plant",1:33)
    pol_stuff <- paste("Pollinator", 1:25)
    dimnames(r) <- list(Plants = plant_stuff,
                        Pollinators = pol_stuff)
    
    pl <- viridis(33)
    pl2 <- plasma(25)
    
    chorddiag(r, type="bipartite", groupColors = c(pl,pl2), showTicks = F, groupnameFontsize = 10)
  })
  
  
  results <- reactive(guim_sim())
  
  
  
  output$slider <- renderUI({
    sliderInput("Guim_time",
                "Timestep", 
                min = 1, 
                max=results()$tau,
                value=1, 
                step=1,
                animate=animationOptions(interval=750,loop=T,playButton = "Play"))
  })
  
  output$Guim_eq <- renderPlot({


   z_w_stars <-  results()$peaks_z
   
   long_sim <- pivot_longer(z_w_stars,c(z,z_star),names_to = "origin", values_to = "data")
   
   gif2 <- ggplot(data = long_sim[long_sim$time==input$Guim_time,], aes(x=index, y=data))+geom_point(size=3,aes(colour=origin, shape=type))+
     geom_line(aes(group=index))+ylim(0,10)
   gif2
   
  })
  
  
})
