library(deSolve)
library(ggplot2)
library(reshape2)

data <- data <- read_excel("C:/Users/user/OneDrive - London School of Hygiene and Tropical Medicine/Preethi_MSc_CHIK/data.xlsx")

age <- 1:100
sero_data <- as.data.frame(matrix(NA, nrow = 100, ncol =2))
sero_data[,1] <- 1:100
sero_data[,2] <- 1 - exp(-data[1,]$lambda*age)
colnames(sero_data) <- c("age", "seropositive")
sero_age70 <- sero_data[70,]
cum_sero_epi1 <- sero_age70 / (70/8)
N <- data[1,]$`Total Population`
N_cum_inc <- cum_sero_epi1 * N

sir_model <- function(time, state, parameters) { 
  with(as.list(c(state, parameters)), {
    N <- S + I + R
    dS <- -parameters['beta'] * (I / N) * S
    dI <- parameters['beta'] * (I / N) * S - parameters['gamma'] * I
    dR <- parameters['gamma'] * I
    return(list(c(dS, dI, dR)))
  })
}

beta_calibration <- function(beta, data, cum_inc = 9.8) {
  
  parameters <- c(beta = beta, gamma = 1/7) # 7 days duration of infectious
  
  # Initial values
  initial_values <- c(S = data[1,]$`Initial Susceptible` - 1,
                      I = 1, 
                      R = data[1,]$`Initial Recovered`)
  
  # Time
  times <- seq(from = 0, to = 180, by = 1) # assume 180 days outbreak with 1 day each step
  
  # Run the model
  output <- as.data.frame(ode(y = initial_values, 
                              times = times, 
                              func = sir_model,
                              parms = parameters))
  
  # Calculate the final cumulative incidence
  final_incidence <- (output[nrow(output), "R"] - initial_values["R"]) / 
    (initial_values["S"] + initial_values["I"] + initial_values["R"]) * 100
  
  return(abs(final_incidence - cum_inc)) # final incidence should be close to 9.8% (one single outbreak = 9.8% assuming 8 yrs of inter-epidemic period)
}

# Use optim to find the beta value that makes cumulative incidence close to 9.8%
result <- optimize(f = beta_calibration,
                   interval = c(0.01, 5), 
                   data = data,
                   cum_inc = 9.8)

(best_beta <- result$minimum)

# Run the model with the optimized beta
parameters <- c(beta = best_beta, gamma = 1/7)

times <- seq(from = 0, to = 180, by = 1)

initial_values <- c(S = data[1,]$`Initial Susceptible` - 1,
                    I = 1, 
                    R = data[1,]$`Initial Recovered`)

output <- as.data.frame(ode(y = initial_values, 
                            times = times, 
                            func = sir_model,
                            parms = parameters))

# Plotting the results
melt(as.data.frame(output), id = "time") %>% 
  ggplot(aes(x=time, y=value, color=variable, shape=variable)) +
  geom_line() +
  geom_jitter(show.legend = FALSE) +
  scale_color_manual(values = c("blue","red","green")) + 
  scale_shape_manual(values = c(0, 4, 1)) +
  xlab("Time (days)") +
  ylab("Number of people") +
  #labs(title=paste("gamma of", parameters['gamma'], 
  #                 " beta of", parameters['beta']),
  #     color = "Compartment") +
  theme_light()