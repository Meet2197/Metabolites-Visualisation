# Packages are loaded below. DeSolve provides support for ODE calculation and Plotly is used for innovative graphs.  
library(deSolve)
library(reshape2)
library(nloptr)

# the selected columns you want to work with in metabolites dataset. 

"LH,MH,LV,VSV,I,SL,CE_VL_HDL,TH,VC,CEH,Total_T,Total_C,HC,TV,CE_VL_VLDL"

# this columns are selected sub-Metabolites from Metabolites dataframe from Metabolites.R file. This file have dataframe names metabolites and df3. We took merged data of df3. 
# Mentioned Metabolites have association with exogenous lipo-protein pathway in liver. 

selected_columns <- metabolites %>%
  select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P')

# the mean values for selected columns

ODE_stats <- selected_columns %>%
  summarise(across(.fns = mean))

# Sub-metabolites calculations are depicted below. calculations of sub metabolites has been provided.

H = round((ODE_stats$LH + ODE_stats$MH) , 2)
V = round((ODE_stats$LV + ODE_stats$VSV) , 2)
I = round((ODE_stats$I) , 2)
L = round((ODE_stats$SL) , 2)
C_V = round((ODE_stats$VC) , 2) 
T_V = round((ODE_stats$TGV) , 2)
E_V = round((ODE_stats$CEV) , 2) 
T_H = round((ODE_stats$TGH) , 2)
E_H = round((ODE_stats$CEH) , 2)
C_H = round((ODE_stats$HC) , 2)
T_I = round((ODE_stats$TGI) , 2)
E_I = round((ODE_stats$CEI) , 2)
CI = round((ODE_stats$CI) , 2) 
C = round((ODE_stats$TLCL) , 2) 
T = round((ODE_stats$TLTG) , 2)  
E = round((ODE_stats$TLEC) , 2)
B = round((ODE_stats$AB) , 2)
A1 = round((ODE_stats$A1) , 2)
P = round((ODE_stats$P) , 2)

# equilibrium states of the ODE. 
par <- c(k1=1, k2=2, k4=3, k6=4, k7=5, k8=6, k9=7, alpha=0.6, beta=0.5, gamma=0.7, C, T, E, B, P)

# Define the function to calculate equilibrium
calculate_equilibrium <- function(parameters) {
  with(as.list(parameters), {
    V <- k1 * (C + T + E) * B / k2
    I <- k1 * (C + T + E) * B * alpha / k4
    L <- beta * k1 * (C + T + E) * B * alpha / k6  
    H <- k9 / k8 / gamma 
    A <- k9 / gamma / k7 / (C + P)
    equilibrium_state <- c(V=V, I=I, L=L, H=H, A=A)
    return(equilibrium_state)})} 

# Calculate equilibrium values
equilibrium_values <- calculate_equilibrium(par)

# Create a data frame with the equilibrium values
equilibrium_df <- data.frame(V = equilibrium_values[1], I = equilibrium_values[2],L = equilibrium_values[3],
                             H = equilibrium_values[4],A = equilibrium_values[5])
print(equilibrium_df)

# Initial conditions from equilibrium values
initial_state <- equilibrium_values # initial state of the system using the equilibrium values obtained earlier

# Initial guess for parameters
parameters_guess <- c(1, 1, 1, 1, 1, 1, 1, 1, 1)

# 3.2 we have fit the equilibrium state to the data.
# Modifications :
ode_system <- function(t, state, parameters) { #time t, current state = state, and parameters as inputs 
  with(as.list(c(state, parameters)), {
    # Calculate derivatives
    dVdt <- k1 * (C + T + E) * B - k2 * V
    dIdt <- k2 * alpha * V - k4 * I
    dLdt <- k4 * I * beta - k6 * L
    dHdt <- k7 * A * (C + P) - k8 * H
    dAdt <- k9 - k7 * A * (C + P) + k8 * (1-gamma) * H
    return(list(c(dVdt, dIdt, dLdt, dHdt, dAdt))) #parameters returns the derivatives of state variables with respect to time
  })}

# Initial conditions from equilibrium values
initial_state = calculate_equilibrium(par) # initial state of the system using the equilibrium values obtained earlier

#calculate derivative at initial condition
derivatives_start = ode_system(1,initial_state,par)

# Time points for evaluation
times = seq(0, 0.99, length.out = 9) #  sequence of time points at which the ODE system is evaluated

# Solve the ODE system
solution <- ode(y = initial_state, #  solves the ODE system using the ode function from an appropriate package
                times = times, #passes the initial state, time points, ODE system function, and parameters as arguments to the ode function.
                func = ode_system, 
                parms = par)
print(solution)
# Print the solution
solution_df <- as.data.frame(solution)
 #contains the values of state variables (V, I, L, H, A) at each time point specified in times

# parameters and derivatives test:

test_parameters = c()
derivatives_test = c()

# check derivatives at initial condition for many different random parameter sets
for (i in 1:100){
  par = c(k1=runif(1), k2=runif(1), k4=runif(1), k6=runif(1), k7=runif(1), k8=runif(1), k9=runif(1), alpha=runif(1), beta=runif(1), gamma=runif(1), C=runif(1), T=runif(1), E=runif(1), B=runif(1), P=runif(1))
  test_parameters = rbind(test_parameters,par)
  initial_state = calculate_equilibrium(par) 
  derivatives_start = ode_system(1,initial_state,par)
  derivatives_test = rbind(derivatives_test, derivatives_start[[1]])}

summary(test_parameters) #should follow uniform distribution over interval [0,1]
summary(derivatives_test) #should equal zero up to numerical precision

# 100 different initial guess for parameters guess

# 3.1 Optimisation for constant values
# Define objective wrapper function for optimization
objective_wrapper <- function(params) {
  equilibrium_state <- calculate_equilibrium(c(k1 = params[1], k2 = params[2], k4 = params[3],
                                               k6 = params[4], k7 = params[5], k8 = params[6],
                                               k9 = params[7], alpha = 0.6, beta = 0.5, gamma = 0.7,
                                               C , T , E , B , P))
  sum((equilibrium_state - equilibrium_values)^2)}

# Set optimization options
opts <- list(algorithm = "NLOPT_LN_COBYLA")

# Initialize a list to store the results
optim_results <- list()

# Initial guess for parameters
parameters_guesses <- replicate(100, rep(1, 7), simplify = FALSE) # Initial guess for K1, K2, K4, K6, K7, K8, K9

# Minimize objective function using nloptr for each initial guess
for (i in 1:100) {
  result <- nloptr(x0 = parameters_guesses[[i]], eval_f = objective_wrapper, opts = opts)
  optim_results[[i]] <- result$solution}

# Create a dataframe for the solutions
df_solution <- data.frame(matrix(unlist(optim_results), ncol = 7, byrow = TRUE))
colnames(df_solution) <- c("k1", "k2", "k4", "k6", "k7", "k8", "k9")

# Print the optimized parameters
print(df_solution)

# Iterate over each set of parameters from df_solution and solve the ODE system
solutions_list <- list()  # Initialize a list to store solutions for each set of parameters

for (i in 1:nrow(df_solution)) { # Get parameter values from df_solution
  parameters <- c(k1 = df_solution[i, "k1"], k2 = df_solution[i, "k2"], k4 = df_solution[i, "k4"], 
                  k6 = df_solution[i, "k6"], k7 = df_solution[i, "k7"], k8 = df_solution[i, "k8"], 
                  k9 = df_solution[i, "k9"], alpha = 0.6, beta = 0.5, gamma = 0.7)
  # Initial conditions from equilibrium values
  initial_state <- calculate_equilibrium(parameters)  # Initial state of the system using the equilibrium values obtained earlier
  # Time points for evaluation
  times <- seq(0, 0.99, length.out = 9)  # Sequence of time points at which the ODE system is evaluated
  # Solve the ODE system
  solution <- data.frame(ode(y = initial_state, times = times, func = ode_system, parms = parameters))
  # Store solution in the list
  solutions_list[[i]] <- solution}

# Print solutions
for (i in 1:length(solutions_list)) {
  cat("Solution for parameters set", i, ":\n")
  print(solutions_list[[i]])}