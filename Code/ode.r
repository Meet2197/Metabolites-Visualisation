# Packages are loaded below. DeSolve provides support for ODE calculation and Plotly is used for innovative graphs.  

library(deSolve)
library(plotly)
library(minpack.lm)


# the selected columns you want to work with in metabolites dataset. 

"L_HDL,M_HDL,L_VLDL,VS_VLDL,IDL,S_LDL,CE_VL_HDL,TG_HDL,VLDL_C,CE_HDL,Total_TG,Total_C,HDL_C,TG_VLDL,CE_VL_VLDL"

# this columns are selected sub-Metabolites from Metabolites dataframe from Metabolites.R file. This file have dataframe names metabolites and df3. We took merged data of df3. 
# Mentioned Metabolites have association with exogenous lipo-protein pathway in liver. 

selected_columns <- metabolites %>%
  select('L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Tl_Esterified_C', 'Tl_TG', 'Tl_C', 'Apo_B','Apo_A1','P_HDL')

# the mean values for selected columns

ODE_stats <- selected_columns %>%
  summarise(across(.fns = mean))

# Sub-metabolites calculations are depicted below. calculations of sub metabolites has been provided.

H = round((ODE_stats$L_HDL + ODE_stats$M_HDL) , 2)
V = round((ODE_stats$L_VLDL + ODE_stats$VS_VLDL) , 2)
I = round((ODE_stats$IDL) , 2)
L = round((ODE_stats$S_LDL) , 2)
C_V = round((ODE_stats$VLDL_C) , 2) 
T_V = round((ODE_stats$TG_VLDL) , 2)
E_V = round((ODE_stats$CE_VLDL) , 2) 
T_H = round((ODE_stats$TG_HDL) , 2)
E_H = round((ODE_stats$CE_HDL) , 2)
C_H = round((ODE_stats$HDL_C) , 2)
T_I = round((ODE_stats$TG_IDL) , 2)
E_I = round((ODE_stats$CE_IDL) , 2)
C_I = round((ODE_stats$C_IDL) , 2) 
C = round((ODE_stats$Tl_C) , 2) 
T = round((ODE_stats$Tl_TG) , 2)  
E = round((ODE_stats$Tl_Esterified_C) , 2)
B = round((ODE_stats$Apo_B) , 2)
A1 = round((ODE_stats$Apo_A1) , 2)
P = round((ODE_stats$P_HDL) , 2)

# ODE system df Extracts initial state values from ODE_table. This ODE system dataframe calculates ODE equation written. 
# This equation provides calculations of VLDL, IDL, LDL, HDL, Apo A1, Apo B metabolits for VLDL synthesis.   
 
ode_system <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dVdt <- k1 * (C + T + E) * B - k2 * V
    dIdt <- k3 * 0.60 * V - k4 * I
    dLdt <- k5 * I * 0.5 - k6 * L
    dHdt <- k7 * 0.7 * A1 * (C + P) - k8 * 0.3 * A1 * H
    dAdt <- k9 - k7 * 0.7 * A1 * (C + P) + k8 * (0.3) * H
    dC_Vdt <- 0  # Add equations for the other variables
    dT_Vdt <- 0
    dE_Vdt <- 0
    dC_Hdt <- 0
    dT_Hdt <- 0
    
    return(list(c(dVdt, dIdt, dLdt, dHdt, dAdt, dC_Vdt, dT_Vdt, dE_Vdt, dC_Hdt, dT_Hdt)))
  })
}

# initial conditions(from metabolites dataframe) and parameters(constant values)

initial_state <- c( C , T , E , V , I , L , B , A1 ,H, P)
parameters <- c(k1 = 1.5, k2 = 2.0, k3 = 1.2, k4 = 1.8, k5 = 1.3, k6 = 1.7, k7 = 1.4, k8 = 1.9, k9 = 1.1, k10 = 0)

# time points for evaluation

times <- c(0)

# Set hmax to a non-negative value (e.g., 0.1)

hmax <- 0.1

#  ODE system

solution <- ode(y = initial_state, times = times, func = ode_system, parms = parameters, hmax = hmax)
  
# Initialize Metabolism_df to store derivatives

Metabolism_df <- data.frame(Time = times, dVdt = numeric(length(times)),
                             dIdt = numeric(length(times)), 
                             dLdt = numeric(length(times)), 
                             dHdt = numeric(length(times)), 
                             dAdt = numeric(length(times)))

# derivatives for each time point

for (i in 1:length(times)) {
  derivatives <- tryCatch(
    {
      ode_system(times[i], initial_state, parameters)
    },
    warning = function(w) {
      cat("Warning:", w$message, "\n")
      NULL
    }
  )
  
  if (!inherits(derivatives, "try-error")) {
    Metabolism_df[i, 2:6] <- unlist(derivatives)
  }
}

print(Metabolism_df)

# meatbolism_plot for plotting

meatbolism_plot <- data.frame(Time = times,
                        dVdt = Metabolism_df$dVdt,
                        dIdt = Metabolism_df$dIdt,
                        dLdt = Metabolism_df$dLdt,
                        dHdt = Metabolism_df$dHdt,
                        dAdt = Metabolism_df$dAdt)

# Reshape the meatbolism_plot into long format for ggplot

meatbolism_plot1 <- tidyr::pivot_longer(meatbolism_plot, cols = -Time, names_to = "Variable", values_to = "Value")

# grouped bar plot for Metabolites secretion by using ggplot function. 

ggplot(meatbolism_plot1, aes(x = Variable, y = Value, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(Value, 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(title = " ODE Metabolites Secretion", y = "Value(Mean)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

#  NLS function with ODE system:
nls_function <- function(C, T, E, V, I, L, B, A1, H, P, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10) {
  ode_system <- function(t, state, params) {
    with(as.list(c(state, params)), {
      dVdt <- k1 * (C + T + E) * B - k2 * V
      dIdt <- k3 * 0.60 * V - k4 * I
      dLdt <- k5 * I * 0.5 - k6 * L
      dHdt <- k7 * 0.7 * A1 * (C + P) - k8 * 0.3 * A1 * H
      dAdt <- k9 - k7 * 0.7 * A1 * (C + P) + k8 * (0.3) * H
      dC_Vdt <- 0  # Add equations for the other variables
      dT_Vdt <- 0
      dE_Vdt <- 0
      dC_Hdt <- 0
      dT_Hdt <- 0
      
      return(list(c(dVdt, dIdt, dLdt, dHdt, dAdt, dC_Vdt, dT_Vdt, dE_Vdt, dC_Hdt, dT_Hdt)))
    })
   }
  
  initial_state <- c( C , T , E , V , I , L , B , A1 ,H, P)
  parameters <- c(k1 = 1.0, k2 = 1.0, k3 = 1.0, k4 = 1.0, k5 = 1.0, k6 = 1.0, k7 = 1.0, k8 = 1.0, k9 = 1.0 , k10 = 0)
  
  # time points for evaluation
  
  times <- c(0)
  
  # hmax to a non-negative value
  
  hmax <- 0.1
  
  # ODE system
  
  solution <- ode(y = initial_state, times = times, func = ode_system, parms = parameters, hmax = hmax)
}

# equilibrium states of the ODE

calculate_equilibrium <- function(k1, k2, k4, k6, k7, k8, k9, C, T, E, ApoB, P, HDL) {
  V <- k1 * (C + T + E) * ApoB / k2
  I <- k1 * (C + T + E) * ApoB * 0.60 / k4
  L <- k4 * ApoB * I / (k6 * 0.5)
  H <- k9 / (k8 * 0.7)
  A <- k9 / (k7 * (C + P))
  
  equilibrium_state <- c(V, I, L, H, A)
  # Filter out negative equilibrium states
  positive_equilibrium <- equilibrium_state[equilibrium_state >= 0]
  
  return(positive_equilibrium)
}

# Example usage:
parameters <- list(k1 = 1.0, k2 = 1.0, k4 = 1.0, k6 = 1.0, k7 = 1.0, k8 = 1.0, k9 = 1.0,
                   C, T , E , B, P , H)

equilibrium <- calculate_equilibrium(parameters$k1, parameters$k2, parameters$k4, parameters$k6,
                                     parameters$k7, parameters$k8, parameters$k9,
                                     parameters$C, parameters$T, parameters$E, parameters$B,
                                     parameters$P, parameters$L)
print(equilibrium)

# calculate the difference between measured metabolites and equilibrium state
equilibrium_function <- function(k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, C, T, E, B, A1, P, H) {
  V <- k1 * (C + T + E) * B / k2
  I <- k1 * (C + T + E) * B * 0.60 / k4
  L <- k4 * B * I / (k6 * 0.5)
  H <- k7 * 0.7 * A1 * (C + P) / k8 * 0.3 * A1
  A <- k9 / k8 * (0.7)
  
  return(c(V, I, L, H, A1, P))
}

# Define the objective function
objective_function <- function(parameters, data) {
  predicted_values <- equilibrium_function(parameters[1], parameters[2], parameters[3], parameters[4],
                                           parameters[5], parameters[6], parameters[7], parameters[8],
                                           parameters[9], parameters[10], data[1], data[2], data[3], 
                                           data[4], data[5], data[6])
  # Calculate sum of squared errors
  error <- sum((predicted_values - data)^2)
  
  return(error)
}

# Initial parameter guesses

initial_guess <- c(k1 = 0.5, k2 = 0.5, k3 = 0.5, k4 = 0.5, k5 = 0.5, k6 = 0.5, k7 = 0.5, k8 = 0.5, k9 = 0.5)

# Fit the model using nls
fit <- nls(ODE_stats ~ equilibrium_function(k1, k2, k3, k4, k5, k6, k7, k8, k9, C, T, E, B, A1, P),
           start = initial_guess, algorithm = "port", lower = rep(0, length(initial_guess)),
           upper = rep(Inf, length(initial_guess)), control = list(maxiter = 1000))

# Get the fitted parameters
fitted_parameters <- coef(fit)
