# Packages are loaded below. DeSolve provides support for ODE calculation and Plotly is used for innovative graphs.  

library(deSolve)
library(plotly)

# the selected columns you want to work with in metabolites dataset. 

"L_HDL,M_HDL,L_VLDL,VS_VLDL,IDL,S_LDL,CE_VL_HDL,TG_HDL,VLDL_C,CE_HDL,Total_TG,Total_C,HDL_C,TG_VLDL,CE_VL_VLDL"

# this columns are selected sub-Metabolites from Metabolites dataframe from Metabolites.R file. This file have dataframe names metabolites and df3. We took merged data of df3. 
# Mentioned Metabolites have association with exogenous lipoprotein pathway in liver. 

selected_columns <- df3 %>%
  select('L_HDL', 'M_HDL', 'L_VLDL', 'VS_VLDL', 'IDL', 'S_LDL', 'TG_HDL', 'TG_LDL', 'TG_VLDL', 'TG_IDL', 'VLDL_C', 'HDL_C', 'LDL_C', 'C_IDL', 'CE_VLDL', 'CE_LDL', 'CE_HDL', 'CE_IDL', 'Total_Esterified_C', 'Total_TG', 'Total_C', 'Apo_B','Apo_A1','Phospholipids_in_HDL')

# Calculate the mean values for selected columns

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
C = round((ODE_stats$Total_C) , 2) 
T = round((ODE_stats$Total_TG) , 2)  
E = round((ODE_stats$Total_Esterified_C) , 2)
B = round((ODE_stats$Apo_B) , 2)
A1 = round((ODE_stats$Apo_A1) , 2)
P = round((ODE_stats$Phospholipids_in_HDL) , 2)

# ODE system df Extracts initial state values from ODE_table. This ODE system dataframe calculates ODE equation written. 
# This equation provides calculations of VLDL, IDL, LDL, HDL, Apo A1, Apo B metabolits for VLDL synthesis.   
 
ode_system <- function(t, state, params) {
  with(as.list(c(state, params)), {
    dVdt <- k1 * (C_V + T_V + E_V) * B - k2 * V
    dIdt <- k3 * 0.60 * V - k4 * I
    dLdt <- k5 * I * 0.5 - k6 * L
    dHdt <- k7 * 0.7 * A1 * (C_H + P) - k8 * 0.3 * A1 * H
    dAdt <- k7 * 0.7 * A1 * (C_H + P) + k9 * H 
    dC_Vdt <- 0  # Add equations for the other variables
    dT_Vdt <- 0
    dE_Vdt <- 0
    dC_Hdt <- 0
    dT_Hdt <- 0
    dE_Hdt <- 0
    dBdt <- 0
    dPdt <- 0
    return(list(c(dVdt, dIdt, dLdt, dHdt, dAdt, dC_Vdt, dT_Vdt, dE_Vdt, dC_Hdt, dT_Hdt, dE_Hdt, dBdt, dPdt)))
  })
}

# Set initial conditions(from metabolites dataframe) and parameters(constant values)

initial_state <- c(H , V , I , L , A1, C_V , T_V , E_V , C_H , T_H , E_H , B , P)
parameters <- c(k1 = 1.0, k2 = 1.0, k3 = 1.0, k4 = 1.0, k5 = 1.0, k6 = 1.0, k7 = 1.0, k8 = 1.0, k9 = 1.0 )


# Define the time points for evaluation

times <- c(0)

# Set hmax to a non-negative value (e.g., 0.1)

hmax <- 0.1

# Solve the ODE system

solution <- ode(y = initial_state, times = times, func = ode_system, parms = parameters, hmax = hmax)
  
print(solution)

# Initialize a data frame to store derivatives

Metabolism_df <- data.frame(Time = times, dVdt = numeric(length(times)),
                             dIdt = numeric(length(times)), 
                             dLdt = numeric(length(times)), 
                             dHdt = numeric(length(times)), 
                             dAdt = numeric(length(times)))

# Calculate derivatives for each time point

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

# Create a data frame for plotting

meatbolism_plot <- data.frame(Time = times,
                        dVdt = Metabolism_df$dVdt,
                        dIdt = Metabolism_df$dIdt,
                        dLdt = Metabolism_df$dLdt,
                        dHdt = Metabolism_df$dHdt,
                        dAdt = Metabolism_df$dAdt)

# Reshape the data frame into long format for ggplot

meatbolism_plot1 <- tidyr::pivot_longer(meatbolism_plot, cols = -Time, names_to = "Variable", values_to = "Value")

# Create a grouped bar plot for Metabolites secretion by using ggplot function. 

ggplot(meatbolism_plot1, aes(x = Variable, y = Value, fill = Variable)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_text(aes(label = round(Value, 2)), vjust = -0.5, position = position_dodge(width = 0.9)) +
  labs(title = " ODE Metabolites Secretion", y = "Value(IQR)") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# Particle size distribution

particle_size_distribution <- function(t, state, parameters1) {
  dVdt <- parameters1$k1 * (state$C_V + state$T_V + state$E_V) * state$B - parameters1$k2 * state$V
  dIdt <- parameters1$k3 * 0.57 * dVdt - parameters1$k4 * state$I
  dLdt <- parameters1$k5 * dIdt * 0.5 - parameters1$k6 * state$L * 0.5
  dHdt <- parameters1$k7 * 0.7 * state$A1 * (state$C_H + state$P) - parameters1$k8 * 0.3 * state$A1 * state$H
  dAdt <- parameters1$k7 * 0.7 * state$A1 * (state$C_H + state$P) + parameters1$k9 * dHdt
  
# Coefficients for the particle size distribution

  alpha1 <- parameters1$alpha1
  alpha2 <- parameters1$alpha2
  alpha3 <- parameters1$alpha3
  alpha4 <- parameters1$alpha4
  alpha5 <- parameters1$alpha5
  alpha6 <- parameters1$alpha6
  
  dPSDdt <- alpha1 * 0.57 * dVdt + alpha2 * dIdt * 0.5 + alpha3 * dLdt * 0.5 + alpha4 * 0.7 * state$A1 - alpha5 * 0.3 * state$A1 + alpha6 * 0.7 * state$A1
  
  return(list(c(dVdt, dIdt, dLdt, dHdt, dAdt, dPSDdt)))
}

# Example usage

initial_state <- list(
  V = round((ODE_stats$L_VLDL + ODE_stats$VS_VLDL), 2),
  I = round(ODE_stats$IDL, 2),
  L = round(ODE_stats$S_LDL, 2),
  C_V = round(ODE_stats$VLDL_C, 2),
  T_V = round(ODE_stats$TG_VLDL, 2),
  E_V = round(ODE_stats$CE_VLDL, 2),
  T_H = round(ODE_stats$TG_HDL, 2),
  E_H = round(ODE_stats$CE_HDL, 2),
  C_H = round(ODE_stats$HDL_C, 2),
  T_I = round(ODE_stats$TG_IDL, 2),
  E_I = round(ODE_stats$CE_IDL, 2),
  C_I = round(ODE_stats$C_IDL, 2),
  C = round(ODE_stats$Total_C, 2),
  T = round(ODE_stats$Total_TG, 2),
  E = round(ODE_stats$Total_Esterified_C, 2),
  B = round(ODE_stats$Apo_B, 2),
  A1 = round(ODE_stats$Apo_A1, 2),
  P = round(ODE_stats$Phospholipids_in_HDL, 2)
)

parameters1 <- list(k1 = 1, k2 = 1, k3 = 1, k4 = 1, k5 = 1, k6 = 1, k7 = 1, k8 = 1, k9 = 1, alpha1 = 1, alpha2 = 1, alpha3 = 1, alpha4 = 1, alpha5 = 1, alpha6 = 1)

result <- particle_size_distribution(0, initial_state, parameters1)

dPSDdt_values <- result$dPSDdt

# Create a data frame for plotting

plot_data <- data.frame(PSD = dPSDdt_values)

# Create the time series plot without the time variable

ggplot(data = data.frame(PSD = result$dPSDdt), aes(y = PSD)) +
  geom_line() +
  labs(y = "dPSD/dt") +
  ggtitle("Particle Size Distribution Over Time")
