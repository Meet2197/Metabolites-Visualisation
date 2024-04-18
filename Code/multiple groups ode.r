# selected metabolite values :

metabolites2 <- metabolites %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P')

ODE_stats2 <- No_Meds %>% summarise(across(.fns = mean))

# the mean values for selected columns
H = round((ODE_stats2$LH + ODE_stats2$MH) , 2)
C = round((ODE_stats2$TLCL) , 2) 
T = round((ODE_stats2$TLTG) , 2)  
E = round((ODE_stats2$TLEC) , 2)
P = round((ODE_stats2$P) , 2)
A1 = round((ODE_stats2$A1) , 2)
B = round((ODE_stats2$AB) , 2)


#Different group mean values:

Healthy <- merge(metabolites, Healthy_Patients, by = "eid", all.x = TRUE) %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P')
No_Meds <- intersect(metabolites$eid, No_Medication$eid)
MASLD1 <- merge(metabolites, MASLD, by = "eid", all.x = TRUE) %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P')
MASH1 <- merge(metabolites, MASH, by = "eid", all.x = TRUE) %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P') 
metformin1 <- merge(metabolites, metformin, by = "eid", all.x = TRUE) %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P') 
pioglitazone1 <- merge(metabolites, pioglitazone, by = "eid", all.x = TRUE) %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P') 
ramipril1 <- merge(metabolites, ramipril, by = "eid", all.x = TRUE) %>% select('LH', 'MH', 'LV', 'VSV', 'I', 'SL', 'TGH', 'TGL', 'TGV', 'TGI', 'VC', 'HC', 'LC', 'CI', 'CEV', 'CEL', 'CEH', 'CEI', 'TLEC', 'TLTG', 'TLCL', 'AB','A1','P') 

# different mean values:

No_MedsM <- colMeans(No_Meds, na.rm = TRUE)
MASLD1M <- colMeans(MASLD1, na.rm = TRUE)
MASH1M <- colMeans(MASH1, na.rm = TRUE)
metformin1M <- colMeans(metformin1, na.rm = TRUE)
pioglitazone1M <- colMeans(pioglitazone1, na.rm = TRUE)
ramipril1M <- colMeans(ramipril1, na.rm = TRUE)

# equilibrium constant :

calculate_equilibrium1 <- function(parameters) {
  with(as.list(c(parameters)), {
    V <- k1 * (C + T + E) * B1 / k2
    I <- k1 * (C + T + E) * B1 * alpha / k4
    L <- beta * k1 * (C + T + E) * B * alpha / k6  
    H <- k9 / k8 / gamma 
    A <- k9 / gamma / k7 / (C + P)
    equilibrium_state <- c(V=V, I=I, L=L, H=H, A=A)
    return(equilibrium_state)})} 

# Calculate equilibrium values for No_Meds_mean

par1 <- c(k1=1, k2=2, k4=3, k6=4, k7=5, k8=6, k9=7, alpha=0.6, beta=0.5, gamma=0.7, C, T, E, B, P)

# Calculate equilibrium values
equilibrium_values1 <- calculate_equilibrium1(par1)
equilibrium_values1 <- calculate_equilibrium1(par1, initial_conditions = list(C1 = C1, T1 = T1, E1 = E1, B1 = B1, P1 = P1))

# Create a data frame with the equilibrium values
equilibrium_df1 <- data.frame(V = equilibrium_values1[1], I = equilibrium_values1[2],L = equilibrium_values1[3],
                             H = equilibrium_values1[4],A = equilibrium_values1[5])
print(equilibrium_df1)

# Initial conditions from equilibrium values
initial_state <- equilibrium_values # initial state of the system using the equilibrium values obtained earlier
