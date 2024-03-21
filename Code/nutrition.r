# readcsv :

nutritrion <- read.csv(file = 'C:/Users/User/Desktop/Data/nutritrion.csv')

# setnames of columns :

setnames(nutritrion, old = colnames(nutritrion), new = c('eid', 'meals_typ', 'alcohol_fr','alcohol_cnsm', 'Diet_Var', 'Fruit','Beef',
                                                         'Phys_act','Pr_meat','No_edws','Sleep','oil_typ','Work_sat','Miserable'))

# Nutrition associated analysis :

ALL_Nutrition <- merge(nutritrion, ALL_Death, by.x="eid", by.y="eid", all.x = TRUE, all.y = TRUE)  %>%
  select('eid', 'Age_AC', 'Gender','BMI','Smoking', 'Drinking', 'Qualifications', 'Diabetes','liverdisease','death','meals_typ', 'alcohol_fr',
         'alcohol_cnsm', 'Diet_Var', 'Fruit','Beef','Phys_act','Pr_meat','No_edws','Sleep','oil_typ','Miserable','Work_sat')

ALL_Nutrition[is.na(ALL_Nutrition)] <- 0

# Type of meals:
meals <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$meals_typ %in% c(461, 462, 463, 464, 465), ]
mealcounts <- table(meals$meals)
mostmeals <- head(sort(counts, decreasing = TRUE),5)  # Sort counts in descending order and take top 5 values
maxmeals <- meals[meals$liverdisease == 1 & meals$meals %in% mostmeals, ] # Filter dataframe to include only top 5 values
meal_per <- prop.table(mostmeals) * 100
meals_df <- data.frame(Value = as.character(names(mostmeals)),Count = as.numeric(mostmeals),Percentage = as.numeric(meal_per))

# Diet variation: 

mvar <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$Diet_Var %in% c(1, 2, 3), ] # numbers are type of variation
mvar$Diet_Var <- as.character(mvar$Diet_Var) # change of data class
mvarcounts <- table(mvar$Diet_Var) # counts of variation
mostmvar <- head(sort(mvarcounts, decreasing = TRUE),3)  # Sort counts in descending order and take top 5 values
maxmvar <- mvar[mvar$liverdisease == 1 & mvar$Diet_Var %in% mostmvar, ] # Filter dataframe to include only top 5 values
mvar_per <- prop.table(mostmvar) * 100 # pcentage of meal variation
mvar_df <- data.frame(Value = as.character(names(mostmvar)),Count = as.numeric(mostmvar),Percentage = as.numeric(mvar_per)) # percentage of meal var. 

# beer intake :
beef_int <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$Beef %in% c(0, 1, 2, 3, 4), ] # numbers are type of variation
beef_int$Beef <- as.character(beef_int$Beef) # change of data class
beef_intcounts <- table(beef_int$Beef) # counts of variation
mostbeef_int <- head(sort(beef_intcounts, decreasing = TRUE),5)  # Sort counts in descending order and take top 5 values
maxbeef <- beef_int[beef_int$liverdisease == 1 & beef_int$Beef %in% mostbeef_int, ] # Filter dataframe to include only top 5 values
mbeef_per <- prop.table(mostbeef_int) * 100 # pcentage of meal variation
beef_df <- data.frame(Value = as.character(names(mostbeef_int)),Count = as.numeric(mostbeef_int),
                      Percentage = as.numeric(mbeef_per)) # percentage of meal var. 


# physical activity  :
phy_int <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$Phys_act %in% c(1, 2, 3, 4, 5, -7), ] # numbers are type of variation
phy_int$Phys_act <- as.character(phy_int$Phys_act) # change of data class
Physact_c <- table(phy_int$Phys_act) # counts of variation
mostPhysact <- head(sort(Physact_c, decreasing = TRUE),6)  # Sort counts in descending order and take top 6 values
maxphyact <- phy_int[phy_int$liverdisease == 1 & phy_int$Phys_act %in% mostPhysact, ] # Filter dataframe to include only top 5 values
mphy_per <- prop.table(mostPhysact) * 100 # percentage of meal variation
mphy_df <- data.frame(Value = as.character(names(mostPhysact)),Count = as.numeric(mostPhysact),
                      Percentage = as.numeric(mphy_per)) # percentage of meal var.

# alcohol intake frequency :
alcohol_fr <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$alcohol_fr %in% c(1, 2, 3, 4, 5, 6), ] # numbers are type of variation
alcohol_fr$alcohol_fr <- as.character(alcohol_fr$alcohol_fr) # change of data class
alcohol_c <- table(alcohol_fr$alcohol_fr) # counts of variation
mostalcoholfr <- head(sort(alcohol_c, decreasing = TRUE),6)  # Sort counts in descending order and take top 6 values
maxalcoholfr <- alcohol_fr[alcohol_fr$liverdisease == 1 & alcohol_fr$alcohol_fr %in% mostalcoholfr, ] # Filter dataframe to include only top 5 values
malchfr_per <- prop.table(mostalcoholfr) * 100 # percentage of meal variation
malchfr_df <- data.frame(Value = as.character(names(mostalcoholfr)),Count = as.numeric(mostalcoholfr),
                      Percentage = as.numeric(malchfr_per)) # percentage of meal var.

# Processed meat intake frequency :
prmeat <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$Pr_meat %in% c(1, 2, 3, 4, 5, 6), ] # numbers are type of variation
prmeat$Pr_meat <- as.character(prmeat$Pr_meat) # change of data class
prmeat_c <- table(prmeat$Pr_meat) # counts of variation
mostmeatpr <- head(sort(prmeat_c, decreasing = TRUE),6)  # Sort counts in descending order and take top 6 values
maxmeat <- prmeat[prmeat$liverdisease == 1 & prmeat$Pr_meat %in% mostmeatpr, ] # Filter dataframe to include only top 5 values
maxmeat_per <- prop.table(mostmeatpr) * 100 # percentage of meal variation
maxmeat_df <- data.frame(Value = as.character(names(mostmeatpr)),Count = as.numeric(mostmeatpr),
                         Percentage = as.numeric(maxmeat_per)) # percentage of meal var.

# No_edws : Never eat eggs, dairy, wheat, sugar
edws <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$No_edws %in% c(1, 2, 3, 4, 5 ), ] # numbers are type of variation
edws$No_edws <- as.character(edws$No_edws) # change of data class
predws<- table(edws$No_edws) # counts of variation
mostedws <- head(sort(predws, decreasing = TRUE),5)  # Sort counts in descending order and take top 5 values
maxedws <- edws[edws$liverdisease == 1 & edws$No_edws %in% mostedws, ] # Filter dataframe to include only top 5 values
maxedws_pr <- prop.table(mostedws) * 100 # percentage of meal variation
maxmeat_df <- data.frame(Value = as.character(names(mostedws)),Count = as.numeric(mostedws),
                         Percentage = as.numeric(maxedws_pr)) # percentage of meal var.

# Type of meals
meals <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$meals_typ %in% c(461, 462, 463, 464, 465), ]
mealcounts <- table(meals$meals)
mostmeals <- head(sort(counts, decreasing = TRUE),5)  # Sort counts in descending order and take top 5 values
maxmeals <- meals[meals$liverdisease == 1 & meals$meals %in% mostmeals, ] # Filter dataframe to include only top 5 values
meal_per <- prop.table(mostmeals) * 100
meals_df <- data.frame(Value = as.character(names(mostmeals)),Count = as.numeric(mostmeals),Percentage = as.numeric(meal_per))

# sleep 
mean_sleep <- mean(metformin_MASLD$eid, na.rm = TRUE)
mean_sleep <- mean(subset(ALL_Nutrition$Sleep, ALL_Nutrition$liverdisease == 1))
mean_sleep2 <- mean(subset(ALL_Nutrition$Sleep, ALL_Nutrition$liverdisease == 0))

# Type of fat/oil used in cooking

oiltyp <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$oil_typ %in% c(353, 354, 355, 359, 362, 365, 371, 373), ]
mealcounts <- table(oiltyp$oil_typ)
mostoiltyp <- head(sort(mealcounts, decreasing = TRUE),8)  # Sort counts in descending order and take top 5 values
maxoiltyp <- oiltyp[oiltyp$liverdisease == 1 & oiltyp$oil_typ %in% mostoiltyp, ] # Filter dataframe to include only top 5 values
oiltyp_per <- prop.table(mostoiltyp) * 100
oiltyp_df <- data.frame(Value = as.character(names(mostoiltyp)),Count = as.numeric(mostoiltyp),Percentage = as.numeric(oiltyp_per))

# Work_sat pschological cause :

Worksat <- ALL_Nutrition[ALL_Nutrition$liverdisease == 1 & ALL_Nutrition$Work_sat %in% c(1, 2, 3, 4, 5, 6, 7), ]
Workcounts <- table(Worksat$Work_sat)
mostWork <- head(sort(Workcounts, decreasing = TRUE),7)  # Sort counts in descending order and take top 5 values
maxWorksat <- Worksat[Worksat$liverdisease == 1 & Worksat$Work_sat %in% mostWork, ] # Filter dataframe to include only top 5 values
worksat_per <- prop.table(mostWork) * 100
worksat_df <- data.frame(Value = as.character(names(mostWork)),Count = as.numeric(mostWork),Percentage = as.numeric(worksat_per))

