library(sf)
library(dplyr)
library(spdep)
library(splm)
library(ggplot2)
library(plm)
library(readr)

world <- st_read("ArcGIS/Projects/APCOMP209_EDA/R/World_Countries.gpkg")
cardio_COMPLETE <- read_csv("ArcGIS/Projects/APCOMP209_EDA/R/cardio_COMPLETE.csv")
cardio_COMPLETE$ISO_CC <- cardio_COMPLETE$COUNTRY_REGION

countries_ischemic <- world %>% left_join(cardio_COMPLETE, by = 'ISO_CC', relationship = "many-to-many")
countries_ischemic <- countries_ischemic[which(!is.na(countries_ischemic$Year)),]

# get unique coutry names and group by ISO-3 digit
unique_countries <- countries_ischemic %>%
group_by(ISO_CC) %>%
summarize(.groups = "drop")

# set rownames to ISO_CC s.t. neighbor list retains names
rownames(unique_countries) <- unique_countries$ISO_CC

# keep geometry and model from standard df
df_panel <- st_drop_geometry(countries_ischemic)

#UPDATE 4/17 LOG and SCALE, skewed variables from ipynb
vars_to_transform <- c('IschemicHeartDisease_All','PopDensity', 'GDPperCapita', 'UnemploymentRate',
    'FruitsVeggiesAvailability', 'Hospitalsper100k', 'HospitalBedper100k',
    'AvgLengthOfStay')



# remove duplicates
vars <- c("IschemicHeartDisease_All", "PopDensity", "UnemploymentRate", "GDPperCapita", 
          "SocialBenefits", "HDI", "CalorieAvailability", "FatAvailability", 
          "FruitsVeggiesAvailability", "Hospitalsper100k", "HospitalBedper100k", 
          "Physiciansper100k", "AvgLengthOfStay", "HealthExpenditurePctGDP", "GINIcoefficient", 
          "CirculatoryDisease_65plus", "IschemicHeartDisease_0to64", "Pop65plus", 
          "OverweightPrevalence")

df_panel <- df_panel %>% 
  mutate(across(all_of(vars_to_transform), log))

# Then, scale (z-score standardize) all variables of interest.
# For variables that were log-transformed, scaling is applied on the log‐values.
df_panel <- df_panel %>% 
  mutate(across(all_of(vars), ~ as.vector(scale(.x))))


df_panel_unique <- df_panel %>%
  group_by(ISO_CC, Year) %>%
  summarize(across(all_of(vars), ~ mean(.x, na.rm = TRUE)), .groups = "drop")

# Specify the years of interest
years <- 1985:2020

# Subset the data to only include the specified years
df_subset <- df_panel_unique[df_panel_unique$Year %in% years, ]

# Determine the number of years in the subset
n_years <- length(years)

# Identify countries with data in every year of interest
country_counts <- table(df_subset$ISO_CC)
balanced_countries <- names(country_counts[country_counts == n_years])

# Subset the panel data to only include balanced countries
df_balanced <- df_subset[df_subset$ISO_CC %in% balanced_countries, ]
cat("Number of countries in balanced panel:", length(balanced_countries), "\n")

### SMOOTHING PREP

library(gt)
data(countrypops, package = "gt")
str(countrypops)  # Inspect the structure to check column names

# Create a mapping of ISO codes to total population
# (Adjust column names based on the actual structure, here we assume 'iso3' and 'pop_est')
iso_pop_mapping <- countrypops[, c("country_code_3", "population")]
print(iso_pop_mapping)

# Filter for countries with a total population below 50,000
small_pop <- iso_pop_mapping[iso_pop_mapping$population < 50000, ]

# Remove duplicates by grouping by the ISO code and taking the row with the maximum population
small_pop_unique <- small_pop %>%
  group_by(country_code_3) %>%
  slice_max(population, with_ties = FALSE) %>%
  ungroup()

cat("Unique countries with total population below 50,000 (most recent estimates):\n")
smoothing_threshold_countries <- small_pop_unique[small_pop_unique$country_code_3 %in% countries_ischemic$ISO_CC,]
print(small_pop_unique[small_pop_unique$country_code_3 %in% countries_ischemic$ISO_CC,])


# --- EMPIRICAL BAYES SMOOTHING --- (not spatially applied here)


# Merge the small-population mapping with our balanced panel to add the population info
df_balanced_sm <- df_balanced %>%
  left_join(smoothing_threshold_countries, by = c("ISO_CC" = "country_code_3"))
# (Countries not present in the small-pop list will have NA in the 'population' column.)

# Set a smoothing constant (λ); here, we use 50,000
lambda <- 50000

# Compute global mean rates from the balanced panel (for each variable to be smoothed)
global_IHD <- mean(df_balanced$IschemicHeartDisease_All, na.rm = TRUE)
global_HBed <- mean(df_balanced$HospitalBedper100k, na.rm = TRUE)

# For each record in df_balanced_sm, if the country is small (population not NA),
# compute a weight, then calculate a smoothed value as a weighted average of the observed value and the global mean.
df_balanced_sm <- df_balanced_sm %>%
  mutate(
    weight = ifelse(!is.na(population), population / (population + lambda), 1),
    IHD_EB = ifelse(!is.na(population),
                    weight * IschemicHeartDisease_All + (1 - weight) * global_IHD,
                    IschemicHeartDisease_All),
    HospitalBed_EB = ifelse(!is.na(population),
                            weight * HospitalBedper100k + (1 - weight) * global_HBed,
                            HospitalBedper100k)
  )

# Document decisions:
# - For countries with total population < 50,000 (as per the mapping), we replace the observed rate (or value)
#   with a smoothed value.
# - The weight is computed as w = population / (population + λ) with λ set to 50,000.
# - Thus, if a country has a very small population, w will be small and the smoothed value will move closer to the global mean.
# - For countries not flagged as small, the weight is 1 and the observed value remains unchanged.
cat("\nA summary of smoothed values (first few rows):\n")
print(head(df_balanced_sm[, c("ISO_CC", "Year", "IschemicHeartDisease_All", "IHD_EB",
                              "HospitalBedper100k", "HospitalBed_EB")]))



# --- SPATIAL WEIGHTS MATRIX ---

# 2. Create a spatial weights matrix based on the balanced set of countries
# Subset the country polygons to only those in the balanced panel
unique_countries_balanced <- unique_countries[unique_countries$ISO_CC %in% balanced_countries, ]

# Compute centroids and coordinates
centroids <- st_centroid(unique_countries_balanced)
coords <- st_coordinates(centroids)

# Determine the number of neighbors (at most 4, but must be less than the number of available countries)
k_value <- min(4, nrow(coords) - 1)
knn_obj <- knearneigh(coords, k = k_value)
knn_nb <- knn2nb(knn_obj)

# Assign names to the neighbors list using ISO country codes
names(knn_nb) <- unique_countries_balanced$ISO_CC

# Create a spatial weights list with row-standardized weights
listw <- nb2listw(knn_nb, style = "W")

cat("\nSummary of spatial weights matrix:\n")
print(summary(listw))

neighbors_count <- sapply(knn_nb, length)
avg_neighbors <- mean(neighbors_count)
cat("\nAverage number of neighbors:", avg_neighbors, "\n")


# 3. Convert the balanced data to a panel data object
panel_data <- pdata.frame(df_balanced, index = c("ISO_CC", "Year"))


# --- GLOBAL MORAN'S I ON DEPENDENT VAR ---

vars_of_interest <- c("PopDensity", "UnemploymentRate", "GDPperCapita", 
                      "SocialBenefits", "HDI", "CalorieAvailability", 
                      "FatAvailability", "FruitsVeggiesAvailability", 
                      "Hospitalsper100k", "HospitalBedper100k", "Physiciansper100k", 
                      "AvgLengthOfStay", "HealthExpenditurePctGDP", "GINIcoefficient", 
                      "CirculatoryDisease_65plus", "IschemicHeartDisease_0to64", 
                      "Pop65plus", "OverweightPrevalence", "IschemicHeartDisease_All")

# Aggregate: compute the mean for each variable by country.
df_country <- df_balanced %>%
  group_by(ISO_CC) %>%
  summarize(across(all_of(vars_of_interest), mean, na.rm = TRUE), .groups = "drop")

avg_pop_mapping <- countrypops %>%
  filter(year %in% years) %>%
  group_by(country_code_3) %>%
  summarize(TotalPopulation = mean(population, na.rm = TRUE), .groups = "drop")

# -------------------------------
# 3. Join the Average Population to the Aggregated Country Data
# -------------------------------
# Join on ISO_CC (from your data) to country_code_3 (from countrypops)
df_country <- left_join(df_country, avg_pop_mapping, by = c("ISO_CC" = "country_code_3"))

# For subsequent Moran's I, create a vector that includes all variables of interest plus TotalPopulation.
vars_moran <- c(vars_of_interest, "TotalPopulation")

# -------------------------------
# 4. Reorder the Aggregated Data to Match the Spatial Weights List Ordering
# -------------------------------
# Here we assume the names in the 'listw' neighbor list correspond to ISO_CC.
ordered_df <- df_country[match(names(listw$neighbours), df_country$ISO_CC), ]

# Check for missing matches.
if(any(is.na(ordered_df$ISO_CC))) {
  warning("Some countries in the spatial weights list are not present in df_country. Please verify your join.")
}

# -------------------------------
# 5. Run Global Moran’s I for Each Variable and Store the Results
# -------------------------------
result_list <- lapply(vars_moran, function(var) {
  vec <- ordered_df[[var]]
  # Run Global Moran's I test.
  test <- moran.test(vec, listw)
  # Extract the required statistics.
  stat   <- as.numeric(test$estimate["Moran I statistic"])
  expect <- as.numeric(test$estimate["Expectation"])
  var_val<- as.numeric(test$estimate["Variance"])
  z_val  <- as.numeric(test$statistic)  # Z-score
  p_val  <- test$p.value               # p-value
  return(c("Moran I statistic" = stat,
           "Expectation" = expect,
           "Variance" = var_val,
           "Z-score" = z_val,
           "P-value" = p_val))
})

# Convert the results list to a data frame.
result_df <- as.data.frame(do.call(cbind, result_list))
colnames(result_df) <- vars_moran
rownames(result_df) <- c("Moran I statistic", "Expectation", "Variance", "Z-score", "P-value")






# 4. Fit a spatial panel model using spml (fixed effects with spatial lag)
#include factor(Year) to control for time effects.
# consult prior literature for Year as ind. variable => treated as cat. variable
# before spatial models, perform OLS models, check skewness, normalization dist. of variables, 
# decide if necessary to log transform variables and scale (Z-score or min-max)
# for spatial regression models, begin with spatial lag, spatial error, GWR, 
#   GTWR, must base on OLS models, check Moran's I of dep. var., check Moran's I for residuals and error
#       see if Moran's I residuals significantly clustered for each country, explore GWR at this point
#          will give attribute table for coef for each country (check Moran's I for coef of each country)
#   3 types of spatial regression: different for which variables have significant clustering
#   Prefer Moran's I highest and z-score but could be different.

# Some of the predictors need to be transformed
# e.g. log transform + z-score normalization for GDP
# check multicollinearity and VIF

### ---- STEPS 2/28 -----
# 1. Select important variables (may include 1-2 potential variables you have not decided if we want to include or not), show the distribution
# 2. Based on the distribution, try log transformation for those highly skewed variables or try log transformation for all variables (including dep.) if necessary
### => update Dr. Chen on intended distribution 
# 3. Z-score normalize all variables, including dependent variables (b/c better OLS coef interpretation )
# 4. Check Variance Inflation Factor (VIF) for multi-collinearity among variables, exclude variables with VIR>5
# 5. Run OLS model, play around with the variables you are not certain about to see how they affect your OLS results.


# 6. Now we will use those variables that are significantly correlated for further spatial analysis


### ---- SPEAK AFTER step 6 to decide spat model

# 7.  Check global moran's I for depedent variables, error/residual of OLS models and compare them
# 8. Run GWR model and check global moran's Is for coefficients for the variables
# 9. Compare all global moran's I and decide which spatial models we want to use
# 10. If we go with GWR, we will further check spatial nonstationarity of variables by comparing OLS and GWR variables
important_predictors <- c('GDPperCapita', 'FruitsVeggiesAvailability', 'HospitalBedper100k', 'Pop65plus', 'UnemploymentRate', 'Physiciansper100k')

# FIRST, fit OLS and run spatial autocorrelation for residuals to determine type of spatial model


formula_str <- paste("IschemicHeartDisease_All ~", paste(important_predictors, collapse = " + "))
formula_model <- as.formula(formula_str)

# ---- UPDATE 5/7 - 
# residual diagnostics to determine spatial lag/error
#  covariate coefficients come from within-country variation over time
# lm treats (country, year) as independent observations

formula_str <- paste("IschemicHeartDisease_All ~", paste(important_predictors, collapse = " + "))
ols_cs <- lm(formula_str, data = panel_data)

# 2. Extract and aggregate residuals by country
resid_df <- data.frame(
  ISO_CC   = attr(panel_data, "index")[, "ISO_CC"],
  residual = residuals(ols_cs),
  stringsAsFactors = FALSE
)

mean_resid <- resid_df %>%
  group_by(ISO_CC) %>%
  summarize(resid = mean(residual, na.rm = TRUE))

# 3. Re‐order so that names(mean_resid_vec) == names(listw$neighbours)
mean_resid_vec <- setNames(mean_resid$resid, mean_resid$ISO_CC)
mean_resid_vec <- mean_resid_vec[ names(listw$neighbours) ]

# 4. Build an intercept‐only lm whose rownames match the weight matrix
mean_resid_df <- data.frame(resid = mean_resid_vec)
row.names(mean_resid_df) <- names(mean_resid_vec)
lm_agg <- lm(resid ~ 1, data = mean_resid_df)

# 5. Run LM and robust‐LM tests to see if you need a lag or an error term
lm_tests <- lm.RStests(
  model       = lm_agg,
  listw       = listw,
  test        = c("LMerr","LMlag"),
  zero.policy = TRUE
)
print(lm_tests)

# in practice one almost always fits either a pure spatial‐lag or a pure spatial‐error mode
# unless there is a very clear substantive reason to do otherwise (LeSage & Pace 2009, Ch. 7)
# “if both basic Lagrange multiplier tests reject, compare the robust versions and include only the term with the smaller p‐value” (Florax et al. 2003, p. 568)
# via Monte Carlo that blindly estimating both a ρ and a λ without theoretical backing leads to inflated size and poor power
# within removes cross-sectional noise
# lose N-1 df compared to pooled model b/c intercept for each country

# ---- UPDATE 5/7 - 
# first, rescale Year to start at zero to interpret the intercept)
panel_data <- panel_data %>%
  # create a numeric year column:
  mutate(YearNum = as.integer(as.character(Year))) %>%
  # then compute “time since first year”:
  group_by(ISO_CC) %>%
  mutate(time = YearNum - min(YearNum)) %>%
  ungroup()

# include continuous time var
formula_str_time <- paste0(
  "IschemicHeartDisease_All ~ time + ",
  paste(important_predictors, collapse = " + ")
)

model_lag <- spml(
  formula = as.formula(formula_str_time),
  listw = listw,
  data = panel_data,
  model = "within",   # pooled effects specification (based on ARIMA anal. similar dec. over time) => reducation r2
  lag = TRUE   ,       # include spatial lag of the dependent variable
  spatial.error = "none"        # exclude spatial error term
)

# Display the summary of the model
summary(model_lag)

# --- Extract Residuals ---
model_resid <- residuals(model_lag)
fitted_values <- fitted(model_lag)

# Create a data frame pairing each observation (country-year) with its residual
resid_df <- data.frame(
  ISO_CC   = as.character(index(panel_data, "individual")),
  Year     = index(panel_data, "time"),
  residual = model_resid
)
agg_resid <- resid_df %>%
  group_by(ISO_CC) %>%
  summarize(mean_resid = mean(residual, na.rm = TRUE))

unique_countries_balanced <- unique_countries_balanced %>%
  # if it isn’t already character, convert this too:
  mutate(ISO_CC = as.character(ISO_CC)) %>%
  left_join(agg_resid, by = "ISO_CC")


# Match the order of residuals with the order in the spatial weights list
# keep only those countries that actually appear in your listw
common_iso <- intersect(names(listw$neighbours), agg_resid$ISO_CC)

# then pull the means in the exact same order as the weights object
ordered_resid <- agg_resid %>%
  filter(ISO_CC %in% common_iso) %>%
  arrange(match(ISO_CC, names(listw$neighbours))) %>%
  pull(mean_resid)


# Run Moran's I test on aggregated residuals
#moran_test <- moran.test(ordered_resid, listw, na.action = na.exclude)
#print(moran_test)


# --- Moran Scatter Plot ---
#moran.plot(model_resid, listw, main = "Moran Scatterplot of Residuals")

# --- Residual vs. Fitted Values Plot ---
plot(fitted_values, model_resid, 
     main = "Residuals vs Fitted Values", 
     xlab = "Fitted Values", 
     ylab = "Residuals",
     pch = 20, col = "blue")
abline(h = 0, col = "red", lwd = 2)

# --- QQ Plot for Residuals ---
qqnorm(model_resid, main = "Normal Q-Q Plot of Residuals", pch = 20, col = "blue")
qqline(model_resid, col = "red", lwd = 2)

# --- Spatial Map of Residuals ---
# Merge the residuals with the spatial polygons of the balanced countries.
# Assumes 'unique_countries_balanced' has a column 'ISO_CC' that matches panel_data.
unique_countries_balanced$residuals <- model_resid[match(unique_countries_balanced$ISO_CC, names(model_resid))]

#library(ggplot2)
#ggplot(unique_countries_balanced) +
 # geom_sf(aes(fill = residuals)) +
  #scale_fill_viridis_c(option = "viridis") +
  #labs(title = "Spatial Distribution of Model Residuals",
   #    fill = "Residuals") +
  #theme_minimal()



### ADDITIONAL DIAGNOSTICS
# Extract the observed dependent variable
y <- panel_data$IschemicHeartDisease_All

# Compute Sum of Squared Residuals (SSR)
SSR <- sum(residuals(model_lag)^2)

# Compute Total Sum of Squares (TSS) relative to the overall mean of y
TSS <- sum((y - mean(y))^2)

# Compute the "pseudo" R²
R2 <- 1 - SSR / TSS

# Determine the number of observations and estimated parameters.
# Note: if your model includes an intercept, subtract one from the parameter count.
n <- length(y)
p <- length(coef(model_lag)) - 1  # subtract 1 if an intercept is included

# Compute the adjusted R²
adjR2 <- 1 - ((n - 1) / (n - p - 1)) * (1 - R2)

# Print the results
cat("R-squared:", R2, "\n")
cat("Adjusted R-squared:", adjR2, "\n")

res  <- residuals(model_lag)
fit  <- fitted(model_lag)
n    <- length(res)

mse_spatial   <- mean(res^2)
mae_spatial   <- mean(abs(res))
rmse_spatial  <- sqrt(mse_spatial)
cat("mse", mse_spatial, "rmse:", rmse_spatial, "mae", mae_spatial)


# -----UPDATE 5/7

# moving forward with spatial lag for theoretical interpretability, high df with cont. time term, and lack of spatial autocorrelation among residuals.
# drop1 anal based on likelihood ratio test

library(purrr)

# extract the terms to test
full_form    <- as.formula(formula_str_time)
terms_to_drop <- attr(terms(full_form), "term.labels")

# store full model’s logLik once
full <- model_lag
ll_full <- full$logLik

drop_tests <- map_dfr(terms_to_drop, function(trm) {
  # build a new RHS omitting trm
  rhs <- setdiff(terms_to_drop, trm)
  new_form <- as.formula(
    paste("IschemicHeartDisease_All ~", paste(rhs, collapse = " + "))
  )
  # refit with spml()
  reduced <- spml(
    formula       = new_form,
    data          = panel_data,
    listw         = listw,
    model         = "within",
    lag           = TRUE,
    spatial.error = "none"
  )
  ll_red  <- reduced$logLik
  
  # compute LRT & AIC by hand
  LRstat  <- 2*(ll_full - ll_red)
  df_diff <- length(coef(full)) - length(coef(reduced))
  p_val   <- pchisq(LRstat, df_diff, lower.tail = FALSE)
  AIC_red <- -2*ll_red + 2*length(coef(reduced))
  
  tibble(
    dropped    = trm,
    df_diff    = df_diff,
    LRstat     = LRstat,
    p.value    = p_val,
    AIC_reduced= AIC_red
  )
})

print(drop_tests)

# note: within for fixed effects by country versus pooling reduced AIC by half even with more est. params

rhs_terms    <- attr(terms(full_form), "term.labels")
reduced_rhs  <- setdiff(rhs_terms, "UnemploymentRate")
reduced_form <- as.formula(
  paste("IschemicHeartDisease_All ~", paste(reduced_rhs, collapse = " + "))
)

model_nounemp <- spml(
  formula       = reduced_form,
  data          = panel_data,
  listw         = listw,
  model         = "within",
  lag           = TRUE,
  spatial.error = "none"
)

summary(model_nounemp)

resid_nu <- residuals(model_nounemp)
fit_nu   <- fitted(model_nounemp)

y         <- panel_data$IschemicHeartDisease_All
SSR_nu    <- sum(resid_nu^2, na.rm=TRUE)
TSS       <- sum((y - mean(y))^2)
R2_nu     <- 1 - SSR_nu/TSS
n_obs     <- length(resid_nu)
p_nu      <- length(coef(model_nounemp)) - 1    # subtract intercept if any
adjR2_nu  <- 1 - ((n_obs - 1)/(n_obs - p_nu - 1))*(1 - R2_nu)

cat("Pseudo-R² (no Unemp):",    round(R2_nu,3),  "\n")
cat("Adj. R² (no Unemp):",      round(adjR2_nu,3),"\n")

# 6. Compute basic error metrics
mse_nu  <- mean(resid_nu^2,   na.rm=TRUE)
mae_nu  <- mean(abs(resid_nu),na.rm=TRUE)
rmse_nu <- sqrt(mse_nu)

cat("MSE:",   round(mse_nu,3),
    "  RMSE:", round(rmse_nu,3),
    "  MAE:", round(mae_nu,3), "\n")

# Moran’s I on country‐mean residuals
resid_df_nu <- data.frame(
  ISO_CC   = as.character(index(panel_data, "individual")),
  residual = resid_nu
) %>%
  group_by(ISO_CC) %>%
  summarize(mean_resid = mean(residual, na.rm=TRUE))

common_iso <- intersect(names(listw$neighbours), resid_df_nu$ISO_CC)
ord_resid  <- resid_df_nu %>%
  filter(ISO_CC %in% common_iso) %>%
  arrange(match(ISO_CC, names(listw$neighbours))) %>%
  pull(mean_resid)

#moran_nu <- moran.test(ord_resid, listw, na.action = na.exclude)
#print(moran_nu)

# Re-plot Residual vs. Fitted and Q–Q
plot(fit_nu, resid_nu,
     main = "Residuals vs Fitted (no Unemp)",
     xlab = "Fitted", ylab = "Residuals",
     pch = 20, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)

qqnorm(resid_nu, main = "Q–Q Plot of Residuals (no Unemp)", pch=20)
qqline(resid_nu, col="red", lwd=2)

# drop1 test again 

nu_form    <- formula(model_nounemp)
terms_nu   <- attr(terms(nu_form), "term.labels")
ll_nu_full <- model_nounemp$logLik

drop_tests_nu <- map_dfr(terms_nu, function(trm) {
  rhs       <- setdiff(terms_nu, trm)
  new_form  <- as.formula(paste("IschemicHeartDisease_All ~", paste(rhs, collapse = " + ")))
  reduced_nu <- spml(
    formula       = new_form,
    data          = panel_data,
    listw         = listw,
    model         = "within",
    lag           = TRUE,
    spatial.error = "none"
  )
  ll_nu_red <- reduced_nu$logLik
  LRstat    <- 2 * (ll_nu_full - ll_nu_red)
  df_diff   <- length(coef(model_nounemp)) - length(coef(reduced_nu))
  p_val     <- pchisq(LRstat, df_diff, lower.tail = FALSE)
  AIC_red   <- -2 * ll_nu_red + 2 * length(coef(reduced_nu))
  tibble(
    dropped     = trm,
    df_diff     = df_diff,
    LRstat      = LRstat,
    p.value     = p_val,
    AIC_reduced = AIC_red
  )
})

print(drop_tests_nu)

# DROP HospitalBedper100k

rhs_terms2    <- attr(terms(formula(model_nounemp)), "term.labels")
reduced_rhs2  <- setdiff(rhs_terms2, "HospitalBedper100k")
reduced_form2 <- as.formula(
  paste("IschemicHeartDisease_All ~", paste(reduced_rhs2, collapse = " + "))
)

model_nobh <- spml(
  formula       = reduced_form2,
  data          = panel_data,
  listw         = listw,
  model         = "within",
  lag           = TRUE,
  spatial.error = "none"
)

summary(model_nobh)

resid_nobh <- residuals(model_nobh)
fit_nobh   <- fitted(model_nobh)

# Re-plot Residual vs. Fitted and Q–Q
plot(fit_nobh, resid_nobh,
     main = "Residuals vs Fitted",
     xlab = "Fitted", ylab = "Residuals",
     pch = 20, col = "darkgreen")
abline(h = 0, col = "red", lwd = 2)

qqnorm(resid_nu, main = "Q–Q Plot of Residuals", pch=20)
qqline(resid_nu, col="red", lwd=2)

y           <- panel_data$IschemicHeartDisease_All
SSR_nobh    <- sum(resid_nobh^2,   na.rm = TRUE)
TSS         <- sum((y - mean(y))^2)
R2_nobh     <- 1 - SSR_nobh / TSS
n_obs2      <- length(resid_nobh)
p_nobh      <- length(coef(model_nobh)) - 1
adjR2_nobh  <- 1 - ((n_obs2 - 1) / (n_obs2 - p_nobh - 1)) * (1 - R2_nobh)

cat("Pseudo-R² (no Unemp, no HB):", round(R2_nobh, 3), "\n")
cat("Adj. R² (no Unemp, no HB):", round(adjR2_nobh, 3), "\n")

mse_nobh  <- mean(resid_nobh^2,    na.rm = TRUE)
mae_nobh  <- mean(abs(resid_nobh), na.rm = TRUE)
rmse_nobh <- sqrt(mse_nobh)

cat("MSE:", round(mse_nobh, 3),
    " RMSE:", round(rmse_nobh, 3),
    " MAE:", round(mae_nobh, 3), "\n")

resid_df_nobh <- data.frame(
  ISO_CC   = as.character(index(panel_data, "individual")),
  residual = resid_nobh
) %>%
  group_by(ISO_CC) %>%
  summarize(mean_resid = mean(residual, na.rm = TRUE))

common_iso2 <- intersect(names(listw$neighbours), resid_df_nobh$ISO_CC)
ord_resid2  <- resid_df_nobh %>%
  filter(ISO_CC %in% common_iso2) %>%
  arrange(match(ISO_CC, names(listw$neighbours))) %>%
  pull(mean_resid)

trs2           <- attr(terms(formula(model_nobh)), "term.labels")
ll_nobh_full   <- model_nobh$logLik

drop_tests_nobh <- map_dfr(trs2, function(trm) {
  rhs2      <- setdiff(trs2, trm)
  new_form2 <- as.formula(
    paste("IschemicHeartDisease_All ~", paste(rhs2, collapse = " + "))
  )
  reduced2  <- spml(
    formula       = new_form2,
    data          = panel_data,
    listw         = listw,
    model         = "within",
    lag           = TRUE,
    spatial.error = "none"
  )
  ll2       <- reduced2$logLik
  LR2       <- 2 * (ll_nobh_full - ll2)
  df2       <- length(coef(model_nobh)) - length(coef(reduced2))
  p2        <- pchisq(LR2, df2, lower.tail = FALSE)
  AIC2      <- -2 * ll2 + 2 * length(coef(reduced2))
  tibble(
    dropped     = trm,
    df_diff     = df2,
    LRstat      = LR2,
    p.value     = p2,
    AIC_reduced = AIC2
  )
})

print(drop_tests_nobh)

#### UPDATE 5/8 ---- export for visualization in ArcGIS, Moran's I on residuals

# 1. raw residuals (in the same order as panel_data)
resid_nobh <- residuals(model_nobh)

# 2. observed Y in the same order
y_obs <- panel_data$IschemicHeartDisease_All

# 3. fitted = observed − residual
fitted_nobh <- y_obs - resid_nobh

# 4. pull out the ISO_CC and Year vectors
iso <- as.character(panel_data$ISO_CC)
yr  <- as.integer(as.character(panel_data$Year))

# 5. assemble into a data.frame
df_nobh_obs <- data.frame(
  ISO_CC   = iso,
  Year     = yr,
  residual = resid_nobh,
  fitted   = fitted_nobh
)

# inspect
head(df_nobh_obs)

# average each country's residuals

average_country_residuals <- df_nobh_obs %>%
  group_by(ISO_CC) %>%
  summarize(
    avg_residual = mean(residual, na.rm = TRUE),
    avg_fitted   = mean(fitted,   na.rm = TRUE),
    std_residual = abs(mean(residual)),
    .groups = "drop"
  )

# 2. Export to CSV in your working directory
write_csv(average_country_residuals, "average_country_residuals.csv")

## DECADE AVERAGES


# 1. add a decade column
df_nobh_obs <- df_nobh_obs %>%
  mutate(decade = case_when(
    Year >= 1980 & Year < 1990 ~ "1980s",
    Year >= 1990 & Year < 2000 ~ "1990s",
    Year >= 2000 & Year < 2010 ~ "2000s",
    Year >= 2010 & Year < 2020 ~ "2010s",
    TRUE                       ~ NA_character_
  ))

# 2. helper function to average for a given decade
avg_by_decade <- function(dd) {
  df_nobh_obs %>%
    filter(decade == dd) %>%
    group_by(ISO_CC) %>%
    summarize(
      avg_residual  = mean(residual,   na.rm = TRUE),
      avg_fitted    = mean(fitted,     na.rm = TRUE),
      std_residual  = abs(mean(residual, na.rm = TRUE)),  # per your original definition
      .groups       = "drop"
    )
}

# 3. create the four data‐frames
average_1980s <- avg_by_decade("1980s")
average_1990s <- avg_by_decade("1990s")
average_2000s <- avg_by_decade("2000s")
average_2010s <- avg_by_decade("2010s")

# inspect
head(average_1980s)
head(average_1990s)
head(average_2000s)
head(average_2010s)

write_csv(average_1980s, "resid_aver_1980s.csv")
write_csv(average_1990s, "resid_aver_1990s.csv")
write_csv(average_2000s, "resid_aver_2000s.csv")
write_csv(average_2010s, "resid_aver_2010s.csv")

#final Moran's I resid

resid_vec <- setNames(
  average_country_residuals$avg_residual,
  average_country_residuals$ISO_CC
)

# 2. reorder to match the ordering in your listw neighbour list
resid_ordered <- resid_vec[ names(listw$neighbours) ]

# 3. run Global Moran’s I
moran_res <- moran.test(
  x       = resid_ordered,
  listw   = listw,
  zero.policy = TRUE
)

# 4. inspect the results
print(moran_res)

# 1. Get the OLS residual for each row in your panel:
resid_all <- residuals(ols_cs)

# 2. Pull out the ISO code for each row:
iso_all   <- as.character(panel_data[["ISO_CC"]])

# 3. Build a data.frame and average by ISO_CC:
mean_resid_df <- tibble(
  ISO_CC = iso_all,
  resid  = resid_all
) %>%
  group_by(ISO_CC) %>%
  summarize(avg_resid = mean(resid, na.rm = TRUE), .groups="drop")

# 4. Reorder to match the order in your spatial weights:
#    names(listw$neighbours) is the exact country order in `listw`
resid_ols <- mean_resid_df$avg_resid[
  match(names(listw$neighbours), mean_resid_df$ISO_CC)
]

# 5. Run Global Moran’s I:
moran_ols <- moran.test(
  x           = resid_ols,
  listw       = listw,
  zero.policy = TRUE
)
print(moran_ols)

# (Optional) Local Moran’s I for each country:
local_ols <- localmoran(
  x           = resid_ols,
  listw       = listw,
  zero.policy = TRUE
)
head(local_ols)


#### DATA VIZ BARPLOT IHD VS POP
library(forcats)
library(dplyr)
library(ggplot2)

# 1. Prep: grab only the two metrics and turn ISO into a factor ordered by population
plot_df <- df_country %>%
  select(ISO_CC, TotalPopulation, RawIHD = IschemicHeartDisease_All) %>%
  mutate(ISO_CC = fct_reorder(ISO_CC, TotalPopulation))

# 2. Compute a scale factor so that the IHD‐line will sit nicely on the same plot
scale_factor <- max(plot_df$TotalPopulation, na.rm = TRUE) / 
  max(plot_df$RawIHD,       na.rm = TRUE)

# 3. Make the combo plot
ggplot(plot_df, aes(x = ISO_CC)) +
  # bars for population
  geom_col(aes(y = TotalPopulation), fill = "steelblue") +
  
  # line + points for IHD (scaled up)
  geom_line(aes(y = RawIHD * scale_factor, group = 1), 
            color = "firebrick", size = 1) +
  geom_point(aes(y = RawIHD * scale_factor), 
             color = "firebrick", size = 2) +
  
  # dual y‐axis: left = population, right = back‐transformed IHD
  scale_y_continuous(
    name      = "Total Population",
    sec.axis  = sec_axis(
      ~ . / scale_factor, 
      name = "Mean IHD (z-score)"
    )
  ) +
  
  # labels & theme
  labs(
    x     = "Country (ISO-3)",
    title = "Countries Ordered by Population\nBars = Population | Line = Mean IHD"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )



# 1) Pull out the stats and set up colors (unchanged)
vars_to_plot <- c(important_predictors, "IschemicHeartDisease_All")
moran_I      <- as.numeric( result_df["Moran I statistic", vars_to_plot] )
names(moran_I) <- vars_to_plot
null_I       <- as.numeric( result_df["Expectation", vars_to_plot][1] )
bar_colors   <- ifelse(names(moran_I) == "IschemicHeartDisease_All",
                       "red", "steelblue")

# 2) Make room on the bottom for angled labels
old_mar <- par("mar")
par(mar = c(10, 4, 4, 2) + 0.1)  

# 3) Draw the bars, suppressing the default x axis
bp <- barplot(
  height = moran_I,
  col    = bar_colors,
  ylim   = c(-0.1, max(moran_I, null_I) + 0.1),  # y goes from -1 up
  xaxt   = "n",                                # turn off x axis
  main   = "Global Moran's I Clustering\nacross Selected Variables",
  ylab   = "Moran's I statistic"
)

# 4) Add a custom y-axis (optional: choose your own tick spacing)
axis(2, at = seq(-1, max(moran_I, null_I) + 0.1, by = 0.2))

# 5) Add the null‐expectation line
abline(h = null_I, col = "darkgrey", lty = 2, lwd = 2)

# 6) Draw angled x‐labels at 45°
text(
  x     = bp, 
  y     = par("usr")[3] - 0.03 * diff(par("usr")[3:4]),  # slightly below the bars
  labels= names(moran_I),
  srt   = 45, 
  adj   = 1,
  xpd   = TRUE
)

# 7) Legend (unchanged)
legend("topright",
       legend = c("Predictors",
                  "IschemicHeartDisease_All",
                  sprintf("Null E[I] = %.3f", null_I)),
       fill   = c("steelblue", "red", NA),
       border = NA,
       bty    = "n",
       lty    = c(NA,         NA,     2),
       col    = c(NA,         NA,     "darkgrey")
)

# 8) Restore original margins
par(mar = old_mar)


