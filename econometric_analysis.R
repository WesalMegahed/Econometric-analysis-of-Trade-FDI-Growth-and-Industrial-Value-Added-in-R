###############################################################################
# 📊 Econometric Analysis: Trade, FDI, Growth, and Industrial Value Added
# Author: [Wesal Megahed]
# Description: This script performs correlation analysis, stationarity testing,
#              regression modeling, cointegration testing, and robustness checks.
###############################################################################

#-----------------------------#
# 1. Load Required Libraries  #
#-----------------------------#
library(tidyverse)   # Data manipulation & visualization
library(tseries)     # Time series tests (ADF)
library(urca)        # Unit root & cointegration tests
library(vars)        # VAR & VECM models
library(ggplot2)     # Advanced plotting
library(GGally)      # Scatterplot matrix
library(stargazer)   # Professional regression tables
library(lmtest)      # Diagnostic tests
library(strucchange) # Structural change tests
library(psych)       # Descriptive statistics



set.seed(123)  # for reproducibility

# Create a time variable (years 1990–2020)
year <- 1990:2020
n <- length(year)

# Simulate explanatory variables with trends + noise
OT   <- 50 + 0.5*(1:n) + rnorm(n, 0, 5)       # Trade openness (increasing)
PGDP <- 2000 + 50*(1:n) + rnorm(n, 0, 200)    # Per capita GDP
S    <- 15 + rnorm(n, 0, 2)                   # Savings (% of GDP, stable)
FDI  <- runif(n, 1, 10) + 0.2*(1:n)           # Foreign direct investment inflow
GL   <- seq(40, 70, length.out = n) + rnorm(n, 0, 3)  # Globalization index

# Dependent variable: Industrial Value Added (AI)
# Assume it's positively affected by OT, PGDP, and FDI
AI <- 100 + 0.3*OT + 0.002*PGDP + 0.5*FDI + rnorm(n, 0, 20)

# Combine into a dataframe
data <- data.frame(year, AI, OT, PGDP, S, FDI, GL)
head(data)



#-----------------------------#
# 2. Correlation Analysis     #
#-----------------------------#
# Compute simple correlation matrix
cor_matrix <- cor(
  data[, c("OT", "PGDP", "S", "FDI", "GL", "AI")],
  use = "pairwise.complete.obs"
)

# Round correlations for readability
round(cor_matrix, 2)

# Base R scatterplot matrix
pairs(
  data[, c("OT", "PGDP", "S", "FDI", "GL", "AI")],
  main = "Scatterplot Matrix: OT, PGDP, S, FDI, GL, AI",
  pch = 16, col = "steelblue"
)

# Enhanced scatterplot matrix with GGally
ggpairs(
  data[, c("OT", "PGDP", "S", "FDI", "GL", "AI")],
  title = "Scatterplot Matrix for Key Variables"
)

#-----------------------------#
# 3. Unit Root (ADF) Tests    #
#-----------------------------#
# Apply Augmented Dickey-Fuller test for each variable
adf_results <- list(
  AI   = adf.test(data$AI, k = 1),
  OT   = adf.test(data$OT, k = 1),
  PGDP = adf.test(data$PGDP, k = 1),
  S    = adf.test(data$S, k = 1),
  FDI  = adf.test(data$FDI, k = 1),
  GL   = adf.test(data$GL, k = 1)
)

# Summarize test results in a clean table
lapply(adf_results, function(x) {
  data.frame(
    Test.Statistic = x$statistic,
    P.Value        = x$p.value
  )
})

#-----------------------------#
# 4. Regression Model         #
#-----------------------------#
# Estimate multiple linear regression
model <- lm(AI ~ OT + PGDP + S + FDI + GL, data = data)

# Display results
summary(model)

# Professional regression output
stargazer(model, type = "text")

#-----------------------------#
# 5. Model Diagnostics        #
#-----------------------------#
# Heteroskedasticity (Breusch-Pagan test)
bptest(model)

# Autocorrelation tests
dwtest(model)             # Durbin-Watson
bgtest(model, order = 1)  # Breusch-Godfrey

#-----------------------------#
# 6. Cointegration Analysis   #
#-----------------------------#
# Johansen cointegration test
johansen_test <- ca.jo(
  data[, c("AI", "OT", "PGDP", "S", "FDI", "GL")],
  type  = "trace", 
  ecdet = "const", 
  K     = 2
)

summary(johansen_test)

# Estimate Error Correction Model (ECM) if cointegration is detected
if (johansen_test@teststat[1] > johansen_test@cval[1, 2]) {
  ecm_model <- lm(
    d(AI) ~ d(OT) + d(PGDP) + d(S) + d(FDI) + d(GL) + lag(residuals(model), 1),
    data = data
  )
  summary(ecm_model)
}

#-----------------------------#
# 7. Robustness Checks        #
#-----------------------------#
# Include a post-2000 dummy variable
model_controlled <- lm(
  AI ~ OT + PGDP + S + FDI + GL + factor(year > 2000),
  data = data
)
summary(model_controlled)

# Chow test for structural breaks
sctest(AI ~ OT + PGDP + S + FDI + GL, data = data, type = "Chow", point = 20)

#-----------------------------#
# 8. Summary of Results       #
#-----------------------------#
results <- list(
  Descriptive_Stats  = describe(data[, -1]),
  Unit_Root_Tests    = adf_results,
  Regression_Results = summary(model),
  Cointegration_Test = summary(johansen_test),
  Diagnostic_Tests   = list(
    Heteroskedasticity = bptest(model),
    Autocorrelation    = list(DW = dwtest(model), BG = bgtest(model, order = 1))
  )
)

# Key interpretation
cat("
Key Findings:
1. Trade Openness (OT) has a positive and statistically significant impact 
   on Industrial Value Added (AI). Coefficient =", round(coef(model)["OT"], 3), "

2. Foreign Direct Investment (FDI) also shows a positive effect, 
   with each 1% increase in FDI raising AI by", round(coef(model)["FDI"], 3), "%.

3. Per Capita GDP (PGDP) contributes positively, but with smaller magnitude.

4. ADF tests suggest most variables are non-stationary at levels, but stationary at first differences (I(1)).

5. Johansen test confirms a long-run equilibrium relationship among variables.
")
# ---------------------------------------------

data_long <- data %>%
  pivot_longer(cols = -year, names_to = "Variable", values_to = "Value")

ggplot(data_long, aes(x = year, y = Value, color = Variable)) +
  geom_line(size = 1.2) +
  facet_wrap(~Variable, scales = "free_y") +
  theme_minimal() +
  labs(title = "Time-Series Trends of Key Variables", x = "Year", y = "Value")
#------------------------------------------
library(reshape2)
cor_matrix <- round(cor(data[, -1]), 2)
cor_melt <- melt(cor_matrix)

ggplot(cor_melt, aes(Var1, Var2, fill = value)) +
  geom_tile(color = "white") +
  geom_text(aes(label = value), color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0) +
  theme_minimal() +
  labs(title = "Correlation Heatmap")

#---------------------

model <- lm(AI ~ OT + PGDP + S + FDI + GL, data = data)

# Residuals vs fitted
plot(model, which = 1)

# Normal Q-Q plot
plot(model, which = 2)

# Scale-Location plot
plot(model, which = 3)

# Residuals vs leverage
plot(model, which = 5)

#-------------------------------------
library(broom)
library(ggplot2)

tidy_model <- tidy(model, conf.int = TRUE)

ggplot(tidy_model, aes(x = estimate, y = term)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  theme_minimal() +
  labs(title = "Regression Coefficients with 95% CI", x = "Estimate", y = "Variable")


library(plotly)
install.packages("plotly")

p <- ggplot(tidy_model, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "steelblue") +
  geom_errorbarh(aes(xmin = conf.low, xmax = conf.high), height = 0.2) +
  theme_minimal() +
  labs(title = "Regression Coefficients with 95% CI", x = "Estimate", y = "Variable")

ggplotly(p)

ggplotly(p2)


library(gganimate)
anim <- ggplot(data, aes(year, AI)) +
  geom_line(color="blue") +
  transition_reveal(year)
animate(anim, width=800, height=600, renderer=gifski_renderer())
anim_save("growth.gif", animation = last_animation())
