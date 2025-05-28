library(forecast)
library(tseries)
library(dplyr)
library(lmtest)


# Importation
chemin <- "C:/Users/USER/OneDrive/Documents/Projet SC (LAGZOULI Riyad - OUAJA Med Nabil)/defaillances_paris.txt"
data <- read.table(chemin, header = TRUE)

# Conversion en série temporelle
ts_data <- ts(data$Nombre_de_faillites, 
              start = c(1995, 1), 
              frequency = 4)

# Statistiques descriptives de base
summary(data$Nombre_de_faillites)
cat("Écart-type :", sd(data$Nombre_de_faillites), "\n")
cat("Variance :",var(data$Nombre_de_faillites), "\n")

# Histogramme
hist(data$Nombre_de_faillites,
     main = "Histogramme du nombre de faillites",
     xlab = "Nombre de faillites",
     col = "skyblue",
     border = "white")

# Boîte à moustaches
boxplot(data$Nombre_de_faillites,
        main = "Boîte à moustaches",
        ylab = "Nombre de faillites",
        col = "orange")

# Division  (1995-2016 / 2017-2021)
train <- window(ts_data, end = c(2011, 4))
test <- window(ts_data, start = c(2012, 1))

plot(ts_data, 
     main = "Faillites trimestrielles à Paris (1995-2021)",
     ylab = "Nombre de faillites", 
     xlab = "Année",
     col = "darkblue",
     xaxt = "n")
axis(1, at = seq(1995, 2021, by = 1), labels = seq(1995, 2021, by = 1), las = 2)

# Test ADF
adf.test(train)
# Test KPSS
kpss.test(train)

train_diff <- diff(train, differences = 1)
adf.test(train_diff)  
kpss.test(train_diff) 

acf(train_diff, lag.max = 36, main = "ACF après Différenciation")
pacf(train_diff, lag.max = 36, main = "PACF après Différenciation")

grid <- expand.grid(
  p = 0:3,      # AR non-saisonnier (0 à 3)
  d = 1,        # Différenciation régulière fixée
  q = 0:1,      # MA non-saisonnier (0 ou 1)
  P = 0,        # SAR saisonnier fixé à 0
  D = 1,        # Différenciation saisonnière fixée
  Q = 0:1       # SMA saisonnier (0 ou 1)
)
head(grid)

evaluate_sarima <- function(train, grid) {
  results <- data.frame()
  
  for (i in 1:nrow(grid)) {
    order <- c(grid$p[i], grid$d[i], grid$q[i])
    seasonal <- c(grid$P[i], grid$D[i], grid$Q[i])
    
    # Estimation avec gestion des erreurs
    model <- tryCatch({
      Arima(train, 
            order = order,
            seasonal = list(order = seasonal, period = 4),
            method = "ML")  # Maximum Likelihood pour comparer AIC
    }, error = function(e) NULL)
    
    if (!is.null(model)) {
      # Tests sur les résidus
      residuals <- residuals(model)
      box_test <- Box.test(residuals, lag = 24, type = "Ljung-Box")
      
      # Stockage des résultats
      results <- rbind(results, data.frame(
        Model = sprintf("(%d,%d,%d)(%d,%d,%d)[4]", 
                        order[1], order[2], order[3],
                        seasonal[1], seasonal[2], seasonal[3]),
        AIC = model$aic,
        BIC = AIC(model, k = log(length(train)))
      ))
    }
  }
  
  return(results)
}
evaluate_sarima(train, grid)

# Estimation du modèle
model <- arima(train, 
               order = c(0, 1, 1),
               seasonal = list(order = c(0, 1, 1), period = 4))

coeftest(model)

#Calcul des résidus 
residus <- residuals(model)

# Test de Box-Pierce 
box_test <- Box.test(residus, type = "Box-Pierce")

# Affichage des résultats
cat("=== Test de Box-Pierce ===\n")
cat("Statistique X-squared:", box_test$statistic, "\n")
cat("p-value:", box_test$p.value, "\n")

# Préparation des prévisions 
forecast_values <- forecast(model, h = length(test))

plot(forecast_values, main = "Prévisions vs Réelles")
lines(test, col = "red")
legend("topleft", 
       col = c("black", "blue", "red"), lty = 1)

# Calcul des mesures d'erreur
mape=mean(abs(1-forecast_values$mean/test ))*100
mape
rmse=sqrt(mean((test - forecast_values$mean)^2))
rmse

# Test de normalité des erreurs
forecast_errors <- test - forecast_values$mean
shapiro.test(forecast_errors)