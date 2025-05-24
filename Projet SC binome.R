library(forecast)
library(tseries)

# Importation
chemin <- "defaillances_paris.txt"
data <- read.table(chemin, header = TRUE)

# Conversion en série temporelle
ts_data <- ts(data$Nombre_de_faillites, 
              start = c(1995, 1), 
              frequency = 4)

# Division  (1995-2016 / 2017-2021)
train <- window(ts_data, end = c(2011, 4))
test <- window(ts_data, start = c(2012, 1))

# Paramétrage de l'axe des x pour afficher chaque année
plot(ts_data, 
     main = "Faillites trimestrielles à Paris (1995-2021)",
     ylab = "Nombre de faillites", 
     xlab = "Année",
     col = "darkblue",
     xaxt = "n")  # Désactive l'axe x automatique

# Ajout manuel de l'axe x avec toutes les années
axis(1, at = seq(1995, 2021, by = 1), labels = seq(1995, 2021, by = 1), las = 2)


# Test ADF
adf.test(train)

# Test KPSS
kpss.test(train)

train_diff <- diff(train, differences = 1)

adf.test(train_diff)  # Doit donner p-value < 0.05
kpss.test(train_diff) # Doit donner p-value > 0.05


acf(train_diff, lag.max = 36, main = "ACF après Différenciation")
pacf(train_diff, lag.max = 36, main = "PACF après Différenciation")

auto.arima(train, d = 1, seasonal = TRUE)


# 1. Estimation du modèle
model <- arima(train, 
               order = c(0, 1, 1),
               seasonal = list(order = c(2, 1, 1), period = 4))

# 2. Fonction pour afficher les résultats avec tests de Student
summary_sarima <- function(model) {
  # Calcul des statistiques
  estimates <- model$coef
  std_errors <- sqrt(diag(model$var.coef))
  t_values <- estimates/std_errors
  p_values <- 2*(1-pnorm(abs(t_values)))
  
  # Création du data.frame (sans virgule superflue)
  coef_table <- data.frame(
    Paramètre = names(estimates),
    Estimation = round(estimates, 4),
    Ecart_Type = round(std_errors, 4),
    t_value = round(t_values, 2),
    p_value = format.pval(p_values, eps = 0.0001),
    Signif = ifelse(p_values < 0.001, "***", 
                    ifelse(p_values < 0.01, "**",
                           ifelse(p_values < 0.05, "*", "")))
  )
  
  # Affichage
  cat("=== Modèle SARIMA(0,1,1)(2,1,1)[4] ===\n\n")
  print(coef_table, row.names = FALSE)
  cat("\nSigma^2 =", model$sigma2, 
      "\nLogLikelihood =", model$loglik,
      "\nAIC =", model$aic)
}

# 3. Application
summary_sarima(model)

model_simple <- arima(train, 
                      order = c(0,1,1),
                      seasonal = list(order = c(0,1,1), period = 4))



# 1. Calcul des résidus du modèle
residus <- residuals(model_simple)

# 2. Test de Box-Pierce avec différents lags (24 par défaut pour 6 ans)
box_test <- Box.test(residus, lag = 24, type = "Box-Pierce")

# 3. Affichage des résultats
cat("=== Test de Box-Pierce ===\n")
cat("Statistique X-squared:", box_test$statistic, "\n")
cat("p-value:", box_test$p.value, "\n")
cat("Lags utilisés:", box_test$parameter, "\n\n")

# 4. Interprétation
if(box_test$p.value > 0.05) {
  cat("Conclusion: Les résidus sont non corrélés (bruit blanc) au seuil de 5%")
} else {
  cat("Conclusion: Les résidus présentent une autocorrélation (modèle à revoir)")
}




# 1. Préparation des prévisions avec le modèle simplifié (recommandé)
forecast_values <- forecast(model_simple, h = length(test))

plot(forecast_values, main = "Prévisions vs Réelles")
lines(test, col = "red")
legend("topleft", 
       col = c("black", "blue", "red"), lty = 1)



# Calcul des mesures d'erreur
accuracy(forecast_values, test)

# Visualisation améliorée
plot(forecast_values, 
     main = "Prévisions vs Observations Réelles (2017-2021)",
     xlab = "Année", 
     ylab = "Nombre de faillites",
     col = "blue")
lines(test, col = "red")
legend("topleft", 
       legend = c("Prévisions", "IC 95%", "Observations réelles"),
       col = c("blue", "lightblue", "red"),
       lty = c(1, 1, 1), 
       cex = 0.8)

# Analyse des résidus de prévision
forecast_errors <- test - forecast_values$mean
plot(forecast_errors, main = "Erreurs de prévision")
abline(h = 0, col = "red")

# Test de normalité des erreurs
shapiro.test(forecast_errors)