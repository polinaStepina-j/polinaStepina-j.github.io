library(sf)#для работы с картами
library(tidyverse)#для работы с таблицами данных
library(readxl)#для xls файлов
library(spdep)
library(robustHD) #для стандартизации переменных - standardized()
library(tmap)
library(rgeoda)
#дополнительно используется библиотека dplyr
options(scipen = 999)
#ЗАГРУЗКА ДАННЫХ ПО СУБЪЕКТАМ РФ
map <- st_read('RFregions_2019.gpkg')#загружаемая карта представлена в границах 2019г.
class(map)
##таблица с данными
data <- read_xlsx("DATA.xlsx")
class(data)
##объединяем данные по oktmo (обратить внимание на Архангельскую и Тюменскую области)
MAP <- right_join(map, data, by  ="oktmo")#если ставить первой таблицу данных то тип sf теряется
class(MAP)


#######1) Ферзь 1го порядка (бинарная, нормированная по строкам)######
nb_queen1 = poly2nb(MAP) 
WbinW_qSpdep0 <- nb2listw(nb_queen1, style = "W", zero.policy=TRUE)#взвешивание построчно
nb_queen1[[which(MAP$NAME_RUS == "Калининградская область", arr.ind = FALSE, useNames = TRUE)]] <- c(which(MAP$NAME_RUS == "Смоленская область", arr.ind = FALSE, useNames = TRUE),
                                                                                                     which(MAP$NAME_RUS == "Город Санкт-Петербург город федерального значения", arr.ind = FALSE, useNames = TRUE))
nb_queen1[[which(MAP$NAME_RUS == "Сахалинская область", arr.ind = FALSE, useNames = TRUE)]] <- which(MAP$NAME_RUS == "Хабаровский край", arr.ind = FALSE, useNames = TRUE)
nb_queen1[[which(MAP$NAME_RUS == "Краснодарский край", arr.ind = FALSE, useNames = TRUE)]] <- sort(c(nb_queen1[[which(MAP$NAME_RUS == "Краснодарский край", arr.ind = FALSE, useNames = TRUE)]],
                                                                                                     which(MAP$NAME_RUS == "Республика Крым", arr.ind = FALSE, useNames = TRUE)))
nb_queen1 <- make.sym.nb(nb_queen1)
coords = st_coordinates(st_centroid(MAP)) 
plot(st_geometry(MAP), border = "gray50")
plot(nb_queen1, coords, pch = 19, cex = 0.5, add = TRUE)
title(main = "Соседи первого порядка (правило ферзя):
      черные линии по центроидам")
WbinW_qSpdep1 <- nb2listw(nb_queen1, style = "W", zero.policy=TRUE)
######2) Соседство по методу k-ближайших соседей ######
nb_knn_5 = knn2nb(knearneigh(coords, k = 5))
plot(st_geometry(MAP$geom), border = "grey70")
plot(nb_knn_5, coords, pch = 19, cex = 0.5, add = TRUE)
title(main = paste("Соседство по методу k-ближайших соседей (k = 5)"))
dist_k5 <- nbdists(nb_knn_5, coords)
inv_dist <- lapply(dist_k5, function(x) 1/x)
W_k5_inv <- nb2listw(nb_knn_5, glist = inv_dist, style = "W", zero.policy = TRUE)
###### 3) Ядро Епанечникова с адаптивной полосой пропускания ######
#Берем наших 5 соседей
nb_k5_ss <- include.self(nb_knn_5) 
dist_k5_ss <- nbdists(nb_k5_ss, coords)
#Считаем адаптивный радиус (максимальное расстояние до 5-го соседа для каждой точки)
max_dist <- sapply(dist_k5_ss, max) 
#Применяем формулу ядра (индивидуально для каждой строки)
epa_weights <- mapply(function(x, h) 0.75 * (1 - (x/h)^2), 
                      dist_k5_ss, max_dist, SIMPLIFY = FALSE)

#Создаем listw (ОБЯЗАТЕЛЬНО style="W" для корректности моделей)
W_kernel_manual <- nb2listw(nb_k5_ss, glist = epa_weights, style = "W")

f <- GRGDP_21 ~ IR_21 + U_21 + VR_21
####Теперь строим модели####
library(spatialreg)
library(lmtest)

# --- 1. Классическая модель OLS ---
model_ols <- lm(f, data = MAP)
summary(model_ols)

# --- 2. Модели на матрице ФЕРЗЬ (WbinW_qSpdep1) ---
sar_q <- lagsarlm(f, data = MAP, listw = WbinW_qSpdep1)
sem_q <- errorsarlm(f, data = MAP, listw = WbinW_qSpdep1)

# --- 3. Модели на матрице K=5 ближайших соседей (обратные расстояния)(W_k5_inv) ---
sar_k5 <- lagsarlm(f, data = MAP, listw = W_k5_inv)
sem_k5 <- errorsarlm(f, data = MAP, listw = W_k5_inv)

# --- 4. Модели Ядро Епанечникова с адаптивной полосой пропускания (W_kernel_manual) ---
# Добавляем zero.policy=TRUE, так как в ядерных весах структура сложнее
sar_kern <- lagsarlm(f, data = MAP, listw = W_kernel_manual, zero.policy = TRUE)
sem_kern <- errorsarlm(f, data = MAP, listw = W_kernel_manual, zero.policy = TRUE)
# Результаты для весов Ферзя
summary(sar_q)
summary(sem_q)

# Результаты для весов K=5
summary(sar_k5)
summary(sem_k5)

# Результаты для Ядерных весов
summary(sar_kern)
summary(sem_kern)
#####Хотим представить результаты в одной таблице:
install.packages("stargazer")
library(stargazer)

# Создаем сводную таблицу для всех моделей
stargazer(model_ols, sar_q, sem_q, sar_k5, sem_k5, sar_kern, sem_kern, 
          type = "text", 
          title = "Сравнение результатов моделей (Вариант 1)",
          column.labels = c("OLS", "SAR-Q", "SEM-Q", "SAR-K5", "SEM-K5", "SAR-Kern", "SEM-Kern"),
          dep.var.labels = "GRGDP_21",
          omit.stat = c("f", "ser")) # Убираем лишнюю статистику для компактности

# ПУНКТ 3: СРАВНЕНИЕ РЕЗУЛЬТАТОВ
#Значимость коэффициентов:
#VR_21: Это самый стабильный и значимый фактор во всех моделях
#U_21:значим во всех моделях, но на разных уровнях значимости
#IR_21: В данной выборке фактор оказался незначимым 

# Сравнение информационных критериев AIC и BIC
# Для корректного расчета BIC используем формулу: AIC - 2*k + log(n)*k
n_obs <- nrow(MAP)
# Количество параметров в пространственных моделях = 6 (Intercept, 3 фактора, sigma, rho/lambda)
k_spatial <- 6 
k_ols <- 5

# Создаем сводную таблицу качества моделей
comparison_metrics <- data.frame(
  Model = c("OLS", "SAR-Queen", "SEM-Queen", "SAR-K5", "SEM-K5", "SAR-Kernel", "SEM-Kernel"),
  AIC = c(AIC(model_ols), AIC(sar_q), AIC(sem_q), AIC(sar_k5), AIC(sem_k5), AIC(sar_kern), AIC(sem_kern)),
  LogLik = c(logLik(model_ols), logLik(sar_q), logLik(sem_q), logLik(sar_k5), logLik(sem_k5), logLik(sar_kern), logLik(sem_kern))
)

# Добавляем расчет BIC для каждой модели
comparison_metrics$BIC <- comparison_metrics$AIC - 2*ifelse(comparison_metrics$Model=="OLS", k_ols, k_spatial) + 
  log(n_obs)*ifelse(comparison_metrics$Model=="OLS", k_ols, k_spatial)

print("--- Сравнение информационных критериев (AIC и BIC) ---")
print(comparison_metrics)

#ВЫВОДЫ ПО ПУНКТУ 3
# Анализ информационных критериев (AIC и BIC):
# - Согласно полученной таблице, наилучшее качество подгонки демонстрируют модели 
#   на матрице «Ферзь» (Queen). Их AIC (440-442) и BIC (455-457) значительно ниже, 
#   чем у других моделей.
# - наилучшей является модель SEM-Queen, 
#   у которой минимальный AIC = 440.38 и минимальный BIC = 455.04.

# Анализ Log-Likelihood (Логарифмическая функция правдоподобия):
# - Максимальное значение LogLik (-214.19) также наблюдается у модели SEM-Queen. 
#   Это математически подтверждает, что данная спецификация лучше всего 
#   соответствует имеющимся данным по ВРП регионов.


#ПУНКТ 4 и 5: Диагностика остатков и выбор модели 
# --- 1. Диагностика OLS ---
# Проверяем гетероскедастичность
bptest(model_ols) 
# Проверяем пространственную автокорреляцию (через спец. функцию для OLS)
lm.morantest(model_ols, WbinW_qSpdep1)

# --- 2. Диагностика пространственных моделей (через Тест Морана) ---
# Мы проверяем остатки (residuals) каждой модели. 
# Если p-value > 0.05, значит модель хорошая и зависимости в остатках нет.

# Для матрицы Ферзя (Queen)
moran.test(residuals(sar_q), WbinW_qSpdep1)
moran.test(residuals(sem_q), WbinW_qSpdep1)

# Для матрицы K=5
moran.test(residuals(sar_k5), W_k5_inv)
moran.test(residuals(sem_k5), W_k5_inv)

# --- 3. Выбор модели по Элхорсту ---
lm_tests <- lm.LMtests(model_ols, WbinW_qSpdep1, test = "all")
summary(lm_tests)

# ВЫВОДЫ ПО ПУНКТУ 4: ДИАГНОСТИКА ОСТАТКОВ
# 1. Тест Бройша-Пагана (bptest) для модели OLS показал p-value = 0.335.
#    Так как p > 0.05, мы принимаем гипотезу о гомоскедастичности (постоянстве дисперсии).
#    Это означает, что оценки параметров надежны и не требуют корректировки на гетероскедастичность.

# 2. Тест Морана для остатков OLS (p-value = 0.00001569) значим.
#    Это подтверждает, что в классической модели осталась сильная пространственная 
#    автокорреляция, и использование обычного МНК (OLS) приводит к смещенным результатам.

# 3. Диагностика пространственных моделей (SAR и SEM) на разных матрицах весов:
#    - Для всех моделей (SAR-Q, SEM-Q, SAR-K5, SEM-K5) p-value теста Морана > 0.05.
#    - В частности, для лучшей модели SEM-Queen p-value = 0.3586.
#    ВЫВОД: Пространственные модели успешно устранили зависимость в остатках. 


# ВЫВОДЫ ПО ПУНКТУ 5: ОБОСНОВАНИЕ ВЫБОРА МОДЕЛИ

# Согласно комбинированному подходу Элхорста, выбор был сделан следующим образом:

# 1. ВЫБОР ВЕСОВ: Сравнение информационных критериев (Пункт 3) показало явное 
#    преимущество матрицы «Ферзь» (AIC = 440.38) над матрицами расстояний (AIC > 450).

# 2. РЕЗУЛЬТАТЫ LM-ТЕСТОВ:
#    - Базовые тесты RSerr (LM-error) и RSlag (LM-lag) оба высокозначимы (p < 0.001).
#    - В такой ситуации Элхорст рекомендует смотреть на робастные (Robust) версии.
#    - Тест adjRSerr (Robust LM-error) значим на уровне 10% (p = 0.09), в то время как 
#      тест adjRSlag (Robust LM-lag)  незначим (p = 0.55).
#    Это является аргументом в пользу модели пространственной ошибки (SEM).

# Модель SEM-Queen признана оптимальной. 
#    Она обладает минимальным AIC (440.383), устраняет пространственную 
#    автокорреляцию в остатках и показывает высокую значимость ключевых факторов: 
#    VR_21 и U_21.

# Модель SEM-Queen лучше всего описывает индекс ВРП регионов РФ, учитывая 
# внешние шоки и ненаблюдаемые факторы, общие для соседних территорий.
