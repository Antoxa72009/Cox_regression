# --- Автоматичне визначення теки, де лежить цей скрипт ---
if (!require("rstudioapi")) install.packages("rstudioapi")
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# --- Пакети ---
if (!require("survival")) install.packages("survival")
if (!require("survminer")) install.packages("survminer")

library(survival)
library(survminer)

# --- Зчитування даних ---
df <- read.csv("SR_2.csv", sep = ",", header = TRUE)

# --- Підготовка змінних ---
df$event_bin <- ifelse(df$event == 2, 1, 0)  # 1 = подія, 0 = цензура
df$stage     <- factor(df$stage, ordered = TRUE)
df$gleason   <- factor(df$gleason, ordered = TRUE)
df$T         <- factor(df$T, ordered = TRUE)
df$epr       <- factor(df$epr)
df$radiation <- factor(df$radiation)
df$OE        <- factor(df$OE)
df$CD        <- factor(ifelse(df$CD %in% c(4, "4", "4.0"), 1, 0))

# --- Модель Кокса ---
fit <- coxph(Surv(time_months, event_bin) ~ gleason + epr + stage + T + CD + radiation + OE, data = df)
summary(fit)

# --- Kaplan–Meier (всі пацієнти) ---
km_all <- survfit(Surv(time_months, event_bin) ~ 1, data = df)
plot_all <- ggsurvplot(km_all,
                       conf.int = TRUE,
                       risk.table = TRUE,
                       xlab = "Місяці спостереження",
                       ylab = "Ймовірність без рецидиву",
                       ggtheme = theme_minimal())

# --- Kaplan–Meier по стадіях ---
km_stage <- survfit(Surv(time_months, event_bin) ~ stage, data = df)
plot_stage <- ggsurvplot(km_stage,
                         conf.int = FALSE,
                         risk.table = TRUE,
                         pval = TRUE,
                         legend.title = "Стадія",
                         xlab = "Місяці спостереження",
                         ylab = "Ймовірність без рецидиву",
                         ggtheme = theme_bw())

# --- Кумулятивна функція ризику ---
plot_cumhaz <- ggsurvplot(km_stage,
                          fun = "cumhaz",
                          legend.title = "Стадія",
                          xlab = "Місяці",
                          ylab = "Кумулятивний ризик рецидиву",
                          ggtheme = theme_light())

# --- Перевірка пропорційності ризиків ---
test.ph <- cox.zph(fit)
print(test.ph)
plot(test.ph)

# --- Збереження графіків ---
ggsave("survival_all.png", plot = print(plot_all), width = 8, height = 6, dpi = 300)
ggsave("survival_by_stage.png", plot = print(plot_stage), width = 8, height = 6, dpi = 300)
ggsave("cumhazard.png", plot = print(plot_cumhaz), width = 8, height = 6, dpi = 300)

cat("✅ Аналіз завершено. Результати збережено у поточній теці.\n")

