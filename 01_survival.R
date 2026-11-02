# Heart failure survival analysis
# Data: Heart Failure Clinical Records (UCI), 299 patients
# Methods: Kaplan-Meier, log-rank test, Cox regression

library(dplyr)
library(survival)

# ---- Load the data and check it ----
hf <- read.csv("heart_failure_clinical_records_dataset.csv")

dim(hf)          # 299 patients, 13 columns
names(hf)
sum(is.na(hf))   # 0 missing values

# DEATH_EVENT: 1 = died during follow-up, 0 = alive at the end (censored)
table(hf$DEATH_EVENT)   # 203 alive, 96 died

# ---- Kaplan-Meier: survival for all patients ----

# Surv() joins follow-up time (days) with whether the death happened
km_all <- survfit(Surv(time, DEATH_EVENT) ~ 1, data = hf)

km_all   # no median shown, since more than half of the patients survived

# Share of patients still alive at chosen days
summary(km_all, times = c(50, 100, 150, 200, 250))

# ---- Plot the overall survival curve ----

plot(km_all,
     xlab = "Follow-up time (days)",
     ylab = "Survival probability",
     main = "Kaplan-Meier survival, all patients")

# Save the same plot as a picture
dir.create("output", showWarnings = FALSE)
png("output/km_overall.png", width = 800, height = 600)
plot(km_all,
     xlab = "Follow-up time (days)",
     ylab = "Survival probability",
     main = "Kaplan-Meier survival, all patients")
dev.off()

# ---- Compare survival: low vs normal ejection fraction ----

# Two groups: ejection fraction below 40% (reduced) and 40% or above
hf <- hf %>%
  mutate(ef_group = ifelse(ejection_fraction < 40, "EF below 40", "EF 40 or above"))

# Group sizes and deaths in each group
table(hf$ef_group)
table(hf$ef_group, hf$DEATH_EVENT)

# Kaplan-Meier curves for the two groups
km_ef <- survfit(Surv(time, DEATH_EVENT) ~ ef_group, data = hf)
km_ef

# Log-rank test: do the two groups survive differently?
survdiff(Surv(time, DEATH_EVENT) ~ ef_group, data = hf)

# ---- Plot: survival by ejection fraction group ----

plot(km_ef,
     col = c("blue", "red"),
     lwd = 2,
     xlab = "Follow-up time (days)",
     ylab = "Survival probability",
     main = "Survival by ejection fraction group")

legend("bottomleft",
       legend = c("EF 40 or above", "EF below 40"),
       col = c("blue", "red"), lwd = 2)

# Save the same plot as a picture
png("output/km_by_ef.png", width = 800, height = 600)
plot(km_ef, col = c("blue", "red"), lwd = 2,
     xlab = "Follow-up time (days)", ylab = "Survival probability",
     main = "Survival by ejection fraction group")
legend("bottomleft", legend = c("EF 40 or above", "EF below 40"),
       col = c("blue", "red"), lwd = 2)
dev.off()

# ---- Cox proportional hazards regression ----

# Several factors in one model; exp(coef) is the hazard ratio
# HR above 1 means higher risk of death, below 1 means lower risk
cox_model <- coxph(Surv(time, DEATH_EVENT) ~ age + ejection_fraction +
                     serum_creatinine + serum_sodium + anaemia +
                     high_blood_pressure + diabetes + sex + smoking,
                   data = hf)

summary(cox_model)

# ---- Save the Cox results as a table ----

cox_table <- data.frame(
  factor   = names(coef(cox_model)),
  HR       = round(exp(coef(cox_model)), 2),
  lower_95 = round(exp(confint(cox_model))[, 1], 2),
  upper_95 = round(exp(confint(cox_model))[, 2], 2),
  p_value  = ifelse(summary(cox_model)$coefficients[, "Pr(>|z|)"] < 0.001,
                    "<0.001",
                    round(summary(cox_model)$coefficients[, "Pr(>|z|)"], 3))
)

cox_table
write.csv(cox_table, "output/cox_results.csv", row.names = FALSE)

# ---- Check the proportional hazards assumption ----

# Tests whether each factor's effect on risk stays steady over time.
# p below 0.05 is a warning for that factor.
zph <- cox.zph(cox_model)
zph   # GLOBAL p = 0.376; only ejection_fraction is below 0.05 (p = 0.036)

# Plot for ejection fraction, to see if its effect changes over time
plot(zph, var = "ejection_fraction")
