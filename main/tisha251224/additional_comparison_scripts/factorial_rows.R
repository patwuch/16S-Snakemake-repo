# The factorial's report rows, shared by the LinDA and MaAsLin3 scripts.
# The three diets are peers: each row is the direct coefficient of a fit
# whose reference levels make it one, so no row depends on a chosen baseline.
DIETS <- c("Cellulose", "Inulin", "RS")
CONDS <- c("Control", "IR")

# Report rows: the fit (diet ref, condition ref) and coefficient that estimate
# each one directly.
ROWS <- data.frame(
  row    = c("Inulin vs Cellulose", "RS vs Cellulose", "RS vs Inulin",
             "Inulin vs Cellulose", "RS vs Cellulose", "RS vs Inulin",
             "IR within Cellulose", "IR within Inulin", "IR within RS",
             "Interaction: Inulin vs Cellulose", "Interaction: RS vs Cellulose",
             "Interaction: RS vs Inulin"),
  stratum = c(rep("Control", 3), rep("IR", 3), rep("", 6)),
  family = c(rep("diet_Control", 3), rep("diet_IR", 3),
             "IR_Cellulose", "IR_Inulin", "IR_RS", rep("interaction", 3)),
  d_ref = c("Cellulose", "Cellulose", "Inulin", "Cellulose", "Cellulose", "Inulin",
            "Cellulose", "Inulin", "RS", "Cellulose", "Cellulose", "Inulin"),
  c_ref = c(rep("Control", 3), rep("IR", 3), rep("Control", 6)),
  coef  = c("DietInulin", "DietRS", "DietRS", "DietInulin", "DietRS", "DietRS",
            rep("ConditionIR", 3),
            "DietInulin:ConditionIR", "DietRS:ConditionIR", "DietRS:ConditionIR"),
  stringsAsFactors = FALSE)

# The cells (diet, condition) a row's coefficient compares: two for a diet
# pair or irradiation within a diet, four for an interaction.
row_cells <- function(k) {
  d <- ROWS$d_ref[k]; cf <- ROWS$coef[k]
  other <- sub(":.*", "", sub("^Diet", "", cf))
  if (ROWS$family[k] == "interaction") return(expand.grid(diet = c(other, d), cond = CONDS, stringsAsFactors = FALSE))
  if (grepl("^Condition", cf)) return(data.frame(diet = d, cond = CONDS, stringsAsFactors = FALSE))
  data.frame(diet = c(other, d), cond = ROWS$c_ref[k], stringsAsFactors = FALSE)
}
