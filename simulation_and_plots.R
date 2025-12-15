library(rvinecopulib); library(ggplot2)
library(dplyr); library(tidyr)

set.seed(123); n <- 5000

structure <- cvine_structure(c(1, 2, 3))
indep_par <- numeric(0)

# core sim function
simulate_vine <- function(family_xz, rot_xz, par_xz, 
                          family_yz, rot_yz, par_yz,
                          family_xy_z, rot_xy_z, par_xy_z) {
  
  bicop_xz <- bicop_dist(family_xz, rot_xz, par_xz)
  bicop_yz <- bicop_dist(family_yz, rot_yz, par_yz)
  bicop_xy_z <- bicop_dist(family_xy_z, rot_xy_z, par_xy_z)
  
  vine <- vinecop_dist(list(list(bicop_xz, bicop_yz), list(bicop_xy_z)), structure)
  data <- rvinecop(n, vine)
  
  X <- data[, 1]; Y <- data[, 2]; Z <- data[, 3]
  U_X <- hbicop(cbind(X, Z), cond_var = 2, bicop_xz)
  U_Y <- hbicop(cbind(Y, Z), cond_var = 2, bicop_yz)
  
  data.frame(X = X, Y = Y, Z = Z, U_X = U_X, U_Y = U_Y)
}

# scenario configurations (1-10)
configs <- list(
  list(scen = 1, fam = "frank",   pars = list(c(1,1,1), c(3,3,1), c(5,5,1))),
  list(scen = 1, fam = "gumbel",  pars = list(c(1.5,1.5,1.5), c(2,2,1.5), c(3,3,1.5))),
  list(scen = 1, fam = "clayton", pars = list(c(1,1,1), c(3,3,1), c(5,5,1))),
  
  list(scen = 2, fam = "frank",   pars = list(c(1,1,NA), c(3,3,NA), c(5,5,NA)), cond_indep = TRUE),
  list(scen = 2, fam = "gumbel",  pars = list(c(1.5,1.5,NA), c(2,2,NA), c(3,3,NA)), cond_indep = TRUE),
  list(scen = 2, fam = "clayton", pars = list(c(1,1,NA), c(3,3,NA), c(5,5,NA)), cond_indep = TRUE),
  
  list(scen = 3, fam = "frank",      pars = list(c(-1,-1,NA), c(-3,-3,NA), c(-5,-5,NA)), cond_indep = TRUE),
  list(scen = 3, fam = "gumbel_90",  pars = list(c(1.5,1.5,NA), c(2,2,NA), c(3,3,NA)), rot = 90, cond_indep = TRUE),
  list(scen = 3, fam = "clayton_90", pars = list(c(1,1,NA), c(3,3,NA), c(5,5,NA)), rot = 90, cond_indep = TRUE),
  
  list(scen = 4, fam = "frank",    pars = list(c(1,-1,NA), c(3,-3,NA), c(5,-5,NA)), cond_indep = TRUE),
  list(scen = 4, fam = "gaussian", pars = list(c(0.3,-0.3,NA), c(0.5,-0.5,NA), c(0.7,-0.7,NA)), cond_indep = TRUE),
  
  list(scen = 5, fam = "frank",    pars = list(c(1,1,-1), c(1,1,-3), c(1,1,-5))),
  list(scen = 5, fam = "gaussian", pars = list(c(0.2,0.2,-0.2), c(0.2,0.2,-0.5), c(0.2,0.2,-0.7))),
  list(scen = 5, fam = "frank",    pars = list(c(5,5,-1), c(5,5,-0.5)), trial = TRUE),
  list(scen = 5, fam = "gaussian", pars = list(c(0.7,0.7,-0.2), c(0.7,0.7,-0.1)), trial = TRUE),
  
  list(scen = 6, fam = "frank",      pars = list(c(-1,-1,1), c(-3,-3,1), c(-5,-5,1))),
  list(scen = 6, fam = "gumbel_90",  pars = list(c(1.5,1.5,1.5), c(2,2,1.5), c(3,3,1.5)), rot_xz = 90, rot_yz = 90),
  list(scen = 6, fam = "clayton_90", pars = list(c(1,1,1), c(3,3,1), c(5,5,1)), rot_xz = 90, rot_yz = 90),
  
  list(scen = 7, fam = "frank",      pars = list(c(-1,-1,-1), c(-1,-1,-3), c(-1,-1,-5))),
  list(scen = 7, fam = "gumbel_90",  pars = list(c(1.5,1.5,1.5), c(1.5,1.5,2), c(1.5,1.5,3)), rot = 90),
  list(scen = 7, fam = "clayton_90", pars = list(c(1,1,1), c(1,1,3), c(1,1,5)), rot = 90),
  
  list(scen = 8, fam = "frank",    pars = list(c(1,-1,1), c(3,-3,1), c(5,-5,1))),
  list(scen = 8, fam = "gaussian", pars = list(c(0.2,-0.2,0.2), c(0.5,-0.5,0.2), c(0.7,-0.7,0.2))),
  
  list(scen = 9, fam = "frank",    pars = list(c(1,-1,-1), c(3,-3,-1), c(5,-5,-1))),
  list(scen = 9, fam = "gaussian", pars = list(c(0.2,-0.2,-0.2), c(0.5,-0.5,-0.2), c(0.7,-0.7,-0.2))),
  
  list(scen = 10, fam = "frank", pars = list(c(NA,NA,-1), c(NA,NA,1)), marginal_indep = TRUE),
  list(scen = 10, fam = "indep", pars = list(c(NA,NA,NA)), marginal_indep = TRUE, cond_indep = TRUE)
)

# run all simulations (scenarios 1-10)
all_data <- data.frame()
results <- data.frame()
case_counter <- list()

for (cfg in configs) {
  base_fam <- gsub("_90", "", cfg$fam)
  rot <- if (!is.null(cfg$rot)) cfg$rot else 0
  rot_xz <- if (!is.null(cfg$rot_xz)) cfg$rot_xz else rot
  rot_yz <- if (!is.null(cfg$rot_yz)) cfg$rot_yz else rot
  
  for (i in seq_along(cfg$pars)) {
    p <- cfg$pars[[i]]
    
    if (isTRUE(cfg$marginal_indep)) {
      fam_xz <- "indep"; par_xz <- indep_par; r_xz <- 0
      fam_yz <- "indep"; par_yz <- indep_par; r_yz <- 0
    } else {
      fam_xz <- base_fam; par_xz <- p[1]; r_xz <- rot_xz
      fam_yz <- base_fam; par_yz <- p[2]; r_yz <- rot_yz
    }
    
    if (isTRUE(cfg$cond_indep)) {
      fam_xy_z <- "indep"; par_xy_z <- indep_par; r_xy_z <- 0
    } else {
      fam_xy_z <- base_fam; par_xy_z <- p[3]; r_xy_z <- rot
    }
    
    df <- simulate_vine(fam_xz, r_xz, par_xz, fam_yz, r_yz, par_yz, fam_xy_z, r_xy_z, par_xy_z)
    
    scen_key <- as.character(cfg$scen)
    if (is.null(case_counter[[scen_key]])) case_counter[[scen_key]] <- 0
    case_counter[[scen_key]] <- case_counter[[scen_key]] + 1
    case_letter <- letters[case_counter[[scen_key]]]
    case_suffix <- if (isTRUE(cfg$trial)) "*" else ""
    
    display_fam <- paste0(toupper(substr(cfg$fam, 1, 1)), substr(cfg$fam, 2, nchar(cfg$fam)))
    par_xz_str <- ifelse(is.na(p[1]), "ø", p[1])
    par_yz_str <- ifelse(is.na(p[2]), "ø", p[2])
    par_xy_z_str <- ifelse(isTRUE(cfg$cond_indep), "ø", p[3])
    label <- sprintf("%s(%s,%s,%s)", display_fam, par_xz_str, par_yz_str, par_xy_z_str)
    
    rho_xy <- cor(df$X, df$Y, method = "spearman")
    rho_uv <- cor(df$U_X, df$U_Y, method = "spearman")
    tau_xy <- cor(df$X, df$Y, method = "kendall")
    tau_uv <- cor(df$U_X, df$U_Y, method = "kendall")
    
    case_id <- paste0(scen_key, case_letter, case_suffix)
    
    df$Scenario <- cfg$scen
    df$Case <- case_id
    df$Label <- label
    df$Case_Label <- paste0(case_id, ": ", label)
    df$Family <- display_fam
    df$rho_XY <- rho_xy
    df$rho_UV <- rho_uv
    df$tau_XY <- tau_xy
    df$tau_UV <- tau_uv
    
    all_data <- rbind(all_data, df)
    
    results <- rbind(results, data.frame(
      Scenario = cfg$scen,
      Case = case_id,
      Family = display_fam,
      Parameters = label,
      rho_XY = round(rho_xy, 3),
      rho_UV = round(rho_uv, 3),
      tau_XY = round(tau_xy, 3),
      tau_UV = round(tau_uv, 3)
    ))
  }
}

# scenario 11: violation of the simplifying assumption
simulate_nonsimplified <- function(par_func, family = "frank") {
  bicop_xz <- bicop_dist("frank", 0, 1)
  bicop_yz <- bicop_dist("frank", 0, 1)
  
  Z <- runif(n)
  X <- Y <- U_X <- U_Y <- numeric(n)
  
  for (i in 1:n) {
    param_z <- par_func(Z[i])
    cond_cop <- bicop_dist(family, 0, param_z)
    ST <- rbicop(1, cond_cop)
    X[i] <- hbicop(c(ST[1,1], Z[i]), cond_var = 2, bicop_xz, inverse = TRUE)
    Y[i] <- hbicop(c(ST[1,2], Z[i]), cond_var = 2, bicop_yz, inverse = TRUE)
    U_X[i] <- hbicop(c(X[i], Z[i]), cond_var = 2, bicop_xz)
    U_Y[i] <- hbicop(c(Y[i], Z[i]), cond_var = 2, bicop_yz)
  }
  
  data.frame(X = X, Y = Y, Z = Z, U_X = U_X, U_Y = U_Y)
}

nonsimplified_configs <- list(
  list(case = "11a", family = "frank",    par_func = function(z) exp(z),  label = "Frank(exp(z))"),
  list(case = "11b", family = "frank",    par_func = function(z) -exp(z), label = "Frank(-exp(z))"),
  list(case = "11c", family = "gaussian", par_func = function(z) 1 - 2*z, label = "Gaussian(1-2z)")
)

data_11 <- data.frame()
results_11 <- data.frame()

for (cfg in nonsimplified_configs) {
  df <- simulate_nonsimplified(cfg$par_func, cfg$family)
  
  rho_xy <- cor(df$X, df$Y, method = "spearman")
  rho_uv <- cor(df$U_X, df$U_Y, method = "spearman")
  tau_xy <- cor(df$X, df$Y, method = "kendall")
  tau_uv <- cor(df$U_X, df$U_Y, method = "kendall")
  
  df$Scenario <- 11
  df$Case <- cfg$case
  df$Label <- cfg$label
  df$Case_Label <- paste0(cfg$case, ": ", cfg$label)
  df$Family <- cfg$family
  df$rho_XY <- rho_xy
  df$rho_UV <- rho_uv
  df$tau_XY <- tau_xy
  df$tau_UV <- tau_uv
  
  data_11 <- rbind(data_11, df)
  
  results_11 <- rbind(results_11, data.frame(
    Scenario = 11,
    Case = cfg$case,
    Family = cfg$family,
    Parameters = cfg$label,
    rho_XY = round(rho_xy, 3),
    rho_UV = round(rho_uv, 3),
    tau_XY = round(tau_xy, 3),
    tau_UV = round(tau_uv, 3)
  ))
}

# combine all data
all_data <- rbind(all_data, data_11)
results <- rbind(results, results_11)

scenario_info <- data.frame(
  Scenario = 1:11,
  Title = c(
    "Scenario 1: Positive Confounding, Positive Conditional (+,+,+)",
    "Scenario 2: Positive Confounding, Conditional Independence (+,+,ø)",
    "Scenario 3: Negative Confounding, Conditional Independence (-,-,ø)",
    "Scenario 4: Opposite Confounding, Conditional Independence (+,-,ø)",
    "Scenario 5: Positive Confounding, Negative Conditional (+,+,-)",
    "Scenario 6: Negative Confounding, Positive Conditional (-,-,+)",
    "Scenario 7: Negative Confounding, Negative Conditional (-,-,-)",
    "Scenario 8: Opposite Confounding, Positive Conditional (+,-,+)",
    "Scenario 9: Opposite Confounding, Negative Conditional (+,-,-)",
    "Scenario 10: No Confounding (ø,ø,*)",
    "Scenario 11: Simplifying Assumption Violated"
  ),
  Short = c("1: (+,+,+)", "2: (+,+,ø)", "3: (-,-,ø)", "4: (+,-,ø)", "5: (+,+,-)",
            "6: (-,-,+)", "7: (-,-,-)", "8: (+,-,+)", "9: (+,-,-)", "10: (ø,ø,*)",
            "11: Non-simplified")
)

# SCATTER PLOTS
create_scenario_scatter <- function(data, scenario_num, scenario_title) {
  
  scen_data <- data %>% filter(Scenario == scenario_num)
  
  case_order <- scen_data %>%
    distinct(Case, Case_Label) %>%
    arrange(Case) %>%
    pull(Case_Label)
  
  plot_data <- scen_data %>%
    select(Case_Label, X, Y, U_X, U_Y, Z, rho_XY, rho_UV, tau_XY, tau_UV) %>%
    pivot_longer(cols = c(X, U_X), names_to = "x_type", values_to = "x_val") %>%
    pivot_longer(cols = c(Y, U_Y), names_to = "y_type", values_to = "y_val") %>%
    filter((x_type == "X" & y_type == "Y") | (x_type == "U_X" & y_type == "U_Y")) %>%
    mutate(
      Type = ifelse(x_type == "X", "Marginal (X, Y)", "Partial (U_X, U_Y)"),
      Type = factor(Type, levels = c("Marginal (X, Y)", "Partial (U_X, U_Y)")),
      Case_Label = factor(Case_Label, levels = case_order)
    )
  
  annot_data <- scen_data %>%
    distinct(Case_Label, rho_XY, rho_UV, tau_XY, tau_UV) %>%
    pivot_longer(
      cols = c(rho_XY, rho_UV, tau_XY, tau_UV),
      names_to = "measure",
      values_to = "value"
    ) %>%
    mutate(
      Type = ifelse(grepl("UV", measure), "Partial (U_X, U_Y)", "Marginal (X, Y)"),
      Type = factor(Type, levels = c("Marginal (X, Y)", "Partial (U_X, U_Y)")),
      Coef = ifelse(grepl("rho", measure), "rho", "tau"),
      Case_Label = factor(Case_Label, levels = case_order)
    ) %>%
    pivot_wider(names_from = Coef, values_from = value) %>%
    mutate(label = sprintf("ρ=%.2f, τ=%.2f", rho, tau))
  
  ggplot(plot_data, aes(x = x_val, y = y_val)) +
    geom_point(aes(color = Z), alpha = 0.4, size = 0.6) +
    geom_smooth(method = "lm", se = FALSE, color = "black", linetype = "dashed", linewidth = 0.8) +
    geom_text(data = annot_data, aes(x = 0.95, y = 0.05, label = label),
              hjust = 1, vjust = 0, size = 2.5, color = "black") +
    scale_color_viridis_c(option = "plasma", name = "Z") +
    facet_grid(Case_Label ~ Type, scales = "fixed") +
    labs(title = scenario_title, x = NULL, y = NULL) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) +
    theme_minimal(base_size = 9) +
    theme(
      strip.text = element_text(face = "bold", size = 8),
      strip.text.y = element_text(angle = 0, hjust = 0),
      legend.position = "right",
      plot.title = element_text(face = "bold", size = 10, hjust = 0.5),
      panel.grid.minor = element_blank()
    )
}

# scatter plots for each scenario
for (s in 1:11) {
  title <- scenario_info$Title[scenario_info$Scenario == s]
  p <- create_scenario_scatter(all_data, s, title)
  print(p)
}

# DUMBBELL PLOT
results_dumbbell <- results %>%
  pivot_longer(
    cols = c(rho_XY, rho_UV, tau_XY, tau_UV),
    names_to = "measure",
    values_to = "value"
  ) %>%
  mutate(
    Type = ifelse(grepl("UV", measure), "Partial", "Marginal"),
    Coefficient = ifelse(grepl("rho", measure), "Spearman (ρ)", "Kendall (τ)")
  ) %>%
  pivot_wider(
    id_cols = c(Scenario, Case, Family, Parameters, Coefficient),
    names_from = Type, 
    values_from = value
  ) %>%
  mutate(
    Scenario_label = factor(
      scenario_info$Short[match(Scenario, scenario_info$Scenario)],
      levels = scenario_info$Short
    ),
    Case_Label = paste0(Case, ": ", Parameters)
  )

results_dumbbell <- results_dumbbell %>%
  arrange(Scenario, Case) %>%
  group_by(Scenario) %>%
  mutate(Case_Label = factor(Case_Label, levels = rev(unique(Case_Label)))) %>%
  ungroup()

ggplot(results_dumbbell, aes(y = Case_Label)) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_segment(aes(x = Marginal, xend = Partial, yend = Case_Label, linetype = Coefficient), 
               color = "grey50", linewidth = 0.5) +
  geom_point(aes(x = Marginal, color = "Marginal", shape = Coefficient), size = 2.5) +
  geom_point(aes(x = Partial, color = "Partial", shape = Coefficient), size = 2.5) +
  facet_wrap(~ Scenario_label, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("Marginal" = "steelblue", "Partial" = "coral")) +
  scale_shape_manual(values = c("Spearman (ρ)" = 16, "Kendall (τ)" = 17)) +
  scale_linetype_manual(values = c("Spearman (ρ)" = "solid", "Kendall (τ)" = "dotted")) +
  labs(
    x = "Correlation", 
    y = NULL, 
    color = "Type",
    shape = "Coefficient",
    linetype = "Coefficient",
    title = "Marginal vs Partial Correlations Across Scenarios"
  ) +
  theme_minimal(base_size = 9) +
  theme(
    strip.text = element_text(face = "bold"),
    legend.position = "bottom",
    legend.box = "horizontal",
    plot.title = element_text(face = "bold", hjust = 0.5),
    axis.text.y = element_text(size = 7)
  ) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2),
    linetype = guide_legend(order = 2)
  )