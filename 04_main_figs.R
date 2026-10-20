# =============================================================================
# Main-text figures, built from the flattened tables written by 03.

# Run after 01_packages_data_functions.R, which attaches the packages and defines all.dat / the model parameters.
# =============================================================================

all.dat = read_rds("./data/all_tick_dat.rds")

# =================================================
# Fig 2  Seasonality of H. longicornis life stages
# =================================================

all.dat %>% 
  pivot_longer(3:5,
               names_to = "stage",
               values_to = "count") %>% 
  mutate(date = as.Date(t - 1, origin = "2020-01-01"),
         month = as.numeric(format(date, "%m"))) %>% 
  group_by(clean.var, Lat, stage, month) %>% 
  summarise(count = sum(count), .groups = "drop") %>% 
  group_by(clean.var, Lat, stage) %>% 
  mutate(count = scales::rescale(count)) %>% 
  ungroup() %>% 
  mutate(stage = factor(stage, levels = c("Nymph", "Adult", "Larva"))) %>% 
  ggplot(aes(month, reorder(clean.var, Lat), fill = count)) +
  geom_tile(color = "grey50", linewidth = 0.2) +
  facet_wrap(~stage, ncol = 3) +
  scale_fill_gradientn(colors = c("#f7fbff", "#6baed6", "#08306b"),
                       name = "Std. density") +
  scale_x_continuous(expand = c(0, 0),
                     breaks = 1:12) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = "Month", y = "") +
  theme_bw() +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 15),
        axis.text.x = element_text(size = 11.1),
        axis.text.y = element_text(size = 12.1),
        axis.title = element_text(size = 13)) +
  coord_fixed()

# ggsave("./Visualization/Fig2.pdf", height = 8, width = 11.5)

# =================================================
# Fig 3  Distribution of R0 across the 27 scenarios
# =================================================
LOOP_SPLITS = str_c("./output/all.loop.slim_", 1:3, ".rds")


# R0 of every posterior draw in all 27 scenarios - used only by the violin plot
R0.all = map_dfr(LOOP_SPLITS, ~ read_rds(.x) %>%
                   select(sce, ratio, kappa, R0.origin))

# the loop table restricted to the baseline dilution rate kappa = 0.15
all.res = map_dfr(LOOP_SPLITS, ~ read_rds(.x) %>% filter(kappa == 0.15))

all.res$Var %>% unique()

R0.all %>%
  mutate(
    sce = factor(sce, 
                 levels = c("low", "medium", "high"),
                 labels = c("Low", "Medium", "High")),
    ratio = factor(ratio,
                   levels = c("200:10:1", "100:10:1", "52:5:1")),
    kappa = case_when(kappa == 0.15 ~ "~kappa~'= 0.15'",
                      kappa == 0.05 ~ "~kappa~'= 0.05'",
                      kappa == 0.5 ~ "~kappa~'= 0.5'"),
    kappa = factor(kappa, levels = c("~kappa~'= 0.15'",
                                     "~kappa~'= 0.05'",
                                     "~kappa~'= 0.5'"))) %>%
  ggplot(aes(x = sce, y = R0.origin, fill = ratio)) +
  geom_violin(col = "grey",
              linewidth = 1,
              alpha = 0.5,
              trim = FALSE,
              position = position_dodge(width = 0.8)) +
  geom_hline(yintercept = 1, lty = "dotted") +
  geom_boxplot(position = position_dodge(width = 0.8),
               linewidth = 0.65,
               notch = TRUE,
               alpha = 0.5,
               width = 0.2,
               col = "black",
               outlier.shape = NA) +
  facet_wrap(~ kappa, 
             labeller = label_parsed,
             nrow = 1) +
  scale_fill_manual(values = c("#4DBBD5", "#E64B35", "#00A087"), name = "L:N:A") +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 16)) +
  labs(x = NULL, y = expression(R[0]))

# ggsave("./Visualization/Fig3.pdf", height = 4.2, width = 10)


# ====================================================
# Fig 4  Elasticity of R0
#   A  relative contribution of the three transmission pathways
#   B  49 individual matrix elements, mean with 2.5 and 97.5 percentiles.
# ====================================================
all.eij.res = read_rds("./output/all.eij.wide.rds")

ela.raw =  all.eij.res %>% 
  filter(lna == "100:10:1", kappa == 0.15) %>%
  as_tibble() %>%
  pivot_longer(e11:e77,
               names_to = "Eij",
               values_to = "elasticity") %>%
  mutate(Var = ifelse(Var == "new Jersey 2020_2020", "new Jersey 2020", Var),
         Var = str_trim(Var))

ela.sum = ela.raw %>%
  group_by(sce, Eij) %>%
  summarise(mean = mean(elasticity),
            low  = quantile(elasticity, 0.025),
            high = quantile(elasticity, 0.975),
            .groups = "drop")

fmt = function(x) sub("^(-?)0\\.", "\\1.", sprintf("%.2f", x))

heat.dat = ela.sum %>%
  mutate(
    row = factor(str_sub(Eij, 2, 2), levels = 7:1),
    col = factor(str_sub(Eij, 3, 3), levels = 1:7),
    sce = factor(sce,
                 levels = c("low", "medium", "high"),
                 labels = c("Low", "Medium", "High")),
    fill = case_when(
      mean == 0 ~ NA_real_,
      round(mean, 2) == 0 ~ 0,
      TRUE ~ mean),
    label = case_when(
      mean == 0 ~ NA_character_,
      round(mean, 2) == 0 ~ "0",
      TRUE ~ paste0(fmt(mean), "\n(", fmt(low), ", ", fmt(high), ")")))


P1 = ggplot(heat.dat, aes(col, row, fill = fill)) +
  geom_tile(colour = "white",
            linewidth = 0.6,
            alpha = 0.75) +
  geom_text(aes(label = label),
            size = 3.2,
            lineheight = 0.9) +
  scale_fill_gradientn(
    colours = c("#ffffff", "#eef5fc", "#dceef9", "#c6e2f3", "#a9d2ec", "#8ec3e6"),
    trans = "sqrt", na.value = "grey92",
    name = "Elasticity") +
  facet_wrap(~ sce, nrow = 3) +
  theme_classic() +
  scale_x_discrete(labels = c("TIE", "TIL", "TIN", "TIA", "HIL", "HIN", "HIA")) +
  scale_y_discrete(labels = rev(c("TIE", "TIL", "TIN", "TIA", "HIL", "HIN", "HIA"))) +
  theme(strip.background = element_blank(),
        strip.text = element_text(size = 16),
        axis.text = element_text(size = 12.5),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 12.5),
        legend.position = "bottom",
        legend.margin = margin(t = -30, r = 0, b = 0, l = 0),
        legend.box.spacing = grid::unit(0, "pt"),
        legend.box.margin = margin(0, 0, 0, 0),
        plot.margin = margin(t = 5, r = 20, b = -5, l = 5),
        panel.border = element_rect(colour = "black",
                                    fill = NA,
                                    linewidth = 0.8),
        axis.ticks = element_blank()) +
  labs(x = NULL, y = "")


path_ela = ela.raw %>%
  mutate(T3 = case_when(
    Eij %in% c("e11", "e12", "e13", "e14") ~ "Transovarial",
    Eij %in% c("e21", "e22", "e23",
               "e31", "e32", "e33",
               "e41", "e42", "e43") ~ "Non-systemic",
    TRUE ~ "Systemic")) %>%
  group_by(Var, sce, j, T3) %>%
  summarise(elasticity = sum(elasticity), .groups = "drop") %>%
  group_by(sce, T3) %>%
  summarise(mean = mean(elasticity),
            low  = quantile(elasticity, 0.025),
            high = quantile(elasticity, 0.975),
            .groups = "drop")

P2 = path_ela %>%
  mutate(sce = factor(sce, levels = c("low", "medium", "high"),
                      labels = c("Low", "Medium", "High")),
         T3 = factor(T3, levels = c("Transovarial", "Non-systemic", "Systemic"))) %>%
  ggplot(aes(x = T3, y = mean, fill = T3)) +
  geom_col(width = 0.7, color = "black", linewidth = 0.3) +
  geom_errorbar(aes(ymin = low, ymax = high),
                width = 0.1, linewidth = 0.5) +
  facet_wrap(~ sce, nrow = 3) +
  scale_fill_npg() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
  theme_classic(base_size = 12) +
  theme(legend.position = "none",
        strip.background = element_blank(),
        strip.text  = element_text(size = 16),
        panel.border = element_rect(color = "black",
                                    fill = NA,
                                    linewidth = 0.8),
        axis.text    = element_text(size = 14),
        axis.text.x  = element_text(size = 13, 
                                    angle = 35, 
                                    hjust = 1),
        axis.title   = element_text(size = 16),
        plot.margin  = margin(t = 5, r = 5, b = -5, l = 5)) +
  labs(x = NULL,
       y = "Relative contribution")


fig = cowplot::plot_grid(P2, P1, 
                         ncol = 2,
                         rel_widths = c(1, 1.4),
                         align = "hv",
                         labels = c("A", "B"))
fig

# ggsave("./Visualization/Fig4.pdf", fig, height = 8.5, width = 10)


# ===================================================
# Fig 5  Elasticities of transmission loops
# Loop value = (lower + upper)/2, i.e. the Adams midpoint of the elasticity decomposition. Loops are ranked by their mean pooled over the three abundance levels and only those above 0.025 are kept
# ===================================================
top_loops = all.res %>%
  filter(kappa == 0.15 & ratio == "100:10:1") %>% 
  mutate(loop.elas = (loop.elas.low + loop.elas.high)/2) %>% 
  group_by(loop) %>%
  summarise(global_mean = mean(loop.elas, na.rm = TRUE)) %>%
  filter(global_mean > 0.025) %>% 
  pull(loop)

loop_100_10_1 = all.res %>%
  filter(kappa == 0.15 & ratio == "100:10:1") %>% 
  mutate(loop.elas = (loop.elas.low + loop.elas.high)/2) %>% 
  group_by(loop, Var, sce) %>% 
  summarise(mean = mean(loop.elas, na.rm = TRUE))  %>%
  filter(loop %in% top_loops) %>% 
  mutate(sce = factor(sce, levels = c("low", "medium", "high"), labels = c("Low", "Medium", "High")))


loop_100_10_1 %>% 
  left_join(distinct(all.dat[c(1, 6)]), by = c("Var" = "clean.var")) %>% 
  ggplot(aes(reorder(loop, -mean), reorder(Var, Lat), fill = mean)) +
  geom_tile(color = "white", linewidth = 0.2) +
  geom_text(aes(label = ifelse(is.na(mean), "", 
                               round(mean, 2))),
            size = 3,         
            color = "black") +   
  scale_fill_gradientn(
    colors = c("white", "#FDB863", "#E31A1C"),
    na.value = "transparent",   # missing site-year -> leave the cell empty
    name = "Elasticity"
  ) +
  facet_wrap(~sce, ncol = 3,
             scales = "free_x") +
  # coord_fixed() +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5,
                               hjust = 1, size = 14),
    axis.text.y = element_text(size = 11),
    panel.border = element_rect(color = "black", 
                                fill = NA, 
                                linewidth = 0.8),
    panel.grid = element_blank(),     
    legend.position = "right",
    strip.text = element_text(size = 16)) +
  labs(x = "", y = "")

# ggsave("./Visualization/Fig5.pdf", height = 9.5, width = 11.5)


# ===================================================
# Phenological overlap indices (used by Fig 6 and Fig 7)
# ===================================================
# Tick activity is re-weighted to the 100:10:1 larval:nymph:adult ratio and normalised into a daily activity profile a_L(t), a_N(t), a_A(t) that sums to 1 over the 365-day cycle. Each O_* index is the area under the pointwise
# minimum of two such profiles (trapezoid rule), i.e. the fraction of the year in which two stages are simultaneously active:
# self-overlap:  O_LL, O_NN, O_AA (profile squared)
# cross-stage:  O_NL, O_LA, O_NA, O_LNA   (two or three stages)

all_ts  = read_rds("./data/all_tick_dat.rds") %>% 
  rename("Var" = clean.var)

compute_shape = function(df) {
  df %>%
    mutate(total = sum(Nymph + Larva + Adult),
           sum_N = sum(Nymph), sum_L = sum(Larva), sum_A = sum(Adult),
           wN = 10/111, wL = 100/111, wA = 1/111,
           target_N = total*wN, target_L = total*wL, target_A = total*wA,
           kN = target_N/sum_N, kL = target_L/sum_L, kA = target_A/sum_A,
           Nymph = Nymph*kN, Larva = Larva*kL, Adult = Adult*kA) %>%
    mutate(a_N = Nymph/sum(Nymph),
           a_L = Larva/sum(Larva),
           a_A = Adult/sum(Adult)) %>%
    select(Var, t, a_L, a_N, a_A)
}

shape.dat = all_ts %>%
  group_by(Var) %>%
  slice(1:365) %>%
  group_split() %>%
  map_dfr(compute_shape)

# trapezoid integration over the 365 daily values (no endpoints)
trap = function(x) sum((x[-1] + x[-length(x)]) / 2)

ind = shape.dat %>%
  group_by(Var) %>%
  summarise(
    O_NL  = trap(pmin(a_L, a_N)),
    O_LA  = trap(pmin(a_L, a_A)),
    O_NA  = trap(pmin(a_N, a_A)),
    O_LNA = trap(pmin(a_L, pmin(a_N, a_A))),
    O_LL  = trap(a_L^2),
    O_NN  = trap(a_N^2),
    O_AA  = trap(a_A^2),
    .groups = "drop"
  )

ind_cols = c("O_LL", "O_NN", "O_AA", "O_NL", "O_LA", "O_NA", "O_LNA")
ind_labs = c("O[LL]", "O[NN]", "O[AA]", "O[L*N]", "O[LA]", "O[N*A]", "O[LNA]")

SELF_IND = c("O_LL", "O_NN", "O_AA") 
CROSS_IND = c("O_NL", "O_LA", "O_NA", "O_LNA")

# baseline cell only: every location-year x abundance level, all draws
raw = all.res %>%
  filter(ratio == "100:10:1", kappa == 0.15) %>%
  mutate(Var = ifelse(Var == "new Jersey 2020_2020", "new Jersey 2020", Var),
         Var = str_trim(Var))

# one R0 per site-year and abundance level (median over the 1000 draws)
R0_df = raw %>%
  group_by(Var, sce) %>%
  summarise(R0 = median(R0.origin, na.rm = TRUE), .groups = "drop")

top4 = c("1->2->1", "1->2->3->1", "1->3->1", "1->4->1")

loop_df = raw %>%
  filter(loop %in% top4) %>%
  mutate(loop.elas = (loop.elas.low + loop.elas.high) / 2) %>%
  group_by(loop, Var, sce) %>%
  summarise(elas = mean(loop.elas, na.rm = TRUE), .groups = "drop") %>%
  mutate(loop = factor(loop, levels = top4))

R0_df   = R0_df   %>% left_join(ind, by = "Var")
loop_df = loop_df %>% left_join(ind, by = "Var")
stopifnot(all(ind_cols %in% names(R0_df)), all(ind_cols %in% names(loop_df)))


# Controls used when testing a cross-stage overlap: such an overlap always involves two (or three) stages, so it is adjusted for the concentration of each of them.
CTRL_MAP = list(
  O_NL  = c("O_LL", "O_NN"),
  O_LA  = c("O_LL", "O_AA"),
  O_NA  = c("O_NN", "O_AA"),
  O_LNA = c("O_LL", "O_NN", "O_AA"),
  O_LL  = character(0),
  O_NN  = character(0),
  O_AA  = character(0)
)

# Partial correlation via added-variable residuals: regress both y and x on the controls, correlate the two residual vectors, and use df = n - 2 - k, with k the number of controls. Returns NA for self-overlaps (k = 0).
partial_cor = function(dat, y, x, ctrl) {
  if (length(ctrl) == 0) return(tibble(r = NA_real_, p = NA_real_))
  fy = as.formula(paste(y, "~", paste(ctrl, collapse = " + ")))
  fx = as.formula(paste(x, "~", paste(ctrl, collapse = " + ")))
  ry = residuals(lm(fy, data = dat))
  rx = residuals(lm(fx, data = dat))
  dfree = nrow(dat) - 2 - length(ctrl)
  r = unname(cor(rx, ry))
  t = r * sqrt(dfree / (1 - r^2))
  p = 2 * pt(-abs(t), dfree)
  tibble(r = r, p = p)
}

# For one response (R0 or one loop elasticity) and each of the three abundance levels, compute both the marginal and the partial correlation with all 7 overlap indices. r_final keeps the marginal value for self-overlaps and the adjusted value for cross-stage overlaps.

build_one = function(d, yvar) {
  d = d %>% mutate(sce = as.character(sce))
  expand_grid(sce = unique(d$sce), ind = ind_cols) %>%
    mutate(
      marg = map2(sce, ind, ~ {
        g = filter(d, sce == .x)
        ct = cor.test(g[[.y]], g[[yvar]])
        tibble(r = unname(ct$estimate), p = ct$p.value)
      }),
      part = map2(sce, ind, ~ {
        g = filter(d, sce == .x)
        partial_cor(g, yvar, .y, CTRL_MAP[[.y]])
      })
    ) %>%
    unnest_wider(marg, names_sep = "_") %>%
    unnest_wider(part, names_sep = "_") %>%
    rename(r_marg = marg_r, p_marg = marg_p,
           r_part = part_r, p_part = part_p) %>%
    # reported coefficient: marginal for self-overlaps, adjusted otherwise
    mutate(is_self  = ind %in% SELF_IND,
           r_final  = ifelse(is_self, r_marg, r_part),
           p_final  = ifelse(is_self, p_marg, p_part),
           ind      = factor(ind, levels = ind_cols, labels = ind_labs))
}


loop_corr = loop_df %>%
  group_by(loop) %>%
  group_modify(~ build_one(.x, "elas")) %>%
  ungroup()
stopifnot("loop" %in% names(loop_corr))

mk_lab = function(d, rcol, pcol) {
  r = d[[rcol]]; p = d[[pcol]]
  d %>% mutate(
    sig = case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ ""),
    lab = ifelse(is.na(r), "", paste0(sprintf("%.2f", r), sig))
  )
}

loop_final = mk_lab(loop_corr, "r_final", "p_final")

sce_levels = c("low", "medium", "high")
sce_labs   = c("Low", "Medium", "High")


# ===================================================
# Fig 6  R0 against phenological overlap indices
#   A  self-overlaps (O_LL, O_NN, O_AA): marginal correlation, raw scale.
#   B  cross-stage overlaps (O_NL, O_LA, O_NA, O_LNA): partial correlation,
# ===================================================

ov.base = R0_df %>% select(Var, sce, R0) %>% left_join(ind, by = "Var")

ov.dat = expand_grid(sce = unique(ov.base$sce), idx = ind_cols) %>%
  mutate(res = map2(sce, idx, ~ {
    g    = ov.base %>% filter(sce == .x) %>% drop_na()  
    ctrl = CTRL_MAP[[.y]]
    if (length(ctrl) == 0) {                      # self-overlap: raw values + Pearson
      tibble(Var = g$Var,
             x = g[[.y]], y = g$R0,
             r = unname(cor(g[[.y]], g$R0)),
             p = cor.test(g[[.y]], g$R0)$p.value)
    } else {                                      # cross-stage: residuals
      fy    = as.formula(paste("R0 ~", paste(ctrl, collapse = " + ")))
      fx    = as.formula(paste(.y, " ~", paste(ctrl, collapse = " + ")))
      ry    = residuals(lm(fy, data = g))
      rx    = residuals(lm(fx, data = g))
      r     = unname(cor(rx, ry))
      dfree = nrow(g) - 2 - length(ctrl)
      t     = r * sqrt(dfree / (1 - r^2))
      tibble(Var = g$Var, x = rx, y = ry, r = r, p = 2 * pt(-abs(t), dfree))
    }
  })) %>%
  unnest(res) %>%
  mutate(idx_raw = idx,
         idx     = factor(idx, levels = ind_cols, labels = ind_labs),
         sce     = factor(sce, levels = c("low", "medium", "high"),
                          labels = c("Low", "Medium", "High")))


ov.lab = ov.dat %>%
  distinct(sce, idx, idx_raw, r, p) %>%
  mutate(star = case_when(p < 0.001 ~ "***",
                          p < 0.01  ~ "**",
                          p < 0.05  ~ "*",
                          TRUE      ~ ""),
         lab  = case_when(
           p < 1e-4 ~ paste0("r = ", sprintf("%.2f", r), ", p < 0.0001", star),
           p < 1e-3 ~ paste0("r = ", sprintf("%.2f", r), ", p < 0.001",  star),
           TRUE     ~ paste0("r = ", sprintf("%.2f", r), ", p = ",
                             sprintf("%.3g", p), star)))


# p-value colour scale: pale = not significant, dark = significant
blues = RColorBrewer::brewer.pal(9, "Blues")[9:4]

ov_col = scale_color_gradientn(colours = blues, limits = c(0, 1), trans = "sqrt",
                               breaks = c(0.05, 0.5, 1),
                               labels = c("0.05", "0.50", "1.00"),
                               name   = "p value",
                               guide  = guide_colourbar(barwidth  = unit(0.38, "cm"),
                                                        barheight = unit(4.2, "cm"),
                                                        ticks = FALSE))
ov_fill = scale_fill_gradientn(colours = blues, limits = c(0, 1), trans = "sqrt",
                               name = "p value", guide = "none")

ov_theme = theme_bw(base_size = 11) +
  theme(panel.grid = element_blank(),
        panel.spacing.x = unit(0.5, "cm"),
        panel.spacing.y = unit(0.3, "cm"),
        axis.title = element_text(size = 14),
        axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 11),
        strip.text = element_text(size = 13),
        strip.background = element_rect(fill = NA, colour = "black"),
        legend.position = "right",
        legend.title = element_text(size = 12),
        legend.text  = element_text(size = 10))


# panel A: self-overlaps, marginal correlation, raw scale
dA   = ov.dat %>% filter(idx_raw %in% SELF_IND)
labA = ov.lab %>% filter(idx_raw %in% SELF_IND)

pA = ggplot(dA, aes(x, y)) +
  geom_point(aes(color = p), size = 2) +
  geom_smooth(aes(color = p, fill = p), method = "lm", formula = y ~ x,
              se = TRUE, linewidth = 0.7, alpha = 0.15) +
  ov_col + ov_fill +
  geom_text(data = labA, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.06, vjust = 1.5, size = 2.2, inherit.aes = FALSE) +
  facet_grid(sce ~ idx, scales = "free",
             labeller = labeller(idx = label_parsed, sce = label_value)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3),
                     expand = expansion(mult = 0.12)) +
  labs(x = "Activity concentration",
       y = expression(R[0])) +
  ov_theme

# panel B: cross-stage overlaps, partial correlation, both axes are residuals
dB   = ov.dat %>% filter(idx_raw %in% CROSS_IND)
labB = ov.lab %>% filter(idx_raw %in% CROSS_IND)

pB = ggplot(dB, aes(x, y)) +
  geom_point(aes(color = p), size = 2) +
  geom_smooth(aes(color = p, fill = p), method = "lm", formula = y ~ x,
              se = TRUE, linewidth = 0.7, alpha = 0.15) +
  ov_col + ov_fill +
  geom_text(data = labB, aes(x = -Inf, y = Inf, label = lab),
            hjust = -0.06, vjust = 1.5, size = 2.2, inherit.aes = FALSE) +
  facet_grid(sce ~ idx, scales = "free",
             labeller = labeller(idx = label_parsed, sce = label_value)) +
  scale_x_continuous(breaks = scales::breaks_pretty(n = 3),
                     expand = expansion(mult = 0.12)) +
  labs(x = "Overlap index residual",
       y = expression(R[0] ~ "residual")) +
  ov_theme

# combine: keep a single p-value colour bar to the right of panel B
ov_leg = cowplot::get_legend(pB)

ov_combo = cowplot::plot_grid(
  pA + theme(legend.position = "none"),
  pB + theme(legend.position = "none"),
  ov_leg,
  nrow = 1,
  rel_widths = c(3, 4, 0.5),
  align = "h", axis = "b",
  labels = c("A", "B", ""), label_size = 15)

ov_combo

# ggsave("./Visualization/Fig6.pdf", ov_combo, width = 13, height = 5.2)


# ===================================================
# Fig 7  Loop elasticities and their relation to overlap indices
# A  violin + box plot of the loops above the 0.025 cut-off
# B  Pearson r between each loop elasticity and the 7 overlap indices.
# ===================================================
heat_theme = theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 11),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(size = 14),
        strip.background = element_rect(fill = NA, colour = "black"),
        legend.position = "right")

fill_scale = scale_fill_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
                                  midpoint = 0, limits = c(-1, 1),
                                  name = "Pearson r",
                                  breaks = c(-1, -0.5, 0, 0.5, 1))


loop_cols = c("#A6CEE3", "#A20056FF", "#FDBF6F", "#CAB2D6",  
              "#FFFF99", "#FB9A99", "#FCCDE5", "#80B1D3",  
              "#B3DE69")


PB1 = loop_100_10_1 %>% 
  mutate(loop = factor(loop, levels = levels(top_loops))) %>% 
  ggplot(aes(x = loop, y = mean, fill = loop)) +   
  geom_violin(col = "grey",
              linewidth = 0.7,
              alpha = 0.9,
              trim = FALSE) +
  geom_boxplot(linewidth = 0.5,
               alpha = 0.5,
               width = 0.15,
               col = "black",
               outlier.shape = NA,
               position = position_dodge(width = 0.9)) +
  scale_fill_manual(values = loop_cols) +     
  facet_wrap(~ sce, ncol = 3,
             scales = "free_y") +
  scale_y_sqrt() +
  theme_classic() +
  theme(legend.position = "none",
        axis.title = element_text(size = 15),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 90, 
                                   hjust = 1),
        strip.background = element_rect(fill = NA, colour = "black"),
        strip.text = element_text(size = 15)) +
  labs(x = "", y = "Loop elasticity")

print(PB1)

SEP = 3.5

PB2 = loop_final %>%
  mutate(sce = factor(sce, levels = sce_levels, labels = sce_labs)) %>%
  mutate(loop = factor(loop, levels = rev(top_loops))) %>% 
  ggplot(aes(ind, loop, fill = r_final)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_vline(xintercept = SEP, colour = "grey30", linewidth = 0.7, linetype = "dashed") +
  geom_text(aes(label = lab), size = 3.5) +
  fill_scale +
  facet_wrap(~ sce, nrow = 1) +
  scale_x_discrete(labels = function(x) parse(text = x)) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 12),
        panel.grid = element_blank(),
        strip.text = element_text(size = 15),
        strip.background = element_rect(fill = NA, colour = "black"),
        legend.position = "bottom",
        plot.margin = margin(t = -15, b =0, r = 0, l = 0)) + 
  labs(x = NULL, y = NULL)

plot_grid(PB1, PB2,
          ncol = 1, 
          rel_heights = c(1.4, 1),
          align = "v",
          labels = c("A", "B"),
          label_size = 17)

# ggsave("./Visualization/Fig.7.pdf", height = 8, width = 12.5)