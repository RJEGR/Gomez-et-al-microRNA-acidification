# =============================================================================
# chemistries_TA_DIC.R  (curated copy for microRNA-red-abalone-pub)
#
# Seawater carbonate-chemistry statistics and figures (Supplementary Tables S1-S2).
#
# PROVENANCE: the carbonate-system variables (DIC, TA, CO3, Ara = Omega_ara, pCO2,
# pH) in `quimicas.csv` were computed with CO2SYS v2.0 from the measured Total
# Alkalinity and Dissolved Inorganic Carbon. This script performs only the
# downstream statistical analysis on those values (outliers, normality,
# homoscedasticity, t-tests / ANOVA, salinity normalisation after Friis et al.
# 2003, pH regression and correlation matrices).
#
# The original exploratory `seacarb` block (which recomputed the carbonate system
# inside R) has been REMOVED: CO2SYS v2.0 is the method reported in the manuscript.
# =============================================================================

# analyze chemistry ----

# test https://datatab.net/tutorial/get-started

# Se considera hacer dos 'cruzas' de comparaciones
# 1) evaluar que no existan (Ho) diferencias entre los valores de DIC y TA entre dias de muestreo en los grupos experimentales I y II y control. Para ello evaluar intercuartiles (IC) e identificar que valores si salen o no?
# Evaluar que existan (Ho) diferencias significativas entre los valores de DIC y TA entre los grupos experimentales I y II vs control

# revisar figura 1 de History of Seawater Carbonate Chemistry, Atmospheric CO2, and Ocean Acidification Richard E. Zeebe = La invasión y liberación de CO2 hacia/desde el océano cambia solo la DIC, mientras que la fotosíntesis y la respiración también cambian ligeramente la TA debido a la absorción y liberación de nitrato. La formación de CaCO3 disminuye el DIC y el TA en una proporción de 1:2 y, contradictoriamente, aumenta el [CO2], aunque la concentración de carbono inorgánico total ha disminuido. La disolución de CaCO3 tiene el efecto inverso.

rm(list = ls())

options(stringsAsFactors = FALSE)

library(tidyverse)
library(rstatix)

# quimicas.csv (CO2SYS v2.0 output) lives in the repo data folder
path <- "data/water_chemistry/"
if (!dir.exists(path)) path <- "../data/water_chemistry/"   # run from repo root or module dir

ggsavepath <- paste0(getwd(), '/Figures')

# stats.R provides helper functions (e.g. IC = 95% CI half-width). Optional:
if (file.exists(file.path(getwd(), "stats.R"))) source(file.path(getwd(), "stats.R"))
if (!exists("IC")) IC <- function(x) qt(0.975, df = length(x) - 1) * sd(x) / sqrt(length(x))


head(df <- read.csv(paste0(path, 'quimicas.csv')))

ylabs = expression(CO[3]^{-2}~(µmol~Kg^{-1}))
xlabs = expression(Omega["ara"])

# Lnames <- c(expression(Omega["ara"]),
#   expression(CO[3]^{-2}~(µmol~Kg^{-1})),
#   expression(pCO[2]~(uatm)),
#   expression(O[2]~(µmol~Kg^{-1})), "pH")


vars <- c('DIC', 'TA', 'CO3', 'Ara')

df %>% select(ID, name, vars) %>%
  pivot_longer(cols = all_of(vars), names_to = 'var') %>%
  mutate(var = factor(var, levels = vars)) -> df_longer


# 0) Test if outlier (TRUE) -----

df_longer %>% 
  ggplot(aes(y = value, x = name)) + 
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  stat_boxplot(geom ='errorbar', width = 0.07) +
  geom_boxplot(width = 0.3, outlier.alpha = 0) +
  geom_point(alpha = 0.5) +
  # stat_summary(fun=mean, geom="point", shape=23, size=1, color = 'red') +
  facet_grid(var ~., scales = 'free_y') 


df_longer %>% group_by(name, var) %>% 
  mutate(is.outlier = is_outlier(value), 
    is.extreme = is_extreme(value)) -> df_longer

# df_longer %>% group_by(name, var, is.extreme) %>% tally() %>% 
#   filter(is.extreme %in% TRUE)
df_longer %>% group_by(name, var, is.outlier) %>% tally() %>%
  filter(is.outlier %in% TRUE)

# after remove outliers, continue: 

# Find if difference between DIC and TA across experiments and control
# 1) test if gaussianity:

df_longer %>% filter(is.outlier == FALSE) %>% 
  group_by(name, var) %>% rstatix::shapiro_test(value)

# From the output, if the p-value > 0.05 implying that the distribution of the data are NOT significantly different from a normal distribution. In other words, we can assume the normality. Otherwise, if the p-value < 0.05 implying that the distribution of the data are  significantly different from a normal distribution. Therefore, we can assume not gaussianity. In other words, we argue a non parametric test to priori and posteriri analysis.

# 1) Homocelasticidad (TRUE) ----
# Before doing parametric or not test, lets to analyze homogeneity of variance across experimental Levene’s test:

df_longer %>%
  # filter(is.outlier == FALSE) %>%
  group_by(var) %>%
  levene_test(value ~ as.factor(name))

# If the p-value or significance is above  0.05, we can assume homogeneous variance based on the available data. Ex. The p value of 0.149 is greater than 0.05 so the null hypothesis is maintained and there is no difference between the variances

# check https://datatab.net/tutorial/levene-test


# Continue w/ a parametric analysis.

# Priori statistical test ----
# T test ----

df_longer %>%
  filter(is.outlier == FALSE) %>%
  group_by(var) %>%
  rstatix::pairwise_t_test(value ~ name, ref.group = 'Control') -> t_test_stats

# or tukey_hsd, porque este test ???


df_longer %>%
  filter(is.outlier == FALSE) %>%
  group_by(var) %>%
  rstatix::tukey_hsd(value ~ name)

# En los dias que yo evalue las quimicas, no se identificaron diferencias, sin embargo los datos de pH si senalan diferencias entre los sistemas experimetales. Excepto la TA del experimento II vs control. Esto puede deberse a los valores puntuales de las quimicas de agua


# Posteriori test ----

# https://www.datanovia.com/en/lessons/anova-in-r/

# ANOVA assumes that the variance of the residuals is equal for all groups. This can be checked using the 

df_longer %>%
  filter(is.outlier == FALSE) %>%
  group_by(var) %>%
  # ungroup() %>%
  rstatix::anova_test(value ~ name) -> res_aov

res_aov

# Normalize nTA and nDIC ----
# Some analyses in the following sections required normalization of data to a constant salinity. Various studies have used a simple technique to normalize data (e.g., normalization to a constant salinity of 35, such as described by Millero et al. 1998), but Friis et al. (2003) pointed out that this technique can lead to erroneous results if a nonzero y-intercept is present, as was the case in this study. We followed the approach of Friis et al. (2003) according to the following:

df %>%
  group_by(name) %>%
  mutate(nDIC = mean(Sal)*(DIC/Sal), nTA = mean(Sal)*(TA/Sal)) %>%
  select(name, nTA, nDIC, TA, DIC) %>%
  ggplot(aes(color = name)) +
  geom_point(aes(nTA, TA))
  # geom_point(aes(nDIC, DIC))

# boxplot(df$Sal)
# mean(df$Sal)

# adjust or model the DIC and TA based on the pH measures made during the experiment
out_stats <- read_rds(paste0(getwd(), '/pH_aLLdatasets_stats.rds'))

sum(df$date %in% out_stats$date) / length(df$date)

out_stats %>% ungroup() %>% distinct(name)

fill_na <- function(x,fill = 8, na.rm = FALSE) ifelse(is.na(x), fill, x)

df %>% 
  left_join(out_stats %>% select(date, name, a), 
    by = c("date", "name")) %>%
  mutate_at(c('a'), fill_na) -> df_w_pH

# Continue ----

# in comparison to pH average from the time series dataset

out_stats %>% 
  ungroup() %>%
  select(date, name, a) %>%
  pivot_wider(names_from = name, values_from = a) 

level_key <- structure(Labels, names = c('Canal-2', 'Canal-3', 'Canal-4'))

# como no hay pH de control al inicio, no vas a tener join con los datos de df
# sum(dviz$date %in% df$date)

out_stats %>% ungroup() %>%
  # mutate(name = recode_factor(name, !!!level_key)) %>%
  select(date, name, a) %>%
  inner_join(df) %>%
  ggplot(aes(a, Ara)) +
  geom_point(aes(color = name)) +
  geom_smooth(method = "lm", linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", p.accuracy = 0.001)


# Previz ----
cols <- c('DIC', 'TA', 'Ara', 'CO3')

df_w_pH %>%
  select(name, cols, a) %>%
  pivot_longer(cols = all_of(cols), names_to = 'Chem') %>%
  mutate(Chem = factor(Chem, levels = cols)) %>%
  ggplot(aes(value, a, color = name)) +
  geom_smooth(method = "lm", 
    linetype="dashed", size = 0.5, alpha=0.5, 
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", 
    p.accuracy = 0.001) +
  geom_point(size = 2.5, alpha = 0.5) +
  geom_rug(aes(color = name)) +
  facet_grid(~Chem, scales = 'free_x') +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'top') +
  labs(y = 'pH (Average)', x = '', color = '') -> psave

psave

# or cor matrix ----
# var_sel <- c('DIC', 'TA', 'Sal', 'Ara', 'CO3', 'pH_Predicted', a)

var_select <- df_w_pH %>% select_if(is.numeric) %>% names()

cor_mat_group <- function(data, var_select, var_group) {
  
  var_group <- data %>% distinct_at(var_group) %>% pull(var_group)
  
  out <- list()
  
  for(i in var_group) {
    var_group <- i
    
    data %>%
      filter(name %in% var_group) %>%
      select(all_of(var_select)) %>%
      cor_mat() -> cor_mat
    
    dist.method <- 'euclidean'
    linkage.method <- 'complete'
    
    m <- cor_mat %>% select_if(is.numeric)
    m <- m %>% mutate(across(where(is.numeric), replace_na, 0))
    
    hclust_mat <- hclust(dist(t(m), method = dist.method), 
      method = linkage.method)
    
    # hc_sam_ord <- hclust_mat$labels[hclust_mat$order]
    
    hc_sam_ord <- data.frame(var2 = hclust_mat$labels,
      order = hclust_mat$order)
    
    cor_mat %>%
      mutate(across(where(is.numeric), replace_na, 0)) %>%
      cor_reorder() %>%
      # pull_lower_triangle() %>%
      cor_gather() %>%
      left_join(hc_sam_ord) %>%
      mutate(g = var_group) %>%
      as_tibble() -> cor_mat
    
    out[[var_group]] <- cor_mat
  }
  
  out <- do.call(rbind, out)
  
  return(out)

  
}

cor_mat_group(df_w_pH, var_select, 'name') -> cor_mat

# cor_mat %>%
#   filter(g %in% 'Control') %>%
#   select(-g) %>%
#   cor_reorder() %>%
#   pull_lower_triangle() %>%
#   cor_plot(label = TRUE)

values <- corrplot::COL2("RdBu", n = 200)

# values <- structure(values, names = cor_vals)

# no veo la manera de ordenar la matriz

cor_mat %>% distinct(var1) %>% arrange(var1) %>% pull ()-> arrangeL

cor_mat %>%
  mutate(start = 'ns') %>%
  mutate(start = ifelse(p < 0.05, '*', 
    ifelse(p < 0.01, '**', start))) %>%
  # group_by(g) %>%
  mutate(var1 = factor(var1, levels = arrangeL)) %>%
  mutate(var2 = factor(var2, levels = arrangeL)) %>%
  # mutate(var2 = fct_reorder(var2, desc(order))) %>%
  # mutate(var1 = fct_reorder(var1, match(var1, var2))) %>%
  ggplot(aes(var1, var2)) +
  facet_wrap(~g) +
  geom_tile(aes(fill = cor), color = 'black') +
  scale_fill_gradientn(colours = values) +
  # geom_point(aes(color = cor), size = 5) +
  scale_color_gradientn(colours = values) +
  geom_text(aes(label = start), color = 'black')  +
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45,
      hjust = 1, vjust = 1, size = 10)) 
  
  # geom_raster(aes(var1, var2, fill = cor))

# ggsave(psave, filename = 'CO3_chem_scatterp.png', path = ggsavepath, 
  # width = 6 , height = 2.5)


df %>%
  group_by(name) %>%
  # summarise(m = mean(Sal)) %>%
  mutate(nDIC = mean(Sal)*(DIC/Sal), nTA = mean(Sal)*(TA/Sal))%>%
  ggplot(aes(nDIC, nTA, group = name)) +
  # ggplot(aes(DIC, TA, group = name)) +
  geom_smooth(method = "lm", linetype="dashed", color = 'black', se = F) +
  ggpubr::stat_cor(method = "pearson", cor.coef.name = "R", 
    p.accuracy = 0.001, label.y = 2300) +
  # ggrepel::geom_label_repel(aes(label = ID, color = pH)) +
  geom_point(size = 2.5, alpha = 0.5, aes(color = name)) +
  theme_bw(base_family = "GillSans", base_size = 10) +
  theme(legend.position = 'top') +
  labs(color = '') +
  # scale_color_manual('', values = structure(getPalette[-1], names =  Labels)) +
  facet_grid(~name) -> psave

ggsave(psave, filename = 'nTA_nDIC_scatterp.png', path = ggsavepath, 
  width = 6, height = 2.5)

#


# if normality and homogenity is true:

df_longer %>%
  filter(is.outlier == 'FALSE') %>%
  group_by(var, name) %>%
  summarise(
    a = mean(value), sd = sd(value), IC = IC(value),
    upper = a+sd, lower = a-sd, n = n()) -> df_stats

df_stats %>%
  ggplot(aes(y = a, x = name)) + 
  facet_grid( var ~., scales = 'free_y') +
  theme_bw(base_family = "GillSans", base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, 
    hjust = 1, vjust = 1, size = 10)) +
  geom_point(alpha = 0.5) +
  geom_text(aes(label = n), size = 3, vjust = -2, family = "GillSans") +
  # geom_line(aes(group = name), orientation = "x", linetype = 'dashed') +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2)

# Linear regression to predict DIC and TA (omit) -----
# Minimos cuadrados

# https://rpubs.com/joser/RegresionSimple
#  tidyverse-friendly introductory linear regression

# https://moderndive.github.io/moderndive/index.html
# remotes::install_github("moderndive/moderndive")
library(moderndive)

fits <- lme4::lmList(DIC ~ a | name, data=df_w_pH)

get_regression_table(regresion)

lm_df <- function(data, Intercept = 'DIC', val = 'a', new_set) {
  
  # new_set <- out_stats$a # vector values to be modeled
  
  form <- paste0(Intercept, "~", val)
  
  regresion <- lm(formula(form), data)
  
  # summary(regresion)
  
  # En este ejemplo la ecuación de la recta de mínimos cuadrados es:
  # y = 803.2 (Intercept) + 156.7x
  
  new_set <- data.frame(new_set)
  
  names(new_set) <- val
  
  summary(predicted <- predict(regresion, new_set))
  
  # Intervalos de confianza de la respuesta media:
  # ic es una matriz con tres columnas: la primera es la prediccion, las otras dos son los extremos del intervalo
  
  ic_conf <- predict(regresion, new_set, interval = 'confidence')
  
  # confint(regresion)
  # 
  # plot(df_w_pH$a, df_w_pH$DIC)
  # abline(regresion)
  # lines(new.set$a, ic[, 2], lty = 2)
  # lines(new.set$a, ic[, 3], lty = 2)
  # 
  # Intervalos de prediccion
  
  ic_interv <- predict(regresion, new_set, interval = 'prediction')
  
  # lines(new.set$a, ic[, 2], lty = 2, col = 'red')
  # lines(new.set$a, ic[, 3], lty = 2, col = 'red')
  # 
  # La tabla de análisis de la varianza de los errores se obtiene con el comando anova:
 
  # anova(regresion)
  
  # residuos <- rstandard(regresion)
  # valores.ajustados <- fitted(regresion)
  # plot(valores.ajustados, residuos)
  
  # No se observa ningún patrón especial, por lo que tanto la homocedasticidad como la linealidad resultan hipótesis razonables.
  
  # La hipótesis de normalidad se suele comprobar mediante un QQ plot de los residuos. El siguiente código sirve para obtenerlo:
  
  # qqnorm(residuos)
  # qqline(residuos)
  
  # Dado que los puntos están bastante alineados, la normalidad también parece aceptable.
  # summary(ic_conf)
  # summary(ic_interv)
  # 
  colnames(ic_conf) <- paste(colnames(ic_conf), 'conf', sep = '-')
  colnames(ic_interv) <- paste(colnames(ic_interv), 'inter', sep = '-')
  
  df_predicted <- new_set %>%
    mutate(predicted = c(predicted), 
      as.data.frame(ic_conf), 
      as.data.frame(ic_interv), g = Intercept)
 
  return(df_predicted)
  
}

groups <- df_w_pH %>% distinct(name) %>% pull

out <- list()

for(i in groups) {
  j <- i
  
  new_set <- out_stats %>% filter(name %in% j) %>% pull(a)
  data <- df_w_pH %>% filter(name %in% j)
  out[[j]] <- lm_df(data, Intercept = 'DIC', val = 'a', new_set)
  
}


# tienes que separar por grupos experimentales los datos!!!!

df_predicted %>%
  ggplot(aes(a, predicted_dic, color = name)) +
  geom_smooth(method = "lm", 
    linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", 
    cor.coef.name = "R", p.accuracy = 0.001) +
  geom_point()

df_w_pH %>%
  ggplot(aes(a, DIC, color = name)) +
  geom_smooth(method = "lm", 
    linetype="dashed", size = 0.5, alpha=0.5,
    se = F, na.rm = TRUE) +
  ggpubr::stat_cor(method = "pearson", 
    cor.coef.name = "R", p.accuracy = 0.001) +
  geom_point() 

