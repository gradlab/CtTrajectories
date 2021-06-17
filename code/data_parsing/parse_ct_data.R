# =============================================================================
# Import
# =============================================================================

library(tidyverse) 

ct_dat_clean <- read.csv("data/ct_dat_clean.csv", stringsAsFactors=FALSE) %>% as_tibble()

# Visualize a regression of Yale values on Florida values: 
fig_yale_florida_regression <- ct_dat_clean %>% 
	filter(!is.na(CT.T1) & !is.na(N1_CT_Value)) %>% 
	ggplot(aes(x=CT.T1, y=N1_CT_Value)) + 
		geom_abline(intercept=0, slope=1, linetype="dashed", col="grey") + 
		geom_point() + 
		geom_line(stat="smooth", method="lm", fullrange=TRUE) + 
		theme_minimal() + 
		theme(text=element_text(size=14)) + 
		xlim(min(ct_dat_clean$CT.T1, na.rm=TRUE), max(ct_dat_clean$CT.T1, na.rm=TRUE)) + 
		labs(x="Florida Ct (T1)", y="Yale Ct (N1)")

# write_csv(ct_dat_clean %>% 
# 	filter(!is.na(CT.T1) & !is.na(N1_CT_Value)), "figure_data/FigS15/figs15.csv")

# Run the regression:
ct_reg_fit <- ct_dat_clean %>% 
	filter(!is.na(CT.T1) & !is.na(N1_CT_Value)) %>%
	lm(N1_CT_Value ~ CT.T1, data=.) 
ct_reg_fit_b0 <- ct_reg_fit$coef[[1]]
ct_reg_fit_b1 <- ct_reg_fit$coef[[2]]

# Store the residuals: 
ct_resid_df <- ct_dat_clean %>% 
	filter(!is.na(CT.T1) & !is.na(N1_CT_Value)) %>%
	mutate(N1_CT_Value_resid = N1_CT_Value - (ct_reg_fit_b0 + ct_reg_fit_b1*CT.T1)) %>%
	select(CT.T1, N1_CT_Value, N1_CT_Value_resid)

# Store the sd of the residuals: 
orlando_residual_sd <- sd(ct_resid_df$N1_CT_Value_resid)

# Visualize the residuals: 
fig_yale_florida_residuals <- ct_resid_df %>%
	ggplot(aes(x=CT.T1, y=N1_CT_Value_resid)) + 
		geom_point() + 
		theme_minimal()  + 
		theme(text=element_text(size=14)) + 
		labs(x="Florida Ct (T1)", y="Yale Ct (N1) residual")

# write_csv(ct_resid_df, file="figure_data/FigS16/figs16.csv")

# QQ plot of the residuals: 
fig_yale_florida_qq <- ct_resid_df %>%
	mutate(N1_CT_Value_resid_std = (N1_CT_Value_resid - mean(N1_CT_Value_resid))/sd(N1_CT_Value_resid))  %>% 
	ggplot(aes(sample=N1_CT_Value_resid_std)) + 
		stat_qq() + 
		stat_qq_line() + 
		theme_minimal() + 
		theme(text=element_text(size=14)) + 
		labs(x="Normal(0,1)", y="Standardized Yale Ct (N1) residual")

# write_csv(ct_resid_df %>%
# 	mutate(N1_CT_Value_resid_std = (N1_CT_Value_resid - mean(N1_CT_Value_resid))/sd(N1_CT_Value_resid)), file="figure_data/FigS17/figs17.csv")
