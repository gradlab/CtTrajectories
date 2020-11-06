m_conv <- -3.609714286
b_conv <- 40.93733333


convert_Ct_logGEML <- function(Ct, m_conv=-3.609714286, b_conv=40.93733333){
	out <- (Ct-b_conv)/m_conv * log10(10) + log10(250)
	return(out) 
}

10^convert_Ct_logGEML(34.55)

ct_dat_clean_new <- ct_dat_clean %>% 
	mutate(log10_GEperML = (CT.Mean-b_conv)/m_conv * log10(10) + log10(250)) 

write.csv(ct_dat_clean_new, file="data/ct_dat_clean_new.csv", row.names=FALSE) 

	%>%
	select(Person.ID, CT.Mean, log10_GEperML, GEperML) %>% 
	print(n=100) 
