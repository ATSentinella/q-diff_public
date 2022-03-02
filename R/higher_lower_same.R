#Runs a t-test and checks if sample 2 is "greater", "less" 
# or ~ "notdifferent" (not significantly different)
# than sample 1
higher_lower_same <- function(m1, m2, s1, s2, n1, n2) {
  
  se <- sqrt((s1^2/n1) + (s2^2/n2))
  
  t <- (m1 - m2)/se 
  
  # welch-satterthwaite df
  df <- ((s1^2/n1 + s2^2/n2)^2) / ((s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1))
  
  p.value <- 2*pt(-abs(t), df)
  
  p.value <- ifelse(is.na(p.value), 0, p.value)
  
  mean.difference <- m1 - m2
  
  if (is.na(mean.difference)) return(NA)
  
  if (p.value > 0.05) return("notdifferent")
  
  if (mean.difference == 0) return("notdifferent")
  
  if (mean.difference > 0) return("greater")
  
  if (mean.difference < 0) return("less")
  
}