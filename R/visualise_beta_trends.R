##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##'
##' @title

visualise_beta_trends <- function() {

  #How do each of the measures vary from 0 to 1
  #i.e. location 1 has a allele frequency of (p = 0)
  # and location 2  has a allele frequency of (p = 0 to 1)
  #So differentiation should range from 0 (no differentiation)
  #                                  to 1 (complete differentiation)
  
  p.start = 0.1
  p.end = 0.9
  n.samples = 10
  n.loci = 1000
  #Make a variable which is the distance between each value of d (distance)
  increment.size = 0.01
  
  
  #Create a table with 
  data <- tibble(
    n.samples = n.samples,
    i = increment.size,
    p.start = p.start,
    p.end = p.end,
    d = seq(0, 1, increment.size)) %>% #Location of population along distance (0 to 1)
    mutate(p = p.start + d*(p.end - p.start))
  
  #p.binom:
  #Take n.samples from a binomial distribution around allele proportion p
  #Divide by n.samples to get 'measured' allele frequency
  #e.g. n.samples = 5, allele proportion p = 0.1
  # rbinom(1,5,0.1) / 5
  # Can only be 0, 0.2, 0.4, 0.6, 0.8, 1
  # But would more likely be 0/0.2  
  # repeat over n.loci
  data <- data %>% 
    .[rep(1:nrow(.), times = n.loci),] %>% #replicate over multiple loci with same p
    rowwise() %>% #allows for mutate to work row by row (rather than as a vector)
    mutate(p.binom = rbinom(1, n.samples, p)/n.samples) %>% #p.binom
    ungroup() #stops rowwise operations

  data <- data %>%
    
    mutate(H0b = get.Hq.beta.per.locus(p.start, p.binom, 0),
           H1b =  get.Hq.beta.per.locus(p.start, p.binom, 1),
           H2b.JOST = get.Hq.beta.per.locus(p.start, p.binom, 2),
           H2b.GST = get.Hq.beta.per.locus(p.start, p.binom, 2, "Gst"),
           D0b = get.Dq.beta.per.locus(p.start, p.binom, 0),
           D1b = get.Dq.beta.per.locus(p.start, p.binom, 1),
           D2b = get.Dq.beta.per.locus(p.start, p.binom, 2),
           AFD = get.AFD(p.start, p.binom)
    )
  
  data_summary <- data %>%
    group_by(d) %>%
    summarise_each(list(mean = mean, sd = sd)) %>%
    as_tibble()%>%
    mutate(H0b = get.Hq.beta.per.locus(p.start_mean, p_mean, 0),
           H1b =  get.Hq.beta.per.locus(p.start_mean, p_mean, 1),
           H2b.JOST = get.Hq.beta.per.locus(p.start_mean, p_mean, 2),
           H2b.GST = get.Hq.beta.per.locus(p.start_mean, p_mean, 2, "Gst"),
           D0b = get.Dq.beta.per.locus(p.start_mean, p_mean, 0),
           D1b = get.Dq.beta.per.locus(p.start_mean, p_mean, 1),
           D2b = get.Dq.beta.per.locus(p.start_mean, p_mean, 2),
           AFD = get.AFD(p.start_mean, p_mean)
    )
  
  ggplot(data_summary, aes(x = d))+
    geom_point(shape = 3, aes(y = H2b.GST_mean), colour = "darkgreen")+
    geom_point(shape = 3, aes(y = H1b_mean), colour = "lightblue")+
    geom_point(shape = 3, aes(y = H0b_mean), colour = "red")+
    geom_point(shape = 3, aes(y = H2b.JOST_mean), colour = "purple")+
    geom_point(shape = 3, aes(y = AFD_mean), colour = "orange")+
    geom_point(shape = 3, aes(y = D1b_mean-1), colour = "blue")+
    geom_point(aes(y = D1b-1), colour = "blue")+
    geom_point(aes(y = H2b.GST), colour = "darkgreen")+
    geom_point(aes(y = H1b), colour = "lightblue")+
    geom_point(aes(y = H0b), colour = "red")+
    geom_point(aes(y = H2b.JOST), colour = "purple")+
    geom_point(aes(y = AFD), colour = "orange") 
  
  ggplot(data_summary)+
    geom_point(aes(x = H2b.GST, y = H2b.GST_mean), colour = "darkgreen")+
    geom_point(aes(x = D1b-1, y = D1b_mean-1), colour = "green")+
    geom_point(aes(x = D1b-1, y = D1b_mean-1), colour = "blue")+
    geom_point(aes(x = H0b, y = H0b_mean), colour = "pink")+
    geom_point(aes(x = D0b-1, y = D0b_mean-1), colour = "red")+
    geom_point(aes(x = H1b, y = H1b_mean), colour = "lightblue")+
    geom_point(aes(x = H2b.JOST, y = H2b.JOST_mean), colour = "purple")+
    geom_point(aes(x = AFD, y = AFD_mean), colour = "orange")+
    geom_abline(slope=1)+
    xlab("True differentiation (0 to 1)")+
    ylab("Measured differentiation (0 to 1)")
  
  
  
}
