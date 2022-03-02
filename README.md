q-diff Notebook
================

This is the overview of the “q-diff” project for the paper:
"Detecting steps in spatial genetic data: which diversity measures are best?"

We aim to compare the q-profile with commonly used measures e.g. GST (q
= 2) using a simplified one allele system.

-----

Core Packages used:

  - drake: project management and saves time for reanalysis
  - tidyverse: for data manipulation
  
Scripts:

[Run all analyses](R/plan.R)

We use the package ‘drake’ to make a “plan” of all analyses

It detects changes in files and scripts and only reruns what will be
different (i.e. it will skip analyses that don't need to be rerun).

  - make(plan) runs all analyses
  - readd(target) to extract the data included in any target

[q diversity functions](R/q_diversity_functions.R)

The core diversity calculation functions. More detail in paper.

One allele simulations

  - get.one.allele.data(step.k, n.samples)
    
      - Enter
          - strength of step (k) from 0 to 1 (0.5 at d = 0.5)
              - 0: linear
              - 1: gentle curve
              - 10: steep curve
              - 100: solid step
          - number of samples (n)
              - Number of sample showing th same pattern
              - For calculating standard error
      - Return
          - Table of simulated allele frequencies
              - step strength (k)
              - number of sample (n)
              - distance (d) from 0 to 1, increments of 0.01
              - alelle proportion (p) from 0 to 1
                  - calculated by: qbeta(d, 1/(1 + k), 1/(1 + k)
                  - qbeta function makes a beta distribution
                  - NOT related to q-profile or beta diversity
              - standard error for p (se.p)
              - upper esitmate of p (p.u), p + se.p (maximum value 1)
              - lower estimate of p (p,l), p - se.p (minimum value 0)
          - Plot of allele frequencies (p) over distance (d) with
            standard errors (p.u, p.l)
