#The plan - runs all of the code for this project

the_plan <- drake_plan(
    
    #Default variables for step sensitivity tests
    default_values = list(
      step = c(0, 1, 5, 50),  #4 levels of step
      n.samples = 20, #default
      p.start = c(0, 0.1, 0.3), #Different allele proportion treatments
      p.end = c(0.2, 0.5, 0.9, 1), #Different allele proportion treatments
      n.loci = 1000, #default
      n.pops = 10, #default
      n_reps = 100 #replication number
    ),
    
    #Run simulations of all combinations of default values and varied n.samples
    step_sensitivities_samples = measure_step_detection_sensitivity(
      step = default_values$step,
      n.samples = seq(2, 30, 2), #varied input
      p.start = default_values$p.start,
      p.end = default_values$p.end,
      n.loci = default_values$n.loci,
      n.pops = default_values$n.pops,
      n_reps = default_values$n_reps
    ),
    
    # Run simulations of all combinations of default values and varied n.loci
    step_sensitivities_loci = measure_step_detection_sensitivity(
      step = default_values$step,
      n.samples = default_values$n.samples,
      p.start = default_values$p.start,
      p.end = default_values$p.end,
      n.loci = seq(100, 2000, 100), #varied input
      n.pops = default_values$n.pops,
      n_reps = default_values$n_reps
    ), 
    
    # Run simulations of all combinations of default values and varied n.pops
    step_sensitivities_pops = measure_step_detection_sensitivity(
      step = default_values$step,
      n.samples = default_values$n.samples,
      p.start = default_values$p.start,
      p.end = default_values$p.end,
      n.loci = default_values$n.loci,
      n.pops = seq(4, 15, 1), #varied input
      n_reps = default_values$n_reps
    ),
    
    #Merge results into a single data frame
    merged_results = merge_results(step_sensitivities_samples,
                                   step_sensitivities_loci,
                                   step_sensitivities_pops),
    
    #Turn results into tables for publication, including all formatting
    #Standard treatment results for q = 0, q = 1,q = 2 and BC)
    results_table = get_results_table(merged_results), 
    
    #Percentage of all step sensitivities for best 6 measures
    relative_sensitivities = extract_relative_senstivities(merged_results),
    
    #Summary of measures' properties, all treatments
    properties_table = make_properties_summary_table(merged_results),
    
    #Creates supplementary figures for all tested measures
    individual_figures = generate_sensitivity_figures(merged_results),
    
    #Creates supplementary figures comparing the best 6 measures
    comparison_figures = create_comparison_figures(),
    
    alpha_peaks = vis_alpha_peaks()
    
  )

