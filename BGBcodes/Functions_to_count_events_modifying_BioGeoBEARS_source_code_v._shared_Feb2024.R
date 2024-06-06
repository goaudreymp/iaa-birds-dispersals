# New version of event counter function. February 2024
# Editing some source code in BioGeoBEARS functions
library(cladoRcpp); library(BioGeoBEARS); library(reshape2); library(picante); library(betapart); library(ape); library(phytools)


count_ana_clado_eventsM <- function (clado_events_tables, ana_events_tables, areanames, 
                                     actual_names, timeperiod = NULL, verbose=TRUE) 
{
  defaults = "\n\tareanames = c(\"A\", \"B\", \"C\", \"D\", \"E\", \"F\", \"G\")\n\tactual_names = c(\"SAm\", \"CAm\", \"Car\", \"NAm\", \"AF\", \"EU\", \"OZ\")\n\tclado_events_table = clado_events_tables[[1]]\n\t"
  TF = nchar(areanames) > 1
  if (any(TF)) {
    txt = paste0("STOP ERROR: all 'areanames' must be only 1 character long for the counting of events in stochastic mapping to work. You have violations for these 'areanames': ", 
                 paste(areanames[TF], sep = ", "), ".")
    cat("\n\n")
    cat(txt)
    cat("\n\n")
    stop(txt)
  }
  if (length(areanames) != length(actual_names)) {
    txt = paste0("STOP ERROR in count_ana_clado_events(): areanames and actual_names must have the same length. Instead, length(areanames)=", 
                 length(areanames), ", length(actual_names)=", length(actual_names), 
                 ".")
    cat("\n\n")
    cat(txt)
    cat("\n\n")
    stop(txt)
  }
  if (is.null(timeperiod) == FALSE) {
    if (length(timeperiod) != 2) {
      errortxt = paste("\n\nERROR in count_events_huge_tables(): timeperiod must either be NULL or be 2 numbers (min and max of time bin, in time before present).\n\n", 
                       sep = "")
      cat(errortxt)
      stop(errortxt)
    }
    minT = min(timeperiod)
    maxT = max(timeperiod)
    original_clado_events_tables = clado_events_tables
    original_ana_events_tables = ana_events_tables
    for (i in 1:length(original_clado_events_tables)) { # This loop trims the clado and anagenetic events tables by the timeperiod selected
      clado_events_table = clado_events_tables[[i]]
      ana_events_table = ana_events_tables[[i]]
      TF1 = clado_events_table$time_bp < maxT
      TF2 = clado_events_table$time_bp >= minT
      TF = (TF1 + TF2) == 2
      clado_events_table = clado_events_table[TF, ]
      if(is.data.frame(ana_events_table)){ # Conditional added by Octavio to avoid error:
        # Error in ana_events_table$abs_event_time : 
        #   $ operator is invalid for atomic vectors
        TF1 = ana_events_table$abs_event_time < maxT
        TF2 = ana_events_table$abs_event_time >= minT
        TF = (TF1 + TF2) == 2
        ana_events_table = ana_events_table[TF, ]
      }
      clado_events_tables[[i]] = clado_events_table
      ana_events_tables[[i]] = ana_events_table
    }
    sapply(clado_events_tables, nrow)
    sapply(ana_events_tables, nrow)
    sapply(original_clado_events_tables, nrow)
    sapply(original_ana_events_tables, nrow)
  }
  else {
    minT = 0
    maxT = 1e+50
  }
  dimdata = c(length(areanames), length(areanames))
  dims = c(dimdata, length(clado_events_tables)) # dimensions to define a cubic matrix with dimensions of number of areas x areas x maps. This will be used to define the cubic matrices ending in "_counts_cube"
  founder_counts_cube = array(NA, dim = dims) 
  founder_totals_list = rep(NA, length(clado_events_tables))
  dimdata = c(length(areanames), length(areanames))
  dims = c(dimdata, length(ana_events_tables))
  ana_counts_cube = array(NA, dim = dims)
  ana_totals_list = rep(NA, length(ana_events_tables))
  a_counts_cube = array(NA, dim = dims)
  a_totals_list = rep(NA, length(ana_events_tables))
  d_counts_cube = array(NA, dim = dims)
  d_totals_list = rep(NA, length(ana_events_tables))
  e_counts_rectangle = matrix(NA, nrow = length(ana_events_tables), 
                              ncol = length(areanames))
  e_totals_list = rep(NA, length(ana_events_tables))
  vicariance_clado_events_table_bigTable = NULL
  subsetSymp_clado_events_table_bigTable = NULL
  sympatry_clado_events_table_bigTable = NULL
  vicariance_totals_list = rep(NA, length(clado_events_tables))
  subsetSymp_totals_list = rep(NA, length(clado_events_tables))
  sympatry_totals_list = rep(NA, length(clado_events_tables))
  for (i in 1:length(clado_events_tables)) { 
    clado_events_table = clado_events_tables[[i]]
    if(nrow(clado_events_table[which(clado_events_table$clado_event_type!=""),])>1){ # conditional added here by Octavio
      clado_events_table = uniquify_clado_events(clado_events_table)
      # This gives the following error when there are less than two cladogenetic events in a particular time slice
      # Tips appear in the cladogenetic events table with clado_event_type column as "",
      # Those rows are later trimmed at the start of uniquify_clado_events, creating the same problem
      # It gives the following error:
      # Error in names(nchar_LRdf) <- c("L", "R") : 
      #   'names' attribute [2] must be the same length as the vector [1]
    }  # else { if(verbose) print(paste("WARNING!! Map", i, "has", nrow(clado_events_table), "cladogenetic event(s)") ) } # added by Octavio
    
    vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == 
                                                         "vicariance (v)", ]
    vicariance_clado_events_table
    vic_table = table(vicariance_clado_events_table$clado_event_txt)
    names_vic_table = names(vic_table)
    vic_table
    if (length(vic_table) > 0) {
      vic_table2 = cbind(i, names(vic_table), vic_table)
    }
    else {
      vic_table2 = NULL
    }
    vicariance_clado_events_table_bigTable = rbind(vicariance_clado_events_table_bigTable, 
                                                   vic_table2)
    subsetSymp_clado_events_table = clado_events_table[clado_events_table$clado_event_type == 
                                                         "subset (s)", ]
    subsetSymp_clado_events_table
    sub_table = table(subsetSymp_clado_events_table$clado_event_txt)
    names_sub_table = names(sub_table)
    sub_table
    if (length(sub_table) > 0) {
      sub_table2 = cbind(i, names(sub_table), sub_table)
    }
    else {
      sub_table2 = NULL
    }
    subsetSymp_clado_events_table_bigTable = rbind(subsetSymp_clado_events_table_bigTable, 
                                                   sub_table2)
    sympatry_clado_events_table = clado_events_table[clado_events_table$clado_event_type == 
                                                       "sympatry (y)", ]
    sympatry_clado_events_table
    sym_table = table(sympatry_clado_events_table$clado_event_txt)
    names_sym_table = names(sym_table)
    sym_table
    if (length(sym_table) > 0) {
      sym_table2 = cbind(i, names(sym_table), sym_table)
    }
    else {
      sym_table2 = NULL
    }
    sympatry_clado_events_table_bigTable = rbind(sympatry_clado_events_table_bigTable, 
                                                 sym_table2)
  }
  unique_vic_events = NULL
  unique_sub_events = NULL
  unique_sym_events = NULL
  if (is.null(vicariance_clado_events_table_bigTable) == FALSE) {
    unique_vic_events = sort(unique(vicariance_clado_events_table_bigTable[, 
                                                                           2]))
    unique_vic_events
  }
  if (is.null(subsetSymp_clado_events_table_bigTable) == FALSE) {
    unique_sub_events = sort(unique(subsetSymp_clado_events_table_bigTable[, 
                                                                           2]))
    unique_sub_events
  }
  if (is.null(sympatry_clado_events_table_bigTable) == FALSE) {
    unique_sym_events = sort(unique(sympatry_clado_events_table_bigTable[, 
                                                                         2]))
    unique_sym_events
  }
  unique_vic_counts = matrix(NA, nrow = length(clado_events_tables), 
                             ncol = length(unique_vic_events))
  unique_sub_counts = matrix(NA, nrow = length(clado_events_tables), 
                             ncol = length(unique_sub_events))
  unique_sym_counts = matrix(NA, nrow = length(clado_events_tables), 
                             ncol = length(unique_sym_events))
  # cat("\ncount_ana_clado_events() is counting events over Biogeographical Stochastic Map (BSM) #:\n")
  for (i in 1:length(clado_events_tables)) { # when there is only one cladogenetic event per map, the original function was also stuck in this loop after fixing the previous one. Now with the conditional below it runs
    # This loop goes map by map, filling the _counts_cube cubic matrices with the 2x2 matrices corresponding to area to area events for each map i.
    # cat(i, " ", sep = "")
    clado_events_table = clado_events_tables[[i]]
    ana_events_table = ana_events_tables[[i]]
    founder_counts_df = count_clado_dispersal_events(clado_events_table, 
                                                     areanames = areanames, actual_names = actual_names)
    # count_clado_dispersal_events() creates a area x area matrix with the cladogenetic dispersal events (founder events j)
    # it goes area by area looking at the columns clado_dispersals_to and clado_dispersals_from in the clado_events_table 
    # and with that it reconstructs the events in the area(from) x area(to) matrix
    # The column clado_dispersals_from is created by simulate_source_areas_clado()
    # In the case of ancestral area with 2 or more areas, one of the areas is chosen using a dispersal matrix dmat 
    # This dispersal matrix dmat (and its sliced times if it has) is gotten from the state_list object within the input of the BioGeoBEARS result object (lines 3-5 in simulate_source_areas_ana_clado) 
    
    
    
    if (is.null(founder_counts_df) == TRUE) {
      founder_counts_df = matrix(data = 0, nrow = length(areanames), 
                                 ncol = length(areanames))
    }
    founder_counts_cube[, , i] = as.matrix(founder_counts_df)
    founder_totals_list[i] = sum(founder_counts_df)
    if ((length(ana_events_table) > 1) || (is.na(ana_events_table) == 
                                           FALSE)) {
      ana_disp_counts = count_ana_dispersal_events(ana_events_table, 
                                                   areanames = areanames, actual_names = actual_names)
      # count_ana_dispersal_events() create 4 area x area matrices with the anagenetic dispersal events:
      e_counts_df = ana_disp_counts$e_counts_df # e, range contraction
      a_counts_df = ana_disp_counts$a_counts_df # a, range switching
      d_counts_df = ana_disp_counts$d_counts_df # d, range expansion
      ana_counts_df = ana_disp_counts$counts_df # 
      
      # count_ana_dispersal_events() goes area by area looking at the columns dispersals_to and ana_dispersals_from in the ana_events_table 
      # and with that it reconstructs the events in the area(from) x area(to) matrix
      # The column ana_dispersals_from is created by simulate_source_areas_ana() [lines 53-82]
      
    }
    if ((length(ana_events_table) > 1) || (is.na(ana_events_table) == 
                                           FALSE)) {
      if (!is.null(ana_counts_df)) {
        ana_counts_cube[, , i] = as.matrix(ana_counts_df)
        ana_totals_list[i] = sum(ana_counts_df)
      }
      else {
        ana_counts_cube[, , i] = 0
        ana_totals_list[i] = 0
      }
      if (!is.null(a_counts_df)) {
        a_counts_cube[, , i] = as.matrix(a_counts_df)
        a_totals_list[i] = sum(a_counts_df)
      }
      else {
        a_counts_cube[, , i] = 0
        a_totals_list[i] = 0
      }
      if (!is.null(d_counts_df)) {
        d_counts_cube[, , i] = as.matrix(d_counts_df)
        d_totals_list[i] = sum(d_counts_df)
      }
      else {
        d_counts_cube[, , i] = 0
        d_totals_list[i] = 0
      }
      if (!is.null(e_counts_df)) {
        e_counts_rectangle[i, ] = as.matrix(e_counts_df)
        e_totals_list[i] = sum(e_counts_df)
      }
      else {
        e_counts_rectangle[i, ] = 0
        e_totals_list[i] = 0
      }
    }
    if(nrow(clado_events_table[which(clado_events_table$clado_event_type!=""),])>1){ # conditional added by Octavio to avoid error
      clado_events_table = uniquify_clado_events(clado_events_table)
    } 
    
    vicariance_clado_events_table = clado_events_table[clado_events_table$clado_event_type == 
                                                         "vicariance (v)", ]
    vicariance_clado_events_table
    if (nrow(vicariance_clado_events_table) > 0) {
      vic_table = table(vicariance_clado_events_table$clado_event_txt)
      names_vic_table = names(vic_table)
      event_in_uniqs_TF = unique_vic_events %in% names_vic_table
      unique_vic_counts[i, event_in_uniqs_TF] = vic_table
      vicariance_totals_list[i] = sum(unique_vic_counts[i, 
                                                        event_in_uniqs_TF])
    }
    else {
      vic_table = NA
      vicariance_totals_list[i] = 0
    }
    subsetSymp_clado_events_table = clado_events_table[clado_events_table$clado_event_type == 
                                                         "subset (s)", ]
    subsetSymp_clado_events_table
    if (nrow(subsetSymp_clado_events_table) > 0) {
      sub_table = table(subsetSymp_clado_events_table$clado_event_txt)
      names_sub_table = names(sub_table)
      event_in_uniqs_TF = unique_sub_events %in% names_sub_table
      unique_sub_counts[i, event_in_uniqs_TF] = sub_table
      subsetSymp_totals_list[i] = sum(unique_sub_counts[i, 
                                                        event_in_uniqs_TF])
    }
    else {
      sub_table = NA
      subsetSymp_totals_list[i] = 0
    }
    sub_table
    sympatry_clado_events_table = clado_events_table[clado_events_table$clado_event_type == 
                                                       "sympatry (y)", ]
    sympatry_clado_events_table
    if (nrow(sympatry_clado_events_table) > 0) {
      sym_table = table(sympatry_clado_events_table$clado_event_txt)
      names_sym_table = names(sym_table)
      event_in_uniqs_TF = unique_sym_events %in% names_sym_table
      unique_sym_counts[i, event_in_uniqs_TF] = sym_table
      sympatry_totals_list[i] = sum(unique_sym_counts[i, 
                                                      event_in_uniqs_TF])
    }
    else {
      sym_table = NA
      sympatry_totals_list[i] = 0
    }
  }
  unique_vic_counts = adf2(unique_vic_counts)
  names(unique_vic_counts) = unique_vic_events
  if (ncol(unique_vic_counts) > 0) {
    unique_vic_counts[is.na(unique_vic_counts)] = 0
  }
  unique_sub_counts = adf2(unique_sub_counts)
  names(unique_sub_counts) = unique_sub_events
  if (ncol(unique_sub_counts) > 0) {
    unique_sub_counts[is.na(unique_sub_counts)] = 0
  }
  unique_sym_counts = adf2(unique_sym_counts)
  names(unique_sym_counts) = unique_sym_events
  if (ncol(unique_sym_counts) > 0) {
    unique_sym_counts[is.na(unique_sym_counts)] = 0
  }
  head(unique_vic_counts)
  head(unique_sub_counts)
  head(unique_sym_counts)
  round(apply(X = unique_vic_counts, MARGIN = 2, FUN = mean), 
        1)
  round(apply(X = unique_vic_counts, MARGIN = 2, FUN = sd), 
        1)
  round(apply(X = unique_sub_counts, MARGIN = 2, FUN = mean), 
        1)
  round(apply(X = unique_sub_counts, MARGIN = 2, FUN = sd), 
        1)
  round(apply(X = unique_sym_counts, MARGIN = 2, FUN = mean), 
        1)
  round(apply(X = unique_sym_counts, MARGIN = 2, FUN = sd), 
        1)
  vic_counts_means = apply(X = unique_vic_counts, MARGIN = 2, 
                           FUN = mean)
  vic_counts_sds = apply(X = unique_vic_counts, MARGIN = 2, 
                         FUN = sd)
  vic_counts_sums = apply(X = unique_vic_counts, MARGIN = 2, 
                          FUN = sum)
  sub_counts_means = apply(X = unique_sub_counts, MARGIN = 2, 
                           FUN = mean)
  sub_counts_sds = apply(X = unique_sub_counts, MARGIN = 2, 
                         FUN = sd)
  sub_counts_sums = apply(X = unique_sub_counts, MARGIN = 2, 
                          FUN = sum)
  sym_counts_means = apply(X = unique_sym_counts, MARGIN = 2, 
                           FUN = mean)
  sym_counts_sds = apply(X = unique_sym_counts, MARGIN = 2, 
                         FUN = sd)
  sym_counts_sums = apply(X = unique_sym_counts, MARGIN = 2, 
                          FUN = sum)
  sum(sym_counts_means)
  sum(sub_counts_means)
  sum(vic_counts_means)
  counts_list = NULL
  founder_counts_cube[is.na(founder_counts_cube)] = 0
  ana_counts_cube[is.na(ana_counts_cube)] = 0
  a_counts_cube[is.na(a_counts_cube)] = 0
  d_counts_cube[is.na(d_counts_cube)] = 0
  e_counts_rectangle[is.na(e_counts_rectangle)] = 0
  unique_sub_counts[is.na(unique_sub_counts)] = 0
  unique_vic_counts[is.na(unique_vic_counts)] = 0
  unique_sym_counts[is.na(unique_sym_counts)] = 0
  founder_totals_list[is.na(founder_totals_list)] = 0
  ana_totals_list[is.na(ana_totals_list)] = 0
  a_totals_list[is.na(a_totals_list)] = 0
  d_totals_list[is.na(d_totals_list)] = 0
  e_totals_list[is.na(e_totals_list)] = 0
  subsetSymp_totals_list[is.na(subsetSymp_totals_list)] = 0
  vicariance_totals_list[is.na(vicariance_totals_list)] = 0
  sympatry_totals_list[is.na(sympatry_totals_list)] = 0
  all_dispersals_counts_cube = founder_counts_cube
  all_dispersals_totals_list = founder_totals_list
  anagenetic_dispersals_counts_cube = d_counts_cube
  anagenetic_dispersals_totals_list = d_totals_list
  for (i in 1:length(all_dispersals_totals_list)) {
    all_dispersals_counts_cube[, , i] = all_dispersals_counts_cube[, 
                                                                   , i] + a_counts_cube[, , i] + d_counts_cube[, , 
                                                                                                               i]
    all_dispersals_totals_list[i] = all_dispersals_totals_list[i] + 
      a_totals_list[i] + d_totals_list[i]
    anagenetic_dispersals_counts_cube[, , i] = anagenetic_dispersals_counts_cube[, 
                                                                                 , i] + a_counts_cube[, , i]
    anagenetic_dispersals_totals_list[i] = anagenetic_dispersals_totals_list[i] + 
      a_totals_list[i]
  }
  counts_list$all_dispersals_counts_cube = all_dispersals_counts_cube
  counts_list$anagenetic_dispersals_counts_cube = anagenetic_dispersals_counts_cube
  counts_list$founder_counts_cube = founder_counts_cube
  counts_list$ana_counts_cube = ana_counts_cube
  counts_list$a_counts_cube = a_counts_cube
  counts_list$d_counts_cube = d_counts_cube
  counts_list$e_counts_rectangle = e_counts_rectangle
  counts_list$unique_sub_counts = unique_sub_counts
  counts_list$unique_vic_counts = unique_vic_counts
  counts_list$unique_sym_counts = unique_sym_counts
  counts_list$all_dispersals_totals_list = all_dispersals_totals_list
  counts_list$anagenetic_dispersals_totals_list = anagenetic_dispersals_totals_list
  counts_list$founder_totals_list = founder_totals_list
  counts_list$ana_totals_list = ana_totals_list
  counts_list$a_totals_list = a_totals_list
  counts_list$d_totals_list = d_totals_list
  counts_list$e_totals_list = e_totals_list
  counts_list$subsetSymp_totals_list = subsetSymp_totals_list
  counts_list$vicariance_totals_list = vicariance_totals_list
  counts_list$sympatry_totals_list = sympatry_totals_list
  clado_totals_list = founder_totals_list + subsetSymp_totals_list + 
    vicariance_totals_list + sympatry_totals_list
  all_totals_list = founder_totals_list + subsetSymp_totals_list + 
    vicariance_totals_list + sympatry_totals_list + ana_totals_list
  counts_list$clado_totals_list = clado_totals_list
  counts_list$all_totals_list = all_totals_list
  clado_totals_list
  all_totals_list
  means = c(mean(founder_totals_list), mean(a_totals_list), 
            mean(d_totals_list), mean(e_totals_list), mean(subsetSymp_totals_list), 
            mean(vicariance_totals_list), mean(sympatry_totals_list), 
            mean(all_dispersals_totals_list), mean(anagenetic_dispersals_totals_list), 
            mean(ana_totals_list), mean(clado_totals_list), mean(all_totals_list))
  sds = c(sd(founder_totals_list), sd(a_totals_list), sd(d_totals_list), 
          sd(e_totals_list), sd(subsetSymp_totals_list), sd(vicariance_totals_list), 
          sd(sympatry_totals_list), sd(all_dispersals_totals_list), 
          sd(anagenetic_dispersals_totals_list), sd(ana_totals_list), 
          sd(clado_totals_list), sd(all_totals_list))
  sums = c(sum(founder_totals_list), sum(a_totals_list), sum(d_totals_list), 
           sum(e_totals_list), sum(subsetSymp_totals_list), sum(vicariance_totals_list), 
           sum(sympatry_totals_list), sum(all_dispersals_totals_list), 
           sum(anagenetic_dispersals_totals_list), sum(ana_totals_list), 
           sum(clado_totals_list), sum(all_totals_list))
  all_dispersals_counts_fromto_means = adf2(apply(X = all_dispersals_counts_cube, 
                                                  MARGIN = c(1, 2), FUN = mean))
  ana_dispersals_counts_fromto_means = adf2(apply(X = anagenetic_dispersals_counts_cube, 
                                                  MARGIN = c(1, 2), FUN = mean))
  founder_counts_fromto_means = adf2(apply(X = founder_counts_cube, 
                                           MARGIN = c(1, 2), FUN = mean))
  ana_events_counts_fromto_means = adf2(apply(X = ana_counts_cube, 
                                              MARGIN = c(1, 2), FUN = mean))
  a_counts_fromto_means = adf2(apply(X = a_counts_cube, MARGIN = c(1, 
                                                                   2), FUN = mean))
  d_counts_fromto_means = adf2(apply(X = d_counts_cube, MARGIN = c(1, 
                                                                   2), FUN = mean))
  all_dispersals_counts_fromto_sds = adf2(apply(X = all_dispersals_counts_cube, 
                                                MARGIN = c(1, 2), FUN = sd))
  ana_dispersals_counts_fromto_sds = adf2(apply(X = anagenetic_dispersals_counts_cube, 
                                                MARGIN = c(1, 2), FUN = sd))
  founder_counts_fromto_sds = adf2(apply(X = founder_counts_cube, 
                                         MARGIN = c(1, 2), FUN = sd))
  ana_counts_fromto_sds = adf2(apply(X = ana_counts_cube, 
                                     MARGIN = c(1, 2), FUN = sd))
  a_counts_fromto_sds = adf2(apply(X = a_counts_cube, MARGIN = c(1, 
                                                                 2), FUN = sd))
  d_counts_fromto_sds = adf2(apply(X = d_counts_cube, MARGIN = c(1, 
                                                                 2), FUN = sd))
  names(all_dispersals_counts_fromto_means) = actual_names
  row.names(all_dispersals_counts_fromto_means) = actual_names
  names(ana_dispersals_counts_fromto_means) = actual_names
  row.names(ana_dispersals_counts_fromto_means) = actual_names
  names(founder_counts_fromto_means) = actual_names
  row.names(founder_counts_fromto_means) = actual_names
  names(a_counts_fromto_means) = actual_names
  row.names(a_counts_fromto_means) = actual_names
  names(d_counts_fromto_means) = actual_names
  row.names(d_counts_fromto_means) = actual_names
  names(all_dispersals_counts_fromto_sds) = actual_names
  row.names(all_dispersals_counts_fromto_sds) = actual_names
  names(ana_dispersals_counts_fromto_sds) = actual_names
  row.names(ana_dispersals_counts_fromto_sds) = actual_names
  names(founder_counts_fromto_sds) = actual_names
  row.names(founder_counts_fromto_sds) = actual_names
  names(a_counts_fromto_sds) = actual_names
  row.names(a_counts_fromto_sds) = actual_names
  names(d_counts_fromto_sds) = actual_names
  row.names(d_counts_fromto_sds) = actual_names
  counts_list$all_dispersals_counts_fromto_means = all_dispersals_counts_fromto_means
  counts_list$ana_dispersals_counts_fromto_means = ana_dispersals_counts_fromto_means
  counts_list$founder_counts_fromto_means = founder_counts_fromto_means
  counts_list$a_counts_fromto_means = a_counts_fromto_means
  counts_list$d_counts_fromto_means = d_counts_fromto_means
  counts_list$all_dispersals_counts_fromto_sds = all_dispersals_counts_fromto_sds
  counts_list$ana_dispersals_counts_fromto_sds = ana_dispersals_counts_fromto_sds
  counts_list$founder_counts_fromto_sds = founder_counts_fromto_sds
  counts_list$a_counts_fromto_sds = a_counts_fromto_sds
  counts_list$d_counts_fromto_sds = d_counts_fromto_sds
  {
    # cat("\n\nRange-switching dispersal (all observed 'a' dispersals):\n")
    # cat("\nmeans:\n")
    # print(conditional_format_table(counts_list$a_counts_fromto_means))
    # cat("\nstandard deviations:\n")
    # print(conditional_format_table(counts_list$a_counts_fromto_sds))
    # cat("\n\nRange-expansion dispersal (all observed 'd' dispersals):\n")
    # cat("\nmeans:\n")
    # print(conditional_format_table(counts_list$d_counts_fromto_means))
    # cat("\nstandard deviations:\n")
    # print(conditional_format_table(counts_list$d_counts_fromto_sds))
    # cat("\n\nAnagenetic dispersal (mean of all observed anagenetic 'a' or 'd' dispersals):\n")
    # cat("\nmeans:\n")
    # print(conditional_format_table(counts_list$ana_dispersals_counts_fromto_means))
    # cat("\nstandard deviations:\n")
    # print(conditional_format_table(counts_list$ana_dispersals_counts_fromto_sds))
    # cat("\n\nCladogenetic dispersal (mean of all observed jump 'j' dispersals):\n")
    # cat("\nmeans:\n")
    # print(conditional_format_table(counts_list$founder_counts_fromto_means))
    # cat("\nstandard deviations:\n")
    # print(conditional_format_table(counts_list$founder_counts_fromto_sds))
    # cat("\n\nALL dispersal (mean of all observed anagenetic 'a', 'd' dispersals, PLUS cladogenetic founder/jump dispersal):\n")
    # cat("\nmeans:\n")
    # print(conditional_format_table(counts_list$all_dispersals_counts_fromto_means))
    # cat("\nstandard deviations:\n")
    # print(conditional_format_table(counts_list$all_dispersals_counts_fromto_sds))
  } # chunk of printing code
  summary_counts_BSMs = rbind(means, sds, sums)
  summary_counts_BSMs = adf(summary_counts_BSMs)
  names(summary_counts_BSMs) = c("founder", "a", "d", "e", 
                                 "subset", "vicariance", "sympatry", "ALL_disp", "ana_disp", 
                                 "all_ana", "all_clado", "total_events")
  row.names(summary_counts_BSMs) = c("means", "stdevs", "sums")
  counts_list$summary_counts_BSMs = summary_counts_BSMs
  if(verbose) { cat("\nSummary of event counts from ", length(clado_events_tables), 
                    " BSMs:\n", sep = "")
    print(conditional_format_table(summary_counts_BSMs)) }
  
  extract = "\n\tall_dispersals_counts_cube = counts_list$all_dispersals_counts_cube\n\tanagenetic_dispersals_counts_cube = counts_list$anagenetic_dispersals_counts_cube\n\tfounder_counts_cube = counts_list$founder_counts_cube\n\tana_counts_cube = counts_list$ana_counts_cube\n\ta_counts_cube = counts_list$a_counts_cube\n\td_counts_cube = counts_list$d_counts_cube\n\te_counts_rectangle = counts_list$e_counts_rectangle\n\tunique_sub_counts = counts_list$unique_sub_counts\n\tunique_vic_counts = counts_list$unique_vic_counts\n\tunique_sym_counts = counts_list$unique_sym_counts\n\n\tall_dispersals_totals_list = counts_list$all_dispersals_totals_list\n\tanagenetic_dispersals_totals_list = counts_list$anagenetic_dispersals_totals_list\n\tfounder_totals_list = counts_list$founder_totals_list\n\tana_totals_list = counts_list$ana_totals_list\n\ta_totals_list = counts_list$a_totals_list\n\td_totals_list = counts_list$d_totals_list\n\te_totals_list = counts_list$e_totals_list\n\tsubsetSymp_totals_list = counts_list$subsetSymp_totals_list\n\tvicariance_totals_list = counts_list$vicariance_totals_list\n\tsympatry_totals_list = counts_list$sympatry_totals_list\n\t\n\n\tsummary_counts_BSMs = counts_list$summary_counts_BSMs\n\t"
  return(counts_list)
}

EventCounter <- function(BGBresults, BSM, statenames=NULL,
                         time.slice.size=10, oldest.time.slice=100, verbose=TRUE, 
                         results.per.map=FALSE){ # the results.per.map take a lot of memory space, so I ended up putting it as optional, as it is something that I would normally not want to have
  
  # 0. Preparation of objects before t loop 
  { # FOLD THE CODE HERE
    if(file.exists(BGBresults$inputs$geogfn)) { areas <- colnames(read_PHYLIP_data(BGBresults$inputs$geogfn) )
    }  else {
      if(!is.null(statenames)){
        areas <- unlist(statenames)[which(nchar(statenames)==1)]
        areas <- areas[- which(areas == "_")]
        print("Retrieving area names from statenames object, assuming that each area is represented by a single caracter. If this is not the case, the function may crash or give wrong output")
      } else {
        stop(paste0("The file ", BGBresults$inputs$geogfn, " could not be found and there is no statenames object provided, so the areas cannot be named"))
      }
    }
    
    if(is.na(BGBresults[["inputs"]][["max_range_size"]])) {
      BGBresults[["inputs"]][["max_range_size"]] <- length(areas)
    }
    
    if(!"states_list" %in% names(BGBresults$inputs)) { 
      print("The BGBresults object does not have states_list. This is esential for later assigning states to each branch")
      print(" Therefore, the analysis continues assuming that all possible states (combinations of areas, depending on the BGBresults$inputs$max_range_size) are allowed. If it is not the case, please provide a states_list in BGBresults$inputs")
      BGBresults$inputs$states_list <- rcpp_areas_list_to_states_list(areas=areas, 
                                                                      maxareas=BGBresults$inputs$max_range_size, 
                                                                      include_null_range=BGBresults$inputs$include_null_range)
    }
    names_of_states <- sapply(BGBresults$inputs$states_list, FUN = function(x) {y <- paste(areas[x+1], collapse = ""); return(y)}) ### CHECK THE REDUNDANCY OF names_of_states WITH statenames
    
    if("clado_dispersal_from" %in% names(BSM[[1]][[1]]) & 
       "ana_dispersal_from" %in% names(BSM[[2]][[1]])  ) {
      BSMs_w_sourceAreas <- BSM
      if(verbose) print( "The stochastic maps provided had already the simulated source areas." )
    } else {
      BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res = BGBresults, clado_events_tables = BSM$RES_clado_events_tables,
                                                           ana_events_tables = BSM$RES_ana_events_tables, areanames = areas)  
    }
    
    time.slices <- seq(0, oldest.time.slice+time.slice.size, by=time.slice.size) # we add one time slice more at the end to calculate the diversities at the beginning of the last required period
    rootage <- max(nodeHeights(read.tree(BGBresults$inputs$trfn)))
    if(rootage < oldest.time.slice){ 
      print(paste("Age of tree (", round(rootage,1), "MY) is more modern than the specified oldest time slice (", oldest.time.slice-time.slice.size, "—", oldest.time.slice,
                  "MY). Resetting oldest.time.slice to", max(time.slices[-which(time.slices > (rootage+time.slice.size) )]), "MY" ))
      time.slices <- time.slices[-which(time.slices > (rootage+time.slice.size) )]
      oldest.time.slice <- time.slices[length(time.slices)]
    }
    
    # section the tree, to count branch lengths later
    BGB_run_object <- define_BioGeoBEARS_run(trfn = BGBresults$inputs$trfn)
    BGB_run_object <- readfiles_BioGeoBEARS_run(BGB_run_object)
    if(rootage > time.slices[length(time.slices)]) {  BGB_run_object$timeperiods <- c(time.slices[-1], oldest.time.slice+time.slice.size)
    } else { BGB_run_object$timeperiods <- time.slices[-1] }
    invisible(capture.output( #https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message
      tree_sectioning <- section_the_tree(BGB_run_object, make_master_table = T, plot_pieces = F, cut_fossils = F)   ))
    tree_sections_list <- tree_sectioning$tree_sections_list
    master_table <- tree_sectioning$master_table
    rm(tree_sectioning, BGB_run_object)
    
    time.slices.nh <- rootage - time.slices
    # In time-sliced BioGeoBEARS clado_events_tables are missing the node_ht column. Therefore we need to create it.
    if("node_ht" %in% names(BSMs_w_sourceAreas$clado_events_tables[[1]]) == F |
       "node_ht" %in% names(BSMs_w_sourceAreas$ana_events_tables[[1]]) == F){ 
      include.node_ht <- function(event_table, root.age){ event_table$node_ht <- root.age - event_table$time_bp; return(event_table)}  
      BSMs_w_sourceAreas$clado_events_tables <- lapply(BSMs_w_sourceAreas$clado_events_tables, include.node_ht, root.age = rootage) 
      BSMs_w_sourceAreas$ana_events_tables <- lapply(BSMs_w_sourceAreas$ana_events_tables, include.node_ht, root.age = rootage )
    }
    
    # structure to save the Events Through Time
    eventsCube <- array(NA, dim = c(length(areas), length(areas), (length(time.slices)-1) ),
                        dimnames = list(from=c(areas), to=c(areas) , time= time.slices[1:length(time.slices)-1]))
    eventsCube.states <- array(NA, dim = c(length(names_of_states), length(names_of_states), (length(time.slices)-1) ),
                               dimnames = list(from=c(names_of_states), to=c(names_of_states) , time= time.slices[1:length(time.slices)-1]))
    ETT.directional <- melt(eventsCube)[-4]
    ETT.directionalStates <- melt(eventsCube.states)[-4]
    ETT.beta.div <- melt(eventsCube)[-4] # these would be non directional distances between areas
    ETT.withinAreas <- melt(eventsCube[1,,])[-4]
    ETT.global <- data.frame(time = time.slices[1:length(time.slices)-1])
    Global.Disp.Rates.states.uncorr <- Global.Disp.Rates.states <- Global.Disp.Rates <- list()
    if(results.per.map) results.per.map2 <- list()
    
    rm(eventsCube)
    
    if("SimulationNumber" %in% names(BSM[[1]][[1]]) &  "Map" %in% names(BSM[[2]][[1]])  ){ # in case we have simulations
      Nsimulations <- max(sapply(BSM$clado_events_tables, function(x) x$SimulationNumber))
      MapsPerSimulation <- max(sapply(BSM$clado_events_tables, function(x) x$Map))
      print(paste0("The BSM (Biogeography Stochastic Map) object is a list of ", Nsimulations*MapsPerSimulation ," maps from ", Nsimulations, " simulations (", MapsPerSimulation, " maps each simulation)."))
      FirstMapEachSimulation <- seq(1, length(BSM$clado_events_tables), by=MapsPerSimulation)
      if(length(FirstMapEachSimulation)!=Nsimulations) stop("The index of simulations has a different length than the number of simulations")
    } else {
      Nsimulations <- 1
      FirstMapEachSimulation <- 1
      MapsPerSimulation <- length(BSM[[1]])
    }
    
  }
  
  ## Start of t loop ###
  for(t in 1:(length(time.slices)-1)){
    time.slice <- time.slices[t:(t+1)]
    time.slice.nh <- time.slices.nh[t:(t+1)]
    
    if(verbose) cat("\nTIME BIN", t, ":", time.slice[1], "-", time.slice[2], "MY.")
    if(verbose & any(BGBresults$inputs$timeperiods > time.slice[1] & BGBresults$inputs$timeperiods < time.slice[2])  ) cat(" Includes different time slices in the BioGeoBEARS model.")
    
    time <- time.slices[t]
    
    if(!any(c (sapply(BSMs_w_sourceAreas$clado_events_tables, 
                      FUN = function(x) any(x$time_bp > time.slice[1] & x$time_bp < time.slice[2]) ) ) )){
      if(verbose) cat("  No cladogenetic events. ")
    }
    if(any(is.na(BSMs_w_sourceAreas$ana_events_tables))) {
      if(verbose) cat(" No anagenetic events.")
    } else {
      if(!any(sapply(BSMs_w_sourceAreas$ana_events_tables, 
                     FUN = function(x) any(x$abs_event_time > time.slice[1] & x$abs_event_time < time.slice[2]) ))) {
        if(verbose) cat(" And no anagenetic events")
      } 
    }
    counts_list.t = count_ana_clado_eventsM(clado_events_tables = BSMs_w_sourceAreas$clado_events_tables,
                                            ana_events_tables = BSMs_w_sourceAreas$ana_events_tables,
                                            areanames= areas, actual_names = areas, timeperiod = time.slice, verbose = verbose)
    
    # 1. Additional calculations per map ####
    # Add to this object the number of lineages per area for each map and the dispersal rates (disp. events / branch length of source area). 
    # This involves doing several loops through the maps, which are coded independently for every process for clarity
    
    ## 1.1. Total branch length of each area per time bin ####
    {# The time spent in each area (within the bin), should be used to account for the dispersal rates as events/lineages/time
      # The result should be a rectangle of areas x map
      
      # get all the pieces of the tree (subbranches and subtrees) in that timeslice
      pieces_in_time_slice <- master_table[which(master_table$time_top == time.slice[1] & master_table$time_bot == time.slice[2]),]
      # exclude the root if it does not have a branch length
      if(any(pieces_in_time_slice$node.type == "root" & is.na(pieces_in_time_slice$edge.length) )) {
        pieces_in_time_slice <- pieces_in_time_slice[-which(pieces_in_time_slice$node.type == "root" & is.na(pieces_in_time_slice$edge.length)),]
      } 
      
      rounding.error.allowed <- 0.00001
      BL.perArea.inTSlice <- matrix(0, nrow=dim(counts_list.t$all_dispersals_counts_cube)[3], 
                                    ncol=length(areas), dimnames = list(c(), areas))
      BL.perState.inTSlice <- matrix(0, nrow=dim(counts_list.t$all_dispersals_counts_cube)[3], 
                                     ncol=length(BGBresults$inputs$states_list), dimnames = list(c(), names_of_states))
      
      for(mapI in 1:dim(counts_list.t$all_dispersals_counts_cube)[3]){
        mappedBranches <- BSMs_w_sourceAreas$clado_events_tables[[mapI]][which(BSMs_w_sourceAreas$clado_events_tables[[mapI]]$node %in% pieces_in_time_slice$node),]
        if("stratum" %in% names(mappedBranches)){ 
          if(any(mappedBranches$time_bot < time.slice[1] | mappedBranches$time_top > time.slice[2] ) ){
            mappedBranches <- mappedBranches[-which(mappedBranches$time_bot < time.slice[1] |
                                                      mappedBranches$time_top > time.slice[2] ),]
          }
        }
        
        for(branchJ in 1:nrow(mappedBranches)) {
          # Calculate length of the branch within the time slice. This includes new branches after speciation events
          nhBranch <- c(mappedBranches$node_ht[branchJ], mappedBranches$node_ht[branchJ]-mappedBranches$edge.length[branchJ])
          if("stratum" %in% names(mappedBranches)){ 
            nhTop <- master_table$time_bp[which(master_table$node.type=="root")] - mappedBranches$time_top[branchJ]
            nhBot <- master_table$time_bp[which(master_table$node.type=="root")] - mappedBranches$time_bot[branchJ]
            if(nhBot < 0) {nhBot = 0}
            nhBranch <- c( ifelse( nhTop<=(nhBranch[1]+rounding.error.allowed), nhTop, nhBranch[1]),
                           ifelse(nhBot > nhBranch[2], nhBot, nhBranch[2]) )
          }
          nhBranchWithinSlice <- c( ifelse( time.slice.nh[1]<=(nhBranch[1]+rounding.error.allowed), time.slice.nh[1], nhBranch[1]),
                                    ifelse(time.slice.nh[2] > nhBranch[2], time.slice.nh[2], nhBranch[2]) )
          if(nhBranchWithinSlice[1]-nhBranchWithinSlice[2] > time.slice.size) { stop("The piece of that branch is larger than the time.slice.size")}
          if(mappedBranches$sampled_states_AT_nodes[branchJ] == mappedBranches$sampled_states_AT_brbots[branchJ]) {# if there is no anagenetic change through that branch
            states.tt <- data.frame(state_num_1based = mappedBranches$sampled_states_AT_nodes[branchJ],
                                    state_text = names_of_states[mappedBranches$sampled_states_AT_nodes[branchJ]], # verified that this matches the biogeography matrix when applied to the tips of the tree (most recent timeslice)
                                    branch.length = round(nhBranchWithinSlice[1] - nhBranchWithinSlice[2],5))
          } else { # in case there are changes through that branch, There might be more than one.
            mappedAnaChanges <- BSMs_w_sourceAreas$ana_events_tables[[mapI]][which(BSMs_w_sourceAreas$ana_events_tables[[mapI]]$node == mappedBranches$node[branchJ]),]
            if(all((rootage - mappedAnaChanges$abs_event_time) > time.slice.nh[1]) ) { 
              # In case changes happen after that time slice. 
              states.tt <- data.frame(state_num_1based = mappedAnaChanges$current_rangenum_1based[1],
                                      state_text = mappedAnaChanges$current_rangetxt[1],             
                                      branch.length = round(nhBranchWithinSlice[1] - nhBranchWithinSlice[2],5))
            } else {
              # Select only the changes that happened within the time slice or before
              mappedAnaChanges <- mappedAnaChanges[which( (rootage - mappedAnaChanges$abs_event_time) < time.slice.nh[1]),]
              if(length(which((rootage - mappedAnaChanges$abs_event_time) < time.slice.nh[2]))>1) { # if there were more than 1 anagenetic events that happened before that time slice, keep only the most recent one
                mappedAnaChanges <- mappedAnaChanges[-c(1:(length(which((rootage - mappedAnaChanges$abs_event_time) < time.slice.nh[2]))-1)),] }
              if(any( (rootage - mappedAnaChanges$abs_event_time) > time.slice.nh[2])){ # if there is any changes within the time slice
                states.tt <- data.frame(startTime = c(nhBranchWithinSlice[2], rootage-mappedAnaChanges$abs_event_time),
                                        endTime = c(rootage-mappedAnaChanges$abs_event_time, nhBranchWithinSlice[1]),
                                        state_text = c(mappedAnaChanges$current_rangetxt, mappedAnaChanges$new_rangetxt[nrow(mappedAnaChanges)]),
                                        state_num_1based = c(mappedAnaChanges$current_rangenum_1based, mappedAnaChanges$new_rangenum_1based[nrow(mappedAnaChanges)]))
                states.tt$branch.length <- states.tt$endTime - states.tt$startTime
              } else { # if there are no changes within the time slice. We still account for the changes before
                states.tt <- data.frame(state_num_1based = mappedAnaChanges$new_rangenum_1based,
                                        state_text = mappedAnaChanges$new_rangetxt,
                                        branch.length = round(nhBranchWithinSlice[1] - nhBranchWithinSlice[2],5))
              } # end of else
            } # end of else
          } # end of else
          
          if(sum(states.tt$branch.length)>time.slice.size) stop(paste0("The sum of the branch lengths for branch ", branchJ, " at map ", mapI, " in time ", t, " is bigger than the time.slice.size. FIX THIS!!!"))
          
          # use states.tt to sum the right amount of branch lengths to each state and area
          # It is good that this is summarized for both areas and states (for later getting probabilities of ancestry by branch section)
          for(k in 1:nrow(states.tt)){
            BL.perState.inTSlice[mapI, states.tt$state_num_1based[k] ] <- BL.perState.inTSlice[mapI, states.tt$state_num_1based[k] ] + states.tt$branch.length[k]
            
            areasAffected <- unlist(strsplit(states.tt$state_text[k],split = ""))
            BL.perArea.inTSlice[mapI, which(areas %in% areasAffected) ] <- BL.perArea.inTSlice[mapI, which(areas %in% areasAffected) ] + states.tt$branch.length[k]
          } # end of k loop
          rm(states.tt)
        } # end for branchJ loop
        if(mapI!=1) if(round(sum(BL.perState.inTSlice[mapI,])-sum(BL.perState.inTSlice[mapI-1,]),3)!=0) print(paste0("BL problem with map ", mapI, " at time ", t))
      } # end for mapI loop
      
      counts_list.t$branch.lengths.per.state <- BL.perState.inTSlice 
      counts_list.t$branch.lengths.per.area <- BL.perArea.inTSlice
      
      if(#length(unique(round(rowSums(BL.perState.inTSlice),3) ))>1
        any(round(rowSums(BL.perState.inTSlice)[-1] - rowSums(BL.perState.inTSlice)[-dim(counts_list.t$all_dispersals_counts_cube)[3]], 3)!=0)
      ) stop(paste0("The sums of all the branch lengths for every state are different between maps at time slice ", t, ". FIX THIS!!!"))
      rowSums(BL.perArea.inTSlice)
      
      rm(BL.perArea.inTSlice, BL.perState.inTSlice, mappedBranches, mapI, nhBot, nhTop, nhBranch, nhBranchWithinSlice, mappedAnaChanges, k, areasAffected, branchJ)
    }
    
    ## 1.2. Sympatric speciation (both y and s) ####
    {names.subCounts <- substr(sub(".*->", "", names(counts_list.t$unique_sub_counts)), 1, 1)
    Gains_subset_symp <- data.frame(matrix(0, nrow = length(BSMs_w_sourceAreas$clado_events_tables), ncol = length(areas),
                                           dimnames = list(c(),areas)))
    Gains_y_symp <- Gains_subset_symp
    for(a in 1:length(areas)){
      if(length(which(names.subCounts==areas[a])) > 1){
        Gains_subset_symp[,a] <- rowSums(counts_list.t$unique_sub_counts[,which(names.subCounts==areas[a])])
      }
      if(length(which(names.subCounts==areas[a])) == 1){
        Gains_subset_symp[,a] <- counts_list.t$unique_sub_counts[,which(names.subCounts==areas[a])]
      }
      if(length (which(substr(names(counts_list.t$unique_sym_counts), 1, 1) == areas[a])) ==1 ) {
        Gains_y_symp[,a] <- counts_list.t$unique_sym_counts[, which(substr(names(counts_list.t$unique_sym_counts), 1, 1) == areas[a] ) ]
      }
    }
    Gains_all_symp <- Gains_subset_symp + Gains_y_symp   # Interestingly, there is less variation when counting all sympatric speciation events. When separating subset and y sympatric speciation, the s>
    counts_list.t$y_counts_rectangle <- Gains_y_symp
    counts_list.t$s_counts_rectangle <- Gains_subset_symp
    counts_list.t$all_symp_counts_rectangle <- Gains_all_symp
    rm(Gains_subset_symp, Gains_y_symp, Gains_all_symp)
    }
    
    ## 1.3. Immigration and emigration ####
    # we should get matrices of n areas x n maps
    counts_list.t$immigration_counts_rectangle <- apply(counts_list.t$all_dispersals_counts_cube, c(3,2), sum) # dimension 2 is "to"
    counts_list.t$emigration_counts_rectangle <- apply(counts_list.t$all_dispersals_counts_cube, c(3,1), sum) # dimension 1 is "from"
    
    counts_list.t$netMig_rectangle <- counts_list.t$emigration_counts_rectangle - counts_list.t$immigration_counts_rectangle # Net migration (source - sink)
    
    w.emigration <- counts_list.t$emigration_counts_rectangle/counts_list.t$all_dispersals_totals_list # This operation goes row-wise
    w.immigration <- counts_list.t$immigration_counts_rectangle/counts_list.t$all_dispersals_totals_list # This operation goes row-wise
    # There is a problem here, as x/0 = Inf and 0/0 = NaN
    w.emigration[is.nan(w.emigration)] <- 0
    w.immigration[is.nan(w.immigration)] <- 0
    counts_list.t$w.netMig_rectangle <- w.emigration - w.immigration # Weighted Net migration (w.source - w.sink)
    rm(w.emigration, w.immigration)
    
    # counts_list.t$outindisp.ratio_rectangle <- counts_list.t$emigration_counts_rectangle / counts_list.t$immigration_counts_rectangle # Outward:Inward dispersal ratio
    # There is a problem here, as x/0 = Inf and 0/0 = NaN
    
    ## 1.4. Taxonomic and phylogenetic diversity ####
    {# Both at the start and end of each period
      # In this case we need the branch lengths until the root of the tree, so we will need to slice the tree again.
      
      existingPhyDiv <- existingLineages <- matrix(NA, ncol = length(areas), nrow = 0, dimnames = list(c(), areas) )
      counts_list.t$SpBetaDiv.nes_cube <- counts_list.t$SpBetaDiv.sim_cube <- counts_list.t$SpBetaDiv.sor_cube <-
        counts_list.t$PhyBetaDiv.nes_cube <- counts_list.t$PhyBetaDiv.sim_cube <- counts_list.t$PhyBetaDiv.sor_cube <- 
        array(NA, dim = c(length(areas),  length(areas), MapsPerSimulation*Nsimulations))
      
      # Diversities at end of each time slice
      if(t == 1){ # There are two ways of obtaining the current lineages per area
        for(map1 in FirstMapEachSimulation){
          # 1.4.1. Current diversities from biogeography matrix
          if(is.null(statenames)){
            print(paste0("Calculating current lineages per area from biogeography matrix specified in BGBresults$inputs$geogfn: ", BGBresults$inputs$geogfn))
            biogeoMatrix <- read_PHYLIP_data(BGBresults$inputs$geogfn) # This is not feasible if we do not keep the biogeography input files (for example in simulations), but we can recover this information from the tips of the stochastic maps and the customised states_list object with statenames
            
            # phylogenetic diversity
            biogeo <- t(biogeoMatrix)
            biogeo <- apply(biogeo, 2, as.numeric)
            row.names(biogeo) <- areas #names(read_PHYLIP_data(BGBresults$inputs$geogfn))
          } else { # when there are statenames
            # 1.4.2. Current diversities from tips of stochastic maps (using the statenames to recognize them) 
            # The following lines are an alternative in case we do not have the biogeography input file(s), but we have the states_list object with names
            if(map1 == 1) print("Calculating current lineages per area from stochastic maps and given state names")
            tip.states <- BSMs_w_sourceAreas$clado_events_tables[[map1]]$sampled_states_AT_nodes[which(BSMs_w_sourceAreas$clado_events_tables[[map1]]$node.type=="tip")] # 1based
            tip.states.letter <- unlist(statenames)[as.numeric(as.character(tip.states))]
            tip.areas <- strsplit(tip.states.letter, "")
            tip.areas.count <- data.frame(table(unlist(tip.areas)))
            
            biogeo <- matrix(data = 0, nrow = length(areas), ncol = length(which(BSMs_w_sourceAreas$clado_events_tables[[map1]]$node.type=="tip")),
                             dimnames = c(list(areas),
                                          list(BSMs_w_sourceAreas$clado_events_tables[[map1]]$label[which(BSMs_w_sourceAreas$clado_events_tables[[map1]]$node.type=="tip")])  ))
            for(sp in 1:ncol(biogeo)){
              biogeo[which(rownames(biogeo) %in% tip.areas[[sp]]), sp] <- 1
            }
            
          } # end of else is.null(statenames)
          
          Div <- pd(samp = biogeo, tree = read.tree(BGBresults$inputs$trfn), include.root = T)
          existLin <- as.data.frame(t(Div$SR))
          existPD <- as.data.frame(t(Div$PD))
          names(existLin) <- names(existPD) <- areas
          
          existLin <- existLin[rep(1, MapsPerSimulation),]
          existingLineages <- rbind(existingLineages, existLin)
          existPD <- existPD[rep(1, MapsPerSimulation),]
          existingPhyDiv <- rbind(existingPhyDiv, existPD)
          
          PhyBetaDiv <- phylo.beta.pair(biogeo,read.tree(BGBresults$inputs$trfn))
          SpBetaDiv <- beta.pair(biogeo)
          
          counts_list.t$PhyBetaDiv.sor_cube[,,map1:(map1+MapsPerSimulation-1)] <- array(rep(as.matrix(PhyBetaDiv$phylo.beta.sor), MapsPerSimulation), 
                                                                                        dim = c(length(areas),  length(areas), MapsPerSimulation) ) 
          counts_list.t$PhyBetaDiv.sim_cube[,,map1:(map1+MapsPerSimulation-1)] <- array(rep(as.matrix(PhyBetaDiv$phylo.beta.sim), MapsPerSimulation), 
                                                                                        dim = c(length(areas),  length(areas), MapsPerSimulation) )
          counts_list.t$PhyBetaDiv.nes_cube[,,map1:(map1+MapsPerSimulation-1)] <- array(rep(as.matrix(PhyBetaDiv$phylo.beta.sne), MapsPerSimulation), 
                                                                                        dim = c(length(areas),  length(areas), MapsPerSimulation) )
          counts_list.t$SpBetaDiv.sor_cube[,,map1:(map1+MapsPerSimulation-1)] <- array(rep(as.matrix(SpBetaDiv$beta.sor), MapsPerSimulation), 
                                                                                       dim = c(length(areas),  length(areas), MapsPerSimulation) )
          counts_list.t$SpBetaDiv.sim_cube[,,map1:(map1+MapsPerSimulation-1)] <- array(rep(as.matrix(SpBetaDiv$beta.sim), MapsPerSimulation), 
                                                                                       dim = c(length(areas),  length(areas), MapsPerSimulation) )
          counts_list.t$SpBetaDiv.nes_cube[,,map1:(map1+MapsPerSimulation-1)] <- array(rep(as.matrix(SpBetaDiv$beta.sne), MapsPerSimulation), 
                                                                                       dim = c(length(areas),  length(areas), MapsPerSimulation) )
        } # end of map1 loop
        
        counts_list.t$existingLineages <- existingLineages
        counts_list.t$existingPhyDiv <- existingPhyDiv
        
      } else { # when t is not 1
        
        # section the tree in two parts, to count branch lengths until the root of the tree for calculating phylogenetic diversity in each community
        BGB_run_object <- define_BioGeoBEARS_run(trfn = BGBresults$inputs$trfn)
        BGB_run_object <- readfiles_BioGeoBEARS_run(BGB_run_object)
        BGB_run_object$timeperiods <- c(time.slice[1], time.slices[length(time.slices)])
        invisible(capture.output( #https://stackoverflow.com/questions/34208564/how-to-hide-or-disable-in-function-printed-message
          tree_sectioning <- section_the_tree(BGB_run_object, make_master_table = T, plot_pieces = F)   ))
        new.tree <- tree_sectioning$tree_sections_list[[2]]$return_pieces_list[[1]]
        master_table2 <- tree_sectioning$master_table
        rm(tree_sectioning, BGB_run_object)
        new.tips <- master_table2$node[which(master_table2$time_top==time.slice[1] & master_table2$SUBnode.type == "tip") ]
        
        if(!all(new.tree$tip.label %in% master_table2$SUBlabel[which(master_table2$time_top==time.slice[1] & master_table2$SUBnode.type == "tip")]) ) stop("The tip labels of the new.tree do not match the SUBlabels in master_table2") 
        
        counts_list.t$existingPhyDiv <- counts_list.t$existingLineages <- matrix(NA, ncol = length(areas), nrow = dim(counts_list.t$all_dispersals_counts_cube)[3], dimnames = list(c(), areas) )
        counts_list.t$PhyBetaDiv.nes_cube <- counts_list.t$SpBetaDiv.nes_cube <- counts_list.t$PhyBetaDiv.sim_cube <- 
          counts_list.t$SpBetaDiv.sim_cube <- counts_list.t$PhyBetaDiv.sor_cube <- counts_list.t$SpBetaDiv.sor_cube <-
          array(NA, dim = dim(counts_list.t$all_dispersals_counts_cube) )
        
        for(mapD in 1:dim(counts_list.t$all_dispersals_counts_cube)[3]){
          tip.values <- matrix(0, nrow=length(areas), ncol= length(new.tree$tip.label), #length(new.tips), 
                               dimnames = list(c(areas), c(new.tree$tip.label
                                                           #master_table2$label[which(master_table2$time_bot==time.slice[1])]
                               )))
          mapD_clado <- BSMs_w_sourceAreas$clado_events_tables[[mapD]]# we create a copy of this so that we can modify in the time-sliced cases. It gets overwritten for every map and time slice 
          if("stratum" %in% names(BSMs_w_sourceAreas$clado_events_tables[[mapD]])){ # if the BioGeoBEARS analysis was time-sliced
            mapD_clado <- mapD_clado[which(mapD_clado$time_top <= time.slice[1] &
                                             mapD_clado$time_bot >= time.slice[2]),]
            if(nrow(mapD_clado)==0) { 
              mapD_clado <- BSMs_w_sourceAreas$clado_events_tables[[mapD]]
            }
          }
          
          for(b in 1:length(new.tips)){
            if(all(mapD_clado$anagenetic_events_txt_below_node[which(mapD_clado$node %in% new.tips[b] ) ] %in% c("none", NA) )){ # if there are no anagenetic events
              state_text = names_of_states[mapD_clado$sampled_states_AT_nodes[which(mapD_clado$node %in% new.tips[b]) ] ]
            } else { # if there are anagenetic changes
              anaChanges <- BSMs_w_sourceAreas$ana_events_tables[[mapD]][which(BSMs_w_sourceAreas$ana_events_tables[[mapD]]$node %in% new.tips[b] ), ]
              anaChanges <- anaChanges[which(anaChanges$node_ht-anaChanges$event_time < time.slice.nh[1]),] # keep only the anagenetic events that happen before timeslice
              if(nrow(anaChanges)==0) {
                state_text <- names_of_states[mapD_clado$sampled_states_AT_nodes[which(mapD_clado$node %in% new.tips[b] )] ]
                if(length(state_text) > 1) {
                  state_text <- names_of_states[mapD_clado$sampled_states_AT_nodes[which(mapD_clado$node %in% new.tips[b] &
                                                                                           mapD_clado$anagenetic_events_txt_below_node %in% c("none", NA))] ]
                  if(length(state_text)==0 ){
                    state_text <- names_of_states[mapD_clado$sampled_states_AT_nodes[which(mapD_clado$node %in% new.tips[b] )] ] [-1]
                  } 
                  
                }
              } else { # if the anagenetic changes happened in that t time-slice
                anaChanges <- anaChanges[nrow(anaChanges),]
                state_text = anaChanges$new_rangetxt
              } # end of else when the ana events are in that t time-slice
            } # end of else, when there are anagenetic changes
            if(length(state_text)>1){ # This was added because in the structure of clado_events_tables in time-sliced analysis, there are multiple rows per branch
              if(all(state_text[-1]==state_text[1])) state_text <- state_text[1] else {
                
                stop(paste("There is a problem with multiple different values in state_text for time slice", t, ", map =", mapD, ", new.tip b =", b)) }
            }
            
            areasAffected <- unique(unlist(strsplit(state_text,split = "")))
            tip.values[which(areas %in% areasAffected), b ] <- 1
          } # end of b loop, now with the tip.values matrix complete, calculate diversities for all the communities, and save them:
          failing_tips <- which(colSums(tip.values)==0)
          if(any(colSums(tip.values)==0)) stop(paste("There are", length(which(colSums(tip.values)==0)), "out of", 
                                                     ncol(tip.values), "tip labels of new.tree that did not have any assigned area"))
          AlphaDiv <- pd(tip.values, new.tree, include.root = T)
          PhyBetaDiv <- phylo.beta.pair(tip.values,new.tree)
          SpBetaDiv <- beta.pair(tip.values)
          
          PhyBetaDiv <- lapply(PhyBetaDiv, as.matrix)
          SpBetaDiv <- lapply(SpBetaDiv, as.matrix)
          
          counts_list.t$existingPhyDiv[mapD,] <- AlphaDiv$PD
          counts_list.t$existingLineages[mapD,] <- AlphaDiv$SR
          counts_list.t$PhyBetaDiv.sor_cube[,,mapD] <- PhyBetaDiv$phylo.beta.sor
          counts_list.t$PhyBetaDiv.sim_cube[,,mapD] <- PhyBetaDiv$phylo.beta.sim
          counts_list.t$PhyBetaDiv.nes_cube[,,mapD] <- PhyBetaDiv$phylo.beta.sne
          counts_list.t$SpBetaDiv.sor_cube[,,mapD] <- SpBetaDiv$beta.sor
          counts_list.t$SpBetaDiv.sim_cube[,,mapD] <- SpBetaDiv$beta.sim
          counts_list.t$SpBetaDiv.nes_cube[,,mapD] <- SpBetaDiv$beta.sne   
        } # end of mapD loop
      } # end of else when t is not 1
      
      if(any(c(is.nan(counts_list.t$SpBetaDiv.sim_cube), is.nan(counts_list.t$SpBetaDiv.nes_cube),
               is.nan(counts_list.t$PhyBetaDiv.sim_cube), is.nan(counts_list.t$PhyBetaDiv.nes_cube) )) ) {
        if(verbose) print(paste0("There are NaN in some of the Simpson and nestedness dissimilarities for time ", t, 
                                 ": ", time.slice[1], "-", time.slice[2], " MY, because areas ",
                                 paste(colnames(counts_list.t$existingLineages)[ unique(which(counts_list.t$existingLineages==0, arr.ind = T)[,2]) ], collapse = " and "),
                                 " sometimes have 0 species. Omitting ",
                                 length(unique(which(counts_list.t$existingLineages == 0, arr.ind = T)[,1]) ),
                                 #length(which(rowSums(counts_list.t$existingLineages == 0) >= 2)) ,
                                 " maps out of ",  nrow(counts_list.t$existingLineages), 
                                 ", to calculate the average values of Simpson and nestedness dissimilarities between those areas") )
      }
      rm(existingLineages, existingPhyDiv, AlphaDiv, PhyBetaDiv, SpBetaDiv, areasAffected, tip.values, failing.tips, anaChanges, mapD, mapD_clado, newtips)
    }
    
    
    
    
    ## 1.5 Dispersal, extirpation and sympatric speciation rates ####
    {# For dispersal we need a cube of n areas x n areas x n maps
      
      counts_list.t$w.disp0_cube <- counts_list.t$w.disp_cube  <- array(NA, dim=dim(counts_list.t$all_dispersals_counts_cube))
      
      
      # For the sympatric speciation and extirpation rates we need a rectangle of n areas x n maps
      counts_list.t$w.symp_rectangle <- counts_list.t$w.extirp_rectangle <- matrix(NA, ncol = length(areas), nrow = dim(counts_list.t$all_dispersals_counts_cube)[3])
      
      counts_list.t$global.disp.rate.states <- counts_list.t$global.disp.rate <- 
        rep(NA, dim(counts_list.t$all_dispersals_counts_cube)[3]) # Global dispersal in that time period. Vectors of the length of the number of maps (for both areas and states)
      
      
      for (mapd in 1:dim(counts_list.t$w.disp_cube)[3]){
        
        total.bl.states <- sum(as.numeric(counts_list.t$branch.lengths.per.state[mapd,]))
        
        counts_list.t$w.disp_cube[, , mapd] <- sweep(counts_list.t$all_dispersals_counts_cube[, , mapd], MARGIN = 1, 
                                                     as.numeric(counts_list.t$branch.lengths.per.area[mapd,]),
                                                     function(x,y) ifelse(y==0, NA, x/y )) # When the denominator is 0, x/0= Inf; and when the numerator and denominator are 0, 0/0 = NaN. Therefore instead of just using the function '/', we use an ifelse() so that no rates are calculated in those cases (before we were assuming those were 0s)
        
        counts_list.t$w.disp0_cube[, , mapd] <- sweep(counts_list.t$all_dispersals_counts_cube[, , mapd], MARGIN = 1, 
                                                      as.numeric(counts_list.t$branch.lengths.per.area[mapd,]),
                                                      function(x,y) ifelse(y==0, 0, x/y )) # When the denominator is 0, x/0= Inf; and when the numerator and denominator are 0, 0/0 = NaN. Therefore instead of just using the function '/', we use an ifelse() so that rates are assumed to be 0.
        
        counts_list.t$w.extirp_rectangle[mapd,] <- apply(as.data.frame(counts_list.t$e_counts_rectangle[mapd,]), MARGIN = 2,
                                                         function(x) ifelse(as.numeric(counts_list.t$branch.lengths.per.area[mapd,])==0, NA,
                                                                            x/as.numeric(counts_list.t$branch.lengths.per.area[mapd,]) ))  # When the denominator is 0, x/0= Inf; and when the numerator and denominator are 0, 0/0 = NaN. Therefore instead of just using the function '/', we use an ifelse() so that no rates are calculated in those cases (before we were assuming those were 0s)
        
        counts_list.t$w.symp_rectangle[mapd,] <- apply(counts_list.t$all_symp_counts_rectangle[mapd,], MARGIN = 1,
                                                       function(x) ifelse(as.numeric(counts_list.t$branch.lengths.per.area[mapd,])==0, NA,
                                                                          x/as.numeric(counts_list.t$branch.lengths.per.area[mapd,]) ))  # When the denominator is 0, x/0= Inf; and when the numerator and denominator are 0, 0/0 = NaN. Therefore instead of just using the function '/', we use an ifelse() so that no rates are calculated in those cases (before we were assuming those were 0s)
        
        total.disp.events <- sum(counts_list.t$all_dispersals_counts_cube[, , mapd])
        total.disp.events.states <- sum(counts_list.t$all_dispersals_counts_cube[, , mapd])
        
        counts_list.t$global.disp.rate[mapd] <- total.disp.events/total.bl.states
        counts_list.t$global.disp.rate.states[mapd] <- total.disp.events.states/total.bl.states
        
      } # end of mapd loop
      
      if(any(is.na(counts_list.t$w.disp_cube)) & verbose) {
        print(paste0("There are ", length(unique(which(counts_list.t$branch.lengths.per.area == 0, arr.ind = T)[,1])) , " maps, in time ", t,  
                     ": ", time.slice[1], "-", time.slice[2], " MY, with some NAs as dispersal rates from area(s) ",
                     paste(colnames(counts_list.t$branch.lengths.per.area)[ unique(which(counts_list.t$branch.lengths.per.area==0, arr.ind = T)[,2]) ], collapse = " and "),
                     " because they sometimes have 0 branch length"
        ))
      }
    }
    
    # END of calculations per map #
    
    # Saving calculations per map #
    if(results.per.map){ results.per.map2[[t]] <- counts_list.t
    } 
    
    # 2. Directional between area events ####
    {# Dispersal, which includes anagenetic dispersal(d), range switching (a, not included in the traditional DEC model, but in others such as DECX) and cladogenetic dispersal (j, founder events)
      # BGBresults$outputs@params_table$est[1] # rate of anagenetic dispersals d
      
      rows.ETTd.t <- ((t-1)*length(areas)*length(areas) + 1):(t*length(areas)*length(areas))
      
      # To avoid message from melt() function, include id.vars=NULL  # https://stackoverflow.com/questions/48774207/in-r-using-melt-how-can-i-hide-warning-messages
      ETT.directional$all.dispersals[rows.ETTd.t] <- melt(counts_list.t$all_dispersals_counts_fromto_means, id.vars=NULL)$value
      ETT.directional$all.dispersals.sd[rows.ETTd.t] <- melt(counts_list.t$all_dispersals_counts_fromto_sds, id.vars=NULL)$value
      ETT.directional$all.dispersals.UQ[rows.ETTd.t] <- apply(counts_list.t$all_dispersals_counts_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$all.dispersals.LQ[rows.ETTd.t] <- apply(counts_list.t$all_dispersals_counts_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      
      ETT.directional$ana.dispersals[rows.ETTd.t] <- melt(counts_list.t$ana_dispersals_counts_fromto_means, id.vars=NULL)$value   # includes a and d
      ETT.directional$ana.dispersals.sd[rows.ETTd.t] <- melt(counts_list.t$ana_dispersals_counts_fromto_sds, id.vars=NULL)$value
      ETT.directional$ana.dispersals.UQ[rows.ETTd.t] <- apply(counts_list.t$ana_counts_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$ana.dispersals.LQ[rows.ETTd.t] <- apply(counts_list.t$ana_counts_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      
      ETT.directional$founder[rows.ETTd.t] <- melt(counts_list.t$founder_counts_fromto_means, id.vars=NULL)$value  # j is founder event (cladogenetic dispersal)
      ETT.directional$founder.sd[rows.ETTd.t] <- melt(counts_list.t$founder_counts_fromto_sds, id.vars=NULL)$value
      ETT.directional$founder.UQ[rows.ETTd.t] <- apply(counts_list.t$founder_counts_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$founder.LQ[rows.ETTd.t] <- apply(counts_list.t$founder_counts_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      
      ETT.directional$a[rows.ETTd.t] <- melt(counts_list.t$a_counts_fromto_means, id.vars=NULL)$value   # a is anagenetic range switching (the same species loses one are as it colonises another one)
      ETT.directional$a.sd[rows.ETTd.t] <- melt(counts_list.t$a_counts_fromto_sds, id.vars=NULL)$value
      ETT.directional$a.UQ[rows.ETTd.t] <- apply(counts_list.t$a_counts_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$a.LQ[rows.ETTd.t] <- apply(counts_list.t$a_counts_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      
      ETT.directional$d[rows.ETTd.t] <- melt(counts_list.t$d_counts_fromto_means, id.vars=NULL)$value  # d is anagenetic range expansion (the same species colonises one area, expanding its range)
      ETT.directional$d.sd[rows.ETTd.t] <- melt(counts_list.t$d_counts_fromto_sds, id.vars=NULL)$value
      ETT.directional$d.UQ[rows.ETTd.t] <- apply(counts_list.t$d_counts_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$d.LQ[rows.ETTd.t] <- apply(counts_list.t$d_counts_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      
      ETT.directional$w.disp.rates[rows.ETTd.t] <- apply(counts_list.t$w.disp_cube, c(1,2), mean, na.rm=T)  # This includes all dispersal event types, divided by the branch length of their area of origin in that time slice
      ETT.directional$w.disp.rates.sd[rows.ETTd.t] <- apply(counts_list.t$w.disp_cube, c(1,2), sd, na.rm=T)
      ETT.directional$w.disp.rates.LQ[rows.ETTd.t] <- apply(counts_list.t$w.disp_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$w.disp.rates.UQ[rows.ETTd.t] <- apply(counts_list.t$w.disp_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      ETT.directional$w.disp.rates.N[rows.ETTd.t] <- apply(counts_list.t$w.disp_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each direction of dispersal (because in some cases it varies between maps, given that the denominator in some cases can be 0)
      
      ETT.directional$w.disp.rates0[rows.ETTd.t] <- apply(counts_list.t$w.disp0_cube, c(1,2), mean, na.rm=T)  # This includes all dispersal event types, divided by the branch length of their area of origin in that time slice. Assuming it is 0 when branch length is 0.
      ETT.directional$w.disp.rates0.sd[rows.ETTd.t] <- apply(counts_list.t$w.disp0_cube, c(1,2), sd, na.rm=T)
      ETT.directional$w.disp.rates0.LQ[rows.ETTd.t] <- apply(counts_list.t$w.disp0_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
      ETT.directional$w.disp.rates0.UQ[rows.ETTd.t] <- apply(counts_list.t$w.disp0_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
      ETT.directional$w.disp.rates0.N[rows.ETTd.t] <- apply(counts_list.t$w.disp0_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each direction of dispersal (because in some cases it varies between maps, given that the denominator in some cases can be 0)
      
    }
    
    
    # 3. Beta diversities between areas ####
    {ETT.beta.div$sor[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sor_cube, c(1,2), mean, na.rm=T)
    ETT.beta.div$sor.sd[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sor_cube, c(1,2), sd, na.rm=T)
    ETT.beta.div$sor.LQ[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sor_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
    ETT.beta.div$sor.UQ[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sor_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
    ETT.beta.div$sor.N[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sor_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each pairwise comparison (because in some cases it varies between maps, given that the denominator in some cases can be 0)
    
    ETT.beta.div$sim[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sim_cube, c(1,2), mean, na.rm =T) # The pairwise comparisons in the maps that produce NaNs are excluded by using na.rm
    ETT.beta.div$sim.sd[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sim_cube, c(1,2), sd, na.rm =T)
    ETT.beta.div$sim.LQ[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sim_cube, c(1,2), quantile, probs = 0.025, na.rm = T) 
    ETT.beta.div$sim.UQ[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sim_cube, c(1,2), quantile, probs = 0.975, na.rm = T)
    ETT.beta.div$sim.N[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.sim_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each pairwise comparison (because in some cases it varies between maps, given that the denominator in some cases can be 0)
    
    ETT.beta.div$nes[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.nes_cube, c(1,2), mean, na.rm =T)
    ETT.beta.div$nes.sd[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.nes_cube, c(1,2), sd, na.rm =T)
    ETT.beta.div$nes.LQ[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.nes_cube, c(1,2), quantile, probs = 0.025, na.rm = T)
    ETT.beta.div$nes.UQ[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.nes_cube, c(1,2), quantile, probs = 0.975, na.rm = T)
    ETT.beta.div$nes.N[rows.ETTd.t] <- apply(counts_list.t$SpBetaDiv.nes_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each pairwise comparison (because in some cases it varies between maps, given that the denominator in some cases can be 0)
    
    ETT.beta.div$phylo.sor[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sor_cube, c(1,2), mean, na.rm=T)
    ETT.beta.div$phylo.sor.sd[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sor_cube, c(1,2), sd, na.rm=T)
    ETT.beta.div$phylo.sor.LQ[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sor_cube, c(1,2), quantile, probs = 0.025, na.rm=T)
    ETT.beta.div$phylo.sor.UQ[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sor_cube, c(1,2), quantile, probs = 0.975, na.rm=T)
    ETT.beta.div$phylo.sor.N[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sor_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each pairwise comparison (because in some cases it varies between maps, given that the denominator in some cases can be 0)
    
    ETT.beta.div$phylo.sim[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sim_cube, c(1,2), mean, na.rm =T)
    ETT.beta.div$phylo.sim.sd[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sim_cube, c(1,2), sd, na.rm =T)
    ETT.beta.div$phylo.sim.LQ[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sim_cube, c(1,2), quantile, probs = 0.025, na.rm = T)
    ETT.beta.div$phylo.sim.UQ[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sim_cube, c(1,2), quantile, probs = 0.975, na.rm = T)
    ETT.beta.div$phylo.sim.N[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.sim_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each pairwise comparison (because in some cases it varies between maps, given that the denominator in some cases can be 0)
    
    ETT.beta.div$phylo.nes[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.nes_cube, c(1,2), mean, na.rm =T)
    ETT.beta.div$phylo.nes.sd[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.nes_cube, c(1,2), sd, na.rm =T)
    ETT.beta.div$phylo.nes.LQ[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.nes_cube, c(1,2), quantile, probs = 0.025, na.rm = T)
    ETT.beta.div$phylo.nes.UQ[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.nes_cube, c(1,2), quantile, probs = 0.975, na.rm = T)
    ETT.beta.div$phylo.nes.N[rows.ETTd.t] <- apply(counts_list.t$PhyBetaDiv.nes_cube, c(1,2), FUN = function(x) length(which(!is.na(x) )) ) # This is to track the sample size for each pairwise comparison (because in some cases it varies between maps, given that the denominator in some cases can be 0)
    }
    
    
    # 4. Within Area events ####
    {rows.ETTw.t <- ((t-1)*length(areas) + 1):(t*length(areas))
    
    # 4.1. Extirpation
    # e is extinction (only partial, as extirpation, range contraction by loosing one area in the range)
    ETT.withinAreas$e[rows.ETTw.t] <- colMeans(counts_list.t$e_counts_rectangle, na.rm=T)
    ETT.withinAreas$e.LQ[rows.ETTw.t] <- apply(counts_list.t$e_counts_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$e.UQ[rows.ETTw.t] <- apply(counts_list.t$e_counts_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$e.rates[rows.ETTw.t] <- colMeans(counts_list.t$w.extirp_rectangle, na.rm=T)
    ETT.withinAreas$e.rates.LQ[rows.ETTw.t] <- apply(counts_list.t$w.extirp_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$e.rates.UQ[rows.ETTw.t] <- apply(counts_list.t$w.extirp_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    # 4.2. Sympatric speciation
    # There's two variants of sympatric speciation
    #       s is subset sympatric speciation (increases species diversity locally)
    #       y is sympatric speciation (increases species diversity locally)
    
    ETT.withinAreas$sympatryY[rows.ETTw.t] <- colMeans(counts_list.t$y_counts_rectangle, na.rm=T)
    ETT.withinAreas$sympatryY.LQ[rows.ETTw.t] <- apply(counts_list.t$y_counts_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$sympatryY.UQ[rows.ETTw.t] <- apply(counts_list.t$y_counts_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$sympatryS[rows.ETTw.t] <- colMeans(counts_list.t$s_counts_rectangle, na.rm=T)
    ETT.withinAreas$sympatryS.LQ[rows.ETTw.t] <- apply(counts_list.t$s_counts_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$sympatryS.UQ[rows.ETTw.t] <- apply(counts_list.t$s_counts_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$sympatryAll[rows.ETTw.t] <- colMeans(counts_list.t$all_symp_counts_rectangle, na.rm=T)
    ETT.withinAreas$sympatryAll.LQ[rows.ETTw.t] <- apply(counts_list.t$all_symp_counts_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$sympatryAll.UQ[rows.ETTw.t] <- apply(counts_list.t$all_symp_counts_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$sympatric.speciation.rate[rows.ETTw.t] <- colMeans(counts_list.t$w.symp_rectangle, na.rm=T)
    ETT.withinAreas$sympatric.speciation.rate.LQ[rows.ETTw.t] <- apply(counts_list.t$w.symp_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$sympatric.speciation.rate.UQ[rows.ETTw.t] <- apply(counts_list.t$w.symp_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    # 4.3. Immigration and emigration
    ETT.withinAreas$immigration[rows.ETTw.t] <- colMeans(counts_list.t$immigration_counts_rectangle, na.rm=T)[match(ETT.withinAreas$to[rows.ETTw.t], areas )]
    ETT.withinAreas$immigration.LQ[rows.ETTw.t] <- apply(counts_list.t$immigration_counts_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$immigration.UQ[rows.ETTw.t] <- apply(counts_list.t$immigration_counts_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$emigration[rows.ETTw.t] <- colMeans(counts_list.t$emigration_counts_rectangle, na.rm=T)[match(ETT.withinAreas$to[rows.ETTw.t], areas )]
    ETT.withinAreas$emigration.LQ[rows.ETTw.t] <- apply(counts_list.t$emigration_counts_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$emigration.UQ[rows.ETTw.t] <- apply(counts_list.t$emigration_counts_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$netMig[rows.ETTw.t] <- colMeans(counts_list.t$netMig_rectangle, na.rm=T)[match(ETT.withinAreas$to[rows.ETTw.t], areas )]
    ETT.withinAreas$netMig.LQ[rows.ETTw.t] <- apply(counts_list.t$netMig_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$netMig.UQ[rows.ETTw.t] <- apply(counts_list.t$netMig_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    ETT.withinAreas$w.netMig[rows.ETTw.t] <- colMeans(counts_list.t$w.netMig_rectangle, na.rm=T)[match(ETT.withinAreas$to[rows.ETTw.t], areas )]
    ETT.withinAreas$w.netMig.LQ[rows.ETTw.t] <- apply(counts_list.t$w.netMig_rectangle, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$w.netMig.UQ[rows.ETTw.t] <- apply(counts_list.t$w.netMig_rectangle, 2, quantile, probs = 0.975, na.rm=T)
    
    # 4.4. Accumulation of lineages
    ETT.withinAreas$AccuLineages[rows.ETTw.t] <- colMeans(counts_list.t$existingLineages, na.rm=T)
    ETT.withinAreas$AccuLineages.LQ[rows.ETTw.t] <- apply(counts_list.t$existingLineages, 2, quantile, probs = 0.025, na.rm=T)
    ETT.withinAreas$AccuLineages.UQ[rows.ETTw.t] <- apply(counts_list.t$existingLineages, 2, quantile, probs = 0.975, na.rm=T)
    
    # 4.5. Accumulation of phylogenetic alpha diversity
    ETT.withinAreas$AccuPhyDiv[rows.ETTw.t] <- colMeans(counts_list.t$existingPhyDiv, na.rm=T)
    ETT.withinAreas$AccuPhyDiv.LQ[rows.ETTw.t] <- apply(counts_list.t$existingPhyDiv, 2, quantile, probs = 0.025)
    ETT.withinAreas$AccuPhyDiv.UQ[rows.ETTw.t] <- apply(counts_list.t$existingPhyDiv, 2, quantile, probs = 0.975)
    
    ETT.withinAreas$branchLengths[rows.ETTw.t] <- colMeans(counts_list.t$branch.lengths.per.area, na.rm=T)
    ETT.withinAreas$branchLengths.LQ[rows.ETTw.t] <- apply(counts_list.t$branch.lengths.per.area, 2, quantile, probs = 0.025)
    ETT.withinAreas$branchLengths.UQ[rows.ETTw.t] <- apply(counts_list.t$branch.lengths.per.area, 2, quantile, probs = 0.975)
    }
    
    # 5. Calculations of global dispersal parameter through time ####
    {ETT.global$disp.rate[t] <- mean(counts_list.t$global.disp.rate, na.rm=T)
    ETT.global$disp.rate.LQ[t] <- quantile(counts_list.t$global.disp.rate, probs = 0.025, na.rm=T)
    ETT.global$disp.rate.UQ[t] <- quantile(counts_list.t$global.disp.rate, probs = 0.975, na.rm=T)
    ETT.global$disp.rate.states[t] <- mean(counts_list.t$global.disp.rate.states, na.rm=T)
    ETT.global$disp.rate.states.LQ[t] <- quantile(counts_list.t$global.disp.rate.states, probs = 0.025, na.rm=T)
    ETT.global$disp.rate.states.UQ[t] <- quantile(counts_list.t$global.disp.rate.states, probs = 0.975, na.rm=T)
    
    Global.Disp.Rates[[t]] <- counts_list.t$global.disp.rate
    Global.Disp.Rates.states[[t]] <- counts_list.t$global.disp.rate.states
    }
    #} # end of else
  } # end of t for loop
  # 6. Calculate dispersal parameter globally between all areas/states for the whole age of the tree ####
  { GLOBAL <- data.frame(d.areas.mean= mean(unlist(Global.Disp.Rates), na.rm=T),
                         d.areas.LQ= quantile(unlist(Global.Disp.Rates), probs = 0.025, na.rm=T),
                         d.areas.UQ= quantile(unlist(Global.Disp.Rates), probs = 0.975, na.rm=T),
                         d.states.mean= mean(unlist(Global.Disp.Rates.states), na.rm=T),
                         d.states.LQ= quantile(unlist(Global.Disp.Rates.states), probs = 0.025, na.rm=T),
                         d.states.UQ= quantile(unlist(Global.Disp.Rates.states), probs = 0.975, na.rm=T) )
  }
  if(results.per.map ==FALSE) {results.per.map2 <- NULL}
  
  if(verbose) cat("\nEvents counted through all the time slices. DONE!\n")
  EventsThroughTime <- list(directional= ETT.directional, beta.div= ETT.beta.div, withinAreas=ETT.withinAreas, global.tt = ETT.global, 
                            global.across.time = GLOBAL, results.per.map = results.per.map2)
  
  return(EventsThroughTime)
}

# # Plot for Net Migration ####
plot_w.netMig <- function(ETT, Area, color, add=F){
  if(!add)  {plot(w.netMig~I(time*-1), data=ETT$withinAreas[which(ETT$withinAreas$to==Area),], type="l", xlab="Million years ago (MYA)", ylab="Net Migration",
                  col=color, lwd=2, ylim=c(-1,1))} else {
                    lines(w.netMig~I(time*-1), data=ETT$withinAreas[which(ETT$withinAreas$to==Area),], type="l", xlab="Million years ago (MYA)", ylab="Net Migration",
                          col=color, lwd=2)
                  }
  xPoints <- c(ETT$withinAreas$time[which(ETT$withinAreas$to == Area)]*-1,
               rev(ETT$withinAreas$time[which(ETT$withinAreas$to == Area)])*-1)

  yPoints <- c(ETT$withinAreas$w.netMig.LQ[which(ETT$withinAreas$to == Area)],
               rev(ETT$withinAreas$w.netMig.UQ[which(ETT$withinAreas$to == Area)]))
  if(any(is.na(yPoints))) yPoints[which(is.na(yPoints))] <- 0
  RGB<- col2rgb(color)
  alpha.color <- rgb(RGB[1]/255,RGB[2]/255,RGB[3]/255,.5)
  polygon(xPoints, yPoints, col=alpha.color, border=NA)
  }

  
##EventCounter ####
####
resDECj$inputs$geogfn <- paste0(resDECj$inputs$geogfn)
resDECj$inputs$trfn <- paste0(resDECj$inputs$trfn)
resDECj$inputs$timesfn <- paste0(resDECj$inputs$timesfn)
resDECj$inputs$dispersal_multipliers_fn <- paste0(resDECj$inputs$dispersal_multipliers_fn)
resDECj$inputs$distsfn <- paste0(resDECj$inputs$distsfn)
resDECj$inputs$areas_allowed_fn <- paste0(resDECj$inputs$areas_allowed_fn)
#resDECj$inputs$states_list <- states_list_0based

ETT_DECj <- EventCounter(BGBresults =  resDECj,
                         BSM =  BSM_output,
                         time.slice.size = 5, oldest.time.slice = 60, results.per.map = T)

save(ETT_DECj, file="ETT_DECj_5.Rdata")
load("ETT_DEC_j_5.Rdata")
  
  
}
par(mfrow=c(1,1))

hex = c("#d73027", "#ffa600", "#ff7c43", "#f95d6a", "#d45087", "#a05195","#4575b4","#4575b4", "#2f4b7c","#003f5c")

# How to aggregate areas ####
# ETT_DECj$results.per.map[[1]]$immigration_counts_rectangle
# ETT_DECj$results.per.map[[1]]$emigration_counts_rectangle
## Adding group category of Sahul/Sunda/Island to the results ####
ETT_DECj$withinAreas <- 
  ETT_DECj$withinAreas %>% mutate(group = case_when(to %in% c("A", "N") ~ "Sahul",
                                                    to %in% c("E", "M", "U", "L") ~ "Island",
                                                    TRUE ~ "Sunda"))
ETT_DECj$directional <- 
  ETT_DECj$directional %>% mutate(groupfrom = case_when(from %in% c("A", "N") ~ "Sahul",
                                                        from %in% c("E", "M", "U", "L") ~ "Island",
                                                        TRUE ~ "Sunda"),
                                  groupto = case_when(to %in% c("A", "N") ~ "Sahul",
                                                      to %in% c("E", "M", "U", "L") ~ "Island",
                                                      TRUE ~ "Sunda"))

#Aggregating 
## weighted net mig aggregated ####
ETT_wNetMig_agg <- as.data.frame(NULL)
for (i in 1:12){
  wNetMig_agg <- 
    as.data.frame(rbind(
      quantile(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,1:2]), probs = c(0.25, 0.75, 0.5)),
      quantile(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,3], probs = c(0.25, 0.75, 0.5)),
      quantile(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,4:6]), probs = c(0.25, 0.75, 0.5)),
      quantile(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,7:10]), probs = c(0.25, 0.75, 0.5))))
  wNetMig_agg$time <- (i*5-5)
  wNetMig_agg$label <- c("Sahul", "EMN", "Wallacea", "Sunda")
  
  ETT_wNetMig_agg <- rbind(ETT_wNetMig_agg, wNetMig_agg)
}

ETT_wNetMig_agg <- as.data.frame(NULL)
for (i in 1:12){
  wNetMig_agg <- 
    as.data.frame(rbind(
    c(t.test(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,1:2]))$conf.int,
      mean(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,1:2]))),
    c(t.test((ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,3]))$conf.int,
      mean(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,3])),
    c(t.test(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,4:6]))$conf.int,
      mean(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,4:6]))),
    c(t.test(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,7:10]))$conf.int,
      mean(rowSums(ETT_DECj$results.per.map[[i]]$w.netMig_rectangle[,7:10])))))
  wNetMig_agg$time <- (i*5-5)
  wNetMig_agg$label <- c("Sahul", "EMN", "Wallacea", "Sunda")
  
  ETT_wNetMig_agg <- rbind(ETT_wNetMig_agg, wNetMig_agg)
}

colnames(ETT_wNetMig_agg) <- c("w.netMig.LCI", "w.netMig.UCI", "w.netMig", "time", "group")
ETT_wNetMig_agg[is.nan(ETT_wNetMig_agg$w.netMig.LCI) | is.nan(ETT_wNetMig_agg$w.netMig.UCI), c("w.netMig.LCI", "w.netMig.UCI")] <- 0

netmig3 <- ggplot(ETT_wNetMig_agg, aes(x=time, y = w.netMig, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = w.netMig.LCI, ymax = w.netMig.UCI), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  ggtitle("a) across three regions") +
  theme_classic()+
  xlab("Time (myr)")+
  ylab("Net Migration") +
  scale_y_continuous(limits=c(-0.5,0.6)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0), limits = c(35,0)) +
  scale_color_manual(values = c("#ff7c43", "#d73027", "#003f5c", "#a05195")) +
  scale_fill_manual(values = c("#ff7c43", "#d73027", "#003f5c", "#a05195")) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

## weighted net mig not agg ####
# ETT_wNetMig_CI <- as.data.frame(NULL)
# for(i in 1:10){
#   for (t in 1:12){
#     wNetMig_CI <- as.data.frame(rbind((t.test((ETT_DECj$results.per.map[[t]]$w.netMig_rectangle[,i]))$conf.int)),
#                                       mean(ETT_DECj$results.per.map[[t]]$w.netMig_rectangle[,i]))
#     wNetMig_CI$time <- (t*5-5)
#     wNetMig_CI$label <- i
#     
#     ETT_wNetMig_CI <- rbind(ETT_wNetMig_CI, wNetMig_CI)
#   }
# }
# 
# mapping <- c("AUS", "NG", "EMN", "MAU", "SUL", "LSU", "BOR", "JAV", "SUM", "SEA")
# ETT_wNetMig_CI$label <- mapping[ETT_wNetMig_CI$label]
# 
# colnames(ETT_wNetMig_CI) <- c("w.netMig.UCI", "w.netMig.LCI", "time", "area")
# ETT_wNetMig_CI$w.netMig.mean <- as.numeric(rownames(ETT_wNetMig_CI))
# rownames(ETT_wNetMig_CI) <- NULL
# ETT_wNetMig_CI[is.nan(ETT_wNetMig_CI$w.netMig.LCI) | is.nan(ETT_wNetMig_CI$w.netMig.UCI), c("w.netMig.LCI", "w.netMig.UCI")] <- 0
# ETT_wNetMig_CI$w.netMig.mean[ETT_wNetMig_CI$w.netMig.mean > 1] <- 0
# 
# ETT_DECj$withinAreas$w.netMig.LCI <- ETT_wNetMig_CI$w.netMig.LCI
# ETT_DECj$withinAreas$w.netMig.UCI <- ETT_wNetMig_CI$w.netMig.UCI
# ETT_DECj$withinAreas$w.netMig.mean <- ETT_wNetMig_CI$w.netMig.mean

ggplot(ETT_DECj$withinAreas, aes(x = time, y = emigration, color = to, fill = to)) +
  geom_line(linewidth = 2) 

netmig10 <- ggplot(ETT_DECj$withinAreas, aes(x = time, y = w.netMig, color = to, fill = to)) +
  geom_line() +
  geom_ribbon(aes(ymin = w.netMig.LQ, ymax = w.netMig.UQ), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_classic()+
  ggtitle("b) across ten areas") +
  xlab("Time (myr)")+
  ylab("Net Migration") +
  scale_y_continuous(limits=c(-1,1)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0), limits = c(35,0)) +
  scale_color_manual(values = hex) +
  scale_fill_manual(values = hex) +
  facet_wrap(~group) +
  theme(panel.grid.minor = element_blank(), legend.position = "none")

grid.arrange(netmig3, netmig10)


### Plotting sp richness, phylogenetic richness, and branch lengths ####
ETT_DECj_sumagg <- aggregate(cbind(AccuLineages, AccuLineages.LQ, AccuLineages.UQ,
                                   AccuPhyDiv, AccuPhyDiv.LQ, AccuPhyDiv.UQ,
                                   branchLengths, branchLengths.LQ, branchLengths.UQ) ~ group + time, data = ETT_DECj$withinAreas, FUN = sum)

# ETT_DECj_3 <- EventCounter(BGBresults =  resDECj,
#                          BSM =  BSM_output,
#                          time.slice.size = 3, oldest.time.slice = 60)


LTT_3_plot <- 
  ggplot(ETT_DECj_sumagg, aes(x = time, y = AccuLineages, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = AccuLineages.LQ, ymax = AccuLineages.UQ), alpha = 0.3) +
  theme_minimal()+
  xlab("Time (mya)")+
  ylab("Species Richness") +
  scale_y_continuous(limits=c(0,1200)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "top")

LTT_10_plot <- 
  ggplot(ETT_DECj$withinAreas, aes(x = time, y = AccuLineages, color = to, fill = to)) +
  geom_line() +
  geom_ribbon(aes(ymin = AccuLineages.LQ, ymax = AccuLineages.UQ), alpha = 0.3) +
  theme_minimal()+
  xlab("Time (mya)")+
  ylab("Species Richness") +
  scale_y_continuous(limits=c(0,450)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "top")

PDT_3_plot <- 
  ggplot(ETT_DECj_sumagg, aes(x = time, y = AccuPhyDiv, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = AccuPhyDiv.LQ, ymax = AccuPhyDiv.UQ), alpha = 0.3) +
  theme_minimal()+
  xlab("Time (mya)")+
  ylab("Phylogenetic Richness") +
  scale_y_continuous(limits=c(0,15000)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "none")

PDT_10_plot <- 
  ggplot(ETT_DECj$withinAreas, aes(x = time, y = AccuPhyDiv, color = to, fill = to)) +
  geom_line() +
  geom_ribbon(aes(ymin = AccuPhyDiv.LQ, ymax = AccuPhyDiv.UQ), alpha = 0.3) +
  theme_minimal()+
  xlab("Time (mya)")+
  ylab("Phylogenetic Richness") +
  scale_y_continuous(limits=c(0,5000)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "none")

BL_3_plot <- 
  ggplot(ETT_DECj_sumagg, aes(x = time, y = branchLengths, color = group, fill = group)) +
  geom_line() +
  geom_ribbon(aes(ymin = branchLengths.LQ, ymax = branchLengths.UQ), alpha = 0.3) +
  theme_minimal()+
  xlab("Time (mya)")+
  ylab("Branch Lengths") +
  scale_y_continuous(limits=c(0,5000)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "none")

BL_10_plot <- 
  ggplot(ETT_DECj$withinAreas, aes(x = time, y = branchLengths, color = to, fill = to)) +
  geom_line() +
  geom_ribbon(aes(ymin = branchLengths.LQ, ymax = branchLengths.UQ), alpha = 0.3) +
  theme_minimal()+
  xlab("Time (mya)")+
  ylab("Branch Lengths") +
  scale_y_continuous(limits=c(0,2000)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "none")

grid.arrange(LTT_3_plot, LTT_10_plot, PDT_3_plot, PDT_10_plot, BL_3_plot, BL_10_plot, ncol = 2)

ETT_DECj_dispagg <- aggregate(cbind(immigration, immigration.LQ, immigration.UQ,
                        emigration, emigration.LQ, emigration.UQ) ~ group + time, data = ETT_DECj$withinAreas, FUN = sum)

ETT_DECj_dispagg$ratio <- ETT_DECj_dispagg$emigration/ETT_DECj_dispagg$immigration
ETT_DECj_dispagg$ratio.LQ <- ETT_DECj_dispagg$emigration.LQ/ETT_DECj_dispagg$immigration.LQ
ETT_DECj_dispagg$ratio.UQ <- ETT_DECj_dispagg$emigration.UQ/ETT_DECj_dispagg$immigration.UQ

ETT_DECj_dispagg$ratio[is.nan(ETT_DECj_dispagg$ratio) | is.infinite(ETT_DECj_dispagg$ratio)] <- 0
ETT_DECj_dispagg$ratio.LQ[is.nan(ETT_DECj_dispagg$ratio.LQ) | is.infinite(ETT_DECj_dispagg$ratio.LQ)] <- 0
ETT_DECj_dispagg$ratio.UQ[is.nan(ETT_DECj_dispagg$ratio.UQ) | is.infinite(ETT_DECj_dispagg$ratio.UQ)] <- 0

ETT_DECj$withinAreas$ratio <- ETT_DECj$withinAreas$emigration/ETT_DECj$withinAreas$immigration
ETT_DECj$withinAreas$ratio.LQ <- ETT_DECj$withinAreas$emigration.LQ/ETT_DECj$withinAreas$immigration.LQ
ETT_DECj$withinAreas$ratio.UQ <- ETT_DECj$withinAreas$emigration.UQ/ETT_DECj$withinAreas$immigration.UQ

ETT_DECj$withinAreas$ratio[is.nan(ETT_DECj$withinAreas$ratio) | is.infinite(ETT_DECj$withinAreas$ratio)] <- 0
ETT_DECj$withinAreas$ratio.LQ[is.nan(ETT_DECj$withinAreas$ratio.LQ) | is.infinite(ETT_DECj$withinAreas$ratio.LQ)] <- 0
ETT_DECj$withinAreas$ratio.UQ[is.nan(ETT_DECj$withinAreas$ratio.UQ) | is.infinite(ETT_DECj$withinAreas$ratio.UQ)] <- 0

### Plotting net disp over time in each areas and each group ####
ggplot(ETT_DECj$withinAreas, aes(x = time, y = (emigration - immigration), color = to, fill = to)) +
  geom_line() +
  geom_ribbon(aes(ymin = (emigration.LQ - immigration.LQ), ymax = emigration.UQ - immigration.UQ), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal()+
  xlab("Time (myr)")+
  ylab("Net lineages (source - sink)") +
  #scale_y_continuous(limits=c(-80,180)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "right")

ggplot(ETT_DECj$withinAreas, aes(x = time, y = (ratio), color = to, fill = to)) +
  geom_line(lwd = 2) +
  geom_ribbon(aes(ymin = ratio.LQ, ymax = ratio.UQ ), alpha = 0.3) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_minimal()+
  xlab("Time (myr)")+
  ylab("Outward:Inward Dispersal Ratio") +
  #scale_y_continuous(limits=c(-80,180)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "right")

ggplot(ETT_DECj_dispagg, aes(x = time, y = (emigration-immigration), color = group, fill = group)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = (emigration.LQ-immigration.LQ), ymax = emigration.UQ-immigration.UQ), alpha = 0.3) +
  geom_hline(yintercept = 0, lty = 2) +
  theme_minimal()+
  xlab("Time (myr)")+
  ylab("Net lineages (source - sink)") +
#  scale_y_continuous(limits=c(-20,20)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "right")

ggplot(ETT_DECj_dispagg, aes(x = time, y = (ratio), color = group, fill = group)) +
  geom_line(lwd = 1) +
  geom_ribbon(aes(ymin = (ratio.LQ), ymax = ratio.UQ), alpha = 0.3) +
  geom_hline(yintercept = 1, lty = 2) +
  theme_minimal()+
  xlab("Time (myr)")+
  ylab("Outward:Inward Dispersal Ratio") +
  #scale_y_continuous(limits=c(0,3)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "right")



# Network diagram for dispersal connectivity ####


# Outward vs inward dispersal ####

# Sympatry ####
sympatryrate <- ggplot(ETT_DECj$withinAreas, aes(x = time, y=sympatric.speciation.rate, group = to, color = to, fill = to)) +
  geom_line() +
  geom_ribbon(aes(ymin = sympatric.speciation.rate.LQ, ymax = sympatric.speciation.rate.UQ), alpha = 0.3) +
  theme_classic()+
  ggtitle("d) Sympatric speciation rate") +
  xlab("Time (myr)")+
  ylab("mean and quartile rate") +
  scale_y_continuous(limits=c(0,0.5), expand = c(0,0)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0), limits = c(35,0)) +
  scale_color_manual(values = hex) +
  scale_fill_manual(values = hex) +
  facet_wrap(~group) +
  theme(legend.position = "none", axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

cladogenesis <- ggplot(ETT_DECj$withinAreas, aes(x=time, group = to, color = to, fill = to)) +
  geom_line(aes(y = log(sympatryAll))) +
  geom_ribbon(aes(ymin = log(sympatryAll.LQ), ymax = log(sympatryAll.UQ)), alpha = 0.3) +
  theme_classic()+
  ggtitle("d) Endemic cladogenesis") +
  xlab("Time (myr)")+
  ylab("mean and quartile event (log)") +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0), limits = c(35,0)) +
  scale_color_manual(values = hex) +
  scale_fill_manual(values = hex) +
  facet_wrap(~group) +
  theme(legend.position = "none", axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))

ETT_DECj$directional %>% group_by(time, from) %>% summarise(alldisp = sum(all.dispersals))

# Events counts over time ####
merge <- merge(ETT_DECj$withinAreas, 
               (ETT_DECj$directional %>% group_by(time, from) %>% summarise(alldisp = sum(all.dispersals)))[, c("time","from","alldisp")],
               by.x = c("time", "to"), by.y = c("time", "from"), all.x = TRUE)
ggplot(merge, aes(x = time)) +
  geom_bar(aes(y = sympatryAll + alldisp), stat = "identity") +
  theme_classic()+
  ggtitle("") +
  xlab("Time (myr)")+
  ylab("total event counts") +
  scale_y_continuous(expand = c(0,0), limits = c(0, 1300)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme(legend.position = "none", axis.title = element_text(size = 12), 
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10))
  

# Sympatry slopes ####
ETT_DECj_symagg <- aggregate(cbind(sympatryAll, sympatryAll.LQ, sympatryAll.UQ, 
                                   AccuLineages, AccuLineages.LQ, AccuLineages.UQ) ~ group + time, 
                             data = ETT_DECj$withinAreas, FUN = sum)

ggplot(ETT_DECj_symagg, aes(x = time, y = (sympatryAll), color = group, fill = group)) +
  geom_line(lwd=1) +
  geom_ribbon(aes(ymin = sympatryAll.LQ, ymax = sympatryAll.UQ, alpha = 0.3)) +
  theme_minimal()+
  xlab("Time (myr)")+
  ylab("Within area diversification") +
  #scale_y_continuous(limits=c(0,3)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "right")

calculate_slope <- function(data) {
  slopes <- diff(data$sympatryAll) / diff(data$time)
  return(slopes)
}

calculate_slope.LQ <- function(data) {
  slopes <- diff(data$sympatryAll.LQ) / diff(data$time)
  return(slopes)
}

calculate_slope.UQ <- function(data) {
  slopes <- diff(data$sympatryAll.UQ) / diff(data$time)
  return(slopes)
}

A.df <- ETT_DECj_symagg %>% 
  filter(group == "Island") 
B.df <- ETT_DECj_symagg %>% 
  filter(group == "Sunda")
C.df <- ETT_DECj_symagg %>% 
  filter(group == "Sahul")
slopes <- as.data.frame(cbind(calculate_slope(A.df), calculate_slope(B.df), calculate_slope(C.df),
                              calculate_slope.LQ(A.df), calculate_slope.LQ(B.df), calculate_slope.LQ(C.df),
                              calculate_slope.UQ(A.df), calculate_slope.UQ(B.df), calculate_slope.UQ(C.df)))
colnames(slopes) <- c("Island", "Sunda", "Sahul", "Island.LQ", "Sunda.LQ", "Sahul.LQ", "Island.UQ", "Sunda.UQ", "Sahul.UQ")
slopes$time <- c(0,5,10,15,20,25,30,35,40,45,50)

ggplot(slopes, aes(x=time)) +
  geom_line(aes(y = -1*Island), col = "blue") +
  #geom_ribbon(aes(ymin = -1*Island.LQ, ymax = -1*Island.UQ), fill = "blue", alpha = 0.3) +
  geom_line(aes(y= -1*Sunda), col = "green") +
  #geom_ribbon(aes(ymin = -1*Sunda.LQ, ymax = -1*Sunda.UQ), fill = "green", alpha = 0.3) +
  geom_line(aes(y= -1*Sahul), col = "red") +
  #geom_ribbon(aes(ymin = -1*Sahul.LQ, ymax = -1*Sahul.UQ), fill = "red", alpha = 0.3) +
  theme_minimal()+
  xlab("Time (myr)")+
  ylab("Rate of change in within area diversification events")+
  #scale_y_continuous(limits=c(0,3)) +
  scale_x_reverse(breaks = c(0,10,20,30,40,50,55), expand = c(0,0)) +
  theme( panel.grid.minor = element_blank(), legend.position = "right")
