# Function to calculate PET-PEESE
PET.PEESE <- function(data) {
  mod <- list()
  fit_PET <- lm(yi ~ sqrt(vi), 
                weights = 1/vi, 
                data = data)
  
  pet_p <- coef(summary(fit_PET))["(Intercept)", "Pr(>|t|)"] # pet p-value < .10 -> peese
  
  if(pet_p >= .1) {
    mod$b <- coef(summary(fit_PET))["(Intercept)", "Estimate"] # pet estimate
    mod$ci.lb <- confint(fit_PET)["(Intercept)", "2.5 %"] 
    mod$ci.ub<- confint(fit_PET)["(Intercept)", "97.5 %"] 
    mod$pval <- pet_p
    mod$type <- "PET"
    
  }else{
    
    fit_PEESE <- lm(yi ~ vi, 
                    weights = 1/vi, 
                    data = data)
    
    mod$pval <- coef(summary(fit_PEESE))["(Intercept)", "Pr(>|t|)"] # pet p-value < .10 -> peese
    mod$b  <- coef(summary(fit_PEESE))["(Intercept)", "Estimate"] # peese estimate
    mod$ci.lb <- confint(fit_PEESE)["(Intercept)", "2.5 %"] 
    mod$ci.ub <- confint(fit_PEESE)["(Intercept)", "97.5 %"] 
    mod$type <- "PEESE"
    
  }
  return(mod)
}

# Function to create flat violin plots for raincloud plots
geom_flat_violin <- function(mapping = NULL, data = NULL, stat = "ydensity",
                             position = "dodge", trim = TRUE, scale = "width",
                             show.legend = NA, inherit.aes = TRUE, ...) {
  layer(
    data = data,
    mapping = mapping,
    stat = stat,
    geom = GeomFlatViolin,
    position = position,
    show.legend = show.legend,
    inherit.aes = inherit.aes,
    params = list(
      trim = trim,
      scale = scale,
      ...
    )
  )
}

# Function to calculate UWLS
calculate_uwls  <- function(dat) {
  mod <- list()
  d = dat$yi
  sed = sqrt(dat$vi) 
  t = d/sed
  Precision=1/sed
  reg_uwls = lm(t ~ 0 + Precision)
  mod$pval <- coef(summary(reg_uwls))["Precision", "Pr(>|t|)"] # pet p-value < .10 -> peese
  mod$b  <- coef(summary(reg_uwls))["Precision", "Estimate"] # peese estimate
  mod$ci.lb <- confint(reg_uwls)["Precision", "2.5 %"] 
  mod$ci.ub <- confint(reg_uwls)["Precision", "97.5 %"] 
  
  return(mod)
}

# Function to calculate WAAP
calculate_waap <- function(dat) {
  mod <- list()
  d = dat$yi
  sed = sqrt(dat$vi) 
  k = length(d) #number of studies
  t = d/sed
  Precision=1/sed
  reg_uwls = lm(t ~ 0 + Precision)
  UWLS <- as.numeric(reg_uwls$coefficients)
  powered<-sed<abs(UWLS)/2.8
  if(sum(powered)<2) {
    mod$pval  <- NA
    mod$b     <- NA
    mod$ci.lb <- NA
    mod$ci.ub <- NA
  } else{
    reg_waap=lm(t[powered]~ 0 + Precision[powered])
    mod$pval  <- coef(summary(reg_waap))["Precision[powered]", "Pr(>|t|)"] # pet p-value < .10 -> peese
    mod$b     <- coef(summary(reg_waap))["Precision[powered]", "Estimate"] # peese estimate
    mod$ci.lb <- confint(reg_waap)["Precision[powered]", "2.5 %"] 
    mod$ci.ub <- confint(reg_waap)["Precision[powered]", "97.5 %"] 
  }
  return(mod)
}

# Plot Vibration of Effects
plot_VoE <- function(specs, cutoff) {
  
  specs <- specs %>% mutate(p = ifelse(p < 0.0000000001, 0.0000000001, p))
  
  specs$density = densCols(specs$mean, 
                           specs$p, 
                           colramp = colorRampPalette(rev(rainbow(10, end = 4/6))))
  
  k <- nrow(specs)
  p1  <- ggplot(specs,
                aes(x=mean,y=p)) +
    geom_jitter(aes(colour = density), 
                size=1,
                width = .001
    ) +
    scale_color_identity() +
    theme_bw()+
    geom_hline(aes(fill=specs$p),
               yintercept = 0.05, 
               linetype=2)+
    geom_vline(aes(fill=specs$p),
               xintercept =quantile(specs$mean,(.10), na.rm = T),
               color='red', 
               linetype=2)+
    geom_vline(aes(fill=specs$p),
               xintercept =quantile(specs$mean,(.90), na.rm = T),
               color='red', 
               linetype=2)+
    xlab('ES')+ ylab("p-value")
  
  p1f <- p1+
    scale_y_continuous(trans="log",
                       breaks=c(0, 0.0000000001, 0.000000001, 0.00000001, 0.0000001, 0.000001, 0.00001, 0.0001, 0.001,0.01,0.05,0.10,0.50,1),
                       labels=c("<1e-11", "<1e-10","1e-09", "1e-08", "1e-07", "1e-06", "1e-05", "1e-04","0.001","0.01","0.05","0.10","0.50","1"),
                       limits = c(0.0000000001, 1))+
    #limits = c(min(specs$p)-.000000000000000000000000001, 1))+
    scale_x_continuous(breaks = c(-0.5,-0.25,0,0.25,0.5,0.75,1,1.25,1.5,1.75,2),
                       limits = c(-0.7, 2.1))+
    theme(panel.border = element_blank(), 
          plot.title=element_text(size=10)) + 
    ggtitle(paste(k, "meta-analyses with at least", cutoff, "studies."))
  
  print(p1f)
}

# Summarise which and how factors with mean, median, etc.
summarise_factors <- function(data = data, 
                              wf) {
  data %>% 
    group_by({{wf}}) %>% 
    dplyr::summarise(
      Mean = mean(mean),
      LB = mean(lb),
      UB = mean(ub),
      N = n(),
      Median = median(mean),
      Min = min(mean),
      Max = max(mean)) %>% 
    arrange(desc(Mean))
}

### `create_summary_table()` of which factors and filter for subgroups.
create_summary_table <- function(data,
                                 variable_selection,
                                 filter_variable, 
                                 value) {
  
  data  %>%  
    filter({{filter_variable}} == value) %>% 
    select({{variable_selection}}, everything()) %>% 
    group_by({{variable_selection}}) %>% 
    dplyr::summarise(                   
      # "95% CI LB" = mean(lb, na.rm = T),
      "Mean" = mean(mean),
      "Median" = median(mean),
      "Min" = min(mean, na.rm = T),
      "Max" = max(mean, na.rm = T),
      #"95% CI UB" = mean(ub, na.rm = T),
      "k" = n()) %>% 
    arrange(desc(Median)) %>% 
    mutate_if(is.numeric, round, 2)
}

# Creates Rain Cloud Plots of Which Factors
create_rain_plot <- function(data,
                             method,
                             group_variable, 
                             filter_variable,
                             value) {
  set.seed(42)
  if (method == "all") {
    
    data_plot <- data %>% 
      dplyr::select(g = mean, everything()) %>% 
      dplyr::mutate(
        {{group_variable}} := reorder(factor({{group_variable}}), g, median),
        vector = as.integer({{group_variable}}))
    
  } else {
    data_plot <- data %>% 
      dplyr::filter({{filter_variable}} == value) %>% 
      dplyr::select(g = mean, everything()) %>% 
      dplyr::mutate(
        {{group_variable}} := reorder(factor({{group_variable}}), g, median),
        vector = as.integer({{group_variable}}))
  }
  data_plot %>% 
    ggplot(aes(x=reorder({{group_variable}}, g, median),
               y=g, 
               fill = {{group_variable}}, 
               colour = {{group_variable}}))+
    geom_flat_violin(position = position_nudge(x = .25, y = 0),
                     adjust =2, 
                     trim = TRUE,
                     alpha = 0.4,
                     colour = NA) +
    geom_point(position = position_jitter(width = .25), size = .25) +
    geom_boxplot(aes(x = vector + 0.25, 
                     y = g),
                 size = 0.2,
                 alpha = 0.3, width = .1, colour = "BLACK") +
    labs(title = paste(value),
         x = "",
         y = "SMD") +
    coord_flip() +
    theme_cowplot(font_size = 12)+
    guides(fill = FALSE, colour = FALSE) +
    scale_colour_brewer(palette = "Dark2")+
    scale_fill_brewer(palette = "Dark2")
}


# `plot_descriptive_spec_curve()`
## method = "focus" grays out everything except value of variable
## method = "select" creates only a multiverse for the selected value

plot_descriptive_spec_curve  <- function(data = data_specifications,
                                         variable, 
                                         value, 
                                         method = "all") {
  
  if (method == "select") {
    specifications_manual_studies_added <- data %>% 
      filter({{variable}} == value) 
    
  } else if (method =="all") {
    specifications_manual_studies_added <- data
    
  }  else if (method == "focus") {
    
    specifications_manual_studies <- data %>% 
      filter({{variable}} == value) 
    
    specifications_manual_studies_added <- data
  }
  
  wf_1 <- c("website", "mobile", "total_wf_1")
  wf_2 <- c("minimal to no support", 
            "automated encouragement", 
            "human encouragement",
            "guided",
            "total_wf_2")
  wf_3 <- c("other group",
            "med",
            "ppd",
            "young",
            "adul",
            "total_wf_3")
  wf_4 <- c("cbt-based",
            "not-cbt-based",
            "total_wf_4")
  wf_5 <- c("cau", "other ctr", "wl", "total_wf_5")
  
  wf_6 <- c("cut",
            "mood",
            "mdd",
            "total_wf_6")
  wf_7 <- c("include_best", "exclude_worst", "total_wf_7")
  wf_8 <- c("follow-up", "post")
  
  ma_method    <- c("3-level", 
                    "rve", 
                    "reml", 
                    "fe",
                    "p-uniform", 
                    "pet-peese", 
                    "pet-peese (corrected)",
                    "uwls", 
                    "waap")
  
  number_which_how_factors <- 9
  
  x_rank <- rank(specifications_manual_studies_added$mean, 
                 ties.method = "random")
  
  yvar <- rep(factor(rev(c(
    wf_1,
    wf_2,
    wf_3,
    wf_4,
    wf_5,
    wf_6,
    wf_7,
    wf_8,
    ma_method )), 
    levels = rev(c(
      wf_1,
      wf_2,
      wf_3,
      wf_4,
      wf_5,
      wf_6,
      wf_7,
      wf_8,
      ma_method ))), 
    times = nrow(specifications_manual_studies_added))
  
  xvar <- rep(x_rank, 
              each = length(levels(yvar)))
  spec <- NULL
  
  specifications_manual_studies_added <- specifications_manual_studies_added %>%
    select(wf_1:wf_8, ma_method, mean, lb, ub, p, k, set)
  
  for(i in 1:nrow(specifications_manual_studies_added)) {
    id <- as.numeric(levels(yvar) %in% 
                       as.character(unlist(
                         specifications_manual_studies_added[i, 1:number_which_how_factors])))  
    spec <- c(spec, id)
  }
  
  plotdata <- data.frame(xvar, 
                         yvar, 
                         spec)
  
  ylabels <- rev(c(
    # wf_1 <- c("website", "mobile", "total_wf_1")
    "Tech: Website",
    "Tech: Mobile",    
    "Tech: All",  
    
    # wf_2 <- c("minimal to no support", 
    #        "automated encouragement", 
    #        "human encouragement",
    #        "guided",
    #        "total_wf_2")
    
    "Guidance: Minimal to no Support",  
    "Guidance: Automated Encouragement",  
    "Guidance: Human Encouragement", 
    "Guidance: Guided",  
    "Guidance: All",  
    
    # wf_3 <- c("other group",
    #        "med",
    #        "ppd",
    #        "young",
    #        "adul",
    #        "total_wf_3")
    "Group: Other",
    "Group: Medical",
    "Group: PPD",
    "Group: Young",
    "Group: Adults",
    "Group: All",
    
    # wf_4 <- c("cbt-based",
    #        "not-cbt-based",
    #        "total_wf_4")
    
    "Intervention: CBT-Based",
    "Intervention: Not CBT-Based",  
    "Intervention: All", 
    
    # wf_5 <- c("cau", "other ctr", "wl", "total_wf_5")
    
    
    
    "Control: CAU",
    "Control: Other",  
    "Control: Wait List", 
    "Control: All", 
    
    # wf_6 <- c("cut",
    #        "mood",
    #        "mdd",
    #        "total_wf_6")
    "Diagnoses: Cutoff",
    "Diagnoses: Mood Disorder",
    "Diagnoses: MDD",
    "Diagnoses: All",
    
    # wf_7 <- c("exclude_worst", "total_wf_7")
    "ROB: Include Best",  
    "ROB: Excluded Worst",  
    "ROB: All ROB",
    
    # wf_8 <- c("follow-up", "post")
    "Time: Follow Up",
    "Time: Post Intervention",
    
    "Method: 3-Level",
    "Method: RVE",
    "Method: REML",
    "Method: FE",
    "Method: p-uniform",
    "Method: PET-PEESE",
    "Method: PET-PEESE (corrected)",
    "Method: UWLS",
    "Method: WAAP"
  ))
  
  plotdata$k <- rep(specifications_manual_studies_added$k, 
                    each = length(levels(yvar)))  
  
  plotdata$fill <- as.factor(plotdata$k * plotdata$spec)
  
  fill_quantiles <- quantile(plotdata$k, c(.10, .20, .30, .40, .50, .60, .70, .80, .90))
  
  plotdata_rel <- plotdata %>% 
    dplyr::mutate(fill_manual = case_when(
      spec == 0 ~                                      0,           # white
      k  != 0               & k <= fill_quantiles[1] ~ 1, #"#5E4FA2", 
      k > fill_quantiles[1] & k <= fill_quantiles[2] ~ 2, #"#3288BD",
      k > fill_quantiles[2] & k <= fill_quantiles[3] ~ 3, #"#66C2A5",
      k > fill_quantiles[3] & k <= fill_quantiles[4] ~ 4,  #"#ABDDA4",
      k > fill_quantiles[4] & k <= fill_quantiles[5] ~ 5,  #"#E6F598",
      k > fill_quantiles[5] & k <= fill_quantiles[6] ~ 6,  #"#FEE08B",
      k > fill_quantiles[6] & k <= fill_quantiles[7] ~ 7,  #"#FDAE61",
      k > fill_quantiles[7] & k <= fill_quantiles[8] ~ 8,  #"#F46D43",
      k > fill_quantiles[8] & k <= fill_quantiles[9] ~ 9,  #"#D53E4F",
      k > fill_quantiles[9]                          ~ 10), #"#9E0142"),   # dark red
      fill_manual = as.factor(fill_manual))
  
  plotdata_rel$fill_manual %>% levels() 
  # specify the desired order of levels
  
  length_of_each_factor <- c(
    length(ma_method),
    length(wf_8) + length(ma_method),
    length(wf_7) + length(wf_8) + length(ma_method),
    length(wf_6) + length(wf_7) + length(wf_8) + length(ma_method),
    length(wf_5) + length(wf_6) + length(wf_7) + length(wf_8) + length(ma_method),
    length(wf_4) + length(wf_5) + length(wf_6) + length(wf_7) + length(wf_8) + length(ma_method),
    length(wf_3) + length(wf_4) + length(wf_5) + length(wf_6) + length(wf_7) + length(wf_8) + length(ma_method),
    length(wf_2) + length(wf_3) + length(wf_4) + length(wf_5) + length(wf_6) + length(wf_7) + length(wf_8) + length(ma_method),
    length(wf_1) + length(wf_2) + length(wf_3) + length(wf_4) + length(wf_5) + length(wf_6) + length(wf_7) + length(wf_8) + length(ma_method)
  )
  
  if (method == "focus") {
    
    plotdata_grey <- plotdata_rel %>% 
      mutate(id = rep(1:(nrow(plotdata_rel)/length(levels(yvar))),
                      each = length(levels(yvar))),
             fill_manual_num = as.numeric(fill_manual)) %>% 
      group_by(id) %>%
      mutate(
        fill_grey_replace = case_when(
          any(yvar == value & spec == 1) ~ fill_manual_num,
          any(yvar != value & spec == 1) ~ 99),
        fill_grey = ifelse(spec == 0, 0, fill_grey_replace)) %>% 
      ungroup() %>% 
      mutate(fill_manual = as.factor(fill_grey))
    
    plotdata <- plotdata_grey
    
    cols <- RColorBrewer::brewer.pal(min(11, length(levels(plotdata$fill_manual))), "Spectral")
    
    cols <- cols[floor(seq(from = length(cols), to = 0, # change from to to reverse color coding!
                           length.out = length(levels(plotdata$fill_manual))))] # change - 1 for direction #change
    #cols <- c(cols, "grey")
    cols[length(cols)] <- "grey"
    
  } else if (method  %in% c("select", "all")) {
    plotdata <- plotdata_rel
    cols <- RColorBrewer::brewer.pal(min(11, length(levels(plotdata$fill_manual)) - 1), "Spectral")
    
    cols <- cols[floor(seq(from = length(cols), to = 0, # change from to to reverse color coding!
                           length.out = length(levels(plotdata$fill_manual))))] # change - 1 for direction change
  }
  
  tile_plot <- ggplot(data = plotdata, 
                      aes(x = xvar, 
                          y = as.factor(yvar), 
                          fill =  fill_manual)) +
    geom_raster() + 
    geom_hline(yintercept = length_of_each_factor + 0.5) +  # Change lines here here
    scale_x_continuous(position = "bottom") +
    scale_y_discrete(labels = ylabels) +
    scale_fill_manual(
      values = c("white", cols)) + # 
    labs(x = "Specification number", 
         y = "Which/How factors") +
    coord_cartesian(
      expand = F, 
      xlim = c(0, nrow(specifications_manual_studies_added))) +
    theme_bw() + 
    theme(legend.position = "none",
          axis.text.y = element_text(colour = "black", size = 8),
          axis.text.x = element_text(colour = "black"),
          axis.ticks = element_line(colour = "black"),
          plot.margin = margin(t = 5.5, 
                               r = 5.5, 
                               b = 5.5, 
                               l = 5.5, 
                               unit = "pt"))
  
  fill_quantiles <- quantile(specifications_manual_studies_added$k, c(.10, .20, .30, .40, .50, .60, .70, .80, .90))
  
  specifications_rel <- specifications_manual_studies_added %>% 
    mutate(fill_manual = case_when(
      k  != 0               & k <= fill_quantiles[1] ~  1, #"#5E4FA2",  
      k > fill_quantiles[1] & k <= fill_quantiles[2] ~  2, #"#3288BD",
      k > fill_quantiles[2] & k <= fill_quantiles[3] ~  3 , #"#66C2A5",
      k > fill_quantiles[3] & k <= fill_quantiles[4] ~  4, #"#ABDDA4",
      k > fill_quantiles[4] & k <= fill_quantiles[5] ~  5, #"#E6F598",
      k > fill_quantiles[5] & k <= fill_quantiles[6] ~  6, #"#FEE08B",
      k > fill_quantiles[6] & k <= fill_quantiles[7] ~  7, #"#FDAE61",
      k > fill_quantiles[7] & k <= fill_quantiles[8] ~  8, #"#F46D43",
      k > fill_quantiles[8] & k <= fill_quantiles[9]  ~ 9, #"#F46D43",
      k > fill_quantiles[9]                          ~  10), #"#9E0142"),   # dark red
      fill_manual = as.factor(fill_manual))
  
  specifications_rel$xvar <- x_rank
  yrng <- range(c(0, specifications_manual_studies_added$lb, specifications_manual_studies_added$ub))
  ylimit <- c(-.75, 2)
  y_breaks_forest <- seq(-.5, 2, 0.25)
  
  y_labels_forest <- format(y_breaks_forest, nsmall = 2)
  y_breaks_forest <- c(ylimit[1], y_breaks_forest)
  y_labels_forest <- c(ylabels[which.max(nchar(ylabels))], y_labels_forest)
  
  if (method %in% c("select", "focus")) {
    
    specifications_grey <- specifications_rel %>% 
      mutate(fill_manual = case_when(
        {{variable}} != value ~ "grey",
        TRUE ~ as.character(fill_manual)),   
        fill_manual = as.factor(fill_manual))
    
    specifications <- specifications_grey
    
  } else if (method %in% c("all")){
    specifications <- specifications_rel
  }
  
  spec_curve_plot <- specifications %>% 
    ggplot(aes(x = xvar, 
               y = mean))+ 
    geom_errorbar(aes(ymin = lb, 
                      ymax = ub,
                      col = as.factor(fill_manual)), 
                  size = .4,
                  alpha = ifelse(specifications$fill_manual == "grey", 0.3, 1)
    ) +
    geom_line(col = "black", size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", size = 0.7) +
    geom_hline(yintercept = .24, linetype = "dotted", size = 0.7, color = "red") +
    scale_x_continuous(name = "") +
    scale_y_continuous(name = expression(paste("Summary effect (", italic("g"),")")),
                       breaks = y_breaks_forest, labels = y_labels_forest) + 
    scale_color_manual(values = cols) +
    coord_cartesian(ylim = ylimit, 
                    xlim = c(0, nrow(specifications_manual_studies_added)), 
                    expand = FALSE) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y = element_text(colour = c("white", rep("black", times = length(y_labels_forest) - 1))),
          axis.ticks.y = element_line(colour = c("white", rep("black", times = length(y_breaks_forest) - 1))),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          panel.grid.major.y = element_line(),
          panel.grid.minor.y = element_blank(),
          plot.margin = margin(t = 5.5, r = 5.5, b = -15, l = 5.5, unit = "pt"))
  
  
  plot_grid(spec_curve_plot,
            tile_plot,
            labels = ifelse(method == "all", "", paste(str_to_sentence(value))),
            label_x = -0.01,
            label_size = 10,
            ncol = 1,
            align = "v",
            rel_heights = c(4,5))
}