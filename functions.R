run_jags <- function(data, file, em, N, thin){
    jags.data <- list(y.13 = data$d13C,
                      y.14 = data$D14C,
                      age = data$age,
                      n = nrow(data),
                      em13.m = filter(em, Variable == "mean", Isotope == "d13C")$value,
                      em14.m = filter(em, Variable == "mean", Isotope == "D14C")$value,
                      em13.s = filter(em, Variable == "s.d.", Isotope == "d13C")$value,
                      em14.s = filter(em, Variable == "s.d.", Isotope == "D14C")$value)
    jm <- jags.model(file,
                     data = jags.data,
                     n.chains = 1, quiet=TRUE)
    update(jm, 50000)
    coda.samples(jm, c("p", "m.13", "m.14", "mu.13", "mu.14" , "scale.13", "scale.14"),
                 n.iter = N, thin = thin) %>% 
        as.matrix() %>% 
        as_tibble()
}


plot_p <- function(out_young, out_old, data_young, data_old){
    p_y <- get_p(out_young, data_young) %>% 
        group_by(age, source) %>% 
        summarise(mean = mean(value), 
                  median = median(value), 
                  q25 = quantile(value, .25),
                  q75 = quantile(value, .75),
                  q05 = quantile(value, .05),
                  q95 = quantile(value, .95)) %>% 
        filter(age < 4000) %>% 
        mutate(period = "a) Late Holocene")
    p_o <- get_p(out_old, data_old) %>% 
        group_by(age, source) %>% 
        summarise(mean = mean(value), 
                  median = median(value), 
                  q25 = quantile(value, .25),
                  q75 = quantile(value, .75),
                  q05 = quantile(value, .05),
                  q95 = quantile(value, .95)) %>% 
        mutate(period = "c) B/A-YD-HTM")
    ggplot() + 
        geom_ribbon(data = p_o, aes(x = age, ymax = q75, ymin = q25, fill = source), alpha = 0.4) +
        geom_ribbon(data = p_y, aes(x = age, ymax = q75, ymin = q25, fill = source), alpha = 0.4) +
        geom_line(data = p_o, aes(x = age, y = median, color = source)) +
        geom_line(data = p_y, aes(x = age, y = median, color = source)) +
        expand_limits(y = c(0,1)) + 
        ylab("Source fraction") + 
        scale_color_manual(values=c("olivedrab3", "darkorange1", "cadetblue3"), guide = guide_legend(title = NULL)) +
        scale_fill_manual(values=c("olivedrab3", "darkorange1", "cadetblue3"), guide = guide_legend(title = NULL)) + 
        xlab("Cal yr BP") + 
        scale_y_continuous(expand = c(0, 0)) +
        theme_bw() +
        scale_x_continuous(breaks = seq(0, 13000, by = 500)) +
        theme(strip.background = element_blank(), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
              legend.position="none",
              panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
        geom_rug(data = bind_rows(data_young, data_old), aes(x = age)) 
}

plot_p2 <- function(out_young, out_old, data_young, data_old){
  bind_rows(get_p(out_young, data_young), 
            get_p(out_old, data_old)) %>% 
    filter(age %in% filter(bind_rows(data_young, data_old), !is.na(D14C))$age) %>% 
    mutate(age = round(age)) %>% 
    mutate(label = fct_reorder(paste(age, "BP"), age)) %>% 
    ggplot(aes(x=value, fill = source, color = source))+geom_density(alpha = .4) + 
    facet_wrap(~ label) +
    scale_fill_manual(values=c("olivedrab3", "darkorange1", "cadetblue3"), 
                      guide = guide_legend(title = NULL)) +
    scale_color_manual(values=c("olivedrab3", "darkorange1", "cadetblue3"), 
                       guide = guide_legend(title = NULL)) + 
    theme_bw() +
    theme(legend.position=c(0.7, 0.05), axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.text.y=element_blank()) +
    ylab("") + xlab("")
}


get_p <- function(coda_out, data){
    coda_out %>% 
        select(starts_with("p[")) %>% 
        gather() %>% 
        mutate(source = c("ICD", "Active Layer", "Marine OC")[as.numeric(str_sub(key, -2, -2))],
               age = data$age[as.numeric(str_sub(key, 3, -4))])
}



