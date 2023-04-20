### INFO: Some basics to work with RNAseq
### DATE: 04.10.2021
### AUTHOR: Artem Baranovskii



write_gmt <- function(l.gs, path.out, experiment_name = NA) {
  if (!str_detect(path.out, "\\.gmt")) {
    path.out <- paste0(path.out, ".gmt")
  } 
  if (length(experiment_name) < length(l.gs)) {
    experiment_name <- rep(experiment_name, length(l.gs))
  }
  if (is.null(names(l.gs))) {
    names(l.gs) <- rep(NA, length(l.gs))
  }
  l.gs <- purrr::pmap(list(experiment_name,
                           names(l.gs), 
                           l.gs),
                      ~ c(..1, ..2, ..3) %>% paste0(., collapse = "\t")
  ) %>% 
    as.character()
  readr::write_lines(x = l.gs, file = path.out)
}

read_gmt <- function(path.in) {
  if (!str_detect(path.in, "\\.gmt")) {
    stop("Provide path to the .gmt file")
  } 
  l.gs <- purrr::map(read_lines(path.in),
                     ~ str_split_fixed(.x, "\t", Inf)) %>% 
    purrr::map(., ~ .x[which(str_length(.x) > 0)])
  names(l.gs) <- purrr::map(l.gs, ~ .x[[2]])
  experiment_name <- l.gs[[1]][[1]]
  l.gs <- purrr::map(l.gs, ~ .x[3:length(.x)])
  attr(l.gs, "experiment_name") <- experiment_name
  return(l.gs)
}


plot_ma <- function(de_res, de = TRUE, ..., comps = c("",""), 
                    lg2cts_cutoff = 5.5, lg2fc_cutoff = 1, padj_cutoff = 0.05, y_cutoff = c(-7, 7), 
                    empty = FALSE, custom_palette = c("#1565C0", "#b92b27", "#333333")) {
  #
  de_res <- de_res %>% mutate(de_test = "Not significant", 
                              de_test = ifelse(log2(baseMean) > lg2cts_cutoff & log2FoldChange > lg2fc_cutoff & padj < padj_cutoff,
                                               paste0("Enriched"), 
                                               ifelse(log2(baseMean) > lg2cts_cutoff & log2FoldChange < -lg2fc_cutoff & padj < padj_cutoff,
                                                      paste0("Depleted"), 
                                                      "Not significant")), 
                              de_test = ifelse(is.na(de_test), "Not significant", de_test))
  stat_tab <- de_res %>% janitor::tabyl(de_test)
  #
  if (empty == TRUE) {
    #
    ggobj <- ggplot(data = de_res,
                    aes(x = log2(baseMean), 
                        y = log2FoldChange, 
                        color = de_test, 
                        size = de_test, 
                        alpha = de_test)) + 
      geom_point() + 
      geom_hline(yintercept = c(-lg2fc_cutoff, 0, lg2fc_cutoff), 
                 color = c("grey30", "red", "grey30"), 
                 linetype = c("solid", "dashed", "solid"), 
                 size = c(1, 0.5, 1)) +
      scale_size_manual(values = c(1, 1, 0.5)) + 
      scale_alpha_manual(values = c(0.6, 0.6, 0.4)) +
      scale_color_manual(values = custom_palette) +
      ylim(y_cutoff) +
      labs(title = if ((comps == c("", ""))[1]) "" else paste0(comps[1], " vs ", comps[2])) +
      theme_bw() +
      labs(x = "", y = "") +
      theme(axis.text = element_blank(), legend.position = "none")
  } else {
    #
    ggobj <- ggplot(data = de_res,
                    aes(x = log2(baseMean), 
                        y = log2FoldChange, 
                        color = de_test, 
                        size = de_test, 
                        alpha = de_test)) + 
      geom_point() + 
      geom_hline(yintercept = c(-lg2fc_cutoff, 0, lg2fc_cutoff), 
                 color = c("grey30", "red", "grey30"), 
                 linetype = c("solid", "dashed", "solid"), 
                 size = c(1, 0.5, 1)) +
      geom_text(data = tibble(label = c(paste0(stat_tab[2, 1], ": ", stat_tab[2, 2]), 
                                        paste0(stat_tab[3, 1], ": ", stat_tab[3, 2]), 
                                        paste0(stat_tab[1, 1], ": ", stat_tab[1, 2])), 
                              x_coord = c(log2(max(de_res[["baseMean"]])) - 4, 
                                          log2(max(de_res[["baseMean"]])) - 4.5, 
                                          log2(max(de_res[["baseMean"]])) - 4), 
                              y_coord = c((max(de_res[["log2FoldChange"]]) + 2), 
                                          (max(de_res[["log2FoldChange"]]) + 2) - 1, 
                                          (max(de_res[["log2FoldChange"]]) + 2) - 2), 
                              de_test = c("Enriched", "Not significant", "Depleted")), 
                aes(x = x_coord, 
                    y = y_coord,
                    label = label, 
                    color = de_test),
                size = 4) +
      scale_size_manual(values = c(1, 1, 0.5)) + 
      scale_alpha_manual(values = c(0.6, 0.6, 0.4)) +
      scale_color_manual(values = custom_palette) +
      ylim(y_cutoff) +
      labs(title = if ((comps == c("", ""))[1]) "" else paste0(comps[1], " vs ", comps[2])) +
      theme_bw() +
      theme(legend.position = "none")
  }
  return(ggobj)
}


## Helper funtion to plot ECDFs split by percentiles of tAIs in ORFeome experiment
plot_tai_ecdf <- function(t.orfm.res, seqby, are_filt = FALSE) {
  pal.loc <- c("#82b21e", "#666666", "#006beb")
  
  if (are_filt) {
    t.tmp <- t.orfm.res %>% 
      filter(!str_detect(PARENT_GENE_SYMBOL, "^HB"),
             are <= 2) %>% 
      mutate(tai.m_perc = fact_percentile(tai.m, seq_by = seqby, fact = F)) 
  } else {
    t.tmp <- t.orfm.res %>% 
      filter(!str_detect(PARENT_GENE_SYMBOL, "^HB")) %>% 
      mutate(tai.m_perc = fact_percentile(tai.m, seq_by = seqby, fact = F)) 
  }
  
  
  ks.r <- ks.test(t.tmp[t.tmp$tai.m_perc == paste0(seqby * 100, "%"), ]$lg2FC_NvsS, 
                  t.tmp[t.tmp$tai.m_perc == "100%", ]$lg2FC_NvsS)
  
  ggplot(data = t.tmp %>% mutate(tai.m_perc = "General\npopulation"), 
         aes(x = lg2FC_NvsS)) + 
    geom_vline(xintercept = 0) + 
    geom_hline(yintercept = 0.5) +
    stat_ecdf(size = 0.75, geom = "step", color = "grey50", alpha = 0.5) +
    stat_ecdf(data = t.tmp %>% filter(tai.m_perc %in% c(paste0(seqby * 100, "%"), "100%")),
              aes(x = lg2FC_NvsS, color = tai.m_perc), 
              geom = "point", size = 0.3) +
    geom_text(data = tibble(x = 2, y = 0.25, 
                            txt = paste0("D = ", round(ks.r$statistic, 3), 
                                         "\np-value = ", round(ks.r$p.value, 3))),
              aes(x, y, label = txt)) +
    scale_color_manual(values = alpha(setNames(pal.loc[c(3, 1)], c(paste0(seqby * 100, "%"), "100%")), 0.67)) +
    scale_x_continuous(limits = c(-3.5, 3.5)) +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    labs(color = "Optimality\npercentile", y = "Cumulative fraction", x = "lg2FC Neurite vs Soma") +
    scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
    coord_cartesian(xlim = c(-2.5, 2.5)) +
    #theme_bw(base_size = 14) + 
    theme(legend.position = c(0.25, 0.75))
}


