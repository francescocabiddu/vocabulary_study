# load libraries
lib <- c("magrittr", "tidyverse", "janitor")
lapply(lib, require, character.only = TRUE)
rm(lib)


#### log PP ####
# thirds
pp_qu_subtlex_us_fun <- function(df1_var, df2_var) {
  subtlex_us_pp %>%
    filter(Input %in% (df1_var %>% 
                         unlist() %>%
                         str_replace_all("_", ".") %>%
                         str_replace_all("UU", "AH") %>%
                         c(df2_var %>% 
                             unlist() %>%
                             str_replace_all("_", ".") %>%
                             str_replace_all("UU", "AH"))
    )) %>%
    .$unsLBPAV %>%
    quantile(probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
}

pp_qu_subtlex_us <- pp_qu_subtlex_us_fun(mot_uni_on_subtlex_us$phon, chi_uni_on_subtlex_us$phon)

quartiles_pp_sets <- function(df, pp_qu = pp_qu_subtlex_us) {
  df %<>%
    mutate(ph_pr_qu1 = sapply(ph_pr, function(x) {
      length(x[which(x < pp_qu[2])])
    }),
    ph_pr_qu2 = sapply(ph_pr, function(x) {
      length(x[which(x >= pp_qu[2] & x < pp_qu[3])])
    }),
    ph_pr_qu3 = sapply(ph_pr, function(x) {
      length(x[which(x >= pp_qu[3])])
    })
    ) %>%
    group_by(baby) %>%
    mutate(ph_pr_qu1_cum = c(cumsum(ph_pr_qu1)),
           ph_pr_qu2_cum = c(cumsum(ph_pr_qu2)),
           ph_pr_qu3_cum = c(cumsum(ph_pr_qu3))) %>%
    select(baby:ph_pr_qu3_cum) %>%
    (function(x) {
      percentages <- x %>%
        select(ph_pr_qu1_cum:ph_pr_qu3_cum) %>%
        adorn_percentages() %>%
        .[,-1] %>%
        (function(y) {
          colnames(y) <- gsub("_cum$", "_perc", colnames(y))
          y
        }) %>%
        mutate(ph_pr_qu1_perc = ph_pr_qu1_perc*100,
               ph_pr_qu2_perc = ph_pr_qu2_perc*100,
               ph_pr_qu3_perc = ph_pr_qu3_perc*100)
      
      x %>%
        ungroup() %>%
        cbind(percentages)
    })
}

chi_uni_on_subtlex_us %<>%
  select(baby:ph_pr) %>%
  quartiles_pp_sets()
mot_uni_on_subtlex_us %<>%
  select(baby:ph_pr) %>%
  quartiles_pp_sets()
mod_on_subtlex_us %<>%
  select(baby:ph_pr) %>%
  quartiles_pp_sets()

#### log pp - phonemic lengths ####
# pp at phonemic length
pp_at_length <- function(list_dfs, df_labels, var1_qu, var2_qu) {
  for (df in df_labels) {
    list_dfs[[df]] <- list_dfs[[df]] %>%
      select(baby:phon_plu, ph_pr) %>%
      quartiles_pp_sets(pp_qu = pp_qu_subtlex_us_fun(var1_qu, 
                                                 var2_qu))
  }
  
  list_dfs
} 

log_pp_lengths %<>%
  pp_at_length(c("chi_uni_len2", "chi_uni_len3", "chi_uni_len4", "chi_uni_len5", "chi_uni_len6",
                 "mod_len2", "mod_len3", "mod_len4", "mod_len5", "mod_len6"),
               mot_uni_len23456$phon_plu, 
               chi_uni_len23456$phon_plu)

log_pp_lengths %<>%
  pp_at_length(c("chi_uni_len2b",
               "mod_len2b"),
               mot_uni_len2$phon_plu,
               chi_uni_len2$phon_plu)

log_pp_lengths[["mod_len3b"]] %<>%
  .[, c(1:4,22)] 

log_pp_lengths %<>%
  pp_at_length(c("chi_uni_len3b",
                 "mod_len3b"),
               mot_uni_len3$phon_plu, 
               chi_uni_len3$phon_plu)

log_pp_lengths %<>%
  pp_at_length(c("chi_uni_len4b",
                 "mod_len4b"),
               mot_uni_len4$phon_plu, 
               chi_uni_len4$phon_plu)

log_pp_lengths %<>%
  pp_at_length(c("chi_uni_len5b",
                 "mod_len5b"),
               mot_uni_len5$phon_plu, 
               chi_uni_len5$phon_plu)

log_pp_lengths %<>%
  pp_at_length(c("chi_uni_len6b",
                 "mod_len6b"),
               mot_uni_len6$phon_plu, 
               chi_uni_len6$phon_plu)



#### ND ####
nei_qu_fun <- function(df1_var, df2_var) {
  online %>%
    filter(phon %in% c(df1_var %>% unlist(),
                       df2_var  %>% unlist())) %>%
    .$nei %>%
    quantile(probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
}

nei_qu <- nei_qu_fun(mot_nei_on$phon, chi_nei_on$phon)

quartiles_nei_sets <- function(df, nd_qu = nei_qu) {
  df %<>%
    mutate(nei_qu1 = sapply(nei, function(x) {
      length(x[which(x < nd_qu[2])])
    }),
    nei_qu2 = sapply(nei, function(x) {
      length(x[which(x >= nd_qu[2] & x < nd_qu[3])])
    }),
    nei_qu3 = sapply(nei, function(x) {
      length(x[which(x >= nd_qu[3])])
    })
    ) %>%
    group_by(baby) %>%
    mutate(nei_qu1_cum = c(cumsum(nei_qu1)),
           nei_qu2_cum = c(cumsum(nei_qu2)),
           nei_qu3_cum = c(cumsum(nei_qu3))) %>%
    select(baby:nei_qu3_cum) %>%
    (function(x) {
      percentages <- x %>%
        select(nei_qu1_cum:nei_qu3_cum) %>%
        adorn_percentages() %>%
        .[,-1] %>%
        (function(y) {
          colnames(y) <- gsub("_cum$", "_perc", colnames(y))
          y
        }) %>%
        mutate(nei_qu1_perc = nei_qu1_perc*100,
               nei_qu2_perc = nei_qu2_perc*100,
               nei_qu3_perc = nei_qu3_perc*100)
      
      x %>%
        ungroup() %>%
        cbind(percentages)
    })
}

chi_nei_on %<>%
  select(baby:nei) %>%
  quartiles_nei_sets()
mot_nei_on %<>%
  select(baby:nei) %>%
  quartiles_nei_sets()
mod_nei_on %<>%
  select(baby:nei) %>%
  quartiles_nei_sets()



#### ND - phonemic lengths ####
# pp at phonemic length
nd_at_length <- function(list_dfs, df_labels, var1_qu, var2_qu) {
  for (df in df_labels) {
    list_dfs[[df]] <- df %>%
      get() %>%
      select(baby:phon_plu, nei) %>%
      quartiles_nei_sets(nd_qu = nei_qu_fun(var1_qu, 
                                            var2_qu))
  }
  
  list_dfs
} 

nd_lengths <- list()

nd_lengths %<>%
  nd_at_length(c("chi_uni_len2", "chi_uni_len3", "chi_uni_len4", "chi_uni_len5", "chi_uni_len6",
                 "mod_len2", "mod_len3", "mod_len4", "mod_len5", "mod_len6"),
               mot_uni_len23456$phon_plu, 
               chi_uni_len23456$phon_plu)

nd_lengths %<>%
  nd_at_length(c("chi_uni_len2b",
                 "mod_len2b"),
               mot_uni_len2$phon_plu,
               chi_uni_len2$phon_plu)

nd_lengths %<>%
  nd_at_length(c("chi_uni_len3b",
                 "mod_len3b"),
               mot_uni_len3$phon_plu, 
               chi_uni_len3$phon_plu)

nd_lengths %<>%
  nd_at_length(c("chi_uni_len4b",
                 "mod_len4b"),
               mot_uni_len4$phon_plu, 
               chi_uni_len4$phon_plu)

nd_lengths %<>%
  nd_at_length(c("chi_uni_len5b",
                 "mod_len5b"),
               mot_uni_len5$phon_plu, 
               chi_uni_len5$phon_plu)

nd_lengths %<>%
  nd_at_length(c("chi_uni_len6b",
                 "mod_len6b"),
               mot_uni_len6$phon_plu, 
               chi_uni_len6$phon_plu)

#### General frequency ####
iphod_freq <- iphod %>%
  group_by(UnTrn) %>%
  summarise(SFreq = sum(SFreq)) %>% 
  ungroup()

assign_freq <- function(df, var = phon, general = TRUE) {
  var <- enquo(var)
  
  if (general == TRUE) {
    df %>%
      mutate(freq = sapply(!!var, function(x) {
        x %>%
          str_replace_all("_", ".") %>%
          str_replace_all("UU", "AH") %>%
          (function(y) {
            tibble(UnTrn = y) %>%
              left_join(., iphod_freq, "UnTrn") %>%
              .$SFreq
          })
      }))
  } else {
    df %>%
      mutate(freq = sapply(!!var, function(x) {
        x %>%
          (function(y) {
            tibble(UnTrn = y) %>%
              left_join(., mot_freq, "UnTrn") %>%
              .$SFreq
          })
      }))
  }

}

chi_uni_on_freq <- assign_freq(chi_uni_on) %>%
  select(baby:phon, freq)
mot_uni_on_freq <- assign_freq(mot_uni_on) %>%
  select(baby:phon, freq)
mod_on_freq <- assign_freq(mod_on) %>%
  select(baby:phon, freq)

freq_qu_fun <- function(df1_var, df2_var, general = TRUE) {
  if (general == TRUE) {
    iphod_freq %>%
      filter(UnTrn %in% (df1_var %>% 
                           unlist() %>%
                           str_replace_all("_", ".") %>%
                           str_replace_all("UU", "AH") %>%
                           c(df2_var %>% 
                               unlist() %>%
                               str_replace_all("_", ".") %>%
                               str_replace_all("UU", "AH"))
      )) %>%
      .$SFreq %>%
      quantile(probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  } else {
    mot_freq %>%
      filter(UnTrn %in% (df1_var %>% 
                           unlist() %>%
                           c(df2_var %>% 
                               unlist())
      )) %>%
      .$SFreq %>%
      quantile(probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
  }
  
}

freq_qus <- freq_qu_fun(mot_uni_on_freq$phon, chi_uni_on_freq$phon)

quartiles_freq_sets <- function(df, freq_qu = freq_qus) {
  df %<>%
    mutate(freq_qu1 = sapply(freq, function(x) {
      length(x[which(x < freq_qu[2])])
    }),
    freq_qu2 = sapply(freq, function(x) {
      length(x[which(x >= freq_qu[2] & x < freq_qu[3])])
    }),
    freq_qu3 = sapply(freq, function(x) {
      length(x[which(x >= freq_qu[3])])
    })) %>%
    group_by(baby) %>%
    mutate(freq_qu1_cum = c(cumsum(freq_qu1)),
           freq_qu2_cum = c(cumsum(freq_qu2)),
           freq_qu3_cum = c(cumsum(freq_qu3))) %>%
    select(baby:freq_qu3_cum) %>%
    (function(x) {
      percentages <- x %>%
        select(freq_qu1_cum:freq_qu3_cum) %>%
        adorn_percentages() %>%
        .[,-1] %>%
        (function(y) {
          colnames(y) <- gsub("_cum$", "_perc", colnames(y))
          y
        }) %>%
        mutate(freq_qu1_perc = freq_qu1_perc*100,
               freq_qu2_perc = freq_qu2_perc*100,
               freq_qu3_perc = freq_qu3_perc*100)
      
      x %>%
        ungroup() %>%
        cbind(percentages)
    })
  
}

chi_uni_on_freq <- quartiles_freq_sets(chi_uni_on_freq)
mot_uni_on_freq <- quartiles_freq_sets(mot_uni_on_freq)
mod_on_freq <- quartiles_freq_sets(mod_on_freq)

#### General frequency at each phonemic length ####
freq_lengths <- list()

freq_at_length <- function(list_dfs, df_labels, var1_qu, var2_qu) {
  for (df in df_labels) {
    list_dfs[[df]] <- df %>%
      get() %>%
      assign_freq(var = phon_plu) %>%
      quartiles_freq_sets(freq_qu = freq_qu_fun(var1_qu, 
                                            var2_qu))
  }
  
  list_dfs
} 

freq_lengths %<>%
  freq_at_length(c("chi_uni_len2", "chi_uni_len3", "chi_uni_len4", "chi_uni_len5", "chi_uni_len6",
                 "mod_len2", "mod_len3", "mod_len4", "mod_len5", "mod_len6"),
               mot_uni_len23456$phon_plu, 
               chi_uni_len23456$phon_plu)

freq_lengths %<>%
  freq_at_length(c("chi_uni_len2b",
                   "mod_len2b"),
                 mot_uni_len2$phon_plu, 
                 chi_uni_len2$phon_plu)

freq_lengths %<>%
  freq_at_length(c("chi_uni_len3b",
                   "mod_len3b"),
                 mot_uni_len3$phon_plu, 
                 chi_uni_len3$phon_plu)

freq_lengths %<>%
  freq_at_length(c("chi_uni_len4b",
                   "mod_len4b"),
                 mot_uni_len4$phon_plu, 
                 chi_uni_len4$phon_plu)

freq_lengths %<>%
  freq_at_length(c("chi_uni_len5b",
                   "mod_len5b"),
                 mot_uni_len5$phon_plu, 
                 chi_uni_len5$phon_plu)

freq_lengths %<>%
  freq_at_length(c("chi_uni_len6b",
                   "mod_len6b"),
                 mot_uni_len6$phon_plu, 
                 chi_uni_len6$phon_plu)


#### Maternal frequency ####
mot_freq <- "mot_tokens.txt" %>%
  read_tsv() %>%
  rename(UnTrn = phon) %>%
  group_by(UnTrn) %>%
  summarise(SFreq = n()) %>%
  ungroup()

chi_uni_on_mot_freq <- assign_freq(chi_uni_on, general = FALSE) %>%
  select(baby:phon, freq)
mot_uni_on_mot_freq <- assign_freq(mot_uni_on, general = FALSE) %>%
  select(baby:phon, freq)
mod_on_mot_freq <- assign_freq(mod_on, general = FALSE) %>%
  select(baby:phon, freq)

freq_mot_qus <- freq_qu_fun(mot_uni_on_mot_freq$phon, chi_uni_on_mot_freq$phon, general = FALSE)

chi_uni_on_mot_freq <- quartiles_freq_sets(chi_uni_on_mot_freq, freq_qu = freq_mot_qus)
mot_uni_on_mot_freq <- quartiles_freq_sets(mot_uni_on_mot_freq, freq_qu = freq_mot_qus)
mod_on_mot_freq <- quartiles_freq_sets(mod_on_mot_freq, freq_qu = freq_mot_qus)

#### Maternal frequency at each phonemic length ####
freq_mot_lengths <- list()

freq_mot_at_length <- function(list_dfs, df_labels, var1_qu, var2_qu) {
  for (df in df_labels) {
    list_dfs[[df]] <- df %>%
      get() %>%
      assign_freq(var = phon_plu, general = FALSE) %>%
      quartiles_freq_sets(freq_qu = freq_qu_fun(var1_qu, 
                                                var2_qu, general = FALSE))
  }
  
  list_dfs
} 

freq_mot_lengths %<>%
  freq_mot_at_length(c("chi_uni_len2", "chi_uni_len3", "chi_uni_len4", "chi_uni_len5", "chi_uni_len6",
                   "mod_len2", "mod_len3", "mod_len4", "mod_len5", "mod_len6"),
                 mot_uni_len23456$phon_plu, 
                 chi_uni_len23456$phon_plu)

freq_mot_lengths %<>%
  freq_mot_at_length(c("chi_uni_len2b",
                   "mod_len2b"),
                 mot_uni_len2$phon_plu, 
                 chi_uni_len2$phon_plu)

freq_mot_lengths %<>%
  freq_mot_at_length(c("chi_uni_len3b",
                   "mod_len3b"),
                 mot_uni_len3$phon_plu, 
                 chi_uni_len3$phon_plu)

freq_mot_lengths %<>%
  freq_mot_at_length(c("chi_uni_len4b",
                   "mod_len4b"),
                 mot_uni_len4$phon_plu, 
                 chi_uni_len4$phon_plu)

freq_mot_lengths %<>%
  freq_mot_at_length(c("chi_uni_len5b",
                   "mod_len5b"),
                 mot_uni_len5$phon_plu, 
                 chi_uni_len5$phon_plu)

freq_mot_lengths %<>%
  freq_mot_at_length(c("chi_uni_len6b",
                   "mod_len6b"),
                 mot_uni_len6$phon_plu, 
                 chi_uni_len6$phon_plu)

#### ND by general frequency ####
freq_nei_qus <- freq_qu_fun(mot_nei_on$phon, chi_nei_on$phon)

select_freq <- function(df, tertile, qus, value = TRUE) {
  prova <- df %>%
    select(baby:nei) %>%
    assign_freq(general = value) %>%
    mutate(freq_q1 = sapply(freq, function(x) {
      y <- x < qus[2]
      y[is.na(y)] <- FALSE
      y
    }),
    freq_q2 = sapply(freq, function(x) {
      y <- x >= qus[2] & x < qus[3]
      y[is.na(y)] <- FALSE
      y
    }),
    freq_q3 = sapply(freq, function(x) {
      y <- x >= qus[3]
      y[is.na(y)] <- FALSE
      y
    }))
  
  for (i in seq_along(prova$baby)) {
    prova$phon[[i]] <- prova$phon[[i]][prova[[tertile]][[i]]]
    prova$nei[[i]] <- prova$nei[[i]][prova[[tertile]][[i]]]
    prova$freq[[i]] <- prova$freq[[i]][prova[[tertile]][[i]]]
  }
  
  prova %>%
    select(baby:nei)
} 

quartiles_nei_sets <- function(df, qu = nei_qu) {
  df %<>%
    mutate(nei_qu1 = sapply(nei, function(x) {
      length(x[which(x < qu[2])])
    }),
    nei_qu2 = sapply(nei, function(x) {
      length(x[which(x >= qu[2] & x < qu[3])])
    }),
    nei_qu3 = sapply(nei, function(x) {
      length(x[which(x >= qu[3])])
    })
    ) %>%
    group_by(baby) %>%
    mutate(nei_qu1_cum = c(cumsum(nei_qu1)),
           nei_qu2_cum = c(cumsum(nei_qu2)),
           nei_qu3_cum = c(cumsum(nei_qu3))) %>%
    select(baby:nei_qu3_cum) %>%
    (function(x) {
      percentages <- x %>%
        select(nei_qu1_cum:nei_qu3_cum) %>%
        adorn_percentages() %>%
        .[,-1] %>%
        (function(y) {
          colnames(y) <- gsub("_cum$", "_perc", colnames(y))
          y
        }) %>%
        mutate(nei_qu1_perc = nei_qu1_perc*100,
               nei_qu2_perc = nei_qu2_perc*100,
               nei_qu3_perc = nei_qu3_perc*100)
      
      x %>%
        ungroup() %>%
        cbind(percentages)
    })
}

nei_freqs <- list()
nei_freqs[["chi_nei_freq1"]] <- select_freq(chi_nei_on, "freq_q1", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["chi_nei_freq2"]] <- select_freq(chi_nei_on, "freq_q2", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["chi_nei_freq3"]] <- select_freq(chi_nei_on, "freq_q3", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["mot_nei_freq1"]] <- select_freq(mot_nei_on, "freq_q1", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["mot_nei_freq2"]] <- select_freq(mot_nei_on, "freq_q2", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["mot_nei_freq3"]] <- select_freq(mot_nei_on, "freq_q3", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["mod_nei_freq1"]] <- select_freq(mod_nei_on, "freq_q1", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["mod_nei_freq2"]] <- select_freq(mod_nei_on, "freq_q2", freq_nei_qus) %>%
  quartiles_nei_sets()
nei_freqs[["mod_nei_freq3"]] <- select_freq(mod_nei_on, "freq_q3", freq_nei_qus) %>%
  quartiles_nei_sets()

#### ND by maternal frequency ####
freq_nei_mot_qus <- freq_qu_fun(mot_nei_on$phon, chi_nei_on$phon, general = FALSE)

nei_mot_freqs <- list()
nei_mot_freqs[["chi_nei_freq1"]] <- select_freq(chi_nei_on, "freq_q1", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["chi_nei_freq2"]] <- select_freq(chi_nei_on, "freq_q2", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["chi_nei_freq3"]] <- select_freq(chi_nei_on, "freq_q3", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["mot_nei_freq1"]] <- select_freq(mot_nei_on, "freq_q1", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["mot_nei_freq2"]] <- select_freq(mot_nei_on, "freq_q2", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["mot_nei_freq3"]] <- select_freq(mot_nei_on, "freq_q3", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["mod_nei_freq1"]] <- select_freq(mod_nei_on, "freq_q1", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["mod_nei_freq2"]] <- select_freq(mod_nei_on, "freq_q2", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()
nei_mot_freqs[["mod_nei_freq3"]] <- select_freq(mod_nei_on, "freq_q3", 
                                                freq_nei_mot_qus, value = FALSE) %>%
  quartiles_nei_sets()

