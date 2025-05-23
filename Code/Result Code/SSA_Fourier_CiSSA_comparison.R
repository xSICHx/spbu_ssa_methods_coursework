library(Rssa)
library(signal)
library(gsignal)
library(xtable)
source("methods/eossa_new.R")
source("methods/CiSSA.R")
source("methods/time_series_functions.R")

# ------------------------------------
# Зависимость ошибки от длины окна
good_experiment_cissa <- function(
    L_vec,         
    n,             
    x,             
    t_series,      
    groups_freqs,  
    y = NULL,      
    extend_flags = FALSE
) {
  
  if (is.null(y)) {
    y <- Reduce("+", t_series)
  }
  y_real <- Reduce("+", t_series)
  
  K <- length(t_series)
  
  # Определяем имена компонент по names(groups_freqs)
  grp_names <- names(groups_freqs)
  if (is.null(grp_names) || any(grp_names == "")) {
    # если имён нет или есть пустые — дефолтные
    grp_names <- paste0("comp", seq_len(K))
  }
  # Все колонки результата
  all_cols <- c("L", "method", grp_names, "total")
  
  # Заготовка пустого data.frame
  results <- data.frame(matrix(NA, nrow = 0, ncol = length(all_cols)),
                        stringsAsFactors = FALSE)
  colnames(results) <- all_cols
  
  # вспомогательная MSE
  mse <- function(a, b) mean((a - b)^2)
  
  # Главный цикл по L и двум флагам extend
  for (L in L_vec) {
    for (extend_flag in extend_flags) {
      # 1) делаем CiSSA
      cissa_obj <- circulant_SSA(y, L, extend_flag = extend_flag)
      rec_list  <- grouping_cissa(cissa_obj, groups = groups_freqs)$t_series
      
      # 2) считаем MSE для каждой компоненты
      comp_errs <- numeric(K)
      for (i in seq_len(K)) {
        comp_errs[i] <- mse(t_series[[i]], rec_list[[i]])
      }
      # 3) MSE для суммарного
      total_err <- mse(y_real, Reduce("+", rec_list[1:length(rec_list) -1]))
      
      # 4) имя метода
      method_name <- if (extend_flag) "CiSSA_ext" else "CiSSA"
      
      # 5) формируем строку и добавляем в results
      row <- c(L, method_name, comp_errs, total_err)
      results <- rbind(results, setNames(as.list(row), all_cols))
    }
  }
  
  # Приводим типы
  results$L      <- as.integer(results$L)
  results$method <- as.character(results$method)
  for (nm in grp_names) {
    results[[nm]] <- as.numeric(results[[nm]])
  }
  results$total  <- as.numeric(results$total)
  
  rownames(results) <- NULL
  return(results)
}




plot_cissa_combined <- function(df_errors) {
  # Проверяем и загружаем пакеты
  for (pkg in c("ggplot2", "tidyr", "dplyr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf("Установите пакет %s: install.packages('%s')", pkg, pkg))
    }
  }
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  
  # Переводим в длинный формат
  df_long <- df_errors %>%
    pivot_longer(
      cols = -c(L, method),
      names_to  = "series",
      values_to = "mse"
    )
  
  # Получаем все уникальные L
  L_values <- sort(unique(df_long$L))
  
  # Строим график
  p <- ggplot(df_long, aes(x = L, y = mse, color = series, linetype = method)) +
    geom_line(size = 1) +
    geom_point(size = 2) +
    scale_x_continuous(breaks = L_values) +
    labs(
      title = "CiSSA: MSE по всем компонентам и total",
      x     = "Окно L",
      y     = "MSE",
      color = "Серия",
      linetype = "Метод"
    ) +
    theme_minimal() +
    theme(
      plot.title   = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    )
  
  # print(p)
  invisible(p)
}


library(dplyr)

run_cissa_experiments <- function(
    n_experiment,   
    L_vec,          
    n,              
    x,              
    t_series,       
    groups_freqs,   
    extend_flags = c(FALSE, TRUE)  # можно запускать сразу обе опции extend = FALSE и TRUE
) {
  all_runs <- list()
  
  for (run_id in seq_len(n_experiment)) {
    set.seed(1232 + run_id)
    
    x_local <- 0:(n)
    y1 <- sin(2 * pi / 12 * x_local)
    y2 <- 0.5 * cos(2 * pi / 3 * x_local)
    y_noisy <- y1 + y2 + rnorm(n + 1, mean = 0, sd = 0.1)
    
    t_series_clean <- list(y1 = y1, y2 = y2)
    
    df_this_run <- good_experiment_cissa(
      L_vec       = L_vec,
      n           = n,
      x           = x_local,
      t_series    = t_series_clean,
      groups_freqs= groups_freqs,
      y           = y_noisy,
      extend_flags= extend_flags
    )
    
    df_this_run$run_id <- run_id
    
    all_runs[[run_id]] <- df_this_run
  }
  
  # Склеим все результаты по вертикали:
  df_all <- bind_rows(all_runs)
  
  # усредняем
  df_mean <- df_all %>%
    select(-run_id) %>%      
    group_by(L, method) %>%
    summarise_all(mean) %>%  
    ungroup()
  
  return(df_mean)
}


library(dplyr)

run_cissa_experiments_exp <- function(
    n_experiment,     
    L_vec,            
    n,                
    groups_freqs,     
    extend_flags = c(FALSE, TRUE)  # флаги extend = FALSE и TRUE
) {
  all_runs <- list()
  
  for (run_id in seq_len(n_experiment)) {
    
    
    x_local <- 0:(n)
    y1 <- exp(x_local / 200) * sin(2 * pi / 12 * x_local)
    y2 <- exp(x_local / 100) * cos(2 * pi / 3 * x_local)
    
    y_noisy <- y1 + y2 + rnorm(n + 1, mean = 0, sd = 1)
    
    t_series_clean <- list(y1 = y1, y2 = y2)
    
    df_this_run <- good_experiment_cissa(
      L_vec        = L_vec,
      n            = n,
      x            = x_local,
      t_series     = t_series_clean,
      groups_freqs = groups_freqs,
      y            = y_noisy
      # extend_flags = extend_flags
    )
    
    
    df_this_run$run_id <- run_id
    all_runs[[run_id]] <- df_this_run
  }
  
  
  df_all <- bind_rows(all_runs)
  
  df_mean <- df_all %>%
    select(-run_id) %>%
    group_by(L, method) %>%
    summarise_all(mean) %>%
    ungroup()
  
  return(df_mean)
}

#----------------------------------


#---------------------------------
# Функции для численных экспериментов



good_experiment <- function(
    t_series,
    
    noisy_y = NULL,
    
    
    ssa_params = list(),      
    eessa_params = list(),    
    fft_params = list(),         
    fft_ext_params = list(),     
    cissa_params = list(),       
    cissa_ext_params = list(),   
    
    groups_ssa = list(),
    eossa_nested_groups = list(),
    groups_ssa_e = list(),
    groups_fft = list(),
    groups_fft_ext = list(),
    groups_cissa = list(),
    groups_cissa_ext = list()
) {
  if (length(t_series) == 0) {
    stop("t_series не задан.")
  }
  m <- length(t_series)
  
  runs <- list()
  if (length(ssa_params) > 0) {
    for (pr in ssa_params) {
      runs[[length(runs) + 1]] <- list(
        method = "SSA",
        L = pr$L,
        n = pr$n
      )
    }
  }
  if (length(eessa_params) > 0) {
    for (pr in eessa_params) {
      runs[[length(runs) + 1]] <- list(
        method = "E-SSA",
        L = pr$L,
        n = pr$n
      )
    }
  }
  if (length(fft_params) > 0) {
    for (pr in fft_params) {
      runs[[length(runs) + 1]] <- list(
        method = "FFT",
        L = NA,
        n = pr$n
      )
    }
  }
  if (length(fft_ext_params) > 0) {
    for (pr in fft_ext_params) {
      runs[[length(runs) + 1]] <- list(
        method = "FFT-ext",
        L = NA,
        n = pr$n
      )
    }
  }
  if (length(cissa_params) > 0) {
    for (pr in cissa_params) {
      runs[[length(runs) + 1]] <- list(
        method = "CISSA",
        L = pr$L,
        n = pr$n
      )
    }
  }
  if (length(cissa_ext_params) > 0) {
    for (pr in cissa_ext_params) {
      runs[[length(runs) + 1]] <- list(
        method = "CISSA-ext",
        L = pr$L,
        n = pr$n
      )
    }
  }
  if (length(runs) == 0) {
    stop("Ни один параметр для методов не задан.")
  }
  
  all_n <- vapply(runs, function(r) r$n, numeric(1))
  max_n <- max(all_n, na.rm = TRUE)
  
  if (!is.null(noisy_y)) {
    if (length(noisy_y) < max_n) {
      stop("Длина noisy_y меньше максимального n = ", max_n)
    }
    y_full <- noisy_y
  } else {
    y_clean <- Reduce("+", t_series)
    if (length(y_clean) < max_n) {
      stop("Длина суммы t_series меньше максимального n = ", max_n)
    }
    y_full <- y_clean
  }
  
  mse_fun <- function(actual, pred) {
    mean((actual - pred)^2)
  }
  
  results <- vector("list", length(runs))
  
  for (idx in seq_along(runs)) {
    run <- runs[[idx]]
    method <- run$method
    L_val   <- run$L
    n_val   <- run$n
    
    y_partial <- y_full[1:n_val]
    
    rec_list <- vector("list", m)
    
    if (method == "SSA") {
      s_obj <- ssa(y_partial, L_val)
      rec   <- reconstruct(s_obj, groups = groups_ssa)
      for (i in seq_len(m)) {
        rec_list[[i]] <- rec[[i]]
      }
      
    } else if (method == "E-SSA") {
      s_obj <- ssa(y_partial, L_val)
      e_obj <- eossa_new(s_obj, nested.groups = eossa_nested_groups, clust_type = "distance")
      g_e   <- grouping.auto(e_obj, base = "eigen", freq.bins = groups_ssa_e, threshold = 0.5)
      rec   <- reconstruct(e_obj, groups = g_e)
      for (i in seq_len(m)) {
        rec_list[[i]] <- rec[[i]]
      }
      
    } else if (method == "FFT") {
      x_fft   <- seq_len(n_val) - 1
      rec_all  <- reconstruct_fft(x_fft, y_partial)
      rec_grp  <- grouping_cissa(rec_all, groups = groups_fft)$t_series
      for (i in seq_len(m)) {
        rec_list[[i]] <- rec_grp[[i]]
      }
      
    } else if (method == "FFT-ext") {
      x_fft    <- seq_len(n_val) - 1
      rec_all  <- reconstruct_fft(x_fft, y_partial, extend_flag = TRUE)
      rec_grp  <- grouping_cissa(rec_all, groups = groups_fft_ext)$t_series
      for (i in seq_len(m)) {
        rec_list[[i]] <- rec_grp[[i]]
      }
      
    } else if (method == "CISSA") {
      rec_all  <- circulant_SSA(y_partial, L_val)
      rec_grp  <- grouping_cissa(rec_all, groups = groups_cissa)$t_series
      for (i in seq_len(m)) {
        rec_list[[i]] <- rec_grp[[i]]
      }
      
    } else if (method == "CISSA-ext") {
      rec_all  <- circulant_SSA(y_partial, L_val, extend_flag = TRUE)
      rec_grp  <- grouping_cissa(rec_all, groups = groups_cissa_ext)$t_series
      for (i in seq_len(m)) {
        rec_list[[i]] <- rec_grp[[i]]
      }
      
    } else {
      stop("Неизвестный метод: ", method)
    }
    
    series_errors <- numeric(m)
    for (i in seq_len(m)) {
      err <- mse_fun(
        t_series[[i]][1:n_val],
        rec_list[[i]]
      )
      series_errors[i] <- err
    }
    
    rec_sum  <- Reduce("+", rec_list[1:length(t_series)])
    y_clean <- Reduce("+", t_series)[1:n_val]
    total_err <- mse_fun(y_clean, rec_sum)
    
    if (!is.na(L_val)) {
      param_str <- paste0("L=", L_val, ", n=", n_val)
    } else {
      param_str <- paste0("n=", n_val)
    }
    row <- list(
      Method = method,
      Param  = param_str
    )
    for (i in seq_len(m)) {
      row[[paste0("Series", i)]] <- series_errors[i]
    }
    row[["Total"]] <- total_err
    
    results[[idx]] <- row
  }
  
  
  df <- do.call(rbind, lapply(results, function(r) {
    as.data.frame(r, stringsAsFactors = FALSE)
  }))
  
  
  num_cols <- setdiff(names(df), c("Method", "Param"))
  for (col in num_cols) {
    df[[col]] <- as.numeric(df[[col]])
  }
  
  return(df)
}



df_to_latex_print <- function(df) {
  df_fmt <- df
  num_cols <- setdiff(names(df_fmt), c("Method", "Param"))
  df_fmt[num_cols] <- lapply(df_fmt[num_cols], function(col) {
    formatC(col, format = "e", digits = 1)
  })
  table_latex <- xtable(df_fmt, caption = "Example Table")
  print(table_latex, include.rownames = FALSE)
}




run_good_experiments <- function(
    n_experiments,     
    t_series,          
    ssa_params = list(),      
    eessa_params = list(),    
    fft_params = list(),         
    fft_ext_params = list(),     
    cissa_params = list(),       
    cissa_ext_params = list(),   
    groups_ssa = list(),
    eossa_nested_groups = list(),
    groups_ssa_e = list(),
    groups_fft = list(),
    groups_fft_ext = list(),
    groups_cissa = list(),
    groups_cissa_ext = list()
) {
  
  all_n <- c(
    vapply(ssa_params,   function(pr) pr$n, numeric(1)),
    vapply(eessa_params, function(pr) pr$n, numeric(1)),
    vapply(fft_params,   function(pr) pr$n, numeric(1)),
    vapply(fft_ext_params, function(pr) pr$n, numeric(1)),
    vapply(cissa_params, function(pr) pr$n, numeric(1)),
    vapply(cissa_ext_params, function(pr) pr$n, numeric(1))
  )
  max_n <- max(all_n, na.rm = TRUE)
  
  y_clean_full <- Reduce("+", t_series)
  if (length(y_clean_full) < max_n) {
    stop("Длина суммы t_series меньше максимального n = ", max_n)
  }
  
  all_runs <- vector("list", n_experiments)
  
  for (run_id in seq_len(n_experiments)) {
    set.seed(100 + run_id)  
    noise_full <- rnorm(max_n, mean = 0, sd = 0.1)
    
    y_noisy_full <- y_clean_full[1:max_n] + noise_full
    
    df_run <- good_experiment(
      t_series            = t_series,
      noisy_y             = y_noisy_full,
      ssa_params          = ssa_params,
      eessa_params        = eessa_params,
      fft_params          = fft_params,
      fft_ext_params      = fft_ext_params,
      cissa_params        = cissa_params,
      cissa_ext_params    = cissa_ext_params,
      groups_ssa          = groups_ssa,
      eossa_nested_groups = eossa_nested_groups,
      groups_ssa_e        = groups_ssa_e,
      groups_fft          = groups_fft,
      groups_fft_ext      = groups_fft_ext,
      groups_cissa        = groups_cissa,
      groups_cissa_ext    = groups_cissa_ext
    )
    
    df_run$run_id <- run_id
    all_runs[[run_id]] <- df_run
  }
  
  df_all <- bind_rows(all_runs)
  
  num_cols <- setdiff(names(df_all), c("Method", "Param", "run_id"))
  
  df_mean <- df_all %>%
    select(-run_id) %>%
    group_by(Method, Param) %>%
    summarise(across(all_of(num_cols), mean), .groups = "drop")
  
  return(df_mean)
}

