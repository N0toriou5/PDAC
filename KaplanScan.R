library(survival)
library(survminer)

kaplan_scan <- function(expr_mat, surv_data, gene,
                        method = c("scan", "median", "mean"),
                        min_group_size = 0.3,
                        palette = c("black", "red")) {
  
  method <- match.arg(method)
  
  # --- 1. CONTROLLI E PREPARAZIONE DATI ---
  if (!gene %in% rownames(expr_mat)) stop(paste("Gene", gene, "non trovato nella matrice."))
  
  # Standardizziamo il nome colonna: se arriva "patient", lo trattiamo come "sample"
  if ("patient" %in% colnames(surv_data)) {
    colnames(surv_data)[colnames(surv_data) == "patient"] <- "sample"
  }
  
  if (!all(c("sample", "time", "status") %in% colnames(surv_data))) {
    stop("surv_data deve contenere colonne: sample, time, status")
  }
  
  # Estrazione espressione e merge
  expr <- as.numeric(expr_mat[gene, ])
  samples_expr <- colnames(expr_mat)
  df_expr <- data.frame(sample = samples_expr, expr = expr, stringsAsFactors = FALSE)
  
  # Merge mantenendo solo i campioni comuni
  df <- merge(df_expr, surv_data, by = "sample")
  df <- df[complete.cases(df$expr, df$time, df$status), ]
  
  # Pulizia status (assicura numerico 0/1)
  df$status <- as.numeric(as.character(df$status))
  
  if (nrow(df) < 10) stop("Troppi pochi campioni dopo l'allineamento (<10).")
  
  # --- 2. DEFINIZIONE CUTOFF (SCAN vs FIXED) ---
  n <- nrow(df)
  min_n <- max(1, ceiling(n * min_group_size))
  
  cutoff_used <- NA
  best_p <- 1
  
  # A. METODO SCAN (Maxstat)
  if (method == "scan") {
    cutoffs <- sort(unique(df$expr))
    # Escludi estremi per rispettare min_group_size
    cutoffs <- cutoffs[cutoffs > quantile(df$expr, min_group_size) & 
                         cutoffs < quantile(df$expr, 1 - min_group_size)]
    
    if(length(cutoffs) == 0) stop("Nessun cutoff valido rispetta il min_group_size.")
    
    best_cut <- NA
    
    for (c in cutoffs) {
      tmp_group <- ifelse(df$expr >= c, "High", "Low")
      # Calcolo rapido Chi-Quadro
      sd <- survdiff(Surv(time, status) ~ tmp_group, data = df)
      ptmp <- 1 - pchisq(sd$chisq, length(sd$n) - 1)
      
      if (!is.na(ptmp) && ptmp < best_p) {
        best_p <- ptmp
        best_cut <- c
      }
    }
    
    if (is.na(best_cut)) stop("Scansione fallita (nessun p-value calcolabile).")
    cutoff_used <- best_cut
    
    # B. METODO MEDIANA / MEDIA
  } else if (method == "median") {
    cutoff_used <- median(df$expr, na.rm = TRUE)
  } else if (method == "mean") {
    cutoff_used <- mean(df$expr, na.rm = TRUE)
  }
  
  # --- 3. CREAZIONE GRUPPI FINALI ---
  df$Expression <- factor(ifelse(df$expr >= cutoff_used, "High", "Low"),
                          levels = c("Low", "High"))
  
  # Ricalcolo P-value Log-Rank finale per il gruppo scelto
  diff <- survdiff(Surv(time, status) ~ Expression, data = df)
  p_val_logrank <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
  
  # --- 4. CALCOLO COX HAZARD RATIO ---
  cox_model <- coxph(Surv(time, status) ~ Expression, data = df)
  cox_sum <- summary(cox_model)
  
  HR_val  <- round(cox_sum$coefficients[2], 2)
  CI_low  <- round(cox_sum$conf.int[3], 2)
  CI_high <- round(cox_sum$conf.int[4], 2)
  
  # Testo HR formattato
  hr_text <- paste0("HR = ", HR_val, " (", CI_low, " - ", CI_high, ")")
  
  # Testo P-value formattato
  p_label_custom <- paste0("Logrank P = ", format.pval(p_val_logrank, digits = 3, eps = 0.001))
  
  # --- 5. PLOTTING (STILE AGGIORNATO) ---
  fit_km <- survfit(Surv(time, status) ~ Expression, data = df)
  
  # Coordinate per il p-value (basso a sinistra)
  max_time <- max(df$time, na.rm=TRUE)
  coord_x <- max_time * 0.05
  coord_y <- 0.15 
  
  ggs <- ggsurvplot(
    fit_km,
    data = df,
    
    # Gestione P-value custom
    pval = p_label_custom,
    pval.size = 4,
    pval.coord = c(coord_x, coord_y),
    pval.method = FALSE,
    
    # Stile
    palette = palette,
    ggtheme = theme_classic(),
    xlab = "Time (Months)",
    legend.title = "Expression",
    legend.labs = c("Low", "High"),
    
    # Titoli e HR
    risk.table = TRUE,
    title = paste0("Gene: ", gene, " (Cutoff: ", round(cutoff_used,2), ")\n", hr_text),
    tables.height = 0.2
  )
  
  # --- 6. OUTPUT ---
  return(list(
    gene = gene,
    method = method,
    cutoff_value = cutoff_used,
    p_value = p_val_logrank,
    hr_text = hr_text,
    stats = data.frame(HR=HR_val, CI_L=CI_low, CI_H=CI_high, P=p_val_logrank),
    plot = ggs,
    data = df
  ))
}
