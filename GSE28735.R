################################################################################
### ESTRAZIONE E INTEGRAZIONE: GSE28735
### Input: Microarray Affymetrix
### Filtri: Tissue == "T"
### Survival: Survival_month (Time), Cancer_death (Status 0/1)
################################################################################
################################################################################
### ESTRAZIONE CORRETTA: GSE28735
### Target: survival_month:ch1, cancer_deat:ch1, tissue:ch1
################################################################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
})

gse_id <- "GSE28735"
list_expr_28735 <- list()
list_surv_28735 <- list()

message(paste0(">>> Elaborazione ", gse_id, "..."))

# --- 1. DOWNLOAD ---
gset_list <- tryCatch(getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE), error = function(e) NULL)
# Fallback
if(is.null(gset_list)) gset_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

gset <- gset_list[[1]]
pd <- pData(gset)

message(paste0("   Campioni totali grezzi: ", nrow(pd)))

# --- 2. FILTRO TESSUTO (Tissue:ch1 == "T") ---
# Cerchiamo la colonna che contiene "tissue" e "ch1"
# R trasforma "tissue:ch1" in "tissue.ch1" o simili, usiamo grep per sicurezza
col_tissue <- grep("tissue.*ch1", colnames(pd), ignore.case = TRUE, value = TRUE)

if(length(col_tissue) > 0) {
  message(paste0("   Colonna Tessuto trovata: ", col_tissue[1]))
  
  # Applichiamo filtro: Tieni solo chi ha "T"
  # trimws pulisce eventuali spazi vuoti
  is_tumor <- trimws(as.character(pd[[col_tissue[1]]])) == "T"
  tumor_ids <- rownames(pd)[is_tumor]
  
  if(length(tumor_ids) == 0) stop("Errore: Nessun campione 'T' trovato.")
  
  message(paste0("   Campioni Tumorali selezionati: ", length(tumor_ids)))
  
  # Sottocampionamento
  pd_tumor <- pd[tumor_ids, ]
  
} else {
  # Debug: stampa colonne se fallisce
  print(colnames(pd))
  stop("Errore: Colonna 'tissue:ch1' non trovata.")
}

# --- 3. MAPPING SONDE -> GENE SYMBOLS ---
message("   Mapping Sonde -> Geni...")

# Recupero annotazioni
if (is.null(annotation(gset)) || annotation(gset) == "") {
  tbl <- fData(gset)
} else {
  gpl <- tryCatch(getGEO(annotation(gset), AnnotGPL = TRUE), error = function(e) NULL)
  if (is.null(gpl)) tbl <- fData(gset) else tbl <- Table(gpl)
}

sym_col <- grep("Gene.?Symbol|Symbol|GENE_SYMBOL", names(tbl), ignore.case = TRUE, value = TRUE)[1]

if(!is.na(sym_col)) {
  mapping <- tbl[, c("ID", sym_col)]
  colnames(mapping) <- c("ID", "Symbol")
  mapping$Symbol <- gsub(" ///.*", "", as.character(mapping$Symbol))
  mapping <- mapping[!is.na(mapping$Symbol) & mapping$Symbol != "" & mapping$Symbol != "---", ]
  
  raw_expr <- exprs(gset)
  common_probes <- intersect(rownames(raw_expr), mapping$ID)
  
  ex_sub <- raw_expr[common_probes, ]
  mapping_sub <- mapping[match(common_probes, mapping$ID), ]
  
  final_expr <- rowsum(ex_sub, group = mapping_sub$Symbol)
  count <- table(mapping_sub$Symbol)
  final_expr <- final_expr / as.vector(count[rownames(final_expr)])
  
  # Tieni solo colonne tumorali
  final_expr_tumor <- final_expr[, intersect(colnames(final_expr), tumor_ids)]
  
} else {
  stop("Errore: Colonna Gene Symbol non trovata.")
}

# --- 4. ESTRAZIONE SURVIVAL (Time & Status) ---
message("   Parsing Survival...")

# Cerchiamo le colonne specifiche indicate
# Nota: usiamo "deat" perché hai indicato "cancer_deat" (possibile refuso nel dataset originale)
col_time_name <- grep("survival_month.*ch1", colnames(pd_tumor), ignore.case = TRUE, value = TRUE)
col_stat_name <- grep("cancer_death.*ch1", colnames(pd_tumor), ignore.case = TRUE, value = TRUE)

if(length(col_time_name) > 0 && length(col_stat_name) > 0) {
  
  message(paste0("   Time Col:   ", col_time_name[1]))
  message(paste0("   Status Col: ", col_stat_name[1]))
  
  t_vals <- as.numeric(as.character(pd_tumor[[col_time_name[1]]]))
  s_vals <- as.numeric(as.character(pd_tumor[[col_stat_name[1]]]))
  
  # Creazione DataFrame
  surv_df <- data.frame(
    sample = rownames(pd_tumor),
    time = t_vals,
    status = s_vals,
    dataset = gse_id,
    stringsAsFactors = FALSE
  )
  
  # Rimuovi NA
  surv_df <- surv_df[!is.na(surv_df$time) & !is.na(surv_df$status), ]
  
} else {
  message("Colonne disponibili:")
  print(colnames(pd_tumor))
  stop("Errore: Colonne survival_month:ch1 o cancer_deat:ch1 non trovate.")
}

# --- 5. INTEGRAZIONE ---

common_ids <- intersect(colnames(final_expr_tumor), surv_df$sample)

if(length(common_ids) > 5) {
  
  list_expr_28735[[gse_id]] <- final_expr_tumor[, common_ids]
  list_surv_28735[[gse_id]] <- surv_df[surv_df$sample %in% common_ids, ]
  
  message(paste0("\n>>> SUCCESSO: ", gse_id, " elaborato."))
  message(paste0("    Campioni Finali: ", length(common_ids)))
  message(paste0("    Eventi (Morti): ", sum(list_surv_28735[[gse_id]]$status == 1)))
  
} else {
  stop("Errore intersezione finale vuota.")
}
save(list_expr_28735,list_surv_28735,file='./GSE28735.rda')
# --- COMANDO PER UNIRE ---
# list_expr <- c(list_expr, list_expr_28735)
# list_surv <- c(list_surv, list_surv_28735)