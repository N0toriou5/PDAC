################################################################################
### Extraction and format of PDAC dataset from GEO: GSE21501
################################################################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
  library(stringr)
})

gse_id <- "GSE21501"
message(paste0(">>> Elaborazione ", gse_id, "..."))

# --- 1. DOWNLOAD ---
gset_list <- tryCatch(getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE), error = function(e) NULL)
if(is.null(gset_list)) gset_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)

gset <- gset_list[[1]]
pd <- pData(gset)
fd <- fData(gset)

# --- 2. MAPPING SONDE -> GENI ---
message("   Mapping Sonde -> Geni...")

sym_col <- grep("GeneName|GENE_SYMBOL|Symbol", colnames(fd), ignore.case = TRUE, value = TRUE)[1]

if(is.na(sym_col)) {
  # Fallback if GeneName not found
  sym_col <- grep("ORF", colnames(fd), ignore.case = TRUE, value = TRUE)[1]
}

if(!is.na(sym_col)) {
  message(paste0("   Usando colonna annotazione: ", sym_col))
  
  mapping <- fd[, c("ID", sym_col)]
  colnames(mapping) <- c("ID", "Symbol")
  mapping$Symbol <- as.character(mapping$Symbol)
  mapping <- mapping[!is.na(mapping$Symbol) & mapping$Symbol != "", ]
  
  raw_expr <- exprs(gset)
  common_probes <- intersect(rownames(raw_expr), mapping$ID)
  
  ex_sub <- raw_expr[common_probes, ]
  mapping_sub <- mapping[match(common_probes, mapping$ID), ]
  
  # Aggregazione (Media)
  final_expr <- rowsum(ex_sub, group = mapping_sub$Symbol)
  count <- table(mapping_sub$Symbol)
  final_expr <- final_expr / as.vector(count[rownames(final_expr)])
  
} else {
  stop("Errore: Impossibile trovare la colonna dei Gene Symbols nelle annotazioni.")
}

# --- 3. FILTRO CAMPIONI E ESTRAZIONE SURVIVAL ---
message("   Filtraggio Campioni e Parsing Survival...")

# Identifichiamo le colonne target nei metadati
col_risk   <- grep("characteristics_ch2\\.5", colnames(pd), value = TRUE)
col_time   <- grep("characteristics_ch2($|\\.)", colnames(pd), value = TRUE) # Cerca ch2 esatto (o seguito da punto)
col_event  <- grep("characteristics_ch2\\.1", colnames(pd), value = TRUE)

col_time <- col_time[!col_time %in% c(col_risk, col_event)]
if(length(col_time) > 1) col_time <- col_time[1] 

message(paste0("   Colonna Risk:  ", col_risk))
message(paste0("   Colonna Time:  ", col_time))
message(paste0("   Colonna Event: ", col_event))

if(length(col_risk)>0 && length(col_time)>0 && length(col_event)>0) {
  
  # A. FILTRO RIGHE (Risk Group)
  # Teniamo solo High Risk o Low Risk
  keep_idx <- grep("risk group: (High|Low) Risk", pd[[col_risk]], ignore.case = TRUE)
  pd_sub <- pd[keep_idx, ]
  
  if(nrow(pd_sub) == 0) stop("Nessun campione trovato col filtro 'risk group'.")
  
  # B. ESTRAZIONE DATI
  # Pulizia stringhe: "os time: 24.5" -> 24.5
  # Pulizia stringhe: "os event: 1" -> 1
  
  raw_time <- as.character(pd_sub[[col_time]])
  raw_event <- as.character(pd_sub[[col_event]])
  
  # Regex per estrarre numeri (anche decimali)
  clean_time <- as.numeric(str_extract(raw_time, "[0-9\\.]+"))
  # Regex per estrarre 0 o 1
  clean_event <- as.numeric(str_extract(raw_event, "[0-9]"))
  
  surv_df <- data.frame(
    sample = rownames(pd_sub),
    time = clean_time,
    status = clean_event,
    dataset = gse_id,
    stringsAsFactors = FALSE
  )
  
  # Rimuovi NA
  surv_df <- surv_df[!is.na(surv_df$time) & !is.na(surv_df$status), ]
  
  # --- 4. ALLINEAMENTO E INSERIMENTO NELLE LISTE ---
  
  # Intersezione ID
  common_ids <- intersect(colnames(final_expr), surv_df$sample)
  
  if(length(common_ids) > 5) {
    
    # Crea lista temporanea per questo dataset
    list_expr_21501 <- list()
    list_surv_21501 <- list()
    
    list_expr_21501[[gse_id]] <- final_expr[, common_ids]
    list_surv_21501[[gse_id]] <- surv_df[surv_df$sample %in% common_ids, ]
    
    message(paste0("\n>>> SUCCESSO: ", gse_id, " elaborato."))
    message(paste0("    Campioni Finali: ", length(common_ids)))
    message(paste0("    Eventi (Morti): ", sum(list_surv_21501[[gse_id]]$status == 1)))
    
    # UNIONE ALLE LISTE PRINCIPALI
    list_expr <- c(list_expr, list_expr_21501)
    list_surv <- c(list_surv, list_surv_21501)
    message(">>> Dataset aggiunto alle liste 'list_expr' e 'list_surv'.")
    
  } else {
    stop("Errore: Nessun campione in comune tra espressione e survival.")
  }
  
} else {
  stop("Errore: Una delle colonne specificate non è stata trovata nei metadati.")
}

save(list_expr_21501,list_surv_21501,file='./GSE21501.rda')
