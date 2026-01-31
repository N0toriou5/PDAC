################################################################################
### FIX GSE71729: Mapping corretto usando 'GENE_NAME'
################################################################################

suppressPackageStartupMessages({
  library(GEOquery)
  library(Biobase)
  library(dplyr)
})

gse_id <- "GSE71729"
list_expr_71729 <- list()
list_surv_71729 <- list()

message(paste0(">>> Elaborazione FIX per: ", gse_id))

# 1. Download
gset_list <- tryCatch(getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = TRUE), error = function(e) NULL)
if(is.null(gset_list)) gset_list <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)
gset <- gset_list[[1]]
pd <- pData(gset)

# 2. Filtro Tumori 
# Usiamo source_name_ch2 
tumor_ids <- rownames(pd)[pd$source_name_ch2 == "Pancreas_Primary"]
message(paste0("   Campioni Tumorali Selezionati: ", length(tumor_ids)))

if(length(tumor_ids) == 0) stop("Errore filtro campioni.")

# 3. MAPPING SONDE (IL FIX È QUI)
message("   Mapping Sonde -> Geni (Usando colonna 'GENE_NAME')...")
tbl <- fData(gset)

# Qui forziamo l'uso di GENE_NAME 
if("GENE_NAME" %in% colnames(tbl)) {
  mapping <- tbl[, c("ID", "GENE_NAME")]
  colnames(mapping) <- c("ID", "Symbol") # Rinominiamo per standardizzare
} else {
  stop("Errore: Colonna GENE_NAME non trovata in fData.")
}

# Pulizia Mapping
mapping$Symbol <- as.character(mapping$Symbol)
mapping <- mapping[mapping$Symbol != "" & !is.na(mapping$Symbol), ]

# Creazione Matrice Espressione
raw_expr <- exprs(gset)
common_probes <- intersect(rownames(raw_expr), mapping$ID)
ex_sub <- raw_expr[common_probes, ]
mapping_sub <- mapping[match(common_probes, mapping$ID), ]

# Aggregazione per media
final_expr <- rowsum(ex_sub, group = mapping_sub$Symbol)
count <- table(mapping_sub$Symbol)
final_expr <- final_expr / as.vector(count[rownames(final_expr)])

# Filtra solo i campioni tumorali nella matrice di espressione
final_expr_tumor <- final_expr[, intersect(colnames(final_expr), tumor_ids)]

message(paste0("   Matrice Espressione Pronta: ", nrow(final_expr_tumor), " geni x ", ncol(final_expr_tumor), " campioni"))

# 4. ESTRAZIONE SURVIVAL
message("   Estrazione Survival...")
pd_tumor <- pd[colnames(final_expr_tumor), ]

# Colonne specifiche per GSE71729
# Usiamo grep per sicurezza sui nomi
c_time <- grep("survival_months", colnames(pd_tumor), ignore.case=TRUE, value=TRUE)
c_stat <- grep("death_event", colnames(pd_tumor), ignore.case=TRUE, value=TRUE)

if(length(c_time) > 0 && length(c_stat) > 0) {
  t_vals <- as.numeric(as.character(pd_tumor[[c_time[1]]]))
  s_vals <- as.numeric(as.character(pd_tumor[[c_stat[1]]]))
  
  surv_df <- data.frame(
    sample = rownames(pd_tumor),
    time = t_vals,
    status = s_vals,
    dataset = gse_id,
    stringsAsFactors = FALSE
  )
  
  # Rimuovi NA
  surv_df <- surv_df[!is.na(surv_df$time) & !is.na(surv_df$status), ]
  
  # 5. SALVATAGGIO NELLA LISTA
  common_final <- intersect(colnames(final_expr_tumor), surv_df$sample)
  
  if(length(common_final) > 0) {
    list_expr_71729[[gse_id]] <- final_expr_tumor[, common_final]
    list_surv_71729[[gse_id]] <- surv_df[surv_df$sample %in% common_final, ]
    
    message(">>> SUCCESSO: GSE71729 aggiunto correttamente.")
    message(paste0("    Campioni finali: ", length(common_final)))
    message(paste0("    Eventi (Morti): ", sum(list_surv_71729[[gse_id]]$status == 1)))
  } else {
    stop("Errore: Nessun campione in comune tra espressione e survival.")
  }
  
} else {
  message("Nomi colonne disponibili:")
  print(colnames(pd_tumor))
  stop("Errore: Colonne 'survival_months' o 'death_event' non trovate.")
}
save(list_expr_71729,list_surv_71729,file='./GSE71729.rda')
