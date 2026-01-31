################################################################################
### ESTRAZIONE E PROCESSING RNA-SEQ: GSE79668
### 1. Download Raw Counts
### 2. Conversione Entrez ID -> Gene Symbol
### 3. Normalizzazione log2(CPM + 1)
### 4. Estrazione Survival
################################################################################
suppressPackageStartupMessages({
  library(GEOquery)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(AnnotationDbi)
  library(org.Hs.eg.db) 
})

gse_id <- "GSE79668"
list_expr_79668 <- list()
list_surv_79668 <- list()

message(paste0(">>> Elaborazione RNA-Seq: ", gse_id, "..."))

# ------------------------------------------------------------------------------
# 1. DOWNLOAD RAW COUNTS (Metodo GEO2R)
# ------------------------------------------------------------------------------
message("   Scaricamento e lettura tabella dei conteggi...")

url_path <- paste0("https://www.ncbi.nlm.nih.gov/geo/download/?format=file&type=rnaseq_counts",
                   "&acc=GSE79668",
                   "&file=GSE79668_raw_counts_GRCh38.p13_NCBI.tsv.gz")

# Leggiamo la tabella
# row.names = 1 imposta la prima colonna (EntrezID) come nomi di riga
raw_counts <- as.matrix(data.table::fread(url_path, header=TRUE), rownames=1)

message(paste0("   Dimensioni Raw Counts: ", nrow(raw_counts), " geni x ", ncol(raw_counts), " campioni"))

# Filtro geni poco espressi (come suggerito da GEO2R)
keep <- rowSums(raw_counts >= 10) >= 2
raw_counts <- raw_counts[keep, ]

# ------------------------------------------------------------------------------
# 2. CONVERSIONE ENTREZ ID -> GENE SYMBOL
# ------------------------------------------------------------------------------
message("   Conversione Entrez ID -> Gene Symbol...")

entrez_ids <- rownames(raw_counts)

# Mappatura usando il database org.Hs.eg.db
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = entrez_ids,
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

# Rimuovi geni che non hanno un simbolo mappato (NA)
valid_genes <- !is.na(gene_symbols)
raw_counts_sub <- raw_counts[valid_genes, ]
gene_symbols_sub <- gene_symbols[valid_genes]

# Aggregazione (se ci sono duplicati di simboli, prendiamo la media o somma)
# Per RNA-seq sommare i count dei duplicati è la prassi, ma qui usiamo rowsum per semplicità
# Prima assegniamo i nomi
rownames(raw_counts_sub) <- gene_symbols_sub

# Gestione duplicati (somma dei counts per lo stesso simbolo)
final_counts <- rowsum(raw_counts_sub, group = rownames(raw_counts_sub))

message(paste0("   Geni annotati (Symbols): ", nrow(final_counts)))

# ------------------------------------------------------------------------------
# 3. NORMALIZZAZIONE (log2 CPM)
# ------------------------------------------------------------------------------
message("   Normalizzazione log2(CPM + 1) per compatibilità Microarray...")

# Calcolo CPM manuale
# CPM = (counts / total_counts_in_sample) * 1.000.000
#cpm <- apply(final_counts, 2, function(x) (x / sum(x)) * 1e6)

# Trasformazione Logaritmica
#log_cpm <- log2(cpm + 1)

# Trasformazione VST
library(DESeq2)
log_cpm <- vst(final_counts, blind = TRUE, nsub = 1000, fitType = "parametric")


# ------------------------------------------------------------------------------
# 4. DOWNLOAD METADATI E PARSING SURVIVAL
# ------------------------------------------------------------------------------
message("   Scaricamento Metadati (pData)...")

# Scarichiamo solo i metadati tramite getGEO
gset_meta <- getGEO(gse_id, GSEMatrix = TRUE, AnnotGPL = FALSE)[[1]]
pd <- pData(gset_meta)

message("   Parsing Survival...")

# A. STATUS (Patient survival status)
col_stat_name <- grep("Patient.*survival.*status", colnames(pd), ignore.case=TRUE, value=TRUE)
if(length(col_stat_name) == 0) stop("Colonna status non trovata.")

raw_status <- as.character(pd[[col_stat_name[1]]])
clean_status <- rep(NA, length(raw_status))
clean_status[grep("Alive", raw_status, ignore.case=TRUE)] <- 0
clean_status[grep("Dead", raw_status, ignore.case=TRUE)] <- 1

# B. TIME (Characteristics -> survival time (days))
char_cols <- grep("characteristics", colnames(pd), value=TRUE)
target_col_time <- NULL
for(col in char_cols) {
  if(any(grepl("survival time", pd[[col]], ignore.case=TRUE))) {
    target_col_time <- col; break
  }
}

if(is.null(target_col_time)) stop("Colonna survival time non trovata.")

raw_time_str <- as.character(pd[[target_col_time]])
# Estrazione numero (giorni)
extracted_days <- as.numeric(str_extract(raw_time_str, "(?<=: )\\s*[0-9]+"))
if(all(is.na(extracted_days))) extracted_days <- as.numeric(str_extract(raw_time_str, "[0-9]+"))

# Conversione in Mesi
clean_time_months <- extracted_days / 30.42

# ------------------------------------------------------------------------------
# 5. ALLINEAMENTO E OUTPUT
# ------------------------------------------------------------------------------

surv_df <- data.frame(
  sample = rownames(pd),
  time = clean_time_months,
  status = clean_status,
  dataset = gse_id,
  stringsAsFactors = FALSE
)
surv_df <- surv_df[!is.na(surv_df$time) & !is.na(surv_df$status), ]

# L'espressione RNA-seq ha nomi colonna tipo "GSM2100...", controlliamo l'intersezione
common_ids <- intersect(colnames(log_cpm), surv_df$sample)

if(length(common_ids) > 5) {
  
  # Salvataggio nelle liste
  list_expr_79668[[gse_id]] <- log_cpm[, common_ids]
  list_surv_79668[[gse_id]] <- surv_df[surv_df$sample %in% common_ids, ]
  
  message(paste0("\n>>> SUCCESSO: ", gse_id, " (RNA-seq) elaborato."))
  message(paste0("    Campioni finali: ", length(common_ids)))
  message(paste0("    Eventi (Morti): ", sum(list_surv_79668[[gse_id]]$status == 1)))
  
} else {
  stop("Errore: I nomi dei campioni nel file counts non corrispondono ai metadati.")
}
save(list_expr_79668,list_surv_79668,file='./GSE79668.rda')
