################################################################################
### 0. SETUP E LIBRERIE
################################################################################
setwd('/set/working/directory')

suppressPackageStartupMessages({
  library(GEOquery)
  library(limma)
  library(sva)
  library(survival)
  library(survminer)
  library(dplyr)
  library(tibble)
  library(ggplot2)
  library(tidyr)
})

# Dataset PDAC gestibili insieme
target_gse_ids <- c("GSE183795", "GSE62452", "GSE57495", "GSE85916")

gene_target <- "TGFB1"

################################################################################
### 1. FUNZIONI DI UTILITÀ (PARSING, MAPPING E NORMALIZZAZIONE)
################################################################################

# --- A. Mapping Probes -> Gene Symbols ---
map_probes_to_symbols <- function(gset) {
  gpl_id <- annotation(gset)
  gpl <- tryCatch(getGEO(gpl_id, AnnotGPL = TRUE), error = function(e) NULL)
  
  if (is.null(gpl)) tbl <- fData(gset) else tbl <- Table(gpl)
  
  sym_col <- grep("Gene.?Symbol|Symbol|GENE_SYMBOL", names(tbl), ignore.case = TRUE, value = TRUE)[1]
  if (is.na(sym_col)) return(NULL)
  
  mapping <- tbl[, c("ID", sym_col)]
  colnames(mapping) <- c("ID", "Symbol")
  mapping$Symbol <- gsub(" ///.*", "", as.character(mapping$Symbol))
  mapping <- mapping[mapping$Symbol != "" & !is.na(mapping$Symbol), ]
  
  ex <- exprs(gset)
  common_probes <- intersect(rownames(ex), mapping$ID)
  ex <- ex[common_probes, ]
  mapping <- mapping[match(common_probes, mapping$ID), ]
  
  # Aggregazione per media
  ex_collapsed <- rowsum(ex, group = mapping$Symbol)
  count <- table(mapping$Symbol)
  ex_collapsed <- ex_collapsed / as.vector(count[rownames(ex_collapsed)])
  
  return(ex_collapsed)
}

# --- B. Filtro Tumori (PDAC) ---
filter_tumor_samples <- function(gset, expr_mat) {
  pd <- pData(gset)
  # Keywords specifiche per includere tumori ed escludere normale/adiacente
  keywords_tumor <- "tumor|cancer|carcinoma|pdac|adenocarcinoma|primary"
  keywords_normal <- "normal|adjacent|control|non-tumor|healthy"
  
  text_meta <- paste(pd$title, pd$source_name_ch1, pd$characteristics_ch1)
  is_tumor <- grepl(keywords_tumor, text_meta, ignore.case = TRUE) & 
    !grepl(keywords_normal, text_meta, ignore.case = TRUE)
  
  tumor_ids <- rownames(pd)[is_tumor]
  valid_ids <- intersect(tumor_ids, colnames(expr_mat))
  
  if (length(valid_ids) < 5) {
    warning("Meno di 5 campioni tumorali in ", annotation(gset))
    return(NULL)
  }
  return(expr_mat[, valid_ids])
}

# C. Parsing Survival con NORMALIZZAZIONE MESI 
extract_survival_robust <- function(gset) {
  pd <- pData(gset)
  gse_name <- gset@experimentData@name
  if(is.null(gse_name) || gse_name == "") gse_name <- "Unknown"
  
  clean_num <- function(x) suppressWarnings(as.numeric(as.character(x)))
  
  time_vals <- NULL; status_vals <- NULL; time_unit <- "unknown"
  
  cols <- colnames(pd)
  
  # 1. Cerca colonne "time" e "status" standard
  # Nota: aggiunto "deaths" (specifico per GSE85916) e "os.year" (se fosse colonna diretta)
  t_col <- grep("os\\.time|survival\\.time|overall.*survival|survival.*month|time.*month|futime|survival.*day|time.*day|os\\.year", cols, ignore.case=TRUE, value=TRUE)
  s_col <- grep("os\\.event|survival\\.status|event|status|death|deaths", cols, ignore.case=TRUE, value=TRUE)
  
  # Rimuovi colonne che contengono sia time che status ambiguamente, ma tieni os.year se esiste
  t_col <- t_col[!grepl("status|event|death", t_col, ignore.case=TRUE)]
  
  if(length(t_col) > 0 && length(s_col) > 0) {
    # Priorità alla prima colonna trovata
    t_cand <- clean_num(pd[[t_col[1]]])
    s_raw <- as.character(pd[[s_col[2]]])
    s_cand <- rep(NA, length(s_raw))
    
    # Parsing Status: Gestione robusta 0/1 e Alive/Dead
    s_cand[grep("alive|0|no|cens", s_raw, ignore.case=TRUE)] <- 0
    s_cand[grep("dead|1|yes|event|death", s_raw, ignore.case=TRUE)] <- 1
    
    if (mean(is.na(t_cand)) < 0.5 && mean(is.na(s_cand)) < 0.5) {
      time_vals <- t_cand
      status_vals <- s_cand
      
      # Check label colonna per unità
      if (grepl("day", t_col[1], ignore.case=TRUE)) time_unit <- "days"
      if (grepl("month", t_col[1], ignore.case=TRUE)) time_unit <- "months"
      if (grepl("year", t_col[1], ignore.case=TRUE)) time_unit <- "years"
    }
  }
  
  # 2. Parsing "Characteristics" (Se colonne falliscono o se GSE85916 ha info qui)
  # Nota: GSE85916 ha spesso 'os.year' dentro characteristics_ch1
  if (is.null(time_vals)) {
    char_cols <- grep("characteristics", cols, ignore.case=TRUE, value=TRUE)
    t_vec <- rep(NA, nrow(pd)); s_vec <- rep(NA, nrow(pd))
    
    found_time <- FALSE
    
    for (cc in char_cols) {
      vals <- as.character(pd[[cc]])
      
      #  Parsing Time
      # Caso specifico: os.year (GSE85916)
      if (any(grepl("os\\.year", vals, ignore.case=TRUE))) {
        clean <- gsub(".*os\\.year:\\s*", "", vals, ignore.case=TRUE)
        # Rimuovi eventuale testo successivo se presente (es. split per ;)
        clean <- gsub(";.*", "", clean) 
        nums <- clean_num(clean)
        if (mean(is.na(nums)) < 0.5) {
          t_vec <- nums
          time_unit <- "years" # Forza unità anni
          found_time <- TRUE
        }
      }
      # Caso generico: survival time
      else if (any(grepl("survival.*:|os.*:|month.*:|time.*:", vals, ignore.case=TRUE))) {
        if (any(grepl("day", vals, ignore.case=TRUE))) time_unit <- "days"
        if (any(grepl("year", vals, ignore.case=TRUE))) time_unit <- "years"
        
        clean <- gsub(".*:\\s*", "", vals)
        nums <- clean_num(clean)
        if (mean(is.na(nums)) < 0.5) {
          t_vec <- nums
          found_time <- TRUE
        }
      }
      
      # Parsing Status
      if (any(grepl("status.*:|event.*:|dead|alive", vals, ignore.case=TRUE))) {
        clean <- gsub(".*:\\s*", "", vals)
        s_tmp <- rep(NA, length(clean))
        s_tmp[grep("alive|0|no", clean, ignore.case=TRUE)] <- 0
        s_tmp[grep("dead|1|yes", clean, ignore.case=TRUE)] <- 1
        if (mean(is.na(s_tmp)) < 0.5) s_vec <- s_tmp
      }
    }
    
    # Se abbiamo trovato il tempo in Characteristics ma lo status era in una colonna (es. Deaths)
    # proviamo a recuperare lo status dalle colonne standard se s_vec è vuoto
    if (found_time && all(is.na(s_vec)) && length(s_col) > 0) {
      s_raw <- as.character(pd[[s_col[1]]])
      s_cand <- rep(NA, length(s_raw))
      s_cand[grep("alive|0|no|cens", s_raw, ignore.case=TRUE)] <- 0
      s_cand[grep("dead|1|yes|event|death", s_raw, ignore.case=TRUE)] <- 1
      if (mean(is.na(s_cand)) < 0.5) s_vec <- s_cand
    }
    
    if (!all(is.na(t_vec)) && !all(is.na(s_vec))) {
      time_vals <- t_vec; status_vals <- s_vec
    }
  }
  
  # --- NORMALIZZAZIONE A MESI ---
  if (!is.null(time_vals)) {
    
    # A. Unità ANNI -> Mesi (Nuovo fix per GSE85916)
    if (time_unit == "years") {
      message("   -> Unità rilevata: Anni. Converto in Mesi (x12).")
      time_vals <- time_vals * 12
    }
    # B. Unità GIORNI -> Mesi
    else if (time_unit == "days") {
      message("   -> Unità rilevata: Giorni. Converto in Mesi (/30.42).")
      time_vals <- time_vals / 30.42
    } 
    # C. Euristica: se mediana > 60 e unità ignota -> sono Giorni
    else if (time_unit == "unknown" && median(time_vals, na.rm=TRUE) > 60) {
      message("   -> Mediana temporale alta (", round(median(time_vals, na.rm=TRUE),1), "). Assumo GIORNI e converto in MESI.")
      time_vals <- time_vals / 30.42
    }
    
    return(data.frame(sample=rownames(pd), time=time_vals, status=status_vals, dataset=gse_name))
  }
  return(NULL)
}

################################################################################
### 2. PIPELINE (DOWNLOAD -> PARSE -> INTEGRATE)
################################################################################

list_expr <- list()
list_surv <- list()

for (gse in target_gse_ids) {
  message("\n>>> Processando: ", gse)
  
  # Download
  gset_list <- tryCatch(getGEO(gse, GSEMatrix = TRUE, AnnotGPL = FALSE), error=function(e) NULL)
  if(is.null(gset_list)) { message("Errore download ", gse); next }
  gset <- gset_list[[1]]
  
  # Mappatura
  mat_mapped <- map_probes_to_symbols(gset)
  if(is.null(mat_mapped)) next
  
  # Filtro Tumori
  mat_tumor <- filter_tumor_samples(gset, mat_mapped)
  if(is.null(mat_tumor)) next
  
  # Survival (con nuova logica per GSE85916)
  surv_data <- extract_survival_robust(gset)
  
  if (!is.null(surv_data)) {
    common_s <- intersect(colnames(mat_tumor), surv_data$sample)
    if (length(common_s) > 10) {
      list_expr[[gse]] <- mat_tumor[, common_s]
      list_surv[[gse]] <- surv_data[surv_data$sample %in% common_s, ]
      message("   -> OK. ", length(common_s), " pazienti aggiunti.")
    } else {
      message("   -> Skip: Pochi campioni in comune dopo filtro survival.")
    }
  } else {
    message("   -> Skip: Dati survival assenti o non parsabili.")
  }
}

### INTRODUCE ALSO TCGA DATA ####
packages<-c("corto","DESeq2","ggforce","stringr","matrixStats","TCGAbiolinks","SummarizedExperiment",
            "recount","TCGAutils","biomaRt","limma","factoextra","survival","dbparser","sva","xlsx"
)
for(p in packages){
  if (!p %in% rownames(installed.packages())){
    BiocManager::install(p)
  }
  library(p,character.only=TRUE)
}
source("squisher.R") # load the squisher function, R script is included in this repository (main branch)
source("ensembl2symbol.R")
source("KaplanScan.R")
## find delete_list elements: doi: 10.1158/1078-0432.CCR-18-0290
tcga_ids <- c(
  "TCGA-H6-8124-11", "TCGA-L1-A7W4-01", "TCGA-YB-A89D-11",
  "TCGA-F2-6880-01", "TCGA-F2-7273-01", "TCGA-F2-7276-01",
  "TCGA-H8-A6C1-01", "TCGA-HZ-7920-01", "TCGA-HZ-7923-01",
  "TCGA-H6-A45N-11", "TCGA-HZ-7924-01", "TCGA-IB-AAUV-01",
  "TCGA-IB-AAUW-01", "TCGA-RL-AAAS-01", "TCGA-US-A77J-01",
  "TCGA-HV-A5A3-11", "TCGA-FB-A7DR-01", "TCGA-HZ-8638-01",
  "TCGA-HZ-7289-01", "TCGA-2L-AAQM-01", "TCGA-3A-A9IO-01",
  "TCGA-3A-A9IJ-01", "TCGA-3A-A9IR-01", "TCGA-3A-A9IL-01",
  "TCGA-3A-A9IS-01", "TCGA-3A-A9IN-01", "TCGA-3A-A9IV-01",
  "TCGA-HZ-7918-01", "TCGA-FB-AAPP-01", "TCGA-HV-A7OP-01",
  "TCGA-HZ-A9TJ-06", "TCGA-2J-AABP-01", "TCGA-IB-7654-01"
)

###---- LOAD TCGA PAAD
query <- GDCquery(
  project = "TCGA-PAAD",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

paad_matrix <- assay(data, 'unstranded')


# Crea un'unica espressione regolare con OR (|)
pattern <- paste(tcga_ids, collapse = "|")

# Trova le colonne che contengono almeno uno degli ID
colonne_da_selezionare <- grep(pattern, colnames(paad_matrix), value = TRUE, invert = TRUE)

# Seleziona solo quelle colonne dalla matrice
paad_matrix <- paad_matrix[, colonne_da_selezionare, drop = FALSE]

dim(paad_matrix) # 60660   150

rawcounts <- paad_matrix
# convert ENsembl gene names to HUGO Symbols

ensgenes<-sub('\\..*',"",rownames(rawcounts))
idx<-which(duplicated(ensgenes))
ifelse(length(idx)>0,rawcounts<-rawcounts[-(idx),],rawcounts<-rawcounts)
ensgenes<-sub('\\..*',"",rownames(rawcounts))
idx<-which(duplicated(ensgenes)) # empty
rownames(rawcounts)<-ensgenes

# create a convlist for the squishy thing
symbol_list<-as.character(ens2sym(ensgenes))
names(symbol_list)<-ensgenes
input<-as.matrix(rawcounts[,1:ncol(rawcounts)])
rawmat<-squish(input,symbol_list,method="sum")
#cpm <- apply(rawmat, 2, function(x) (x / sum(x)) * 1e6)

# Trasformazione Logaritmica
#expmat <- log2(cpm + 1)

# VST Transformation
expmat <- vst(rawmat, blind = TRUE, nsub = 1000, fitType = "parametric")

###---- LOAD SURVIVAL DATA
clinical <- GDCquery_clinic("TCGA-PAAD", type = "clinical")

surv_data <- clinical %>%
  dplyr::transmute(
    sample = submitter_id,
    time = ifelse(is.na(days_to_death), days_to_last_follow_up, days_to_death),
    status = ifelse(vital_status == "Dead", 1, 0),
    dataset = 'TCGA'
  ) %>%
  dplyr::filter(!is.na(time))

surv_data$time <- round(surv_data$time/30, digits=0)
colnames(expmat) <- substr(colnames(expmat), 1, 12)
common <- intersect(colnames(expmat), surv_data$sample)

expmat <- expmat[, common]
surv_data <- surv_data[match(common, surv_data$sample), ]

list_expr$'TCGA-PAAD'<- expmat
list_surv$'TCGA-PAAD' <- surv_data

load('./GSE71729.rda')
# # Unione (Assumendo che list_expr e list_surv esistano già con gli altri dataset)
list_expr <- c(list_expr, list_expr_71729)
list_surv <- c(list_surv, list_surv_71729)
load('./GSE21501.rda')
list_expr <- c(list_expr, list_expr_21501)
list_surv <- c(list_surv, list_surv_21501)
load('./GSE79668.rda')
list_expr <- c(list_expr, list_expr_79668)
list_surv <- c(list_surv, list_surv_79668)
#load('./GSE28735.rda')
#list_expr <- c(list_expr, list_expr_28735)
#list_surv <- c(list_surv, list_surv_28735)

################################################################################
### 3. INTEGRAZIONE E BATCH CORRECTION
################################################################################

if (length(list_expr) < 2) stop("Troppi pochi dataset validi (minimo 2).")

common_genes <- Reduce(intersect, lapply(list_expr, rownames))
message("\nIntegrazione su ", length(common_genes), " geni e ", length(list_expr), " studi.")

expr_list_sub <- lapply(list_expr, function(x) x[common_genes, ])
expr_merged <- do.call(cbind, expr_list_sub)

# Batch vector
batch_info <- data.frame(
  sample = colnames(expr_merged),
  batch = rep(names(expr_list_sub), times = sapply(expr_list_sub, ncol))
)

# PCA Pre
pca_pre <- prcomp(t(expr_merged), scale. = TRUE)
p1 <- ggplot(data.frame(pca_pre$x, batch=batch_info$batch), aes(PC1, PC2, color=batch)) +
  geom_point(alpha=0.6) + theme_bw() + ggtitle("PCA - Pre Batch Correction")
print(p1)

# ComBat
message("Esecuzione ComBat...")
colnames(expr_merged) <- make.unique(colnames(expr_merged)) 
expr_corrected <- ComBat(dat = expr_merged, batch = batch_info$batch, par.prior = TRUE)

# PCA Post
pca_post <- prcomp(t(expr_corrected), scale. = TRUE)
p2 <- ggplot(data.frame(pca_post$x, batch=batch_info$batch), aes(PC1, PC2, color=batch)) +
  geom_point(alpha=0.6) + theme_bw() + ggtitle("PCA - Post Batch Correction")
print(p2)

################################################################################
### 4. ANALISI SURVIVAL E MODELLISTICA COX ESTESA
################################################################################

# 1. Preparazione Dati
surv_merged <- do.call(rbind, list_surv)
rownames(surv_merged) <- surv_merged$sample
# Allinea con i campioni corretti (nel caso ComBat abbia cambiato ordine o nomi)
surv_final <- surv_merged[surv_merged$sample %in% colnames(expr_corrected), ]

if (gene_target %in% rownames(expr_corrected)) {
  
  message(">>> Preparazione modelli per: ", gene_target)
  
  # Aggiungi espressione
  expr_vals <- expr_corrected[gene_target, surv_final$sample]
  surv_final$gene_expr <- as.numeric(expr_vals)
  
  # Definisci Gruppi (High/Low) basati sulla Mediana
  med_val <- median(surv_final$gene_expr, na.rm=TRUE)
  surv_final$Expression <- ifelse(surv_final$gene_expr >= med_val, "High", "Low")
  surv_final$Expression <- factor(surv_final$Expression, levels = c("Low", "High"))
  
  # --- A. GRAFICO KAPLAN-MEIER ---
  fit_km <- survfit(Surv(time, status) ~ Expression, data = surv_final)
  
  p_surv <- ggsurvplot(
    fit_km, data = surv_final,
    pval = TRUE, risk.table = TRUE,
    palette = c("dodgerblue", "firebrick"),
    title = paste("Kaplan-Meier Curve -", gene_target, "(Integrated)"),
    xlab = "Time (Months)",
    legend.labs = c("Low Expr", "High Expr"),
    ggtheme = theme_classic(),
    break.time.by = 12
  )
  print(p_surv)
  
  # --- B. MODELLISTICA COX STRATIFICATA ---
  message("\n[Modello Cox] Continuo - STRATIFICATO per Dataset")
  cox_strat <- coxph(Surv(time, status) ~ gene_expr + strata(dataset), data = surv_final)
  print(summary(cox_strat)$coefficients)
  
} else {
  stop("Gene target non trovato.")
}

################################################################################
### 4. ANALISI SURVIVAL E MODELLISTICA COX (Con Hazard Ratio)
################################################################################

# Merge dei dati survival
surv_merged <- do.call(rbind, list_surv)
rownames(surv_merged) <- surv_merged$sample

# Allineamento con dati di espressione
surv_final <- surv_merged[surv_merged$sample %in% colnames(expr_corrected), ]

# Controllo finale status
message("\n>>> Controllo Distribuzione Status Finale:")
print(table(surv_final$status))

if (gene_target %in% rownames(expr_corrected)) {
  
  message(">>> Analisi per Gene: ", gene_target)
  
  # 1. Prepara Dati
  expr_vals <- expr_corrected[gene_target, surv_final$sample]
  surv_final$gene_expr <- as.numeric(expr_vals)
  
  # Cutoff Mediana
  med_val <- median(surv_final$gene_expr, na.rm=TRUE)
  surv_final$Expression <- ifelse(surv_final$gene_expr >= med_val, "High", "Low")
  surv_final$Expression <- factor(surv_final$Expression, levels = c("Low", "High"))
  
  # 2. Calcolo Hazard Ratio (Cox Univariato)
  cox_model <- coxph(Surv(time, status) ~ Expression, data = surv_final)
  cox_res <- summary(cox_model)
  
  HR <- round(cox_res$coefficients[2], 2)        # Exp(coef)
  CI_low <- round(cox_res$conf.int[3], 2)        # Lower .95
  CI_high <- round(cox_res$conf.int[4], 2)       # Upper .95
  P_cox <- signif(cox_res$coefficients[5], 3)    # P-value Cox
  
  hr_text <- paste0("HR = ", HR, " (", CI_low, " - ", CI_high, ")")
  message("   ", hr_text, " | Cox p = ", P_cox)
  
  # 3. Calcolo P-value Log-Rank Esplicito (per personalizzare l'etichetta)
  # Necessario per scrivere "Logrank P" invece del default
  diff <- survdiff(Surv(time, status) ~ Expression, data = surv_final)
  p_val_logrank <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
  p_label_custom <- paste0("log-rank p = ", format.pval(p_val_logrank, digits = 3, eps = 0.001))
  
  # 4. Plot Kaplan Meier Modificato
  fit_km <- survfit(Surv(time, status) ~ Expression, data = surv_final)
  
  # Calcolo posizione ottimale (opzionale, ma utile)
  # Mettiamo il testo a x = 10% del tempo massimo e y = 0.20 (in basso)
  max_time <- max(surv_final$time, na.rm=TRUE)
  coord_x <- 0  # Leggermente staccato dall'asse Y
  coord_y <- 0.10             # In basso (scala 0-1)
  
  p <- ggsurvplot(
    fit_km, 
    data = surv_final,
    
    # --- GESTIONE P-VALUE ---
    pval = p_label_custom,      
    pval.size = 3,              
    pval.method = FALSE,
    
    # Aggiungi questa riga per spostarlo: c(x, y)
    # x = coordinata temporale, y = probabilità (0-1)
    pval.coord = c(coord_x, coord_y), 
    
    # --- ALTRE IMPOSTAZIONI ---
    legend.title = "Expression", 
    risk.table = TRUE,
    title = paste0("Gene: ", gene_target, " - Cutoff: ",round(med_val,2),"\n", hr_text),
    legend.labs = c("Low", "High"),
    palette = c("black", "red"),
    xlab = "Time (Months)",
    ggtheme = theme_classic()
  )
  png(filename = paste0(gene_target, '_KM.png'),w=1500, h=2000,res = 300)
  print(p)
  dev.off()
} else {
  stop("Gene non trovato nella matrice corretta.")
}

################################################################################
### 5. DIFFERENTIAL EXPRESSION ANALYSIS (LIMMA): HIGH vs LOW
################################################################################

message("\n>>> Avvio Analisi Differenziale: ", gene_target, " High vs Low")

# 1. Allineamento rigoroso Campioni
# Assicuriamoci che l'ordine delle righe nei metadati corrisponda esattamente 
# all'ordine delle colonne nella matrice di espressione
surv_final <- surv_final[colnames(expr_corrected), ]

# Verifica (deve essere TRUE)
if(!all(rownames(surv_final) == colnames(expr_corrected))) stop("Errore allineamento campioni!")

# 2. Creazione Design Matrix
# Usiamo la colonna 'Expression' creata in precedenza (High/Low)
# Model matrix senza intercetta (0 + ...) per definire esplicitamente i due gruppi
design <- model.matrix(~ 0 + Expression, data = surv_final)
colnames(design) <- c("Low", "High") # Rinominiamo per chiarezza (levels alfabetici: High, Low -> verifica)
# NOTA: R ordina i livelli alfabeticamente. Se levels(surv_final$Expression) è c("Low", "High"),
# model.matrix produrrà col1=ExpressionLow, col2=ExpressionHigh.
# Per sicurezza, controlliamo i nomi originali e assegniamo correttamente:
colnames(design) <- gsub("Expression", "", colnames(design))

# 3. Fit del Modello Lineare
fit <- lmFit(expr_corrected, design)

# 4. Definizione del Contrasto (High - Low)
# Vogliamo sapere cosa cambia nel gruppo High RISPETTO al gruppo Low
cont.matrix <- makeContrasts(Diff = High - Low, levels = design)

# 5. Calcolo coefficienti e statistica bayesiana
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)

# 6. Estrazione Risultati (Top Table)
# n=Inf per ottenere tutti i geni
de_results <- topTable(fit2, adjust.method = "fdr", number = Inf)

# Anteprima
message("Top 10 Geni differenzialmente espressi:")
print(head(de_results, 10))

# Salvataggio su file
write.csv(de_results, paste0("DE_Results_", gene_target, "_HighVsLow.csv"))

################################################################################
### 6. VISUALIZZAZIONE: VOLCANO PLOT
################################################################################

# Preparazione dati per ggplot
de_results$Symbol <- rownames(de_results)
de_results$DiffExpressed <- "NO"

# Imposta threshold
logFC_cutoff <- 0.58 # circa 1.5 fold change
adjP_cutoff <- 0.05

de_results$DiffExpressed[de_results$logFC > logFC_cutoff & de_results$adj.P.Val < adjP_cutoff] <- "UP"
de_results$DiffExpressed[de_results$logFC < -logFC_cutoff & de_results$adj.P.Val < adjP_cutoff] <- "DOWN"

# Etichette per i top geni (opzionale, i primi 20 per significatività)
de_results$Label <- NA
top_genes <- head(order(de_results$adj.P.Val), 20)
de_results$Label[top_genes] <- de_results$Symbol[top_genes]

# Creazione Plot
p_volcano <- ggplot(de_results, aes(x = logFC, y = -log10(adj.P.Val), col = DiffExpressed, label = Label)) +
  geom_point(alpha = 0.6, size = 1.5) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  geom_vline(xintercept = c(-logFC_cutoff, logFC_cutoff), col = "black", linetype = "dashed") +
  geom_hline(yintercept = -log10(adjP_cutoff), col = "black", linetype = "dashed") +
  theme_bw() +
  labs(title = paste("Volcano Plot:", gene_target, "High vs Low"),
       subtitle = paste("Cutoff: logFC >", logFC_cutoff, "& FDR <", adjP_cutoff),
       x = "log2 Fold Change",
       y = "-log10 FDR") +
  theme(legend.position = "top")

# Aggiungi nomi ai geni (se ggrepel è installato è meglio, altrimenti geom_text standard)
if ("ggrepel" %in% installed.packages()) {
  library(ggrepel)
  p_volcano <- p_volcano + geom_text_repel(max.overlaps = 15, size=3)
} else {
  p_volcano <- p_volcano + geom_text(size=3, check_overlap = TRUE, vjust=1.5)
}

# Salva Volcano
png(filename = paste0(gene_target, '_Volcano.png'), w=2000, h=1800, res = 300)
print(p_volcano)
dev.off()

message(">>> Analisi completata. Salvati CSV risultati e Volcano Plot.")

source('./KaplanScan.R')
# # Assicurati che surv_data abbia le colonne corrette
surv_input <- surv_final[, c("sample", "time", "status")]

# Esegui la scansione
res_scan <- kaplan_scan(
  expr_mat = expr_corrected,
  surv_data = surv_input,
  gene = "TGFBR1",
  method = "mean"  # Oppure "median"
)

# Stampa il grafico
print(res_scan$plot)
