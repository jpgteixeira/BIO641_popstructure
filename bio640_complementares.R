# =============================================================================
# BIO641 - Genética de Populações
# Estrutura populacional dos painéis CornFed Dent e Flint (GEO GSE50558)
#
# Análises complementares ao STRUCTURE
#   1. Códigos de genótipo (AA, AB, BB, sem chamada) por tipo de linha
#   2. Funil do controle de qualidade (números dos logs do PLINK)
#   3. Distribuição da MAF
#   4. Decaimento do desequilíbrio de ligação (LD) com a distância
#   5. PCA (análise de componentes principais)
#   6. Estatísticas F de Wright (F_IS, F_ST, F_IT) e F_ST de Weir & Cockerham
#
# Pré-requisito: seções 1-8 do bio640_popstructure.R (dados após o CQ e
# conjuntos podados por LD). Saída: results/complementares/
#
# Rodar a partir da raiz do projeto (BIO641_popstructure/).
# =============================================================================

library(data.table)
source("funcoes_popstructure.R")

dir_proc   <- file.path("data", "processed")
dir_str    <- file.path("results", "structure")
dir_saida  <- file.path("results", "complementares")
dir.create(dir_saida, showWarnings = FALSE, recursive = TRUE)

arq_series <- file.path("data", "GSE50558_series_matrix.txt.gz")
prefixo_cq <- file.path(dir_proc, "cornfed_cq")

amostras <- fread(file.path(dir_proc, "amostras.tsv"))
fam_cq   <- fread(paste0(prefixo_cq, ".fam"), header = FALSE)
bim_cq   <- fread(paste0(prefixo_cq, ".bim"), header = FALSE)
amostras[, passou_cq := id %in% fam_cq$V2]

# Grupo de cada amostra para as contagens: a família (DH) ou o pool do
# parental. Os demais níveis (pool, tipo, total) são somas desses grupos.
amostras[, grupo := fifelse(tipo == "DH", familia, paste0("Par_", pool))]

# -----------------------------------------------------------------------------
# 1. Contagem dos genótipos por SNP e por grupo
# -----------------------------------------------------------------------------
# O PLINK recebeu os heterozigotos (AB) como dado faltante (seção 4 do script
# principal). Para as estatísticas F precisamos deles, então voltamos ao
# Series Matrix. A matriz completa (56 mil SNPs x 2.290 amostras, em texto)
# ocupa vários GB; lemos 250 colunas por vez e guardamos só as contagens de
# AA, AB e BB de cada grupo.

# (a última linha, "!series_matrix_table_end", marca o fim da tabela)
ler_bloco <- function(colunas) {
  g <- suppressWarnings(fread(arq_series, skip = "\"ID_REF\"", header = TRUE,
                              quote = "\"", select = c("ID_REF", colunas)))
  g[ID_REF != "!series_matrix_table_end"]
}
ids_snp <- ler_bloco(character(0))$ID_REF

grupos <- sort(unique(amostras$grupo))
contar <- function() {
  matrix(0L, nrow = length(ids_snp), ncol = length(grupos),
         dimnames = list(NULL, grupos))
}
n_AA <- contar(); n_AB <- contar(); n_BB <- contar()
codigos_tipo <- list()

blocos <- split(amostras$gsm, ceiling(seq_len(nrow(amostras)) / 250))
for (b in blocos) {
  g <- ler_bloco(b)
  stopifnot(identical(g$ID_REF, ids_snp))
  for (j in b) {
    a <- amostras[gsm == j]
    x <- g[[j]]
    if (a$passou_cq) {
      n_AA[, a$grupo] <- n_AA[, a$grupo] + (x == "AA")
      n_AB[, a$grupo] <- n_AB[, a$grupo] + (x == "AB")
      n_BB[, a$grupo] <- n_BB[, a$grupo] + (x == "BB")
    }
    codigos_tipo[[j]] <- data.table(tipo = a$tipo, codigo = x)[, .N,
                                                               by = .(tipo, codigo)]
  }
  rm(g); invisible(gc())
}

# Frequência de cada código por tipo de linha (todas as amostras e SNPs)
codigos <- rbindlist(codigos_tipo)[, .(N = sum(N)), by = .(tipo, codigo)]
codigos[, prop := N / sum(N), by = tipo]
codigos
fwrite(codigos, file.path(dir_saida, "codigos_genotipo.tsv"), sep = "\t")

# Daqui em diante, só os SNPs que passaram no CQ
manter <- ids_snp %in% bim_cq$V2
n_AA <- n_AA[manter, ]; n_AB <- n_AB[manter, ]; n_BB <- n_BB[manter, ]

# -----------------------------------------------------------------------------
# 2. Funil do controle de qualidade
# -----------------------------------------------------------------------------
# Números tirados dos logs do PLINK (data/processed/*.log e
# results/structure/<conjunto>/<conjunto>.log).

log_cq <- readLines(paste0(prefixo_cq, ".log"))
num_log <- function(linhas, padrao) {
  as.integer(sub("^\\D*(\\d+).*", "\\1", grep(padrao, linhas, value = TRUE)[1]))
}
n_total   <- length(ids_snp)
n_mapeado <- num_log(log_cq, "variants loaded")
n_gentrain <- length(readLines(file.path(dir_proc, "snps_gentrain_baixo.txt")))
n_geno    <- num_log(log_cq, "removed due to missing genotype")
n_maf     <- num_log(log_cq, "removed due to minor allele")
n_final   <- num_log(log_cq, "pass filters and QC")

funil <- data.table(
  etapa = c("Chip (MaizeSNP50)", "Com posição no genoma", "GenTrain ≥ 0,7",
            "Call rate ≥ 90%", "MAF ≥ 1%"),
  snps  = c(n_total, n_mapeado, n_mapeado - n_gentrain,
            n_mapeado - n_gentrain - n_geno, n_final)
)
stopifnot(funil[.N, snps] == n_mapeado - n_gentrain - n_geno - n_maf)

for (nome in c("global", "dent", "flint")) {
  n <- length(readLines(file.path(dir_str, nome, paste0(nome, ".prune.in"))))
  funil <- rbind(funil, data.table(etapa = paste("Poda por LD:", nome), snps = n))
}
funil
fwrite(funil, file.path(dir_saida, "funil_cq.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 3. Distribuição da MAF
# -----------------------------------------------------------------------------
# MAF (frequência do alelo menos comum) de cada SNP após o CQ, no conjunto
# todo e dentro de cada pool. Um SNP comum no total pode ser raro (ou fixo)
# dentro de um pool: por isso os filtros foram refeitos por conjunto.

maf_de <- function(cols) {
  nA <- rowSums(n_AA[, cols, drop = FALSE]) * 2 + rowSums(n_AB[, cols, drop = FALSE])
  nB <- rowSums(n_BB[, cols, drop = FALSE]) * 2 + rowSums(n_AB[, cols, drop = FALSE])
  p <- nB / (nA + nB)
  pmin(p, 1 - p)
}
cols_dent  <- grep("^CFD|^Par_Dent", grupos, value = TRUE)
cols_flint <- grep("^CFF|^Par_Flint", grupos, value = TRUE)

maf <- rbind(data.table(conjunto = "global", maf = maf_de(grupos)),
             data.table(conjunto = "dent",   maf = maf_de(cols_dent)),
             data.table(conjunto = "flint",  maf = maf_de(cols_flint)))
maf[, .(mediana = median(maf, na.rm = TRUE),
        monomorficos = sum(maf == 0, na.rm = TRUE)), by = conjunto]
fwrite(maf, file.path(dir_saida, "maf.tsv.gz"), sep = "\t")

# -----------------------------------------------------------------------------
# 4. Decaimento do LD
# -----------------------------------------------------------------------------
# r² entre pares de SNPs no mesmo cromossomo, até 10 Mb de distância, em três
# grupos: DH dent, DH flint e os 23 parentais. Nas DH houve uma única meiose
# desde o cruzamento central x fundadora, então os blocos herdados são longos
# e o LD cai devagar. É isso que faz a poda por LD remover tantos SNPs.

ld_grupo <- function(rotulo, keep) {
  pref <- file.path(dir_saida, paste0("ld_", gsub(" ", "_", rotulo)))
  arq_keep <- paste0(pref, "_keep.txt")
  fwrite(keep[, .(familia, id)], arq_keep, sep = " ", col.names = FALSE)
  rodar_plink("--bfile", prefixo_cq, "--keep", arq_keep,
              "--r2", "--ld-window", "200", "--ld-window-kb", "10000",
              "--ld-window-r2", "0", "--out", pref)
  ld <- fread(paste0(pref, ".ld"), select = c("BP_A", "BP_B", "R2"))
  file.remove(paste0(pref, ".ld"))
  # média de r² em faixas de 250 kb
  ld[, faixa := floor(abs(BP_B - BP_A) / 250e3) * 0.25 + 0.125]
  ld[, .(r2 = mean(R2, na.rm = TRUE), pares = .N), by = faixa][
    , grupo := rotulo][order(faixa)]
}

a_cq <- amostras[(passou_cq)]
decaimento <- rbind(
  ld_grupo("DH dent",    a_cq[tipo == "DH" & pool == "Dent"]),
  ld_grupo("DH flint",   a_cq[tipo == "DH" & pool == "Flint"]),
  ld_grupo("Parentais",  a_cq[tipo == "Parental"])
)
setnames(decaimento, "faixa", "distancia_mb")
decaimento[distancia_mb < 1]
fwrite(decaimento, file.path(dir_saida, "ld_decaimento.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 5. PCA
# -----------------------------------------------------------------------------
# Todas as amostras que passaram no CQ, com os SNPs podados do conjunto
# global. A PCA não supõe nenhum modelo genético: só resume a variação em
# poucos eixos. É uma conferência independente do STRUCTURE.

pref_pca <- file.path(dir_saida, "pca")
rodar_plink("--bfile", prefixo_cq,
            "--extract", file.path(dir_str, "global", "global.prune.in"),
            "--pca", "10", "header", "tabs", "--out", pref_pca)

pca <- fread(paste0(pref_pca, ".eigenvec"))
autovalores <- scan(paste0(pref_pca, ".eigenval"), quiet = TRUE)
cat("\nVariância explicada pelos 10 primeiros PCs (% do total desses 10):\n")
round(100 * autovalores / sum(autovalores), 1)

# -----------------------------------------------------------------------------
# 6. Estatísticas F
# -----------------------------------------------------------------------------
# Wright dividiu a endogamia em três níveis (indivíduo I, subpopulação S,
# total T), comparando heterozigosidades:
#   H_O: heterozigosidade observada (fração de AB)
#   H_S: esperada sob Hardy-Weinberg dentro de cada subpopulação, 2p(1-p),
#        média entre subpopulações
#   H_T: esperada no total, 2 p̄ (1 - p̄), com p̄ a média das frequências
#   F_IS = 1 - H_O/H_S   endogamia dentro das subpopulações
#   F_ST = 1 - H_S/H_T   diferenciação entre subpopulações
#   F_IT = 1 - H_O/H_T   endogamia total;  (1 - F_IT) = (1 - F_IS)(1 - F_ST)
# H_S < H_T é o efeito Wahlund: juntar subpopulações diferentes produz
# "falta" de heterozigotos mesmo sem endogamia.
#
# Dois níveis de subpopulação: os pools (dent x flint) e as famílias DH
# dentro de cada pool. Os valores são razões das médias entre os SNPs.

estatisticas_f <- function(subpops) {
  # subpops: lista nomeada; cada elemento são as colunas (grupos) somadas
  por_sub <- lapply(subpops, function(cols) {
    aa <- rowSums(n_AA[, cols, drop = FALSE])
    ab <- rowSums(n_AB[, cols, drop = FALSE])
    bb <- rowSums(n_BB[, cols, drop = FALSE])
    n <- aa + ab + bb
    p <- (2 * bb + ab) / (2 * n)
    list(p = p, ho = ab / n, hs = 2 * p * (1 - p))
  })
  p  <- sapply(por_sub, `[[`, "p")
  ho <- rowMeans(sapply(por_sub, `[[`, "ho"))
  hs <- rowMeans(sapply(por_sub, `[[`, "hs"))
  pm <- rowMeans(p)
  ht <- 2 * pm * (1 - pm)
  ok <- is.finite(ho) & is.finite(hs) & ht > 0
  H_O <- mean(ho[ok]); H_S <- mean(hs[ok]); H_T <- mean(ht[ok])
  data.table(H_O = H_O, H_S = H_S, H_T = H_T,
             F_IS = 1 - H_O / H_S, F_ST = 1 - H_S / H_T, F_IT = 1 - H_O / H_T,
             snps = sum(ok))
}

fam_dent  <- grep("^CFD", grupos, value = TRUE)
fam_flint <- grep("^CFF", grupos, value = TRUE)

estat_f <- rbind(
  data.table(nivel = "Pools (dent x flint)",
             estatisticas_f(list(Dent = cols_dent, Flint = cols_flint))),
  data.table(nivel = "Famílias dent",
             estatisticas_f(setNames(as.list(fam_dent), fam_dent))),
  data.table(nivel = "Famílias flint",
             estatisticas_f(setNames(as.list(fam_flint), fam_flint)))
)
estat_f[, lapply(.SD, function(x) if (is.numeric(x)) round(x, 3) else x)]
fwrite(estat_f, file.path(dir_saida, "estatisticas_f.tsv"), sep = "\t")

# F_ST de Weir & Cockerham (1984) pelo PLINK, entre os pools, como
# conferência. Este estimador corrige pelo tamanho amostral; usa os dados do
# CQ, em que AB já é faltante (o que não afeta F_ST entre pools).
arq_pool <- file.path(dir_saida, "pools.txt")
fwrite(a_cq[, .(familia, id, pool)], arq_pool, sep = " ", col.names = FALSE)
saida_fst <- rodar_plink("--bfile", prefixo_cq, "--within", arq_pool,
                         "--fst", "--out", file.path(dir_saida, "fst_pools"))
grep("Fst estimate", sub("^.*%", "", saida_fst), value = TRUE)
