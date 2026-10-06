# =============================================================================
# BIO641 - Genética de Populações
# Estrutura populacional dos painéis CornFed Dent e Flint (GEO GSE50558)
#
# Parte 1: leitura dos dados do GEO, montagem dos arquivos do PLINK e
#          controle de qualidade (CQ)
#
# Fluxo geral:
#   R (preparar dados) -> PLINK (CQ, poda por LD) -> STRUCTURE -> R (pophelper)
#
# Rodar a partir da raiz do projeto (BIO641_popstructure/).
# =============================================================================

library(data.table)

# -----------------------------------------------------------------------------
# 0. Caminhos e programas externos
# -----------------------------------------------------------------------------

dir_data <- "data"
dir_proc <- file.path(dir_data, "processed")   # arquivos gerados (ignorados no git)
dir.create(dir_proc, showWarnings = FALSE, recursive = TRUE)

arq_series <- file.path(dir_data, "GSE50558_series_matrix.txt.gz")
arq_gpl    <- file.path(dir_data, "GPL17677_annotation.txt.gz")
arq_pais   <- file.path(dir_data, "GSE50558_Parental_matrix_GEO.txt.gz")

plink     <- "/usr/local/bin/plink"
structure <- "/home/phenylalanine/structure/structure_linux_console/console/structure"

# O PLINK assume genoma humano por padrão. O milho tem 10 cromossomos
# autossômicos e nenhum cromossomo sexual.
plink_especie <- c("--chr-set", "10", "no-xy")

# Função auxiliar: roda o PLINK, mostra o log e para se der erro
rodar_plink <- function(...) {
  args <- c(plink_especie, ...)
  cat("\n$ plink", paste(args, collapse = " "), "\n")
  saida <- system2(plink, args, stdout = TRUE, stderr = TRUE)
  status <- attr(saida, "status")
  if (!is.null(status) && status != 0) {
    cat(saida, sep = "\n")
    stop("PLINK terminou com erro.")
  }
  # mostra só as linhas com contagens, que são as informativas
  # (removendo os indicadores de progresso "1%2%3%...")
  saida_limpa <- sub("^.*%", "", saida)
  cat(grep("variants|samples|people|pass|removed|loaded", saida_limpa,
           value = TRUE), sep = "\n")
  invisible(saida)
}

# -----------------------------------------------------------------------------
# 1. Metadados das amostras (cabeçalho do Series Matrix)
# -----------------------------------------------------------------------------
# As linhas que começam com "!Sample_" descrevem as 2.290 amostras (uma coluna
# por amostra). Interessam: o título (nome da linha), o status (parental ou DH)
# e o pool (Dent ou Flint).

cab <- readLines(arq_series, n = 80)
cab <- cab[startsWith(cab, "!Sample_")]
campos <- strsplit(cab, "\t", fixed = TRUE)
campos <- lapply(campos, function(x) gsub('"', "", x))
chave  <- vapply(campos, `[`, "", 1)

titulo <- campos[[which(chave == "!Sample_title")]][-1]
gsm    <- campos[[which(chave == "!Sample_geo_accession")]][-1]

# As características vêm em várias linhas "chave: valor", e a ordem muda entre
# amostras. Por isso buscamos o valor pela chave, e não pela posição.
carac <- do.call(rbind, lapply(campos[chave == "!Sample_characteristics_ch1"],
                               `[`, -1))
pegar_carac <- function(nome) {
  apply(carac, 2, function(col) {
    v <- col[startsWith(col, paste0(nome, ":"))]
    if (length(v) == 0) NA_character_ else sub("^[^:]*: ", "", v[1])
  })
}

amostras <- data.table(
  gsm    = gsm,
  titulo = titulo,
  status = pegar_carac("status"),
  pool   = pegar_carac("maize population"),
  info   = pegar_carac("population information")
)

# Nome da linha: "Parental line, B73" -> "B73"; "DH line, CFD01-003" -> "CFD01-003"
amostras[, id := sub("^.*, ", "", titulo)]
amostras[, tipo := fifelse(startsWith(status, "DH"), "DH", "Parental")]
# Família: para DH é o prefixo (CFD01); para os parentais usamos "Parental"
amostras[, familia := fifelse(tipo == "DH", sub("-.*$", "", id), "Parental")]

# Correção: no GEO, o parental EC169 está marcado como Flint, mas ele é Dent
# (Tabela S1 de Bauer et al. 2013; fundador da família dent CFD05 = F353 x EC169).
amostras[id == "EC169", pool := "Dent"]

amostras[, .N, by = .(tipo, pool)]
amostras[tipo == "DH", .N, by = .(pool, familia)][order(familia)]

# Conferência: o pool declarado bate com o prefixo da família?
# (CFD = Dent, CFF = Flint)
amostras[tipo == "DH", pool_prefixo := fifelse(startsWith(familia, "CFD"),
                                               "Dent", "Flint")]
inconsistentes <- amostras[tipo == "DH" & pool != pool_prefixo]
cat("\nLinhas DH com pool diferente do prefixo da família:", nrow(inconsistentes), "\n")
print(inconsistentes[, .(gsm, id, pool, familia)])

# -----------------------------------------------------------------------------
# 2. Genótipos (tabela do Series Matrix)
# -----------------------------------------------------------------------------
# Uma linha por SNP e uma coluna por amostra (identificada pelo GSM).
# Códigos: AA, BB (homozigotos), AB (heterozigoto), NC (sem chamada).

geno <- fread(arq_series, skip = "\"ID_REF\"", header = TRUE, quote = "\"")
geno <- geno[ID_REF != "!series_matrix_table_end"]
stopifnot(identical(names(geno)[-1], amostras$gsm))
dim(geno)   # 56.110 SNPs x (1 + 2.290 amostras)

# Frequência dos códigos de genótipo por tipo de amostra
contar_codigos <- function(cols) {
  table(unlist(geno[, ..cols], use.names = FALSE))
}
cat("\nCódigos nos parentais:\n");  print(contar_codigos(amostras[tipo == "Parental", gsm]))
cat("\nCódigos nas linhas DH:\n");  print(contar_codigos(amostras[tipo == "DH", gsm]))

# -----------------------------------------------------------------------------
# 3. Anotação dos SNPs (plataforma GPL17677) e GenTrain score
# -----------------------------------------------------------------------------
# A anotação traz cromossomo e posição na referência B73 v2.
# Cromossomos 0 e 99 = SNPs sem posição conhecida.

gpl <- fread(arq_gpl, skip = "ID\tIlmnID", sep = "\t", fill = TRUE,
             select = c("ID", "B73_v2_chr", "B73_v2_BP", "SNP"))
gpl <- gpl[ID != "!platform_table_end"]
setnames(gpl, c("snp", "chr", "pos", "alelos"))

# O GenTrain score usado pelo artigo (cluster file do projeto) está nos arquivos
# "supplementary"; basta o dos parentais, pois o valor é o mesmo para todas as
# amostras de um SNP.
gentrain <- fread(arq_pais, select = c("ID_REF", "B73.Score"))
setnames(gentrain, c("snp", "gentrain"))

snps <- merge(gpl, gentrain, by = "snp", all = TRUE)
snps <- snps[match(geno$ID_REF, snp)]   # mesma ordem da matriz de genótipos
stopifnot(identical(snps$snp, geno$ID_REF))

snps[, .N, by = chr][order(chr)]
summary(snps$gentrain)

# -----------------------------------------------------------------------------
# 4. Arquivos do PLINK no formato transposto (.tped / .tfam)
# -----------------------------------------------------------------------------
# O .tped tem uma linha por SNP, igual à nossa matriz, o que evita transpor
# 2.290 x 56.110 genótipos no R:
#   cromossomo  id_snp  posição_cM  posição_bp  alelo1 alelo2 (amostra 1) ...
# Codificamos os alelos como "A" e "B" (os mesmos do chip) e faltante como "0".
# SNPs sem posição (chr 0 e 99) ficam de fora.
#
# Heterozigotos (AB) viram dado faltante: todas as amostras são linhas
# endogâmicas (DH ou parentais), então um AB é erro de chamada ou heterozigose
# residual, e não informação útil sobre estrutura.

mapeados <- snps$chr %in% 1:10
cat("\nSNPs com posição no genoma:", sum(mapeados), "de", nrow(snps), "\n")

conv <- c(AA = "A A", BB = "B B", AB = "0 0", NC = "0 0")
tped <- geno[mapeados]
for (j in amostras$gsm) set(tped, j = j, value = unname(conv[tped[[j]]]))
tped[, ID_REF := NULL]
tped <- cbind(snps[mapeados, .(chr, snp, cm = 0, pos)], tped)

# .tfam: família, indivíduo, pai, mãe, sexo, fenótipo
# Família = CFDxx/CFFxx para as DH e "Parental" para os 23 parentais
tfam <- amostras[, .(fid = familia, iid = id, pai = 0, mae = 0, sexo = 0,
                     fenotipo = -9)]

prefixo_bruto <- file.path(dir_proc, "cornfed_bruto")
fwrite(tped, paste0(prefixo_bruto, ".tped"), sep = " ", col.names = FALSE,
       quote = FALSE)
fwrite(tfam, paste0(prefixo_bruto, ".tfam"), sep = " ", col.names = FALSE,
       quote = FALSE)
rm(tped); invisible(gc())

# Converte para o formato binário (.bed/.bim/.fam), mais compacto e rápido.
rodar_plink("--tfile", prefixo_bruto,
            "--make-bed", "--out", prefixo_bruto)
file.remove(paste0(prefixo_bruto, c(".tped", ".tfam")))

# Tabela de amostras para uso posterior (pool, família, tipo)
fwrite(amostras[, .(gsm, id, tipo, pool, familia)],
       file.path(dir_proc, "amostras.tsv"), sep = "\t")

# -----------------------------------------------------------------------------
# 5. Controle de qualidade (critérios de Lehermeier et al. 2014)
# -----------------------------------------------------------------------------
#   - GenTrain score < 0,7       -> SNP removido (feito com a lista abaixo)
#   - call rate < 0,9 por SNP    -> --geno 0.1
#   - MAF < 0,01                 -> --maf 0.01
#   - > 10% de faltantes por indivíduo -> removido, mas só para as linhas DH
#
# Os parentais são poucos e servem de referência para nomear os clusters, então
# não são removidos por dados faltantes (o UH304 tem ~11% de faltantes e seria
# perdido com --mind 0.1).

snps_gentrain_baixo <- snps[mapeados & gentrain < 0.7, snp]
cat("\nSNPs mapeados com GenTrain < 0,7:", length(snps_gentrain_baixo), "\n")
arq_excluir <- file.path(dir_proc, "snps_gentrain_baixo.txt")
writeLines(snps_gentrain_baixo, arq_excluir)

# Taxa de dados faltantes por amostra (após remover GenTrain baixo)
prefixo_faltantes <- file.path(dir_proc, "faltantes")
rodar_plink("--bfile", prefixo_bruto, "--exclude", arq_excluir,
            "--missing", "--out", prefixo_faltantes)
imiss <- fread(paste0(prefixo_faltantes, ".imiss"))
cat("\nFaltantes por amostra (F_MISS):\n"); print(summary(imiss$F_MISS))
print(imiss[F_MISS > 0.1, .(FID, IID, F_MISS)])

arq_remover <- file.path(dir_proc, "dh_muitos_faltantes.txt")
fwrite(imiss[F_MISS > 0.1 & FID != "Parental", .(FID, IID)], arq_remover,
       sep = " ", col.names = FALSE)

prefixo_cq <- file.path(dir_proc, "cornfed_cq")
rodar_plink("--bfile", prefixo_bruto,
            "--exclude", arq_excluir,
            "--remove", arq_remover,
            "--geno", "0.1",
            "--maf", "0.01",
            "--make-bed", "--out", prefixo_cq)

# Resumo após o CQ
bim_cq <- fread(paste0(prefixo_cq, ".bim"), header = FALSE)
fam_cq <- fread(paste0(prefixo_cq, ".fam"), header = FALSE)
cat("\nApós o CQ:", nrow(bim_cq), "SNPs e", nrow(fam_cq), "amostras\n")
cat("Referência (Lehermeier et al. 2014): 34.116 SNPs\n")

# Amostras removidas por excesso de dados faltantes
removidas <- amostras[!id %in% fam_cq$V2, .(id, tipo, pool, familia)]
cat("\nAmostras removidas no CQ:", nrow(removidas), "\n")
print(removidas[, .N, by = .(tipo, familia)])

# =============================================================================
# Parte 2: subamostragem, conjuntos de análise, poda por LD e STRUCTURE
# =============================================================================

# -----------------------------------------------------------------------------
# 6. Subamostragem: 30 linhas DH por família + todos os parentais
# -----------------------------------------------------------------------------
# O tempo do STRUCTURE cresce com o número de indivíduos. Com 30 DH por família
# mantemos as 24 famílias (inclusive os controles CFD01, CFF01 e CFF02) e
# reduzimos o conjunto de 2.290 para 743 amostras.

n_por_familia <- 30
set.seed(641)

amostras_cq <- amostras[id %in% fam_cq$V2]
sub_dh <- amostras_cq[tipo == "DH", .SD[sample(.N, min(.N, n_por_familia))],
                      by = familia]
subamostra <- rbind(amostras_cq[tipo == "Parental"], sub_dh, use.names = TRUE)
subamostra[, .N, by = .(tipo, pool)]

# -----------------------------------------------------------------------------
# 7. Três conjuntos de análise
# -----------------------------------------------------------------------------
#   global: todas as famílias + 23 parentais
#   dent:   famílias CFD + parentais dent + UH007 (pai da CFD01)
#   flint:  famílias CFF + parentais flint + F353 (pai da CFF01) e B73 (pai da CFF02)
# Os parentais de fora do pool entram para podermos identificar a contribuição
# deles nas famílias "misturadas".

conjuntos <- list(
  global = subamostra,
  dent   = subamostra[startsWith(familia, "CFD") |
               (tipo == "Parental" & (pool == "Dent" | id == "UH007"))],
  flint  = subamostra[startsWith(familia, "CFF") |
               (tipo == "Parental" & (pool == "Flint" | id %in% c("F353", "B73")))]
)
sapply(conjuntos, nrow)

# -----------------------------------------------------------------------------
# 8. Filtros por conjunto e exportação para o STRUCTURE (PLINK)
# -----------------------------------------------------------------------------
# Para cada conjunto:
#   - MAF e call rate recalculados dentro do conjunto (um SNP polimórfico entre
#     dent e flint pode ser fixo dentro de um pool);
#   - poda por LD: em janelas de 50 SNPs (passo de 5), de cada par com r² > 0,2
#     um SNP é removido. O STRUCTURE (modelo sem ligação) assume que os locos
#     são independentes dentro de cada população;
#   - exportação no formato do STRUCTURE (--recode structure): uma linha por
#     indivíduo, alelos codificados como 1/2 e faltante como 0.

dir_str <- file.path("results", "structure")
dir.create(dir_str, showWarnings = FALSE, recursive = TRUE)

info_conjuntos <- list()
for (nome in names(conjuntos)) {
  d <- file.path(dir_str, nome)
  dir.create(d, showWarnings = FALSE)
  pref <- file.path(d, nome)

  fwrite(conjuntos[[nome]][, .(familia, id)], paste0(pref, "_keep.txt"),
         sep = " ", col.names = FALSE)

  rodar_plink("--bfile", prefixo_cq, "--keep", paste0(pref, "_keep.txt"),
              "--maf", "0.01", "--geno", "0.1",
              "--indep-pairwise", "50", "5", "0.2",
              "--out", pref)

  rodar_plink("--bfile", prefixo_cq, "--keep", paste0(pref, "_keep.txt"),
              "--extract", paste0(pref, ".prune.in"),
              "--recode", "structure", "--out", pref)

  # O PLINK reordena as amostras; guardamos a ordem do arquivo para o pophelper
  linhas <- readLines(paste0(pref, ".recode.strct_in"))
  ids <- sub(" .*", "", linhas[-(1:2)])
  ordem <- conjuntos[[nome]][match(ids, id), .(id, tipo, pool, familia)]
  fwrite(ordem, paste0(pref, "_ordem.tsv"), sep = "\t")

  info_conjuntos[[nome]] <- list(arquivo = paste0(pref, ".recode.strct_in"),
                                 n_ind = length(ids),
                                 n_loci = length(strsplit(linhas[1], " ")[[1]]))
}
rbindlist(lapply(info_conjuntos, as.data.table), idcol = "conjunto")

# -----------------------------------------------------------------------------
# 9. Arquivos de parâmetros do STRUCTURE
# -----------------------------------------------------------------------------
# Partimos dos modelos que vêm com o programa e trocamos só o necessário.
#   mainparams: formato do arquivo (gerado pelo PLINK) e tamanho da MCMC
#   extraparams: modelo com mistura (NOADMIX 0), frequências alélicas
#                correlacionadas (FREQSCORR 1), como em Caniato et al. 2011,
#                e RANDOMIZE 0 para que a semente (-D) seja respeitada

definir <- function(linhas, parametro, valor) {
  padrao <- paste0("^#define\\s+", parametro, "\\s+\\S+")
  stopifnot(sum(grepl(padrao, linhas)) == 1)
  sub(padrao, paste("#define", parametro, valor), linhas)
}

dir_modelos <- dirname(structure)

escrever_parametros <- function(destino, burnin, numreps) {
  mp <- readLines(file.path(dir_modelos, "mainparams"))
  mp <- definir(mp, "BURNIN", burnin)
  mp <- definir(mp, "NUMREPS", numreps)
  mp <- definir(mp, "PLOIDY", 2)
  mp <- definir(mp, "MISSING", 0)
  mp <- definir(mp, "ONEROWPERIND", 1)
  mp <- definir(mp, "LABEL", 1)
  mp <- definir(mp, "POPDATA", 1)
  mp <- definir(mp, "POPFLAG", 0)
  mp <- definir(mp, "LOCDATA", 0)
  mp <- definir(mp, "PHENOTYPE", 0)
  mp <- definir(mp, "EXTRACOLS", 0)
  mp <- definir(mp, "MARKERNAMES", 1)
  mp <- definir(mp, "MAPDISTANCES", 1)
  writeLines(mp, file.path(destino, "mainparams"))

  ep <- readLines(file.path(dir_modelos, "extraparams"))
  ep <- definir(ep, "NOADMIX", 0)
  ep <- definir(ep, "LINKAGE", 0)
  ep <- definir(ep, "USEPOPINFO", 0)
  ep <- definir(ep, "FREQSCORR", 1)
  ep <- definir(ep, "INFERALPHA", 1)
  ep <- definir(ep, "RANDOMIZE", 0)
  writeLines(ep, file.path(destino, "extraparams"))
}

# -----------------------------------------------------------------------------
# 10. Execução do STRUCTURE (fase exploratória, local)
# -----------------------------------------------------------------------------
# Exploração: 10.000 de burn-in + 20.000 iterações, 3 réplicas por K.
# Cada execução usa um núcleo; rodamos várias em paralelo.
# Execuções que já terminaram (arquivo *_f existe) são puladas, então dá para
# interromper e retomar.

fase      <- "exploracao"
burnin    <- 10000
numreps   <- 20000
replicas  <- 1:3
faixa_K   <- list(global = 1:6, dent = 1:12, flint = 1:14)
n_nucleos <- 16

dir_fase <- file.path(dir_str, fase)
dir.create(dir_fase, showWarnings = FALSE)
escrever_parametros(dir_fase, burnin, numreps)

execucoes <- rbindlist(lapply(names(faixa_K), function(nome)
  CJ(conjunto = nome, K = faixa_K[[nome]], rep = replicas)))
# semente distinta e reprodutível para cada execução
execucoes[, semente := 1e5 * match(conjunto, names(faixa_K)) + 100 * K + rep]
execucoes[, saida := file.path(dir_fase, conjunto,
                               sprintf("%s_K%02d_r%d", conjunto, K, rep))]
# execuções mais longas primeiro, para distribuir melhor entre os núcleos
execucoes[, custo := info_conjuntos[conjunto][[1]]$n_ind, by = conjunto]
setorder(execucoes, -custo, -K)

rodar_structure <- function(i) {
  e <- execucoes[i]
  if (file.exists(paste0(e$saida, "_f"))) return("pulada")
  dir.create(dirname(e$saida), showWarnings = FALSE)
  inf <- info_conjuntos[[e$conjunto]]
  log <- system2(structure,
                 c("-m", file.path(dir_fase, "mainparams"),
                   "-e", file.path(dir_fase, "extraparams"),
                   "-i", inf$arquivo, "-o", e$saida,
                   "-N", inf$n_ind, "-L", inf$n_loci,
                   "-K", e$K, "-D", e$semente),
                 stdout = paste0(e$saida, ".log"), stderr = paste0(e$saida, ".log"))
  if (log != 0) "erro" else "ok"
}

cat("\nExecuções do STRUCTURE:", nrow(execucoes), "\n")
inicio <- Sys.time()
status <- parallel::mclapply(seq_len(nrow(execucoes)), rodar_structure,
                             mc.cores = n_nucleos, mc.preschedule = FALSE)
execucoes[, status := unlist(status)]
print(execucoes[, .N, by = .(conjunto, status)])
print(Sys.time() - inicio)
