# =============================================================================
# BIO641 - Genética de Populações
# Regenera os arquivos de results/ a partir de um clone do repositório
#
# Uso (da raiz do projeto):
#   Rscript tools/gerar_resultados.R
#
# O STRUCTURE NÃO é executado aqui: as execuções foram feitas no cluster da
# UFV (levariam semanas num computador comum). Os resultados delas estão numa
# Release do GitHub, que este script baixa. A partir deles são geradas as
# tabelas e figuras usadas na apresentação.
#
# Etapas:
#   1. pacotes do R (instala os que faltarem)
#   2. resultados do STRUCTURE do cluster (Release do GitHub) -> data/
#   3. análise dos resultados do STRUCTURE (bio640_resultados.R)
#        -> results/structure/producao/analise/
#   4. dados do GEO + controle de qualidade + análises complementares
#      (bio640_popstructure.R e bio640_complementares.R)
#        -> results/complementares/
#      Precisa do PLINK 1.9. Sem ele, esta etapa é pulada com um aviso.
#
# Cada script roda num processo R separado, como se fosse chamado à mão.
# =============================================================================

# Permite rodar também de dentro de tools/
if (!file.exists("bio640_popstructure.R") && file.exists("../bio640_popstructure.R")) {
  setwd("..")
}
stopifnot("Rode a partir da raiz do projeto (BIO641_popstructure/)" =
            file.exists("bio640_popstructure.R"))

options(timeout = max(1200, getOption("timeout")))   # downloads grandes
inicio <- Sys.time()

etapa <- function(texto) cat("\n", strrep("=", 70), "\n ", texto, "\n",
                             strrep("=", 70), "\n", sep = "")

# Roda um script do projeto num processo R novo; para tudo se ele falhar.
# pdf(NULL) evita que gráficos abertos sem destino criem um Rplots.pdf.
rodar_script <- function(arquivo) {
  cat("\n$ Rscript", arquivo, "\n")
  status <- system2(file.path(R.home("bin"), "Rscript"),
                    c("-e", shQuote(sprintf("pdf(NULL); source('%s', echo = FALSE)",
                                            arquivo))))
  if (status != 0) stop("O script ", arquivo, " terminou com erro.")
}

baixar <- function(urls, destino, md5 = NULL) {
  for (url in urls) {
    cat("Baixando", basename(destino), "de", url, "\n")
    ok <- tryCatch({
      download.file(url, destino, mode = "wb", quiet = TRUE,
                    headers = c(Accept = "application/octet-stream"))
      TRUE
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (ok && file.exists(destino) && file.size(destino) > 0) break
  }
  if (!file.exists(destino) || file.size(destino) == 0) {
    stop("Não foi possível baixar ", basename(destino))
  }
  # confere se o arquivo chegou inteiro
  if (!is.null(md5) && unname(tools::md5sum(destino)) != md5) {
    stop("MD5 não confere: ", destino)
  }
}

# -----------------------------------------------------------------------------
# 1. Pacotes do R
# -----------------------------------------------------------------------------
etapa("1. Pacotes do R")

repos <- "https://cloud.r-project.org"
cran <- c("data.table", "ggplot2", "patchwork", "ggrepel", "remotes")
faltam <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(faltam) > 0) {
  cat("Instalando:", faltam, "\n")
  install.packages(faltam, repos = repos)
}
if (!requireNamespace("pophelper", quietly = TRUE)) {
  cat("Instalando pophelper (GitHub)\n")
  remotes::install_github("royfrancis/pophelper", upgrade = "never")
}
cat("Pacotes OK\n")

# -----------------------------------------------------------------------------
# 2. Resultados do STRUCTURE do cluster (Release do GitHub)
# -----------------------------------------------------------------------------
etapa("2. Resultados do STRUCTURE do cluster")

dir.create("data", showWarnings = FALSE)
dir_exportacao <- file.path("data", "exportacao_structure_R_20261015")

if (dir.exists(file.path(dir_exportacao, "resultados_completos"))) {
  cat("Já existe:", dir_exportacao, "\n")
} else {
  arq_tar <- file.path("data", "exportacao_structure_R_20261015.tar.gz")
  baixar(c(
    # link da Release; o da API é a alternativa se o primeiro falhar
    paste0("https://github.com/jpgteixeira/BIO641_popstructure/releases/download/",
           "structure-cluster-20261015/exportacao_structure_R_20261015.tar.gz"),
    "https://api.github.com/repos/jpgteixeira/BIO641_popstructure/releases/assets/657758906"),
    arq_tar,
    md5 = "a079ddc7867c665e1bdf2cb3f879dd01")
  untar(arq_tar, exdir = "data")
  file.remove(arq_tar)
  cat("Descompactado em", dir_exportacao, "\n")
}

# -----------------------------------------------------------------------------
# 3. Análise dos resultados do STRUCTURE
# -----------------------------------------------------------------------------
etapa("3. Análise dos resultados do STRUCTURE (pophelper)")

rodar_script("bio640_resultados.R")

# -----------------------------------------------------------------------------
# 4. Dados do GEO, controle de qualidade e análises complementares
# -----------------------------------------------------------------------------
etapa("4. Dados do GEO, controle de qualidade e análises complementares")

source("funcoes_popstructure.R")   # localiza o plink (variável PLINK ou PATH)
tem_plink <- nzchar(plink) && file.exists(plink)

if (!tem_plink) {
  cat("\nAVISO: PLINK 1.9 não encontrado; esta etapa foi pulada.\n",
      "Instale (https://www.cog-genomics.org/plink/1.9/) e informe o caminho\n",
      "no PATH ou na variável PLINK (por exemplo, no ~/.Renviron:\n",
      "PLINK=/caminho/para/plink). As tabelas de results/complementares/\n",
      "já versionadas no git continuam valendo.\n", sep = "")
} else {
  cat("PLINK:", plink, "\n")

  geo <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50558"
  arquivos_geo <- c(
    GSE50558_series_matrix.txt.gz =
      file.path(geo, "matrix", "GSE50558_series_matrix.txt.gz"),
    GSE50558_Parental_matrix_GEO.txt.gz =
      file.path(geo, "suppl", "GSE50558_Parental_matrix_GEO.txt.gz"))
  for (nome in names(arquivos_geo)) {
    destino <- file.path("data", nome)
    if (!file.exists(destino)) baixar(arquivos_geo[[nome]], destino)
  }
  # A anotação da plataforma vem em texto; é comprimida aqui
  arq_gpl <- file.path("data", "GPL17677_annotation.txt.gz")
  if (!file.exists(arq_gpl)) {
    tmp <- tempfile(fileext = ".txt")
    baixar(paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?",
                  "acc=GPL17677&targ=self&form=text&view=data"), tmp)
    con <- gzfile(arq_gpl, "w")
    writeLines(readLines(tmp), con)
    close(con)
  }

  # Partes 1 e 2: leitura do GEO, CQ, subamostra e poda por LD
  # (o STRUCTURE fica desligado: executar_structure = FALSE)
  rodar_script("bio640_popstructure.R")
  # PCA, LD, MAF, funil do CQ e estatísticas F
  rodar_script("bio640_complementares.R")
}

# -----------------------------------------------------------------------------
# Fim
# -----------------------------------------------------------------------------
etapa("Pronto")
gerados <- c(list.files(file.path("results", "structure", "producao", "analise"),
                        full.names = TRUE),
             list.files(file.path("results", "complementares"),
                        pattern = "\\.(tsv|tsv\\.gz|eigenvec)$", full.names = TRUE))
cat(gerados, sep = "\n")
cat("\nTempo total:", format(round(Sys.time() - inicio, 1)), "\n")
cat("Para gerar a apresentação: cd presentation && quarto render apresentacao.qmd\n")
