# =============================================================================
# BIO641 - Funções compartilhadas pelos scripts de estrutura populacional
#   bio640_popstructure.R  (análise local)
#   bio640_cluster.R       (pacote para o cluster da UFV)
# =============================================================================

plink     <- "/usr/local/bin/plink"
structure <- "/home/phenylalanine/structure/structure_linux_console/console/structure"

# O PLINK assume genoma humano por padrão. O milho tem 10 cromossomos
# autossômicos e nenhum cromossomo sexual.
plink_especie <- c("--chr-set", "10", "no-xy")

# Roda o PLINK, mostra as linhas informativas do log e para se der erro
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

# Troca o valor de um parâmetro "#define NOME valor" num arquivo de parâmetros
# do STRUCTURE (o parâmetro precisa existir exatamente uma vez)
definir <- function(linhas, parametro, valor) {
  padrao <- paste0("^#define\\s+", parametro, "\\s+\\S+")
  stopifnot(sum(grepl(padrao, linhas)) == 1)
  # format(): sem isso, 100000 vira "1e+05", que o STRUCTURE leria como 1
  valor <- format(valor, scientific = FALSE)
  sub(padrao, paste("#define", parametro, valor), linhas)
}

# Escreve mainparams e extraparams em `destino`, partindo dos modelos que vêm
# com o STRUCTURE.
#   mainparams: formato do arquivo gerado pelo PLINK e tamanho da MCMC
#   extraparams: modelo com mistura (NOADMIX 0), frequências alélicas
#                correlacionadas (FREQSCORR 1), como em Caniato et al. 2011,
#                e RANDOMIZE 0 para que a semente (-D) seja respeitada
escrever_parametros <- function(destino, burnin, numreps,
                                dir_modelos = dirname(structure)) {
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

# Define os três conjuntos de análise a partir de uma tabela de amostras
# (colunas id, tipo, pool, familia):
#   global: todas as famílias + 23 parentais
#   dent:   famílias CFD + parentais dent + UH007 (pai da CFD01)
#   flint:  famílias CFF + parentais flint + F353 (pai da CFF01) e B73 (pai da CFF02)
# Os parentais de fora do pool entram para podermos identificar a contribuição
# deles nas famílias "misturadas".
definir_conjuntos <- function(a) {
  list(
    global = a,
    dent   = a[startsWith(familia, "CFD") |
               (tipo == "Parental" & (pool == "Dent" | id == "UH007"))],
    flint  = a[startsWith(familia, "CFF") |
               (tipo == "Parental" & (pool == "Flint" | id %in% c("F353", "B73")))]
  )
}

# Para um conjunto: filtros de MAF e call rate recalculados dentro do conjunto,
# poda por LD (r² > 0,2 em janelas de 50 SNPs, passo 5) e exportação no formato
# do STRUCTURE. Devolve o arquivo, o número de indivíduos e de locos, e grava a
# ordem das amostras no arquivo (necessária para rotular os gráficos).
preparar_conjunto <- function(bfile, amostras_conj, pref) {
  dir.create(dirname(pref), showWarnings = FALSE, recursive = TRUE)
  fwrite(amostras_conj[, .(familia, id)], paste0(pref, "_keep.txt"),
         sep = " ", col.names = FALSE)

  rodar_plink("--bfile", bfile, "--keep", paste0(pref, "_keep.txt"),
              "--maf", "0.01", "--geno", "0.1",
              "--indep-pairwise", "50", "5", "0.2",
              "--out", pref)

  rodar_plink("--bfile", bfile, "--keep", paste0(pref, "_keep.txt"),
              "--extract", paste0(pref, ".prune.in"),
              "--recode", "structure", "--out", pref)

  # O PLINK reordena as amostras; guardamos a ordem do arquivo
  linhas <- readLines(paste0(pref, ".recode.strct_in"))
  ids <- sub(" .*", "", linhas[-(1:2)])
  ordem <- amostras_conj[match(ids, id), .(id, tipo, pool, familia)]
  fwrite(ordem, paste0(pref, "_ordem.tsv"), sep = "\t")

  list(arquivo = paste0(pref, ".recode.strct_in"),
       n_ind   = length(ids),
       n_loci  = length(strsplit(linhas[1], " ")[[1]]))
}
