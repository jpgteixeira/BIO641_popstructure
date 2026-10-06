# =============================================================================
# BIO641 - Genética de Populações
# Pacote para rodar o STRUCTURE no cluster da UFV (conjunto completo)
#
# Pré-requisito: ter rodado as seções 1-5 do bio640_popstructure.R, que geram
#   data/processed/cornfed_cq.{bed,bim,fam}   (dados após o CQ)
#   data/processed/amostras.tsv               (pool e família de cada amostra)
#
# Saída: cluster/pacote/, pronto para enviar ao cluster (ver docs/cluster_ufv.md)
#
# Rodar a partir da raiz do projeto (BIO641_popstructure/).
# =============================================================================

library(data.table)
source("funcoes_popstructure.R")

# -----------------------------------------------------------------------------
# 1. Configuração da análise final
# -----------------------------------------------------------------------------
# MCMC: 50.000 de burn-in + 100.000 iterações, 5 réplicas por K.
# Faixas de K provisórias: ajustar depois da exploração local.
# Walltime com folga sobre o tempo medido localmente (global ~15 h, pools
# ~4-6 h por execução), dentro do máximo de 72 h da fila qtime.

burnin   <- 50000
numreps  <- 100000
replicas <- 1:5
faixa_K  <- list(global = 1:8, dent = 1:13, flint = 1:14)
walltime <- c(global = "48:00:00", dent = "24:00:00", flint = "24:00:00")

prefixo_cq <- file.path("data", "processed", "cornfed_cq")
dir_pacote <- file.path("cluster", "pacote")

# -----------------------------------------------------------------------------
# 2. Amostras: todas as que passaram no CQ (sem subamostragem)
# -----------------------------------------------------------------------------

amostras <- fread(file.path("data", "processed", "amostras.tsv"))
fam_cq   <- fread(paste0(prefixo_cq, ".fam"), header = FALSE)
amostras <- amostras[id %in% fam_cq$V2]
amostras[, .N, by = .(tipo, pool)]

conjuntos <- definir_conjuntos(amostras)
sapply(conjuntos, nrow)

# -----------------------------------------------------------------------------
# 3. Arquivos de entrada do STRUCTURE (filtros e poda por LD em cada conjunto)
# -----------------------------------------------------------------------------

for (d in c("bin", "dados", "params", "resultados"))
  dir.create(file.path(dir_pacote, d), showWarnings = FALSE, recursive = TRUE)

info <- lapply(names(conjuntos), function(nome)
  preparar_conjunto(prefixo_cq, conjuntos[[nome]],
                    file.path(dir_pacote, "dados", nome)))
names(info) <- names(conjuntos)
rbindlist(lapply(info, as.data.table), idcol = "conjunto")

# Os arquivos intermediários do PLINK (.prune.in, .log, ...) não precisam ir
# para o cluster; ficam só o arquivo do STRUCTURE e a ordem das amostras.
intermediarios <- list.files(file.path(dir_pacote, "dados"), full.names = TRUE)
intermediarios <- intermediarios[!grepl("\\.recode\\.strct_in$|_ordem\\.tsv$",
                                        intermediarios)]
invisible(file.remove(intermediarios))

# -----------------------------------------------------------------------------
# 4. Parâmetros do STRUCTURE
# -----------------------------------------------------------------------------

escrever_parametros(file.path(dir_pacote, "params"), burnin, numreps)

# -----------------------------------------------------------------------------
# 5. Tabela de execuções (uma linha = um job no cluster)
# -----------------------------------------------------------------------------
# Sementes diferentes das da exploração local (que começam em 1e5).

execucoes <- rbindlist(lapply(names(faixa_K), function(nome)
  CJ(conjunto = nome, K = faixa_K[[nome]], rep = replicas)))
execucoes[, semente := 1e6 * match(conjunto, names(faixa_K)) + 100 * K + rep]
execucoes[, n_ind := info[conjunto][[1]]$n_ind, by = conjunto]
execucoes[, n_loci := info[conjunto][[1]]$n_loci, by = conjunto]
execucoes[, walltime := walltime[conjunto]]
# execuções mais longas primeiro: entram antes na fila
setorder(execucoes, -n_ind, -K, rep)

fwrite(execucoes[, .(conjunto, K, rep, semente = format(semente, scientific = FALSE),
                     n_ind, n_loci, walltime)],
       file.path(dir_pacote, "execucoes.tsv"), sep = "\t")
execucoes[, .N, by = conjunto]

# -----------------------------------------------------------------------------
# 6. Scripts do cluster
# -----------------------------------------------------------------------------

scripts <- c("structure_job.pbs", "submeter.sh", "compilar_structure.sh")
file.copy(file.path("cluster", scripts), file.path(dir_pacote, scripts),
          overwrite = TRUE)
Sys.chmod(file.path(dir_pacote, c("submeter.sh", "compilar_structure.sh")), "755")

cat("\nPacote pronto em", dir_pacote, "\n")
cat("Jobs:", nrow(execucoes), "\n")
cat("Lembrete: editar o e-mail (#PBS -M) em structure_job.pbs antes de enviar.\n")
