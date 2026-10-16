# =============================================================================
# BIO641 - Genética de Populações
# Estrutura populacional dos painéis CornFed Dent e Flint (GEO GSE50558)
#
# Parte 3: análise dos resultados do STRUCTURE com o pophelper
#   - conferência das execuções e da estabilidade das cadeias
#   - ln P(D|K) e ΔK de Evanno
#   - alinhamento das réplicas e média por K
#   - gráficos de barras ordenados por família e dos parentais
#   - ancestralidade média por família
#
# Pré-requisito: execuções do STRUCTURE já terminadas, seja a exploração local
# (seção 10 do bio640_popstructure.R) ou a produção no cluster da UFV
# (exportada em data/exportacao_structure_R_20261015/). Este script só lê os
# arquivos *_f; não roda o STRUCTURE de novo.
#
# Rodar a partir da raiz do projeto (BIO641_popstructure/).
# =============================================================================

library(data.table)
library(pophelper)   # remotes::install_github("royfrancis/pophelper")
source("funcoes_popstructure.R")   # ler_traco()

# -----------------------------------------------------------------------------
# 0. Qual fase analisar
# -----------------------------------------------------------------------------
# "exploracao": 743 amostras (30 DH por família + parentais), 3 réplicas,
#               10 mil de burn-in + 20 mil iterações (rodada local)
# "producao":   todas as 2.269 amostras, 5 réplicas, 50 mil + 100 mil
#               (cluster da UFV). O cluster levaria ~14 dias para terminar
#               tudo; usamos o que ficou pronto: dent e flint completos e o
#               global sem 2 réplicas de K = 6 e 4 de K = 7.

fase <- "producao"

dir_str <- file.path("results", "structure")
dir_exportacao <- file.path("data", "exportacao_structure_R_20261015")

if (fase == "exploracao") {
  dir_resultados <- file.path(dir_str, "exploracao")      # <conjunto>/*_f
  dir_ordem      <- dir_str                               # <conjunto>/<conjunto>_ordem.tsv
} else {
  dir_resultados <- file.path(dir_exportacao, "resultados_completos")
  dir_ordem      <- file.path(dir_exportacao, "metadados", "dados")
}
dir_saida <- file.path(dir_str, fase, "analise")
dir.create(dir_saida, showWarnings = FALSE, recursive = TRUE)

conjuntos <- c("global", "dent", "flint")

# Ordem das amostras no arquivo de entrada do STRUCTURE (gravada na preparação
# dos dados). O *_f traz só a posição de cada indivíduo; é esse arquivo que diz
# quem é quem.
arq_ordem <- function(nome) {
  if (fase == "exploracao") {
    file.path(dir_ordem, nome, paste0(nome, "_ordem.tsv"))
  } else {
    file.path(dir_ordem, paste0(nome, "_ordem.tsv"))
  }
}

# -----------------------------------------------------------------------------
# 1. Leitura dos resultados
# -----------------------------------------------------------------------------
# readQ() lê cada arquivo *_f e guarda:
#   - a matriz Q: uma linha por indivíduo, uma coluna por grupo (cluster),
#     com a fração do genoma atribuída a cada grupo;
#   - os atributos da execução: K, burn-in, iterações, ln P(D|K) ("elpd"),
#     média ("mvll") e variância ("vll") do ln da verossimilhança.
# tabulateQ() junta esses atributos numa tabela, uma linha por execução.

q <- list()
execucoes <- list()
for (nome in conjuntos) {
  arqs <- list.files(file.path(dir_resultados, nome), pattern = "_f$",
                     full.names = TRUE)
  q[[nome]] <- readQ(arqs)
  execucoes[[nome]] <- as.data.table(tabulateQ(q[[nome]]))
}
execucoes <- rbindlist(execucoes, idcol = "conjunto")
execucoes[, rep := as.integer(sub(".*_r(\\d+)_f$", "\\1", file))]

# Quantas execuções por conjunto e K (deve ser o número de réplicas em todos;
# na produção faltam réplicas de K = 6 e 7 no global)
dcast(execucoes, conjunto ~ k, value.var = "rep", fun.aggregate = length)

# -----------------------------------------------------------------------------
# 2. Estabilidade das cadeias (convergência)
# -----------------------------------------------------------------------------
# Depois do burn-in, o ln da verossimilhança deve oscilar em torno de um valor
# estável. Se a cadeia ainda estiver "pulando" entre soluções, a variância
# (vll) fica muito maior do que nas outras réplicas do mesmo K. Isso importa
# porque o STRUCTURE calcula ln P(D|K) = média - variância / 2: uma variância
# enorme derruba o ln P(D|K) e distorce o ΔK de Evanno.
#
# Critério: réplica instável se vll > 3 x a mediana das réplicas do mesmo K.
# Com 3 ou mais réplicas, a mediana não é afetada por uma única réplica ruim.
# (Com uma réplica só, como o K = 7 global da produção, não há o que comparar.)

execucoes[, instavel := vll > 3 * median(vll), by = .(conjunto, k)]
execucoes[(instavel), .(conjunto, K = k, rep, elpd, vll)]

fwrite(execucoes, file.path(dir_saida, "execucoes.tsv"), sep = "\t")

# As réplicas instáveis ficam de fora daqui em diante
for (nome in conjuntos) {
  fora <- execucoes[conjunto == nome & (instavel), file]
  q[[nome]] <- q[[nome]][!names(q[[nome]]) %in% fora]
}

# -----------------------------------------------------------------------------
# 3. ln P(D|K) e ΔK de Evanno
# -----------------------------------------------------------------------------
# Método de Evanno et al. (2005), a partir da média e do desvio-padrão de
# ln P(D|K) entre as réplicas:
#   L'(K)  = L(K) - L(K-1)                ganho ao acrescentar um grupo
#   L''(K) = |L(K+1) - 2 L(K) + L(K-1)|   mudança desse ganho
#   ΔK     = L''(K) / dp(L(K))
# O ΔK é alto onde a curva de ln P(D|K) "faz um cotovelo". Ele não existe para
# o primeiro e o último K testados, nem para um K com uma única réplica (sem
# desvio-padrão; o pophelper avisa, mas calcula os demais).
#
# Cuidado na interpretação: o ΔK tende a apontar o nível mais alto da
# hierarquia (aqui, dent x flint em K = 2). Os picos secundários e os gráficos
# de barras é que mostram a estrutura dentro dos pools.

evanno <- list()
for (nome in conjuntos) {
  resumo <- summariseQ(tabulateQ(q[[nome]]))
  evanno[[nome]] <- as.data.table(evannoMethodStructure(
    resumo, exportplot = TRUE, exportpath = dir_saida,
    outputfilename = paste0("evanno_", nome), imgtype = "png",
    width = 16, height = 12, basesize = 8))
}
evanno <- rbindlist(evanno, idcol = "conjunto")
fwrite(evanno, file.path(dir_saida, "evanno.tsv"), sep = "\t")

# Os três maiores ΔK de cada conjunto
evanno[!is.na(deltaK)][order(-deltaK), head(.SD, 3), by = conjunto][
  , .(conjunto, K = k, elpdmean, elpdsd, deltaK = round(deltaK, 1))]

# -----------------------------------------------------------------------------
# 4. Alinhamento das réplicas e média por K
# -----------------------------------------------------------------------------
# Os grupos do STRUCTURE não têm nome: o "grupo 1" de uma réplica pode ser o
# "grupo 3" de outra (troca de rótulos, ou label switching). Antes de tirar a
# média das réplicas, é preciso alinhá-las.
#   a) dentro de cada K: alignK(type = "within"), semelhante ao CLUMPP;
#   b) modo principal: réplicas podem chegar a soluções diferentes (no dent,
#      K = 2, três réplicas dividem as fundadoras de um jeito e duas de
#      outro). A média de soluções diferentes "borra" o gráfico. Por isso
#      ficam só as réplicas parecidas com a de maior ln P(D|K), como faz o
#      CLUMPAK ("major mode");
#   c) média das réplicas do modo principal;
#   d) entre K: cada K é alinhado ao K anterior, para que um mesmo grupo
#      mantenha a mesma cor nos gráficos.
#
# Para (b) não usamos o mergeQ() do pophelper: na versão 2.3.1 ele falha no
# R 4.6 (erro em sortQ). A média é simples de fazer diretamente.
# Para (c) não usamos alignK(type = "across"): ele troca as cores de grupos
# que não mudaram de um K para o seguinte.

# Reordena as colunas de `q` para que fiquem o mais parecidas possível com as
# de `ref` (o K anterior). A semelhança entre um grupo novo e um antigo é a
# soma, entre os indivíduos, da menor das duas frações (sobreposição).
# Os pares são escolhidos do mais parecido para o menos parecido; o grupo
# que sobra (o "novo" grupo deste K) vai para a última coluna.
alinhar_ao_anterior <- function(q, ref) {
  sobrep <- sapply(seq_len(ncol(ref)), function(j)
    sapply(seq_len(ncol(q)), function(i) sum(pmin(q[[i]], ref[[j]]))))
  ordem <- rep(NA_integer_, ncol(ref))
  for (passo in seq_len(ncol(ref))) {
    pos <- which(sobrep == max(sobrep), arr.ind = TRUE)[1, ]
    ordem[pos[2]] <- pos[1]
    sobrep[pos[1], ] <- -Inf
    sobrep[, pos[2]] <- -Inf
  }
  ordem <- c(ordem, setdiff(seq_len(ncol(q)), ordem))
  saida <- q[, ordem, drop = FALSE]
  names(saida) <- paste0("Cluster", seq_len(ncol(q)))
  saida
}

# Semelhança entre duas matrizes Q já alinhadas: fração do genoma atribuída
# ao mesmo grupo, em média entre os indivíduos (1 = idênticas).
semelhanca <- function(a, b) mean(rowSums(pmin(as.matrix(a), as.matrix(b))))

# Limite para duas réplicas estarem no mesmo modo
limite_modo <- 0.9

modos <- list()   # registro: quantas réplicas ficaram no modo principal

media_por_K <- function(qlista, nome) {
  k <- sapply(qlista, ncol)
  medias <- lapply(split(qlista, k), function(reps) {
    reps <- alignK(as.qlist(reps), type = "within")
    elpd <- sapply(reps, attr, "elpd")
    melhor <- reps[[which.max(elpd)]]
    sem <- sapply(reps, semelhanca, melhor)
    principal <- reps[sem >= limite_modo]
    modos[[length(modos) + 1]] <<- data.table(
      conjunto = nome, K = ncol(melhor), replicas = length(reps),
      no_modo_principal = length(principal),
      semelhanca_minima = round(min(sem), 3))
    Reduce(`+`, principal) / length(principal)
  })
  medias <- medias[order(as.integer(names(medias)))]
  for (i in seq_along(medias)[-1]) {
    medias[[i]] <- alinhar_ao_anterior(medias[[i]], medias[[i - 1]])
  }
  names(medias) <- paste0("K", names(medias))
  as.qlist(medias)
}

q_media <- lapply(conjuntos, function(nome) media_por_K(q[[nome]], nome))
names(q_media) <- conjuntos

modos <- rbindlist(modos)
modos[no_modo_principal < replicas]   # os K com mais de uma solução
fwrite(modos, file.path(dir_saida, "modos.tsv"), sep = "\t")

# Matrizes Q médias em formato longo: uma linha por indivíduo, K e grupo
amostras <- lapply(conjuntos, function(nome) fread(arq_ordem(nome)))
names(amostras) <- conjuntos

q_longo <- rbindlist(lapply(conjuntos, function(nome) {
  rbindlist(lapply(names(q_media[[nome]]), function(kk) {
    d <- cbind(amostras[[nome]], q_media[[nome]][[kk]])
    d <- melt(d, id.vars = names(amostras[[nome]]),
              variable.name = "grupo", value.name = "fracao")
    d[, K := as.integer(sub("K", "", kk))]
  }))
}), idcol = "conjunto")
q_longo[, conjunto := conjuntos[conjunto]]
fwrite(q_longo, file.path(dir_saida, "q_medio.tsv.gz"), sep = "\t")

# -----------------------------------------------------------------------------
# 5. Gráficos de barras
# -----------------------------------------------------------------------------
# Cada barra é um indivíduo; as cores são as frações de cada grupo. As DH
# ficam agrupadas por família e, dentro da família, ordenadas pelo grupo
# majoritário (sortind = "all"). Um painel por K, de K = 2 ao maior testado.
# Os parentais vão num gráfico separado, com o nome de cada linha, porque são
# poucos (23) e se perdem no meio das DH.
#
# Paleta fixa: o pophelper troca de paleta quando K passa de 12, e aí a mesma
# cor deixaria de indicar o mesmo grupo entre gráficos. Como os K estão
# alinhados (seção 4), o grupo n tem sempre a cor n.
cores <- c("#1F78B4", "#E31A1C", "#33A02C", "#FF7F00", "#6A3D9A", "#FFD92F",
           "#A6CEE3", "#FB9A99", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#B15928",
           "#1B9E77", "#E7298A", "#666666")

# Seleciona as linhas `linhas` de cada matriz Q da lista
selecionar <- function(qlista, linhas, nomes = NULL) {
  as.qlist(lapply(qlista, function(x) {
    x <- x[linhas, , drop = FALSE]
    if (!is.null(nomes)) rownames(x) <- nomes
    x
  }))
}

for (nome in conjuntos) {
  a <- amostras[[nome]]
  qk <- q_media[[nome]][-1]   # sem K = 1, que é uma barra única
  rot_K <- sub("K", "K = ", names(qk))

  # Só as DH, agrupadas por família
  dh <- which(a$tipo == "DH")
  plotQ(selecionar(qk, dh), imgoutput = "join", clustercol = cores,
        grplab = data.frame(familia = a$familia[dh], stringsAsFactors = FALSE),
        ordergrp = TRUE, sortind = "all", sharedindlab = FALSE,
        splab = rot_K, splabsize = 6,
        grplabangle = 90, grplabsize = 1.8, grplabheight = 3,
        grplabpos = 0.8, grplabjust = 1,
        width = 30, height = 2.2, basesize = 7,
        exportpath = dir_saida, outputfilename = paste0("barras_", nome))

  # Só os parentais, com o nome de cada linha
  par <- which(a$tipo == "Parental")
  plotQ(selecionar(qk, par, a$id[par]), imgoutput = "join", clustercol = cores,
        grplab = data.frame(pool = a$pool[par], stringsAsFactors = FALSE),
        ordergrp = TRUE, sharedindlab = TRUE, showindlab = TRUE,
        useindlab = TRUE, indlabsize = 5, indlabheight = 0.6,
        splab = rot_K, splabsize = 6, grplabsize = 2.5,
        width = 12, height = 2.2, basesize = 7,
        exportpath = dir_saida, outputfilename = paste0("parentais_", nome))
}

# -----------------------------------------------------------------------------
# 6. Ancestralidade média por família
# -----------------------------------------------------------------------------
# Média das frações de cada grupo dentro de cada família (e dos parentais de
# cada pool). Permite ler numericamente o que os gráficos mostram.

anc_familia <- q_longo[, .(fracao = mean(fracao)),
                       by = .(conjunto, K,
                              familia = ifelse(tipo == "Parental",
                                               paste("Par.", pool), familia),
                              grupo)]
fwrite(anc_familia, file.path(dir_saida, "ancestralidade_familias.tsv"),
       sep = "\t")

# Conferência com o esperado (docs/estrutura_dos_dados.md): no conjunto
# global com K = 2, os grupos devem separar dent e flint, e as famílias com
# pais dos dois pools (CFD01 = F353 x UH007, CFF01 = UH007 x F353,
# CFF02 = UH007 x B73) devem ter cerca de 50% de cada grupo.
dcast(anc_familia[conjunto == "global" & K == 2],
      familia ~ grupo, value.var = "fracao")[
  , lapply(.SD, function(x) if (is.numeric(x)) round(x, 2) else x)]

# -----------------------------------------------------------------------------
# 7. Traços da MCMC
# -----------------------------------------------------------------------------
# O .log de cada execução tem, a cada 100 iterações, o alfa e (depois do
# burn-in) o ln da verossimilhança. Juntamos tudo numa tabela pequena, usada
# nos gráficos de convergência da apresentação (os .log da produção somam ~110 MB
# e não vão para o git).

tracos <- rbindlist(lapply(conjuntos, function(nome) {
  logs <- list.files(file.path(dir_resultados, nome), pattern = "_r\\d+\\.log$",
                     full.names = TRUE)
  rbindlist(lapply(logs, function(arq) {
    ler_traco(arq)[, `:=`(
      conjunto = nome,
      K   = as.integer(sub(".*_K(\\d+)_r.*", "\\1", basename(arq))),
      rep = as.integer(sub(".*_r(\\d+)\\.log$", "\\1", basename(arq))))]
  }))
}))
# só as execuções que terminaram (as interrompidas não têm *_f)
tracos <- tracos[execucoes[, .(conjunto, K = k, rep)], on = .(conjunto, K, rep),
                 nomatch = NULL]
setcolorder(tracos, c("conjunto", "K", "rep"))
fwrite(tracos, file.path(dir_saida, "tracos.tsv.gz"), sep = "\t")

cat("\nResultados em", dir_saida, "\n")
