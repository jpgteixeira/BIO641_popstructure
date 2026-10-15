# =============================================================================
# Configuração comum dos gráficos da apresentação
# Lido pelo chunk "setup" de apresentacao.qmd (renderizado a partir de
# presentation/, por isso os caminhos começam com "..").
# =============================================================================

library(data.table)
library(ggplot2)
library(patchwork)
library(ggrepel)

dir_res  <- file.path("..", "results")
dir_expl <- file.path(dir_res, "structure", "exploracao")
dir_comp <- file.path(dir_res, "complementares")

# Letras grandes em todos os gráficos
tema <- theme_minimal(base_size = 24) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = rel(1)),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))
theme_set(tema)

# Vírgula como separador decimal nos eixos e rótulos
options(OutDec = ",")

cor_pool <- c(Dent = "#E31A1C", Flint = "#1F78B4")
cor_alelo <- c(A = "#F2B701", a = "#6A3D9A")

# Mesma paleta do bio640_resultados.R: o grupo n tem sempre a cor n
cores <- c("#1F78B4", "#E31A1C", "#33A02C", "#FF7F00", "#6A3D9A", "#FFD92F",
           "#A6CEE3", "#FB9A99", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#B15928",
           "#1B9E77", "#E7298A", "#666666")
names(cores) <- paste0("Cluster", seq_along(cores))

amostras <- fread(file.path("..", "data", "processed", "amostras.tsv"))
q_medio  <- fread(file.path(dir_expl, "analise", "q_medio.tsv.gz"))

# Gráfico de barras do STRUCTURE com ggplot: uma barra por indivíduo, um painel
# por K (linhas) e por família (colunas). Dentro da família, os indivíduos
# ficam na mesma ordem em todos os K (ordenados pelo grupo do primeiro K).
barras <- function(conj, Ks, so_dh = TRUE, rotulo_familia = TRUE) {
  d <- q_medio[conjunto == conj & K %in% Ks]
  if (so_dh) d <- d[tipo == "DH"]
  # grupo majoritário de cada indivíduo no primeiro K
  ordem <- d[K == min(Ks), .SD[which.max(fracao)], by = .(familia, id)]
  ordem <- ordem[order(familia, grupo, -fracao)][, x := seq_len(.N), by = familia]
  d <- merge(d, ordem[, .(id, x)], by = "id")
  d[, painel_K := factor(paste("K =", K), levels = paste("K =", sort(Ks)))]
  d[, familia := sub("^CF", "", familia)]
  ggplot(d, aes(x, fracao, fill = grupo)) +
    geom_col(width = 1) +
    facet_grid(painel_K ~ familia, scales = "free_x", space = "free_x",
               switch = "x") +
    scale_fill_manual(values = cores, guide = "none") +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme(panel.spacing.x = unit(1.5, "pt"), panel.grid = element_blank(),
          strip.text.x = element_text(angle = 90, size = rel(0.7),
                                      hjust = 1),
          strip.text.y = element_text(angle = 0, size = rel(0.8)))
}

# Só os parentais, com o nome de cada linha
barras_parentais <- function(conj, Ks) {
  d <- q_medio[conjunto == conj & K %in% Ks & tipo == "Parental"]
  d[, painel_K := factor(paste("K =", K), levels = paste("K =", sort(Ks)))]
  ggplot(d, aes(id, fracao, fill = grupo)) +
    geom_col(width = 0.95) +
    facet_grid(painel_K ~ pool, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = cores, guide = "none") +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = rel(0.75)),
          strip.text.y = element_text(angle = 0, size = rel(0.8)))
}

# Traço da MCMC a partir do .log de uma execução do STRUCTURE: a cada 100
# iterações o programa imprime alfa e, depois do burn-in, o ln da
# verossimilhança. As linhas do traço são as que têm Lambda = 1.00.
ler_traco <- function(arq) {
  l <- grep("^\\s*\\d+:\\s+1\\.00\\s", readLines(arq), value = TRUE)
  campos <- strsplit(trimws(l), "\\s+")
  data.table(
    iteracao = as.integer(sub(":", "", sapply(campos, `[`, 1))),
    alpha    = as.numeric(sapply(campos, `[`, 3)),
    # no burn-in a última coluna é "--"; depois, ln da verossimilhança e um
    # contador, nessa ordem
    lnlike   = sapply(campos, function(x)
      if (x[length(x)] == "--") NA_real_ else as.numeric(x[length(x) - 1]))
  )[, fase := fifelse(is.na(lnlike), "burn-in", "amostragem")]
}

# Diagrama de caixas e setas desenhado com ggplot (fica igual no HTML e no PDF)
#   caixas: data.table com id, x, y, rotulo e cor (preenchimento)
#   setas:  data.table com de, para e rotulo (texto sobre a seta, opcional)
diagrama <- function(caixas, setas, largura = 2.4, altura = 1.3,
                     tamanho = 8) {
  s <- merge(setas, caixas[, .(de = id, x0 = x, y0 = y)], by = "de")
  s <- merge(s, caixas[, .(para = id, x1 = x, y1 = y)], by = "para")
  # a seta começa e termina na borda das caixas, não no centro
  s[, `:=`(dx = x1 - x0, dy = y1 - y0)]
  s[, t := pmin((largura / 2 + 0.08) / abs(dx), (altura / 2 + 0.08) / abs(dy))]
  s[, `:=`(xa = x0 + t * dx, ya = y0 + t * dy, xb = x1 - t * dx, yb = y1 - t * dy)]
  if (!"rotulo" %in% names(s)) s[, rotulo := ""]
  ggplot() +
    geom_segment(data = s, aes(xa, ya, xend = xb, yend = yb), linewidth = 1.5,
                 colour = "grey35",
                 arrow = arrow(length = unit(0.35, "cm"), type = "closed")) +
    geom_label(data = s[rotulo != ""], aes((xa + xb) / 2, (ya + yb) / 2,
                                           label = rotulo),
               size = tamanho * 0.75, label.size = 0, lineheight = 0.85,
               colour = "grey25") +
    geom_tile(data = caixas, aes(x, y, fill = I(cor)), width = largura,
              height = altura, colour = "grey30", linewidth = 0.8) +
    geom_text(data = caixas, aes(x, y, label = rotulo), size = tamanho,
              lineheight = 0.9, colour = ifelse(caixas$cor %in% c("white", "grey92"),
                                                "black", "white"),
              fontface = "bold") +
    coord_cartesian(clip = "off") +
    theme_void()
}

# Números grandes com legenda, lado a lado
numeros <- function(valores, legendas) {
  d <- data.table(x = seq_along(valores), valor = valores, legenda = legendas)
  ggplot(d) +
    geom_text(aes(x, 1, label = valor), size = 26, fontface = "bold",
              colour = "#1F3A5F") +
    geom_text(aes(x, 0.35, label = legenda), size = 9, colour = "grey30",
              lineheight = 0.9) +
    scale_y_continuous(limits = c(0, 1.4)) +
    scale_x_continuous(limits = c(0.5, length(valores) + 0.5)) +
    theme_void()
}

# No PDF (beamer) o slide é baixo: limita a altura das figuras para sobrar
# espaço para o título e a frase de destaque
if (knitr::is_latex_output()) {
  knitr::opts_chunk$set(out.height = "0.66\\textheight", out.width = NULL)
}
