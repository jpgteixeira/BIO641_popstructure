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

# Resultados do cluster da UFV (produção), analisados pelo bio640_resultados.R
# (com fase <- "producao"). A apresentação só lê tabelas de results/, que são
# versionadas no git: dá para renderizar sem os dados brutos.
dir_prod   <- file.path(dir_res, "structure", "producao", "analise")

# Letras grandes em todos os gráficos
tema <- theme_minimal(base_size = 24) +
  theme(panel.grid.minor = element_blank(),
        plot.title = element_text(face = "bold", size = rel(1)),
        legend.position = "bottom",
        strip.text = element_text(face = "bold"))
theme_set(tema)

# Apresentação em inglês: ponto como separador decimal (padrão do R)
options(OutDec = ".")

cor_pool <- c(Dent = "#E31A1C", Flint = "#1F78B4")
cor_alelo <- c(A = "#F2B701", a = "#6A3D9A")

# Mesma paleta do bio640_resultados.R: o grupo n tem sempre a cor n
cores <- c("#1F78B4", "#E31A1C", "#33A02C", "#FF7F00", "#6A3D9A", "#FFD92F",
           "#A6CEE3", "#FB9A99", "#B2DF8A", "#FDBF6F", "#CAB2D6", "#B15928",
           "#1B9E77", "#E7298A", "#666666")
names(cores) <- paste0("Cluster", seq_along(cores))

q_medio  <- fread(file.path(dir_prod, "q_medio.tsv.gz"))
# Tabela de amostras (id, tipo, pool, família): o conjunto global com K = 1
# tem todas as 2.269 amostras, uma linha cada
amostras <- unique(q_medio[conjunto == "global" & K == 1,
                           .(id, tipo, pool, familia)])

# Gráfico de barras do STRUCTURE com ggplot: uma barra por indivíduo, um painel
# por K (linhas) e por família (colunas). Dentro da família, os indivíduos
# ficam na mesma ordem em todos os K (ordenados pelo grupo do primeiro K).
# K_destaque: um K para marcar como o "K provável" (pico do ΔK); a linha dele
# ganha moldura e o rótulo "ΔK peak"
barras <- function(conj, Ks, so_dh = TRUE, rotulo_familia = TRUE,
                   K_destaque = NULL) {
  d <- q_medio[conjunto == conj & K %in% Ks]
  if (so_dh) d <- d[tipo == "DH"]
  # grupo majoritário de cada indivíduo no primeiro K
  ordem <- d[K == min(Ks), .SD[which.max(fracao)], by = .(familia, id)]
  ordem <- ordem[order(familia, grupo, -fracao)][, x := seq_len(.N), by = familia]
  d <- merge(d, ordem[, .(id, x)], by = "id")
  rotulo_K <- function(k) ifelse(k %in% K_destaque,
                                 paste0("K = ", k, "\nΔK peak"), paste("K =", k))
  d[, painel_K := factor(rotulo_K(K), levels = rotulo_K(sort(Ks)))]
  d[, familia := sub("^CF", "", familia)]
  moldura <- unique(d[K %in% K_destaque, .(painel_K, familia)])
  # Rótulos dos K escritos ao lado da última família (em vez da faixa de
  # faceta), para o K destacado poder ficar em vermelho
  rotulos <- unique(d[, .(painel_K, K)])[, familia := max(d$familia)]
  rotulos[, destaque := K %in% K_destaque]
  ggplot(d, aes(x, fracao, fill = grupo)) +
    # contorno fino da mesma cor: tapa as frestas brancas que o PNG mostra
    # entre barras coladas (efeito da suavização das bordas, antialiasing)
    geom_col(aes(colour = grupo), width = 1, linewidth = 0.3) +
    # moldura só quando há um K destacado
    (if (nrow(moldura) > 0)
      geom_rect(data = moldura, aes(xmin = -Inf, xmax = Inf, ymin = -Inf,
                                    ymax = Inf), inherit.aes = FALSE,
                fill = NA, colour = "grey10", linewidth = 1.6)) +
    geom_text(data = rotulos, aes(x = Inf, y = 0.5, label = painel_K),
              inherit.aes = FALSE, hjust = -0.15, size = 6.5,
              fontface = "bold", lineheight = 0.9,
              colour = ifelse(rotulos$destaque, "#B22222", "grey10")) +
    facet_grid(painel_K ~ familia, scales = "free_x", space = "free_x",
               switch = "x") +
    scale_fill_manual(values = cores, guide = "none") +
    scale_colour_manual(values = cores, guide = "none") +
    scale_y_continuous(breaks = NULL, expand = c(0, 0)) +
    scale_x_continuous(breaks = NULL, expand = c(0, 0)) +
    coord_cartesian(clip = "off") +
    labs(x = NULL, y = NULL) +
    theme(panel.spacing.x = unit(1.5, "pt"), panel.grid = element_blank(),
          strip.text.x = element_text(angle = 90, size = rel(0.7),
                                      hjust = 1),
          strip.text.y = element_blank(),
          plot.margin = margin(5, 120, 5, 5))
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


# Diagrama de caixas e setas desenhado com ggplot (fica igual no HTML e no PDF)
#   caixas: data.table com id, x, y, rotulo e cor (preenchimento)
#   setas:  data.table com de, para e rotulo (texto sobre a seta, opcional);
#           colunas lx e ly (opcionais) dão a posição do rótulo; sem elas, o
#           rótulo fica no meio da seta
diagrama <- function(caixas, setas, largura = 2.4, altura = 1.3,
                     tamanho = 8, tamanho_rotulo = tamanho * 0.75) {
  s <- merge(setas, caixas[, .(de = id, x0 = x, y0 = y)], by = "de")
  s <- merge(s, caixas[, .(para = id, x1 = x, y1 = y)], by = "para")
  # a seta começa e termina na borda das caixas, não no centro
  s[, `:=`(dx = x1 - x0, dy = y1 - y0)]
  s[, t := pmin((largura / 2 + 0.08) / abs(dx), (altura / 2 + 0.08) / abs(dy))]
  s[, `:=`(xa = x0 + t * dx, ya = y0 + t * dy, xb = x1 - t * dx, yb = y1 - t * dy)]
  if (!"rotulo" %in% names(s)) s[, rotulo := ""]
  if (!"lx" %in% names(s)) s[, lx := (xa + xb) / 2]
  if (!"ly" %in% names(s)) s[, ly := (ya + yb) / 2]
  ggplot() +
    geom_segment(data = s, aes(xa, ya, xend = xb, yend = yb), linewidth = 1.5,
                 colour = "grey35",
                 arrow = arrow(length = unit(0.35, "cm"), type = "closed")) +
    geom_tile(data = caixas, aes(x, y, fill = I(cor)), width = largura,
              height = altura, colour = "grey30", linewidth = 0.8) +
    geom_text(data = caixas, aes(x, y, label = rotulo), size = tamanho,
              lineheight = 0.9, colour = ifelse(caixas$cor %in% c("white", "grey92"),
                                                "black", "white"),
              fontface = "bold") +
    # rótulos das setas por último, para nunca ficarem escondidos
    geom_label(data = s[rotulo != ""], aes(lx, ly, label = rotulo),
               size = tamanho_rotulo, label.size = 0, lineheight = 0.85,
               colour = "grey20", fontface = "bold") +
    coord_cartesian(clip = "off") +
    theme_void()
}

# Tema das ilustrações sem eixos (mesmo título dos demais gráficos)
# Tamanhos iguais aos do `tema`, para os títulos ficarem alinhados quando uma
# ilustração fica ao lado de um gráfico comum (slide "Why does it matter?")
tema_vazio <- theme_void(base_size = 24) +
  theme(plot.title = element_text(face = "bold", size = rel(1)),
        plot.subtitle = element_text(colour = "grey30", size = rel(0.9)))

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
