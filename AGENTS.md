# Instruções para assistentes de IA

Leia primeiro o [README.md](README.md): ele descreve o objetivo, os dados, o
fluxo de análise, as decisões já tomadas e o que falta fazer. Os documentos em
`docs/` trazem os detalhes dos dados e do cluster.

## Contexto

- Trabalho acadêmico de Genética de Populações (UFV). Os autores estão
  **aprendendo** a ferramenta: explique o raciocínio além de entregar o código.
- Idioma: **português** em textos, comentários de código e nomes de variáveis.
- O professor exige o **STRUCTURE**. ADMIXTURE ou LEA/sNMF só como comparação.

## Convenções de código

- R com `data.table`; o PLINK é chamado de dentro do R (`rodar_plink()` em
  `funcoes_popstructure.R`), para que todo o fluxo fique num script só.
- Rodar os scripts a partir da raiz do projeto.
- Scripts organizados em seções numeradas com comentários explicativos, porque
  vão virar um `.qmd`.
- Funções compartilhadas ficam em `funcoes_popstructure.R`. Não duplicar.
- PLINK sempre com `--chr-set 10 no-xy` (milho).
- Nos parâmetros do STRUCTURE, números sempre sem notação científica
  (`format(x, scientific = FALSE)`): `1e+05` seria lido como 1.
- `RANDOMIZE 0` e semente explícita (`-D`) em toda execução.

## Cuidados

- `data/` e `results/` não são versionados. Os dados vêm do GEO (ver README).
- Os caminhos de `plink` e `structure` em `funcoes_popstructure.R` (e na
  seção 0 do `bio640_popstructure.R`) são da máquina do autor do repositório; ajuste para a sua.
- O STRUCTURE roda em 1 núcleo por execução e é lento: estime o tempo antes de
  rodar (veja os números no README) e prefira subamostras para testes.
- Não reverter as decisões listadas no README sem discutir: cada uma tem um
  motivo registrado.
