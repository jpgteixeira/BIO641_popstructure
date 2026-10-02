# Estrutura dos dados — BIO641, Estrutura de Populações

Visão geral do banco de dados que vamos usar para a análise com o STRUCTURE.

## Origem

Os dados vêm do projeto europeu **CornFed**, com populações multiparentais de milho (*Zea mays* L.):

- **Bauer et al. 2013**, *Genome Biology* 14:R103 — desenho das populações e genotipagem (estudo de taxa de recombinação).
- **Lehermeier et al. 2014**, *Genetics* 198(1):3 — predição genômica; fornece os fenótipos (Supporting File S1).
- Genótipos públicos no NCBI GEO: **GSE50558**.

## Artigos de referência do trabalho

| Artigo | Assunto | Papel no trabalho |
|---|---|---|
| Caniato et al. 2011, *PLoS ONE* | Estrutura populacional e tolerância ao alumínio em sorgo (38 SSRs) | Modelo de uso do STRUCTURE: K = 1–13, burn-in 100.000, MCMC 1.000.000, 5 réplicas, modelo *admixture* com frequências correlacionadas, ΔK de Evanno |
| Yu et al. 2006, *Nature Genetics* 38:203 | Modelo misto unificado Q + K | Base conceitual: Q (estrutura, vinda do STRUCTURE) × K (parentesco) |
| Lehermeier et al. 2014, *Genetics* | Populações multiparentais de milho | Fonte dos dados; usou PCA, parentesco e LD, mas **não** STRUCTURE |

## Desenho experimental

O desenho é do tipo NAM (*nested association mapping*), com dois *pools* heteróticos:

| | Dent | Flint |
|---|---|---|
| Linha central | **F353** | **UH007** |
| Prefixo das famílias | CFD | CFF |
| Testador (fenótipo) | UH007 | F353 |

- Cada família vem do cruzamento **linha central × linha fundadora**.
- Das F1 foram geradas **linhas duplo-haploides (DH)**, que são totalmente homozigotas.
- Em média, cada linha DH tem **50% do genoma da linha central e 50% da fundadora**.
- Assim, as famílias de um mesmo *pool* são **meio-irmãs**, porque todas compartilham a linha central.

### Famílias e pedigrees (segundo os arquivos de fenótipo)

| Dent | Cruzamento | Flint | Cruzamento |
|---|---|---|---|
| CFD02 | F353 × B73 | CFF02 | UH007 × **B73** (fundadora dent!) |
| CFD03 | F353 × D06 | CFF03 | UH007 × D152 |
| CFD04 | F353 × D09 | CFF04 | UH007 × EC49A |
| CFD05 | F353 × EC169 | CFF05 | UH007 × EP44 |
| CFD06 | F353 × F252 | CFF06 | UH007 × EZ5 |
| CFD07 | F353 × F618 | CFF07 | UH007 × F03802 |
| CFD08 | F353 × F98902 | CFF08 | UH007 × F2 |
| CFD09 | F353 × Mo17 | CFF09 | UH007 × F283 |
| CFD10 | F353 × UH250 | CFF10 | UH007 × F64 |
| CFD11 | F353 × UH304 | CFF12 | UH007 × UH006 |
| CFD12 | F353 × W117 | CFF13 | UH007 × UH009 |
| | | CFF15 | UH007 × DK105 |

### Linhas parentais (Tabela S1 de Bauer et al. 2013 — as mesmas 23 do GEO)

- **Dent (11):** B73, D06, D09, EC169, F252, **F353** (central), F618, Mo17, UH250, UH304, W117
- **Flint (12):** D152, DK105, EC49A, EP44, EZ5, F03802, F2, F64, F283, UH006, **UH007** (central), UH009

## Arquivos

```
data/
├── GSE50558_Parental_matrix_GEO.txt.gz      14 MB   23 parentais   — genótipos
├── GSE50558_CFD_matrix_GEO.txt.gz          530 MB   1.005 DH dent  — genótipos
├── GSE50558_CFF_matrix_GEO.txt.gz          657 MB   1.262 DH flint — genótipos
├── GSE50558_*_intensities_GEO.txt.gz       ~1,7 GB  intensidades brutas (não necessárias)
├── gb-2013-14-9-r103-S1.xlsx                        Tabela S1 de Bauer et al. 2013 (parentais)
├── gb-2013-14-9-r103-S4.xlsx                        Suplemento S4 de Bauer et al. 2013
└── genetics.114.161943-17/
    ├── SupportingFileS1_Readme.txt
    ├── PhenotypicDataDent.csv                       4.800 parcelas
    └── PhenotypicDataFlint.csv                      7.680 parcelas
```

### Arquivos de genótipo (`*_matrix_GEO.txt.gz`)

- Chip **Illumina MaizeSNP50**, com **56.110 SNPs** (uma linha por SNP).
- A primeira coluna (`ID_REF`) é o nome do SNP (ex.: `abph1.15`). **Não há posição cromossômica nesses arquivos**; ela está na anotação da plataforma (GPL) no GEO ou no *manifest* do chip.
- Para cada amostra existem 5 colunas:

| Coluna | Conteúdo |
|---|---|
| `<amostra>.GType` | Genótipo: `AA`, `AB`, `BB` ou `NC` (sem chamada) |
| `<amostra>.Top Alleles` | Alelos em nucleotídeos (fita TOP da Illumina), ex.: `GG`; `--` = sem chamada |
| `<amostra>.Score` | GenTrain score do SNP (qualidade do *cluster*; igual para todas as amostras) |
| `<amostra>.Theta` | Ângulo normalizado da intensidade (0 ≈ AA, 1 ≈ BB) |
| `<amostra>.R` | Intensidade total normalizada |

Nos 23 parentais, cerca de **16% das chamadas são `NC`** e apenas **0,8% são heterozigotas (`AB`)**, como se espera de linhas endogâmicas.

### Amostras genotipadas por família (GEO)

| Dent | n | Flint | n |
|---|---|---|---|
| CFD01 | 86 | CFF01 | 99 |
| CFD02 | 73 | CFF02 | 120 |
| CFD03 | 103 | CFF03 | 112 |
| CFD04 | 105 | CFF04 | 53 |
| CFD05 | 77 | CFF05 | 34 |
| CFD06 | 105 | CFF06 | 50 |
| CFD07 | 108 | CFF07 | 129 |
| CFD09 | 63 | CFF08 | 77 |
| CFD10 | 99 | CFF09 | 134 |
| CFD11 | 86 | CFF10 | 108 |
| CFD12 | 100 | CFF12 | 114 |
| | | CFF13 | 117 |
| | | CFF15 | 115 |
| **Total** | **1.005** | **Total** | **1.262** |

### Arquivos de fenótipo (`PhenotypicData*.csv`)

Uma linha por parcela de campo (*testcross*). As colunas principais são `Genotype` (código da linha DH, igual ao GEO), `Population` (família), `Pedigree`, `Tester`, `LOC` (local), `Rep` e as características DMY, DMC, PH, DtTAS, DtSILK e NBPL.

**Para o STRUCTURE, o fenótipo não é usado.** Esses arquivos servem só para mapear família → fundadora na interpretação dos resultados.

### Controle de qualidade aplicado no artigo (Lehermeier et al. 2014)

Foram removidos os SNPs com GenTrain < 0,7, *call frequency* < 0,9, MAF < 0,01 ou mais de 10% de dados faltantes. Restaram **34.116 SNPs**. O artigo analisou 841 linhas DH dent e 811 flint.

## Pontos de atenção

1. **CFF02 (UH007 × B73)** — uma fundadora dent dentro do painel flint. É um controle natural: o STRUCTURE deve atribuir ~50% de ancestralidade dent a essas linhas.
2. **Divergências entre GEO e fenótipos**:
   - **CFD01** e **CFF01** estão no GEO, mas não nos fenótipos, então o pedigree delas não consta nos nossos arquivos.
   - **CFD08** (F353 × F98902) está nos fenótipos, mas não no GEO, e F98902 não está entre os parentais genotipados.
3. **Premissas do STRUCTURE** (equilíbrio de Hardy-Weinberg e equilíbrio de ligação dentro dos *clusters*) são violadas:
   - As linhas DH não vêm de cruzamento aleatório e são homozigotas.
   - As famílias são aparentadas (meio-irmãs).
   - Houve só uma rodada de recombinação, então o LD se estende por blocos longos.
4. **Ploidia**: as linhas DH podem ser codificadas como diploides homozigotas ou como haploides (`PLOIDY=1`).
5. **Custo computacional**: ~2.300 indivíduos × ~34 mil SNPs é inviável no STRUCTURE com os parâmetros do Caniato. Será preciso selecionar SNPs (espaçados no genoma) e possivelmente subamostrar linhas.

## Expectativa de resultado

| Nível | O que esperamos ver |
|---|---|
| Dent × Flint | K = 2 separando os *pools* (diferenciação histórica, estrutura "clássica") |
| Dentro de cada *pool* | *Clusters* que reproduzem as famílias/fundadoras (estrutura de pedigree, ou seja, parentesco) |
| CFF02 | Ancestralidade mista, ~50% dent |

A pergunta central do trabalho: **o que o STRUCTURE detecta quando a "população" foi construída por melhoristas?**
