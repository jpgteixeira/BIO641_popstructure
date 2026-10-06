# BIO641 — Estrutura populacional em milho com o STRUCTURE

Trabalho da disciplina **BIO641 – Genética de Populações (UFV)**. O objetivo é
calcular a estrutura populacional dos painéis de milho **CornFed Dent e Flint**
(GEO [GSE50558](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE50558))
com o **STRUCTURE**, como pedido pelo professor, e apresentar o tema
"estrutura de populações" com uma parte prática.

> **Pergunta central:** o que o STRUCTURE detecta quando a "população" foi
> construída por melhoristas?

## Resumo do material

- **Desenho do tipo NAM** com dois *pools* heteróticos. No **Dent**, a linha
  central F353 foi cruzada com 10–11 fundadoras (famílias CFD). No **Flint**,
  a central UH007 foi cruzada com 11–12 fundadoras (famílias CFF).
- Cada família é formada por **linhas duplo-haploides (DH)**, que são
  homozigotas e têm em média 50% do genoma da linha central e 50% da
  fundadora.
- **Genotipagem** com o chip Illumina MaizeSNP50 (56.110 SNPs), em 2.267 DH e
  nos 23 parentais.
- **Controles de mistura conhecida:**
  - **CFD01 e CFF01** são F353 × UH007 (dent × flint);
  - **CFF02** é UH007 × B73, com uma fundadora dent no painel flint.

  Na análise global, essas famílias devem aparecer com cerca de 50% de cada
  *pool*.

Detalhes em [docs/estrutura_dos_dados.md](docs/estrutura_dos_dados.md).

### Artigos de referência

| Artigo | Papel |
|---|---|
| Caniato et al. 2011, *PLoS ONE* (sorgo, SSR) | Modelo de uso do STRUCTURE: modelo *admixture*, frequências correlacionadas, ΔK de Evanno |
| Yu et al. 2006, *Nature Genetics* 38:203 | Base conceitual: matriz Q (estrutura) × matriz K (parentesco) |
| Lehermeier et al. 2014, *Genetics* 198:3 | Origem dos dados; critérios de controle de qualidade |
| Bauer et al. 2013, *Genome Biology* 14:R103 | Desenho das populações (Tabela S1: parentais) |

## Fluxo de análise

```
R (data.table)        PLINK 1.9                      STRUCTURE 2.3.4         R (pophelper)
ler GEO          →    controle de qualidade     →    MCMC para cada K   →    ΔK de Evanno,
montar .tped          poda por LD                    e cada réplica          alinhamento das
                      exportar formato STRUCTURE                             réplicas, gráficos
```

**Três análises hierárquicas:** **global** (dent + flint + parentais),
**dent** e **flint**. O ΔK tende a captar o nível mais alto da estrutura
(dent × flint), então os *pools* precisam de análises separadas para revelar a
estrutura interna de cada um.

**Duas fases:**
1. **Exploração local**: 30 DH por família + parentais (743 amostras),
   10 mil de burn-in + 20 mil iterações, 3 réplicas.
2. **Produção no cluster da UFV**: todas as 2.269 amostras, 50 mil de burn-in +
   100 mil iterações, 5 réplicas. Ver [docs/cluster_ufv.md](docs/cluster_ufv.md).

## Organização do repositório

```
bio640_popstructure.R    script principal (análise local)
                           Parte 1 (seções 1-5): GEO -> PLINK -> controle de qualidade
                           Parte 2 (seções 6-10): subamostra, 3 conjuntos, poda por LD,
                                                  parâmetros e execução do STRUCTURE
                           Parte 3 (a fazer): pophelper (ΔK, gráficos de barras)
bio640_cluster.R         monta o pacote para o cluster a partir dos dados após o CQ
funcoes_popstructure.R   funções compartilhadas (PLINK, parâmetros, conjuntos)
                           -> contém os CAMINHOS do plink e do structure
cluster/                 modelos dos scripts PBS (OpenPBS/UFV) e de compilação
docs/
  estrutura_dos_dados.md   descrição dos dados e pontos de atenção
  cluster_ufv.md           guia para rodar no cluster da UFV
data/                    (não versionado) dados brutos e processados
results/                 (não versionado) saídas do PLINK e do STRUCTURE
```

Os nomes `bio640_*` são intencionais (ainda que a disciplina seja BIO641).

## Como reproduzir

### 1. Programas

| Programa | Versão usada | Observação |
|---|---|---|
| R | 4.6.1 | pacotes `data.table`, `ggplot2`, `pophelper` (`remotes::install_github("royfrancis/pophelper")`) |
| PLINK | 1.9 | sempre com `--chr-set 10 no-xy` (milho) |
| STRUCTURE | 2.3.4 (console) | [site do Pritchard Lab](https://web.stanford.edu/group/pritchardlab/structure.html) |

**Ajuste os caminhos** de `plink` e `structure` no início de
`funcoes_popstructure.R` (e na seção 0 do `bio640_popstructure.R`). Os
arquivos `mainparams`/`extraparams` são gerados a partir dos modelos que
acompanham o STRUCTURE, na mesma pasta do executável.

### 2. Dados (pasta `data/`, não versionada)

| Arquivo | Origem |
|---|---|
| `GSE50558_series_matrix.txt.gz` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50558/matrix/> (12 MB, **fonte principal dos genótipos**) |
| `GPL17677_annotation.txt.gz` | tabela da plataforma, com cromossomo e posição: `curl -sL "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL17677&targ=self&form=text&view=data" \| gzip > data/GPL17677_annotation.txt.gz` |
| `GSE50558_Parental_matrix_GEO.txt.gz` | <https://ftp.ncbi.nlm.nih.gov/geo/series/GSE50nnn/GSE50558/suppl/> (usado só para o GenTrain score) |
| `genetics.114.161943-17/` | Supporting File S1 de Lehermeier et al. 2014 (fenótipos e pedigrees) |
| `gb-2013-14-9-r103-S1.xlsx`, `-S4.xlsx` | Arquivos adicionais de Bauer et al. 2013 |

Os arquivos `*_CFD_matrix`, `*_CFF_matrix` e `*_intensities` do
*supplementary* **não são necessários**: o Series Matrix tem os mesmos
genótipos (conferido: 0 diferenças nos parentais).

### 3. Execução

```bash
Rscript bio640_popstructure.R    # Partes 1 e 2; a exploração leva ~2 h com 16 núcleos
Rscript bio640_cluster.R         # gera cluster/pacote/ para o cluster
```

Rodar da raiz do projeto. As execuções do STRUCTURE já concluídas são
puladas, então dá para interromper e retomar.

## Decisões tomadas (e por quê)

| Decisão | Motivo |
|---|---|
| Genótipos do **Series Matrix**, posições da **GPL17677** | Um arquivo de 12 MB com as 2.290 amostras e os metadados |
| **EC169 reclassificado como Dent** | No GEO consta Flint, mas é fundador da CFD05 e Dent na Tabela S1 |
| Controle de qualidade como em Lehermeier 2014: GenTrain ≥ 0,7, call rate ≥ 0,9, MAF ≥ 0,01 | Comparabilidade com o artigo; restaram 33.430 SNPs |
| SNPs dos cromossomos 0 e 99 removidos | Sem posição no genoma, não dá para controlar o LD |
| **Heterozigotos → dado faltante** | Todas as linhas são endogâmicas; AB é erro ou heterozigose residual |
| Filtro de dados faltantes por amostra (> 10%) **só nas DH** | O parental UH304 tem ~12% de faltantes, mas é necessário para rotular os clusters |
| 21 DH removidas | Mais de 10% de faltantes depois da troca dos AB; várias com numeração 3xx, possivelmente não são DH verdadeiras |
| **Diploide** (`PLOIDY 2`) | Prática comum para linhas endogâmicas. Ressalva: cada DH conta como duas cópias "independentes", o que deixa as estimativas de Q mais confiantes do que deveriam |
| MAF, call rate e **poda por LD recalculados em cada conjunto** (r² > 0,2, janelas de 50 SNPs, passo 5) | Um SNP polimórfico entre *pools* pode ser fixo dentro de um *pool* |
| Modelo *admixture*, `FREQSCORR 1`, `INFERALPHA 1` | Como em Caniato et al. 2011 |
| `RANDOMIZE 0` e sementes fixas | O padrão (`RANDOMIZE 1`) ignora a semente e impede reproduzir os resultados |
| Subamostra de 30 DH por família na exploração | O tempo do STRUCTURE cresce com o número de indivíduos; mantém todas as 24 famílias e os controles |
| Menos iterações que Caniato (1,1 milhão) | Com milhares de SNPs a cadeia converge mais rápido; a convergência será verificada |

### Números da exploração (subamostra)

| Conjunto | Indivíduos | SNPs após a poda | K testado |
|---|---|---|---|
| global | 743 | 2.255 | 1–6 |
| dent | 342 | 780 | 1–12 |
| flint | 404 | 720 | 1–14 |

Dentro de cada *pool* sobram poucos SNPs porque as DH passaram por uma única
meiose: o LD se estende por blocos longos.

## Situação atual

- [x] Dados entendidos e documentados
- [x] Controle de qualidade (Parte 1)
- [x] Subamostra, conjuntos, poda por LD e exportação (Parte 2)
- [ ] **Exploração local do STRUCTURE em andamento** (96 execuções)
- [ ] Parte 3: pophelper. Curva ln P(D|K), ΔK de Evanno, alinhamento das réplicas, gráficos de barras ordenados por família, verificação de convergência
- [ ] Refatorar o `bio640_popstructure.R` para usar o `funcoes_popstructure.R` (hoje parte do código está repetido)
- [ ] Ajustar as faixas de K e rodar a produção no cluster da UFV
- [ ] Análises complementares: PCA (`plink --pca`), F_ST entre *pools* e famílias (efeito Wahlund: H_T × H_S)
- [ ] `.qmd` explicando cada etapa, para o relatório e a apresentação

## Resultados esperados

| Nível | Expectativa |
|---|---|
| Dent × Flint | K = 2 separando os *pools* (diferenciação histórica, F_ST alto) |
| Dentro de cada *pool* | *Clusters* que reproduzem as famílias e fundadoras, ou seja, estrutura de pedigree (parentesco), e não populações no sentido clássico |
| CFD01, CFF01, CFF02 | Ancestralidade mista, cerca de 50% |
