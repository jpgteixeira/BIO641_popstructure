# STRUCTURE no cluster da UFV — conjunto completo

Guia para rodar a análise final (todas as 2.269 amostras que passaram no controle de qualidade) no cluster da UFV.
A preparação dos dados e a análise dos resultados continuam no computador
local; o cluster só executa o STRUCTURE.

```
LOCAL                         CLUSTER (head node)              LOCAL
bio640_cluster.R        →     compilar_structure.sh      →     bio640_popstructure.R
  gera o pacote               submeter.sh (qsub)                 Parte 3: pophelper
  (dados + parâmetros)        jobs nos nós de cálculo            (ΔK, gráficos)
```

Referências: <https://dct.ufv.br/uso-do-ambiente/> e <https://dct.ufv.br/politicas/>.

## 1. O que roda no cluster

| Conjunto | Indivíduos | SNPs (após LD) | K | Réplicas | Jobs | Tempo por job (estimado)* |
|---|---|---|---|---|---|---|
| global | 2.269 | 2.226 | 1–8 | 5 | 40 | até ~15 h |
| dent | 1.014 | 764 | 1–13 | 5 | 65 | até ~4 h |
| flint | 1.258 | 799 | 1–14 | 5 | 70 | até ~6 h |
| **Total** | | | | | **175** | |

\* MCMC com 50.000 de burn-in + 100.000 iterações, medido na máquina local.
Os nós do cluster podem ser mais lentos; por isso o walltime pedido tem folga
(48 h no global e 24 h nos *pools*), dentro do máximo de 72 h da fila `qtime`.
As faixas de K serão ajustadas depois da exploração local.

Cada job é **uma execução do STRUCTURE** (um conjunto, um K, uma réplica) e
usa **1 processador e ~100 MB de RAM** (o STRUCTURE não é paralelo).

## 2. Limites do cluster que afetam o plano

| Limite (por usuário) | Valor | Consequência |
|---|---|---|
| Jobs na fila `qtime` | 48 | Não dá para submeter os 175 de uma vez |
| CPUs em fila + execução | 96 | No máximo ~96 execuções simultâneas |
| Walltime na `qtime` | 72 h | Cabe a execução mais longa com folga |
| Acesso interativo aos nós | proibido | Tudo via `qsub` |
| Processadores usados > pedidos | job encerrado | Pedimos 1 CPU e o STRUCTURE usa 1 |

O `submeter.sh` mantém no máximo 45 jobs ativos. Quando alguns terminarem,
basta rodá-lo de novo: ele pula as execuções concluídas, em fila ou rodando.

## 3. Conteúdo do pacote

Gerado pelo `bio640_cluster.R` na pasta `cluster/pacote/` (27 MB). Ele parte dos
dados após o controle de qualidade (seções 1–5 do `bio640_popstructure.R`) e usa as
funções de `funcoes_popstructure.R`, as mesmas da análise local:

```bash
Rscript bio640_cluster.R
```


```
pacote/
├── bin/                     (vazio; o STRUCTURE é compilado lá)
├── dados/
│   ├── global.recode.strct_in
│   ├── dent.recode.strct_in
│   ├── flint.recode.strct_in
│   └── *_ordem.tsv          ordem das amostras em cada arquivo (para os gráficos)
├── params/
│   ├── mainparams           burn-in, iterações, formato do arquivo
│   └── extraparams          modelo admixture, FREQSCORR, RANDOMIZE 0
├── execucoes.tsv            uma linha por job: conjunto, K, réplica, semente, ...
├── structure_job.pbs        script PBS de uma execução
├── submeter.sh              submete os jobs respeitando os limites
├── compilar_structure.sh    compila o STRUCTURE no cluster
└── resultados/              saídas do STRUCTURE (preenchida pelos jobs)
```

Antes de enviar, edite a linha `#PBS -M <seu_email>@ufv.br` do
`structure_job.pbs` com o seu e-mail institucional. O script pede aviso por
e-mail só quando um job é abortado (`-m a`), para não receber 175 mensagens.

## 4. Passo a passo

O endereço do cluster e o seu usuário estão no
[Sistema Cluster](https://www3.dti.ufv.br/dti/cluster/autenticacao)
→ sua conta → **Detalhes** → *Detalhes da Solicitação*. Abaixo,
`<usuario>` e `<host>` devem ser trocados por esses valores.

### 4.1 Enviar o pacote (máquina local)

```bash
rsync -av cluster/pacote/ <usuario>@<host>:~/bio641_structure/
```

### 4.2 Preparar o STRUCTURE (head node)

```bash
ssh <usuario>@<host>
cd ~/bio641_structure
module avail structure          # se existir um módulo, dá para usá-lo
bash compilar_structure.sh      # senão, compila (alguns segundos)
```

O executável de Linux distribuído no site é de 32 bits e pode não rodar nos
nós; compilado no cluster ele fica com 64 bits.

### 4.3 Testar com um conjunto pequeno

```bash
bash submeter.sh dent           # 45 primeiros jobs do dent
qstat -anu <usuario>            # devem aparecer como Q (fila) e depois R
jinfo <job_id>                  # CPU e RAM usados vs. pedidos
```

Quando o primeiro job terminar, confira `resultados/dent/*_f` e o `.log`
correspondente. Se tudo estiver certo, siga.

### 4.4 Submeter o restante

```bash
bash submeter.sh                # repetir até "Execuções ainda não concluídas: 0"
```

Comandos úteis: `qstat -anu <usuario>` (seus jobs), `qdel <job_id>` (cancelar),
`qstat -f <job_id>` (detalhes). Os arquivos `str_*.o<id>` e `str_*.e<id>` são
a saída padrão e de erro de cada job.

### 4.5 Trazer os resultados (máquina local)

```bash
rsync -av <usuario>@<host>:~/bio641_structure/resultados/ results/structure/producao/
```

A Parte 3 do `bio640_popstructure.R` lê essa pasta com o pophelper.

## 5. Boas práticas (políticas da UFV)

- Não rodar análises no head node; só compilar, editar e submeter.
- Backup é responsabilidade do usuário: traga os resultados para a máquina
  local assim que terminarem.
- Apague do cluster os arquivos que não forem mais necessários.
