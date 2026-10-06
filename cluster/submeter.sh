#!/bin/bash
#
# Submete ao PBS uma execução do STRUCTURE para cada linha de execucoes.tsv.
# Rodar no head node, de dentro da pasta do pacote:
#
#   bash submeter.sh            # submete o que falta, respeitando o limite
#   bash submeter.sh dent       # só um conjunto
#
# Limites do cluster da UFV (dct.ufv.br/politicas):
#   - fila qtime: até 48 jobs por usuário;
#   - até 96 CPUs por usuário, somando jobs em fila e em execução.
# Por isso o script submete no máximo MAX_JOBS jobs ativos. Quando alguns
# terminarem, rode de novo: execuções concluídas, em fila ou rodando são puladas.
#
# Colunas de execucoes.tsv:
#   conjunto  K  rep  semente  n_ind  n_loci  walltime

set -euo pipefail
filtro="${1:-}"
MAX_JOBS=45

ativos=$(qselect -u "$USER" -s QRH | wc -l)
echo "Jobs ativos agora: $ativos (limite usado por este script: $MAX_JOBS)"

while IFS=$'\t' read -r conj k rep semente nind nloc walltime; do
  [ -n "$filtro" ] && [ "$conj" != "$filtro" ] && continue

  nome="str_${conj}_K${k}_r${rep}"
  saida="resultados/${conj}/${conj}_K$(printf "%02d" "$k")_r${rep}_f"

  [ -f "$saida" ] && continue                                   # concluída
  [ -n "$(qselect -u "$USER" -N "$nome" -s QRH)" ] && continue  # em fila/rodando

  if [ "$ativos" -ge "$MAX_JOBS" ]; then
    echo "Limite de $MAX_JOBS jobs atingido. Rode este script de novo mais tarde."
    break
  fi

  id=$(qsub -N "$nome" \
            -l walltime="${walltime}" \
            -v CONJ="$conj",K="$k",REP="$rep",SEMENTE="$semente",NIND="$nind",NLOC="$nloc" \
            structure_job.pbs)
  echo "$id  $nome"
  ativos=$((ativos + 1))
done < <(tail -n +2 execucoes.tsv)

restantes=0
while IFS=$'\t' read -r conj k rep _; do
  [ -f "resultados/${conj}/${conj}_K$(printf "%02d" "$k")_r${rep}_f" ] || restantes=$((restantes + 1))
done < <(tail -n +2 execucoes.tsv)
echo "Execuções ainda não concluídas: $restantes"
