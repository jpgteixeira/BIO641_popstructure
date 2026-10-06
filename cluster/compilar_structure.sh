#!/bin/bash
#
# Compila o STRUCTURE 2.3.4 (versão de linha de comando) a partir do código-fonte
# e coloca o executável em bin/. Rodar uma vez no head node, dentro do pacote.
#
# Por que compilar? O executável distribuído para Linux é de 32 bits e depende
# das bibliotecas i386 do sistema, que podem não existir nos nós do cluster.
# Compilado lá, ele fica com 64 bits e usa as bibliotecas do próprio cluster.
#
# Compilar é uma tarefa curta (segundos) e não é análise, por isso pode ser
# feita no head node.

set -euo pipefail

URL=https://web.stanford.edu/group/pritchardlab/structure_software/release_versions/v2.3.4/structure_kernel_source.tar.gz

source /etc/profile.d/modules.sh
module load gcc 2>/dev/null || echo "Aviso: módulo gcc não carregado; usando o gcc do sistema."
gcc --version | head -1

mkdir -p bin fonte
cd fonte
wget -q -O structure_kernel_source.tar.gz "$URL"
tar -xzf structure_kernel_source.tar.gz
cd structure_kernel_src

# Compiladores gcc >= 10 recusam variáveis globais definidas em mais de um
# arquivo (padrão -fno-common); o código do STRUCTURE é antigo e precisa de
# -fcommon.
make clean >/dev/null 2>&1 || true
make CC="gcc -fcommon"

cp structure ../../bin/structure
cd ../..
file bin/structure
echo "STRUCTURE compilado em bin/structure"
