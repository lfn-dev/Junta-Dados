# Junta-Dados
Programa para juntar dados climaticos em pontos de grade (binário).
No laboratório, foi utilziado para juntar dados do GPCC com dados observados.

---

Este é um arquivo com o máximo de informações sobre o programa reunidas, para uma referência rápida ver `Resumo.md` e para uma explicação matemática dos métodos utilizados: `Interpolação.tex` ou `Interpolação.pdf`.

# Funcionamento

O programa lê dois arquivos (.ctl) no formato [GrADS Data Descriptor File](http://cola.gmu.edu/grads/gadoc/descriptorfile.html),
e cria um novo arquivo binário com seu respectivo .ctl que corresponde à junção de ambos.

## Compilação
Para gerar o executável, basta utilizar o comando `make` no diretório `src/` que contém o arquivo `Makefile`.

    make all

De forma alternativa, para gerar um executável de depuração (veja seção [`Modo Depuração`](#Debug) ) utilizar:

    make debug

Caso os arquivos binários foram escritos com precisão dupla (`double`) ao invés de `float`, é preciso gerar o executável compatível com o comando:

    make double

Para excluir os arquivos objetos gerados durante a compilção:

    make clean

Para excluir tanto os arquivos objeto quanto os executáveis:

    make purge
    
## Executando
Para executar o programa com valores padrão:

    ./compose arquivo_princial.ctl arquivo_secundario.ctl prefixo_saída
    
É possível modificar os valores padrão com opções na linha de comando. Para ver todas as opções disponĩveis:

    ./compose -h
    ./compose --help

### Formatação dos arquivos de Entrada
Para executar a junção de dados é necessário formatar os arquivos `.ctl` de entrada.
Os arquivos precisam ser de apenas 1 variável e 1 nível, contendo os seguintes campos na ordem:

> dset ^file.bin
> title Título
> undef `<value>`
> xdef `<int_value>` linear `<value>` `<value>`
> ydef `<int_value>` linear `<value>` `<value>`
> zdef `<int_value>` linear 1 1
> tdef `<int_value>` linear `<date>` `<step>`
> vars 1
> <name>  0  777.7  cx_chuva
> endvars

Onde `<int_value>` são valores inteiros e `<value>` são valores em ponto flutuante. Veja [GrADS Data Descriptor File](http://cola.gmu.edu/grads/gadoc/descriptorfile.html) para mais informações.

