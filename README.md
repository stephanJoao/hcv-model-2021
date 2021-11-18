# hcv_model

Modelo matemático que representa o tratamento contra hepatite C.

Na pasta differential_evolution temos a pasta depend, que possui a função de custo(em python), o modelo(em C++) e um arquivo makefile. O makefile é chamado pela função de custo(em python) e é responsável por executar o modelo(em C++). No arquivo parametros_DE é escrito, pela cost.py, os parâmetros selecionados pela evolução diferencial e esses parâmetros são lidos pelo HCV_model.h. Após terminar de executar, o modelo escreve os resultados no arquivo saida e retorna o comando para a função cost. Quando a DE termina de executar para um determinado paciente ela escreve o conjunto de parâmetros que mais bem se ajustou a curva no arquivo relatorio_DE e salva o gráfico na pasta figs. O nome de cada arquivo gerado é composto pelos parâmetros utilizados pela DE. 

O arquivo plotDados_e_Modelo.py executa o modelo e compara com os dados experimentais. Este arquivo não tem acesso automático aos parâmetros gerados pelo differential_evolution.py portanto após executar o differential_evolution.py, por exemplo, para o pacienteB06 é necessário pegar os parâmetros em relatorio_DE e escrever eles na sua respectiva estrutura dentro do arquivo plotDados_e_Modelo.py. 

Qualquer mudança feita nos parâmetros que serão ajustados pelo differential_evolution.py deve ser feita nos arquivos cost.py, HCV_model.h e plotDados_e_Modelo.py de forma manual.

# Como usar

De dentro da pasta depend, execute o código differential_evolution.py. É necessário rodar de dentro da pasta depend pois o código em python faz uma chamada de sistema pelo comando "make run". Para que essa chamada dê certo é necessário ter um arquivo makefile, e ele está na pasta depend.

No terminal:

cd differential_evolution/depend

python3 ../differential_evolution.py
