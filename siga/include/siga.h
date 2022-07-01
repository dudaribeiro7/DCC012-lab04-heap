#ifndef SIGA_H
#define SIGA_H

#include <iostream>
#include <fstream>
#include <string>

#include "config.h"
#include "estudante.h"

using namespace std;


typedef enum {BUBBLESORT, INSERTIONSORT, SELECTIONSORT, MERGESORT, QUICKSORT} sorting_method;

class Siga
{
    
    public:
        //Construtor: Cria ou abre um dataset binário.
        Siga(string dataset_name);

        // Cadastra um estudante no dataset.
        void CadastraEstudante(Estudante est);

        // Pesquisa um estudante pelo número de matrícula.
        int  PesquisaPorMatricula(int matricula);

        int PesquisaPorNome(const char* nome);

        // Busca por um estudante pela matricula e imprime seus dados.
        void ImprimeEstudantePorMatricula(int matricula);

        // Busca por um estudante pelo nome e altera seus dados.
        void AlteraCadastroEstudante(Estudante est);

        // Imprime todos os estudantes do dataset em um arquivo texto.
        void ExportDataSetToCSV(string arquivo_txt);

        // Obtem o numero de estudantes armazenados no dataset.
        int  ObterNumeroEstudantesArmazenados();

        // Ler um arquivo CSV e importa os dados para o dataset.
        int ImportCSVData(string arquivo_csv);

        // Ordena os estudantes do dataset por nome e imprime em um arquivo texto.
        void SalvarCSVOrdenadoPorNome(string arquivo_txt, sorting_method method);

        // Ordena os estudantes do dataset por curs, seguido peo nome e imprime em um arquivo texto.
        void SalvarCSVOrdenadoPorCursoENome(string arquivo_txt);


        void clean();

        // Destrutor: fecha o dataset binário.
        ~Siga();

        // Função auxiliar para Extrair um vetor de estudante do dataset.
        Estudante* ExtractAllEstudantsFromDataSet();
         
    private:

        void ExportaVectorToCSV(Estudante* estudantes, int n_estudantes, string arquivo_txt);
    
        // atributos
        string   nome_dataset;
        fstream  file_stream;
        int      n_estudantes;
};


#endif /* SIGA_H */
