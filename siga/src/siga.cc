#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>

using namespace std;

#include "siga.h"
#include "sort.h"


Siga::Siga(string arquivo_dados_estudante)
{

    this->nome_dataset = arquivo_dados_estudante;
    this->file_stream.open(this->nome_dataset, ios::in | ios::out | ios::binary );

    if(this->file_stream.is_open())
    {
        cout << "SIGA: Inicializado com sucesso" << endl;
    }
    else
    {
        // arquivo não existe, cria um novo dataset.
        this->file_stream.open(this->nome_dataset, ios::out | ios::trunc);
        this->file_stream.close();

        this->file_stream.open(this->nome_dataset, ios::in | ios::out | ios::binary );
        if(!this->file_stream.is_open()) {
            cout << "SIGA: Erro ao abrir arquivo" << endl;
            exit(-1);
        }
        
    }

     this->n_estudantes = 0;

    if ( this->file_stream.is_open() ){
        this->file_stream.seekg(0, this->file_stream.end);
        int length = this->file_stream.tellg();
        this->file_stream.seekg(0, this->file_stream.beg);

        this->n_estudantes = length / ESTUDANTE_SIZE;
    }
    cout <<  this->n_estudantes << " registros de estudantes" << endl;
}


int  Siga::PesquisaPorMatricula(int matricula)
{
    //Posicionando o cursor no inicio do arquivo
    this->file_stream.seekg(0, this->file_stream.beg);

    //Declarando Estudante
    Estudante est = Estudante();

    for ( int i = 0; i < this->n_estudantes; i++){
        //Ler do arquivo binário
        est.LerDoArquivoBinario(this->file_stream);
        //Verifica se a matricula é a mesma
        if ( est.ObterMatricula() == matricula )
            return i;
    }

    //Retornando o cursor pro fim de leitura
    this->file_stream.seekg(0, this->file_stream.end);

    //Caso não encontrar a matrícula retorna -1
    return -1;
}

void Siga::clean()
{
    this->n_estudantes = 0;
}

int  Siga::PesquisaPorNome(const char* nome)
{
    //Posicionando o cursor no inicio do arquivo
    this->file_stream.seekg(0, this->file_stream.beg);

    //Declarando Estudante
    Estudante est = Estudante();

    for ( int i = 0; i < this->n_estudantes; i++){
        //Ler do arquivo binário
        est.LerDoArquivoBinario(this->file_stream);
        //Verifica se a matricula é a mesma
        if (strcmp(est.ObterNome(), nome) == 0)
            return i;
    }

    //Retornando o cursor pro fim de leitura
    this->file_stream.seekg(0, this->file_stream.end);

    //Caso não encontrar a matrícula retorna -1
    return -1;
}
        
void Siga::CadastraEstudante(Estudante est)
{
    //Verifica se est já está cadastrado
    if( PesquisaPorNome(est.ObterNome() ) != -1 )
        return;

    //Colocando o Cursor no Fim
    this->file_stream.seekp(0, this->file_stream.end);

    //Adicionando o aluno
    est.EscreverNoArquivoBinario(this->file_stream);

    //Incrementando o número de estudantes
    this->n_estudantes = this->n_estudantes + 1;

}
        
void Siga::ImprimeEstudantePorMatricula(int matricula)
{
    int index = PesquisaPorMatricula(matricula);
    if( index == -1){
        cout << "Estudante não cadastrado" << endl;
    } else {
        this->file_stream.seekg((index * ESTUDANTE_SIZE), this->file_stream.beg );
        Estudante est = Estudante();

        est.LerDoArquivoBinario(this->file_stream);

        est.Imprimir();
    }
}
        
void Siga::ExportDataSetToCSV(string arquivo_txt)
{
    //Cursor no inicio
    this->file_stream.seekg(0, this->file_stream.beg);

    //Verifica se o arquivo está vazio
    if(this->n_estudantes == 0)
        return;

    //Abre o arquivo
    fstream arq;
    arq.open(arquivo_txt, ios::in | ios::trunc);
    if(!arq.is_open()){
        cout << "FALHA AO ABRIR O ARQUIVO" << endl;
        exit(-1);
        //return;
    }

    //Cria um vetor de alunos
    Estudante estudante;
        
    //Salva a primeira linha
    arq << "Matricula;Nome;Ano de Ingresso;Curso;IRA" << endl;

    //Le o vetor de Estudantes
    for( int i = 0; i < this->n_estudantes; i++){
        estudante.LerDoArquivoBinario(this->file_stream);
        arq << estudante.ObterMatricula() << ";";
        arq << estudante.ObterNome() << ";";
        arq << estudante.ObterAnoIngresso() << ";";
        arq << estudante.ObterCurso() << ";";
        arq << estudante.ObterIRA() << endl;
    }

    //Fechando o arquivo
    arq.close();


}
        
        
void Siga::AlteraCadastroEstudante(Estudante est)
{
    //Pesquisa o aluno
    int index = PesquisaPorMatricula( est.ObterMatricula() );

    //Caso nao encontrado imprime e retorna
    if ( index == -1 ){
        cout << "Estudante não cadastrado" << endl;
        exit(-1);
        //return;
    }

    //Coloca o cursor no estudante que será reescrito
    this->file_stream.seekp(index * ESTUDANTE_SIZE, this->file_stream.beg);

    //Reescreve os dados dele
    est.EscreverNoArquivoBinario(this->file_stream);

    cout << "Dados Alterados" << endl;
}
        
Siga::~Siga()
{
    this->file_stream.close();
}

int Siga::ObterNumeroEstudantesArmazenados()
{
    return this->n_estudantes;
}


int Siga::ImportCSVData(string file)
{
     cout << "Importando dados do arquivo " << endl;
    int n_importados = this->n_estudantes;

    ifstream arq;
    arq.open(file, ios::in);
    if(!arq.is_open()){
        cout << "FALHA AO ABRIR O ARQUIVO" << endl;
        return -1;
    }
    string line;
    getline(arq, line); // header line
    while(getline(arq, line)){
        stringstream ss(line);
        string token;
        vector<string> tokens;
        while(getline(ss, token, ';')){
            tokens.push_back(token);
        }
        if(tokens.size() == 5){
            Estudante est;
            est.TrocarMatricula(stoi(tokens[0]));
            est.TrocarNome(tokens[1].c_str());
            est.TrocarAnoIngresso(stoi(tokens[2]));
            est.TrocarCurso(stoi(tokens[3]));
            est.TrocarIRA(stof(tokens[4]));
            this->CadastraEstudante(est);
        }
    }
    n_importados = this->n_estudantes - n_importados;
 
    // Imprima o numero de estudantes importados:
    cout << "Importacao concluida: " << n_importados << " novos alunos cadastrados" << endl;
    return n_importados;
     
}

void Siga::SalvarCSVOrdenadoPorNome(string arquivo_txt, sorting_method method)
{
    // Iremos aplicar a ordenação na memoria, para isso faca:
    // 1. Ler todos os dados do arquivo binário colocando-os em um vetor
    // 2. Ordene o vetor usando o método de ordenação passado como parametro:
    //    - Ordenar por nome
    Estudante* estudantes = ExtractAllEstudantsFromDataSet();

    switch(method){
        case BUBBLESORT:
            bubble_sort<Estudante>(estudantes, this->n_estudantes, CompareByNome);
            break;
        case INSERTIONSORT:
            insertion_sort<Estudante>(estudantes, this->n_estudantes, CompareByNome);
            break;
        case SELECTIONSORT:
            selection_sort<Estudante>(estudantes, this->n_estudantes, CompareByNome);
            break;
        case MERGESORT:
            merge_sort<Estudante>(estudantes, this->n_estudantes, CompareByNome);
            break;
        case QUICKSORT:
            quick_sort<Estudante>(estudantes, this->n_estudantes, CompareByNome);
            break;
        default:
           cout << "metodo de ordenação não encontrado" << endl;
            break;
    }
    ExportaVectorToCSV(estudantes, this->n_estudantes, arquivo_txt);


    delete [] estudantes;

}

void Siga::SalvarCSVOrdenadoPorCursoENome(std::string arquivo_txt)
{

    Estudante* estudantes = ExtractAllEstudantsFromDataSet();
    // TODO: observe que quicksort foi usado primeiro pois ele não é estável
    insertion_sort<Estudante>(estudantes, this->n_estudantes, CompareByNome);
    insertion_sort<Estudante>(estudantes, this->n_estudantes, CompareByCurso);

    ExportaVectorToCSV(estudantes, this->n_estudantes, arquivo_txt);
    delete [] estudantes;

}

Estudante* Siga::ExtractAllEstudantsFromDataSet()
{
    cout <<" Extraindo " << this->n_estudantes << " registros de estudantes" << endl;
    //TODO: Extrair todos em estudantes cadastrados do arquivo.
    // 1. Alocar um vetor de estudante.
    Estudante* est = new Estudante[this->n_estudantes];
    // 2. Ler os registros de estudantes do arquivo binário, armazenando-os no vetor.
    this->file_stream.seekg(0, this->file_stream.beg);
    for( int i = 0; i < this->n_estudantes; i++){
        est[i].LerDoArquivoBinario(this->file_stream);
    }
    // 3. Retornar o vetor.
    return est;
}

void Siga::ExportaVectorToCSV(Estudante *E, int n, string file)
{
    ofstream arq;
    arq.open(file);
    if(arq.is_open())
    {
        //arq << "matricula;nome;ano de ingresso;curso;IRA" << endl;
        //Le o vetor de Estudantes
        for( int i = 0; i < this->n_estudantes; i++)
        {        
            arq << E[i].ObterMatricula() << ";";
            arq << E[i].ObterNome() << ";";
            arq << E[i].ObterAnoIngresso() << ";";
            arq << E[i].ObterCurso() << ";";
            arq << E[i].ObterIRA() << endl   ; 
        }
        arq.close();
    }
    else
    {
        cout << "FALHA AO ABRIR O ARQUIVO" << endl;
    }
}


