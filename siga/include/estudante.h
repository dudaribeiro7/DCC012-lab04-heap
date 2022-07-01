#ifndef ESTUDANTE_H
#define ESTUDANTE_H

#include <fstream>
#include <cstring>

using namespace std;

#define STR_MAX_SIZE 64
#define ESTUDANTE_SIZE (sizeof(char)*STR_MAX_SIZE + 3*sizeof(int) + sizeof(float))

class Estudante
{
    public:
        Estudante();
        Estudante(const char* nome, int matricula, int ano_ingresso, int curso, float ira);
        Estudante(const Estudante& e);

        int         ObterMatricula();
        int         ObterAnoIngresso();
        int         ObterCurso();
        const char* ObterNome();
        float       ObterIRA();

        void TrocarMatricula(int matricula);
        void TrocarNome(const char* nome);
        void TrocarIRA(float ira);
        void TrocarAnoIngresso(int ano_ingresso);
        void TrocarCurso(int curso);

        void Imprimir();

        void LerDoArquivoBinario(fstream &file);
        void EscreverNoArquivoBinario(fstream &file);
        

        //friend ostream& operator<<(ostream& os, const Estudante& e);

    private:
        char   _nome[STR_MAX_SIZE];
        int    _matricula;
        int    _ano_ingresso;
        int    _codigo_curso;
        float  _ira;
};


bool CompareByMatricula(Estudante& a, Estudante& b);

bool CompareByNome(Estudante& a, Estudante& b);

bool CompareByIRA(Estudante& a, Estudante& b);

bool CompareByAnoIngresso(Estudante& a, Estudante& b);

bool CompareByCurso(Estudante& a, Estudante& b);


#endif /* ESTUDANTE_H */
