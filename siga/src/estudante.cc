#include <cstring>
#include <iostream>
#include <fstream>

#include "estudante.h"

using namespace std;


Estudante::Estudante(): _nome(""), _matricula(0), _ano_ingresso(0), _codigo_curso(0), _ira(0.0)
{};

Estudante::Estudante(const char* nome, int matricula, int ano_ingresso, int curso, float ira)
{    
     strcpy(_nome,nome);
    _matricula    = matricula;
    _ano_ingresso = ano_ingresso;
    _codigo_curso = curso;
    _ira          = ira;
}

Estudante::Estudante(const Estudante& e)
{
    strcpy(_nome,e._nome);
    _matricula    = e._matricula;
    _ano_ingresso = e._ano_ingresso;
    _codigo_curso = e._codigo_curso;
    _ira          = e._ira;
}


int Estudante::ObterMatricula()
{
    return this->_matricula;
}
        
const char* Estudante::ObterNome()
{
    return this->_nome;

}     
float Estudante::ObterIRA()
{
    return this->_ira;
}

int Estudante::ObterAnoIngresso()
{
    return this->_ano_ingresso;
}

int Estudante::ObterCurso()
{
    return this->_codigo_curso;
}

void Estudante::TrocarMatricula(int matricula)
{
    this->_matricula = matricula;
}

void Estudante::TrocarNome(const char* nome)
{
    strcpy(_nome,nome);
}
        
void Estudante::TrocarIRA(float ira)
{
    this->_ira = ira;
}

void Estudante::TrocarAnoIngresso(int ano_ingresso)
{
    this->_ano_ingresso = ano_ingresso;
}

void Estudante::TrocarCurso(int curso)
{
    this->_codigo_curso = curso;
}


void Estudante::Imprimir()
{
    cout <<"Nome: "            << this->_nome << endl;
    cout <<"Matricula: "       << this->_matricula << endl;
    cout <<"Ano de Ingresso: " << this->_ano_ingresso << endl;
    cout <<"Curso: "           << this->_codigo_curso << endl;
    cout <<"IRA: "             << this->_ira << endl;
}

void Estudante::LerDoArquivoBinario(fstream &file)
{
    //CHECK: veja como é feita a leitura do arquivo binário!
    file.read((char*)&_nome, sizeof(_nome));
    file.read((char*)&_matricula, sizeof(_matricula));
    file.read((char*)&_ano_ingresso, sizeof(_ano_ingresso));
    file.read((char*)&_codigo_curso, sizeof(_codigo_curso));
    file.read((char*)&_ira, sizeof(_ira));
}

void Estudante::EscreverNoArquivoBinario(fstream &file)
{
    //CHECK: veja como é feita a escrita do arquivo binário!
    file.write((char*)&_nome, sizeof(_nome));
    file.write((char*)&_matricula, sizeof(_matricula));
    file.write((char*)&_ano_ingresso, sizeof(_ano_ingresso));
    file.write((char*)&_codigo_curso, sizeof(_codigo_curso));
    file.write((char*)&_ira, sizeof(_ira));
}


bool CompareByMatricula(Estudante& a, Estudante& b)
{
    return a.ObterMatricula() < b.ObterMatricula();
}

bool CompareByNome(Estudante& a, Estudante& b)
{
    return (strcmp(a.ObterNome(),b.ObterNome())<0);
}

bool CompareByIRA(Estudante& a, Estudante& b)
{
    return a.ObterIRA() < b.ObterIRA();
}

bool CompareByAnoIngresso(Estudante& a, Estudante& b)
{
    return a.ObterAnoIngresso() < b.ObterAnoIngresso();
}

bool CompareByCurso(Estudante& a, Estudante& b)
{
    return a.ObterCurso() < b.ObterCurso();
}