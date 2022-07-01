#include <iostream>
#include <chrono>

#include "siga.h"
#include "leftist_heap.h"
#include "sort.h" 

int main(int argc, char* argv[])
{
 
    string input = INPUT_DIR + string("alunos.csv");
    Siga siga("dataset.bin");
    siga.ImportCSVData(input);

    int n_estudantes = siga.ObterNumeroEstudantesArmazenados();

    int  mid = n_estudantes / 2;
    Estudante* Estudantes = siga.ExtractAllEstudantsFromDataSet();

    LeftistHeap<Estudante> heapA(CompareByIRA);
    LeftistHeap<Estudante> heapB(CompareByIRA);

    for(int i = 0; i < mid; i++)
    {
        heapA.Insert(Estudantes[i]);
    }

    for(int i=mid; i < n_estudantes; i++)
    {
        heapB.Insert(Estudantes[i]);
    }

    Estudante MaiorA = heapA.FindMax();
    Estudante MaiorB = heapB.FindMax();
    MaiorA.Imprimir();
    MaiorB.Imprimir();

    if(MaiorA.ObterMatricula() == MaiorB.ObterMatricula())
    {
        cout << "Teste 6 falhou" << endl;
        return -1;
    }
    
    MaiorA = heapA.RemoveMax();
    MaiorA = heapA.FindMax();
    MaiorA.Imprimir();
    if(MaiorA.ObterMatricula() != 200490)
    {
        cout << "Teste 6 falhou" << endl;
        return -1;
    }
    MaiorB = heapB.RemoveMax();
    MaiorB = heapB.FindMax();
    MaiorB.Imprimir();
    if(MaiorB.ObterMatricula() != 200154)
    {
        cout << "Teste 6 falhou" << endl;
        return -1;
    }


    heapA.Merge(heapB);
    MaiorA = heapA.FindMax();
    MaiorA.Imprimir();
    if(MaiorA.ObterMatricula() != 200154)
    {
        cout << "Teste 6 falhou" << endl;
        return -1;
    }


    delete [] Estudantes;

    return 0;
}

  
    
    
