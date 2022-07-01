#include "siga.h"
#include "leftist_heap.h"
#include "sort.h"


// Nesse teste, extraimos os alunos com maior IRA usando uma 
// heap esquerdista.

int main(int argc, char *argv[])
{
    string input = INPUT_DIR + string("alunos.csv");
    Siga siga("dataset.bin");
    siga.ImportCSVData(input);

    int n_estudantes = siga.ObterNumeroEstudantesArmazenados();

    Estudante *Estudantes = siga.ExtractAllEstudantsFromDataSet();

    // Construtor receve como parametro uma função de comparacao
    // Se não for passada nenhuma função de comparacao, a heap
    // emprega a comparacao < entre os elementos
    LeftistHeap<Estudante> heap(CompareByIRA);

    auto t1 = Clock::now();
    for (int i = 0; i < n_estudantes; i++)
    {
        heap.Insert(Estudantes[i]);
    }
    Estudante MaiorIra = heap.RemoveMax();
    auto t2 = Clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    cout << "Tempo gasto para inserir os elementos na heap: " << diff.count() << " segundos" << endl;

    if (MaiorIra.ObterMatricula() != 200770)
    {
        cout << "Teste falhou" << endl;
        return -1;
    }

    MaiorIra.Imprimir();

    Estudante SegundoMaiorIra = heap.RemoveMax();

    if (SegundoMaiorIra.ObterMatricula() != 200789)
    {
        cout << "Teste falhou" << endl;
        return -1;
    }

    SegundoMaiorIra.Imprimir();

    delete[] Estudantes;

    return 0;
}