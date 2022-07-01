#ifndef SORT_H
#define SORT_H

#include <chrono>
#include <iostream>
#include <vector>
using namespace std;

#include "metrics.h"


// Use essa função para movimentar dados 
template <typename T>
void troca(T& a, T& b)
{
    T tmp = a;
    a     = b;
    b     = tmp;
}


template <typename T>
void bubble_sort_internal(T* array, int size, bool (*compare)(T&, T&), PerformanceMetrics *p)  
{
    // bubble sort
    for (int i = size-2; i >= 0; i--)
    {
        for (int j = 0; j <= i; j++)
        {
            p->n_comp++;   // incrementa o número de comparações

            if (compare(array[j+1], array[j]))
            {
                p->n_mov+=3; // troca realiza tres movimentacoes de dados
                troca<T>(array[j], array[j + 1]);
            }
        }
    }
}



template <typename T>
void insertion_sort_internal(T* array, int size,  bool (*compare)(T&, T&), PerformanceMetrics *p)
{
    // Implementação do Insert Sort
    // TODO: Insira as métricas de performance  OK
    for (int i = 1; i < size; i++)
    {
        T key = array[i];
        int j = i - 1;
        while (j >= 0 && compare(key, array[j]))
        {
            //incrementa o número de comparações:
            p->n_comp+=2; 

            array[j + 1] = array[j];
            j = j - 1;

            //2 movimentações de dados realizadas acima:
            p->n_mov+=2;
        }
        //incrementa o número de comparações para a última feita no comando while:
        p->n_comp+=2;

        array[j + 1] = key;

        //3 movimentações de dados realizadas nesse for:
        p->n_mov+=3;
    }
}


template <typename T>
void selection_sort_internal(T* array, int size, bool (*compare)(T&, T&), PerformanceMetrics *p)
{
    // TODO: Insira as métricas de performance  OK
    for (int i = 0; i < size; i++)
    {
        int min = i;
        for (int j = i + 1; j < size; j++)
        {
            //incrementa o número de comparações:
            p->n_comp++;  
            if (compare(array[j], array[min])){
                min = j;

                //1 movimentação de dados realizada acima:
                p->n_mov++;
            }
        }
        troca(array[i], array[min]);

        //"troca" realiza 3 movimentações de dados e 1 movimentação foi feita nesse for::
        p->n_mov+=4;
    }
}

template <typename T>
void merge(T* array, int start, int middle, int end, bool (*compare)(T&, T&), PerformanceMetrics *p){
    int i = start;      //1
    int j = middle;     //2
    int k = 0;          //3

    //3 movimentações de dados realizadas acima:
    p->n_mov+=3;

    T* aux = new T[end - start];
    while((i < middle) && (j < end)){
        //incrementa o numero de comparações feitas pelo while:
        p->n_comp+=2;

        //incrementa o numero de comparações feita pelo if:
        p->n_comp++;

        if(compare(array[i], array[j])){
            aux[k] = array[i];  //1
            i++;                //2

            //2 movimentações de dados realizadas acima:
            p->n_mov+=2;
        }
        else{
            aux[k] = array[j];  //1
            j++;                //2

            //2 movimentações de dados realizadas acima:
            p->n_mov+=2;
        }
        k++;

        //1 movimentação de dados realizada acima:
        p->n_mov++;
    }
    //incrementa o numero de comparações feitas pelo while qdo ele acaba:
    p->n_comp+=2;

    while(i < middle){
        //incrementa o numero de comparações feitas pelo while:
        p->n_comp++;

        aux[k] = array[i];  //1
        i++;                //2
        k++;                //3

        //3 movimentações de dados realizadas acima:
        p->n_mov+=3;
    }
    //incrementa o numero de comparações feitas pelo while qdo ele acaba:
    p->n_comp++;

    while(j < end){
        //incrementa o numero de comparações feitas pelo while:
        p->n_comp++;

        aux[k] = array[j];  //1
        j++;                //2
        k++;                //3

        //3 movimentações de dados realizadas acima:
        p->n_mov+=3;
    }
    //incrementa o numero de comparações feitas pelo while qdo ele acaba:
    p->n_comp++;

    for(i = start; i < end; i++){
        array[i] = aux[i-start];

        //1 movimentação de dados realizada acima:
        p->n_mov++;
    }
}

template <typename T>
void mergesort_internal(T* array, int start, int end, bool (*compare)(T&, T&), PerformanceMetrics *p)
{
    // TODO: Implementar merge sort aqui    OK

    int middle;
    if(start < (end-1)){
        middle = (start + end) / 2;
        //1 movimentação de dados realizada acima:
        p->n_mov++;
        mergesort_internal(array, start, middle, compare, p);
        mergesort_internal(array, middle, end, compare, p);
        merge(array, start, middle, end, compare, p);
        //as movimentações de "merge" já estão sendo contabilizadas na própria função.
    }
    return;
}

template <typename T>
int particionamento(T* array, int low, int high, bool (*compare)(T&, T&), PerformanceMetrics *p){
    T pivo = array[(high + low) / 2];       //1
    int i = low - 1;                        //2
    int j = high + 1;                       //3

    //3 movimentações de dados realizadas acima:
    p->n_mov+=3;

    do
    {
        do
        {
            //incrementa o numero de comparações feitas pelo while:
            p->n_comp++;

            i++;

            //1 movimentação de dados realizada acima:
            p->n_mov++;

        } while (compare(array[i], pivo));

        do
        {
            //incrementa o numero de comparações feitas pelo while:
            p->n_comp++;

            j--;

            //1 movimentação de dados realizada acima:
            p->n_mov++;

        } while (compare(pivo, array[j]));

        //incrementa o numero de comparações feitas pelo if:
        p->n_comp++;

        if(i >= j){
            return j;
        }

        troca(array[i], array[j]);

        //"troca" realiza 3 movimentações de dados:
        p->n_mov+=3;

    } while (true);
}

template <typename T>
void quicksort_internal(T* array, int low, int high, bool (*compare)(T&, T&), PerformanceMetrics *p)
{
    // TODO: Implementação do quicksort     OK

    int middle;
    //incrementa o numero de comparações feitas pelo if:
    p->n_comp++;
    if(low < high){
        middle = particionamento(array, low, high, compare, p);
        quicksort_internal(array, low, middle, compare, p);
        quicksort_internal(array, middle+1, high, compare, p);

        //1 movimentação de dados realizada acima:
        p->n_mov++;
        //as movimentações de "particionamento" já estão sendo contabilizadas na própria função.
    }
    return;
}


template <typename T>
void bubble_sort(T* array, int size, bool (*compare)(T&, T&))  
{
    PerformanceMetrics p;
    SetUpPerformanceMetrics(&p);
    auto t1 = Clock::now();
    bubble_sort_internal<T>(array, size, compare, &p);
    auto t2 = Clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    PerformanceMetricsCPUTime(&p, diff.count());
    cout << "Bubble Sort:" << endl; 
    PerformanceMetricsPrint(&p);
}

template <typename T>
void insertion_sort(T* array, int size, bool (*compare)(T&, T&))  
{
    PerformanceMetrics p;
    SetUpPerformanceMetrics(&p);
    auto t1 = Clock::now();
    insertion_sort_internal(array, size, compare, &p);
    auto t2 = Clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    PerformanceMetricsCPUTime(&p, diff.count());
    cout << "Insertion Sort:" << endl; 
    PerformanceMetricsPrint(&p);
}

template <typename T>
void selection_sort(T* array, int size, bool (*compare)(T&, T&))  
{
    PerformanceMetrics p;
    SetUpPerformanceMetrics(&p);
    auto t1 = Clock::now();
    selection_sort_internal(array, size, compare, &p);
    auto t2 = Clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    PerformanceMetricsCPUTime(&p, diff.count());
    cout << "Selection Sort:" << endl; 
    PerformanceMetricsPrint(&p);
}

template <typename T>
void merge_sort(T* array, int size, bool (*compare)(T&, T&))  
{
    PerformanceMetrics p;
    SetUpPerformanceMetrics(&p);
    auto t1 = Clock::now();
    mergesort_internal(array, 0, size-1, compare, &p);
    auto t2 = Clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    PerformanceMetricsCPUTime(&p, diff.count());
    cout << "Merge Sort:" << endl; 
    PerformanceMetricsPrint(&p);
}

template <typename T>
void quick_sort(T* array, int size, bool (*compare)(T&, T&))  
{
    PerformanceMetrics p;
    SetUpPerformanceMetrics(&p);
    auto t1 = Clock::now();
    quicksort_internal(array, 0, size-1, compare, &p);
    auto t2 = Clock::now();
    std::chrono::duration<double> diff = t2 - t1;
    PerformanceMetricsCPUTime(&p, diff.count());
    cout << "Quick Sort:" << endl; 
    PerformanceMetricsPrint(&p);
}

#endif /* SORT_H */
