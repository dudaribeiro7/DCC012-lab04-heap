#ifndef LEFTIST_HEAP_H
#define LEFTIST_HEAP_H

#include <iostream>
#include <string>

#include "leftist_heap_node.h"



template<typename T>
class LeftistHeap {
    public:
        /**
         * @brief Constroi uma nova heap esquerdista vazia.
         *        
         */
        LeftistHeap();

        // Cria uma heap vazia e usa a função compare para comparar os elementos.
        LeftistHeap(bool (*compare)(T&, T&));

        void Insert(T);     // Insere chave na heap
        T&   FindMax();     // Seleciona o maior elemento da heap
        T    RemoveMax();   // Remove o maior elemento da heap
        void Merge(LeftistHeap<T>&); // Intercala duas heaps esquerdistas
        void Print();                // Imprime a heap
        bool isEmpty();              // Verifica se a heap está vazia
        bool isLeftistHeap();        // Verifica se a heap é uma heap esquerdista
        ~LeftistHeap();
    private:

        // Funções Auxiliares
        // troca os filhos do nó
        void SwapChildren(LeftistHeapNode<T>* node); 

        // libera memoria alocada 
        void DeallocateMemory(LeftistHeapNode<T>*);  

        // auxilio para a impressão.
        void print(LeftistHeapNode<T>* node, const std::string& prefix, bool isLeft, std::ostream& out);
        
        // calcula o caminho minimo de um nó para a raiz.
        int  ComputeDist(LeftistHeapNode<T>*);

        // realiza o merge dados duas arvores esquerdistas.
        LeftistHeapNode<T>* merge(LeftistHeapNode<T>* left, LeftistHeapNode<T>* right);
        // auxiliar no merge
        LeftistHeapNode<T>* merge_auxiliar(LeftistHeapNode<T>* left, LeftistHeapNode<T>* right);    
        
        // auxiliar na 
        bool checkLeftistHeap(LeftistHeapNode<T>* node);
        
        // Ponteiro para a função de comparação
        bool (*_compare)(T&, T&); 

        // Atributos
        // Ponteiro para a raiz da heap
        LeftistHeapNode<T>* root;
};

/*
 *  Funções auxiliares
 */

template<typename T>
bool compare_function(T& a, T& b)
{
    return a < b;
}


template<typename T>
void LeftistHeap<T>::SwapChildren(LeftistHeapNode<T>* node)
{
    LeftistHeapNode<T>* temp = node->getLeftChild();
    node->setLeftChild(node->getRightChild());
    node->setRightChild(temp);
}

template<typename T>
int LeftistHeap<T>::ComputeDist(LeftistHeapNode<T>* node)
{
    if (node == nullptr)
        return -1;
    else if(node->getLeftChild() != nullptr && node->getRightChild() != nullptr)
            return 1 + std::min(ComputeDist(node->getLeftChild()), ComputeDist(node->getRightChild()));
         else 
            return 0;
}

// template<typename T>
// void troca(LeftistHeapNode<T>* A, LeftistHeapNode<T>* B){
//     LeftistHeapNode<T> aux = LeftistHeapNode<T>(A->getData(), A->getLeftChild(), A->getRightChild());
//     *A = LeftistHeapNode<T>(B->getData(), B->getLeftChild(), B->getRightChild());
//     *B = LeftistHeapNode<T>(aux.getData(), aux.getLeftChild(), aux.getRightChild());
// }

template<typename T>
LeftistHeapNode<T>* LeftistHeap<T>::merge_auxiliar(LeftistHeapNode<T>* A, LeftistHeapNode<T>* B)
{
    if(A->getLeftChild() == NULL){
        A->setLeftChild(B);
    }
    else{
        A->setRightChild(merge(A->getRightChild(), B));
        if(A->getRightChild()->GetDist() > A->getLeftChild()->GetDist()){
            SwapChildren(A);
        }
        if(A->getRightChild() == NULL){
            A->SetDist(0);
        }
        else{
            A->SetDist(1 + A->getRightChild()->GetDist());
        }
    }
    
    return A;
}

template<typename T>
LeftistHeapNode<T>* LeftistHeap<T>::merge(LeftistHeapNode<T>* A, LeftistHeapNode<T>* B)
{
    //TODO: IMPLEMENTAR MERGE   OK
   
    if(A == NULL){
        return B;
    }
    if(B == NULL){
        return A;
    }
    if(this->_compare(A->getData(), B->getData())){
        return merge_auxiliar(B, A);
    }
    else{
        return merge_auxiliar(A, B);
    }
}

template <typename T>
bool LeftistHeap<T>::checkLeftistHeap(LeftistHeapNode<T>* node)
{
    if (node == nullptr)
        return true;

    int distL = ComputeDist(node->getLeftChild());
    int distR = ComputeDist(node->getRightChild());
    if( distL < distR)
        return false;
    else
        return checkLeftistHeap(node->getLeftChild()) && checkLeftistHeap(node->getRightChild());

}

template<typename T>
LeftistHeap<T>::LeftistHeap() : root(nullptr) {
    this->_compare = compare_function;
}

template<typename T>
LeftistHeap<T>::LeftistHeap(bool (*compare)(T&, T&)) : root(nullptr) {
    this->_compare = compare;
}

template<typename T>
void LeftistHeap<T>::Insert(T data) 
{
    //TODO: IMPLEMENTAR INSERÇÃO    OK

    //cria uma heap com único nó:
    LeftistHeapNode<T> nova_heap = LeftistHeapNode<T>(data);
    //mescla com a heap original:
    this->root = merge(this->root, &nova_heap);
}

template<typename T>
T LeftistHeap<T>::RemoveMax()
{
    //TODO: IMPLEMENTAR REMOÇÃO     OK

    LeftistHeapNode<T>* r = this->root;
    //remove a raíz da heap e mescla os filhos:
    this->root = merge(this->root->getLeftChild(), this->root->getRightChild());
    //retorna o elemento removido:
    return r->getData();
}

template<typename T>
T& LeftistHeap<T>::FindMax()
{
    //TODO: IMPLEMENTAR SELEÇÃO     OK

    //retorna a primeira posição da heap, a raíz:
    return this->root->getData();
}

template<typename T>
void LeftistHeap<T>::Merge(LeftistHeap<T>& other)
{
    if(this == &other)
        return;
    root = merge(root, other.root);
    other.root = nullptr;
}

template<typename T>
LeftistHeap<T>::~LeftistHeap()
{
    DeallocateMemory(root);
}

template<typename T>
void LeftistHeap<T>::DeallocateMemory(LeftistHeapNode<T>* node)
{
    if(node == nullptr)
        return;
    DeallocateMemory(node->getLeftChild());
    DeallocateMemory(node->getRightChild());
    delete node;
}

template<typename T>
void LeftistHeap<T>::print(LeftistHeapNode<T>* node, const std::string& prefix, bool isLeft, std::ostream& out)
{
    if(node == nullptr)
        return;
    
    out << prefix;
    out << (isLeft ? "├──" : "└──" );

    // print the value of the node
    out <<"["<<node->getData()<<","<<node->GetDist()<<"]" << std::endl;

    print(node->getLeftChild() , prefix + (isLeft ? "│   " : "    "), true, out);
    print(node->getRightChild(), prefix + (isLeft ? "│   " : "    "), false, out);
}

template<typename T>
void LeftistHeap<T>::Print()
{
    cout << root << "" << false << endl;
    // print(root, "", false,  std::cout);
    // std::cout << std::endl;
}

template<typename T>
bool LeftistHeap<T>::isEmpty()
{
    return root == nullptr;
}


template<typename T>
bool LeftistHeap<T>::isLeftistHeap()
{
    return checkLeftistHeap(root);
}

#endif
