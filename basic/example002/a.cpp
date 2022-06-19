
#include "a.h"

template<class T>
auto A<T>::Add(auto a)->decltype(a+x){
    return a+x;
}

// list types of implementations
template
auto A<double>::Add(double a)->decltype(a+x);

template
auto A<int>::Add(int a)->decltype(a+x);

