#include <iostream>
#include "a.h"

int main(){

  A<double> a1{.x=10};
  A<int> a2{.x=10};
  A<int> a3{.x=10};

  auto b = a1.Add(1.0); // from header, compiler knows b is double

  std::cout<<a1.Add(20.0)<<std::endl; // 30
  std::cout<<a2.Add(20)<<std::endl; //30
 // std::cout<<a3.Add(20.0)<<std::endl; // Error: double A<int>::Add<double>(double) not found!

  return 0;
}