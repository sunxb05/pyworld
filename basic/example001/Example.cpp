
#include <iostream>
#include "Shape/Box.h"

using namespace std;
using namespace Example::Shape;


int main(){
    
    Box box(5);

    cout<< box.GetArea()<<endl;
    return 0;
}
