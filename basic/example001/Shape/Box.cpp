#include <iostream>
#include "Box.h"
namespace Example::Shape{

Box::Box(double side_):side(side_){};

double Box::GetArea(){
    return side*side;
}
 
}