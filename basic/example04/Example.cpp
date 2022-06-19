#include <iostream>
#include <vector>

using namespace std;

// Create metafunction
template <class T>
struct TypeDecider
{
    using inType = T;
    using outType = string;
};

template <>
struct TypeDecider<int>
{
    using inType = int;
    using outType = string;
};

template <>
struct TypeDecider<double>
{
    using inType = double;
    using outType = char;
};
// end of metafunction

// Define a function that does the same job
template<class T>
auto constexpr GetOutType(){
    if constexpr (is_same_v<T, int>)
        return string{};
    else if (is_same_v<T, double>)
        return char{};

}

// A class example created using meta function TypeDecider
template <class T>
class App
{
public:
    typename TypeDecider<T>::outType m;
    decltype(GetOutType<T>()) n;
};

// A function example using the metafunction, TypeDecider.
// We also, use if constexpr() to decide types within function body.
template <class T>
void f(vector<T> v)
{
    typename TypeDecider<T>::outType x;
    cout << x.size()<<"\n";

    // deciding the type and do stuff
    // Note this is a compile-time if,
    // so if a condition is not met, the
    // subsequent lines won't be in the 
    // compiled program.
    if constexpr (is_same_v<T, int>)
        cout << "type is int\n";
}

int main()
{
    App<double> app1; // has member => char m, char n
    App<int> app2; // has member => string m, string n

    vector<int> v1; // has local variable => string x
    f(v1);

    vector<double> v2;
    //f(v2); this will give error because char doesn't have function size().

}
