template<class T>
struct A
{
    T x;
    auto Add(auto a)->decltype(a+x);
};

