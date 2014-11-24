namespace std
{
}
//test mixed element access member functions
#include <vector>

int main()
{
int i;
std::vector< int  , class std::allocator< int  >  > v1(100);
for (i = 1; i <= v1 . size() - 1; i += 1) {
v1 . at(i) = v1[i] + 1;
}
return 0;
}
