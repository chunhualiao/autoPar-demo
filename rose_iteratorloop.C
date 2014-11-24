namespace std
{
}
/* Test case contributed by Jeff Keasler
 * 1/26/2009
 *
 * */
#include <vector>

void foo()
{
std::vector< int ,class std::allocator< int  > > ::iterator i_nom_3;
std::vector< int  , class std::allocator< int  >  > v(10,7);
v . push_back(0);
for (i_nom_3 = v . begin();  * i_nom_3 != 0; i_nom_3 ++ ) {
}
}
