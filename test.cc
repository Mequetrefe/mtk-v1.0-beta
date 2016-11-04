#include <iostream>

int main(void)
{
int ii;
int jj;
int kk=8;
int dim_null=3;
for (jj=0;jj<kk+1;jj++)
{
  for(ii=0; ii<dim_null; ii++)
  {
    std::cout<< (ii+kk-dim_null)+jj*kk << " " ;
    std::cout << (dim_null-ii-1)+jj*dim_null<<" ";
  }
  std::cout << std::endl;
}
return 0;
}
