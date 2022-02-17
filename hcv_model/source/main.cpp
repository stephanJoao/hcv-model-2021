#include <iostream>
#include <unistd.h>

#include "HCV_Model.hpp"

using namespace std;

int main()
{
  	HCV_Model* model = new HCV_Model("../diferential_evolution/output/parametros_DE.txt");
  	
	char cwd[1024];
  	getcwd(cwd, sizeof(cwd));
  	
	model->solve();
	return 0;
}