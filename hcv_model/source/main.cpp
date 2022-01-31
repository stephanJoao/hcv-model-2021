#include <iostream>
#include <unistd.h>

#include "HCV_Model.hpp"

using namespace std;

int main()
{
  	HCV_Model* model = new HCV_Model();
  	char cwd[1024];
  	getcwd(cwd, sizeof(cwd));
  	// cout << cwd<< endl;
	model->solve();
	return 0;
}