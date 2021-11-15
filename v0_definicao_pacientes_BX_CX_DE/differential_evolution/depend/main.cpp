#include "HCV_Model.h"
#include <iostream>
#include <unistd.h>
using namespace std;

int main(){
  HCV_Model* model = new HCV_Model();
  char cwd[1024];
  getcwd(cwd, sizeof(cwd));
  // cout << cwd<< endl;
  model->solve();//solve
  return 0;
}