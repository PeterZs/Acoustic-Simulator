#include <iostream>
#include <map>
#include <sstream>

#include "simulator.hpp"
#include "constList.h"

using namespace std;

int main()
{
	string vtk_file,config_file;
	int op_type; 


	cin >> vtk_file >> config_file;
	cin >> op_type; //0 for muffler simulation, 1 for experiment. The difference is the method of reading config.

	Simulator simulator(vtk_file,config_file,op_type);
	
	simulator.simulate();

	return 0;
}

