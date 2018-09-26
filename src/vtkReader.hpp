#ifndef _VTKREADER_H_
#define _VTKREADER_H_

#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <string>
#include <vector>
#include <fstream>

#include "basicClass.hpp"

using namespace std;
using namespace Eigen;

class VtkReader
{
public:
	int readData(const string& fileName, Matrix<double,Dynamic,3>& n_set, vector<Voxel>& v_set, int& n_num,int& v_num)
	{
		int current_state=READ_NODES;
		string temp;

		fstream fin(fileName.c_str(), ios::in);
		
		if (!fin.is_open())
		{
			cout << "open vtk failed." << endl;
			return ERR_FILE_OPEN;
		}
		bool isEnd=false;
		for (int i = 0; i < 9; i++)
		{
			fin >> temp;
		} 

		while (!isEnd && (fin >> temp))
		{
			switch (current_state)
			{
				case READ_NODES:
					if (temp=="POINTS")
					{
						fin >> n_num;
						n_set=Matrix<double, Dynamic, 3>::Zero(n_num,3);
						fin >> temp;
						for (int i=0;i<n_num;i++)
						{
							for (int j=0;j<3;j++)
							{
								double num;
								fin >> num;
								n_set(i,j)=num;
							}
						}
						current_state=READ_VOXELS;
					}
					else return ERR_FILE_FORMAT;
				break;
				case READ_VOXELS:
					if (temp=="CELLS")
					{
						int v_tot_num;
						fin >> v_num;
						fin >> v_tot_num;
						if (v_tot_num==5*v_num)
						{
							v_set.reserve(v_num);
							for (int i=0;i<v_num;i++)
							{
								fin >> temp;
								Voxel v;
								for (int j=0;j<4;j++) fin >> v[j];
								v.sort();
								v_set.push_back(v);
							}
							isEnd=true;
						}
						else return ERR_FILE_FORMAT;
					}
					else return ERR_FILE_FORMAT;
				break;
			}
		}
		fin.close();
		return 0;
	}
};

#endif 