#include <fstream>
#include <iostream>
#include <string>

using namespace std;

int main()
{
	string buffer;
	string vol_file;
	string vtk_file;
	cin>>vol_file>>vtk_file;

	fstream fin(vol_file.c_str(),ios::in);
	if (!fin.is_open())
	{
		cout<<"open vol file failed."<<endl;
		return 0;
	}

	fstream fout(vtk_file.c_str(),ios::out);
	fout<<"# vtk DataFile Version 2.0" <<endl;
	fout<<"TET"<<endl;
	fout<<"ASCII"<<endl<<endl;
	fout<<"DATASET UNSTRUCTURED_GRID"<<endl;

	while(fin>>buffer)
	{
		if (buffer=="points")
		{
			int n_num;
			fin>>n_num;
			fout<<"POINTS "<<n_num<<" float"<<endl;
			for (int i=0;i<n_num;i++)
			{
				float x,y,z;
				fin>>x>>y>>z;
				fout<<x<<" "<<y<<" "<<z<<endl;
			}
		}
	}

	fin.clear();
	fin.seekg(0,ios::beg);
	while (fin>>buffer)
	{
		if (buffer=="volumeelements")
		{
			int v_num;
			fin>>v_num;
			fout<<"CELLS "<<v_num<<" "<<v_num*5<<endl;
			for (int i=0;i<v_num;i++)
			{
				int temp,num,v1,v2,v3,v4;
				fin>>temp>>num>>v1>>v2>>v3>>v4;
				fout<<num<<" "<<v1-1<<" "<<v2-1<<" "<<v3-1<<" "<<v4-1<<endl;
			}
			fout<<"CELL_TYPES "<<v_num<<endl;
			for (int i=0;i<v_num;i++)
			{
				fout<<"10"<<endl;
			}

		}
	}

	fin.close();
	fout.close();

	return 0;
}