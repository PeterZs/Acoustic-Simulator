#ifndef _CONFIGREADER_H_
#define _CONFIGREADER_H_

#include "basicClass.hpp"
#include <string>
#include <vector>
#include <fstream>

using namespace std;

class ConfigReader
{
public:
	int readConfig(const string& fileName,double& freq_start, double& freq_end, double& freq_step,
				   double& speed, double& rho,double& p0_real,double& p0_img,Plain& inlet,Plain& outlet)
	{
		string command;
		char temp[256];
		double chamber_radius[10];
		double chamber_length[10];
		double max_mesh_size;

		fstream fin(fileName.c_str(), ios::in);
		if (!fin.is_open()) return ERR_FILE_OPEN;
		int chamber_num;

		while (1)
		{
			fin >> command;
			if (command=="chamber_num:") 
			{
				fin>>chamber_num;
			}
			else if (command=="chamber_radius:")
			{
				for (int i=0;i<chamber_num;i++)
				{
					fin>>chamber_radius[i];
				}
			}
			else if (command=="chamber_length:")
			{
				double length,total_chamber_length=0.0;
				for (int i=0;i<chamber_num;i++)
				{
					fin>>chamber_length[i];
					total_chamber_length+=chamber_length[i];
				}
				inlet[0]=0;inlet[1]=1;inlet[2]=0;inlet[3]=0; //plain y=0
				outlet[0]=0;outlet[1]=1;outlet[2]=0;outlet[3]=-total_chamber_length;
			}
			else if (command=="max_mesh_size:")
			{
				fin>>max_mesh_size;
			}
			else if (command=="freq_start:")
			{
				fin>>freq_start;
			}
			else if (command=="freq_end:")
			{
				fin>>freq_end;
			}
			else if (command=="freq_step:")
			{
				fin>>freq_step;
			}
			else if (command=="speed:")
			{
				fin>>speed;
			}
			else if (command=="rho:")
			{
				fin>>rho;
			}
			else if (command=="p0_real:")
			{
				fin>>p0_real;
			}
			else if (command=="p0_img:")
			{
				fin>>p0_img;
			}
			else if (command=="#")
			{
				fin.getline(temp,256);
			}
			else if (command=="end")
			{
				break;
			}
			else return ERR_FILE_FORMAT;
		}
		return 0;	
    }
	//the difference is the format of the config file.
    int readConfig_ex(const string& fileName,double& freq_start, double& freq_end, double& freq_step,
  				  double& speed, double& rho,double& p0_real,double& p0_img,Plain& inlet, Plain& outlet)
    {
    	string command;
		char temp[256];

		fstream fin(fileName.c_str(), ios::in);
		if (!fin.is_open()) return ERR_FILE_OPEN;

		inlet[0]=0;inlet[1]=1;inlet[2]=0;inlet[3]=0;
		outlet[0]=0;outlet[1]=1;outlet[2]=0;outlet[3]=100;
		while (1)
		{
			fin >> command;

			if (command=="freq_start:")
			{
				fin>>freq_start;
			}
			else if (command=="freq_end:")
			{
				fin>>freq_end;
			}
			else if (command=="freq_step:")
			{
				fin>>freq_step;
			}
			else if (command=="speed:")
			{
				fin>>speed;
			}
			else if (command=="rho:")
			{
				fin>>rho;
			}
			else if (command=="p0_real:")
			{
				fin>>p0_real;
			}
			else if (command=="p0_img:")
			{
				fin>>p0_img;
			}
			else if (command=="#")
			{
				fin.getline(temp,256);
			}
			else if (command=="end")
			{
				break;
			}
			else return ERR_FILE_FORMAT;
		}
		return 0;	
    }
};

#endif