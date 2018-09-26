#ifndef _DATAANALYSIS_H_
#define _DATAANALYSIS_H_

#include "constList.h"
#include "basicClass.hpp"
#include "utilityTool.hpp"
#include <Eigen/Dense>
#include <fstream>
#include <iostream>
#include <map>
#include <set>
#include <vector>

using namespace std;
using namespace Eigen;

class dataAnalysis
{
private:
	map<Face,vector<VoxelIndex>> f_temp_map;
	const UtilityTool* tools;

private:	
	void getFaceTempMap(const vector<Voxel>& v_set);

public:
	dataAnalysis(){}
	~dataAnalysis(){}
	void setData(const UtilityTool *_tools, const vector<Voxel>& v_set);

	void getFaceSet(const vector<Voxel>& v_set, vector<Face>& f_set);
	void getFaceNearNode(const vector<Voxel>& v_set,const vector<Face>& f_set,vector<set<NodeIndex>>& face_near_n);
	void getFaceNearVoxel(const vector<Voxel>& v_set,const vector<Face>& f_set,vector<set<VoxelIndex>>& face_near_v);
	void getOneRingVoxel(const int& n_num,const vector<Voxel>& v_set,vector<set<VoxelIndex>>& ring_1_v);
	void getOneRingNode(const int& n_num,const vector<Voxel>& v_set,vector<set<NodeIndex>>& ring_1_n);
	void getOneRingFace(const int& n_num,const vector<set<NodeIndex>>& face_near_n,vector<set<FaceIndex>>& ring_1_f);
	void getFaceType(const Plain& inlet,const Plain& outlet,const vector<Face>& f_set,const vector<set<NodeIndex>>& face_near_n,vector<ElementType>& f_type);
	void getNodeType(const int& n_num,const vector<Face>& f_set,const vector<ElementType>& f_type,vector<ElementType>& n_type);
	void getVolumeSet(const int& v_num, const vector<Voxel>& v_set, vector<double>& volume_set);
	void getGradientSet(const vector<Voxel>& v_set, const vector<double>& volume_set, const Matrix<double,Dynamic,3>& n_set, vector<vector<Vector3d>>& gradient_set);
};

void dataAnalysis::getGradientSet(const vector<Voxel>& v_set, const vector<double>& volume_set, const Matrix<double,Dynamic,3>& n_set, vector<vector<Vector3d>>& gradient_set)
{
	for (int i=0 ; i<v_set.size();i++)
	{
		vector<Vector3d> temp; 
		gradient_set.push_back(temp);
	}

	for (int i=0;i<v_set.size();i++) 
	{
		for (int j=0;j<VOXEL_NODE_NUM;j++)   
		{
			int a=(j+1)%VOXEL_NODE_NUM;
			int b=(j+2)%VOXEL_NODE_NUM;
			int c=(j+3)%VOXEL_NODE_NUM;
			Vector3d e1=n_set.row(v_set[i][a])-n_set.row(v_set[i][b]);
			Vector3d e2=n_set.row(v_set[i][a])-n_set.row(v_set[i][c]);
			Vector3d cross_res=e1.cross(e2);
			double area = 0.5*cross_res.norm();
			if (area < 0) area = -area;
			cross_res.normalize();

			Vector3d e=n_set.row(v_set[i][j])-n_set.row(v_set[i][a]); 
			if (e.dot(cross_res)<0) cross_res=-cross_res; 

			Vector3d gradient=cross_res*area/(3.0*volume_set[i]); 
			gradient_set[i].push_back(gradient);  
		}
	}
}

void dataAnalysis::getVolumeSet(const int& v_num, const vector<Voxel>& v_set, vector<double>& volume_set)
{
	volume_set.reserve(v_num);

	for (int i=0;i<v_num;i++)
	{
		double volume=tools->volume(v_set[i]);
		volume_set.push_back(volume);
	}
}

void dataAnalysis::setData(const UtilityTool *_tools,const vector<Voxel>& v_set)
{
	tools=_tools;
	getFaceTempMap(v_set);
}

void dataAnalysis::getFaceTempMap(const vector<Voxel>& v_set)
{
	for (int i=0;i<v_set.size();i++)
	{
		for (int j=0;j<4;j++)
		{
			Face f;
			for (int k=0;k<3;k++)
				f[k]=v_set[i][(j+1+k)%4];
			f.sort();
			if (f_temp_map.find(f)==f_temp_map.end())
			{
				vector<VoxelIndex> temp;
				temp.push_back(i);
				f_temp_map[f]=temp;
			}
			else f_temp_map[f].push_back(i);
		} //4 faces for each voxel
	}
}

void dataAnalysis::getFaceSet(const vector<Voxel>& v_set, vector<Face>& f_set)
{
	map<Face,vector<VoxelIndex>>::iterator it;
	for (it=f_temp_map.begin();it!=f_temp_map.end();it++)
	{
		f_set.push_back(it->first);
	}
}

void dataAnalysis::getFaceNearNode(const vector<Voxel>& v_set,const vector<Face>& f_set,vector<set<NodeIndex>>& face_near_n)
{
	face_near_n.reserve(f_set.size());
	for (int i=0;i<f_set.size();i++)
	{
		set<NodeIndex> temp;
		face_near_n.push_back(temp);
	}//init

	for (int i=0;i<f_set.size();i++)
	{
		for (int j=0;j<f_temp_map[f_set[i]].size();j++)
		{
			VoxelIndex v_index=f_temp_map[f_set[i]][j];
			Voxel v=v_set[v_index];
			for (int k=0;k<4;k++)
			{
				if (!f_set[i].contains(v[k]))
				{
					face_near_n[i].insert(v[k]);
					break;
				}
			}
		}

	}
}

void dataAnalysis::getFaceNearVoxel(const vector<Voxel>& v_set,const vector<Face>& f_set,vector<set<VoxelIndex>>& face_near_v)
{
	face_near_v.reserve(f_set.size());
	for (int i=0;i<f_set.size();i++)
	{
		set<VoxelIndex> temp;
		face_near_v.push_back(temp);
	}

	for (int i=0;i<f_set.size();i++)
	{
		for (int j=0;j<f_temp_map[f_set[i]].size();j++)
		{
			face_near_v[i].insert(f_temp_map[f_set[i]][j]);
		}
	}
}

void dataAnalysis::getOneRingVoxel(const int& n_num,const vector<Voxel>& v_set,vector<set<VoxelIndex>>& ring_1_v)
{
	ring_1_v.reserve(n_num);
	for (int i=0;i<n_num;i++)
	{
		set<VoxelIndex> temp;
		ring_1_v.push_back(temp);
	}

	for (int i=0;i<v_set.size();i++)
	{
		for (int j=0;j<4;j++)
			ring_1_v[v_set[i][j]].insert(i);
	}
}

void dataAnalysis::getOneRingNode(const int& n_num,const vector<Voxel>& v_set,vector<set<NodeIndex>>& ring_1_n)
{
	ring_1_n.reserve(n_num);
	for (int i=0;i<n_num;i++)
	{
		set<NodeIndex> temp;
		ring_1_n.push_back(temp);
	}

	for (int i=0;i<v_set.size();i++)
	{
		for (int j=0;j<4;j++)
		{
			for (int k=0;k<3;k++)
				ring_1_n[v_set[i][j]].insert(v_set[i][(j+1+k)%4]);
		}
	}
}

void dataAnalysis::getOneRingFace(const int& n_num,const vector<set<NodeIndex>>& face_near_n,vector<set<FaceIndex>>& ring_1_f)
{
	ring_1_f.reserve(n_num);
	for (int i=0;i<n_num;i++)
	{
		set<FaceIndex> temp;
		ring_1_f.push_back(temp);
	}

	for (int i=0;i<face_near_n.size();i++)
	{
		set<NodeIndex>::iterator it;
		for (it=face_near_n[i].begin();it!=face_near_n[i].end();it++)
		{
			ring_1_f[*it].insert(i);
		}
	}
}

void dataAnalysis::getFaceType(const Plain& inlet,const Plain& outlet,const vector<Face>& f_set,const vector<set<NodeIndex>>& face_near_n,vector<ElementType>& f_type)
{
	f_type.reserve(face_near_n.size());

	for (int i=0;i<face_near_n.size();i++)
	{
		ElementType temp;
		f_type.push_back(temp);
	}

	/*------------------------------------------------------------ */
	/*
	Plain outlet1, outlet2, outlet3, outlet4, outlet5, outlet6;
	outlet1[0] = 0; outlet1[1] = 1; outlet1[2] = 0; outlet1[3] = -3; //y=3;
	outlet2[0] = 0; outlet2[1] = 1; outlet2[2] = 0; outlet2[3] = -10; //y=10;
	outlet3[0] = 1; outlet3[1] = 0; outlet3[2] = 0; outlet3[3] = -3; //x=3;
	outlet4[0] = 1; outlet4[1] = 0; outlet4[2] = 0; outlet4[3] = 3;// x=-3;
	outlet5[0] = 0; outlet5[1] = 0; outlet5[2] = 1; outlet5[3] = -3;//z=3;
	outlet6[0] = 0; outlet6[1] = 0; outlet6[2] = 1; outlet6[3] = 3;//z=-3;*/
	/*--------------------------------------------------------------*/

	int hard_num=0, inlet_num=0, outlet_num=0, inner_num=0;

	for (int i=0;i<face_near_n.size();i++)
	{
		switch(face_near_n[i].size())
		{
			case 1:
			{
				f_type[i]=HARD;
				hard_num++;

				bool isInlet=true;
				for (int j=0;j<FACE_NODE_NUM;j++)
				{
					if (!tools->in_plain(f_set[i][j],inlet))
					{
						isInlet=false;
						break;
					}
				}
				if (isInlet)
				{
					f_type[i]=INLET;
					hard_num--;
					inlet_num++;
					break;
				}

				bool isOutlet=true;
				for (int j=0;j<FACE_NODE_NUM;j++)
				{
					if (!tools->in_plain(f_set[i][j],outlet))
					/*if (!tools->in_plain(f_set[i][j], outlet1) && !tools->in_plain(f_set[i][j], outlet2) &&
						!tools->in_plain(f_set[i][j], outlet3) && !tools->in_plain(f_set[i][j], outlet4) &&
						!tools->in_plain(f_set[i][j], outlet5) && !tools->in_plain(f_set[i][j], outlet6) )*/
					{
						isOutlet=false;
						break;
					}
				}
				if (isOutlet)
				{
					f_type[i]=OUTLET;
					hard_num--;
					outlet_num++;
					break;
				}
				break;
			}
			case 2:
			{
				f_type[i]=INNER;
				inner_num++;
				break;
			}
		}
	}
}

void dataAnalysis::getNodeType(const int& n_num,const vector<Face>& f_set,const vector<ElementType>& f_type,vector<ElementType>& n_type)
{
	n_type.reserve(n_num);
	for (int i=0;i<n_num;i++)
	{
		ElementType temp=INNER;
		n_type.push_back(temp);
	}

	for (int i=0;i<f_type.size();i++) 
	{
		for (int j=0;j<FACE_NODE_NUM;j++)
		{
			if (f_type[i]>n_type[f_set[i][j]])
			{
				n_type[f_set[i][j]]=f_type[i];
			}
		}
	}
}

#endif