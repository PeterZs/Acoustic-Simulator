#ifndef _UTILITYTOOL_H_
#define _UTILITYTOOL_H_

#include "basicClass.hpp"
#include <Eigen/Dense>

using namespace Eigen;

class UtilityTool
{
private:
	const Matrix<double,Dynamic,3> *n_set;
public:
	UtilityTool();
	void setTool(const Matrix<double,Dynamic,3>* _n_set);
	~UtilityTool(){}

	Vector3d norm_vec(const Face& f) const; // return the normal vector of the face 
	Vector3d norm_vec(const Face& f, const NodeIndex& n) const;
	double area(const Face& f) const;  //return the area of the face
	double volume(const Voxel& v) const;  //return the volume of the voxel
	double cos(const Face& f, const NodeIndex& n)const; //return cos(a), where a is the angle at point n in the face
	double cot(const Face& f, const NodeIndex& n)const;
	bool in_plain(const NodeIndex& n, const Plain& p)const;  
};

UtilityTool::UtilityTool()
{
	n_set=NULL;
}

void UtilityTool::setTool(const Matrix<double, Dynamic, 3>* _n_set)
{
	n_set=_n_set;
}

Vector3d UtilityTool::norm_vec(const Face& f) const
{
	Vector3d e1=n_set->row(f[0])-n_set->row(f[1]);
	Vector3d e2=n_set->row(f[0])-n_set->row(f[2]);
	Vector3d norm_res=e1.cross(e2);
	return norm_res;
}
Vector3d UtilityTool::norm_vec(const Face& f, const NodeIndex& n)const
{
	Vector3d norm_res=norm_vec(f);
	Vector3d e=n_set->row(f[1])-n_set->row(n);
	if (e.dot(norm_res)<0)
	{
		norm_res=-norm_res;
	}
	return norm_res;
}

double UtilityTool::area(const Face& f) const
{
	double area_res=0.5*norm_vec(f).norm();
	if (area_res<0) area_res=-area_res;
	return area_res;
}

double UtilityTool::volume(const Voxel& v) const
{
	Matrix3d vol;
	vol.col(0)=n_set->row(v[0])-n_set->row(v[1]);
	vol.col(1)=n_set->row(v[0])-n_set->row(v[2]);
	vol.col(2)=n_set->row(v[0])-n_set->row(v[3]);

	double volume_res=vol.determinant()/6.0;
	if (volume_res<0) volume_res=-volume_res;
	return volume_res;
}

double UtilityTool::cos(const Face& f,const NodeIndex& n) const
{
	double cos_res=2;
	for (int i=0;i<FACE_NODE_NUM;i++)
	{
		if (f[i]==n)
		{
			Vector3d e1=n_set->row(f[i+1]%FACE_NODE_NUM)-n_set->row(f[i]);
			Vector3d e2=n_set->row(f[i+2]%FACE_NODE_NUM)-n_set->row(f[i]);
			cos_res=e1.dot(e2)/sqrt(e1.dot(e1)*e2.dot(e2));
			break;
		}
	}
	return cos_res;
}
double UtilityTool::cot(const Face& f,const NodeIndex& n)const
{
	double cosine=cos(f,n);
	double sine=sqrt(1-cosine*cosine);
	double cot_res=1.0*cosine/sine;
	return cot_res;
}

bool UtilityTool::in_plain(const NodeIndex& n,const Plain& P)const
{
	double judgement=0;
	for (int i=0;i<3;i++)
	{
		judgement+=n_set->coeff(n,i)*P[i];
	}
	judgement+=P[3];

	if (judgement<TOLERANCE&&judgement>-TOLERANCE) return true;
	else return false;
} 

#endif