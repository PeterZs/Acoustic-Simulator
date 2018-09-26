#ifndef _BASICCLASS_H_
#define _BASICCLASS_H_

#include "constList.h"

class Complex
{
public:
	double real;
	double img;

public:
	Complex(){}
    Complex(double r, double i):real(r),img(i){}
    ~Complex(){}

    Complex operator=(const Complex& op)
    {
    	real = op.real;
    	img  = op.img;
    	return *this;
    }

    Complex operator+(const Complex& op)
    {
    	real += op.real;
    	img  += op.img;
    	return *this;
    }

    Complex operator-(const Complex& op)
    {
    	real -= op.real;
    	img  -= op.img;
    	return *this; 
    }

    Complex operator*(const Complex& op)
    {
    	real = real * op.real - img * op.img;
    	img  = real * op.img + img * op.real;
    	return *this;
    }

    Complex operator/(const Complex& op)
    {
    	double denominator = op.real * op.real + op.img * op.img;
    	real = (real * op.real + img * op.img) * 1.0 / denominator;
    	img  = (img  * op.real - real* op.img) * 1.0 / denominator;
    	return *this;
    }
};

template <class T,int n>
class MultiNode
{
private:
	T node[n];
public:
	MultiNode(){}
	MultiNode(const MultiNode& op)
	{
		for (int i=0;i<n;i++) node[i]=op.node[i];
	}
    ~MultiNode(){}

    T& operator[](int i)
    {
    	return node[i];
    }

	const T& operator[](int i) const
	{
		return node[i];
	}

    void operator=(const MultiNode& op)
    {
    	for (int i=0;i<n;i++) node[i]=op.node[i];
    }
	
	bool operator==(const MultiNode& op)
	{
		for (int i=0;i<n;i++) if (node[i]!=op.node[i]) return false;
		return true;
	}

	bool operator<(const MultiNode& op) const
	{
		for (int i=0;i<n;i++)
		{
			if (node[i]<op.node[i]) return true;
			if (node[i]>op.node[i]) return false;
		}
		return false;
	}

	bool contains(const T& op) const 
	{
		for (int i=0;i<n;i++)
		{
			if (node[i]==op) return true;
		}
		return false;
	}

	void sort()
	{
		for (int i=0;i<n;i++)
			for (int j=i+1;j<n;j++)
			{
				if (node[i]>node[j])
				{
					T temp=node[i];
					node[i]=node[j];
					node[j]=temp;
				}
			}
	}
};

typedef int NodeIndex;
typedef int EdgeIndex;
typedef int FaceIndex;
typedef int VoxelIndex;

typedef MultiNode<NodeIndex,EDGE_NODE_NUM> Edge;
typedef MultiNode<NodeIndex,FACE_NODE_NUM> Face;
typedef MultiNode<NodeIndex,VOXEL_NODE_NUM> Voxel;
typedef MultiNode<double,4> Plain;

#endif