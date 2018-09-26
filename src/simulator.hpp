#ifndef _SIMULATOR_H_
#define _SIMULATOR_H_

#include "simulation.h"

using namespace std;
using namespace Eigen;

void Simulator::load_vtk(const string vtk_file)
{
	vtk_reader.readData(vtk_file, n_set, v_set, n_num, v_num);
}

void Simulator::setup_utilitytools()
{
	tools.setTool(&n_set);
}

void Simulator::load_config(const string config_file, const int op_type)
{
	if (op_type == 0)
		config_reader.readConfig(config_file, freq_start, freq_end, freq_step,
			speed, rho, p0_real, p0_img, inlet, outlet);
	else if (op_type == 1) //without outlet because the whole room is closed
	{
		config_reader.readConfig_ex(config_file, freq_start, freq_end, freq_step,
			speed, rho, p0_real, p0_img, inlet, outlet);

	}
}
void Simulator::setup_data()
{
	data_analysis.setData(&tools, v_set);
	data_analysis.getFaceSet(v_set, f_set);
	data_analysis.getFaceNearNode(v_set, f_set, face_near_n);
	data_analysis.getFaceNearVoxel(v_set, f_set, face_near_v);
	data_analysis.getOneRingVoxel(n_num, v_set, ring_1_v);
	data_analysis.getOneRingNode(n_num, v_set, ring_1_n);
	data_analysis.getOneRingFace(n_num, face_near_n, ring_1_f);
	data_analysis.getFaceType(inlet, outlet, f_set, face_near_n, f_type);//outlet
	data_analysis.getNodeType(n_num, f_set, f_type, n_type);
	data_analysis.getVolumeSet(v_num, v_set, volume_set);
	data_analysis.getGradientSet(v_set, volume_set, n_set, gradient_set);
	
}
Simulator::Simulator(const string vtk_file, const string config_file, const int op_type)
{
	load_vtk(vtk_file);
	load_config(config_file,op_type);
	setup_utilitytools();
	setup_data();
	evaluateSize();
}

void Simulator::evaluateSize()
{
	this->equation_num = 0;
		
	for (int i = 0; i < n_num; i++)
	{
		if (n_type[i]==INNER||n_type[i]==HARD) this->equation_num += 2; //Helmholtz

	}

	for (int i = 0; i < f_set.size(); i++)
	{
		if (f_type[i] == HARD) this->equation_num += 2;
		if (f_type[i] == OUTLET) this->equation_num += 2;	
		if (f_type[i] == INLET) this->equation_num += 2;
	}
}

void Simulator::simulate() //simulate and output
{
	int index = 0;
	for (double freq = freq_start; freq <= freq_end; freq += freq_step,index++)
	{
		simulate(freq);
		for (int twice = 0; twice < 2; twice++)
		{
			const string fileName = to_string(index) + "_" + to_string(twice) + ".vtk";
			fstream fout(fileName.c_str(), ios::out);
			fout << "# vtk DataFile Version 2.0" << endl;
			fout << "TET" << endl;
			fout << "ASCII" << endl << endl;
			fout << "DATASET UNSTRUCTURED_GRID" << endl;
			fout << "POINTS " << n_num << " float" << endl;
			for (int i = 0; i < n_num; i++)
			{
				for (int j = 0; j < 3; j++) fout << n_set(i, j) << " ";
				fout << endl;
			}
			fout << "CELLS " << v_set.size() << " " << v_set.size()*(VOXEL_NODE_NUM + 1) << endl;
			for (int i = 0; i < v_set.size(); i++)
			{
				fout << VOXEL_NODE_NUM << " ";
				for (int j = 0; j < VOXEL_NODE_NUM; j++)
					fout << v_set[i][j] << " ";
				fout << endl;
			}
			fout << "CELL_TYPES " << v_set.size() << endl;
			for (int i = 0; i < v_set.size(); i++)
			{
				fout << "10" << endl;
			}
			fout << "POINT_DATA " << n_num << endl;
			fout << "SCALARS point_scalar_data float 1" << endl;
			fout << "LOOKUP_TABLE data_color" << endl;

			VectorXd real, img;
			real.resize(n_num);
			img.resize(n_num);
			for (int i = 0; i < n_num; i++)
			{
				real(i) = p[i].real;
				img(i) = p[i].img;
			}

			double max_real = 1;
			double min_real = -1;
			double max_img = 1;
			double min_img = -1;


			for (int i = 0; i < n_num; i++)
			{
				if (twice == 0)
				{
					fout << (real(i) - min_real) / (max_real - min_real) << endl;
				}
				if (twice == 1) 
				{ 
					fout << (img(i) - min_img) / (max_img - min_img) << endl; 
				}
			}

			fout << "LOOKUP_TABLE data_color 101" << endl;
			for (int i = 0; i <= 100; i++)
			{
				fout << 1.0*i / 100 << " " << 1.0 - 1.0*i / 100 << " 0.0 1.0" << endl;
			}

			fout.close();
		}
	}
}

void Simulator::simulate(double frequency)
{
	frequency = 2.0 * PI * frequency;
	SparseMatrix<double> A(Dynamic,Dynamic); 
	A.resize(equation_num,2*n_num);
	VectorXd x(2*n_num);      
	VectorXd y(equation_num); 	
	for (int i=0;i<equation_num;i++) y(i)=0;
	
	int equation_count=0;
	for (int i=0;i<n_num;i++)
	{
		if (n_type[i] == INNER || n_type[i] == HARD)
		{
			double *equation1 = (double*)malloc(sizeof(double) * 2 * n_num);
			memset(equation1, 0, sizeof(double) * 2 * n_num);
			double *equation2 = (double*)malloc(sizeof(double) * 2 * n_num);
			memset(equation2, 0, sizeof(double) * 2 * n_num);
			double tot_volume = 0;
			set<VoxelIndex>::iterator it;

			for (it = ring_1_v[i].begin(); it != ring_1_v[i].end(); it++)
			{
				tot_volume += volume_set[*it];
			}
			for (it = ring_1_v[i].begin(); it != ring_1_v[i].end(); it++)
			{
				int current_index; 
				for (int j = 0; j < VOXEL_NODE_NUM; j++)
				{
					if (v_set[*it][j] == i)
					{
						current_index = j;
						break;
					}
				}
				int a = v_set[*it][(current_index + 1) % VOXEL_NODE_NUM];
				int b = v_set[*it][(current_index + 2) % VOXEL_NODE_NUM];
				int c = v_set[*it][(current_index + 3) % VOXEL_NODE_NUM];

				Vector3d e1 = n_set.row(a) - n_set.row(b);
				Vector3d e2 = n_set.row(a) - n_set.row(c);
				Vector3d face_norm = e1.cross(e2);  

				double area = face_norm.norm()*0.5;
				if (area < 0) area = -area;
				face_norm.normalize();
				Vector3d e = n_set.row(a) - n_set.row(i);
				if (e.dot(face_norm) < 0) face_norm = -face_norm;

				for (int k = 0; k < VOXEL_NODE_NUM; k++)
				{
					int index = v_set[*it][k];
					equation1[index] += gradient_set[*it][k].dot(face_norm)*area / tot_volume;
					equation2[index+n_num] += gradient_set[*it][k].dot(face_norm)*area / tot_volume;
				}
			}
			equation1[i] += ((frequency*frequency)*1.0 / (speed*speed));
			equation2[i +n_num] += ((frequency*frequency)*1.0 / (speed*speed));

			for (int j = 0; j < 2 * n_num; j++)
			{
				if (equation1[j] != 0)
					A.insert(equation_count, j) = equation1[j];
			}
			equation_count++;
			for (int j = 0; j < 2 * n_num; j++)
			{
				if (equation2[j] != 0)
					A.insert(equation_count, j) = equation2[j];
			}
			equation_count++;
		}
		else continue;

	}	
	//HARD
	for (int i=0;i<f_set.size();i++)
	{
		if (f_type[i]==HARD)
		{
			double *equation1 = (double*)malloc(sizeof(double) * 2 * n_num);
			memset(equation1, 0, sizeof(double) * 2 * n_num);
			double *equation2 = (double*)malloc(sizeof(double) * 2 * n_num);
			memset(equation2, 0, sizeof(double) * 2 * n_num);
				
			set<NodeIndex>::iterator it = face_near_n[i].begin();
			Vector3d e1=n_set.row(f_set[i][0])-n_set.row(f_set[i][1]);
			Vector3d e2=n_set.row(f_set[i][0])-n_set.row(f_set[i][2]);
			Vector3d cross_res=e1.cross(e2); 
			cross_res.normalize();
			Vector3d e = n_set.row(f_set[i][0]) - n_set.row(*it);
			if (e.dot(cross_res)<0) cross_res = -cross_res; 
			
			set<VoxelIndex>::iterator iv = face_near_v[i].begin();
			for (int k=0;k<VOXEL_NODE_NUM;k++)
			{
				int index=v_set[*iv][k];
				equation1[index] += (gradient_set[*iv][k].dot(cross_res));
				equation2[index+n_num] += (gradient_set[*iv][k].dot(cross_res));
			}
			for (int j=0;j<2*n_num;j++)
			{
				if (equation1[j]!=0)
				A.insert(equation_count,j) = equation1[j];
			}
			equation_count++;
			for (int j=0;j<2*n_num;j++)
			{
				if (equation2[j]!=0)
				A.insert(equation_count,j) = equation2[j];
			}
			equation_count++;
		}
	}
	//OUTLET 
	for (int i=0;i<f_set.size();i++)
	{
		if (f_type[i]==OUTLET)
		{
				double *equation1 = (double*)malloc(sizeof(double) * 2 * n_num);
				memset(equation1, 0, sizeof(double) * 2 * n_num);
				double *equation2 = (double*)malloc(sizeof(double) * 2 * n_num);
				memset(equation2, 0, sizeof(double) * 2 * n_num);

				set<NodeIndex>::iterator it = face_near_n[i].begin();
				Vector3d e1=n_set.row(f_set[i][0])-n_set.row(f_set[i][1]);
				Vector3d e2=n_set.row(f_set[i][0])-n_set.row(f_set[i][2]);
				Vector3d cross_res=e1.cross(e2);
				cross_res.normalize();	
				Vector3d e=n_set.row(f_set[i][0])-n_set.row(*it);
				if (e.dot(cross_res)<0) cross_res=-cross_res; 
				
				set<VoxelIndex>::iterator iv = face_near_v[i].begin();
				for (int k=0;k<VOXEL_NODE_NUM;k++)
				{
					int index = v_set[*iv][k];
					equation1[index] += (gradient_set[*iv][k].dot(cross_res));
					equation2[index+n_num] += (gradient_set[*iv][k].dot(cross_res));
					//equation1[index+n_num] -= (frequency / speed / 4.0);
					//equation2[index] += (frequency / speed / 4.0);
				}

				for (int k = 0; k < FACE_NODE_NUM; k++)
				{
					int index = f_set[i][k];
					equation1[index+n_num] -= (frequency / speed / 3.0);
					equation2[index] += (frequency / speed / 3.0);
				}
				for (int j=0;j<2*n_num;j++)
				{
					if (equation1[j]!=0)
					A.insert(equation_count,j) = equation1[j];
				}
				equation_count++;
				for (int j=0;j<2*n_num;j++)
				{
					if (equation2[j]!=0)
					A.insert(equation_count,j) = equation2[j];
				}
				equation_count++;
		}
	}
	//INLET
	for (int i=0;i<f_set.size();i++)
	{
		if (f_type[i]==INLET)
		{
				double *equation1 = (double*)malloc(sizeof(double) * 2 * n_num);
				memset(equation1, 0, sizeof(double) * 2 * n_num);
				double *equation2 = (double*)malloc(sizeof(double) * 2 * n_num);
				memset(equation2, 0, sizeof(double) * 2 * n_num);

				set<NodeIndex>::iterator it = face_near_n[i].begin();
				Vector3d e1 = n_set.row(f_set[i][0]) - n_set.row(f_set[i][1]);
				Vector3d e2 = n_set.row(f_set[i][0]) - n_set.row(f_set[i][2]);
				Vector3d cross_res = e1.cross(e2);
				cross_res.normalize();
				Vector3d e = n_set.row(f_set[i][0]) - n_set.row(*it);
				if (e.dot(cross_res) < 0) cross_res = -cross_res;

				set<VoxelIndex>::iterator iv = face_near_v[i].begin();
				for (int k = 0; k < VOXEL_NODE_NUM; k++)
				{
					int index = v_set[*iv][k];
					equation1[index] += (gradient_set[*iv][k].dot(cross_res));
					equation2[index + n_num] += (gradient_set[*iv][k].dot(cross_res));
					//equation1[index+n_num] -= (frequency / speed / 4.0);
					//equation2[index] += (frequency / speed / 4.0);
				}
				
				for (int k = 0; k < FACE_NODE_NUM; k++)
				{
					int index = f_set[i][k];
					equation1[index+n_num] -= (frequency / speed/3.0);
					equation2[index] += (frequency / speed / 3.0);
				}
				for (int j = 0; j<2 * n_num; j++)
				{
					if (equation1[j] != 0)
						A.insert(equation_count, j) = equation1[j];
				}
				y(equation_count) =-2.0* p0_img * frequency / speed;
				equation_count++;

				for (int j = 0; j<2 * n_num; j++)
				{
					if (equation2[j] != 0)
						A.insert(equation_count, j) = equation2[j];
				}
				y(equation_count) = p0_real * (2.0)*frequency / speed;
				equation_count++;
		}
	}
	A.makeCompressed();
	SparseMatrix <double> B = A.transpose()*A;
	SimplicialLDLT<SparseMatrix<double>> solver;
	solver.compute(B);
	assert(solver.info() == Success);
	MatrixXd b = A.transpose()*y;
	x = solver.solve(b);
	assert(solver.info() == Success);

	for (int i = 0; i < n_num; i++)
	{
		Complex temp;
		temp.real = x(i);
		temp.img = x(i+n_num);
		p.push_back(temp);
	}
	
	/*for (int i = 0; i < n_num; i++)
	{
		if (n_type[i] == OUTLET)
		{
			cout << p[i].real << " ";
		}
	}
	cout << endl;


	double average = 0;
	int sum = 0;
	for (int i = 0; i < n_num; i++)
	{
		if (n_type[i] == OUTLET)
		{
			sum++;
			average += p[i].real;
		}
	}
	average = average * 1.0 / sum;
	cout << "average: " << average << endl;
	for (int i = 0; i < n_num; i++)
	{
		if (n_type[i] == OUTLET)
		{
			cout << (p[i].real - average) / p[i].real << " ";
		}
	}*/

}


#endif