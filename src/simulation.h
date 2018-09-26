#ifndef _SIMULATION_H_
#define _SIMULATION_H_

#include <Eigen/Sparse>

#include "constList.h"
#include "basicClass.hpp"
#include "configReader.hpp"
#include "vtkReader.hpp"
#include "utilityTool.hpp"
#include "dataAnalysis.hpp"

class Simulator
{
private:
	//Friend Class
	ConfigReader config_reader;
	VtkReader vtk_reader;
	UtilityTool tools;
	dataAnalysis data_analysis;
	
	//Config Attribute config_reader
	double freq_start,freq_end,freq_step;
	double speed,rho;
	Plain inlet,outlet;
	double p0_real,p0_img;
	
	//Original Data vtk_reader
	int n_num,v_num;
	Eigen::Matrix<double, Dynamic, 3> n_set;
	
	//Topology Structure data_analysis
	std::vector<Voxel> v_set;
	std::vector<Face> f_set;

	std::vector<std::set<VoxelIndex>> ring_1_v;
	std::vector<std::set<NodeIndex>>  ring_1_n;
	std::vector<std::set<FaceIndex>>  ring_1_f;

	std::vector<std::set<NodeIndex>> face_near_n;
	std::vector<std::set<VoxelIndex>> face_near_v;


	std::vector<double> volume_set;
	std::vector<std::vector<Eigen::Vector3d>> gradient_set; 
	std::vector<ElementType> f_type;
	std::vector<ElementType> n_type;

	int equation_num;
	//Final result
	std::vector<Complex> p;

private:
	void load_vtk(const string vtk_file);
	void setup_utilitytools();
	void setup_data();
	void load_config(const string config_file,const int op_type);
	void simulate(double frequency);
	void evaluateSize();

public:
	Simulator(const string vtk_file, const string config_file, const int op_type);
	void simulate();
};

#endif