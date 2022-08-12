#pragma once
#include <iostream>

#include <utility>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <igl/heat_geodesics.h>


typedef double real_type;

typedef Eigen::Matrix<real_type, -1, -1> MatrixXR;
typedef Eigen::MatrixXi MatrixXI;
typedef Eigen::Matrix<real_type, -1, 1> VectorXR;

typedef 
 struct Mesh {
	MatrixXR V;
	MatrixXI F;
	Mesh() {}
	Mesh(MatrixXR& verts, MatrixXI& faces)
		: V(std::move(verts)), F(std::move(faces)) {}
}Mesh;


class SplocsSolver {
private:
	
public :
	std::unique_ptr<Mesh> meanshape_mesh_;
	std::vector<Mesh> meshes_;

	real_type d_min_, d_max_;

	bool is_compiled_;
	int component_num_;// component K
	int num_opt_;

	real_type scale_factor_;
	SplocsSolver();

	std::unique_ptr<MatrixXR> matrix_X_;
	std::unique_ptr<MatrixXR> matrix_W_;
	std::unique_ptr<MatrixXR> matrix_C_;

	//local support matrix
	std::unique_ptr<MatrixXR> matrix_local_sup_;

private:
	//initialize
	void make_X_matrix_by_meshes_();
	void setup_W_C_lsup_(); // setup W C local_support.
	


	//
	void compile();
	


	struct igl::HeatGeodesicsData<real_type> data_;
	void get_local_support(const int idx, MatrixXR& result);
	void precompute_local_support(const Mesh& t);

public:

	void solve(int component_num, real_type d_min, real_type d_max);
	void solve(int component_num);

	void phase1(MatrixXR& C, MatrixXR& W);


	
	void add_triangle_mesh(std::string& name );
	void add_triangle_meshes(std::vector<std::string>& names);

	Mesh* get_mesh_component();
	std::vector<Mesh*> get_mesh_components();
	std::vector<Mesh*> get_weights_of_componets();

};