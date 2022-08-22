#pragma once
#include <iostream>

#include <utility>
#include <vector>
#include <string>
#include <Eigen/Core>
#include <igl/heat_geodesics.h>



#define F_EPS 10e-6

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

	real_type d_min_, d_max_, sparsity_lambda_, rho_;


	bool is_compiled_;
	int component_num_;// component K
	int num_opt_ = 10;
	int num_admm_iterations_;
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
	void get_local_support(const int idx, MatrixXR& result, real_type min, real_type max);
	void precompute_local_support(const Mesh& t);
	void scale_meshes();
	void find_rbm_procrustes();
public:

	void solve(	int component_num = 50, int num_iter_max = 10, int num_admm_iterations = 10,
				real_type d_min = 0.1, real_type d_max = 0.7, 
				real_type sparsity_lambda = 2.0, real_type rho_ = 10.0);

	//void solve(int component_num);

	void phase1(MatrixXR& C, MatrixXR& W);


	
	void add_triangle_mesh(std::string& name );
	void add_triangle_meshes(std::vector<std::string>& names);

	Mesh* get_mesh_component();
	std::vector<Mesh*> get_mesh_components();
	std::vector<Mesh*> get_weights_of_componets();

};