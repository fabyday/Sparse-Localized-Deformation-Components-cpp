#include <iostream>
#include <igl/heat_geodesics.h>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <stdlib.h>
#include "splocs.h"

SplocsSolver::SplocsSolver()
	:matrix_W_(std::make_unique<MatrixXR>()), 
	matrix_X_(std::make_unique<MatrixXR>()),
	matrix_C_(std::make_unique<MatrixXR>()), 
	d_min_(0.0), d_max_(1.0)
{}

void SplocsSolver::make_X_matrix_by_meshes_(){
	//const int mesh_num = meshes_.size();
	//const int v_size = meshes_[0].V.rows();

	//matrix_X_->resize(mesh_num, 3 * v_size);
	//MatrixXR mean_shape;
	//mean_shape.resize(1, matrix_X_->cols());
	//mean_shape.setZero();
	//std::for_each(meshes_.begin(), meshes_.end(), [&mean_shape, this, mesh_num,i = 0](Mesh& m) mutable {
	//	auto x = Eigen::Map<MatrixXR>(m.V.data(), 1, m.V.size());
	//	mean_shape += x / mesh_num;
	//	matrix_X_->row(i++) = x;
	//}

	//);


	//// check page 4
	////normalize part
	//matrix_X_->rowwise() -= mean_shape;
	////across vertices
	//MatrixXR standard_deviation = matrix_X_->colwise().squaredNorm() / std::sqrt(matrix_X_->rows());
	////across frames.
	//standard_deviation = standard_deviation.rowwise().squaredNorm() / std::sqrt(matrix_X_->cols());
	//matrix_X_->colwise() /= standard_deviation;

}


void SplocsSolver::setup_W_C_lsup_()
{
	const int K = component_num_;
	matrix_W_->resize(matrix_X_->rows(), K);
	matrix_C_->resize(K, matrix_X_->cols());

	//matrix_local_sup_->resize(K, matrix_X_->rows());



}


void SplocsSolver::compile(){
	if (is_compiled_)
		return;

	is_compiled_ = true;
	
	make_X_matrix_by_meshes_();
	setup_W_C_lsup_();

	calc_local_support();


}

void SplocsSolver::calc_local_support()
{
	/*
	std::for_each(meshes_.begin(), meshes_.end(), [this](Mesh& m) {
		igl::HeatGeodesicsData<real_type> data;
		real_type t = std::pow(igl::avg_edge_length(m.V, m.F), 2);
		if (!igl::heat_geodesics_precompute(m.V, m.F, t, data)) {
			std::cerr << "Error: heat_geodesics_precompute failed." << std::endl;
			exit(EXIT_FAILURE);
		}
		
		MatrixXR D;
		MatrixXI vids = MatrixXI::LinSpaced(Eigen::Sequential, 0, m.V.rows());
		igl::heat_geodesics_solve(data, vids, D);
		*(this->matrix_local_sup_) = std::move(D);
	});*/


}




void SplocsSolver::solve(int component_num, real_type d_min, real_type d_max)
{
	if (component_num_ != component_num) {
		is_compiled_ = false;
		component_num_ = component_num;
	}
	d_min_ = d_min;
	d_max_ = d_max;


	compile();


}

void SplocsSolver::solve(int component_num)
{
	solve(component_num, d_min_, d_max_);
}

void SplocsSolver::add_triangle_mesh(std::string& name)
{
	static int s = 0;
	MatrixXR V;	
	MatrixXI F;
	igl::read_triangle_mesh(name, V, F);
	this->meshes_.emplace_back(V, F);
	std::string ss = "test";
	char sst[1];
	igl::writePLY(ss+std::to_string(s)+".ply", this->meshes_[0].V, this->meshes_[0].F);
}

void SplocsSolver::add_triangle_meshes(std::vector<std::string>& names)
{
	for (auto name : names) {
		add_triangle_mesh(name);
	}
}

Mesh* SplocsSolver::get_mesh_component()
{
	return nullptr;
}

std::vector<Mesh*> SplocsSolver::get_mesh_components()
{
	return std::vector<Mesh*>();
}

std::vector<Mesh*> SplocsSolver::get_weights_of_componets()
{
	return std::vector<Mesh*>();
}
