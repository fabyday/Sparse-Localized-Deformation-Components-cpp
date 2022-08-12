#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <stdlib.h>
#include "splocs.h"

SplocsSolver::SplocsSolver()
	:matrix_W_(std::make_unique<MatrixXR>()), 
	matrix_X_(std::make_unique<MatrixXR>()),
	matrix_C_(std::make_unique<MatrixXR>()), 
	meanshape_mesh_(std::make_unique<Mesh>()), 
	d_min_(0.0), d_max_(1.0)
{}

void SplocsSolver::make_X_matrix_by_meshes_(){
	const int mesh_num = meshes_.size();
	const int v_size = meshes_[0].V.rows();
	matrix_X_->resize(mesh_num, 3 * v_size);
	VectorXR mean_shape;
	mean_shape.resize( matrix_X_->cols());
	mean_shape.setZero();
	std::for_each(meshes_.begin(), meshes_.end(), [&mean_shape, this, mesh_num, i = 0](Mesh& m) mutable {
		auto x = Eigen::Map<MatrixXR>(m.V.data(), m.V.size(), 1);
		mean_shape += x / mesh_num;
		matrix_X_->row(i++) = x.transpose();
	}

	);


	// check page 4
	//normalize part
	matrix_X_->rowwise() -= mean_shape.transpose();
	meanshape_mesh_->V = std::move(Eigen::Map<MatrixXR>(mean_shape.data(), v_size, 3));
	meanshape_mesh_->F = meshes_[0].F;


	Eigen::Map<VectorXR> x_vec(matrix_X_->data(), matrix_X_->size());

	//std
	scale_factor_ = sqrt((x_vec.array() - x_vec.mean()).square().sum() / (x_vec.size() - 1));
	*matrix_X_ /= scale_factor_;



	////across vertices
	//MatrixXR standard_deviation = matrix_X_->colwise().squaredNorm() / std::sqrt(matrix_X_->rows());
	////across frames.
	//real_type frame_standard_deviation = (standard_deviation.rowwise().squaredNorm() / std::sqrt(matrix_X_->cols()))(0,0);
	//matrix_X_->colwise() /= frame_standard_deviation;
	//*matrix_X_ /= frame_standard_deviation;


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



}

void SplocsSolver::get_local_support(const int idx, MatrixXR& result)
{

	Eigen::Matrix<real_type, -1, 1> D;
	Eigen::VectorXi vid(1, 1); vid(0) = idx;

	igl::heat_geodesics_solve(data_, vid, D);

	result = std::move(D.array().min(d_min_).max(d_max_));
	result = (result.array() - d_min_) / (d_max_- d_min_);

}

void SplocsSolver::precompute_local_support(const Mesh& m)
{
	real_type t = std::pow(igl::avg_edge_length(m.V, m.F), 2);
	if (!igl::heat_geodesics_precompute(m.V, m.F, t, data_)) {
		std::cerr << "Error: heat_geodesics_precompute failed." << std::endl;
		exit(EXIT_FAILURE);
	}
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

	phase1(*matrix_C_, *matrix_W_);


}

void SplocsSolver::solve(int component_num)
{
	solve(component_num, d_min_, d_max_);
}



		

MatrixXR proj(const MatrixXR& src) {
	return src/src.cwiseAbs().maxCoeff();
}



void save_component_to_mesh(const MatrixXR& component) {

}


void outerProduct(MatrixXR& result, const MatrixXR& t1, const MatrixXR& t2) {
	 
}



void SplocsSolver::phase1(MatrixXR& C, MatrixXR& W){
	
	const int v_size = meshes_[0].V.rows();
	const int frame_size = meshes_.size();
	//copy
	MatrixXR residual_X = *matrix_X_;

	precompute_local_support(meshes_[0]);
	C.resize(component_num_, v_size*3);
	W.resize(frame_size, component_num_);
	



	for (int k = 0; k < component_num_; k++) {

		using vec = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
		vec magnitude(v_size);
		magnitude.setZero();
		for (int f = 0; f < frame_size; f++) {
			Eigen::Map<MatrixXR, 0>
				tmp(residual_X.block(f, 0, 1, v_size).data(),
					v_size, 3);
			magnitude += tmp.rowwise().squaredNorm();
		}


		int max_v_idx;
		magnitude.maxCoeff(&max_v_idx);
		







		Eigen::JacobiSVD<MatrixXR> svd(matrix_X_->
			block(0, max_v_idx, frame_size, 3),
			Eigen::ComputeThinU | Eigen::ComputeThinV);
		MatrixXR U = svd.matrixU();
		MatrixXR S = svd.singularValues().asDiagonal();//sigma 
		MatrixXR Vt = svd.matrixV().transpose();

		
		
		//MatrixXR wk = S(0, 0) * V.col(0);
		MatrixXR wk = U.col(0) * S(0, 0);
		const MatrixXR wk_proj = proj(wk);
		const MatrixXR wk_proj_negative = proj(-wk);
		
		wk = wk_proj.norm() > wk_proj_negative.norm() ? wk_proj : wk_proj_negative;



		MatrixXR local_sup;
		get_local_support(max_v_idx, local_sup);

		MatrixXR s = 1 - local_sup.array();

		MatrixXR ck = wk.transpose() * (residual_X);
		Eigen::Map<MatrixXR> ck_3(ck.data(), 3, ck.cols() / 3);
		ck_3.array() *= s.transpose().replicate(ck_3.rows(), 1).array();

		C.row(k) = ck;
		W.col(k) = wk;

		MatrixXR tmp;
		tmp.noalias() = wk * ck; //outer product
		residual_X -= tmp;

		Eigen::Map<MatrixXR> qq(C.row(k).data(), 3, C.cols()/3);

		
	}

	MatrixXR tmpC = C * scale_factor_;
	for (int i = 0; i < C.rows(); i++){
	Mesh t;
	t.F = meanshape_mesh_->F;
	t.V = Eigen::Map<MatrixXR>(tmpC.row(i).data(), meanshape_mesh_->V.rows(), meanshape_mesh_->V.cols()) + meanshape_mesh_->V;
	igl::writeOBJ("test" + std::to_string(i) + ".obj", t.V, t.F);
	}
	std::cout << "first" << std::endl;

	// prepare auxiluary variables
	MatrixXR Lambda(component_num_, meanshape_mesh_->V.rows()); // K x N
	MatrixXR U(component_num_, meanshape_mesh_->V.size()); // K x N * 3



	// main global opt
	num_opt_ = 10;
	for (int i = 0; i < i < num_opt_; i++) {

		MatrixXR& Rflat = residual_X;
		for (int k = 0; k < C.rows(); k++) {
			Eigen::Block<MatrixXR, -1, 1, true> Ck = C.col(k);
			real_type Ck_norm = Ck.dot(Ck);
			if (Ck_norm <= 1.e-8) {
				W.col(k).array() = 0.0;
				continue;
			}


		}



	}





}

void SplocsSolver::add_triangle_mesh(std::string& name)
{
	static int s = 0;
	MatrixXR V;	
	MatrixXI F;
	igl::read_triangle_mesh(name, V, F);
	this->meshes_.emplace_back(V, F);
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
