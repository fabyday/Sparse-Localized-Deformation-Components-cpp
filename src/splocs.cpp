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
	// normalize part
	matrix_X_->rowwise() -= mean_shape.transpose();
	meanshape_mesh_->V = std::move(Eigen::Map<MatrixXR>(mean_shape.data(), v_size, 3));
	meanshape_mesh_->F = meshes_[0].F;


	Eigen::Map<VectorXR> x_vec(matrix_X_->data(), matrix_X_->size());

	std::cout << x_vec.array().block<5,1>(0,0).transpose() << std::endl << std::endl;
	std::cout << (x_vec.mean()) << std::endl<<std::endl;
	std::cout << (x_vec.array() - x_vec.mean()).block<5, 1>(0, 0).transpose() << std::endl<<std::endl;
	std::cout << ( ( x_vec.array() - x_vec.mean()).square()).block<5, 1>(0, 0).transpose() << std::endl<<std::endl;
	std::cout << ( ( x_vec.array() - x_vec.mean()).square().sum()) << std::endl<<std::endl;
	std::cout << ( ( x_vec.array() - x_vec.mean()).square().sum())/x_vec.size() << std::endl << std::endl;
	std::cout << ( ( x_vec.array() - x_vec.mean()).square().mean()) << std::endl << std::endl;
	std::cout << x_vec.size() << std::endl << std::endl;
	//std
	scale_factor_ = sqrt((x_vec.array() - x_vec.mean()).square().sum() / (x_vec.size()));
	std::cout << 1/ scale_factor_ << std::endl << std::endl;
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


void SplocsSolver::compile() {
	if (is_compiled_)
		return;

	is_compiled_ = true;
	write_mesh_debug(std::string("org"));
	scale_meshes();
	write_mesh_debug(std::string("scale"));
	find_rbm_procrustes();
	write_mesh_debug(std::string("rot"));
	make_X_matrix_by_meshes_();
	setup_W_C_lsup_();



}

void SplocsSolver::get_local_support(const int idx, MatrixXR& result, real_type min, real_type max)
{




	Eigen::Matrix<real_type, -1, 1> D;
	Eigen::VectorXi vid(1, 1); vid(0) = idx;

	igl::heat_geodesics_solve(data_, vid, D);
	result = D.array().min(max).max(min);
	result = (result.array() - min) / (max- min);

}

void SplocsSolver::precompute_local_support(const Mesh& m)
{
	//real_type t = std::pow(igl::avg_edge_length(m.V, m.F), 2);
	real_type t = 10*std::pow(igl::avg_edge_length(m.V, m.F), 2);
	if (!igl::heat_geodesics_precompute(m.V, m.F, t, data_)) {
		std::cerr << "Error: heat_geodesics_precompute failed." << std::endl;
		exit(EXIT_FAILURE);
	}
}
void SplocsSolver::scale_meshes() {
	Mesh& ref = meshes_[0];
	VectorXR t0 = ref.V.colwise().mean();
	MatrixXR mean_Vs;
	mean_Vs.resizeLike(ref.V);
	mean_Vs.setZero();
	real_type v_scale = 0;
	for (int i = 0; i < meshes_.size(); i++) {
		Mesh& m = meshes_[i];
		mean_Vs += m.V / meshes_.size();

		real_type tmp_scale = (m.V.colwise().maxCoeff() - m.V.colwise().minCoeff()).cwiseAbs().maxCoeff();
		v_scale = v_scale > tmp_scale ? v_scale : tmp_scale;
	}
	VectorXR mean_v = mean_Vs.colwise().mean();
	for (int i = 0; i < meshes_.size(); i++) {
		meshes_[i].V.array().rowwise() -= mean_v.transpose().array();
		meshes_[i].V /= v_scale;
	} 



}
void SplocsSolver::find_rbm_procrustes()
{


	Mesh& ref = meshes_[0];
	VectorXR t0 = ref.V.colwise().mean();
	MatrixXR frompts_local = ref.V.array().rowwise() - t0.transpose().array();
	for (int i = 0; i < meshes_.size(); i++) {
		Mesh& m = meshes_[i];
		VectorXR t1 = m.V.colwise().mean();
		MatrixXR topts_local = m.V.array().rowwise() - t1.transpose().array();

		MatrixXR M = topts_local.transpose()* frompts_local;
		Eigen::JacobiSVD<MatrixXR> svd(M.transpose(),
			Eigen::ComputeThinU | Eigen::ComputeThinV);
		MatrixXR U = svd.matrixU();
		MatrixXR S = svd.singularValues().asDiagonal();//sigma 
		MatrixXR Vt = svd.matrixV().transpose();
		MatrixXR R = U * Vt;
		if (R.determinant() < 0.0) {
			R *= -1;
		}
		MatrixXR T0 = MatrixXR::Identity(4, 4);
		T0.block<3, 3>(0, 0) = R;
		T0.block<3, 1>(0, 3) = t1 - R*t0;
		
		MatrixXR m_tmp(m.V.rows(), 4); 
		const int rows = m.V.rows();
		m_tmp.block(0,0, rows, 3) = std::move(m.V);
		m_tmp.block(0, 3, rows, 1).array() = 1.0;
		m.V = ( T0* m_tmp.transpose()).transpose().block(0, 0, rows, 3);
		
	}
}




void SplocsSolver::solve(int component_num , int num_iter_max, int num_admm_iterations,
						real_type d_min, real_type d_max,
						real_type sparsity_lambda, real_type rho)
{
	if (component_num_ != component_num) {
		is_compiled_ = false;
		component_num_ = component_num;
	}
	d_min_ = d_min;
	d_max_ = d_max;
	sparsity_lambda_ = sparsity_lambda;
	num_opt_ = num_iter_max;
	num_admm_iterations_ = num_admm_iterations;
	rho_ = rho;

	compile();

	phase1(*matrix_C_, *matrix_W_);


}




//
//void SplocsSolver::solve(int component_num)
//{
//	solve(component_num, d_min_, d_max_);
//}



		
MatrixXR proj(const MatrixXR& src) {
	MatrixXR x = src.cwiseMax(0.0);
	real_type max_x = x.maxCoeff();
	if (max_x < F_EPS) {
		return x;
	}
	else {
		return x / max_x;
	}



	//return src.array() / src.cwiseAbs().maxCoeff();
}



void save_component_to_mesh(const MatrixXR& component) {

}


static void outerProduct(MatrixXR& result, const MatrixXR& t1, const MatrixXR& t2) {

	// t1, t2 col-vector
	// 
	//result;
	Eigen::Map<MatrixXR> t1v(const_cast<double*>(t1.data()), t1.size(),1);
	Eigen::Map<MatrixXR> t2v(const_cast<double*>(t2.data()), t2.size(),1);
	result = t1v * t2v.transpose();

	//result.noalias() = t1 * t2; //outer product
}


static MatrixXR l1l2(const MatrixXR& Lambda, const MatrixXR& x, real_type beta) {

	const int v_size = x.cols() / 3;
	MatrixXR x_len(x.rows(), v_size);
	typedef Eigen::Map<MatrixXR, 0, Eigen::OuterStride<Eigen::Dynamic>> DMap;


	for (int i = 0; i < x.rows(); i++) {
		DMap tmp(const_cast<real_type*>(x.row(i).data()), v_size, 3, Eigen::OuterStride(v_size));
		x_len.row(i) = tmp.rowwise().squaredNorm();
		//Eigen::Map<Eigen::Matrix<real_type, -1, -1, Eigen::RowMajor>> tmp(const_cast<real_type*>(x.row(i).data()), x.row(i).size() / 3, 3);
		//x_len.row(i) = tmp.rowwise().squaredNorm();
	}

	MatrixXR shrinkage = (1 - beta * Lambda.array() / x_len.array()).max(0.0).matrix();
	
	MatrixXR res;
	res.resizeLike(x);
	for (int i = 0; i < x.rows(); i++) {
		DMap tmp(const_cast<real_type*>(x.row(i).data()), v_size, 3, Eigen::OuterStride(v_size));
		Eigen::Map<VectorXR> tmp_shrinkage(shrinkage.row(i).data(), shrinkage.cols(), 1);
		MatrixXR&& st = tmp.array().colwise()* tmp_shrinkage.array();
		res.row(i) = Eigen::Map<MatrixXR>(st.data(), 1, st.size());
		//Eigen::Map<Eigen::Matrix<real_type, -1, -1, Eigen::RowMajor>> tmp(const_cast<real_type*>(x.row(i).data()), x.row(i).size() / 3, 3);
		//auto ts = Eigen::Map<MatrixXR>(shrinkage.row(i).data(), shrinkage.row(i).size(), 1);// .array().matrix();
		//MatrixXR st = (tmp.array() * ts.replicate(1,3).array()).matrix();
		//res.row(i) = Eigen::Map<MatrixXR>(const_cast<real_type*>( st.data() ), 1, res.cols());
	}


	return res;
	//return x * shrinkage;

}



void SplocsSolver::phase1(MatrixXR& C, MatrixXR& W){
	
	// matrix X : memory layout [-X-Y-Z-]
	//							[   .   ]
	//							[   .   ]
	//							[-X-Y-Z-]


	const int v_size = meshes_[0].V.rows();
	const int frame_size = meshes_.size();
	//copy
	MatrixXR residual_X = *matrix_X_;

	precompute_local_support(*meanshape_mesh_);
	C.resize(component_num_, v_size*3);
	W.resize(frame_size, component_num_);

	for (int k = 0; k < component_num_; k++) {

		using vec = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
		vec magnitude(v_size);
		magnitude.setZero();
		for (int f = 0; f < frame_size; f++) {
			Eigen::Map<Eigen::Matrix<real_type, -1,-1>, 0>
				tmp(residual_X.row(f).data(),
					v_size, 3);
			//std::cout << tmp.row(0) << std::endl;
			//std::cout <<  tmp(0,0)*tmp(0,0)<<" " << tmp(0, 1) * tmp(0, 1) << " " << tmp(0, 2) * tmp(0, 2) << std::endl;
			//std::cout << tmp.rowwise().squaredNorm()(0) << std::endl;
			magnitude += tmp.rowwise().squaredNorm();
		}

		int max_v_idx;
		magnitude.maxCoeff(&max_v_idx);


		/*
		Eigen::Map<MatrixXR, 0, Eigen::OuterStride<Eigen::Dynamic>> 
			svd_data(residual_X.data(), frame_size, 3, Eigen::OuterStride(v_size));*/
		// see End of eigen_test.cpp
		Eigen::Map<MatrixXR, 0, Eigen::OuterStride<Eigen::Dynamic>>
			svd_data(residual_X.col(max_v_idx).data(), frame_size, 3, Eigen::OuterStride(frame_size*v_size));

  		Eigen::JacobiSVD<MatrixXR> svd(svd_data.transpose(),
			Eigen::ComputeThinU | Eigen::ComputeThinV);
		MatrixXR U = svd.matrixU();
		MatrixXR S = svd.singularValues().asDiagonal();//sigma 
		MatrixXR Vt = svd.matrixV().transpose();
		//std::cout << U << std::endl;
		//std::cout << S << std::endl;
		//std::cout << Vt << std::endl;
		
		
		//MatrixXR wk = U.col(0) * S(0, 0);
		MatrixXR wk = Vt.row(0) * S(0, 0);
		const MatrixXR wk_proj = proj(wk);
		const MatrixXR wk_proj_negative = proj(-wk);
		
		wk = wk_proj.norm() > wk_proj_negative.norm() ? wk_proj : wk_proj_negative;

		MatrixXR local_sup;
		get_local_support(max_v_idx, local_sup, d_min_, d_max_);

		MatrixXR s = 1 - local_sup.array();

		//MatrixXR ck = wk.transpose() * (residual_X);
		MatrixXR ck = wk * (residual_X);

		//std::cout << ck.block<1, 10> (0,0) << std::endl;
		Eigen::Map<MatrixXR> ck_3(ck.data(), v_size, 3);
		//std::cout << s.maxCoeff() << std::endl;
		//ck_3.array() *= s.transpose().replicate(ck_3.rows(), 1).array() / (wk*(wk.transpose().eval()))(0,0);
		s.array()  /= (wk * (wk.transpose().eval()))(0, 0);
		Eigen::Map<VectorXR> s_tmp(s.data(), s.size());
		 ck_3= ck_3.array().colwise() *s_tmp.array();
		 //std::cout << ck.maxCoeff() << std::endl;
		//C.row(k) = ck;
		//std::cout << ck.block<1, 10> (0,0) << std::endl;
		C.row(k) = ck;
		//W.col(k) = wk;
		W.col(k) = wk.transpose();

		MatrixXR tmp;
		//tmp.noalias() = wk * ck; //outer product
		outerProduct(tmp, wk, ck);
		residual_X -= tmp;

		//Eigen::Map<MatrixXR> qq(C.row(k).data(), 3, C.cols()/3);
	}
  /*	MatrixXR tmpC = C * scale_factor_;
	for (int i = 0; i < C.rows(); i++){
	Mesh t;
	t.F = meanshape_mesh_->F;
	t.V = Eigen::Map<MatrixXR>(tmpC.row(i).data(), meanshape_mesh_->V.rows(), meanshape_mesh_->V.cols()) + meanshape_mesh_->V;
	igl::writeOBJ("test" + std::to_string(i) + ".obj", t.V, t.F);
	}
	std::cout << "first" << std::endl;*/

	// prepare auxiluary variables
	MatrixXR Lambda(component_num_, meanshape_mesh_->V.rows()); // K x N
	MatrixXR U(component_num_, meanshape_mesh_->V.size()); // K x N * 3

	// main global opt
	for (int i = 0; i < num_opt_; i++) {

		MatrixXR& Rflat = residual_X;
		for (int k = 0; k < C.rows(); k++) {
			Eigen::Block<MatrixXR, 1, -1, false> Ck = C.row(k);
			//std::cout << Ck.rows() << ", " << Ck.cols() << std::endl;
			real_type Ck_norm = Ck.dot(Ck);
			if (Ck_norm <= F_EPS) {
				W.col(k).array() = 0.0;
				continue;
			}
			// block coordinate descent update
			MatrixXR tmp;
			outerProduct(tmp, W.col(k), Ck);
			Rflat += tmp;
			MatrixXR opt = Rflat* Ck.transpose() / Ck_norm;
			W.col(k) = proj(opt);
			tmp.setZero();
			outerProduct(tmp, W.col(k), Ck);
			Rflat -= tmp;
			// update spatially varying regularization strength
		}
		
		for (int k = 0; k < C.rows(); k++) {
			//Eigen::Block<MatrixXR, 1, -1, false> Ck = 
			Eigen::Map<MatrixXR> Ck(C.row(k).data(), v_size, 3);
			int idx;
			Ck.array().pow(2.0).matrix().rowwise().sum().maxCoeff(&idx);
			MatrixXR support_map;
			get_local_support(idx, support_map, d_min_, d_max_);
			Lambda.row(k) = sparsity_lambda_ * support_map.transpose();
		}
		//copy C
		MatrixXR Z = C;
		MatrixXR G = W.transpose().eval() * W;
		


		//MatrixXR c = W.transpose() * Eigen::Map<MatrixXR>(matrix_X_->data(), matrix_X_->rows(), matrix_X_->cols());
		MatrixXR c = W.transpose() * *matrix_X_;
		Eigen::LLT<MatrixXR> cho_slvr;
		cho_slvr.compute(G + rho_ * MatrixXR::Identity(G.rows(), G.cols()));
		for (int admm_it = 0; admm_it < num_admm_iterations_; admm_it++) {
			C = cho_slvr.solve(c + rho_ * (Z - U));
			Z = l1l2(Lambda, C + U, 1. / rho_);
			U = U + C - Z;
		}
		C = Z;
		residual_X = *(matrix_X_) - W * C;

		
		MatrixXR C_norm(C.rows(), static_cast<int>(C.cols() / 3));
		//Eigen::Map<MatrixXR, 0, Eigen::Stride<1, 3>>   sss;
		for (int i = 0; i < C.rows(); i++) {
			Eigen::Map<Eigen::Matrix<real_type,-1,-1, Eigen::RowMajor>> tmp(C.row(i).data(), C.row(i).size() / 3, 3);
			C_norm.row(i) = tmp.rowwise().norm();
		}

		real_type sparsity = (Lambda.array() * C_norm.array()).sum();
		real_type e = (residual_X.array().pow(2.0).sum() + sparsity);
		std::cout << "energy : " << e << std::endl;

	}

	// undo scaling
	C *= scale_factor_;
	for (int i = 0; i < C.rows(); i++) {
		Mesh tt;
		tt.V = Eigen::Map<Eigen::Matrix<real_type, -1,-1, Eigen::RowMajor>>(C.row(i).data(), C.row(0).size()/3, 3) + meanshape_mesh_->V;
		tt.F = meanshape_mesh_->F;
		igl::write_triangle_mesh("testtest"+std::to_string(i) + ".obj", tt.V, tt.F);
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
