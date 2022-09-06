#include <iostream>
#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/avg_edge_length.h>
#include <stdlib.h>
#include "splocs.h"

//
//// currently, it' just testing.
//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Surface_mesh.h>
//#include <CGAL/Polygon_mesh_processing/remesh.h>
//#include <CGAL/Polygon_mesh_processing/border.h>
//#include <CGAL/Polygon_mesh_processing/IO/polygon_mesh_io.h>
//#include <boost/iterator/function_output_iterator.hpp>
//#include <CGAL/Surface_mesh/Surface_mesh.h>
//#include <CGAL/boost/graph/graph_traits_Surface_mesh.h>


//22995
#include "awesome_header.h"

SplocsSolver::SplocsSolver()
	:matrix_W_(std::make_unique<MatrixXR>()), 
	matrix_X_(std::make_unique<MatrixXR>()),
	matrix_C_(std::make_unique<MatrixXR>()), 
	meanshape_mesh_(std::make_unique<Mesh>()), 
	d_min_(0.0), d_max_(1.0)
{}

void SplocsSolver::make_X_matrix_by_meshes_(){

	typedef Eigen::Map<MatrixXR, 0, Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>> DMap;

	const int mesh_num = meshes_.size();
	const int v_size = meshes_[0].V.rows();
	matrix_X_->resize(mesh_num, 3 * v_size);
	matrix_X_->setZero();
	VectorXR mean_shape;
	mean_shape.resize( matrix_X_->cols());
	mean_shape.setZero();
	std::for_each(meshes_.begin(), meshes_.end(), [&mean_shape, this,v_size, mesh_num, i = 0](Mesh& m) mutable {
		Eigen::Map<MatrixXR> x(m.V.data(), m.V.size(), 1);
		mean_shape += x;
		matrix_X_->row(i++) = x.transpose();
	});
	mean_shape /= mesh_num;
	//MatrixXR vv;
	//MatrixXI ff;
	//igl::read_triangle_mesh("D:/lab/2022/external_code/Mesh_DNN/splocs/xmean.obj", vv, ff);
	//std::cout << (meanshape_mesh_->V - vv).sum() << std::endl;
	//std::cout << (meanshape_mesh_->F - ff).sum() << std::endl;

	// check page 4
	// normalize part


 	matrix_X_->rowwise() -= mean_shape.transpose();

	for (int i = 0; i < matrix_X_->rows(); i++) {
		std::cout << matrix_X_->row(i).mean() << std::endl;;
	}


	meanshape_mesh_->V = std::move(Eigen::Map<MatrixXR>(mean_shape.data(), v_size, 3));
	meanshape_mesh_->F = meshes_[0].F;

	Eigen::Map<VectorXR> x_vec(matrix_X_->data(), matrix_X_->size());


	//std
	scale_factor_ = sqrt((x_vec.array() - x_vec.mean()).square().sum() / (x_vec.size()));
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
}

//
//void mesh_refinement(std::vector<Mesh>& meshes) {
//
//	typedef CGAL::Exact_predicates_inexact_constructions_kernel   K;
//	typedef CGAL::Surface_mesh<K::Point_3>                        CGAL_Mesh;
//	typedef boost::graph_traits<CGAL_Mesh>::halfedge_descriptor        halfedge_descriptor;
//	typedef boost::graph_traits<CGAL_Mesh>::edge_descriptor            edge_descriptor;
//	typedef CGAL_Mesh::Vertex_index vertex_descriptor;
//
//	namespace PMP = CGAL::Polygon_mesh_processing;
//	struct halfedge2edge
//	{
//		halfedge2edge(const CGAL_Mesh& m, std::vector<edge_descriptor>& edges)
//			: m_mesh(m), m_edges(edges)
//		{}
//		void operator()(const halfedge_descriptor& h) const
//		{
//			m_edges.push_back(edge(h, m_mesh));
//		}
//		const CGAL_Mesh& m_mesh;
//		std::vector<edge_descriptor>& m_edges;
//	};
//	
//
//	CGAL_Mesh m;
//	std::vector<K::Point_3> points;
//	std::vector<std::vector<std::size_t> > m_face;
//	Mesh m_mesh;
//	CGAL_Mesh mm;
//	mm.reserve(meshes[0].V.rows(), meshes[0].F.rows() * 3, meshes[0].F.rows());
//	
//	for (int i = 0; i < meshes.size(); i++) {
//		mm.clear_without_removing_property_maps();
//		for (int j = 0; j<meshes[i].V.rows(); j++) {
//			auto& t = meshes[i].V;
//			mm.add_vertex(K::Point_3(t(j, 0), t(j,1), t(j,2)));
//		}
//
//		for (int j = 0; j<meshes[i].F.rows(); j++) {
//			CGAL_Mesh::Vertex_index  a(meshes[i].F(j,0)),
//									b(meshes[i].F(j,1)),
//									c(meshes[i].F(j,2));
//			mm.add_face(a, b, c);
//		}
//
//
//		//for (vertex_descriptor vd : mm.vertices()) {
//		//	meshes[i].V(vd.idx(), 0) = mm.point(vd).x();
//		//	meshes[i].V(vd.idx(), 1) = mm.point(vd).y();
//		//	meshes[i].V(vd.idx(), 2) = mm.point(vd).z();
//		//}
//		//igl::write_triangle_mesh("fu_test.obj", meshes[i].V, meshes[i].F);
//		//std::cout << "end it" << std::endl;
//
//		typedef CGAL_Mesh::Face_range faces;
//
//		double target_edge_length = 0.04;
//		unsigned int nb_iter = 3;
//		std::cout << "Split border...";
//		std::vector<edge_descriptor> border;
//		PMP::border_halfedges(mm.faces(), mm, boost::make_function_output_iterator(halfedge2edge(mm, border)));
//		PMP::split_long_edges(border, target_edge_length, mm);
//		std::cout << "done." << std::endl;
//		PMP::isotropic_remeshing(mm.faces(), target_edge_length, mm,
//			CGAL::parameters::number_of_iterations(nb_iter)
//			.protect_constraints(true));
//
//
//		//MatrixXR tm(mm.number_of_vertices(),3);
//		//MatrixXI tf(mm.number_of_faces(),3);
//		//for (vertex_descriptor vd : mm.vertices()) {
//		//	tm(vd.idx(), 0) =  mm.point(vd).x();
//		//	tm(vd.idx(), 1) =  mm.point(vd).y();
//		//	tm(vd.idx(), 2) =  mm.point(vd).z();
//		//}
//
//
//		CGAL::IO::write_OBJ("test_reg"+std::to_string(i) + ".obj", mm);
//	}
//	
//}

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
	result = D.array().min(max).max( min );
	result = (result.array() - min) / ( max - min );
}

void SplocsSolver::precompute_local_support(const Mesh& m)
{

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
		mean_Vs += m.V ;

		real_type tmp_scale = (m.V.colwise().maxCoeff() - m.V.colwise().minCoeff()).cwiseAbs().maxCoeff();
		v_scale = v_scale > tmp_scale ? v_scale : tmp_scale;
	}
	mean_Vs /= meshes_.size();
	VectorXR mean_v = mean_Vs.colwise().mean();
	for (int i = 0; i < meshes_.size(); i++) {
		meshes_[i].V.rowwise() -= mean_v.transpose();
		meshes_[i].V /= v_scale;
	} 
	



}
void SplocsSolver::find_rbm_procrustes()
{


	Mesh& ref = meshes_[0];
	VectorXR t1  = ref.V.colwise().mean();
	MatrixXR topts_local  = ref.V.rowwise() - t1.transpose();
	for (int i = 0; i < meshes_.size(); i++) {
		Mesh& m = meshes_[i];
		VectorXR t0 = m.V.colwise().mean();
		MatrixXR frompts_local = m.V.rowwise() - t0.transpose();
		MatrixXR M = topts_local.transpose()* frompts_local;
		Eigen::JacobiSVD<MatrixXR> svd(M,
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

		const int rows = m.V.rows();
		MatrixXR m_tmp(rows, 4);

		m_tmp.block(0,0, rows, 3) = (m.V);
		m_tmp.block(0, 3, rows, 1).array() = 1.0;
		//m.V = ( T0* m_tmp.transpose()).transpose().eval().block(0, 0, rows, 3);
		m.V = ( m_tmp* T0.transpose()).block(0, 0, rows, 3);
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




		
MatrixXR proj(const MatrixXR& src) {
	MatrixXR x = src.cwiseMax(0.0);
	real_type max_x = x.maxCoeff();
	//if (max_x < F_EPS) {
	if (max_x == 0.0) {
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
	//result = t1v * t2v.transpose();
	result.noalias() = t1v * t2v.transpose();
	

	//result.noalias() = t1 * t2; //outer product
}


static MatrixXR l1l2(const MatrixXR& Lambda, const MatrixXR& x, real_type beta) {

	const int v_size = x.cols() / 3;
	MatrixXR x_len(x.rows(), v_size);
	for (int i = 0; i < x.rows(); i++) {
		DMap tmp = get_matrix_row_to_3dim( const_cast<MatrixXR&>(x), i);
		x_len.row(i) = tmp.rowwise().norm();
	}
	MatrixXR shrinkage = (1 - beta * Lambda.array() / x_len.array()).max(0.0);
	MatrixXR res;
	res.resizeLike(x);
	//res.setZero();
	for (int i = 0; i < x.rows(); i++) {
		DMap reshaped_x = get_matrix_row_to_3dim(*const_cast<MatrixXR*>(&x), i);
		MatrixXR mat = (reshaped_x.array().colwise() * shrinkage.row(i).transpose().array());
		res.row(i) = Eigen::Map<MatrixXR> ( mat.data(), 1, res.cols());
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
	//trace_22995(*matrix_X_, 0)
	//trace(*matrix_X_, 0,0)
	for (int k = 0; k < component_num_; k++) {
		using vec = Eigen::Matrix<real_type, Eigen::Dynamic, 1>;
		vec magnitude(v_size);
		magnitude.setZero();
		for (int f = 0; f < frame_size; f++) {
			DMap tmp = get_matrix_row_to_3dim(residual_X, f);
			magnitude += tmp.rowwise().squaredNorm();
			//magnitude += tmp.rowwise().norm();
		}

		int max_v_idx;
		real_type res = magnitude.maxCoeff(&max_v_idx);
		// see End of eigen_test.cpp
		Eigen::Map<MatrixXR, 0, Eigen::OuterStride<Eigen::Dynamic>>
			svd_data(residual_X.col(max_v_idx).data(), frame_size, 3, Eigen::OuterStride(frame_size*v_size));
  		Eigen::JacobiSVD<MatrixXR> svd(svd_data.transpose(),
			Eigen::ComputeThinU | Eigen::ComputeThinV);
		MatrixXR U = svd.matrixU();
		MatrixXR S = svd.singularValues().asDiagonal();//sigma 
		MatrixXR Vt = svd.matrixV().transpose();
		
		MatrixXR wk = Vt.row(0) * S(0, 0);
		const MatrixXR wk_proj = proj(wk);
		const MatrixXR wk_proj_negative = proj(-wk);
		
		wk = wk_proj.norm() > wk_proj_negative.norm() ? wk_proj : wk_proj_negative;
		MatrixXR local_sup;
		get_local_support(max_v_idx, local_sup, d_min_, d_max_);

		MatrixXR s = 1 - local_sup.array();
		
		//s.setOnes();//test

		MatrixXR ck = wk * (residual_X);
		Eigen::Map<MatrixXR> ck_3(ck.data(), v_size, 3);
		s.array()  /= (wk * (wk.transpose().eval()))(0, 0);
		Eigen::Map<VectorXR> s_tmp(s.data(), s.size());
		ck_3= ck_3.array().colwise() *s_tmp.array();

		C.row(k) = ck;
		W.col(k) = wk.transpose();

		MatrixXR tmp;
		outerProduct(tmp, wk, ck);
		residual_X -= tmp;
	}
	//trace_22995(C, 0)
	//trace(C, 0, 0);


	// prepare auxiluary variables
	MatrixXR Lambda(component_num_, meanshape_mesh_->V.rows()); // K x N
	MatrixXR U(component_num_, meanshape_mesh_->V.size()); // K x N * 3
	U.setZero();

	// main global opt
	for (int i = 0; i < num_opt_; i++) {
		MatrixXR& Rflat = residual_X;
		for (int k = 0; k < C.rows(); k++) {
			Eigen::Block<MatrixXR, 1, -1, false> Ck = C.row(k);
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
			outerProduct(tmp, W.col(k), Ck);
			Rflat -= tmp;
			// update spatially varying regularization strength
		}
		for (int k = 0; k < C.rows(); k++) {
			DMap Ck = get_matrix_row_to_3dim(C, k);
			int idx;
			Ck.rowwise().squaredNorm().maxCoeff(&idx);
			MatrixXR support_map;
			get_local_support(idx, support_map, d_min_, d_max_);
			Lambda.row(k) = sparsity_lambda_ * support_map.transpose();
		}
		//copy C
		MatrixXR Z = C;
		MatrixXR G = W.transpose().eval() * W;

		MatrixXR c = W.transpose() * (*matrix_X_);
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
		C_norm.setZero();
		for (int k = 0; k < C.rows(); k++) {
			DMap Ck =get_matrix_row_to_3dim(C, k);
			C_norm.row(k) = Ck.rowwise().norm().transpose();
		}
		real_type sparsity = (Lambda.array() * C_norm.array()).sum();
		real_type e = (residual_X.array().pow(2.0).sum() + sparsity);
		std::cout <<"=====" << i <<"=====" << std::endl;
		std::cout << "energy : " << e << std::endl;

	}
	// undo scaling
	C *= scale_factor_;
	std::cout << "save ..." << std::endl;
	for (int i = 0; i < C.rows(); i++) {
		Mesh tt; 
		tt.V = get_matrix_row_to_3dim(C, i) + meanshape_mesh_->V;
		//tt.V = Eigen::Map<Eigen::Matrix<real_type, -1,-1, Eigen::RowMajor>>(C.row(i).data(), C.row(0).size()/3, 3) + meanshape_mesh_->V;
		tt.F = meanshape_mesh_->F;
		igl::write_triangle_mesh("component_"+std::to_string(i) + ".obj", tt.V, tt.F);
		std::cout << i << " saved..." << std::endl;

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

inline DMap get_matrix_row_to_3dim(MatrixXR& mat, const int row_idx)
{
	const int total_rows = mat.rows();
	const int v_size = mat.cols()/3;
	assert(mat.cols() % 3 == 0 && " it's not 3 colums");
	return DMap(mat.row(row_idx).data(),
		v_size, 3,
		Eigen::Stride<Eigen::Dynamic, Eigen::Dynamic>
		((total_rows * v_size), total_rows)
	);
}
