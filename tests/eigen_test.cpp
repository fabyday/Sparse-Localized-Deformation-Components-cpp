#include <Eigen/Core>


#include <iostream>



int main() {
	using namespace Eigen;
	typedef Eigen::Matrix<int, -1, -1, Eigen::RowMajor> ROWMAT;
	typedef Eigen::Matrix<int, -1, -1, Eigen::ColMajor> COLMAT;
	std::vector<int> s{ 1,2,3,4,5,6,7,8,9 };
	std::vector<int> ss{ 1,2,3, 1,2,3, 1,2,3, 1,2,3 };
	COLMAT t(2,6); 
	t << 111, 112, 113, 121, 122, 123,
		211,212,213,221,222,223;

	std::vector<int> sss{ 11,11, 12,12, 13,13, //
						21,21, 22,22 , 23, 23};
	Map<ROWMAT> q(s.data(), 3, 3);
	Map<COLMAT> r(s.data(), 3, 3);
	Map<COLMAT, 0, Eigen::Stride<1,2>> qq(s.data(), 3, 3);
	Map<COLMAT, 0, Eigen::Stride<1,3>> qqq(ss.data(), 4, 3);
	Map<ROWMAT, 0, Eigen::Stride<6,2>> qqqq(t.data(), 2, 3);

	std::cout << t<< std::endl;
	std::cout << qqqq << std::endl;
	std::cout << "===="<< std::endl;
	COLMAT m1(3,1), m2(3,3),m3;
	m1 << 1, 2, 3;
	m2 << 4, 5, 6 ,7,8,9,10,11,12;

	m3 = m1.array().replicate(1,3) * m2.array();

	VectorXi ttm(3);

	ttm << 1, 2, 3;


	m2.colwise() -= ttm;
	
	std::cout << m3 << std::endl <<std::endl;
	std::cout << m2.colwise().sum() << std::endl << std::endl;
	std::cout <<m2<< std::endl << std::endl;
	std::cout << m2.rowwise().sum() << std::endl << std::endl;


	COLMAT sx(2, 3);
	sx << 1, 2, 3, 4, 5, 6;
	std::cout << Eigen::Map<ROWMAT>(sx.data(), 2, 3) << std::endl;
	std::cout << Eigen::Map<ROWMAT>(sx.transpose().data(), 2, 3) << std::endl;







}