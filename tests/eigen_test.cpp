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
	Map<ROWMAT, 0, Eigen::Stride<2,2>> qqqq(t.data(), 4, 3);

	std::cout << t<< std::endl;
	std::cout << qqqq << std::endl;
}