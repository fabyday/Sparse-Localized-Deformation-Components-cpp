#include <splocs.h>
#include <igl/opengl/glfw/Viewer.h>




int main() {

	SplocsSolver slvr;


	slvr.add_triangle_meshes(std::vector<std::string> { 
		face205"aaa.ply", 
		face205"face_mesh_000306.ply" , 
		face205"face_mesh_000759.ply"});
	slvr.solve(50);





	return 0;
}