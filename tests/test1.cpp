#include <splocs.h>
#include <igl/opengl/glfw/Viewer.h>




int main() {

	SplocsSolver slvr;


	//slvr.add_triangle_meshes(std::vector<std::string> { 
	//	face205"face-reference.obj", 
	//	face205"face-01-anger.obj" ,
	//	});

	slvr.add_triangle_meshes(std::vector<std::string> {
		navitan"cat_000001.obj",
		navitan"cat_000012.obj",
		navitan"cat_000017.obj",
		navitan"cat_000030.obj",
	});
	slvr.solve(50);





	return 0;
}