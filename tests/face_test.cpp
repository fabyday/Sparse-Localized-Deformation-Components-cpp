#include <splocs.h>
#include <igl/opengl/glfw/Viewer.h>




int main() {

	SplocsSolver slvr;


	slvr.add_triangle_meshes(std::vector<std::string> { 
		//face205"face-reference.obj", 
		face205"face-02-cry.obj" ,
		face205"face-03-fury.obj" ,
		face205"face-05-laugh.obj" ,
		face205"face-06-rage.obj" ,
		face205"face-07-sad.obj" ,
		face205"face-08-smile.obj" ,
		face205"face-09-surprise.obj" ,
		});

	slvr.solve(50);





	return 0;
}