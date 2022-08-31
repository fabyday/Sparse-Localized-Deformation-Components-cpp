#include <splocs.h>
#include <igl/opengl/glfw/Viewer.h>




int main() {

	SplocsSolver slvr;

	slvr.add_triangle_meshes(std::vector<std::string> {
			big_face"head-01-anger.obj",
			big_face"head-02-cry.obj",
			big_face"head-03-fury.obj",
			big_face"head-04-grin.obj",
			big_face"head-05-laugh.obj",
			big_face"head-06-rage.obj",
			big_face"head-07-sad.obj",
			big_face"head-08-smile.obj",
			big_face"head-09-surprise.obj",
			big_face"head-reference.obj",
	});

	//slvr.solve(11);
	slvr.solve(50);





	return 0;
}	