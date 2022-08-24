#include <splocs.h>
#include <igl/opengl/glfw/Viewer.h>




int main() {

	SplocsSolver slvr;

	slvr.add_triangle_meshes(std::vector<std::string> { 
		face205"face-01-anger.ply",
		face205"face-02-cry.ply",
		face205"face-03-fury.ply",
		face205"face-04-grin.ply",
		face205"face-05-laugh.ply",
		face205"face-06-rage.ply",
		face205"face-07-sad.ply",
		face205"face-08-smile.ply",
		face205"face-09-surprise.ply",
		face205"face-reference.ply"
		});

	slvr.solve(50);





	return 0;
}