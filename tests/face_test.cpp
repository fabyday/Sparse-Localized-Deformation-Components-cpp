#include <splocs.h>
#include <igl/opengl/glfw/Viewer.h>


int main() {

	SplocsSolver slvr;
	//add mesh
	slvr.add_triangle_meshes(std::vector<std::string> { 
		face205"face-01-anger.ply",
		//face205"face-02-cry.ply",
		//face205"face-03-fury.ply",
		face205"face-04-grin.ply",
		face205"face-05-laugh.ply",
		face205"face-06-rage.ply",
		face205"face-07-sad.ply",
		face205"face-08-smile.ply",
		face205"face-09-surprise.ply",
		//face205"face-reference.ply"
		});

	//(Optional) select verbose level for energy term. default value is false.
	slvr.set_info_level(true);
	//solving
	slvr.solve(50);
	//get component
	std::vector<Mesh> components = slvr.get_mesh_components();
	//saving
	for (int i = 0; i < components.size(); i++) {
		std::cout << i << "-th components saved..." << std::endl;
		components[i].save_as("component_" + std::to_string(i) + ".obj");
	}


	return 0;
}