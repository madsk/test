#include "driver_state.h"
#include <cstring>
#include <vector>
#include <limits>
#include <algorithm>


// float min(float a, float b, float c)
// {
// 	float min;
// 	(a <= b) ? min = a : min = b;
// 	if (min > c) { min = c; }
// 	return min;
// }

// float max(float a, float b, float c)
// {
// 	float max;
// 	(a >= b) ? max = a : max = b;
// 	if (max < c) { max = c; }
// 	return max;
// }

driver_state::driver_state()
{
}

driver_state::~driver_state()
{
    delete [] image_color;
    delete [] image_depth;
}

// This function should allocate and initialize the arrays that store color and
// depth.  This is not done during the constructor since the width and height
// are not known when this class is constructed.
void initialize_render(driver_state& state, int width, int height)
{
    state.image_width=width;
    state.image_height=height;
    state.image_color=0;
    state.image_depth=0;
	state.image_depth = new float[width*height];
    //std::cout<<"TODO: allocate and initialize state.image_color and state.image_depth."<<std::endl;

	state.image_color = new pixel[height*width]; 
	
	for(int i = 0; i < width*height; i++) {
			state.image_color[i] = make_pixel(0, 0, 0); //initialize to black
			state.image_depth[i] = std::numeric_limits<float>::max(); //1
			//state.image_depth[i] = 1; //1
	}
	
/*Implement initialize_render in driver_state.cpp by allocating the memory for color_image. Initialize all
the pixels in color_image to black. We will not be using depth_image until we implement z-buffer, so
you can ignore it for now. Make sure your code compiles and run without issues on valgrind. You can
compile the code using scons and you can run on test 00.txt using ./driver -i tests/00.txt.*/
}

// This function will be called to render the data that has been stored in this class.
// Valid values of type are:
//   render_type::triangle - Each group of three vertices corresponds to a triangle.
//   render_type::indexed -  Each group of three indices in index_data corresponds
//                           to a triangle.  These numbers are indices into vertex_data.
//   render_type::fan -      The vertices are to be interpreted as a triangle fan.
//   render_type::strip -    The vertices are to be interpreted as a triangle strip.
void render(driver_state& state, render_type type)
{
    //std::cout<<"TODO: implement rendering."<<std::endl;
	
	switch(type) {
		
		default:
		break;
		case render_type:: indexed:
		break;
		case render_type:: fan:
		break;
		case render_type:: strip:
		break;
		
		case render_type:: triangle: {
	
			data_geometry g_array[3];
			data_vertex v_array[3];
			const data_geometry* g[3];
			int k = 0;

			//int num_v = state.num_vertices / 3;
			
			//std::cout << state.num_vertices << std::endl;

			//int test = state.num_triangles;
			for(int i = 0; i < state.num_vertices / 3; i++) { //works
				for(int j = 0; j < 3; j++) {
					v_array[j].data = &state.vertex_data[k];
					g_array[j].data = v_array[j].data;

					state.vertex_shader(v_array[j], g_array[j], state.uniform_data);
					g[j] = &g_array[j];

					k += state.floats_per_vertex;
					//std::cout << j << std::endl;
				}
				//rasterize_triangle(state, g);
				clip_triangle(state, g, 0); //infinite loop
			}
			

			// for(int i = 0; i < state.num_vertices / 3; i++) { //alt way, also works
			// 	for(int j = 0; k < state.num_vertices * state.floats_per_vertex; j++, k += state.floats_per_vertex) {
			// 		v_array[j].data = &state.vertex_data[k];
			// 		g_array[j].data = v_array[j].data;
			// 	}
			// 	for(unsigned j = 0; j < 3; j++) {
			// 		state.vertex_shader(v_array[j], g_array[j], state.uniform_data);
			// 		g[j] = &g_array[j];
			// 	}

			// 	rasterize_triangle(state, g);
			// }

			} //end case
			break;
	}
}

// This function clips a triangle (defined by the three vertices in the "in" array).
// It will be called recursively, once for each clipping face (face=0, 1, ..., 5) to
// clip against each of the clipping faces in turn.  When face=6, clip_triangle should
// simply pass the call on to rasterize_triangle.
void clip_triangle(driver_state& state, const data_geometry* in[3],int face)
{
	//new data_geometry[3]; //new
	//if f(p) < 0 is inside plane f(p) > 0 is outside plane
	//implicit eq for plane through point q and normal n: f(p) = n.D = 0 (D is (p-q)
	// solving for t, we get: t = n . a + D / n . (a-b)
	//Implement clipping against the near and far planes. (test 10)
	//Note that full perspective-correct clipping is not required to pass this test.
	//image_depth to determine if near or far?

    if(face==6)
    {
        rasterize_triangle(state, in);
        return;
    }

    //std::cout<<"TODO: implement clipping. (The current code passes the triangle through without clipping them.)"<<std::endl;
	//near, far, top, right, left, bottom
	
	vec4 a = in[0] -> gl_Position; //each contains points x, y, z, w
	vec4 b = in[1] -> gl_Position; //0, 1, 2, 3
	vec4 c = in[2] -> gl_Position;

	//inside clipping plane
	if(b[2] >= -b[3] && c[2] >= -c[3] && a[2] >= -a[3]) { //bot & top = y plane, near & far = z plane

	//std::cout << "test" << std::endl;
	//float a1;
	//float b2;

	//a1 = (-1*b[3]-b[2]) / a[2]+a[3]-b[3]-b[2]); 
	//case for everything is inside, no calculation necessary
	//leave as is
	
	//vec4 p = a1 * a[2] + (1-a1)*b[2];

	//p(α) = αa + (1 − α)b //two points only, only one line segment
	//< -w not is_unsigned >= -w near plane
	//
	//x < w is right
	}
    clip_triangle(state,in,face+1);
}

//i made new data geometry objects and copied the ones from in to the new data geometry objects,
//then changed the gl position on the interpolated vertices. that should work right?

// Rasterize the triangle defined by the three vertices in the "in" array.  This
// function is responsible for rasterization, interpolation of data to
// fragments, calling the fragment shader, and z-buffering.

void rasterize_triangle(driver_state& state, const data_geometry* in[3])
{
    //std::cout<<"TODO: implement rasterization"<<std::endl;
	//divide by w position, 
	//x, y, z, w
	//NDC -> divide by w
	//barycentric coordinate -> triangle intersection with a pixel
	float Ax = ((in[0]->gl_Position[0] / in[0]->gl_Position[3]) * (state.image_width / 2)) + (((state.image_width / 2) - 0.5));
	float Ay = ((in[0]->gl_Position[1] / in[0]->gl_Position[3]) * (state.image_height / 2)) + (((state.image_height / 2) - 0.5)); 
	float Bx = ((in[1]->gl_Position[0] / in[1]->gl_Position[3]) * (state.image_width / 2)) + (((state.image_width / 2) - 0.5)); 
	float By = ((in[1]->gl_Position[1] / in[1]->gl_Position[3]) * (state.image_height / 2)) + (((state.image_height / 2) - 0.5));
	float Cx = ((in[2]->gl_Position[0] / in[2]->gl_Position[3]) * (state.image_width / 2)) + (((state.image_width / 2) - 0.5)); 
	float Cy = ((in[2]->gl_Position[1] / in[2]->gl_Position[3]) * (state.image_height / 2)) + (((state.image_height / 2) - 0.5)); 

	//get z coordinate for z buffering, one for each vertex
	float z0 = ((in[0]->gl_Position[2] / in[0]->gl_Position[3]) * (state.image_width / 2)) + (((state.image_width / 2) - 0.5));
	float z1 = ((in[1]->gl_Position[2] / in[1]->gl_Position[3]) * (state.image_width / 2)) + (((state.image_width / 2) - 0.5));
	float z2 = ((in[2]->gl_Position[2] / in[2]->gl_Position[3]) * (state.image_width / 2)) + (((state.image_width / 2) - 0.5));
	
	float w0 = in[0]->gl_Position[3] * (state.image_width / 2);
	float w1 = in[1]->gl_Position[3] * (state.image_width / 2);
	float w2 = in[2]->gl_Position[3] * (state.image_width / 2);

	vec2 v0 = {Ax, Ay};
	vec2 v1 = {Bx, By};
	vec2 v2 = {Cx, Cy};
	
	float area_abc = (v2[0] - v0[0]) * (v1[1] - v0[1]) - (v2[1] - v0[1]) * (v1[0] - v0[0]);
	//vec3 area_abc = {(v2[0] - v0[0]) * (v1[1] - v0[1]) - (v2[1] - v0[1]) * (v1[0] - v0[0]), 0, 0};
	//for (int w = minX; w < maxX + 1; w++) {                                      // Iterate trhough all pixels in the square
        //for (int h = minY; h < maxY + 1; h++) {
	
	for(int i = 0; i < state.image_width; i++) {
		for(int j = 0; j < state.image_height; j++) {

			//two vertices and a point, area of subtriangle w/n triangle
			vec2 p = {i, j};
			
			//barycentric coordinates - points where the intersection occured
			float alpha = ((p[0] - v1[0]) * (v2[1] - v1[1]) - (p[1] - v1[1]) * (v2[0] - v1[0])) / area_abc;
			float beta = ((p[0] - v2[0]) * (v0[1] - v2[1]) - (p[1] - v2[1]) * (v0[0] - v2[0])) / area_abc;
			float gamma = ((p[0] - v0[0]) * (v1[1] - v0[1]) - (p[1] - v0[1]) * (v1[0] - v0[0])) / area_abc;
			//vec3 b_prime = {area_abc[0] / area[0], area[2]/ area[0], 0};
			//b_prime[2] = 1 - b_prime[0] - b_prime[1];

			float new_alpha;
			float new_beta;
			float new_gamma;
			
			if(alpha >= 0 && beta >= 0 && gamma >= 0) { //to make sure that they are inside the triangle
			float test2 = alpha * z0 + beta * z1 + gamma * z2; //new
			
				if(test2 < state.image_depth[i + j * state.image_width] /*&& test2 < 1 && test2 > -1*/) { //new
					data_fragment data_frag; //new
					data_frag.data = new float[state.floats_per_vertex]; //iterate each float and find interpolating rule of each
					data_output data_out;

					for(int k = 0; k < state.floats_per_vertex; k++) { //new
							if(state.interp_rules[k] == interp_type::flat) { //how do we know which rule it is
								//which rule that the element in data geo matches with k
								data_frag.data[k] = in[0]->data[k]; //data ptr is first float in v data? 3rd element in first vertex data
							}
							else if(state.interp_rules[k] == interp_type::noperspective) {
								data_frag.data[k] = alpha * in[0]->data[k] + beta * in[1]->data[k] + gamma * in[2]->data[k];
							}
							else if(state.interp_rules[k] == interp_type::smooth) { //test 19
								float m = (alpha / w0) + (beta / w1) + (gamma / w2);
								new_alpha = (alpha / w0) / m;
								new_beta = (beta / w1) / m;
								new_gamma = (gamma / w2) / m;
								data_frag.data[k] = (new_alpha * in[0]->data[k] + new_beta * in[1]->data[k]) + (new_gamma * in[2] -> data[k]);
							}
					} //1-9 17-18 21
				//z buffering if triangles overlap know which one is colored over the over
				//image_depth = 0 initialize to large number for z buffer set to num_limits large number
				state.fragment_shader(data_frag, data_out, state.uniform_data);
				state.image_color[i + j * state.image_width] = make_pixel(data_out.output_color[0]*255, data_out.output_color[1]*255, data_out.output_color[2]*255);
				state.image_depth[i + j * state.image_width] = test2;
				}
			//}
			//}
			}
		}
	}	
}














