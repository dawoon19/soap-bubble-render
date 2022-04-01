#include "camera.h"

#include <iostream>
#include <sstream>
#include <fstream>

#include "CGL/misc.h"
#include "CGL/vector2D.h"
#include "CGL/vector3D.h"

using std::cout;
using std::endl;
using std::max;
using std::min;
using std::ifstream;
using std::ofstream;

namespace CGL {

using Collada::CameraInfo;

Ray Camera::generate_ray_for_thin_lens(double x, double y, double rndR, double rndTheta) const {

  // TODO Project 3-2: Part 4
  // compute position and direction of ray from the input sensor sample coordinate.
  // Note: use rndR and rndTheta to uniformly sample a unit disk.
  
    // Red segment
    float hFov_rad = hFov * (PI / 180);
    float vFov_rad = vFov * (PI / 180);

    float max_x = tan(0.5 * hFov_rad);
    float max_y = tan(0.5 * vFov_rad);

    float offset_x_img = x - 0.5;
    float offset_y_img = y - 0.5;

    float scale_x = max_x / 0.5;
    float scale_y = max_y / 0.5;

    float x_cam = scale_x * offset_x_img;
    float y_cam = scale_y * offset_y_img;

    Vector3D d = Vector3D(x_cam, y_cam, -1);
    d.normalize();

    Vector3D pLens = Vector3D(lensRadius * sqrt(rndR) * cos(rndTheta),
                            lensRadius * sqrt(rndR) * sin(rndTheta),0);


    double t = dot(Vector3D(0, 0, -focalDistance), Vector3D(0, 0, 1)) / dot(d, Vector3D(0, 0, 1));
    Vector3D pFocus = t * d;

    Vector3D blue_sgmt = pFocus - pLens;
    blue_sgmt.normalize();

    Ray finalRay = Ray(pLens, blue_sgmt);
    finalRay.o = c2w * finalRay.o;
    finalRay.o += pos;

    finalRay.d = c2w * finalRay.d;
    finalRay.d.normalize();

    finalRay.min_t = nClip;
    finalRay.max_t = fClip;

    return finalRay;
}


} // namespace CGL
