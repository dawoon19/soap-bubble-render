#include "bsdf.h"

#include <algorithm>
#include <iostream>
#include <utility>

#include "application/visual_debugger.h"

using std::max;
using std::min;
using std::swap;

namespace CGL {

// Mirror BSDF //

Vector3D MirrorBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D MirrorBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Implement MirrorBSDF
    reflect(wo, wi);
    *pdf = 1;
    return reflectance / abs_cos_theta(*wi);
}

void MirrorBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Mirror BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    ImGui::TreePop();
  }
}

// Microfacet BSDF //

double MicrofacetBSDF::G(const Vector3D wo, const Vector3D wi) {
  return 1.0 / (1.0 + Lambda(wi) + Lambda(wo));
}

double MicrofacetBSDF::D(const Vector3D h) {
  // TODO Project 3-2: Part 2
  // Compute Beckmann normal distribution function (NDF) here.
  // You will need the roughness alpha.
    Vector3D n = Vector3D(0,0,1);
    double tan_h = cross(h, n).norm() / dot(h,n);
    double cos_h = dot(h, n) / (h.norm() * n.norm());
    
    return exp(-pow(tan_h,2) / pow(alpha,2) ) / (PI * pow(alpha, 2) * pow(cos_h,4));
//  return 1.0;
}

Vector3D MicrofacetBSDF::F(const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Compute Fresnel term for reflection on dielectric-conductor interface.
  // You will need both eta and etaK, both of which are Vector3D.
//    double n = 0.093390;
//    double n = 0.35678;
//    double k = 8.9358;
//    double k = 8.2614;
//    Vector3D cobalt_n = Vector3D(2.1849,2.05,1.7925);
//    Vector3D cobalt_k = Vector3D(4.0971,3.82,3.3775);
//
//    Vector3D gold_n = Vector3D(0.21646,0.42833,1.3284);
//    Vector3D gold_k = Vector3D(3.2390,2.4599,1.8661);
//
////    Vector3D eta = n * Vector3D(614, 549, 466);
////    Vector3D etaK = k * Vector3D(614, 549, 466);
//    Vector3D eta = gold_n;
//    Vector3D etaK = gold_k;
    
    double cosine = dot(wi, Vector3D(0, 0, 1));
    
    Vector3D R_s = (eta*eta + k*k - 2 * eta * cosine + cosine * cosine) / (eta*eta + k*k + 2 * eta * cosine + cosine * cosine);
    Vector3D R_p = ((eta * eta + k*k)* cosine*cosine - 2 *eta*cosine + 1) / ((eta * eta + k*k)* cosine*cosine + 2 *eta*cosine + 1);
    
    return (R_s + R_p) / 2;
}

Vector3D MicrofacetBSDF::f(const Vector3D wo, const Vector3D wi) {
  // TODO Project 3-2: Part 2
  // Implement microfacet model here.
    Vector3D n = Vector3D(0,0,1);
    
    if (dot(n, wi) <= 0 || dot(n,wo) <= 0) {
        return 0;
    }
    
    Vector3D h = (wo + wi) / (wo + wi).norm();
    
    return (F(wi) * G(wo, wi) * D(h)) / (4 * dot(n,wo) * dot(n,wi));
//  return Vector3D();
}

Vector3D MicrofacetBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 2
  // *Importance* sample Beckmann normal distribution function (NDF) here.
  // Note: You should fill in the sampled direction *wi and the corresponding *pdf,
  //       and return the sampled BRDF value.

    alpha = 0.5;
    Vector2D r = sampler.get_sample();
    double r1 = r.x;
    double r2 = r.y;
    
    double theta_h = atan(sqrt(-pow(alpha,2) * log(1 - r1) ));
    double phi_h = 2 * PI * r2;
    
    double p_theta = (2 * sin(theta_h)) / (pow(alpha,2) * pow(cos(theta_h), 3)) * exp(-pow(tan(theta_h),2) / pow(alpha,2));
    double p_phi = 1 / (2 * PI);
    
    Vector3D h = Vector3D(sin(theta_h) * cos(phi_h), sin(theta_h) * sin(phi_h), cos(theta_h));
    h.normalize();
    
//    *wi = 2 * h * dot(wo, h) - wo;
//    *wi = wo - 2 * dot(wo, h) * h;
    
    if (dot(*wi, Vector3D(0,0,1)) <= 0) {
        *pdf = 0;
        return 0;
    }
    
    double p_w_h = p_theta * p_phi / sin(theta_h);
    *pdf = p_w_h / (4 * dot(*wi,h));
    
  *wi = cosineHemisphereSampler.get_sample(pdf);
  return MicrofacetBSDF::f(wo, *wi);
}

void MicrofacetBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Micofacet BSDF"))
  {
    DragDouble3("eta", &eta[0], 0.005);
    DragDouble3("K", &k[0], 0.005);
    DragDouble("alpha", &alpha, 0.005);
    ImGui::TreePop();
  }
}

// Refraction BSDF //

Vector3D RefractionBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D RefractionBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {
  // TODO Project 3-2: Part 1
  // Implement RefractionBSDF
    *pdf = 1;
    double n = 1.0 / this->ior;
    //exiting
    if (wo.z < 0) {
        n = this->ior;
    }
    if (!refract(wo, wi, this->ior)) {
        //*pdf = 1;
        //return transmittance / abs_cos_theta(*wi) / (n * n);
        return Vector3D();
    }
      //return Vector3D();
      return transmittance / abs_cos_theta(*wi) / (n * n);
}

void RefractionBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

// Glass BSDF //

Vector3D GlassBSDF::f(const Vector3D wo, const Vector3D wi) {
  return Vector3D();
}

Vector3D GlassBSDF::sample_f(const Vector3D wo, Vector3D* wi, double* pdf) {

  // TODO Project 3-2: Part 1
  // Compute Fresnel coefficient and either reflect or refract based on it.

  // compute Fresnel coefficient and use it as the probability of reflection
  // - Fundamentals of Computer Graphics page 305
    if (!refract(wo, wi, ior)) {
        reflect(wo, wi);
        *pdf = 1;
        return reflectance / abs_cos_theta(*wi);
    }
    else {
        double R0 = pow((1.0 - ior) / (1.0 + ior), 2);
        double R = R0 + (1 - R0) * pow((1.0 - abs_cos_theta(wo)), 5);
        if (coin_flip(R)) {
            reflect(wo, wi);
            *pdf = R;
            return R * reflectance / abs_cos_theta(*wi);
        }
        else {
            //refraction of wo to *wi
            double n = 1.0 / ior;
            int dir = -1;
            if (wo.z < 0) {
                n = ior;
                dir = 1;
            }
            double cos2 = 1 - (n * n) * (1 - (wo.z * wo.z));
            if (cos2 < 0) {

            }
            *wi = Vector3D(-n * wo.x, -n * wo.y, dir * sqrt(cos2));

            *pdf = 1 - R;
            return (1 - R) * transmittance / abs_cos_theta(*wi) / (n * n);
        }
    }

  return Vector3D();
}

void GlassBSDF::render_debugger_node()
{
  if (ImGui::TreeNode(this, "Refraction BSDF"))
  {
    DragDouble3("Reflectance", &reflectance[0], 0.005);
    DragDouble3("Transmittance", &transmittance[0], 0.005);
    DragDouble("ior", &ior, 0.005);
    ImGui::TreePop();
  }
}

void BSDF::reflect(const Vector3D wo, Vector3D* wi) {

  // TODO Project 3-2: Part 1
  // Implement reflection of wo about normal (0,0,1) and store result in wi.
    *wi = Vector3D(-wo.x, -wo.y, wo.z);
}

bool BSDF::refract(const Vector3D wo, Vector3D* wi, double ior) {

  // TODO Project 3-2: Part 1
  // Use Snell's Law to refract wo surface and store result ray in wi.
  // Return false if refraction does not occur due to total internal reflection
  // and true otherwise. When dot(wo,n) is positive, then wo corresponds to a
  // ray entering the surface through vacuum.

  //entering 
  double n = 1.0 / ior;
  int dir = -1;
  //exiting
  if (wo.z < 0) {
      n = ior;
      dir = 1;
  }
  double cos2 = 1 - (n * n) * (1 - (wo.z * wo.z));
  if (cos2 < 0) {
      return false;
  }
  *wi = Vector3D( -n * wo.x, -n * wo.y, dir * sqrt(cos2));
  return true;
}

} // namespace CGL
