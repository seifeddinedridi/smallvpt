#define _USE_MATH_DEFINES
#include <math.h>   // smallpt, a Path Tracer by Kevin Beason, 2008
#include <stdlib.h> // Make : g++ -O3 -fopenmp smallpt.cpp -o smallpt
#include <stdio.h>  //        Remove "-fopenmp" for g++ version < 4.2
#include <algorithm>
#pragma warning(disable: 4244) // Disable double to float warning
namespace XORShift {
	// XOR shift PRNG
	unsigned int x = 123456789;
	unsigned int y = 362436069;
	unsigned int z = 521288629;
	unsigned int w = 88675123; 

	inline float frand()
	{ 
		unsigned int t;

		t = x ^ (x << 11);
		x = y; y = z; z = w;
		return (w = (w ^ (w >> 19)) ^ (t ^ (t >> 8))) * (1.0f / 4294967295.0f); 
	}
}
struct Vec {        // Usage: time ./smallpt 5000 && xv image.ppm
	double x, y, z;                  // position, also color (r,g,b)
	Vec(double x_=0, double y_=0, double z_=0){ x=x_; y=y_; z=z_; }
	Vec operator+(const Vec &b) const { return Vec(x+b.x,y+b.y,z+b.z); }
	Vec operator-(const Vec &b) const { return Vec(x-b.x,y-b.y,z-b.z); }
	Vec operator*(double b) const { return Vec(x*b,y*b,z*b); }
	Vec mult(const Vec &b) const { return Vec(x*b.x,y*b.y,z*b.z); }
	Vec& norm(){ return *this = *this * (1/sqrt(x*x+y*y+z*z)); }
	float length() {return sqrt(x*x+y*y+z*z); }
	double dot(const Vec &b) const { return x*b.x+y*b.y+z*b.z; } // cross:
	Vec operator%(Vec&b){return Vec(y*b.z-z*b.y,z*b.x-x*b.z,x*b.y-y*b.x);}
};
struct Ray { Vec o, d; Ray() {} Ray(Vec o_, Vec d_) : o(o_), d(d_) {} };
enum Refl_t { DIFF, SPEC, REFR };  // material types, used in radiance()
bool intersectRayAABB(const Ray &r, const Vec &pmin, const Vec &pmax, float *tnear, float *tfar) {
	float tmin, tmax, tymin, tymax, tzmin, tzmax;
	Vec inv_direction = Vec(1/r.d.x, 1/r.d.y, 1/r.d.z);
	int sign[3];
	sign[0] = (inv_direction.x < 0);
	sign[1] = (inv_direction.y < 0);
	sign[2] = (inv_direction.z < 0);
	Vec parameters[2] = {pmin, pmax};
	tmin = (parameters[sign[0]].x - r.o.x) * inv_direction.x;
	tmax = (parameters[1-sign[0]].x - r.o.x) * inv_direction.x;
	tymin = (parameters[sign[1]].y - r.o.y) * inv_direction.y;
	tymax = (parameters[1-sign[1]].y - r.o.y) * inv_direction.y;
	if ( (tmin > tymax) || (tymin > tmax) ) 
		return false;
	if (tymin > tmin)
		tmin = tymin;
	if (tymax < tmax)
		tmax = tymax;
	tzmin = (parameters[sign[2]].z - r.o.z) * inv_direction.z;
	tzmax = (parameters[1-sign[2]].z - r.o.z) * inv_direction.z;
	if ( (tmin > tzmax) || (tzmin > tmax) ) 
		return false;
	if (tzmin > tmin)
		tmin = tzmin;
	if (tzmax < tmax)
		tmax = tzmax;
	*tnear = std::max(tmin, std::max(tymin,tzmin));
	*tfar = std::min(tmax, std::min(tymax, tzmax));
	return true;
}
struct Sphere {
	double rad;       // radius
	Vec p, e, c;      // position, emission, color
	Refl_t refl;      // reflection type (DIFFuse, SPECular, REFRactive)
	Sphere(double rad_, Vec p_, Vec e_, Vec c_, Refl_t refl_):
	rad(rad_), p(p_), e(e_), c(c_), refl(refl_) {}
	double intersect(const Ray &r) const { // returns distance, 0 if nohit
		Vec op = p-r.o; // Solve t^2*d.d + 2*t*(o-p).d + (o-p).(o-p)-R^2 = 0
		double t, eps=1e-4, b=op.dot(r.d), det=b*b-op.dot(op)+rad*rad;
		if (det<0) return 0; else det=sqrt(det);
		return (t=b-det)>eps ? t : ((t=b+det)>eps ? t : 0);
	}
};
Sphere spheres[] = {//Scene: radius, position, emission, color, material 
	Sphere(1e5, Vec( 1e5+1,40.8,81.6), Vec(),Vec(.75,.25,.25),DIFF),//Left
	Sphere(1e5, Vec(-1e5+99,40.8,81.6),Vec(),Vec(.25,.25,.75),DIFF),//Right
	Sphere(1e5, Vec(50,40.8, 1e5),     Vec(),Vec(.75,.75,.75),DIFF),//Back
	Sphere(1e5, Vec(50,40.8,-1e5+170), Vec(),Vec(),           DIFF),//Frnt
	Sphere(1e5, Vec(50, 1e5, 81.6),    Vec(),Vec(.75,.75,.75),DIFF),//Botm
	Sphere(1e5, Vec(50,-1e5+81.6,81.6),Vec(),Vec(.75,.75,.75),DIFF),//Top
	Sphere(16.5,Vec(27,16.5,47),       Vec(),Vec(1,1,1)*.75, DIFF),//Mirr
	Sphere(16.5,Vec(73,16.5,78),       Vec(),Vec(1,1,1)*.75, DIFF),//Glas
	Sphere(600, Vec(50,681.6-.17,81.6),Vec(12,12,12),  Vec(), DIFF) //Lite
};
const int lightId = 8;
struct HomogeneousMedium {
	HomogeneousMedium(double sigs, double siga): sigma_s(sigs), sigma_a(siga), sigma_t(sigma_s + sigma_a) {}
	double sigma_s, sigma_a, sigma_t;
} medium(0.04, 0.001);
inline double clamp(double x){ return x<0 ? 0 : x>1 ? 1 : x; }
inline int toInt(double x){ return int(pow(clamp(x),1/2.2)*255+.5); }
inline bool intersect(const Ray &r, double &t, int &id){
	double n=sizeof(spheres)/sizeof(Sphere), d, inf=t=1e20;
	for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){t=d;id=i;}
	return t<inf;
}
inline bool intersectP(const Ray &r, double t) {
	double n=sizeof(spheres)/sizeof(Sphere), d;
	for(int i=int(n);i--;) if((d=spheres[i].intersect(r))&&d<t){return true;}
	return false;
}
inline double sampleSegment(double epsilon, float sigma, float smax) {
	return -log(1.0 - epsilon * (1.0 - exp(-sigma * smax))) / sigma;
}
inline Vec sampleSphere(double e1, double e2) {
	double z = 1.0 - 2.0 * e1, xx = sqrt(1.0 - z * z);
	return Vec(cos(2.0 * M_PI * e2) * xx, sin(2.0 * M_PI * e2) * xx, z);
}
inline Vec singleScatter(const Ray &r, double tmax) {
	// Sample a point along the ray's extents
	double s = sampleSegment(XORShift::frand(), medium.sigma_s, tmax);
	Vec x = r.o + r.d * s;
	Vec dir = sampleSphere(XORShift::frand(), XORShift::frand()); // Sample a direction ~ uniform phase function
	double tLight;
	int id = 0;
	if (intersect(Ray(x, dir), tLight, id) && id == lightId) // make sure that the closest intersection is a light source
		return spheres[lightId].e * exp(-medium.sigma_a * tLight) * (1.0 - exp(-medium.sigma_s * tmax));
	return Vec();
}
inline float multipleScatter(const Ray &r, Ray *sRay, double tin, float tout) {
	// Sample a point along the ray's extents
	double s = sampleSegment(XORShift::frand(), medium.sigma_s, tout-tin);
	Vec x = r.o + r.d *(tin < 0 ? 0 : tin) + r.d * s;
	Vec dir = sampleSphere(XORShift::frand(), XORShift::frand()); // Sample a direction ~ uniform phase function
	if (sRay)	*sRay = Ray(x, dir);
	return (1.0 - exp(-medium.sigma_s * (tout-tin)));
}/*
Vec radiance(const Ray &r, int depth) {
	if (++depth == 20)
		return Vec();
	float tnear, tfar;
	bool intrsctmd = intersectRayAABB(r, Vec(30, 40, 20), Vec(70, 80, 100), &tnear, &tfar);
	double t;                               // distance to intersection
	int id=0;                               // id of intersected object
	if (!intersect(r, t, id)) {
		if (!intrsctmd) return Vec();
		Ray sRay;
		float ms = multipleScatter(r, &sRay, tnear, tfar);
		return radiance(sRay, depth) * ms;
	}
	const Sphere &obj = spheres[id];        // the hit object
	if (intrsctmd)
		return obj.e  * expf(-medium.sigma_a * (t-tnear));
	else
		return obj.e;
}*/
Vec radiance(const Ray &r, int depth) {
	double t;                               // distance to intersection
	int id=0;                               // id of intersected object
	float tnear, tfar;
	bool intrsctmd = intersectRayAABB(r, Vec(30, 30, 0), Vec(80, 100, 100), &tnear, &tfar);
	if (!intersect(r, t, id)) {
		if (++depth >= 20 || !intrsctmd) return Vec();
		Ray sRay;
		return radiance(sRay, depth) * multipleScatter(r, &sRay, tnear, tfar);
	}
	const Sphere &obj = spheres[id];        // the hit object
	Vec x=r.o+r.d*t, n=(x-obj.p).norm(), nl=n.dot(r.d)<0?n:n*-1, f=obj.c;
	double p = f.x>f.y && f.x>f.z ? f.x : f.y>f.z ? f.y : f.z; // max refl
	Ray sRay;
	float ms = 0;
	if (intrsctmd)
		ms = multipleScatter(r, &sRay, tnear, tfar);
	if (++depth>5) if (XORShift::frand()<p) {f=f*(1/p);ms = ms *(1/p);} else return Vec(); //R.R.
	float scaleBy = 1.0f;
	Vec Le = obj.e;
	if (intrsctmd)
	{
		scaleBy = 2.0f;
		float dist = (t > tfar ? tfar - tnear : t - tnear);
		f = f * expf(-medium.sigma_a * dist); // Absorption
		Le = obj.e * expf(-medium.sigma_a * dist);
		// Sample surface or volume?
		if (XORShift::frand() <= 0.5f && (n.dot(nl)>0)) // no scattering inside a dielectric
		{
			return radiance(sRay, depth) * ms * 2.0f;
		}
	}
	if (obj.refl == DIFF) {                  // Ideal DIFFUSE reflection
		double r1=2*M_PI*XORShift::frand(), r2=XORShift::frand(), r2s=sqrt(r2);
		Vec w=nl, u=((fabs(w.x)>.1?Vec(0,1):Vec(1))%w).norm(), v=w%u;
		Vec d = (u*cos(r1)*r2s + v*sin(r1)*r2s + w*sqrt(1-r2)).norm();
		return (Le + f.mult(radiance(Ray(x,d),depth))) * scaleBy;
	} else if (obj.refl == SPEC)            // Ideal SPECULAR reflection
		return (Le + f.mult(radiance(Ray(x,r.d-n*2*n.dot(r.d)),depth))) * scaleBy;
	Ray reflRay(x, r.d-n*2*n.dot(r.d));     // Ideal dielectric REFRACTION
	bool into = n.dot(nl)>0;                // Ray from outside going in?
	double nc=1, nt=1.5, nnt=into?nc/nt:nt/nc, ddn=r.d.dot(nl), cos2t;
	if ((cos2t=1-nnt*nnt*(1-ddn*ddn))<0)    // Total internal reflection
		return (Le + f.mult(radiance(reflRay,depth))) * scaleBy;
	Vec tdir = (r.d*nnt - n*((into?1:-1)*(ddn*nnt+sqrt(cos2t)))).norm();
	double a=nt-nc, b=nt+nc, R0=a*a/(b*b), c = 1-(into?-ddn:tdir.dot(n));
	double Re=R0+(1-R0)*c*c*c*c*c,Tr=1-Re,P=.25+.5*Re,RP=Re/P,TP=Tr/(1-P);
	return (Le + f.mult(depth>2 ? (XORShift::frand()<P ?   // Russian roulette
		radiance(reflRay,depth)*RP:radiance(Ray(x,tdir),depth)*TP) :
	radiance(reflRay,depth)*Re+radiance(Ray(x,tdir),depth)*Tr)) * scaleBy;
}
int main(int argc, char *argv[]) {
	int w=400, h=400, samps = argc==2 ? atoi(argv[1])/4 : 1; // # samples
	Ray cam(Vec(50,52,345.6), Vec(0,-0.042612,-1).norm()); // cam pos, dir
	Vec cx=Vec(w*.5135/h), cy=(cx%cam.d).norm()*.5135, r, *c=new Vec[w*h];
#pragma omp parallel for schedule(dynamic, 1) private(r)       // OpenMP
	for (int y=0; y<h; y++){                       // Loop over image rows
		fprintf(stderr,"\rRendering (%d spp) %5.2f%%",samps*4,100.*y/(h-1));
		for (unsigned short x=0; x<w; x++)   // Loop cols
			for (int sy=0, i=(h-y-1)*w+x; sy<2; sy++)     // 2x2 subpixel rows
				for (int sx=0; sx<2; sx++, r=Vec()){        // 2x2 subpixel cols
					for (int s=0; s<samps; s++){
						double r1=2*XORShift::frand(), dx=r1<1 ? sqrt(r1)-1: 1-sqrt(2-r1);
						double r2=2*XORShift::frand(), dy=r2<1 ? sqrt(r2)-1: 1-sqrt(2-r2);
						Vec d = cx*( ( (sx+.5 + dx)/2 + x)/w - .5) +
							cy*( ( (sy+.5 + dy)/2 + y)/h - .5) + cam.d;
						r = r + radiance(Ray(cam.o+d*140,d.norm()),0)*(1./samps);
					} // Camera rays are pushed ^^^^^ forward to start in interior
					c[i] = c[i] + Vec(clamp(r.x),clamp(r.y),clamp(r.z))*.25;
				}
	}
	FILE *f = fopen("..\\Renders\\image.ppm", "w"); // Write image to PPM file.
	fprintf(f, "P3\n%d %d\n%d\n", w, h, 255);
	for (int i=0; i<w*h; i++)
		fprintf(f,"%d %d %d ", toInt(c[i].x), toInt(c[i].y), toInt(c[i].z));
}
