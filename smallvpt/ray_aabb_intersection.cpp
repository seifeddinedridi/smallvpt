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