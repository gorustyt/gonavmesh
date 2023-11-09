package recast

import "math"

func dtDistancePtSegSqr2D(pt, p, q []float64) (t float64, res float64) {
	pqx := q[0] - p[0]
	pqz := q[2] - p[2]
	dx := pt[0] - p[0]
	dz := pt[2] - p[2]
	d := pqx*pqx + pqz*pqz
	t = pqx*dx + pqz*dz
	if d > 0 {
		t /= d
	}
	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}
	dx = p[0] + t*pqx - pt[0]
	dz = p[2] + t*pqz - pt[2]
	return t, dx*dx + dz*dz
}

func dtCalcPolyCenter(idx []int, nidx int, verts []float64) (tc []float64) {
	tc[0] = 0.0
	tc[1] = 0.0
	tc[2] = 0.0
	for j := 0; j < nidx; j++ {
		v := rcGetVert(verts, idx[j])
		tc[0] += v[0]
		tc[1] += v[1]
		tc[2] += v[2]
	}
	s := 1.0 / float64(nidx)
	tc[0] *= s
	tc[1] *= s
	tc[2] *= s
	return tc
}

func dtClosestHeightPointTriangle(p, a, b, c []float64) (h float64, ok bool) {
	EPS := 1e-6

	v0 := dtVsub(c, a)
	v1 := dtVsub(b, a)
	v2 := dtVsub(p, a)

	// Compute scaled barycentric coordinates
	denom := v0[0]*v1[2] - v0[2]*v1[0]
	if math.Abs(denom) < EPS {
		return h, false
	}
	u := v1[2]*v2[0] - v1[0]*v2[2]
	v := v0[0]*v2[2] - v0[2]*v2[0]

	if denom < 0 {
		denom = -denom
		u = -u
		v = -v
	}

	// If point lies inside the triangle, return interpolated ycoord.
	if u >= 0.0 && v >= 0.0 && (u+v) <= denom {
		h = a[1] + (v0[1]*u+v1[1]*v)/denom
		return h, true
	}
	return h, false
}
