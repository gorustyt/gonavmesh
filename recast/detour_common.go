package recast

import (
	"gonavamesh/common"
	"math"
)

func DtDistancePtSegSqr2D(pt, p, q []float64) (t float64, res float64) {
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
	v0 := make([]float64, 3)
	v1 := make([]float64, 3)
	v2 := make([]float64, 3)
	common.Vsub(v0, c, a)
	common.Vsub(v1, b, a)
	common.Vsub(v2, p, a)

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

func dtOppositeTile(side int) int { return (side + 4) & 0x7 }

// / Determines if two axis-aligned bounding boxes overlap.
// /  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
// /  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
// /  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
// /  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
// / @return True if the two AABB's overlap.
// / @see dtOverlapBounds
func dtOverlapQuantBounds(amin [3]int, amax [3]int, bmin [3]int, bmax [3]int) bool {
	overlap := true
	if amin[0] > bmax[0] || amax[0] < bmin[0] {
		overlap = false
	}
	if amin[1] > bmax[1] || amax[1] < bmin[1] {
		overlap = false
	}
	if amin[2] > bmax[2] || amax[2] < bmin[2] {
		overlap = false
	}

	return overlap
}

// / Determines if two axis-aligned bounding boxes overlap.
// /  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
// /  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
// /  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
// /  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
// / @return True if the two AABB's overlap.
// / @see dtOverlapQuantBounds
func dtOverlapBounds(amin, amax, bmin, bmax []float64) bool {
	overlap := true
	if amin[0] > bmax[0] || amax[0] < bmin[0] {
		overlap = false
	}
	if amin[1] > bmax[1] || amax[1] < bmin[1] {
		overlap = false
	}
	if amin[2] > bmax[2] || amax[2] < bmin[2] {
		overlap = false
	}
	return overlap
}

// 在凸多边形里面随机一个点
// Returns a random point in a convex polygon.
// Adapted from Graphics Gems article.
func dtRandomPointInConvexPoly(pts []float64, npts int, areas []float64, s, t float64) (pt []float64) {
	pt = make([]float64, 3)
	// Calc triangle araes
	areasum := 0.0
	for i := 2; i < npts; i++ {
		areas[i] = common.TriArea2D(rcGetVert(pts, 0), rcGetVert(pts, i-1), rcGetVert(pts, i))
		areasum += common.Max(0.001, areas[i])
	}
	// Find sub triangle weighted by area.
	thr := s * areasum
	acc := 0.0
	u := 1.0
	tri := npts - 1
	for i := 2; i < npts; i++ {
		dacc := areas[i]
		if thr >= acc && thr < (acc+dacc) {
			u = (thr - acc) / dacc
			tri = i
			break
		}
		acc += dacc
	}

	v := math.Sqrt(t)

	a := 1 - v
	b := (1 - u) * v
	c := u * v
	pa := rcGetVert(pts, 0)
	pb := rcGetVert(pts, tri-1)
	pc := rcGetVert(pts, tri)

	pt[0] = a*pa[0] + b*pb[0] + c*pc[0]
	pt[1] = a*pa[1] + b*pb[1] + c*pc[1]
	pt[2] = a*pa[2] + b*pb[2] + c*pc[2]
	return pt
}

func dtDistancePtPolyEdgesSqr(pt, verts []float64, nverts int, ed, et []float64) (c bool) {
	// TODO: Replace pnpoly with triArea2D tests?
	i := 0
	j := nverts - 1
	for i < nverts {
		vi := rcGetVert(verts, i)
		vj := rcGetVert(verts, j)
		if ((vi[2] > pt[2]) != (vj[2] > pt[2])) && (pt[0] < (vj[0]-vi[0])*(pt[2]-vi[2])/(vj[2]-vi[2])+vi[0]) {
			c = !c
		}
		_, ed[j] = DtDistancePtSegSqr2D(pt, vj, vi)
		j = i
		i++
	}
	return c
}
func vperpXZ(a, b []float64) float64 { return a[0]*b[2] - a[2]*b[0] }
func dtIntersectSegSeg2D(ap, aq, bp, bq []float64) (s, t float64, ok bool) {
	u := make([]float64, 3)
	v := make([]float64, 3)
	w := make([]float64, 3)
	common.Vsub(u, aq, ap)
	common.Vsub(v, bq, bp)
	common.Vsub(w, ap, bp)
	d := vperpXZ(u, v)
	if math.Abs(d) < 1e-6 {
		return s, t, false
	}
	s = vperpXZ(v, w) / d
	t = vperpXZ(u, w) / d
	return s, t, true
}
func dtAssertTrue(ok bool) {
	if !ok {
		panic("")
	}
}

func dtIntersectSegmentPoly2D(p0, p1, verts []float64, nverts int) (tmin, tmax float64, segMin, segMax int, ok bool) {
	EPS := 0.000001

	tmin = 0
	tmax = 1
	segMin = -1
	segMax = -1
	dir := make([]float64, 3)
	common.Vsub(dir, p1, p0)

	i := 0
	j := nverts - 1
	for i < nverts {
		edge := make([]float64, 3)
		diff := make([]float64, 3)
		common.Vsub(edge, rcGetVert(verts, i), rcGetVert(verts, j))
		common.Vsub(diff, p0, rcGetVert(verts, j))
		n := common.Vperp2D(edge, diff)
		d := common.Vperp2D(dir, edge)
		if math.Abs(d) < EPS {
			// S is nearly parallel to this edge
			if n < 0 {
				return tmin, tmax, segMin, segMax, false
			} else {
				j = i
				i++
				continue
			}

		}
		t := n / d
		if d < 0 {
			// segment S is entering across this edge
			if t > tmin {
				tmin = t
				segMin = j
				// S enters after leaving polygon
				if tmin > tmax {
					return tmin, tmax, segMin, segMax, false
				}

			}
		} else {
			// segment S is leaving across this edge
			if t < tmax {
				tmax = t
				segMax = j
				// S leaves before entering polygon
				if tmax < tmin {
					return tmin, tmax, segMin, segMax, false
				}

			}
		}
		j = i
		i++
	}

	return tmin, tmax, segMin, segMax, true
}

// / @par
// /
// / All vertices are projected onto the xz-plane, so the y-values are ignored.
func dtOverlapPolyPoly2D(polya []float64, npolya int,
	polyb []float64, npolyb int) bool {
	eps := 1e-4
	i := 0
	j := npolya - 1
	for i < npolya {
		va := polya[j*3 : j*3+3]
		vb := polya[i*3 : i*3+3]
		n := [3]float64{vb[2] - va[2], 0, -(vb[0] - va[0])}
		amin, amax := projectPoly(n[:], polya, npolya)
		bmin, bmax := projectPoly(n[:], polyb, npolyb)
		if !overlapRange(amin, amax, bmin, bmax, eps) {
			// Found separating axis
			return false
		}
		j = i
		i++
	}
	i = 0
	j = npolyb - 1
	for i < npolyb {
		va := polyb[j*3 : j*3+3]
		vb := polyb[i*3 : i*3+3]
		var n = [3]float64{vb[2] - va[2], 0, -(vb[0] - va[0])}
		amin, amax := projectPoly(n[:], polya, npolya)
		bmin, bmax := projectPoly(n[:], polyb, npolyb)
		if !overlapRange(amin, amax, bmin, bmax, eps) {
			// Found separating axis
			return false
		}
		j = i
		i++
	}
	return true
}

func projectPoly(axis, poly []float64, npoly int) (rmin, rmax float64) {
	rmax = common.Vdot2D(axis, poly[0:3])
	rmin = rmax
	for i := 1; i < npoly; i++ {
		d := common.Vdot2D(axis, poly[i*3:i*3+3])
		rmin = common.Min(rmin, d)
		rmax = common.Max(rmax, d)
	}
	return rmin, rmax
}

func overlapRange(amin, amax, bmin, bmax, eps float64) bool {
	if (amin+eps) > bmax || (amax-eps) < bmin {
		return false
	}
	return true
}
