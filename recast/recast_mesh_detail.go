package recast

import (
	"gonavamesh/common"
	"math"
)

const (
	RC_UNSET_HEIGHT = 0xffff
)

// / Polygon touches multiple regions.
// / If a polygon has this region ID it was merged with or created
// / from polygons of different regions during the polymesh
// / build step that removes redundant border vertices.
// / (Used during the polymesh and detail polymesh build processes)
// / @see RcPolyMesh::regs
const RC_MULTIPLE_REGS = 0

type rcHeightPatch struct {
	data                      []int
	xmin, ymin, width, height int
}

func vdot2(a, b []float64) float64 {
	return a[0]*b[0] + a[2]*b[2]
}

func vdistSq2(p, q []float64) float64 {
	dx := q[0] - p[0]
	dy := q[2] - p[2]
	return dx*dx + dy*dy
}

func vdist2(p, q []float64) float64 {
	return math.Sqrt(vdistSq2(p, q))
}

func circumCircle(p1, p2, p3, c []float64, r *float64) bool {
	EPS := 1e-6
	// Calculate the circle relative to p1, to avoid some precision issues.
	v1 := []float64{0, 0, 0}
	v2 := make([]float64, 3)
	v3 := make([]float64, 3)
	common.Vsub(v2, p2, p1)
	common.Vsub(v3, p3, p1)

	cp := vcross2(v1, v2, v3)
	if math.Abs(cp) > EPS {
		v1Sq := vdot2(v1, v1)
		v2Sq := vdot2(v2, v2)
		v3Sq := vdot2(v3, v3)
		c[0] = (v1Sq*(v2[2]-v3[2]) + v2Sq*(v3[2]-v1[2]) + v3Sq*(v1[2]-v2[2])) / (2 * cp)
		c[1] = 0
		c[2] = (v1Sq*(v3[0]-v2[0]) + v2Sq*(v1[0]-v3[0]) + v3Sq*(v2[0]-v1[0])) / (2 * cp)
		*r = vdist2(c, v1)
		common.Vadd(c, c, p1)
		return true
	}

	copy(c, p1)
	*r = 0
	return false
}

func vcross2(p1, p2, p3 []float64) float64 {
	u1 := p2[0] - p1[0]
	v1 := p2[2] - p1[2]
	u2 := p3[0] - p1[0]
	v2 := p3[2] - p1[2]
	return u1*v2 - v1*u2
}

func distPtTri(p, a, b, c []float64) float64 {
	v0 := make([]float64, 3)
	v1 := make([]float64, 3)
	v2 := make([]float64, 3)
	common.Vsub(v0, c, a)
	common.Vsub(v1, b, a)
	common.Vsub(v2, p, a)

	dot00 := vdot2(v0, v0)
	dot01 := vdot2(v0, v1)
	dot02 := vdot2(v0, v2)
	dot11 := vdot2(v1, v1)
	dot12 := vdot2(v1, v2)

	// Compute barycentric coordinates
	invDenom := 1.0 / (dot00*dot11 - dot01*dot01)
	u := (dot11*dot02 - dot01*dot12) * invDenom
	v := (dot00*dot12 - dot01*dot02) * invDenom

	// If point lies inside the triangle, return interpolated y-coord.
	EPS := 1e-4
	if u >= -EPS && v >= -EPS && (u+v) <= 1+EPS {
		y := a[1] + v0[1]*u + v1[1]*v
		return math.Abs(y - p[1])
	}
	return math.MaxFloat64
}

func distancePtSeg(pt, p, q []float64) float64 {
	pqx := q[0] - p[0]
	pqy := q[1] - p[1]
	pqz := q[2] - p[2]
	dx := pt[0] - p[0]
	dy := pt[1] - p[1]
	dz := pt[2] - p[2]
	d := pqx*pqx + pqy*pqy + pqz*pqz
	t := pqx*dx + pqy*dy + pqz*dz
	if d > 0 {
		t /= d
	}

	if t < 0 {
		t = 0
	} else if t > 1 {
		t = 1
	}

	dx = p[0] + t*pqx - pt[0]
	dy = p[1] + t*pqy - pt[1]
	dz = p[2] + t*pqz - pt[2]

	return dx*dx + dy*dy + dz*dz
}

func distancePtSeg2d(pt, p, q []float64) float64 {
	pqx := q[0] - p[0]
	pqz := q[2] - p[2]
	dx := pt[0] - p[0]
	dz := pt[2] - p[2]
	d := pqx*pqx + pqz*pqz
	t := pqx*dx + pqz*dz
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

	return dx*dx + dz*dz
}

func distToTriMesh(p, verts []float64, nverts int, tris []int, ntris int) float64 {
	dmin := math.MaxFloat64
	for i := 0; i < ntris; i++ {
		va := rcGetVert(verts, tris[i*4+0])
		vb := rcGetVert(verts, tris[i*4+1])
		vc := rcGetVert(verts, tris[i*4+2])
		d := distPtTri(p, va, vb, vc)
		if d < dmin {
			dmin = d
		}
	}
	if dmin == math.MaxFloat64 {
		return -1
	}
	return dmin
}

func distToPoly(nvert int, verts []float64, p []float64) float64 {

	dmin := math.MaxFloat64
	var i, j int
	i = 0
	j = nvert - 1
	c := false
	for i < nvert {
		vi := rcGetVert(verts, i)
		vj := rcGetVert(verts, j)
		if ((vi[2] > p[2]) != (vj[2] > p[2])) &&
			(p[0] < (vj[0]-vi[0])*(p[2]-vi[2])/(vj[2]-vi[2])+vi[0]) {
			c = !c
		}

		dmin = common.Min(dmin, distancePtSeg2d(p, vj, vi))
		j = i
		i++
	}
	if c {
		return -dmin
	}
	return dmin
}
func getHeight(fx, fy, fz,
	cs, ics, ch float64,
	radius int, hp *rcHeightPatch) int {
	ix := int(math.Floor(fx*ics + 0.01))
	iz := int(math.Floor(fz*ics + 0.01))
	ix = common.Clamp(ix-hp.xmin, 0, hp.width-1)
	iz = common.Clamp(iz-hp.ymin, 0, hp.height-1)
	h := hp.data[ix+iz*hp.width]
	if h == RC_UNSET_HEIGHT {
		// Special case when data might be bad.
		// Walk adjacent cells in a spiral up to 'radius', and look
		// for a pixel which has a valid height.
		x := 1
		z := 0
		dx := 1
		dz := 0
		maxSize := radius*2 + 1
		maxIter := maxSize*maxSize - 1

		nextRingIterStart := 8
		nextRingIters := 16

		dmin := math.MaxFloat64
		for i := 0; i < maxIter; i++ {
			nx := ix + x
			nz := iz + z

			if nx >= 0 && nz >= 0 && nx < hp.width && nz < hp.height {
				nh := hp.data[nx+nz*hp.width]
				if nh != RC_UNSET_HEIGHT {
					d := math.Abs(float64(nh)*ch - fy)
					if d < dmin {
						h = nh
						dmin = d
					}
				}
			}

			// We are searching in a grid which looks approximately like this:
			//  __________
			// |2 ______ 2|
			// | |1 __ 1| |
			// | | |__| | |
			// | |______| |
			// |__________|
			// We want to find the best height as close to the center cell as possible. This means that
			// if we find a height in one of the neighbor cells to the center, we don't want to
			// expand further out than the 8 neighbors - we want to limit our search to the closest
			// of these "rings", but the best height in the ring.
			// For example, the center is just 1 cell. We checked that at the entrance to the function.
			// The next "ring" contains 8 cells (marked 1 above). Those are all the neighbors to the center cell.
			// The next one again contains 16 cells (marked 2). In general each ring has 8 additional cells, which
			// can be thought of as adding 2 cells around the "center" of each side when we expand the ring.
			// Here we detect if we are about to enter the next ring, and if we are and we have found
			// a height, we abort the search.
			if i+1 == nextRingIterStart {
				if h != RC_UNSET_HEIGHT {
					break
				}

				nextRingIterStart += nextRingIters
				nextRingIters += 8
			}

			if (x == z) || ((x < 0) && (x == -z)) || ((x > 0) && (x == 1-z)) {
				tmp := dx
				dx = -dz
				dz = tmp
			}
			x += dx
			z += dz
		}
	}
	return h
}

const (
	EV_UNDEF = -1
	EV_HULL  = -2
)

func findEdge(edges []int, nedges int, s, t int) int {
	for i := 0; i < nedges; i++ {
		e := rcGetVert4(edges, i)
		if (e[0] == s && e[1] == t) || (e[0] == t && e[1] == s) {
			return i
		}

	}
	return EV_UNDEF
}

func addEdge(edges []int, nedges *int, maxEdges int, s, t, l, r int) int {
	if *nedges >= maxEdges {
		return EV_UNDEF
	}

	// Add edge if not already in the triangulation.
	e := findEdge(edges, *nedges, s, t)
	if e == EV_UNDEF {
		edge := rcGetVert4(edges, *nedges)
		edge[0] = s
		edge[1] = t
		edge[2] = l
		edge[3] = r
		old := *nedges
		*nedges++
		return old
	} else {
		return EV_UNDEF
	}
}

func updateLeftFace(e []int, s, t, f int) {
	if e[0] == s && e[1] == t && e[2] == EV_UNDEF {
		e[2] = f
	} else if e[1] == s && e[0] == t && e[3] == EV_UNDEF {
		e[3] = f
	}

}

func overlapSegSeg2d(a, b, c, d []float64) int {
	a1 := vcross2(a, b, d)
	a2 := vcross2(a, b, c)
	if a1*a2 < 0.0 {
		a3 := vcross2(c, d, a)
		a4 := a3 + a2 - a1
		if a3*a4 < 0.0 {
			return 1
		}

	}
	return 0
}

func overlapEdges(pts []float64, edges []int, nedges, s1, t1 int) bool {
	for i := 0; i < nedges; i++ {
		s0 := edges[i*4+0]
		t0 := edges[i*4+1]
		// Same or connected edges do not overlap.
		if s0 == s1 || s0 == t1 || t0 == s1 || t0 == t1 {
			continue
		}

		if overlapSegSeg2d(rcGetVert(pts, s0), rcGetVert(pts, t0), rcGetVert(pts, s1), rcGetVert(pts, t1)) != 0 {
			return true
		}

	}
	return false
}

func completeFacet(pts []float64, npts int, edges []int, nedges *int, maxEdges int, nfaces *int, e int) {
	EPS := 1e-5

	edge := rcGetVert4(edges, e)

	// Cache s and t.
	var s, t int
	if edge[2] == EV_UNDEF {
		s = edge[0]
		t = edge[1]
	} else if edge[3] == EV_UNDEF {
		s = edge[1]
		t = edge[0]
	} else {
		// Edge already completed.
		return
	}

	// Find best point on left of edge.
	pt := npts
	c := []float64{0, 0, 0}
	r := float64(-1)
	for u := 0; u < npts; u++ {
		if u == s || u == t {
			continue
		}
		if vcross2(rcGetVert(pts, s), rcGetVert(pts, t), rcGetVert(pts, u)) > EPS {
			if r < 0 {
				// The circle is not updated yet, do it now.
				pt = u
				circumCircle(rcGetVert(pts, s), rcGetVert(pts, t), rcGetVert(pts, u), c, &r)
				continue
			}
			d := vdist2(c, rcGetVert(pts, u))
			tol := 0.001
			if d > r*(1+tol) {
				// Outside current circumcircle, skip.
				continue
			} else if d < r*(1-tol) {
				// Inside safe circumcircle, update circle.
				pt = u
				circumCircle(rcGetVert(pts, s), rcGetVert(pts, t), rcGetVert(pts, u), c, &r)
			} else {
				// Inside epsilon circum circle, do extra tests to make sure the edge is valid.
				// s-u and t-u cannot overlap with s-pt nor t-pt if they exists.
				if overlapEdges(pts, edges, *nedges, s, u) {
					continue
				}

				if overlapEdges(pts, edges, *nedges, t, u) {
					continue
				}

				// Edge is valid.
				pt = u
				circumCircle(rcGetVert(pts, s), rcGetVert(pts, t), rcGetVert(pts, u), c, &r)
			}
		}
	}

	// Add new triangle or update edge info if s-t is on hull.
	if pt < npts {
		// Update face information of edge being completed.
		updateLeftFace(rcGetVert4(edges, e), s, t, *nfaces)

		// Add new edge or update face info of old edge.
		e = findEdge(edges, *nedges, pt, s)
		if e == EV_UNDEF {
			addEdge(edges, nedges, maxEdges, pt, s, *nfaces, EV_UNDEF)
		} else {
			updateLeftFace(rcGetVert4(edges, e), pt, s, *nfaces)
		}

		// Add new edge or update face info of old edge.
		e = findEdge(edges, *nedges, t, pt)
		if e == EV_UNDEF {
			addEdge(edges, nedges, maxEdges, t, pt, *nfaces, EV_UNDEF)
		} else {
			updateLeftFace(rcGetVert4(edges, e), t, pt, *nfaces)
		}
		*nfaces++
	} else {
		updateLeftFace(rcGetVert4(edges, e), s, t, EV_HULL)
	}
}

func delaunayHull(npts int, pts []float64,
	nhull int, hull []int,
	tris Stack[int], edges Stack[int]) {
	nfaces := 0
	nedges := 0
	maxEdges := npts * 10
	edges.Resize(maxEdges * 4)

	i := 0
	j := nhull - 1
	for i < nhull {
		addEdge(edges.Data(), &nedges, maxEdges, hull[j], hull[i], EV_HULL, EV_UNDEF)
		j = i
		i++
	}

	currentEdge := 0
	for currentEdge < nedges {
		if edges.Index(currentEdge*4+2) == EV_UNDEF {
			completeFacet(pts, npts, edges.Data(), &nedges, maxEdges, &nfaces, currentEdge)
		}

		if edges.Index(currentEdge*4+3) == EV_UNDEF {
			completeFacet(pts, npts, edges.Data(), &nedges, maxEdges, &nfaces, currentEdge)
		}

		currentEdge++
	}

	// Create tris
	tris.Resize(nfaces * 4)
	for i := 0; i < nfaces*4; i++ {
		tris.SetByIndex(i, -1)
	}

	for i := 0; i < nedges; i++ {
		e := rcGetVert4(edges.Data(), i)
		if e[3] >= 0 {
			// Left face
			t := rcGetVert4(tris.Slice(0, tris.Len()), e[3])
			if t[0] == -1 {
				t[0] = e[0]
				t[1] = e[1]
			} else if t[0] == e[1] {
				t[2] = e[0]
			} else if t[1] == e[0] {
				t[2] = e[1]
			}

		}
		if e[2] >= 0 {
			// Right
			t := rcGetVert4(tris.Data(), e[2])
			if t[0] == -1 {
				t[0] = e[1]
				t[1] = e[0]
			} else if t[0] == e[0] {
				t[2] = e[1]
			} else if t[1] == e[1] {
				t[2] = e[0]
			}

		}
	}

	for i := 0; i < tris.Len()/4; i++ {
		t := rcGetVert4(tris.Data(), i)
		if t[0] == -1 || t[1] == -1 || t[2] == -1 {
			t[0] = tris.Index(tris.Len() - 4)
			t[1] = tris.Index(tris.Len() - 3)
			t[2] = tris.Index(tris.Len() - 2)
			t[3] = tris.Index(tris.Len() - 1)
			tris.Resize(tris.Len() - 4)
			i--
		}
	}
}

// Calculate minimum extend of the polygon.
func polyMinExtent(verts []float64, nverts int) float64 {
	minDist := math.MaxFloat64
	for i := 0; i < nverts; i++ {
		ni := (i + 1) % nverts
		p1 := rcGetVert(verts, i)
		p2 := rcGetVert(verts, ni)
		maxEdgeDist := float64(0)
		for j := 0; j < nverts; j++ {
			if j == i || j == ni {
				continue
			}
			d := distancePtSeg2d(rcGetVert(verts, j), p1, p2)
			maxEdgeDist = common.Max(maxEdgeDist, d)
		}
		minDist = common.Min(minDist, maxEdgeDist)
	}
	return common.Sqrt(minDist)
}
func triangulateHull(nverts int, verts []float64, nhull int, hull []int, nin int, tris Stack[int]) {
	start := 0
	left := 1
	right := nhull - 1

	// Start from an ear with shortest perimeter.
	// This tends to favor well formed triangles as starting point.
	dmin := math.MaxFloat64
	for i := 0; i < nhull; i++ {
		if hull[i] >= nin {
			continue
		} // Ears are triangles with original vertices as middle vertex while others are actually line segments on edges
		pi := prev(i, nhull)
		ni := next(i, nhull)
		pv := rcGetVert(verts, hull[pi])
		cv := rcGetVert(verts, hull[i])
		nv := rcGetVert(verts, hull[ni])
		d := vdist2(pv, cv) + vdist2(cv, nv) + vdist2(nv, pv)
		if d < dmin {
			start = i
			left = ni
			right = pi
			dmin = d
		}
	}

	// Add first triangle
	tris.Push(hull[start])
	tris.Push(hull[left])
	tris.Push(hull[right])
	tris.Push(0)

	// Triangulate the polygon by moving left or right,
	// depending on which triangle has shorter perimeter.
	// This heuristic was chose empirically, since it seems
	// handle tessellated straight edges well.
	for next(left, nhull) != right {
		// Check to see if se should advance left or right.
		nleft := next(left, nhull)
		nright := prev(right, nhull)

		cvleft := rcGetVert(verts, hull[left])
		nvleft := rcGetVert(verts, hull[nleft])
		cvright := rcGetVert(verts, hull[right])
		nvright := rcGetVert(verts, hull[nright])
		dleft := vdist2(cvleft, nvleft) + vdist2(nvleft, cvright)
		dright := vdist2(cvright, nvright) + vdist2(cvleft, nvright)

		if dleft < dright {
			tris.Push(hull[left])
			tris.Push(hull[nleft])
			tris.Push(hull[right])
			tris.Push(0)
			left = nleft
		} else {
			tris.Push(hull[left])
			tris.Push(hull[nright])
			tris.Push(hull[right])
			tris.Push(0)
			right = nright
		}
	}
}

func getJitterX(i int) float64 {
	return (float64((i*0x8da6b343)&0xffff) / 65535.0 * 2.0) - 1.0
}

func getJitterY(i int) float64 {
	return (float64((i*0xd8163841)&0xffff) / 65535.0 * 2.0) - 1.0
}

func buildPolyDetail(in []float64, nin int,
	sampleDist, sampleMaxError float64,
	heightSearchRadius int, chf *RcCompactHeightfield,
	hp *rcHeightPatch, verts []float64, nverts *int,
	tris, edges, samples Stack[int]) bool {
	MAX_VERTS := 127
	MAX_TRIS := 255 // Max tris for delaunay is 2n-2-k (n=num verts, k=num hull verts).
	MAX_VERTS_PER_EDGE := 32
	edge := make([]float64, (MAX_VERTS_PER_EDGE+1)*3)
	hull := make([]int, MAX_VERTS)
	nhull := 0

	*nverts = nin

	for i := 0; i < nin; i++ {
		copy(rcGetVert(verts, i), rcGetVert(in, i))
	}

	edges.Clear()
	tris.Clear()

	cs := chf.Cs
	ics := 1.0 / cs

	// Calculate minimum extents of the polygon based on input data.
	minExtent := polyMinExtent(verts, *nverts)

	// Tessellate outlines.
	// This is done in separate pass in order to ensure
	// seamless height values across the ply boundaries.
	if sampleDist > 0 {
		i := 0
		j := nin - 1
		for i < nin {
			vj := rcGetVert(in, j)
			vi := rcGetVert(in, i)
			swapped := false
			// Make sure the segments are always handled in same order
			// using lexological sort or else there will be seams.
			if math.Floor(vj[0]-vi[0]) < 1e-6 {
				if vj[2] > vi[2] {
					vj, vi = vi, vj
					swapped = true
				}
			} else {
				if vj[0] > vi[0] {
					vj, vi = vi, vj
					swapped = true
				}
			}
			// Create samples along the edge.
			dx := vi[0] - vj[0]
			dy := vi[1] - vj[1]
			dz := vi[2] - vj[2]
			d := math.Sqrt(dx*dx + dz*dz)
			nn := int(1 + math.Floor(d/sampleDist))
			if nn >= MAX_VERTS_PER_EDGE {
				nn = MAX_VERTS_PER_EDGE - 1
			}
			if *nverts+nn >= MAX_VERTS {
				nn = MAX_VERTS - 1 - *nverts
			}

			for k := 0; k <= nn; k++ {
				u := float64(k) / float64(nn)
				pos := rcGetVert(edge, k)
				pos[0] = vj[0] + dx*u
				pos[1] = vj[1] + dy*u
				pos[2] = vj[2] + dz*u
				pos[1] = float64(getHeight(pos[0], pos[1], pos[2], cs, ics, chf.Ch, heightSearchRadius, hp)) * chf.Ch
			}
			// Simplify samples.
			idx := make([]int, MAX_VERTS_PER_EDGE)
			idx[0] = 0
			idx[1] = nn
			nidx := 2
			for k := 0; k < nidx-1; {
				a := idx[k]
				b := idx[k+1]
				va := rcGetVert(edge, a)
				vb := rcGetVert(edge, b)
				// Find maximum deviation along the segment.
				maxd := float64(0)
				maxi := -1
				for m := a + 1; m < b; m++ {
					dev := distancePtSeg(rcGetVert(edge, m), va, vb)
					if dev > maxd {
						maxd = dev
						maxi = m
					}
				}
				// If the max deviation is larger than accepted error,
				// add new point, else continue to next segment.
				if maxi != -1 && maxd > common.Sqr(sampleMaxError) {
					for m := nidx; m > k; m-- {
						idx[m] = idx[m-1]
					}

					idx[k+1] = maxi
					nidx++
				} else {
					k++
				}
			}

			hull[nhull] = j
			nhull++
			// Add new vertices.
			if swapped {
				for k := nidx - 2; k > 0; k-- {
					copy(rcGetVert(verts, *nverts), rcGetVert(edge, idx[k]*3))
					hull[nhull] = *nverts
					nhull++
					*nverts++
				}
			} else {
				for k := 1; k < nidx-1; k++ {
					copy(rcGetVert(verts, *nverts), rcGetVert(edge, idx[k]))
					hull[nhull] = *nverts
					nhull++
					*nverts++
				}
			}
			j = i
			i++
		}
	}

	// If the polygon minimum extent is small (sliver or small triangle), do not try to add internal points.
	if minExtent < sampleDist*2 {
		triangulateHull(*nverts, verts, nhull, hull, nin, tris)
		return true
	}

	// Tessellate the base mesh.
	// We're using the triangulateHull instead of delaunayHull as it tends to
	// create a bit better triangulation for long thin triangles when there
	// are no internal points.
	triangulateHull(*nverts, verts, nhull, hull, nin, tris)

	if tris.Len() == 0 {
		// Could not triangulate the poly, make sure there is some valid data there.
		return true
	}

	if sampleDist > 0 {
		// Create sample locations in a grid.
		bmin := make([]float64, 3)
		bmax := make([]float64, 3)
		copy(bmin, in)
		copy(bmax, in)
		for i := 1; i < nin; i++ {
			common.Vmin(bmin, rcGetVert(in, i))
			common.Vmax(bmax, rcGetVert(in, i))
		}
		x0 := int(math.Floor(bmin[0] / sampleDist))
		x1 := int(math.Ceil(bmax[0] / sampleDist))
		z0 := int(math.Floor(bmin[2] / sampleDist))
		z1 := int(math.Ceil(bmax[2] / sampleDist))
		samples.Clear()
		for z := z0; z < z1; z++ {
			for x := x0; x < x1; x++ {
				pt := make([]float64, 3)
				pt[0] = float64(x) * sampleDist
				pt[1] = (bmax[1] + bmin[1]) * 0.5
				pt[2] = float64(z) * sampleDist
				// Make sure the samples are not too close to the edges.
				if distToPoly(nin, in, pt) > -sampleDist/2 {
					continue
				}
				samples.Push(x)
				samples.Push(getHeight(pt[0], pt[1], pt[2], cs, ics, chf.Ch, heightSearchRadius, hp))
				samples.Push(z)
				samples.Push(0) // Not added
			}
		}

		// Add the samples starting from the one that has the most
		// error. The procedure stops when all samples are added
		// or when the max error is within treshold.
		nsamples := samples.Len() / 4
		for iter := 0; iter < nsamples; iter++ {
			if *nverts >= MAX_VERTS {
				break
			}

			// Find sample with most error.
			bestpt := []float64{0, 0, 0}
			bestd := float64(0)
			besti := -1
			for i := 0; i < nsamples; i++ {
				s := samples.Slice(i*4, samples.Len())
				if s[3] != 0 {
					continue
				} // skip added.
				pt := make([]float64, 3)
				// The sample location is jittered to get rid of some bad triangulations
				// which are cause by symmetrical data from the grid structure.
				pt[0] = float64(s[0])*sampleDist + getJitterX(i)*float64(cs)*0.1
				pt[1] = float64(s[1]) * chf.Ch
				pt[2] = float64(s[2])*sampleDist + getJitterY(i)*cs*0.1
				d := distToTriMesh(pt, verts, *nverts, tris.Data(), tris.Len()/4)
				if d < 0 {
					continue
				} // did not hit the mesh.
				if d > bestd {
					bestd = d
					besti = i
					copy(bestpt, pt)
				}
			}
			// If the max error is within accepted threshold, stop tesselating.
			if bestd <= sampleMaxError || besti == -1 {
				break
			}

			// Mark sample as added.
			samples.SetByIndex(besti*4+3, 1)
			// Add the new sample point.
			copy(verts[*nverts*3:*nverts*3+3], bestpt)
			*nverts++

			// Create new triangulation.
			// TODO: Incremental add instead of full rebuild.
			edges.Clear()
			tris.Clear()
			delaunayHull(*nverts, verts, nhull, hull, tris, edges)
		}
	}

	ntris := tris.Len() / 4
	if ntris > MAX_TRIS {
		tris.Resize(MAX_TRIS * 4)
	}

	return true
}

func seedArrayWithPolyCenter(chf *RcCompactHeightfield, poly []int, npoly int, verts []int, bs int, hp *rcHeightPatch, array Stack[int]) {
	// Note: Reads to the compact heightfield are offset by border size (bs)
	// since border size offset is already removed from the polymesh vertices.

	offset :=
		[18]int{ //9*2
			0, 0, -1, -1, 0, -1, 1, -1, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0,
		}

	// Find cell closest to a poly vertex
	startCellX := 0
	startCellY := 0
	startSpanIndex := -1
	dmin := RC_UNSET_HEIGHT
	for j := 0; j < npoly && dmin > 0; j++ {
		for k := 0; k < 9 && dmin > 0; k++ {
			ax := verts[poly[j]*3+0] + offset[k*2+0]
			ay := verts[poly[j]*3+1]
			az := verts[poly[j]*3+2] + offset[k*2+1]
			if ax < hp.xmin || ax >= hp.xmin+hp.width ||
				az < hp.ymin || az >= hp.ymin+hp.height {
				continue
			}

			c := chf.Cells[(ax+bs)+(az+bs)*chf.Width]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni && dmin > 0; i++ {
				s := chf.Spans[i]
				d := common.Abs(ay - s.Y)
				if d < dmin {
					startCellX = ax
					startCellY = az
					startSpanIndex = i
					dmin = d
				}
			}
		}
	}

	if startSpanIndex != -1 {
		panic("")
	}
	// Find center of the polygon
	pcx := 0
	pcy := 0
	for j := 0; j < npoly; j++ {
		pcx += verts[poly[j]*3+0]
		pcy += verts[poly[j]*3+2]
	}
	pcx /= npoly
	pcy /= npoly

	// Use seeds array as a stack for DFS
	array.Clear()
	array.Push(startCellX)
	array.Push(startCellY)
	array.Push(startSpanIndex)

	dirs := []int{0, 1, 2, 3}
	hp.data = make([]int, hp.width*hp.height)
	// DFS to move to the center. Note that we need a DFS here and can not just move
	// directly towards the center without recording intermediate nodes, even though the polygons
	// are convex. In very rare we can get stuck due to contour simplification if we do not
	// record nodes.
	cx := -1
	cy := -1
	ci := -1
	for true {
		if array.Len() < 3 {
			break
		}

		ci = array.Pop()
		cy = array.Pop()
		cx = array.Pop()

		if cx == pcx && cy == pcy {
			break
		}

		// If we are already at the correct X-position, prefer direction
		// directly towards the center in the Y-axis; otherwise prefer
		// direction in the X-axis
		var directDir int
		if cx == pcx {
			y := -1
			if pcy > cy {
				y = 1
			}
			directDir = common.GetDirForOffset(0, y)
		} else {
			x := -1
			if pcx > cx {
				x = 1
			}
			directDir = common.GetDirForOffset(x, 0)
		}

		// Push the direct dir last so we start with this on next iteration
		dirs[directDir], dirs[3] = dirs[3], dirs[directDir]

		cs := chf.Spans[ci]
		for i := 0; i < 4; i++ {
			dir := dirs[i]
			if rcGetCon(cs, dir) == RC_NOT_CONNECTED {
				continue
			}

			newX := cx + common.GetDirOffsetX(dir)
			newY := cy + common.GetDirOffsetY(dir)

			hpx := newX - hp.xmin
			hpy := newY - hp.ymin
			if hpx < 0 || hpx >= hp.width || hpy < 0 || hpy >= hp.height {
				continue
			}

			if hp.data[hpx+hpy*hp.width] != 0 {
				continue
			}

			hp.data[hpx+hpy*hp.width] = 1
			array.Push(newX)
			array.Push(newY)
			array.Push(chf.Cells[(newX+bs)+(newY+bs)*chf.Width].Index + rcGetCon(cs, dir))
		}

		dirs[directDir], dirs[3] = dirs[3], dirs[directDir]
	}

	array.Clear()
	// getHeightData seeds are given in coordinates with borders
	array.Push(cx + bs)
	array.Push(cy + bs)
	array.Push(ci)
	for i := range hp.data {
		hp.data[i] = 0xff
	}

	cs := chf.Spans[ci]
	hp.data[cx-hp.xmin+(cy-hp.ymin)*hp.width] = cs.Y
}

func push3(queue Stack[int], v1, v2, v3 int) {
	queue.Resize(queue.Len() + 3)
	queue.SetByIndex(queue.Len()-3, v1)
	queue.SetByIndex(queue.Len()-2, v2)
	queue.SetByIndex(queue.Len()-1, v3)
}

func etHeightData(chf *RcCompactHeightfield,
	poly []int, npoly int,
	verts []int, bs int,
	hp *rcHeightPatch, queue Stack[int],
	region int) {
	// Note: Reads to the compact heightfield are offset by border size (bs)
	// since border size offset is already removed from the polymesh vertices.

	queue.Clear()
	if hp.data == nil {
		hp.data = make([]int, hp.width*hp.height)
	}
	for i := range hp.data {
		hp.data[i] = 0xff
	}
	// Set all heights to RC_UNSET_HEIGHT.

	empty := true

	// We cannot sample from this poly if it was created from polys
	// of different regions. If it was then it could potentially be overlapping
	// with polys of that region and the heights sampled here could be wrong.
	if region != RC_MULTIPLE_REGS {
		// Copy the height from the same region, and mark region borders
		// as seed points to fill the rest.
		for hy := 0; hy < hp.height; hy++ {
			y := hp.ymin + hy + bs
			for hx := 0; hx < hp.width; hx++ {
				x := hp.xmin + hx + bs
				c := chf.Cells[x+y*chf.Width]
				i := c.Index
				ni := (c.Index + c.Count)
				for ; i < ni; i++ {
					s := chf.Spans[i]
					if s.Reg == region {
						// Store height
						hp.data[hx+hy*hp.width] = s.Y
						empty = false

						// If any of the neighbours is not in same region,
						// add the current location as flood fill start
						border := false
						for dir := 0; dir < 4; dir++ {
							if rcGetCon(s, dir) != RC_NOT_CONNECTED {
								ax := x + common.GetDirOffsetX(dir)
								ay := y + common.GetDirOffsetY(dir)
								ai := chf.Cells[ax+ay*chf.Width].Index + rcGetCon(s, dir)
								as := chf.Spans[ai]
								if as.Reg != region {
									border = true
									break
								}
							}
						}
						if border {
							push3(queue, x, y, i)
						}

						break
					}
				}
			}
		}
	}

	// if the polygon does not contain any points from the current region (rare, but happens)
	// or if it could potentially be overlapping polygons of the same region,
	// then use the center as the seed point.
	if empty {
		seedArrayWithPolyCenter(chf, poly, npoly, verts, bs, hp, queue)
	}

	RETRACT_SIZE := 256
	head := 0

	// We assume the seed is centered in the polygon, so a BFS to collect
	// height data will ensure we do not move onto overlapping polygons and
	// sample wrong heights.
	for head*3 < queue.Len() {
		cx := queue.Index(head*3 + 0)
		cy := queue.Index(head*3 + 1)
		ci := queue.Index(head*3 + 2)
		head++
		if head >= RETRACT_SIZE {
			head = 0
			if queue.Len() > RETRACT_SIZE*3 {
				copy(queue.Data(), queue.Data()[RETRACT_SIZE*3:queue.Len()-RETRACT_SIZE*3])
			}

			queue.Resize(queue.Len() - RETRACT_SIZE*3)
		}

		cs := chf.Spans[ci]
		for dir := 0; dir < 4; dir++ {
			if rcGetCon(cs, dir) == RC_NOT_CONNECTED {
				continue
			}

			ax := cx + common.GetDirOffsetX(dir)
			ay := cy + common.GetDirOffsetY(dir)
			hx := ax - hp.xmin - bs
			hy := ay - hp.ymin - bs

			if hx >= hp.width || hy >= hp.height {
				continue
			}

			if hp.data[hx+hy*hp.width] != RC_UNSET_HEIGHT {
				continue
			}

			ai := chf.Cells[ax+ay*chf.Width].Index + rcGetCon(cs, dir)
			as := chf.Spans[ai]

			hp.data[hx+hy*hp.width] = as.Y

			push3(queue, ax, ay, ai)
		}
	}
}
func getHeightData(chf *RcCompactHeightfield,
	poly []int, npoly int,
	verts []int, bs int,
	hp *rcHeightPatch, queue Stack[int],
	region int) {
	// Note: Reads to the compact heightfield are offset by border size (bs)
	// since border size offset is already removed from the polymesh vertices.

	queue.Clear()
	// Set all heights to RC_UNSET_HEIGHT.
	if hp.data == nil {
		hp.data = make([]int, hp.width*hp.height)
	}
	for i := range hp.data {
		hp.data[i] = 0xff
	}
	empty := true

	// We cannot sample from this poly if it was created from polys
	// of different regions. If it was then it could potentially be overlapping
	// with polys of that region and the heights sampled here could be wrong.
	if region != RC_MULTIPLE_REGS {
		// Copy the height from the same region, and mark region borders
		// as seed points to fill the rest.
		for hy := 0; hy < hp.height; hy++ {
			y := hp.ymin + hy + bs
			for hx := 0; hx < hp.width; hx++ {
				x := hp.xmin + hx + bs
				c := chf.Cells[x+y*chf.Width]
				i := c.Index
				ni := (c.Index + c.Count)
				for ; i < ni; i++ {
					s := chf.Spans[i]
					if s.Reg == region {
						// Store height
						hp.data[hx+hy*hp.width] = s.Y
						empty = false

						// If any of the neighbours is not in same region,
						// add the current location as flood fill start
						border := false
						for dir := 0; dir < 4; dir++ {
							if rcGetCon(s, dir) != RC_NOT_CONNECTED {
								ax := x + common.GetDirOffsetX(dir)
								ay := y + common.GetDirOffsetY(dir)
								ai := chf.Cells[ax+ay*chf.Width].Index + rcGetCon(s, dir)
								as := chf.Spans[ai]
								if as.Reg != region {
									border = true
									break
								}
							}
						}
						if border {
							push3(queue, x, y, i)
						}

						break
					}
				}
			}
		}
	}

	// if the polygon does not contain any points from the current region (rare, but happens)
	// or if it could potentially be overlapping polygons of the same region,
	// then use the center as the seed point.
	if empty {
		seedArrayWithPolyCenter(chf, poly, npoly, verts, bs, hp, queue)
	}

	RETRACT_SIZE := 256
	head := 0

	// We assume the seed is centered in the polygon, so a BFS to collect
	// height data will ensure we do not move onto overlapping polygons and
	// sample wrong heights.
	for head*3 < queue.Len() {
		cx := queue.Index(head*3 + 0)
		cy := queue.Index(head*3 + 1)
		ci := queue.Index(head*3 + 2)
		head++
		if head >= RETRACT_SIZE {
			head = 0
			if queue.Len() > RETRACT_SIZE*3 {
				copy(queue.Data(), queue.Data()[RETRACT_SIZE*3:RETRACT_SIZE*3+queue.Len()-RETRACT_SIZE*3])
			}

			queue.Resize(queue.Len() - RETRACT_SIZE*3)
		}

		cs := chf.Spans[ci]
		for dir := 0; dir < 4; dir++ {
			if rcGetCon(cs, dir) == RC_NOT_CONNECTED {
				continue
			}

			ax := cx + common.GetDirOffsetX(dir)
			ay := cy + common.GetDirOffsetY(dir)
			hx := ax - hp.xmin - bs
			hy := ay - hp.ymin - bs

			if hx >= hp.width || hy >= hp.height {
				continue
			}

			if hp.data[hx+hy*hp.width] != RC_UNSET_HEIGHT {
				continue
			}

			ai := chf.Cells[ax+ay*chf.Width].Index + rcGetCon(cs, dir)
			as := chf.Spans[ai]

			hp.data[hx+hy*hp.width] = as.Y

			push3(queue, ax, ay, ai)
		}
	}
}

func getEdgeFlags(va, vb []float64,
	vpoly []float64, npoly int) int {
	// The flag returned by this function matches dtDetailTriEdgeFlags in Detour.
	// Figure out if edge (va,vb) is part of the polygon boundary.
	thrSqr := common.Sqr(0.001)
	i := 0
	j := npoly - 1
	for i < npoly {

		if distancePtSeg2d(va, rcGetVert(vpoly, j), rcGetVert(vpoly, i)) < thrSqr &&
			distancePtSeg2d(vb, rcGetVert(vpoly, j), rcGetVert(vpoly, i)) < thrSqr {
			return 1
		}
		j = i
		i++
	}
	return 0
}

func getTriFlags(va, vb, vc []float64,
	vpoly []float64, npoly int) int {
	flags := 0
	flags |= getEdgeFlags(va, vb, vpoly, npoly) << 0
	flags |= getEdgeFlags(vb, vc, vpoly, npoly) << 2
	flags |= getEdgeFlags(vc, va, vpoly, npoly) << 4
	return flags
}

// / Contains triangle meshes that represent detailed height data associated
// / with the polygons in its associated polygon mesh object.
// / @ingroup recast
type RcPolyMeshDetail struct {
	Meshes  []int     ///< The sub-mesh data. [Size: 4*#nmeshes]
	Verts   []float64 ///< The mesh vertices. [Size: 3*#nverts]
	Tris    []int     ///< The mesh triangles. [Size: 4*#ntris]
	Nmeshes int       ///< The number of sub-meshes defined by #meshes.
	Nverts  int       ///< The number of vertices in #verts.
	Ntris   int       ///< The number of triangles in #tris.
	// Explicitly-disabled copy constructor and copy assignment operator.
}

// / @par
// /
// / See the #RcConfig documentation for more information on the configuration parameters.
// /
// / @see rcAllocPolyMeshDetail, RcPolyMesh, RcCompactHeightfield, RcPolyMeshDetail, RcConfig
func RcBuildPolyMeshDetail(mesh *RcPolyMesh, chf *RcCompactHeightfield, sampleDist float64, sampleMaxError float64,
	dmesh *RcPolyMeshDetail) bool {
	if mesh.Nverts == 0 || mesh.Npolys == 0 {
		return true
	}

	nvp := mesh.Nvp
	cs := mesh.Cs
	ch := mesh.Ch
	orig := mesh.Bmin
	borderSize := mesh.BorderSize
	heightSearchRadius := int(common.Max(1, math.Ceil(mesh.MaxEdgeError)))
	c := func() int { return 0 }
	edges := NewStackArray(c, 64)
	tris := NewStackArray(c, 512)
	arr := NewStackArray(c, 512)
	samples := NewStackArray(c, 512)
	verts := make([]float64, 256*3)
	var hp rcHeightPatch
	nPolyVerts := 0
	maxhw := 0
	maxhh := 0
	bounds := make([]int, mesh.Npolys*4)

	poly := make([]float64, nvp*3)
	// Find max size for a polygon area.
	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		xmin := &bounds[i*4+0]
		xmax := &bounds[i*4+1]
		ymin := &bounds[i*4+2]
		ymax := &bounds[i*4+3]
		*xmin = chf.Width
		*xmax = 0
		*ymin = chf.Height
		*ymax = 0
		for j := 0; j < nvp; j++ {
			if p[j] == RC_MESH_NULL_IDX {
				break
			}
			v := rcGetVert(mesh.Verts, p[j])
			*xmin = common.Min(*xmin, v[0])
			*xmax = common.Max(*xmax, v[0])
			*ymin = common.Min(*ymin, v[2])
			*ymax = common.Max(*ymax, v[2])
			nPolyVerts++
		}
		*xmin = common.Max(0, *xmin-1)
		*xmax = common.Min(chf.Width, *xmax+1)
		*ymin = common.Max(0, *ymin-1)
		*ymax = common.Min(chf.Height, *ymax+1)
		if *xmin >= *xmax || *ymin >= *ymax {
			continue
		}
		maxhw = common.Max(maxhw, *xmax-*xmin)
		maxhh = common.Max(maxhh, *ymax-*ymin)
	}

	hp.data = make([]int, maxhw*maxhh)

	dmesh.Nmeshes = mesh.Npolys
	dmesh.Nverts = 0
	dmesh.Ntris = 0

	dmesh.Meshes = make([]int, dmesh.Nmeshes*4)

	vcap := nPolyVerts + nPolyVerts/2
	tcap := vcap * 2

	dmesh.Nverts = 0

	dmesh.Verts = make([]float64, vcap*3)
	dmesh.Ntris = 0

	dmesh.Tris = make([]int, tcap*4)

	for i := 0; i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]

		// Store polygon vertices for processing.
		npoly := 0
		for j := 0; j < nvp; j++ {
			if p[j] == RC_MESH_NULL_IDX {
				break
			}
			v := rcGetVert(mesh.Verts, p[j])
			poly[j*3+0] = float64(v[0]) * cs
			poly[j*3+1] = float64(v[1]) * ch
			poly[j*3+2] = float64(v[2]) * cs
			npoly++
		}

		// Get the height data from the area of the polygon.
		hp.xmin = bounds[i*4+0]
		hp.ymin = bounds[i*4+2]
		hp.width = bounds[i*4+1] - bounds[i*4+0]
		hp.height = bounds[i*4+3] - bounds[i*4+2]
		getHeightData(chf, p, npoly, mesh.Verts, borderSize, &hp, arr, mesh.Regs[i])

		// Build detail mesh.
		nverts := 0
		if !buildPolyDetail(poly, npoly,
			sampleDist, sampleMaxError,
			heightSearchRadius, chf, &hp,
			verts, &nverts, tris,
			edges, samples) {
			return false
		}

		// Move detail verts to world space.
		for j := 0; j < nverts; j++ {
			verts[j*3+0] += orig[0]
			verts[j*3+1] += orig[1] + chf.Ch // Is this offset necessary?
			verts[j*3+2] += orig[2]
		}
		// Offset poly too, will be used to flag checking.
		for j := 0; j < npoly; j++ {
			poly[j*3+0] += orig[0]
			poly[j*3+1] += orig[1]
			poly[j*3+2] += orig[2]
		}

		// Store detail submesh.
		ntris := tris.Len() / 4

		dmesh.Meshes[i*4+0] = dmesh.Nverts
		dmesh.Meshes[i*4+1] = nverts
		dmesh.Meshes[i*4+2] = dmesh.Ntris
		dmesh.Meshes[i*4+3] = ntris

		// Store vertices, allocate more memory if necessary.
		if dmesh.Nverts+nverts > vcap {
			for dmesh.Nverts+nverts > vcap {
				vcap += 256
			}

			newv := make([]float64, vcap*3)
			if dmesh.Nverts != 0 {
				copy(newv, dmesh.Verts[:3*dmesh.Nverts])
			}
			dmesh.Verts = newv
		}
		for j := 0; j < nverts; j++ {
			dmesh.Verts[dmesh.Nverts*3+0] = verts[j*3+0]
			dmesh.Verts[dmesh.Nverts*3+1] = verts[j*3+1]
			dmesh.Verts[dmesh.Nverts*3+2] = verts[j*3+2]
			dmesh.Nverts++
		}

		// Store triangles, allocate more memory if necessary.
		if dmesh.Ntris+ntris > tcap {
			for dmesh.Ntris+ntris > tcap {
				tcap += 256
			}

			newt := make([]int, tcap*4)

			if dmesh.Ntris != 0 {
				copy(newt, dmesh.Tris[:4*dmesh.Ntris])
			}
			dmesh.Tris = newt
		}
		for j := 0; j < ntris; j++ {
			t := tris.Slice(j*4, tris.Len())
			dmesh.Tris[dmesh.Ntris*4+0] = t[0]
			dmesh.Tris[dmesh.Ntris*4+1] = t[1]
			dmesh.Tris[dmesh.Ntris*4+2] = t[2]
			dmesh.Tris[dmesh.Ntris*4+3] = getTriFlags(rcGetVert(verts, t[0]), rcGetVert(verts, t[1]), rcGetVert(verts, t[2]), poly, npoly)
			dmesh.Ntris++
		}
	}

	return true
}

// / @see rcAllocPolyMeshDetail, RcPolyMeshDetail
func rcMergePolyMeshDetails(meshes []*RcPolyMeshDetail, nmeshes int, mesh *RcPolyMeshDetail) bool {

	maxVerts := 0
	maxTris := 0
	maxMeshes := 0

	for i := 0; i < nmeshes; i++ {
		if meshes[i] == nil {
			continue
		}
		maxVerts += meshes[i].Nverts
		maxTris += meshes[i].Ntris
		maxMeshes += meshes[i].Nmeshes
	}

	mesh.Nmeshes = 0
	mesh.Meshes = make([]int, maxMeshes*4)
	mesh.Ntris = 0
	mesh.Tris = make([]int, maxTris*4)
	mesh.Nverts = 0

	mesh.Verts = make([]float64, maxVerts*3)
	// Merge datas.
	for i := 0; i < nmeshes; i++ {
		dm := meshes[i]
		if dm == nil {
			continue
		}
		for j := 0; j < dm.Nmeshes; j++ {
			dst := rcGetVert4(mesh.Meshes, mesh.Nmeshes)
			src := rcGetVert4(dm.Meshes, j)
			dst[0] = mesh.Nverts + src[0]
			dst[1] = src[1]
			dst[2] = mesh.Ntris + src[2]
			dst[3] = src[3]
			mesh.Nmeshes++
		}

		for k := 0; k < dm.Nverts; k++ {
			copy(rcGetVert(mesh.Verts, mesh.Nverts), rcGetVert(dm.Verts, k))
			mesh.Nverts++
		}
		for k := 0; k < dm.Ntris; k++ {
			mesh.Tris[mesh.Ntris*4+0] = dm.Tris[k*4+0]
			mesh.Tris[mesh.Ntris*4+1] = dm.Tris[k*4+1]
			mesh.Tris[mesh.Ntris*4+2] = dm.Tris[k*4+2]
			mesh.Tris[mesh.Ntris*4+3] = dm.Tris[k*4+3]
			mesh.Ntris++
		}
	}

	return true
}
