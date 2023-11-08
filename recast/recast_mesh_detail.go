package recast

import "math"

const (
	RC_UNSET_HEIGHT = 0xffff
)

type rcHeightPatch struct {
	data                      []int
	xmin, ymin, width, height int
}

func vdot2(a, b *Point) float64 {
	return a.X*b.X + a.Y*b.Y
}

func vdistSq2(p, q *Point) float64 {
	dx := q.X - p.X
	dy := q.Y - p.Y
	return dx*dx + dy*dy
}

func vdist2(p, q *Point) float64 {
	return math.Sqrt(vdistSq2(p, q))
}

func vcross2(p1, p2, p3 *Point) float64 {
	u1 := p2.X - p1.X
	v1 := p2.Y - p1.Y
	u2 := p3.X - p1.X
	v2 := p3.Y - p1.Y
	return u1*v2 - v1*u2
}

// 外接圆
func circumCircle(p1, p2, p3 *Point, r float64) bool {
	var c *Point
	EPS := 1e-6
	// Calculate the circle relative to p1, to avoid some precision issues.
	v1 := &Point{}
	v2 := rcVsub(p2, p1)
	v3 := rcVsub(p3, p1)

	cp := vcross2(v1, v2, v3)
	if math.Abs(cp) > EPS {
		v1Sq := vdot2(v1, v1)
		v2Sq := vdot2(v2, v2)
		v3Sq := vdot2(v3, v3)
		c.X = (v1Sq*(v2.Z-v3.Z) + v2Sq*(v3.Z-v1.Z) + v3Sq*(v1.Z-v2.Z)) / (2 * cp)
		c.Y = 0
		c.Z = (v1Sq*(v3.X-v2.X) + v2Sq*(v1.X-v3.X) + v3Sq*(v2.X-v1.X)) / (2 * cp)
		r = vdist2(c, v1)
		c = rcVadd(c, p1)
		return true
	}
	c = p1.Copy()
	r = 0
	return false
}

func distPtTri(p, a, b, c *Point) float64 {
	v0 := rcVsub(c, a)
	v1 := rcVsub(b, a)
	v2 := rcVsub(p, a)

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
		y := a.Y + v0.Y*u + v1.Y*v
		return math.Abs(y - p.Y)
	}
	return math.MaxFloat64
}

func distancePtSeg(pt, p, q *Point) float64 {
	pqx := q.X - p.X
	pqy := q.Y - p.Y
	pqz := q.Z - p.Z
	dx := pt.X - p.X
	dy := pt.Y - p.Y
	dz := pt.Z - p.Z
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

	dx = p.X + t*pqx - pt.X
	dy = p.Y + t*pqy - pt.Y
	dz = p.Z + t*pqz - pt.Z

	return dx*dx + dy*dy + dz*dz
}

func distancePtSeg2d(pt, p, q *Point) float64 {
	pqx := q.X - p.X
	pqz := q.Z - p.Z
	dx := pt.X - p.X
	dz := pt.Z - p.Z
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

	dx = p.X + t*pqx - pt.X
	dz = p.Z + t*pqz - pt.Z

	return dx*dx + dz*dz
}

func distToTriMesh(p *Point, tris [][]*Point) float64 {
	dmin := math.MaxFloat64
	for i := 0; i < len(tris); i++ {
		va := tris[i][0]
		vb := tris[i][1]
		vc := tris[i][2]
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

func distToPoly(verts []*Point, p *Point) float64 {

	dmin := math.MaxFloat64
	c := false
	i := 0
	j := len(verts) - 1
	for i < len(verts) {
		vi := verts[i*3]
		vj := verts[j*3]
		if ((vi.Y > p.Y) != (vj.Y > p.Y)) && (p.X < (vj.X-vi.X)*(p.Y-vi.Y)/(vj.Y-vi.Y)+vi.X) {
			c = true
		}
		dmin = math.Min(dmin, distancePtSeg2d(p, vj, vi))
		j = i
		i++

	}
	if c {
		return -dmin
	}
	return dmin

}

func getHeight(fx, fy, fz float64, ics, ch float64, radius int, hp *rcHeightPatch) int {
	ix := int(math.Floor(fx*ics + 0.01))
	iz := int(math.Floor(fz*ics + 0.01))
	ix = rcClamp(ix-hp.xmin, 0, hp.width-1)
	iz = rcClamp(iz-hp.ymin, 0, hp.height-1)
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

type Edge struct {
	s int
	t int
	l int
	r int
}

func findEdge(edges []*Edge, s, t int) int {
	for i := 0; i < len(edges); i++ {
		e := edges[i]
		if (e.s == s && e.t == t) || (e.s == t && e.t == s) {
			return i
		}

	}
	return EV_UNDEF
}

func addEdge(edges []*Edge, s, t, l, r int) []*Edge {
	// Add edge if not already in the triangulation.
	e := findEdge(edges, s, t)
	if e == EV_UNDEF {
		var edge Edge
		edge.s = s
		edge.t = t
		edge.l = l
		edge.r = r
		edges = append(edges, &edge)
	}
	return edges
}

func updateLeftFace(e *Edge, s, t, f int) {
	if e.s == s && e.t == t && e.l == EV_UNDEF {
		e.l = f
	} else if e.t == s && e.s == t && e.r == EV_UNDEF {
		e.r = f
	}
}

func overlapSegSeg2d(a, b, c, d *Point) int {
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

func polyMinExtent(verts []*Point) float64 {
	minDist := math.MaxFloat64
	for i := 0; i < len(verts); i++ {
		ni := (i + 1) % len(verts)
		p1 := verts[i]
		p2 := verts[ni]
		maxEdgeDist := float64(0)
		for j := 0; j < len(verts); j++ {
			if j == i || j == ni {
				continue
			}
			d := distancePtSeg2d(verts[j], p1, p2)
			maxEdgeDist = math.Max(maxEdgeDist, d)
		}
		minDist = math.Min(minDist, maxEdgeDist)
	}
	return math.Sqrt(minDist)
}

func overlapEdges(pts []*Point, edges []*Edge, s1, t1 int) bool {
	for i := 0; i < len(edges); i++ {
		s0 := edges[i].s
		t0 := edges[i].t
		// Same or connected edges do not overlap.
		if s0 == s1 || s0 == t1 || t0 == s1 || t0 == t1 {
			continue
		}

		if overlapSegSeg2d(pts[s0], pts[t0], pts[s1], pts[t1]) == 1 {
			return true
		}
	}
	return false
}

func getJitterX(i int) float64 {
	return (float64((i*0x8da6b343)&0xffff) / 65535.0 * 2.0) - 1.0
}

func getJitterY(i int) float64 {
	return (float64((i*0xd8163841)&0xffff) / 65535.0 * 2.0) - 1.0
}

func triangulateHull(verts []*Point, hull []int, nin int) {
	nhull := len(hull)
	var tris [][]*Point
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
		pv := verts[hull[pi]]
		cv := verts[hull[i]]
		nv := verts[hull[ni]]
		d := vdist2(pv, cv) + vdist2(cv, nv) + vdist2(nv, pv)
		if d < dmin {
			start = i
			left = ni
			right = pi
			dmin = d
		}
	}

	// Add first triangle
	tris = append(tris, []*Point{verts[hull[start]], verts[hull[left]], verts[hull[right]]})
	// Triangulate the polygon by moving left or right,
	// depending on which triangle has shorter perimeter.
	// This heuristic was chose empirically, since it seems
	// handle tessellated straight edges well.
	for next(left, nhull) != right {
		// Check to see if se should advance left or right.
		nleft := next(left, nhull)
		nright := prev(right, nhull)

		cvleft := verts[hull[left]]
		nvleft := verts[hull[nleft]]
		cvright := verts[hull[right]]
		nvright := verts[hull[nright]]
		dleft := vdist2(cvleft, nvleft) + vdist2(nvleft, cvright)
		dright := vdist2(cvright, nvright) + vdist2(cvleft, nvright)

		if dleft < dright {
			tris = append(tris, []*Point{verts[hull[start]], verts[hull[left]], verts[hull[right]]})
			left = nleft
		} else {
			tris = append(tris, []*Point{verts[hull[start]], verts[hull[left]], verts[hull[right]]})
			right = nright
		}
	}
}
