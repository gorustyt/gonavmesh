package recast

import (
	"github.com/gorustyt/gonavmesh/common"
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
	data                      []uint16
	xmin, ymin, width, height int32
}

func getHeight(fx, fy, fz,
	cs, ics, ch float32,
	radius int32, hp *rcHeightPatch) uint16 {
	ix := int32(math.Floor(float64(fx*ics) + 0.01))
	iz := int32(math.Floor(float64(fz*ics) + 0.01))
	ix = common.Clamp(ix-hp.xmin, 0, hp.width-1)
	iz = common.Clamp(iz-hp.ymin, 0, hp.height-1)
	h := hp.data[ix+iz*hp.width]
	if h == RC_UNSET_HEIGHT {
		// Special case when data might be bad.
		// Walk adjacent cells in a spiral up to 'radius', and look
		// for a pixel which has a valid height.
		x := int32(1)
		z := int32(0)
		dx := int32(1)
		dz := int32(0)
		maxSize := radius*2 + 1
		maxIter := maxSize*maxSize - 1

		nextRingIterStart := int32(8)
		nextRingIters := int32(16)

		dmin := float32(math.MaxFloat32)
		for i := int32(0); i < maxIter; i++ {
			nx := ix + x
			nz := iz + z

			if nx >= 0 && nz >= 0 && nx < hp.width && nz < hp.height {
				nh := hp.data[nx+nz*hp.width]
				if nh != RC_UNSET_HEIGHT {
					d := float32(math.Abs(float64(nh)*float64(ch) - float64(fy)))
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

func findEdge(edges []int32, nedges int32, s, t int32) int32 {
	for i := int32(0); i < nedges; i++ {
		e := common.GetVert4(edges, i)
		if (e[0] == s && e[1] == t) || (e[0] == t && e[1] == s) {
			return i
		}

	}
	return EV_UNDEF
}

func addEdge(edges []int32, nedges *int32, maxEdges int32, s, t, l, r int32) int32 {
	if *nedges >= maxEdges {
		return EV_UNDEF
	}

	// Add edge if not already in the triangulation.
	e := findEdge(edges, *nedges, s, t)
	if e == EV_UNDEF {
		edge := common.GetVert4(edges, *nedges)
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

func updateLeftFace(e []int32, s, t, f int32) {
	if e[0] == s && e[1] == t && e[2] == EV_UNDEF {
		e[2] = f
	} else if e[1] == s && e[0] == t && e[3] == EV_UNDEF {
		e[3] = f
	}

}

func overlapSegSeg2d(a, b, c, d []float32) int32 {
	a1 := common.Vcross2(a, b, d)
	a2 := common.Vcross2(a, b, c)
	if a1*a2 < 0.0 {
		a3 := common.Vcross2(c, d, a)
		a4 := a3 + a2 - a1
		if a3*a4 < 0.0 {
			return 1
		}

	}
	return 0
}

func overlapEdges(pts []float32, edges []int32, nedges, s1, t1 int32) bool {
	for i := int32(0); i < nedges; i++ {
		s0 := edges[i*4+0]
		t0 := edges[i*4+1]
		// Same or connected edges do not overlap.
		if s0 == s1 || s0 == t1 || t0 == s1 || t0 == t1 {
			continue
		}

		if overlapSegSeg2d(common.GetVert3(pts, s0), common.GetVert3(pts, t0), common.GetVert3(pts, s1), common.GetVert3(pts, t1)) != 0 {
			return true
		}

	}
	return false
}

func completeFacet(pts []float32, npts int32, edges []int32, nedges *int32, maxEdges int32, nfaces *int32, e int32) {
	EPS := float32(1e-5)

	edge := common.GetVert4(edges, e)

	// Cache s and t.
	var s, t int32
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
	c := []float32{0, 0, 0}
	r := float32(-1)
	for u := int32(0); u < npts; u++ {
		if u == s || u == t {
			continue
		}
		if common.Vcross2(common.GetVert3(pts, s), common.GetVert3(pts, t), common.GetVert3(pts, u)) > EPS {
			if r < 0 {
				// The circle is not updated yet, do it now.
				pt = u
				common.CircumCircle(common.GetVert3(pts, s), common.GetVert3(pts, t), common.GetVert3(pts, u), c, &r)
				continue
			}
			d := common.Vdist2(c, common.GetVert3(pts, u))
			tol := float32(0.001)
			if d > r*(1+tol) {
				// Outside current circumcircle, skip.
				continue
			} else if d < r*(1-tol) {
				// Inside safe circumcircle, update circle.
				pt = u
				common.CircumCircle(common.GetVert3(pts, s), common.GetVert3(pts, t), common.GetVert3(pts, u), c, &r)
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
				common.CircumCircle(common.GetVert3(pts, s), common.GetVert3(pts, t), common.GetVert3(pts, u), c, &r)
			}
		}
	}

	// Add new triangle or update edge info if s-t is on hull.
	if pt < npts {
		// Update face information of edge being completed.
		updateLeftFace(common.GetVert4(edges, e), s, t, *nfaces)

		// Add new edge or update face info of old edge.
		e = findEdge(edges, *nedges, pt, s)
		if e == EV_UNDEF {
			addEdge(edges, nedges, maxEdges, pt, s, *nfaces, EV_UNDEF)
		} else {
			updateLeftFace(common.GetVert4(edges, e), pt, s, *nfaces)
		}

		// Add new edge or update face info of old edge.
		e = findEdge(edges, *nedges, t, pt)
		if e == EV_UNDEF {
			addEdge(edges, nedges, maxEdges, t, pt, *nfaces, EV_UNDEF)
		} else {
			updateLeftFace(common.GetVert4(edges, e), t, pt, *nfaces)
		}
		*nfaces++
	} else {
		updateLeftFace(common.GetVert4(edges, e), s, t, EV_HULL)
	}
}

func delaunayHull(npts int32, pts []float32,
	nhull int32, hull []int32,
	tris Stack[int32], edges Stack[int32]) {
	nfaces := int32(0)
	nedges := int32(0)
	maxEdges := npts * 10
	edges.Resize(int(maxEdges * 4))

	i := int32(0)
	j := nhull - 1
	for i < nhull {
		addEdge(edges.Data(), &nedges, maxEdges, hull[j], hull[i], EV_HULL, EV_UNDEF)
		j = i
		i++
	}

	currentEdge := int32(0)
	for currentEdge < nedges {
		if edges.Index(int(currentEdge*4+2)) == EV_UNDEF {
			completeFacet(pts, npts, edges.Data(), &nedges, maxEdges, &nfaces, currentEdge)
		}

		if edges.Index(int(currentEdge*4+3)) == EV_UNDEF {
			completeFacet(pts, npts, edges.Data(), &nedges, maxEdges, &nfaces, currentEdge)
		}

		currentEdge++
	}

	// Create tris
	tris.Resize(int(nfaces * 4))
	for i := 0; i < int(nfaces*4); i++ {
		tris.SetByIndex(i, -1)
	}

	for i := 0; i < int(nedges); i++ {
		e := common.GetVert4(edges.Data(), i)
		if e[3] >= 0 {
			// Left face
			t := common.GetVert4(tris.Slice(0, tris.Len()), e[3])
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
			t := common.GetVert4(tris.Data(), e[2])
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
		t := common.GetVert4(tris.Data(), i)
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
func polyMinExtent(verts []float32, nverts int32) float32 {
	minDist := float32(math.MaxFloat32)
	for i := int32(0); i < nverts; i++ {
		ni := (i + 1) % nverts
		p1 := common.GetVert3(verts, i)
		p2 := common.GetVert3(verts, ni)
		maxEdgeDist := float32(0)
		for j := int32(0); j < nverts; j++ {
			if j == i || j == ni {
				continue
			}
			d := common.DistancePtSeg2d(common.GetVert3(verts, j), p1, p2)
			maxEdgeDist = common.Max(maxEdgeDist, d)
		}
		minDist = common.Min(minDist, maxEdgeDist)
	}
	return float32(common.Sqrt(float64(minDist)))
}
func triangulateHull(nverts int32, verts []float32, nhull int32, hull []int32, nin int32, tris Stack[int32]) {
	start := int32(0)
	left := int32(1)
	right := nhull - 1

	// Start from an ear with shortest perimeter.
	// This tends to favor well formed triangles as starting point.
	dmin := float32(math.MaxFloat32)
	for i := int32(0); i < nhull; i++ {
		if hull[i] >= nin {
			continue
		} // Ears are triangles with original vertices as middle vertex while others are actually line segments on edges
		pi := common.Prev(i, nhull)
		ni := common.Next(i, nhull)
		pv := common.GetVert3(verts, hull[pi])
		cv := common.GetVert3(verts, hull[i])
		nv := common.GetVert3(verts, hull[ni])
		d := common.Vdist2(pv, cv) + common.Vdist2(cv, nv) + common.Vdist2(nv, pv)
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
	for common.Next(left, nhull) != right {
		// Check to see if se should advance left or right.
		nleft := common.Next(left, nhull)
		nright := common.Prev(right, nhull)

		cvleft := common.GetVert3(verts, hull[left])
		nvleft := common.GetVert3(verts, hull[nleft])
		cvright := common.GetVert3(verts, hull[right])
		nvright := common.GetVert3(verts, hull[nright])
		dleft := common.Vdist2(cvleft, nvleft) + common.Vdist2(nvleft, cvright)
		dright := common.Vdist2(cvright, nvright) + common.Vdist2(cvleft, nvright)

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

func getJitterX(i int) float32 {
	return (float32((i*0x8da6b343)&0xffff) / 65535.0 * 2.0) - 1.0
}

func getJitterY(i int) float32 {
	return (float32((i*0xd8163841)&0xffff) / 65535.0 * 2.0) - 1.0
}

func buildPolyDetail(in []float32, nin int32,
	sampleDist, sampleMaxError float32,
	heightSearchRadius int32, chf *RcCompactHeightfield,
	hp *rcHeightPatch, verts []float32, nverts *int32,
	tris, edges, samples Stack[int32]) bool {
	MAX_VERTS := int32(127)
	MAX_TRIS := int32(255) // Max tris for delaunay is 2n-2-k (n=num verts, k=num hull verts).
	MAX_VERTS_PER_EDGE := int32(32)
	edge := make([]float32, (MAX_VERTS_PER_EDGE+1)*3)
	hull := make([]int32, MAX_VERTS)
	nhull := int32(0)

	*nverts = nin

	for i := int32(0); i < nin; i++ {
		copy(common.GetVert3(verts, i), common.GetVert3(in, i))
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
		i := int32(0)
		j := nin - 1
		for i < nin {
			vj := common.GetVert3(in, j)
			vi := common.GetVert3(in, i)
			swapped := false
			// Make sure the segments are always handled in same order
			// using lexological sort or else there will be seams.
			if math.Floor(float64(vj[0]-vi[0])) < 1e-6 {
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
			d := math.Sqrt(float64(dx*dx + dz*dz))
			nn := int32(1 + math.Floor(d/float64(sampleDist)))
			if nn >= MAX_VERTS_PER_EDGE {
				nn = MAX_VERTS_PER_EDGE - 1
			}
			if *nverts+nn >= MAX_VERTS {
				nn = MAX_VERTS - 1 - *nverts
			}

			for k := int32(0); k <= nn; k++ {
				u := float32(k) / float32(nn)
				pos := common.GetVert3(edge, k)
				pos[0] = vj[0] + dx*u
				pos[1] = vj[1] + dy*u
				pos[2] = vj[2] + dz*u
				pos[1] = float32(getHeight(pos[0], pos[1], pos[2], cs, ics, chf.Ch, heightSearchRadius, hp)) * chf.Ch
			}
			// Simplify samples.
			idx := make([]int32, MAX_VERTS_PER_EDGE)
			idx[0] = 0
			idx[1] = nn
			nidx := 2
			for k := 0; k < nidx-1; {
				a := idx[k]
				b := idx[k+1]
				va := common.GetVert3(edge, a)
				vb := common.GetVert3(edge, b)
				// Find maximum deviation along the segment.
				maxd := float32(0)
				maxi := int32(-1)
				for m := a + 1; m < b; m++ {
					dev := common.DistancePtSeg(common.GetVert3(edge, m), va, vb)
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
					copy(common.GetVert3(verts, *nverts), common.GetVert3(edge, idx[k]*3))
					hull[nhull] = *nverts
					nhull++
					*nverts++
				}
			} else {
				for k := 1; k < nidx-1; k++ {
					copy(common.GetVert3(verts, *nverts), common.GetVert3(edge, idx[k]))
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
		bmin := make([]float32, 3)
		bmax := make([]float32, 3)
		copy(bmin, in)
		copy(bmax, in)
		for i := int32(1); i < nin; i++ {
			common.Vmin(bmin, common.GetVert3(in, i))
			common.Vmax(bmax, common.GetVert3(in, i))
		}
		x0 := int32(math.Floor(float64(bmin[0] / sampleDist)))
		x1 := int32(math.Ceil(float64(bmax[0] / sampleDist)))
		z0 := int32(math.Floor(float64(bmin[2] / sampleDist)))
		z1 := int32(math.Ceil(float64(bmax[2] / sampleDist)))
		samples.Clear()
		for z := z0; z < z1; z++ {
			for x := x0; x < x1; x++ {
				pt := make([]float32, 3)
				pt[0] = float32(x) * sampleDist
				pt[1] = (bmax[1] + bmin[1]) * 0.5
				pt[2] = float32(z) * sampleDist
				// Make sure the samples are not too close to the edges.
				if common.DistToPoly(nin, in, pt) > -sampleDist/2 {
					continue
				}
				samples.Push(x)
				samples.Push(int32(getHeight(pt[0], pt[1], pt[2], cs, ics, chf.Ch, heightSearchRadius, hp)))
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
			bestpt := []float32{0, 0, 0}
			bestd := float32(0)
			besti := -1
			for i := 0; i < nsamples; i++ {
				s := samples.Slice(i*4, samples.Len())
				if s[3] != 0 {
					continue
				} // skip added.
				pt := make([]float32, 3)
				// The sample location is jittered to get rid of some bad triangulations
				// which are cause by symmetrical data from the grid structure.
				pt[0] = float32(s[0])*sampleDist + getJitterX(i)*float32(cs)*0.1
				pt[1] = float32(s[1]) * chf.Ch
				pt[2] = float32(s[2])*sampleDist + getJitterY(i)*cs*0.1
				d := common.DistToTriMesh(pt, verts, *nverts, tris.Data(), int32(tris.Len()/4))
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

	ntris := int32(tris.Len() / 4)
	if ntris > MAX_TRIS {
		tris.Resize(int(MAX_TRIS * 4))
	}

	return true
}

func seedArrayWithPolyCenter(chf *RcCompactHeightfield, poly []uint16, npoly int32, verts []uint16, bs int32, hp *rcHeightPatch, array Stack[int32]) {
	// Note: Reads to the compact heightfield are offset by border size (bs)
	// since border size offset is already removed from the polymesh vertices.

	offset :=
		[18]int32{ //9*2
			0, 0, -1, -1, 0, -1, 1, -1, 1, 0, 1, 1, 0, 1, -1, 1, -1, 0,
		}

	// Find cell closest to a poly vertex
	startCellX := int32(0)
	startCellY := int32(0)
	startSpanIndex := int32(-1)
	dmin := int32(RC_UNSET_HEIGHT)
	for j := int32(0); j < npoly && dmin > 0; j++ {
		for k := 0; k < 9 && dmin > 0; k++ {
			ax := int32(verts[poly[j]*3+0]) + offset[k*2+0]
			ay := int32(verts[poly[j]*3+1])
			az := int32(verts[poly[j]*3+2]) + offset[k*2+1]
			if ax < hp.xmin || ax >= hp.xmin+hp.width ||
				az < hp.ymin || az >= hp.ymin+hp.height {
				continue
			}

			c := chf.Cells[(ax+bs)+(az+bs)*chf.Width]
			i := c.Index
			ni := (c.Index + c.Count)
			for ; i < ni && dmin > 0; i++ {
				s := chf.Spans[i]
				d := common.Abs(ay - int32(s.Y))
				if d < dmin {
					startCellX = ax
					startCellY = az
					startSpanIndex = int32(i)
					dmin = d
				}
			}
		}
	}

	if startSpanIndex != -1 {
		panic("")
	}
	// Find center of the polygon
	pcx := int32(0)
	pcy := int32(0)
	for j := int32(0); j < npoly; j++ {
		pcx += int32(verts[poly[j]*3+0])
		pcy += int32(verts[poly[j]*3+2])
	}
	pcx /= npoly
	pcy /= npoly

	// Use seeds array as a stack for DFS
	array.Clear()
	array.Push(startCellX)
	array.Push(startCellY)
	array.Push(startSpanIndex)

	dirs := []int32{0, 1, 2, 3}
	hp.data = make([]uint16, hp.width*hp.height)
	// DFS to move to the center. Note that we need a DFS here and can not just move
	// directly towards the center without recording intermediate nodes, even though the polygons
	// are convex. In very rare we can get stuck due to contour simplification if we do not
	// record nodes.
	cx := int32(-1)
	cy := int32(-1)
	ci := int32(-1)
	for {
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
		var directDir int32
		if cx == pcx {
			y := int32(-1)
			if pcy > cy {
				y = 1
			}
			directDir = common.GetDirForOffset(0, y)
		} else {
			x := int32(-1)
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
			array.Push(int32(chf.Cells[(newX+bs)+(newY+bs)*chf.Width].Index) + rcGetCon(cs, dir))
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

func push3(queue Stack[int32], v1, v2, v3 int32) {
	queue.Resize(queue.Len() + 3)
	queue.SetByIndex(queue.Len()-3, v1)
	queue.SetByIndex(queue.Len()-2, v2)
	queue.SetByIndex(queue.Len()-1, v3)
}

func getHeightData(chf *RcCompactHeightfield,
	poly []uint16, npoly int32,
	verts []uint16, bs int32,
	hp *rcHeightPatch, queue Stack[int32],
	region int32) {
	// Note: Reads to the compact heightfield are offset by border size (bs)
	// since border size offset is already removed from the polymesh vertices.

	queue.Clear()
	// Set all heights to RC_UNSET_HEIGHT.
	if hp.data == nil {
		hp.data = make([]uint16, hp.width*hp.height)
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
		for hy := int32(0); hy < hp.height; hy++ {
			y := hp.ymin + hy + bs
			for hx := int32(0); hx < hp.width; hx++ {
				x := hp.xmin + hx + bs
				c := chf.Cells[x+y*chf.Width]
				i := c.Index
				ni := (c.Index + c.Count)
				for ; i < ni; i++ {
					s := chf.Spans[i]
					if int32(s.Reg) == region {
						// Store height
						hp.data[hx+hy*hp.width] = s.Y
						empty = false

						// If any of the neighbours is not in same region,
						// add the current location as flood fill start
						border := false
						for dir := int32(0); dir < 4; dir++ {
							if rcGetCon(s, dir) != RC_NOT_CONNECTED {
								ax := x + common.GetDirOffsetX(dir)
								ay := y + common.GetDirOffsetY(dir)
								ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(s, dir)
								as := chf.Spans[ai]
								if int32(as.Reg) != region {
									border = true
									break
								}
							}
						}
						if border {
							push3(queue, x, y, int32(i))
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
		for dir := int32(0); dir < 4; dir++ {
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

			ai := int32(chf.Cells[ax+ay*chf.Width].Index) + rcGetCon(cs, dir)
			as := chf.Spans[ai]

			hp.data[hx+hy*hp.width] = as.Y

			push3(queue, ax, ay, ai)
		}
	}
}

func getEdgeFlags(va, vb []float32,
	vpoly []float32, npoly int32) int32 {
	// The flag returned by this function matches dtDetailTriEdgeFlags in Detour.
	// Figure out if edge (va,vb) is part of the polygon boundary.
	thrSqr := float32(common.Sqr(0.001))
	i := int32(0)
	j := npoly - 1
	for i < npoly {

		if common.DistancePtSeg2d(va, common.GetVert3(vpoly, j), common.GetVert3(vpoly, i)) < thrSqr &&
			common.DistancePtSeg2d(vb, common.GetVert3(vpoly, j), common.GetVert3(vpoly, i)) < thrSqr {
			return 1
		}
		j = i
		i++
	}
	return 0
}

func getTriFlags(va, vb, vc []float32,
	vpoly []float32, npoly int32) int32 {
	flags := int32(0)
	flags |= getEdgeFlags(va, vb, vpoly, npoly) << 0
	flags |= getEdgeFlags(vb, vc, vpoly, npoly) << 2
	flags |= getEdgeFlags(vc, va, vpoly, npoly) << 4
	return flags
}

// / Contains triangle meshes that represent detailed height data associated
// / with the polygons in its associated polygon mesh object.
// / @ingroup recast
type RcPolyMeshDetail struct {
	Meshes  []uint32  ///< The sub-mesh data. [Size: 4*#nmeshes]
	Verts   []float32 ///< The mesh vertices. [Size: 3*#nverts]
	Tris    []uint8   ///< The mesh triangles. [Size: 4*#ntris]
	Nmeshes int32     ///< The number of sub-meshes defined by #meshes.
	Nverts  int32     ///< The number of vertices in #verts.
	Ntris   int32     ///< The number of triangles in #tris.
	// Explicitly-disabled copy constructor and copy assignment operator.
}

// / @par
// /
// / See the #RcConfig documentation for more information on the configuration parameters.
// /
// / @see rcAllocPolyMeshDetail, RcPolyMesh, RcCompactHeightfield, RcPolyMeshDetail, RcConfig
func RcBuildPolyMeshDetail(mesh *RcPolyMesh, chf *RcCompactHeightfield, sampleDist float32, sampleMaxError float32,
	dmesh *RcPolyMeshDetail) bool {
	if mesh.Nverts == 0 || mesh.Npolys == 0 {
		return true
	}

	nvp := mesh.Nvp
	cs := mesh.Cs
	ch := mesh.Ch
	orig := mesh.Bmin
	borderSize := mesh.BorderSize
	heightSearchRadius := int32(common.Max(1, math.Ceil(float64(mesh.MaxEdgeError))))
	c := func() int32 { return 0 }
	edges := NewStackArray(c, 64)
	tris := NewStackArray(func() int32 { return 0 }, 512)
	arr := NewStackArray(c, 512)
	samples := NewStackArray(c, 512)
	verts := make([]float32, 256*3)
	var hp rcHeightPatch
	nPolyVerts := 0
	maxhw := int32(0)
	maxhh := int32(0)
	bounds := make([]int32, mesh.Npolys*4)

	poly := make([]float32, nvp*3)
	// Find max size for a polygon area.
	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		xmin := &bounds[i*4+0]
		xmax := &bounds[i*4+1]
		ymin := &bounds[i*4+2]
		ymax := &bounds[i*4+3]
		*xmin = chf.Width
		*xmax = 0
		*ymin = chf.Height
		*ymax = 0
		for j := int32(0); j < nvp; j++ {
			if p[j] == RC_MESH_NULL_IDX {
				break
			}
			v := common.GetVert3(mesh.Verts, p[j])
			*xmin = min(*xmin, int32(v[0]))
			*xmax = max(*xmax, int32(v[0]))
			*ymin = min(*ymin, int32(v[2]))
			*ymax = max(*ymax, int32(v[2]))
			nPolyVerts++
		}
		*xmin = max(0, *xmin-1)
		*xmax = min(chf.Width, *xmax+1)
		*ymin = max(0, *ymin-1)
		*ymax = min(chf.Height, *ymax+1)
		if *xmin >= *xmax || *ymin >= *ymax {
			continue
		}
		maxhw = max(maxhw, *xmax-*xmin)
		maxhh = max(maxhh, *ymax-*ymin)
	}

	hp.data = make([]uint16, maxhw*maxhh)

	dmesh.Nmeshes = mesh.Npolys
	dmesh.Nverts = 0
	dmesh.Ntris = 0

	dmesh.Meshes = make([]uint32, dmesh.Nmeshes*4)

	vcap := nPolyVerts + nPolyVerts/2
	tcap := vcap * 2

	dmesh.Nverts = 0

	dmesh.Verts = make([]float32, vcap*3)
	dmesh.Ntris = 0

	dmesh.Tris = make([]uint8, tcap*4)

	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]

		// Store polygon vertices for processing.
		npoly := int32(0)
		for j := int32(0); j < nvp; j++ {
			if p[j] == RC_MESH_NULL_IDX {
				break
			}
			v := common.GetVert3(mesh.Verts, p[j])
			poly[j*3+0] = float32(v[0]) * cs
			poly[j*3+1] = float32(v[1]) * ch
			poly[j*3+2] = float32(v[2]) * cs
			npoly++
		}

		// Get the height data from the area of the polygon.
		hp.xmin = bounds[i*4+0]
		hp.ymin = bounds[i*4+2]
		hp.width = bounds[i*4+1] - bounds[i*4+0]
		hp.height = bounds[i*4+3] - bounds[i*4+2]
		getHeightData(chf, p, npoly, mesh.Verts, borderSize, &hp, arr, int32(mesh.Regs[i]))

		// Build detail mesh.
		nverts := int32(0)
		if !buildPolyDetail(poly, npoly,
			sampleDist, sampleMaxError,
			heightSearchRadius, chf, &hp,
			verts, &nverts, tris,
			edges, samples) {
			return false
		}

		// Move detail verts to world space.
		for j := int32(0); j < nverts; j++ {
			verts[j*3+0] += orig[0]
			verts[j*3+1] += orig[1] + chf.Ch // Is this offset necessary?
			verts[j*3+2] += orig[2]
		}
		// Offset poly too, will be used to flag checking.
		for j := int32(0); j < npoly; j++ {
			poly[j*3+0] += orig[0]
			poly[j*3+1] += orig[1]
			poly[j*3+2] += orig[2]
		}

		// Store detail submesh.
		ntris := tris.Len() / 4

		dmesh.Meshes[i*4+0] = uint32(dmesh.Nverts)
		dmesh.Meshes[i*4+1] = uint32(nverts)
		dmesh.Meshes[i*4+2] = uint32(dmesh.Ntris)
		dmesh.Meshes[i*4+3] = uint32(ntris)

		// Store vertices, allocate more memory if necessary.
		if int(dmesh.Nverts)+int(nverts) > vcap {
			for int(dmesh.Nverts)+int(nverts) > vcap {
				vcap += 256
			}

			newv := make([]float32, vcap*3)
			if dmesh.Nverts != 0 {
				copy(newv, dmesh.Verts[:3*dmesh.Nverts])
			}
			dmesh.Verts = newv
		}
		for j := int32(0); j < nverts; j++ {
			dmesh.Verts[dmesh.Nverts*3+0] = verts[j*3+0]
			dmesh.Verts[dmesh.Nverts*3+1] = verts[j*3+1]
			dmesh.Verts[dmesh.Nverts*3+2] = verts[j*3+2]
			dmesh.Nverts++
		}

		// Store triangles, allocate more memory if necessary.
		if int(dmesh.Ntris)+ntris > tcap {
			for int(dmesh.Ntris)+ntris > tcap {
				tcap += 256
			}

			newt := make([]uint8, tcap*4)

			if dmesh.Ntris != 0 {
				copy(newt, dmesh.Tris[:4*dmesh.Ntris])
			}
			dmesh.Tris = newt
		}
		for j := 0; j < ntris; j++ {
			t := tris.Slice(j*4, tris.Len())
			dmesh.Tris[dmesh.Ntris*4+0] = uint8(t[0])
			dmesh.Tris[dmesh.Ntris*4+1] = uint8(t[1])
			dmesh.Tris[dmesh.Ntris*4+2] = uint8(t[2])
			dmesh.Tris[dmesh.Ntris*4+3] = uint8(getTriFlags(common.GetVert3(verts, t[0]), common.GetVert3(verts, t[1]), common.GetVert3(verts, t[2]), poly, npoly))
			dmesh.Ntris++
		}
	}

	return true
}

// / @see rcAllocPolyMeshDetail, RcPolyMeshDetail
func rcMergePolyMeshDetails(meshes []*RcPolyMeshDetail, nmeshes int32, mesh *RcPolyMeshDetail) bool {

	maxVerts := int32(0)
	maxTris := int32(0)
	maxMeshes := int32(0)

	for i := int32(0); i < nmeshes; i++ {
		if meshes[i] == nil {
			continue
		}
		maxVerts += meshes[i].Nverts
		maxTris += meshes[i].Ntris
		maxMeshes += meshes[i].Nmeshes
	}

	mesh.Nmeshes = 0
	mesh.Meshes = make([]uint32, maxMeshes*4)
	mesh.Ntris = 0
	mesh.Tris = make([]uint8, maxTris*4)
	mesh.Nverts = 0

	mesh.Verts = make([]float32, maxVerts*3)
	// Merge datas.
	for i := int32(0); i < nmeshes; i++ {
		dm := meshes[i]
		if dm == nil {
			continue
		}
		for j := int32(0); j < dm.Nmeshes; j++ {
			dst := common.GetVert4(mesh.Meshes, mesh.Nmeshes)
			src := common.GetVert4(dm.Meshes, j)
			dst[0] = uint32(mesh.Nverts) + src[0]
			dst[1] = src[1]
			dst[2] = uint32(mesh.Ntris) + src[2]
			dst[3] = src[3]
			mesh.Nmeshes++
		}

		for k := int32(0); k < dm.Nverts; k++ {
			copy(common.GetVert3(mesh.Verts, mesh.Nverts), common.GetVert3(dm.Verts, k))
			mesh.Nverts++
		}
		for k := int32(0); k < dm.Ntris; k++ {
			mesh.Tris[mesh.Ntris*4+0] = dm.Tris[k*4+0]
			mesh.Tris[mesh.Ntris*4+1] = dm.Tris[k*4+1]
			mesh.Tris[mesh.Ntris*4+2] = dm.Tris[k*4+2]
			mesh.Tris[mesh.Ntris*4+3] = dm.Tris[k*4+3]
			mesh.Ntris++
		}
	}

	return true
}
