package recast

import (
	"log"
	"math"
	"reflect"
)

const (
	VERTEX_BUCKET_COUNT = (1 << 12)
	RC_MESH_NULL_IDX    = 0xffff
)

type rcEdge struct {
	vert     [2]int
	polyEdge [2]int
	poly     [2]int
}

func buildMeshAdjacency(polys []int, npolys, nverts, vertsPerPoly int) bool {
	// Based on code by Eric Lengyel from:
	// https://web.archive.org/web/20080704083314/http://www.terathon.com/code/edges.php

	maxEdgeCount := npolys * vertsPerPoly
	firstEdge := make([]int, nverts+maxEdgeCount)
	nextEdge := firstEdge[nverts:]
	edgeCount := 0
	edges := make([]*rcEdge, maxEdgeCount)
	for i := range edges {
		edges[i] = &rcEdge{}
	}
	for i := 0; i < nverts; i++ {
		firstEdge[i] = RC_MESH_NULL_IDX
	}

	for i := 0; i < npolys; i++ {
		t := rcGetVert2(polys, i*vertsPerPoly)
		for j := 0; j < vertsPerPoly; j++ {
			if t[j] == RC_MESH_NULL_IDX {
				break
			}
			v0 := t[j]
			v1 := t[j+1]
			if j+1 >= vertsPerPoly || t[j+1] == RC_MESH_NULL_IDX {
				v1 = t[0]
			}
			if v0 < v1 {
				edge := edges[edgeCount]
				edge.vert[0] = v0
				edge.vert[1] = v1
				edge.poly[0] = i
				edge.polyEdge[0] = j
				edge.poly[1] = i
				edge.polyEdge[1] = 0
				// Insert edge
				nextEdge[edgeCount] = firstEdge[v0]
				firstEdge[v0] = edgeCount
				edgeCount++
			}
		}
	}

	for i := 0; i < npolys; i++ {
		t := rcGetVert2(polys, i*vertsPerPoly)
		for j := 0; j < vertsPerPoly; j++ {
			if t[j] == RC_MESH_NULL_IDX {
				break
			}
			v0 := t[j]
			v1 := t[j+1]
			if j+1 >= vertsPerPoly || t[j+1] == RC_MESH_NULL_IDX {
				v1 = t[0]
			}
			if v0 > v1 {
				for e := firstEdge[v1]; e != RC_MESH_NULL_IDX; e = nextEdge[e] {
					edge := edges[e]
					if edge.vert[1] == v0 && edge.poly[0] == edge.poly[1] {
						edge.poly[1] = i
						edge.polyEdge[1] = j
						break
					}
				}
			}
		}
	}

	// Store adjacency
	for i := 0; i < edgeCount; i++ {
		e := edges[i]
		if e.poly[0] != e.poly[1] {
			p0 := rcGetVert2(polys, e.poly[0]*vertsPerPoly)
			p1 := rcGetVert2(polys, e.poly[1]*vertsPerPoly)
			p0[vertsPerPoly+e.polyEdge[0]] = e.poly[1]
			p1[vertsPerPoly+e.polyEdge[1]] = e.poly[0]
		}
	}
	return true
}

func computeVertexHash(x, y, z int) int {
	h1 := 0x8da6b343 // Large multiplicative constants;
	h2 := 0xd8163841 // here arbitrarily chosen primes
	h3 := 0xcb1ab31f
	n := h1*x + h2*y + h3*z
	return (n & (VERTEX_BUCKET_COUNT - 1))
}

func addVertex(x int, y int, z int,
	verts []int, firstVert []int, nextVert []int) (nv int, t int) {
	bucket := computeVertexHash(x, 0, z)
	i := firstVert[bucket]

	for i != -1 {
		v := rcGetVert(verts, i)
		if v[0] == x && (rcAbs(v[1]-y) <= 2) && v[2] == z {
			return nv, i
		}

		i = nextVert[i] // next
	}

	// Could not find, create new.
	i = nv
	nv++
	v := rcGetVert(verts, i)
	v[0] = x
	v[1] = y
	v[2] = z
	nextVert[i] = firstVert[bucket]
	firstVert[bucket] = i

	return nv, i
}

// Last time I checked the if version got compiled using cmov, which was a lot faster than module (with idiv).
func prev(i, n int) int {
	if i-1 >= 0 {
		return i - 1
	}
	return n - 1
}
func next(i, n int) int {
	if i+1 < n {
		return i + 1
	}
	return 0
}

func area2(a, b, c []int) int {
	return (b[0]-a[0])*(c[2]-a[2]) - (c[0]-a[0])*(b[2]-a[2])
}

// Returns true iff c is strictly to the left of the directed
// line through a to b.
func left(a, b, c []int) bool {
	return area2(a, b, c) < 0
}

func leftOn(a, b, c []int) bool {
	return area2(a, b, c) <= 0
}

func collinear(a, b, c []int) bool {
	return area2(a, b, c) == 0
}

// Exclusive or: true iff exactly one argument is true.
// The arguments are negated to ensure that they are 0/1
// values.  Then the bitwise Xor operator may apply.
// (This idea is due to Michael Baldwin.)
func xorb(x, y bool) bool {
	if (x && !y) || (!x && y) {
		return true
	}
	return false
}

// Returns true iff ab properly intersects cd: they share
// a point interior to both segments.  The properness of the
// intersection is ensured by using strict leftness.
func intersectProp(a, b, c, d []int) bool {
	// Eliminate improper cases.
	if collinear(a, b, c) || collinear(a, b, d) ||
		collinear(c, d, a) || collinear(c, d, b) {
		return false
	}

	return xorb(left(a, b, c), left(a, b, d)) && xorb(left(c, d, a), left(c, d, b))
}

// Returns T iff (a,b,c) are collinear and point c lies
// on the closed segement ab.
func between(a, b, c []int) bool {
	if !collinear(a, b, c) {
		return false
	}

	// If ab not vertical, check betweenness on x; else on y.
	if a[0] != b[0] {
		return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]))
	}

	return ((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]))
}

// Returns true iff segments ab and cd intersect, properly or improperly.
func intersect(a, b, c, d []int) bool {
	if intersectProp(a, b, c, d) {
		return true
	}

	if between(a, b, c) || between(a, b, d) ||
		between(c, d, a) || between(c, d, b) {
		return true
	}

	return false
}

func vequal(a, b []int) bool {
	return a[0] == b[0] && a[2] == b[2]
}

// Returns T iff (v_i, v_j) is a proper internal *or* external
// diagonal of P, *ignoring edges incident to v_i and v_j*.
func diagonalie(i, j, n int, verts []int, indices []int) bool {
	d0 := rcGetVert4(verts, (indices[i] & 0x0fffffff))
	d1 := rcGetVert4(verts, (indices[j] & 0x0fffffff))

	// For each edge (k,k+1) of P
	for k := 0; k < n; k++ {
		k1 := next(k, n)
		// Skip edges incident to i or j
		if !((k == i) || (k1 == i) || (k == j) || (k1 == j)) {
			p0 := rcGetVert4(verts, (indices[k] & 0x0fffffff))
			p1 := rcGetVert4(verts, (indices[k1] & 0x0fffffff))

			if vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1) {
				continue
			}

			if intersect(d0, d1, p0, p1) {
				return false
			}

		}
	}
	return true
}

// Returns true iff the diagonal (i,j) is strictly internal to the
// polygon P in the neighborhood of the i endpoint.
func inCone(i, j, n int, verts []int, indices []int) bool {
	pi := rcGetVert4(verts, (indices[i] & 0x0fffffff))
	pj := rcGetVert4(verts, (indices[j] & 0x0fffffff))
	pi1 := rcGetVert4(verts, (indices[next(i, n)] & 0x0fffffff))
	pin1 := rcGetVert4(verts, (indices[prev(i, n)] & 0x0fffffff))

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if leftOn(pin1, pi, pi1) {
		return left(pi, pj, pin1) && left(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1))
}

// Returns T iff (v_i, v_j) is a proper internal
// diagonal of P.
func diagonal(i, j, n int, verts []int, indices []int) bool {
	return inCone(i, j, n, verts, indices) && diagonalie(i, j, n, verts, indices)
}

func diagonalieLoose(i, j, n int, verts []int, indices []int) bool {
	d0 := rcGetVert4(verts, (indices[i] & 0x0fffffff))
	d1 := rcGetVert4(verts, (indices[j] & 0x0fffffff))

	// For each edge (k,k+1) of P
	for k := 0; k < n; k++ {
		k1 := next(k, n)
		// Skip edges incident to i or j
		if !((k == i) || (k1 == i) || (k == j) || (k1 == j)) {
			p0 := rcGetVert(verts, (indices[k] & 0x0fffffff))
			p1 := rcGetVert(verts, (indices[k1] & 0x0fffffff))

			if vequal(d0, p0) || vequal(d1, p0) || vequal(d0, p1) || vequal(d1, p1) {
				continue
			}

			if intersectProp(d0, d1, p0, p1) {
				return false
			}

		}
	}
	return true
}

func inConeLoose(i, j, n int, verts []int, indices []int) bool {
	pi := rcGetVert4(verts, (indices[i] & 0x0fffffff))
	pj := rcGetVert4(verts, (indices[j] & 0x0fffffff))
	pi1 := rcGetVert4(verts, (indices[next(i, n)] & 0x0fffffff))
	pin1 := rcGetVert4(verts, (indices[prev(i, n)] & 0x0fffffff))

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if leftOn(pin1, pi, pi1) {
		return leftOn(pi, pj, pin1) && leftOn(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1))
}

func diagonalLoose(i, j, n int, verts []int, indices []int) bool {
	return inConeLoose(i, j, n, verts, indices) && diagonalieLoose(i, j, n, verts, indices)
}

func countPolyVerts(p []int, nvp int) int {
	for i := 0; i < nvp; i++ {
		if p[i] == RC_MESH_NULL_IDX {
			return i
		}
	}

	return nvp
}

func uleft(a, b, c []int) bool {
	return (b[0]-a[0])*(c[2]-a[2])-(c[0]-a[0])*(b[2]-a[2]) < 0
}

func getPolyMergeValue(pa, pb []int, verts []int, ea, eb int, nvp int) int {
	na := countPolyVerts(pa, nvp)
	nb := countPolyVerts(pb, nvp)

	// If the merged polygon would be too big, do not merge.
	if na+nb-2 > nvp {
		return -1
	}

	// Check if the polygons share an edge.
	ea = -1
	eb = -1

	for i := 0; i < na; i++ {
		va0 := pa[i]
		va1 := pa[(i+1)%na]
		if va0 > va1 {
			va0, va1 = va1, va0
		}

		for j := 0; j < nb; j++ {
			vb0 := pb[j]
			vb1 := pb[(j+1)%nb]
			if vb0 > vb1 {
				vb0, vb1 = vb1, vb0
			}

			if va0 == vb0 && va1 == vb1 {
				ea = i
				eb = j
				break
			}
		}
	}

	// No common edge, cannot merge.
	if ea == -1 || eb == -1 {
		return -1
	}

	// Check to see if the merged polygon would be convex.

	va := pa[(ea+na-1)%na]
	vb := pa[ea]
	vc := pb[(eb+2)%nb]

	if !uleft(rcGetVert(verts, va), rcGetVert(verts, vb), rcGetVert(verts, vc)) {
		return -1
	}

	va = pb[(eb+nb-1)%nb]
	vb = pb[eb]
	vc = pa[(ea+2)%na]
	if !uleft(rcGetVert(verts, va), rcGetVert(verts, vb), rcGetVert(verts, vc)) {
		return -1
	}

	va = pa[ea]
	vb = pa[(ea+1)%na]

	dx := verts[va*3+0] - verts[vb*3+0]
	dy := verts[va*3+2] - verts[vb*3+2]

	return dx*dx + dy*dy
}

func mergePolyVerts(pa, pb []int, ea, eb int, tmp []int, nvp int) {
	na := countPolyVerts(pa, nvp)
	nb := countPolyVerts(pb, nvp)

	for i := range tmp {
		if i < nvp {
			tmp[i] = 0xff
		}
	}
	// Merge polygons.
	n := 0
	// Add pa
	for i := 0; i < na-1; i++ {
		tmp[n] = pa[(ea+1+i)%na]
		n++
	}

	// Add pb
	for i := 0; i < nb-1; i++ {
		tmp[n] = pb[(eb+1+i)%nb]
		n++
	}

	copy(pa, tmp[:nvp])
}

func pushFront(v int, arr []int, an *int) {
	*an++
	for i := *an - 1; i > 0; i-- {
		arr[i] = arr[i-1]
	}
	arr[0] = v
}

func pushBack(v int, arr []int, an *int) {
	arr[*an] = v
	*an++
}

func canRemoveVertex(mesh *rcPolyMesh, rem int) bool {
	nvp := mesh.nvp

	// Count number of polygons to remove.
	numTouchedVerts := 0
	numRemainingEdges := 0
	for i := 0; i < mesh.npolys; i++ {
		p := mesh.polys[i*nvp*2 : i*nvp*2+2]
		nv := countPolyVerts(p, nvp)
		numRemoved := 0
		numVerts := 0
		for j := 0; j < nv; j++ {
			if p[j] == rem {
				numTouchedVerts++
				numRemoved++
			}
			numVerts++
		}
		if numRemoved > 0 {
			numRemainingEdges += numVerts - (numRemoved + 1)
		}
	}

	// There would be too few edges remaining to create a polygon.
	// This can happen for example when a tip of a triangle is marked
	// as deletion, but there are no other polys that share the vertex.
	// In this case, the vertex should not be removed.
	if numRemainingEdges <= 2 {
		return false
	}

	// Find edges which share the removed vertex.
	maxEdges := numTouchedVerts * 2
	nedges := 0
	edges := make([]int, maxEdges*3)
	for i := 0; i < mesh.npolys; i++ {
		p := mesh.polys[i*nvp*2 : i*nvp*2+2]
		nv := countPolyVerts(p, nvp)

		// Collect edges which touches the removed vertex.
		j := 0
		k := nv - 1
		for j < nv {
			if p[j] == rem || p[k] == rem {
				// Arrange edge so that a=rem.
				a := p[j]
				b := p[k]
				if b == rem {
					a, b = b, a
				}

				// Check if the edge exists
				exists := false
				for m := 0; m < nedges; m++ {
					e := edges[m*3 : m*3+3]
					if e[1] == b {
						// Exists, increment vertex share count.
						e[2]++
						exists = true
					}
				}
				// Add new edge.
				if !exists {
					e := edges[nedges*3 : nedges*3+3]
					e[0] = a
					e[1] = b
					e[2] = 1
					nedges++
				}
			}
		}
		k = j
		j++
	}

	// There should be no more than 2 open edges.
	// This catches the case that two non-adjacent polygons
	// share the removed vertex. In that case, do not remove the vertex.
	numOpenEdges := 0
	for i := 0; i < nedges; i++ {
		if edges[i*3+2] < 2 {
			numOpenEdges++
		}

	}
	if numOpenEdges > 2 {
		return false
	}

	return true
}

func triangulate(n int, verts, indices []int, tris []int) int {
	ntris := 0
	dst := 0

	// The last bit of the index is used to indicate if the vertex can be removed.
	for i := 0; i < n; i++ {
		i1 := next(i, n)
		i2 := next(i1, n)
		if diagonal(i, i2, n, verts, indices) {
			indices[i1] |= 0x80000000
		}

	}

	for n > 3 {
		minLen := -1
		mini := -1
		for i := 0; i < n; i++ {
			i1 := next(i, n)
			if indices[i1]&0x80000000 > 0 {
				p0 := rcGetVert(verts, (indices[i] & 0x0fffffff))
				p2 := rcGetVert(verts, (indices[next(i1, n)] & 0x0fffffff))

				dx := p2[0] - p0[0]
				dy := p2[2] - p0[2]
				length := dx*dx + dy*dy

				if minLen < 0 || length < minLen {
					minLen = length
					mini = i
				}
			}
		}

		if mini == -1 {
			// We might get here because the contour has overlapping segments, like this:
			//
			//  A o-o=====o---o B
			//   /  |C   D|    \.
			//  o   o     o     o
			//  :   :     :     :
			// We'll try to recover by loosing up the inCone test a bit so that a diagonal
			// like A-B or C-D can be found and we can continue.
			minLen = -1
			mini = -1
			for i := 0; i < n; i++ {
				i1 := next(i, n)
				i2 := next(i1, n)
				if diagonalLoose(i, i2, n, verts, indices) {
					p0 := rcGetVert4(verts, (indices[i] & 0x0fffffff))
					p2 := rcGetVert4(verts, (indices[next(i2, n)] & 0x0fffffff))
					dx := p2[0] - p0[0]
					dy := p2[2] - p0[2]
					length := dx*dx + dy*dy

					if minLen < 0 || length < minLen {
						minLen = length
						mini = i
					}
				}
			}
			if mini == -1 {
				// The contour is messed up. This sometimes happens
				// if the contour simplification is too aggressive.
				return -ntris
			}
		}

		i := mini
		i1 := next(i, n)
		i2 := next(i1, n)

		tris[dst] = indices[i] & 0x0fffffff
		dst++
		tris[dst] = indices[i1] & 0x0fffffff
		dst++
		tris[dst] = indices[i2] & 0x0fffffff
		dst++
		ntris++

		// Removes P[i1] by copying P[i+1]...P[n-1] left one index.
		n--
		for k := i1; k < n; k++ {
			indices[k] = indices[k+1]
		}

		if i1 >= n {
			i1 = 0
		}
		i = prev(i1, n)

		// Update diagonal flags.
		if diagonal(prev(i, n), i1, n, verts, indices) {
			indices[i] |= 0x80000000
		} else {
			indices[i] &= 0x0fffffff
		}

		if diagonal(i, next(i1, n), n, verts, indices) {
			indices[i1] |= 0x80000000
		} else {
			indices[i1] &= 0x0fffffff
		}

	}

	// Append the remaining triangle.
	tris[dst] = indices[0] & 0x0fffffff
	dst++
	tris[dst] = indices[1] & 0x0fffffff
	dst++
	tris[dst] = indices[2] & 0x0fffffff
	dst++
	ntris++

	return ntris
}
func removeVertex(mesh *rcPolyMesh, rem int, maxTris int) bool {
	nvp := mesh.nvp

	// Count number of polygons to remove.
	numRemovedVerts := 0
	for i := 0; i < mesh.npolys; i++ {
		p := mesh.polys[i*nvp*2:]
		nv := countPolyVerts(p, nvp)
		for j := 0; j < nv; j++ {
			if p[j] == rem {
				numRemovedVerts++
			}

		}
	}

	nedges := 0
	edges := make([]int, numRemovedVerts*nvp*4)

	nhole := 0
	hole := make([]int, numRemovedVerts*nvp)
	nhreg := 0
	hreg := make([]int, numRemovedVerts*nvp)
	nharea := 0
	harea := make([]int, numRemovedVerts*nvp)

	for i := 0; i < mesh.npolys; i++ {
		p := mesh.polys[i*nvp*2:]
		nv := countPolyVerts(p, nvp)
		hasRem := false
		for j := 0; j < nv; j++ {
			if p[j] == rem {
				hasRem = true
			}
		}

		if hasRem {
			// Collect edges which does not touch the removed vertex.
			j := 0
			k := nv - 1
			for j < nv {
				if p[j] != rem && p[k] != rem {
					e := rcGetVert4(edges, nedges)
					e[0] = p[k]
					e[1] = p[j]
					e[2] = mesh.regs[i]
					e[3] = mesh.areas[i]
					nedges++
				}
				k = j
				j++
			}
			// Remove the polygon.
			p2 := mesh.polys[(mesh.npolys-1)*nvp*2:]
			if !reflect.DeepEqual(p, p2) {
				copy(p, p2[:nvp])
			}
			for k := nvp; i < nvp+nvp; i++ {
				p[k] = 0xff
			}
			mesh.regs[i] = mesh.regs[mesh.npolys-1]
			mesh.areas[i] = mesh.areas[mesh.npolys-1]
			mesh.npolys--
			i--
		}
	}

	// Remove vertex.
	for i := rem; i < mesh.nverts-1; i++ {
		mesh.verts[i*3+0] = mesh.verts[(i+1)*3+0]
		mesh.verts[i*3+1] = mesh.verts[(i+1)*3+1]
		mesh.verts[i*3+2] = mesh.verts[(i+1)*3+2]
	}
	mesh.nverts--

	// Adjust indices to match the removed vertex layout.
	for i := 0; i < mesh.npolys; i++ {
		p := mesh.polys[i*nvp*2:]
		nv := countPolyVerts(p, nvp)
		for j := 0; j < nv; j++ {
			if p[j] > rem {
				p[j]--
			}
		}

	}
	for i := 0; i < nedges; i++ {
		if edges[i*4+0] > rem {
			edges[i*4+0]--
		}
		if edges[i*4+1] > rem {
			edges[i*4+1]--
		}
	}

	if nedges == 0 {
		return true
	}

	// Start with one vertex, keep appending connected
	// segments to the start and end of the hole.
	pushBack(edges[0], hole, &nhole)
	pushBack(edges[2], hreg, &nhreg)
	pushBack(edges[3], harea, &nharea)

	for nedges > 0 {
		match := false

		for i := 0; i < nedges; i++ {
			ea := edges[i*4+0]
			eb := edges[i*4+1]
			r := edges[i*4+2]
			a := edges[i*4+3]
			add := false
			if hole[0] == eb {
				// The segment matches the beginning of the hole boundary.
				pushFront(ea, hole, &nhole)
				pushFront(r, hreg, &nhreg)
				pushFront(a, harea, &nharea)
				add = true
			} else if hole[nhole-1] == ea {
				// The segment matches the end of the hole boundary.
				pushBack(eb, hole, &nhole)
				pushBack(r, hreg, &nhreg)
				pushBack(a, harea, &nharea)
				add = true
			}
			if add {
				// The edge segment was added, remove it.
				edges[i*4+0] = edges[(nedges-1)*4+0]
				edges[i*4+1] = edges[(nedges-1)*4+1]
				edges[i*4+2] = edges[(nedges-1)*4+2]
				edges[i*4+3] = edges[(nedges-1)*4+3]
				nedges--
				match = true
				i--
			}
		}

		if !match {
			break
		}

	}

	tris := make([]int, nhole*3)
	tverts := make([]int, nhole*4)
	thole := make([]int, nhole)
	// Generate temp vertex array for triangulation.
	for i := 0; i < nhole; i++ {
		pi := hole[i]
		tverts[i*4+0] = mesh.verts[pi*3+0]
		tverts[i*4+1] = mesh.verts[pi*3+1]
		tverts[i*4+2] = mesh.verts[pi*3+2]
		tverts[i*4+3] = 0
		thole[i] = i
	}

	// Triangulate the hole.
	ntris := triangulate(nhole, tverts, thole, tris)
	if ntris < 0 {
		ntris = -ntris
		log.Printf("removeVertex: triangulate() returned bad results.")
	}

	// Merge the hole triangles back to polygons.
	polys := make([]int, (ntris+1)*nvp)
	for i := range polys {
		polys[i] = 0xff
	}
	pregs := make([]int, ntris)
	pareas := make([]int, ntris)
	tmpPoly := polys[ntris*nvp:]
	// Build initial polygons.
	npolys := 0
	for j := 0; j < ntris; j++ {
		t := rcGetVert(tris, j)
		if t[0] != t[1] && t[0] != t[2] && t[1] != t[2] {
			polys[npolys*nvp+0] = hole[t[0]]
			polys[npolys*nvp+1] = hole[t[1]]
			polys[npolys*nvp+2] = hole[t[2]]

			// If this polygon covers multiple region types then
			// mark it as such
			if hreg[t[0]] != hreg[t[1]] || hreg[t[1]] != hreg[t[2]] {
				pregs[npolys] = RC_MULTIPLE_REGS
			} else {
				pregs[npolys] = hreg[t[0]]
			}

			pareas[npolys] = harea[t[0]]
			npolys++
		}
	}
	if npolys == 0 {
		return true
	}

	// Merge polygons.
	if nvp > 3 {
		for {
			// Find best polygons to merge.
			var bestMergeVal = 0
			var bestPa = 0
			var bestPb = 0
			var bestEa = 0
			var bestEb = 0

			for j := 0; j < npolys-1; j++ {
				pj := polys[j*nvp:]
				for k := j + 1; k < npolys; k++ {
					pk := polys[k*nvp:]
					var ea, eb int
					v := getPolyMergeValue(pj, pk, mesh.verts, ea, eb, nvp)
					if v > bestMergeVal {
						bestMergeVal = v
						bestPa = j
						bestPb = k
						bestEa = ea
						bestEb = eb
					}
				}
			}

			if bestMergeVal > 0 {
				// Found best, merge.
				pa := polys[bestPa*nvp:]
				pb := polys[bestPb*nvp:]
				mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp)
				if pregs[bestPa] != pregs[bestPb] {
					pregs[bestPa] = RC_MULTIPLE_REGS
				}

				last := polys[(npolys-1)*nvp:]
				if !(reflect.DeepEqual(pb, last)) {
					copy(pb, last[:nvp])
				}

				pregs[bestPb] = pregs[npolys-1]
				pareas[bestPb] = pareas[npolys-1]
				npolys--
			} else {
				// Could not merge any polygons, stop.
				break
			}
		}
	}

	// Store polygons.
	for i := 0; i < npolys; i++ {
		if mesh.npolys >= maxTris {
			break
		}
		p := mesh.polys[mesh.npolys*nvp*2:]
		for i := 0; i < nvp*2; i++ {
			p[i] = 0xff
		}
		for j := 0; j < nvp; j++ {
			p[j] = polys[i*nvp+j]
		}
		mesh.regs[mesh.npolys] = pregs[i]
		mesh.areas[mesh.npolys] = pareas[i]
		mesh.npolys++
		if mesh.npolys > maxTris {
			log.Printf("removeVertex: Too many polygons %d (max:%d).", mesh.npolys, maxTris)
			return false
		}
	}

	return true
}

// / @par
// /
// / @note If the mesh data is to be used to construct a Detour navigation mesh, then the upper
// / limit must be restricted to <= #DT_VERTS_PER_POLYGON.
// /
// / @see rcAllocPolyMesh, rcContourSet, rcPolyMesh, rcConfig
func rcBuildPolyMesh(cset *rcContourSet, nvp int, mesh *rcPolyMesh) bool {
	copy(mesh.bmin, cset.bmin)
	copy(mesh.bmax, cset.bmax)
	mesh.cs = cset.cs
	mesh.ch = cset.ch
	mesh.borderSize = cset.borderSize
	mesh.maxEdgeError = cset.maxError

	maxVertices := 0
	maxTris := 0
	maxVertsPerCont := 0
	for i := 0; i < cset.nconts; i++ {
		// Skip null contours.
		if cset.conts[i].nverts < 3 {
			continue
		}
		maxVertices += cset.conts[i].nverts
		maxTris += cset.conts[i].nverts - 2
		maxVertsPerCont = rcMax(maxVertsPerCont, cset.conts[i].nverts)
	}

	if maxVertices >= 0xfffe {
		log.Printf("rcBuildPolyMesh: Too many vertices %d.", maxVertices)
		return false
	}

	vflags := make([]int, maxVertices)
	mesh.verts = make([]int, maxVertices*3)

	mesh.polys = make([]int, maxTris*nvp*2)
	for i := range mesh.polys {
		mesh.polys[i] = 0xff
	}
	mesh.regs = make([]int, maxTris)
	mesh.areas = make([]int, maxTris)

	mesh.nverts = 0
	mesh.npolys = 0
	mesh.nvp = nvp
	mesh.maxpolys = maxTris
	nextVert := make([]int, maxVertices)
	firstVert := make([]int, VERTEX_BUCKET_COUNT)
	for i := 0; i < VERTEX_BUCKET_COUNT; i++ {
		firstVert[i] = -1
	}

	indices := make([]int, maxVertsPerCont)
	tris := make([]int, maxVertsPerCont*3)
	polys := make([]int, (maxVertsPerCont+1)*nvp)
	tmpPoly := polys[maxVertsPerCont*nvp:]

	for i := 0; i < cset.nconts; i++ {
		cont := cset.conts[i]

		// Skip null contours.
		if cont.nverts < 3 {
			continue
		}

		// Triangulate contour
		for j := 0; j < cont.nverts; j++ {
			indices[j] = j
		}

		ntris := triangulate(cont.nverts, cont.verts, indices, tris)
		if ntris <= 0 {
			// Bad triangulation, should not happen.
			/*			printf("\tconst float bmin[3] = {%ff,%ff,%ff};\n", cset.bmin[0], cset.bmin[1], cset.bmin[2]);
						printf("\tconst float cs = %ff;\n", cset.cs);
						printf("\tconst float ch = %ff;\n", cset.ch);
						printf("\tconst int verts[] = {\n");
						for (int k = 0; k < cont.nverts; ++k)
						{
							const int* v = &cont.verts[k*4];
							printf("\t\t%d,%d,%d,%d,\n", v[0], v[1], v[2], v[3]);
						}
						printf("\t};\n\tconst int nverts = sizeof(verts)/(sizeof(int)*4);\n");*/
			log.Printf("rcBuildPolyMesh: Bad triangulation Contour %d.", i)
			ntris = -ntris
		}

		// Add and merge vertices.
		for j := 0; j < cont.nverts; j++ {
			v := rcGetVert4(cont.verts, j)
			mesh.nverts, indices[j] = addVertex(v[0], v[1], v[2], mesh.verts, firstVert, nextVert)
			if v[3]&RC_BORDER_VERTEX > 0 {
				// This vertex should be removed.
				vflags[indices[j]] = 1
			}
		}

		// Build initial polygons.
		npolys := 0
		for i := range polys {
			polys[i] = 0xff
		}
		for j := 0; j < ntris; j++ {
			t := rcGetVert(tris, j)
			if t[0] != t[1] && t[0] != t[2] && t[1] != t[2] {
				polys[npolys*nvp+0] = indices[t[0]]
				polys[npolys*nvp+1] = indices[t[1]]
				polys[npolys*nvp+2] = indices[t[2]]
				npolys++
			}
		}
		if npolys == 0 {
			continue
		}

		// Merge polygons.
		if nvp > 3 {
			for {
				// Find best polygons to merge.
				bestMergeVal := 0
				bestPa := 0
				bestPb := 0
				bestEa := 0
				bestEb := 0

				for j := 0; j < npolys-1; j++ {
					pj := polys[j*nvp:]
					for k := j + 1; k < npolys; k++ {
						pk := polys[k*nvp:]
						var ea, eb int
						v := getPolyMergeValue(pj, pk, mesh.verts, ea, eb, nvp)
						if v > bestMergeVal {
							bestMergeVal = v
							bestPa = j
							bestPb = k
							bestEa = ea
							bestEb = eb
						}
					}
				}

				if bestMergeVal > 0 {
					// Found best, merge.
					pa := polys[bestPa*nvp:]
					pb := polys[bestPb*nvp:]
					mergePolyVerts(pa, pb, bestEa, bestEb, tmpPoly, nvp)
					lastPoly := polys[(npolys-1)*nvp:]
					if !reflect.DeepEqual(pb, lastPoly) {
						copy(pb, lastPoly[:nvp])
					}

					npolys--
				} else {
					// Could not merge any polygons, stop.
					break
				}
			}
		}

		// Store polygons.
		for j := 0; j < npolys; j++ {
			p := mesh.polys[mesh.npolys*nvp*2:]
			q := polys[j*nvp:]
			for k := 0; k < nvp; k++ {
				p[k] = q[k]
			}

			mesh.regs[mesh.npolys] = cont.reg
			mesh.areas[mesh.npolys] = cont.area
			mesh.npolys++
			if mesh.npolys > maxTris {
				log.Printf("rcBuildPolyMesh: Too many polygons %d (max:%d).", mesh.npolys, maxTris)
				return false
			}
		}
	}

	// Remove edge vertices.
	for i := 0; i < mesh.nverts; i++ {
		if vflags[i] != 0 {
			if !canRemoveVertex(mesh, i) {
				continue
			}

			if !removeVertex(mesh, i, maxTris) {
				// Failed to remove vertex
				log.Printf("rcBuildPolyMesh: Failed to remove edge vertex %d.", i)
				return false
			}
			// Remove vertex
			// Note: mesh.nverts is already decremented inside removeVertex()!
			// Fixup vertex flags
			for j := i; j < mesh.nverts; j++ {
				vflags[j] = vflags[j+1]
			}

			i--
		}
	}

	// Calculate adjacency.
	if !buildMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, nvp) {
		log.Printf("rcBuildPolyMesh: Adjacency failed.")
		return false
	}

	// Find portal edges
	if mesh.borderSize > 0 {
		w := cset.width
		h := cset.height
		for i := 0; i < mesh.npolys; i++ {
			p := mesh.polys[i*2*nvp:]
			for j := 0; j < nvp; j++ {
				if p[j] == RC_MESH_NULL_IDX {
					break
				}
				// Skip connected edges.
				if p[nvp+j] != RC_MESH_NULL_IDX {
					continue
				}

				nj := j + 1
				if nj >= nvp || p[nj] == RC_MESH_NULL_IDX {
					nj = 0
				}
				va := rcGetVert(mesh.verts, p[j])
				vb := rcGetVert(mesh.verts, p[nj])

				if va[0] == 0 && vb[0] == 0 {
					p[nvp+j] = 0x8000 | 0
				} else if va[2] == h && vb[2] == h {
					p[nvp+j] = 0x8000 | 1
				} else if va[0] == w && vb[0] == w {
					p[nvp+j] = 0x8000 | 2
				} else if va[2] == 0 && vb[2] == 0 {
					p[nvp+j] = 0x8000 | 3
				}

			}
		}
	}

	// Just allocate the mesh flags array. The user is resposible to fill it.
	mesh.flags = make([]int, mesh.npolys)
	if mesh.nverts > 0xffff {
		log.Printf("rcBuildPolyMesh: The resulting mesh has too many vertices %d (max %d). Data can be corrupted.", mesh.nverts, 0xffff)
	}
	if mesh.npolys > 0xffff {
		log.Printf("rcBuildPolyMesh: The resulting mesh has too many polygons %d (max %d). Data can be corrupted.", mesh.npolys, 0xffff)
	}

	return true
}

// / @see rcAllocPolyMesh, rcPolyMesh
func rcMergePolyMeshes(meshes []*rcPolyMesh, nmeshes int, mesh *rcPolyMesh) bool {

	if nmeshes == 0 || len(meshes) == 0 {
		return true
	}
	mesh.nvp = meshes[0].nvp
	mesh.cs = meshes[0].cs
	mesh.ch = meshes[0].ch
	copy(mesh.bmin, meshes[0].bmin)
	copy(mesh.bmax, meshes[0].bmax)

	maxVerts := 0
	maxPolys := 0
	maxVertsPerMesh := 0
	for i := 0; i < nmeshes; i++ {
		mesh.bmin = rcVmin(mesh.bmin, meshes[i].bmin)
		mesh.bmax = rcVmax(mesh.bmax, meshes[i].bmax)
		maxVertsPerMesh = rcMax(maxVertsPerMesh, meshes[i].nverts)
		maxVerts += meshes[i].nverts
		maxPolys += meshes[i].npolys
	}

	mesh.nverts = 0
	mesh.verts = make([]int, maxVerts*3)

	mesh.npolys = 0
	mesh.polys = make([]int, maxPolys*2*mesh.nvp)
	for i := range mesh.polys {
		mesh.polys[i] = 0xff
	}

	mesh.regs = make([]int, maxPolys)

	mesh.areas = make([]int, maxPolys)

	mesh.flags = make([]int, maxPolys)

	nextVert := make([]int, maxVerts)

	firstVert := make([]int, VERTEX_BUCKET_COUNT)
	for i := 0; i < VERTEX_BUCKET_COUNT; i++ {
		firstVert[i] = -1
	}

	vremap := make([]int, maxVertsPerMesh)
	for i := 0; i < nmeshes; i++ {
		pmesh := meshes[i]

		ox := int(math.Floor((pmesh.bmin[0]-mesh.bmin[0])/mesh.cs + 0.5))
		oz := int(math.Floor((pmesh.bmin[2]-mesh.bmin[2])/mesh.cs + 0.5))

		isMinX := (ox == 0)
		isMinZ := (oz == 0)
		isMaxX := (math.Floor((mesh.bmax[0]-pmesh.bmax[0])/mesh.cs + 0.5)) == 0
		isMaxZ := (math.Floor((mesh.bmax[2]-pmesh.bmax[2])/mesh.cs + 0.5)) == 0
		isOnBorder := (isMinX || isMinZ || isMaxX || isMaxZ)

		for j := 0; j < pmesh.nverts; j++ {
			v := rcGetVert(pmesh.verts, j)
			mesh.nverts, vremap[j] = addVertex(v[0]+ox, v[1], v[2]+oz, mesh.verts, firstVert, nextVert)
		}

		for j := 0; j < pmesh.npolys; j++ {
			tgt := mesh.polys[mesh.npolys*2*mesh.nvp:]
			src := pmesh.polys[j*2*mesh.nvp:]
			mesh.regs[mesh.npolys] = pmesh.regs[j]
			mesh.areas[mesh.npolys] = pmesh.areas[j]
			mesh.flags[mesh.npolys] = pmesh.flags[j]
			mesh.npolys++
			for k := 0; k < mesh.nvp; k++ {
				if src[k] == RC_MESH_NULL_IDX {
					break
				}
				tgt[k] = vremap[src[k]]
			}

			if isOnBorder {
				for k := mesh.nvp; k < mesh.nvp*2; k++ {
					if src[k]&0x8000 > 0 && src[k] != 0xffff {
						dir := src[k] & 0xf
						switch dir {
						case 0: // Portal x-
							if isMinX {
								tgt[k] = src[k]
							}
							break
						case 1: // Portal z+
							if isMaxZ {
								tgt[k] = src[k]
							}
							break
						case 2: // Portal x+
							if isMaxX {
								tgt[k] = src[k]
							}

							break
						case 3: // Portal z-
							if isMinZ {
								tgt[k] = src[k]
							}

							break
						}
					}
				}
			}
		}
	}

	// Calculate adjacency.
	if !buildMeshAdjacency(mesh.polys, mesh.npolys, mesh.nverts, mesh.nvp) {
		log.Printf("rcMergePolyMeshes: Adjacency failed.")
		return false
	}

	if mesh.nverts > 0xffff {
		log.Printf("rcMergePolyMeshes: The resulting mesh has too many vertices %d (max %d). Data can be corrupted.", mesh.nverts, 0xffff)
	}
	if mesh.npolys > 0xffff {
		log.Printf("rcMergePolyMeshes: The resulting mesh has too many polygons %d (max %d). Data can be corrupted.", mesh.npolys, 0xffff)
	}

	return true
}

func rcCopyPolyMesh(src *rcPolyMesh, dst *rcPolyMesh) bool {
	// Destination must be empty.
	if len(dst.verts) == 0 {
		panic("")
	}
	if len(dst.polys) == 0 {
		panic("")
	}
	if len(dst.regs) == 0 {
		panic("")
	}
	if len(dst.areas) == 0 {
		panic("")
	}
	if len(dst.flags) == 0 {
		panic("")
	}

	dst.nverts = src.nverts
	dst.npolys = src.npolys
	dst.maxpolys = src.npolys
	dst.nvp = src.nvp
	copy(dst.bmin, src.bmin)
	copy(dst.bmax, src.bmax)
	dst.cs = src.cs
	dst.ch = src.ch
	dst.borderSize = src.borderSize
	dst.maxEdgeError = src.maxEdgeError

	dst.verts = make([]int, src.nverts*3)
	copy(dst.verts, src.verts[:src.nverts*3])

	dst.polys = make([]int, src.npolys*2*src.nvp)
	copy(dst.polys, src.polys[:src.npolys*2*src.nvp])

	dst.regs = make([]int, src.npolys)
	copy(dst.regs, src.regs[:src.npolys])

	dst.areas = make([]int, src.npolys)
	copy(dst.areas, src.areas[:src.npolys])

	dst.flags = make([]int, src.npolys)

	copy(dst.flags, src.flags[:src.npolys])

	return true
}
