package recast

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
