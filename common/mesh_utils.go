package common

import "math"

// Last time I checked the if version got compiled using cmov, which was a lot faster than module (with idiv).
func Prev[T IT](i, n T) T {
	if i-1 >= 0 {
		return i - 1
	}
	return n - 1
}
func Next[T IT](i, n T) T {
	if i+1 < n {
		return i + 1
	}
	return 0
}

func Area2[T IT](a, b, c []T) T {
	return (b[0]-a[0])*(c[2]-a[2]) - (c[0]-a[0])*(b[2]-a[2])
}

// Returns true iff c is strictly to the left of the directed
// line through a to b.
func Left[T IT](a, b, c []T) bool {
	return Area2(a, b, c) < 0
}

func LeftOn[T IT](a, b, c []T) bool {
	return Area2(a, b, c) <= 0
}

func Collinear[T IT](a, b, c []T) bool {
	return Area2(a, b, c) == 0
}

// Exclusive or: true iff exactly one argument is true.
// The arguments are negated to ensure that they are 0/1
// values.  Then the bitwise Xor operator may apply.
// (This idea is due to Michael Baldwin.)
func Xorb(x, y bool) bool {
	if (x && !y) || (!x && y) {
		return true
	}
	return false
}

// Returns true iff ab properly intersects cd: they share
// a point interior to both segments.  The properness of the
// intersection is ensured by using strict leftness.
func IntersectProp[T IT](a, b, c, d []T) bool {
	// Eliminate improper cases.
	if Collinear(a, b, c) || Collinear(a, b, d) ||
		Collinear(c, d, a) || Collinear(c, d, b) {
		return false
	}

	return Xorb(Left(a, b, c), Left(a, b, d)) && Xorb(Left(c, d, a), Left(c, d, b))
}

// Returns T iff (a,b,c) are collinear and point c lies
// on the closed segement ab.
func Between[T IT](a, b, c []T) bool {
	if !Collinear(a, b, c) {
		return false
	}

	// If ab not vertical, check betweenness on x; else on y.
	if a[0] != b[0] {
		return ((a[0] <= c[0]) && (c[0] <= b[0])) || ((a[0] >= c[0]) && (c[0] >= b[0]))
	}

	return ((a[2] <= c[2]) && (c[2] <= b[2])) || ((a[2] >= c[2]) && (c[2] >= b[2]))
}

// Returns true iff segments ab and cd intersect, properly or improperly.
func Intersect[T IT](a, b, c, d []T) bool {
	if IntersectProp(a, b, c, d) {
		return true
	}

	if Between(a, b, c) || Between(a, b, d) ||
		Between(c, d, a) || Between(c, d, b) {
		return true
	}

	return false
}

func RcVequal[T IT](a, b []T) bool {
	return a[0] == b[0] && a[2] == b[2]
}

// Returns T iff (v_i, v_j) is a proper internal *or* external
// diagonal of P, *ignoring edges incident to v_i and v_j*.
func Diagonalie[T1 IT, T2 int32 | uint8 | int, T3 int32 | uint16 | int](i, j, n T1, verts []T2, indices []T3) bool {
	d0 := GetVert4(verts, int(indices[int64(i)])&0x0fffffff)
	d1 := GetVert4(verts, int(indices[int64(j)])&0x0fffffff)

	// For each edge (k,k+1) of P
	for k := int64(0); k < int64(n); k++ {
		k1 := Next(k, int64(n))
		// Skip edges incident to i or j
		if !((k == int64(i)) || (k1 == int64(i)) || (k == int64(j)) || (k1 == int64(j))) {
			p0 := GetVert4(verts, int(indices[k])&0x0fffffff)
			p1 := GetVert4(verts, int(indices[k1])&0x0fffffff)

			if RcVequal(d0, p0) || RcVequal(d1, p0) || RcVequal(d0, p1) || RcVequal(d1, p1) {
				continue
			}

			if Intersect(d0, d1, p0, p1) {
				return false
			}

		}
	}
	return true
}

func ContourInCone(i, n int32, verts, pj []int32) bool {
	pi := GetVert4(verts, i)
	pi1 := GetVert4(verts, Next(i, n))
	pin1 := GetVert4(verts, Prev(i, n))

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if LeftOn(pin1, pi, pi1) {
		return Left(pi, pj, pin1) && Left(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(LeftOn(pi, pj, pi1) && LeftOn(pj, pi, pin1))
}

// Returns true iff the diagonal (i,j) is strictly internal to the
// polygon P in the neighborhood of the i endpoint.
func InCone[T1 IT, T2 int32 | uint8 | int, T3 int32 | uint16 | int](i, j, n T1, verts []T2, indices []T3) bool {
	pi := GetVert4(verts, int(indices[int64(i)])&0x0fffffff)
	pj := GetVert4(verts, int(indices[int64(j)])&0x0fffffff)
	pi1 := GetVert4(verts, int(indices[int64(Next(i, n))])&0x0fffffff)
	pin1 := GetVert4(verts, int(indices[int64(Prev(i, n))])&0x0fffffff)

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if LeftOn(pin1, pi, pi1) {
		return Left(pi, pj, pin1) && Left(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(LeftOn(pi, pj, pi1) && LeftOn(pj, pi, pin1))
}

// Returns T iff (v_i, v_j) is a proper internal
// diagonal of P.
func Diagonal[T1 IT, T2 int32 | uint8 | int, T3 int32 | uint16 | int](i, j, n T1, verts []T2, indices []T3) bool {
	return InCone(i, j, n, verts, indices) && Diagonalie(i, j, n, verts, indices)
}

func DiagonalieLoose[T1 IT, T2 int32 | uint8 | int, T3 int32 | uint16 | int](i, j, n T1, verts []T2, indices []T3) bool {
	tmp := 0x0fffffff
	d0 := GetVert4(verts, indices[int64(i)]&T3(tmp))
	d1 := GetVert4(verts, indices[int64(j)]&T3(tmp))

	// For each edge (k,k+1) of P
	for k := int64(0); k < int64(n); k++ {
		k1 := Next(k, int64(n))
		// Skip edges incident to i or j
		if !((k == int64(i)) || (k1 == int64(i)) || (k == int64(j)) || (k1 == int64(j))) {
			p0 := GetVert3(verts, indices[k]&T3(tmp))
			p1 := GetVert3(verts, indices[k1]&T3(tmp))

			if RcVequal(d0, p0) || RcVequal(d1, p0) || RcVequal(d0, p1) || RcVequal(d1, p1) {
				continue
			}

			if IntersectProp(d0, d1, p0, p1) {
				return false
			}

		}
	}
	return true
}

func InConeLoose[T1 IT, T2 int32 | uint8 | int, T3 int32 | uint16 | int](i, j, n T1, verts []T2, indices []T3) bool {
	tmp := 0x0fffffff
	pi := GetVert4(verts, indices[int64(i)]&T3(tmp))
	pj := GetVert4(verts, indices[int64(j)]&T3(tmp))
	pi1 := GetVert4(verts, indices[int64(Next(i, n))]&T3(tmp))
	pin1 := GetVert4(verts, indices[int64(Prev(i, n))]&T3(tmp))

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if LeftOn(pin1, pi, pi1) {
		return LeftOn(pi, pj, pin1) && LeftOn(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(LeftOn(pi, pj, pi1) && LeftOn(pj, pi, pin1))
}

func DiagonalLoose[T1 IT, T2 int32 | uint8 | int, T3 int32 | uint16 | int](i, j, n T1, verts []T2, indices []T3) bool {
	return InConeLoose(i, j, n, verts, indices) && DiagonalieLoose(i, j, n, verts, indices)
}

func Triangulate[T1 int32 | uint8 | int, T2 int32 | uint16 | int](n int32, verts []T1, indices []T2, tris []T2) int32 {
	ntris := int32(0)
	dst := 0
	tmpT1 := 0x0fffffff
	tmpT2 := 0x80000000
	// The last bit of the index is used to indicate if the vertex can be removed.
	for i := int32(0); i < n; i++ {
		i1 := Next(i, n)
		i2 := Next(i1, n)
		if Diagonal[int32, T1, T2](i, i2, n, verts, indices) {
			indices[i1] |= T2(tmpT2)
		}

	}

	for n > 3 {
		minLen := -1
		mini := -1
		for i := int32(0); i < n; i++ {
			i1 := Next(i, n)
			if indices[i1]&T2(tmpT2) > 0 {
				p0 := GetVert3(verts, (indices[i] & T2(tmpT1)))
				p2 := GetVert3(verts, (indices[Next(i1, n)] & T2(tmpT1)))

				dx := p2[0] - p0[0]
				dy := p2[2] - p0[2]
				length := dx*dx + dy*dy

				if minLen < 0 || int(length) < minLen {
					minLen = int(length)
					mini = int(i)
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
			for i := int32(0); i < n; i++ {
				i1 := Next(i, n)
				i2 := Next(i1, n)
				if DiagonalLoose(i, i2, n, verts, indices) {
					p0 := GetVert4(verts, (indices[i] & T2(tmpT1)))
					p2 := GetVert4(verts, (indices[Next(i2, n)] & T2(tmpT1)))
					dx := p2[0] - p0[0]
					dy := p2[2] - p0[2]
					length := dx*dx + dy*dy

					if minLen < 0 || int(length) < int(minLen) {
						minLen = int(length)
						mini = int(i)
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
		i1 := Next(i, int(n))
		i2 := Next(i1, int(n))

		tris[dst] = indices[i] & T2(tmpT1)
		dst++
		tris[dst] = indices[i1] & T2(tmpT1)
		dst++
		tris[dst] = indices[i2] & T2(tmpT1)
		dst++
		ntris++

		// Removes P[i1] by copying P[i+1]...P[n-1] left one index.
		n--
		for k := i1; k < int(n); k++ {
			indices[k] = indices[k+1]
		}

		if i1 >= int(n) {
			i1 = 0
		}
		i = Prev(i1, int(n))

		// Update diagonal flags.
		if Diagonal[int32, T1, T2](Prev(int32(i), n), int32(i1), n, verts, indices) {
			indices[i] |= T2(tmpT2)
		} else {
			indices[i] &= T2(tmpT1)
		}

		if Diagonal[int32, T1, T2](int32(i), Next(int32(i1), n), n, verts, indices) {
			indices[i1] |= T2(tmpT2)
		} else {
			indices[i1] &= T2(tmpT1)
		}

	}

	// Append the remaining triangle.
	tris[dst] = indices[0] & T2(tmpT1)
	dst++
	tris[dst] = indices[1] & T2(tmpT1)
	dst++
	tris[dst] = indices[2] & T2(tmpT1)
	dst++
	ntris++

	return ntris
}

func IntersectSegContour(d0, d1 []int32, i, n int32, verts []int32) bool {
	// For each edge (k,k+1) of P
	for k := int32(0); k < n; k++ {
		k1 := Next(k, n)
		// Skip edges incident to i.
		if i == k || i == k1 {
			continue
		}

		p0 := GetVert4(verts, k)
		p1 := GetVert4(verts, k1)
		if RcVequal(d0, p0) || RcVequal(d1, p0) || RcVequal(d0, p1) || RcVequal(d1, p1) {
			continue
		}

		if Intersect(d0, d1, p0, p1) {
			return true
		}

	}
	return false
}

func CalcAreaOfPolygon2D(verts []int32, nverts int32) int32 {
	area := int32(0)
	i := int32(0)
	j := nverts - 1
	for i < nverts {
		vi := GetVert4(verts, i)
		vj := GetVert4(verts, j)
		area += vi[0]*vj[2] - vj[0]*vi[2]
		j = i
		i++
	}
	return (area + 1) / 2
}

func ContourdistancePtSeg[T int | int32](x, z, px, pz, qx, qz T) float32 {
	pqx := float32(qx - px)
	pqz := float32(qz - pz)
	dx := float32(x - px)
	dz := float32(z - pz)
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

	dx = float32(px) + t*pqx - float32(x)
	dz = float32(pz) + t*pqz - float32(z)

	return dx*dx + dz*dz
}

// TODO (graham): This is duplicated in the ConvexVolumeTool in RecastDemo
// / Checks if a point is contained within a polygon
// /
// / @param[in]	numVerts	Number of vertices in the polygon
// / @param[in]	verts		The polygon vertices
// / @param[in]	point		The point to check
// / @returns true if the point lies within the polygon, false otherwise.
func PointInPoly(numVerts int32, verts []float32, point []float32) bool {
	inPoly := false
	i := int32(0)
	j := numVerts - 1
	for i < numVerts {
		vi := verts[i*3 : i*3+3]
		vj := verts[j*3 : j*3+3]

		if (vi[2] > point[2]) == (vj[2] > point[2]) {
			j = i
			i++
			continue
		}

		if point[0] >= (vj[0]-vi[0])*(point[2]-vi[2])/(vj[2]-vi[2])+vi[0] {
			j = i
			i++
			continue
		}
		inPoly = !inPoly
		j = i
		i++
	}
	return inPoly
}

func Vcross2[T float64 | float32](p1, p2, p3 []T) T {
	u1 := p2[0] - p1[0]
	v1 := p2[2] - p1[2]
	u2 := p3[0] - p1[0]
	v2 := p3[2] - p1[2]
	return u1*v2 - v1*u2
}

func Vdot2[T float64 | float32](a, b []T) T {
	return a[0]*b[0] + a[2]*b[2]
}

func VdistSq2[T float64 | float32](p, q []T) T {
	dx := q[0] - p[0]
	dy := q[2] - p[2]
	return dx*dx + dy*dy
}

func Vdist2[T float64 | float32](p, q []T) T {
	return T(math.Sqrt(float64(VdistSq2(p, q))))
}

func CircumCircle[T float64 | float32](p1, p2, p3, c []T, r *T) bool {
	EPS := 1e-6
	// Calculate the circle relative to p1, to avoid some precision issues.
	v1 := []T{0, 0, 0}
	v2 := make([]T, 3)
	v3 := make([]T, 3)
	Vsub(v2, p2, p1)
	Vsub(v3, p3, p1)

	cp := Vcross2(v1, v2, v3)
	if math.Abs(float64(cp)) > EPS {
		v1Sq := T(Vdot2(v1, v1))
		v2Sq := Vdot2(v2, v2)
		v3Sq := Vdot2(v3, v3)
		c[0] = (v1Sq*(v2[2]-v3[2]) + v2Sq*(v3[2]-v1[2]) + v3Sq*(v1[2]-v2[2])) / (2 * cp)
		c[1] = 0
		c[2] = (v1Sq*(v3[0]-v2[0]) + v2Sq*(v1[0]-v3[0]) + v3Sq*(v2[0]-v1[0])) / (2 * cp)
		*r = Vdist2(c, v1)
		Vadd(c, c, p1)
		return true
	}

	copy(c, p1)
	*r = 0
	return false
}

func DistPtTri[T float64 | float32](p, a, b, c []T) (res T) {
	v0 := make([]T, 3)
	v1 := make([]T, 3)
	v2 := make([]T, 3)
	Vsub(v0, c, a)
	Vsub(v1, b, a)
	Vsub(v2, p, a)

	dot00 := Vdot2(v0, v0)
	dot01 := Vdot2(v0, v1)
	dot02 := Vdot2(v0, v2)
	dot11 := Vdot2(v1, v1)
	dot12 := Vdot2(v1, v2)

	// Compute barycentric coordinates
	invDenom := 1.0 / (dot00*dot11 - dot01*dot01)
	u := (dot11*dot02 - dot01*dot12) * invDenom
	v := (dot00*dot12 - dot01*dot02) * invDenom

	// If point lies inside the triangle, return interpolated y-coord.
	EPS := T(1e-4)
	if u >= -EPS && v >= -EPS && (u+v) <= 1+EPS {
		y := a[1] + v0[1]*u + v1[1]*v
		return T(math.Abs(float64(y - p[1])))
	}
	return GetTFloatMax(res)
}

func DistancePtSeg[T float64 | float32](pt, p, q []T) T {
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

func DistancePtSeg2d[T float64 | float32](pt, p, q []T) T {
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

func DistToTriMesh[T float64 | float32](p, verts []T, nverts int32, tris []int32, ntris int32) (res T) {
	maxValue := GetTFloatMax(res)
	dmin := maxValue
	for i := int32(0); i < ntris; i++ {
		va := GetVert3(verts, tris[i*4+0])
		vb := GetVert3(verts, tris[i*4+1])
		vc := GetVert3(verts, tris[i*4+2])
		d := DistPtTri(p, va, vb, vc)
		if d < dmin {
			dmin = d
		}
	}
	if dmin == maxValue {
		return -1
	}
	return dmin
}

func DistToPoly[T float64 | float32](nvert int32, verts []T, p []T) (res T) {

	dmin := GetTFloatMax(res)

	i := int32(0)
	j := nvert - 1
	c := false
	for i < nvert {
		vi := GetVert3(verts, i)
		vj := GetVert3(verts, j)
		if ((vi[2] > p[2]) != (vj[2] > p[2])) &&
			(p[0] < (vj[0]-vi[0])*(p[2]-vi[2])/(vj[2]-vi[2])+vi[0]) {
			c = !c
		}

		dmin = min(dmin, DistancePtSeg2d(p, vj, vi))
		j = i
		i++
	}
	if c {
		return -dmin
	}
	return dmin
}

func PushFront(v int32, arr []int32, an *int32) {
	*an++
	for i := *an - 1; i > 0; i-- {
		arr[i] = arr[i-1]
	}
	arr[0] = v
}

func PushBack(v int32, arr []int32, an *int32) {
	arr[*an] = v
	*an++
}

func Uleft[T IT](a, b, c []T) bool {
	return (b[0]-a[0])*(c[2]-a[2])-(c[0]-a[0])*(b[2]-a[2]) < 0
}

func ComputeTileHash(x, y, mask int32) int32 {
	h1 := uint32(0x8da6b343) // Large multiplicative constants;
	h2 := uint32(0xd8163841) // here arbitrarily chosen primes
	n := h1*uint32(x) + h2*uint32(y)
	return int32(n & uint32(mask))
}

// / Determines if two axis-aligned bounding boxes overlap.
// /  @param[in]		amin	Minimum bounds of box A. [(x, y, z)]
// /  @param[in]		amax	Maximum bounds of box A. [(x, y, z)]
// /  @param[in]		bmin	Minimum bounds of box B. [(x, y, z)]
// /  @param[in]		bmax	Maximum bounds of box B. [(x, y, z)]
// / @return True if the two AABB's overlap.
// / @see dtOverlapQuantBounds
func DtOverlapBounds(amin, amax, bmin, bmax []float32) bool {
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
