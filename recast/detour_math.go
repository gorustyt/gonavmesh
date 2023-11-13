package recast

import "math"

// / Returns the minimum of two values.
// /  @param[in]		a	Value A
// /  @param[in]		b	Value B
// /  @return The minimum of the two values.
func dtMin[T IT](a, b T) T {
	if a < b {
		return a
	}
	return b
}

// / Returns the maximum of two values.
// /  @param[in]		a	Value A
// /  @param[in]		b	Value B
// /  @return The maximum of the two values.
func dtMax[T IT](a, b T) T {
	if a > b {
		return a
	}
	return b
}

// / Returns the absolute value.
// /  @param[in]		a	The value.
// /  @return The absolute value of the specified value.
func dtAbs[T IT](a T) T {
	if a < 0 {
		return -a
	}
	return a
}

// / Returns the square of the value.
// /  @param[in]		a	The value.
func dtSqr[T IT](a T) T { return a * a }

// / Performs a vector addition. (@p v1 + @p v2)
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		v1		The base vector. [(x, y, z)]
// /  @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
func dtVadd(v1, v2 []float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[0] + v2[0]
	res[1] = v1[1] + v2[1]
	res[2] = v1[2] + v2[2]
	return res
}

// / Performs a vector subtraction. (@p v1 - @p v2)
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		v1		The base vector. [(x, y, z)]
// /  @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
func dtVsub(v1, v2 []float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[0] - v2[0]
	res[1] = v1[1] - v2[1]
	res[2] = v1[2] - v2[2]
	return res
}

// / Clamps the value to the specified range.
// /  @param[in]		v	The value to clamp.
// /  @param[in]		mn	The minimum permitted return value.
// /  @param[in]		mx	The maximum permitted return value.
// /  @return The value, clamped to the specified range.
func dtClamp[T IT](v, mn, mx T) T {
	if v < mn {
		return mn
	}
	if v > mx {
		return mx
	}
	return v
}

// / Performs a linear interpolation between two vectors. (@p v1 toward @p v2)
// /  @param[out]	dest	The result vector. [(x, y, x)]
// /  @param[in]		v1		The starting vector.
// /  @param[in]		v2		The destination vector.
// /	 @param[in]		t		The interpolation factor. [Limits: 0 <= value <= 1.0]
func dtVlerp(v1, v2 []float64, t float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[0] + (v2[0]-v1[0])*t
	res[1] = v1[1] + (v2[1]-v1[1])*t
	res[2] = v1[2] + (v2[2]-v1[2])*t
	return res
}

// / Performs a vector copy.
// /  @param[out]	dest	The result. [(x, y, z)]
// /  @param[in]		a		The vector to copy. [(x, y, z)]
func dtVcopy(a []float64) []float64 {
	res := make([]float64, 3)
	res[0] = a[0]
	res[1] = a[1]
	res[2] = a[2]
	return res
}

// / Derives the square of the scalar length of the vector. (len * len)
// /  @param[in]		v The vector. [(x, y, z)]
// / @return The square of the scalar length of the vector.
func dtVlenSqr(v []float64) float64 {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
}

// / Selects the minimum value of each element from the specified vectors.
// /  @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
// /  @param[in]	v	A vector. [(x, y, z)]
func dtVmin(mn, v []float64) []float64 {
	res := make([]float64, 3)
	res[0] = dtMin(mn[0], v[0])
	res[1] = dtMin(mn[1], v[1])
	res[2] = dtMin(mn[2], v[2])
	return res
}

// / Selects the maximum value of each element from the specified vectors.
// /  @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
// /  @param[in]		v	A vector. [(x, y, z)]
func dtVmax(mx, v []float64) []float64 {
	res := make([]float64, 3)
	res[0] = dtMax(mx[0], v[0])
	res[1] = dtMax(mx[1], v[1])
	res[2] = dtMax(mx[2], v[2])
	return res
}

/// @}
/// @name Computational geometry helper functions.
/// @{

// / Derives the signed xz-plane area of the triangle ABC, or the relationship of line AB to point C.
// /  @param[in]		a		Vertex A. [(x, y, z)]
// /  @param[in]		b		Vertex B. [(x, y, z)]
// /  @param[in]		c		Vertex C. [(x, y, z)]
// / @return The signed xz-plane area of the triangle.
func dtTriArea2D(a, b, c []float64) float64 {
	abx := b[0] - a[0]
	abz := b[2] - a[2]
	acx := c[0] - a[0]
	acz := c[2] - a[2]
	return acx*abz - abx*acz
}

// / Checks that the specified vector's components are all finite.
// /  @param[in]		v	A point. [(x, y, z)]
// / @return True if all of the point's components are finite, i.e. not NaN
// / or any of the infinities.
func dtVisfinite(v []float64) bool {
	return !(dtIsFinite(v[0]) && dtIsFinite(v[1]) && dtIsFinite(v[2]))
}
func dtIsFinite(v float64) bool {
	return math.IsInf(v, -1) || math.IsInf(v, 1) || math.IsNaN(v)
}

// / Checks that the specified vector's 2D components are finite.
// /  @param[in]		v	A point. [(x, y, z)]
func dtVisfinite2D(v []float64) bool {
	return dtIsFinite(v[0]) && dtIsFinite(v[2])
}

// / Returns the distance between two points.
// /  @param[in]		v1	A point. [(x, y, z)]
// /  @param[in]		v2	A point. [(x, y, z)]
// / @return The distance between the two points.
func dtVdist(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return math.Sqrt(dx*dx + dy*dy + dz*dz)
}

// / Normalizes the vector.
// /  @param[in,out]	v	The vector to normalize. [(x, y, z)]
func dtVnormalize(v []float64) {
	d := 1.0 / math.Sqrt(dtSqr(v[0])+dtSqr(v[1])+dtSqr(v[2]))
	v[0] *= d
	v[1] *= d
	v[2] *= d
	return
}

// / Performs a 'sloppy' colocation check of the specified points.
// /  @param[in]		p0	A point. [(x, y, z)]
// /  @param[in]		p1	A point. [(x, y, z)]
// / @return True if the points are considered to be at the same location.
// /
// / Basically, this function will return true if the specified points are
// / close enough to eachother to be considered colocated.
func dtVequal(p0 []float64, p1 []float64) bool {
	thr := dtSqr(1.0 / 16384.0)
	d := dtVdistSqr(p0, p1)
	return d < thr
}

// / Returns the square of the distance between two points.
// /  @param[in]		v1	A point. [(x, y, z)]
// /  @param[in]		v2	A point. [(x, y, z)]
// / @return The square of the distance between the two points.
func dtVdistSqr(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return dx*dx + dy*dy + dz*dz
}

// / Scales the vector by the specified value. (@p v * @p t)
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		v		The vector to scale. [(x, y, z)]
// /  @param[in]		t		The scaling factor.
func dtVscale(v []float64, t float64) []float64 {
	res := make([]float64, 3)
	res[0] = v[0] * t
	res[1] = v[1] * t
	res[2] = v[2] * t
	return res
}

// / Derives the xz-plane 2D perp product of the two vectors. (uz*vx - ux*vz)
// /  @param[in]		u		The LHV vector [(x, y, z)]
// /  @param[in]		v		The RHV vector [(x, y, z)]
// / @return The perp dot product on the xz-plane.
// /
// / The vectors are projected onto the xz-plane, so the y-values are ignored.
func dtVperp2D(u, v []float64) float64 {
	return u[2]*v[0] - u[0]*v[2]
}

// / Sets the vector elements to the specified values.
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		x		The x-value of the vector.
// /  @param[in]		y		The y-value of the vector.
// /  @param[in]		z		The z-value of the vector.
func dtVset(dest []float64, x, y, z float64) {
	dest[0] = x
	dest[1] = y
	dest[2] = z
}

// / Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		v1		The base vector. [(x, y, z)]
// /  @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
// /  @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
func dtVmad(dest, v1, v2 []float64, s float64) {
	dest[0] = v1[0] + v2[0]*s
	dest[1] = v1[1] + v2[1]*s
	dest[2] = v1[2] + v2[2]*s
}

// / Derives the dot product of two vectors on the xz-plane. (@p u . @p v)
// /  @param[in]		u		A vector [(x, y, z)]
// /  @param[in]		v		A vector [(x, y, z)]
// / @return The dot product on the xz-plane.
// /
// / The vectors are projected onto the xz-plane, so the y-values are ignored.
func dtVdot2D(u, v []float64) float64 {
	return u[0]*v[0] + u[2]*v[2]
}

func dtAlign4(x int) int { return (x + 3) & ^3 }
