package recast

import "math"

type IT interface {
	~int | ~int8 | ~int16 | ~int32 | ~int64 |
		~uint | ~uint8 | ~uint16 | ~uint32 | ~uint64 |
		~float32 | ~float64
}

// / Gets the standard width (x-axis) offset for the specified direction.
// / @param[in]		direction		The direction. [Limits: 0 <= value < 4]
// / @return The width offset to apply to the current cell position to move in the direction.
func rcGetDirOffsetX(direction int) int {
	offset := [4]int{-1, 0, 1, 0}
	return offset[direction&0x03]
}

// TODO (graham): Rename this to rcGetDirOffsetZ
// / Gets the standard height (z-axis) offset for the specified direction.
// / @param[in]		direction		The direction. [Limits: 0 <= value < 4]
// / @return The height offset to apply to the current cell position to move in the direction.
func rcGetDirOffsetY(direction int) int {
	offset := [4]int{0, 1, 0, -1}
	return offset[direction&0x03]
}

// / Returns the minimum of two values.
// / @param[in]		a	Value A
// / @param[in]		b	Value B
// / @return The minimum of the two values.
func rcMin[T IT](a, b T) T {
	if a < b {
		return a
	}
	return b
}

// / Returns the maximum of two values.
// / @param[in]		a	Value A
// / @param[in]		b	Value B
// / @return The maximum of the two values.
func rcMax[T IT](a, b T) T {
	if a > b {
		return a
	}
	return b
}

// / Returns the absolute value.
// / @param[in]		a	The value.
// / @return The absolute value of the specified value.
func rcAbs[T IT](a T) T {
	if a < 0 {
		return -a
	}
	return a
}

// / Returns the square of the value.
// / @param[in]		a	The value.
// / @return The square of the value.
func rcSqr[T IT](a T) T {
	return a * a
}

func rcSqrt(x float64) float64 {
	return math.Sqrt(x)
}

/// @}
/// @name Vector helper functions.
/// @{

// / Derives the cross product of two vectors. (@p v1 x @p v2)
// / @param[out]		dest	The cross product. [(x, y, z)]
// / @param[in]		v1		A Vector [(x, y, z)]
// / @param[in]		v2		A vector [(x, y, z)]
func rcVcross(v1, v2 []float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[1]*v2[2] - v1[2]*v2[1]
	res[1] = v1[2]*v2[0] - v1[0]*v2[2]
	res[2] = v1[0]*v2[1] - v1[1]*v2[0]
	return res
}

// / Derives the dot product of two vectors. (@p v1 . @p v2)
// / @param[in]		v1	A Vector [(x, y, z)]
// / @param[in]		v2	A vector [(x, y, z)]
// / @return The dot product.
func rcVdot(v1, v2 []float64) float64 {
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}

// / Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
// / @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
func rcVmad(v1, v2 []float64, s float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[0] + v2[0]*s
	res[1] = v1[1] + v2[1]*s
	res[2] = v1[2] + v2[2]*s
	return res
}

// / Performs a vector addition. (@p v1 + @p v2)
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
func rcVadd(v1, v2 []float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[0] + v2[0]
	res[1] = v1[1] + v2[1]
	res[2] = v1[2] + v2[2]
	return res
}

// / Performs a vector subtraction. (@p v1 - @p v2)
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
func rcVsub(v1, v2 []float64) []float64 {
	res := make([]float64, 3)
	res[0] = v1[0] - v2[0]
	res[1] = v1[1] - v2[1]
	res[2] = v1[2] - v2[2]
	return res
}

// / Selects the minimum value of each element from the specified vectors.
// / @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
// / @param[in]		v	A vector. [(x, y, z)]
func rcVmin(mn, v []float64) []float64 {
	res := make([]float64, 3)
	res[0] = rcMin(mn[0], v[0])
	res[1] = rcMin(mn[1], v[1])
	res[2] = rcMin(mn[2], v[2])
	return res
}

// / Selects the maximum value of each element from the specified vectors.
// / @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
// / @param[in]		v	A vector. [(x, y, z)]
func rcVmax(mx, v []float64) []float64 {
	res := make([]float64, 3)
	res[0] = rcMax(mx[0], v[0])
	res[1] = rcMax(mx[1], v[1])
	res[2] = rcMax(mx[2], v[2])
	return res
}

// / Performs a vector copy.
// / @param[out]		dest	The result. [(x, y, z)]
// / @param[in]		v		The vector to copy. [(x, y, z)]
func rcVcopy(v []float64) []float64 {
	res := make([]float64, 3)
	res[0] = v[0]
	res[1] = v[1]
	res[2] = v[2]
	return res
}

// / Returns the distance between two points.
// / @param[in]		v1	A point. [(x, y, z)]
// / @param[in]		v2	A point. [(x, y, z)]
// / @return The distance between the two points.
func rcVdist(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return rcSqrt(dx*dx + dy*dy + dz*dz)
}

// / Returns the square of the distance between two points.
// / @param[in]		v1	A point. [(x, y, z)]
// / @param[in]		v2	A point. [(x, y, z)]
// / @return The square of the distance between the two points.
func rcVdistSqr(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return dx*dx + dy*dy + dz*dz
}

// / Normalizes the vector.
// / @param[in,out]	v	The vector to normalize. [(x, y, z)]
func rcVnormalize(v []float64) []float64 {
	res := make([]float64, 3)
	d := 1.0 / rcSqrt(rcSqr(v[0])+rcSqr(v[1])+rcSqr(v[2]))
	res[0] *= d
	res[1] *= d
	res[2] *= d
	return res
}

func rcGetVert[T IT](verts []T, index int) []T {
	return verts[index*3 : index*3+3]
}

func rcGetVert2[T IT](verts []T, index int) []T {
	return verts[index*2 : index*2+2]
}
func rcGetVert4[T IT](verts []T, index int) []T {
	return verts[index*4 : index*4+4]
}

func rcGetTris(tris []int, index int) []int {
	return tris[index*4 : index*4+4]
}
