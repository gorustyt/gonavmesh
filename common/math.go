package common

import (
	"cmp"
	"github.com/go-gl/mathgl/mgl32"
	"math"
)

type IT interface {
	~int | ~int8 | ~int16 | ~int32 | ~int64 |
		~uint | ~uint8 | ~uint16 | ~uint32 | ~uint64 |
		~float32 | ~float64
}

// / Returns the square of the value.
// / @param[in]		a	The value.
// / @return The square of the value.
func Sqr[T IT](a T) T {
	return a * a
}

func Sqrt(x float64) float64 {
	return math.Sqrt(x)
}

// / Returns the absolute value.
// / @param[in]		a	The value.
// / @return The absolute value of the specified value.
func Abs[T IT](a T) T {
	if a < 0 {
		return -a
	}
	return a
}

// / Performs a vector addition. (@p v1 + @p v2)
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
func Vadd(res []float64, v1, v2 []float64) {
	res[0] = v1[0] + v2[0]
	res[1] = v1[1] + v2[1]
	res[2] = v1[2] + v2[2]
}

// / Performs a vector subtraction. (@p v1 - @p v2)
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
func Vsub(res, v1, v2 []float64) {
	res[0] = v1[0] - v2[0]
	res[1] = v1[1] - v2[1]
	res[2] = v1[2] - v2[2]
}

// / Selects the minimum value of each element from the specified vectors.
// / @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
// / @param[in]		v	A vector. [(x, y, z)]
func Vmin(mn, v []float64) {
	mn[0] = Min(mn[0], v[0])
	mn[1] = Min(mn[1], v[1])
	mn[2] = Min(mn[2], v[2])
}

// / Selects the maximum value of each element from the specified vectors.
// / @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
// / @param[in]		v	A vector. [(x, y, z)]
func Vmax(mx, v []float64) {
	mx[0] = Max(mx[0], v[0])
	mx[1] = Max(mx[1], v[1])
	mx[2] = Max(mx[2], v[2])
}

// / Returns the distance between two points.
// / @param[in]		v1	A point. [(x, y, z)]
// / @param[in]		v2	A point. [(x, y, z)]
// / @return The distance between the two points.
func Vdist(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return Sqrt(dx*dx + dy*dy + dz*dz)
}

// / Returns the square of the distance between two points.
// / @param[in]		v1	A point. [(x, y, z)]
// / @param[in]		v2	A point. [(x, y, z)]
// / @return The square of the distance between the two points.
func VdistSqr(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return dx*dx + dy*dy + dz*dz
}

// / Normalizes the vector.
// / @param[in,out]	v	The vector to normalize. [(x, y, z)]
func Vnormalize(v []float64) { //向量的单位化
	d := 1.0 / Sqrt(Sqr(v[0])+Sqr(v[1])+Sqr(v[2]))
	v[0] *= d
	v[1] *= d
	v[2] *= d
}

// / Returns the minimum of two values.
// / @param[in]		a	Value A
// / @param[in]		b	Value B
// / @return The minimum of the two values.
func Min[T IT](a, b T) T {
	if a < b {
		return a
	}
	return b
}

// / Returns the maximum of two values.
// / @param[in]		a	Value A
// / @param[in]		b	Value B
// / @return The maximum of the two values.
func Max[T IT](a, b T) T {
	if a > b {
		return a
	}
	return b
}

// / Performs a 'sloppy' colocation check of the specified points.
// /  @param[in]		p0	A point. [(x, y, z)]
// /  @param[in]		p1	A point. [(x, y, z)]
// / @return True if the points are considered to be at the same location.
// /
// / Basically, this function will return true if the specified points are
// / close enough to eachother to be considered colocated.
func Vequal(p0 []float64, p1 []float64) bool {
	thr := Sqr(1.0 / 16384.0)
	d := VdistSqr(p0, p1)
	return d < thr
}

// / Scales the vector by the specified value. (@p v * @p t)
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		v		The vector to scale. [(x, y, z)]
// /  @param[in]		t		The scaling factor.
func Vscale(res []float64, v []float64, t float64) {
	res[0] = v[0] * t
	res[1] = v[1] * t
	res[2] = v[2] * t
}

/// @}
/// @name Vector helper functions.
/// @{

// / Derives the cross product of two vectors. (@p v1 x @p v2)
// / @param[out]		dest	The cross product. [(x, y, z)]
// / @param[in]		v1		A Vector [(x, y, z)]
// / @param[in]		v2		A vector [(x, y, z)]
func Vcross(res []float64, v1, v2 []float64) { //求向量的叉集
	res[0] = v1[1]*v2[2] - v1[2]*v2[1]
	res[1] = v1[2]*v2[0] - v1[0]*v2[2]
	res[2] = v1[0]*v2[1] - v1[1]*v2[0]
}

// / Derives the dot product of two vectors. (@p v1 . @p v2)
// / @param[in]		v1	A Vector [(x, y, z)]
// / @param[in]		v2	A vector [(x, y, z)]
// / @return The dot product.
func Vdot(v1, v2 []float64) float64 { //求向量的点积
	return v1[0]*v2[0] + v1[1]*v2[1] + v1[2]*v2[2]
}

// / Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
// / @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
func Vmad(res []float64, v1, v2 []float64, s float64) {
	res[0] = v1[0] + v2[0]*s
	res[1] = v1[1] + v2[1]*s
	res[2] = v1[2] + v2[2]*s
}

// / Clamps the value to the specified range.
// / @param[in]		value			The value to clamp.
// / @param[in]		minInclusive	The minimum permitted return value.
// / @param[in]		maxInclusive	The maximum permitted return value.
// / @return The value, clamped to the specified range.
func Clamp[T cmp.Ordered](value, minInclusive, maxInclusive T) T {
	if value < minInclusive {
		return minInclusive
	}
	if value > maxInclusive {
		return maxInclusive
	}
	return value
}

// / Sets the vector elements to the specified values.
// /  @param[out]	dest	The result vector. [(x, y, z)]
// /  @param[in]		x		The x-value of the vector.
// /  @param[in]		y		The y-value of the vector.
// /  @param[in]		z		The z-value of the vector.
func Vset(dest []float64, x, y, z float64) {
	dest[0] = x
	dest[1] = y
	dest[2] = z
}

// / Checks that the specified vector's 2D components are finite.
// /  @param[in]		v	A point. [(x, y, z)]
func Visfinite2D(v []float64) bool {
	return IsFinite(v[0]) && IsFinite(v[2])
}

// / Derives the xz-plane 2D perp product of the two vectors. (uz*vx - ux*vz)
// /  @param[in]		u		The LHV vector [(x, y, z)]
// /  @param[in]		v		The RHV vector [(x, y, z)]
// / @return The perp dot product on the xz-plane.
// /
// / The vectors are projected onto the xz-plane, so the y-values are ignored.
func Vperp2D(u, v []float64) float64 {
	return u[2]*v[0] - u[0]*v[2]
}

// / Derives the dot product of two vectors on the xz-plane. (@p u . @p v)
// /  @param[in]		u		A vector [(x, y, z)]
// /  @param[in]		v		A vector [(x, y, z)]
// / @return The dot product on the xz-plane.
// /
// / The vectors are projected onto the xz-plane, so the y-values are ignored.
func Vdot2D(u, v []float64) float64 {
	return u[0]*v[0] + u[2]*v[2]
}

// / Derives the square of the distance between the specified points on the xz-plane.
// /  @param[in]		v1	A point. [(x, y, z)]
// /  @param[in]		v2	A point. [(x, y, z)]
// / @return The square of the distance between the point on the xz-plane.
func Vdist2DSqr(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dz := v2[2] - v1[2]
	return dx*dx + dz*dz
}

// / Derives the distance between the specified points on the xz-plane.
// /  @param[in]		v1	A point. [(x, y, z)]
// /  @param[in]		v2	A point. [(x, y, z)]
// / @return The distance between the point on the xz-plane.
// /
// / The vectors are projected onto the xz-plane, so the y-values are ignored.
func Vdist2D(v1, v2 []float64) float64 {
	dx := v2[0] - v1[0]
	dz := v2[2] - v1[2]
	return math.Sqrt(dx*dx + dz*dz)
}

// / Derives the scalar length of the vector.
// /  @param[in]		v The vector. [(x, y, z)]
// / @return The scalar length of the vector.
func Vlen(v []float64) float64 {
	return math.Sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])
}

// / Performs a linear interpolation between two vectors. (@p v1 toward @p v2)
// /  @param[out]	dest	The result vector. [(x, y, x)]
// /  @param[in]		v1		The starting vector.
// /  @param[in]		v2		The destination vector.
// /	 @param[in]		t		The interpolation factor. [Limits: 0 <= value <= 1.0]
func Vlerp(dest []float64, v1, v2 []float64, t float64) []float64 {
	dest[0] = v1[0] + (v2[0]-v1[0])*t
	dest[1] = v1[1] + (v2[1]-v1[1])*t
	dest[2] = v1[2] + (v2[2]-v1[2])*t
	return dest
}

// / Derives the square of the scalar length of the vector. (len * len)
// /  @param[in]		v The vector. [(x, y, z)]
// / @return The square of the scalar length of the vector.
func VlenSqr(v []float64) float64 {
	return v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
}

/// @}
/// @name Computational geometry helper functions.
/// @{

// / Derives the signed xz-plane area of the triangle ABC, or the relationship of line AB to point C.
// /  @param[in]		a		Vertex A. [(x, y, z)]
// /  @param[in]		b		Vertex B. [(x, y, z)]
// /  @param[in]		c		Vertex C. [(x, y, z)]
// / @return The signed xz-plane area of the triangle.
func TriArea2D(a, b, c []float64) float64 {
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
func Visfinite(v []float64) bool {
	return !(IsFinite(v[0]) && IsFinite(v[1]) && IsFinite(v[2]))
}
func IsFinite(v float64) bool {
	return math.IsInf(v, -1) || math.IsInf(v, 1) || math.IsNaN(v)
}

// / Gets the direction for the specified offset. One of x and y should be 0.
// / @param[in]		offsetX		The x offset. [Limits: -1 <= value <= 1]
// / @param[in]		offsetZ		The z offset. [Limits: -1 <= value <= 1]
// / @return The direction that represents the offset.
func GetDirForOffset(offsetX, offsetZ int) int {
	dirs := []int{3, 0, -1, 2, 1}
	return dirs[((offsetZ+1)<<1)+offsetX]
}

// / Gets the standard width (x-axis) offset for the specified direction.
// / @param[in]		direction		The direction. [Limits: 0 <= value < 4]
// / @return The width offset to apply to the current cell position to move in the direction.
func GetDirOffsetX(direction int) int {
	offset := [4]int{-1, 0, 1, 0}
	return offset[direction&0x03]
}

// TODO (graham): Rename this to rcGetDirOffsetZ
// / Gets the standard height (z-axis) offset for the specified direction.
// / @param[in]		direction		The direction. [Limits: 0 <= value < 4]
// / @return The height offset to apply to the current cell position to move in the direction.
func GetDirOffsetY(direction int) int {
	offset := [4]int{0, 1, 0, -1}
	return offset[direction&0x03]
}

func NextPow2(v int) int {
	v--
	v |= v >> 1
	v |= v >> 2
	v |= v >> 4
	v |= v >> 8
	v |= v >> 16
	v++
	return v
}

func Ilog2(v int) int {
	getBool := func(b bool) int {
		if b {
			return 1
		}
		return 0
	}
	var r int
	var shift int
	r = getBool((v > 0xffff)) << 4
	v >>= r
	shift = getBool((v > 0xff)) << 3
	v >>= shift
	r |= shift
	shift = getBool((v > 0xf)) << 2
	v >>= shift
	r |= shift
	shift = getBool((v > 0x3)) << 1
	v >>= shift
	r |= shift
	r |= (v >> 1)
	return r
}

func GetVs3[T IT](verts []T, index int) []T {
	return verts[index*3 : index*3+3]
}

func GluProject(obj []float64, projection, modelview []float64, view []int) []float64 {
	o := mgl32.Vec3{}
	for i := range obj {
		o[i] = float32(obj[i])
	}
	p := mgl32.Mat4{}
	for i := range projection {
		p[i] = float32(projection[i])
	}
	v := mgl32.Mat4{}
	for i := range modelview {
		v[i] = float32(modelview[i])
	}
	res := mgl32.Project(o, v, p, view[0], view[1], view[2], view[3])
	return []float64{float64(res[0]), float64(res[1]), float64(res[2])} //TODO
}
