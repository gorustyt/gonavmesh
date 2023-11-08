package recast

import (
	"cmp"
	"math"
)

// 叉乘，向量a->b * a->c
// area = 0→abc共线
// area<0→c在ab左侧
// area>0→c在ab右边
func area2(a, b, c *Point) int {
	return int((b.X-a.X)*(c.Y-a.Y) - (c.X-a.X)*(b.Y-a.Y))
}

func left(a, b, c *Point) bool {
	return area2(a, b, c) < 0
}

func leftOn(a, b, c *Point) bool {
	return area2(a, b, c) <= 0
}

// a,b,c是在一条直线上
func collinear(a, b, c *Point) bool {
	return area2(a, b, c) == 0
}

// 判断c是否在a,b线段上,垂直时需要判断y坐标，
func between(a, b, c *Point) bool {
	if !collinear(a, b, c) { //不共线
		return false
	}

	//不垂直直接判断x坐标
	if a.X != b.X {
		return ((a.X <= c.X) && (c.X <= b.X)) || ((a.X >= c.X) && (c.X >= b.X))
	} else {
		return ((a.Y <= c.Y) && (c.Y <= b.Y)) || ((a.Y >= c.Y) && (c.Y >= b.Y))
	}
}

// 判断两个点是否相等
func vequal(a, b *Point) bool {
	return a.X == b.X && a.Y == b.Y
}

// 判断线段相交时，如果用斜截式联立方程组求根的方式，一者斜截式有斜率90度、两直线平行等特殊情况；
// 二者方程组的解很可能不是整数，出现浮点计算误差影响结果。因此换用别的方式。
// 如果两条线段ab和cd内相交，那么cd就被ab所确定的直线切割为两部分。
// 所以对直线ab来说，点c和d一个在左，一个在右。类似地，对直线cd来说，点a和b也是一个在左，一个在右。如果这两个条件都满足，
// 就可以判定两条线段相交。即：
func intersectProp(a, b, c, d *Point) bool {
	// Eliminate improper cases.
	if collinear(a, b, c) || collinear(a, b, d) || collinear(c, d, a) || collinear(c, d, b) {
		return false
	}

	return xorb(left(a, b, c), left(a, b, d)) && xorb(left(c, d, a), left(c, d, b))
}

// 只要有一个参数为true才为true
func xorb(x, y bool) bool {
	return (x && !y) || (!x && y)
}

// 判断ab 和cd两个线段是否相交
func intersect(a, b, c, d *Point) bool {
	if intersectProp(a, b, c, d) { //相交
		return true
	}

	if between(a, b, c) || between(a, b, d) || //三点共线
		between(c, d, a) || between(c, d, b) {
		return true
	}

	return false
}

// 判断线段ab是不是一个多边形的对角线。
// 对角线要满足以下条件：
// 1- 端点是多边形的顶点；
// 2- 除了端点之外，不与多边形的任何一条边有交集；
// 3- 在多边形内部。
// 第一步，先不考虑条件（3），条件1、2实现起来不难。
// 基本的思路就是针对给定的线段ab，遍历多边形的所有边：
// 如果边的顶点分别是c、d，且c、d和a、b都不重合，并且ab与cd相交，
// 那么按照定义ab就一定不是多边形的顶点。
func diagonalie(i, j, n int, verts []*Point, indices []int) bool {
	d0 := verts[(indices[i]&0x0fffffff)*4]
	d1 := verts[(indices[j]&0x0fffffff)*4]
	// For each edge (k,k+1) of P
	for k := 0; k < n; k++ {
		k1 := next(k, n)
		// Skip edges incident to i or j
		if !((k == i) || (k1 == i) || (k == j) || (k1 == j)) {
			p0 := verts[(indices[k]&0x0fffffff)*4]
			p1 := verts[(indices[k1]&0x0fffffff)*4]
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

// 如何判断一个对角线ab是否在多边形内部？
// 假设要检查的对角线是线段ab，a是多边形的一个顶点。和a相邻的两个顶点是a0和a1。
// 按照反时针的顺序，三个顶点的次序是a0,a,a1。
// 分成两种情况讨论：
// 一，角a小于等于180度（这时候a称为一个“凸顶点”）
// 二，角a大于180度（这时候a称为一个“凹顶点”）
// 如果是凸顶点，ab在多边形内部可以等价为“a0在有向线段ab左边，a1在有向线段ba左边”。
//
// 如果是凹顶点，上面的判断条件就不对了，可能a0,a1在ab的同一侧而ab仍然在多边形内。
// 对凹顶点a，将多边形的“内部”和“外部”互换，a就成了一个凸顶点。
// 此时，“ab在多边形内”等价为“以下命题不成立：ab在新多边形内，或与角a的任何一边共线”。
// 即“以下命题不成立：a1在ab左侧或在ab上，且a0在ba左侧或在ba上”。
// 注意这里的条件都包含了“在直线上”的情形，因为按定义，
// 对角线和多边形边界的交集，只能是对角线的两个端点。
//
// 区分一个顶点a是凸顶点还是凹顶点，是通过判断a0,a,a1的关系。
// a是凸顶点，当且仅当a0在aa1左边或在aa1上。注意如果a0,a,a1共线，a的内角是180度。
// 按定义这种情况a被当作一个凸顶点。
func inCone(i, j, n int, verts []*Point, indices []int) bool {
	pi := verts[(indices[i]&0x0fffffff)*4]
	pj := verts[(indices[j]&0x0fffffff)*4]
	pi1 := verts[(indices[next(i, n)]&0x0fffffff)*4]
	pin1 := verts[(indices[prev(i, n)]&0x0fffffff)*4]
	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if leftOn(pin1, pi, pi1) {
		return left(pi, pj, pin1) && left(pj, pi, pi1)
	}
	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1))
}

// 是一个凸点，并且是对角线
func diagonal(i, j, n int, verts []*Point, indices []int) bool {
	return inCone(i, j, n, verts, indices) && diagonalie(i, j, n, verts, indices)
}

func diagonalieLoose(i, j, n int, verts []*Point, indices []int) bool {
	d0 := verts[(indices[i]&0x0fffffff)*4]
	d1 := verts[(indices[j]&0x0fffffff)*4]

	// For each edge (k,k+1) of P
	for k := 0; k < n; k++ {
		k1 := next(k, n)
		// Skip edges incident to i or j
		if !((k == i) || (k1 == i) || (k == j) || (k1 == j)) {
			p0 := verts[(indices[k]&0x0fffffff)*4]
			p1 := verts[(indices[k1]&0x0fffffff)*4]

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

func inConeLoose(i, j, n int, verts []*Point, indices []int) bool {
	pi := verts[(indices[i]&0x0fffffff)*4]
	pj := verts[(indices[j]&0x0fffffff)*4]
	pi1 := verts[(indices[next(i, n)]&0x0fffffff)*4]
	pin1 := verts[(indices[prev(i, n)]&0x0fffffff)*4]

	// If P[i] is a convex vertex [ i+1 left or on (i-1,i) ].
	if leftOn(pin1, pi, pi1) {
		return leftOn(pi, pj, pin1) && leftOn(pj, pi, pi1)
	}

	// Assume (i-1,i,i+1) not collinear.
	// else P[i] is reflex.
	return !(leftOn(pi, pj, pi1) && leftOn(pj, pi, pin1))
}

func diagonalLoose(i, j, n int, verts []*Point, indices []int) bool {
	return inConeLoose(i, j, n, verts, indices) && diagonalieLoose(i, j, n, verts, indices)
}

/// @}
/// @name Vector helper functions.
/// @{

// / Derives the cross product of two vectors. (@p v1 x @p v2)
// / @param[out]		dest	The cross product. [(x, y, z)]
// / @param[in]		v1		A Vector [(x, y, z)]
// / @param[in]		v2		A vector [(x, y, z)]
func rcVcross(v1, v2 *Point) *Point {
	var p Point
	p.X = v1.Y*v2.Z - v1.Z*v2.Y
	p.Y = v1.Z*v2.X - v1.X*v2.Z
	p.Z = v1.X*v2.Y - v1.Y*v2.X
	return &p
}

// / Derives the dot product of two vectors. (@p v1 . @p v2)
// / @param[in]		v1	A Vector [(x, y, z)]
// / @param[in]		v2	A vector [(x, y, z)]
// / @return The dot product.
func rcVdot(v1, v2 *Point) float64 {
	return v1.X*v2.X + v1.Y*v2.Y + v1.Z*v2.Z
}

// / Performs a scaled vector addition. (@p v1 + (@p v2 * @p s))
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to scale and add to @p v1. [(x, y, z)]
// / @param[in]		s		The amount to scale @p v2 by before adding to @p v1.
func rcVmad(v1, v2 *Point, s float64) *Point {
	var p Point
	p.X = v1.X + v2.X*s
	p.Y = v1.Y + v2.Y*s
	p.Z = v1.Z + v2.Z*s
	return &p
}

// / Performs a vector addition. (@p v1 + @p v2)
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to add to @p v1. [(x, y, z)]
func rcVadd(v1, v2 *Point) *Point {
	var p Point
	p.X = v1.X + v2.X
	p.Y = v1.Y + v2.Y
	p.Z = v1.Z + v2.Z
	return &p
}

// / Performs a vector subtraction. (@p v1 - @p v2)
// / @param[out]		dest	The result vector. [(x, y, z)]
// / @param[in]		v1		The base vector. [(x, y, z)]
// / @param[in]		v2		The vector to subtract from @p v1. [(x, y, z)]
func rcVsub(v1, v2 *Point) *Point {
	var p Point
	p.X = v1.X - v2.X
	p.Y = v1.Y - v2.Y
	p.Z = v1.Z - v2.Z
	return &p
}

// / Selects the minimum value of each element from the specified vectors.
// / @param[in,out]	mn	A vector.  (Will be updated with the result.) [(x, y, z)]
// / @param[in]		v	A vector. [(x, y, z)]
func rcVmin(mn, v *Point) *Point {
	var p Point
	mn.X = math.Min(mn.X, v.X)
	mn.Y = math.Min(mn.Y, v.Y)
	mn.Z = math.Min(mn.Z, v.Z)
	return &p
}

// / Selects the maximum value of each element from the specified vectors.
// / @param[in,out]	mx	A vector.  (Will be updated with the result.) [(x, y, z)]
// / @param[in]		v	A vector. [(x, y, z)]
func rcVmax(mx, v *Point) *Point {
	var p Point
	p.X = math.Max(mx.X, v.X)
	p.Y = math.Max(mx.Y, v.Y)
	p.Z = math.Max(mx.Z, v.Z)
	return &p
}

// / Clamps the value to the specified range.
// / @param[in]		value			The value to clamp.
// / @param[in]		minInclusive	The minimum permitted return value.
// / @param[in]		maxInclusive	The maximum permitted return value.
// / @return The value, clamped to the specified range.
func rcClamp[T cmp.Ordered](value, minInclusive, maxInclusive T) T {
	if value < minInclusive {
		return minInclusive
	} else {
		if value > maxInclusive {
			return maxInclusive
		} else {
			return value
		}
	}
}
