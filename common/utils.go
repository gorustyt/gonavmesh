package common

import "github.com/go-gl/mathgl/mgl32"

type Vec3 = mgl32.Vec3
type Vec4 = mgl32.Vec4
type Vec2 = mgl32.Vec2

type IT interface {
	~int | ~int8 | ~int16 | ~int32 | ~int64 |
		~uint | ~uint8 | ~uint16 | ~uint32 | ~uint64 |
		~float32 | ~float64
}
type IIndex interface {
	~int | ~int8 | ~int16 | ~int32 | ~uint | ~uint8 | ~uint16 | ~uint32
}

func GetVert3[T IT, T1 IIndex](verts []T, index T1) []T {
	return verts[index*3 : index*3+3]
}
func GetVert2[T IT, T1 IIndex](verts []T, index T1) []T {
	return verts[index*2 : index*2+2]
}
func GetVert4[T IT, T1 IIndex](verts []T, index T1) []T {
	return verts[index*4 : index*4+4]
}

func DoWhile(do func() (stop bool), while func() bool) {
	if do() {
		return
	}
	for while() {
		if do() {
			return
		}
	}
}

func SliceTToSlice[T1, T2 IT](v1 []T1) (v2 []T2) {
	for _, v := range v1 {
		v2 = append(v2, T2(v))
	}
	return v2
}
