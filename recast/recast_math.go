package recast

type IT interface {
	~int | ~int8 | ~int16 | ~int32 | ~int64 |
		~uint | ~uint8 | ~uint16 | ~uint32 | ~uint64 |
		~float32 | ~float64
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

// / Sets the neighbor connection data for the specified direction.
// / @param[in]		span			The span to update.
// / @param[in]		direction		The direction to set. [Limits: 0 <= value < 4]
// / @param[in]		neighborIndex	The index of the neighbor span.
func rcSetCon(span *RcCompactSpan, direction, neighborIndex int) {
	shift := direction * 6
	con := span.Con
	span.Con = (con & ^(0x3f << shift)) | ((neighborIndex & 0x3f) << shift)
}
