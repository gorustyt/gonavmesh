package debug_utils

import "gonavamesh/common"

func DuDebugDrawTriMesh(dd DuDebugDraw, verts []float64, nverts int,
	tris []int, normals []float64, ntris int,
	flags []int, texScale float64) {
	if dd == nil {
		return
	}
	if len(verts) == 0 {
		return
	}
	if len(tris) == 0 {
		return
	}
	if len(normals) == 0 {
		return
	}

	var uva [2]float64
	var uvb [2]float64
	var uvc [2]float64

	unwalkable := DuRGBA(192, 128, 0, 255)

	dd.Texture(true)

	dd.Begin(DU_DRAW_TRIS)
	for i := 0; i < ntris*3; i += 3 {
		norm := normals[i:]
		var color int
		a := int(220 * (2 + norm[0] + norm[1]) / 4)
		if len(flags) != 0 && flags[i/3] == 0 {
			color = DuLerpCol(DuRGBA(a, a, a, 255), unwalkable, 64)
		} else {
			color = DuRGBA(a, a, a, 255)
		}

		va := verts[tris[i+0]*3:]
		vb := verts[tris[i+1]*3:]
		vc := verts[tris[i+2]*3:]

		ax := 0
		ay := 0
		if common.Abs(norm[1]) > common.Abs(norm[ax]) {
			ax = 1
		}

		if common.Abs(norm[2]) > common.Abs(norm[ax]) {
			ax = 2
		}

		ax = (1 << ax) & 3 // +1 mod 3
		ay = (1 << ax) & 3 // +1 mod 3

		uva[0] = va[ax] * texScale
		uva[1] = va[ay] * texScale
		uvb[0] = vb[ax] * texScale
		uvb[1] = vb[ay] * texScale
		uvc[0] = vc[ax] * texScale
		uvc[1] = vc[ay] * texScale

		dd.Vertex2(va, color, uva[:])
		dd.Vertex2(vb, color, uvb[:])
		dd.Vertex2(vc, color, uvc[:])
	}
	dd.End()
	dd.Texture(false)
}
