package debug_utils

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/recast"
	"math"
)

func DuDebugDrawTriMesh(verts []float32, nverts int,
	tris []int, normals []float32, ntris int,
	flags []int, texScale float32) (res []float32) {

	if len(verts) == 0 {
		return
	}
	if len(tris) == 0 {
		return
	}
	if len(normals) == 0 {
		return
	}

	var uva [2]float32
	var uvb [2]float32
	var uvc [2]float32

	unwalkable := DuRGBA(192, 128, 0, 255)
	for i := 0; i < ntris*3; i += 3 {
		var color []float32
		norm := normals[i:]
		a := int(220 * (2 + norm[0] + norm[1]) / 4)
		if len(flags) != 0 && flags[i/3] == 0 {
			color = NormalizeDuLerpCol(DuRGBA(a, a, a, 255), unwalkable, 64)
		} else {
			color = NormalizeRgba(a, a, a, 255)
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
		res = append(res, va...)
		res = append(res, color...)
		res = append(res, uva[:]...)
		res = append(res, vb...)
		res = append(res, color...)
		res = append(res, uvb[:]...)
		res = append(res, vc...)
		res = append(res, color...)
		res = append(res, uvc[:]...)
	}
	return res
}

func DuDebugDrawTriMeshSlope(verts []float32, nverts int,
	tris []int, normals []float32, ntris int,
	walkableSlopeAngle float32, texScale float32) (res []float32) {
	if len(verts) == 0 {
		return
	}
	if len(tris) == 0 {
		return
	}
	if len(normals) == 0 {
		return
	}

	walkableThr := float32(math.Cos(float64(walkableSlopeAngle / 180.0 * math.Pi)))

	uva := make([]float32, 2)
	uvb := make([]float32, 2)
	uvc := make([]float32, 2)
	unwalkable := DuRGBA(192, 128, 0, 255)
	for i := 0; i < ntris*3; i += 3 {
		norm := normals[i:]
		var color []float32
		a := int(220 * (2 + norm[0] + norm[1]) / 4)
		if norm[1] < walkableThr {
			color = NormalizeDuLerpCol(DuRGBA(a, a, a, 255), unwalkable, 64)
		} else {
			color = NormalizeRgba(a, a, a, 255)
		}

		va := verts[tris[i+0]*3:]
		vb := verts[tris[i+1]*3:]
		vc := verts[tris[i+2]*3:]

		ax := 0
		ay := 0
		if math.Abs(float64(norm[1])) > math.Abs(float64(norm[ax])) {
			ax = 1
		}

		if math.Abs(float64(norm[2])) > math.Abs(float64(norm[ax])) {
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

		res = append(res, va...)
		res = append(res, color...)
		res = append(res, uva[:]...)
		res = append(res, vb...)
		res = append(res, color...)
		res = append(res, uvb[:]...)
		res = append(res, vc...)
		res = append(res, color...)
		res = append(res, uvc[:]...)
	}
	return res
}

func DuDebugDrawHeightfieldSolid(dd DuDebugDraw, hf *recast.RcHeightfield) {
	if dd == nil {
		return
	}

	orig := hf.Bmin
	cs := float32(hf.Cs)
	ch := float32(hf.Ch)

	w := int(hf.Width)
	h := int(hf.Height)

	fcol := make([]int, 6)
	DuCalcBoxColors(fcol, DuRGBA(255, 255, 255, 255), DuRGBA(255, 255, 255, 255))

	dd.Begin(DU_DRAW_QUADS)

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			fx := float32(orig[0]) + float32(x)*cs
			fz := float32(orig[2]) + float32(y)*cs
			s := hf.Spans[x+y*w]
			for s != nil {
				DuAppendBox(dd, fx, float32(orig[1])+float32(s.Smin)*ch, fz, fx+cs, float32(
					orig[1])+float32(s.Smax)*ch, fz+cs, fcol)
				s = s.Next
			}
		}
	}
	dd.End()
}

func DuDebugDrawHeightfieldWalkable(dd DuDebugDraw, hf *recast.RcHeightfield) {
	if dd == nil {
		return
	}

	orig := hf.Bmin
	cs := float32(hf.Cs)
	ch := float32(hf.Ch)

	w := int(hf.Width)
	h := int(hf.Height)

	fcol := make([]int, 6)
	DuCalcBoxColors(fcol, DuRGBA(255, 255, 255, 255), DuRGBA(217, 217, 217, 255))

	dd.Begin(DU_DRAW_QUADS)

	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			fx := float32(orig[0]) + float32(x)*cs
			fz := float32(orig[2]) + float32(y)*cs
			s := hf.Spans[x+y*w]
			for s != nil {
				if s.Area == recast.RC_WALKABLE_AREA {
					fcol[0] = DuRGBA(64, 128, 160, 255)
				} else if s.Area == recast.RC_NULL_AREA {
					fcol[0] = DuRGBA(64, 64, 64, 255)
				} else {
					fcol[0] = duMultCol(dd.AreaToCol(int(s.Area)), 200)
				}

				DuAppendBox(dd, fx, float32(orig[1])+float32(s.Smin)*ch, fz, fx+cs, float32(orig[1])+float32(s.Smax)*ch, fz+cs, fcol)
				s = s.Next
			}
		}
	}

	dd.End()
}

func DuDebugDrawCompactHeightfieldSolid(dd DuDebugDraw, chf *recast.RcCompactHeightfield) {
	if dd == nil {
		return
	}

	cs := float32(chf.Cs)
	ch := float32(chf.Ch)

	dd.Begin(DU_DRAW_QUADS)

	for y := int32(0); y < chf.Height; y++ {
		for x := int32(0); x < chf.Width; x++ {
			fx := float32(
				chf.Bmin[0]) + float32(x)*cs
			fz := float32(chf.Bmin[2]) + float32(y)*cs
			c := chf.Cells[x+y*chf.Width]
			i := c.Index
			ni := c.Index + c.Count
			for ; i < ni; i++ {
				s := chf.Spans[i]

				area := chf.Areas[i]
				var color int
				if area == recast.RC_WALKABLE_AREA {
					color = DuRGBA(0, 192, 255, 64)
				} else if area == recast.RC_NULL_AREA {
					color = DuRGBA(0, 0, 0, 64)
				} else {
					color = dd.AreaToCol(int(area))
				}

				fy := float32(chf.Bmin[1]) + float32(s.Y+1)*ch
				dd.Vertex1(fx, fy, fz, color)
				dd.Vertex1(fx, fy, fz+cs, color)
				dd.Vertex1(fx+cs, fy, fz+cs, color)
				dd.Vertex1(fx+cs, fy, fz, color)
			}
		}
	}
	dd.End()
}

func DuDebugDrawCompactHeightfieldRegions(dd DuDebugDraw, chf *recast.RcCompactHeightfield) {
	if dd == nil {
		return
	}

	cs := float32(chf.Cs)
	ch := float32(chf.Ch)

	dd.Begin(DU_DRAW_QUADS)

	for y := int32(0); y < chf.Height; y++ {
		for x := int32(0); x < chf.Width; x++ {
			fx := float32(chf.Bmin[0]) + float32(x)*cs
			fz := float32(chf.Bmin[2]) + float32(y)*cs
			c := chf.Cells[x+y*chf.Width]

			i := c.Index
			ni := c.Index + c.Count
			for ; i < ni; i++ {
				s := chf.Spans[i]
				fy := float32(chf.Bmin[1]) + float32(s.Y)*ch
				var color int
				if s.Reg != 0 {
					color = DuIntToCol(int(s.Reg), 192)
				} else {

					color = DuRGBA(0, 0, 0, 64)
				}

				dd.Vertex1(fx, fy, fz, color)
				dd.Vertex1(fx, fy, fz+cs, color)
				dd.Vertex1(fx+cs, fy, fz+cs, color)
				dd.Vertex1(fx+cs, fy, fz, color)
			}
		}
	}

	dd.End()
}

func DuDebugDrawCompactHeightfieldDistance(dd DuDebugDraw, chf *recast.RcCompactHeightfield) {
	if dd == nil {
		return
	}
	if len(chf.Dist) == 0 {
		return
	}

	cs := float32(chf.Cs)
	ch := float32(chf.Ch)

	maxd := chf.MaxDistance
	if maxd < 1.0 {
		maxd = 1
	}
	dscale := 255.0 / maxd

	dd.Begin(DU_DRAW_QUADS)

	for y := int32(0); y < chf.Height; y++ {
		for x := int32(0); x < chf.Width; x++ {
			fx := float32(chf.Bmin[0]) + float32(x)*cs
			fz := float32(chf.Bmin[2]) + float32(y)*cs
			c := chf.Cells[x+y*chf.Width]

			i := c.Index
			ni := c.Index + c.Count
			for ; i < ni; i++ {
				s := chf.Spans[i]
				fy := float32(chf.Bmin[1]) + float32(s.Y+1)*ch
				cd := int(chf.Dist[i] * dscale)
				color := DuRGBA(cd, cd, cd, 255)
				dd.Vertex1(fx, fy, fz, color)
				dd.Vertex1(fx, fy, fz+cs, color)
				dd.Vertex1(fx+cs, fy, fz+cs, color)
				dd.Vertex1(fx+cs, fy, fz, color)
			}
		}
	}
	dd.End()
}

func drawLayerPortals(dd DuDebugDraw, layer *recast.RcHeightfieldLayer) {
	cs := float32(layer.Cs)
	ch := float32(layer.Ch)
	w := int(layer.Width)
	h := int(layer.Height)

	pcol := DuRGBA(255, 255, 255, 255)

	segs := [4 * 4]int{0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0}

	// Layer portals
	dd.Begin(DU_DRAW_LINES, 2.0)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			idx := x + y*w
			lh := layer.Heights[idx]
			if lh == 255 {
				continue
			}

			for dir := 0; dir < 4; dir++ {
				if layer.Cons[idx]&(1<<(dir+4)) != 0 {
					seg := segs[dir*4:]
					ax := float32(layer.Bmin[0]) + float32(x+seg[0])*cs
					ay := float32(layer.Bmin[1]) + float32(lh+2)*ch
					az := float32(layer.Bmin[2]) + float32(y+seg[1])*cs
					bx := float32(layer.Bmin[0]) + float32(x+seg[2])*cs
					by := float32(layer.Bmin[1]) + float32(lh+2)*ch
					bz := float32(layer.Bmin[2]) + float32(y+seg[3])*cs
					dd.Vertex1(ax, ay, az, pcol)
					dd.Vertex1(bx, by, bz, pcol)
				}
			}
		}
	}
	dd.End()
}

func DuDebugDrawHeightfieldLayer(dd DuDebugDraw, layer *recast.RcHeightfieldLayer, idx int) {
	cs := float32(layer.Cs)
	ch := float32(layer.Ch)
	w := int(layer.Width)
	h := int(layer.Height)

	color := DuIntToCol(idx+1, 255)

	// Layer bounds
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	bmin[0] = float32(layer.Bmin[0]) + float32(layer.Minx)*cs
	bmin[1] = float32(layer.Bmin[1])
	bmin[2] = float32(layer.Bmin[2]) + float32(layer.Miny)*cs
	bmax[0] = float32(layer.Bmin[0]) + float32(layer.Maxx+1)*cs
	bmax[1] = float32(layer.Bmax[1])
	bmax[2] = float32(layer.Bmin[2]) + float32(layer.Maxy+1)*cs
	DuDebugDrawBoxWire(dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], DuTransCol(color, 128), 2.0)

	// Layer height
	dd.Begin(DU_DRAW_QUADS)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			lidx := x + y*w
			lh := layer.Heights[lidx]
			if h == 0xff {
				continue
			}
			area := layer.Areas[lidx]

			var col int
			if area == recast.RC_WALKABLE_AREA {
				col = DuLerpCol(color, DuRGBA(0, 192, 255, 64), 32)
			} else if area == recast.RC_NULL_AREA {
				col = DuLerpCol(color, DuRGBA(0, 0, 0, 64), 32)
			} else {
				col = DuLerpCol(color, dd.AreaToCol(int(area)), 32)
			}

			fx := float32(layer.Bmin[0]) + float32(x)*cs
			fy := float32(layer.Bmin[1]) + float32(lh+1)*ch
			fz := float32(layer.Bmin[2]) + float32(y)*cs

			dd.Vertex1(fx, fy, fz, col)
			dd.Vertex1(fx, fy, fz+cs, col)
			dd.Vertex1(fx+cs, fy, fz+cs, col)
			dd.Vertex1(fx+cs, fy, fz, col)
		}
	}
	dd.End()

	// Portals
	drawLayerPortals(dd, layer)
}

func DuDebugDrawHeightfieldLayers(dd DuDebugDraw, lset *recast.RcHeightfieldLayerSet) {
	if dd == nil {
		return
	}
	for i := int32(0); i < lset.Nlayers; i++ {
		DuDebugDrawHeightfieldLayer(dd, lset.Layers[i], int(i))
	}

}

/*
void duDebugDrawLayerContours(duDebugDraw* dd, const struct rcLayerContourSet& lcset)
{
	if (!dd) return;

	const float* orig = lcset.bmin;
	const float cs = lcset.cs;
	const float ch = lcset.ch;

	const unsigned char a = 255;// (unsigned char)(alpha*255.0f);

	const int offs[2*4] = {-1,0, 0,1, 1,0, 0,-1};

	dd->begin(DU_DRAW_LINES, 2.0f);

	for (int i = 0; i < lcset.nconts; ++i)
	{
		const rcLayerContour& c = lcset.conts[i];
		unsigned int color = 0;

		color = duIntToCol(i, a);

		for (int j = 0; j < c.nverts; ++j)
		{
			const int k = (j+1) % c.nverts;
			const unsigned char* va = &c.verts[j*4];
			const unsigned char* vb = &c.verts[k*4];
			const float ax = orig[0] + va[0]*cs;
			const float ay = orig[1] + (va[1]+1+(i&1))*ch;
			const float az = orig[2] + va[2]*cs;
			const float bx = orig[0] + vb[0]*cs;
			const float by = orig[1] + (vb[1]+1+(i&1))*ch;
			const float bz = orig[2] + vb[2]*cs;
			unsigned int col = color;
			if ((va[3] & 0xf) != 0xf)
			{
				col = duRGBA(255,255,255,128);
				int d = va[3] & 0xf;

				const float cx = (ax+bx)*0.5f;
				const float cy = (ay+by)*0.5f;
				const float cz = (az+bz)*0.5f;

				const float dx = cx + offs[d*2+0]*2*cs;
				const float dy = cy;
				const float dz = cz + offs[d*2+1]*2*cs;

				dd->vertex(cx,cy,cz,duRGBA(255,0,0,255));
				dd->vertex(dx,dy,dz,duRGBA(255,0,0,255));
			}

			duAppendArrow(dd, ax,ay,az, bx,by,bz, 0.0f, cs*0.5f, col);
		}
	}
	dd->end();

	dd->begin(DU_DRAW_POINTS, 4.0f);

	for (int i = 0; i < lcset.nconts; ++i)
	{
		const rcLayerContour& c = lcset.conts[i];
		unsigned int color = 0;

		for (int j = 0; j < c.nverts; ++j)
		{
			const unsigned char* va = &c.verts[j*4];

			color = duDarkenCol(duIntToCol(i, a));
			if (va[3] & 0x80)
				color = duRGBA(255,0,0,255);

			float fx = orig[0] + va[0]*cs;
			float fy = orig[1] + (va[1]+1+(i&1))*ch;
			float fz = orig[2] + va[2]*cs;
			dd->vertex(fx,fy,fz, color);
		}
	}
	dd->end();
}

void duDebugDrawLayerPolyMesh(duDebugDraw* dd, const struct rcLayerPolyMesh& lmesh)
{
	if (!dd) return;

	const int nvp = lmesh.nvp;
	const float cs = lmesh.cs;
	const float ch = lmesh.ch;
	const float* orig = lmesh.bmin;

	const int offs[2*4] = {-1,0, 0,1, 1,0, 0,-1};

	dd->begin(DU_DRAW_TRIS);

	for (int i = 0; i < lmesh.npolys; ++i)
	{
		const unsigned short* p = &lmesh.polys[i*nvp*2];

		unsigned int color;
		if (lmesh.areas[i] == RC_WALKABLE_AREA)
			color = duRGBA(0,192,255,64);
		else if (lmesh.areas[i] == RC_NULL_AREA)
			color = duRGBA(0,0,0,64);
		else
			color = duIntToCol(lmesh.areas[i], 255);

		unsigned short vi[3];
		for (int j = 2; j < nvp; ++j)
		{
			if (p[j] == RC_MESH_NULL_IDX) break;
			vi[0] = p[0];
			vi[1] = p[j-1];
			vi[2] = p[j];
			for (int k = 0; k < 3; ++k)
			{
				const unsigned short* v = &lmesh.verts[vi[k]*3];
				const float x = orig[0] + v[0]*cs;
				const float y = orig[1] + (v[1]+1)*ch;
				const float z = orig[2] + v[2]*cs;
				dd->vertex(x,y,z, color);
			}
		}
	}
	dd->end();

	// Draw neighbours edges
	const unsigned int coln = duRGBA(0,48,64,32);
	dd->begin(DU_DRAW_LINES, 1.5f);
	for (int i = 0; i < lmesh.npolys; ++i)
	{
		const unsigned short* p = &lmesh.polys[i*nvp*2];
		for (int j = 0; j < nvp; ++j)
		{
			if (p[j] == RC_MESH_NULL_IDX) break;
			if (p[nvp+j] & 0x8000) continue;
			const int nj = (j+1 >= nvp || p[j+1] == RC_MESH_NULL_IDX) ? 0 : j+1;
			int vi[2] = {p[j], p[nj]};

			for (int k = 0; k < 2; ++k)
			{
				const unsigned short* v = &lmesh.verts[vi[k]*3];
				const float x = orig[0] + v[0]*cs;
				const float y = orig[1] + (v[1]+1)*ch + 0.1f;
				const float z = orig[2] + v[2]*cs;
				dd->vertex(x, y, z, coln);
			}
		}
	}
	dd->end();

	// Draw boundary edges
	const unsigned int colb = duRGBA(0,48,64,220);
	dd->begin(DU_DRAW_LINES, 2.5f);
	for (int i = 0; i < lmesh.npolys; ++i)
	{
		const unsigned short* p = &lmesh.polys[i*nvp*2];
		for (int j = 0; j < nvp; ++j)
		{
			if (p[j] == RC_MESH_NULL_IDX) break;
			if ((p[nvp+j] & 0x8000) == 0) continue;
			const int nj = (j+1 >= nvp || p[j+1] == RC_MESH_NULL_IDX) ? 0 : j+1;
			int vi[2] = {p[j], p[nj]};

			unsigned int col = colb;
			if ((p[nvp+j] & 0xf) != 0xf)
			{
				const unsigned short* va = &lmesh.verts[vi[0]*3];
				const unsigned short* vb = &lmesh.verts[vi[1]*3];

				const float ax = orig[0] + va[0]*cs;
				const float ay = orig[1] + (va[1]+1+(i&1))*ch;
				const float az = orig[2] + va[2]*cs;
				const float bx = orig[0] + vb[0]*cs;
				const float by = orig[1] + (vb[1]+1+(i&1))*ch;
				const float bz = orig[2] + vb[2]*cs;

				const float cx = (ax+bx)*0.5f;
				const float cy = (ay+by)*0.5f;
				const float cz = (az+bz)*0.5f;

				int d = p[nvp+j] & 0xf;

				const float dx = cx + offs[d*2+0]*2*cs;
				const float dy = cy;
				const float dz = cz + offs[d*2+1]*2*cs;

				dd->vertex(cx,cy,cz,duRGBA(255,0,0,255));
				dd->vertex(dx,dy,dz,duRGBA(255,0,0,255));

				col = duRGBA(255,255,255,128);
			}

			for (int k = 0; k < 2; ++k)
			{
				const unsigned short* v = &lmesh.verts[vi[k]*3];
				const float x = orig[0] + v[0]*cs;
				const float y = orig[1] + (v[1]+1)*ch + 0.1f;
				const float z = orig[2] + v[2]*cs;
				dd->vertex(x, y, z, col);
			}
		}
	}
	dd->end();

	dd->begin(DU_DRAW_POINTS, 3.0f);
	const unsigned int colv = duRGBA(0,0,0,220);
	for (int i = 0; i < lmesh.nverts; ++i)
	{
		const unsigned short* v = &lmesh.verts[i*3];
		const float x = orig[0] + v[0]*cs;
		const float y = orig[1] + (v[1]+1)*ch + 0.1f;
		const float z = orig[2] + v[2]*cs;
		dd->vertex(x,y,z, colv);
	}
	dd->end();
}
*/

func GetContourCenter(cont *recast.RcContour, orig []float32, cs, ch float32, center []float32) {
	center[0] = 0
	center[1] = 0
	center[2] = 0
	if cont.Nverts == 0 {
		return
	}

	for i := int32(0); i < cont.Nverts; i++ {
		v := cont.Verts[i*4:]
		center[0] += float32(v[0])
		center[1] += float32(v[1])
		center[2] += float32(v[2])
	}
	s := 1.0 / float32(cont.Nverts)
	center[0] *= s * cs
	center[1] *= s * ch
	center[2] *= s * cs
	center[0] += orig[0]
	center[1] += orig[1] + 4*ch
	center[2] += orig[2]
}

func findContourFromSet(cset *recast.RcContourSet, reg int) *recast.RcContour {
	for i := int32(0); i < cset.Nconts; i++ {
		if int(cset.Conts[i].Reg) == reg {
			return cset.Conts[i]
		}

	}
	return nil
}

func DuDebugDrawRegionConnections(dd DuDebugDraw, cset *recast.RcContourSet, alphas ...float32) {
	alpha := float32(1.0)
	if len(alphas) > 0 {
		alpha = alphas[0]
	}
	if dd == nil {
		return
	}

	orig := cset.Bmin
	cs := float32(cset.Cs)
	ch := float32(cset.Ch)

	// Draw centers
	pos := make([]float32, 3)
	pos2 := make([]float32, 3)

	color := DuRGBA(0, 0, 0, 196)

	dd.Begin(DU_DRAW_LINES, 2.0)

	for i := int32(0); i < cset.Nconts; i++ {
		cont := cset.Conts[i]
		GetContourCenter(cont, common.SliceTToSlice[float32, float32](orig), cs, ch, pos)
		for j := int32(0); j < cont.Nverts; j++ {
			v := cont.Verts[j*4:]
			if v[3] == 0 || v[3] < int32(cont.Reg) {
				continue
			}
			cont2 := findContourFromSet(cset, int(v[3]))
			if cont2 != nil {
				GetContourCenter(cont2, common.SliceTToSlice[float32, float32](orig), cs, ch, pos2)
				DuAppendArc(dd, pos[0], pos[1], pos[2], pos2[0], pos2[1], pos2[2], 0.25, 0.6, 0.6, color)
			}
		}
	}

	dd.End()

	a := int(alpha * 255.0)

	dd.Begin(DU_DRAW_POINTS, 7.0)

	for i := int32(0); i < cset.Nconts; i++ {
		cont := cset.Conts[i]
		col := DuDarkenCol(DuIntToCol(int(cont.Reg), a))
		GetContourCenter(cont, common.SliceTToSlice[float32, float32](orig), cs, ch, pos)
		dd.Vertex(pos, col)
	}
	dd.End()
}

func DuDebugDrawRawContours(dd DuDebugDraw, cset *recast.RcContourSet, alphas ...float32) {
	alpha := float32(1.0)
	if len(alphas) > 0 {
		alpha = alphas[0]
	}
	if dd == nil {
		return
	}

	orig := cset.Bmin
	cs := float32(cset.Cs)
	ch := float32(cset.Ch)

	a := int(alpha * 255.0)

	dd.Begin(DU_DRAW_LINES, 2.0)

	for i := int32(0); i < cset.Nconts; i++ {
		c := cset.Conts[i]
		color := DuIntToCol(int(c.Reg), a)

		for j := int32(0); j < c.Nrverts; j++ {
			v := c.Rverts[j*4:]
			fx := float32(orig[0]) + float32(v[0])*cs
			fy := float32(orig[1]) + float32(v[1]+1+(i&1))*ch
			fz := float32(orig[2]) + float32(v[2])*cs
			dd.Vertex1(fx, fy, fz, color)
			if j > 0 {
				dd.Vertex1(fx, fy, fz, color)
			}

		}
		// Loop last segment.
		v := c.Rverts
		fx := float32(orig[0]) + float32(v[0])*cs
		fy := float32(orig[1]) + float32(v[1]+1+(i&1))*ch
		fz := float32(orig[2]) + float32(v[2])*cs
		dd.Vertex1(fx, fy, fz, color)
	}
	dd.End()

	dd.Begin(DU_DRAW_POINTS, 2.0)

	for i := int32(0); i < cset.Nconts; i++ {
		c := cset.Conts[i]
		color := DuDarkenCol(DuIntToCol(int(c.Reg), a))

		for j := int32(0); j < c.Nrverts; j++ {
			v := c.Rverts[j*4:]
			off := float32(0.0)
			colv := color
			if v[3]&recast.RC_BORDER_VERTEX != 0 {
				colv = DuRGBA(255, 255, 255, a)
				off = ch * 2
			}

			fx := orig[0] + float32(v[0])*cs
			fy := orig[1] + float32(v[1]+1+(i&1))*ch + off
			fz := orig[2] + float32(v[2])*cs
			dd.Vertex1(fx, fy, fz, colv)
		}
	}
	dd.End()
}

func DuDebugDrawContours(dd DuDebugDraw, cset *recast.RcContourSet, alphas ...float32) {
	alpha := float32(1.0)
	if len(alphas) > 0 {
		alpha = alphas[0]
	}
	if dd == nil {
		return
	}

	orig := cset.Bmin
	cs := float32(cset.Cs)
	ch := float32(cset.Ch)

	a := int(alpha * 255.0)

	dd.Begin(DU_DRAW_LINES, 2.5)

	for i := int32(0); i < cset.Nconts; i++ {
		c := cset.Conts[i]
		if c.Nverts == 0 {
			continue
		}

		color := DuIntToCol(int(c.Reg), a)
		bcolor := DuLerpCol(color, DuRGBA(255, 255, 255, a), 128)
		j := int32(0)
		k := c.Nverts - 1
		for j < c.Nverts {
			va := c.Verts[k*4:]
			vb := c.Verts[j*4:]
			col := color
			if va[3]&recast.RC_AREA_BORDER != 0 {
				col = bcolor
			}
			fx := float32(orig[0]) + float32(va[0])*cs
			fy := float32(orig[1]) + float32(va[1]+1+(i&1))*ch
			fz := float32(orig[2]) + float32(va[2])*cs
			dd.Vertex1(fx, fy, fz, col)
			fx = float32(orig[0]) + float32(vb[0])*cs
			fy = float32(orig[1]) + float32(vb[1]+1+(i&1))*ch
			fz = float32(orig[2]) + float32(vb[2])*cs
			dd.Vertex1(fx, fy, fz, col)
			k = j
			j++
		}
	}
	dd.End()

	dd.Begin(DU_DRAW_POINTS, 3.0)

	for i := int32(0); i < cset.Nconts; i++ {
		c := cset.Conts[i]
		color := DuDarkenCol(DuIntToCol(int(c.Reg), a))
		for j := int32(0); j < c.Nverts; j++ {
			v := c.Verts[j*4:]
			off := float32(0.0)
			colv := color
			if v[3]&recast.RC_BORDER_VERTEX != 0 {
				colv = DuRGBA(255, 255, 255, a)
				off = ch * 2
			}

			fx := float32(orig[0]) + float32(v[0])*cs
			fy := float32(orig[1]) + float32(v[1]+1+(i&1))*ch + off
			fz := float32(orig[2]) + float32(v[2])*cs
			dd.Vertex1(fx, fy, fz, colv)
		}
	}
	dd.End()
}

func DuDebugDrawPolyMesh(dd DuDebugDraw, mesh *recast.RcPolyMesh) {
	if dd == nil {
		return
	}

	nvp := mesh.Nvp
	cs := float32(mesh.Cs)
	ch := float32(mesh.Ch)
	orig := mesh.Bmin

	dd.Begin(DU_DRAW_TRIS)

	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		area := mesh.Areas[i]

		var color int
		if area == recast.RC_WALKABLE_AREA {
			color = DuRGBA(0, 192, 255, 64)
		} else if area == recast.RC_NULL_AREA {
			color = DuRGBA(0, 0, 0, 64)
		} else {
			color = dd.AreaToCol(int(area))
		}

		vi := make([]int, 3)
		for j := int32(2); j < nvp; j++ {
			if p[j] == recast.RC_MESH_NULL_IDX {
				break
			}
			vi[0] = int(p[0])
			vi[1] = int(p[j-1])
			vi[2] = int(p[j])
			for k := 0; k < 3; k++ {
				v := mesh.Verts[vi[k]*3:]
				x := float32(orig[0]) + float32(v[0])*cs
				y := float32(orig[1]) + float32(v[1]+1)*ch
				z := float32(orig[2]) + float32(v[2])*cs
				dd.Vertex1(x, y, z, color)
			}
		}
	}
	dd.End()

	// Draw neighbours edges
	coln := DuRGBA(0, 48, 64, 32)
	dd.Begin(DU_DRAW_LINES, 1.5)
	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		for j := int32(0); j < nvp; j++ {
			if p[j] == recast.RC_MESH_NULL_IDX {
				break
			}
			if p[nvp+j]&0x8000 != 0 {
				continue
			}
			nj := j + 1
			if j+1 >= nvp || p[j+1] == recast.RC_MESH_NULL_IDX {
				nj = 0
			}
			vi := []int{int(p[j]), int(p[nj])}

			for k := 0; k < 2; k++ {
				v := mesh.Verts[vi[k]*3:]
				x := float32(orig[0]) + float32(v[0])*cs
				y := float32(orig[1]) + float32(v[1]+1)*ch + 0.1
				z := float32(orig[2]) + float32(v[2])*cs
				dd.Vertex1(x, y, z, coln)
			}
		}
	}
	dd.End()

	// Draw boundary edges
	colb := DuRGBA(0, 48, 64, 220)
	dd.Begin(DU_DRAW_LINES, 2.5)
	for i := int32(0); i < mesh.Npolys; i++ {
		p := mesh.Polys[i*nvp*2:]
		for j := int32(0); j < nvp; j++ {
			if p[j] == recast.RC_MESH_NULL_IDX {
				break
			}
			if (p[nvp+j] & 0x8000) == 0 {
				continue
			}
			nj := j + 1
			if j+1 >= nvp || p[j+1] == recast.RC_MESH_NULL_IDX {
				nj = 0
			}
			vi := [2]int{int(p[j]), int(p[nj])}

			col := colb
			if (p[nvp+j] & 0xf) != 0xf {
				col = DuRGBA(255, 255, 255, 128)
			}

			for k := 0; k < 2; k++ {
				v := mesh.Verts[vi[k]*3:]
				x := float32(orig[0]) + float32(v[0])*cs
				y := float32(orig[1]) + float32(v[1]+1)*ch + 0.1
				z := float32(orig[2]) + float32(v[2])*cs
				dd.Vertex1(x, y, z, col)
			}
		}
	}
	dd.End()

	dd.Begin(DU_DRAW_POINTS, 3.0)
	colv := DuRGBA(0, 0, 0, 220)
	for i := int32(0); i < mesh.Nverts; i++ {
		v := mesh.Verts[i*3:]
		x := float32(orig[0]) + float32(v[0])*cs
		y := float32(orig[1]) + float32(v[1]+1)*ch + 0.1
		z := float32(orig[2]) + float32(v[2])*cs
		dd.Vertex1(x, y, z, colv)
	}
	dd.End()
}

func DuDebugDrawPolyMeshDetail(dd DuDebugDraw, dmesh *recast.RcPolyMeshDetail) {
	if dd == nil {
		return
	}

	dd.Begin(DU_DRAW_TRIS)

	for i := int32(0); i < dmesh.Nmeshes; i++ {
		m := dmesh.Meshes[i*4:]
		bverts := m[0]
		btris := m[2]
		ntris := m[3]
		verts := dmesh.Verts[bverts*3:]
		tris := dmesh.Tris[btris*4:]

		color := DuIntToCol(int(i), 192)

		for j := uint32(0); j < ntris; j++ {
			dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, tris[j*4+0])), color)
			dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, tris[j*4+1])), color)
			dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, tris[j*4+2])), color)
		}
	}
	dd.End()

	// Internal edges.
	dd.Begin(DU_DRAW_LINES, 1.0)
	coli := DuRGBA(0, 0, 0, 64)
	for i := int32(0); i < dmesh.Nmeshes; i++ {
		m := dmesh.Meshes[i*4:]
		bverts := m[0]
		btris := m[2]
		ntris := m[3]
		verts := dmesh.Verts[bverts*3:]
		tris := dmesh.Tris[btris*4:]

		for j := uint32(0); j < ntris; j++ {
			t := tris[j*4:]
			k := 0
			kp := 2
			for k < 3 {
				ef := (t[3] >> (kp * 2)) & 0x3
				if ef == 0 {
					// Internal edge
					if t[kp] < t[k] {
						dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, t[kp])), coli)
						dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, t[k])), coli)
					}
				}
				kp = k
				k++
			}
		}
	}
	dd.End()

	// External edges.
	dd.Begin(DU_DRAW_LINES, 2.0)
	cole := DuRGBA(0, 0, 0, 64)
	for i := int32(0); i < dmesh.Nmeshes; i++ {
		m := dmesh.Meshes[i*4:]
		bverts := m[0]
		btris := m[2]
		ntris := m[3]
		verts := dmesh.Verts[bverts*3:]
		tris := dmesh.Tris[btris*4:]

		for j := uint32(0); j < ntris; j++ {
			t := tris[j*4:]
			k := 0
			kp := 2
			for k < 3 {
				ef := (t[3] >> (kp * 2)) & 0x3
				if ef != 0 {
					// Ext edge
					dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, t[kp])), cole)
					dd.Vertex(common.SliceTToSlice[float32, float32](common.GetVert3(verts, t[k])), cole)
				}
				kp = k
				k++
			}
		}
	}
	dd.End()

	dd.Begin(DU_DRAW_POINTS, 3.0)
	colv := DuRGBA(0, 0, 0, 64)
	for i := int32(0); i < dmesh.Nmeshes; i++ {
		m := dmesh.Meshes[i*4:]
		bverts := m[0]
		nverts := m[1]
		verts := dmesh.Verts[bverts*3:]
		for j := uint32(0); j < nverts; j++ {
			dd.Vertex(common.SliceTToSlice[float32, float32](verts[j*3:]), colv)
		}

	}
	dd.End()
}
