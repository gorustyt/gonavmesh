package debug_utils

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/detour"
	"github.com/gorustyt/gonavmesh/detour_tile_cache"
)

const (
	DU_DRAWNAVMESH_OFFMESHCONS = 0x01
	DU_DRAWNAVMESH_CLOSEDLIST  = 0x02
	DU_DRAWNAVMESH_COLOR_TILES = 0x04
)

func distancePtLine2d(pt, p, q []float32) float32 {
	pqx := q[0] - p[0]
	pqz := q[2] - p[2]
	dx := pt[0] - p[0]
	dz := pt[2] - p[2]
	d := pqx*pqx + pqz*pqz
	t := pqx*dx + pqz*dz
	if d != 0 {
		t /= d
	}
	dx = p[0] + t*pqx - pt[0]
	dz = p[2] + t*pqz - pt[2]
	return dx*dx + dz*dz
}

func drawPolyBoundaries(dd DuDebugDraw, tile *detour.DtMeshTile, col Colorb, linew float32, inner bool) {
	const thr = 0.01 * 0.01

	dd.Begin(DU_DRAW_LINES, linew)

	for i := int32(0); i < tile.Header.PolyCount; i++ {
		p := tile.Polys[i]

		if p.GetType() == detour.DT_POLYTYPE_OFFMESH_CONNECTION {
			continue
		}

		pd := tile.DetailMeshes[i]
		j := 0
		nj := int(p.VertCount)
		for ; j < int(nj); j++ {
			c := col
			if inner {
				if p.Neis[j] == 0 {
					continue
				}
				if p.Neis[j]&detour.DT_EXT_LINK != 0 {
					con := false
					for k := p.FirstLink; k != detour.DT_NULL_LINK; k = tile.Links[k].Next {
						if int(tile.Links[k].Edge) == j {
							con = true
							break
						}
					}
					if con {
						c = DuRGBA(255, 255, 255, 48)
					} else {
						c = DuRGBA(0, 0, 0, 48)
					}

				} else {
					c = DuRGBA(0, 48, 64, 32)
				}

			} else {
				if p.Neis[j] != 0 {
					continue
				}
			}

			v0 := tile.Verts[p.Verts[j]*3:]
			v1 := tile.Verts[p.Verts[(j+1)%nj]*3:]

			// Draw detail mesh edges which align with the actual poly edge.
			// This is really slow.
			for k := 0; k < int(pd.TriCount); k++ {
				t := tile.DetailTris[(int(pd.TriBase)+k)*4:]
				tv := make([][]float32, 3)
				for m := 0; m < 3; m++ {
					if t[m] < p.VertCount {
						tv[m] = tile.Verts[p.Verts[t[m]]*3:]
					} else {
						tv[m] = tile.DetailVerts[(int(pd.VertBase)+(int(t[m])-int(p.VertCount)))*3:]
					}

				}
				m := 0
				n := 2
				for m < 3 {
					if (detour.DtGetDetailTriEdgeFlags(t[3], int32(n)) & detour.DT_DETAIL_EDGE_BOUNDARY) == 0 {
						n = m
						m++
						continue
					}

					if distancePtLine2d(tv[n], v0, v1) < thr && distancePtLine2d(tv[m], v0, v1) < thr {
						dd.Vertex(common.SliceTToSlice[float32, float32](tv[n]), c)
						dd.Vertex(common.SliceTToSlice[float32, float32](tv[m]), c)
					}
					n = m
					m++
				}
			}
		}
	}
	dd.End()
}

func DrawMeshTile(dd DuDebugDraw, mesh detour.IDtNavMesh, query detour.NavMeshQuery, tile *detour.DtMeshTile, flags int) {
	base := mesh.GetPolyRefBase(tile)

	tileNum := mesh.DecodePolyIdTile(base)
	tileColor := DuIntToCol(int(tileNum), 128)

	dd.DepthMask(false)

	dd.Begin(DU_DRAW_TRIS)
	for i := int32(0); i < tile.Header.PolyCount; i++ {
		p := tile.Polys[i]
		if p.GetType() == detour.DT_POLYTYPE_OFFMESH_CONNECTION {
			continue
		} // Skip off-mesh links.

		pd := tile.DetailMeshes[i]

		var col Colorb
		if query != nil && query.IsInClosedList(base|detour.DtPolyRef(i)) {
			col = DuRGBA(255, 196, 0, 64)
		} else {
			if flags&DU_DRAWNAVMESH_COLOR_TILES != 0 {
				col = tileColor
			} else {
				col = DuTransCol(dd.AreaToCol(int(p.GetArea())), 64)
			}

		}

		for j := uint32(0); j < uint32(pd.TriCount); j++ {
			t := tile.DetailTris[(pd.TriBase+j)*4:]
			for k := 0; k < 3; k++ {
				if t[k] < p.VertCount {
					dd.Vertex(common.SliceTToSlice[float32, float32](tile.Verts[p.Verts[t[k]]*3:]), col)
				} else {
					dd.Vertex(common.SliceTToSlice[float32, float32](tile.DetailVerts[(pd.VertBase+uint32(t[k])-uint32(p.VertCount))*3:]), col)
				}

			}
		}
	}
	dd.End()

	// Draw inter poly boundaries
	drawPolyBoundaries(dd, tile, DuRGBA(0, 48, 64, 32), 1.5, true)

	// Draw outer poly boundaries
	drawPolyBoundaries(dd, tile, DuRGBA(0, 48, 64, 220), 2.5, false)

	if flags&DU_DRAWNAVMESH_OFFMESHCONS != 0 {
		dd.Begin(DU_DRAW_LINES, 2.0)
		for i := int32(0); i < tile.Header.PolyCount; i++ {
			p := tile.Polys[i]
			if p.GetType() != detour.DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			} // Skip regular polys.

			var col, col2 Colorb
			if query != nil && query.IsInClosedList(base|detour.DtPolyRef(i)) {
				col = DuRGBA(255, 196, 0, 220)
			} else {
				col = DuDarkenCol(DuTransCol(dd.AreaToCol(int(p.GetArea())), 220))
			}

			con := tile.OffMeshCons[i-tile.Header.OffMeshBase]
			va := common.SliceTToSlice[float32, float32](tile.Verts[p.Verts[0]*3:])
			vb := common.SliceTToSlice[float32, float32](tile.Verts[p.Verts[1]*3:])

			// Check to see if start and end end-points have links.
			startSet := false
			endSet := false
			for k := p.FirstLink; k != detour.DT_NULL_LINK; k = tile.Links[k].Next {
				if tile.Links[k].Edge == 0 {
					startSet = true
				}

				if tile.Links[k].Edge == 1 {
					endSet = true
				}

			}

			// End points and their on-mesh locations.
			dd.Vertex1(va[0], va[1], va[2], col)
			dd.Vertex1(float32(con.Pos[0]), float32(con.Pos[1]), float32(con.Pos[2]), col)
			col2 = DuRGBA(220, 32, 16, 196)
			if startSet {
				col2 = col
			}
			DuAppendCircle(dd, float32(con.Pos[0]), float32(con.Pos[1]+0.1), float32(con.Pos[2]), float32(con.Rad), col2)

			dd.Vertex1(vb[0], vb[1], vb[2], col)
			dd.Vertex1(float32(con.Pos[3]), float32(con.Pos[4]), float32(con.Pos[5]), col)
			col2 = DuRGBA(220, 32, 16, 196)
			if endSet {
				col2 = col
			}
			DuAppendCircle(dd, float32(con.Pos[3]), float32(con.Pos[4]+0.1), float32(con.Pos[5]), float32(con.Rad), col2)

			// End point vertices.
			dd.Vertex1(float32(con.Pos[0]), float32(con.Pos[1]), float32(con.Pos[2]), DuRGBA(0, 48, 64, 196))
			dd.Vertex1(float32(con.Pos[0]), float32(con.Pos[1]+0.2), float32(con.Pos[2]), DuRGBA(0, 48, 64, 196))

			dd.Vertex1(float32(con.Pos[3]), float32(con.Pos[4]), float32(con.Pos[5]), DuRGBA(0, 48, 64, 196))
			dd.Vertex1(float32(con.Pos[3]), float32(con.Pos[4]+0.2), float32(con.Pos[5]), DuRGBA(0, 48, 64, 196))

			// Connection arc.
			a := float32(0.0)
			if con.Flags&1 != 0 {
				a = 0.6
			}
			DuAppendArc(dd, float32(con.Pos[0]), float32(con.Pos[1]), float32(con.Pos[2]), float32(con.Pos[3]), float32(con.Pos[4]), float32(con.Pos[5]), 0.25,
				a, 0.6, col)
		}
		dd.End()
	}

	vcol := DuRGBA(0, 0, 0, 196)
	dd.Begin(DU_DRAW_POINTS, 3.0)
	for i := int32(0); i < tile.Header.VertCount; i++ {
		v := tile.Verts[i*3:]
		dd.Vertex1(float32(v[0]), float32(v[1]), float32(v[2]), vcol)
	}
	dd.End()

	dd.DepthMask(true)
}

func DuDebugDrawNavMesh(dd DuDebugDraw, mesh detour.IDtNavMesh, flags int) {
	if dd == nil {
		return
	}

	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile.Header == nil {
			continue
		}
		DrawMeshTile(dd, mesh, nil, tile, flags)
	}
}

func DuDebugDrawNavMeshWithClosedList(dd DuDebugDraw, mesh detour.IDtNavMesh, query detour.NavMeshQuery, flags int) {
	if dd == nil {
		return
	}

	var q detour.NavMeshQuery
	if flags&DU_DRAWNAVMESH_CLOSEDLIST != 0 {
		q = query
	}

	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile.Header == nil {
			continue
		}
		DrawMeshTile(dd, mesh, q, tile, flags)
	}
}

func DuDebugDrawNavMeshNodes(dd DuDebugDraw, query detour.NavMeshQuery) {
	if dd == nil {
		return
	}

	pool := query.GetNodePool()
	if pool != nil {
		off := float32(0.5)
		dd.Begin(DU_DRAW_POINTS, 4.0)
		for i := int32(0); i < pool.GetHashSize(); i++ {
			for j := pool.GetFirst(int(i)); j != detour.DT_NULL_IDX; j = pool.GetNext(int(j)) {
				node := pool.GetNodeAtIdx(int(j) + 1)
				if node == nil {
					continue
				}
				dd.Vertex1(float32(node.Pos[0]), float32(node.Pos[1])+off, float32(node.Pos[2]), DuRGBA(255, 192, 0, 255))
			}
		}
		dd.End()

		dd.Begin(DU_DRAW_LINES, 2.0)
		for i := int32(0); i < pool.GetHashSize(); i++ {
			for j := pool.GetFirst(int(i)); j != detour.DT_NULL_IDX; j = pool.GetNext(int(j)) {
				node := pool.GetNodeAtIdx(int(j) + 1)
				if node == nil {
					continue
				}
				if node.Pidx == 0 {
					continue
				}
				parent := pool.GetNodeAtIdx(int(node.Pidx))
				if parent == nil {
					continue
				}
				dd.Vertex1(float32(node.Pos[0]), float32(node.Pos[1])+off, float32(node.Pos[2]), DuRGBA(255, 192, 0, 128))
				dd.Vertex1(float32(parent.Pos[0]), float32(parent.Pos[1])+off, float32(parent.Pos[2]), DuRGBA(255, 192, 0, 128))
			}
		}
		dd.End()
	}
}

func DrawMeshTileBVTree(dd DuDebugDraw, tile *detour.DtMeshTile) {
	// Draw BV nodes.
	cs := 1.0 / float32(tile.Header.BvQuantFactor)
	dd.Begin(DU_DRAW_LINES, 1.0)
	for i := 0; i < int(tile.Header.BvNodeCount); i++ {
		n := tile.BvTree[i]
		if n.I < 0 {
			continue
		} // Leaf indices are positive.

		DuAppendBoxWire(dd, float32(tile.Header.Bmin[0])+float32(n.Bmin[0])*cs,
			float32(tile.Header.Bmin[1])+float32(n.Bmin[1])*cs,
			float32(tile.Header.Bmin[2])+float32(n.Bmin[2])*cs,
			float32(tile.Header.Bmin[0])+float32(n.Bmin[0])*cs,
			float32(tile.Header.Bmin[1])+float32(n.Bmin[1])*cs,
			float32(tile.Header.Bmin[2])+float32(n.Bmin[2])*cs,
			DuRGBA(255, 255, 255, 128))
	}
	dd.End()
}

func DuDebugDrawNavMeshBVTree(dd DuDebugDraw, mesh detour.IDtNavMesh) {
	if dd == nil {
		return
	}

	for i := 0; i < int(mesh.GetMaxTiles()); i++ {
		tile := mesh.GetTile(i)
		if tile.Header == nil {
			continue
		}
		DrawMeshTileBVTree(dd, tile)
	}
}

func drawMeshTilePortal(dd DuDebugDraw, tile *detour.DtMeshTile) {
	// Draw portals
	padx := float32(0.04)
	pady := tile.Header.WalkableClimb

	dd.Begin(DU_DRAW_LINES, 2.0)

	for side := 0; side < 8; side++ {
		m := detour.DT_EXT_LINK | side

		for i := 0; i < int(tile.Header.PolyCount); i++ {
			poly := tile.Polys[i]

			// Create new links.
			nv := int(poly.VertCount)
			for j := 0; j < nv; j++ {
				// Skip edges which do not point to the right side.
				if int(poly.Neis[j]) != m {
					continue
				}

				// Create new links
				va := tile.Verts[poly.Verts[j]*3:]
				vb := tile.Verts[poly.Verts[(j+1)%nv]*3:]

				if side == 0 || side == 4 {
					col := DuRGBA(128, 0, 128, 128)
					if side == 0 {
						col = DuRGBA(128, 0, 0, 128)
					}
					x := float32(va[0]) + padx
					if side == 0 {
						x = float32(va[0]) - padx
					}
					dd.Vertex1(x, float32(va[1]-pady), float32(va[2]), col)
					dd.Vertex1(x, float32(va[1]+pady), float32(va[2]), col)

					dd.Vertex1(x, float32(va[1]+pady), float32(va[2]), col)
					dd.Vertex1(x, float32(vb[1]+pady), float32(vb[2]), col)

					dd.Vertex1(x, float32(vb[1]+pady), float32(vb[2]), col)
					dd.Vertex1(x, float32(vb[1]-pady), float32(vb[2]), col)

					dd.Vertex1(x, float32(vb[1]-pady), float32(vb[2]), col)
					dd.Vertex1(x, float32(va[1]-pady), float32(va[2]), col)
				} else if side == 2 || side == 6 {
					col := DuRGBA(0, 128, 128, 128)
					if side == 2 {
						col = DuRGBA(0, 128, 0, 128)
					}
					z := va[2] + padx
					if side == 2 {
						z = va[2] - padx
					}
					dd.Vertex1(float32(va[0]), float32(va[1]-pady), z, col)
					dd.Vertex1(float32(va[0]), float32(va[1]+pady), z, col)

					dd.Vertex1(float32(va[0]), float32(va[1]+pady), z, col)
					dd.Vertex1(float32(vb[0]), float32(vb[1]+pady), z, col)

					dd.Vertex1(float32(vb[0]), float32(vb[1]+pady), z, col)
					dd.Vertex1(float32(vb[0]), float32(vb[1]-pady), z, col)

					dd.Vertex1(float32(vb[0]), float32(vb[1]-pady), z, col)
					dd.Vertex1(float32(va[0]), float32(va[1]-pady), z, col)
				}

			}
		}
	}

	dd.End()
}

func DuDebugDrawNavMeshPortals(dd DuDebugDraw, mesh detour.IDtNavMesh) {
	if dd == nil {
		return
	}

	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile.Header == nil {
			continue
		}
		drawMeshTilePortal(dd, tile)
	}
}

func DuDebugDrawNavMeshPolysWithFlags(dd DuDebugDraw, mesh detour.IDtNavMesh,
	polyFlags int, col Colorb) {
	if dd == nil {
		return
	}

	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile.Header == nil {
			continue
		}
		base := mesh.GetPolyRefBase(tile)

		for j := int32(0); j < tile.Header.PolyCount; j++ {
			p := tile.Polys[j]
			if (int(p.Flags) & polyFlags) == 0 {
				continue
			}
			DuDebugDrawNavMeshPoly(dd, mesh, base|detour.DtPolyRef(j), col)
		}
	}
}

func DuDebugDrawNavMeshPoly(dd DuDebugDraw, mesh detour.IDtNavMesh, ref detour.DtPolyRef, col Colorb) {
	if dd == nil {
		return
	}

	tile, poly, status := mesh.GetTileAndPolyByRef(ref)
	if status.DtStatusFailed() {
		return
	}

	dd.DepthMask(false)

	c := DuTransCol(col, 64)
	ip := tile.GetIndexByPloy(poly)
	if poly.GetType() == detour.DT_POLYTYPE_OFFMESH_CONNECTION {
		con := tile.OffMeshCons[ip-int(tile.Header.OffMeshBase)]

		dd.Begin(DU_DRAW_LINES, 2.0)

		// Connection arc.
		as0 := 00.0
		if con.Flags&1 != 0 {
			as0 = 0.6
		}
		DuAppendArc(dd, float32(con.Pos[0]), float32(con.Pos[1]), float32(con.Pos[2]), float32(con.Pos[3]), float32(con.Pos[4]), float32(con.Pos[5]), 0.25,
			float32(as0), 0.6, c)

		dd.End()
	} else {
		pd := tile.DetailMeshes[ip]

		dd.Begin(DU_DRAW_TRIS)
		for i := uint32(0); i < uint32(pd.TriCount); i++ {
			t := tile.DetailTris[(pd.TriBase+i)*4:]
			for j := 0; j < 3; j++ {
				if t[j] < poly.VertCount {
					dd.Vertex(common.SliceTToSlice[float32, float32](tile.Verts[poly.Verts[t[j]]*3:]), c)
				} else {
					dd.Vertex(common.SliceTToSlice[float32, float32](tile.DetailVerts[(pd.VertBase+uint32(t[j])-uint32(poly.VertCount))*3:]), c)
				}

			}
		}
		dd.End()
	}

	dd.DepthMask(true)

}

func DebugDrawTileCachePortals(dd DuDebugDraw, layer *detour_tile_cache.DtTileCacheLayer, cs, ch float32) {
	w := layer.Header.Width
	h := layer.Header.Height
	bmin := layer.Header.Bmin

	// Portals
	pcol := DuRGBA(255, 255, 255, 255)

	segs := [4 * 4]int{0, 0, 0, 1, 0, 1, 1, 1, 1, 1, 1, 0, 1, 0, 0, 0}

	// Layer portals
	dd.Begin(DU_DRAW_LINES, 2.0)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			idx := x + y*w
			lh := layer.Heights[idx]
			if lh == 0xff {
				continue
			}

			for dir := 0; dir < 4; dir++ {
				if layer.Cons[idx]&(1<<(dir+4)) != 0 {
					seg := segs[dir*4:]
					ax := float32(bmin[0]) + float32(x+seg[0])*cs
					ay := float32(bmin[1]) + float32(lh+2)*ch
					az := float32(bmin[2]) + float32(y+seg[1])*cs
					bx := float32(bmin[0]) + float32(x+seg[2])*cs
					by := float32(bmin[1]) + float32(lh+2)*ch
					bz := float32(bmin[2]) + float32(y+seg[3])*cs
					dd.Vertex1(ax, ay, az, pcol)
					dd.Vertex1(bx, by, bz, pcol)
				}
			}
		}
	}
	dd.End()
}

func DuDebugDrawTileCacheLayerAreas(dd DuDebugDraw, layer *detour_tile_cache.DtTileCacheLayer, cs, ch float32) {
	w := layer.Header.Width
	h := layer.Header.Height
	bmin := layer.Header.Bmin
	bmax := layer.Header.Bmax
	idx := layer.Header.Tlayer

	color := DuIntToCol(idx+1, 255)

	// Layer bounds
	lbmin := make([]float32, 3)
	lbmax := make([]float32, 3)
	lbmin[0] = float32(bmin[0]) + float32(layer.Header.Minx)*cs
	lbmin[1] = float32(bmin[1])
	lbmin[2] = float32(bmin[2]) + float32(layer.Header.Miny)*cs
	lbmax[0] = float32(bmin[0]) + float32(layer.Header.Maxx+1)*cs
	lbmax[1] = float32(bmax[1])
	lbmax[2] = float32(bmin[2]) + float32(layer.Header.Maxy+1)*cs
	DuDebugDrawBoxWire(dd, lbmin[0], lbmin[1], lbmin[2], lbmax[0], lbmax[1], lbmax[2], DuTransCol(color, 128), 2.0)

	// Layer height
	dd.Begin(DU_DRAW_QUADS)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			lidx := x + y*w
			lh := layer.Heights[lidx]
			if lh == 0xff {
				continue
			}

			area := layer.Areas[lidx]
			var col Colorb
			if area == 63 {
				col = DuLerpCol(color, DuRGBA(0, 192, 255, 64), 32)
			} else if area == 0 {
				col = DuLerpCol(color, DuRGBA(0, 0, 0, 64), 32)
			} else {
				col = DuLerpCol(color, dd.AreaToCol(area), 32)
			}

			fx := float32(bmin[0]) + float32(x)*cs
			fy := float32(bmin[1]) + float32(lh+1)*ch
			fz := float32(bmin[2]) + float32(y)*cs

			dd.Vertex1(float32(fx), fy, float32(fz), col)
			dd.Vertex1(float32(fx), float32(fy), float32(fz+cs), col)
			dd.Vertex1(float32(fx+cs), float32(fy), float32(fz+cs), col)
			dd.Vertex1(float32(fx+cs), float32(fy), float32(fz), col)
		}
	}
	dd.End()

	DebugDrawTileCachePortals(dd, layer, cs, ch)
}

func DuDebugDrawTileCacheLayerRegions(dd DuDebugDraw, layer *detour_tile_cache.DtTileCacheLayer, cs, ch float32) {
	w := layer.Header.Width
	h := layer.Header.Height
	bmin := layer.Header.Bmin
	bmax := layer.Header.Bmax
	idx := layer.Header.Tlayer

	color := DuIntToCol(idx+1, 255)

	// Layer bounds
	lbmin := make([]float32, 3)
	lbmax := make([]float32, 3)
	lbmin[0] = float32(bmin[0]) + float32(layer.Header.Minx)*cs
	lbmin[1] = float32(bmin[1])
	lbmin[2] = float32(bmin[2]) + float32(layer.Header.Miny)*cs
	lbmax[0] = float32(bmin[0]) + float32(layer.Header.Maxx+1)*cs
	lbmax[1] = float32(bmax[1])
	lbmax[2] = float32(bmin[2]) + float32(layer.Header.Maxy+1)*cs
	DuDebugDrawBoxWire(dd, lbmin[0], lbmin[1], lbmin[2], lbmax[0], lbmax[1], lbmax[2], DuTransCol(color, 128), 2.0)

	// Layer height
	dd.Begin(DU_DRAW_QUADS)
	for y := 0; y < h; y++ {
		for x := 0; x < w; x++ {
			lidx := x + y*w
			lh := layer.Heights[lidx]
			if lh == 0xff {
				continue
			}
			reg := layer.Regs[lidx]

			col := DuLerpCol(color, DuIntToCol(reg, 255), 192)

			fx := float32(bmin[0]) + float32(x)*cs
			fy := float32(bmin[1]) + float32(lh+1)*ch
			fz := float32(bmin[2]) + float32(y)*cs

			dd.Vertex1(fx, fy, fz, col)
			dd.Vertex1(fx, fy, fz+cs, col)
			dd.Vertex1(fx+cs, fy, fz+cs, col)
			dd.Vertex1(fx+cs, fy, fz, col)
		}
	}
	dd.End()

	DebugDrawTileCachePortals(dd, layer, cs, ch)
}

/*struct dtTileCacheContour
{
	int nverts;
	unsigned char* verts;
	unsigned char reg;
	unsigned char area;
};

struct dtTileCacheContourSet
{
	int nconts;
	dtTileCacheContour* conts;
};*/

func DuDebugDrawTileCacheContours(dd DuDebugDraw, lcset *detour_tile_cache.DtTileCacheContourSet,
	orig []float32, cs, ch float32) {
	if dd == nil {
		return
	}

	a := 255 // (unsigned char)(alpha*255.0f);

	offs := [2 * 4]float32{-1, 0, 0, 1, 1, 0, 0, -1}

	dd.Begin(DU_DRAW_LINES, 2.0)

	for i := 0; i < lcset.Nconts; i++ {
		c := lcset.Conts[i]

		color := DuIntToCol(i, a)

		for j := 0; j < c.Nverts; j++ {
			k := (j + 1) % c.Nverts
			va := c.Verts[j*4:]
			vb := c.Verts[k*4:]
			ax := orig[0] + float32(va[0])*cs
			ay := orig[1] + float32(va[1]+1+(i&1))*ch
			az := orig[2] + float32(va[2])*cs
			bx := orig[0] + float32(vb[0])*cs
			by := orig[1] + float32(vb[1]+1+(i&1))*ch
			bz := orig[2] + float32(vb[2])*cs
			col := color
			if (va[3] & 0xf) != 0xf {
				// Portal segment
				col = DuRGBA(255, 255, 255, 128)
				d := va[3] & 0xf

				cx := (ax + bx) * 0.5
				cy := (ay + by) * 0.5
				cz := (az + bz) * 0.5

				dx := cx + offs[d*2+0]*2*cs
				dy := cy
				dz := cz + offs[d*2+1]*2*cs

				dd.Vertex1(cx, cy, cz, DuRGBA(255, 0, 0, 255))
				dd.Vertex1(dx, dy, dz, DuRGBA(255, 0, 0, 255))
			}

			DuAppendArrow(dd, ax, ay, az, bx, by, bz, 0.0, cs*0.5, col)
		}
	}
	dd.End()

	dd.Begin(DU_DRAW_POINTS, 4.0)

	for i := 0; i < lcset.Nconts; i++ {
		c := lcset.Conts[i]
		var color Colorb

		for j := 0; j < c.Nverts; j++ {
			va := c.Verts[j*4:]

			color = DuDarkenCol(DuIntToCol(i, a))
			if va[3]&0x80 != 0 {
				// Border vertex
				color = DuRGBA(255, 0, 0, 255)
			}

			fx := orig[0] + float32(va[0])*cs
			fy := orig[1] + float32(va[1]+1+(i&1))*ch
			fz := orig[2] + float32(va[2])*cs
			dd.Vertex1(fx, fy, fz, color)
		}
	}
	dd.End()
}

func DuDebugDrawTileCachePolyMesh(dd DuDebugDraw, lmesh *detour_tile_cache.DtTileCachePolyMesh,
	orig []float32, cs, ch float32) {
	if dd == nil {
		return
	}

	nvp := lmesh.Nvp

	offs := [2 * 4]float32{-1, 0, 0, 1, 1, 0, 0, -1}

	dd.Begin(DU_DRAW_TRIS)

	for i := 0; i < lmesh.Npolys; i++ {
		p := lmesh.Polys[i*nvp*2:]
		area := lmesh.Areas[i]

		var color Colorb
		if area == detour_tile_cache.DT_TILECACHE_WALKABLE_AREA {
			color = DuRGBA(0, 192, 255, 64)
		} else if area == detour_tile_cache.DT_TILECACHE_NULL_AREA {
			color = DuRGBA(0, 0, 0, 64)
		} else {
			color = dd.AreaToCol(area)
		}

		vi := make([]int, 3)
		for j := 2; j < nvp; j++ {
			if p[j] == detour_tile_cache.DT_TILECACHE_NULL_IDX {
				break
			}
			vi[0] = p[0]
			vi[1] = p[j-1]
			vi[2] = p[j]
			for k := 0; k < 3; k++ {
				v := lmesh.Verts[vi[k]*3:]
				x := orig[0] + float32(v[0])*cs
				y := orig[1] + float32(v[1]+1)*ch
				z := orig[2] + float32(v[2])*cs
				dd.Vertex1(x, y, z, color)
			}
		}
	}
	dd.End()

	// Draw neighbours edges
	coln := DuRGBA(0, 48, 64, 32)
	dd.Begin(DU_DRAW_LINES, 1.5)
	for i := 0; i < lmesh.Npolys; i++ {
		p := lmesh.Polys[i*nvp*2:]
		for j := 0; j < nvp; j++ {
			if p[j] == detour_tile_cache.DT_TILECACHE_NULL_IDX {
				break
			}
			if p[nvp+j]&0x8000 != 0 {
				continue
			}
			nj := j + 1
			if j+1 >= nvp || p[j+1] == detour_tile_cache.DT_TILECACHE_NULL_IDX {
				nj = 0
			}
			vi := [2]int{p[j], p[nj]}

			for k := 0; k < 2; k++ {
				v := lmesh.Verts[vi[k]*3:]
				x := orig[0] + float32(v[0])*cs
				y := orig[1] + float32(v[1]+1)*ch + 0.1
				z := orig[2] + float32(v[2])*cs
				dd.Vertex1(x, y, z, coln)
			}
		}
	}
	dd.End()

	// Draw boundary edges
	colb := DuRGBA(0, 48, 64, 220)
	dd.Begin(DU_DRAW_LINES, 2.5)
	for i := 0; i < lmesh.Npolys; i++ {
		p := lmesh.Polys[i*nvp*2:]
		for j := 0; j < nvp; j++ {
			if p[j] == detour_tile_cache.DT_TILECACHE_NULL_IDX {
				break
			}
			if (p[nvp+j] & 0x8000) == 0 {
				continue
			}
			nj := j + 1
			if j+1 >= nvp || p[j+1] == detour_tile_cache.DT_TILECACHE_NULL_IDX {
				nj = 0
			}
			vi := [2]int{p[j], p[nj]}

			col := colb
			if (p[nvp+j] & 0xf) != 0xf {
				va := lmesh.Verts[vi[0]*3:]
				vb := lmesh.Verts[vi[1]*3:]

				ax := orig[0] + float32(va[0])*cs
				ay := orig[1] + float32(va[1]+1+(i&1))*ch
				az := orig[2] + float32(va[2])*cs
				bx := orig[0] + float32(vb[0])*cs
				by := orig[1] + float32(vb[1]+1+(i&1))*ch
				bz := orig[2] + float32(vb[2])*cs

				cx := (ax + bx) * 0.5
				cy := (ay + by) * 0.5
				cz := (az + bz) * 0.5

				d := p[nvp+j] & 0xf

				dx := cx + offs[d*2+0]*2*cs
				dy := cy
				dz := cz + offs[d*2+1]*2*cs

				dd.Vertex1(cx, cy, cz, DuRGBA(255, 0, 0, 255))
				dd.Vertex1(dx, dy, dz, DuRGBA(255, 0, 0, 255))

				col = DuRGBA(255, 255, 255, 128)
			}

			for k := 0; k < 2; k++ {
				v := lmesh.Verts[vi[k]*3:]
				x := orig[0] + float32(v[0])*cs
				y := orig[1] + float32(v[1]+1)*ch + 0.1
				z := orig[2] + float32(v[2])*cs
				dd.Vertex1(x, y, z, col)
			}
		}
	}
	dd.End()

	dd.Begin(DU_DRAW_POINTS, 3.0)
	colv := DuRGBA(0, 0, 0, 220)
	for i := 0; i < lmesh.Nverts; i++ {
		v := lmesh.Verts[i*3:]
		x := orig[0] + float32(v[0])*cs
		y := orig[1] + float32(v[1]+1)*ch + 0.1
		z := orig[2] + float32(v[2])*cs
		dd.Vertex1(x, y, z, colv)
	}
	dd.End()
}
