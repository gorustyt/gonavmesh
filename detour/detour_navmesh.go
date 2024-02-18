package detour

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/recast"
	"math"
)

/// @{
/// @name Initialization and Tile Management

// / Initializes the navigation mesh for tiled use.
// /  @param[in]	params		Initialization parameters.
// / @return The status flags for the operation.
func NewDtNavMeshWithParams(params *NavMeshParams) (m IDtNavMesh, status DtStatus) {
	mesh := &DtNavMesh{}
	mesh.m_orig = params.Orig
	mesh.m_tileWidth = params.TileWidth
	mesh.m_tileHeight = params.TileHeight
	mesh.m_params = params

	// Init tiles
	mesh.m_maxTiles = params.MaxTiles
	mesh.m_tileLutSize = int32(common.NextPow2(uint32(params.MaxTiles) / 4))
	if mesh.m_tileLutSize == 0 {
		mesh.m_tileLutSize = 1
	}
	mesh.m_tileLutMask = mesh.m_tileLutSize - 1
	var m_nextFree *DtMeshTile
	for i := mesh.m_maxTiles - 1; i >= 0; i-- {
		mesh.m_tiles[i].salt = 1
		mesh.m_tiles[i].Next = m_nextFree
		m_nextFree = mesh.m_tiles[i]
	}

	// Init ID generator values.
	if DT_POLYREF64 == 1 {
		mesh.m_tileBits = common.Ilog2(common.NextPow2(uint32(params.MaxTiles)))
		mesh.m_polyBits = common.Ilog2(common.NextPow2(uint32(params.MaxPolys)))
		// Only allow 31 salt bits, since the salt mask is calculated using 32bit uint and it will overflow.
		mesh.m_saltBits = common.Min(31, 32-mesh.m_tileBits-mesh.m_polyBits)

		if mesh.m_saltBits < 10 {
			return mesh, DT_FAILURE | DT_INVALID_PARAM
		}

	}
	return mesh, DT_SUCCESS
}

// / Initializes the navigation mesh for single tile use.
// /  @param[in]	data		Data of the new tile. (See: #dtCreateNavMeshData)
// /  @param[in]	dataSize	The data size of the new tile.
// /  @param[in]	flags		The tile flags. (See: #dtTileFlags)
// / @return The status flags for the operation.
// /  @see dtCreateNavMeshData
func NewDtNavMesh(data *recast.NavMeshData, flags int32) (m IDtNavMesh, result DtTileRef, status DtStatus) {
	// Make sure the data is in right format.
	header := data.Header
	if header.Magic != DT_NAVMESH_MAGIC {
		return m, result, DT_FAILURE | DT_WRONG_MAGIC
	}

	if header.Version != DT_NAVMESH_VERSION {
		return m, result, DT_FAILURE | DT_WRONG_VERSION
	}

	var params NavMeshParams
	params.Orig = header.Bmin
	params.TileWidth = header.Bmax[0] - header.Bmin[0]
	params.TileHeight = header.Bmax[2] - header.Bmin[2]
	params.MaxTiles = 1
	params.MaxPolys = header.PolyCount

	mesh, status := NewDtNavMeshWithParams(&params)
	if status.DtStatusFailed() {
		return mesh, result, status
	}
	result, status = mesh.AddTile(data, flags, 0)
	return mesh, result, status
}

// / @par
// /
// / @note The parameters are created automatically when the single tile
// / initialization is performed.
func (mesh *DtNavMesh) GetParams() *NavMeshParams {
	return mesh.m_params
}
func overlapSlabs(amin, amax, bmin, bmax []float32, px, py float32) bool {
	// Check for horizontal overlap.
	// The segment is shrunken a little so that slabs which touch
	// at end points are not connected.
	minx := common.Max(amin[0]+px, bmin[0]+px)
	maxx := common.Min(amax[0]-px, bmax[0]-px)
	if minx > maxx {
		return false
	}
	// Check vertical overlap.
	ad := (amax[1] - amin[1]) / (amax[0] - amin[0])
	ak := amin[1] - ad*amin[0]
	bd := (bmax[1] - bmin[1]) / (bmax[0] - bmin[0])
	bk := bmin[1] - bd*bmin[0]
	aminy := ad*minx + ak
	amaxy := ad*maxx + ak
	bminy := bd*minx + bk
	bmaxy := bd*maxx + bk
	dmin := bminy - aminy
	dmax := bmaxy - amaxy

	// Crossing segments always overlap.
	if dmin*dmax < 0 {
		return true
	}

	// Check for overlap at endpoints.
	thr := common.Sqr(py * 2)
	if dmin*dmin <= thr || dmax*dmax <= thr {
		return true
	}

	return false
}

func getSlabCoord(va []float32, side int32) float32 {
	if side == 0 || side == 4 {
		return va[0]
	} else if side == 2 || side == 6 {
		return va[2]
	}
	return 0
}

func calcSlabEndPoints(va, vb []float32, side int32) ([]float32, []float32) {
	bmin, bmax := make([]float32, 3), make([]float32, 3)
	if side == 0 || side == 4 {
		if va[2] < vb[2] {
			bmin[0] = va[2]
			bmin[1] = va[1]
			bmax[0] = vb[2]
			bmax[1] = vb[1]
		} else {
			bmin[0] = vb[2]
			bmin[1] = vb[1]
			bmax[0] = va[2]
			bmax[1] = va[1]
		}
	} else if side == 2 || side == 6 {
		if va[0] < vb[0] {
			bmin[0] = va[0]
			bmin[1] = va[1]
			bmax[0] = vb[0]
			bmax[1] = vb[1]
		} else {
			bmin[0] = vb[0]
			bmin[1] = vb[1]
			bmax[0] = va[0]
			bmax[1] = va[1]
		}
	}
	return bmin, bmax
}

func allocLink(tile *DtMeshTile) uint32 {
	if tile.linksFreeList == DT_NULL_LINK {
		return DT_NULL_LINK
	}

	link := tile.linksFreeList
	tile.linksFreeList = tile.Links[link].Next
	return link
}

func freeLink(tile *DtMeshTile, link uint32) {
	tile.Links[link].Next = tile.linksFreeList
	tile.linksFreeList = link
}

// ////////////////////////////////////////////////////////////////////////////////////////
func (mesh *DtNavMesh) FindConnectingPolys(va, vb []float32, tile *DtMeshTile, side int32, maxcon int32) (con []DtPolyRef, conarea []float32, n int32) {
	con = make([]DtPolyRef, 4)
	conarea = make([]float32, 4*2)
	if tile == nil {
		return
	}

	amin, amax := calcSlabEndPoints(va, vb, side)
	apos := getSlabCoord(va, side)

	// Remove links pointing to 'side' and compact the links array.

	m := uint16(DT_EXT_LINK | side)
	base := mesh.GetPolyRefBase(tile)

	for i := int32(0); i < tile.Header.PolyCount; i++ {
		poly := tile.Polys[i]
		nv := int32(poly.VertCount)
		for j := int32(0); j < nv; j++ {
			// Skip edges which do not point to the right side.
			if poly.Neis[j] != m {
				continue
			}

			vc := common.GetVert3(tile.Verts, poly.Verts[j])
			vd := common.GetVert3(tile.Verts, (j+1)%nv)
			bpos := getSlabCoord(vc, side)

			// Segments are not close enough.
			if common.Abs(apos-bpos) > 0.01 {
				continue
			}

			// Check if the segments touch.
			bmin, bmax := calcSlabEndPoints(vc, vd, side)

			if !overlapSlabs(amin, amax, bmin, bmax, 0.01, tile.Header.WalkableClimb) {
				continue
			}

			// Add return value.
			if n < maxcon {
				conarea[n*2+0] = common.Max(amin[0], bmin[0])
				conarea[n*2+1] = common.Min(amax[0], bmax[0])
				con[n] = DtPolyRef(int32(base) | i)
				n++
			}
			break
		}
	}
	return
}

func (mesh *DtNavMesh) GetPolyRefBase(tile *DtMeshTile) DtPolyRef {
	if tile == nil {
		return 0
	}
	it := mesh.getTileIndex(tile)
	return mesh.EncodePolyId(tile.salt, uint32(it), 0)
}

func (mesh *DtNavMesh) unconnectLinks(tile *DtMeshTile, target *DtMeshTile) {
	if tile == nil || target == nil {
		return
	}

	targetNum := mesh.DecodePolyIdTile(DtPolyRef(mesh.GetTileRef(target)))

	for i := int32(0); i < tile.Header.PolyCount; i++ {
		poly := tile.Polys[i]
		j := poly.FirstLink
		pj := uint32(DT_NULL_LINK)
		for j != DT_NULL_LINK {
			if mesh.DecodePolyIdTile(tile.Links[j].Ref) == targetNum {
				// Remove link.
				nj := tile.Links[j].Next
				if pj == DT_NULL_LINK {
					poly.FirstLink = nj
				} else {
					tile.Links[pj].Next = nj
				}

				freeLink(tile, j)
				j = nj
			} else {
				// Advance
				pj = j
				j = tile.Links[j].Next
			}
		}
	}
}
func (mesh *DtNavMesh) getTileIndex(tile *DtMeshTile) int {
	for index, v := range mesh.m_tiles {
		if v == tile {
			return index
		}
	}
	return -1
}
func (mesh *DtNavMesh) GetTileRef(tile *DtMeshTile) DtTileRef {
	if tile == nil {
		return DtTileRef(0)
	}
	it := mesh.getTileIndex(tile)
	return DtTileRef(mesh.EncodePolyId(tile.salt, uint32(it), 0))
}

func (mesh *DtNavMesh) GetTileRefAt(x, y, layer int32) DtTileRef {
	// Find tile based on hash.
	h := common.ComputeTileHash(x, y, mesh.m_tileLutMask)
	tile := mesh.m_posLookup[h]
	for tile != nil {
		if tile.Header != nil &&
			tile.Header.X == x &&
			tile.Header.Y == y &&
			tile.Header.Layer == layer {
			return mesh.GetTileRef(tile)
		}
		tile = tile.Next
	}
	return 0
}

func (mesh *DtNavMesh) connectExtLinks(tile *DtMeshTile, target *DtMeshTile, side int32) {
	if tile == nil {
		return
	}

	// Connect border links.
	for i := int32(0); i < tile.Header.PolyCount; i++ {
		poly := tile.Polys[i]

		// Create new links.
		//		unsigned short m = DT_EXT_LINK | (unsigned short)side;
		nv := int32(poly.VertCount)
		for j := int32(0); j < nv; j++ {
			// Skip non-portal edges.
			if (poly.Neis[j] & DT_EXT_LINK) == 0 {
				continue
			}

			dir := int32(poly.Neis[j] & 0xff)
			if side != -1 && dir != side {
				continue
			}

			// Create new links
			va := common.GetVert3(tile.Verts, poly.Verts[j])
			vb := common.GetVert3(tile.Verts, poly.Verts[(j+1)%nv])

			nei, neia, nnei := mesh.FindConnectingPolys(va, vb, target, dtOppositeTile(dir), 4)
			for k := int32(0); k < nnei; k++ {
				idx := allocLink(tile)
				if idx != DT_NULL_LINK {
					link := tile.Links[idx]
					link.Ref = nei[k]
					link.Edge = uint8(j)
					link.Side = uint8(dir)
					link.Next = poly.FirstLink
					poly.FirstLink = idx

					// Compress portal limits to a byte value.
					if dir == 0 || dir == 4 {
						tmin := (neia[k*2+0] - va[2]) / (vb[2] - va[2])
						tmax := (neia[k*2+1] - va[2]) / (vb[2] - va[2])
						if tmin > tmax {
							tmin, tmax = tmax, tmin
						}

						link.Bmin = uint8(math.Round(float64(common.Clamp(tmin, 0.0, 1.0) * 255.0)))
						link.Bmax = uint8(math.Round(float64(common.Clamp(tmax, 0.0, 1.0) * 255.0)))
					} else if dir == 2 || dir == 6 {
						tmin := (neia[k*2+0] - va[0]) / (vb[0] - va[0])
						tmax := (neia[k*2+1] - va[0]) / (vb[0] - va[0])
						if tmin > tmax {
							tmin, tmax = tmax, tmin
						}
						link.Bmin = uint8(math.Round(float64(common.Clamp(tmin, 0.0, 1.0) * 255.0)))
						link.Bmax = uint8(common.Clamp(tmax, 0.0, 1.0) * 255.0)
					}
				}
			}
		}
	}
}

func (mesh *DtNavMesh) ConnectExtOffMeshLinks(tile *DtMeshTile, target *DtMeshTile, side int32) {
	if tile == nil {
		return
	}

	// Connect off-mesh links.
	// We are interested on links which land from target tile to this tile.
	oppositeSide := dtOppositeTile(side)
	if side == -1 {
		oppositeSide = 0xff
	}

	for i := int32(0); i < target.Header.OffMeshConCount; i++ {
		targetCon := target.OffMeshCons[i]
		if int32(targetCon.Side) != oppositeSide {
			continue
		}

		targetPoly := target.Polys[targetCon.Poly]
		// Skip off-mesh connections which start location could not be connected at all.
		if targetPoly.FirstLink == DT_NULL_LINK {
			continue
		}

		var halfExtents = []float32{targetCon.Rad, target.Header.WalkableClimb, targetCon.Rad}

		// Find polygon to connect to.
		p := targetCon.Pos[3:]

		nearestPt, ref := mesh.FindNearestPolyInTile(tile, p, halfExtents)
		if ref == 0 {
			continue
		}

		// findNearestPoly may return too optimistic results, further check to make sure.
		if common.Sqr(nearestPt[0]-p[0])+common.Sqr(nearestPt[2]-p[2]) > common.Sqr(targetCon.Rad) {
			continue
		}

		// Make sure the location is on current mesh.
		copy(target.Verts[targetPoly.Verts[1]*3:targetPoly.Verts[1]*3+3], nearestPt)
		// Link off-mesh connection to target poly.
		idx := allocLink(target)
		if idx != DT_NULL_LINK {
			link := target.Links[idx]
			link.Ref = ref
			link.Edge = 1
			link.Side = uint8(oppositeSide)
			link.Bmax = 0
			link.Bmin = link.Bmax
			// Add to linked list.
			link.Next = targetPoly.FirstLink
			targetPoly.FirstLink = idx
		}

		// Link target poly to off-mesh connection.
		if targetCon.Flags&DT_OFFMESH_CON_BIDIR > 0 {
			tidx := allocLink(tile)
			if tidx != DT_NULL_LINK {
				landPolyIdx := mesh.DecodePolyIdPoly(ref)
				landPoly := tile.Polys[landPolyIdx]
				link := tile.Links[tidx]
				link.Ref = mesh.GetPolyRefBase(target) | DtPolyRef(targetCon.Poly)
				link.Edge = 0xff
				link.Side = uint8(side)
				if side == -1 {
					side = 0xff
					link.Side = uint8(side)
				}
				link.Bmax = 0
				link.Bmin = link.Bmax
				// Add to linked list.
				link.Next = landPoly.FirstLink
				landPoly.FirstLink = tidx
			}
		}
	}

}

func (mesh *DtNavMesh) FindNearestPolyInTile(tile *DtMeshTile, center []float32, halfExtents []float32) (nearestPt []float32, nearest DtPolyRef) {
	nearestPt = make([]float32, 3)
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	common.Vsub(bmin, center, halfExtents)
	common.Vadd(bmax, center, halfExtents)

	// Get nearby polygons from proximity grid.
	polys, polyCount := mesh.QueryPolygonsInTile(tile, bmin, bmax, 128)

	// Find nearest polygon amongst the nearby polygons.

	nearestDistanceSqr := float32(math.MaxFloat32)
	for i := int32(0); i < polyCount; i++ {
		ref := polys[i]
		var d float32
		closestPtPoly := make([]float32, 3)
		posOverPoly := false
		mesh.ClosestPointOnPoly(ref, center, closestPtPoly, &posOverPoly)
		diff := make([]float32, 3)
		// If a point is directly over a polygon and closer than
		// climb height, favor that instead of straight line nearest point.
		common.Vsub(diff, center, closestPtPoly)
		if posOverPoly {
			d = common.Abs(diff[1]) - tile.Header.WalkableClimb
			d = 0
			if d > 0 {
				d = d * d
			}
		} else {
			d = common.VlenSqr(diff)
		}

		if d < nearestDistanceSqr {
			copy(nearestPt, closestPtPoly)
			nearestDistanceSqr = d
			nearest = ref
		}
	}

	return nearestPt, nearest
}

func (mesh *DtNavMesh) QueryPolygonsInTile(tile *DtMeshTile, qmin []float32, qmax []float32, maxPolys int32) (polys []DtPolyRef, n int32) {
	polys = make([]DtPolyRef, 128)
	if tile.BvTree != nil {
		node := int32(0)
		end := tile.Header.BvNodeCount
		tbmin := tile.Header.Bmin
		tbmax := tile.Header.Bmax
		qfac := tile.Header.BvQuantFactor

		// Calculate quantized box
		var (
			bmin [3]uint16
			bmax [3]uint16
		)
		// dtClamp query box to world box.
		minx := common.Clamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0]
		miny := common.Clamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1]
		minz := common.Clamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2]
		maxx := common.Clamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0]
		maxy := common.Clamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1]
		maxz := common.Clamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2]
		// Quantize
		bmin[0] = uint16(qfac*minx) & 0xfffe
		bmin[1] = uint16(qfac*miny) & 0xfffe
		bmin[2] = uint16(qfac*minz) & 0xfffe
		bmax[0] = uint16(qfac*maxx+1) | 1
		bmax[1] = uint16(qfac*maxy+1) | 1
		bmax[2] = uint16(qfac*maxz+1) | 1

		// Traverse tree
		base := mesh.GetPolyRefBase(tile)
		for node < end {
			overlap := dtOverlapQuantBounds(bmin, bmax, tile.BvTree[node].Bmin, tile.BvTree[node].Bmax)
			isLeafNode := tile.BvTree[node].I >= 0

			if isLeafNode && overlap {
				if n < maxPolys {
					polys[n] = DtPolyRef(int32(base) | tile.BvTree[node].I)
				}
				n++
			}

			if overlap || isLeafNode {
				node++
			} else {
				escapeIndex := -tile.BvTree[node].I
				node += escapeIndex
			}
		}

		return polys, n
	} else {
		bmin := make([]float32, 3)
		bmax := make([]float32, 3)
		base := mesh.GetPolyRefBase(tile)
		for i := int32(0); i < tile.Header.PolyCount; i++ {
			p := tile.Polys[i]
			// Do not return off-mesh connection polygons.
			if p.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Calc polygon bounds.
			v := common.GetVert3(tile.Verts, p.Verts[0])

			copy(bmin, v)
			copy(bmax, v)
			for j := int32(1); j < int32(p.VertCount); j++ {
				v = common.GetVert3(tile.Verts, p.Verts[j])
				common.Vmin(bmin, v)
				common.Vmax(bmax, v)
			}
			if common.DtOverlapBounds(qmin, qmax, bmin, bmax) {
				if n < maxPolys {
					polys[n] = DtPolyRef(int32(base) | i)
					n++
				}

			}
		}
		return polys, n
	}
}

func (mesh *DtNavMesh) ClosestPointOnPoly(ref DtPolyRef, pos []float32, closest []float32, posOverPoly *bool) {
	tile, poly := mesh.GetTileAndPolyByRefUnsafe(ref)
	copy(closest, pos)
	if height, ok := mesh.GetPolyHeight(tile, poly, pos); ok {
		closest[1] = height
		*posOverPoly = true
		return
	}
	*posOverPoly = false
	// Off-mesh connections don't have detail polygons.
	if poly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		v0 := common.GetVert3(tile.Verts, poly.Verts[0])
		v1 := common.GetVert3(tile.Verts, poly.Verts[1])
		t, _ := DtDistancePtSegSqr2D(pos, v0, v1)
		common.Vlerp(closest, v0, v1, t)
		return
	}

	// Outside poly that is not an offmesh connection.
	closestPointOnDetailEdges(true, tile, poly, pos, closest)
	return
}
func getPolyIndexByTitle(tile *DtMeshTile, poly *DtPoly) int {
	for index, v := range tile.Polys {
		if v == poly {
			return index
		}
	}
	return -1
}
func (mesh *DtNavMesh) GetPolyHeight(tile *DtMeshTile, poly *DtPoly, pos []float32) (height float32, ok bool) {
	// Off-mesh connections do not have detail polys and getting height
	// over them does not make sense.
	if poly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		return height, false
	}

	ip := getPolyIndexByTitle(tile, poly)
	pd := tile.DetailMeshes[ip]

	verts := make([]float32, DT_VERTS_PER_POLYGON*3)
	nv := int32(poly.VertCount)
	for i := int32(0); i < nv; i++ {
		copy(verts[i*3:i*3+3], common.GetVert3(tile.Verts, poly.Verts[i]))
	}
	if !dtPointInPolygon(pos, verts, nv) {
		return height, false
	}
	// Find height at the location.
	for j := uint8(0); j < pd.TriCount; j++ {
		t := common.GetVert4(tile.DetailTris, pd.TriBase+uint32(j))
		var v [3][]float32
		for k := 0; k < 3; k++ {
			if t[k] < poly.VertCount {
				v[k] = common.GetVert3(tile.Verts, poly.Verts[t[k]])
			} else {
				v[k] = common.GetVert3[float32, uint32](tile.Verts, pd.VertBase+uint32(t[k]-poly.VertCount))
			}
			if h, ok := dtClosestHeightPointTriangle(pos, v[0], v[1], v[2]); ok {
				return h, true
			}
		}
	}
	// If all triangle checks failed above (can happen with degenerate triangles
	// or larger floating point values) the point is on an edge, so just select
	// closest. This should almost never happen so the extra iteration here is
	// ok.
	closest := make([]float32, 3)
	closestPointOnDetailEdges(false, tile, poly, pos, closest)
	height = closest[1]
	return height, true

}

// / @par
// /
// / @warning Only use this function if it is known that the provided polygon
// / reference is valid. This function is faster than #getTileAndPolyByRef, but
// / it does not validate the reference.
func (mesh *DtNavMesh) GetTileAndPolyByRefUnsafe(ref DtPolyRef) (tile *DtMeshTile, poly *DtPoly) {
	_, it, ip := mesh.DecodePolyId(ref)
	tile = mesh.m_tiles[it]
	poly = mesh.m_tiles[it].Polys[ip]
	return
}
func (mesh *DtNavMesh) GetTileAndPolyByRef(ref DtPolyRef) (tile *DtMeshTile, poly *DtPoly, status DtStatus) {
	if ref == 0 {
		return tile, poly, DT_FAILURE
	}

	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return tile, poly, DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return tile, poly, DT_FAILURE | DT_INVALID_PARAM
	}
	if int64(ip) >= int64(mesh.m_tiles[it].Header.PolyCount) {
		return tile, poly, DT_FAILURE | DT_INVALID_PARAM
	}
	tile = mesh.m_tiles[it]
	poly = mesh.m_tiles[it].Polys[ip]
	return tile, poly, DT_SUCCESS
}

func closestPointOnDetailEdges(onlyBoundary bool, tile *DtMeshTile, poly *DtPoly, pos []float32, closest []float32) {
	ip := getPolyIndexByTitle(tile, poly)
	pd := tile.DetailMeshes[ip]
	dmin := float32(math.MaxFloat32)
	tmin := float32(0)
	var pmin []float32
	var pmax []float32

	for i := uint8(0); i < pd.TriCount; i++ {
		tris := common.GetVert4(tile.DetailTris, pd.TriBase+uint32(i))
		ANY_BOUNDARY_EDGE := (DT_DETAIL_EDGE_BOUNDARY << 0) | (DT_DETAIL_EDGE_BOUNDARY << 2) | (DT_DETAIL_EDGE_BOUNDARY << 4)
		if onlyBoundary && (int(tris[3])&ANY_BOUNDARY_EDGE) == 0 {
			continue
		}
		var v [3][]float32
		for j := 0; j < 3; j++ {
			if tris[j] < poly.VertCount {
				v[j] = common.GetVert3(tile.Verts, poly.Verts[tris[j]])
			} else {
				v[j] = common.GetVert3[float32, uint32](tile.Verts, pd.VertBase+uint32(tris[j]-poly.VertCount))
			}
		}
		k := int32(0)
		j := int32(2)
		for k < 3 {
			if (DtGetDetailTriEdgeFlags(tris[3], j)&DT_DETAIL_EDGE_BOUNDARY) == 0 && (onlyBoundary || tris[j] < tris[k]) {
				// Only looking at boundary edges and this is internal, or
				// this is an inner edge that we will see again or have already seen.
				j = k
				k++
				continue
			}

			t, d := DtDistancePtSegSqr2D(pos, v[j], v[k])
			if d < dmin {
				dmin = d
				tmin = t
				pmin = v[j]
				pmax = v[k]
			}
			j = k
			k++
		}
	}

	common.Vlerp(closest, pmin, pmax, tmin)
}

// / @par
// /
// / All points are projected onto the xz-plane, so the y-values are ignored.
func dtPointInPolygon(pt, verts []float32, nverts int32) bool {
	// TODO: Replace pnpoly with triArea2D tests?
	i := int32(0)
	j := nverts - 1
	c := false
	for i < nverts {
		vi := common.GetVert3(verts, i)
		vj := common.GetVert3(verts, j)
		if ((vi[2] > pt[2]) != (vj[2] > pt[2])) && (pt[0] < (vj[0]-vi[0])*(pt[2]-vi[2])/(vj[2]-vi[2])+vi[0]) {
			c = !c
		}

		j = i
		i++
	}
	return c
}

type dtTileState struct {
	magic   int32     // Magic number, used to identify the data.
	version int32     // Data version number.
	ref     DtTileRef // Tile ref at the time of storing the data.
}

type dtPolyState struct {
	flags uint16 // Flags (see dtPolyFlags).
	area  uint8  // Area ID of the polygon.
}

// / @par
// /
// / Tile state includes non-structural data such as polygon flags, area ids, etc.
// / @note This function does not impact the tile's #DtTileRef and #DtPolyRef's.
// / @see #storeTileState
func (mesh *DtNavMesh) RestoreTileState(tile *DtMeshTile, tileState *dtTileState, polyStates []*dtPolyState, maxDataSize int32) DtStatus {
	// Make sure there is enough space to store the state.
	// Check that the restore is possible.
	if tileState.magic != DT_NAVMESH_STATE_MAGIC {
		return DT_FAILURE | DT_WRONG_MAGIC
	}

	if tileState.version != DT_NAVMESH_STATE_VERSION {
		return DT_FAILURE | DT_WRONG_VERSION
	}

	if tileState.ref != mesh.GetTileRef(tile) {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	// Restore per poly state.
	for i := int32(0); i < tile.Header.PolyCount; i++ {
		p := tile.Polys[i]
		s := polyStates[i]
		p.Flags = s.flags
		p.SetArea(s.area)
	}

	return DT_SUCCESS
}

// / @par
// /
// / The add operation will fail if the data is in the wrong format, the allocated tile
// / space is full, or there is a tile already at the specified reference.
// /
// / The lastRef parameter is used to restore a tile with the same tile
// / reference it had previously used.  In this case the #DtPolyRef's for the
// / tile will be restored to the same values they were before the tile was
// / removed.
// /
// / The nav mesh assumes exclusive access to the data passed and will make
// / changes to the dynamic portion of the data. For that reason the data
// / should not be reused in other nav meshes until the tile has been successfully
// / removed from this nav mesh.
// /
// / @see dtCreateNavMeshData, #removeTile
func (mesh *DtNavMesh) AddTile(data *recast.NavMeshData, flags int32, lastRef DtTileRef) (result DtTileRef, status DtStatus) {
	// Make sure the data is in right format.
	header := data.Header
	if header.Magic != DT_NAVMESH_MAGIC {
		return result, DT_FAILURE | DT_WRONG_MAGIC
	}

	if header.Version != DT_NAVMESH_VERSION {
		return result, DT_FAILURE | DT_WRONG_VERSION
	}

	if DT_POLYREF64 == 1 {
		// Do not allow adding more polygons than specified in the NavMesh's maxPolys constraint.
		// Otherwise, the poly ID cannot be represented with the given number of bits.
		if mesh.m_polyBits < common.Ilog2(common.NextPow2(uint32(header.PolyCount))) {
			return result, DT_FAILURE | DT_INVALID_PARAM
		}

	}
	// Make sure the location is free.
	if mesh.GetTileAt(header.X, header.Y, header.Layer) != nil {
		return result, DT_FAILURE | DT_ALREADY_OCCUPIED
	}

	var tile *DtMeshTile
	if lastRef == 0 {
		if mesh.m_nextFree != nil {
			tile = mesh.m_nextFree
			mesh.m_nextFree = tile.Next
			tile.Next = nil
		}
	} else {
		// Try to relocate the tile to specific index with same salt.
		tileIndex := mesh.DecodePolyIdTile(DtPolyRef(lastRef))
		if int64(tileIndex) >= int64(mesh.m_maxTiles) {
			return result, DT_FAILURE | DT_OUT_OF_MEMORY
		}

		// Try to find the specific tile id from the free list.
		target := mesh.m_tiles[tileIndex]
		var prev *DtMeshTile
		tile = mesh.m_nextFree
		for tile != nil && tile != target {
			prev = tile
			tile = tile.Next
		}
		// Could not find the correct location.
		if tile != target {
			return result, DT_FAILURE | DT_OUT_OF_MEMORY
		}

		// Remove from freelist
		if prev == nil {
			mesh.m_nextFree = tile.Next
		} else {
			prev.Next = tile.Next
		}

		// Restore salt.
		tile.salt = mesh.DecodePolyIdSalt(DtPolyRef(lastRef))
	}

	// Make sure we could allocate a tile.
	if tile == nil {
		return result, DT_FAILURE | DT_OUT_OF_MEMORY
	}

	// Insert tile into the position lut.
	h := common.ComputeTileHash(header.X, header.Y, mesh.m_tileLutMask)
	tile.Next = mesh.m_posLookup[h]
	mesh.m_posLookup[h] = tile

	// Patch header pointers.
	tile.Verts = data.NavVerts
	tile.Polys = data.NavPolys
	tile.Links = data.Links
	tile.DetailMeshes = data.NavDMeshes
	tile.DetailVerts = data.NavDVerts
	tile.DetailTris = data.NavDTris
	tile.BvTree = data.NavBvtree
	tile.OffMeshCons = data.OffMeshCons

	// If there are no items in the bvtree, reset the tree pointer.
	// Build links freelist
	tile.linksFreeList = 0
	tile.Links[header.MaxLinkCount-1].Next = DT_NULL_LINK
	for i := int32(0); i < header.MaxLinkCount-1; i++ {
		tile.Links[i].Next = uint32(i) + 1
	}

	// Init tile.
	tile.Header = header
	tile.Data = data
	tile.Flags = flags

	mesh.connectIntLinks(tile)

	// Base off-mesh connections to their starting polygons and connect connections inside the tile.
	mesh.baseOffMeshLinks(tile)
	mesh.ConnectExtOffMeshLinks(tile, tile, -1)

	// Create connections with neighbour tiles.
	MAX_NEIS := int32(32)

	// Connect with layers in current tile.
	neis, nneis := mesh.GetTilesAt(header.X, header.Y, MAX_NEIS)
	for j := int32(0); j < nneis; j++ {
		if neis[j] == tile {
			continue
		}

		mesh.connectExtLinks(tile, neis[j], -1)
		mesh.connectExtLinks(neis[j], tile, -1)
		mesh.ConnectExtOffMeshLinks(tile, neis[j], -1)
		mesh.ConnectExtOffMeshLinks(neis[j], tile, -1)
	}

	// Connect with neighbour tiles.
	for i := int32(0); i < 8; i++ {
		neis, nneis := mesh.getNeighbourTilesAt(header.X, header.Y, i, MAX_NEIS)
		for j := int32(0); j < nneis; j++ {
			mesh.connectExtLinks(tile, neis[j], i)
			mesh.connectExtLinks(neis[j], tile, dtOppositeTile(i))
			mesh.ConnectExtOffMeshLinks(tile, neis[j], i)
			mesh.ConnectExtOffMeshLinks(neis[j], tile, dtOppositeTile(i))
		}
	}

	result = mesh.GetTileRef(tile)

	return result, DT_SUCCESS
}

func (mesh *DtNavMesh) baseOffMeshLinks(tile *DtMeshTile) {
	if tile == nil {
		return
	}

	base := mesh.GetPolyRefBase(tile)

	// Base off-mesh connection start points.
	for i := int32(0); i < tile.Header.OffMeshConCount; i++ {
		con := tile.OffMeshCons[i]
		poly := tile.Polys[con.Poly]

		halfExtents := [3]float32{con.Rad, tile.Header.WalkableClimb, con.Rad}

		// Find polygon to connect to.
		p := con.Pos[0:3] // First vertex
		nearestPt, ref := mesh.FindNearestPolyInTile(tile, p, halfExtents[:])
		if ref == 0 {
			continue
		}
		// findNearestPoly may return too optimistic results, further check to make sure.
		if common.Sqr(nearestPt[0]-p[0])+common.Sqr(nearestPt[2]-p[2]) > common.Sqr(con.Rad) {
			continue
		}

		// Make sure the location is on current mesh.
		copy(tile.Verts[poly.Verts[0]*3:poly.Verts[0]*3+3], nearestPt)

		// Link off-mesh connection to target poly.
		idx := allocLink(tile)
		if idx != DT_NULL_LINK {
			link := tile.Links[idx]
			link.Ref = ref
			link.Edge = 0
			link.Side = 0xff
			link.Bmax = 0
			link.Bmin = link.Bmax
			// Add to linked list.
			link.Next = poly.FirstLink
			poly.FirstLink = idx
		}

		// Start end-point is always connect back to off-mesh connection.
		tidx := allocLink(tile)
		if tidx != DT_NULL_LINK {
			landPolyIdx := mesh.DecodePolyIdPoly(ref)
			landPoly := tile.Polys[landPolyIdx]
			link := tile.Links[tidx]
			link.Ref = base | DtPolyRef(con.Poly)
			link.Edge = 0xff
			link.Side = 0xff
			link.Bmax = 0
			link.Bmin = link.Bmax
			// Add to linked list.
			link.Next = landPoly.FirstLink
			landPoly.FirstLink = tidx
		}
	}
}

func (mesh *DtNavMesh) connectIntLinks(tile *DtMeshTile) {
	if tile == nil {
		return
	}

	base := mesh.GetPolyRefBase(tile)

	for i := int32(0); i < tile.Header.PolyCount; i++ {
		poly := tile.Polys[i]
		poly.FirstLink = DT_NULL_LINK

		if poly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
			continue
		}

		// Build edge links backwards so that the links will be
		// in the linked list from lowest index to highest.
		for j := poly.VertCount - 1; j >= 0; j-- {
			// Skip hard and non-internal edges.
			if poly.Neis[j] == 0 || (poly.Neis[j]&DT_EXT_LINK) > 0 {
				continue
			}

			idx := allocLink(tile)
			if idx != DT_NULL_LINK {
				link := tile.Links[idx]
				link.Ref = base | (DtPolyRef)(poly.Neis[j]-1)
				link.Edge = j
				link.Side = 0xff
				link.Bmax = 0
				link.Bmin = link.Bmax
				// Add to linked list.
				link.Next = poly.FirstLink
				poly.FirstLink = idx
			}
		}
	}
}

func (mesh *DtNavMesh) GetTileAt(x, y, layer int32) *DtMeshTile {
	// Find tile based on hash.
	h := common.ComputeTileHash(x, y, mesh.m_tileLutMask)
	tile := mesh.m_posLookup[h]
	for tile != nil {
		if tile.Header != nil && tile.Header.X == x && tile.Header.Y == y && tile.Header.Layer == layer {
			return tile
		}
		tile = tile.Next
	}
	return nil
}

// / @par
// /
// / This function returns the data for the tile so that, if desired,
// / it can be added back to the navigation mesh at a later point.
// /
// / @see #addTile
func (mesh *DtNavMesh) RemoveTile(ref DtTileRef) (data *recast.NavMeshData, status DtStatus) {
	if ref == 0 {
		return data, DT_FAILURE | DT_INVALID_PARAM
	}

	tileIndex := mesh.DecodePolyIdTile((DtPolyRef(ref)))
	tileSalt := mesh.DecodePolyIdSalt((DtPolyRef(ref)))
	if int64(tileIndex) >= int64(mesh.m_maxTiles) {
		return data, DT_FAILURE | DT_INVALID_PARAM
	}

	tile := mesh.m_tiles[tileIndex]
	if tile.salt != tileSalt {
		return data, DT_FAILURE | DT_INVALID_PARAM
	}

	// Remove tile from hash lookup.
	h := common.ComputeTileHash(tile.Header.X, tile.Header.Y, mesh.m_tileLutMask)
	var prev *DtMeshTile
	cur := mesh.m_posLookup[h]
	for cur != nil {
		if cur == tile {
			if prev != nil {
				prev.Next = cur.Next
			} else {
				mesh.m_posLookup[h] = cur.Next
			}

			break
		}
		prev = cur
		cur = cur.Next
	}
	MAX_NEIS := int32(32)
	// Remove connections to neighbour tiles.
	// Disconnect from other layers in current tile.
	neis, nneis := mesh.GetTilesAt(tile.Header.X, tile.Header.Y, MAX_NEIS)
	for j := int32(0); j < nneis; j++ {
		if neis[j] == tile {
			continue
		}
		mesh.unconnectLinks(neis[j], tile)
	}

	// Disconnect from neighbour tiles.
	for i := int32(0); i < 8; i++ {
		neis, nneis = mesh.getNeighbourTilesAt(tile.Header.X, tile.Header.Y, i, MAX_NEIS)
		for j := int32(0); j < nneis; j++ {
			mesh.unconnectLinks(neis[j], tile)
		}

	}

	// Reset tile.
	if tile.Flags&DT_TILE_FREE_DATA > 0 {
		// Owns data
		tile.Data = nil
	} else if data != nil {
		data = tile.Data
	}

	tile.Header = nil
	tile.Flags = 0
	tile.linksFreeList = 0
	tile.Polys = tile.Polys[:]
	tile.Verts = tile.Verts[:]
	tile.Links = tile.Links[:]
	tile.DetailMeshes = tile.DetailMeshes[:]
	tile.DetailVerts = tile.DetailVerts[:]
	tile.DetailTris = tile.DetailTris[:]
	tile.BvTree = tile.BvTree[:]
	tile.OffMeshCons = tile.OffMeshCons[:]

	// Update salt, salt should never be zero.
	if DT_POLYREF64 == 1 {
		tile.salt = (tile.salt + 1) & ((1 << DT_SALT_BITS) - 1)
	} else {
		tile.salt = (tile.salt + 1) & ((1 << mesh.m_saltBits) - 1)
	}

	if tile.salt == 0 {
		tile.salt++
	}

	// Add to free list.
	tile.Next = mesh.m_nextFree
	mesh.m_nextFree = tile

	return data, DT_SUCCESS
}

func (mesh *DtNavMesh) getNeighbourTilesAt(x, y, side, maxTiles int32) ([]*DtMeshTile, int32) {
	nx := x
	ny := y
	switch side {
	case 0:
		nx++
		break
	case 1:
		nx++
		ny++
		break
	case 2:
		ny++
		break
	case 3:
		nx--
		ny++
		break
	case 4:
		nx--
		break
	case 5:
		nx--
		ny--
		break
	case 6:
		ny--
		break
	case 7:
		nx++
		ny--
		break
	}

	return mesh.GetTilesAt(nx, ny, maxTiles)
}

func (mesh *DtNavMesh) GetTilesAt(x, y int32, maxTiles int32) ([]*DtMeshTile, int32) {
	MAX_NEIS := 32
	tiles := make([]*DtMeshTile, MAX_NEIS)
	n := int32(0)
	// Find tile based on hash.
	h := common.ComputeTileHash(x, y, mesh.m_tileLutMask)
	tile := mesh.m_posLookup[h]
	for tile != nil {
		if tile.Header != nil && tile.Header.X == x && tile.Header.Y == y {
			if n < maxTiles {
				tiles[n] = tile
				n++
			}

		}
		tile = tile.Next
	}

	return tiles, n
}
func (mesh *DtNavMesh) IsValidPolyRef(ref DtPolyRef) bool {
	if ref == 0 {
		return false
	}
	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return false
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return false
	}
	if int64(ip) >= int64(mesh.m_tiles[it].Header.PolyCount) {
		return false
	}
	return true
}

// / @par
// /
// / Off-mesh connections are stored in the navigation mesh as special 2-vertex
// / polygons with a single edge. At least one of the vertices is expected to be
// / inside a normal polygon. So an off-mesh connection is "entered" from a
// / normal polygon at one of its endpoints. This is the polygon identified by
// / the prevRef parameter.
func (mesh *DtNavMesh) GetOffMeshConnectionPolyEndPoints(prevRef DtPolyRef, polyRef DtPolyRef, startPos []float32, endPos []float32) DtStatus {

	if polyRef == 0 {
		return DT_FAILURE
	}

	// Get current polygon
	salt, it, ip := mesh.DecodePolyId(polyRef)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if int64(ip) >= int64(tile.Header.PolyCount) {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.Polys[ip]

	// Make sure that the current poly is indeed off-mesh link.
	if poly.GetType() != DT_POLYTYPE_OFFMESH_CONNECTION {
		return DT_FAILURE
	}

	// Figure out which way to hand out the vertices.
	idx0 := 0
	idx1 := 1

	// Find link that points to first vertex.
	for i := poly.FirstLink; i != DT_NULL_LINK; i = tile.Links[i].Next {
		if tile.Links[i].Edge == 0 {
			if tile.Links[i].Ref != prevRef {
				idx0 = 1
				idx1 = 0
			}
			break
		}
	}

	copy(startPos, common.GetVert3(tile.Verts, poly.Verts[idx0]))
	copy(endPos, common.GetVert3(tile.Verts, poly.Verts[idx1]))

	return DT_SUCCESS
}

func (mesh *DtNavMesh) GetOffMeshConnectionByRef(ref DtPolyRef) *DtOffMeshConnection {

	if ref == 0 {
		return nil
	}

	// Get current polygon
	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return nil
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return nil
	}
	tile := mesh.m_tiles[it]
	if int64(ip) >= int64(tile.Header.PolyCount) {
		return nil
	}
	poly := tile.Polys[ip]

	// Make sure that the current poly is indeed off-mesh link.
	if poly.GetType() != DT_POLYTYPE_OFFMESH_CONNECTION {
		return nil
	}

	idx := int(ip) - int(tile.Header.OffMeshBase)
	if int64(idx) < int64(tile.Header.OffMeshConCount) {
		panic("not expect")
	}
	return tile.OffMeshCons[idx]
}

func (mesh *DtNavMesh) SetPolyFlags(ref DtPolyRef, flags uint16) DtStatus {
	if ref == 0 {
		return DT_FAILURE
	}
	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if int64(ip) >= int64(tile.Header.PolyCount) {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.Polys[ip]
	// Change flags.
	poly.Flags = flags
	return DT_SUCCESS
}

func (mesh *DtNavMesh) GetPolyFlags(ref DtPolyRef) (resultFlags uint16, status DtStatus) {
	if ref == 0 {
		return resultFlags, DT_FAILURE
	}
	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return resultFlags, DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return resultFlags, DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if int64(ip) >= int64(tile.Header.PolyCount) {
		return resultFlags, DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.Polys[ip]

	resultFlags = poly.Flags

	return resultFlags, DT_SUCCESS
}

func (mesh *DtNavMesh) SetPolyArea(ref DtPolyRef, area uint8) DtStatus {
	if ref == 0 {
		return DT_FAILURE
	}
	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if int64(ip) >= int64(tile.Header.PolyCount) {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.Polys[ip]

	poly.SetArea(area)

	return DT_SUCCESS
}
func (mesh *DtNavMesh) GetPolyArea(ref DtPolyRef) (resultArea uint8, status DtStatus) {
	if ref == 0 {
		return resultArea, DT_FAILURE
	}
	salt, it, ip := mesh.DecodePolyId(ref)
	if int64(it) >= int64(mesh.m_maxTiles) {
		return resultArea, DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].Header == nil {
		return resultArea, DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if int64(ip) >= int64(tile.Header.PolyCount) {
		return resultArea, DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.Polys[ip]

	resultArea = poly.GetArea()

	return resultArea, DT_SUCCESS
}
