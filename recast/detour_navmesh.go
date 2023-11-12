package recast

import "math"

func (mesh *dtNavMesh) init1(params *dtNavMeshParams) dtStatus {
	mesh.m_orig = params.orig
	mesh.m_tileWidth = params.tileWidth
	mesh.m_tileHeight = params.tileHeight
	mesh.m_params = params

	// Init tiles
	mesh.m_maxTiles = params.maxTiles
	mesh.m_tileLutSize = dtNextPow2(params.maxTiles / 4)
	if mesh.m_tileLutSize == 0 {
		mesh.m_tileLutSize = 1
	}
	mesh.m_tileLutMask = mesh.m_tileLutSize - 1
	var m_nextFree *dtMeshTile
	for i := mesh.m_maxTiles - 1; i >= 0; i-- {
		mesh.m_tiles[i].salt = 1
		mesh.m_tiles[i].next = m_nextFree
		m_nextFree = mesh.m_tiles[i]
	}

	// Init ID generator values.
	if DT_POLYREF64 == 1 {
		mesh.m_tileBits = dtIlog2(dtNextPow2(params.maxTiles))
		mesh.m_polyBits = dtIlog2(dtNextPow2(params.maxPolys))
		// Only allow 31 salt bits, since the salt mask is calculated using 32bit uint and it will overflow.
		mesh.m_saltBits = dtMin(31, 32-mesh.m_tileBits-mesh.m_polyBits)

		if mesh.m_saltBits < 10 {
			return DT_FAILURE | DT_INVALID_PARAM
		}

	}
	return DT_SUCCESS
}

func (mesh *dtNavMesh) init3(header *dtMeshHeader, titleData *dtMeshTile, dataSize int, flags int) (result dtTileRef, status dtStatus) {
	// Make sure the data is in right format.
	if header.magic != DT_NAVMESH_MAGIC {
		return result, DT_FAILURE | DT_WRONG_MAGIC
	}

	if header.version != DT_NAVMESH_VERSION {
		return result, DT_FAILURE | DT_WRONG_VERSION
	}

	var params dtNavMeshParams
	params.orig = header.bmin
	params.tileWidth = header.bmax[0] - header.bmin[0]
	params.tileHeight = header.bmax[2] - header.bmin[2]
	params.maxTiles = 1
	params.maxPolys = header.polyCount

	status = mesh.init1(&params)
	if status.dtStatusFailed() {
		return result, status
	}

	return mesh.addTile(header, titleData, dataSize, flags, 0)
}

// / @par
// /
// / @note The parameters are created automatically when the single tile
// / initialization is performed.
func (mesh *dtNavMesh) getParams() *dtNavMeshParams {
	return mesh.m_params
}
func overlapSlabs(amin, amax, bmin, bmax []float64, px, py float64) bool {
	// Check for horizontal overlap.
	// The segment is shrunken a little so that slabs which touch
	// at end points are not connected.
	minx := dtMax(amin[0]+px, bmin[0]+px)
	maxx := dtMin(amax[0]-px, bmax[0]-px)
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
	thr := dtSqr(py * 2)
	if dmin*dmin <= thr || dmax*dmax <= thr {
		return true
	}

	return false
}

func getSlabCoord(va []float64, side int) float64 {
	if side == 0 || side == 4 {
		return va[0]
	} else if side == 2 || side == 6 {
		return va[2]
	}
	return 0
}

func calcSlabEndPoints(va, vb []float64, side int) ([]float64, []float64) {
	bmin, bmax := make([]float64, 3), make([]float64, 3)
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

func computeTileHash(x, y, mask int) int {
	h1 := 0x8da6b343 // Large multiplicative constants;
	h2 := 0xd8163841 // here arbitrarily chosen primes
	n := h1*x + h2*y
	return n & mask
}

func allocLink(tile *dtMeshTile) int {
	if tile.linksFreeList == DT_NULL_LINK {
		return DT_NULL_LINK
	}

	link := tile.linksFreeList
	tile.linksFreeList = tile.links[link].next
	return link
}

func freeLink(tile *dtMeshTile, link int) {
	tile.links[link].next = tile.linksFreeList
	tile.linksFreeList = link
}

// ////////////////////////////////////////////////////////////////////////////////////////
func (mesh *dtNavMesh) findConnectingPolys(va, vb []float64, tile *dtMeshTile, side int, maxcon int) (con []dtPolyRef, conarea []float64, n int) {
	con = make([]dtPolyRef, 4)
	conarea = make([]float64, 4*2)
	if tile == nil {
		return
	}

	amin, amax := calcSlabEndPoints(va, vb, side)
	apos := getSlabCoord(va, side)

	// Remove links pointing to 'side' and compact the links array.

	m := DT_EXT_LINK | side
	base := mesh.getPolyRefBase(tile)

	for i := 0; i < tile.header.polyCount; i++ {
		poly := tile.polys[i]
		nv := poly.vertCount
		for j := 0; j < nv; j++ {
			// Skip edges which do not point to the right side.
			if poly.neis[j] != m {
				continue
			}

			vc := rcGetVert(tile.verts, poly.verts[j])
			vd := rcGetVert(tile.verts, (j+1)%nv)
			bpos := getSlabCoord(vc, side)

			// Segments are not close enough.
			if dtAbs(apos-bpos) > 0.01 {
				continue
			}

			// Check if the segments touch.
			bmin, bmax := calcSlabEndPoints(vc, vd, side)

			if !overlapSlabs(amin, amax, bmin, bmax, 0.01, tile.header.walkableClimb) {
				continue
			}

			// Add return value.
			if n < maxcon {
				conarea[n*2+0] = dtMax(amin[0], bmin[0])
				conarea[n*2+1] = dtMin(amax[0], bmax[0])
				con[n] = dtPolyRef(int(base) | i)
				n++
			}
			break
		}
	}
	return
}

func (mesh *dtNavMesh) getPolyRefBase(tile *dtMeshTile) dtPolyRef {
	if tile == nil {
		return 0
	}
	it := mesh.getTileIndex(tile)
	return mesh.encodePolyId(tile.salt, it, 0)
}

func (mesh *dtNavMesh) unconnectLinks(tile *dtMeshTile, target *dtMeshTile) {
	if tile == nil || target == nil {
		return
	}

	targetNum := mesh.decodePolyIdTile(dtPolyRef(mesh.getTileRef(target)))

	for i := 0; i < tile.header.polyCount; i++ {
		poly := tile.polys[i]
		j := poly.firstLink
		pj := DT_NULL_LINK
		for j != DT_NULL_LINK {
			if mesh.decodePolyIdTile(tile.links[j].ref) == targetNum {
				// Remove link.
				nj := tile.links[j].next
				if pj == DT_NULL_LINK {
					poly.firstLink = nj
				} else {
					tile.links[pj].next = nj
				}

				freeLink(tile, j)
				j = nj
			} else {
				// Advance
				pj = j
				j = tile.links[j].next
			}
		}
	}
}
func (mesh *dtNavMesh) getTileIndex(tile *dtMeshTile) int {
	for index, v := range mesh.m_tiles {
		if v == tile {
			return index
		}
	}
	return -1
}
func (mesh *dtNavMesh) getTileRef(tile *dtMeshTile) dtTileRef {
	if tile == nil {
		return dtTileRef(0)
	}
	it := mesh.getTileIndex(tile)
	return dtTileRef(mesh.encodePolyId(tile.salt, it, 0))
}

func (mesh *dtNavMesh) getTileRefAt(x, y, layer int) dtTileRef {
	// Find tile based on hash.
	h := computeTileHash(x, y, mesh.m_tileLutMask)
	tile := mesh.m_posLookup[h]
	for tile != nil {
		if tile.header != nil &&
			tile.header.x == x &&
			tile.header.y == y &&
			tile.header.layer == layer {
			return mesh.getTileRef(tile)
		}
		tile = tile.next
	}
	return 0
}

func (mesh *dtNavMesh) connectExtLinks(tile *dtMeshTile, target *dtMeshTile, side int) {
	if tile == nil {
		return
	}

	// Connect border links.
	for i := 0; i < tile.header.polyCount; i++ {
		poly := tile.polys[i]

		// Create new links.
		//		unsigned short m = DT_EXT_LINK | (unsigned short)side;
		nv := poly.vertCount
		for j := 0; j < nv; j++ {
			// Skip non-portal edges.
			if (poly.neis[j] & DT_EXT_LINK) == 0 {
				continue
			}

			dir := poly.neis[j] & 0xff
			if side != -1 && dir != side {
				continue
			}

			// Create new links
			va := rcGetVert(tile.verts, poly.verts[j])
			vb := rcGetVert(tile.verts, poly.verts[(j+1)%nv])

			nei, neia, nnei := mesh.findConnectingPolys(va, vb, target, dtOppositeTile(dir), 4)
			for k := 0; k < nnei; k++ {
				idx := allocLink(tile)
				if idx != DT_NULL_LINK {
					link := tile.links[idx]
					link.ref = nei[k]
					link.edge = j
					link.side = dir
					link.next = poly.firstLink
					poly.firstLink = idx

					// Compress portal limits to a byte value.
					if dir == 0 || dir == 4 {
						tmin := (neia[k*2+0] - va[2]) / (vb[2] - va[2])
						tmax := (neia[k*2+1] - va[2]) / (vb[2] - va[2])
						if tmin > tmax {
							tmin, tmax = tmax, tmin
						}

						link.bmin = int(math.Round(dtClamp(tmin, 0.0, 1.0) * 255.0))
						link.bmax = int(math.Round(dtClamp(tmax, 0.0, 1.0) * 255.0))
					} else if dir == 2 || dir == 6 {
						tmin := (neia[k*2+0] - va[0]) / (vb[0] - va[0])
						tmax := (neia[k*2+1] - va[0]) / (vb[0] - va[0])
						if tmin > tmax {
							tmin, tmax = tmax, tmin
						}
						link.bmin = int(math.Round(dtClamp(tmin, 0.0, 1.0) * 255.0))
						link.bmax = int(dtClamp(tmax, 0.0, 1.0) * 255.0)
					}
				}
			}
		}
	}
}

func (mesh *dtNavMesh) connectExtOffMeshLinks(tile *dtMeshTile, target *dtMeshTile, side int) {
	if tile == nil {
		return
	}

	// Connect off-mesh links.
	// We are interested on links which land from target tile to this tile.
	oppositeSide := dtOppositeTile(side)
	if side == -1 {
		oppositeSide = 0xff
	}

	for i := 0; i < target.header.offMeshConCount; i++ {
		targetCon := target.offMeshCons[i]
		if targetCon.side != oppositeSide {
			continue
		}

		targetPoly := target.polys[targetCon.poly]
		// Skip off-mesh connections which start location could not be connected at all.
		if targetPoly.firstLink == DT_NULL_LINK {
			continue
		}

		var halfExtents = []float64{targetCon.rad, target.header.walkableClimb, targetCon.rad}

		// Find polygon to connect to.
		p := targetCon.pos[3:]

		nearestPt, ref := mesh.findNearestPolyInTile(tile, p, halfExtents)
		if ref == 0 {
			continue
		}

		// findNearestPoly may return too optimistic results, further check to make sure.
		if dtSqr(nearestPt[0]-p[0])+dtSqr(nearestPt[2]-p[2]) > dtSqr(targetCon.rad) {
			continue
		}

		// Make sure the location is on current mesh.
		copy(target.verts[targetPoly.verts[1]*3:targetPoly.verts[1]*3+3], nearestPt)
		// Link off-mesh connection to target poly.
		idx := allocLink(target)
		if idx != DT_NULL_LINK {
			link := target.links[idx]
			link.ref = ref
			link.edge = 1
			link.side = oppositeSide
			link.bmax = 0
			link.bmin = link.bmax
			// Add to linked list.
			link.next = targetPoly.firstLink
			targetPoly.firstLink = idx
		}

		// Link target poly to off-mesh connection.
		if targetCon.flags&DT_OFFMESH_CON_BIDIR > 0 {
			tidx := allocLink(tile)
			if tidx != DT_NULL_LINK {
				landPolyIdx := mesh.decodePolyIdPoly(ref)
				landPoly := tile.polys[landPolyIdx]
				link := tile.links[tidx]
				link.ref = mesh.getPolyRefBase(target) | dtPolyRef(targetCon.poly)
				link.edge = 0xff
				link.side = side
				if side == -1 {
					side = 0xff
					link.side = side
				}
				link.bmax = 0
				link.bmin = link.bmax
				// Add to linked list.
				link.next = landPoly.firstLink
				landPoly.firstLink = tidx
			}
		}
	}

}

func (mesh *dtNavMesh) findNearestPolyInTile(tile *dtMeshTile, center []float64, halfExtents []float64) (nearestPt []float64, nearest dtPolyRef) {
	nearestPt = make([]float64, 3)
	bmin := dtVsub(center, halfExtents)
	bmax := dtVadd(center, halfExtents)

	// Get nearby polygons from proximity grid.
	polys, polyCount := mesh.queryPolygonsInTile(tile, bmin, bmax, 128)

	// Find nearest polygon amongst the nearby polygons.

	nearestDistanceSqr := math.MaxFloat64
	for i := 0; i < polyCount; i++ {
		ref := polys[i]
		var d float64
		closestPtPoly, posOverPoly := mesh.closestPointOnPoly(ref, center)

		// If a point is directly over a polygon and closer than
		// climb height, favor that instead of straight line nearest point.
		diff := dtVsub(center, closestPtPoly)
		if posOverPoly {
			d = dtAbs(diff[1]) - tile.header.walkableClimb
			d = 0
			if d > 0 {
				d = d * d
			}
		} else {
			d = dtVlenSqr(diff)
		}

		if d < nearestDistanceSqr {
			nearestPt = dtVcopy(closestPtPoly)
			nearestDistanceSqr = d
			nearest = ref
		}
	}

	return nearestPt, nearest
}

func (mesh *dtNavMesh) queryPolygonsInTile(tile *dtMeshTile, qmin []float64, qmax []float64, maxPolys int) (polys []dtPolyRef, n int) {
	polys = make([]dtPolyRef, 128)
	if tile.bvTree != nil {
		node := 0
		end := tile.header.bvNodeCount
		tbmin := tile.header.bmin
		tbmax := tile.header.bmax
		qfac := tile.header.bvQuantFactor

		// Calculate quantized box
		var (
			bmin [3]int
			bmax [3]int
		)
		// dtClamp query box to world box.
		minx := dtClamp(qmin[0], tbmin[0], tbmax[0]) - tbmin[0]
		miny := dtClamp(qmin[1], tbmin[1], tbmax[1]) - tbmin[1]
		minz := dtClamp(qmin[2], tbmin[2], tbmax[2]) - tbmin[2]
		maxx := dtClamp(qmax[0], tbmin[0], tbmax[0]) - tbmin[0]
		maxy := dtClamp(qmax[1], tbmin[1], tbmax[1]) - tbmin[1]
		maxz := dtClamp(qmax[2], tbmin[2], tbmax[2]) - tbmin[2]
		// Quantize
		bmin[0] = int(qfac*minx) & 0xfffe
		bmin[1] = int(qfac*miny) & 0xfffe
		bmin[2] = int(qfac*minz) & 0xfffe
		bmax[0] = int(qfac*maxx+1) | 1
		bmax[1] = int(qfac*maxy+1) | 1
		bmax[2] = int(qfac*maxz+1) | 1

		// Traverse tree
		base := mesh.getPolyRefBase(tile)
		for node < end {
			overlap := dtOverlapQuantBounds(bmin, bmax, tile.bvTree[node].bmin, tile.bvTree[node].bmax)
			isLeafNode := tile.bvTree[node].i >= 0

			if isLeafNode && overlap {
				if n < maxPolys {
					polys[n] = dtPolyRef(int(base) | tile.bvTree[node].i)
				}
				n++
			}

			if overlap || isLeafNode {
				node++
			} else {
				escapeIndex := -tile.bvTree[node].i
				node += escapeIndex
			}
		}

		return polys, n
	} else {
		base := mesh.getPolyRefBase(tile)
		for i := 0; i < tile.header.polyCount; i++ {
			p := tile.polys[i]
			// Do not return off-mesh connection polygons.
			if p.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Calc polygon bounds.
			v := rcGetVert(tile.verts, p.verts[0])
			bmin := dtVcopy(v)
			bmax := dtVcopy(v)
			for j := 1; j < p.vertCount; j++ {
				v = rcGetVert(tile.verts, p.verts[j])
				bmin = dtVmin(bmin, v)
				bmax = dtVmax(bmax, v)
			}
			if dtOverlapBounds(qmin, qmax, bmin, bmax) {
				if n < maxPolys {
					polys[n] = dtPolyRef(int(base) | i)
					n++
				}

			}
		}
		return polys, n
	}
}

func (mesh *dtNavMesh) closestPointOnPoly(ref dtPolyRef, pos []float64) (closest []float64, posOverPoly bool) {
	closest = make([]float64, 3)
	tile, poly := mesh.getTileAndPolyByRefUnsafe(ref)
	closest = dtVcopy(pos)
	if height, ok := mesh.getPolyHeight(tile, poly, pos); ok {
		closest[1] = height
		if posOverPoly {
			posOverPoly = true
		}
		return
	}

	if posOverPoly {
		posOverPoly = false
	}

	// Off-mesh connections don't have detail polygons.
	if poly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		v0 := rcGetVert(tile.verts, poly.verts[0])
		v1 := rcGetVert(tile.verts, poly.verts[1])
		t, _ := dtDistancePtSegSqr2D(pos, v0, v1)
		closest = dtVlerp(v0, v1, t)
		return
	}

	// Outside poly that is not an offmesh connection.
	closest = closestPointOnDetailEdges(true, tile, poly, pos)
	return
}
func getPolyIndexByTitle(tile *dtMeshTile, poly *dtPoly) int {
	for index, v := range tile.polys {
		if v == poly {
			return index
		}
	}
	return -1
}
func (mesh *dtNavMesh) getPolyHeight(tile *dtMeshTile, poly *dtPoly, pos []float64) (height float64, ok bool) {
	// Off-mesh connections do not have detail polys and getting height
	// over them does not make sense.
	if poly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		return height, false
	}

	ip := getPolyIndexByTitle(tile, poly)
	pd := tile.detailMeshes[ip]

	verts := make([]float64, DT_VERTS_PER_POLYGON*3)
	nv := poly.vertCount
	for i := 0; i < nv; i++ {
		copy(verts[i*3:i*3+3], rcGetVert(tile.verts, poly.verts[i]))
	}
	if !dtPointInPolygon(pos, verts, nv) {
		return height, false
	}
	// Find height at the location.
	for j := 0; j < pd.triCount; j++ {
		t := rcGetTris(tile.detailTris, pd.triBase+j)
		var v [3][]float64
		for k := 0; k < 3; k++ {
			if t[k] < poly.vertCount {
				v[k] = rcGetVert(tile.verts, poly.verts[t[k]])
			} else {
				v[k] = rcGetVert(tile.verts, pd.vertBase+(t[k]-poly.vertCount))
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
	closest := closestPointOnDetailEdges(false, tile, poly, pos)
	height = closest[1]
	return height, true

}

// / @par
// /
// / @warning Only use this function if it is known that the provided polygon
// / reference is valid. This function is faster than #getTileAndPolyByRef, but
// / it does not validate the reference.
func (mesh *dtNavMesh) getTileAndPolyByRefUnsafe(ref dtPolyRef) (tile *dtMeshTile, poly *dtPoly) {
	_, it, ip := mesh.decodePolyId(ref)
	tile = mesh.m_tiles[it]
	poly = mesh.m_tiles[it].polys[ip]
	return
}
func (mesh *dtNavMesh) getTileAndPolyByRef(ref dtPolyRef) (tile *dtMeshTile, poly *dtPoly, status dtStatus) {
	if ref == 0 {
		return tile, poly, DT_FAILURE
	}

	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return tile, poly, DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return tile, poly, DT_FAILURE | DT_INVALID_PARAM
	}
	if ip >= mesh.m_tiles[it].header.polyCount {
		return tile, poly, DT_FAILURE | DT_INVALID_PARAM
	}
	tile = mesh.m_tiles[it]
	poly = mesh.m_tiles[it].polys[ip]
	return tile, poly, DT_SUCCESS
}

func closestPointOnDetailEdges(onlyBoundary bool, tile *dtMeshTile, poly *dtPoly, pos []float64) (closest []float64) {
	closest = make([]float64, 3)
	ip := getPolyIndexByTitle(tile, poly)
	pd := tile.detailMeshes[ip]
	dmin := math.MaxFloat64
	tmin := float64(0)
	pmin := []float64{}
	pmax := []float64{}

	for i := 0; i < pd.triCount; i++ {
		tris := rcGetTris(tile.detailTris, pd.triBase+i)
		ANY_BOUNDARY_EDGE := (DT_DETAIL_EDGE_BOUNDARY << 0) | (DT_DETAIL_EDGE_BOUNDARY << 2) | (DT_DETAIL_EDGE_BOUNDARY << 4)
		if onlyBoundary && (tris[3]&ANY_BOUNDARY_EDGE) == 0 {
			continue
		}
		var v [3][]float64
		for j := 0; j < 3; j++ {
			if tris[j] < poly.vertCount {
				v[j] = rcGetVert(tile.verts, poly.verts[tris[j]])
			} else {
				v[j] = rcGetVert(tile.verts, pd.vertBase+(tris[j]-poly.vertCount))
			}
		}
		k := 0
		j := 2
		for k < 3 {
			if (dtGetDetailTriEdgeFlags(tris[3], j)&DT_DETAIL_EDGE_BOUNDARY) == 0 && (onlyBoundary || tris[j] < tris[k]) {
				// Only looking at boundary edges and this is internal, or
				// this is an inner edge that we will see again or have already seen.
				continue
			}

			t, d := dtDistancePtSegSqr2D(pos, v[j], v[k])
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

	closest = dtVlerp(pmin, pmax, tmin)
	return closest
}

// / @par
// /
// / All points are projected onto the xz-plane, so the y-values are ignored.
func dtPointInPolygon(pt, verts []float64, nverts int) bool {
	// TODO: Replace pnpoly with triArea2D tests?
	i := 0
	j := nverts - 1
	c := false
	for i < nverts {
		vi := rcGetVert(verts, i)
		vj := rcGetVert(verts, j)
		if ((vi[2] > pt[2]) != (vj[2] > pt[2])) && (pt[0] < (vj[0]-vi[0])*(pt[2]-vi[2])/(vj[2]-vi[2])+vi[0]) {
			c = !c
		}

		j = i
		i++
	}
	return c
}

type dtTileState struct {
	magic   int       // Magic number, used to identify the data.
	version int       // Data version number.
	ref     dtTileRef // Tile ref at the time of storing the data.
}

type dtPolyState struct {
	flags int // Flags (see dtPolyFlags).
	area  int // Area ID of the polygon.
}

// / @par
// /
// / Tile state includes non-structural data such as polygon flags, area ids, etc.
// / @note This function does not impact the tile's #dtTileRef and #dtPolyRef's.
// / @see #storeTileState
func (mesh *dtNavMesh) restoreTileState(tile *dtMeshTile, tileState *dtTileState, polyStates []*dtPolyState, maxDataSize int) dtStatus {
	// Make sure there is enough space to store the state.
	// Check that the restore is possible.
	if tileState.magic != DT_NAVMESH_STATE_MAGIC {
		return DT_FAILURE | DT_WRONG_MAGIC
	}

	if tileState.version != DT_NAVMESH_STATE_VERSION {
		return DT_FAILURE | DT_WRONG_VERSION
	}

	if tileState.ref != mesh.getTileRef(tile) {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	// Restore per poly state.
	for i := 0; i < tile.header.polyCount; i++ {
		p := tile.polys[i]
		s := polyStates[i]
		p.flags = s.flags
		p.setArea(s.area)
	}

	return DT_SUCCESS
}

// / @par
// /
// / The add operation will fail if the data is in the wrong format, the allocated tile
// / space is full, or there is a tile already at the specified reference.
// /
// / The lastRef parameter is used to restore a tile with the same tile
// / reference it had previously used.  In this case the #dtPolyRef's for the
// / tile will be restored to the same values they were before the tile was
// / removed.
// /
// / The nav mesh assumes exclusive access to the data passed and will make
// / changes to the dynamic portion of the data. For that reason the data
// / should not be reused in other nav meshes until the tile has been successfully
// / removed from this nav mesh.
// /
// / @see dtCreateNavMeshData, #removeTile
func (mesh *dtNavMesh) addTile(header *dtMeshHeader, titleData *dtMeshTile, dataSize int, flags int, lastRef dtTileRef) (result dtTileRef, status dtStatus) {
	// Make sure the data is in right format.

	if header.magic != DT_NAVMESH_MAGIC {
		return result, DT_FAILURE | DT_WRONG_MAGIC
	}

	if header.version != DT_NAVMESH_VERSION {
		return result, DT_FAILURE | DT_WRONG_VERSION
	}

	if DT_POLYREF64 == 1 {
		// Do not allow adding more polygons than specified in the NavMesh's maxPolys constraint.
		// Otherwise, the poly ID cannot be represented with the given number of bits.
		if mesh.m_polyBits < dtIlog2(dtNextPow2(header.polyCount)) {
			return result, DT_FAILURE | DT_INVALID_PARAM
		}

	}
	// Make sure the location is free.
	if mesh.getTileAt(header.x, header.y, header.layer) != nil {
		return result, DT_FAILURE | DT_ALREADY_OCCUPIED
	}

	var tile *dtMeshTile
	if lastRef == 0 {
		if mesh.m_nextFree != nil {
			tile = mesh.m_nextFree
			mesh.m_nextFree = tile.next
			tile.next = nil
		}
	} else {
		// Try to relocate the tile to specific index with same salt.
		tileIndex := mesh.decodePolyIdTile(dtPolyRef(lastRef))
		if tileIndex >= mesh.m_maxTiles {
			return result, DT_FAILURE | DT_OUT_OF_MEMORY
		}

		// Try to find the specific tile id from the free list.
		target := mesh.m_tiles[tileIndex]
		var prev *dtMeshTile
		tile = mesh.m_nextFree
		for tile != nil && tile != target {
			prev = tile
			tile = tile.next
		}
		// Could not find the correct location.
		if tile != target {
			return result, DT_FAILURE | DT_OUT_OF_MEMORY
		}

		// Remove from freelist
		if prev == nil {
			mesh.m_nextFree = tile.next
		} else {
			prev.next = tile.next
		}

		// Restore salt.
		tile.salt = mesh.decodePolyIdSalt(dtPolyRef(lastRef))
	}

	// Make sure we could allocate a tile.
	if tile == nil {
		return result, DT_FAILURE | DT_OUT_OF_MEMORY
	}

	// Insert tile into the position lut.
	h := computeTileHash(header.x, header.y, mesh.m_tileLutMask)
	tile.next = mesh.m_posLookup[h]
	mesh.m_posLookup[h] = tile

	// Patch header pointers.
	tile.verts = titleData.verts
	tile.polys = titleData.polys
	tile.links = titleData.links
	tile.detailMeshes = titleData.detailMeshes
	tile.detailVerts = titleData.detailVerts
	tile.detailTris = titleData.detailTris
	tile.bvTree = titleData.bvTree
	tile.offMeshCons = titleData.offMeshCons

	// If there are no items in the bvtree, reset the tree pointer.
	// Build links freelist
	tile.linksFreeList = 0
	tile.links[header.maxLinkCount-1].next = DT_NULL_LINK
	for i := 0; i < header.maxLinkCount-1; i++ {
		tile.links[i].next = i + 1
	}

	// Init tile.
	tile.header = header
	//tile.data = data
	tile.dataSize = dataSize
	tile.flags = flags

	mesh.connectIntLinks(tile)

	// Base off-mesh connections to their starting polygons and connect connections inside the tile.
	mesh.baseOffMeshLinks(tile)
	mesh.connectExtOffMeshLinks(tile, tile, -1)

	// Create connections with neighbour tiles.
	MAX_NEIS := 32

	// Connect with layers in current tile.
	neis, nneis := mesh.getTilesAt(header.x, header.y, MAX_NEIS)
	for j := 0; j < nneis; j++ {
		if neis[j] == tile {
			continue
		}

		mesh.connectExtLinks(tile, neis[j], -1)
		mesh.connectExtLinks(neis[j], tile, -1)
		mesh.connectExtOffMeshLinks(tile, neis[j], -1)
		mesh.connectExtOffMeshLinks(neis[j], tile, -1)
	}

	// Connect with neighbour tiles.
	for i := 0; i < 8; i++ {
		neis, nneis := mesh.getNeighbourTilesAt(header.x, header.y, i, MAX_NEIS)
		for j := 0; j < nneis; j++ {
			mesh.connectExtLinks(tile, neis[j], i)
			mesh.connectExtLinks(neis[j], tile, dtOppositeTile(i))
			mesh.connectExtOffMeshLinks(tile, neis[j], i)
			mesh.connectExtOffMeshLinks(neis[j], tile, dtOppositeTile(i))
		}
	}

	result = mesh.getTileRef(tile)

	return result, DT_SUCCESS
}

func (mesh *dtNavMesh) baseOffMeshLinks(tile *dtMeshTile) {
	if tile == nil {
		return
	}

	base := mesh.getPolyRefBase(tile)

	// Base off-mesh connection start points.
	for i := 0; i < tile.header.offMeshConCount; i++ {
		con := tile.offMeshCons[i]
		poly := tile.polys[con.poly]

		halfExtents := [3]float64{con.rad, tile.header.walkableClimb, con.rad}

		// Find polygon to connect to.
		p := con.pos[0:3] // First vertex
		nearestPt, ref := mesh.findNearestPolyInTile(tile, p, halfExtents[:])
		if ref == 0 {
			continue
		}
		// findNearestPoly may return too optimistic results, further check to make sure.
		if dtSqr(nearestPt[0]-p[0])+dtSqr(nearestPt[2]-p[2]) > dtSqr(con.rad) {
			continue
		}

		// Make sure the location is on current mesh.
		copy(tile.verts[poly.verts[0]*3:poly.verts[0]*3+3], nearestPt)

		// Link off-mesh connection to target poly.
		idx := allocLink(tile)
		if idx != DT_NULL_LINK {
			link := tile.links[idx]
			link.ref = ref
			link.edge = 0
			link.side = 0xff
			link.bmax = 0
			link.bmin = link.bmax
			// Add to linked list.
			link.next = poly.firstLink
			poly.firstLink = idx
		}

		// Start end-point is always connect back to off-mesh connection.
		tidx := allocLink(tile)
		if tidx != DT_NULL_LINK {
			landPolyIdx := mesh.decodePolyIdPoly(ref)
			landPoly := tile.polys[landPolyIdx]
			link := tile.links[tidx]
			link.ref = base | dtPolyRef(con.poly)
			link.edge = 0xff
			link.side = 0xff
			link.bmax = 0
			link.bmin = link.bmax
			// Add to linked list.
			link.next = landPoly.firstLink
			landPoly.firstLink = tidx
		}
	}
}

func (mesh *dtNavMesh) connectIntLinks(tile *dtMeshTile) {
	if tile == nil {
		return
	}

	base := mesh.getPolyRefBase(tile)

	for i := 0; i < tile.header.polyCount; i++ {
		poly := tile.polys[i]
		poly.firstLink = DT_NULL_LINK

		if poly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
			continue
		}

		// Build edge links backwards so that the links will be
		// in the linked list from lowest index to highest.
		for j := poly.vertCount - 1; j >= 0; j-- {
			// Skip hard and non-internal edges.
			if poly.neis[j] == 0 || (poly.neis[j]&DT_EXT_LINK) > 0 {
				continue
			}

			idx := allocLink(tile)
			if idx != DT_NULL_LINK {
				link := tile.links[idx]
				link.ref = base | (dtPolyRef)(poly.neis[j]-1)
				link.edge = j
				link.side = 0xff
				link.bmax = 0
				link.bmin = link.bmax
				// Add to linked list.
				link.next = poly.firstLink
				poly.firstLink = idx
			}
		}
	}
}

func (mesh *dtNavMesh) getTileAt(x, y, layer int) *dtMeshTile {
	// Find tile based on hash.
	h := computeTileHash(x, y, mesh.m_tileLutMask)
	tile := mesh.m_posLookup[h]
	for tile != nil {
		if tile.header != nil && tile.header.x == x && tile.header.y == y && tile.header.layer == layer {
			return tile
		}
		tile = tile.next
	}
	return nil
}

// / @par
// /
// / This function returns the data for the tile so that, if desired,
// / it can be added back to the navigation mesh at a later point.
// /
// / @see #addTile
func (mesh *dtNavMesh) removeTile(ref dtTileRef) (data []int, dataSize int, status dtStatus) {
	if ref == 0 {
		return data, dataSize, DT_FAILURE | DT_INVALID_PARAM
	}

	tileIndex := mesh.decodePolyIdTile((dtPolyRef(ref)))
	tileSalt := mesh.decodePolyIdSalt((dtPolyRef(ref)))
	if tileIndex >= mesh.m_maxTiles {
		return data, dataSize, DT_FAILURE | DT_INVALID_PARAM
	}

	tile := mesh.m_tiles[tileIndex]
	if tile.salt != tileSalt {
		return data, dataSize, DT_FAILURE | DT_INVALID_PARAM
	}

	// Remove tile from hash lookup.
	h := computeTileHash(tile.header.x, tile.header.y, mesh.m_tileLutMask)
	var prev *dtMeshTile
	cur := mesh.m_posLookup[h]
	for cur != nil {
		if cur == tile {
			if prev != nil {
				prev.next = cur.next
			} else {
				mesh.m_posLookup[h] = cur.next
			}

			break
		}
		prev = cur
		cur = cur.next
	}
	MAX_NEIS := 32
	// Remove connections to neighbour tiles.
	// Disconnect from other layers in current tile.
	neis, nneis := mesh.getTilesAt(tile.header.x, tile.header.y, MAX_NEIS)
	for j := 0; j < nneis; j++ {
		if neis[j] == tile {
			continue
		}
		mesh.unconnectLinks(neis[j], tile)
	}

	// Disconnect from neighbour tiles.
	for i := 0; i < 8; i++ {
		neis, nneis = mesh.getNeighbourTilesAt(tile.header.x, tile.header.y, i, MAX_NEIS)
		for j := 0; j < nneis; j++ {
			mesh.unconnectLinks(neis[j], tile)
		}

	}

	// Reset tile.
	if tile.flags&DT_TILE_FREE_DATA > 0 {
		// Owns data
		tile.data = tile.data[:]
		tile.dataSize = 0
		if len(data) > 0 {
			data = []int{}
		}
		if dataSize > 0 {
			dataSize = 0
		}
	} else {
		if len(data) > 0 {
			data = tile.data
		}
		if dataSize > 0 {
			dataSize = tile.dataSize
		}
	}

	tile.header = nil
	tile.flags = 0
	tile.linksFreeList = 0
	tile.polys = tile.polys[:]
	tile.verts = tile.verts[:]
	tile.links = tile.links[:]
	tile.detailMeshes = tile.detailMeshes[:]
	tile.detailVerts = tile.detailVerts[:]
	tile.detailTris = tile.detailTris[:]
	tile.bvTree = tile.bvTree[:]
	tile.offMeshCons = tile.offMeshCons[:]

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
	tile.next = mesh.m_nextFree
	mesh.m_nextFree = tile

	return data, dataSize, DT_SUCCESS
}

func (mesh *dtNavMesh) getNeighbourTilesAt(x, y, side, maxTiles int) ([]*dtMeshTile, int) {
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

	return mesh.getTilesAt(nx, ny, maxTiles)
}

func (mesh *dtNavMesh) getTilesAt(x, y int, maxTiles int) ([]*dtMeshTile, int) {
	MAX_NEIS := 32
	tiles := make([]*dtMeshTile, MAX_NEIS)
	n := 0
	// Find tile based on hash.
	h := computeTileHash(x, y, mesh.m_tileLutMask)
	tile := mesh.m_posLookup[h]
	for tile != nil {
		if tile.header != nil && tile.header.x == x && tile.header.y == y {
			if n < maxTiles {
				tiles[n] = tile
				n++
			}

		}
		tile = tile.next
	}

	return tiles, n
}
func (mesh *dtNavMesh) isValidPolyRef(ref dtPolyRef) bool {
	if ref == 0 {
		return false
	}
	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return false
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return false
	}
	if ip >= mesh.m_tiles[it].header.polyCount {
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
func (mesh *dtNavMesh) getOffMeshConnectionPolyEndPoints(prevRef dtPolyRef, polyRef dtPolyRef, startPos []float64, endPos []float64) dtStatus {

	if polyRef == 0 {
		return DT_FAILURE
	}

	// Get current polygon
	salt, it, ip := mesh.decodePolyId(polyRef)
	if it >= mesh.m_maxTiles {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if ip >= tile.header.polyCount {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.polys[ip]

	// Make sure that the current poly is indeed off-mesh link.
	if poly.getType() != DT_POLYTYPE_OFFMESH_CONNECTION {
		return DT_FAILURE
	}

	// Figure out which way to hand out the vertices.
	idx0 := 0
	idx1 := 1

	// Find link that points to first vertex.
	for i := poly.firstLink; i != DT_NULL_LINK; i = tile.links[i].next {
		if tile.links[i].edge == 0 {
			if tile.links[i].ref != prevRef {
				idx0 = 1
				idx1 = 0
			}
			break
		}
	}

	startPos = dtVcopy(rcGetVert(tile.verts, poly.verts[idx0]))
	endPos = dtVcopy(rcGetVert(tile.verts, poly.verts[idx1]))

	return DT_SUCCESS
}

func (mesh *dtNavMesh) getOffMeshConnectionByRef(ref dtPolyRef) *dtOffMeshConnection {

	if ref == 0 {
		return nil
	}

	// Get current polygon
	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return nil
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return nil
	}
	tile := mesh.m_tiles[it]
	if ip >= tile.header.polyCount {
		return nil
	}
	poly := tile.polys[ip]

	// Make sure that the current poly is indeed off-mesh link.
	if poly.getType() != DT_POLYTYPE_OFFMESH_CONNECTION {
		return nil
	}

	idx := ip - tile.header.offMeshBase
	if idx < tile.header.offMeshConCount {
		panic("not expect")
	}
	return tile.offMeshCons[idx]
}

func (mesh *dtNavMesh) setPolyFlags(ref dtPolyRef, flags int) dtStatus {
	if ref == 0 {
		return DT_FAILURE
	}
	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if ip >= tile.header.polyCount {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.polys[ip]
	// Change flags.
	poly.flags = flags
	return DT_SUCCESS
}

func (mesh *dtNavMesh) getPolyFlags(ref dtPolyRef) (resultFlags int, status dtStatus) {
	if ref == 0 {
		return resultFlags, DT_FAILURE
	}
	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return resultFlags, DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return resultFlags, DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if ip >= tile.header.polyCount {
		return resultFlags, DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.polys[ip]

	resultFlags = poly.flags

	return resultFlags, DT_SUCCESS
}

func (mesh *dtNavMesh) setPolyArea(ref dtPolyRef, area int) dtStatus {
	if ref == 0 {
		return DT_FAILURE
	}
	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if ip >= tile.header.polyCount {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.polys[ip]

	poly.setArea(area)

	return DT_SUCCESS
}
func (mesh *dtNavMesh) getPolyArea(ref dtPolyRef) (resultArea int, status dtStatus) {
	if ref == 0 {
		return resultArea, DT_FAILURE
	}
	salt, it, ip := mesh.decodePolyId(ref)
	if it >= mesh.m_maxTiles {
		return resultArea, DT_FAILURE | DT_INVALID_PARAM
	}
	if mesh.m_tiles[it].salt != salt || mesh.m_tiles[it].header == nil {
		return resultArea, DT_FAILURE | DT_INVALID_PARAM
	}
	tile := mesh.m_tiles[it]
	if ip >= tile.header.polyCount {
		return resultArea, DT_FAILURE | DT_INVALID_PARAM
	}
	poly := tile.polys[ip]

	resultArea = poly.getArea()

	return resultArea, DT_SUCCESS
}
