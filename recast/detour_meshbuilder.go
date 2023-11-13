package recast

import (
	"math"
	"sort"
)

// / Represents the source data used to build an navigation mesh tile.
// / @ingroup detour
type dtNavMeshCreateParams struct {

	/// @name Polygon Mesh Attributes
	/// Used to create the base navigation graph.
	/// See #rcPolyMesh for details related to these attributes.
	/// @{

	verts     []int ///< The polygon mesh vertices. [(x, y, z) * #vertCount] [Unit: vx]
	vertCount int   ///< The number vertices in the polygon mesh. [Limit: >= 3]
	polys     []int ///< The polygon data. [Size: #polyCount * 2 * #nvp]
	polyFlags []int ///< The user defined flags assigned to each polygon. [Size: #polyCount]
	polyAreas []int ///< The user defined area ids assigned to each polygon. [Size: #polyCount]
	polyCount int   ///< Number of polygons in the mesh. [Limit: >= 1]
	nvp       int   ///< Number maximum number of vertices per polygon. [Limit: >= 3]

	/// @}
	/// @name Height Detail Attributes (Optional)
	/// See #rcPolyMeshDetail for details related to these attributes.
	/// @{

	detailMeshes     []int     ///< The height detail sub-mesh data. [Size: 4 * #polyCount]
	detailVerts      []float64 ///< The detail mesh vertices. [Size: 3 * #detailVertsCount] [Unit: wu]
	detailVertsCount int       ///< The number of vertices in the detail mesh.
	detailTris       []int     ///< The detail mesh triangles. [Size: 4 * #detailTriCount]
	detailTriCount   int       ///< The number of triangles in the detail mesh.

	/// @}
	/// @name Off-Mesh Connections Attributes (Optional)
	/// Used to define a custom point-to-point edge within the navigation graph, an
	/// off-mesh connection is a user defined traversable connection made up to two vertices,
	/// at least one of which resides within a navigation mesh polygon.
	/// @{

	/// Off-mesh connection vertices. [(ax, ay, az, bx, by, bz) * #offMeshConCount] [Unit: wu]
	offMeshConVerts []float64
	/// Off-mesh connection radii. [Size: #offMeshConCount] [Unit: wu]
	offMeshConRad []float64
	/// User defined flags assigned to the off-mesh connections. [Size: #offMeshConCount]
	offMeshConFlags []int
	/// User defined area ids assigned to the off-mesh connections. [Size: #offMeshConCount]
	offMeshConAreas []int
	/// The permitted travel direction of the off-mesh connections. [Size: #offMeshConCount]
	///
	/// 0 = Travel only from endpoint A to endpoint B.<br/>
	/// #DT_OFFMESH_CON_BIDIR = Bidirectional travel.
	offMeshConDir []int
	/// The user defined ids of the off-mesh connection. [Size: #offMeshConCount]
	offMeshConUserID []int
	/// The number of off-mesh connections. [Limit: >= 0]
	offMeshConCount int

	/// @}
	/// @name Tile Attributes
	/// @note The tile grid/layer data can be left at zero if the destination is a single tile mesh.
	/// @{

	userId    int        ///< The user defined id of the tile.
	tileX     int        ///< The tile's x-grid location within the multi-tile destination mesh. (Along the x-axis.)
	tileY     int        ///< The tile's y-grid location within the multi-tile destination mesh. (Along the z-axis.)
	tileLayer int        ///< The tile's layer within the layered destination mesh. [Limit: >= 0] (Along the y-axis.)
	bmin      [3]float64 ///< The minimum bounds of the tile. [(x, y, z)] [Unit: wu]
	bmax      [3]float64 ///< The maximum bounds of the tile. [(x, y, z)] [Unit: wu]

	/// @}
	/// @name General Configuration Attributes
	/// @{

	walkableHeight float64 ///< The agent height. [Unit: wu]
	walkableRadius float64 ///< The agent radius. [Unit: wu]
	walkableClimb  float64 ///< The agent maximum traversable ledge. (Up/Down) [Unit: wu]
	cs             float64 ///< The xz-plane cell size of the polygon mesh. [Limit: > 0] [Unit: wu]
	ch             float64 ///< The y-axis cell height of the polygon mesh. [Limit: > 0] [Unit: wu]

	/// True if a bounding volume tree should be built for the tile.
	/// @note The BVTree is not normally needed for layered navigation meshes.
	buildBvTree bool

	/// @}
}

const MESH_NULL_IDX = 0xffff

type BVItem struct {
	bmin [3]int
	bmax [3]int
	i    int
}

func compareItemX(va, vb *BVItem) int {
	if va.bmin[0] < vb.bmin[0] {
		return -1
	}

	if va.bmin[0] > vb.bmin[0] {
		return 1
	}

	return 0
}

func compareItemY(va, vb *BVItem) int {

	if va.bmin[1] < vb.bmin[1] {
		return -1
	}

	if va.bmin[1] > vb.bmin[1] {
		return 1
	}

	return 0
}

func compareItemZ(va, vb *BVItem) int {

	if va.bmin[2] < vb.bmin[2] {
		return -1
	}

	if va.bmin[2] > vb.bmin[2] {
		return 1
	}

	return 0
}

func calcExtends(items []*BVItem, imin int, imax int) (bmin, bmax []int) {
	bmin, bmax = make([]int, 3), make([]int, 3)
	bmin[0] = items[imin].bmin[0]
	bmin[1] = items[imin].bmin[1]
	bmin[2] = items[imin].bmin[2]

	bmax[0] = items[imin].bmax[0]
	bmax[1] = items[imin].bmax[1]
	bmax[2] = items[imin].bmax[2]

	for i := imin + 1; i < imax; i++ {
		it := items[i]
		if it.bmin[0] < bmin[0] {
			bmin[0] = it.bmin[0]
		}
		if it.bmin[1] < bmin[1] {
			bmin[1] = it.bmin[1]
		}
		if it.bmin[2] < bmin[2] {
			bmin[2] = it.bmin[2]
		}

		if it.bmax[0] > bmax[0] {
			bmax[0] = it.bmax[0]
		}
		if it.bmax[1] > bmax[1] {
			bmax[1] = it.bmax[1]
		}
		if it.bmax[2] > bmax[2] {
			bmax[2] = it.bmax[2]
		}
	}
	return
}
func longestAxis(x, y, z int) int {
	axis := 0
	maxVal := x
	if y > maxVal {
		axis = 1
		maxVal = y
	}
	if z > maxVal {
		axis = 2
	}
	return axis
}

func classifyOffMeshPoint(pt, bmin, bmax []float64) int {
	XP := 1 << 0
	ZP := 1 << 1
	XM := 1 << 2
	ZM := 1 << 3

	outcode := 0
	if pt[0] >= bmax[0] {
		outcode |= XP
	} else {
		outcode |= 0
	}
	if pt[2] >= bmax[2] {
		outcode |= ZP
	} else {
		outcode |= 0
	}

	if pt[0] < bmin[0] {
		outcode |= XM
	} else {
		outcode |= 0
	}
	if pt[2] < bmin[2] {
		outcode |= ZM
	} else {
		outcode |= 0
	}

	switch outcode {
	case XP:
		return 0
	case XP | ZP:
		return 1
	case ZP:
		return 2
	case XM | ZP:
		return 3
	case XM:
		return 4
	case XM | ZM:
		return 5
	case ZM:
		return 6
	case XP | ZM:
		return 7
	}

	return 0xff
}

func subdivide(items []*BVItem, nitems int, imin, imax int, curNode *int, nodes []*dtBVNode) {
	inum := imax - imin
	icur := curNode

	node := nodes[*curNode]
	*curNode++
	if inum == 1 {
		// Leaf
		node.bmin[0] = items[imin].bmin[0]
		node.bmin[1] = items[imin].bmin[1]
		node.bmin[2] = items[imin].bmin[2]

		node.bmax[0] = items[imin].bmax[0]
		node.bmax[1] = items[imin].bmax[1]
		node.bmax[2] = items[imin].bmax[2]

		node.i = items[imin].i
	} else {
		// Split
		bmin, bmax := calcExtends(items, imin, imax)
		copy(node.bmin[:], bmin)
		copy(node.bmax[:], bmax)
		axis := longestAxis(node.bmax[0]-node.bmin[0],
			node.bmax[1]-node.bmin[1],
			node.bmax[2]-node.bmin[2])

		if axis == 0 {
			sort.Slice(items[imin:inum], func(i, j int) bool {
				return compareItemX(items[imin+i], items[imin+j]) == -1
			})
			// Sort along x-axis

		} else if axis == 1 {
			// Sort along y-axis
			sort.Slice(items[imin:inum], func(i, j int) bool {
				return compareItemY(items[imin+i], items[imin+j]) == -1
			})
		} else {
			// Sort along z-axis
			sort.Slice(items[imin:inum], func(i, j int) bool {
				return compareItemZ(items[imin+i], items[imin+j]) == -1
			})
		}

		isplit := imin + inum/2

		// Left
		subdivide(items, nitems, imin, isplit, curNode, nodes)
		// Right
		subdivide(items, nitems, isplit, imax, curNode, nodes)

		iescape := *curNode - *icur
		// Negative index means escape.
		node.i = -iescape
	}
}

func createBVTree(params *dtNavMeshCreateParams, nodes []*dtBVNode) int {
	// Build tree
	quantFactor := 1 / params.cs
	var items = make([]*BVItem, params.polyCount)
	for i := 0; i < params.polyCount; i++ {
		it := items[i]
		it.i = i
		// Calc polygon bounds. Use detail meshes if available.
		if len(params.detailMeshes) > 0 {
			vb := params.detailMeshes[i*4+0]
			ndv := params.detailMeshes[i*4+1]
			var bmin [3]float64
			var bmax [3]float64

			dv := params.detailVerts[vb*3 : vb*3+3]
			copy(bmin[:], dv[:])
			copy(bmax[:], dv[:])

			for j := 1; j < ndv; j++ {
				copy(bmin[:], dv[j*3:j*3+3])
				copy(bmax[:], dv[j*3:j+3+3])
			}

			// BV-tree uses cs for all dimensions
			it.bmin[0] = int(dtClamp(((bmin[0] - params.bmin[0]) * quantFactor), float64(0), float64(0xffff)))
			it.bmin[1] = int(dtClamp(((bmin[1] - params.bmin[1]) * quantFactor), float64(0), float64(0xffff)))
			it.bmin[2] = int(dtClamp(((bmin[2] - params.bmin[2]) * quantFactor), float64(0), float64(0xffff)))

			it.bmax[0] = int(dtClamp(((bmax[0] - params.bmin[0]) * quantFactor), float64(0), float64(0xffff)))
			it.bmax[1] = int(dtClamp(((bmax[1] - params.bmin[1]) * quantFactor), float64(0), float64(0xffff)))
			it.bmax[2] = int(dtClamp(((bmax[2] - params.bmin[2]) * quantFactor), float64(0), float64(0xffff)))
		} else {
			p := params.polys[i*params.nvp*2 : i*params.nvp*2+2]
			it.bmax[0] = params.verts[p[0]*3+0]
			it.bmin[0] = it.bmax[0]
			it.bmax[1] = params.verts[p[0]*3+1]
			it.bmin[1] = it.bmax[1]
			it.bmax[2] = params.verts[p[0]*3+2]
			it.bmin[2] = it.bmax[2]

			for j := 1; j < params.nvp; j++ {
				if p[j] == MESH_NULL_IDX {
					break
				}
				x := params.verts[p[j]*3+0]
				y := params.verts[p[j]*3+1]
				z := params.verts[p[j]*3+2]

				if x < it.bmin[0] {
					it.bmin[0] = x
				}
				if y < it.bmin[1] {
					it.bmin[1] = y
				}
				if z < it.bmin[2] {
					it.bmin[2] = z
				}

				if x > it.bmax[0] {
					it.bmax[0] = x
				}
				if y > it.bmax[1] {
					it.bmax[1] = y
				}
				if z > it.bmax[2] {
					it.bmax[2] = z
				}
			}
			// Remap y
			it.bmin[1] = int(math.Floor(float64(it.bmin[1]) * params.ch / params.cs))
			it.bmax[1] = int(math.Ceil(float64(it.bmax[1]) * params.ch / params.cs))
		}
	}

	curNode := 0
	subdivide(items, params.polyCount, 0, params.polyCount, &curNode, nodes)
	return curNode
}

// / @par
// /
// / The output data array is allocated using the detour allocator (dtAlloc()).  The method
// / used to free the memory will be determined by how the tile is added to the navigation
// / mesh.
// /
// / @see dtNavMesh, dtNavMesh::addTile()
func dtCreateNavMeshData(params *dtNavMeshCreateParams) bool {
	if params.nvp > DT_VERTS_PER_POLYGON {
		return false
	}

	if params.vertCount >= 0xffff {
		return false
	}

	if params.vertCount == 0 || len(params.verts) == 0 {
		return false
	}

	if params.polyCount == 0 || len(params.polys) == 0 {
		return false
	}
	nvp := params.nvp
	offMeshConClass := make([]int, params.offMeshConCount*2)
	// Classify off-mesh connection points. We store only the connections
	// whose start point is inside the tile.
	storedOffMeshConCount := 0
	offMeshConLinkCount := 0

	if params.offMeshConCount > 0 {
		// Find tight heigh bounds, used for culling out off-mesh start locations.
		hmin := math.MaxFloat64
		hmax := math.SmallestNonzeroFloat64

		if len(params.detailVerts) > 0 && params.detailVertsCount > 0 {
			for i := 0; i < params.detailVertsCount; i++ {
				h := params.detailVerts[i*3+1]
				hmin = dtMin(hmin, h)
				hmax = dtMax(hmax, h)
			}
		} else {
			for i := 0; i < params.vertCount; i++ {
				iv := params.verts[i*3 : i*3+3]
				h := params.bmin[1] + float64(iv[1])*params.ch
				hmin = dtMin(hmin, h)
				hmax = dtMax(hmax, h)
			}
		}
		hmin -= params.walkableClimb
		hmax += params.walkableClimb
		var bmin, bmax [3]float64
		copy(bmin[:], params.bmin[:])
		copy(bmax[:], params.bmax[:])
		bmin[1] = hmin
		bmax[1] = hmax

		for i := 0; i < params.offMeshConCount; i++ {
			p0 := rcGetVert(params.offMeshConVerts, (i*2 + 0))
			p1 := rcGetVert(params.offMeshConVerts, (i*2 + 1))
			offMeshConClass[i*2+0] = classifyOffMeshPoint(p0, bmin[:], bmax[:])
			offMeshConClass[i*2+1] = classifyOffMeshPoint(p1, bmin[:], bmax[:])

			// Zero out off-mesh start positions which are not even potentially touching the mesh.
			if offMeshConClass[i*2+0] == 0xff {
				if p0[1] < bmin[1] || p0[1] > bmax[1] {
					offMeshConClass[i*2+0] = 0
				}

			}

			// Cound how many links should be allocated for off-mesh connections.
			if offMeshConClass[i*2+0] == 0xff {
				offMeshConLinkCount++
			}

			if offMeshConClass[i*2+1] == 0xff {
				offMeshConLinkCount++
			}

			if offMeshConClass[i*2+0] == 0xff {
				storedOffMeshConCount++
			}

		}
	}

	// Off-mesh connections are stored as polygons, adjust values.
	totPolyCount := params.polyCount + storedOffMeshConCount
	totVertCount := params.vertCount + storedOffMeshConCount*2

	// Find portal edges which are at tile borders.
	edgeCount := 0
	portalCount := 0
	for i := 0; i < params.polyCount; i++ {
		p := params.polys[i*2*nvp:]
		for j := 0; j < nvp; j++ {
			if p[j] == MESH_NULL_IDX {
				break
			}
			edgeCount++

			if p[nvp+j]&0x8000 > 0 {
				dir := p[nvp+j] & 0xf
				if dir != 0xf {
					portalCount++
				}

			}
		}
	}

	maxLinkCount := edgeCount + portalCount*2 + offMeshConLinkCount*2

	// Find unique detail vertices.
	uniqueDetailVertCount := 0
	detailTriCount := 0
	if len(params.detailMeshes) > 0 {
		// Has detail mesh, count unique detail vertex count and use input detail tri count.
		detailTriCount = params.detailTriCount
		for i := 0; i < params.polyCount; i++ {
			p := params.polys[i*nvp*2:]
			ndv := params.detailMeshes[i*4+1]
			nv := 0
			for j := 0; j < nvp; j++ {
				if p[j] == MESH_NULL_IDX {
					break
				}
				nv++
			}
			ndv -= nv
			uniqueDetailVertCount += ndv
		}
	} else {
		// No input detail mesh, build detail mesh from nav polys.
		uniqueDetailVertCount = 0 // No extra detail verts.
		detailTriCount = 0
		for i := 0; i < params.polyCount; i++ {
			p := params.polys[i*nvp*2:]
			nv := 0
			for j := 0; j < nvp; j++ {
				if p[j] == MESH_NULL_IDX {
					break
				}
				nv++
			}
			detailTriCount += nv - 2
		}
	}
	var header dtMeshHeader
	navVerts := make([]float64, 3*totVertCount)
	navPolys := make([]*dtPoly, totPolyCount)
	for index := range navPolys {
		navPolys[index] = &dtPoly{}
	}
	navDMeshes := make([]*dtPolyDetail, params.polyCount)
	for index := range navDMeshes {
		navDMeshes[index] = &dtPolyDetail{}
	}
	navDTris := make([]int, 4*detailTriCount)
	offMeshCons := make([]*dtOffMeshConnection, storedOffMeshConCount)
	for index := range offMeshCons {
		offMeshCons[index] = &dtOffMeshConnection{}
	}
	navDVerts := make([]float64, 3*uniqueDetailVertCount)
	navBvtreeSize := 0
	if params.buildBvTree {
		navBvtreeSize = params.polyCount * 2
	}
	navBvtree := make([]*dtBVNode, navBvtreeSize)
	// Store header
	header.magic = DT_NAVMESH_MAGIC
	header.version = DT_NAVMESH_VERSION
	header.x = params.tileX
	header.y = params.tileY
	header.layer = params.tileLayer
	header.userId = params.userId
	header.polyCount = totPolyCount
	header.vertCount = totVertCount
	header.maxLinkCount = maxLinkCount
	copy(header.bmin[:], params.bmin[:])
	copy(header.bmax[:], params.bmax[:])
	header.detailMeshCount = params.polyCount
	header.detailVertCount = uniqueDetailVertCount
	header.detailTriCount = detailTriCount
	header.bvQuantFactor = 1.0 / params.cs
	header.offMeshBase = params.polyCount
	header.walkableHeight = params.walkableHeight
	header.walkableRadius = params.walkableRadius
	header.walkableClimb = params.walkableClimb
	header.offMeshConCount = storedOffMeshConCount
	if params.buildBvTree {
		header.bvNodeCount = params.polyCount * 2
	}

	offMeshVertsBase := params.vertCount
	offMeshPolyBase := params.polyCount
	// Mesh vertices
	for i := 0; i < params.vertCount; i++ {
		iv := params.verts[i*3 : i*3+3]
		v := rcGetVert(navVerts, i)
		v[0] = params.bmin[0] + float64(iv[0])*params.cs
		v[1] = params.bmin[1] + float64(iv[1])*params.ch
		v[2] = params.bmin[2] + float64(iv[2])*params.cs
	}
	// Off-mesh link vertices.
	n := 0
	for i := 0; i < params.offMeshConCount; i++ {
		// Only store connections which start from this tile.
		if offMeshConClass[i*2+0] == 0xff {
			linkv := params.offMeshConVerts[i*2 : i*2+3]
			v := rcGetVert(navVerts, (offMeshVertsBase + n*2))
			copy(v[:3], linkv[:3])
			copy(v[3:], linkv[3:])
			n++
		}
	}
	// Store polygons
	// Mesh polys
	src := params.polys
	for i := 0; i < params.polyCount; i++ {
		p := navPolys[i]
		p.vertCount = 0
		p.flags = params.polyFlags[i]
		p.setArea(params.polyAreas[i])
		p.setType(DT_POLYTYPE_GROUND)
		for j := 0; j < nvp; j++ {
			if src[j] == MESH_NULL_IDX {
				break
			}
			p.verts[j] = src[j]
			if src[nvp+j]&0x8000 > 0 {
				// Border or portal edge.
				dir := src[nvp+j] & 0xf
				if dir == 0xf {
					p.neis[j] = 0
				} else if dir == 0 { // Border
					p.neis[j] = DT_EXT_LINK | 4
				} else if dir == 1 { // Portal x-
					p.neis[j] = DT_EXT_LINK | 2
				} else if dir == 2 { // Portal z+
					p.neis[j] = DT_EXT_LINK | 0
				} else if dir == 3 { // Portal x+
					p.neis[j] = DT_EXT_LINK | 6
				}
			} else { // Portal z-
				// Normal connection
				p.neis[j] = src[nvp+j] + 1
			}

			p.vertCount++
		}
		src = src[nvp*2:]
	}
	// Off-mesh connection vertices.
	n = 0
	for i := 0; i < params.offMeshConCount; i++ {
		// Only store connections which start from this tile.
		if offMeshConClass[i*2+0] == 0xff {
			p := navPolys[offMeshPolyBase+n]
			p.vertCount = 2
			p.verts[0] = (offMeshVertsBase + n*2 + 0)
			p.verts[1] = (offMeshVertsBase + n*2 + 1)
			p.flags = params.offMeshConFlags[i]
			p.setArea(params.offMeshConAreas[i])
			p.setType(DT_POLYTYPE_OFFMESH_CONNECTION)
			n++
		}
	}

	// Store detail meshes and vertices.
	// The nav polygon vertices are stored as the first vertices on each mesh.
	// We compress the mesh data by skipping them and using the navmesh coordinates.
	if len(params.detailMeshes) > 0 {
		vbase := 0
		for i := 0; i < params.polyCount; i++ {
			dtl := navDMeshes[i]
			vb := params.detailMeshes[i*4+0]
			ndv := params.detailMeshes[i*4+1]
			nv := navPolys[i].vertCount
			dtl.vertBase = vbase
			dtl.vertCount = (ndv - nv)
			dtl.triBase = params.detailMeshes[i*4+2]
			dtl.triCount = params.detailMeshes[i*4+3]
			// Copy vertices except the first 'nv' verts which are equal to nav poly verts.
			if ndv-nv > 0 {
				copy(navDVerts[vbase*3:], params.detailVerts[(vb+nv)*3:(vb+nv)*3+3*(ndv-nv)])
				vbase += (ndv - nv)
			}
		}
		// Store triangles.

		copy(navDTris, params.detailTris[:4*params.detailTriCount])
	} else {
		// Create dummy detail mesh by triangulating polys.
		tbase := 0
		for i := 0; i < params.polyCount; i++ {
			dtl := navDMeshes[i]
			nv := navPolys[i].vertCount
			dtl.vertBase = 0
			dtl.vertCount = 0
			dtl.triBase = tbase
			dtl.triCount = (nv - 2)
			// Triangulate polygon (local indices).
			for j := 2; j < nv; j++ {
				t := navDTris[tbase*4 : tbase*4+4]
				t[0] = 0
				t[1] = (j - 1)
				t[2] = j
				// Bit for each edge that belongs to poly boundary.
				t[3] = (1 << 2)
				if j == 2 {
					t[3] |= (1 << 0)
				}
				if j == nv-1 {
					t[3] |= (1 << 4)
				}
				tbase++
			}
		}
	}

	// Store and create BVtree.
	if params.buildBvTree {
		createBVTree(params, navBvtree)
	}

	// Store Off-Mesh connections.
	n = 0
	for i := 0; i < params.offMeshConCount; i++ {
		// Only store connections which start from this tile.
		if offMeshConClass[i*2+0] == 0xff {
			con := offMeshCons[n]
			con.poly = (offMeshPolyBase + n)
			// Copy connection end-points.
			endPts := params.offMeshConVerts[i*2*3 : i*2*3]
			copy(con.pos[:3], endPts[:3])
			copy(con.pos[3:], endPts[3:])
			con.rad = params.offMeshConRad[i]
			if params.offMeshConDir[i] > 0 {
				con.flags = DT_OFFMESH_CON_BIDIR
			} else {
				con.flags = 0
			}

			con.side = offMeshConClass[i*2+1]
			if len(params.offMeshConUserID) > 0 {
				con.userId = params.offMeshConUserID[i]
			}

			n++
		}
	}
	return true

}
