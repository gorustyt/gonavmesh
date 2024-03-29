package detour

import (
	"github.com/gorustyt/gonavmesh/common"
	"math"
	"sort"
)

// / Represents the source data used to build an navigation mesh tile.
// / @ingroup detour
type DtNavMeshCreateParams struct {

	/// @name Polygon Mesh Attributes
	/// Used to create the base navigation graph.
	/// See #RcPolyMesh for details related to these attributes.
	/// @{

	Verts     []uint16 ///< The polygon mesh vertices. [(x, y, z) * #vertCount] [Unit: vx]
	VertCount int32    ///< The number vertices in the polygon mesh. [Limit: >= 3]
	Polys     []uint16 ///< The polygon data. [Size: #polyCount * 2 * #nvp]
	PolyFlags []uint16 ///< The user defined flags assigned to each polygon. [Size: #polyCount]
	PolyAreas []uint8  ///< The user defined area ids assigned to each polygon. [Size: #polyCount]
	PolyCount int32    ///< Number of polygons in the mesh. [Limit: >= 1]
	Nvp       int32    ///< Number maximum number of vertices per polygon. [Limit: >= 3]

	/// @}
	/// @name Height Detail Attributes (Optional)
	/// See #RcPolyMeshDetail for details related to these attributes.
	/// @{

	DetailMeshes     []uint32  ///< The height detail sub-mesh data. [Size: 4 * #polyCount]
	DetailVerts      []float32 ///< The detail mesh vertices. [Size: 3 * #detailVertsCount] [Unit: wu]
	DetailVertsCount int32     ///< The number of vertices in the detail mesh.
	DetailTris       []uint8   ///< The detail mesh triangles. [Size: 4 * #detailTriCount]
	DetailTriCount   int32     ///< The number of triangles in the detail mesh.

	/// @}
	/// @name Off-Mesh Connections Attributes (Optional)
	/// Used to define a custom point-to-point edge within the navigation graph, an
	/// off-mesh connection is a user defined traversable connection made up to two vertices,
	/// at least one of which resides within a navigation mesh polygon.
	/// @{

	/// Off-mesh connection vertices. [(ax, ay, az, bx, by, bz) * #offMeshConCount] [Unit: wu]
	OffMeshConVerts []float32
	/// Off-mesh connection radii. [Size: #offMeshConCount] [Unit: wu]
	OffMeshConRad []float32
	/// User defined flags assigned to the off-mesh connections. [Size: #offMeshConCount]
	OffMeshConFlags []uint16
	/// User defined area ids assigned to the off-mesh connections. [Size: #offMeshConCount]
	OffMeshConAreas []uint8
	/// The permitted travel direction of the off-mesh connections. [Size: #offMeshConCount]
	///
	/// 0 = Travel only from endpoint A to endpoint B.<br/>
	/// #DT_OFFMESH_CON_BIDIR = Bidirectional travel.
	OffMeshConDir []uint8
	/// The user defined ids of the off-mesh connection. [Size: #offMeshConCount]
	OffMeshConUserID []uint32
	/// The number of off-mesh connections. [Limit: >= 0]
	OffMeshConCount int32

	/// @}
	/// @name Tile Attributes
	/// @note The tile grid/layer data can be left at zero if the destination is a single tile mesh.
	/// @{

	UserId    uint32     ///< The user defined id of the tile.
	TileX     int32      ///< The tile's x-grid location within the multi-tile destination mesh. (Along the x-axis.)
	TileY     int32      ///< The tile's y-grid location within the multi-tile destination mesh. (Along the z-axis.)
	TileLayer int32      ///< The tile's layer within the layered destination mesh. [Limit: >= 0] (Along the y-axis.)
	Bmin      [3]float32 ///< The minimum bounds of the tile. [(x, y, z)] [Unit: wu]
	Bmax      [3]float32 ///< The maximum bounds of the tile. [(x, y, z)] [Unit: wu]

	/// @}
	/// @name General Configuration Attributes
	/// @{

	WalkableHeight float32 ///< The agent height. [Unit: wu]
	WalkableRadius float32 ///< The agent radius. [Unit: wu]
	WalkableClimb  float32 ///< The agent maximum traversable ledge. (Up/Down) [Unit: wu]
	Cs             float32 ///< The xz-plane cell size of the polygon mesh. [Limit: > 0] [Unit: wu]
	Ch             float32 ///< The y-axis cell height of the polygon mesh. [Limit: > 0] [Unit: wu]

	/// True if a bounding volume tree should be built for the tile.
	/// @note The BVTree is not normally needed for layered navigation meshes.
	BuildBvTree bool

	/// @}
}

const MESH_NULL_IDX = 0xffff

type BVItem struct {
	bmin [3]uint16
	bmax [3]uint16
	i    int32
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

func calcExtends(items []*BVItem, imin int32, imax int32) (bmin, bmax []uint16) {
	bmin, bmax = make([]uint16, 3), make([]uint16, 3)
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
func longestAxis(x, y, z uint16) int32 {
	axis := int32(0)
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

func classifyOffMeshPoint(pt, bmin, bmax []float32) uint8 {
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

func subdivide(items []*BVItem, nitems int32, imin, imax int32, curNode *int32, nodes []*DtBVNode) {
	inum := imax - imin
	icur := curNode

	node := nodes[*curNode]
	*curNode++
	if inum == 1 {
		// Leaf
		node.Bmin[0] = items[imin].bmin[0]
		node.Bmin[1] = items[imin].bmin[1]
		node.Bmin[2] = items[imin].bmin[2]

		node.Bmax[0] = items[imin].bmax[0]
		node.Bmax[1] = items[imin].bmax[1]
		node.Bmax[2] = items[imin].bmax[2]

		node.I = items[imin].i
	} else {
		// Split
		bmin, bmax := calcExtends(items, imin, imax)
		copy(node.Bmin[:], bmin)
		copy(node.Bmax[:], bmax)
		axis := longestAxis(node.Bmax[0]-node.Bmin[0],
			node.Bmax[1]-node.Bmin[1],
			node.Bmax[2]-node.Bmin[2])

		if axis == 0 {
			sort.Slice(items[imin:inum], func(i, j int) bool {
				return compareItemX(items[imin+int32(i)], items[imin+int32(j)]) == -1
			})
			// Sort along x-axis

		} else if axis == 1 {
			// Sort along y-axis
			sort.Slice(items[imin:inum], func(i, j int) bool {
				return compareItemY(items[imin+int32(i)], items[imin+int32(j)]) == -1
			})
		} else {
			// Sort along z-axis
			sort.Slice(items[imin:inum], func(i, j int) bool {
				return compareItemZ(items[imin+int32(i)], items[imin+int32(j)]) == -1
			})
		}

		isplit := imin + inum/2

		// Left
		subdivide(items, nitems, imin, isplit, curNode, nodes)
		// Right
		subdivide(items, nitems, isplit, imax, curNode, nodes)

		iescape := *curNode - *icur
		// Negative index means escape.
		node.I = -iescape
	}
}

func createBVTree(params *DtNavMeshCreateParams, nodes []*DtBVNode) int32 {
	// Build tree
	quantFactor := 1 / params.Cs
	var items = make([]*BVItem, params.PolyCount)
	for i := int32(0); i < params.PolyCount; i++ {
		it := items[i]
		it.i = i
		// Calc polygon bounds. Use detail meshes if available.
		if len(params.DetailMeshes) > 0 {
			vb := params.DetailMeshes[i*4+0]
			ndv := params.DetailMeshes[i*4+1]
			var bmin [3]float32
			var bmax [3]float32

			dv := params.DetailVerts[vb*3 : vb*3+3]
			copy(bmin[:], dv[:])
			copy(bmax[:], dv[:])

			for j := uint32(1); j < ndv; j++ {
				copy(bmin[:], dv[j*3:j*3+3])
				copy(bmax[:], dv[j*3:j+3+3])
			}

			// BV-tree uses cs for all dimensions
			it.bmin[0] = uint16(common.Clamp((bmin[0]-params.Bmin[0])*quantFactor, float32(0), float32(0xffff)))
			it.bmin[1] = uint16(common.Clamp((bmin[1]-params.Bmin[1])*quantFactor, float32(0), float32(0xffff)))
			it.bmin[2] = uint16(common.Clamp((bmin[2]-params.Bmin[2])*quantFactor, float32(0), float32(0xffff)))

			it.bmax[0] = uint16(common.Clamp((bmax[0]-params.Bmin[0])*quantFactor, float32(0), float32(0xffff)))
			it.bmax[1] = uint16(common.Clamp((bmax[1]-params.Bmin[1])*quantFactor, float32(0), float32(0xffff)))
			it.bmax[2] = uint16(common.Clamp((bmax[2]-params.Bmin[2])*quantFactor, float32(0), float32(0xffff)))
		} else {
			p := params.Polys[i*params.Nvp*2 : i*params.Nvp*2+2]
			it.bmax[0] = params.Verts[p[0]*3+0]
			it.bmin[0] = it.bmax[0]
			it.bmax[1] = params.Verts[p[0]*3+1]
			it.bmin[1] = it.bmax[1]
			it.bmax[2] = params.Verts[p[0]*3+2]
			it.bmin[2] = it.bmax[2]

			for j := int32(1); j < params.Nvp; j++ {
				if p[j] == MESH_NULL_IDX {
					break
				}
				x := params.Verts[p[j]*3+0]
				y := params.Verts[p[j]*3+1]
				z := params.Verts[p[j]*3+2]

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
			it.bmin[1] = uint16(math.Floor(float64(it.bmin[1]) * float64(params.Ch) / float64(params.Cs)))
			it.bmax[1] = uint16(math.Ceil(float64(it.bmax[1]) * float64(params.Ch) / float64(params.Cs)))
		}
	}

	curNode := int32(0)
	subdivide(items, params.PolyCount, 0, params.PolyCount, &curNode, nodes)
	return curNode
}

// / @par
// /
// / The output data array is allocated using the detour allocator (dtAlloc()).  The method
// / used to free the memory will be determined by how the tile is added to the navigation
// / mesh.
// /
// / @see DtNavMesh, DtNavMesh::addTile()
func DtCreateNavMeshData(params *DtNavMeshCreateParams) (outData *NavMeshData, ok bool) {
	if params.Nvp > DT_VERTS_PER_POLYGON {
		return outData, false
	}

	if params.VertCount >= 0xffff {
		return outData, false
	}

	if params.VertCount == 0 || len(params.Verts) == 0 {
		return outData, false
	}

	if params.PolyCount == 0 || len(params.Polys) == 0 {
		return outData, false
	}
	nvp := params.Nvp
	offMeshConClass := make([]uint8, params.OffMeshConCount*2)
	// Classify off-mesh connection points. We store only the connections
	// whose start point is inside the tile.
	storedOffMeshConCount := int32(0)
	offMeshConLinkCount := int32(0)

	if params.OffMeshConCount > 0 {
		// Find tight heigh bounds, used for culling out off-mesh start locations.
		hmin := float32(math.MaxFloat32)
		hmax := float32(math.SmallestNonzeroFloat32)

		if len(params.DetailVerts) > 0 && params.DetailVertsCount > 0 {
			for i := int32(0); i < params.DetailVertsCount; i++ {
				h := params.DetailVerts[i*3+1]
				hmin = min(hmin, h)
				hmax = max(hmax, h)
			}
		} else {
			for i := int32(0); i < params.VertCount; i++ {
				iv := params.Verts[i*3 : i*3+3]
				h := params.Bmin[1] + float32(iv[1])*params.Ch
				hmin = min(hmin, h)
				hmax = max(hmax, h)
			}
		}
		hmin -= params.WalkableClimb
		hmax += params.WalkableClimb
		var bmin, bmax [3]float32
		bmin = params.Bmin
		bmax = params.Bmax
		bmin[1] = hmin
		bmax[1] = hmax

		for i := int32(0); i < params.OffMeshConCount; i++ {
			p0 := common.GetVert3(params.OffMeshConVerts, i*2+0)
			p1 := common.GetVert3(params.OffMeshConVerts, i*2+1)
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
	totPolyCount := params.PolyCount + storedOffMeshConCount
	totVertCount := params.VertCount + storedOffMeshConCount*2

	// Find portal edges which are at tile borders.
	edgeCount := int32(0)
	portalCount := int32(0)
	for i := int32(0); i < params.PolyCount; i++ {
		p := params.Polys[i*2*nvp:]
		for j := int32(0); j < nvp; j++ {
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
	uniqueDetailVertCount := int32(0)
	detailTriCount := int32(0)
	if len(params.DetailMeshes) > 0 {
		// Has detail mesh, count unique detail vertex count and use input detail tri count.
		detailTriCount = params.DetailTriCount
		for i := int32(0); i < params.PolyCount; i++ {
			p := params.Polys[i*nvp*2:]
			ndv := int32(params.DetailMeshes[i*4+1])
			nv := int32(0)
			for j := int32(0); j < nvp; j++ {
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
		uniqueDetailVertCount = int32(0) // No extra detail verts.
		detailTriCount = int32(0)
		for i := int32(0); i < params.PolyCount; i++ {
			p := params.Polys[i*nvp*2:]
			nv := int32(0)
			for j := int32(0); j < nvp; j++ {
				if p[j] == MESH_NULL_IDX {
					break
				}
				nv++
			}
			detailTriCount += nv - 2
		}
	}
	var header DtMeshHeader
	navVerts := make([]float32, 3*totVertCount)
	navPolys := make([]*DtPoly, totPolyCount)
	for index := range navPolys {
		navPolys[index] = &DtPoly{}
	}
	navDMeshes := make([]*DtPolyDetail, params.PolyCount)
	for index := range navDMeshes {
		navDMeshes[index] = &DtPolyDetail{}
	}
	navDTris := make([]uint8, 4*detailTriCount)
	offMeshCons := make([]*DtOffMeshConnection, storedOffMeshConCount)
	for index := range offMeshCons {
		offMeshCons[index] = &DtOffMeshConnection{}
	}
	navDVerts := make([]float32, 3*uniqueDetailVertCount)
	navBvtreeSize := int32(0)
	if params.BuildBvTree {
		navBvtreeSize = params.PolyCount * 2
	}
	navBvtree := make([]*DtBVNode, navBvtreeSize)
	data := &NavMeshData{
		Header:      &header,
		NavVerts:    navVerts,
		NavPolys:    navPolys,
		NavDMeshes:  navDMeshes,
		NavDVerts:   navDVerts,
		NavBvtree:   navBvtree,
		NavDTris:    navDTris,
		OffMeshCons: offMeshCons,
	}
	// Store header
	header.Magic = DT_NAVMESH_MAGIC
	header.Version = DT_NAVMESH_VERSION
	header.X = params.TileX
	header.Y = params.TileY
	header.Layer = params.TileLayer
	header.UserId = params.UserId
	header.PolyCount = totPolyCount
	header.VertCount = totVertCount
	header.MaxLinkCount = maxLinkCount
	header.Bmin = params.Bmin
	header.Bmax = params.Bmax
	header.DetailMeshCount = params.PolyCount
	header.DetailVertCount = uniqueDetailVertCount
	header.DetailTriCount = detailTriCount
	header.BvQuantFactor = 1.0 / params.Cs
	header.OffMeshBase = params.PolyCount
	header.WalkableHeight = params.WalkableHeight
	header.WalkableRadius = params.WalkableRadius
	header.WalkableClimb = params.WalkableClimb
	header.OffMeshConCount = storedOffMeshConCount
	if params.BuildBvTree {
		header.BvNodeCount = params.PolyCount * 2
	}

	offMeshVertsBase := params.VertCount
	offMeshPolyBase := params.PolyCount
	// Mesh vertices
	for i := int32(0); i < params.VertCount; i++ {
		iv := params.Verts[i*3 : i*3+3]
		v := common.GetVert3(navVerts, i)
		v[0] = params.Bmin[0] + float32(iv[0])*params.Cs
		v[1] = params.Bmin[1] + float32(iv[1])*params.Ch
		v[2] = params.Bmin[2] + float32(iv[2])*params.Cs
	}
	// Off-mesh link vertices.
	n := int32(0)
	for i := int32(0); i < params.OffMeshConCount; i++ {
		// Only store connections which start from this tile.
		if offMeshConClass[i*2+0] == 0xff {
			linkv := params.OffMeshConVerts[i*2 : i*2+3]
			v := common.GetVert3(navVerts, offMeshVertsBase+n*2)
			copy(v[:3], linkv[:3])
			copy(v[3:], linkv[3:])
			n++
		}
	}
	// Store polygons
	// Mesh polys
	src := params.Polys
	for i := int32(0); i < params.PolyCount; i++ {
		p := navPolys[i]
		p.VertCount = 0
		p.Flags = params.PolyFlags[i]
		p.SetArea(params.PolyAreas[i])
		p.SetType(DT_POLYTYPE_GROUND)
		for j := int32(0); j < nvp; j++ {
			if src[j] == MESH_NULL_IDX {
				break
			}
			p.Verts[j] = src[j]
			if src[nvp+j]&0x8000 > 0 {
				// Border or portal edge.
				dir := src[nvp+j] & 0xf
				if dir == 0xf {
					p.Neis[j] = 0
				} else if dir == 0 { // Border
					p.Neis[j] = DT_EXT_LINK | 4
				} else if dir == 1 { // Portal x-
					p.Neis[j] = DT_EXT_LINK | 2
				} else if dir == 2 { // Portal z+
					p.Neis[j] = DT_EXT_LINK | 0
				} else if dir == 3 { // Portal x+
					p.Neis[j] = DT_EXT_LINK | 6
				}
			} else { // Portal z-
				// Normal connection
				p.Neis[j] = src[nvp+j] + 1
			}

			p.VertCount++
		}
		src = src[nvp*2:]
	}
	// Off-mesh connection vertices.
	n = int32(0)
	for i := int32(0); i < params.OffMeshConCount; i++ {
		// Only store connections which start from this tile.
		if offMeshConClass[i*2+0] == 0xff {
			p := navPolys[offMeshPolyBase+n]
			p.VertCount = 2
			p.Verts[0] = uint16(offMeshVertsBase + n*2 + 0)
			p.Verts[1] = uint16(offMeshVertsBase + n*2 + 1)
			p.Flags = params.OffMeshConFlags[i]
			p.SetArea(params.OffMeshConAreas[i])
			p.SetType(DT_POLYTYPE_OFFMESH_CONNECTION)
			n++
		}
	}

	// Store detail meshes and vertices.
	// The nav polygon vertices are stored as the first vertices on each mesh.
	// We compress the mesh data by skipping them and using the navmesh coordinates.
	if len(params.DetailMeshes) > 0 {
		vbase := uint8(0)
		for i := int32(0); i < params.PolyCount; i++ {
			dtl := navDMeshes[i]
			vb := int32(params.DetailMeshes[i*4+0])
			ndv := int32(params.DetailMeshes[i*4+1])
			nv := int32(navPolys[i].VertCount)
			dtl.VertBase = uint32(vbase)
			dtl.VertCount = uint8(ndv - nv)
			dtl.TriBase = params.DetailMeshes[i*4+2]
			dtl.TriCount = uint8(params.DetailMeshes[i*4+3])
			// Copy vertices except the first 'nv' verts which are equal to nav poly verts.
			if ndv-nv > 0 {
				copy(navDVerts[vbase*3:], params.DetailVerts[(vb+nv)*3:(vb+nv)*3+3*(ndv-nv)])
				vbase += uint8(ndv - nv)
			}
		}
		// Store triangles.

		copy(navDTris, params.DetailTris[:4*params.DetailTriCount])
	} else {
		// Create dummy detail mesh by triangulating polys.
		tbase := int32(0)
		for i := int32(0); i < params.PolyCount; i++ {
			dtl := navDMeshes[i]
			nv := int32(navPolys[i].VertCount)
			dtl.VertBase = 0
			dtl.VertCount = 0
			dtl.TriBase = uint32(tbase)
			dtl.TriCount = uint8(nv - 2)
			// Triangulate polygon (local indices).
			for j := int32(2); j < nv; j++ {
				t := navDTris[tbase*4 : tbase*4+4]
				t[0] = 0
				t[1] = uint8(j - 1)
				t[2] = uint8(j)
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
	if params.BuildBvTree {
		createBVTree(params, navBvtree)
	}

	// Store Off-Mesh connections.
	n = 0
	for i := int32(0); i < params.OffMeshConCount; i++ {
		// Only store connections which start from this tile.
		if offMeshConClass[i*2+0] == 0xff {
			con := offMeshCons[n]
			con.Poly = uint16(offMeshPolyBase + n)
			// Copy connection end-points.
			endPts := params.OffMeshConVerts[i*2*3 : i*2*3]
			copy(con.Pos[:3], endPts[:3])
			copy(con.Pos[3:], endPts[3:])
			con.Rad = params.OffMeshConRad[i]
			if params.OffMeshConDir[i] > 0 {
				con.Flags = DT_OFFMESH_CON_BIDIR
			} else {
				con.Flags = 0
			}

			con.Side = offMeshConClass[i*2+1]
			if len(params.OffMeshConUserID) > 0 {
				con.UserId = params.OffMeshConUserID[i]
			}

			n++
		}
	}

	return data, true

}
