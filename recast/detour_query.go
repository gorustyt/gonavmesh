package recast

import (
	"math"
	"math/rand"
)

const (
	DT_FINDPATH_ANY_ANGLE = 0x02 ///< use raycasts during pathfind to "shortcut" (raycast still consider costs)
)
const (
	H_SCALE = 0.999 // Search heuristic scale.

	/// Vertex flags returned by dtNavMeshQuery::findStraightPath.
	DT_STRAIGHTPATH_START              = 0x01 ///< The vertex is the start position in the path.
	DT_STRAIGHTPATH_END                = 0x02 ///< The vertex is the end position in the path.
	DT_STRAIGHTPATH_OFFMESH_CONNECTION = 0x04 ///< The vertex is the start of an off-mesh connection.

)
const (
	/// Options for dtNavMeshQuery::findStraightPath.
	DT_STRAIGHTPATH_AREA_CROSSINGS = 0x01 ///< Add a vertex at every polygon edge crossing where area changes.
	DT_STRAIGHTPATH_ALL_CROSSINGS  = 0x02 ///< Add a vertex at every polygon edge crossing.
)

var (
	DT_VIRTUAL_QUERYFILTER = 1
)

type dtQueryData struct {
	status           dtStatus
	lastBestNode     *dtNode
	lastBestNodeCost float64
	startRef, endRef dtPolyRef
	startPos         [3]float64
	endPos           [3]float64
	filter           *dtQueryFilter
	options          int
	raycastLimitSqr  float64
}

type dtNavMeshQuery struct {
	m_nav          *dtNavMesh ///< Pointer to navmesh data.
	m_nodePool     map[dtPolyRef]*dtNode
	m_tinyNodePool map[dtPolyRef]*dtNode
	m_openList     NodeQueue[*dtNode]
	m_query        dtQueryData
}

func NewDtNavMeshQuery() *dtNavMeshQuery {

	return &dtNavMeshQuery{m_nodePool: map[dtPolyRef]*dtNode{},
		m_tinyNodePool: map[dtPolyRef]*dtNode{},
		m_openList: NewNodeQueue(func(t1, t2 *dtNode) bool {
			return t1.total < t2.total //花费最小
		})}
}

// / Returns random location on navmesh.
// / Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
// /  @param[in]		filter			The polygon filter to apply to the query.
// /  @param[in]		frand			Function returning a random number [0..1).
// /  @param[out]	randomRef		The reference id of the random location.
// /  @param[out]	randomPt		The random location.
// / @returns The status flags for the query.
func (query *dtNavMeshQuery) findRandomPoint(filter *dtQueryFilter) (randomRef dtPolyRef, randomPt []float64, status dtStatus) {
	// Randomly pick one tile. Assume that all tiles cover roughly the same area.
	tsum := 0.0
	var tile *dtMeshTile
	for i := 0; i < query.m_nav.getMaxTiles(); i++ {
		t := query.m_nav.getTile(i)
		if t == nil || t.header == nil {
			continue
		}

		// Choose random tile using reservoir sampling.
		area := 1.0 // Could be tile area too.
		tsum += area
		u := rand.Float64()
		if u*tsum <= area {
			tile = t
		}

	}
	if tile == nil {
		return randomRef, randomPt, DT_FAILURE
	}

	var poly *dtPoly
	var polyRef dtPolyRef
	// Randomly pick one polygon weighted by polygon area.
	base := query.m_nav.getPolyRefBase(tile)

	areaSum := 0.0
	for i := 0; i < tile.header.polyCount; i++ {
		p := tile.polys[i]
		// Do not return off-mesh connection polygons.
		if p.getType() != DT_POLYTYPE_GROUND {
			continue
		}

		// Must pass filter
		ref := base | dtPolyRef(i)
		if !filter.passFilter(p) {
			continue
		}

		// Calc area of the polygon.
		polyArea := 0.0
		for j := 2; j < p.vertCount; j++ {
			va := rcGetVert(tile.verts, p.verts[0])
			vb := rcGetVert(tile.verts, p.verts[j-1])
			vc := rcGetVert(tile.verts, p.verts[j])
			polyArea += dtTriArea2D(va, vb, vc)
		}

		// Choose random polygon weighted by area, using reservoir sampling.
		areaSum += polyArea
		u := rand.Float64()
		if u*areaSum <= polyArea {
			poly = p
			polyRef = ref
		}
	}

	if poly == nil {
		return randomRef, randomPt, DT_FAILURE
	}

	// Randomly pick point on polygon.
	verts := make([]float64, 3*DT_VERTS_PER_POLYGON)
	areas := make([]float64, DT_VERTS_PER_POLYGON)
	copy(verts[0:3], rcGetVert(tile.verts, poly.verts[0]))
	for j := 1; j < poly.vertCount; j++ {
		copy(verts[j:3], rcGetVert(tile.verts, poly.verts[j]))
	}

	s := rand.Float64()
	t := rand.Float64()

	pt := dtRandomPointInConvexPoly(verts, poly.vertCount, areas, s, t)

	randomPt, _, _ = query.closestPointOnPoly(polyRef, pt)
	randomRef = polyRef

	return randomRef, randomPt, DT_SUCCESS
}

// ////////////////////////////////////////////////////////////////////////////////////////
// 查找目标点与指定polygon最近的点
// / @par
// /
// / Uses the detail polygons to find the surface height. (Most accurate.)
// /
// / @p pos does not have to be within the bounds of the polygon or navigation mesh.
// /
// / See closestPointOnPolyBoundary() for a limited but faster option.
// /
func (query *dtNavMeshQuery) closestPointOnPoly(ref dtPolyRef, pos []float64) (closest []float64, posOverPoly bool, status dtStatus) {
	if query.m_nav == nil {
		panic("")
	}
	if query.m_nav.isValidPolyRef(ref) || len(pos) == 0 || !dtVisfinite(pos) {
		return closest, posOverPoly, DT_FAILURE | DT_INVALID_PARAM
	}

	closest, posOverPoly = query.m_nav.closestPointOnPoly(ref, pos)
	return closest, posOverPoly, DT_SUCCESS
}

// / @par
// /
// / Much faster than closestPointOnPoly().
// /
// / If the provided position lies within the polygon's xz-bounds (above or below),
// / then @p pos and @p closest will be equal.
// /
// / The height of @p closest will be the polygon boundary.  The height detail is not used.
// /
// / @p pos does not have to be within the bounds of the polybon or the navigation mesh.
// /
func (query *dtNavMeshQuery) closestPointOnPolyBoundary(ref dtPolyRef, pos []float64) (closest []float64, status dtStatus) {
	if query.m_nav == nil {
		panic("")
	}
	closest = make([]float64, 3)
	tile, poly, status := query.m_nav.getTileAndPolyByRef(ref)
	if status.dtStatusFailed() {
		return closest, DT_FAILURE | DT_INVALID_PARAM
	}

	if len(pos) == 0 || !dtVisfinite(pos) {
		return closest, DT_FAILURE | DT_INVALID_PARAM
	}

	// Collect vertices.
	var verts [DT_VERTS_PER_POLYGON * 3]float64
	var edged [DT_VERTS_PER_POLYGON]float64
	var edget [DT_VERTS_PER_POLYGON]float64
	nv := 0
	for i := 0; i < poly.vertCount; i++ {
		copy(verts[nv*3:nv*3+3], rcGetVert(tile.verts, poly.verts[i]))
		nv++
	}

	inside := dtDistancePtPolyEdgesSqr(pos, verts[:], nv, edged[:], edget[:])
	if inside {
		// Point is inside the polygon, return the point.
		copy(closest, pos)
	} else {
		// Point is outside the polygon, dtClamp to nearest edge.
		dmin := edged[0]
		imin := 0
		for i := 1; i < nv; i++ {
			if edged[i] < dmin {
				dmin = edged[i]
				imin = i
			}
		}
		va := rcGetVert(verts[:], imin)
		vb := rcGetVert(verts[:], (imin+1)%nv)
		closest = dtVlerp(va, vb, edget[imin])
	}

	return closest, DT_SUCCESS
}

type dtQueryFilter struct {
	m_areaCost     [DT_MAX_AREAS]float64 ///< Cost per area type. (Used by default implementation.)
	m_includeFlags int                   ///< Flags for polygons that can be visited. (Used by default implementation.)
	m_excludeFlags int                   ///< Flags for polygons that should not be visited. (Used by default implementation.)
}

func (filter *dtQueryFilter) getCost(pa, pb []float64, curPoly *dtPoly) float64 {
	return dtVdist(pa, pb) * filter.m_areaCost[curPoly.getArea()]
}
func (filter *dtQueryFilter) passFilter(poly *dtPoly) bool {
	return (poly.flags&filter.m_includeFlags) != 0 && (poly.flags&filter.m_excludeFlags) == 0
}

func (query *dtNavMeshQuery) findRandomPointAroundCircle(startRef dtPolyRef, centerPos []float64, maxRadius float64,
	filter *dtQueryFilter) (randomRef dtPolyRef, randomPt []float64, status dtStatus) {
	if query.m_nav == nil {
		return
	}
	query.m_openList.Reset()
	query.m_nodePool = map[dtPolyRef]*dtNode{}
	// Validate input
	if !query.m_nav.isValidPolyRef(startRef) || len(centerPos) == 0 || !dtVisfinite(centerPos) || maxRadius < 0 || !dtIsFinite(maxRadius) || filter == nil {
		return randomRef, randomPt, DT_FAILURE | DT_INVALID_PARAM
	}

	_, startPoly := query.m_nav.getTileAndPolyByRefUnsafe(startRef)
	if !filter.passFilter(startPoly) {
		return randomRef, randomPt, DT_FAILURE | DT_INVALID_PARAM
	}
	var startNode dtNode
	copy(startNode.pos[:], centerPos)
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = 0
	startNode.id = startRef
	startNode.flags = DT_NODE_OPEN
	query.m_openList.Offer(&startNode)
	query.m_nodePool[startRef] = &startNode
	radiusSqr := dtSqr(maxRadius)
	areaSum := 0.0

	var randomTile *dtMeshTile
	var randomPoly *dtPoly
	var randomPolyRef dtPolyRef

	for !query.m_openList.Empty() {
		bestNode := query.m_openList.Poll()
		bestNode.flags &= ^DT_NODE_OPEN
		bestNode.flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.id
		bestTile, bestPoly := query.m_nav.getTileAndPolyByRefUnsafe(bestRef)

		// Place random locations on on ground.
		if bestPoly.getType() == DT_POLYTYPE_GROUND {
			// Calc area of the polygon.
			polyArea := 0.0
			for j := 2; j < bestPoly.vertCount; j++ {
				va := rcGetVert(bestTile.verts, bestPoly.verts[0])
				vb := rcGetVert(bestTile.verts, bestPoly.verts[j-1])
				vc := rcGetVert(bestTile.verts, bestPoly.verts[j])
				polyArea += dtTriArea2D(va, vb, vc)
			}
			// Choose random polygon weighted by area, using reservoir sampling.
			areaSum += polyArea
			u := rand.Float64()
			if u*areaSum <= polyArea {
				randomTile = bestTile
				randomPoly = bestPoly
				randomPolyRef = bestRef
			}
		}

		// Get parent poly and tile.
		var parentRef dtPolyRef
		if bestNode.pidx != nil {
			parentRef = bestNode.pidx.id
		}

		if parentRef > 0 {
			_, _ = query.m_nav.getTileAndPolyByRefUnsafe(parentRef)
		}

		for i := bestPoly.firstLink; i != DT_NULL_LINK; i = bestTile.links[i].next {
			link := bestTile.links[i]
			neighbourRef := link.ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef == 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour
			neighbourTile, neighbourPoly := query.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.

			va, vb, status := query.getPortalPoints1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile)
			if status.dtStatusFailed() {
				continue
			}

			// If the circle is not touching the next polygon, skip it.

			_, distSqr := dtDistancePtSegSqr2D(centerPos, va, vb)
			if distSqr > radiusSqr {
				continue
			}

			neighbourNode := query.m_nodePool[neighbourRef]
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.flags == 0 {
				copy(neighbourNode.pos[:], dtVlerp(va, vb, 0.5))
			}

			total := bestNode.total + dtVdist(bestNode.pos[:], neighbourNode.pos[:])

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_OPEN > 0) && total >= neighbourNode.total {
				continue
			}

			neighbourNode.id = neighbourRef
			neighbourNode.flags = neighbourNode.flags & ^DT_NODE_CLOSED
			neighbourNode.pidx = bestNode
			neighbourNode.total = total

			if neighbourNode.flags&DT_NODE_OPEN > 0 {
				query.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.flags = DT_NODE_OPEN
				query.m_openList.Offer(neighbourNode)
			}
		}
	}

	if randomPoly == nil {
		return randomRef, randomPt, DT_FAILURE
	}

	// Randomly pick point on polygon.
	var verts [3 * DT_VERTS_PER_POLYGON]float64
	var areas [DT_VERTS_PER_POLYGON]float64
	copy(verts[0:3], rcGetVert(randomTile.verts, randomPoly.verts[0]))
	for j := 1; j < randomPoly.vertCount; j++ {
		copy(verts[j*3:j*3+3], rcGetVert(randomTile.verts, randomPoly.verts[j]))
	}

	s := rand.Float64()
	t := rand.Float64()

	pt := dtRandomPointInConvexPoly(verts[:], randomPoly.vertCount, areas[:], s, t)

	randomPt, _, _ = query.closestPointOnPoly(randomPolyRef, pt)

	randomRef = randomPolyRef

	return randomRef, randomPt, status
}

func (query *dtNavMeshQuery) getPortalPoints(from dtPolyRef, to dtPolyRef) (left, right []float64, fromType int, toType int, status dtStatus) {
	left = make([]float64, 3)
	right = make([]float64, 3)
	if query.m_nav == nil {
		panic("")
	}

	fromTile, fromPoly, status := query.m_nav.getTileAndPolyByRef(from)
	if status.dtStatusFailed() {
		return left, right, fromType, toType, DT_FAILURE | DT_INVALID_PARAM
	}

	fromType = fromPoly.getType()

	toTile, toPoly, status := query.m_nav.getTileAndPolyByRef(to)
	if status.dtStatusFailed() {
		return left, right, fromType, toType, DT_FAILURE | DT_INVALID_PARAM
	}

	toType = toPoly.getType()
	left, right, status = query.getPortalPoints1(from, fromPoly, fromTile, to, toPoly, toTile)
	return left, right, fromType, toType, status
}

// Returns portal points between two polygons.
func (query *dtNavMeshQuery) getPortalPoints1(from dtPolyRef, fromPoly *dtPoly, fromTile *dtMeshTile,
	to dtPolyRef, toPoly *dtPoly, toTile *dtMeshTile) (left, right []float64, status dtStatus) {
	// Find the link that points to the 'to' polygon.
	// Find the link that points to the 'to' polygon.
	left = make([]float64, 3)
	right = make([]float64, 3)
	var link *dtLink
	for i := fromPoly.firstLink; i != DT_NULL_LINK; i = fromTile.links[i].next {
		if fromTile.links[i].ref == to {
			link = fromTile.links[i]
			break
		}
	}
	if link == nil {
		return left, right, DT_FAILURE | DT_INVALID_PARAM
	}

	// Handle off-mesh connections.
	if fromPoly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		// Find link that points to first vertex.
		for i := fromPoly.firstLink; i != DT_NULL_LINK; i = fromTile.links[i].next {
			if fromTile.links[i].ref == to {
				v := fromTile.links[i].edge
				copy(left, rcGetVert(fromTile.verts, fromPoly.verts[v]))
				copy(right, rcGetVert(fromTile.verts, fromPoly.verts[v]))
				return left, right, DT_SUCCESS
			}
		}
		return left, right, DT_FAILURE | DT_INVALID_PARAM
	}

	if toPoly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		for i := toPoly.firstLink; i != DT_NULL_LINK; i = toTile.links[i].next {
			if toTile.links[i].ref == from {
				v := toTile.links[i].edge
				copy(left, rcGetVert(toTile.verts, toPoly.verts[v]))
				copy(right, rcGetVert(toTile.verts, toPoly.verts[v]))
				return left, right, DT_SUCCESS
			}
		}
		return left, right, DT_FAILURE | DT_INVALID_PARAM
	}

	// Find portal vertices.
	v0 := fromPoly.verts[link.edge]
	v1 := fromPoly.verts[(link.edge+1)%fromPoly.vertCount]
	copy(left, rcGetVert(fromTile.verts, v0))
	copy(right, rcGetVert(fromTile.verts, v1))
	// If the link is at tile boundary, dtClamp the vertices to
	// the link width.
	if link.side != 0xff {
		// Unpack portal limits.
		if link.bmin != 0 || link.bmax != 255 {
			s := 1.0 / 255.0
			tmin := float64(link.bmin) * s
			tmax := float64(link.bmax) * s
			left = dtVlerp(rcGetVert(fromTile.verts, v0), rcGetVert(fromTile.verts, v1), tmin)
			right = dtVlerp(rcGetVert(fromTile.verts, v0), rcGetVert(fromTile.verts, v1), tmax)
		}
	}

	return left, right, DT_SUCCESS
}

func (q *dtNavMeshQuery) queryPolygonsInTile(tile *dtMeshTile, qmin, qmax []float64, filter *dtQueryFilter, query dtPolyQuery) {
	if q.m_nav == nil {
		panic("")
	}
	batchSize := 32
	var polyRefs [32]dtPolyRef
	var polys [32]*dtPoly
	n := 0

	if tile.bvTree != nil {
		node := 0
		end := tile.header.bvNodeCount
		tbmin := tile.header.bmin
		tbmax := tile.header.bmax
		qfac := tile.header.bvQuantFactor

		// Calculate quantized box
		var bmin [3]int
		var bmax [3]int
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
		base := q.m_nav.getPolyRefBase(tile)
		for node < end {

			overlap := dtOverlapQuantBounds(bmin, bmax, tile.bvTree[node].bmin, tile.bvTree[node].bmax)
			isLeafNode := tile.bvTree[node].i >= 0

			if isLeafNode && overlap {
				ref := base | dtPolyRef(tile.bvTree[node].i)
				if filter.passFilter(tile.polys[tile.bvTree[node].i]) {
					polyRefs[n] = ref
					polys[n] = tile.polys[tile.bvTree[node].i]

					if n == batchSize-1 {
						query.process(tile, polyRefs[:], batchSize)
						n = 0
					} else {
						n++
					}
				}
			}

			if overlap || isLeafNode {
				node++
			} else {
				escapeIndex := -tile.bvTree[node].i
				node += escapeIndex
			}
		}
	} else {
		bmin := make([]float64, 3)
		bmax := make([]float64, 3)
		base := q.m_nav.getPolyRefBase(tile)
		for i := 0; i < tile.header.polyCount; i++ {
			p := tile.polys[i]
			// Do not return off-mesh connection polygons.
			if p.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Must pass filter
			ref := base | dtPolyRef(i)
			if filter != nil && filter.passFilter(p) {
				continue
			}

			// Calc polygon bounds.
			v := rcGetVert(tile.verts, p.verts[0])
			copy(bmin, v)
			copy(bmax, v)
			for j := 1; j < p.vertCount; j++ {
				v := rcGetVert(tile.verts, p.verts[j])
				bmin = dtVmin(bmin, v)
				bmax = dtVmax(bmax, v)
			}
			if dtOverlapBounds(qmin, qmax, bmin, bmax) {
				polyRefs[n] = ref
				polys[n] = p

				if n == batchSize-1 {
					query.process(tile, polyRefs[:], batchSize)
					n = 0
				} else {
					n++
				}
			}
		}
	}

	// Process the last polygons that didn't make a full batch.
	if n > 0 {
		query.process(tile, polyRefs[:], n)
	}

}

// / Provides custom polygon query behavior.
// / Used by dtNavMeshQuery::queryPolygons.
// / @ingroup detour
type dtPolyQuery interface {
	/// Called for each batch of unique polygons touched by the search area in dtNavMeshQuery::queryPolygons.
	/// This can be called multiple times for a single query.
	process(tile *dtMeshTile, refs []dtPolyRef, count int)
}

type dtFindNearestPolyQuery struct {
	m_query              *dtNavMeshQuery
	m_center             []float64
	m_nearestDistanceSqr float64
	m_nearestRef         dtPolyRef
	m_nearestPoint       [3]float64
	m_overPoly           bool
}

func (query *dtFindNearestPolyQuery) process(tile *dtMeshTile, refs []dtPolyRef, count int) {

	for i := 0; i < count; i++ {
		ref := refs[i]
		var d float64
		closestPtPoly, posOverPoly, _ := query.m_query.closestPointOnPoly(ref, query.m_center)

		// If a point is directly over a polygon and closer than
		// climb height, favor that instead of straight line nearest point.
		diff := dtVsub(query.m_center, closestPtPoly)
		if posOverPoly {
			d = dtAbs(diff[1]) - tile.header.walkableClimb
			if d > 0 {
				d = d * d
			} else {
				d = 0
			}
		} else {
			d = dtVlenSqr(diff)
		}

		if d < query.m_nearestDistanceSqr {
			copy(query.m_nearestPoint[:], closestPtPoly)
			query.m_nearestDistanceSqr = d
			query.m_nearestRef = ref
			query.m_overPoly = posOverPoly
		}
	}
}

type dtCollectPolysQuery struct {
	m_polys        []dtPolyRef
	m_maxPolys     int
	m_numCollected int
	m_overflow     bool
}

func (query *dtCollectPolysQuery) process(refs []dtPolyRef, count int) {
	numLeft := query.m_maxPolys - query.m_numCollected
	toCopy := count
	if toCopy > numLeft {
		query.m_overflow = true
		toCopy = numLeft
	}
	query.m_polys = append(query.m_polys, refs[:toCopy]...)
	query.m_numCollected += toCopy
}

type dtSegInterval struct {
	ref        dtPolyRef
	tmin, tmax int
}

func insertInterval(ints []*dtSegInterval, maxInts int, tmin, tmax int, ref dtPolyRef) (nints int) {
	if nints+1 > maxInts {
		return
	}
	// Find insertion point.
	idx := 0
	for idx < nints {
		if tmax <= ints[idx].tmin {
			break
		}
		idx++
	}
	// Move current results.
	if nints-idx > 0 {
		ints = append(ints)
		copy(ints[idx+1:], ints[idx:nints-idx]) //TODO
	}
	// Store
	ints[idx].ref = ref
	ints[idx].tmin = tmin
	ints[idx].tmax = tmax
	nints++
	return
}

// / @par
// /
// / @p hitPos is not adjusted using the height detail data.
// /
// / @p hitDist will equal the search radius if there is no wall within the
// / radius. In this case the values of @p hitPos and @p hitNormal are
// / undefined.
// /
// / The normal will become unpredicable if @p hitDist is a very small number.
// /
func (q *dtNavMeshQuery) findDistanceToWall(startRef dtPolyRef, centerPos []float64, maxRadius float64,
	filter *dtQueryFilter, hitPos []float64, hitNormal []float64) (hitDist float64, status dtStatus) {

	// Validate input
	if !q.m_nav.isValidPolyRef(startRef) ||
		len(centerPos) == 0 || !dtVisfinite(centerPos) || maxRadius < 0 || !dtIsFinite(maxRadius) || filter == nil {
		return hitDist, DT_FAILURE | DT_INVALID_PARAM
	}

	q.m_nodePool = map[dtPolyRef]*dtNode{}
	q.m_openList.Reset()

	startNode := q.m_nodePool[startRef]
	copy(startNode.pos[:], centerPos)
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = 0
	startNode.id = startRef
	startNode.flags = DT_NODE_OPEN
	q.m_openList.Offer(startNode)

	radiusSqr := dtSqr(maxRadius)

	for !q.m_openList.Empty() {
		bestNode := q.m_openList.Poll()
		bestNode.flags &= ^DT_NODE_OPEN
		bestNode.flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.id
		var bestTile *dtMeshTile
		var bestPoly *dtPoly
		bestTile, bestPoly = q.m_nav.getTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef dtPolyRef

		if bestNode.pidx != nil {
			parentRef = bestNode.pidx.id
		}

		if parentRef != 0 {
			q.m_nav.getTileAndPolyByRefUnsafe(parentRef)
		}

		i := 0
		j := bestPoly.vertCount - 1
		// Hit test walls.
		for i < bestPoly.vertCount {
			// Skip non-solid edges.
			if bestPoly.neis[j]&DT_EXT_LINK > 0 {
				// Tile border.
				solid := true
				for k := bestPoly.firstLink; k != DT_NULL_LINK; k = bestTile.links[k].next {
					link := bestTile.links[k]
					if link.edge == j {
						if link.ref != 0 {
							_, neiPoly := q.m_nav.getTileAndPolyByRefUnsafe(link.ref)
							if filter.passFilter(neiPoly) {
								solid = false
							}

						}
						break
					}
				}
				if !solid {
					continue
				}
			} else if bestPoly.neis[j] > 0 {
				// Internal edge
				idx := (bestPoly.neis[j] - 1)
				if filter.passFilter(bestTile.polys[idx]) {
					continue
				}

			}

			// Calc distance to the edge.
			vj := rcGetVert(bestTile.verts, bestPoly.verts[j])
			vi := rcGetVert(bestTile.verts, bestPoly.verts[j])

			tseg, distSqr := dtDistancePtSegSqr2D(centerPos, vj, vi)

			// Edge is too far, skip.
			if distSqr > radiusSqr {
				continue
			}

			// Hit wall, update radius.
			radiusSqr = distSqr
			// Calculate hit pos.
			hitPos[0] = vj[0] + (vi[0]-vj[0])*tseg
			hitPos[1] = vj[1] + (vi[1]-vj[1])*tseg
			hitPos[2] = vj[2] + (vi[2]-vj[2])*tseg
			j = i
			i++
		}

		for i := bestPoly.firstLink; i != DT_NULL_LINK; i = bestTile.links[i].next {
			link := bestTile.links[i]
			neighbourRef := link.ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour.

			neighbourTile, neighbourPoly := q.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			// Skip off-mesh connections.
			if neighbourPoly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Calc distance to the edge.
			va := rcGetVert(bestTile.verts, bestPoly.verts[link.edge])
			vb := rcGetVert(bestTile.verts, (link.edge+1)%bestPoly.vertCount)

			_, distSqr := dtDistancePtSegSqr2D(centerPos, va, vb)

			// If the circle is not touching the next polygon, skip it.
			if distSqr > radiusSqr {
				continue
			}

			if !filter.passFilter(neighbourPoly) {
				continue
			}

			neighbourNode := q.m_nodePool[neighbourRef]
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.flags == 0 {
				pos, _ := q.getEdgeMidPoint1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile)
				copy(neighbourNode.pos[:], pos)
			}

			total := bestNode.total + dtVdist(bestNode.pos[:], neighbourNode.pos[:])

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_OPEN > 0) && total >= neighbourNode.total {
				continue
			}

			neighbourNode.id = neighbourRef
			neighbourNode.flags = (neighbourNode.flags & ^DT_NODE_CLOSED)
			neighbourNode.pidx = bestNode
			neighbourNode.total = total

			if neighbourNode.flags&DT_NODE_OPEN > 0 {
				q.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.flags |= DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}
		}
	}

	// Calc hit normal.
	hitNormal = dtVsub(centerPos, hitPos)
	dtVnormalize(hitNormal)

	hitDist = math.Sqrt(radiusSqr)

	return hitDist, status
}

// / @par
// /
// / The closed list is the list of polygons that were fully evaluated during
// / the last navigation graph search. (A* or Dijkstra)
// /
func (q *dtNavMeshQuery) isInClosedList(ref dtPolyRef) bool {
	if q.m_nodePool == nil {
		return false
	}
	n := q.m_nodePool[ref]
	if n.flags&DT_NODE_CLOSED > 0 {
		return true
	}
	return false
}

// Returns edge mid point between two polygons.
func (q *dtNavMeshQuery) getEdgeMidPoint(from, to dtPolyRef) (mid []float64, status dtStatus) {
	mid = make([]float64, 3)
	left, right, _, _, status := q.getPortalPoints(from, to)
	if status.dtStatusFailed() {
		return mid, DT_FAILURE | DT_INVALID_PARAM
	}

	mid[0] = (left[0] + right[0]) * 0.5
	mid[1] = (left[1] + right[1]) * 0.5
	mid[2] = (left[2] + right[2]) * 0.5
	return mid, DT_SUCCESS
}

func (q *dtNavMeshQuery) getEdgeMidPoint1(from dtPolyRef, fromPoly *dtPoly, fromTile *dtMeshTile,
	to dtPolyRef, toPoly *dtPoly, toTile *dtMeshTile) (mid []float64, status dtStatus) {
	mid = make([]float64, 3)

	left, right, status := q.getPortalPoints1(from, fromPoly, fromTile, to, toPoly, toTile)
	if status.dtStatusFailed() {
		return mid, DT_FAILURE | DT_INVALID_PARAM
	}

	mid[0] = (left[0] + right[0]) * 0.5
	mid[1] = (left[1] + right[1]) * 0.5
	mid[2] = (left[2] + right[2]) * 0.5
	return mid, DT_SUCCESS
}

func (q *dtNavMeshQuery) appendVertex(pos []float64, flags int, ref dtPolyRef,
	straightPath []float64, straightPathFlags []int, straightPathRefs []dtPolyRef, straightPathCount, maxStraightPath int) (int, dtStatus) {
	if (straightPathCount) > 0 && dtVequal(straightPath[((straightPathCount)-1)*3:((straightPathCount)-1)*3+3], pos) {
		// The vertices are equal, update flags and poly.
		if len(straightPathFlags) > 0 {
			straightPathFlags[(straightPathCount)-1] = flags
		}

		if len(straightPathFlags) > 0 {
			straightPathRefs[(straightPathCount)-1] = ref
		}

	} else {
		// Append new vertex.
		copy(straightPath[(straightPathCount)*3:(straightPathCount)*3+3], pos)
		if len(straightPathFlags) > 0 {
			straightPathFlags[(straightPathCount)] = flags
		}

		if len(straightPathRefs) > 0 {
			straightPathRefs[(straightPathCount)] = ref
		}

		(straightPathCount)++

		// If there is no space to append more vertices, return.
		if (straightPathCount) >= maxStraightPath {
			return straightPathCount, DT_SUCCESS | DT_BUFFER_TOO_SMALL
		}

		// If reached end of path, return.
		if flags == DT_STRAIGHTPATH_END {
			return straightPathCount, DT_SUCCESS
		}
	}
	return straightPathCount, DT_IN_PROGRESS
}

// / @par
// /
// / This method peforms what is often called 'string pulling'.
// /
// / The start position is clamped to the first polygon in the path, and the
// / end position is clamped to the last. So the start and end positions should
// / normally be within or very near the first and last polygons respectively.
// /
// / The returned polygon references represent the reference id of the polygon
// / that is entered at the associated path position. The reference id associated
// / with the end point will always be zero.  This allows, for example, matching
// / off-mesh link points to their representative polygons.
// /
// / If the provided result buffers are too small for the entire result set,
// / they will be filled as far as possible from the start toward the end
// / position.
// /
func (q *dtNavMeshQuery) findStraightPath(startPos, endPos []float64,
	path []dtPolyRef, pathSize int,
	straightPath []float64, straightPathFlags []int, straightPathRefs []dtPolyRef,
	maxStraightPath int, options int) (straightPathCount int, status dtStatus) {
	if q.m_nav == nil {
		panic("")
	}
	//if (!startPos || !dtVisfinite(startPos) ||//TODO
	//!endPos || !dtVisfinite(endPos) ||
	//!path || pathSize <= 0 || !path[0] ||
	//maxStraightPath <= 0)
	//{
	//return DT_FAILURE | DT_INVALID_PARAM;
	//}

	// TODO: Should this be callers responsibility?

	closestStartPos, _ := q.closestPointOnPolyBoundary(path[0], startPos)
	if status.dtStatusFailed() {
		return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	closestEndPos, status := q.closestPointOnPolyBoundary(path[pathSize-1], endPos)
	if status.dtStatusFailed() {
		return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	// Add start point.
	straightPathCount, status = q.appendVertex(closestStartPos, DT_STRAIGHTPATH_START, path[0],
		straightPath, straightPathFlags, straightPathRefs,
		straightPathCount, maxStraightPath)
	if status != DT_IN_PROGRESS {
		return straightPathCount, status
	}

	if pathSize > 1 {
		var portalApex, portalLeft, portalRight [3]float64
		copy(portalApex[:], closestStartPos)
		copy(portalLeft[:], portalApex[:])
		copy(portalRight[:], portalApex[:])
		var apexIndex = 0
		var leftIndex = 0
		var rightIndex = 0

		var leftPolyType = 0
		var rightPolyType = 0

		var leftPolyRef dtPolyRef = path[0]
		var rightPolyRef dtPolyRef = path[0]

		for i := 0; i < pathSize; i++ {
			var left, right []float64
			var toType int
			if i+1 < pathSize {

				left, right, _, toType, status = q.getPortalPoints(path[i], path[i+1])
				// Next portal.
				if status.dtStatusFailed() {
					// Failed to get portal points, in practice this means that path[i+1] is invalid polygon.
					// Clamp the end point to path[i], and return the path so far.
					closestEndPos, status = q.closestPointOnPolyBoundary(path[i], endPos)
					if status.dtStatusFailed() {
						// This should only happen when the first polygon is invalid.
						return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
					}

					// Apeend portals along the current straight path segment.
					if options&(DT_STRAIGHTPATH_AREA_CROSSINGS|DT_STRAIGHTPATH_ALL_CROSSINGS) > 0 {
						// Ignore status return value as we're just about to return anyway.
						q.appendPortals(apexIndex, i, closestEndPos, path,
							straightPath, straightPathFlags, straightPathRefs,
							straightPathCount, maxStraightPath, options)
					}

					// Ignore status return value as we're just about to return anyway.
					straightPathCount, _ = q.appendVertex(closestEndPos, 0, path[i],
						straightPath, straightPathFlags, straightPathRefs,
						straightPathCount, maxStraightPath)
					tmp := 0
					if straightPathCount >= maxStraightPath {
						tmp = DT_BUFFER_TOO_SMALL
					}
					return straightPathCount, DT_SUCCESS | DT_PARTIAL_RESULT | dtStatus(tmp)
				}

				// If starting really close the portal, advance.
				if i == 0 {
					t, _ := dtDistancePtSegSqr2D(portalApex[:], left, right)
					if t < dtSqr(0.001) {
						continue
					}

				}
			} else {
				// End of the path.
				copy(left[:], closestEndPos)
				copy(right[:], closestEndPos)

				toType = DT_POLYTYPE_GROUND
			}

			// Right vertex.
			if dtTriArea2D(portalApex[:], portalRight[:], right) <= 0.0 {
				if dtVequal(portalApex[:], portalRight[:]) || dtTriArea2D(portalApex[:], portalLeft[:], right) > 0.0 {
					copy(portalRight[:], right)
					if i+1 < pathSize {
						rightPolyRef = path[i+1]
					} else {
						rightPolyRef = 0
					}
					rightPolyType = toType
					rightIndex = i
				} else {
					// Append portals along the current straight path segment.
					if options&(DT_STRAIGHTPATH_AREA_CROSSINGS|DT_STRAIGHTPATH_ALL_CROSSINGS) > 0 {
						straightPathCount, status = q.appendPortals(apexIndex, leftIndex, portalLeft[:], path,
							straightPath, straightPathFlags, straightPathRefs,
							straightPathCount, maxStraightPath, options)
						if status != DT_IN_PROGRESS {
							return straightPathCount, status
						}

					}

					copy(portalApex[:], portalLeft[:])
					apexIndex = leftIndex

					var flags = 0
					if leftPolyRef == 0 {
						flags = DT_STRAIGHTPATH_END
					} else if leftPolyType == DT_POLYTYPE_OFFMESH_CONNECTION {
						flags = DT_STRAIGHTPATH_OFFMESH_CONNECTION
					}

					ref := leftPolyRef

					// Append or update vertex
					straightPathCount, status = q.appendVertex(portalApex[:], flags, ref,
						straightPath, straightPathFlags, straightPathRefs,
						straightPathCount, maxStraightPath)
					if status != DT_IN_PROGRESS {
						return straightPathCount, status
					}

					copy(portalLeft[:], portalApex[:])
					copy(portalRight[:], portalApex[:])
					leftIndex = apexIndex
					rightIndex = apexIndex

					// Restart
					i = apexIndex

					continue
				}
			}

			// Left vertex.
			if dtTriArea2D(portalApex[:], portalLeft[:], left) >= 0.0 {
				if dtVequal(portalApex[:], portalLeft[:]) || dtTriArea2D(portalApex[:], portalRight[:], left) < 0.0 {
					copy(portalLeft[:], left)
					leftPolyRef = 0
					if i+1 < pathSize {
						leftPolyRef = path[i+1]
					}
					leftPolyType = toType
					leftIndex = i
				} else {
					// Append portals along the current straight path segment.
					if options&(DT_STRAIGHTPATH_AREA_CROSSINGS|DT_STRAIGHTPATH_ALL_CROSSINGS) > 0 {
						straightPathCount, status = q.appendPortals(apexIndex, rightIndex, portalRight[:], path,
							straightPath, straightPathFlags, straightPathRefs,
							straightPathCount, maxStraightPath, options)
						if status != DT_IN_PROGRESS {
							return straightPathCount, status
						}

					}

					copy(portalApex[:], portalRight[:])
					apexIndex = rightIndex

					var flags = 0
					if rightPolyRef == 0 {
						flags = DT_STRAIGHTPATH_END
					} else if rightPolyType == DT_POLYTYPE_OFFMESH_CONNECTION {
						flags = DT_STRAIGHTPATH_OFFMESH_CONNECTION
					}

					ref := rightPolyRef

					// Append or update vertex
					straightPathCount, status = q.appendVertex(portalApex[:], flags, ref,
						straightPath, straightPathFlags, straightPathRefs,
						straightPathCount, maxStraightPath)
					if status != DT_IN_PROGRESS {
						return straightPathCount, status
					}

					copy(portalLeft[:], portalApex[:])
					copy(portalRight[:], portalApex[:])
					leftIndex = apexIndex
					rightIndex = apexIndex

					// Restart
					i = apexIndex

					continue
				}
			}
		}

		// Append portals along the current straight path segment.
		if options&(DT_STRAIGHTPATH_AREA_CROSSINGS|DT_STRAIGHTPATH_ALL_CROSSINGS) > 0 {
			straightPathCount, status = q.appendPortals(apexIndex, pathSize-1, closestEndPos, path,
				straightPath, straightPathFlags, straightPathRefs,
				straightPathCount, maxStraightPath, options)
			if status != DT_IN_PROGRESS {
				return straightPathCount, status
			}

		}
	}

	// Ignore status return value as we're just about to return anyway.
	straightPathCount, _ = q.appendVertex(closestEndPos, DT_STRAIGHTPATH_END, 0,
		straightPath, straightPathFlags, straightPathRefs,
		straightPathCount, maxStraightPath)
	tmp := 0
	if straightPathCount >= maxStraightPath {
		tmp = DT_BUFFER_TOO_SMALL
	}
	return straightPathCount, DT_SUCCESS | dtStatus(tmp)
}

func (q *dtNavMeshQuery) appendPortals(startIdx int, endIdx int, endPos []float64, path []dtPolyRef,
	straightPath []float64, straightPathFlags []int, straightPathRefs []dtPolyRef,
	straightPathCount int, maxStraightPath int, options int) (int, dtStatus) {
	startPos := straightPath[(straightPathCount-1)*3 : (straightPathCount-1)*3+3]
	// Append or update last vertex
	for i := startIdx; i < endIdx; i++ {
		// Calculate portal
		from := path[i]
		fromTile, fromPoly, status := q.m_nav.getTileAndPolyByRef(from)
		if status.dtStatusFailed() {
			return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
		}

		to := path[i+1]
		toTile, toPoly, status := q.m_nav.getTileAndPolyByRef(to)
		if status.dtStatusFailed() {
			return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
		}

		left, right, status := q.getPortalPoints1(from, fromPoly, fromTile, to, toPoly, toTile)
		if status.dtStatusFailed() {
			break
		}

		if options&DT_STRAIGHTPATH_AREA_CROSSINGS > 0 {
			// Skip intersection if only area crossings are requested.
			if fromPoly.getArea() == toPoly.getArea() {
				continue
			}

		}

		// Append intersection
		_, t, ok := dtIntersectSegSeg2D(startPos, endPos, left, right)
		if ok {

			pt := dtVlerp(left, right, t)

			straightPathCount, status = q.appendVertex(pt, 0, path[i+1],
				straightPath, straightPathFlags, straightPathRefs,
				straightPathCount, maxStraightPath)
			if status != DT_IN_PROGRESS {
				return straightPathCount, status
			}

		}
	}
	return straightPathCount, DT_IN_PROGRESS
}

// / @par
// /
// / This method is optimized for small delta movement and a small number of
// / polygons. If used for too great a distance, the result set will form an
// / incomplete path.
// /
// / @p resultPos will equal the @p endPos if the end is reached.
// / Otherwise the closest reachable position will be returned.
// /
// / @p resultPos is not projected onto the surface of the navigation
// / mesh. Use #getPolyHeight if this is needed.
// /
// / This method treats the end position in the same manner as
// / the #raycast method. (As a 2D point.) See that method's documentation
// / for details.
// /
// / If the @p visited array is too small to hold the entire result set, it will
// / be filled as far as possible from the start position toward the end
// / position.
// /
func (q *dtNavMeshQuery) moveAlongSurface(startRef dtPolyRef, startPos, endPos []float64,
	filter *dtQueryFilter, maxVisitedSize int) (resultPos []float64, visited []dtPolyRef, visitedCount int, status dtStatus) {
	dtAssert(q.m_nav)
	dtAssert(q.m_tinyNodePool)
	//if (!m_nav->isValidPolyRef(startRef) ||
	// !startPos || !dtVisfinite(startPos) ||
	// !endPos || !dtVisfinite(endPos) ||
	// !filter || !resultPos || !visited ||
	// maxVisitedSize <= 0)
	//{
	// return DT_FAILURE | DT_INVALID_PARAM;
	//}

	status = DT_SUCCESS
	MAX_STACK := 48
	var stack [48]*dtNode
	nstack := 0

	q.m_tinyNodePool = map[dtPolyRef]*dtNode{}
	var startNode dtNode
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = 0
	startNode.id = startRef
	startNode.flags = DT_NODE_CLOSED
	q.m_tinyNodePool[startRef] = &startNode
	stack[nstack] = &startNode
	nstack++

	bestDist := math.MaxFloat64
	bestPos := dtVcopy(startPos)
	var bestNode *dtNode
	// Search constraints
	searchPos := dtVlerp(startPos, endPos, 0.5)
	searchRadSqr := dtSqr(dtVdist(startPos, endPos)/2.0 + 0.001)

	var verts [DT_VERTS_PER_POLYGON * 3]float64

	for nstack > 0 {
		// Pop front.
		curNode := stack[0]
		for i := 0; i < nstack-1; i++ {
			stack[i] = stack[i+1]
		}

		nstack--

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		curRef := curNode.id
		curTile, curPoly := q.m_nav.getTileAndPolyByRefUnsafe(curRef)

		// Collect vertices.
		nverts := curPoly.vertCount
		for i := 0; i < nverts; i++ {
			copy(verts[i*3:i*3+3], rcGetVert(curTile.verts, curPoly.verts[i]))
		}

		// If target is inside the poly, stop search.
		if dtPointInPolygon(endPos, verts[:], nverts) {
			bestNode = curNode
			copy(bestPos, endPos)
			break
		}
		i := 0
		j := curPoly.vertCount - 1
		// Find wall edges and find nearest point inside the walls.
		for i < curPoly.vertCount {
			// Find links to neighbours.
			MAX_NEIS := 8
			nneis := 0
			var neis [8]dtPolyRef

			if curPoly.neis[j]&DT_EXT_LINK > 0 {
				// Tile border.
				for k := curPoly.firstLink; k != DT_NULL_LINK; k = curTile.links[k].next {
					link := curTile.links[k]
					if link.edge == j {
						if link.ref != 0 {
							_, neiPoly := q.m_nav.getTileAndPolyByRefUnsafe(link.ref)
							if filter.passFilter(neiPoly) {
								if nneis < MAX_NEIS {
									neis[nneis] = link.ref
									nneis++
								}

							}
						}
					}

				}
			} else if curPoly.neis[j] > 0 {
				idx := (curPoly.neis[j] - 1)
				ref := q.m_nav.getPolyRefBase(curTile) | dtPolyRef(idx)
				if filter.passFilter(curTile.polys[idx]) {
					// Internal edge, encode id.
					neis[nneis] = ref
					nneis++
				}
			}

			if nneis == 0 {
				// Wall edge, calc distance.
				vj := rcGetVert(verts[:], j)
				vi := rcGetVert(verts[:], i)
				tseg, distSqr := dtDistancePtSegSqr2D(endPos, vj, vi)
				if distSqr < bestDist {
					// Update nearest distance.
					bestPos = dtVlerp(vj, vi, tseg)
					bestDist = distSqr
					bestNode = curNode
				}
			} else {
				for k := 0; k < nneis; k++ {
					// Skip if no node can be allocated.
					neighbourNode := q.m_tinyNodePool[neis[k]]
					if neighbourNode == nil {
						continue
					}

					// Skip if already visited.
					if neighbourNode.flags&DT_NODE_CLOSED > 0 {
						continue
					}

					// Skip the link if it is too far from search constraint.
					// TODO: Maybe should use getPortalPoints(), but this one is way faster.
					vj := rcGetVert(verts[:], j)
					vi := rcGetVert(verts[:], i)
					_, distSqr := dtDistancePtSegSqr2D(searchPos, vj, vi)
					if distSqr > searchRadSqr {
						continue
					}

					// Mark as the node as visited and push to queue.
					if nstack < MAX_STACK {
						neighbourNode.pidx = curNode
						neighbourNode.flags |= DT_NODE_CLOSED
						stack[nstack] = neighbourNode
						nstack++
					}

				}
			}

			j = i
			i++
		}
	}

	n := 0
	if bestNode != nil {
		// Reverse the path.
		var prev *dtNode
		node := bestNode
		for node != nil {
			next := node.pidx
			node.pidx = prev
			prev = node
			node = next
		}
		// Store result
		node = prev
		for node != nil {
			visited[n] = node.id
			n++
			if n >= maxVisitedSize {
				status |= DT_BUFFER_TOO_SMALL
				break
			}
			node = node.pidx
		}
	}

	copy(resultPos, bestPos)

	visitedCount = n

	return resultPos, visited, visitedCount, status
}

/// @par
///
/// At least one result array must be provided.
///
/// The order of the result set is from least to highest cost to reach the polygon.
///
/// A common use case for this method is to perform Dijkstra searches.
/// Candidate polygons are found by searching the graph beginning at the start polygon.
///
/// If a polygon is not found via the graph search, even if it intersects the
/// search circle, it will not be included in the result set. For example:
///
/// polyA is the start polygon.
/// polyB shares an edge with polyA. (Is adjacent.)
/// polyC shares an edge with polyB, but not with polyA
/// Even if the search circle overlaps polyC, it will not be included in the
/// result set unless polyB is also in the set.
///
/// The value of the center point is used as the start position for cost
/// calculations. It is not projected onto the surface of the mesh, so its
/// y-value will effect the costs.
///
/// Intersection tests occur in 2D. All polygons and the search circle are
/// projected onto the xz-plane. So the y-value of the center point does not
/// effect intersection tests.
///
/// If the result arrays are to small to hold the entire result set, they will be
/// filled to capacity.
///

func (q *dtNavMeshQuery) findPolysAroundCircle(startRef dtPolyRef, centerPos []float64, radius float64,
	filter *dtQueryFilter,
	resultCost []float64, resultRef []dtPolyRef, resultParent []dtPolyRef, maxResult int) (resultCount int, status dtStatus) {
	dtAssert(q.m_nav)
	dtAssert(q.m_nodePool)
	dtAssert(q.m_openList)

	//if (!m_nav->isValidPolyRef(startRef) ||
	//!centerPos || !dtVisfinite(centerPos) ||
	//radius < 0 || !dtMathIsfinite(radius) ||
	//!filter || maxResult < 0)
	//{
	//return DT_FAILURE | DT_INVALID_PARAM;
	//}

	q.m_nodePool = map[dtPolyRef]*dtNode{}
	q.m_openList.Reset()
	var startNode dtNode
	copy(startNode.pos[:], centerPos)
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = 0
	startNode.id = startRef
	startNode.flags = DT_NODE_OPEN
	q.m_nodePool[startRef] = &startNode
	q.m_openList.Offer(&startNode)

	status = DT_SUCCESS
	n := 0
	radiusSqr := dtSqr(radius)

	for !q.m_openList.Empty() {
		bestNode := q.m_openList.Poll()
		bestNode.flags &= ^DT_NODE_OPEN
		bestNode.flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.id
		bestTile, bestPoly := q.m_nav.getTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef dtPolyRef
		if bestNode.pidx != nil {
			parentRef = bestNode.pidx.id
		}

		if parentRef > 0 {
			q.m_nav.getTileAndPolyByRefUnsafe(parentRef)
		}

		if n < maxResult {
			if len(resultRef) > 0 {
				resultRef[n] = bestRef
			}

			if len(resultParent) > 0 {
				resultParent[n] = parentRef
			}

			if len(resultCost) > 0 {
				resultCost[n] = bestNode.total
			}

			n++
		} else {
			status |= DT_BUFFER_TOO_SMALL
		}

		for i := bestPoly.firstLink; i != DT_NULL_LINK; i = bestTile.links[i].next {
			link := bestTile.links[i]
			neighbourRef := link.ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour

			neighbourTile, neighbourPoly := q.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.

			va, vb, status := q.getPortalPoints1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile)
			if status.dtStatusFailed() {
				continue
			}

			// If the circle is not touching the next polygon, skip it.

			_, distSqr := dtDistancePtSegSqr2D(centerPos, va, vb)
			if distSqr > radiusSqr {
				continue
			}

			neighbourNode := q.m_nodePool[neighbourRef]
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.flags == 0 {
				copy(neighbourNode.pos[:], dtVlerp(va, vb, 0.5))
			}

			cost := filter.getCost(bestNode.pos[:], neighbourNode.pos[:], bestPoly)

			total := bestNode.total + cost

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_OPEN > 0) && total >= neighbourNode.total {
				continue
			}

			neighbourNode.id = neighbourRef
			neighbourNode.pidx = bestNode
			neighbourNode.total = total

			if neighbourNode.flags&DT_NODE_OPEN > 0 {
				q.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.flags = DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}
		}
	}

	resultCount = n

	return resultCount, status
}

// / @par
// /
// / The order of the result set is from least to highest cost.
// /
// / At least one result array must be provided.
// /
// / A common use case for this method is to perform Dijkstra searches.
// / Candidate polygons are found by searching the graph beginning at the start
// / polygon.
// /
// / The same intersection test restrictions that apply to findPolysAroundCircle()
// / method apply to this method.
// /
// / The 3D centroid of the search polygon is used as the start position for cost
// / calculations.
// /
// / Intersection tests occur in 2D. All polygons are projected onto the
// / xz-plane. So the y-values of the vertices do not effect intersection tests.
// /
// / If the result arrays are is too small to hold the entire result set, they will
// / be filled to capacity.
// /
func (q *dtNavMeshQuery) findPolysAroundShape(startRef dtPolyRef, verts []float64, nverts int,
	filter *dtQueryFilter,
	resultRef []dtPolyRef, resultParent []dtPolyRef, resultCost []float64, maxResult int) (resultCount int, status dtStatus) {
	dtAssert(q.m_nav)
	dtAssert(q.m_nodePool)
	dtAssert(q.m_openList)

	if !q.m_nav.isValidPolyRef(startRef) ||
		len(verts) == 0 || nverts < 3 ||
		filter == nil || maxResult < 0 {
		return resultCount, DT_FAILURE | DT_INVALID_PARAM
	}

	// Validate input
	if startRef == 0 || !q.m_nav.isValidPolyRef(startRef) {
		return resultCount, DT_FAILURE | DT_INVALID_PARAM
	}

	q.m_nodePool = map[dtPolyRef]*dtNode{}
	q.m_openList.Reset()

	var centerPos = []float64{0, 0, 0}
	for i := 0; i < nverts; i++ {
		centerPos = dtVadd(centerPos, rcGetVert(verts, i))

	}

	centerPos = dtVscale(centerPos, 1.0/float64(nverts))

	var startNode dtNode
	copy(startNode.pos[:], centerPos)
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = 0
	startNode.id = startRef
	startNode.flags = DT_NODE_OPEN
	q.m_nodePool[startRef] = &startNode
	q.m_openList.Offer(&startNode)

	status = DT_SUCCESS

	n := 0

	for !q.m_openList.Empty() {
		bestNode := q.m_openList.Poll()
		bestNode.flags &= ^DT_NODE_OPEN
		bestNode.flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.id

		bestTile, bestPoly := q.m_nav.getTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef dtPolyRef

		if bestNode.pidx != nil {
			parentRef = bestNode.pidx.id
		}

		if parentRef != 0 {
			q.m_nav.getTileAndPolyByRefUnsafe(parentRef)
		}

		if n < maxResult {
			if len(resultRef) > 0 {
				resultRef[n] = bestRef
			}

			if len(resultParent) > 0 {
				resultParent[n] = parentRef
			}

			if len(resultCost) > 0 {
				resultCost[n] = bestNode.total
			}

			n++
		} else {
			status |= DT_BUFFER_TOO_SMALL
		}

		for i := bestPoly.firstLink; i != DT_NULL_LINK; i = bestTile.links[i].next {
			link := bestTile.links[i]
			neighbourRef := link.ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour

			neighbourTile, neighbourPoly := q.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.

			va, vb, status := q.getPortalPoints1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile)
			if status.dtStatusFailed() {
				continue
			}

			// If the poly is not touching the edge to the next polygon, skip the connection it.

			tmin, tmax, _, _, ok := dtIntersectSegmentPoly2D(va, vb, verts, nverts)
			if !ok {
				continue
			}

			if tmin > 1.0 || tmax < 0.0 {
				continue
			}

			neighbourNode := q.m_nodePool[neighbourRef]
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.flags == 0 {
				copy(neighbourNode.pos[:], dtVlerp(va, vb, 0.5))

			}

			cost := filter.getCost(bestNode.pos[:], neighbourNode.pos[:], bestPoly)

			total := bestNode.total + cost

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_OPEN > 0) && total >= neighbourNode.total {
				continue
			}

			neighbourNode.id = neighbourRef
			neighbourNode.pidx = bestNode
			neighbourNode.total = total

			if neighbourNode.flags&DT_NODE_OPEN > 0 {
				q.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.flags = DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}
		}
	}

	resultCount = n

	return resultCount, status
}

func (q *dtNavMeshQuery) getPathFromDijkstraSearch(endRef dtPolyRef, path []dtPolyRef, maxPath int) (pathCount int, status dtStatus) {
	if !q.m_nav.isValidPolyRef(endRef) || len(path) == 0 || maxPath < 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	endNode := q.m_nodePool[endRef]
	if endNode != nil || (endNode.flags&DT_NODE_CLOSED) == 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	pathCount, status = q.getPathToNode(endNode, path, maxPath)
	return
}

func (q *dtNavMeshQuery) getPathToNode(endNode *dtNode, path []dtPolyRef, maxPath int) (pathCount int, status dtStatus) {
	// Find the length of the entire path.
	curNode := endNode
	var length = 0
	for curNode != nil {
		length++
		curNode = curNode.pidx
	}

	// If the path cannot be fully stored then advance to the last node we will be able to store.
	curNode = endNode
	var writeCount int
	for writeCount = length; writeCount > maxPath; writeCount-- {
		dtAssert(curNode)

		curNode = curNode.pidx
	}

	// Write path
	for i := writeCount - 1; i >= 0; i-- {
		dtAssert(curNode)

		path[i] = curNode.id
		curNode = curNode.pidx
	}

	dtAssertTrue(curNode != nil)

	pathCount = dtMin(length, maxPath)

	if length > maxPath {
		return pathCount, DT_SUCCESS | DT_BUFFER_TOO_SMALL
	}

	return pathCount, DT_SUCCESS
}

func (q *dtNavMeshQuery) finalizeSlicedFindPath(path []dtPolyRef, maxPath int) (pathCount int, status dtStatus) {

	if len(path) == 0 || maxPath <= 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	//if (dtStatusFailed(q.m_query.status)) {//TODO
	//// Reset query.
	//memset(&m_query, 0, sizeof(dtQueryData));
	//return pathCount,DT_FAILURE;
	//}

	n := 0

	if q.m_query.startRef == q.m_query.endRef {
		// Special case: the search starts and ends at same poly.
		path[n] = q.m_query.startRef
		n++
	} else {
		// Reverse the path.
		dtAssert(q.m_query.lastBestNode)

		if q.m_query.lastBestNode.id != q.m_query.endRef {
			q.m_query.status |= DT_PARTIAL_RESULT
		}

		var prev *dtNode
		node := q.m_query.lastBestNode
		prevRay := 0
		for node != nil {
			next := node.pidx
			node.pidx = prev
			prev = node
			nextRay := node.flags & DT_NODE_PARENT_DETACHED                // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)
			node.flags = (node.flags & ^DT_NODE_PARENT_DETACHED) | prevRay // and store it in the reversed path's node
			prevRay = nextRay
			node = next
		}

		// Store path
		node = prev
		for node != nil {
			next := node.pidx
			status = 0
			if node.flags&DT_NODE_PARENT_DETACHED > 0 {
				var normal [3]float64
				_, m, status := q.raycast(node.id, node.pos[:], next.pos[:], q.m_query.filter, normal[:], path[n:], maxPath-n)
				n += m
				// raycast ends on poly boundary and the path might include the next poly boundary.
				if path[n-1] == next.id {
					n-- // remove to avoid duplicates}

				} else {
					path[n] = node.id
					n++
					if n >= maxPath {
						status = DT_BUFFER_TOO_SMALL
					}

				}

				if status&DT_STATUS_DETAIL_MASK > 0 {
					q.m_query.status |= status & DT_STATUS_DETAIL_MASK
					break
				}
				node = next
			}

		}
	}
	details := q.m_query.status & DT_STATUS_DETAIL_MASK

	// Reset query.
	q.m_query = dtQueryData{}

	pathCount = n

	return pathCount, DT_SUCCESS | details
}

// / Provides information about raycast hit
// / filled by dtNavMeshQuery::raycast
// / @ingroup detour
type dtRaycastHit struct {

	/// The hit parameter. (FLT_MAX if no wall hit.)
	t float64

	/// hitNormal	The normal of the nearest wall hit. [(x, y, z)]
	hitNormal [3]float64

	/// The index of the edge on the final polygon where the wall was hit.
	hitEdgeIndex int

	/// Pointer to an array of reference ids of the visited polygons. [opt]
	path []dtPolyRef

	/// The number of visited polygons. [opt]
	pathCount int

	/// The maximum number of polygons the @p path array can hold.
	maxPath int

	///  The cost of the path until hit.
	pathCost float64
}

// / @par
// /
// / This method is meant to be used for quick, short distance checks.
// /
// / If the path array is too small to hold the result, it will be filled as
// / far as possible from the start postion toward the end position.
// /
// / <b>Using the Hit Parameter (t)</b>
// /
// / If the hit parameter is a very high value (FLT_MAX), then the ray has hit
// / the end position. In this case the path represents a valid corridor to the
// / end position and the value of @p hitNormal is undefined.
// /
// / If the hit parameter is zero, then the start position is on the wall that
// / was hit and the value of @p hitNormal is undefined.
// /
// / If 0 < t < 1.0 then the following applies:
// /
// / @code
// / distanceToHitBorder = distanceToEndPosition * t
// / hitPoint = startPos + (endPos - startPos) * t
// / @endcode
// /
// / <b>Use Case Restriction</b>
// /
// / The raycast ignores the y-value of the end position. (2D check.) This
// / places significant limits on how it can be used. For example:
// /
// / Consider a scene where there is a main floor with a second floor balcony
// / that hangs over the main floor. So the first floor mesh extends below the
// / balcony mesh. The start position is somewhere on the first floor. The end
// / position is on the balcony.
// /
// / The raycast will search toward the end position along the first floor mesh.
// / If it reaches the end position's xz-coordinates it will indicate FLT_MAX
// / (no wall hit), meaning it reached the end position. This is one example of why
// / this method is meant for short distance checks.
// /
func (q *dtNavMeshQuery) raycast(startRef dtPolyRef, startPos, endPos []float64,
	filter *dtQueryFilter, hitNormal []float64, path []dtPolyRef, maxPath int) (t float64, pathCount int, status dtStatus) {
	var hit dtRaycastHit
	hit.path = path
	hit.maxPath = maxPath

	status = q.raycast1(startRef, startPos, endPos, filter, 0, &hit, 0)

	t = hit.t
	if len(hitNormal) > 0 {
		copy(hitNormal, hit.hitNormal[:])
	}
	pathCount = hit.pathCount
	return t, pathCount, status
}

// / @par
// /
// / This method is meant to be used for quick, short distance checks.
// /
// / If the path array is too small to hold the result, it will be filled as
// / far as possible from the start postion toward the end position.
// /
// / <b>Using the Hit Parameter t of RaycastHit</b>
// /
// / If the hit parameter is a very high value (FLT_MAX), then the ray has hit
// / the end position. In this case the path represents a valid corridor to the
// / end position and the value of @p hitNormal is undefined.
// /
// / If the hit parameter is zero, then the start position is on the wall that
// / was hit and the value of @p hitNormal is undefined.
// /
// / If 0 < t < 1.0 then the following applies:
// /
// / @code
// / distanceToHitBorder = distanceToEndPosition * t
// / hitPoint = startPos + (endPos - startPos) * t
// / @endcode
// /
// / <b>Use Case Restriction</b>
// /
// / The raycast ignores the y-value of the end position. (2D check.) This
// / places significant limits on how it can be used. For example:
// /
// / Consider a scene where there is a main floor with a second floor balcony
// / that hangs over the main floor. So the first floor mesh extends below the
// / balcony mesh. The start position is somewhere on the first floor. The end
// / position is on the balcony.
// /
// / The raycast will search toward the end position along the first floor mesh.
// / If it reaches the end position's xz-coordinates it will indicate FLT_MAX
// / (no wall hit), meaning it reached the end position. This is one example of why
// / this method is meant for short distance checks.
// /
func (q *dtNavMeshQuery) raycast1(startRef dtPolyRef, startPos, endPos []float64,
	filter *dtQueryFilter, options int, hit *dtRaycastHit, prevRef dtPolyRef) (status dtStatus) {
	dtAssert(q.m_nav)

	if hit == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	hit.t = 0
	hit.pathCount = 0
	hit.pathCost = 0

	// Validate input
	if !q.m_nav.isValidPolyRef(startRef) || len(startPos) == 0 ||
		len(endPos) == 0 || filter == nil || (prevRef == 0 && !q.m_nav.isValidPolyRef(prevRef)) {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	var dir []float64
	var curPos, lastPos [3]float64
	var verts [DT_VERTS_PER_POLYGON*3 + 3]float64
	n := 0

	copy(curPos[:], startPos)
	dir = dtVsub(endPos[:], startPos)
	dtVset(hit.hitNormal[:], 0, 0, 0)

	status = DT_SUCCESS

	var prevTile, tile, nextTile *dtMeshTile
	var prevPoly, poly, nextPoly *dtPoly
	var curRef dtPolyRef

	// The API input has been checked already, skip checking internal data.
	curRef = startRef
	tile = nil
	poly = nil
	tile, poly = q.m_nav.getTileAndPolyByRefUnsafe(curRef)
	prevTile = tile
	prevPoly = poly
	nextTile = prevTile
	nextPoly = prevPoly
	if prevRef > 0 {
		prevTile, prevPoly = q.m_nav.getTileAndPolyByRefUnsafe(prevRef)
	}

	for curRef != 0 {
		// Cast ray against current polygon.

		// Collect vertices.
		nv := 0
		for i := 0; i < poly.vertCount; i++ {
			copy(verts[nv*3:nv*3+3], rcGetVert(tile.verts, poly.verts[i]))
			nv++
		}

		_, tmax, _, segMax, ok := dtIntersectSegmentPoly2D(startPos, endPos, verts[:], nv)
		if !ok {
			// Could not hit the polygon, keep the old t and report hit.
			hit.pathCount = n
			return status
		}

		hit.hitEdgeIndex = segMax

		// Keep track of furthest t so far.
		if tmax > hit.t {
			hit.t = tmax
		}

		// Store visited polygons.
		if n < hit.maxPath {
			hit.path[n] = curRef
			n++
		} else {
			status |= DT_BUFFER_TOO_SMALL
		}

		// Ray end is completely inside the polygon.
		if segMax == -1 {
			hit.t = math.MaxFloat64
			hit.pathCount = n

			// add the cost
			if options&DT_RAYCAST_USE_COSTS > 0 {
				hit.pathCost += filter.getCost(curPos[:], endPos, poly)
			}

			return status
		}

		// Follow neighbours.
		var nextRef dtPolyRef

		for i := poly.firstLink; i != DT_NULL_LINK; i = tile.links[i].next {
			link := tile.links[i]

			// Find link which contains this edge.
			if link.edge != segMax {
				continue
			}

			// Get pointer to the next polygon.

			nextTile, nextPoly = q.m_nav.getTileAndPolyByRefUnsafe(link.ref)

			// Skip off-mesh connections.
			if nextPoly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Skip links based on filter.
			if !filter.passFilter(nextPoly) {
				continue
			}

			// If the link is internal, just return the ref.
			if link.side == 0xff {
				nextRef = link.ref
				break
			}

			// If the link is at tile boundary,

			// Check if the link spans the whole edge, and accept.
			if link.bmin == 0 && link.bmax == 255 {
				nextRef = link.ref
				break
			}

			// Check for partial edge links.
			v0 := poly.verts[link.edge]
			v1 := poly.verts[(link.edge+1)%poly.vertCount]
			left := rcGetVert(tile.verts, v0)
			right := rcGetVert(tile.verts, v1)

			// Check that the intersection lies inside the link portal.
			if link.side == 0 || link.side == 4 {
				// Calculate link size.
				s := 1.0 / 255.0
				lmin := left[2] + (right[2]-left[2])*(float64(link.bmin)*s)
				lmax := left[2] + (right[2]-left[2])*(float64(link.bmax)*s)
				if lmin > lmax {
					lmin, lmax = lmax, lmin
				}

				// Find Z intersection.
				z := startPos[2] + (endPos[2]-startPos[2])*tmax
				if z >= lmin && z <= lmax {
					nextRef = link.ref
					break
				}
			} else if link.side == 2 || link.side == 6 {
				// Calculate link size.
				s := 1.0 / 255.0
				lmin := left[0] + (right[0]-left[0])*(float64(link.bmin)*s)
				lmax := left[0] + (right[0]-left[0])*(float64(link.bmax)*s)
				if lmin > lmax {
					lmin, lmax = lmax, lmin
				}

				// Find X intersection.
				x := startPos[0] + (endPos[0]-startPos[0])*tmax
				if x >= lmin && x <= lmax {
					nextRef = link.ref
					break
				}
			}
		}

		// add the cost
		if options&DT_RAYCAST_USE_COSTS > 0 {
			// compute the intersection point at the furthest end of the polygon
			// and correct the height (since the raycast moves in 2d)
			copy(lastPos[:], curPos[:])
			dtVmad(curPos[:], startPos, dir, hit.t)
			e1 := rcGetVert(verts[:], segMax)
			e2 := rcGetVert(verts[:], (segMax+1)%nv)

			eDir := dtVsub(e2, e1)
			diff := dtVsub(curPos[:], e1)
			s := diff[2] / eDir[2]
			if dtSqr(eDir[0]) > dtSqr(eDir[2]) {
				s = diff[0] / eDir[0]
			}
			curPos[1] = e1[1] + eDir[1]*s

			hit.pathCost += filter.getCost(lastPos[:], curPos[:], poly)
		}

		if nextRef == 0 {
			// No neighbour, we hit a wall.

			// Calculate hit normal.
			a := segMax
			b := 0
			if segMax+1 < nv {
				b = segMax + 1
			}
			va := rcGetVert(verts[:], a)
			vb := rcGetVert(verts[:], b)
			dx := vb[0] - va[0]
			dz := vb[2] - va[2]
			hit.hitNormal[0] = dz
			hit.hitNormal[1] = 0
			hit.hitNormal[2] = -dx
			dtVnormalize(hit.hitNormal[:])

			hit.pathCount = n
			return status
		}

		// No hit, advance to neighbour polygon.
		prevRef = curRef
		curRef = nextRef
		prevTile = tile
		tile = nextTile
		prevPoly = poly
		poly = nextPoly
	}

	hit.pathCount = n

	return status
}

// / @par
// /
// / @warning Calling any non-slice methods before calling finalizeSlicedFindPath()
// / or finalizeSlicedFindPathPartial() may result in corrupted data!
// /
// / The @p filter pointer is stored and used for the duration of the sliced
// / path query.
// /
func (q *dtNavMeshQuery) initSlicedFindPath(startRef, endRef dtPolyRef,
	startPos, endPos []float64, filter *dtQueryFilter, options int) dtStatus {
	dtAssert(q.m_nav)
	dtAssert(q.m_nodePool)
	dtAssert(q.m_openList)

	// Init path state.
	q.m_query = dtQueryData{}
	q.m_query.status = DT_FAILURE
	q.m_query.startRef = startRef
	q.m_query.endRef = endRef
	if len(startPos) > 0 {
		copy(q.m_query.startPos[:], startPos)
	}

	if (len(endPos)) > 0 {
		copy(q.m_query.endPos[:], endPos)
	}

	q.m_query.filter = filter
	q.m_query.options = options
	q.m_query.raycastLimitSqr = math.MaxFloat64

	// Validate input
	if !q.m_nav.isValidPolyRef(startRef) || !q.m_nav.isValidPolyRef(endRef) ||
		len(startPos) == 0 ||
		len(endPos) == 0 || filter == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	// trade quality with performance?
	if options&DT_FINDPATH_ANY_ANGLE > 0 {
		// limiting to several times the character radius yields nice results. It is not sensitive
		// so it is enough to compute it from the first tile.
		tile := q.m_nav.getTileByRef(dtTileRef(startRef))
		agentRadius := tile.header.walkableRadius
		q.m_query.raycastLimitSqr = dtSqr(float64(agentRadius * DT_RAY_CAST_LIMIT_PROPORTIONS))
	}

	if startRef == endRef {
		q.m_query.status = DT_SUCCESS
		return DT_SUCCESS
	}

	q.m_nodePool = map[dtPolyRef]*dtNode{}
	q.m_openList.Reset()
	var startNode dtNode

	copy(startNode.pos[:], startPos)
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = dtVdist(startPos, endPos) * H_SCALE
	startNode.id = startRef
	startNode.flags = DT_NODE_OPEN
	q.m_nodePool[startRef] = &startNode
	q.m_openList.Offer(&startNode)

	q.m_query.status = DT_IN_PROGRESS
	q.m_query.lastBestNode = &startNode
	q.m_query.lastBestNodeCost = startNode.total

	return q.m_query.status
}

func (q *dtNavMeshQuery) updateSlicedFindPath(maxIter int) (doneIters int, status dtStatus) {
	if !q.m_query.status.dtStatusInProgress() {
		return doneIters, q.m_query.status
	}

	// Make sure the request is still valid.
	if !q.m_nav.isValidPolyRef(q.m_query.startRef) || !q.m_nav.isValidPolyRef(q.m_query.endRef) {
		q.m_query.status = DT_FAILURE
		return doneIters, DT_FAILURE
	}

	var rayHit dtRaycastHit
	rayHit.maxPath = 0

	iter := 0
	for iter < maxIter && !q.m_openList.Empty() {
		iter++

		// Remove node from open list and put it in closed list.
		bestNode := q.m_openList.Poll()
		bestNode.flags &= ^DT_NODE_OPEN
		bestNode.flags |= DT_NODE_CLOSED

		// Reached the goal, stop searching.
		if bestNode.id == q.m_query.endRef {
			q.m_query.lastBestNode = bestNode
			details := q.m_query.status & DT_STATUS_DETAIL_MASK
			q.m_query.status = DT_SUCCESS | details
			doneIters = iter

			return doneIters, q.m_query.status
		}

		// Get current poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.id
		bestTile, bestPoly, status := q.m_nav.getTileAndPolyByRef(bestRef)
		if status.dtStatusFailed() {
			// The polygon has disappeared during the sliced query, fail.
			q.m_query.status = DT_FAILURE
			doneIters = iter

			return doneIters, q.m_query.status
		}

		// Get parent and grand parent poly and tile.
		var parentRef, grandpaRef dtPolyRef
		var parentNode *dtNode
		if bestNode.pidx != nil {
			parentNode = bestNode.pidx
			parentRef = parentNode.id
			if parentNode.pidx != nil {
				grandpaRef = parentNode.pidx.id
			}

		}
		if parentRef > 0 {
			_, _, status := q.m_nav.getTileAndPolyByRef(parentRef)
			invalidParent := status.dtStatusFailed()
			if invalidParent || (grandpaRef > 0 && !q.m_nav.isValidPolyRef(grandpaRef)) {
				// The polygon has disappeared during the sliced query, fail.
				q.m_query.status = DT_FAILURE

				doneIters = iter
				return doneIters, q.m_query.status
			}
		}

		// decide whether to test raycast to previous nodes
		tryLOS := false
		if q.m_query.options&DT_FINDPATH_ANY_ANGLE > 0 {
			if (parentRef != 0) && (dtVdistSqr(parentNode.pos[:], bestNode.pos[:]) < q.m_query.raycastLimitSqr) {
				tryLOS = true
			}

		}

		for i := bestPoly.firstLink; i != DT_NULL_LINK; i = bestTile.links[i].next {
			neighbourRef := bestTile.links[i].ref

			// Skip invalid ids and do not expand back to where we came from.
			if neighbourRef == 0 || neighbourRef == parentRef {
				continue
			}

			// Get neighbour poly and tile.
			// The API input has been checked already, skip checking internal data.

			neighbourTile, neighbourPoly := q.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			if !q.m_query.filter.passFilter(neighbourPoly) {
				continue
			}

			// get the neighbor node
			neighbourNode := q.m_nodePool[neighbourRef]
			if neighbourNode == nil {
				q.m_query.status |= DT_OUT_OF_NODES
				continue
			}

			// do not expand to nodes that were already visited from the same parent
			if neighbourNode.pidx != nil && neighbourNode.pidx == bestNode.pidx {
				continue
			}

			// If the node is visited the first time, calculate node position.
			if neighbourNode.flags == 0 {
				pos, _ := q.getEdgeMidPoint1(bestRef, bestPoly, bestTile,
					neighbourRef, neighbourPoly, neighbourTile,
				)
				copy(neighbourNode.pos[:], pos)
			}

			// Calculate cost and heuristic.
			cost := float64(0)
			heuristic := float64(0)

			// raycast parent
			foundShortCut := false
			rayHit.t = 0
			rayHit.pathCost = rayHit.t
			if tryLOS {
				q.raycast1(parentRef, parentNode.pos[:], neighbourNode.pos[:], q.m_query.filter, DT_RAYCAST_USE_COSTS, &rayHit, grandpaRef)
				foundShortCut = rayHit.t >= 1.0
			}

			// update move cost
			if foundShortCut {
				// shortcut found using raycast. Using shorter cost instead
				cost = parentNode.cost + rayHit.pathCost
			} else {
				// No shortcut found.
				curCost := q.m_query.filter.getCost(bestNode.pos[:], neighbourNode.pos[:], bestPoly)
				cost = bestNode.cost + curCost
			}

			// Special case for last node.
			if neighbourRef == q.m_query.endRef {
				endCost := q.m_query.filter.getCost(neighbourNode.pos[:], q.m_query.endPos[:], neighbourPoly)

				cost = cost + endCost
				heuristic = 0
			} else {
				heuristic = dtVdist(neighbourNode.pos[:], q.m_query.endPos[:]) * H_SCALE
			}

			total := cost + heuristic

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_OPEN > 0) && total >= neighbourNode.total {
				continue
			}

			// The node is already visited and process, and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_CLOSED > 0) && total >= neighbourNode.total {
				continue
			}

			// Add or update the node.
			neighbourNode.pidx = bestNode
			if foundShortCut {
				neighbourNode.pidx = bestNode.pidx
			}

			neighbourNode.id = neighbourRef
			neighbourNode.flags = (neighbourNode.flags & ^(DT_NODE_CLOSED | DT_NODE_PARENT_DETACHED))
			neighbourNode.cost = cost
			neighbourNode.total = total
			if foundShortCut {
				neighbourNode.flags = (neighbourNode.flags | DT_NODE_PARENT_DETACHED)
			}

			if neighbourNode.flags&DT_NODE_OPEN > 0 {
				// Already in open, update node location.
				q.m_openList.Update(neighbourNode)
			} else {
				// Put the node in open list.
				neighbourNode.flags |= DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}

			// Update nearest node to target so far.
			if heuristic < q.m_query.lastBestNodeCost {
				q.m_query.lastBestNodeCost = heuristic
				q.m_query.lastBestNode = neighbourNode
			}
		}
	}

	// Exhausted all nodes, but could not find path.
	if q.m_openList.Empty() {
		details := q.m_query.status & DT_STATUS_DETAIL_MASK
		q.m_query.status = DT_SUCCESS | details
	}

	doneIters = iter

	return doneIters, q.m_query.status
}

// / @par
// /
// / If the end polygon cannot be reached through the navigation graph,
// / the last polygon in the path will be the nearest the end polygon.
// /
// / If the path array is to small to hold the full result, it will be filled as
// / far as possible from the start polygon toward the end polygon.
// /
// / The start and end positions are used to calculate traversal costs.
// / (The y-values impact the result.)
// /
func (q *dtNavMeshQuery) findPath(startRef, endRef dtPolyRef,
	startPos, endPos []float64, filter *dtQueryFilter, path []dtPolyRef, maxPath int) (pathCount int, status dtStatus) {
	dtAssert(q.m_nav)
	dtAssert(q.m_nodePool)
	dtAssert(q.m_openList)
	// Validate input
	if !q.m_nav.isValidPolyRef(startRef) || !q.m_nav.isValidPolyRef(endRef) ||
		len(startPos) == 0 ||
		len(endPos) == 0 ||
		filter == nil || len(path) == 0 || maxPath <= 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	if startRef == endRef {
		path[0] = startRef
		pathCount = 1
		return pathCount, DT_SUCCESS
	}

	q.m_nodePool = map[dtPolyRef]*dtNode{}
	q.m_openList.Reset()
	var startNode dtNode
	copy(startNode.pos[:], startPos)
	startNode.pidx = nil
	startNode.cost = 0
	startNode.total = dtVdist(startPos, endPos) * H_SCALE
	startNode.id = startRef
	startNode.flags = DT_NODE_OPEN
	q.m_nodePool[startRef] = &startNode
	q.m_openList.Offer(&startNode)

	lastBestNode := &startNode
	lastBestNodeCost := startNode.total

	outOfNodes := false

	for !q.m_openList.Empty() {
		// Remove node from open list and put it in closed list.
		bestNode := q.m_openList.Poll()
		bestNode.flags &= ^DT_NODE_OPEN
		bestNode.flags |= DT_NODE_CLOSED

		// Reached the goal, stop searching.
		if bestNode.id == endRef {
			lastBestNode = bestNode
			break
		}

		// Get current poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.id
		bestTile, bestPoly := q.m_nav.getTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef dtPolyRef

		if bestNode.pidx != nil {
			parentRef = bestNode.pidx.id
		}

		if parentRef > 0 {
			q.m_nav.getTileAndPolyByRefUnsafe(parentRef)
		}

		for i := bestPoly.firstLink; i != DT_NULL_LINK; i = bestTile.links[i].next {
			neighbourRef := bestTile.links[i].ref

			// Skip invalid ids and do not expand back to where we came from.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Get neighbour poly and tile.
			// The API input has been checked already, skip checking internal data.

			neighbourTile, neighbourPoly := q.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// deal explicitly with crossing tile boundaries
			crossSide := 0
			if bestTile.links[i].side != 0xff {
				crossSide = bestTile.links[i].side >> 1
			}
			_ = crossSide
			//TDOO

			// get the node
			neighbourNode := q.m_nodePool[neighbourRef]
			if neighbourNode == nil {
				outOfNodes = true
				continue
			}

			// If the node is visited the first time, calculate node position.
			if neighbourNode.flags == 0 {
				pos, _ := q.getEdgeMidPoint1(bestRef, bestPoly, bestTile,
					neighbourRef, neighbourPoly, neighbourTile)
				copy(neighbourNode.pos[:], pos)

			}

			// Calculate cost and heuristic.
			cost := float64(0)
			heuristic := float64(0)

			// Special case for last node.
			if neighbourRef == endRef {
				curCost := filter.getCost(bestNode.pos[:], neighbourNode.pos[:], bestPoly)
				// Cost
				endCost := filter.getCost(neighbourNode.pos[:], endPos, neighbourPoly)

				cost = bestNode.cost + curCost + endCost
				heuristic = 0
			} else {
				// Cost
				curCost := filter.getCost(bestNode.pos[:], neighbourNode.pos[:], bestPoly)
				cost = bestNode.cost + curCost
				heuristic = dtVdist(neighbourNode.pos[:], endPos) * H_SCALE
			}

			total := cost + heuristic

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_OPEN > 0) && total >= neighbourNode.total {
				continue
			}

			// The node is already visited and process, and the new result is worse, skip.
			if (neighbourNode.flags&DT_NODE_CLOSED > 0) && total >= neighbourNode.total {
				continue
			}

			// Add or update the node.
			neighbourNode.pidx = bestNode
			neighbourNode.id = neighbourRef
			neighbourNode.flags = (neighbourNode.flags & ^DT_NODE_CLOSED)
			neighbourNode.cost = cost
			neighbourNode.total = total

			if neighbourNode.flags&DT_NODE_OPEN > 0 {
				// Already in open, update node location.
				q.m_openList.Update(neighbourNode)
			} else {
				// Put the node in open list.
				neighbourNode.flags |= DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}

			// Update nearest node to target so far.
			if heuristic < lastBestNodeCost {
				lastBestNodeCost = heuristic
				lastBestNode = neighbourNode
			}
		}
	}

	pathCount, status = q.getPathToNode(lastBestNode, path, maxPath)

	if lastBestNode.id != endRef {
		status |= DT_PARTIAL_RESULT
	}

	if outOfNodes {
		status |= DT_OUT_OF_NODES
	}

	return pathCount, status
}

// / @par
// /
// / The query will be invoked with batches of polygons. Polygons passed
// / to the query have bounding boxes that overlap with the center and halfExtents
// / passed to this function. The dtPolyQuery::process function is invoked multiple
// / times until all overlapping polygons have been processed.
// /
func (q *dtNavMeshQuery) queryPolygons(center []float64, halfExtents []float64,
	filter *dtQueryFilter, query dtPolyQuery) dtStatus {
	dtAssert(q.m_nav)

	if len(center) == 0 ||
		len(halfExtents) == 0 || !dtVisfinite(halfExtents) ||
		filter == nil || query == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	bmin := dtVsub(center, halfExtents)
	bmax := dtVadd(center, halfExtents)

	// Find tiles the query touches.

	minx, miny := q.m_nav.calcTileLoc(bmin)
	maxx, maxy := q.m_nav.calcTileLoc(bmax)

	MAX_NEIS := 32

	for y := miny; y <= maxy; y++ {
		for x := minx; x <= maxx; x++ {
			neis, nneis := q.m_nav.getTilesAt(x, y, MAX_NEIS)
			for j := 0; j < nneis; j++ {
				q.queryPolygonsInTile(neis[j], bmin, bmax, filter, query)
			}
		}
	}

	return DT_SUCCESS
}

func (q *dtNavMeshQuery) finalizeSlicedFindPathPartial(existing []dtPolyRef, existingSize int,
	path []dtPolyRef, maxPath int) (pathCount int, status dtStatus) {

	if len(existing) == 0 || existingSize <= 0 || len(path) == 0 || maxPath <= 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	if q.m_query.status.dtStatusFailed() {
		// Reset query.
		q.m_query = dtQueryData{}
		return pathCount, DT_FAILURE
	}

	n := 0

	if q.m_query.startRef == q.m_query.endRef {
		// Special case: the search starts and ends at same poly.
		path[n] = q.m_query.startRef
		n++
	} else {
		// Find furthest existing node that was visited.
		var prev *dtNode
		var node *dtNode
		for i := existingSize - 1; i >= 0; i-- {
			node := q.m_nodePool[existing[i]]
			if node != nil {
				break
			}

		}

		if node == nil {
			q.m_query.status |= DT_PARTIAL_RESULT
			dtAssert(q.m_query.lastBestNode)
			node = q.m_query.lastBestNode
		}

		// Reverse the path.
		prevRay := 0
		for node != nil {
			next := node.pidx
			node.pidx = prev
			prev = node
			nextRay := node.flags & DT_NODE_PARENT_DETACHED                // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)
			node.flags = (node.flags & ^DT_NODE_PARENT_DETACHED) | prevRay // and store it in the reversed path's node
			prevRay = nextRay
			node = next
		}

		// Store path
		node = prev
		for node != nil {
			next := node.pidx
			status = 0
			if node.flags&DT_NODE_PARENT_DETACHED > 0 {
				var normal [3]float64
				_, m, _ := q.raycast(node.id, node.pos[:], next.pos[:], q.m_query.filter, normal[:], path[n:], maxPath-n)
				n += m
				// raycast ends on poly boundary and the path might include the next poly boundary.
				if path[n-1] == next.id {
					n--
				} // remove to avoid duplicates

			} else {
				path[n] = node.id
				n++
				if n >= maxPath {
					status = DT_BUFFER_TOO_SMALL
				}

			}

			if status&DT_STATUS_DETAIL_MASK > 0 {
				q.m_query.status |= status & DT_STATUS_DETAIL_MASK
				break
			}
			node = next
		}

	}

	details := q.m_query.status & DT_STATUS_DETAIL_MASK

	// Reset query.
	q.m_query = dtQueryData{}

	pathCount = n

	return pathCount, DT_SUCCESS | details
}

// / @par
// /
// / This method is optimized for a small search radius and small number of result
// / polygons.
// /
// / Candidate polygons are found by searching the navigation graph beginning at
// / the start polygon.
// /
// / The same intersection test restrictions that apply to the findPolysAroundCircle
// / mehtod applies to this method.
// /
// / The value of the center point is used as the start point for cost calculations.
// / It is not projected onto the surface of the mesh, so its y-value will effect
// / the costs.
// /
// / Intersection tests occur in 2D. All polygons and the search circle are
// / projected onto the xz-plane. So the y-value of the center point does not
// / effect intersection tests.
// /
// / If the result arrays are is too small to hold the entire result set, they will
// / be filled to capacity.
// /
func (q *dtNavMeshQuery) findLocalNeighbourhood(startRef dtPolyRef, centerPos []float64, radius float64,
	filter *dtQueryFilter,
	resultRef []dtPolyRef, resultParent []dtPolyRef, maxResult int) (resultCount int, status dtStatus) {
	dtAssert(q.m_nav)
	dtAssert(q.m_tinyNodePool)
	if !q.m_nav.isValidPolyRef(startRef) ||
		len(centerPos) == 0 ||
		radius < 0 || !dtIsFinite(radius) ||
		filter == nil || maxResult < 0 {
		return resultCount, DT_FAILURE | DT_INVALID_PARAM
	}

	MAX_STACK := 48
	var stack [48]*dtNode
	nstack := 0

	q.m_tinyNodePool = map[dtPolyRef]*dtNode{}

	startNode := q.m_tinyNodePool[startRef]
	startNode.pidx = nil
	startNode.id = startRef
	startNode.flags = DT_NODE_CLOSED
	stack[nstack] = startNode
	nstack++

	radiusSqr := dtSqr(radius)

	var pa [DT_VERTS_PER_POLYGON * 3]float64
	var pb [DT_VERTS_PER_POLYGON * 3]float64

	status = DT_SUCCESS

	n := 0
	if n < maxResult {
		resultRef[n] = startNode.id
		if len(resultParent) > 0 {
			resultParent[n] = 0
		}

		n++
	} else {
		status |= DT_BUFFER_TOO_SMALL
	}

	for nstack > 0 {
		// Pop front.
		curNode := stack[0]
		for i := 0; i < nstack-1; i++ {
			stack[i] = stack[i+1]
		}

		nstack--

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		curRef := curNode.id

		curTile, curPoly := q.m_nav.getTileAndPolyByRefUnsafe(curRef)

		for i := curPoly.firstLink; i != DT_NULL_LINK; i = curTile.links[i].next {
			link := curTile.links[i]
			neighbourRef := link.ref
			// Skip invalid neighbours.
			if neighbourRef == 0 {
				continue
			}

			// Skip if cannot alloca more nodes.
			neighbourNode := q.m_tinyNodePool[neighbourRef]
			if neighbourNode == nil {
				continue
			}

			// Skip visited.
			if neighbourNode.flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Expand to neighbour

			neighbourTile, neighbourPoly := q.m_nav.getTileAndPolyByRefUnsafe(neighbourRef)

			// Skip off-mesh connections.
			if neighbourPoly.getType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.
			va, vb, status := q.getPortalPoints1(curRef, curPoly, curTile, neighbourRef, neighbourPoly, neighbourTile)
			if status.dtStatusFailed() {
				continue
			}

			// If the circle is not touching the next polygon, skip it.
			_, distSqr := dtDistancePtSegSqr2D(centerPos, va, vb)
			if distSqr > radiusSqr {
				continue
			}

			// Mark node visited, this is done before the overlap test so that
			// we will not visit the poly again if the test fails.
			neighbourNode.flags |= DT_NODE_CLOSED
			neighbourNode.pidx = curNode

			// Check that the polygon does not collide with existing polygons.

			// Collect vertices of the neighbour poly.
			npa := neighbourPoly.vertCount
			for k := 0; k < npa; k++ {
				copy(pa[k*3:k*3+3], rcGetVert(neighbourTile.verts, neighbourPoly.verts[k]))
			}

			overlap := false
			for j := 0; j < n; j++ {
				pastRef := resultRef[j]

				// Connected polys do not overlap.
				connected := false
				for k := curPoly.firstLink; k != DT_NULL_LINK; k = curTile.links[k].next {
					if curTile.links[k].ref == pastRef {
						connected = true
						break
					}
				}
				if connected {
					continue
				}

				// Potentially overlapping.

				pastTile, pastPoly := q.m_nav.getTileAndPolyByRefUnsafe(pastRef)

				// Get vertices and test overlap
				npb := pastPoly.vertCount

				for k := 0; k < npb; k++ {
					copy(pb[k*3:k*3+3], rcGetVert(pastTile.verts, pastPoly.verts[k]))
				}

				if dtOverlapPolyPoly2D(pa[:], npa, pb[:], npb) {
					overlap = true
					break
				}
			}
			if overlap {
				continue
			}

			// This poly is fine, store and advance to the poly.
			if n < maxResult {
				{
					resultRef[n] = neighbourRef
				}

				if len(resultParent) > 0 {
					resultParent[n] = curRef
				}

				n++
			} else {
				status |= DT_BUFFER_TOO_SMALL
			}

			if nstack < MAX_STACK {
				stack[nstack] = neighbourNode
				nstack++
			}
		}
	}

	resultCount = n

	return resultCount, status
}
