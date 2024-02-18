package detour

import (
	"github.com/gorustyt/gonavmesh/common"
	"math"
	"math/rand"
)

const (
	DT_FINDPATH_ANY_ANGLE = 0x02 ///< use raycasts during pathfind to "shortcut" (raycast still consider costs)
)
const (
	H_SCALE = 0.999 // Search heuristic scale.

	/// Vertex flags returned by DtNavMeshQuery::findStraightPath.
	DT_STRAIGHTPATH_START              = 0x01 ///< The vertex is the start position in the path.
	DT_STRAIGHTPATH_END                = 0x02 ///< The vertex is the end position in the path.
	DT_STRAIGHTPATH_OFFMESH_CONNECTION = 0x04 ///< The vertex is the start of an off-mesh connection.

)
const (
	/// Options for DtNavMeshQuery::findStraightPath.
	DT_STRAIGHTPATH_AREA_CROSSINGS = 0x01 ///< Add a vertex at every polygon edge crossing where area changes.
	DT_STRAIGHTPATH_ALL_CROSSINGS  = 0x02 ///< Add a vertex at every polygon edge crossing.
)

var (
	DT_VIRTUAL_QUERYFILTER = 1
)

type dtQueryData struct {
	status           DtStatus
	lastBestNode     *DtNode
	lastBestNodeCost float32
	startRef, endRef DtPolyRef
	startPos         [3]float32
	endPos           [3]float32
	filter           *DtQueryFilter
	options          int32
	raycastLimitSqr  float32
}

type NavMeshQuery interface {
	/// @name Standard Pathfinding Functions
	/// @{

	/// Finds a path from the start polygon to the end polygon.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		endRef		The reference id of the end polygon.
	///  @param[in]		startPos	A position within the start polygon. [(x, y, z)]
	///  @param[in]		endPos		A position within the end polygon. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.)
	///  							[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	The number of polygons returned in the @p path array.
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold. [Limit: >= 1]
	FindPath(startRef, endRef DtPolyRef,
		startPos, endPos []float32, filter *DtQueryFilter, path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus)
	/// Finds the straight path from the start to the end position within the polygon corridor.
	///  @param[in]		startPos			Path start position. [(x, y, z)]
	///  @param[in]		endPos				Path end position. [(x, y, z)]
	///  @param[in]		path				An array of polygon references that represent the path corridor.
	///  @param[in]		pathSize			The number of polygons in the @p path array.
	///  @param[out]	straightPath		Points describing the straight path. [(x, y, z) * @p straightPathCount].
	///  @param[out]	straightPathFlags	Flags describing each point. (See: #dtStraightPathFlags) [opt]
	///  @param[out]	straightPathRefs	The reference id of the polygon that is being entered at each point. [opt]
	///  @param[out]	straightPathCount	The number of points in the straight path.
	///  @param[in]		maxStraightPath		The maximum number of points the straight path arrays can hold.  [Limit: > 0]
	///  @param[in]		options				Query options. (see: #dtStraightPathOptions)
	/// @returns The status flags for the query.
	FindStraightPath(startPos, endPos []float32,
		path []DtPolyRef, pathSize int32,
		straightPath []float32, straightPathFlags []int32, straightPathRefs []DtPolyRef,
		maxStraightPath int32, options ...int32) (straightPathCount int32, status DtStatus)
	/// Initializes a sliced path query.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		endRef		The reference id of the end polygon.
	///  @param[in]		startPos	A position within the start polygon. [(x, y, z)]
	///  @param[in]		endPos		A position within the end polygon. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[in]		options		query options (see: #dtFindPathOptions)
	/// @returns The status flags for the query.
	InitSlicedFindPath(startRef, endRef DtPolyRef,
		startPos, endPos []float32, filter *DtQueryFilter, options int32) DtStatus
	/// Updates an in-progress sliced path query.
	///  @param[in]		maxIter		The maximum number of iterations to perform.
	///  @param[out]	doneIters	The actual number of iterations completed. [opt]
	/// @returns The status flags for the query.
	UpdateSlicedFindPath(maxIter int32) (doneIters int32, status DtStatus)
	/// Finalizes and returns the results of a sliced path query.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.)
	///  							[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	The number of polygons returned in the @p path array.
	///  @param[in]		maxPath		The max number of polygons the path array can hold. [Limit: >= 1]
	/// @returns The status flags for the query.
	FinalizeSlicedFindPath(path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus)
	/// Finalizes and returns the results of an incomplete sliced path query, returning the path to the furthest
	/// polygon on the existing path that was visited during the search.
	///  @param[in]		existing		An array of polygon references for the existing path.
	///  @param[in]		existingSize	The number of polygon in the @p existing array.
	///  @param[out]	path			An ordered list of polygon references representing the path. (Start to end.)
	///  								[(polyRef) * @p pathCount]
	///  @param[out]	pathCount		The number of polygons returned in the @p path array.
	///  @param[in]		maxPath			The max number of polygons the @p path array can hold. [Limit: >= 1]
	/// @returns The status flags for the query.
	FinalizeSlicedFindPathPartial(existing []DtPolyRef, existingSize int32,
		path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus)
	///@}
	/// @name Dijkstra Search Functions
	/// @{

	/// Finds the polygons along the navigation graph that touch the specified circle.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		radius			The radius of the search circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the circle. [opt]
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result.
	///  								Zero if a result polygon has no parent. [opt]
	///  @param[out]	resultCost		The search cost from @p centerPos to the polygon. [opt]
	///  @param[out]	resultCount		The number of polygons found. [opt]
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.
	FindPolysAroundCircle(startRef DtPolyRef, centerPos []float32, radius float32,
		filter *DtQueryFilter, resultRef []DtPolyRef, resultParent []DtPolyRef, resultCost []float32, resultCount *int32, maxResult int32) (status DtStatus)

	/// Finds the polygons along the naviation graph that touch the specified convex polygon.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		verts			The vertices describing the convex polygon. (CCW)
	///  								[(x, y, z) * @p nverts]
	///  @param[in]		nverts			The number of vertices in the polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the search polygon. [opt]
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result. Zero if a
	///  								result polygon has no parent. [opt]
	///  @param[out]	resultCost		The search cost from the centroid point to the polygon. [opt]
	///  @param[out]	resultCount		The number of polygons found.
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.
	FindPolysAroundShape(startRef DtPolyRef, verts []float32, nverts int32,
		filter *DtQueryFilter,
		resultRef []DtPolyRef, resultParent []DtPolyRef, resultCost []float32, maxResult int32) (resultCount int32, status DtStatus)
	/// Gets a path from the explored nodes in the previous search.
	///  @param[in]		endRef		The reference id of the end polygon.
	///  @param[out]	path		An ordered list of polygon references representing the path. (Start to end.)
	///  							[(polyRef) * @p pathCount]
	///  @param[out]	pathCount	The number of polygons returned in the @p path array.
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold. [Limit: >= 0]
	///  @returns		The status flags. Returns DT_FAILURE | DT_INVALID_PARAM if any parameter is wrong, or if
	///  				@p endRef was not explored in the previous search. Returns DT_SUCCESS | DT_BUFFER_TOO_SMALL
	///  				if @p path cannot contain the entire path. In this case it is filled to capacity with a partial path.
	///  				Otherwise returns DT_SUCCESS.
	///  @remarks		The result of this function depends on the state of the query object. For that reason it should only
	///  				be used immediately after one of the two Dijkstra searches, findPolysAroundCircle or findPolysAroundShape.
	GetPathFromDijkstraSearch(endRef DtPolyRef, path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus)
	/// Finds polygons that overlap the search box.
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents		The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[in]		query		The query. Polygons found will be batched together and passed to this query.
	QueryPolygons(center []float32, halfExtents []float32,
		filter *DtQueryFilter, query dtPolyQuery) DtStatus
	/// Finds polygons that overlap the search box.
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents		The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	polys		The reference ids of the polygons that overlap the query box.
	///  @param[out]	polyCount	The number of polygons in the search result.
	///  @param[in]		maxPolys	The maximum number of polygons the search result can hold.
	/// @returns The status flags for the query.
	QueryPolygons1(center, halfExtents []float32,
		filter *DtQueryFilter,
		polys []DtPolyRef, maxPolys int32) (polyCount int32, status DtStatus)
	/// Finds the non-overlapping navigation polygons in the local neighbourhood around the center position.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the query circle. [(x, y, z)]
	///  @param[in]		radius			The radius of the query circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultRef		The reference ids of the polygons touched by the circle.
	///  @param[out]	resultParent	The reference ids of the parent polygons for each result.
	///  								Zero if a result polygon has no parent. [opt]
	///  @param[out]	resultCount		The number of polygons found.
	///  @param[in]		maxResult		The maximum number of polygons the result arrays can hold.
	/// @returns The status flags for the query.
	FindLocalNeighbourhood(startRef DtPolyRef, centerPos []float32, radius float32,
		filter *DtQueryFilter,
		resultRef []DtPolyRef, resultParent []DtPolyRef, maxResult int32) (resultCount int32, status DtStatus)
	/// Moves from the start to the end position constrained to the navigation mesh.
	///  @param[in]		startRef		The reference id of the start polygon.
	///  @param[in]		startPos		A position of the mover within the start polygon. [(x, y, x)]
	///  @param[in]		endPos			The desired end position of the mover. [(x, y, z)]
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	resultPos		The result position of the mover. [(x, y, z)]
	///  @param[out]	visited			The reference ids of the polygons visited during the move.
	///  @param[out]	visitedCount	The number of polygons visited during the move.
	///  @param[in]		maxVisitedSize	The maximum number of polygons the @p visited array can hold.
	/// @returns The status flags for the query.
	MoveAlongSurface(startRef DtPolyRef, startPos, endPos []float32,
		filter *DtQueryFilter, resultPos []float32, visited []DtPolyRef, visitedCount *int32, maxVisitedSize int32) (status DtStatus)

	/// Finds the polygon nearest to the specified center point.
	/// [opt] means the specified parameter can be a null pointer, in that case the output parameter will not be set.
	///
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents	The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	nearestRef	The reference id of the nearest polygon. Will be set to 0 if no polygon is found.
	///  @param[out]	nearestPt	The nearest point on the polygon. Unchanged if no polygon is found. [opt] [(x, y, z)]
	///  @param[out]	isOverPoly 	Set to true if the point's X/Z coordinate lies inside the polygon, false otherwise. Unchanged if no polygon is found. [opt]
	/// @returns The status flags for the query.
	FindNearestPoly1(center, halfExtents []float32,
		filter *DtQueryFilter, nearestPt []float32) (nearestRef DtPolyRef, isOverPoly bool, status DtStatus)
	/// @}
	/// @name Local Query Functions
	///@{

	/// Finds the polygon nearest to the specified center point.
	/// [opt] means the specified parameter can be a null pointer, in that case the output parameter will not be set.
	///
	///  @param[in]		center		The center of the search box. [(x, y, z)]
	///  @param[in]		halfExtents	The search distance along each axis. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	nearestRef	The reference id of the nearest polygon. Will be set to 0 if no polygon is found.
	///  @param[out]	nearestPt	The nearest point on the polygon. Unchanged if no polygon is found. [opt] [(x, y, z)]
	/// @returns The status flags for the query.
	FindNearestPoly(center, halfExtents []float32, filter *DtQueryFilter, nearestPt []float32) (nearestRef DtPolyRef, status DtStatus)
	/// Casts a 'walkability' ray along the surface of the navigation mesh from
	/// the start position toward the end position.
	/// @note A wrapper around raycast(..., RaycastHit*). Retained for backward compatibility.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		startPos	A position within the start polygon representing
	///  							the start of the ray. [(x, y, z)]
	///  @param[in]		endPos		The position to cast the ray toward. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[out]	t			The hit parameter. (FLT_MAX if no wall hit.)
	///  @param[out]	hitNormal	The normal of the nearest wall hit. [(x, y, z)]
	///  @param[out]	path		The reference ids of the visited polygons. [opt]
	///  @param[out]	pathCount	The number of visited polygons. [opt]
	///  @param[in]		maxPath		The maximum number of polygons the @p path array can hold.
	/// @returns The status flags for the query.
	Raycast(startRef DtPolyRef, startPos, endPos []float32,
		filter *DtQueryFilter, t *float32, hitNormal []float32, path []DtPolyRef, pathCount *int32, maxPath int32) (status DtStatus)

	/// Casts a 'walkability' ray along the surface of the navigation mesh from
	/// the start position toward the end position.
	///  @param[in]		startRef	The reference id of the start polygon.
	///  @param[in]		startPos	A position within the start polygon representing
	///  							the start of the ray. [(x, y, z)]
	///  @param[in]		endPos		The position to cast the ray toward. [(x, y, z)]
	///  @param[in]		filter		The polygon filter to apply to the query.
	///  @param[in]		options		govern how the raycast behaves. See dtRaycastOptions
	///  @param[out]	hit			Pointer to a raycast hit structure which will be filled by the results.
	///  @param[in]		prevRef		parent of start ref. Used during for cost calculation [opt]
	/// @returns The status flags for the query.
	Raycast1(startRef DtPolyRef, startPos, endPos []float32,
		filter *DtQueryFilter, options int32, hit *dtRaycastHit, prevRef ...DtPolyRef) (status DtStatus)
	/// Finds the distance from the specified position to the nearest polygon wall.
	///  @param[in]		startRef		The reference id of the polygon containing @p centerPos.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		maxRadius		The radius of the search circle.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	hitDist			The distance to the nearest wall from @p centerPos.
	///  @param[out]	hitPos			The nearest position on the wall that was hit. [(x, y, z)]
	///  @param[out]	hitNormal		The normalized ray formed from the wall point to the
	///  								source point. [(x, y, z)]
	/// @returns The status flags for the query.
	FindDistanceToWall(startRef DtPolyRef, centerPos []float32, maxRadius float32,
		filter *DtQueryFilter, hitDist *float32, hitPos []float32, hitNormal []float32) (status DtStatus)
	/// Returns random location on navmesh.
	/// Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[in]		frand			Function returning a random number [0..1).
	///  @param[out]	randomRef		The reference id of the random location.
	///  @param[out]	randomPt		The random location.
	/// @returns The status flags for the query.
	FindRandomPoint(filter *DtQueryFilter) (randomRef DtPolyRef, randomPt []float32, status DtStatus)
	/// Returns the segments for the specified polygon, optionally including portals.
	///  @param[in]		ref				The reference id of the polygon.
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[out]	segmentVerts	The segments. [(ax, ay, az, bx, by, bz) * segmentCount]
	///  @param[out]	segmentRefs		The reference ids of each segment's neighbor polygon.
	///  								Or zero if the segment is a wall. [opt] [(parentRef) * @p segmentCount]
	///  @param[out]	segmentCount	The number of segments returned.
	///  @param[in]		maxSegments		The maximum number of segments the result arrays can hold.
	/// @returns The status flags for the query.
	GetPolyWallSegments(ref DtPolyRef, filter *DtQueryFilter,
		segmentVerts []float32, segmentRefs []DtPolyRef,
		maxSegments int32) (segmentCount int32, status DtStatus)

	/// Returns random location on navmesh within the reach of specified location.
	/// Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
	/// The location is not exactly constrained by the circle, but it limits the visited polygons.
	///  @param[in]		startRef		The reference id of the polygon where the search starts.
	///  @param[in]		centerPos		The center of the search circle. [(x, y, z)]
	///  @param[in]		maxRadius		The radius of the search circle. [Units: wu]
	///  @param[in]		filter			The polygon filter to apply to the query.
	///  @param[in]		frand			Function returning a random number [0..1).
	///  @param[out]	randomRef		The reference id of the random location.
	///  @param[out]	randomPt		The random location. [(x, y, z)]
	/// @returns The status flags for the query.
	FindRandomPointAroundCircle(startRef DtPolyRef, centerPos []float32, maxRadius float32,
		filter *DtQueryFilter) (randomRef DtPolyRef, randomPt []float32, status DtStatus)
	/// Finds the closest point on the specified polygon.
	///  @param[in]		ref			The reference id of the polygon.
	///  @param[in]		pos			The position to check. [(x, y, z)]
	///  @param[out]	closest		The closest point on the polygon. [(x, y, z)]
	///  @param[out]	posOverPoly	True of the position is over the polygon.
	/// @returns The status flags for the query.
	ClosestPointOnPoly(ref DtPolyRef, pos []float32, closest []float32, posOverPoly *bool) (status DtStatus)
	/// Returns a point on the boundary closest to the source point if the source point is outside the
	/// polygon's xz-bounds.
	///  @param[in]		ref			The reference id to the polygon.
	///  @param[in]		pos			The position to check. [(x, y, z)]
	///  @param[out]	closest		The closest point. [(x, y, z)]
	/// @returns The status flags for the query.
	ClosestPointOnPolyBoundary(ref DtPolyRef, pos []float32) (closest []float32, status DtStatus)

	/// Gets the height of the polygon at the provided position using the height detail. (Most accurate.)
	///  @param[in]		ref			The reference id of the polygon.
	///  @param[in]		pos			A position within the xz-bounds of the polygon. [(x, y, z)]
	///  @param[out]	height		The height at the surface of the polygon.
	/// @returns The status flags for the query.
	GetPolyHeight(ref DtPolyRef, pos []float32) (height float32, status DtStatus)
	/// @}
	/// @name Miscellaneous Functions
	/// @{

	/// Returns true if the polygon reference is valid and passes the filter restrictions.
	///  @param[in]		ref			The polygon reference to check.
	///  @param[in]		filter		The filter to apply.
	IsValidPolyRef(ref DtPolyRef, filter *DtQueryFilter) bool
	/// Returns true if the polygon reference is in the closed list.
	///  @param[in]		ref		The reference id of the polygon to check.
	/// @returns True if the polygon is in closed list.
	IsInClosedList(ref DtPolyRef) bool

	GetAttachedNavMesh() IDtNavMesh

	GetNodePool() *DtNodePool
}

type DtNavMeshQuery struct {
	m_nav          IDtNavMesh ///< Pointer to navmesh data.
	m_nodePool     *DtNodePool
	m_tinyNodePool *DtNodePool
	m_openList     NodeQueue[*DtNode]
	m_query        *dtQueryData
}

// / Initializes the query object.
// /  @param[in]		nav			Pointer to the DtNavMesh object to use for all queries.
// /  @param[in]		maxNodes	Maximum number of search nodes. [Limits: 0 < value <= 65535]
// / @returns The status flags for the query.
func NewDtNavMeshQuery(nav IDtNavMesh, maxNodes int32) NavMeshQuery {
	query := &DtNavMeshQuery{
		m_nodePool:     NewDtNodePool(maxNodes, int32(common.NextPow2(uint32(maxNodes/4)))),
		m_tinyNodePool: NewDtNodePool(64, 32),
		m_nav:          nav,
		m_openList: NewNodeQueue(func(t1, t2 *DtNode) bool {
			return t1.Total < t2.Total //花费最小
		})}

	return query
}

// / Gets the node pool.
// / @returns The node pool.
func (query *DtNavMeshQuery) GetNodePool() *DtNodePool { return query.m_nodePool }

// / Gets the navigation mesh the query object is using.
// / @return The navigation mesh the query object is using.
func (query *DtNavMeshQuery) GetAttachedNavMesh() IDtNavMesh { return query.m_nav }

// / Returns random location on navmesh.
// / Polygons are chosen weighted by area. The search runs in linear related to number of polygon.
// /  @param[in]		filter			The polygon filter to apply to the query.
// /  @param[in]		frand			Function returning a random number [0..1).
// /  @param[out]	randomRef		The reference id of the random location.
// /  @param[out]	randomPt		The random location.
// / @returns The status flags for the query.
func (query *DtNavMeshQuery) FindRandomPoint(filter *DtQueryFilter) (randomRef DtPolyRef, randomPt []float32, status DtStatus) {
	// Randomly pick one tile. Assume that all tiles cover roughly the same area.
	tsum := 0.0
	var tile *DtMeshTile
	for i := int32(0); i < query.m_nav.GetMaxTiles(); i++ {
		t := query.m_nav.GetTile(int(i))
		if t == nil || t.Header == nil {
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

	var poly *DtPoly
	var polyRef DtPolyRef
	// Randomly pick one polygon weighted by polygon area.
	base := query.m_nav.GetPolyRefBase(tile)

	areaSum := float32(0.0)
	for i := int32(0); i < tile.Header.PolyCount; i++ {
		p := tile.Polys[i]
		// Do not return off-mesh connection polygons.
		if p.GetType() != DT_POLYTYPE_GROUND {
			continue
		}

		// Must pass filter
		ref := base | DtPolyRef(i)
		if !filter.passFilter(p) {
			continue
		}

		// Calc area of the polygon.
		polyArea := float32(0.0)
		for j := uint8(2); j < p.VertCount; j++ {
			va := common.GetVert3(tile.Verts, p.Verts[0])
			vb := common.GetVert3(tile.Verts, p.Verts[j-1])
			vc := common.GetVert3(tile.Verts, p.Verts[j])
			polyArea += common.TriArea2D(va, vb, vc)
		}

		// Choose random polygon weighted by area, using reservoir sampling.
		areaSum += polyArea
		u := rand.Float32()
		if u*areaSum <= polyArea {
			poly = p
			polyRef = ref
		}
	}

	if poly == nil {
		return randomRef, randomPt, DT_FAILURE
	}

	// Randomly pick point on polygon.
	verts := make([]float32, 3*DT_VERTS_PER_POLYGON)
	areas := make([]float32, DT_VERTS_PER_POLYGON)
	copy(verts[0:3], common.GetVert3(tile.Verts, poly.Verts[0]))
	for j := uint8(1); j < poly.VertCount; j++ {
		copy(verts[j:3], common.GetVert3(tile.Verts, poly.Verts[j]))
	}

	s := rand.Float32()
	t := rand.Float32()

	pt := dtRandomPointInConvexPoly(verts, int32(poly.VertCount), areas, s, t)
	var tmp bool
	query.ClosestPointOnPoly(polyRef, pt, randomPt, &tmp)
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
func (query *DtNavMeshQuery) ClosestPointOnPoly(ref DtPolyRef, pos []float32, closest []float32, posOverPoly *bool) (status DtStatus) {
	if query.m_nav == nil {
		panic("")
	}
	if query.m_nav.IsValidPolyRef(ref) || len(pos) == 0 || !common.Visfinite(pos) {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	query.m_nav.ClosestPointOnPoly(ref, pos, closest, posOverPoly)
	return DT_SUCCESS
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
func (query *DtNavMeshQuery) ClosestPointOnPolyBoundary(ref DtPolyRef, pos []float32) (closest []float32, status DtStatus) {
	if query.m_nav == nil {
		panic("")
	}
	closest = make([]float32, 3)
	tile, poly, status := query.m_nav.GetTileAndPolyByRef(ref)
	if status.DtStatusFailed() {
		return closest, DT_FAILURE | DT_INVALID_PARAM
	}

	if len(pos) == 0 || !common.Visfinite(pos) {
		return closest, DT_FAILURE | DT_INVALID_PARAM
	}

	// Collect vertices.
	var verts [DT_VERTS_PER_POLYGON * 3]float32
	var edged [DT_VERTS_PER_POLYGON]float32
	var edget [DT_VERTS_PER_POLYGON]float32
	nv := 0
	for i := 0; i < int(poly.VertCount); i++ {
		copy(verts[nv*3:nv*3+3], common.GetVert3(tile.Verts, poly.Verts[i]))
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
		va := common.GetVert3(verts[:], imin)
		vb := common.GetVert3(verts[:], (imin+1)%nv)
		common.Vlerp(closest, va, vb, edget[imin])
	}

	return closest, DT_SUCCESS
}

type DtQueryFilter struct {
	m_areaCost     [DT_MAX_AREAS]float32 ///< Cost per area type. (Used by default implementation.)
	m_includeFlags uint16                ///< Flags for polygons that can be visited. (Used by default implementation.)
	m_excludeFlags uint16                ///< Flags for polygons that should not be visited. (Used by default implementation.)
}

/// @name Getters and setters for the default implementation data.
///@{

// / Returns the traversal cost of the area.
// /  @param[in]		i		The id of the area.
// / @returns The traversal cost of the area.
func (filter *DtQueryFilter) GetAreaCost(i int) float32 { return filter.m_areaCost[i] }

// / Sets the traversal cost of the area.
// /  @param[in]		i		The id of the area.
// /  @param[in]		cost	The new cost of traversing the area.
func (filter *DtQueryFilter) SetAreaCost(i int, cost float32) { filter.m_areaCost[i] = cost }

// / Returns the include flags for the filter.
// / Any polygons that include one or more of these flags will be
// / included in the operation.
func (filter *DtQueryFilter) GetIncludeFlags() uint16 { return filter.m_includeFlags }

// / Sets the include flags for the filter.
// / @param[in]		flags	The new flags.
func (filter *DtQueryFilter) SetIncludeFlags(flags uint16) { filter.m_includeFlags = flags }

// / Returns the exclude flags for the filter.
// / Any polygons that include one ore more of these flags will be
// / excluded from the operation.
func (filter *DtQueryFilter) GetExcludeFlags() uint16 { return filter.m_excludeFlags }

// / Sets the exclude flags for the filter.
// / @param[in]		flags		The new flags.
func (filter *DtQueryFilter) SetExcludeFlags(flags uint16) { filter.m_excludeFlags = flags }
func (filter *DtQueryFilter) getCost(pa, pb []float32, curPoly *DtPoly) float32 {
	return common.Vdist(pa, pb) * filter.m_areaCost[curPoly.GetArea()]
}
func (filter *DtQueryFilter) passFilter(poly *DtPoly) bool {
	return (poly.Flags&filter.m_includeFlags) != 0 && (poly.Flags&filter.m_excludeFlags) == 0
}

func (query *DtNavMeshQuery) FindRandomPointAroundCircle(startRef DtPolyRef, centerPos []float32, maxRadius float32,
	filter *DtQueryFilter) (randomRef DtPolyRef, randomPt []float32, status DtStatus) {
	if query.m_nav == nil {
		return
	}
	// Validate input
	if !query.m_nav.IsValidPolyRef(startRef) || len(centerPos) == 0 || !common.Visfinite(centerPos) || maxRadius < 0 || !common.IsFinite(maxRadius) || filter == nil {
		return randomRef, randomPt, DT_FAILURE | DT_INVALID_PARAM
	}
	query.m_openList.Reset()
	query.m_nodePool.Clear()
	_, startPoly := query.m_nav.GetTileAndPolyByRefUnsafe(startRef)
	if !filter.passFilter(startPoly) {
		return randomRef, randomPt, DT_FAILURE | DT_INVALID_PARAM
	}
	startNode := query.m_nodePool.GetNode(startRef)
	copy(startNode.Pos[:], centerPos)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = 0
	startNode.Id = startRef
	startNode.Flags = DT_NODE_OPEN
	query.m_openList.Offer(startNode)
	radiusSqr := common.Sqr(maxRadius)
	areaSum := float32(0.0)

	var randomTile *DtMeshTile
	var randomPoly *DtPoly
	var randomPolyRef DtPolyRef

	for !query.m_openList.Empty() {
		bestNode := query.m_openList.Poll()
		bestNode.Flags &= ^uint32(DT_NODE_OPEN)
		bestNode.Flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.Id
		bestTile, bestPoly := query.m_nav.GetTileAndPolyByRefUnsafe(bestRef)

		// Place random locations on on ground.
		if bestPoly.GetType() == DT_POLYTYPE_GROUND {
			// Calc area of the polygon.
			polyArea := float32(0.0)
			for j := 2; j < int(bestPoly.VertCount); j++ {
				va := common.GetVert3(bestTile.Verts, bestPoly.Verts[0])
				vb := common.GetVert3(bestTile.Verts, bestPoly.Verts[j-1])
				vc := common.GetVert3(bestTile.Verts, bestPoly.Verts[j])
				polyArea += common.TriArea2D(va, vb, vc)
			}
			// Choose random polygon weighted by area, using reservoir sampling.
			areaSum += polyArea
			u := rand.Float32()
			if u*areaSum <= polyArea {
				randomTile = bestTile
				randomPoly = bestPoly
				randomPolyRef = bestRef
			}
		}

		// Get parent poly and tile.
		var parentRef DtPolyRef
		if bestNode.Pidx != 0 {
			parentRef = query.m_nodePool.GetNodeAtIdx(int(bestNode.Pidx)).Id
		}

		if parentRef > 0 {
			_, _ = query.m_nav.GetTileAndPolyByRefUnsafe(parentRef)
		}

		for i := bestPoly.FirstLink; i != DT_NULL_LINK; i = bestTile.Links[i].Next {
			link := bestTile.Links[i]
			neighbourRef := link.Ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef == 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour
			neighbourTile, neighbourPoly := query.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.
			va := make([]float32, 3)
			vb := make([]float32, 3)
			status = query.getPortalPoints1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile, va, vb)
			if status.DtStatusFailed() {
				continue
			}

			// If the circle is not touching the next polygon, skip it.

			_, distSqr := DtDistancePtSegSqr2D(centerPos, va, vb)
			if distSqr > radiusSqr {
				continue
			}

			neighbourNode := query.m_nodePool.GetNode(neighbourRef)
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.Flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.Flags == 0 {
				common.Vlerp(neighbourNode.Pos[:], va, vb, 0.5)
			}

			total := bestNode.Total + common.Vdist(bestNode.Pos[:], neighbourNode.Pos[:])

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_OPEN > 0) && total >= neighbourNode.Total {
				continue
			}

			neighbourNode.Id = neighbourRef
			neighbourNode.Flags = neighbourNode.Flags & ^uint32(DT_NODE_CLOSED)
			neighbourNode.Pidx = uint32(query.m_nodePool.GetNodeIdx(bestNode))
			neighbourNode.Total = total

			if neighbourNode.Flags&DT_NODE_OPEN > 0 {
				query.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.Flags = DT_NODE_OPEN
				query.m_openList.Offer(neighbourNode)
			}
		}
	}

	if randomPoly == nil {
		return randomRef, randomPt, DT_FAILURE
	}

	// Randomly pick point on polygon.
	var verts [3 * DT_VERTS_PER_POLYGON]float32
	var areas [DT_VERTS_PER_POLYGON]float32
	copy(verts[0:3], common.GetVert3(randomTile.Verts, randomPoly.Verts[0]))
	for j := 1; j < int(randomPoly.VertCount); j++ {
		copy(verts[j*3:j*3+3], common.GetVert3(randomTile.Verts, randomPoly.Verts[j]))
	}

	s := rand.Float32()
	t := rand.Float32()

	pt := dtRandomPointInConvexPoly(verts[:], int32(randomPoly.VertCount), areas[:], s, t)
	var tmp bool
	query.ClosestPointOnPoly(randomPolyRef, pt, randomPt, &tmp)

	randomRef = randomPolyRef

	return randomRef, randomPt, status
}

func (query *DtNavMeshQuery) getPortalPoints(from DtPolyRef, to DtPolyRef, left, right []float32) (fromType uint8, toType uint8, status DtStatus) {
	if query.m_nav == nil {
		panic("")
	}

	fromTile, fromPoly, status := query.m_nav.GetTileAndPolyByRef(from)
	if status.DtStatusFailed() {
		return fromType, toType, DT_FAILURE | DT_INVALID_PARAM
	}

	fromType = fromPoly.GetType()

	toTile, toPoly, status := query.m_nav.GetTileAndPolyByRef(to)
	if status.DtStatusFailed() {
		return fromType, toType, DT_FAILURE | DT_INVALID_PARAM
	}

	toType = toPoly.GetType()
	status = query.getPortalPoints1(from, fromPoly, fromTile, to, toPoly, toTile, left, right)
	return fromType, toType, status
}

// Returns portal points between two polygons.
func (query *DtNavMeshQuery) getPortalPoints1(from DtPolyRef, fromPoly *DtPoly, fromTile *DtMeshTile,
	to DtPolyRef, toPoly *DtPoly, toTile *DtMeshTile, left, right []float32) (status DtStatus) {
	// Find the link that points to the 'to' polygon.
	// Find the link that points to the 'to' polygon.
	var link *DtLink
	for i := fromPoly.FirstLink; i != DT_NULL_LINK; i = fromTile.Links[i].Next {
		if fromTile.Links[i].Ref == to {
			link = fromTile.Links[i]
			break
		}
	}
	if link == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	// Handle off-mesh connections.
	if fromPoly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		// Find link that points to first vertex.
		for i := fromPoly.FirstLink; i != DT_NULL_LINK; i = fromTile.Links[i].Next {
			if fromTile.Links[i].Ref == to {
				v := fromTile.Links[i].Edge
				copy(left, common.GetVert3(fromTile.Verts, fromPoly.Verts[v]))
				copy(right, common.GetVert3(fromTile.Verts, fromPoly.Verts[v]))
				return DT_SUCCESS
			}
		}
		return DT_FAILURE | DT_INVALID_PARAM
	}

	if toPoly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		for i := toPoly.FirstLink; i != DT_NULL_LINK; i = toTile.Links[i].Next {
			if toTile.Links[i].Ref == from {
				v := toTile.Links[i].Edge
				copy(left, common.GetVert3(toTile.Verts, toPoly.Verts[v]))
				copy(right, common.GetVert3(toTile.Verts, toPoly.Verts[v]))
				return DT_SUCCESS
			}
		}
		return DT_FAILURE | DT_INVALID_PARAM
	}

	// Find portal vertices.
	v0 := fromPoly.Verts[link.Edge]
	v1 := fromPoly.Verts[(link.Edge+1)%fromPoly.VertCount]
	copy(left, common.GetVert3(fromTile.Verts, v0))
	copy(right, common.GetVert3(fromTile.Verts, v1))
	// If the link is at tile boundary, dtClamp the vertices to
	// the link width.
	if link.Side != 0xff {
		// Unpack portal limits.
		if link.Bmin != 0 || link.Bmax != 255 {
			s := float32(1.0 / 255.0)
			tmin := float32(link.Bmin) * s
			tmax := float32(link.Bmax) * s
			common.Vlerp(left, common.GetVert3(fromTile.Verts, v0), common.GetVert3(fromTile.Verts, v1), tmin)
			common.Vlerp(right, common.GetVert3(fromTile.Verts, v0), common.GetVert3(fromTile.Verts, v1), tmax)
		}
	}

	return DT_SUCCESS
}

func (q *DtNavMeshQuery) QueryPolygonsInTile(tile *DtMeshTile, qmin, qmax []float32, filter *DtQueryFilter, query dtPolyQuery) {
	if q.m_nav == nil {
		panic("")
	}
	batchSize := int32(32)
	var polyRefs [32]DtPolyRef
	var polys [32]*DtPoly
	n := int32(0)

	if tile.BvTree != nil {
		node := 0
		end := int(tile.Header.BvNodeCount)
		tbmin := tile.Header.Bmin
		tbmax := tile.Header.Bmax
		qfac := tile.Header.BvQuantFactor

		// Calculate quantized box
		var bmin [3]uint16
		var bmax [3]uint16
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
		base := q.m_nav.GetPolyRefBase(tile)
		for node < end {

			overlap := dtOverlapQuantBounds(bmin, bmax, tile.BvTree[node].Bmin, tile.BvTree[node].Bmax)
			isLeafNode := tile.BvTree[node].I >= 0

			if isLeafNode && overlap {
				ref := base | DtPolyRef(tile.BvTree[node].I)
				if filter.passFilter(tile.Polys[tile.BvTree[node].I]) {
					polyRefs[n] = ref
					polys[n] = tile.Polys[tile.BvTree[node].I]

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
				escapeIndex := -tile.BvTree[node].I
				node += int(escapeIndex)
			}
		}
	} else {
		bmin := make([]float32, 3)
		bmax := make([]float32, 3)
		base := q.m_nav.GetPolyRefBase(tile)
		for i := int32(0); i < tile.Header.PolyCount; i++ {
			p := tile.Polys[i]
			// Do not return off-mesh connection polygons.
			if p.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Must pass filter
			ref := base | DtPolyRef(i)
			if filter != nil && filter.passFilter(p) {
				continue
			}

			// Calc polygon bounds.
			v := common.GetVert3(tile.Verts, p.Verts[0])
			copy(bmin, v)
			copy(bmax, v)
			for j := 1; j < int(p.VertCount); j++ {
				v := common.GetVert3(tile.Verts, p.Verts[j])
				common.Vmin(bmin, v)
				common.Vmax(bmax, v)
			}
			if common.DtOverlapBounds(qmin, qmax, bmin, bmax) {
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
// / Used by DtNavMeshQuery::queryPolygons.
// / @ingroup detour
type dtPolyQuery interface {
	/// Called for each batch of unique polygons touched by the search area in DtNavMeshQuery::queryPolygons.
	/// This can be called multiple times for a single query.
	process(tile *DtMeshTile, refs []DtPolyRef, count int32)
}

type dtFindNearestPolyQuery struct {
	m_query              *DtNavMeshQuery
	m_center             []float32
	m_nearestDistanceSqr float32
	m_nearestRef         DtPolyRef
	m_nearestPoint       [3]float32
	m_overPoly           bool
}

func newDtFindNearestPolyQuery(query *DtNavMeshQuery, center []float32) *dtFindNearestPolyQuery {
	return &dtFindNearestPolyQuery{
		m_center: center,
		m_query:  query,
	}
}
func (query *dtFindNearestPolyQuery) nearestRef() DtPolyRef   { return query.m_nearestRef }
func (query *dtFindNearestPolyQuery) nearestPoint() []float32 { return query.m_nearestPoint[:] }
func (query *dtFindNearestPolyQuery) isOverPoly() bool        { return query.m_overPoly }
func (query *dtFindNearestPolyQuery) process(tile *DtMeshTile, refs []DtPolyRef, count int32) {

	for i := 0; i < int(count); i++ {
		ref := refs[i]
		var d float32
		closestPtPoly := make([]float32, 3)
		posOverPoly := false
		query.m_query.ClosestPointOnPoly(ref, query.m_center, closestPtPoly, &posOverPoly)

		// If a point is directly over a polygon and closer than
		// climb height, favor that instead of straight line nearest point.
		diff := make([]float32, 3)
		common.Vsub(diff, query.m_center, closestPtPoly)
		if posOverPoly {
			d = common.Abs(diff[1]) - tile.Header.WalkableClimb
			if d > 0 {
				d = d * d
			} else {
				d = 0
			}
		} else {
			d = common.VlenSqr(diff)
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
	m_polys        []DtPolyRef
	m_maxPolys     int32
	m_numCollected int32
	m_overflow     bool
}

func newDtCollectPolysQuery(polys []DtPolyRef, maxPolys int32) *dtCollectPolysQuery {
	return &dtCollectPolysQuery{
		m_polys:    polys,
		m_maxPolys: maxPolys,
	}
}
func (query *dtCollectPolysQuery) numCollected() int32 { return query.m_numCollected }
func (query *dtCollectPolysQuery) overflowed() bool    { return query.m_overflow }

func (query *dtCollectPolysQuery) process(tile *DtMeshTile, refs []DtPolyRef, count int32) {
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
	ref        DtPolyRef
	tmin, tmax int8
}

func insertInterval(ints []*dtSegInterval, maxInts int32, tmin, tmax int8, ref DtPolyRef) (nints int32) {
	if nints+1 > maxInts {
		return
	}
	// Find insertion point.
	idx := int32(0)
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
func (q *DtNavMeshQuery) FindDistanceToWall(startRef DtPolyRef, centerPos []float32, maxRadius float32,
	filter *DtQueryFilter, hitDist *float32, hitPos []float32, hitNormal []float32) (status DtStatus) {

	// Validate input
	if !q.m_nav.IsValidPolyRef(startRef) ||
		len(centerPos) == 0 || !common.Visfinite(centerPos) || maxRadius < 0 || !common.IsFinite(maxRadius) || filter == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	q.m_nodePool.Clear()
	q.m_openList.Reset()

	startNode := q.m_nodePool.GetNode(startRef)
	copy(startNode.Pos[:], centerPos)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = 0
	startNode.Id = startRef
	startNode.Flags = DT_NODE_OPEN
	q.m_openList.Offer(startNode)

	radiusSqr := common.Sqr(maxRadius)

	for !q.m_openList.Empty() {
		bestNode := q.m_openList.Poll()
		bestNode.Flags = bestNode.Flags & (^uint32(DT_NODE_OPEN))
		bestNode.Flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.Id
		var bestTile *DtMeshTile
		var bestPoly *DtPoly
		bestTile, bestPoly = q.m_nav.GetTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef DtPolyRef

		if bestNode.Pidx != 0 {
			parentRef = q.m_nodePool.GetNodeAtIdx(int(bestNode.Pidx)).Id
		}

		if parentRef != 0 {
			q.m_nav.GetTileAndPolyByRefUnsafe(parentRef)
		}

		i := uint8(0)
		j := bestPoly.VertCount - 1
		// Hit test walls.
		for i < bestPoly.VertCount {
			// Skip non-solid edges.
			if bestPoly.Neis[j]&DT_EXT_LINK > 0 {
				// Tile border.
				solid := true
				for k := bestPoly.FirstLink; k != DT_NULL_LINK; k = bestTile.Links[k].Next {
					link := bestTile.Links[k]
					if link.Edge == j {
						if link.Ref != 0 {
							_, neiPoly := q.m_nav.GetTileAndPolyByRefUnsafe(link.Ref)
							if filter.passFilter(neiPoly) {
								solid = false
							}

						}
						break
					}
				}
				if !solid {
					j = i
					i++
					continue
				}
			} else if bestPoly.Neis[j] > 0 {
				// Internal edge
				idx := (bestPoly.Neis[j] - 1)
				if filter.passFilter(bestTile.Polys[idx]) {
					j = i
					i++
					continue
				}

			}

			// Calc distance to the edge.
			vj := common.GetVert3(bestTile.Verts, bestPoly.Verts[j])
			vi := common.GetVert3(bestTile.Verts, bestPoly.Verts[j])

			tseg, distSqr := DtDistancePtSegSqr2D(centerPos, vj, vi)

			// Edge is too far, skip.
			if distSqr > radiusSqr {
				j = i
				i++
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

		for i := bestPoly.FirstLink; i != DT_NULL_LINK; i = bestTile.Links[i].Next {
			link := bestTile.Links[i]
			neighbourRef := link.Ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour.

			neighbourTile, neighbourPoly := q.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			// Skip off-mesh connections.
			if neighbourPoly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Calc distance to the edge.
			va := common.GetVert3(bestTile.Verts, bestPoly.Verts[link.Edge])
			vb := common.GetVert3(bestTile.Verts, (link.Edge+1)%bestPoly.VertCount)

			_, distSqr := DtDistancePtSegSqr2D(centerPos, va, vb)

			// If the circle is not touching the next polygon, skip it.
			if distSqr > radiusSqr {
				continue
			}

			if !filter.passFilter(neighbourPoly) {
				continue
			}

			neighbourNode := q.m_nodePool.GetNode(neighbourRef)
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.Flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.Flags == 0 {
				pos, _ := q.getEdgeMidPoint1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile)
				copy(neighbourNode.Pos[:], pos)
			}

			total := bestNode.Total + float32(common.Vdist(bestNode.Pos[:], neighbourNode.Pos[:]))

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_OPEN > 0) && total >= neighbourNode.Total {
				continue
			}

			neighbourNode.Id = neighbourRef
			neighbourNode.Flags = (neighbourNode.Flags & ^uint32(DT_NODE_CLOSED))
			neighbourNode.Pidx = uint32(q.m_nodePool.GetNodeIdx(bestNode))
			neighbourNode.Total = total

			if neighbourNode.Flags&DT_NODE_OPEN > 0 {
				q.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.Flags |= DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}
		}
	}

	// Calc hit normal.
	common.Vsub(hitNormal, centerPos, hitPos)
	common.Vnormalize(hitNormal)

	*hitDist = float32(math.Sqrt(float64(radiusSqr)))

	return status
}

// / @par
// /
// / The closed list is the list of polygons that were fully evaluated during
// / the last navigation graph search. (A* or Dijkstra)
// /
func (q *DtNavMeshQuery) IsInClosedList(ref DtPolyRef) bool {
	if q.m_nodePool == nil {
		return false
	}

	nodes, n := q.m_nodePool.FindNodes(ref, DT_MAX_STATES_PER_NODE)
	for i := 0; i < n; i++ {
		if nodes[i].Flags&DT_NODE_CLOSED != 0 {
			return true
		}

	}
	return false
}

// Returns edge mid point between two polygons.
func (q *DtNavMeshQuery) getEdgeMidPoint(from, to DtPolyRef) (mid []float32, status DtStatus) {
	mid = make([]float32, 3)
	left := make([]float32, 3)
	right := make([]float32, 3)
	_, _, status = q.getPortalPoints(from, to, left, right)
	if status.DtStatusFailed() {
		return mid, DT_FAILURE | DT_INVALID_PARAM
	}

	mid[0] = (left[0] + right[0]) * 0.5
	mid[1] = (left[1] + right[1]) * 0.5
	mid[2] = (left[2] + right[2]) * 0.5
	return mid, DT_SUCCESS
}

func (q *DtNavMeshQuery) getEdgeMidPoint1(from DtPolyRef, fromPoly *DtPoly, fromTile *DtMeshTile,
	to DtPolyRef, toPoly *DtPoly, toTile *DtMeshTile) (mid []float32, status DtStatus) {
	mid = make([]float32, 3)
	left := make([]float32, 3)
	right := make([]float32, 3)
	status = q.getPortalPoints1(from, fromPoly, fromTile, to, toPoly, toTile, left, right)
	if status.DtStatusFailed() {
		return mid, DT_FAILURE | DT_INVALID_PARAM
	}

	mid[0] = (left[0] + right[0]) * 0.5
	mid[1] = (left[1] + right[1]) * 0.5
	mid[2] = (left[2] + right[2]) * 0.5
	return mid, DT_SUCCESS
}

func (q *DtNavMeshQuery) appendVertex(pos []float32, flags int32, ref DtPolyRef,
	straightPath []float32, straightPathFlags []int32, straightPathRefs []DtPolyRef, straightPathCount, maxStraightPath int32) (int32, DtStatus) {
	if (straightPathCount) > 0 && common.Vequal(straightPath[((straightPathCount)-1)*3:((straightPathCount)-1)*3+3], pos) {
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
func (q *DtNavMeshQuery) FindStraightPath(startPos, endPos []float32,
	path []DtPolyRef, pathSize int32,
	straightPath []float32, straightPathFlags []int32, straightPathRefs []DtPolyRef,
	maxStraightPath int32, optionss ...int32) (straightPathCount int32, status DtStatus) {
	if q.m_nav == nil {
		panic("")
	}
	options := int32(0)
	if len(optionss) > 0 {
		options = optionss[0]
	}
	if len(startPos) == 0 || !common.Visfinite(startPos) || //TODO
		len(endPos) == 0 || !common.Visfinite(endPos) ||
		len(path) == 0 || pathSize <= 0 || path[0] == 0 ||
		maxStraightPath <= 0 {
		return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	// TODO: Should this be callers responsibility?

	closestStartPos, _ := q.ClosestPointOnPolyBoundary(path[0], startPos)
	if status.DtStatusFailed() {
		return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	closestEndPos, status := q.ClosestPointOnPolyBoundary(path[pathSize-1], endPos)
	if status.DtStatusFailed() {
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
		var portalApex, portalLeft, portalRight [3]float32
		copy(portalApex[:], closestStartPos)
		copy(portalLeft[:], portalApex[:])
		copy(portalRight[:], portalApex[:])
		var apexIndex = int32(0)
		var leftIndex = int32(0)
		var rightIndex = int32(0)

		var leftPolyType = uint8(0)
		var rightPolyType = uint8(0)

		var leftPolyRef DtPolyRef = path[0]
		var rightPolyRef DtPolyRef = path[0]

		for i := int32(0); int32(i) < pathSize; i++ {
			left := make([]float32, 3)
			right := make([]float32, 3)
			var toType uint8
			if i+1 < pathSize {

				_, toType, status = q.getPortalPoints(path[i], path[i+1], left, right)
				// Next portal.
				if status.DtStatusFailed() {
					// Failed to get portal points, in practice this means that path[i+1] is invalid polygon.
					// Clamp the end point to path[i], and return the path so far.
					closestEndPos, status = q.ClosestPointOnPolyBoundary(path[i], endPos)
					if status.DtStatusFailed() {
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
					return straightPathCount, DT_SUCCESS | DT_PARTIAL_RESULT | DtStatus(tmp)
				}

				// If starting really close the portal, advance.
				if i == 0 {
					t, _ := DtDistancePtSegSqr2D(portalApex[:], left, right)
					if float64(t) < common.Sqr(0.001) {
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
			if common.TriArea2D(portalApex[:], portalRight[:], right) <= 0.0 {
				if common.Vequal(portalApex[:], portalRight[:]) || common.TriArea2D(portalApex[:], portalLeft[:], right) > 0.0 {
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

					var flags = int32(0)
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
			if common.TriArea2D(portalApex[:], portalLeft[:], left) >= 0.0 {
				if common.Vequal(portalApex[:], portalLeft[:]) || common.TriArea2D(portalApex[:], portalRight[:], left) < 0.0 {
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

					var flags = int32(0)
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
	return straightPathCount, DT_SUCCESS | DtStatus(tmp)
}

func (q *DtNavMeshQuery) appendPortals(startIdx int32, endIdx int32, endPos []float32, path []DtPolyRef,
	straightPath []float32, straightPathFlags []int32, straightPathRefs []DtPolyRef,
	straightPathCount int32, maxStraightPath int32, options int32) (int32, DtStatus) {
	startPos := straightPath[(straightPathCount-1)*3 : (straightPathCount-1)*3+3]
	// Append or update last vertex
	for i := startIdx; i < endIdx; i++ {
		// Calculate portal
		from := path[i]
		fromTile, fromPoly, status := q.m_nav.GetTileAndPolyByRef(from)
		if status.DtStatusFailed() {
			return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
		}

		to := path[i+1]
		toTile, toPoly, status := q.m_nav.GetTileAndPolyByRef(to)
		if status.DtStatusFailed() {
			return straightPathCount, DT_FAILURE | DT_INVALID_PARAM
		}

		left := make([]float32, 3)
		right := make([]float32, 3)
		status = q.getPortalPoints1(from, fromPoly, fromTile, to, toPoly, toTile, left, right)
		if status.DtStatusFailed() {
			break
		}

		if options&DT_STRAIGHTPATH_AREA_CROSSINGS > 0 {
			// Skip intersection if only area crossings are requested.
			if fromPoly.GetArea() == toPoly.GetArea() {
				continue
			}

		}

		// Append intersection
		_, t, ok := dtIntersectSegSeg2D(startPos, endPos, left, right)
		if ok {
			pt := make([]float32, 3)
			common.Vlerp(pt, left, right, t)

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
func (q *DtNavMeshQuery) MoveAlongSurface(startRef DtPolyRef, startPos, endPos []float32,
	filter *DtQueryFilter, resultPos []float32, visited []DtPolyRef, visitedCount *int32, maxVisitedSize int32) (status DtStatus) {
	if !q.m_nav.IsValidPolyRef(startRef) ||
		len(startPos) == 0 || !common.Visfinite(startPos) ||
		len(endPos) == 0 || !common.Visfinite(endPos) ||
		filter == nil ||
		maxVisitedSize <= 0 {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	visited = make([]DtPolyRef, maxVisitedSize)
	resultPos = make([]float32, 3)
	status = DT_SUCCESS
	MAX_STACK := 48
	var stack [48]*DtNode
	nstack := 0

	q.m_tinyNodePool.Clear()
	startNode := q.m_tinyNodePool.GetNode(startRef)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = 0
	startNode.Id = startRef
	startNode.Flags = DT_NODE_CLOSED
	stack[nstack] = startNode
	nstack++

	bestDist := float32(math.MaxFloat32)
	bestPos := make([]float32, 3)
	copy(bestPos, startPos)
	var bestNode *DtNode
	// Search constraints
	searchPos := make([]float32, 3)
	common.Vlerp(searchPos, startPos, endPos, 0.5)
	searchRadSqr := common.Sqr(common.Vdist(startPos, endPos)/2.0 + 0.001)

	var verts [DT_VERTS_PER_POLYGON * 3]float32

	for nstack > 0 {
		// Pop front.
		curNode := stack[0]
		for i := 0; i < nstack-1; i++ {
			stack[i] = stack[i+1]
		}

		nstack--

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		curRef := curNode.Id
		curTile, curPoly := q.m_nav.GetTileAndPolyByRefUnsafe(curRef)

		// Collect vertices.
		nverts := int32(curPoly.VertCount)
		for i := int32(0); i < nverts; i++ {
			copy(verts[i*3:i*3+3], common.GetVert3(curTile.Verts, curPoly.Verts[i]))
		}

		// If target is inside the poly, stop search.
		if dtPointInPolygon(endPos, verts[:], nverts) {
			bestNode = curNode
			copy(bestPos, endPos)
			break
		}
		i := 0
		j := int(curPoly.VertCount - 1)
		// Find wall edges and find nearest point inside the walls.
		for i < int(curPoly.VertCount) {
			// Find links to neighbours.
			MAX_NEIS := 8
			nneis := 0
			var neis [8]DtPolyRef

			if curPoly.Neis[j]&DT_EXT_LINK > 0 {
				// Tile border.
				for k := curPoly.FirstLink; k != DT_NULL_LINK; k = curTile.Links[k].Next {
					link := curTile.Links[k]
					if int(link.Edge) == j {
						if link.Ref != 0 {
							_, neiPoly := q.m_nav.GetTileAndPolyByRefUnsafe(link.Ref)
							if filter.passFilter(neiPoly) {
								if nneis < MAX_NEIS {
									neis[nneis] = link.Ref
									nneis++
								}

							}
						}
					}

				}
			} else if curPoly.Neis[j] > 0 {
				idx := (curPoly.Neis[j] - 1)
				ref := q.m_nav.GetPolyRefBase(curTile) | DtPolyRef(idx)
				if filter.passFilter(curTile.Polys[idx]) {
					// Internal edge, encode id.
					neis[nneis] = ref
					nneis++
				}
			}

			if nneis == 0 {
				// Wall edge, calc distance.
				vj := common.GetVert3(verts[:], j)
				vi := common.GetVert3(verts[:], i)
				tseg, distSqr := DtDistancePtSegSqr2D(endPos, vj, vi)
				if distSqr < bestDist {
					// Update nearest distance.
					common.Vlerp(bestPos, vj, vi, tseg)
					bestDist = distSqr
					bestNode = curNode
				}
			} else {
				for k := 0; k < nneis; k++ {
					// Skip if no node can be allocated.
					neighbourNode := q.m_tinyNodePool.GetNode(neis[k])
					if neighbourNode == nil {
						j = i
						i++
						continue
					}

					// Skip if already visited.
					if neighbourNode.Flags&DT_NODE_CLOSED > 0 {
						j = i
						i++
						continue
					}

					// Skip the link if it is too far from search constraint.
					// TODO: Maybe should use getPortalPoints(), but this one is way faster.
					vj := common.GetVert3(verts[:], j)
					vi := common.GetVert3(verts[:], i)
					_, distSqr := DtDistancePtSegSqr2D(searchPos, vj, vi)
					if distSqr > searchRadSqr {
						j = i
						i++
						continue
					}

					// Mark as the node as visited and push to queue.
					if nstack < MAX_STACK {
						neighbourNode.Pidx = uint32(q.m_tinyNodePool.GetNodeIdx(curNode))
						neighbourNode.Flags |= DT_NODE_CLOSED
						stack[nstack] = neighbourNode
						nstack++
					}

				}
			}

			j = i
			i++
		}
	}

	n := int32(0)
	if bestNode != nil {
		// Reverse the path.
		var prev *DtNode
		node := bestNode
		common.DoWhile(func() (stop bool) {
			next := q.m_tinyNodePool.GetNodeAtIdx(int(node.Pidx))
			node.Pidx = uint32(q.m_tinyNodePool.GetNodeIdx(prev))
			prev = node
			node = next
			return false
		}, func() bool {
			return node != nil
		})
		// Store result
		node = prev
		common.DoWhile(func() (stop bool) {
			visited[n] = node.Id
			n++
			if n >= maxVisitedSize {
				status |= DT_BUFFER_TOO_SMALL
				return true
			}
			node = q.m_tinyNodePool.GetNodeAtIdx(int(node.Pidx))
			return false
		}, func() bool {
			return node != nil
		})
	}

	copy(resultPos, bestPos)

	*visitedCount = n

	return status
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

func (q *DtNavMeshQuery) FindPolysAroundCircle(startRef DtPolyRef, centerPos []float32, radius float32,
	filter *DtQueryFilter, resultRef []DtPolyRef, resultParent []DtPolyRef, resultCost []float32, resultCount *int32, maxResult int32) (status DtStatus) {

	if q.m_nav.IsValidPolyRef(startRef) ||
		len(centerPos) == 0 || !common.Visfinite(centerPos) ||
		radius < 0 || !common.IsFinite(radius) ||
		filter == nil || maxResult < 0 {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	q.m_nodePool.Clear()
	q.m_openList.Reset()
	startNode := q.m_nodePool.GetNode(startRef)
	copy(startNode.Pos[:], centerPos)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = 0
	startNode.Id = startRef
	startNode.Flags = DT_NODE_OPEN
	q.m_openList.Offer(startNode)

	status = DT_SUCCESS
	n := int32(0)
	radiusSqr := common.Sqr(radius)

	for !q.m_openList.Empty() {
		bestNode := q.m_openList.Poll()
		bestNode.Flags &= ^uint32(DT_NODE_OPEN)
		bestNode.Flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.Id
		bestTile, bestPoly := q.m_nav.GetTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef DtPolyRef
		if bestNode.Pidx != 0 {
			parentRef = q.m_nodePool.GetNodeAtIdx(int(bestNode.Pidx)).Id
		}

		if parentRef > 0 {
			q.m_nav.GetTileAndPolyByRefUnsafe(parentRef)
		}

		if n < maxResult {
			if len(resultRef) > 0 {
				resultRef[n] = bestRef
			}

			if len(resultParent) > 0 {
				resultParent[n] = parentRef
			}

			if len(resultCost) > 0 {
				resultCost[n] = bestNode.Total
			}

			n++
		} else {
			status |= DT_BUFFER_TOO_SMALL
		}

		for i := bestPoly.FirstLink; i != DT_NULL_LINK; i = bestTile.Links[i].Next {
			link := bestTile.Links[i]
			neighbourRef := link.Ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour

			neighbourTile, neighbourPoly := q.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.
			va := make([]float32, 3)
			vb := make([]float32, 3)
			status = q.getPortalPoints1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile, va, vb)
			if status.DtStatusFailed() {
				continue
			}

			// If the circle is not touching the next polygon, skip it.

			_, distSqr := DtDistancePtSegSqr2D(centerPos, va, vb)
			if distSqr > radiusSqr {
				continue
			}

			neighbourNode := q.m_nodePool.GetNode(neighbourRef)
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.Flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.Flags == 0 {
				common.Vlerp(neighbourNode.Pos[:], va, vb, 0.5)
			}

			cost := filter.getCost(bestNode.Pos[:], neighbourNode.Pos[:], bestPoly)

			total := bestNode.Total + cost

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_OPEN > 0) && total >= neighbourNode.Total {
				continue
			}

			neighbourNode.Id = neighbourRef
			neighbourNode.Pidx = uint32(q.m_nodePool.GetNodeIdx(bestNode))
			neighbourNode.Total = total

			if neighbourNode.Flags&DT_NODE_OPEN > 0 {
				q.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.Flags = DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}
		}
	}

	*resultCount = n

	return status
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
func (q *DtNavMeshQuery) FindPolysAroundShape(startRef DtPolyRef, verts []float32, nverts int32,
	filter *DtQueryFilter,
	resultRef []DtPolyRef, resultParent []DtPolyRef, resultCost []float32, maxResult int32) (resultCount int32, status DtStatus) {
	if !q.m_nav.IsValidPolyRef(startRef) ||
		len(verts) == 0 || nverts < 3 ||
		filter == nil || maxResult < 0 {
		return resultCount, DT_FAILURE | DT_INVALID_PARAM
	}

	// Validate input
	if startRef == 0 || !q.m_nav.IsValidPolyRef(startRef) {
		return resultCount, DT_FAILURE | DT_INVALID_PARAM
	}

	q.m_nodePool.Clear()
	q.m_openList.Reset()

	var centerPos = []float32{0, 0, 0}
	for i := int32(0); i < nverts; i++ {
		common.Vadd(centerPos, centerPos, common.GetVert3(verts, i))

	}

	common.Vscale(centerPos, centerPos, 1.0/float32(nverts))

	startNode := q.m_nodePool.GetNode(startRef)
	copy(startNode.Pos[:], centerPos)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = 0
	startNode.Id = startRef
	startNode.Flags = DT_NODE_OPEN
	q.m_openList.Offer(startNode)

	status = DT_SUCCESS

	n := int32(0)

	for !q.m_openList.Empty() {
		bestNode := q.m_openList.Poll()
		bestNode.Flags &= ^uint32(DT_NODE_OPEN)
		bestNode.Flags |= DT_NODE_CLOSED

		// Get poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.Id

		bestTile, bestPoly := q.m_nav.GetTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef DtPolyRef

		if bestNode.Pidx != 0 {
			parentRef = q.m_nodePool.GetNodeAtIdx(int(bestNode.Pidx)).Id
		}

		if parentRef != 0 {
			q.m_nav.GetTileAndPolyByRefUnsafe(parentRef)
		}

		if n < maxResult {
			if len(resultRef) > 0 {
				resultRef[n] = bestRef
			}

			if len(resultParent) > 0 {
				resultParent[n] = parentRef
			}

			if len(resultCost) > 0 {
				resultCost[n] = bestNode.Total
			}

			n++
		} else {
			status |= DT_BUFFER_TOO_SMALL
		}

		for i := bestPoly.FirstLink; i != DT_NULL_LINK; i = bestTile.Links[i].Next {
			link := bestTile.Links[i]
			neighbourRef := link.Ref
			// Skip invalid neighbours and do not follow back to parent.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Expand to neighbour

			neighbourTile, neighbourPoly := q.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// Find edge and calc distance to the edge.
			va := make([]float32, 3)
			vb := make([]float32, 3)
			status = q.getPortalPoints1(bestRef, bestPoly, bestTile, neighbourRef, neighbourPoly, neighbourTile, va, vb)
			if status.DtStatusFailed() {
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

			neighbourNode := q.m_nodePool.GetNode(neighbourRef)
			if neighbourNode == nil {
				status |= DT_OUT_OF_NODES
				continue
			}

			if neighbourNode.Flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Cost
			if neighbourNode.Flags == 0 {
				common.Vlerp(neighbourNode.Pos[:], va, vb, 0.5)

			}

			cost := filter.getCost(bestNode.Pos[:], neighbourNode.Pos[:], bestPoly)

			total := bestNode.Total + cost

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_OPEN > 0) && total >= neighbourNode.Total {
				continue
			}

			neighbourNode.Id = neighbourRef
			neighbourNode.Pidx = uint32(q.m_nodePool.GetNodeIdx(bestNode))
			neighbourNode.Total = total

			if neighbourNode.Flags&DT_NODE_OPEN > 0 {
				q.m_openList.Update(neighbourNode)
			} else {
				neighbourNode.Flags = DT_NODE_OPEN
				q.m_openList.Offer(neighbourNode)
			}
		}
	}

	resultCount = n

	return resultCount, status
}

func (q *DtNavMeshQuery) GetPathFromDijkstraSearch(endRef DtPolyRef, path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus) {
	if !q.m_nav.IsValidPolyRef(endRef) || len(path) == 0 || maxPath < 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}
	nodes, n := q.m_nodePool.FindNodes(endRef, 1)
	if n != 1 || (nodes[0].Flags&DT_NODE_CLOSED) == 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}
	endNode := nodes[0]
	pathCount, status = q.getPathToNode(endNode, path, maxPath)
	return
}

func (q *DtNavMeshQuery) getPathToNode(endNode *DtNode, path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus) {
	// Find the length of the entire path.
	curNode := endNode
	var length = 0
	common.DoWhile(func() (stop bool) {
		length++
		curNode = q.m_nodePool.GetNodeAtIdx(int(curNode.Pidx))
		return false
	}, func() bool {
		return curNode != nil
	})
	// If the path cannot be fully stored then advance to the last node we will be able to store.
	curNode = endNode
	var writeCount int
	for writeCount = length; writeCount > int(maxPath); writeCount-- {
		curNode = q.m_nodePool.GetNodeAtIdx(int(curNode.Pidx))
	}

	// Write path
	for i := writeCount - 1; i >= 0; i-- {
		path[i] = curNode.Id
		curNode = q.m_nodePool.GetNodeAtIdx(int(curNode.Pidx))
	}

	dtAssertTrue(curNode != nil)

	pathCount = common.Min(int32(length), maxPath)

	if length > int(maxPath) {
		return pathCount, DT_SUCCESS | DT_BUFFER_TOO_SMALL
	}

	return pathCount, DT_SUCCESS
}

func (q *DtNavMeshQuery) FinalizeSlicedFindPath(path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus) {

	if len(path) == 0 || maxPath <= 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	if q.m_query.status.DtStatusFailed() {
		// Reset query.
		q.m_query = nil
		return pathCount, DT_FAILURE
	}

	n := int32(0)

	if q.m_query.startRef == q.m_query.endRef {
		// Special case: the search starts and ends at same poly.
		path[n] = q.m_query.startRef
		n++
	} else {
		// Reverse the path.
		if q.m_query.lastBestNode.Id != q.m_query.endRef {
			q.m_query.status |= DT_PARTIAL_RESULT
		}

		var prev *DtNode
		node := q.m_query.lastBestNode
		prevRay := uint32(0)
		common.DoWhile(func() bool {
			next := q.m_nodePool.GetNodeAtIdx(int(node.Pidx))
			node.Pidx = uint32(q.m_nodePool.GetNodeIdx(prev))
			prev = node
			nextRay := node.Flags & DT_NODE_PARENT_DETACHED                        // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)
			node.Flags = (node.Flags & ^uint32(DT_NODE_PARENT_DETACHED)) | prevRay // and store it in the reversed path's node
			prevRay = nextRay
			node = next
			return false
		}, func() bool {
			return node != nil
		})

		// Store path
		node = prev

		common.DoWhile(func() bool {
			next := q.m_nodePool.GetNodeAtIdx(int(node.Pidx))
			status = 0
			if node.Flags&DT_NODE_PARENT_DETACHED > 0 {
				var (
					normal [3]float32
					t      float32
					m      int32
				)
				q.Raycast(node.Id, node.Pos[:], next.Pos[:], q.m_query.filter, &t, normal[:], path[n:], &m, maxPath-n)
				n += m
				// raycast ends on poly boundary and the path might include the next poly boundary.
				if path[n-1] == next.Id {
					n-- // remove to avoid duplicates}

				} else {
					path[n] = node.Id
					n++
					if n >= maxPath {
						status = DT_BUFFER_TOO_SMALL
					}

				}

				if status&DT_STATUS_DETAIL_MASK > 0 {
					q.m_query.status |= status & DT_STATUS_DETAIL_MASK
					return true
				}
				node = next

			}
			return false
		}, func() bool {
			return node != nil
		})
	}
	details := q.m_query.status & DT_STATUS_DETAIL_MASK

	// Reset query.
	q.m_query = &dtQueryData{}

	pathCount = n

	return pathCount, DT_SUCCESS | details
}

// / Provides information about raycast hit
// / filled by DtNavMeshQuery::raycast
// / @ingroup detour
type dtRaycastHit struct {

	/// The hit parameter. (FLT_MAX if no wall hit.)
	t float32

	/// hitNormal	The normal of the nearest wall hit. [(x, y, z)]
	hitNormal [3]float32

	/// The index of the edge on the final polygon where the wall was hit.
	hitEdgeIndex int32

	/// Pointer to an array of reference ids of the visited polygons. [opt]
	path []DtPolyRef

	/// The number of visited polygons. [opt]
	pathCount int32

	/// The maximum number of polygons the @p path array can hold.
	maxPath int32

	///  The cost of the path until hit.
	pathCost float32
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
func (q *DtNavMeshQuery) Raycast(startRef DtPolyRef, startPos, endPos []float32,
	filter *DtQueryFilter, t *float32, hitNormal []float32, path []DtPolyRef, pathCount *int32, maxPath int32) (status DtStatus) {
	var hit dtRaycastHit
	hit.path = path
	hit.maxPath = maxPath

	status = q.Raycast1(startRef, startPos, endPos, filter, 0, &hit)

	*t = hit.t
	if len(hitNormal) > 0 {
		copy(hitNormal, hit.hitNormal[:])
	}
	*pathCount = hit.pathCount
	return status
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
func (q *DtNavMeshQuery) Raycast1(startRef DtPolyRef, startPos, endPos []float32,
	filter *DtQueryFilter, options int32, hit *dtRaycastHit, prevRefs ...DtPolyRef) (status DtStatus) {
	var prevRef DtPolyRef
	if len(prevRefs) > 0 {
		prevRef = prevRefs[0]
	}
	if hit == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	hit.t = 0
	hit.pathCount = 0
	hit.pathCost = 0

	// Validate input
	if !q.m_nav.IsValidPolyRef(startRef) || len(startPos) == 0 ||
		len(endPos) == 0 || filter == nil || (prevRef == 0 && !q.m_nav.IsValidPolyRef(prevRef)) {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	var curPos, lastPos [3]float32
	var verts [DT_VERTS_PER_POLYGON*3 + 3]float32
	n := int32(0)
	dir := make([]float32, 3)
	copy(curPos[:], startPos)
	common.Vsub(dir, endPos[:], startPos)
	common.Vset(hit.hitNormal[:], 0, 0, 0)

	status = DT_SUCCESS

	var prevTile, tile, nextTile *DtMeshTile
	var prevPoly, poly, nextPoly *DtPoly
	var curRef DtPolyRef

	// The API input has been checked already, skip checking internal data.
	curRef = startRef
	tile = nil
	poly = nil
	tile, poly = q.m_nav.GetTileAndPolyByRefUnsafe(curRef)
	prevTile = tile
	prevPoly = poly
	nextTile = prevTile
	nextPoly = prevPoly
	if prevRef > 0 {
		prevTile, prevPoly = q.m_nav.GetTileAndPolyByRefUnsafe(prevRef)
	}

	for curRef != 0 {
		// Cast ray against current polygon.

		// Collect vertices.
		nv := int32(0)
		for i := uint8(0); i < poly.VertCount; i++ {
			copy(verts[nv*3:nv*3+3], common.GetVert3(tile.Verts, poly.Verts[i]))
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
			hit.t = float32(math.MaxFloat32)
			hit.pathCount = n

			// add the cost
			if options&DT_RAYCAST_USE_COSTS > 0 {
				hit.pathCost += filter.getCost(curPos[:], endPos, poly)
			}

			return status
		}

		// Follow neighbours.
		var nextRef DtPolyRef

		for i := poly.FirstLink; i != DT_NULL_LINK; i = tile.Links[i].Next {
			link := tile.Links[i]

			// Find link which contains this edge.
			if int32(link.Edge) != segMax {
				continue
			}

			// Get pointer to the next polygon.

			nextTile, nextPoly = q.m_nav.GetTileAndPolyByRefUnsafe(link.Ref)

			// Skip off-mesh connections.
			if nextPoly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Skip links based on filter.
			if !filter.passFilter(nextPoly) {
				continue
			}

			// If the link is internal, just return the ref.
			if link.Side == 0xff {
				nextRef = link.Ref
				break
			}

			// If the link is at tile boundary,

			// Check if the link spans the whole edge, and accept.
			if link.Bmin == 0 && link.Bmax == 255 {
				nextRef = link.Ref
				break
			}

			// Check for partial edge links.
			v0 := poly.Verts[link.Edge]
			v1 := poly.Verts[(link.Edge+1)%poly.VertCount]
			left := common.GetVert3(tile.Verts, v0)
			right := common.GetVert3(tile.Verts, v1)

			// Check that the intersection lies inside the link portal.
			if link.Side == 0 || link.Side == 4 {
				// Calculate link size.
				s := float32(1.0 / 255.0)
				lmin := left[2] + (right[2]-left[2])*(float32(link.Bmin)*s)
				lmax := left[2] + (right[2]-left[2])*(float32(link.Bmax)*s)
				if lmin > lmax {
					lmin, lmax = lmax, lmin
				}

				// Find Z intersection.
				z := startPos[2] + (endPos[2]-startPos[2])*tmax
				if z >= lmin && z <= lmax {
					nextRef = link.Ref
					break
				}
			} else if link.Side == 2 || link.Side == 6 {
				// Calculate link size.
				s := float32(1.0 / 255.0)
				lmin := left[0] + (right[0]-left[0])*(float32(link.Bmin)*s)
				lmax := left[0] + (right[0]-left[0])*(float32(link.Bmax)*s)
				if lmin > lmax {
					lmin, lmax = lmax, lmin
				}

				// Find X intersection.
				x := startPos[0] + (endPos[0]-startPos[0])*tmax
				if x >= lmin && x <= lmax {
					nextRef = link.Ref
					break
				}
			}
		}

		// add the cost
		if options&DT_RAYCAST_USE_COSTS > 0 {
			// compute the intersection point at the furthest end of the polygon
			// and correct the height (since the raycast moves in 2d)
			copy(lastPos[:], curPos[:])
			common.Vmad(curPos[:], startPos, dir, hit.t)
			e1 := common.GetVert3(verts[:], segMax)
			e2 := common.GetVert3(verts[:], (segMax+1)%nv)
			eDir := make([]float32, 3)
			diff := make([]float32, 3)
			common.Vsub(eDir, e2, e1)
			common.Vsub(diff, curPos[:], e1)
			s := diff[2] / eDir[2]
			if common.Sqr(eDir[0]) > common.Sqr(eDir[2]) {
				s = diff[0] / eDir[0]
			}
			curPos[1] = e1[1] + eDir[1]*s

			hit.pathCost += filter.getCost(lastPos[:], curPos[:], poly)
		}

		if nextRef == 0 {
			// No neighbour, we hit a wall.

			// Calculate hit normal.
			a := segMax
			b := int32(0)
			if segMax+1 < nv {
				b = segMax + 1
			}
			va := common.GetVert3(verts[:], a)
			vb := common.GetVert3(verts[:], b)
			dx := vb[0] - va[0]
			dz := vb[2] - va[2]
			hit.hitNormal[0] = dz
			hit.hitNormal[1] = 0
			hit.hitNormal[2] = -dx
			common.Vnormalize(hit.hitNormal[:])

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
func (q *DtNavMeshQuery) InitSlicedFindPath(startRef, endRef DtPolyRef,
	startPos, endPos []float32, filter *DtQueryFilter, options int32) DtStatus {
	// Init path state.
	q.m_query = &dtQueryData{}
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
	q.m_query.raycastLimitSqr = math.MaxFloat32

	// Validate input
	if !q.m_nav.IsValidPolyRef(startRef) || !q.m_nav.IsValidPolyRef(endRef) ||
		len(startPos) == 0 ||
		len(endPos) == 0 || filter == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}

	// trade quality with performance?
	if options&DT_FINDPATH_ANY_ANGLE > 0 {
		// limiting to several times the character radius yields nice results. It is not sensitive
		// so it is enough to compute it from the first tile.
		tile := q.m_nav.GetTileByRef(DtTileRef(startRef))
		agentRadius := tile.Header.WalkableRadius
		q.m_query.raycastLimitSqr = common.Sqr(agentRadius * DT_RAY_CAST_LIMIT_PROPORTIONS)
	}

	if startRef == endRef {
		q.m_query.status = DT_SUCCESS
		return DT_SUCCESS
	}

	q.m_nodePool.Clear()
	q.m_openList.Reset()
	startNode := q.m_nodePool.GetNode(startRef)
	copy(startNode.Pos[:], startPos)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = float32(common.Vdist(startPos, endPos)) * H_SCALE
	startNode.Id = startRef
	startNode.Flags = DT_NODE_OPEN
	q.m_openList.Offer(startNode)

	q.m_query.status = DT_IN_PROGRESS
	q.m_query.lastBestNode = startNode
	q.m_query.lastBestNodeCost = startNode.Total

	return q.m_query.status
}

func (q *DtNavMeshQuery) UpdateSlicedFindPath(maxIter int32) (doneIters int32, status DtStatus) {
	if !q.m_query.status.DtStatusInProgress() {
		return doneIters, q.m_query.status
	}

	// Make sure the request is still valid.
	if !q.m_nav.IsValidPolyRef(q.m_query.startRef) || !q.m_nav.IsValidPolyRef(q.m_query.endRef) {
		q.m_query.status = DT_FAILURE
		return doneIters, DT_FAILURE
	}

	var rayHit dtRaycastHit
	rayHit.maxPath = 0

	iter := int32(0)
	for iter < maxIter && !q.m_openList.Empty() {
		iter++

		// Remove node from open list and put it in closed list.
		bestNode := q.m_openList.Poll()
		bestNode.Flags &= ^uint32(DT_NODE_OPEN)
		bestNode.Flags |= DT_NODE_CLOSED

		// Reached the goal, stop searching.
		if bestNode.Id == q.m_query.endRef {
			q.m_query.lastBestNode = bestNode
			details := q.m_query.status & DT_STATUS_DETAIL_MASK
			q.m_query.status = DT_SUCCESS | details
			doneIters = iter

			return doneIters, q.m_query.status
		}

		// Get current poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.Id
		bestTile, bestPoly, status := q.m_nav.GetTileAndPolyByRef(bestRef)
		if status.DtStatusFailed() {
			// The polygon has disappeared during the sliced query, fail.
			q.m_query.status = DT_FAILURE
			doneIters = iter

			return doneIters, q.m_query.status
		}

		// Get parent and grand parent poly and tile.
		var parentRef, grandpaRef DtPolyRef
		var parentNode *DtNode
		if bestNode.Pidx != 0 {
			parentNode = q.m_nodePool.GetNodeAtIdx(int(bestNode.Pidx))
			parentRef = parentNode.Id
			if parentNode.Pidx != 0 {
				grandpaRef = q.m_nodePool.GetNodeAtIdx(int(parentNode.Pidx)).Id
			}

		}
		if parentRef > 0 {
			_, _, status := q.m_nav.GetTileAndPolyByRef(parentRef)
			invalidParent := status.DtStatusFailed()
			if invalidParent || (grandpaRef > 0 && !q.m_nav.IsValidPolyRef(grandpaRef)) {
				// The polygon has disappeared during the sliced query, fail.
				q.m_query.status = DT_FAILURE

				doneIters = iter
				return doneIters, q.m_query.status
			}
		}

		// decide whether to test raycast to previous nodes
		tryLOS := false
		if q.m_query.options&DT_FINDPATH_ANY_ANGLE > 0 {
			if (parentRef != 0) && (common.VdistSqr(parentNode.Pos[:], bestNode.Pos[:]) < float64(q.m_query.raycastLimitSqr)) {
				tryLOS = true
			}

		}

		for i := bestPoly.FirstLink; i != DT_NULL_LINK; i = bestTile.Links[i].Next {
			neighbourRef := bestTile.Links[i].Ref

			// Skip invalid ids and do not expand back to where we came from.
			if neighbourRef == 0 || neighbourRef == parentRef {
				continue
			}

			// Get neighbour poly and tile.
			// The API input has been checked already, skip checking internal data.

			neighbourTile, neighbourPoly := q.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			if !q.m_query.filter.passFilter(neighbourPoly) {
				continue
			}

			// get the neighbor node
			neighbourNode := q.m_nodePool.GetNode(neighbourRef, 0)
			if neighbourNode == nil {
				q.m_query.status |= DT_OUT_OF_NODES
				continue
			}

			// do not expand to nodes that were already visited from the same parent
			if neighbourNode.Pidx != 0 && neighbourNode.Pidx == bestNode.Pidx {
				continue
			}

			// If the node is visited the first time, calculate node position.
			if neighbourNode.Flags == 0 {
				pos, _ := q.getEdgeMidPoint1(bestRef, bestPoly, bestTile,
					neighbourRef, neighbourPoly, neighbourTile,
				)
				copy(neighbourNode.Pos[:], pos)
			}

			// Calculate cost and heuristic.
			cost := float32(0)
			heuristic := float32(0)

			// raycast parent
			foundShortCut := false
			rayHit.t = 0
			rayHit.pathCost = rayHit.t
			if tryLOS {
				q.Raycast1(parentRef, parentNode.Pos[:], neighbourNode.Pos[:], q.m_query.filter, DT_RAYCAST_USE_COSTS, &rayHit, grandpaRef)
				foundShortCut = rayHit.t >= 1.0
			}

			// update move cost
			if foundShortCut {
				// shortcut found using raycast. Using shorter cost instead
				cost = parentNode.Cost + rayHit.pathCost
			} else {
				// No shortcut found.
				curCost := q.m_query.filter.getCost(bestNode.Pos[:], neighbourNode.Pos[:], bestPoly)
				cost = bestNode.Cost + curCost
			}

			// Special case for last node.
			if neighbourRef == q.m_query.endRef {
				endCost := q.m_query.filter.getCost(neighbourNode.Pos[:], q.m_query.endPos[:], neighbourPoly)

				cost = cost + endCost
				heuristic = 0
			} else {
				heuristic = common.Vdist(neighbourNode.Pos[:], q.m_query.endPos[:]) * H_SCALE
			}

			total := cost + heuristic

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_OPEN > 0) && total >= neighbourNode.Total {
				continue
			}

			// The node is already visited and process, and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_CLOSED > 0) && total >= neighbourNode.Total {
				continue
			}

			// Add or update the node.
			neighbourNode.Pidx = uint32(q.m_nodePool.GetNodeIdx(bestNode))
			if foundShortCut {
				neighbourNode.Pidx = bestNode.Pidx
			}

			neighbourNode.Id = neighbourRef
			neighbourNode.Flags = neighbourNode.Flags & ^uint32(DT_NODE_CLOSED|DT_NODE_PARENT_DETACHED)
			neighbourNode.Cost = cost
			neighbourNode.Total = total
			if foundShortCut {
				neighbourNode.Flags = neighbourNode.Flags | DT_NODE_PARENT_DETACHED
			}

			if neighbourNode.Flags&DT_NODE_OPEN > 0 {
				// Already in open, update node location.
				q.m_openList.Update(neighbourNode)
			} else {
				// Put the node in open list.
				neighbourNode.Flags |= DT_NODE_OPEN
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
func (q *DtNavMeshQuery) FindPath(startRef, endRef DtPolyRef,
	startPos, endPos []float32, filter *DtQueryFilter, path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus) {
	// Validate input
	if !q.m_nav.IsValidPolyRef(startRef) || !q.m_nav.IsValidPolyRef(endRef) ||
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

	q.m_nodePool.Clear()
	q.m_openList.Reset()
	startNode := q.m_nodePool.GetNode(startRef)
	copy(startNode.Pos[:], startPos)
	startNode.Pidx = 0
	startNode.Cost = 0
	startNode.Total = common.Vdist(startPos, endPos) * H_SCALE
	startNode.Id = startRef
	startNode.Flags = DT_NODE_OPEN
	q.m_openList.Offer(startNode)

	lastBestNode := startNode
	lastBestNodeCost := startNode.Total

	outOfNodes := false

	for !q.m_openList.Empty() {
		// Remove node from open list and put it in closed list.
		bestNode := q.m_openList.Poll()
		bestNode.Flags &= ^uint32(DT_NODE_OPEN)
		bestNode.Flags |= DT_NODE_CLOSED

		// Reached the goal, stop searching.
		if bestNode.Id == endRef {
			lastBestNode = bestNode
			break
		}

		// Get current poly and tile.
		// The API input has been checked already, skip checking internal data.
		bestRef := bestNode.Id
		bestTile, bestPoly := q.m_nav.GetTileAndPolyByRefUnsafe(bestRef)

		// Get parent poly and tile.
		var parentRef DtPolyRef

		if bestNode.Pidx != 0 {
			parentRef = q.m_nodePool.GetNodeAtIdx(int(bestNode.Pidx)).Id
		}

		if parentRef > 0 {
			q.m_nav.GetTileAndPolyByRefUnsafe(parentRef)
		}

		for i := bestPoly.FirstLink; i != DT_NULL_LINK; i = bestTile.Links[i].Next {
			neighbourRef := bestTile.Links[i].Ref

			// Skip invalid ids and do not expand back to where we came from.
			if neighbourRef != 0 || neighbourRef == parentRef {
				continue
			}

			// Get neighbour poly and tile.
			// The API input has been checked already, skip checking internal data.

			neighbourTile, neighbourPoly := q.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			if !filter.passFilter(neighbourPoly) {
				continue
			}

			// deal explicitly with crossing tile boundaries
			crossSide := uint32(0)
			if bestTile.Links[i].Side != 0xff {
				crossSide = uint32(bestTile.Links[i].Side) >> 1
			}

			// get the node
			neighbourNode := q.m_nodePool.GetNode(neighbourRef, crossSide)
			if neighbourNode == nil {
				outOfNodes = true
				continue
			}

			// If the node is visited the first time, calculate node position.
			if neighbourNode.Flags == 0 {
				pos, _ := q.getEdgeMidPoint1(bestRef, bestPoly, bestTile,
					neighbourRef, neighbourPoly, neighbourTile)
				copy(neighbourNode.Pos[:], pos)

			}

			// Calculate cost and heuristic.
			cost := float32(0)
			heuristic := float32(0)

			// Special case for last node.
			if neighbourRef == endRef {
				curCost := filter.getCost(bestNode.Pos[:], neighbourNode.Pos[:], bestPoly)
				// Cost
				endCost := filter.getCost(neighbourNode.Pos[:], endPos, neighbourPoly)

				cost = bestNode.Cost + curCost + endCost
				heuristic = 0
			} else {
				// Cost
				curCost := filter.getCost(bestNode.Pos[:], neighbourNode.Pos[:], bestPoly)
				cost = bestNode.Cost + curCost
				heuristic = common.Vdist(neighbourNode.Pos[:], endPos) * H_SCALE
			}

			total := cost + heuristic

			// The node is already in open list and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_OPEN > 0) && total >= neighbourNode.Total {
				continue
			}

			// The node is already visited and process, and the new result is worse, skip.
			if (neighbourNode.Flags&DT_NODE_CLOSED > 0) && total >= neighbourNode.Total {
				continue
			}

			// Add or update the node.
			neighbourNode.Pidx = uint32(q.m_nodePool.GetNodeIdx(bestNode))
			neighbourNode.Id = neighbourRef
			neighbourNode.Flags = neighbourNode.Flags & ^uint32(DT_NODE_CLOSED)
			neighbourNode.Cost = cost
			neighbourNode.Total = total

			if neighbourNode.Flags&DT_NODE_OPEN > 0 {
				// Already in open, update node location.
				q.m_openList.Update(neighbourNode)
			} else {
				// Put the node in open list.
				neighbourNode.Flags |= DT_NODE_OPEN
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

	if lastBestNode.Id != endRef {
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
func (q *DtNavMeshQuery) QueryPolygons(center []float32, halfExtents []float32,
	filter *DtQueryFilter, query dtPolyQuery) DtStatus {

	if len(center) == 0 ||
		len(halfExtents) == 0 || !common.Visfinite(halfExtents) ||
		filter == nil || query == nil {
		return DT_FAILURE | DT_INVALID_PARAM
	}
	bmin := make([]float32, 3)
	bmax := make([]float32, 3)
	common.Vsub(bmin, center, halfExtents)
	common.Vadd(bmax, center, halfExtents)

	// Find tiles the query touches.

	minx, miny := q.m_nav.CalcTileLoc(bmin)
	maxx, maxy := q.m_nav.CalcTileLoc(bmax)

	MAX_NEIS := int32(32)

	for y := miny; y <= maxy; y++ {
		for x := minx; x <= maxx; x++ {
			neis, nneis := q.m_nav.GetTilesAt(x, y, MAX_NEIS)
			for j := int32(0); j < nneis; j++ {
				q.QueryPolygonsInTile(neis[j], bmin, bmax, filter, query)
			}
		}
	}

	return DT_SUCCESS
}

func (q *DtNavMeshQuery) FinalizeSlicedFindPathPartial(existing []DtPolyRef, existingSize int32,
	path []DtPolyRef, maxPath int32) (pathCount int32, status DtStatus) {

	if len(existing) == 0 || existingSize <= 0 || len(path) == 0 || maxPath <= 0 {
		return pathCount, DT_FAILURE | DT_INVALID_PARAM
	}

	if q.m_query.status.DtStatusFailed() {
		// Reset query.
		q.m_query = &dtQueryData{}
		return pathCount, DT_FAILURE
	}

	n := int32(0)

	if q.m_query.startRef == q.m_query.endRef {
		// Special case: the search starts and ends at same poly.
		path[n] = q.m_query.startRef
		n++
	} else {
		// Find furthest existing node that was visited.
		var prev *DtNode
		var node *DtNode
		for i := existingSize - 1; i >= 0; i-- {
			ns, _ := q.m_nodePool.FindNodes(existing[i], 1)
			if len(ns) > 0 {
				node = ns[0]
				break
			}

		}

		if node == nil {
			q.m_query.status |= DT_PARTIAL_RESULT
			node = q.m_query.lastBestNode
		}

		// Reverse the path.
		prevRay := int32(0)
		common.DoWhile(func() (stop bool) {
			next := q.m_nodePool.GetNodeAtIdx(int(node.Pidx))
			node.Pidx = uint32(q.m_nodePool.GetNodeIdx(prev))
			prev = node
			nextRay := int32(node.Flags & DT_NODE_PARENT_DETACHED)                         // keep track of whether parent is not adjacent (i.e. due to raycast shortcut)
			node.Flags = (node.Flags & ^uint32(DT_NODE_PARENT_DETACHED)) | uint32(prevRay) // and store it in the reversed path's node
			prevRay = nextRay
			node = next
			return false
		}, func() bool {
			return node != nil
		})
		// Store path
		node = prev
		common.DoWhile(func() (stop bool) {
			next := q.m_nodePool.GetNodeAtIdx(int(node.Pidx))
			status = 0
			if node.Flags&DT_NODE_PARENT_DETACHED > 0 {
				var (
					normal = make([]float32, 3)
					t      float32
					m      int32
				)
				q.Raycast(node.Id, node.Pos[:], next.Pos[:], q.m_query.filter, &t, normal, path[n:], &m, maxPath-n)
				n += m
				// raycast ends on poly boundary and the path might include the next poly boundary.
				if path[n-1] == next.Id {
					n--
				} // remove to avoid duplicates

			} else {
				path[n] = node.Id
				n++
				if n >= maxPath {
					status = DT_BUFFER_TOO_SMALL
				}

			}

			if status&DT_STATUS_DETAIL_MASK > 0 {
				q.m_query.status |= status & DT_STATUS_DETAIL_MASK
				return true
			}
			node = next
			return false
		}, func() bool {
			return node != nil
		})
	}

	details := q.m_query.status & DT_STATUS_DETAIL_MASK

	// Reset query.
	q.m_query = &dtQueryData{}

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
func (q *DtNavMeshQuery) FindLocalNeighbourhood(startRef DtPolyRef, centerPos []float32, radius float32,
	filter *DtQueryFilter,
	resultRef []DtPolyRef, resultParent []DtPolyRef, maxResult int32) (resultCount int32, status DtStatus) {
	if !q.m_nav.IsValidPolyRef(startRef) ||
		len(centerPos) == 0 ||
		radius < 0 || !common.IsFinite(radius) ||
		filter == nil || maxResult < 0 {
		return resultCount, DT_FAILURE | DT_INVALID_PARAM
	}

	MAX_STACK := 48
	var stack [48]*DtNode
	nstack := 0

	q.m_tinyNodePool.Clear()

	startNode := q.m_tinyNodePool.GetNode(startRef)
	startNode.Pidx = 0
	startNode.Id = startRef
	startNode.Flags = DT_NODE_CLOSED
	stack[nstack] = startNode
	nstack++

	radiusSqr := common.Sqr(radius)

	var pa [DT_VERTS_PER_POLYGON * 3]float32
	var pb [DT_VERTS_PER_POLYGON * 3]float32

	status = DT_SUCCESS

	n := int32(0)
	if n < maxResult {
		resultRef[n] = startNode.Id
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
		curRef := curNode.Id

		curTile, curPoly := q.m_nav.GetTileAndPolyByRefUnsafe(curRef)

		for i := curPoly.FirstLink; i != DT_NULL_LINK; i = curTile.Links[i].Next {
			link := curTile.Links[i]
			neighbourRef := link.Ref
			// Skip invalid neighbours.
			if neighbourRef == 0 {
				continue
			}

			// Skip if cannot alloca more nodes.
			neighbourNode := q.m_tinyNodePool.GetNode(neighbourRef)
			if neighbourNode == nil {
				continue
			}

			// Skip visited.
			if neighbourNode.Flags&DT_NODE_CLOSED > 0 {
				continue
			}

			// Expand to neighbour

			neighbourTile, neighbourPoly := q.m_nav.GetTileAndPolyByRefUnsafe(neighbourRef)

			// Skip off-mesh connections.
			if neighbourPoly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
				continue
			}

			// Do not advance if the polygon is excluded by the filter.
			if !filter.passFilter(neighbourPoly) {
				continue
			}
			va := make([]float32, 3)
			vb := make([]float32, 3)
			// Find edge and calc distance to the edge.
			status = q.getPortalPoints1(curRef, curPoly, curTile, neighbourRef, neighbourPoly, neighbourTile, va, vb)
			if status.DtStatusFailed() {
				continue
			}

			// If the circle is not touching the next polygon, skip it.
			_, distSqr := DtDistancePtSegSqr2D(centerPos, va, vb)
			if distSqr > radiusSqr {
				continue
			}

			// Mark node visited, this is done before the overlap test so that
			// we will not visit the poly again if the test fails.
			neighbourNode.Flags |= DT_NODE_CLOSED
			neighbourNode.Pidx = uint32(q.m_tinyNodePool.GetNodeIdx(curNode))

			// Check that the polygon does not collide with existing polygons.

			// Collect vertices of the neighbour poly.
			npa := int32(neighbourPoly.VertCount)
			for k := int32(0); k < npa; k++ {
				copy(pa[k*3:k*3+3], common.GetVert3(neighbourTile.Verts, neighbourPoly.Verts[k]))
			}

			overlap := false
			for j := int32(0); j < n; j++ {
				pastRef := resultRef[j]

				// Connected polys do not overlap.
				connected := false
				for k := curPoly.FirstLink; k != DT_NULL_LINK; k = curTile.Links[k].Next {
					if curTile.Links[k].Ref == pastRef {
						connected = true
						break
					}
				}
				if connected {
					continue
				}

				// Potentially overlapping.

				pastTile, pastPoly := q.m_nav.GetTileAndPolyByRefUnsafe(pastRef)

				// Get vertices and test overlap
				npb := int32(pastPoly.VertCount)

				for k := int32(0); k < npb; k++ {
					copy(pb[k*3:k*3+3], common.GetVert3(pastTile.Verts, pastPoly.Verts[k]))
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

// / @par
// /
// / @note If the search box does not intersect any polygons the search will
// / return #DT_SUCCESS, but @p nearestRef will be zero. So if in doubt, check
// / @p nearestRef before using @p nearestPt.
// /
func (q *DtNavMeshQuery) FindNearestPoly(center, halfExtents []float32, filter *DtQueryFilter, nearestPt []float32) (nearestRef DtPolyRef, status DtStatus) {
	nearestRef, _, status = q.FindNearestPoly1(center, halfExtents, filter, nearestPt)
	return
}

// If center and nearestPt point to an equal position, isOverPoly will be true;
// however there's also a special case of climb height inside the polygon (see dtFindNearestPolyQuery)
func (q *DtNavMeshQuery) FindNearestPoly1(center, halfExtents []float32,
	filter *DtQueryFilter, nearestPt []float32) (nearestRef DtPolyRef, isOverPoly bool, status DtStatus) {
	// queryPolygons below will check rest of params
	query := newDtFindNearestPolyQuery(q, center)
	status = q.QueryPolygons(center, halfExtents, filter, query)
	if status.DtStatusFailed() {
		return nearestRef, isOverPoly, status
	}

	nearestRef = query.nearestRef()
	// Only override nearestPt if we actually found a poly so the nearest point
	// is valid.
	if len(nearestPt) > 0 {
		copy(nearestPt, query.nearestPoint())
		isOverPoly = query.isOverPoly()

	}

	return nearestRef, isOverPoly, DT_SUCCESS
}

// / @par
// /
// / If no polygons are found, the function will return #DT_SUCCESS with a
// / @p polyCount of zero.
// /
// / If @p polys is too small to hold the entire result set, then the array will
// / be filled to capacity. The method of choosing which polygons from the
// / full set are included in the partial result set is undefined.
// /
func (q *DtNavMeshQuery) QueryPolygons1(center, halfExtents []float32,
	filter *DtQueryFilter,
	polys []DtPolyRef, maxPolys int32) (polyCount int32, status DtStatus) {
	if len(polys) == 0 || maxPolys < 0 {
		return polyCount, DT_FAILURE | DT_INVALID_PARAM
	}

	collector := newDtCollectPolysQuery(polys, maxPolys)
	status = q.QueryPolygons(center, halfExtents, filter, collector)
	if status.DtStatusFailed() {
		return polyCount, status
	}
	polyCount = collector.numCollected()
	if collector.overflowed() {
		return polyCount, DT_SUCCESS | DT_BUFFER_TOO_SMALL
	}
	return polyCount, DT_SUCCESS
}

// / @par
// /
// / If the @p segmentRefs parameter is provided, then all polygon segments will be returned.
// / Otherwise only the wall segments are returned.
// /
// / A segment that is normally a portal will be included in the result set as a
// / wall if the @p filter results in the neighbor polygon becoomming impassable.
// /
// / The @p segmentVerts and @p segmentRefs buffers should normally be sized for the
// / maximum segments per polygon of the source navigation mesh.
// /
func (q *DtNavMeshQuery) GetPolyWallSegments(ref DtPolyRef, filter *DtQueryFilter,
	segmentVerts []float32, segmentRefs []DtPolyRef,
	maxSegments int32) (segmentCount int32, status DtStatus) {

	tile, poly, status := q.m_nav.GetTileAndPolyByRef(ref)
	if status.DtStatusFailed() {
		return segmentCount, DT_FAILURE | DT_INVALID_PARAM
	}

	if filter == nil || len(segmentVerts) == 0 || maxSegments < 0 {
		return segmentCount, DT_FAILURE | DT_INVALID_PARAM
	}

	n := int32(0)
	MAX_INTERVAL := int32(16)
	var ints [16]*dtSegInterval
	var nints int32

	storePortals := len(segmentRefs) != 0

	status = DT_SUCCESS
	i := uint8(0)
	j := poly.VertCount - 1
	for i < poly.VertCount {
		// Skip non-solid edges.
		nints = int32(0)
		if poly.Neis[j]&DT_EXT_LINK > 0 {
			// Tile border.
			for k := poly.FirstLink; k != DT_NULL_LINK; k = tile.Links[k].Next {
				link := tile.Links[k]
				if link.Edge == j {
					if link.Ref != 0 {
						_, neiPoly := q.m_nav.GetTileAndPolyByRefUnsafe(link.Ref)
						if filter.passFilter(neiPoly) {
							nints = insertInterval(ints[:], MAX_INTERVAL, int8(link.Bmin), int8(link.Bmax), link.Ref)
						}
					}
				}
			}
		} else {
			// Internal edge
			neiRef := 0
			if poly.Neis[j] > 0 {
				idx := (poly.Neis[j] - 1)
				neiRef = int(q.m_nav.GetPolyRefBase(tile) | DtPolyRef(idx))
				if !filter.passFilter(tile.Polys[idx]) {
					neiRef = 0
				}

			}

			// If the edge leads to another polygon and portals are not stored, skip.
			if neiRef != 0 && !storePortals {
				j = i
				i++
				continue
			}

			if n < maxSegments {
				vj := common.GetVert3(tile.Verts, poly.Verts[j])
				vi := common.GetVert3(tile.Verts, poly.Verts[i])
				seg := segmentVerts[n*6 : n*6+6]
				copy(seg[:3], vj)
				copy(seg[3:], vi)
				if len(segmentRefs) > 0 {
					segmentRefs[n] = DtPolyRef(neiRef)
				}

				n++
			} else {
				status |= DT_BUFFER_TOO_SMALL
			}
			j = i
			i++
			continue
		}

		// Add sentinels
		nints = insertInterval(ints[:], MAX_INTERVAL, -1, 0, 0)
		nints = insertInterval(ints[:], MAX_INTERVAL, 255, 256, 0)

		// Store segments.
		vj := common.GetVert3(tile.Verts, poly.Verts[j])
		vi := common.GetVert3(tile.Verts, poly.Verts[i])
		for k := int32(1); k < nints; k++ {
			// Portal segment.
			if storePortals && ints[k].ref > 0 {
				tmin := float32(ints[k].tmin) / 255.0
				tmax := float32(ints[k].tmax) / 255.0
				if n < maxSegments {
					seg := segmentVerts[n*6 : n*6+6]
					common.Vlerp(seg[:3], vj, vi, tmin)
					common.Vlerp(seg[3:], vj, vi, tmax)

					if len(segmentRefs) > 0 {
						segmentRefs[n] = ints[k].ref
					}

					n++
				} else {
					status |= DT_BUFFER_TOO_SMALL
				}
			}

			// Wall segment.
			imin := ints[k-1].tmax
			imax := ints[k].tmin
			if imin != imax {
				tmin := float32(imin) / 255.0
				tmax := float32(imax) / 255.0
				if n < maxSegments {
					seg := segmentVerts[n*6 : n*6+6]
					common.Vlerp(seg[:3], vj, vi, tmin)
					common.Vlerp(seg[3:], vj, vi, tmax)

					if len(segmentRefs) > 0 {
						segmentRefs[n] = 0
					}

					n++
				} else {
					status |= DT_BUFFER_TOO_SMALL
				}
			}
		}
		j = i
		i++
	}

	segmentCount = n

	return segmentCount, status
}

// / @par
// /
// / Will return #DT_FAILURE | DT_INVALID_PARAM if the provided position is outside the xz-bounds
// / of the polygon.
// /
func (q *DtNavMeshQuery) GetPolyHeight(ref DtPolyRef, pos []float32) (height float32, status DtStatus) {
	tile, poly, status := q.m_nav.GetTileAndPolyByRef(ref)
	if status.DtStatusFailed() {
		return height, DT_FAILURE | DT_INVALID_PARAM
	}

	if len(pos) == 0 {
		return height, DT_FAILURE | DT_INVALID_PARAM
	}

	// We used to return success for offmesh connections, but the
	// getPolyHeight in DetourNavMesh does not do this, so special
	// case it here.
	if poly.GetType() == DT_POLYTYPE_OFFMESH_CONNECTION {
		v0 := common.GetVert3(tile.Verts, poly.Verts[0])
		v1 := common.GetVert3(tile.Verts, poly.Verts[1])
		t, _ := DtDistancePtSegSqr2D(pos, v0, v1)
		height = v0[1] + (v1[1]-v0[1])*t
		return height, DT_SUCCESS
	}

	height, ok := q.m_nav.GetPolyHeight(tile, poly, pos)
	if ok {
		return height, DT_SUCCESS
	}
	return height, DT_FAILURE | DT_INVALID_PARAM

}

func (q *DtNavMeshQuery) IsValidPolyRef(ref DtPolyRef, filter *DtQueryFilter) bool {

	_, poly, status := q.m_nav.GetTileAndPolyByRef(ref)
	// If cannot get polygon, assume it does not exists and boundary is invalid.
	if status.DtStatusFailed() {
		return false
	}

	// If cannot pass filter, assume flags has changed and boundary is invalid.
	if !filter.passFilter(poly) {
		return false
	}

	return true
}
