package detour_crowd

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/detour"
)

func DtMergeCorridorStartMoved(path []detour.DtPolyRef, npath, maxPath int32,
	visited []detour.DtPolyRef, nvisited int32) int32 {
	furthestPath := int32(-1)
	furthestVisited := int32(-1)

	// Find furthest common polygon.
	for i := npath - 1; i >= 0; i-- {
		found := false
		for j := nvisited - 1; j >= 0; j-- {
			if path[i] == visited[j] {
				furthestPath = i
				furthestVisited = j
				found = true
			}
		}
		if found {
			break
		}

	}

	// If no intersection found just return current path.
	if furthestPath == -1 || furthestVisited == -1 {
		return npath
	}

	// Concatenate paths.

	// Adjust beginning of the buffer to include the visited.
	req := nvisited - furthestVisited
	orig := min(furthestPath+1, npath)
	size := int32(max(0, npath-orig))
	if req+size > maxPath {
		size = maxPath - req
	}

	if size > 0 {
		copy(path[req:], path[orig:orig+size])
	}

	// Store visited
	for i := int32(0); i < req; i++ {
		path[i] = visited[(nvisited-1)-i]
	}

	return req + size
}

func dtMergeCorridorEndMoved(path []detour.DtPolyRef, npath, maxPath int, visited []detour.DtPolyRef, nvisited int) int {
	furthestPath := -1
	furthestVisited := -1

	// Find furthest common polygon.
	for i := 0; i < npath; i++ {
		found := false
		for j := nvisited - 1; j >= 0; j-- {
			if path[i] == visited[j] {
				furthestPath = i
				furthestVisited = j
				found = true
			}
		}
		if found {
			break
		}

	}

	// If no intersection found just return current path.
	if furthestPath == -1 || furthestVisited == -1 {
		return npath
	}

	// Concatenate paths.
	ppos := furthestPath + 1
	vpos := furthestVisited + 1
	count := min(nvisited-vpos, maxPath-ppos)
	if ppos+count <= maxPath {
		panic("")
	}
	if count > 0 {
		copy(path[ppos:], visited[vpos:vpos+count])
	}

	return ppos + count
}

func dtMergeCorridorStartShortcut(path []detour.DtPolyRef, npath, maxPath int32, visited []detour.DtPolyRef, nvisited int32) int32 {
	furthestPath := int32(-1)
	furthestVisited := int32(-1)

	// Find furthest common polygon.
	for i := npath - 1; i >= 0; i-- {
		found := false
		for j := nvisited - 1; j >= 0; j-- {
			if path[i] == visited[j] {
				furthestPath = i
				furthestVisited = j
				found = true
			}
		}
		if found {
			break
		}

	}

	// If no intersection found just return current path.
	if furthestPath == -1 || furthestVisited == -1 {
		return npath
	}

	// Concatenate paths.

	// Adjust beginning of the buffer to include the visited.
	req := furthestVisited
	if req <= 0 {
		return npath
	}

	orig := furthestPath
	size := max(0, npath-orig)
	if req+size > maxPath {
		size = maxPath - req
	}

	if size > 0 {
		copy(path[req:], path[orig:orig+size])
	}

	// Store visited
	for i := int32(0); i < req; i++ {
		path[i] = visited[i]
	}

	return req + size
}

/**
@class dtPathCorridor
@par

The corridor is loaded with a path, usually obtained from a #DtNavMeshQuery::findPath() query. The corridor
is then used to plan local movement, with the corridor automatically updating as needed to deal with inaccurate
agent locomotion.

Example of a common use case:

-# Construct the corridor object and call #init() to allocate its path buffer.
-# Obtain a path from a #DtNavMeshQuery object.
-# Use #reset() to set the agent's current position. (At the beginning of the path.)
-# Use #setCorridor() to load the path and target.
-# Use #findCorners() to plan movement. (This handles dynamic path straightening.)
-# Use #movePosition() to feed agent movement back into the corridor. (The corridor will automatically adjust as needed.)
-# If the target is moving, use #moveTargetPosition() to update the end of the corridor.
   (The corridor will automatically adjust as needed.)
-# Repeat the previous 3 steps to continue to move the agent.

The corridor position and target are always constrained to the navigation mesh.

One of the difficulties in maintaining a path is that floating point errors, locomotion inaccuracies, and/or local
steering can result in the agent crossing the boundary of the path corridor, temporarily invalidating the path.
This class uses local mesh queries to detect and update the corridor as needed to handle these types of issues.

The fact that local mesh queries are used to move the position and target locations results in two beahviors that
need to be considered:

Every time a move function is used there is a chance that the path will become non-optimial. Basically, the further
the target is moved from its original location, and the further the position is moved outside the original corridor,
the more likely the path will become non-optimal. This issue can be addressed by periodically running the
#optimizePathTopology() and #optimizePathVisibility() methods.

All local mesh queries have distance limitations. (Review the #DtNavMeshQuery methods for details.) So the most accurate
use case is to move the position and target in small increments. If a large increment is used, then the corridor
may not be able to accurately find the new location.  Because of this limiation, if a position is moved in a large
increment, then compare the desired and resulting polygon references. If the two do not match, then path replanning
may be needed.  E.g. If you move the target, check #getLastPoly() to see if it is the expected polygon.

*/
// / Represents a dynamic polygon corridor used to plan agent movement.
// / @ingroup crowd, detour
type dtPathCorridor struct {
	m_pos     []float32
	m_target  []float32
	m_path    []detour.DtPolyRef
	m_npath   int32
	m_maxPath int32
}

// / Gets the current position within the corridor. (In the first polygon.)
// / @return The current position within the corridor.
func (d *dtPathCorridor) GetPos() []float32 { return d.m_pos }

// / Gets the current target within the corridor. (In the last polygon.)
// / @return The current target within the corridor.
func (d *dtPathCorridor) getTarget() []float32 { return d.m_target }

// / The polygon reference id of the first polygon in the corridor, the polygon containing the position.
// / @return The polygon reference id of the first polygon in the corridor. (Or zero if there is no path.)
func (d *dtPathCorridor) getFirstPoly() detour.DtPolyRef {
	if d.m_npath > 0 {
		return d.m_path[0]
	}
	return 0
}

// / The polygon reference id of the last polygon in the corridor, the polygon containing the target.
// / @return The polygon reference id of the last polygon in the corridor. (Or zero if there is no path.)
func (d *dtPathCorridor) getLastPoly() detour.DtPolyRef {
	if d.m_npath > 0 {
		return d.m_path[d.m_npath-1]
	}
	return 0
}

// / The corridor's path.
// / @return The corridor's path. [(polyRef) * #getPathCount()]
func (d *dtPathCorridor) GetPath() []detour.DtPolyRef { return d.m_path }

// / The number of polygons in the current corridor path.
// / @return The number of polygons in the current corridor path.
func (d *dtPathCorridor) GetPathCount() int { return int(d.m_npath) }

// / @par
// /
// / @warning Cannot be called more than once.
func newDtPathCorridor(maxPath int32) *dtPathCorridor {
	d := &dtPathCorridor{
		m_pos:    make([]float32, 3),
		m_target: make([]float32, 3),
	}
	d.m_path = make([]detour.DtPolyRef, maxPath)
	d.m_npath = 0
	d.m_maxPath = maxPath
	return d
}

// / @par
// /
// / Essentially, the corridor is set of one polygon in size with the target
// / equal to the position.
func (d *dtPathCorridor) reset(ref detour.DtPolyRef, pos []float32) {
	copy(d.m_pos[:], pos)
	copy(d.m_target[:], pos)
	d.m_path[0] = ref
	d.m_npath = 1
}

/*
*

	@par

	This is the function used to plan local movement within the corridor. One or more corners can be
	detected in order to plan movement. It performs essentially the same function as #DtNavMeshQuery::findStraightPath.

	Due to internal optimizations, the maximum number of corners returned will be (@p maxCorners - 1)
	For example: If the buffers are sized to hold 10 corners, the function will never return more than 9 corners.
	So if 10 corners are needed, the buffers should be sized for 11 corners.

	If the target is within range, it will be the last corner and have a polygon reference id of zero.
*/
func (d *dtPathCorridor) findCorners(cornerVerts []float32, cornerFlags []int32,
	cornerPolys []detour.DtPolyRef, maxCorners int32,
	navquery detour.NavMeshQuery, filter *detour.DtQueryFilter) int32 {
	MIN_TARGET_DIST := 0.01

	ncorners := int32(0)
	ncorners, _ = navquery.FindStraightPath(d.m_pos[:], d.m_target[:], d.m_path, d.m_npath,
		cornerVerts, cornerFlags, cornerPolys, maxCorners, 0)

	// Prune points in the beginning of the path which are too close.
	for ncorners > 0 {
		if (cornerFlags[0]&detour.DT_STRAIGHTPATH_OFFMESH_CONNECTION > 0) || float64(common.Vdist2DSqr(cornerVerts, d.m_pos[:])) > common.Sqr(MIN_TARGET_DIST) {
			break
		}

		ncorners--
		if ncorners > 0 {
			copy(cornerFlags, cornerFlags[1:1+ncorners])
			copy(cornerPolys, cornerPolys[1:1+ncorners])
			copy(cornerVerts, cornerVerts[3:3+3*ncorners])
		}
	}

	// Prune points after an off-mesh connection.
	for i := int32(0); i < ncorners; i++ {
		if cornerFlags[i]&detour.DT_STRAIGHTPATH_OFFMESH_CONNECTION > 0 {
			ncorners = i + 1
			break
		}
	}

	return ncorners
}

/*
*
@par

Inaccurate locomotion or dynamic obstacle avoidance can force the argent position significantly outside the
original corridor. Over time this can result in the formation of a non-optimal corridor. Non-optimal paths can
also form near the corners of tiles.

This function uses an efficient local visibility search to try to optimize the corridor
between the current position and @p next.

The corridor will change only if @p next is visible from the current position and moving directly toward the point
is better than following the existing path.

The more inaccurate the agent movement, the more beneficial this function becomes. Simply adjust the frequency
of the call to match the needs to the agent.

This function is not suitable for long distance searches.
*/
func (d *dtPathCorridor) optimizePathVisibility(next []float32, pathOptimizationRange float32,
	navquery detour.NavMeshQuery, filter *detour.DtQueryFilter) {
	// Clamp the ray to max distance.
	goal := make([]float32, 3)
	copy(goal, next)
	dist := common.Vdist2D(d.m_pos[:], goal)

	// If too close to the goal, do not try to optimize.
	if dist < 0.01 {
		return
	}

	// Overshoot a little. This helps to optimize open fields in tiled meshes.
	dist = min(dist+0.01, pathOptimizationRange)

	// Adjust ray length.
	delta := make([]float32, 3)
	common.Vsub(delta, goal, d.m_pos[:])
	common.Vmad(goal, d.m_pos[:], delta, pathOptimizationRange/dist)

	MAX_RES := int32(32)
	res := make([]detour.DtPolyRef, MAX_RES)

	var (
		t float32
	)
	norm := make([]float32, 3)
	nres := int32(0)

	navquery.Raycast(d.m_path[0], d.m_pos[:], goal, filter, &t, norm, res, &nres, MAX_RES)
	if nres > 1 && t > 0.99 {
		d.m_npath = dtMergeCorridorStartShortcut(d.m_path, d.m_npath, d.m_maxPath, res, int32(nres))
	}
}

/*
*

	@par

	Inaccurate locomotion or dynamic obstacle avoidance can force the agent position significantly outside the
	original corridor. Over time this can result in the formation of a non-optimal corridor. This function will use a
	local area path search to try to re-optimize the corridor.

	The more inaccurate the agent movement, the more beneficial this function becomes. Simply adjust the frequency of
	the call to match the needs to the agent.
*/
func (d *dtPathCorridor) optimizePathTopology(navquery detour.NavMeshQuery, filter *detour.DtQueryFilter) bool {
	if d.m_npath < 3 {
		return false
	}

	MAX_ITER := 32
	MAX_RES := 32

	res := make([]detour.DtPolyRef, MAX_RES)
	navquery.InitSlicedFindPath(d.m_path[0], d.m_path[d.m_npath-1], d.m_pos[:], d.m_target[:], filter, 0)
	navquery.UpdateSlicedFindPath(int32(MAX_ITER))
	nres, status := navquery.FinalizeSlicedFindPathPartial(d.m_path, int32(d.m_npath), res, int32(MAX_RES))
	if status.DtStatusSucceed() && nres > 0 {
		d.m_npath = dtMergeCorridorStartShortcut(d.m_path, d.m_npath, d.m_maxPath, res, int32(nres))
		return true
	}

	return false
}

func (d *dtPathCorridor) moveOverOffmeshConnection(offMeshConRef detour.DtPolyRef, refs []detour.DtPolyRef,
	startPos, endPos []float32,
	navquery detour.NavMeshQuery) bool {
	// Advance the path up to and over the off-mesh connection.
	var prevRef detour.DtPolyRef
	polyRef := d.m_path[0]
	npos := int32(0)
	for npos < d.m_npath && polyRef != offMeshConRef {
		prevRef = polyRef
		polyRef = d.m_path[npos]
		npos++
	}
	if npos == d.m_npath {
		// Could not find offMeshConRef
		return false
	}

	// Prune path
	for i := npos; i < d.m_npath; i++ {
		d.m_path[i-npos] = d.m_path[i]
	}

	d.m_npath -= npos

	refs[0] = prevRef
	refs[1] = polyRef

	nav := navquery.GetAttachedNavMesh()

	status := nav.GetOffMeshConnectionPolyEndPoints(refs[0], refs[1], startPos, endPos)
	if status.DtStatusSucceed() {
		copy(d.m_pos[:], endPos)
		return true
	}

	return false
}

/*
*
@par

Behavior:

- The movement is constrained to the surface of the navigation mesh.
- The corridor is automatically adjusted (shorted or lengthened) in order to remain valid.
- The new position will be located in the adjusted corridor's first polygon.

The expected use case is that the desired position will be 'near' the current corridor. What is considered 'near'
depends on local polygon density, query search half extents, etc.

The resulting position will differ from the desired position if the desired position is not on the navigation mesh,
or it can't be reached using a local search.
*/
func (d *dtPathCorridor) movePosition(npos []float32, navquery detour.NavMeshQuery, filter *detour.DtQueryFilter) bool {
	// Move along navmesh and update new position.
	result := make([]float32, 3)
	MAX_VISITED := 16
	visited := make([]detour.DtPolyRef, MAX_VISITED)
	nvisited := int32(0)
	status := navquery.MoveAlongSurface(d.m_path[0], d.m_pos[:], npos, filter, result, visited, &nvisited, int32(MAX_VISITED))
	if status.DtStatusSucceed() {
		d.m_npath = DtMergeCorridorStartMoved(d.m_path, d.m_npath, d.m_maxPath, visited, nvisited)

		// Adjust the position to stay on top of the navmesh.
		h := d.m_pos[1]
		h, _ = navquery.GetPolyHeight(d.m_path[0], result)
		result[1] = h
		copy(d.m_pos[:], result)
		return true
	}
	return false
}

/*
*
@par

Behavior:

- The movement is constrained to the surface of the navigation mesh.
- The corridor is automatically adjusted (shorted or lengthened) in order to remain valid.
- The new target will be located in the adjusted corridor's last polygon.

The expected use case is that the desired target will be 'near' the current corridor. What is considered 'near' depends on local polygon density, query search half extents, etc.

The resulting target will differ from the desired target if the desired target is not on the navigation mesh, or it can't be reached using a local search.
*/
func (d *dtPathCorridor) moveTargetPosition(npos []float32, navquery *detour.DtNavMeshQuery, filter *detour.DtQueryFilter) bool {
	// Move along navmesh and update new position.
	result := make([]float32, 3)
	MAX_VISITED := 16
	visited := make([]detour.DtPolyRef, MAX_VISITED)
	nvisited := int32(0)
	var status detour.DtStatus
	status = navquery.MoveAlongSurface(d.m_path[d.m_npath-1], d.m_target[:], npos, filter, result, visited, &nvisited,
		int32(MAX_VISITED))
	if status.DtStatusSucceed() {
		d.m_npath = int32(dtMergeCorridorEndMoved(d.m_path, int(d.m_npath), int(d.m_maxPath), visited, int(nvisited)))
		// TODO: should we do that?
		// Adjust the position to stay on top of the navmesh.
		/*	float h = m_target[1];
			navquery->getPolyHeight(m_path[m_npath-1], result, &h);
			result[1] = h;*/

		copy(d.m_target[:], result)

		return true
	}
	return false
}

// / @par
// /
// / The current corridor position is expected to be within the first polygon in the path. The target
// / is expected to be in the last polygon.
// /
// / @warning The size of the path must not exceed the size of corridor's path buffer set during #init().
func (d *dtPathCorridor) setCorridor(target []float32, path []detour.DtPolyRef, npath int) {
	copy(d.m_target[:], target)
	copy(d.m_path, path[:npath])
	d.m_npath = int32(npath)
}

func (d *dtPathCorridor) fixPathStart(safeRef detour.DtPolyRef, safePos []float32) bool {

	copy(d.m_pos[:], safePos)
	if d.m_npath < 3 && d.m_npath > 0 {
		d.m_path[2] = d.m_path[d.m_npath-1]
		d.m_path[0] = safeRef
		d.m_path[1] = 0
		d.m_npath = 3
	} else {
		d.m_path[0] = safeRef
		d.m_path[1] = 0
	}

	return true
}

func (d *dtPathCorridor) trimInvalidPath(safeRef detour.DtPolyRef, safePos []float32,
	navquery *detour.DtNavMeshQuery, filter *detour.DtQueryFilter) bool {
	// Keep valid path as far as possible.
	n := int32(0)
	for n < d.m_npath && navquery.IsValidPolyRef(d.m_path[n], filter) {
		n++
	}

	if n == d.m_npath {
		// All valid, no need to fix.
		return true
	} else if n == 0 {
		// The first polyref is bad, use current safe values.
		copy(d.m_pos[:], safePos)
		d.m_path[0] = safeRef
		d.m_npath = 1
	} else {
		// The path is partially usable.
		d.m_npath = n
	}

	// Clamp target pos to last poly
	tgt := make([]float32, 3)
	copy(tgt, d.m_target[:])
	tgt, _ = navquery.ClosestPointOnPolyBoundary(d.m_path[d.m_npath-1], d.m_target[:])

	return true
}

// / @par
// /
// / The path can be invalidated if there are structural changes to the underlying navigation mesh, or the state of
// / a polygon within the path changes resulting in it being filtered out. (E.g. An exclusion or inclusion flag changes.)
func (d *dtPathCorridor) isValid(maxLookAhead int32, navquery detour.NavMeshQuery, filter *detour.DtQueryFilter) bool {
	// Check that all polygons still pass query filter.
	n := min(d.m_npath, maxLookAhead)
	for i := int32(0); i < n; i++ {
		if !navquery.IsValidPolyRef(d.m_path[i], filter) {
			return false
		}

	}

	return true
}
