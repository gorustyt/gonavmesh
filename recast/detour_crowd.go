package recast

import (
	"math"
	"reflect"
)

// / The maximum number of neighbors that a crowd agent can take into account
// / for steering decisions.
// / @ingroup crowd
const DT_CROWDAGENT_MAX_NEIGHBOURS = 6

// / The maximum number of corners a crowd agent will look ahead in the path.
// / This value is used for sizing the crowd agent corner buffers.
// / Due to the behavior of the crowd manager, the actual number of useful
// / corners will be one less than this number.
// / @ingroup crowd
const DT_CROWDAGENT_MAX_CORNERS = 4

// / The maximum number of crowd avoidance configurations supported by the
// / crowd manager.
// / @ingroup crowd
// / @see dtObstacleAvoidanceParams, dtCrowd::setObstacleAvoidanceParams(), dtCrowd::getObstacleAvoidanceParams(),
// /		 dtCrowdAgentParams::obstacleAvoidanceType
const DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS = 8

// / The maximum number of query filter types supported by the crowd manager.
// / @ingroup crowd
// / @see dtQueryFilter, dtCrowd::getFilter() dtCrowd::getEditableFilter(),
// /		dtCrowdAgentParams::queryFilterType
const DT_CROWD_MAX_QUERY_FILTER_TYPE = 16

// / Provides neighbor data for agents managed by the crowd.
// / @ingroup crowd
// / @see dtCrowdAgent::neis, dtCrowd
type dtCrowdNeighbour struct {
	idx  int     ///< The index of the neighbor in the crowd.
	dist float64 ///< The distance between the current agent and the neighbor.
}

/// The type of navigation mesh polygon the agent is currently traversing.
/// @ingroup crowd

const (
	DT_CROWDAGENT_STATE_INVALID = 0 ///< The agent is not in a valid state.
	DT_CROWDAGENT_STATE_WALKING = 1 ///< The agent is traversing a normal navigation mesh polygon.
	DT_CROWDAGENT_STATE_OFFMESH = 2 ///< The agent is traversing an off-mesh connection.
)

// / Configuration parameters for a crowd agent.
// / @ingroup crowd
type dtCrowdAgentParams struct {
	radius          float64 ///< Agent radius. [Limit: >= 0]
	height          float64 ///< Agent height. [Limit: > 0]
	maxAcceleration float64 ///< Maximum allowed acceleration. [Limit: >= 0]
	maxSpeed        float64 ///< Maximum allowed speed. [Limit: >= 0]

	/// Defines how close a collision element must be before it is considered for steering behaviors. [Limits: > 0]
	collisionQueryRange float64

	pathOptimizationRange float64 ///< The path visibility optimization range. [Limit: > 0]

	/// How aggresive the agent manager should be at avoiding collisions with this agent. [Limit: >= 0]
	separationWeight float64

	/// Flags that impact steering behavior. (See: #UpdateFlags)
	updateFlags int

	/// The index of the avoidance configuration to use for the agent.
	/// [Limits: 0 <= value <= #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
	obstacleAvoidanceType int

	/// The index of the query filter used by this agent.
	queryFilterType int

	/// User defined data attached to the agent.
	//void *userData
}

const (
	DT_CROWDAGENT_TARGET_NONE = iota
	DT_CROWDAGENT_TARGET_FAILED
	DT_CROWDAGENT_TARGET_VALID
	DT_CROWDAGENT_TARGET_REQUESTING
	DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE
	DT_CROWDAGENT_TARGET_WAITING_FOR_PATH
	DT_CROWDAGENT_TARGET_VELOCITY
)

// / Represents an agent managed by a #dtCrowd object.
// / @ingroup crowd
type dtCrowdAgent struct {
	/// True if the agent is active, false if the agent is in an unused slot in the agent pool.
	active bool

	/// The type of mesh polygon the agent is traversing. (See: #CrowdAgentState)
	state int

	/// True if the agent has valid path (targetState == DT_CROWDAGENT_TARGET_VALID) and the path does not lead to the requested position, else false.
	partial bool

	/// The path corridor the agent is using.
	corridor *dtPathCorridor

	/// The local boundary data for the agent.
	boundary *dtLocalBoundary

	/// Time since the agent's path corridor was optimized.
	topologyOptTime float64

	/// The known neighbors of the agent.
	neis [DT_CROWDAGENT_MAX_NEIGHBOURS]*dtCrowdNeighbour

	/// The number of neighbors.
	nneis int

	/// The desired speed.
	desiredSpeed float64

	npos [3]float64 ///< The current agent position. [(x, y, z)]
	disp [3]float64 ///< A temporary value used to accumulate agent displacement during iterative collision resolution. [(x, y, z)]
	dvel [3]float64 ///< The desired velocity of the agent. Based on the current path, calculated from scratch each frame. [(x, y, z)]
	nvel [3]float64 ///< The desired velocity adjusted by obstacle avoidance, calculated from scratch each frame. [(x, y, z)]
	vel  [3]float64 ///< The actual velocity of the agent. The change from nvel -> vel is constrained by max acceleration. [(x, y, z)]

	/// The agent's configuration parameters.
	params *dtCrowdAgentParams

	/// The local path corridor corners for the agent. (Staight path.) [(x, y, z) * #ncorners]
	cornerVerts [DT_CROWDAGENT_MAX_CORNERS * 3]float64

	/// The local path corridor corner flags. (See: #dtStraightPathFlags) [(flags) * #ncorners]
	cornerFlags [DT_CROWDAGENT_MAX_CORNERS]int

	/// The reference id of the polygon being entered at the corner. [(polyRef) * #ncorners]
	cornerPolys [DT_CROWDAGENT_MAX_CORNERS]dtPolyRef

	/// The number of corners.
	ncorners int

	targetState      int            ///< State of the movement request.
	targetRef        dtPolyRef      ///< Target polyref of the movement request.
	targetPos        [3]float64     ///< Target position of the movement request (or velocity in case of DT_CROWDAGENT_TARGET_VELOCITY).
	targetPathqRef   dtPathQueueRef ///< Path finder ref.
	targetReplan     bool           ///< Flag indicating that the current path is being replanned.
	targetReplanTime float64        /// <Time since the agent's target was replanned.
}

type dtCrowdAgentAnimation struct {
	active                    bool
	initPos, startPos, endPos [3]float64
	polyRef                   dtPolyRef
	t, tmax                   float64
}

const (
	DT_CROWD_ANTICIPATE_TURNS   = 1
	DT_CROWD_OBSTACLE_AVOIDANCE = 2
	DT_CROWD_SEPARATION         = 4
	DT_CROWD_OPTIMIZE_VIS       = 8  ///< Use #dtPathCorridor::optimizePathVisibility() to optimize the agent path.
	DT_CROWD_OPTIMIZE_TOPO      = 16 ///< Use dtPathCorridor::optimizePathTopology() to optimize the agent
)

type dtCrowdAgentDebugInfo struct {
	idx              int
	optStart, optEnd [3]float64
	vod              *dtObstacleAvoidanceDebugData
}

// / Provides local steering behaviors for a group of agents.
// / @ingroup crowd
type dtCrowd struct {
	m_maxAgents    int
	m_agents       []*dtCrowdAgent
	m_activeAgents []*dtCrowdAgent
	m_agentAnims   []*dtCrowdAgentAnimation

	m_pathq *dtPathQueue

	m_obstacleQueryParams [DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]*dtObstacleAvoidanceParams
	m_obstacleQuery       *dtObstacleAvoidanceQuery

	m_grid *dtProximityGrid

	m_pathResult    []dtPolyRef
	m_maxPathResult int

	m_agentPlacementHalfExtents [3]float64

	m_filters [DT_CROWD_MAX_QUERY_FILTER_TYPE]*dtQueryFilter

	m_maxAgentRadius float64

	m_velocitySampleCount int

	m_navquery NavMeshQuery
}

func (d *dtCrowd) getAgentIndex(agent *dtCrowdAgent) int {
	for i, v := range d.m_agents {
		if reflect.DeepEqual(v, agent) {
			return i
		}
	}
	return -1
}

const MAX_ITERS_PER_UPDATE = 100

const MAX_PATHQUEUE_NODES = 4096

const MAX_COMMON_NODES = 512

func tween(t, t0, t1 float64) float64 {
	return dtClamp((t-t0)/(t1-t0), 0.0, 1.0)
}

func integrate(ag *dtCrowdAgent, dt float64) {
	// Fake dynamic constraint.
	maxDelta := ag.params.maxAcceleration * dt

	dv := dtVsub(ag.nvel[:], ag.vel[:])
	ds := dtVlen(dv)
	if ds > maxDelta {
		dv = dtVscale(dv, maxDelta/ds)
	}
	copy(ag.vel[:], dtVadd(ag.vel[:], dv))

	// Integrate
	if dtVlen(ag.vel[:]) > 0.0001 {
		dtVmad(ag.npos[:], ag.npos[:], ag.vel[:], dt)
	} else {
		dtVset(ag.vel[:], 0, 0, 0)
	}

}

func overOffmeshConnection(ag *dtCrowdAgent, radius float64) bool {
	if ag.ncorners == 0 {
		return false
	}

	offMeshConnection := false
	if ag.cornerFlags[ag.ncorners-1]&DT_STRAIGHTPATH_OFFMESH_CONNECTION > 0 {
		offMeshConnection = true
	}
	if offMeshConnection {
		distSq := dtVdist2DSqr(ag.npos[:], rcGetVert(ag.cornerVerts[:], ag.ncorners-1))
		if distSq < radius*radius {
			return true
		}

	}

	return false
}

func getDistanceToGoal(ag *dtCrowdAgent, rangef float64) float64 {
	if ag.ncorners == 0 {
		return rangef
	}

	endOfPath := false
	if ag.cornerFlags[ag.ncorners-1]&DT_STRAIGHTPATH_END > 0 {
		endOfPath = true
	}
	if endOfPath {
		return dtMin(dtVdist2D(ag.npos[:], rcGetVert(ag.cornerVerts[:], ag.ncorners-1)), rangef)
	}

	return rangef
}

func calcSmoothSteerDirection(ag *dtCrowdAgent, dir []float64) {
	if ag.ncorners == 0 {
		dtVset(dir, 0, 0, 0)
		return
	}

	ip0 := 0
	ip1 := dtMin(1, ag.ncorners-1)
	p0 := rcGetVert(ag.cornerVerts[:], ip0)
	p1 := rcGetVert(ag.cornerVerts[:], ip1)

	dir0 := dtVsub(p0, ag.npos[:])
	dir1 := dtVsub(p1, ag.npos[:])
	dir0[1] = 0
	dir1[1] = 0

	len0 := dtVlen(dir0)
	len1 := dtVlen(dir1)
	if len1 > 0.001 {
		dir1 = dtVscale(dir1, 1.0/len1)
	}

	dir[0] = dir0[0] - dir1[0]*len0*0.5
	dir[1] = 0
	dir[2] = dir0[2] - dir1[2]*len0*0.5

	dtVnormalize(dir)
}

func calcStraightSteerDirection(ag *dtCrowdAgent, dir []float64) {
	if ag.ncorners == 0 {
		dtVset(dir, 0, 0, 0)
		return
	}
	copy(dir, dtVsub(ag.cornerVerts[:], ag.npos[:]))
	dir[1] = 0
	dtVnormalize(dir)
}

func addNeighbour(idx int, dist float64,
	neis []*dtCrowdNeighbour, nneis int, maxNeis int) int {
	// Insert neighbour based on the distance.
	var nei *dtCrowdNeighbour
	if nneis == 0 {
		nei = neis[nneis]
	} else if dist >= neis[nneis-1].dist {
		if nneis >= maxNeis {
			return nneis
		}

		nei = neis[nneis]
	} else {
		var i int
		for i = 0; i < nneis; i++ {
			if dist <= neis[i].dist {
				break
			}
		}

		tgt := i + 1
		n := dtMin(nneis-i, maxNeis-tgt)

		dtAssert(tgt+n <= maxNeis)

		if n > 0 {
			copy(neis[tgt:], neis[i:i+n])
		}

		nei = neis[i]
	}
	nei.idx = idx
	nei.dist = dist

	return dtMin(nneis+1, maxNeis)
}

func getNeighbours(pos []float64, height float64, rangef float64,
	skip *dtCrowdAgent, result []*dtCrowdNeighbour, maxResult int,
	agents []*dtCrowdAgent, nagents int, grid *dtProximityGrid) int {
	n := 0

	const MAX_NEIS = 32
	ids := make([]int, MAX_NEIS)
	nids := grid.queryItems(pos[0]-rangef, pos[2]-rangef,
		pos[0]+rangef, pos[2]+rangef,
		ids, MAX_NEIS)

	for i := 0; i < nids; i++ {
		ag := agents[ids[i]]

		if ag == skip {
			continue
		}

		// Check for overlap.

		diff := dtVsub(pos, ag.npos[:])
		if math.Abs(diff[1]) >= (height+ag.params.height)/2.0 {
			continue
		}

		diff[1] = 0
		distSqr := dtVlenSqr(diff)
		if distSqr > dtSqr(rangef) {
			continue
		}

		n = addNeighbour(ids[i], distSqr, result, n, maxResult)
	}
	return n
}

func addToOptQueue(newag *dtCrowdAgent, agents []*dtCrowdAgent, nagents, maxAgents int) int {
	// Insert neighbour based on greatest time.
	slot := 0
	if nagents == 0 {
		slot = nagents
	} else if newag.topologyOptTime <= agents[nagents-1].topologyOptTime {
		if nagents >= maxAgents {
			return nagents
		}

		slot = nagents
	} else {
		var i int
		for i := 0; i < nagents; i++ {
			if newag.topologyOptTime >= agents[i].topologyOptTime {
				break
			}

		}

		tgt := i + 1
		n := dtMin(nagents-i, maxAgents-tgt)

		dtAssert(tgt+n <= maxAgents)

		if n > 0 {
			copy(agents[tgt:], agents[i:i+n])
		}

		slot = i
	}

	agents[slot] = newag

	return dtMin(nagents+1, maxAgents)
}

func addToPathQueue(newag *dtCrowdAgent, agents []*dtCrowdAgent, nagents, maxAgents int) int {
	// Insert neighbour based on greatest time.
	slot := 0
	if nagents == 0 {
		slot = nagents
	} else if newag.targetReplanTime <= agents[nagents-1].targetReplanTime {
		if nagents >= maxAgents {
			return nagents
		}

		slot = nagents
	} else {
		var i int
		for i = 0; i < nagents; i++ {
			if newag.targetReplanTime >= agents[i].targetReplanTime {
				break
			}

		}

		tgt := i + 1
		n := dtMin(nagents-i, maxAgents-tgt)

		dtAssert(tgt+n <= maxAgents)

		if n > 0 {
			copy(agents[tgt:], agents[i:i+n])
		}

		slot = i
	}

	agents[slot] = newag

	return dtMin(nagents+1, maxAgents)
}

func (d *dtCrowd) purge() {

	d.m_agents = nil
	d.m_maxAgents = 0

	d.m_activeAgents = nil
	d.m_activeAgents = nil

	d.m_agentAnims = nil

	d.m_pathResult = nil

	d.m_grid = nil

	d.m_obstacleQuery = nil

	d.m_navquery = nil
}

func (d *dtCrowd) setObstacleAvoidanceParams(idx int, params *dtObstacleAvoidanceParams) {
	if idx >= 0 && idx < DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS {
		d.m_obstacleQueryParams[idx] = params
	}

}

func (d *dtCrowd) getObstacleAvoidanceParams(idx int) *dtObstacleAvoidanceParams {
	if idx >= 0 && idx < DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS {
		return d.m_obstacleQueryParams[idx]
	}

	return nil
}

func (d *dtCrowd) getAgentCount() int {
	return d.m_maxAgents
}

// / @par
// /
// / Agents in the pool may not be in use.  Check #dtCrowdAgent.active before using the returned object.
func (d *dtCrowd) getAgent(idx int) *dtCrowdAgent {
	if idx < 0 || idx >= d.m_maxAgents {
		return nil
	}

	return d.m_agents[idx]
}

// /
// / Agents in the pool may not be in use.  Check #dtCrowdAgent.active before using the returned object.
func (d *dtCrowd) getEditableAgent(idx int) *dtCrowdAgent {
	if idx < 0 || idx >= d.m_maxAgents {
		return nil
	}

	return d.m_agents[idx]
}

func (d *dtCrowd) updateAgentParameters(idx int, params *dtCrowdAgentParams) {
	if idx < 0 || idx >= d.m_maxAgents {
		return
	}
	d.m_agents[idx].params = params

}

// / @par
// /
// / The agent's position will be constrained to the surface of the navigation mesh.
func (d *dtCrowd) addAgent(pos []float64, params *dtCrowdAgentParams) int {
	// Find empty slot.
	idx := -1
	for i := 0; i < d.m_maxAgents; i++ {
		if !d.m_agents[i].active {
			idx = i
			break
		}
	}
	if idx == -1 {
		return -1
	}

	ag := d.m_agents[idx]

	d.updateAgentParameters(idx, params)

	// Find nearest position on navmesh and place the agent there.
	nearest := make([]float64, 3)
	var ref dtPolyRef
	copy(nearest, pos)
	ref, status := d.m_navquery.findNearestPoly(pos, d.m_agentPlacementHalfExtents[:], d.m_filters[ag.params.queryFilterType], nearest)
	if status.dtStatusFailed() {
		copy(nearest, pos)
		ref = 0
	}

	ag.corridor.reset(ref, nearest)
	ag.boundary.reset()
	ag.partial = false

	ag.topologyOptTime = 0
	ag.targetReplanTime = 0
	ag.nneis = 0

	dtVset(ag.dvel[:], 0, 0, 0)
	dtVset(ag.nvel[:], 0, 0, 0)
	dtVset(ag.vel[:], 0, 0, 0)
	copy(ag.npos[:], nearest)

	ag.desiredSpeed = 0

	if ref != 0 {
		ag.state = DT_CROWDAGENT_STATE_WALKING
	} else {
		ag.state = DT_CROWDAGENT_STATE_INVALID
	}

	ag.targetState = DT_CROWDAGENT_TARGET_NONE

	ag.active = true

	return idx
}

// / @par
// /
// / The agent is deactivated and will no longer be processed.  Its #dtCrowdAgent object
// / is not removed from the pool.  It is marked as inactive so that it is available for reuse.
func (d *dtCrowd) removeAgent(idx int) {
	if idx >= 0 && idx < d.m_maxAgents {
		d.m_agents[idx].active = false
	}
}

func (d *dtCrowd) requestMoveTargetReplan(idx int, ref dtPolyRef, pos []float64) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.targetRef = ref
	copy(ag.targetPos[:], pos)
	ag.targetPathqRef = DT_PATHQ_INVALID
	ag.targetReplan = true
	if ag.targetRef != 0 {
		ag.targetState = DT_CROWDAGENT_TARGET_REQUESTING
	} else {
		ag.targetState = DT_CROWDAGENT_TARGET_FAILED
	}

	return true
}

// / @par
// /
// / May be called more than once to purge and re-initialize the crowd.
func newDtCrowd(maxAgents int, maxAgentRadius float64, nav *dtNavMesh) *dtCrowd {
	d := &dtCrowd{}
	d.purge()

	d.m_maxAgents = maxAgents
	d.m_maxAgentRadius = maxAgentRadius

	// Larger than agent radius because it is also used for agent recovery.
	dtVset(d.m_agentPlacementHalfExtents[:], d.m_maxAgentRadius*2.0, d.m_maxAgentRadius*1.5, d.m_maxAgentRadius*2.0)

	d.m_grid = newDtProximityGrid(d.m_maxAgents*4, maxAgentRadius*3)
	d.m_obstacleQuery = &dtObstacleAvoidanceQuery{}
	if !d.m_obstacleQuery.init(6, 8) {
		return nil
	}

	// Init obstacle query params.
	for i := 0; i < DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS; i++ {
		params := d.m_obstacleQueryParams[i]
		params.velBias = 0.4
		params.weightDesVel = 2.0
		params.weightCurVel = 0.75
		params.weightSide = 0.75
		params.weightToi = 2.5
		params.horizTime = 2.5
		params.gridSize = 33
		params.adaptiveDivs = 7
		params.adaptiveRings = 2
		params.adaptiveDepth = 5
	}

	// Allocate temp buffer for merging paths.
	d.m_maxPathResult = 256
	d.m_pathResult = make([]dtPolyRef, d.m_maxPathResult)

	if !d.m_pathq.init(d.m_maxPathResult, MAX_PATHQUEUE_NODES, nav) {
		return nil
	}

	d.m_agents = make([]*dtCrowdAgent, d.m_maxAgents)

	d.m_activeAgents = make([]*dtCrowdAgent, d.m_maxAgents)

	d.m_agentAnims = make([]*dtCrowdAgentAnimation, d.m_maxAgents)

	for i := 0; i < d.m_maxAgents; i++ {
		d.m_agents[i] = &dtCrowdAgent{}
		d.m_agents[i].active = false
		d.m_agents[i].corridor = newDtPathCorridor(d.m_maxPathResult)

	}

	for i := 0; i < d.m_maxAgents; i++ {
		d.m_agentAnims[i].active = false
	}

	// The navquery is mostly used for local searches, no need for large node pool.
	d.m_navquery = NewDtNavMeshQuery(nav, MAX_COMMON_NODES)
	return d
}

// / @par
// /
// / This method is used when a new target is set.
// /
// / The position will be constrained to the surface of the navigation mesh.
// /
// / The request will be processed during the next #update().
func (d *dtCrowd) requestMoveTarget(idx int, ref dtPolyRef, pos []float64) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	if ref == 0 {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.targetRef = ref
	copy(ag.targetPos[:], pos)
	ag.targetPathqRef = DT_PATHQ_INVALID
	ag.targetReplan = false
	if ag.targetRef != 0 {
		ag.targetState = DT_CROWDAGENT_TARGET_REQUESTING
	} else {
		ag.targetState = DT_CROWDAGENT_TARGET_FAILED
	}

	return true
}

func (d *dtCrowd) requestMoveVelocity(idx int, vel []float64) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.targetRef = 0
	copy(ag.targetPos[:], vel)
	ag.targetPathqRef = DT_PATHQ_INVALID
	ag.targetReplan = false
	ag.targetState = DT_CROWDAGENT_TARGET_VELOCITY

	return true
}

func (d *dtCrowd) resetMoveTarget(idx int) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.targetRef = 0
	dtVset(ag.targetPos[:], 0, 0, 0)
	dtVset(ag.dvel[:], 0, 0, 0)
	ag.targetPathqRef = DT_PATHQ_INVALID
	ag.targetReplan = false
	ag.targetState = DT_CROWDAGENT_TARGET_NONE

	return true
}

func (d *dtCrowd) getActiveAgents(agents []*dtCrowdAgent, maxAgents int) int {
	n := 0
	for i := 0; i < d.m_maxAgents; i++ {
		if !d.m_agents[i].active {
			continue
		}
		if n < maxAgents {
			agents[n] = d.m_agents[i]
			n++
		}
	}
	return n
}

func (d *dtCrowd) updateMoveRequest(dt float64) {
	PATH_MAX_AGENTS := 8
	queue := make([]*dtCrowdAgent, PATH_MAX_AGENTS)
	nqueue := 0

	// Fire off new requests.
	for i := 0; i < d.m_maxAgents; i++ {
		ag := d.m_agents[i]
		if !ag.active {
			continue
		}

		if ag.state == DT_CROWDAGENT_STATE_INVALID {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_REQUESTING {
			path := ag.corridor.getPath()
			npath := ag.corridor.getPathCount()
			dtAssert(npath)

			MAX_RES := 32
			reqPos := make([]float64, 3)
			reqPath := make([]dtPolyRef, MAX_RES) // The path to the request location
			reqPathCount := 0

			// Quick search towards the goal.
			MAX_ITER := 20
			d.m_navquery.initSlicedFindPath(path[0], ag.targetRef, ag.npos[:], ag.targetPos[:], d.m_filters[ag.params.queryFilterType], 0)
			var status dtStatus
			d.m_navquery.updateSlicedFindPath(MAX_ITER)

			if ag.targetReplan { // && npath > 10)
				// Try to use existing steady path during replan if possible.
				reqPathCount, status = d.m_navquery.finalizeSlicedFindPathPartial(path, npath, reqPath, MAX_RES)
			} else {
				// Try to move towards target when goal changes.
				reqPathCount, status = d.m_navquery.finalizeSlicedFindPath(reqPath, MAX_RES)
			}

			if !status.dtStatusFailed() && reqPathCount > 0 {
				// In progress or succeed.
				if reqPath[reqPathCount-1] != ag.targetRef {
					// Partial path, constrain target position inside the last polygon.
					reqPos, _, status = d.m_navquery.closestPointOnPoly(reqPath[reqPathCount-1], ag.targetPos[:])
					if status.dtStatusFailed() {
						reqPathCount = 0
					}

				} else {
					copy(reqPos, ag.targetPos[:])
				}
			} else {
				reqPathCount = 0
			}

			if reqPathCount == 0 {
				// Could not find path, start the request from current location.
				copy(reqPos, ag.npos[:])
				reqPath[0] = path[0]
				reqPathCount = 1
			}

			ag.corridor.setCorridor(reqPos, reqPath, reqPathCount)
			ag.boundary.reset()
			ag.partial = false

			if reqPath[reqPathCount-1] == ag.targetRef {
				ag.targetState = DT_CROWDAGENT_TARGET_VALID
				ag.targetReplanTime = 0.0
			} else {
				// The path is longer or potentially unreachable, full plan.
				ag.targetState = DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE
			}
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE {
			nqueue = addToPathQueue(ag, queue, nqueue, PATH_MAX_AGENTS)
		}
	}

	for i := 0; i < nqueue; i++ {
		ag := queue[i]
		tmp := ag.corridor.getTarget()
		ag.targetPathqRef = d.m_pathq.request(ag.corridor.getLastPoly(), ag.targetRef,
			tmp[:], ag.targetPos[:], d.m_filters[ag.params.queryFilterType])
		if ag.targetPathqRef != DT_PATHQ_INVALID {
			ag.targetState = DT_CROWDAGENT_TARGET_WAITING_FOR_PATH
		}

	}

	// Update requests.
	d.m_pathq.update(MAX_ITERS_PER_UPDATE)

	var status dtStatus

	// Process path results.
	for i := 0; i < d.m_maxAgents; i++ {
		ag := d.m_agents[i]
		if !ag.active {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_WAITING_FOR_PATH {
			// Poll path queue.
			status = d.m_pathq.getRequestStatus(ag.targetPathqRef)
			if status.dtStatusFailed() {
				// Path find failed, retry if the target location is still valid.
				ag.targetPathqRef = DT_PATHQ_INVALID
				if ag.targetRef != 0 {
					ag.targetState = DT_CROWDAGENT_TARGET_REQUESTING
				} else {
					ag.targetState = DT_CROWDAGENT_TARGET_FAILED
				}

				ag.targetReplanTime = 0.0
			} else if status.dtStatusSucceed() {
				path := ag.corridor.getPath()
				npath := ag.corridor.getPathCount()
				dtAssert(npath)

				// Apply results.
				targetPos := make([]float64, 3)
				copy(targetPos, ag.targetPos[:])

				res := d.m_pathResult
				valid := true
				nres := 0
				status = d.m_pathq.getPathResult(ag.targetPathqRef, res, &nres, d.m_maxPathResult)
				if status.dtStatusFailed() || nres == 0 {
					valid = false
				}

				if status.dtStatusDetail(DT_PARTIAL_RESULT) {
					ag.partial = true
				} else {
					ag.partial = false
				}

				// Merge result and existing path.
				// The agent might have moved whilst the request is
				// being processed, so the path may have changed.
				// We assume that the end of the path is at the same location
				// where the request was issued.

				// The last ref in the old path should be the same as
				// the location where the request was issued..
				if valid && path[npath-1] != res[0] {
					valid = false
				}

				if valid {
					// Put the old path infront of the old path.
					if npath > 1 {
						// Make space for the old path.
						if (npath-1)+nres > d.m_maxPathResult {
							nres = d.m_maxPathResult - (npath - 1)
						}

						copy(res[:npath-1], res[:nres])
						// Copy old path in the beginning.
						copy(res, path[:npath-1])
						nres += npath - 1

						// Remove trackbacks
						for j := 0; j < nres; j++ {
							if j-1 >= 0 && j+1 < nres {
								if res[j-1] == res[j+1] {
									copy(res[j-1:], res[j+1:j+1+nres-(j+1)])
									nres -= 2
									j -= 2
								}
							}
						}

					}

					// Check for partial path.
					if res[nres-1] != ag.targetRef {
						// Partial path, constrain target position inside the last polygon.
						var status dtStatus
						nearest := make([]float64, 3)
						nearest, _, status = d.m_navquery.closestPointOnPoly(res[nres-1], targetPos)
						if status.dtStatusSucceed() {
							copy(targetPos, nearest)
						} else {
							valid = false
						}

					}
				}

				if valid {
					// Set current corridor.
					ag.corridor.setCorridor(targetPos, res, nres)
					// Force to update boundary.
					ag.boundary.reset()
					ag.targetState = DT_CROWDAGENT_TARGET_VALID
				} else {
					// Something went wrong.
					ag.targetState = DT_CROWDAGENT_TARGET_FAILED
				}

				ag.targetReplanTime = 0.0
			}
		}
	}

}

func (d *dtCrowd) updateTopologyOptimization(agents []*dtCrowdAgent, nagents int, dt float64) {
	if nagents == 0 {
		return
	}

	OPT_TIME_THR := 0.5 // seconds
	const OPT_MAX_AGENTS = 1
	queue := [OPT_MAX_AGENTS]*dtCrowdAgent{}
	nqueue := 0

	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		if (ag.params.updateFlags & DT_CROWD_OPTIMIZE_TOPO) == 0 {
			continue
		}

		ag.topologyOptTime += dt
		if ag.topologyOptTime >= OPT_TIME_THR {
			nqueue = addToOptQueue(ag, queue[:], nqueue, OPT_MAX_AGENTS)
		}

	}

	for i := 0; i < nqueue; i++ {
		ag := queue[i]
		ag.corridor.optimizePathTopology(d.m_navquery, d.m_filters[ag.params.queryFilterType])
		ag.topologyOptTime = 0
	}

}

func (d *dtCrowd) checkPathValidity(agents []*dtCrowdAgent, nagents int, dt float64) {
	const CHECK_LOOKAHEAD = 10
	const TARGET_REPLAN_DELAY = 1.0 // seconds

	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		ag.targetReplanTime += dt

		replan := false

		// First check that the current location is valid.
		idx := d.getAgentIndex(ag)
		agentPos := make([]float64, 3)
		agentRef := ag.corridor.getFirstPoly()
		copy(agentPos, ag.npos[:])
		if !d.m_navquery.isValidPolyRef(agentRef, d.m_filters[ag.params.queryFilterType]) {
			// Current location is not valid, try to reposition.
			// TODO: this can snap agents, how to handle that?
			nearest := make([]float64, 3)
			copy(nearest, agentPos)
			agentRef = 0
			agentRef, _ = d.m_navquery.findNearestPoly(ag.npos[:], d.m_agentPlacementHalfExtents[:], d.m_filters[ag.params.queryFilterType], nearest)
			copy(agentPos, nearest)

			if agentRef == 0 {
				// Could not find location in navmesh, set state to invalid.
				ag.corridor.reset(0, agentPos)
				ag.partial = false
				ag.boundary.reset()
				ag.state = DT_CROWDAGENT_STATE_INVALID
				continue
			}

			// Make sure the first polygon is valid, but leave other valid
			// polygons in the path so that replanner can adjust the path better.
			ag.corridor.fixPathStart(agentRef, agentPos)
			//			ag->corridor.trimInvalidPath(agentRef, agentPos, m_navquery, &m_filter);
			ag.boundary.reset()
			copy(ag.npos[:], agentPos)

			replan = true
		}

		// If the agent does not have move target or is controlled by velocity, no need to recover the target nor replan.
		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		// Try to recover move request position.
		if ag.targetState != DT_CROWDAGENT_TARGET_NONE && ag.targetState != DT_CROWDAGENT_TARGET_FAILED {
			if !d.m_navquery.isValidPolyRef(ag.targetRef, d.m_filters[ag.params.queryFilterType]) {
				// Current target is not valid, try to reposition.
				nearest := make([]float64, 3)
				copy(nearest, ag.targetPos[:])
				ag.targetRef = 0
				ag.targetRef, _ = d.m_navquery.findNearestPoly(ag.targetPos[:], d.m_agentPlacementHalfExtents[:], d.m_filters[ag.params.queryFilterType], nearest)
				copy(ag.targetPos[:], nearest)
				replan = true
			}
			if ag.targetRef == 0 {
				// Failed to reposition target, fail moverequest.
				ag.corridor.reset(agentRef, agentPos)
				ag.partial = false
				ag.targetState = DT_CROWDAGENT_TARGET_NONE
			}
		}

		// If nearby corridor is not valid, replan.
		if !ag.corridor.isValid(CHECK_LOOKAHEAD, d.m_navquery, d.m_filters[ag.params.queryFilterType]) {
			// Fix current path.
			//			ag->corridor.trimInvalidPath(agentRef, agentPos, m_navquery, &m_filter);
			//			ag->boundary.reset();
			replan = true
		}

		// If the end of the path is near and it is not the requested location, replan.
		if ag.targetState == DT_CROWDAGENT_TARGET_VALID {
			if ag.targetReplanTime > TARGET_REPLAN_DELAY && ag.corridor.getPathCount() < CHECK_LOOKAHEAD && ag.corridor.getLastPoly() != ag.targetRef {
				replan = true
			}

		}

		// Try to replan path to goal.
		if replan {
			if ag.targetState != DT_CROWDAGENT_TARGET_NONE {
				d.requestMoveTargetReplan(idx, ag.targetRef, ag.targetPos[:])
			}
		}
	}
}

func (d *dtCrowd) update(dt float64, debug *dtCrowdAgentDebugInfo) {
	d.m_velocitySampleCount = 0

	debugIdx := -1
	if debug != nil {
		debugIdx = debug.idx
	}

	agents := d.m_activeAgents
	nagents := d.getActiveAgents(agents, d.m_maxAgents)

	// Check that all agents still have valid paths.
	d.checkPathValidity(agents, nagents, dt)

	// Update async move request and path finder.
	d.updateMoveRequest(dt)

	// Optimize path topology.
	d.updateTopologyOptimization(agents, nagents, dt)

	// Register agents to proximity grid.
	d.m_grid.Clear()
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		p := ag.npos
		r := ag.params.radius
		d.m_grid.addItem(i, p[0]-r, p[2]-r, p[0]+r, p[2]+r)
	}

	// Get nearby navmesh segments and agents to collide with.
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		// Update the collision boundary after certain distance has been passed or
		// if it has become invalid.
		updateThr := ag.params.collisionQueryRange * 0.25
		tmp := ag.boundary.getCenter()
		if dtVdist2DSqr(ag.npos[:], tmp[:]) > dtSqr(updateThr) ||
			!ag.boundary.isValid(d.m_navquery, d.m_filters[ag.params.queryFilterType]) {
			ag.boundary.update(ag.corridor.getFirstPoly(), ag.npos[:], ag.params.collisionQueryRange,
				d.m_navquery, d.m_filters[ag.params.queryFilterType])
		}
		// Query neighbour agents
		ag.nneis = getNeighbours(ag.npos[:], ag.params.height, ag.params.collisionQueryRange,
			ag, ag.neis[:], DT_CROWDAGENT_MAX_NEIGHBOURS,
			agents, nagents, d.m_grid)
		for j := 0; j < ag.nneis; j++ {
			ag.neis[j].idx = d.getAgentIndex(agents[ag.neis[j].idx])
		}

	}

	// Find next corner to steer to.
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		// Find corners for steering
		ag.ncorners = ag.corridor.findCorners(ag.cornerVerts[:], ag.cornerFlags[:], ag.cornerPolys[:],
			DT_CROWDAGENT_MAX_CORNERS, d.m_navquery, d.m_filters[ag.params.queryFilterType])

		// Check to see if the corner after the next corner is directly visible,
		// and short cut to there.
		if (ag.params.updateFlags&DT_CROWD_OPTIMIZE_VIS > 0) && ag.ncorners > 0 {
			target := rcGetVert(ag.cornerVerts[:], dtMin(1, ag.ncorners-1))
			ag.corridor.optimizePathVisibility(target, ag.params.pathOptimizationRange, d.m_navquery, d.m_filters[ag.params.queryFilterType])

			// Copy data for debug purposes.
			if debugIdx == i {
				tmp := ag.corridor.getPos()
				copy(debug.optStart[:], tmp[:])
				copy(debug.optEnd[:], target)
			}
		} else {
			// Copy data for debug purposes.
			if debugIdx == i {
				dtVset(debug.optStart[:], 0, 0, 0)
				dtVset(debug.optEnd[:], 0, 0, 0)
			}
		}
	}

	// Trigger off-mesh connections (depends on corners).
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		// Check
		triggerRadius := ag.params.radius * 2.25
		if overOffmeshConnection(ag, triggerRadius) {
			// Prepare to off-mesh connection.
			idx := d.getAgentIndex(ag)
			anim := d.m_agentAnims[idx]

			// Adjust the path over the off-mesh connection.
			refs := make([]dtPolyRef, 2)
			if ag.corridor.moveOverOffmeshConnection(ag.cornerPolys[ag.ncorners-1], refs,
				anim.startPos[:], anim.endPos[:], d.m_navquery) {
				copy(anim.initPos[:], ag.npos[:])
				anim.polyRef = refs[1]
				anim.active = true
				anim.t = 0.0
				anim.tmax = (dtVdist2D(anim.startPos[:], anim.endPos[:]) / ag.params.maxSpeed) * 0.5

				ag.state = DT_CROWDAGENT_STATE_OFFMESH
				ag.ncorners = 0
				ag.nneis = 0
				continue
			} else {
				// Path validity check will ensure that bad/blocked connections will be replanned.
			}
		}
	}

	// Calculate steering.
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.targetState == DT_CROWDAGENT_TARGET_NONE {
			continue
		}

		dvel := []float64{0, 0, 0}

		if ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			copy(dvel, ag.targetPos[:])
			ag.desiredSpeed = dtVlen(ag.targetPos[:])
		} else {
			// Calculate steering direction.
			if ag.params.updateFlags&DT_CROWD_ANTICIPATE_TURNS > 0 {
				calcSmoothSteerDirection(ag, dvel)
			} else {
				calcStraightSteerDirection(ag, dvel)
			}

			// Calculate speed scale, which tells the agent to slowdown at the end of the path.
			slowDownRadius := ag.params.radius * 2 // TODO: make less hacky.
			speedScale := getDistanceToGoal(ag, slowDownRadius) / slowDownRadius

			ag.desiredSpeed = ag.params.maxSpeed
			copy(dvel, dtVscale(dvel, ag.desiredSpeed*speedScale))

		}

		// Separation
		if ag.params.updateFlags&DT_CROWD_SEPARATION > 0 {
			separationDist := ag.params.collisionQueryRange
			invSeparationDist := 1.0 / separationDist
			separationWeight := ag.params.separationWeight

			w := float64(0)
			disp := []float64{0, 0, 0}

			for j := 0; j < ag.nneis; j++ {
				nei := d.m_agents[ag.neis[j].idx]

				diff := dtVsub(ag.npos[:], nei.npos[:])
				diff[1] = 0

				distSqr := dtVlenSqr(diff)
				if distSqr < 0.00001 {
					continue
				}

				if distSqr > dtSqr(separationDist) {
					continue
				}

				dist := math.Sqrt(distSqr)
				weight := separationWeight * (1.0 - dtSqr(dist*invSeparationDist))

				dtVmad(disp, disp, diff, weight/dist)
				w += 1.0
			}

			if w > 0.0001 {
				// Adjust desired velocity.
				dtVmad(dvel, dvel, disp, 1.0/w)
				// Clamp desired velocity to desired speed.
				speedSqr := dtVlenSqr(dvel)
				desiredSqr := dtSqr(ag.desiredSpeed)
				if speedSqr > desiredSqr {
					copy(dvel, dtVscale(dvel, desiredSqr/speedSqr))
				}

			}
		}

		// Set the desired velocity.
		copy(ag.dvel[:], dvel[:])
	}

	// Velocity planning.
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.params.updateFlags&DT_CROWD_OBSTACLE_AVOIDANCE > 0 {
			d.m_obstacleQuery.reset()

			// Add neighbours as obstacles.
			for j := 0; j < ag.nneis; j++ {
				nei := d.m_agents[ag.neis[j].idx]
				d.m_obstacleQuery.addCircle(nei.npos[:], nei.params.radius, nei.vel[:], nei.dvel[:])
			}

			// Append neighbour segments as obstacles.
			for j := 0; j < ag.boundary.getSegmentCount(); j++ {
				s := ag.boundary.getSegment(j)
				if dtTriArea2D(ag.npos[:], s[:], s[3:]) < 0.0 {
					continue
				}

				d.m_obstacleQuery.addSegment(s[:], s[3:])
			}

			var vod *dtObstacleAvoidanceDebugData
			if debugIdx == i {
				vod = debug.vod
			}

			// Sample new safe velocity.
			adaptive := true
			ns := 0

			params := d.m_obstacleQueryParams[ag.params.obstacleAvoidanceType]

			if adaptive {
				ns = d.m_obstacleQuery.sampleVelocityAdaptive(ag.npos[:], ag.params.radius, ag.desiredSpeed,
					ag.vel[:], ag.dvel[:], ag.nvel[:], params, vod)
			} else {
				ns = d.m_obstacleQuery.sampleVelocityGrid(ag.npos[:], ag.params.radius, ag.desiredSpeed,
					ag.vel[:], ag.dvel[:], ag.nvel[:], params, vod)
			}
			d.m_velocitySampleCount += ns
		} else {
			// If not using velocity planning, new velocity is directly the desired velocity.
			copy(ag.nvel[:], ag.dvel[:])
		}
	}

	// Integrate.
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		integrate(ag, dt)
	}

	// Handle collisions.
	const COLLISION_RESOLVE_FACTOR = 0.7

	for iter := 0; iter < 4; iter++ {
		for i := 0; i < nagents; i++ {
			ag := agents[i]
			idx0 := d.getAgentIndex(ag)

			if ag.state != DT_CROWDAGENT_STATE_WALKING {
				continue
			}

			dtVset(ag.disp[:], 0, 0, 0)

			w := float64(0)

			for j := 0; j < ag.nneis; j++ {
				nei := d.m_agents[ag.neis[j].idx]
				idx1 := d.getAgentIndex(nei)

				diff := dtVsub(ag.npos[:], nei.npos[:])
				diff[1] = 0

				dist := dtVlenSqr(diff)
				if dist > dtSqr(ag.params.radius+nei.params.radius) {
					continue
				}

				dist = math.Sqrt(dist)
				pen := (ag.params.radius + nei.params.radius) - dist
				if dist < 0.0001 {
					// Agents on top of each other, try to choose diverging separation directions.
					if idx0 > idx1 {
						dtVset(diff, -ag.dvel[2], 0, ag.dvel[0])
					} else {
						dtVset(diff, ag.dvel[2], 0, -ag.dvel[0])
					}

					pen = 0.01
				} else {
					pen = (1.0 / dist) * (pen * 0.5) * COLLISION_RESOLVE_FACTOR
				}

				dtVmad(ag.disp[:], ag.disp[:], diff, pen)

				w += 1.0
			}

			if w > 0.0001 {
				iw := 1.0 / w
				copy(ag.disp[:], dtVscale(ag.disp[:], iw))
			}
		}

		for i := 0; i < nagents; i++ {
			ag := agents[i]
			if ag.state != DT_CROWDAGENT_STATE_WALKING {
				continue
			}

			copy(ag.npos[:], dtVadd(ag.npos[:], ag.disp[:]))

		}
	}

	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.state != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		// Move along navmesh.
		ag.corridor.movePosition(ag.npos[:], d.m_navquery, d.m_filters[ag.params.queryFilterType])
		// Get valid constrained position back.
		tmp := ag.corridor.getPos()
		copy(ag.npos[:], tmp[:])

		// If not using path, truncate the corridor to just one poly.
		if ag.targetState == DT_CROWDAGENT_TARGET_NONE || ag.targetState == DT_CROWDAGENT_TARGET_VELOCITY {
			ag.corridor.reset(ag.corridor.getFirstPoly(), ag.npos[:])
			ag.partial = false
		}

	}

	// Update agents using off-mesh connection.
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		idx := d.getAgentIndex(ag)
		anim := d.m_agentAnims[idx]
		if !anim.active {
			continue
		}

		anim.t += dt
		if anim.t > anim.tmax {
			// Reset animation
			anim.active = false
			// Prepare agent for walking.
			ag.state = DT_CROWDAGENT_STATE_WALKING
			continue
		}

		// Update position
		ta := anim.tmax * 0.15
		tb := anim.tmax
		if anim.t < ta {
			u := tween(anim.t, 0.0, ta)
			copy(ag.npos[:], dtVlerp(anim.initPos[:], anim.startPos[:], u))

		} else {
			u := tween(anim.t, ta, tb)
			copy(ag.npos[:], dtVlerp(anim.startPos[:], anim.endPos[:], u))

		}

		// Update velocity.
		dtVset(ag.vel[:], 0, 0, 0)
		dtVset(ag.dvel[:], 0, 0, 0)
	}

}
