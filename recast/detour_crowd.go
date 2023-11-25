package recast

import (
	"gonavamesh/common"
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
// / @see DtObstacleAvoidanceParams, DtCrowd::setObstacleAvoidanceParams(), DtCrowd::getObstacleAvoidanceParams(),
// /		 DtCrowdAgentParams::obstacleAvoidanceType
const DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS = 8

// / The maximum number of query filter types supported by the crowd manager.
// / @ingroup crowd
// / @see DtQueryFilter, DtCrowd::getFilter() DtCrowd::getEditableFilter(),
// /		DtCrowdAgentParams::queryFilterType
const DT_CROWD_MAX_QUERY_FILTER_TYPE = 16

// / Provides neighbor data for agents managed by the crowd.
// / @ingroup crowd
// / @see DtCrowdAgent::neis, DtCrowd
type DtCrowdNeighbour struct {
	Idx  int     ///< The index of the neighbor in the crowd.
	Dist float64 ///< The distance between the current agent and the neighbor.
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
type DtCrowdAgentParams struct {
	Radius          float64 ///< Agent radius. [Limit: >= 0]
	Height          float64 ///< Agent height. [Limit: > 0]
	MaxAcceleration float64 ///< Maximum allowed acceleration. [Limit: >= 0]
	MaxSpeed        float64 ///< Maximum allowed speed. [Limit: >= 0]

	/// Defines how close a collision element must be before it is considered for steering behaviors. [Limits: > 0]
	CollisionQueryRange float64

	PathOptimizationRange float64 ///< The path visibility optimization range. [Limit: > 0]

	/// How aggresive the agent manager should be at avoiding collisions with this agent. [Limit: >= 0]
	SeparationWeight float64

	/// Flags that impact steering behavior. (See: #UpdateFlags)
	UpdateFlags int

	/// The index of the avoidance configuration to use for the agent.
	/// [Limits: 0 <= value <= #DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS]
	ObstacleAvoidanceType int

	/// The index of the query filter used by this agent.
	QueryFilterType int

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

// / Represents an agent managed by a #DtCrowd object.
// / @ingroup crowd
type DtCrowdAgent struct {
	/// True if the agent is active, false if the agent is in an unused slot in the agent pool.
	Active bool

	/// The type of mesh polygon the agent is traversing. (See: #CrowdAgentState)
	State int

	/// True if the agent has valid path (targetState == DT_CROWDAGENT_TARGET_VALID) and the path does not lead to the requested position, else false.
	Partial bool

	/// The path corridor the agent is using.
	Corridor *dtPathCorridor

	/// The local boundary data for the agent.
	Boundary *dtLocalBoundary

	/// Time since the agent's path corridor was optimized.
	TopologyOptTime float64

	/// The known neighbors of the agent.
	Neis [DT_CROWDAGENT_MAX_NEIGHBOURS]*DtCrowdNeighbour

	/// The number of neighbors.
	Nneis int

	/// The desired speed.
	desiredSpeed float64

	Npos [3]float64 ///< The current agent position. [(x, y, z)]
	Disp [3]float64 ///< A temporary value used to accumulate agent displacement during iterative collision resolution. [(x, y, z)]
	Dvel [3]float64 ///< The desired velocity of the agent. Based on the current path, calculated from scratch each frame. [(x, y, z)]
	Nvel [3]float64 ///< The desired velocity adjusted by obstacle avoidance, calculated from scratch each frame. [(x, y, z)]
	Vel  [3]float64 ///< The actual velocity of the agent. The change from nvel -> vel is constrained by max acceleration. [(x, y, z)]

	/// The agent's configuration parameters.
	Params *DtCrowdAgentParams

	/// The local path corridor corners for the agent. (Staight path.) [(x, y, z) * #ncorners]
	CornerVerts [DT_CROWDAGENT_MAX_CORNERS * 3]float64

	/// The local path corridor corner flags. (See: #dtStraightPathFlags) [(flags) * #ncorners]
	CornerFlags [DT_CROWDAGENT_MAX_CORNERS]int

	/// The reference id of the polygon being entered at the corner. [(polyRef) * #ncorners]
	CornerPolys [DT_CROWDAGENT_MAX_CORNERS]DtPolyRef

	/// The number of corners.
	Ncorners int

	TargetState      int            ///< State of the movement request.
	TargetRef        DtPolyRef      ///< Target polyref of the movement request.
	TargetPos        [3]float64     ///< Target position of the movement request (or velocity in case of DT_CROWDAGENT_TARGET_VELOCITY).
	TargetPathqRef   dtPathQueueRef ///< Path finder ref.
	TargetReplan     bool           ///< Flag indicating that the current path is being replanned.
	TargetReplanTime float64        /// <Time since the agent's target was replanned.
}

type DtCrowdAgentAnimation struct {
	Active                    bool
	InitPos, StartPos, EndPos [3]float64
	PolyRef                   DtPolyRef
	T, Tmax                   float64
}

const (
	DT_CROWD_ANTICIPATE_TURNS   = 1
	DT_CROWD_OBSTACLE_AVOIDANCE = 2
	DT_CROWD_SEPARATION         = 4
	DT_CROWD_OPTIMIZE_VIS       = 8  ///< Use #dtPathCorridor::optimizePathVisibility() to optimize the agent path.
	DT_CROWD_OPTIMIZE_TOPO      = 16 ///< Use dtPathCorridor::optimizePathTopology() to optimize the agent
)

type DtCrowdAgentDebugInfo struct {
	Idx              int
	OptStart, OptEnd [3]float64
	Vod              *DtObstacleAvoidanceDebugData
}

// / Provides local steering behaviors for a group of agents.
// / @ingroup crowd
type DtCrowd struct {
	m_maxAgents    int
	m_agents       []*DtCrowdAgent
	m_activeAgents []*DtCrowdAgent
	m_agentAnims   []*DtCrowdAgentAnimation

	m_pathq *dtPathQueue

	m_obstacleQueryParams []*DtObstacleAvoidanceParams
	m_obstacleQuery       *dtObstacleAvoidanceQuery

	m_grid *DtProximityGrid

	m_pathResult    []DtPolyRef
	m_maxPathResult int

	m_agentPlacementHalfExtents []float64

	m_filters []*DtQueryFilter

	m_maxAgentRadius float64

	m_velocitySampleCount int

	m_navquery NavMeshQuery
}

// / @par
// /
// / May be called more than once to purge and re-initialize the crowd.
func NewDtCrowd(maxAgents int, maxAgentRadius float64, nav IDtNavMesh) *DtCrowd {
	d := &DtCrowd{
		m_filters:                   make([]*DtQueryFilter, DT_CROWD_MAX_QUERY_FILTER_TYPE),
		m_agentPlacementHalfExtents: make([]float64, 3),
		m_obstacleQueryParams:       make([]*DtObstacleAvoidanceParams, DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS),
	}
	d.purge()

	d.m_maxAgents = maxAgents
	d.m_maxAgentRadius = maxAgentRadius

	// Larger than agent radius because it is also used for agent recovery.
	common.Vset(d.m_agentPlacementHalfExtents[:], d.m_maxAgentRadius*2.0, d.m_maxAgentRadius*1.5, d.m_maxAgentRadius*2.0)

	d.m_grid = newDtProximityGrid(d.m_maxAgents*4, maxAgentRadius*3)
	d.m_obstacleQuery = &dtObstacleAvoidanceQuery{}
	if !d.m_obstacleQuery.init(6, 8) {
		return nil
	}

	// Init obstacle query params.
	for i := 0; i < DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS; i++ {
		params := d.m_obstacleQueryParams[i]
		params.VelBias = 0.4
		params.WeightDesVel = 2.0
		params.WeightCurVel = 0.75
		params.WeightSide = 0.75
		params.WeightToi = 2.5
		params.HorizTime = 2.5
		params.GridSize = 33
		params.AdaptiveDivs = 7
		params.AdaptiveRings = 2
		params.AdaptiveDepth = 5
	}

	// Allocate temp buffer for merging paths.
	d.m_maxPathResult = 256
	d.m_pathResult = make([]DtPolyRef, d.m_maxPathResult)

	if !d.m_pathq.init(d.m_maxPathResult, MAX_PATHQUEUE_NODES, nav) {
		return nil
	}

	d.m_agents = make([]*DtCrowdAgent, d.m_maxAgents)

	d.m_activeAgents = make([]*DtCrowdAgent, d.m_maxAgents)

	d.m_agentAnims = make([]*DtCrowdAgentAnimation, d.m_maxAgents)

	for i := 0; i < d.m_maxAgents; i++ {
		d.m_agents[i] = &DtCrowdAgent{}
		d.m_agents[i].Active = false
		d.m_agents[i].Corridor = newDtPathCorridor(d.m_maxPathResult)

	}

	for i := 0; i < d.m_maxAgents; i++ {
		d.m_agentAnims[i].Active = false
	}

	// The navquery is mostly used for local searches, no need for large node pool.
	d.m_navquery = NewDtNavMeshQuery(nav, MAX_COMMON_NODES)
	return d
}

// / Gets the filter used by the crowd.
// / @return The filter used by the crowd.
func (d *DtCrowd) GetFilter(i int) *DtQueryFilter {
	if i >= 0 && i < DT_CROWD_MAX_QUERY_FILTER_TYPE {
		return d.m_filters[i]
	}
	return nil
}

// / Gets the filter used by the crowd.
// / @return The filter used by the crowd.
func (d *DtCrowd) GetEditableFilter(i int) *DtQueryFilter {
	if i >= 0 && i < DT_CROWD_MAX_QUERY_FILTER_TYPE {
		return d.m_filters[i]
	}
	return nil
}

// / Gets the search halfExtents [(x, y, z)] used by the crowd for query operations.
// / @return The search halfExtents used by the crowd. [(x, y, z)]
func (d *DtCrowd) GetQueryHalfExtents() []float64 { return d.m_agentPlacementHalfExtents }

// / Same as getQueryHalfExtents. Left to maintain backwards compatibility.
// / @return The search halfExtents used by the crowd. [(x, y, z)]
func (d *DtCrowd) GetQueryExtents() []float64 { return d.m_agentPlacementHalfExtents }

// / Gets the velocity sample count.
// / @return The velocity sample count.
func (d *DtCrowd) GetVelocitySampleCount() int { return d.m_velocitySampleCount }

// / Gets the crowd's proximity grid.
// / @return The crowd's proximity grid.
func (d *DtCrowd) GetGrid() *DtProximityGrid { return d.m_grid }

// / Gets the crowd's path request queue.
// / @return The crowd's path request queue.
func (d *DtCrowd) GetPathQueue() *dtPathQueue { return d.m_pathq }

// / Gets the query object used by the crowd.
func (d *DtCrowd) GetNavMeshQuery() NavMeshQuery { return d.m_navquery }
func (d *DtCrowd) GetAgentIndex(agent *DtCrowdAgent) int {
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
	return common.Clamp((t-t0)/(t1-t0), 0.0, 1.0)
}

func integrate(ag *DtCrowdAgent, dt float64) {
	// Fake dynamic constraint.
	maxDelta := ag.Params.MaxAcceleration * dt

	dv := common.Vsub(ag.Nvel[:], ag.Vel[:])
	ds := common.Vlen(dv)
	if ds > maxDelta {
		dv = common.Vscale(dv, maxDelta/ds)
	}
	copy(ag.Vel[:], common.Vadd(ag.Vel[:], dv))

	// Integrate
	if common.Vlen(ag.Vel[:]) > 0.0001 {
		common.Vmad(ag.Npos[:], ag.Npos[:], ag.Vel[:], dt)
	} else {
		common.Vset(ag.Vel[:], 0, 0, 0)
	}

}

func overOffmeshConnection(ag *DtCrowdAgent, radius float64) bool {
	if ag.Ncorners == 0 {
		return false
	}

	offMeshConnection := false
	if ag.CornerFlags[ag.Ncorners-1]&DT_STRAIGHTPATH_OFFMESH_CONNECTION > 0 {
		offMeshConnection = true
	}
	if offMeshConnection {
		distSq := common.Vdist2DSqr(ag.Npos[:], rcGetVert(ag.CornerVerts[:], ag.Ncorners-1))
		if distSq < radius*radius {
			return true
		}

	}

	return false
}

func getDistanceToGoal(ag *DtCrowdAgent, rangef float64) float64 {
	if ag.Ncorners == 0 {
		return rangef
	}

	endOfPath := false
	if ag.CornerFlags[ag.Ncorners-1]&DT_STRAIGHTPATH_END > 0 {
		endOfPath = true
	}
	if endOfPath {
		return common.Min(common.Vdist2D(ag.Npos[:], rcGetVert(ag.CornerVerts[:], ag.Ncorners-1)), rangef)
	}

	return rangef
}

func calcSmoothSteerDirection(ag *DtCrowdAgent, dir []float64) {
	if ag.Ncorners == 0 {
		common.Vset(dir, 0, 0, 0)
		return
	}

	ip0 := 0
	ip1 := common.Min(1, ag.Ncorners-1)
	p0 := rcGetVert(ag.CornerVerts[:], ip0)
	p1 := rcGetVert(ag.CornerVerts[:], ip1)

	dir0 := common.Vsub(p0, ag.Npos[:])
	dir1 := common.Vsub(p1, ag.Npos[:])
	dir0[1] = 0
	dir1[1] = 0

	len0 := common.Vlen(dir0)
	len1 := common.Vlen(dir1)
	if len1 > 0.001 {
		dir1 = common.Vscale(dir1, 1.0/len1)
	}

	dir[0] = dir0[0] - dir1[0]*len0*0.5
	dir[1] = 0
	dir[2] = dir0[2] - dir1[2]*len0*0.5

	common.Vnormalize(dir)
}

func calcStraightSteerDirection(ag *DtCrowdAgent, dir []float64) {
	if ag.Ncorners == 0 {
		common.Vset(dir, 0, 0, 0)
		return
	}
	copy(dir, common.Vsub(ag.CornerVerts[:], ag.Npos[:]))
	dir[1] = 0
	common.Vnormalize(dir)
}

func addNeighbour(idx int, dist float64,
	neis []*DtCrowdNeighbour, nneis int, maxNeis int) int {
	// Insert neighbour based on the distance.
	var nei *DtCrowdNeighbour
	if nneis == 0 {
		nei = neis[nneis]
	} else if dist >= neis[nneis-1].Dist {
		if nneis >= maxNeis {
			return nneis
		}

		nei = neis[nneis]
	} else {
		var i int
		for i = 0; i < nneis; i++ {
			if dist <= neis[i].Dist {
				break
			}
		}

		tgt := i + 1
		n := common.Min(nneis-i, maxNeis-tgt)

		common.AssertTrue(tgt+n <= maxNeis)

		if n > 0 {
			copy(neis[tgt:], neis[i:i+n])
		}

		nei = neis[i]
	}
	nei.Idx = idx
	nei.Dist = dist

	return common.Min(nneis+1, maxNeis)
}

func getNeighbours(pos []float64, height float64, rangef float64,
	skip *DtCrowdAgent, result []*DtCrowdNeighbour, maxResult int,
	agents []*DtCrowdAgent, nagents int, grid *DtProximityGrid) int {
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

		diff := common.Vsub(pos, ag.Npos[:])
		if math.Abs(diff[1]) >= (height+ag.Params.Height)/2.0 {
			continue
		}

		diff[1] = 0
		distSqr := common.VlenSqr(diff)
		if distSqr > common.Sqr(rangef) {
			continue
		}

		n = addNeighbour(ids[i], distSqr, result, n, maxResult)
	}
	return n
}

func addToOptQueue(newag *DtCrowdAgent, agents []*DtCrowdAgent, nagents, maxAgents int) int {
	// Insert neighbour based on greatest time.
	slot := 0
	if nagents == 0 {
		slot = nagents
	} else if newag.TopologyOptTime <= agents[nagents-1].TopologyOptTime {
		if nagents >= maxAgents {
			return nagents
		}

		slot = nagents
	} else {
		var i int
		for i := 0; i < nagents; i++ {
			if newag.TopologyOptTime >= agents[i].TopologyOptTime {
				break
			}

		}

		tgt := i + 1
		n := common.Min(nagents-i, maxAgents-tgt)

		common.AssertTrue(tgt+n <= maxAgents)

		if n > 0 {
			copy(agents[tgt:], agents[i:i+n])
		}

		slot = i
	}

	agents[slot] = newag

	return common.Min(nagents+1, maxAgents)
}

func addToPathQueue(newag *DtCrowdAgent, agents []*DtCrowdAgent, nagents, maxAgents int) int {
	// Insert neighbour based on greatest time.
	slot := 0
	if nagents == 0 {
		slot = nagents
	} else if newag.TargetReplanTime <= agents[nagents-1].TargetReplanTime {
		if nagents >= maxAgents {
			return nagents
		}

		slot = nagents
	} else {
		var i int
		for i = 0; i < nagents; i++ {
			if newag.TargetReplanTime >= agents[i].TargetReplanTime {
				break
			}

		}

		tgt := i + 1
		n := common.Min(nagents-i, maxAgents-tgt)

		common.AssertTrue(tgt+n <= maxAgents)

		if n > 0 {
			copy(agents[tgt:], agents[i:i+n])
		}

		slot = i
	}

	agents[slot] = newag

	return common.Min(nagents+1, maxAgents)
}

func (d *DtCrowd) purge() {

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

func (d *DtCrowd) SetObstacleAvoidanceParams(idx int, params *DtObstacleAvoidanceParams) {
	if idx >= 0 && idx < DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS {
		d.m_obstacleQueryParams[idx] = params
	}

}

func (d *DtCrowd) GetObstacleAvoidanceParams(idx int) *DtObstacleAvoidanceParams {
	if idx >= 0 && idx < DT_CROWD_MAX_OBSTAVOIDANCE_PARAMS {
		return d.m_obstacleQueryParams[idx]
	}

	return nil
}

func (d *DtCrowd) GetAgentCount() int {
	return d.m_maxAgents
}

// / @par
// /
// / Agents in the pool may not be in use.  Check #DtCrowdAgent.active before using the returned object.
func (d *DtCrowd) GetAgent(idx int) *DtCrowdAgent {
	if idx < 0 || idx >= d.m_maxAgents {
		return nil
	}

	return d.m_agents[idx]
}

// /
// / Agents in the pool may not be in use.  Check #DtCrowdAgent.active before using the returned object.
func (d *DtCrowd) getEditableAgent(idx int) *DtCrowdAgent {
	if idx < 0 || idx >= d.m_maxAgents {
		return nil
	}

	return d.m_agents[idx]
}

func (d *DtCrowd) UpdateAgentParameters(idx int, params *DtCrowdAgentParams) {
	if idx < 0 || idx >= d.m_maxAgents {
		return
	}
	d.m_agents[idx].Params = params

}

// / @par
// /
// / The agent's position will be constrained to the surface of the navigation mesh.
func (d *DtCrowd) AddAgent(pos []float64, params *DtCrowdAgentParams) int {
	// Find empty slot.
	idx := -1
	for i := 0; i < d.m_maxAgents; i++ {
		if !d.m_agents[i].Active {
			idx = i
			break
		}
	}
	if idx == -1 {
		return -1
	}

	ag := d.m_agents[idx]

	d.UpdateAgentParameters(idx, params)

	// Find nearest position on navmesh and place the agent there.
	nearest := make([]float64, 3)
	var ref DtPolyRef
	copy(nearest, pos)
	ref, status := d.m_navquery.FindNearestPoly(pos, d.m_agentPlacementHalfExtents[:], d.m_filters[ag.Params.QueryFilterType], nearest)
	if status.DtStatusFailed() {
		copy(nearest, pos)
		ref = 0
	}

	ag.Corridor.reset(ref, nearest)
	ag.Boundary.reset()
	ag.Partial = false

	ag.TopologyOptTime = 0
	ag.TargetReplanTime = 0
	ag.Nneis = 0

	common.Vset(ag.Dvel[:], 0, 0, 0)
	common.Vset(ag.Nvel[:], 0, 0, 0)
	common.Vset(ag.Vel[:], 0, 0, 0)
	copy(ag.Npos[:], nearest)

	ag.desiredSpeed = 0

	if ref != 0 {
		ag.State = DT_CROWDAGENT_STATE_WALKING
	} else {
		ag.State = DT_CROWDAGENT_STATE_INVALID
	}

	ag.TargetState = DT_CROWDAGENT_TARGET_NONE

	ag.Active = true

	return idx
}

// / @par
// /
// / The agent is deactivated and will no longer be processed.  Its #DtCrowdAgent object
// / is not removed from the pool.  It is marked as inactive so that it is available for reuse.
func (d *DtCrowd) RemoveAgent(idx int) {
	if idx >= 0 && idx < d.m_maxAgents {
		d.m_agents[idx].Active = false
	}
}

func (d *DtCrowd) requestMoveTargetReplan(idx int, ref DtPolyRef, pos []float64) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.TargetRef = ref
	copy(ag.TargetPos[:], pos)
	ag.TargetPathqRef = DT_PATHQ_INVALID
	ag.TargetReplan = true
	if ag.TargetRef != 0 {
		ag.TargetState = DT_CROWDAGENT_TARGET_REQUESTING
	} else {
		ag.TargetState = DT_CROWDAGENT_TARGET_FAILED
	}

	return true
}

// / @par
// /
// / This method is used when a new target is set.
// /
// / The position will be constrained to the surface of the navigation mesh.
// /
// / The request will be processed during the next #update().
func (d *DtCrowd) RequestMoveTarget(idx int, ref DtPolyRef, pos []float64) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	if ref == 0 {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.TargetRef = ref
	copy(ag.TargetPos[:], pos)
	ag.TargetPathqRef = DT_PATHQ_INVALID
	ag.TargetReplan = false
	if ag.TargetRef != 0 {
		ag.TargetState = DT_CROWDAGENT_TARGET_REQUESTING
	} else {
		ag.TargetState = DT_CROWDAGENT_TARGET_FAILED
	}

	return true
}

func (d *DtCrowd) RequestMoveVelocity(idx int, vel []float64) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.TargetRef = 0
	copy(ag.TargetPos[:], vel)
	ag.TargetPathqRef = DT_PATHQ_INVALID
	ag.TargetReplan = false
	ag.TargetState = DT_CROWDAGENT_TARGET_VELOCITY

	return true
}

func (d *DtCrowd) resetMoveTarget(idx int) bool {
	if idx < 0 || idx >= d.m_maxAgents {
		return false
	}

	ag := d.m_agents[idx]

	// Initialize request.
	ag.TargetRef = 0
	common.Vset(ag.TargetPos[:], 0, 0, 0)
	common.Vset(ag.Dvel[:], 0, 0, 0)
	ag.TargetPathqRef = DT_PATHQ_INVALID
	ag.TargetReplan = false
	ag.TargetState = DT_CROWDAGENT_TARGET_NONE

	return true
}

func (d *DtCrowd) getActiveAgents(agents []*DtCrowdAgent, maxAgents int) int {
	n := 0
	for i := 0; i < d.m_maxAgents; i++ {
		if !d.m_agents[i].Active {
			continue
		}
		if n < maxAgents {
			agents[n] = d.m_agents[i]
			n++
		}
	}
	return n
}

func (d *DtCrowd) updateMoveRequest(dt float64) {
	PATH_MAX_AGENTS := 8
	queue := make([]*DtCrowdAgent, PATH_MAX_AGENTS)
	nqueue := 0

	// Fire off new requests.
	for i := 0; i < d.m_maxAgents; i++ {
		ag := d.m_agents[i]
		if !ag.Active {
			continue
		}

		if ag.State == DT_CROWDAGENT_STATE_INVALID {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_REQUESTING {
			path := ag.Corridor.getPath()
			npath := ag.Corridor.getPathCount()

			MAX_RES := 32
			reqPos := make([]float64, 3)
			reqPath := make([]DtPolyRef, MAX_RES) // The path to the request location
			reqPathCount := 0

			// Quick search towards the goal.
			MAX_ITER := 20
			d.m_navquery.InitSlicedFindPath(path[0], ag.TargetRef, ag.Npos[:], ag.TargetPos[:], d.m_filters[ag.Params.QueryFilterType], 0)
			var status DtStatus
			d.m_navquery.UpdateSlicedFindPath(MAX_ITER)

			if ag.TargetReplan { // && npath > 10)
				// Try to use existing steady path during replan if possible.
				reqPathCount, status = d.m_navquery.FinalizeSlicedFindPathPartial(path, npath, reqPath, MAX_RES)
			} else {
				// Try to move towards target when goal changes.
				reqPathCount, status = d.m_navquery.FinalizeSlicedFindPath(reqPath, MAX_RES)
			}

			if !status.DtStatusFailed() && reqPathCount > 0 {
				// In progress or succeed.
				if reqPath[reqPathCount-1] != ag.TargetRef {
					// Partial path, constrain target position inside the last polygon.
					var tmp bool
					status = d.m_navquery.ClosestPointOnPoly(reqPath[reqPathCount-1], ag.TargetPos[:], reqPos, &tmp)
					if status.DtStatusFailed() {
						reqPathCount = 0
					}

				} else {
					copy(reqPos, ag.TargetPos[:])
				}
			} else {
				reqPathCount = 0
			}

			if reqPathCount == 0 {
				// Could not find path, start the request from current location.
				copy(reqPos, ag.Npos[:])
				reqPath[0] = path[0]
				reqPathCount = 1
			}

			ag.Corridor.setCorridor(reqPos, reqPath, reqPathCount)
			ag.Boundary.reset()
			ag.Partial = false

			if reqPath[reqPathCount-1] == ag.TargetRef {
				ag.TargetState = DT_CROWDAGENT_TARGET_VALID
				ag.TargetReplanTime = 0.0
			} else {
				// The path is longer or potentially unreachable, full plan.
				ag.TargetState = DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE
			}
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE {
			nqueue = addToPathQueue(ag, queue, nqueue, PATH_MAX_AGENTS)
		}
	}

	for i := 0; i < nqueue; i++ {
		ag := queue[i]
		tmp := ag.Corridor.getTarget()
		ag.TargetPathqRef = d.m_pathq.request(ag.Corridor.getLastPoly(), ag.TargetRef,
			tmp[:], ag.TargetPos[:], d.m_filters[ag.Params.QueryFilterType])
		if ag.TargetPathqRef != DT_PATHQ_INVALID {
			ag.TargetState = DT_CROWDAGENT_TARGET_WAITING_FOR_PATH
		}

	}

	// Update requests.
	d.m_pathq.update(MAX_ITERS_PER_UPDATE)

	var status DtStatus

	// Process path results.
	for i := 0; i < d.m_maxAgents; i++ {
		ag := d.m_agents[i]
		if !ag.Active {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_WAITING_FOR_PATH {
			// Poll path queue.
			status = d.m_pathq.getRequestStatus(ag.TargetPathqRef)
			if status.DtStatusFailed() {
				// Path find failed, retry if the target location is still valid.
				ag.TargetPathqRef = DT_PATHQ_INVALID
				if ag.TargetRef != 0 {
					ag.TargetState = DT_CROWDAGENT_TARGET_REQUESTING
				} else {
					ag.TargetState = DT_CROWDAGENT_TARGET_FAILED
				}

				ag.TargetReplanTime = 0.0
			} else if status.DtStatusSucceed() {
				path := ag.Corridor.getPath()
				npath := ag.Corridor.getPathCount()
				// Apply results.
				targetPos := make([]float64, 3)
				copy(targetPos, ag.TargetPos[:])

				res := d.m_pathResult
				valid := true
				nres := 0
				status = d.m_pathq.getPathResult(ag.TargetPathqRef, res, &nres, d.m_maxPathResult)
				if status.DtStatusFailed() || nres == 0 {
					valid = false
				}

				if status.DtStatusDetail(DT_PARTIAL_RESULT) {
					ag.Partial = true
				} else {
					ag.Partial = false
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
					if res[nres-1] != ag.TargetRef {
						// Partial path, constrain target position inside the last polygon.
						var status DtStatus
						nearest := make([]float64, 3)
						var tmp bool
						status = d.m_navquery.ClosestPointOnPoly(res[nres-1], targetPos, nearest, &tmp)
						if status.DtStatusSucceed() {
							copy(targetPos, nearest)
						} else {
							valid = false
						}

					}
				}

				if valid {
					// Set current corridor.
					ag.Corridor.setCorridor(targetPos, res, nres)
					// Force to update boundary.
					ag.Boundary.reset()
					ag.TargetState = DT_CROWDAGENT_TARGET_VALID
				} else {
					// Something went wrong.
					ag.TargetState = DT_CROWDAGENT_TARGET_FAILED
				}

				ag.TargetReplanTime = 0.0
			}
		}
	}

}

func (d *DtCrowd) updateTopologyOptimization(agents []*DtCrowdAgent, nagents int, dt float64) {
	if nagents == 0 {
		return
	}

	OPT_TIME_THR := 0.5 // seconds
	const OPT_MAX_AGENTS = 1
	queue := [OPT_MAX_AGENTS]*DtCrowdAgent{}
	nqueue := 0

	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		if (ag.Params.UpdateFlags & DT_CROWD_OPTIMIZE_TOPO) == 0 {
			continue
		}

		ag.TopologyOptTime += dt
		if ag.TopologyOptTime >= OPT_TIME_THR {
			nqueue = addToOptQueue(ag, queue[:], nqueue, OPT_MAX_AGENTS)
		}

	}

	for i := 0; i < nqueue; i++ {
		ag := queue[i]
		ag.Corridor.optimizePathTopology(d.m_navquery, d.m_filters[ag.Params.QueryFilterType])
		ag.TopologyOptTime = 0
	}

}

func (d *DtCrowd) checkPathValidity(agents []*DtCrowdAgent, nagents int, dt float64) {
	const CHECK_LOOKAHEAD = 10
	const TARGET_REPLAN_DELAY = 1.0 // seconds

	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		ag.TargetReplanTime += dt

		replan := false

		// First check that the current location is valid.
		idx := d.GetAgentIndex(ag)
		agentPos := make([]float64, 3)
		agentRef := ag.Corridor.getFirstPoly()
		copy(agentPos, ag.Npos[:])
		if !d.m_navquery.IsValidPolyRef(agentRef, d.m_filters[ag.Params.QueryFilterType]) {
			// Current location is not valid, try to reposition.
			// TODO: this can snap agents, how to handle that?
			nearest := make([]float64, 3)
			copy(nearest, agentPos)
			agentRef = 0
			agentRef, _ = d.m_navquery.FindNearestPoly(ag.Npos[:], d.m_agentPlacementHalfExtents[:], d.m_filters[ag.Params.QueryFilterType], nearest)
			copy(agentPos, nearest)

			if agentRef == 0 {
				// Could not find location in navmesh, set state to invalid.
				ag.Corridor.reset(0, agentPos)
				ag.Partial = false
				ag.Boundary.reset()
				ag.State = DT_CROWDAGENT_STATE_INVALID
				continue
			}

			// Make sure the first polygon is valid, but leave other valid
			// polygons in the path so that replanner can adjust the path better.
			ag.Corridor.fixPathStart(agentRef, agentPos)
			//			ag->corridor.trimInvalidPath(agentRef, agentPos, m_navquery, &m_filter);
			ag.Boundary.reset()
			copy(ag.Npos[:], agentPos)

			replan = true
		}

		// If the agent does not have move target or is controlled by velocity, no need to recover the target nor replan.
		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		// Try to recover move request position.
		if ag.TargetState != DT_CROWDAGENT_TARGET_NONE && ag.TargetState != DT_CROWDAGENT_TARGET_FAILED {
			if !d.m_navquery.IsValidPolyRef(ag.TargetRef, d.m_filters[ag.Params.QueryFilterType]) {
				// Current target is not valid, try to reposition.
				nearest := make([]float64, 3)
				copy(nearest, ag.TargetPos[:])
				ag.TargetRef = 0
				ag.TargetRef, _ = d.m_navquery.FindNearestPoly(ag.TargetPos[:], d.m_agentPlacementHalfExtents[:], d.m_filters[ag.Params.QueryFilterType], nearest)
				copy(ag.TargetPos[:], nearest)
				replan = true
			}
			if ag.TargetRef == 0 {
				// Failed to reposition target, fail moverequest.
				ag.Corridor.reset(agentRef, agentPos)
				ag.Partial = false
				ag.TargetState = DT_CROWDAGENT_TARGET_NONE
			}
		}

		// If nearby corridor is not valid, replan.
		if !ag.Corridor.isValid(CHECK_LOOKAHEAD, d.m_navquery, d.m_filters[ag.Params.QueryFilterType]) {
			// Fix current path.
			//			ag->corridor.trimInvalidPath(agentRef, agentPos, m_navquery, &m_filter);
			//			ag->boundary.reset();
			replan = true
		}

		// If the end of the path is near and it is not the requested location, replan.
		if ag.TargetState == DT_CROWDAGENT_TARGET_VALID {
			if ag.TargetReplanTime > TARGET_REPLAN_DELAY && ag.Corridor.getPathCount() < CHECK_LOOKAHEAD && ag.Corridor.getLastPoly() != ag.TargetRef {
				replan = true
			}

		}

		// Try to replan path to goal.
		if replan {
			if ag.TargetState != DT_CROWDAGENT_TARGET_NONE {
				d.requestMoveTargetReplan(idx, ag.TargetRef, ag.TargetPos[:])
			}
		}
	}
}

func (d *DtCrowd) Update(dt float64, debug *DtCrowdAgentDebugInfo) {
	d.m_velocitySampleCount = 0

	debugIdx := -1
	if debug != nil {
		debugIdx = debug.Idx
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
		p := ag.Npos
		r := ag.Params.Radius
		d.m_grid.addItem(i, p[0]-r, p[2]-r, p[0]+r, p[2]+r)
	}

	// Get nearby navmesh segments and agents to collide with.
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		// Update the collision boundary after certain distance has been passed or
		// if it has become invalid.
		updateThr := ag.Params.CollisionQueryRange * 0.25
		tmp := ag.Boundary.getCenter()
		if common.Vdist2DSqr(ag.Npos[:], tmp[:]) > common.Sqr(updateThr) ||
			!ag.Boundary.isValid(d.m_navquery, d.m_filters[ag.Params.QueryFilterType]) {
			ag.Boundary.update(ag.Corridor.getFirstPoly(), ag.Npos[:], ag.Params.CollisionQueryRange,
				d.m_navquery, d.m_filters[ag.Params.QueryFilterType])
		}
		// Query neighbour agents
		ag.Nneis = getNeighbours(ag.Npos[:], ag.Params.Height, ag.Params.CollisionQueryRange,
			ag, ag.Neis[:], DT_CROWDAGENT_MAX_NEIGHBOURS,
			agents, nagents, d.m_grid)
		for j := 0; j < ag.Nneis; j++ {
			ag.Neis[j].Idx = d.GetAgentIndex(agents[ag.Neis[j].Idx])
		}

	}

	// Find next corner to steer to.
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		// Find corners for steering
		ag.Ncorners = ag.Corridor.findCorners(ag.CornerVerts[:], ag.CornerFlags[:], ag.CornerPolys[:],
			DT_CROWDAGENT_MAX_CORNERS, d.m_navquery, d.m_filters[ag.Params.QueryFilterType])

		// Check to see if the corner after the next corner is directly visible,
		// and short cut to there.
		if (ag.Params.UpdateFlags&DT_CROWD_OPTIMIZE_VIS > 0) && ag.Ncorners > 0 {
			target := rcGetVert(ag.CornerVerts[:], common.Min(1, ag.Ncorners-1))
			ag.Corridor.optimizePathVisibility(target, ag.Params.PathOptimizationRange, d.m_navquery, d.m_filters[ag.Params.QueryFilterType])

			// Copy data for debug purposes.
			if debugIdx == i {
				tmp := ag.Corridor.getPos()
				copy(debug.OptStart[:], tmp[:])
				copy(debug.OptEnd[:], target)
			}
		} else {
			// Copy data for debug purposes.
			if debugIdx == i {
				common.Vset(debug.OptStart[:], 0, 0, 0)
				common.Vset(debug.OptEnd[:], 0, 0, 0)
			}
		}
	}

	// Trigger off-mesh connections (depends on corners).
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			continue
		}

		// Check
		triggerRadius := ag.Params.Radius * 2.25
		if overOffmeshConnection(ag, triggerRadius) {
			// Prepare to off-mesh connection.
			idx := d.GetAgentIndex(ag)
			anim := d.m_agentAnims[idx]

			// Adjust the path over the off-mesh connection.
			refs := make([]DtPolyRef, 2)
			if ag.Corridor.moveOverOffmeshConnection(ag.CornerPolys[ag.Ncorners-1], refs,
				anim.StartPos[:], anim.EndPos[:], d.m_navquery) {
				copy(anim.InitPos[:], ag.Npos[:])
				anim.PolyRef = refs[1]
				anim.Active = true
				anim.T = 0.0
				anim.Tmax = (common.Vdist2D(anim.StartPos[:], anim.EndPos[:]) / ag.Params.MaxSpeed) * 0.5

				ag.State = DT_CROWDAGENT_STATE_OFFMESH
				ag.Ncorners = 0
				ag.Nneis = 0
				continue
			} else {
				// Path validity check will ensure that bad/blocked connections will be replanned.
			}
		}
	}

	// Calculate steering.
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE {
			continue
		}

		dvel := []float64{0, 0, 0}

		if ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			copy(dvel, ag.TargetPos[:])
			ag.desiredSpeed = common.Vlen(ag.TargetPos[:])
		} else {
			// Calculate steering direction.
			if ag.Params.UpdateFlags&DT_CROWD_ANTICIPATE_TURNS > 0 {
				calcSmoothSteerDirection(ag, dvel)
			} else {
				calcStraightSteerDirection(ag, dvel)
			}

			// Calculate speed scale, which tells the agent to slowdown at the end of the path.
			slowDownRadius := ag.Params.Radius * 2 // TODO: make less hacky.
			speedScale := getDistanceToGoal(ag, slowDownRadius) / slowDownRadius

			ag.desiredSpeed = ag.Params.MaxSpeed
			copy(dvel, common.Vscale(dvel, ag.desiredSpeed*speedScale))

		}

		// Separation
		if ag.Params.UpdateFlags&DT_CROWD_SEPARATION > 0 {
			separationDist := ag.Params.CollisionQueryRange
			invSeparationDist := 1.0 / separationDist
			separationWeight := ag.Params.SeparationWeight

			w := float64(0)
			disp := []float64{0, 0, 0}

			for j := 0; j < ag.Nneis; j++ {
				nei := d.m_agents[ag.Neis[j].Idx]

				diff := common.Vsub(ag.Npos[:], nei.Npos[:])
				diff[1] = 0

				distSqr := common.VlenSqr(diff)
				if distSqr < 0.00001 {
					continue
				}

				if distSqr > common.Sqr(separationDist) {
					continue
				}

				dist := math.Sqrt(distSqr)
				weight := separationWeight * (1.0 - common.Sqr(dist*invSeparationDist))

				common.Vmad(disp, disp, diff, weight/dist)
				w += 1.0
			}

			if w > 0.0001 {
				// Adjust desired velocity.
				common.Vmad(dvel, dvel, disp, 1.0/w)
				// Clamp desired velocity to desired speed.
				speedSqr := common.VlenSqr(dvel)
				desiredSqr := common.Sqr(ag.desiredSpeed)
				if speedSqr > desiredSqr {
					copy(dvel, common.Vscale(dvel, desiredSqr/speedSqr))
				}

			}
		}

		// Set the desired velocity.
		copy(ag.Dvel[:], dvel[:])
	}

	// Velocity planning.
	for i := 0; i < nagents; i++ {
		ag := agents[i]

		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		if ag.Params.UpdateFlags&DT_CROWD_OBSTACLE_AVOIDANCE > 0 {
			d.m_obstacleQuery.reset()

			// Add neighbours as obstacles.
			for j := 0; j < ag.Nneis; j++ {
				nei := d.m_agents[ag.Neis[j].Idx]
				d.m_obstacleQuery.addCircle(nei.Npos[:], nei.Params.Radius, nei.Vel[:], nei.Dvel[:])
			}

			// Append neighbour segments as obstacles.
			for j := 0; j < ag.Boundary.getSegmentCount(); j++ {
				s := ag.Boundary.getSegment(j)
				if common.TriArea2D(ag.Npos[:], s[:], s[3:]) < 0.0 {
					continue
				}

				d.m_obstacleQuery.addSegment(s[:], s[3:])
			}

			var vod *DtObstacleAvoidanceDebugData
			if debugIdx == i {
				vod = debug.Vod
			}

			// Sample new safe velocity.
			adaptive := true
			ns := 0

			params := d.m_obstacleQueryParams[ag.Params.ObstacleAvoidanceType]

			if adaptive {
				ns = d.m_obstacleQuery.sampleVelocityAdaptive(ag.Npos[:], ag.Params.Radius, ag.desiredSpeed,
					ag.Vel[:], ag.Dvel[:], ag.Nvel[:], params, vod)
			} else {
				ns = d.m_obstacleQuery.sampleVelocityGrid(ag.Npos[:], ag.Params.Radius, ag.desiredSpeed,
					ag.Vel[:], ag.Dvel[:], ag.Nvel[:], params, vod)
			}
			d.m_velocitySampleCount += ns
		} else {
			// If not using velocity planning, new velocity is directly the desired velocity.
			copy(ag.Nvel[:], ag.Dvel[:])
		}
	}

	// Integrate.
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		integrate(ag, dt)
	}

	// Handle collisions.
	const COLLISION_RESOLVE_FACTOR = 0.7

	for iter := 0; iter < 4; iter++ {
		for i := 0; i < nagents; i++ {
			ag := agents[i]
			idx0 := d.GetAgentIndex(ag)

			if ag.State != DT_CROWDAGENT_STATE_WALKING {
				continue
			}

			common.Vset(ag.Disp[:], 0, 0, 0)

			w := float64(0)

			for j := 0; j < ag.Nneis; j++ {
				nei := d.m_agents[ag.Neis[j].Idx]
				idx1 := d.GetAgentIndex(nei)

				diff := common.Vsub(ag.Npos[:], nei.Npos[:])
				diff[1] = 0

				dist := common.VlenSqr(diff)
				if dist > common.Sqr(ag.Params.Radius+nei.Params.Radius) {
					continue
				}

				dist = math.Sqrt(dist)
				pen := (ag.Params.Radius + nei.Params.Radius) - dist
				if dist < 0.0001 {
					// Agents on top of each other, try to choose diverging separation directions.
					if idx0 > idx1 {
						common.Vset(diff, -ag.Dvel[2], 0, ag.Dvel[0])
					} else {
						common.Vset(diff, ag.Dvel[2], 0, -ag.Dvel[0])
					}

					pen = 0.01
				} else {
					pen = (1.0 / dist) * (pen * 0.5) * COLLISION_RESOLVE_FACTOR
				}

				common.Vmad(ag.Disp[:], ag.Disp[:], diff, pen)

				w += 1.0
			}

			if w > 0.0001 {
				iw := 1.0 / w
				copy(ag.Disp[:], common.Vscale(ag.Disp[:], iw))
			}
		}

		for i := 0; i < nagents; i++ {
			ag := agents[i]
			if ag.State != DT_CROWDAGENT_STATE_WALKING {
				continue
			}

			copy(ag.Npos[:], common.Vadd(ag.Npos[:], ag.Disp[:]))

		}
	}

	for i := 0; i < nagents; i++ {
		ag := agents[i]
		if ag.State != DT_CROWDAGENT_STATE_WALKING {
			continue
		}

		// Move along navmesh.
		ag.Corridor.movePosition(ag.Npos[:], d.m_navquery, d.m_filters[ag.Params.QueryFilterType])
		// Get valid constrained position back.
		tmp := ag.Corridor.getPos()
		copy(ag.Npos[:], tmp[:])

		// If not using path, truncate the corridor to just one poly.
		if ag.TargetState == DT_CROWDAGENT_TARGET_NONE || ag.TargetState == DT_CROWDAGENT_TARGET_VELOCITY {
			ag.Corridor.reset(ag.Corridor.getFirstPoly(), ag.Npos[:])
			ag.Partial = false
		}

	}

	// Update agents using off-mesh connection.
	for i := 0; i < nagents; i++ {
		ag := agents[i]
		idx := d.GetAgentIndex(ag)
		anim := d.m_agentAnims[idx]
		if !anim.Active {
			continue
		}

		anim.T += dt
		if anim.T > anim.Tmax {
			// Reset animation
			anim.Active = false
			// Prepare agent for walking.
			ag.State = DT_CROWDAGENT_STATE_WALKING
			continue
		}

		// Update position
		ta := anim.Tmax * 0.15
		tb := anim.Tmax
		if anim.T < ta {
			u := tween(anim.T, 0.0, ta)
			common.Vlerp(ag.Npos[:], anim.InitPos[:], anim.StartPos[:], u)

		} else {
			u := tween(anim.T, ta, tb)
			common.Vlerp(ag.Npos[:], anim.StartPos[:], anim.EndPos[:], u)

		}

		// Update velocity.
		common.Vset(ag.Vel[:], 0, 0, 0)
		common.Vset(ag.Dvel[:], 0, 0, 0)
	}

}
