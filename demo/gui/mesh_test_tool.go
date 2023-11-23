package gui

import "gonavamesh/recast"

type ToolMode int

const (
	TOOLMODE_PATHFIND_FOLLOW ToolMode = iota
	TOOLMODE_PATHFIND_STRAIGHT
	TOOLMODE_PATHFIND_SLICED
	TOOLMODE_RAYCAST
	TOOLMODE_DISTANCE_TO_WALL
	TOOLMODE_FIND_POLYS_IN_CIRCLE
	TOOLMODE_FIND_POLYS_IN_SHAPE
	TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD
)

const (
	MAX_POLYS        = 256
	MAX_SMOOTH       = 2048
	MAX_STEER_POINTS = 10
	MAX_RAND_POINTS  = 64
)

var (
	menus = map[int]func(){}
)

type MeshTool struct {
	m_sample         *Sample
	m_navMesh        recast.IDtNavMesh
	m_navQuery       recast.NavMeshQuery
	m_filter         *recast.DtQueryFilter
	m_pathFindStatus recast.DtStatus

	m_toolMode ToolMode

	m_startRef recast.DtPolyRef
	m_endRef   recast.DtPolyRef
	m_polys    []recast.DtPolyRef
	m_parent   []recast.DtPolyRef

	m_straightPathOptions int

	m_npolys            int
	m_straightPath      []float64
	m_straightPathFlags []int
	m_straightPathPolys []recast.DtPolyRef
	m_nstraightPath     int
	m_polyPickExt       []float64
	m_smoothPath        []float64
	m_nsmoothPath       int
	m_queryPoly         []float64

	m_randPoints         []float64
	m_nrandPoints        int
	m_randPointsInCircle bool

	m_spos                []float64
	m_epos                []float64
	m_hitPos              []float64
	m_hitNormal           []float64
	m_hitResult           bool
	m_distanceToWall      float64
	m_neighbourhoodRadius float64
	m_randomRadius        float64
	m_sposSet             bool
	m_eposSet             bool

	m_pathIterNum                                     int
	m_pathIterPolys                                   []recast.DtPolyRef
	m_pathIterPolyCount                               int
	m_prevIterPos, m_iterPos, m_steerPos, m_targetPos []float64
	m_steerPoints                                     []float64
	m_steerPointCount                                 int

	gs *guiState
}

func newMeshTool(gs *guiState) *MeshTool {
	m := &MeshTool{
		m_straightPath:      make([]float64, MAX_POLYS*3),
		m_straightPathFlags: make([]int, MAX_POLYS),
		m_straightPathPolys: make([]recast.DtPolyRef, MAX_POLYS),
		m_smoothPath:        make([]float64, MAX_SMOOTH*3),
		m_randPoints:        make([]float64, MAX_RAND_POINTS*3),
		m_queryPoly:         make([]float64, 4*3),
		m_polys:             make([]recast.DtPolyRef, MAX_POLYS),
		m_parent:            make([]recast.DtPolyRef, MAX_POLYS),
		m_spos:              make([]float64, 3),
		m_epos:              make([]float64, 3),
		m_hitPos:            make([]float64, 3),
		m_hitNormal:         make([]float64, 3),
		m_pathIterPolys:     make([]recast.DtPolyRef, MAX_POLYS),
		m_steerPoints:       make([]float64, MAX_STEER_POINTS*3),
		m_prevIterPos:       make([]float64, 3),
		m_iterPos:           make([]float64, 3),
		m_steerPos:          make([]float64, 3),
		m_targetPos:         make([]float64, 3),
		m_polyPickExt:       make([]float64, 3),
		m_toolMode:          TOOLMODE_PATHFIND_FOLLOW,
		m_pathFindStatus:    recast.DT_FAILURE,
		gs:                  gs,
	}
	m.m_filter.SetIncludeFlags(SAMPLE_POLYFLAGS_ALL ^ SAMPLE_POLYFLAGS_DISABLED)
	m.m_filter.SetExcludeFlags(0)

	m.m_polyPickExt[0] = 2
	m.m_polyPickExt[1] = 4
	m.m_polyPickExt[2] = 2

	m.m_neighbourhoodRadius = 2.5
	m.m_randomRadius = 5.0
	return m
}
func (m *MeshTool) init(sample *Sample) {
	m.m_sample = sample
	m.m_navMesh = sample.getNavMesh()
	m.m_navQuery = sample.getNavMeshQuery()
	m.recalc()
	if m.m_navQuery != nil {
		// Change costs.
		m.m_filter.SetAreaCost(SAMPLE_POLYAREA_GROUND, 1.0)
		m.m_filter.SetAreaCost(SAMPLE_POLYAREA_WATER, 10.0)
		m.m_filter.SetAreaCost(SAMPLE_POLYAREA_ROAD, 1.0)
		m.m_filter.SetAreaCost(SAMPLE_POLYAREA_DOOR, 1.0)
		m.m_filter.SetAreaCost(SAMPLE_POLYAREA_GRASS, 2.0)
		m.m_filter.SetAreaCost(SAMPLE_POLYAREA_JUMP, 1.5)
	}
	m.m_neighbourhoodRadius = sample.getAgentRadius() * 20.0
	m.m_randomRadius = sample.getAgentRadius() * 30.0
}

func (m *MeshTool) handleMenu() {
	if m.gs.imguiCheck("Pathfind Follow", m.m_toolMode == TOOLMODE_PATHFIND_FOLLOW) {
		m.m_toolMode = TOOLMODE_PATHFIND_FOLLOW
		m.recalc()
	}
	if m.gs.imguiCheck("Pathfind Straight", m.m_toolMode == TOOLMODE_PATHFIND_STRAIGHT) {
		m.m_toolMode = TOOLMODE_PATHFIND_STRAIGHT
		m.recalc()
	}
	if m.m_toolMode == TOOLMODE_PATHFIND_STRAIGHT {
		m.gs.imguiIndent()
		m.gs.imguiLabel("Vertices at crossings")
		if m.gs.imguiCheck("None", m.m_straightPathOptions == 0) {
			m.m_straightPathOptions = 0
			m.recalc()
		}
		if m.gs.imguiCheck("Area", m.m_straightPathOptions == recast.DT_STRAIGHTPATH_AREA_CROSSINGS) {
			m.m_straightPathOptions = recast.DT_STRAIGHTPATH_AREA_CROSSINGS
			m.recalc()
		}
		if m.gs.imguiCheck("All", m.m_straightPathOptions == recast.DT_STRAIGHTPATH_ALL_CROSSINGS) {
			m.m_straightPathOptions = recast.DT_STRAIGHTPATH_ALL_CROSSINGS
			m.recalc()
		}

		m.gs.imguiUnindent()
	}
	if m.gs.imguiCheck("Pathfind Sliced", m.m_toolMode == TOOLMODE_PATHFIND_SLICED) {
		m.m_toolMode = TOOLMODE_PATHFIND_SLICED
		m.recalc()
	}

	m.gs.imguiSeparator()

	if m.gs.imguiCheck("Distance to Wall", m.m_toolMode == TOOLMODE_DISTANCE_TO_WALL) {
		m.m_toolMode = TOOLMODE_DISTANCE_TO_WALL
		m.recalc()
	}

	m.gs.imguiSeparator()

	if m.gs.imguiCheck("Raycast", m.m_toolMode == TOOLMODE_RAYCAST) {
		m.m_toolMode = TOOLMODE_RAYCAST
		m.recalc()
	}

	m.gs.imguiSeparator()
	if m.gs.imguiCheck("Find Polys in Circle", m.m_toolMode == TOOLMODE_FIND_POLYS_IN_CIRCLE) {
		m.m_toolMode = TOOLMODE_FIND_POLYS_IN_CIRCLE
		m.recalc()
	}
	if m.gs.imguiCheck("Find Polys in Shape", m.m_toolMode == TOOLMODE_FIND_POLYS_IN_SHAPE) {
		m.m_toolMode = TOOLMODE_FIND_POLYS_IN_SHAPE
		m.recalc()
	}

	m.gs.imguiSeparator()

	if m.gs.imguiCheck("Find Local Neighbourhood", m.m_toolMode == TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD) {
		m.m_toolMode = TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD
		m.recalc()
	}

	m.gs.imguiSeparator()

	if m.gs.imguiButton("Set Random Start") {
		var status recast.DtStatus
		m.m_startRef, m.m_spos, status = m.m_navQuery.FindRandomPoint(m.m_filter)
		if status.DtStatusSucceed() {
			m.m_sposSet = true
			m.recalc()
		}
	}
	if m.gs.imguiButton("Set Random End", m.m_sposSet) {
		if m.m_sposSet {
			var status recast.DtStatus
			m.m_startRef, m.m_spos, status = m.m_navQuery.FindRandomPointAroundCircle(m.m_startRef, m.m_spos, m.m_randomRadius, m.m_filter)
			if status.DtStatusSucceed() {
				m.m_eposSet = true
				m.recalc()
			}
		}
	}

	m.gs.imguiSeparator()

	if m.gs.imguiButton("Make Random Points") {
		m.m_randPointsInCircle = false
		m.m_nrandPoints = 0
		for i := 0; i < MAX_RAND_POINTS; i++ {
			_, pt, status := m.m_navQuery.FindRandomPoint(m.m_filter)
			if status.DtStatusSucceed() {
				copy(m.m_randPoints[m.m_nrandPoints*3:], pt)
				m.m_nrandPoints++
			}
		}
	}
	if m.gs.imguiButton("Make Random Points Around", m.m_sposSet) {
		if m.m_sposSet {
			m.m_nrandPoints = 0
			m.m_randPointsInCircle = true
			for i := 0; i < MAX_RAND_POINTS; i++ {
				_, pt, status := m.m_navQuery.FindRandomPointAroundCircle(m.m_startRef, m.m_spos, m.m_randomRadius, m.m_filter)
				if status.DtStatusSucceed() {
					copy(m.m_randPoints[m.m_nrandPoints*3:], pt)
					m.m_nrandPoints++
				}
			}
		}
	}

	m.gs.imguiSeparator()

	m.gs.imguiLabel("Include Flags")

	m.gs.imguiIndent()
	if m.gs.imguiCheck("Walk", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_WALK) != 0) {
		m.m_filter.SetIncludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_WALK)
		m.recalc()
	}
	if m.gs.imguiCheck("Swim", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_SWIM) != 0) {
		m.m_filter.SetIncludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_SWIM)
		m.recalc()
	}
	if m.gs.imguiCheck("Door", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_DOOR) != 0) {
		m.m_filter.SetIncludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_DOOR)
		m.recalc()
	}
	if m.gs.imguiCheck("Jump", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_JUMP) != 0) {
		m.m_filter.SetIncludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_JUMP)
		m.recalc()
	}
	m.gs.imguiUnindent()

	m.gs.imguiSeparator()
	m.gs.imguiLabel("Exclude Flags")

	m.gs.imguiIndent()
	if m.gs.imguiCheck("Walk", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_WALK) != 0) {
		m.m_filter.SetExcludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_WALK)
		m.recalc()
	}
	if m.gs.imguiCheck("Swim", (m.m_filter.GetExcludeFlags()&SAMPLE_POLYFLAGS_SWIM) != 0) {
		m.m_filter.SetExcludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_SWIM)
		m.recalc()
	}
	if m.gs.imguiCheck("Door", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_DOOR) != 0) {
		m.m_filter.SetExcludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_DOOR)
		m.recalc()
	}
	if m.gs.imguiCheck("Jump", (m.m_filter.GetIncludeFlags()&SAMPLE_POLYFLAGS_JUMP) != 0) {
		m.m_filter.SetExcludeFlags(m.m_filter.GetIncludeFlags() ^ SAMPLE_POLYFLAGS_JUMP)
		m.recalc()
	}
	m.gs.imguiUnindent()

	m.gs.imguiSeparator()

}
func (m *MeshTool) handleClick(s, p []float64, shift bool) {
	if shift {
		m.m_sposSet = true
		copy(m.m_spos, p)
	} else {
		m.m_eposSet = true
		copy(m.m_epos, p)
	}
	m.recalc()
}

func (m *MeshTool) handleUpdate(dt float64) {
	if (m.m_toolMode == TOOLMODE_PATHFIND_SLICED) {
		if (m.m_pathFindStatus.DtStatusInProgress()) {
			_, m.m_pathFindStatus = m.m_navQuery.UpdateSlicedFindPath(1);
		}
		if (m.m_pathFindStatus.DtStatusSucceed()) {
			m.m_npolys, _ = m.m_navQuery.FinalizeSlicedFindPath(m.m_polys, MAX_POLYS);
			m.m_nstraightPath = 0;
			if (m.m_npolys != 0) {
				// In case of partial path, make sure the end point is clamped to the last polygon.
				epos := make([]float64, 3)
				copy(epos, m.m_epos);
				if (m.m_polys[m.m_npolys-1] != m.m_endRef) {
					epos, _, _ = m.m_navQuery.ClosestPointOnPoly(m.m_polys[m.m_npolys-1], m.m_epos) //TODO
				}

				m.m_nstraightPath, _ = m.m_navQuery.FindStraightPath(m.m_spos, epos, m.m_polys, m.m_npolys,
					m.m_straightPath, m.m_straightPathFlags,
					m.m_straightPathPolys, MAX_POLYS, recast.DT_STRAIGHTPATH_ALL_CROSSINGS);
			}

			m.m_pathFindStatus = recast.DT_FAILURE;
		}
	}
}

func (m *MeshTool) recalc() {

}
func (m *MeshTool) reset() {
	m.m_startRef = 0
	m.m_endRef = 0
	m.m_npolys = 0
	m.m_nstraightPath = 0
	m.m_nsmoothPath = 0
	m.m_hitPos = make([]float64, len(m.m_hitPos))
	m.m_hitNormal = make([]float64, len(m.m_hitNormal))
	m.m_distanceToWall = 0
}
