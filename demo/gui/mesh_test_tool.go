package gui

import (
	"gonavamesh/common"
	"gonavamesh/debug_utils"
	"gonavamesh/recast"
	"math"
)

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

type NavMeshTesterTool struct {
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

func (m *NavMeshTesterTool) drawAgent(pos []float64, r, h, c float64, col int) {
	dd := m.m_sample.getDebugDraw()

	dd.DepthMask(false)

	// Agent dimensions.
	debug_utils.DuDebugDrawCylinderWire(dd, pos[0]-r, pos[1]+0.02, pos[2]-r, pos[0]+r, pos[1]+h, pos[2]+r, col, 2.0)

	debug_utils.DuDebugDrawCircle(dd, pos[0], pos[1]+c, pos[2], r, debug_utils.DuRGBA(0, 0, 0, 64), 1.0)

	colb := debug_utils.DuRGBA(0, 0, 0, 196)
	dd.Begin(debug_utils.DU_DRAW_LINES)
	dd.Vertex1(pos[0], pos[1]-c, pos[2], colb)
	dd.Vertex1(pos[0], pos[1]+c, pos[2], colb)
	dd.Vertex1(pos[0]-r/2, pos[1]+0.02, pos[2], colb)
	dd.Vertex1(pos[0]+r/2, pos[1]+0.02, pos[2], colb)
	dd.Vertex1(pos[0], pos[1]+0.02, pos[2]-r/2, colb)
	dd.Vertex1(pos[0], pos[1]+0.02, pos[2]+r/2, colb)
	dd.End()

	dd.DepthMask(true)
}

func (m *NavMeshTesterTool) handleRender() {

	dd := m.m_sample.getDebugDraw()

	startCol := debug_utils.DuRGBA(128, 25, 0, 192)
	endCol := debug_utils.DuRGBA(51, 102, 0, 129)
	pathCol := debug_utils.DuRGBA(0, 0, 0, 64)

	agentRadius := m.m_sample.getAgentRadius()
	agentHeight := m.m_sample.getAgentHeight()
	agentClimb := m.m_sample.getAgentClimb()

	dd.DepthMask(false)
	if m.m_sposSet {
		m.drawAgent(m.m_spos, agentRadius, agentHeight, agentClimb, startCol)
	}

	if m.m_eposSet {
		m.drawAgent(m.m_epos, agentRadius, agentHeight, agentClimb, endCol)
	}

	dd.DepthMask(true)

	if m.m_navMesh == nil {
		return
	}

	if m.m_toolMode == TOOLMODE_PATHFIND_FOLLOW {
		debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_startRef, startCol)
		debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_endRef, endCol)

		if m.m_npolys != 0 {
			for i := 0; i < m.m_npolys; i++ {
				if m.m_polys[i] == m.m_startRef || m.m_polys[i] == m.m_endRef {
					continue
				}

				debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_polys[i], pathCol)
			}
		}

		if m.m_nsmoothPath != 0 {
			dd.DepthMask(false)
			spathCol := debug_utils.DuRGBA(0, 0, 0, 220)
			dd.Begin(debug_utils.DU_DRAW_LINES, 3.0)
			for i := 0; i < m.m_nsmoothPath; i++ {
				dd.Vertex1(m.m_smoothPath[i*3], m.m_smoothPath[i*3+1]+0.1, m.m_smoothPath[i*3+2], spathCol)
			}

			dd.End()
			dd.DepthMask(true)
		}

		if m.m_pathIterNum != 0 {
			debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_pathIterPolys[0], debug_utils.DuRGBA(255, 255, 255, 128))

			dd.DepthMask(false)
			dd.Begin(debug_utils.DU_DRAW_LINES, 1.0)

			prevCol := debug_utils.DuRGBA(255, 192, 0, 220)
			curCol := debug_utils.DuRGBA(255, 255, 255, 220)
			steerCol := debug_utils.DuRGBA(0, 192, 255, 220)

			dd.Vertex1(m.m_prevIterPos[0], m.m_prevIterPos[1]-0.3, m.m_prevIterPos[2], prevCol)
			dd.Vertex1(m.m_prevIterPos[0], m.m_prevIterPos[1]+0.3, m.m_prevIterPos[2], prevCol)

			dd.Vertex1(m.m_iterPos[0], m.m_iterPos[1]-0.3, m.m_iterPos[2], curCol)
			dd.Vertex1(m.m_iterPos[0], m.m_iterPos[1]+0.3, m.m_iterPos[2], curCol)

			dd.Vertex1(m.m_prevIterPos[0], m.m_prevIterPos[1]+0.3, m.m_prevIterPos[2], prevCol)
			dd.Vertex1(m.m_iterPos[0], m.m_iterPos[1]+0.3, m.m_iterPos[2], prevCol)

			dd.Vertex1(m.m_prevIterPos[0], m.m_prevIterPos[1]+0.3, m.m_prevIterPos[2], steerCol)
			dd.Vertex1(m.m_steerPos[0], m.m_steerPos[1]+0.3, m.m_steerPos[2], steerCol)

			for i := 0; i < m.m_steerPointCount-1; i++ {
				dd.Vertex1(m.m_steerPoints[i*3+0], m.m_steerPoints[i*3+1]+0.2, m.m_steerPoints[i*3+2], debug_utils.DuDarkenCol(steerCol))
				dd.Vertex1(m.m_steerPoints[(i+1)*3+0], m.m_steerPoints[(i+1)*3+1]+0.2, m.m_steerPoints[(i+1)*3+2], debug_utils.DuDarkenCol(steerCol))
			}

			dd.End()
			dd.DepthMask(true)
		}
	} else if m.m_toolMode == TOOLMODE_PATHFIND_STRAIGHT ||
		m.m_toolMode == TOOLMODE_PATHFIND_SLICED {
		debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_startRef, startCol)
		debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_endRef, endCol)

		if m.m_npolys != 0 {
			for i := 0; i < m.m_npolys; i++ {
				if m.m_polys[i] == m.m_startRef || m.m_polys[i] == m.m_endRef {
					continue
				}

				debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_polys[i], pathCol)
			}
		}

		if m.m_nstraightPath != 0 {
			dd.DepthMask(false)
			spathCol := debug_utils.DuRGBA(64, 16, 0, 220)
			offMeshCol := debug_utils.DuRGBA(128, 96, 0, 220)
			dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
			for i := 0; i < m.m_nstraightPath-1; i++ {
				var col int
				if m.m_straightPathFlags[i]&recast.DT_STRAIGHTPATH_OFFMESH_CONNECTION != 0 {
					col = offMeshCol
				} else {
					col = spathCol
				}

				dd.Vertex1(m.m_straightPath[i*3], m.m_straightPath[i*3+1]+0.4, m.m_straightPath[i*3+2], col)
				dd.Vertex1(m.m_straightPath[(i+1)*3], m.m_straightPath[(i+1)*3+1]+0.4, m.m_straightPath[(i+1)*3+2], col)
			}
			dd.End()
			dd.Begin(debug_utils.DU_DRAW_POINTS, 6.0)
			for i := 0; i < m.m_nstraightPath; i++ {
				var col int
				if m.m_straightPathFlags[i]&recast.DT_STRAIGHTPATH_START != 0 {
					col = startCol
				} else if m.m_straightPathFlags[i]&recast.DT_STRAIGHTPATH_END != 0 {
					col = endCol
				} else if m.m_straightPathFlags[i]&recast.DT_STRAIGHTPATH_OFFMESH_CONNECTION != 0 {
					col = offMeshCol
				} else {
					col = spathCol
				}

				dd.Vertex1(m.m_straightPath[i*3], m.m_straightPath[i*3+1]+0.4, m.m_straightPath[i*3+2], col)
			}
			dd.End()
			dd.DepthMask(true)
		}
	} else if m.m_toolMode == TOOLMODE_RAYCAST {
		debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_startRef, startCol)

		if m.m_nstraightPath != 0 {
			for i := 1; i < m.m_npolys; i++ {
				debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_polys[i], pathCol)
			}

			dd.DepthMask(false)
			spathCol := debug_utils.DuRGBA(240, 240, 240, 220)
			if m.m_hitResult {
				spathCol = debug_utils.DuRGBA(64, 16, 0, 220)
			}
			dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
			for i := 0; i < m.m_nstraightPath-1; i++ {
				dd.Vertex1(m.m_straightPath[i*3], m.m_straightPath[i*3+1]+0.4, m.m_straightPath[i*3+2], spathCol)
				dd.Vertex1(m.m_straightPath[(i+1)*3], m.m_straightPath[(i+1)*3+1]+0.4, m.m_straightPath[(i+1)*3+2], spathCol)
			}
			dd.End()
			dd.Begin(debug_utils.DU_DRAW_POINTS, 4.0)
			for i := 0; i < m.m_nstraightPath; i++ {
				dd.Vertex1(m.m_straightPath[i*3], m.m_straightPath[i*3+1]+0.4, m.m_straightPath[i*3+2], spathCol)

			}
			dd.End()

			if m.m_hitResult {
				hitCol := debug_utils.DuRGBA(0, 0, 0, 128)
				dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
				dd.Vertex1(m.m_hitPos[0], m.m_hitPos[1]+0.4, m.m_hitPos[2], hitCol)
				dd.Vertex1(m.m_hitPos[0]+m.m_hitNormal[0]*agentRadius,
					m.m_hitPos[1]+0.4+m.m_hitNormal[1]*agentRadius,
					m.m_hitPos[2]+m.m_hitNormal[2]*agentRadius, hitCol)
				dd.End()
			}
			dd.DepthMask(true)
		}
	} else if m.m_toolMode == TOOLMODE_DISTANCE_TO_WALL {
		debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_startRef, startCol)
		dd.DepthMask(false)
		debug_utils.DuDebugDrawCircle(dd, m.m_spos[0], m.m_spos[1]+agentHeight/2, m.m_spos[2], m.m_distanceToWall, debug_utils.DuRGBA(64, 16, 0, 220), 2.0)
		dd.Begin(debug_utils.DU_DRAW_LINES, 3.0)
		dd.Vertex1(m.m_hitPos[0], m.m_hitPos[1]+0.02, m.m_hitPos[2], debug_utils.DuRGBA(0, 0, 0, 192))
		dd.Vertex1(m.m_hitPos[0], m.m_hitPos[1]+agentHeight, m.m_hitPos[2], debug_utils.DuRGBA(0, 0, 0, 192))
		dd.End()
		dd.DepthMask(true)
	} else if m.m_toolMode == TOOLMODE_FIND_POLYS_IN_CIRCLE {
		for i := 0; i < m.m_npolys; i++ {
			debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_polys[i], pathCol)
			dd.DepthMask(false)
			if m.m_parent[i] != 0 {
				p0 := make([]float64, 3)
				p1 := make([]float64, 3)
				dd.DepthMask(false)
				getPolyCenter(m.m_navMesh, m.m_parent[i], p0)
				getPolyCenter(m.m_navMesh, m.m_polys[i], p1)
				debug_utils.DuDebugDrawArc(dd, p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], 0.25, 0.0, 0.4, debug_utils.DuRGBA(0, 0, 0, 128), 2.0)
				dd.DepthMask(true)
			}
			dd.DepthMask(true)
		}

		if m.m_sposSet && m.m_eposSet {
			dd.DepthMask(false)
			dx := m.m_epos[0] - m.m_spos[0]
			dz := m.m_epos[2] - m.m_spos[2]
			dist := math.Sqrt(dx*dx + dz*dz)
			debug_utils.DuDebugDrawCircle(dd, m.m_spos[0], m.m_spos[1]+agentHeight/2, m.m_spos[2], dist, debug_utils.DuRGBA(64, 16, 0, 220), 2.0)
			dd.DepthMask(true)
		}
	} else if m.m_toolMode == TOOLMODE_FIND_POLYS_IN_SHAPE {
		for i := 0; i < m.m_npolys; i++ {
			debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_polys[i], pathCol)
			dd.DepthMask(false)
			if m.m_parent[i] != 0 {
				p0 := make([]float64, 3)
				p1 := make([]float64, 3)
				dd.DepthMask(false)
				getPolyCenter(m.m_navMesh, m.m_parent[i], p0)
				getPolyCenter(m.m_navMesh, m.m_polys[i], p1)
				debug_utils.DuDebugDrawArc(dd, p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], 0.25, 0.0, 0.4, debug_utils.DuRGBA(0, 0, 0, 128), 2.0)
				dd.DepthMask(true)
			}
			dd.DepthMask(true)
		}

		if m.m_sposSet && m.m_eposSet {
			dd.DepthMask(false)
			col := debug_utils.DuRGBA(64, 16, 0, 220)
			dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
			i := 0
			j := 3
			for i < 4 {
				p0 := m.m_queryPoly[j*3:]
				p1 := m.m_queryPoly[i*3:]
				dd.Vertex(p0, col)
				dd.Vertex(p1, col)
				j = i
				i++
			}
			dd.End()
			dd.DepthMask(true)
		}
	} else if m.m_toolMode == TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD {
		for i := 0; i < m.m_npolys; i++ {
			debug_utils.DuDebugDrawNavMeshPoly(dd, m.m_navMesh, m.m_polys[i], pathCol)
			dd.DepthMask(false)
			if m.m_parent[i] != 0 {

				p0 := make([]float64, 3)
				p1 := make([]float64, 3)
				dd.DepthMask(false)
				getPolyCenter(m.m_navMesh, m.m_parent[i], p0)
				getPolyCenter(m.m_navMesh, m.m_polys[i], p1)
				debug_utils.DuDebugDrawArc(dd, p0[0], p0[1], p0[2], p1[0], p1[1], p1[2], 0.25, 0.0, 0.4, debug_utils.DuRGBA(0, 0, 0, 128), 2.0)
				dd.DepthMask(true)
			}

			MAX_SEGS := recast.DT_VERTS_PER_POLYGON * 4
			segs := make([]float64, MAX_SEGS*6)
			refs := make([]recast.DtPolyRef, MAX_SEGS)
			nsegs, _ := m.m_navQuery.GetPolyWallSegments(m.m_polys[i], m.m_filter, segs, refs, MAX_SEGS)
			dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
			for j := 0; j < nsegs; j++ {
				s := segs[j*6:]

				// Skip too distant segments.

				_, distSqr := recast.DtDistancePtSegSqr2D(m.m_spos, s, s[3:])
				if distSqr > common.Sqr(m.m_neighbourhoodRadius) {
					continue
				}

				delta := make([]float64, 3)
				norm := make([]float64, 3)
				p0 := make([]float64, 3)
				p1 := make([]float64, 3)
				common.Vsub(delta, s[3:], s)
				common.Vmad(p0, s, delta, 0.5)
				norm[0] = delta[2]
				norm[1] = 0
				norm[2] = -delta[0]
				common.Vnormalize(norm)
				common.Vmad(p1, p0, norm, agentRadius*0.5)

				// Skip backfacing segments.
				if refs[j] != 0 {
					col := debug_utils.DuRGBA(255, 255, 255, 32)
					dd.Vertex1(s[0], s[1]+agentClimb, s[2], col)
					dd.Vertex1(s[3], s[4]+agentClimb, s[5], col)
				} else {
					col := debug_utils.DuRGBA(192, 32, 16, 192)
					if common.TriArea2D(m.m_spos, s, s[3:]) < 0.0 {
						col = debug_utils.DuRGBA(96, 32, 16, 192)
					}

					dd.Vertex1(p0[0], p0[1]+agentClimb, p0[2], col)
					dd.Vertex1(p1[0], p1[1]+agentClimb, p1[2], col)

					dd.Vertex1(s[0], s[1]+agentClimb, s[2], col)
					dd.Vertex1(s[3], s[4]+agentClimb, s[5], col)
				}
			}
			dd.End()

			dd.DepthMask(true)
		}

		if m.m_sposSet {
			dd.DepthMask(false)
			debug_utils.DuDebugDrawCircle(dd, m.m_spos[0], m.m_spos[1]+agentHeight/2, m.m_spos[2], m.m_neighbourhoodRadius, debug_utils.DuRGBA(64, 16, 0, 220), 2.0)
			dd.DepthMask(true)
		}
	}

	if m.m_nrandPoints > 0 {
		dd.Begin(debug_utils.DU_DRAW_POINTS, 6.0)
		for i := 0; i < m.m_nrandPoints; i++ {
			p := m.m_randPoints[i*3:]
			dd.Vertex1(p[0], p[1]+0.1, p[2], debug_utils.DuRGBA(220, 32, 16, 192))
		}
		dd.End()

		if m.m_randPointsInCircle && m.m_sposSet {
			debug_utils.DuDebugDrawCircle(dd, m.m_spos[0], m.m_spos[1]+agentHeight/2, m.m_spos[2], m.m_randomRadius, debug_utils.DuRGBA(64, 16, 0, 220), 2.0)
		}
	}
}

func (m *NavMeshTesterTool) handleRenderOverlay(proj, model []float64, view []int) {

	res := common.GluProject([]float64{m.m_spos[0], m.m_spos[1], m.m_spos[2]},
		model, proj, view)
	x, y := res[0], res[1]
	// Draw start and end point labels
	if m.m_sposSet && len(res) != 0 {
		m.gs.imguiDrawText(int(x), int(y-25), IMGUI_ALIGN_CENTER, "Start", imguiRGBA(0, 0, 0, 220))
	}
	res = common.GluProject([]float64{m.m_epos[0], m.m_epos[1], m.m_epos[2]},
		model, proj, view)
	x, y = res[0], res[1]
	if m.m_eposSet && len(res) != 0 {
		m.gs.imguiDrawText(int(x), int(y-25), IMGUI_ALIGN_CENTER, "End", imguiRGBA(0, 0, 0, 220))
	}

	// Tool help
	h := view[3]
	m.gs.imguiDrawText(280, h-40, IMGUI_ALIGN_LEFT, "LMB+SHIFT: Set start location  LMB: Set end location", imguiRGBA(255, 255, 255, 192))
}
func getPolyCenter(navMesh recast.IDtNavMesh, ref recast.DtPolyRef, center []float64) {
	center[0] = 0
	center[1] = 0
	center[2] = 0

	tile, poly, status := navMesh.GetTileAndPolyByRef(ref)
	if status.DtStatusFailed() {
		return
	}

	for i := 0; i < poly.VertCount; i++ {
		v := tile.Verts[poly.Verts[i]*3:]
		center[0] += v[0]
		center[1] += v[1]
		center[2] += v[2]
	}
	s := 1.0 / float64(poly.VertCount)
	center[0] *= s
	center[1] *= s
	center[2] *= s
}
func inRange(v1, v2 []float64, r, h float64) bool {
	dx := v2[0] - v1[0]
	dy := v2[1] - v1[1]
	dz := v2[2] - v1[2]
	return (dx*dx+dz*dz) < r*r && math.Abs(dy) < h
}

func (m *NavMeshTesterTool) handleToggle() {
	// TODO: merge separate to a path iterator. Use same code in recalc() too.
	if m.m_toolMode != TOOLMODE_PATHFIND_FOLLOW {
		return
	}

	if !m.m_sposSet || !m.m_eposSet || m.m_startRef == 0 || m.m_endRef == 0 {
		return
	}

	STEP_SIZE := 0.5
	SLOP := 0.01

	if m.m_pathIterNum == 0 {
		m.m_npolys, _ = m.m_navQuery.FindPath(m.m_startRef, m.m_endRef, m.m_spos, m.m_epos, m.m_filter, m.m_polys, MAX_POLYS)
		m.m_nsmoothPath = 0

		m.m_pathIterPolyCount = m.m_npolys
		if m.m_pathIterPolyCount != 0 {
			copy(m.m_pathIterPolys, m.m_polys[:m.m_pathIterPolyCount])
		}

		if m.m_pathIterPolyCount != 0 {
			// Iterate over the path to find smooth path on the detail mesh surface.
			tmp := false
			m.m_navQuery.ClosestPointOnPoly(m.m_startRef, m.m_spos, m.m_iterPos, &tmp)
			m.m_navQuery.ClosestPointOnPoly(m.m_pathIterPolys[m.m_pathIterPolyCount-1], m.m_epos, m.m_targetPos, &tmp)

			m.m_nsmoothPath = 0

			copy(common.GetVs3(m.m_smoothPath, m.m_nsmoothPath), m.m_iterPos)
			m.m_nsmoothPath++
		}
	}

	copy(m.m_prevIterPos, m.m_iterPos)

	m.m_pathIterNum++

	if m.m_pathIterPolyCount == 0 {
		return
	}

	if m.m_nsmoothPath >= MAX_SMOOTH {
		return
	}

	// Move towards target a small advancement at a time until target reached or
	// when ran out of memory to store the path.

	// Find location to steer towards.
	steerPos := make([]float64, 3)
	var steerPosFlag int
	var steerPosRef recast.DtPolyRef

	if !getSteerTarget(m.m_navQuery, m.m_iterPos, m.m_targetPos, SLOP,
		m.m_pathIterPolys, m.m_pathIterPolyCount, steerPos, &steerPosFlag, steerPosRef,
		m.m_steerPoints, &m.m_steerPointCount) {
		return
	}

	copy(m.m_steerPos, steerPos)

	endOfPath := false
	if steerPosFlag&recast.DT_STRAIGHTPATH_END != 0 {
		endOfPath = true
	}
	offMeshConnection := false
	if steerPosFlag&recast.DT_STRAIGHTPATH_OFFMESH_CONNECTION != 0 {
		offMeshConnection = true
	}
	// Find movement delta.
	delta := make([]float64, 3)
	common.Vsub(delta, steerPos, m.m_iterPos)
	length := math.Sqrt(common.Vdot(delta, delta))
	// If the steer target is end of path or off-mesh link, do not move past the location.
	if (endOfPath || offMeshConnection) && length < STEP_SIZE {
		length = 1
	} else {
		length = STEP_SIZE / length
	}

	moveTgt := make([]float64, 3)
	common.Vmad(moveTgt, m.m_iterPos, delta, length)

	// Move
	result := make([]float64, 3)
	visited := make([]recast.DtPolyRef, 16)
	nvisited := 0
	m.m_navQuery.MoveAlongSurface(m.m_pathIterPolys[0], m.m_iterPos, moveTgt, m.m_filter,
		result, visited, &nvisited, 16)
	m.m_pathIterPolyCount = recast.DtMergeCorridorStartMoved(m.m_pathIterPolys, m.m_pathIterPolyCount,
		MAX_POLYS, visited, nvisited)
	m.m_pathIterPolyCount = fixupShortcuts(m.m_pathIterPolys, m.m_pathIterPolyCount, m.m_navQuery)
	h, _ := m.m_navQuery.GetPolyHeight(m.m_pathIterPolys[0], result)
	result[1] = h
	copy(m.m_iterPos, result)

	// Handle end of path and off-mesh links when close enough.
	if endOfPath && inRange(m.m_iterPos, steerPos, SLOP, 1.0) {
		// Reached end of path.
		copy(m.m_iterPos, m.m_targetPos)
		if m.m_nsmoothPath < MAX_SMOOTH {
			copy(common.GetVs3(m.m_smoothPath, m.m_nsmoothPath), m.m_iterPos)
			m.m_nsmoothPath++
		}
		return
	} else if offMeshConnection && inRange(m.m_iterPos, steerPos, SLOP, 1.0) {
		// Reached off-mesh connection.
		startPos := make([]float64, 3)
		endPos := make([]float64, 3)
		// Advance the path up to and over the off-mesh connection.
		var prevRef recast.DtPolyRef
		polyRef := m.m_pathIterPolys[0]
		npos := 0
		for npos < m.m_pathIterPolyCount && polyRef != steerPosRef {
			prevRef = polyRef
			polyRef = m.m_pathIterPolys[npos]
			npos++
		}
		for i := npos; i < m.m_pathIterPolyCount; i++ {
			m.m_pathIterPolys[i-npos] = m.m_pathIterPolys[i]
		}

		m.m_pathIterPolyCount -= npos

		// Handle the connection.
		status := m.m_navMesh.GetOffMeshConnectionPolyEndPoints(prevRef, polyRef, startPos, endPos)
		if status.DtStatusSucceed() {
			if m.m_nsmoothPath < MAX_SMOOTH {
				copy(common.GetVs3(m.m_smoothPath, m.m_nsmoothPath), startPos)
				m.m_nsmoothPath++
				// Hack to make the dotted path not visible during off-mesh connection.
				if m.m_nsmoothPath&1 != 0 {
					copy(common.GetVs3(m.m_smoothPath, m.m_nsmoothPath), startPos)
					m.m_nsmoothPath++
				}
			}
			// Move position at the other side of the off-mesh link.
			copy(m.m_iterPos, endPos)

			eh, _ := m.m_navQuery.GetPolyHeight(m.m_pathIterPolys[0], m.m_iterPos)
			m.m_iterPos[1] = eh
		}
	}

	// Store results.
	if m.m_nsmoothPath < MAX_SMOOTH {
		copy(common.GetVs3(m.m_smoothPath, m.m_nsmoothPath), m.m_iterPos)
		m.m_nsmoothPath++
	}
}

func (m *NavMeshTesterTool) handleStep() {

}

func newNavMeshTesterTool(gs *guiState) *NavMeshTesterTool {
	m := &NavMeshTesterTool{
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

func (m *NavMeshTesterTool) init(sample *Sample) {
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

func (m *NavMeshTesterTool) handleMenu() {
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
func (m *NavMeshTesterTool) handleClick(s, p []float64, shift bool) {
	if shift {
		m.m_sposSet = true
		copy(m.m_spos, p)
	} else {
		m.m_eposSet = true
		copy(m.m_epos, p)
	}
	m.recalc()
}

func (m *NavMeshTesterTool) handleUpdate(dt float64) {
	if m.m_toolMode == TOOLMODE_PATHFIND_SLICED {
		if m.m_pathFindStatus.DtStatusInProgress() {
			_, m.m_pathFindStatus = m.m_navQuery.UpdateSlicedFindPath(1)
		}
		if m.m_pathFindStatus.DtStatusSucceed() {
			m.m_npolys, _ = m.m_navQuery.FinalizeSlicedFindPath(m.m_polys, MAX_POLYS)
			m.m_nstraightPath = 0
			if m.m_npolys != 0 {
				// In case of partial path, make sure the end point is clamped to the last polygon.
				epos := make([]float64, 3)
				copy(epos, m.m_epos)
				if m.m_polys[m.m_npolys-1] != m.m_endRef {
					var tmp bool
					m.m_navQuery.ClosestPointOnPoly(m.m_polys[m.m_npolys-1], m.m_epos, epos, &tmp)
				}

				m.m_nstraightPath, _ = m.m_navQuery.FindStraightPath(m.m_spos, epos, m.m_polys, m.m_npolys,
					m.m_straightPath, m.m_straightPathFlags,
					m.m_straightPathPolys, MAX_POLYS, recast.DT_STRAIGHTPATH_ALL_CROSSINGS)
			}

			m.m_pathFindStatus = recast.DT_FAILURE
		}
	}
}

func (m *NavMeshTesterTool) recalc() {

}
func (m *NavMeshTesterTool) Type() int { return TOOL_NAVMESH_TESTER }
func (m *NavMeshTesterTool) reset() {
	m.m_startRef = 0
	m.m_endRef = 0
	m.m_npolys = 0
	m.m_nstraightPath = 0
	m.m_nsmoothPath = 0
	m.m_hitPos = make([]float64, len(m.m_hitPos))
	m.m_hitNormal = make([]float64, len(m.m_hitNormal))
	m.m_distanceToWall = 0
}

// This function checks if the path has a small U-turn, that is,
// a polygon further in the path is adjacent to the first polygon
// in the path. If that happens, a shortcut is taken.
// This can happen if the target (T) location is at tile boundary,
// and we're (S) approaching it parallel to the tile edge.
// The choice at the vertex can be arbitrary,
//
//	+---+---+
//	|:::|:::|
//	+-S-+-T-+
//	|:::|   | <-- the step can end up in here, resulting U-turn path.
//	+---+---+
func fixupShortcuts(path []recast.DtPolyRef, npath int, navQuery recast.NavMeshQuery) int {
	if npath < 3 {
		return npath
	}

	// Get connected polygons
	maxNeis := 16
	neis := make([]recast.DtPolyRef, maxNeis)
	nneis := 0

	tile, poly, status := navQuery.GetAttachedNavMesh().GetTileAndPolyByRef(path[0])
	if status.DtStatusFailed() {
		return npath
	}

	for k := poly.FirstLink; k != recast.DT_NULL_LINK; k = tile.Links[k].Next {
		link := tile.Links[k]
		if link.Ref != 0 {
			if nneis < maxNeis {
				neis[nneis] = link.Ref
				nneis++
			}

		}
	}

	// If any of the neighbour polygons is within the next few polygons
	// in the path, short cut to that polygon directly.
	maxLookAhead := 6
	cut := 0
	for i := common.Min(maxLookAhead, npath) - 1; i > 1 && cut == 0; i-- {
		for j := 0; j < nneis; j++ {
			if path[i] == neis[j] {
				cut = i
				break
			}
		}
	}
	if cut > 1 {
		offset := cut - 1
		npath -= offset
		for i := 1; i < npath; i++ {
			path[i] = path[i+offset]
		}

	}

	return npath
}

func getSteerTarget(navQuery recast.NavMeshQuery, startPos, endPos []float64,
	minTargetDist float64,
	path []recast.DtPolyRef, pathSize int,
	steerPos []float64, steerPosFlag *int, steerPosRef recast.DtPolyRef,
	outPoints []float64, outPointCount *int) bool {
	// Find steer target.
	MAX_STEER_POINTS := 3
	steerPath := make([]float64, MAX_STEER_POINTS*3)
	steerPathFlags := make([]int, MAX_STEER_POINTS)
	steerPathPolys := make([]recast.DtPolyRef, MAX_STEER_POINTS)
	nsteerPath, _ := navQuery.FindStraightPath(startPos, endPos, path, pathSize,
		steerPath, steerPathFlags, steerPathPolys, MAX_STEER_POINTS)
	if nsteerPath == 0 {
		return false
	}

	if len(outPoints) > 0 && outPointCount != nil {
		*outPointCount = nsteerPath
		for i := 0; i < nsteerPath; i++ {
			copy(common.GetVs3(outPoints, i), common.GetVs3(steerPath, i))
		}

	}

	// Find vertex far enough to steer to.
	ns := 0
	for ns < nsteerPath {

		// Stop at Off-Mesh link or when point is further than slop away.
		if (steerPathFlags[ns]&recast.DT_STRAIGHTPATH_OFFMESH_CONNECTION != 0) ||
			!inRange(common.GetVs3(steerPath, ns), startPos, minTargetDist, 1000.0) {
			break
		}

		ns++
	}
	// Failed to find good point to steer to.
	if ns >= nsteerPath {
		return false
	}

	copy(steerPos, common.GetVs3(steerPath, ns))
	steerPos[1] = startPos[1]
	*steerPosFlag = steerPathFlags[ns]
	steerPosRef = steerPathPolys[ns]

	return true
}
