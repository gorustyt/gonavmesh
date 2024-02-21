package mesh

import (
	"fmt"
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/demo/config"
	"github.com/gorustyt/gonavmesh/detour"
	"github.com/gorustyt/gonavmesh/detour_crowd"
	"math"
	"time"
)

func isectSegAABB(sp, sq []float64,
	amin, amax []float64,
	tmin, tmax float64) bool {
	EPS := 1e-6
	d := make([]float64, 3)
	common.Vsub(d, sq, sp)
	tmin = 0               // set to -FLT_MAX to get first hit on line
	tmax = math.MaxFloat64 // set to max distance ray can travel (for segment)

	// For all three slabs
	for i := 0; i < 3; i++ {
		if math.Abs(d[i]) < EPS {
			// Ray is parallel to slab. No hit if origin not within slab
			if sp[i] < amin[i] || sp[i] > amax[i] {
				return false
			}

		} else {
			// Compute intersection t value of ray with near and far plane of slab
			ood := 1.0 / d[i]
			t1 := (amin[i] - sp[i]) * ood
			t2 := (amax[i] - sp[i]) * ood
			// Make t1 be intersection with near plane, t2 with far plane
			if t1 > t2 {
				t1, t2 = t2, t1
			}
			// Compute the intersection of slab intersections intervals
			if t1 > tmin {
				tmin = t1
			}
			if t2 < tmax {
				tmax = t2
			}
			// Exit with no collision as soon as slab intersection becomes empty
			if tmin > tmax {
				return false
			}
		}
	}

	return true
}

func getAgentBounds(ag *detour_crowd.DtCrowdAgent, bmin, bmax []float64) {
	p := common.SliceTToSlice[float32, float64](ag.Npos[:])
	r := float64(ag.Params.Radius)
	h := float64(ag.Params.Height)
	bmin[0] = p[0] - r
	bmin[1] = p[1]
	bmin[2] = p[2] - r
	bmax[0] = p[0] + r
	bmax[1] = p[1] + h
	bmax[2] = p[2] + r
}

const (
	AGENT_MAX_TRAIL = 64
	MAX_AGENTS      = 128
)

type AgentTrail struct {
	trail  [AGENT_MAX_TRAIL * 3]float64
	htrail int
}

type CrowdToolState struct {
	cfg *config.Config

	m_sample    *Sample
	m_nav       detour.IDtNavMesh
	m_crowd     *detour_crowd.DtCrowd
	m_targetPos []float64
	m_targetRef detour.DtPolyRef

	m_agentDebug *detour_crowd.DtCrowdAgentDebugInfo
	m_vod        *detour_crowd.DtObstacleAvoidanceDebugData

	m_trails []AgentTrail

	m_crowdTotalTime   ValueHistory
	m_crowdSampleCount ValueHistory
}

func newCrowdToolState(ctx *Content) *CrowdToolState {
	d := &CrowdToolState{
		m_trails:    make([]AgentTrail, MAX_AGENTS),
		m_targetPos: make([]float64, 3),
		cfg:         ctx.GetConfig(),
	}
	d.m_vod = detour_crowd.NewDtObstacleAvoidanceDebugData(2048)
	d.m_agentDebug.Idx = -1
	d.m_agentDebug.Vod = d.m_vod
	return d
}

func (c *CrowdToolState) init(sample *Sample) {
	if c.m_sample != sample {
		c.m_sample = sample
	}

	nav := c.m_sample.getNavMesh()
	crowd := c.m_sample.getCrowd()
	if crowd == nil {
		crowd = detour_crowd.NewDtCrowd(MAX_AGENTS, float32(c.m_sample.getAgentRadius()), nav)
	}
	if nav != nil && crowd != nil && (c.m_nav != nav || c.m_crowd != crowd) {
		c.m_nav = nav
		c.m_crowd = crowd
		// Make polygons with 'disabled' flag invalid.
		crowd.GetEditableFilter(0).SetExcludeFlags(config.SAMPLE_POLYFLAGS_DISABLED)

		// Setup local avoidance params to different qualities.
		params := crowd.GetObstacleAvoidanceParams(0)
		// Use mostly default settings, copy from dtCrowd.

		// Low (11)
		params.VelBias = 0.5
		params.AdaptiveDivs = 5
		params.AdaptiveRings = 2
		params.AdaptiveDepth = 1
		crowd.SetObstacleAvoidanceParams(0, params)

		// Medium (22)
		params.VelBias = 0.5
		params.AdaptiveDivs = 5
		params.AdaptiveRings = 2
		params.AdaptiveDepth = 2
		crowd.SetObstacleAvoidanceParams(1, params)

		// Good (45)
		params.VelBias = 0.5
		params.AdaptiveDivs = 7
		params.AdaptiveRings = 2
		params.AdaptiveDepth = 3
		crowd.SetObstacleAvoidanceParams(2, params)

		// High (66)
		params.VelBias = 0.5
		params.AdaptiveDivs = 7
		params.AdaptiveRings = 3
		params.AdaptiveDepth = 3

		crowd.SetObstacleAvoidanceParams(3, params)
	}
}
func (c *CrowdToolState) reset() {
}

func (c *CrowdToolState) handleRenderOverlay(proj, model []float64, view []int) {
	res := common.GluProject([]float64{c.m_targetPos[0], c.m_targetPos[1], c.m_targetPos[2]}, model, proj, view)
	x, y := int(res[0]), int(res[1])
	// Draw start and end point labels
	if c.m_targetRef != 0 && len(res) != 0 {
		panic("not impl")
		_ = x
		_ = y
		//c.gs.imguiDrawText(x, y+25, IMGUI_ALIGN_CENTER, "TARGET", imguiRGBA(0, 0, 0, 220))
	}
	if c.cfg.ToolsConfig.GetShowNodes() {
		crowd := c.m_sample.getCrowd()
		if crowd != nil && crowd.GetPathQueue() != nil {
			navquery := crowd.GetPathQueue().GetNavQuery()
			pool := navquery.GetNodePool()
			if pool != nil {
				off := 0.5
				for i := int32(0); i < pool.GetHashSize(); i++ {
					for j := pool.GetFirst(int(i)); j != detour.DT_NULL_IDX; j = pool.GetNext(int(j)) {
						node := pool.GetNodeAtIdx(int(j) + 1)
						if node == nil {
							continue
						}
						res = common.GluProject([]float64{float64(node.Pos[0]), float64(node.Pos[1]) + off, float64(node.Pos[2])}, model, proj, view)
						x, y = int(res[0]), int(res[1])
						if len(res) > 0 {
							heuristic := node.Total // - node.cost;
							label := fmt.Sprintf("%.2f", heuristic)
							panic(fmt.Sprintf("not imple:%v", label))
							//c.gs.imguiDrawText(x, y+15, IMGUI_ALIGN_CENTER, label, imguiRGBA(0, 0, 0, 220))
						}
					}
				}
			}
		}
	}

	if c.cfg.ToolsConfig.GetShowLabels() {
		crowd := c.m_sample.getCrowd()
		if crowd != nil {
			for i := 0; i < crowd.GetAgentCount(); i++ {
				ag := crowd.GetAgent(i)
				if !ag.Active {
					continue
				}
				pos := ag.Npos
				h := ag.Params.Height

				if res := common.GluProject([]float64{float64(pos[0]), float64(pos[1] + h), float64(pos[2])}, model, proj, view); len(res) != 0 {
					x, y = int(res[0]), int(res[1])
					label := fmt.Sprintf("%d", i)
					//c.gs.imguiDrawText(x, y+15, IMGUI_ALIGN_CENTER, label, imguiRGBA(0, 0, 0, 220))
					panic(fmt.Sprintf("not imple:%v", label))
				}
			}
		}
	}
	if c.m_agentDebug.Idx != -1 {
		crowd := c.m_sample.getCrowd()
		if crowd != nil {
			for i := 0; i < crowd.GetAgentCount(); i++ {
				if c.cfg.ToolsConfig.GetShowDetailAll() == false && i != int(c.m_agentDebug.Idx) {
					continue
				}

				ag := crowd.GetAgent(i)
				if !ag.Active {
					continue
				}

				radius := ag.Params.Radius
				if c.cfg.ToolsConfig.GetShowNeis() {
					for j := int32(0); j < ag.Nneis; j++ {
						nei := crowd.GetAgent(ag.Neis[j].Idx)
						if !nei.Active {
							continue
						}

						if res := common.GluProject([]float64{float64(nei.Npos[0]), float64(nei.Npos[1]) + float64(radius), float64(nei.Npos[2])}, model, proj, view); len(res) != 0 {
							x, y = int(res[0]), int(res[1])
							label := fmt.Sprintf("%.3f", ag.Neis[j].Dist)
							panic(fmt.Sprintf("not imple:%v", label))
							//c.gs.imguiDrawText(x, y+15, IMGUI_ALIGN_CENTER, label, imguiRGBA(255, 255, 255, 220))
						}
					}
				}
			}
		}
	}

	if c.cfg.ToolsConfig.GetShowPerfGraph() {
		panic("not impl")
	}

}
func (c *CrowdToolState) handleRender() {
	dd := c.m_sample.getDebugDraw()
	rad := c.m_sample.getAgentRadius()

	nav := c.m_sample.getNavMesh()
	crowd := c.m_sample.getCrowd()
	if nav == nil || crowd == nil {
		return
	}

	if c.cfg.ToolsConfig.GetShowNodes() && crowd.GetPathQueue() != nil {
		navquery := crowd.GetPathQueue().GetNavQuery()
		if navquery != nil {
			debug_utils.DuDebugDrawNavMeshNodes(dd, navquery)
		}

	}

	dd.DepthMask(false)

	// Draw paths
	if c.cfg.ToolsConfig.GetShowPath() {
		for i := 0; i < crowd.GetAgentCount(); i++ {
			if c.cfg.ToolsConfig.GetShowDetailAll() == false && i != int(c.m_agentDebug.Idx) {
				continue
			}

			ag := crowd.GetAgent(i)
			if !ag.Active {
				continue
			}

			path := ag.Corridor.GetPath()
			npath := ag.Corridor.GetPathCount()
			for j := 0; j < npath; j++ {
				debug_utils.DuDebugDrawNavMeshPoly(dd, nav, path[j], debug_utils.DuRGBA(255, 255, 255, 24))
			}

		}
	}

	if c.m_targetRef > 0 {
		debug_utils.DuDebugDrawCross(dd, c.m_targetPos[0], c.m_targetPos[1]+0.1, c.m_targetPos[2], rad, debug_utils.DuRGBA(255, 255, 255, 192), 2.0)
	}

	// Occupancy grid.
	if c.cfg.ToolsConfig.GetShowGrid() {
		gridy := math.SmallestNonzeroFloat64
		for i := 0; i < crowd.GetAgentCount(); i++ {
			ag := crowd.GetAgent(i)
			if !ag.Active {
				continue
			}
			pos := ag.Corridor.GetPos()
			gridy = max(gridy, float64(pos[1]))
		}
		gridy += 1.0

		dd.Begin(debug_utils.DU_DRAW_QUADS)
		grid := crowd.GetGrid()
		bounds := grid.GetBounds()
		cs := float64(grid.GetCellSize())
		for y := bounds[1]; y <= bounds[3]; y++ {
			for x := bounds[0]; x <= bounds[2]; x++ {
				count := grid.GetItemCountAt(int(x), int(y))
				if count == 0 {
					continue
				}
				col := debug_utils.DuRGBA(128, 0, 0, min(count*40, 255))
				dd.Vertex1(float64(x)*cs, gridy, float64(y)*cs, col)
				dd.Vertex1(float64(x)*cs, gridy, float64(y)*cs+cs, col)
				dd.Vertex1(float64(x)*cs+cs, gridy, float64(y)*cs+cs, col)
				dd.Vertex1(float64(x)*cs+cs, gridy, float64(y)*cs, col)
			}
		}
		dd.End()
	}

	// Trail
	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}

		trail := c.m_trails[i]
		pos := ag.Npos

		dd.Begin(debug_utils.DU_DRAW_LINES, 3.0)
		prev := make([]float64, 3)
		preva := 1.0
		copy(prev, common.SliceTToSlice[float32, float64](pos[:]))
		for j := 0; j < AGENT_MAX_TRAIL-1; j++ {
			idx := (trail.htrail + AGENT_MAX_TRAIL - j) % AGENT_MAX_TRAIL
			v := common.GetVert3(trail.trail[:], idx)
			a := float64(1-j) / AGENT_MAX_TRAIL
			dd.Vertex1(prev[0], prev[1]+0.1, prev[2], debug_utils.DuRGBA(0, 0, 0, (int)(128*preva)))
			dd.Vertex1(v[0], v[1]+0.1, v[2], debug_utils.DuRGBA(0, 0, 0, (int)(128*a)))
			preva = a
			copy(prev, v)
		}
		dd.End()

	}

	// Corners & co
	for i := 0; i < crowd.GetAgentCount(); i++ {
		if c.cfg.ToolsConfig.GetShowDetailAll() == false && i != int(c.m_agentDebug.Idx) {
			continue
		}

		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}

		radius := float64(ag.Params.Radius)
		pos := common.SliceTToSlice[float32, float64](ag.Npos[:])

		if c.cfg.ToolsConfig.GetShowCorners() {
			if ag.Ncorners != 0 {
				dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
				for j := int32(0); j < ag.Ncorners; j++ {

					va := common.SliceTToSlice[float32, float64](ag.CornerVerts[(j-1)*3:])
					if j == 0 {
						va = pos[:]
					}
					vb := common.SliceTToSlice[float32, float64](ag.CornerVerts[j*3:])
					dd.Vertex1(va[0], va[1]+radius, va[2], debug_utils.DuRGBA(128, 0, 0, 192))
					dd.Vertex1(vb[0], vb[1]+radius, vb[2], debug_utils.DuRGBA(128, 0, 0, 192))
				}
				if ag.Ncorners > 0 && ag.CornerFlags[ag.Ncorners-1]&detour.DT_STRAIGHTPATH_OFFMESH_CONNECTION > 0 {
					v := common.SliceTToSlice[float32, float64](common.GetVert3(
						ag.CornerVerts[:], (ag.Ncorners - 1)))
					dd.Vertex1(v[0], v[1], v[2], debug_utils.DuRGBA(192, 0, 0, 192))
					dd.Vertex1(v[0], v[1]+radius*2, v[2], debug_utils.DuRGBA(192, 0, 0, 192))
				}

				dd.End()

				if c.cfg.ToolsConfig.GetAnticipateTurns() {
					/*					float dvel[3], pos[3];
										calcSmoothSteerDirection(ag.pos, ag.cornerVerts, ag.ncorners, dvel);
										pos[0] = ag.pos[0] + dvel[0];
										pos[1] = ag.pos[1] + dvel[1];
										pos[2] = ag.pos[2] + dvel[2];

										const float off = ag.radius+0.1f;
										const float* tgt = &ag.cornerVerts[0];
										const float y = ag.pos[1]+off;

										dd.begin(DU_DRAW_LINES, 2.0f);

										dd.vertex(ag.pos[0],y,ag.pos[2], duRGBA(255,0,0,192));
										dd.vertex(pos[0],y,pos[2], duRGBA(255,0,0,192));

										dd.vertex(pos[0],y,pos[2], duRGBA(255,0,0,192));
										dd.vertex(tgt[0],y,tgt[2], duRGBA(255,0,0,192));

										dd.end();*/
				}
			}
		}

		if c.cfg.ToolsConfig.GetShowCollisionSegments() {
			tmp := ag.Boundary.GetCenter()
			center := common.SliceTToSlice[float32, float64](tmp[:])
			debug_utils.DuDebugDrawCross(dd, center[0], center[1]+radius, center[2], 0.2, debug_utils.DuRGBA(192, 0, 128, 255), 2.0)
			debug_utils.DuDebugDrawCircle(dd, center[0], center[1]+radius, center[2], float64(ag.Params.CollisionQueryRange),
				debug_utils.DuRGBA(192, 0, 128, 128), 2.0)

			dd.Begin(debug_utils.DU_DRAW_LINES, 3.0)
			for j := 0; j < ag.Boundary.GetSegmentCount(); j++ {
				tmp := ag.Boundary.GetSegment(j)
				s := common.SliceTToSlice[float32, float64](tmp[:])
				col := debug_utils.DuRGBA(192, 0, 128, 192)
				if common.TriArea2D(pos[:], s[:3], s[3:]) < 0.0 {
					col = debug_utils.DuDarkenCol(col)
				}

				debug_utils.DuAppendArrow(dd, s[0], s[1]+0.2, s[2], s[3], s[4]+0.2, s[5], 0.0, 0.3, col)
			}
			dd.End()
		}

		if c.cfg.ToolsConfig.GetShowNeis() {
			debug_utils.DuDebugDrawCircle(dd, pos[0], pos[1]+radius, pos[2], float64(ag.Params.CollisionQueryRange),
				debug_utils.DuRGBA(0, 192, 128, 128), 2.0)

			dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
			for j := int32(0); j < ag.Nneis; j++ {
				// Get 'n'th active agent.
				// TODO: fix this properly.
				nei := crowd.GetAgent(ag.Neis[j].Idx)
				if nei != nil {
					dd.Vertex1(pos[0], pos[1]+radius, pos[2], debug_utils.DuRGBA(0, 192, 128, 128))
					dd.Vertex1(float64(nei.Npos[0]), float64(nei.Npos[1])+radius, float64(nei.Npos[2]), debug_utils.DuRGBA(0, 192, 128, 128))
				}
			}
			dd.End()
		}

		if c.cfg.ToolsConfig.GetShowOpt() {
			dd.Begin(debug_utils.DU_DRAW_LINES, 2.0)
			dd.Vertex1(float64(c.m_agentDebug.OptStart[0]), float64(c.m_agentDebug.OptStart[1]+0.3), float64(c.m_agentDebug.OptStart[2]), debug_utils.DuRGBA(0, 128, 0, 192))
			dd.Vertex1(float64(c.m_agentDebug.OptEnd[0]), float64(c.m_agentDebug.OptEnd[1]+0.3), float64(c.m_agentDebug.OptEnd[2]), debug_utils.DuRGBA(0, 128, 0, 192))
			dd.End()
		}
	}

	// Agent cylinders.
	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}

		radius := float64(ag.Params.Radius)
		pos := common.SliceTToSlice[float32, float64](ag.Npos[:])

		col := debug_utils.DuRGBA(0, 0, 0, 32)
		if int(c.m_agentDebug.Idx) == i {
			col = debug_utils.DuRGBA(255, 0, 0, 128)
		}

		debug_utils.DuDebugDrawCircle(dd, pos[0], pos[1], pos[2], radius, col, 2.0)
	}

	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}

		height := float64(ag.Params.Height)
		radius := float64(ag.Params.Radius)
		pos := common.SliceTToSlice[float32, float64](ag.Npos[:])

		col := debug_utils.DuRGBA(220, 220, 220, 128)
		if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_REQUESTING || ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE {
			col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(128, 0, 255, 128), 32)
		} else if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH {
			col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(128, 0, 255, 128), 128)
		} else if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_FAILED {
			col = debug_utils.DuRGBA(255, 32, 16, 128)
		} else if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_VELOCITY {
			col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(64, 255, 0, 128), 128)
		}

		debug_utils.DuDebugDrawCylinder(dd, pos[0]-radius, pos[1]+radius*0.1, pos[2]-radius,
			pos[0]+radius, pos[1]+height, pos[2]+radius, col)
	}

	if c.cfg.ToolsConfig.GetShowVO() {
		for i := 0; i < crowd.GetAgentCount(); i++ {
			if c.cfg.ToolsConfig.GetShowDetailAll() == false && i != int(c.m_agentDebug.Idx) {
				continue
			}

			ag := crowd.GetAgent(i)
			if !ag.Active {
				continue
			}

			// Draw detail about agent sela
			vod := c.m_agentDebug.Vod

			dx := float64(ag.Npos[0])
			dy := float64(ag.Npos[1] + ag.Params.Height)
			dz := float64(ag.Npos[2])

			debug_utils.DuDebugDrawCircle(dd, dx, dy, dz, float64(ag.Params.MaxSpeed), debug_utils.DuRGBA(255, 255, 255, 64), 2.0)

			dd.Begin(debug_utils.DU_DRAW_QUADS)
			for j := 0; j < vod.GetSampleCount(); j++ {
				p := common.SliceTToSlice[float32, float64](vod.GetSampleVelocity(j))
				sr := float64(vod.GetSampleSize(j))
				pen := vod.GetSamplePenalty(j)
				pen2 := vod.GetSamplePreferredSidePenalty(j)
				col := debug_utils.DuLerpCol(debug_utils.DuRGBA(255, 255, 255, 220), debug_utils.DuRGBA(128, 96, 0, 220), (int)(pen*255))
				col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(128, 0, 0, 220), (int)(pen2*128))
				dd.Vertex1(dx+p[0]-sr, dy, dz+p[2]-sr, col)
				dd.Vertex1(dx+p[0]-sr, dy, dz+p[2]+sr, col)
				dd.Vertex1(dx+p[0]+sr, dy, dz+p[2]+sr, col)
				dd.Vertex1(dx+p[0]+sr, dy, dz+p[2]-sr, col)
			}
			dd.End()
		}
	}

	// Velocity stuff.
	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}

		radius := float64(ag.Params.Radius)
		height := float64(ag.Params.Height)
		pos := common.SliceTToSlice[float32, float64](ag.Npos[:])
		vel := common.SliceTToSlice[float32, float64](ag.Vel[:])
		dvel := common.SliceTToSlice[float32, float64](ag.Dvel[:])

		col := debug_utils.DuRGBA(220, 220, 220, 192)
		if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_REQUESTING || ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_WAITING_FOR_QUEUE {
			col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(128, 0, 255, 192), 32)
		} else if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_WAITING_FOR_PATH {
			col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(128, 0, 255, 192), 128)
		} else if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_FAILED {
			col = debug_utils.DuRGBA(255, 32, 16, 192)
		} else if ag.TargetState == detour_crowd.DT_CROWDAGENT_TARGET_VELOCITY {
			col = debug_utils.DuLerpCol(col, debug_utils.DuRGBA(64, 255, 0, 192), 128)
		}

		debug_utils.DuDebugDrawCircle(dd, pos[0], pos[1]+height, pos[2], radius, col, 2.0)
		tmp := 1.0
		if int(c.m_agentDebug.Idx) == i {
			tmp = 2.0
		}
		debug_utils.DuDebugDrawArrow(dd, pos[0], pos[1]+height, pos[2],
			pos[0]+dvel[0], pos[1]+height+dvel[1], pos[2]+dvel[2],
			0.0, 0.4, debug_utils.DuRGBA(0, 192, 255, 192), tmp)

		debug_utils.DuDebugDrawArrow(dd, pos[0], pos[1]+height, pos[2],
			pos[0]+vel[0], pos[1]+height+vel[1], pos[2]+vel[2],
			0.0, 0.4, debug_utils.DuRGBA(0, 0, 0, 160), 2.0)
	}

	dd.DepthMask(true)
}

func (c *CrowdToolState) handleUpdate(dt float64) {
	//if c.m_run {
	//	c.updateTick(dt)
	//}

}
func (c *CrowdToolState) updateTick(dt float64) {
	if c.m_sample == nil {
		return
	}
	nav := c.m_sample.getNavMesh()
	crowd := c.m_sample.getCrowd()
	if nav == nil || crowd == nil {
		return
	}

	now := time.Now()

	crowd.Update(float32(dt), c.m_agentDebug)

	endTime := time.Since(now)

	// Update agent trails
	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		trail := c.m_trails[i]
		if !ag.Active {
			continue
		}

		// Update agent movement trail.
		trail.htrail = (trail.htrail + 1) % AGENT_MAX_TRAIL
		copy(common.GetVert3(trail.trail[:], trail.htrail), common.SliceTToSlice[float32, float64](ag.Npos[:]))
	}

	c.m_agentDebug.Vod.NormalizeSamples()
	_ = endTime
	//TODO
	//c.m_crowdSampleCount.addSample(float64(crowd.GetVelocitySampleCount()))
	//c.m_crowdTotalTime.addSample(float64(endTime.Milliseconds()))
}

func (c *CrowdToolState) addAgent(p []float64) {
	if c.m_sample == nil {
		return
	}
	crowd := c.m_sample.getCrowd()

	ap := &detour_crowd.DtCrowdAgentParams{}
	ap.Radius = float32(c.m_sample.getAgentRadius())
	ap.Height = float32(c.m_sample.getAgentHeight())
	ap.MaxAcceleration = 8.0
	ap.MaxSpeed = 3.5
	ap.CollisionQueryRange = ap.Radius * 12.0
	ap.PathOptimizationRange = ap.Radius * 30.0
	ap.UpdateFlags = c.cfg.ToolsConfig.GetExpandOptionsUpdateFlags()

	ap.ObstacleAvoidanceType = int(c.cfg.ToolsConfig.ObstacleAvoidanceType)
	ap.SeparationWeight = float32(c.cfg.ToolsConfig.SeparationWeight)

	idx := crowd.AddAgent(common.SliceTToSlice[float64, float32](p), ap)
	if idx != -1 {
		if c.m_targetRef != 0 {
			crowd.RequestMoveTarget(idx, c.m_targetRef, common.SliceTToSlice[float64, float32](c.m_targetPos))
		}

		// Init trail
		trail := c.m_trails[idx]
		for i := 0; i < AGENT_MAX_TRAIL; i++ {
			copy(common.GetVert3(trail.trail[:], i), p)
		}

		trail.htrail = 0
	}
}
func (c *CrowdToolState) removeAgent(idx int) {
	if c.m_sample == nil {
		return
	}
	crowd := c.m_sample.getCrowd()

	crowd.RemoveAgent(idx)

	if idx == int(c.m_agentDebug.Idx) {
		c.m_agentDebug.Idx = -1
	}

}
func (c *CrowdToolState) hilightAgent(idx int) {
	c.m_agentDebug.Idx = int32(idx)
}

func calcVel(vel, pos, tgt []float64, speed float64) {
	common.Vsub(vel, tgt, pos)
	vel[1] = 0.0
	common.Vnormalize(vel)
	common.Vscale(vel, vel, speed)
}

func (c *CrowdToolState) setMoveTarget(p []float64, adjust bool) {
	if c.m_sample == nil {
		return
	}

	// Find nearest point on navmesh and set move request to that location.
	navquery := c.m_sample.getNavMeshQuery()
	crowd := c.m_sample.getCrowd()
	filter := crowd.GetFilter(0)
	halfExtents := crowd.GetQueryExtents()

	if adjust {
		vel := make([]float64, 3)
		// Request velocity
		if c.m_agentDebug.Idx != -1 {
			ag := crowd.GetAgent(int(c.m_agentDebug.Idx))
			if ag != nil && ag.Active {
				calcVel(vel, common.SliceTToSlice[float32, float64](ag.Npos[:]), p, float64(ag.Params.MaxSpeed))
				crowd.RequestMoveVelocity(int(c.m_agentDebug.Idx), common.SliceTToSlice[float64, float32](vel))
			}
		} else {
			for i := 0; i < crowd.GetAgentCount(); i++ {
				ag := crowd.GetAgent(i)
				if !ag.Active {
					continue
				}
				calcVel(vel, common.SliceTToSlice[float32, float64](ag.Npos[:]), p, float64(ag.Params.MaxSpeed))
				crowd.RequestMoveVelocity(i, common.SliceTToSlice[float64, float32](vel))
			}
		}
	} else {
		c.m_targetRef, _ = navquery.FindNearestPoly(common.SliceTToSlice[float64, float32](p), halfExtents, filter, common.SliceTToSlice[float64, float32](c.m_targetPos))

		if c.m_agentDebug.Idx != -1 {
			ag := crowd.GetAgent(int(c.m_agentDebug.Idx))
			if ag != nil && ag.Active {
				crowd.RequestMoveTarget(int(c.m_agentDebug.Idx), c.m_targetRef, common.SliceTToSlice[float64, float32](c.m_targetPos))
			}

		} else {
			for i := 0; i < crowd.GetAgentCount(); i++ {
				ag := crowd.GetAgent(i)
				if !ag.Active {
					continue
				}
				crowd.RequestMoveTarget(i, c.m_targetRef, common.SliceTToSlice[float64, float32](c.m_targetPos))
			}
		}
	}
}

func (c *CrowdToolState) hitTestAgents(s, p []float64) int {
	if c.m_sample == nil {
		return -1
	}
	crowd := c.m_sample.getCrowd()

	isel := -1
	tsel := math.MaxFloat64

	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}
		bmin := make([]float64, 3)
		bmax := make([]float64, 3)
		getAgentBounds(ag, bmin, bmax)
		var tmin, tmax float64
		if isectSegAABB(s, p, bmin, bmax, tmin, tmax) {
			if tmin > 0 && tmin < tsel {
				isel = i
				tsel = tmin
			}
		}
	}

	return isel
}
func (c *CrowdToolState) updateAgentParams() {
	if c.m_sample == nil {
		return
	}
	crowd := c.m_sample.getCrowd()
	if crowd == nil {
		return
	}

	obstacleAvoidanceType := 0

	obstacleAvoidanceType = int(c.cfg.ToolsConfig.ObstacleAvoidanceType)

	var params detour_crowd.DtCrowdAgentParams

	for i := 0; i < crowd.GetAgentCount(); i++ {
		ag := crowd.GetAgent(i)
		if !ag.Active {
			continue
		}
		params = *ag.Params
		params.UpdateFlags = c.cfg.ToolsConfig.GetExpandOptionsUpdateFlags()
		params.ObstacleAvoidanceType = obstacleAvoidanceType
		params.SeparationWeight = float32(c.cfg.ToolsConfig.SeparationWeight)
		crowd.UpdateAgentParameters(i, &params)
	}
}

type CrowdTool struct {
	m_sample *Sample
	m_state  *CrowdToolState

	ctx *Content
}

func newCrowdTool(ctx *Content) *CrowdTool {
	c := &CrowdTool{ctx: ctx}
	ctx.GetConfig().ToolsConfig.ExpandOptionsOnchange = c.m_state.updateAgentParams
	return c
}

func (c *CrowdTool) Type() int {
	return TOOL_CROWD
}

func (c *CrowdTool) init(sample *Sample) {
	if sample == nil {
		return
	}
	c.m_sample = sample
	c.m_state = sample.getToolState(c.Type()).(*CrowdToolState)
	if c.m_state == nil {
		c.m_state = newCrowdToolState(c.ctx)
		sample.setToolState(c.Type(), c.m_state)
	}
	c.m_state.init(sample)
}

func (c *CrowdTool) handleClick(s []float64, p []float64, shift bool) {
	if c.m_sample == nil {
		return
	}
	if c.m_state == nil {
		return
	}
	geom := c.m_sample.getInputGeom()
	if geom == nil {
		return
	}
	crowd := c.m_sample.getCrowd()
	if crowd == nil {
		return
	}

	if c.ctx.GetConfig().ToolsConfig.ToolModel == config.DescCrowdTool_TOOLMODE_CREATE {
		if shift {
			// Delete
			ahit := c.m_state.hitTestAgents(s, p)
			if ahit != -1 {
				c.m_state.removeAgent(ahit)
			}

		} else {
			// Add
			c.m_state.addAgent(p)
		}
	} else if c.ctx.GetConfig().ToolsConfig.ToolModel == config.DescCrowdTool_TOOLMODE_MOVE_TARGET {
		c.m_state.setMoveTarget(p, shift)
	} else if c.ctx.GetConfig().ToolsConfig.ToolModel == config.DescCrowdTool_TOOLMODE_SELECT {
		// Highlight
		ahit := c.m_state.hitTestAgents(s, p)
		c.m_state.hilightAgent(ahit)
	} else if c.ctx.GetConfig().ToolsConfig.ToolModel == config.DescCrowdTool_TOOLMODE_TOGGLE_POLYS {
		nav := c.m_sample.getNavMesh()
		navquery := c.m_sample.getNavMeshQuery()
		if nav != nil && navquery != nil {
			filter := detour.DTQueryFilter{}
			halfExtents := crowd.GetQueryExtents()
			tgt := make([]float64, 3)
			ref, _ := navquery.FindNearestPoly(common.SliceTToSlice[float64, float32](p), halfExtents, &filter, common.SliceTToSlice[float64, float32](tgt))
			if ref != 0 {
				flags, status := nav.GetPolyFlags(ref)
				if status.DtStatusSucceed() {
					flags ^= config.SAMPLE_POLYFLAGS_DISABLED
					nav.SetPolyFlags(ref, flags)
				}
			}
		}
	}
}

func (c *CrowdTool) handleStep() {
	if c.m_state == nil {
		return
	}

	dt := 1.0 / 20.0
	c.m_state.updateTick(dt)
}
