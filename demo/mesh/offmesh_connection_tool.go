package mesh

import (
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/demo/config"
	"math"
)

type OffMeshConnectionTool struct {
	m_sample    *Sample
	m_hitPos    []float32
	m_hitPosSet bool
	m_bidir     bool
	m_oldFlags  int
	cfg         *config.Config
}

func newOffMeshConnectionTool(ctx *Content) *OffMeshConnectionTool {
	return &OffMeshConnectionTool{
		m_hitPos: make([]float32, 3),
		cfg:      ctx.GetConfig(),
	}
}

func (t *OffMeshConnectionTool) Type() int { return TOOL_OFFMESH_CONNECTION }
func (t *OffMeshConnectionTool) init(sample *Sample) {
	if t.m_sample != sample {
		t.m_sample = sample
		t.m_oldFlags = t.m_sample.getNavMeshDrawFlags()
		t.m_sample.setNavMeshDrawFlags(t.m_oldFlags & ^debug_utils.DU_DRAWNAVMESH_OFFMESHCONS)
	}
}
func (t *OffMeshConnectionTool) reset() { t.m_hitPosSet = false }

func (t *OffMeshConnectionTool) handleClick(s []float32, p []float32, shift bool) {
	if t.m_sample == nil {
		return
	}
	geom := t.m_sample.getInputGeom()
	if geom == nil {
		return
	}

	if shift {
		// Delete
		// Find nearest link end-point
		nearestDist := math.Maxfloat32
		nearestIndex := -1
		verts := geom.getOffMeshConnectionVerts()
		for i := 0; i < geom.getOffMeshConnectionCount()*2; i++ {
			v := common.GetVert3(verts, i)
			d := common.VdistSqr(p, v)
			if d < nearestDist {
				nearestDist = d
				nearestIndex = i / 2 // Each link has two vertices.
			}
		}
		// If end point close enough, delete it.
		if nearestIndex != -1 &&
			math.Sqrt(nearestDist) < t.m_sample.getAgentRadius() {
			geom.deleteOffMeshConnection(nearestIndex)
		}
	} else {
		// Create
		if !t.m_hitPosSet {
			copy(t.m_hitPos, p)
			t.m_hitPosSet = true
		} else {
			area := config.SAMPLE_POLYAREA_JUMP
			flags := config.SAMPLE_POLYFLAGS_JUMP
			tmp := 0
			if t.m_bidir {
				tmp = 1
			}
			geom.addOffMeshConnection(t.m_hitPos, p, t.m_sample.getAgentRadius(), tmp, area, flags)
			t.m_hitPosSet = false
		}
	}
}
func (t *OffMeshConnectionTool) handleRender() {
	dd := t.m_sample.getDebugDraw()
	s := t.m_sample.getAgentRadius()

	if t.m_hitPosSet {
		debug_utils.DuDebugDrawCross(dd, t.m_hitPos[0], t.m_hitPos[1]+0.1, t.m_hitPos[2], s, debug_utils.DuRGBA(0, 0, 0, 128), 2.0)
	}
	geom := t.m_sample.getInputGeom()
	if geom != nil {
		geom.drawOffMeshConnections(dd, true)
	}

}
func (t *OffMeshConnectionTool) handleRenderOverlay(proj, model []float32, view []int) {
	res := common.GluProject([]float32{t.m_hitPos[0], t.m_hitPos[1], t.m_hitPos[2]}, model, proj, view)
	x, y := res[0], res[1]
	//Draw start and end point labels
	if t.m_hitPosSet && len(res) > 0 {
		t.gs.imguiDrawText(int(x), (int)(y-25), IMGUI_ALIGN_CENTER, "Start", imguiRGBA(0, 0, 0, 220))
	}
	// Tool help
	h := view[3]
	if !t.m_hitPosSet {
		t.gs.imguiDrawText(280, h-40, IMGUI_ALIGN_LEFT, "LMB: Create new connection.  SHIFT+LMB: Delete existing connection, click close to start or end point.", imguiRGBA(255, 255, 255, 192))
	} else {
		t.gs.imguiDrawText(280, h-40, IMGUI_ALIGN_LEFT, "LMB: Set connection end point and finish.", imguiRGBA(255, 255, 255, 192))
	}
}
func (t *OffMeshConnectionTool) handleToggle()           {}
func (t *OffMeshConnectionTool) handleStep()             {}
func (t *OffMeshConnectionTool) handleUpdate(dt float32) {}
