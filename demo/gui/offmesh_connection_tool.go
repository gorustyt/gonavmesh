package gui

import (
	"gonavamesh/common"
	"gonavamesh/debug_utils"
	"math"
)

type OffMeshConnectionTool struct {
	m_sample    *Sample
	m_hitPos    []float64
	m_hitPosSet bool
	m_bidir     bool
	m_oldFlags  int
	gs          *guiState
}

func newOffMeshConnectionTool(gs *guiState) *OffMeshConnectionTool {
	return &OffMeshConnectionTool{
		m_hitPos: make([]float64, 3),
		gs:       gs,
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
func (t *OffMeshConnectionTool) handleMenu() {
	if t.gs.imguiCheck("One Way", !t.m_bidir) {
		t.m_bidir = false
	}

	if t.gs.imguiCheck("Bidirectional", t.m_bidir) {
		t.m_bidir = true
	}

}
func (t *OffMeshConnectionTool) handleClick(s []float64, p []float64, shift bool) {
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
		nearestDist := math.MaxFloat64
		nearestIndex := -1
		verts := geom.getOffMeshConnectionVerts()
		for i := 0; i < geom.getOffMeshConnectionCount()*2; i++ {
			v := common.GetVs3(verts, i)
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
			area := SAMPLE_POLYAREA_JUMP
			flags := SAMPLE_POLYFLAGS_JUMP
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
func (t *OffMeshConnectionTool) handleRenderOverlay(proj, model *float64, view []int) {
	//var x, y, z float64

	// Draw start and end point labels
	//if t.m_hitPosSet && mgl32.Project(t.m_hitPos[0], t.m_hitPos[1], t.m_hitPos[2], model, proj, view, &x, &y, &z) {
	//	t.gs.imguiDrawText(int(x), (int)(y-25), IMGUI_ALIGN_CENTER, "Start", imguiRGBA(0, 0, 0, 220))
	//}//TODO
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
func (t *OffMeshConnectionTool) handleUpdate(dt float64) {}
