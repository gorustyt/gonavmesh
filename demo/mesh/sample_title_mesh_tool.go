package mesh

import (
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/demo/config"
)

type meshTitleTool struct {
	m_sample    *Sample_TileMesh
	m_hitPos    [3]float32
	m_hitPosSet bool
	cfg         *config.Config
}

func newMeshTitleTool(ctx *Content) *meshTitleTool {
	t := &meshTitleTool{cfg: ctx.GetConfig()}
	t.cfg.ToolsConfig.OnCreateTilesRemoveAllClick = t.m_sample.removeAllTiles
	t.cfg.ToolsConfig.OnCreateTilesCreateAllClick = t.m_sample.buildAllTiles
	return t
}

func (m *meshTitleTool) init(sample *Sample) {
	m.m_sample = newSampleTileMesh(sample)
}

func (m *meshTitleTool) Type() int { return TOOL_TILE_EDIT }

func (m *meshTitleTool) handleClick(s []float32, p []float32, shift bool) {
	m.m_hitPosSet = true
	copy(m.m_hitPos[:], p)
	if m.m_sample != nil {
		if shift {
			m.m_sample.removeTile(m.m_hitPos[:])
		} else {
			m.m_sample.buildTile(m.m_hitPos[:])
		}

	}
}

func (m *meshTitleTool) handleToggle() {}

func (m *meshTitleTool) handleStep() {}

func (m *meshTitleTool) handleUpdate(dt float32) {}

func (m *meshTitleTool) handleRender() {
	if m.m_hitPosSet {
		s := m.m_sample.getAgentRadius()
		glColor4ub(0, 0, 0, 128)
		glLineWidth(2.0)
		glBegin(GL_LINES)
		glVertex3f(m.m_hitPos[0]-s, m.m_hitPos[1]+0.1, m.m_hitPos[2])
		glVertex3f(m.m_hitPos[0]+s, m.m_hitPos[1]+0.1, m.m_hitPos[2])
		glVertex3f(m.m_hitPos[0], m.m_hitPos[1]-s+0.1, m.m_hitPos[2])
		glVertex3f(m.m_hitPos[0], m.m_hitPos[1]+s+0.1, m.m_hitPos[2])
		glVertex3f(m.m_hitPos[0], m.m_hitPos[1]+0.1, m.m_hitPos[2]-s)
		glVertex3f(m.m_hitPos[0], m.m_hitPos[1]+0.1, m.m_hitPos[2]+s)
		gl.End()
		gl.LineWidth(1.0)
	}
}

func (m *meshTitleTool) handleRenderOverlay(proj, model []float32, view []int) {
	res := common.GluProject([]float32{m.m_hitPos[0], m.m_hitPos[1], m.m_hitPos[2]},
		model, proj, view)
	if m.m_hitPosSet && len(res) > 0 {
		tx := 0
		ty := 0
		x, y := int(res[0]), int(res[1])
		m.m_sample.getTilePos(m.m_hitPos[:], &tx, &ty)
		text := fmt.Sprintf("(%d,%d)", tx, ty)
		m.gs.imguiDrawText(x, y-25, IMGUI_ALIGN_CENTER, text, imguiRGBA(0, 0, 0, 220))
	}

	//Tool help
	h := view[3]
	m.gs.imguiDrawText(280, h-40, IMGUI_ALIGN_LEFT, "LMB: Rebuild hit tile.  Shift+LMB: Clear hit tile.", imguiRGBA(255, 255, 255, 192))
}
