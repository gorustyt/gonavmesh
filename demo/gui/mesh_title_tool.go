package gui

type meshTitleTool struct {
	m_sample    *Sample_TileMesh
	m_hitPos    [3]float64
	m_hitPosSet bool
	gs          *guiState
}

func newMeshTitleTool(gs *guiState) *meshTitleTool {
	t := &meshTitleTool{gs: gs}
	return t
}

func (m *meshTitleTool) init(sample *Sample) {
	m.m_sample = newSampleTileMesh(sample, m.gs, newLogger())
}
func (m *meshTitleTool) Type() int { return TOOL_TILE_EDIT }
func (m *meshTitleTool) reset()    {}

func (m *meshTitleTool) handleMenu() {
	m.gs.imguiLabel("Create Tiles")
	if m.gs.imguiButton("Create All") {
		if m.m_sample != nil {
			m.m_sample.buildAllTiles()
		}
	}
	if m.gs.imguiButton("Remove All") {
		if m.m_sample != nil {
			m.m_sample.removeAllTiles()
		}

	}
}

func (m *meshTitleTool) handleClick(s []float64, p []float64, shift bool) {
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

func (m *meshTitleTool) handleUpdate(dt float64) {}

func (m *meshTitleTool) handleRender() {
	if m.m_hitPosSet {
		//s := m.m_sample.getAgentRadius();
		//glColor4ub(0,0,0,128);
		//glLineWidth(2.0);
		//glBegin(GL_LINES);
		//glVertex3f(m_hitPos[0]-s,m_hitPos[1]+0.1f,m_hitPos[2]);
		//glVertex3f(m_hitPos[0]+s,m_hitPos[1]+0.1f,m_hitPos[2]);
		//glVertex3f(m_hitPos[0],m_hitPos[1]-s+0.1f,m_hitPos[2]);
		//glVertex3f(m_hitPos[0],m_hitPos[1]+s+0.1f,m_hitPos[2]);
		//glVertex3f(m_hitPos[0],m_hitPos[1]+0.1f,m_hitPos[2]-s);
		//glVertex3f(m_hitPos[0],m_hitPos[1]+0.1f,m_hitPos[2]+s);
		//glEnd();
		//glLineWidth(1.0f);
	}
}

func (m *meshTitleTool) handleRenderOverlay(proj, model *float64, view []int) {
	//var  x, y, z float64
	//if (m_hitPosSet && gluProject((GLdouble)m_hitPos[0], (GLdouble)m_hitPos[1], (GLdouble)m_hitPos[2],
	//model, proj, view, &x, &y, &z)){
	//int tx=0, ty=0;
	//m_sample->getTilePos(m_hitPos, tx, ty);
	//char text[32];
	//snprintf(text,32,"(%d,%d)", tx,ty);
	//imguiDrawText((int)x, (int)y-25, IMGUI_ALIGN_CENTER, text, imguiRGBA(0,0,0,220));
	//}

	// Tool help
	//h := view[3];
	//imguiDrawText(280, h-40, IMGUI_ALIGN_LEFT, "LMB: Rebuild hit tile.  Shift+LMB: Clear hit tile.", imguiRGBA(255,255,255,192));
}
