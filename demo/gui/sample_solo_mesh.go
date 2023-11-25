package gui

import (
	"fmt"
	"gonavamesh/recast"
)

type SoloMeshDrawMode int

const (
	SOLOMESH_DRAWMODE_NAVMESH = iota
	SOLOMESH_DRAWMODE_NAVMESH_TRANS
	SOLOMESH_DRAWMODE_NAVMESH_BVTREE
	SOLOMESH_DRAWMODE_NAVMESH_NODES
	SOLOMESH_DRAWMODE_NAVMESH_INVIS
	SOLOMESH_DRAWMODE_MESH
	SOLOMESH_DRAWMODE_VOXELS
	SOLOMESH_DRAWMODE_VOXELS_WALKABLE
	SOLOMESH_DRAWMODE_COMPACT
	SOLOMESH_DRAWMODE_COMPACT_DISTANCE
	SOLOMESH_DRAWMODE_COMPACT_REGIONS
	SOLOMESH_DRAWMODE_REGION_CONNECTIONS
	SOLOMESH_DRAWMODE_RAW_CONTOURS
	SOLOMESH_DRAWMODE_BOTH_CONTOURS
	SOLOMESH_DRAWMODE_CONTOURS
	SOLOMESH_DRAWMODE_POLYMESH
	SOLOMESH_DRAWMODE_POLYMESH_DETAIL
	SOLOMESH_MAX_DRAWMODE
)

type SampleSoloMesh struct {
	*Sample
	m_keepInterResults bool
	m_totalBuildTimeMs float64

	m_triareas []int
	m_solid    *recast.RcHeightfield
	m_chf      *recast.RcCompactHeightfield
	m_cset     *recast.RcContourSet
	m_pmesh    *recast.RcPolyMesh
	m_cfg      *recast.RcConfig
	m_dmesh    *recast.RcPolyMeshDetail
	m_drawMode SoloMeshDrawMode
}

func newSampleSoloMesh(gs *guiState) *SampleSoloMesh {
	s := &SampleSoloMesh{
		m_drawMode:         DRAWMODE_NAVMESH,
		m_keepInterResults: true,
	}
	s.setTool(newMeshTitleTool(gs))
	return s
}

func (s *SampleSoloMesh) cleanup() {

	s.m_triareas = nil
	s.m_solid = nil
	s.m_chf = nil
	s.m_cset = nil
	s.m_pmesh = nil
	s.m_dmesh = nil
	s.m_navMesh = nil
}

func (s *SampleSoloMesh) handleSettings() {
	s.Sample.handleCommonSettings()

	if s.gs.imguiCheck("Keep Itermediate Results", s.m_keepInterResults) {
		s.m_keepInterResults = !s.m_keepInterResults
	}

	s.gs.imguiSeparator()

	s.gs.imguiIndent()
	s.gs.imguiIndent()

	if s.gs.imguiButton("Save") {
		s.Sample.SaveAll("solo_navmesh.bin", s.m_navMesh)
	}

	if s.gs.imguiButton("Load") {
		s.m_navMesh = s.Sample.LoadAll("solo_navmesh.bin")
		s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048)
	}

	s.gs.imguiUnindent()
	s.gs.imguiUnindent()
	msg := fmt.Sprintf("Build Time: %.1fms", s.m_totalBuildTimeMs)
	s.gs.imguiLabel(msg)

	s.gs.imguiSeparator()
}

func (s *SampleSoloMesh) handleTools() {
	tType := s.m_tool.Type()
	if s.m_tool == nil {
		tType = TOOL_NONE
	}

	if s.gs.imguiCheck("Test Navmesh", tType == TOOL_NAVMESH_TESTER) {
		s.setTool(newNavMeshTesterTool(s.gs))
	}
	if s.gs.imguiCheck("Prune Navmesh", tType == TOOL_NAVMESH_PRUNE) {
		s.setTool(newNavMeshPruneTool(s.gs))
	}
	if s.gs.imguiCheck("Create Off-Mesh Connections", tType == TOOL_OFFMESH_CONNECTION) {
		s.setTool(newOffMeshConnectionTool(s.gs))
	}
	if s.gs.imguiCheck("Create Convex Volumes", tType == TOOL_CONVEX_VOLUME) {
		s.setTool(newConvexVolumeTool(s.gs))
	}
	if s.gs.imguiCheck("Create Crowds", tType == TOOL_CROWD) {
		s.setTool(newCrowdTool())
	}

	s.gs.imguiSeparatorLine()

	s.gs.imguiIndent()

	if s.m_tool != nil {
		s.m_tool.handleMenu()
	}

	s.gs.imguiUnindent()

}
