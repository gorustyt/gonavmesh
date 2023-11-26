package gui

import (
	"encoding/binary"
	"fmt"
	"gonavamesh/debug_utils"
	"gonavamesh/recast"
	"io"
	"os"
)

/// Tool types.

const (
	TOOL_NONE = iota
	TOOL_TILE_EDIT
	TOOL_TILE_HIGHLIGHT
	TOOL_TEMP_OBSTACLE
	TOOL_NAVMESH_TESTER
	TOOL_NAVMESH_PRUNE
	TOOL_OFFMESH_CONNECTION
	TOOL_CONVEX_VOLUME
	TOOL_CROWD
	MAX_TOOLS
)

/// These are just sample areas to use consistent values across the samples.
/// The use should specify these base on his needs.

const (
	SAMPLE_POLYAREA_GROUND = iota
	SAMPLE_POLYAREA_WATER
	SAMPLE_POLYAREA_ROAD
	SAMPLE_POLYAREA_DOOR
	SAMPLE_POLYAREA_GRASS
	SAMPLE_POLYAREA_JUMP
)

const (
	SAMPLE_POLYFLAGS_WALK     = 0x01   // Ability to walk (ground, grass, road)
	SAMPLE_POLYFLAGS_SWIM     = 0x02   // Ability to swim (water).
	SAMPLE_POLYFLAGS_DOOR     = 0x04   // Ability to move through doors.
	SAMPLE_POLYFLAGS_JUMP     = 0x08   // Ability to jump.
	SAMPLE_POLYFLAGS_DISABLED = 0x10   // Disabled polygon
	SAMPLE_POLYFLAGS_ALL      = 0xffff // All abilities.
)
const NAVMESHSET_MAGIC = 'M'<<24 | 'S'<<16 | 'E'<<8 | 'T' //'MSET';
const NAVMESHSET_VERSION = 1

const (
	SAMPLE_PARTITION_WATERSHED = iota
	SAMPLE_PARTITION_MONOTONE
	SAMPLE_PARTITION_LAYERS
)

type NavMeshSetHeader struct {
	magic    int
	version  int
	numTiles int
	params   recast.NavMeshParams
}

type NavMeshTileHeader struct {
	tileRef  recast.DtTileRef
	dataSize int
}

type Sample struct {
	m_geom     *InputGeom
	m_navMesh  recast.IDtNavMesh
	m_navQuery recast.NavMeshQuery
	m_crowd    *recast.DtCrowd

	m_navMeshDrawFlags int

	m_cellSize             float64
	m_cellHeight           float64
	m_agentHeight          float64
	m_agentRadius          float64
	m_agentMaxClimb        float64
	m_agentMaxSlope        float64
	m_regionMinSize        float64
	m_regionMergeSize      float64
	m_edgeMaxLen           float64
	m_edgeMaxError         float64
	m_vertsPerPoly         float64
	m_detailSampleDist     float64
	m_detailSampleMaxError float64
	m_partitionType        int

	m_filterLowHangingObstacles    bool
	m_filterLedgeSpans             bool
	m_filterWalkableLowHeightSpans bool

	m_tool       SampleTool
	m_toolStates [MAX_TOOLS]SampleToolState
	m_ctx        *logger
	m_dd         debug_utils.DuDebugDraw

	gs *guiState
}

func (s *Sample) getInputGeom() *InputGeom             { return s.m_geom }
func (s *Sample) getNavMesh() recast.IDtNavMesh        { return s.m_navMesh }
func (s *Sample) getNavMeshQuery() recast.NavMeshQuery { return s.m_navQuery }
func (s *Sample) getCrowd() *recast.DtCrowd            { return s.m_crowd }
func (s *Sample) getAgentRadius() float64              { return s.m_agentRadius }
func (s *Sample) getAgentHeight() float64              { return s.m_agentHeight }
func (s *Sample) getAgentClimb() float64               { return s.m_agentMaxClimb }

func (s *Sample) getNavMeshDrawFlags() int                  { return s.m_navMeshDrawFlags }
func (s *Sample) setNavMeshDrawFlags(flags int)             { s.m_navMeshDrawFlags = flags }
func (s *Sample) getToolState(Type int) SampleToolState     { return s.m_toolStates[Type] }
func (s *Sample) setToolState(Type int, ss SampleToolState) { s.m_toolStates[Type] = ss }

func (s *Sample) getDebugDraw() debug_utils.DuDebugDraw { return s.m_dd }
func newSample(gs *guiState) *Sample {
	return &Sample{
		gs: gs,
	}
}

func (s *Sample) handleSettings() {
}
func (s *Sample) handleRenderOverlay(proj, model float64, view int) {
}
func (s *Sample) setTool(tool SampleTool) {
	s.m_tool = tool
	if tool != nil {
		s.m_tool.init(s)
	}

}

func (s *Sample) handleRender() {
	if s.m_geom == nil {
		return
	}

	// Draw mesh
	debug_utils.DuDebugDrawTriMesh(s.m_dd, s.m_geom.getMesh().getVerts(), s.m_geom.getMesh().getVertCount(),
		s.m_geom.getMesh().getTris(), s.m_geom.getMesh().getNormals(), s.m_geom.getMesh().getTriCount(), []int{}, 1.0)
	// Draw bounds
	bmin := s.m_geom.getMeshBoundsMin()
	bmax := s.m_geom.getMeshBoundsMax()
	debug_utils.DuDebugDrawBoxWire(s.m_dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], debug_utils.DuRGBA(255, 255, 255, 128), 1.0)
}

func (s *Sample) handleMeshChanged(geom *InputGeom) {
	s.m_geom = geom

	buildSettings := geom.getBuildSettings()
	if buildSettings != nil {
		s.m_cellSize = buildSettings.cellSize
		s.m_cellHeight = buildSettings.cellHeight
		s.m_agentHeight = buildSettings.agentHeight
		s.m_agentRadius = buildSettings.agentRadius
		s.m_agentMaxClimb = buildSettings.agentMaxClimb
		s.m_agentMaxSlope = buildSettings.agentMaxSlope
		s.m_regionMinSize = buildSettings.regionMinSize
		s.m_regionMergeSize = buildSettings.regionMergeSize
		s.m_edgeMaxLen = buildSettings.edgeMaxLen
		s.m_edgeMaxError = buildSettings.edgeMaxError
		s.m_vertsPerPoly = buildSettings.vertsPerPoly
		s.m_detailSampleDist = buildSettings.detailSampleDist
		s.m_detailSampleMaxError = buildSettings.detailSampleMaxError
		s.m_partitionType = buildSettings.partitionType
	}
}
func (s *Sample) collectSettings(settings *BuildSettings) {
	settings.cellSize = s.m_cellSize
	settings.cellHeight = s.m_cellHeight
	settings.agentHeight = s.m_agentHeight
	settings.agentRadius = s.m_agentRadius
	settings.agentMaxClimb = s.m_agentMaxClimb
	settings.agentMaxSlope = s.m_agentMaxSlope
	settings.regionMinSize = s.m_regionMinSize
	settings.regionMergeSize = s.m_regionMergeSize
	settings.edgeMaxLen = s.m_edgeMaxLen
	settings.edgeMaxError = s.m_edgeMaxError
	settings.vertsPerPoly = s.m_vertsPerPoly
	settings.detailSampleDist = s.m_detailSampleDist
	settings.detailSampleMaxError = s.m_detailSampleMaxError
	settings.partitionType = s.m_partitionType
}

func (s *Sample) resetCommonSettings() {
	s.m_cellSize = 0.3
	s.m_cellHeight = 0.2
	s.m_agentHeight = 2.0
	s.m_agentRadius = 0.6
	s.m_agentMaxClimb = 0.9
	s.m_agentMaxSlope = 45.0
	s.m_regionMinSize = 8
	s.m_regionMergeSize = 20
	s.m_edgeMaxLen = 12.0
	s.m_edgeMaxError = 1.3
	s.m_vertsPerPoly = 6.0
	s.m_detailSampleDist = 6.0
	s.m_detailSampleMaxError = 1.0
	s.m_partitionType = SAMPLE_PARTITION_WATERSHED
}
func (s *Sample) handleTools() {
}

func (s *Sample) handleDebugMode() {
}
func (s *Sample) handleClick(ss, p []float64, shift bool) {
	if s.m_tool != nil {
		s.m_tool.handleClick(ss, p, shift)
	}
}

func (s *Sample) handleToggle() {
	if s.m_tool != nil {
		s.m_tool.handleToggle()
	}

}

func (s *Sample) handleStep() {
	if s.m_tool != nil {
		s.m_tool.handleStep()
	}

}

func (s *Sample) handleBuild() bool {
	return true
}

func (s *Sample) handleUpdate(dt float64) {
	if s.m_tool != nil {
		s.m_tool.handleUpdate(dt)
	}

	s.updateToolStates(dt)
}

func (s *Sample) updateToolStates(dt float64) {
	for i := 0; i < MAX_TOOLS; i++ {
		if s.m_toolStates[i] != nil {
			s.m_toolStates[i].handleUpdate(dt)
		}

	}
}

func (s *Sample) initToolStates(sample *Sample) {
	for i := 0; i < MAX_TOOLS; i++ {
		if s.m_toolStates[i] != nil {
			s.m_toolStates[i].init(sample)
		}

	}
}

func (s *Sample) resetToolStates() {
	for i := 0; i < MAX_TOOLS; i++ {
		if s.m_toolStates[i] != nil {
			s.m_toolStates[i].reset()
		}

	}
}

func (s *Sample) renderToolStates() {
	for i := 0; i < MAX_TOOLS; i++ {
		if s.m_toolStates[i] != nil {
			s.m_toolStates[i].handleRender()
		}

	}
}
func (s *Sample) renderOverlayToolStates(proj, model []float64, view []int) {
	for i := 0; i < MAX_TOOLS; i++ {
		if s.m_toolStates[i] != nil {
			s.m_toolStates[i].handleRenderOverlay(proj, model, view)
		}

	}
}
func (s *Sample) handleCommonSettings() {
	s.gs.imguiLabel("Rasterization")
	s.gs.imguiSlider("Cell Size", &s.m_cellSize, 0.1, 1.0, 0.01)
	s.gs.imguiSlider("Cell Height", &s.m_cellHeight, 0.1, 1.0, 0.01)

	if s.m_geom != nil {
		bmin := s.m_geom.getNavMeshBoundsMin()
		bmax := s.m_geom.getNavMeshBoundsMax()
		gw := 0
		gh := 0
		recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
		text := fmt.Sprintf("Voxels  %d x %d", gw, gh)
		s.gs.imguiValue(text)
	}

	s.gs.imguiSeparator()
	s.gs.imguiLabel("Agent")
	s.gs.imguiSlider("Height", &s.m_agentHeight, 0.1, 5.0, 0.1)
	s.gs.imguiSlider("Radius", &s.m_agentRadius, 0.0, 5.0, 0.1)
	s.gs.imguiSlider("Max Climb", &s.m_agentMaxClimb, 0.1, 5.0, 0.1)
	s.gs.imguiSlider("Max Slope", &s.m_agentMaxSlope, 0.0, 90.0, 1.0)

	s.gs.imguiSeparator()
	s.gs.imguiLabel("Region")
	s.gs.imguiSlider("Min Region Size", &s.m_regionMinSize, 0.0, 150.0, 1.0)
	s.gs.imguiSlider("Merged Region Size", &s.m_regionMergeSize, 0.0, 150.0, 1.0)

	s.gs.imguiSeparator()
	s.gs.imguiLabel("Partitioning")
	if s.gs.imguiCheck("Watershed", s.m_partitionType == SAMPLE_PARTITION_WATERSHED) {
		s.m_partitionType = SAMPLE_PARTITION_WATERSHED
	}

	if s.gs.imguiCheck("Monotone", s.m_partitionType == SAMPLE_PARTITION_MONOTONE) {
		s.m_partitionType = SAMPLE_PARTITION_MONOTONE
	}

	if s.gs.imguiCheck("Layers", s.m_partitionType == SAMPLE_PARTITION_LAYERS) {
		s.m_partitionType = SAMPLE_PARTITION_LAYERS
	}

	s.gs.imguiSeparator()
	s.gs.imguiLabel("Filtering")
	if s.gs.imguiCheck("Low Hanging Obstacles", s.m_filterLowHangingObstacles) {
		s.m_filterLowHangingObstacles = !s.m_filterLowHangingObstacles
	}

	if s.gs.imguiCheck("Ledge Spans", s.m_filterLedgeSpans) {
		s.m_filterLedgeSpans = !s.m_filterLedgeSpans
	}

	if s.gs.imguiCheck("Walkable Low Height Spans", s.m_filterWalkableLowHeightSpans) {
		s.m_filterWalkableLowHeightSpans = !s.m_filterWalkableLowHeightSpans
	}

	s.gs.imguiSeparator()
	s.gs.imguiLabel("Polygonization")
	s.gs.imguiSlider("Max Edge Length", &s.m_edgeMaxLen, 0.0, 50.0, 1.0)
	s.gs.imguiSlider("Max Edge Error", &s.m_edgeMaxError, 0.1, 3.0, 0.1)
	s.gs.imguiSlider("Verts Per Poly", &s.m_vertsPerPoly, 3.0, 12.0, 1.0)

	s.gs.imguiSeparator()
	s.gs.imguiLabel("Detail Mesh")
	s.gs.imguiSlider("Sample Distance", &s.m_detailSampleDist, 0.0, 16.0, 1.0)
	s.gs.imguiSlider("Max Sample Error", &s.m_detailSampleMaxError, 0.0, 16.0, 1.0)

	s.gs.imguiSeparator()
}
func (s *Sample) LoadAll(p string) *recast.DtNavMesh {
	f, err := os.Open(p)
	defer f.Close()
	if err != nil {
		panic(err)
	}
	readN := func(n int) ([]byte, bool) {
		buf := make([]byte, n)
		n, err := io.ReadFull(f, buf)
		if err == io.EOF {
			return []byte{}, false
		}
		if err != nil {
			panic(err)
		}
		return buf, true
	}
	header := &NavMeshSetHeader{}
	b, ok := readN(4)
	if !ok {
		return nil
	}
	header.magic = int(binary.LittleEndian.Uint32(b))
	if header.magic != NAVMESHSET_MAGIC {
		panic("")
	}
	return nil
}

func (s *Sample) SaveAll(p string, mesh recast.IDtNavMesh) {

}
