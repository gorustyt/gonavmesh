package gui

import (
	"encoding/binary"
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
	m_navQuery *recast.DtNavMeshQuery
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

	m_tool       *SampleTool
	m_toolStates [MAX_TOOLS]*SampleToolState
	m_ctx        *logger
	m_dd         DebugDrawGL
}

func (s *Sample) getInputGeom() *InputGeom             { return s.m_geom }
func (s *Sample) getNavMesh() recast.IDtNavMesh        { return s.m_navMesh }
func (s *Sample) getNavMeshQuery() recast.NavMeshQuery { return s.m_navQuery }
func (s *Sample) getCrowd() *recast.DtCrowd            { return s.m_crowd }
func (s *Sample) getAgentRadius() float64              { return s.m_agentRadius }
func (s *Sample) getAgentHeight() float64              { return s.m_agentHeight }
func (s *Sample) getAgentClimb() float64               { return s.m_agentMaxClimb }

func (s *Sample) getNavMeshDrawFlags() int                   { return s.m_navMeshDrawFlags }
func (s *Sample) setNavMeshDrawFlags(flags int)              { s.m_navMeshDrawFlags = flags }
func (s *Sample) getToolState(Type int) *SampleToolState     { return s.m_toolStates[Type] }
func (s *Sample) setToolState(Type int, ss *SampleToolState) { s.m_toolStates[Type] = ss }

func (s *Sample) getDebugDraw() DebugDrawGL { return s.m_dd }
func newSample() *Sample {
	return &Sample{}
}

func (s *Sample) handleSettings() {
}
func (s *Sample) handleTools() {
}

func (s *Sample) handleDebugMode() {
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

func (s *Sample) SaveAll(p string, mesh *recast.DtNavMesh) {

}
