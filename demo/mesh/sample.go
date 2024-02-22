package mesh

import (
	"fmt"
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/gonavmesh/common/rw"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/detour"
	"github.com/gorustyt/gonavmesh/detour_crowd"
	"github.com/gorustyt/gonavmesh/recast"
	"io"
	"log/slog"
)

const NAVMESHSET_MAGIC = 'M'<<24 | 'S'<<16 | 'E'<<8 | 'T' //'MSET';
const NAVMESHSET_VERSION = 1

/// These are just sample areas to use consistent values across the samples.
/// The use should specify these base on his needs.

type Sample struct {
	m_geom     *InputGeom
	m_navMesh  detour.IDtNavMesh
	m_navQuery detour.NavMeshQuery
	m_crowd    *detour_crowd.DtCrowd

	m_navMeshDrawFlags int

	m_cellSize             float32
	m_cellHeight           float32
	m_agentHeight          float32
	m_agentRadius          float32
	m_agentMaxClimb        float32
	m_agentMaxSlope        float32
	m_regionMinSize        float32
	m_regionMergeSize      float32
	m_edgeMaxLen           float32
	m_edgeMaxError         float32
	m_vertsPerPoly         float32
	m_detailSampleDist     float32
	m_detailSampleMaxError float32
	m_partitionType        int

	m_filterLowHangingObstacles    bool
	m_filterLedgeSpans             bool
	m_filterWalkableLowHeightSpans bool

	m_tool       SampleTool
	m_toolStates [MAX_TOOLS]SampleToolState
	m_dd         debug_utils.DuDebugDraw

	ctx *Content
}

func NewSample(ctx *Content) *Sample {
	return &Sample{ctx: ctx}
}

type NavMeshSetHeader struct {
	Magic    int
	Version  int
	NumTiles int
	Params   *detour.NavMeshParams
}

func NewNavMeshSetHeader() *NavMeshSetHeader {
	return &NavMeshSetHeader{
		Params: &detour.NavMeshParams{},
	}
}

func (n *NavMeshSetHeader) Encode(w *rw.ReaderWriter) {
	w.WriteInt32(n.Magic)
	w.WriteInt32(n.Version)
	w.WriteInt32(n.NumTiles)
	n.Params.ToBin(w)
}

func (n *NavMeshSetHeader) Decode(r *rw.ReaderWriter) {
	n.Magic = int(r.ReadInt32())
	n.Version = int(r.ReadInt32())
	n.NumTiles = int(r.ReadInt32())
	n.Params.FromBin(r)
}

type NavMeshTileHeader struct {
	TileRef  detour.DTTileRef
	DataSize int
}

func (n *NavMeshTileHeader) Decode(r *rw.ReaderWriter) {
	n.TileRef = detour.DTTileRef(r.ReadInt32())
	n.DataSize = int(r.ReadInt32())
}
func (n *NavMeshTileHeader) Encode(w *rw.ReaderWriter) {
	w.WriteInt32(int32(n.TileRef))
	w.WriteInt32(n.DataSize)
}

func (s *Sample) saveAll(writer fyne.URIWriteCloser, mesh detour.IDtNavMesh) {
	if mesh == nil {
		return
	}
	w := rw.NewNavMeshDataBinWriter()
	// Store header.
	header := NewNavMeshSetHeader()
	header.Magic = NAVMESHSET_MAGIC
	header.Version = NAVMESHSET_VERSION
	header.NumTiles = 0
	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile == nil || tile.Header == nil || tile.Data == nil {
			continue
		}
		header.NumTiles++
	}
	header.Params = mesh.GetParams()
	header.Encode(w)

	// Store tiles.
	for i := int32(0); i < mesh.GetMaxTiles(); i++ {
		tile := mesh.GetTile(int(i))
		if tile == nil || tile.Header == nil || tile.Data == nil {
			continue
		}

		var tileHeader NavMeshTileHeader
		tileHeader.TileRef = mesh.GetTileRef(tile)
		dataSize := tile.Data.ToBin(w)
		tileHeader.DataSize = dataSize
		tileHeader.Encode(w)
		_, err := writer.Write(w.GetWriteBytes())
		if err != nil {
			panic(err)
		}
	}
	writer.Close()
}

func (s *Sample) loadAll(reader fyne.URIReadCloser) detour.IDtNavMesh {
	data, err := io.ReadAll(reader)
	if err != nil {
		panic(err)
	}
	r := rw.NewNavMeshDataBinReader(data)
	// Read header.
	header := NewNavMeshSetHeader()
	header.Decode(r)
	if header.Magic != NAVMESHSET_MAGIC {
		slog.Error("header.Magic != NAVMESHSET_MAGIC")
		return nil
	}
	if header.Version != NAVMESHSET_VERSION {
		slog.Error("header.Version != NAVMESHSET_VERSION ")
		return nil
	}

	mesh, _ := detour.NewDtNavMeshWithParams(header.Params)

	// Read tiles.
	for i := 0; i < header.NumTiles; i++ {
		var tileHeader NavMeshTileHeader
		tileHeader.Decode(r)

		if tileHeader.TileRef == 0 || tileHeader.DataSize == 0 {
			break
		}

		meshData := &detour.NavMeshData{}
		err = meshData.FromBin(r)
		if err != nil {
			panic(err)
		}
		mesh.AddTile(meshData, detour.DT_TILE_FREE_DATA, tileHeader.TileRef)
	}
	reader.Close()
	return mesh
}

func (s *Sample) getInputGeom() *InputGeom             { return s.m_geom }
func (s *Sample) getNavMesh() detour.IDtNavMesh        { return s.m_navMesh }
func (s *Sample) getNavMeshQuery() detour.NavMeshQuery { return s.m_navQuery }
func (s *Sample) getCrowd() *detour_crowd.DtCrowd      { return s.m_crowd }
func (s *Sample) getAgentRadius() float32              { return s.m_agentRadius }
func (s *Sample) getAgentHeight() float32              { return s.m_agentHeight }
func (s *Sample) getAgentClimb() float32               { return s.m_agentMaxClimb }

func (s *Sample) getNavMeshDrawFlags() int                  { return s.m_navMeshDrawFlags }
func (s *Sample) setNavMeshDrawFlags(flags int)             { s.m_navMeshDrawFlags = flags }
func (s *Sample) getToolState(Type int) SampleToolState     { return s.m_toolStates[Type] }
func (s *Sample) setToolState(Type int, ss SampleToolState) { s.m_toolStates[Type] = ss }

func (s *Sample) getDebugDraw() debug_utils.DuDebugDraw { return s.m_dd }

func (s *Sample) handleSettings() {
}
func (s *Sample) handleRenderOverlay(proj, model []float32, view []int) {
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

	// Draw rcMeshLoaderObj
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

func (s *Sample) handleDebugMode() {
}
func (s *Sample) handleClick(ss, p []float32, shift bool) {
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

func (s *Sample) handleUpdate(dt float32) {
	if s.m_tool != nil {
		s.m_tool.handleUpdate(dt)
	}

	s.updateToolStates(dt)
}

func (s *Sample) updateToolStates(dt float32) {
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
func (s *Sample) renderOverlayToolStates(proj, model []float32, view []int) {
	for i := 0; i < MAX_TOOLS; i++ {
		if s.m_toolStates[i] != nil {
			s.m_toolStates[i].handleRenderOverlay(proj, model, view)
		}

	}
}

func (s *Sample) handleCommonSettings() {

	if s.m_geom != nil {
		bmin := s.m_geom.getNavMeshBoundsMin()
		bmax := s.m_geom.getNavMeshBoundsMax()
		gw := 0
		gh := 0
		recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
		text := fmt.Sprintf("Voxels  %d x %d", gw, gh)

	}

}
