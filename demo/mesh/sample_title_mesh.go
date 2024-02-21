package mesh

import (
	"fmt"
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/detour"
	"github.com/gorustyt/gonavmesh/recast"
	"log"
	"math"
	"time"
)

type TitleMeshDrawMode int

const (
	TitleMeshDRAWMODE_NAVMESH TitleMeshDrawMode = iota
	TitleMeshDRAWMODE_NAVMESH_TRANS
	TitleMeshDRAWMODE_NAVMESH_BVTREE
	TitleMeshDRAWMODE_NAVMESH_NODES
	TitleMeshDRAWMODE_NAVMESH_PORTALS
	TitleMeshDRAWMODE_NAVMESH_INVIS
	TitleMeshDRAWMODE_MESH
	TitleMeshDRAWMODE_VOXELS
	TitleMeshDRAWMODE_VOXELS_WALKABLE
	TitleMeshDRAWMODE_COMPACT
	TitleMeshDRAWMODE_COMPACT_DISTANCE
	TitleMeshDRAWMODE_COMPACT_REGIONS
	TitleMeshDRAWMODE_REGION_CONNECTIONS
	TitleMeshDRAWMODE_RAW_CONTOURS
	TitleMeshDRAWMODE_BOTH_CONTOURS
	TitleMeshDRAWMODE_CONTOURS
	TitleMeshDRAWMODE_POLYMESH
	TitleMeshDRAWMODE_POLYMESH_DETAIL
	TitleMeshMAX_DRAWMODE
)

type Sample_TileMesh struct {
	*Sample
	m_keepInterResults bool
	m_buildAll         bool
	m_totalBuildTimeMs time.Duration

	m_triareas        []int
	m_solid           *recast.RcHeightfield
	m_chf             *recast.RcCompactHeightfield
	m_cset            *recast.RcContourSet
	m_pmesh           *recast.RcPolyMesh
	m_dmesh           *recast.RcPolyMeshDetail
	m_cfg             *recast.RcConfig
	m_drawMode        TitleMeshDrawMode
	m_maxTiles        int
	m_maxPolysPerTile int
	m_tileSize        float64

	m_tileCol           int
	m_lastBuiltTileBmin [3]float64
	m_lastBuiltTileBmax [3]float64
	m_tileBuildTime     time.Duration
	m_tileMemUsage      float64
	m_tileTriCount      int
}

func newSampleTileMesh(c *Content) *Sample_TileMesh {
	t := &Sample_TileMesh{
		m_drawMode: TitleMeshDRAWMODE_NAVMESH,
		m_tileCol:  debug_utils.DuRGBA(0, 0, 0, 32),
		m_buildAll: true,
		m_tileSize: 32,
	}
	s.setTool(newMeshTitleTool(c))
	return t
}

func (s *Sample_TileMesh) handleSettings() {
	s.Sample.handleCommonSettings()

	if s.gs.imguiCheck("Keep Itermediate Results", s.m_keepInterResults) {
		s.m_keepInterResults = !s.m_keepInterResults
	}

	if s.gs.imguiCheck("Build All Tiles", s.m_buildAll) {
		s.m_buildAll = !s.m_buildAll
	}

	s.gs.imguiLabel("Tiling")
	s.gs.imguiSlider("TileSize", &s.m_tileSize, 16.0, 1024.0, 16.0)

	if s.m_geom != nil {
		gw := 0
		gh := 0
		bmin := s.m_geom.getNavMeshBoundsMin()
		bmax := s.m_geom.getNavMeshBoundsMax()
		recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
		ts := int(s.m_tileSize)
		tw := (gw + ts - 1) / ts
		th := (gh + ts - 1) / ts
		text := fmt.Sprintf("Tiles  %d x %d", tw, th)
		s.gs.imguiValue(text)

		// Max tiles and max polys affect how the tile IDs are caculated.
		// There are 22 bits available for identifying a tile and a polygon.
		tileBits := min(common.Ilog2(common.NextPow2(tw*th)), 14)
		if tileBits > 14 {
			tileBits = 14
		}
		polyBits := 22 - tileBits
		m_maxTiles := 1 << tileBits
		m_maxPolysPerTile := 1 << polyBits
		text = fmt.Sprintf("Max Tiles  %d", m_maxTiles)
		s.gs.imguiValue(text)
		text = fmt.Sprintf("Max Polys  %d", m_maxPolysPerTile)
		s.gs.imguiValue(text)
	} else {
		s.m_maxTiles = 0
		s.m_maxPolysPerTile = 0
	}

	s.gs.imguiSeparator()

	s.gs.imguiIndent()
	s.gs.imguiIndent()

	if s.gs.imguiButton("Save") {
		s.Sample.SaveAll("all_tiles_navmesh.bin", s.m_navMesh)
	}

	if s.gs.imguiButton("Load") {
		s.m_navMesh = s.Sample.LoadAll("all_tiles_navmesh.bin")
		s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048)
	}

	s.gs.imguiUnindent()
	s.gs.imguiUnindent()

	msg := fmt.Sprintf("Build Time: %.1fms", s.m_totalBuildTimeMs)
	s.gs.imguiLabel(msg)

	s.gs.imguiSeparator()

	s.gs.imguiSeparator()
}
func (s *Sample_TileMesh) handleTools() {
	tType := s.m_tool.Type()
	if s.m_tool != nil {
		tType = TOOL_NONE
	}
	if s.gs.imguiCheck("Test Navmesh", tType == TOOL_NAVMESH_TESTER) {
		s.setTool(newNavMeshTesterTool(s.gs))
	}
	if s.gs.imguiCheck("Prune Navmesh", tType == TOOL_NAVMESH_PRUNE) {
		s.setTool(newNavMeshPruneTool(s.gs))
	}
	if s.gs.imguiCheck("Create Tiles", tType == TOOL_TILE_EDIT) {
		s.setTool(newMeshTitleTool(s.gs))
	}
	if s.gs.imguiCheck("Create Off-Mesh Links", tType == TOOL_OFFMESH_CONNECTION) {
		s.setTool(newOffMeshConnectionTool(s.gs))
	}
	if s.gs.imguiCheck("Create Convex Volumes", tType == TOOL_CONVEX_VOLUME) {
		s.setTool(newConvexVolumeTool(s.gs))
	}
	if s.gs.imguiCheck("Create Crowds", tType == TOOL_CROWD) {
		s.setTool(newCrowdTool(s.gs))
	}

	s.gs.imguiSeparatorLine()

	s.gs.imguiIndent()

	if s.m_tool != nil {
		s.m_tool.handleMenu()
	}

	s.gs.imguiUnindent()
}
func (s *Sample_TileMesh) handleRender() {}
func (s *Sample_TileMesh) handleRenderOverlay(proj, model []float64, view []int) {
	// Draw start and end point labels
	res := common.GluProject([]float64{s.m_lastBuiltTileBmin[0] + s.m_lastBuiltTileBmax[0]/2, (s.m_lastBuiltTileBmin[1] + s.m_lastBuiltTileBmax[1]) / 2, (s.m_lastBuiltTileBmin[2] + s.m_lastBuiltTileBmax[2]) / 2}, model, proj, view)
	if s.m_tileBuildTime > 0.0 && len(res) > 0 {
		x, y := int(res[0]), int(res[1])
		text := fmt.Sprintf("%.3fms / %dTris / %.1fkB", s.m_tileBuildTime, s.m_tileTriCount, s.m_tileMemUsage)
		s.gs.imguiDrawText(x, y-25, IMGUI_ALIGN_CENTER, text, imguiRGBA(0, 0, 0, 220))
	}

	if s.m_tool != nil {
		s.m_tool.handleRenderOverlay(proj, model, view)
	}

	s.renderOverlayToolStates(proj, model, view)
}
func (s *Sample_TileMesh) handleMeshChanged(geom *InputGeom) {
	s.Sample.handleMeshChanged(geom)

	buildSettings := geom.getBuildSettings()
	if buildSettings != nil && buildSettings.tileSize > 0 {
		s.m_tileSize = buildSettings.tileSize
	}

	s.cleanup()
	s.m_navMesh = nil

	if s.m_tool != nil {
		s.m_tool.reset()
		s.m_tool.init(s.Sample)
	}
	s.resetToolStates()
	s.initToolStates(s.Sample)
}
func (s *Sample_TileMesh) handleBuild() bool {
	if s.m_geom == nil || s.m_geom.getMesh() == nil {
		log.Printf("buildTiledNavigation: No vertices and triangles.")
		return false
	}

	if s.m_navMesh == nil {
		log.Printf("buildTiledNavigation: Could not allocate navmesh.")
		return false
	}

	var params recast.NavMeshParams
	copy(params.Orig[:], s.m_geom.getNavMeshBoundsMin())
	params.TileWidth = s.m_tileSize * s.m_cellSize
	params.TileHeight = s.m_tileSize * s.m_cellSize
	params.MaxTiles = s.m_maxTiles
	params.MaxPolys = s.m_maxPolysPerTile

	var status detour.DtStatus
	s.m_navMesh, status = recast.NewDtNavMeshWithParams(&params)
	if status.DtStatusFailed() {
		log.Printf("buildTiledNavigation: Could not init navmesh.")
		return false
	}

	s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048)
	if s.m_buildAll {
		s.buildAllTiles()
	}
	if s.m_tool != nil {
		s.m_tool.init(s.Sample)
	}
	s.initToolStates(s.Sample)

	return true
}
func (s *Sample_TileMesh) resetCommonSettings() {

}
func (s *Sample_TileMesh) collectSettings(settings *BuildSettings) {
	s.collectSettings(settings)
	settings.tileSize = s.m_tileSize
}

func (s *Sample_TileMesh) getTilePos(pos []float64, tx, ty *int) {
	if s.m_geom == nil {
		return
	}

	bmin := s.m_geom.getNavMeshBoundsMin()

	ts := s.m_tileSize * s.m_cellSize
	*tx = int((pos[0] - bmin[0]) / ts)
	*ty = int((pos[2] - bmin[2]) / ts)
}

func (s *Sample_TileMesh) buildTile(pos []float64) {
	if s.m_geom == nil {
		return
	}
	if s.m_navMesh == nil {
		return
	}

	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()

	ts := s.m_tileSize * s.m_cellSize
	tx := (int)((pos[0] - bmin[0]) / ts)
	ty := (int)((pos[2] - bmin[2]) / ts)

	s.m_lastBuiltTileBmin[0] = bmin[0] + float64(tx)*ts
	s.m_lastBuiltTileBmin[1] = bmin[1]
	s.m_lastBuiltTileBmin[2] = bmin[2] + float64(ty)*ts

	s.m_lastBuiltTileBmax[0] = bmin[0] + float64(tx+1)*ts
	s.m_lastBuiltTileBmax[1] = bmax[1]
	s.m_lastBuiltTileBmax[2] = bmin[2] + float64(ty+1)*ts

	s.m_tileCol = debug_utils.DuRGBA(255, 255, 255, 64)

	s.l.reset()

	data := s.buildTileMesh(tx, ty, s.m_lastBuiltTileBmin[:], s.m_lastBuiltTileBmax[:])

	// Remove any previous data (navmesh owns and deletes the data).
	s.m_navMesh.RemoveTile(s.m_navMesh.GetTileRefAt(tx, ty, 0))

	// Add tile, or leave the location empty.
	if data != nil {
		// Let the navmesh own the data.
		s.m_navMesh.AddTile(data, detour.Dt_TILE_FREE_DATA, 0)
	}

	s.l.dumpLog("Build Tile (%d,%d):", tx, ty)
}
func (s *Sample_TileMesh) removeTile(pos []float64) {
	if s.m_geom == nil || s.m_navMesh == nil {
		return
	}

	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	gw := 0
	gh := 0
	recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
	ts := int(s.m_tileSize)
	tw := (gw + ts - 1) / ts
	th := (gh + ts - 1) / ts

	for y := 0; y < th; y++ {
		for x := 0; x < tw; x++ {
			s.m_navMesh.RemoveTile(s.m_navMesh.GetTileRefAt(x, y, 0))
		}

	}

}
func (s *Sample_TileMesh) cleanup() {
	s.m_triareas = []int{}
	s.m_solid = nil
	s.m_chf = nil
	s.m_cset = nil
	s.m_pmesh = nil
	s.m_dmesh = nil
}
func (s *Sample_TileMesh) handleDebugMode() {
	// Check which modes are valid.
	valid := make([]bool, TitleMeshMAX_DRAWMODE)
	for i := 0; i < int(TitleMeshMAX_DRAWMODE); i++ {
		valid[i] = false
	}

	if s.m_geom != nil {
		valid[TitleMeshDRAWMODE_NAVMESH] = s.m_navMesh != nil
		valid[TitleMeshDRAWMODE_NAVMESH_TRANS] = s.m_navMesh != nil
		valid[TitleMeshDRAWMODE_NAVMESH_BVTREE] = s.m_navMesh != nil
		valid[TitleMeshDRAWMODE_NAVMESH_NODES] = s.m_navQuery != nil
		valid[TitleMeshDRAWMODE_NAVMESH_PORTALS] = s.m_navMesh != nil
		valid[TitleMeshDRAWMODE_NAVMESH_INVIS] = s.m_navMesh != nil
		valid[TitleMeshDRAWMODE_MESH] = true
		valid[TitleMeshDRAWMODE_VOXELS] = s.m_solid != nil
		valid[TitleMeshDRAWMODE_VOXELS_WALKABLE] = s.m_solid != nil
		valid[TitleMeshDRAWMODE_COMPACT] = s.m_chf != nil
		valid[TitleMeshDRAWMODE_COMPACT_DISTANCE] = s.m_chf != nil
		valid[TitleMeshDRAWMODE_COMPACT_REGIONS] = s.m_chf != nil
		valid[TitleMeshDRAWMODE_REGION_CONNECTIONS] = s.m_cset != nil
		valid[TitleMeshDRAWMODE_RAW_CONTOURS] = s.m_cset != nil
		valid[TitleMeshDRAWMODE_BOTH_CONTOURS] = s.m_cset != nil
		valid[TitleMeshDRAWMODE_CONTOURS] = s.m_cset != nil
		valid[TitleMeshDRAWMODE_POLYMESH] = s.m_pmesh != nil
		valid[TitleMeshDRAWMODE_POLYMESH_DETAIL] = s.m_dmesh != nil
	}

	unavail := 0
	for i := 0; i < int(TitleMeshMAX_DRAWMODE); i++ {
		if !valid[i] {
			unavail++
		}
	}

	if unavail == int(TitleMeshMAX_DRAWMODE) {
		return
	}

	s.gs.imguiLabel("Draw")
	if s.gs.imguiCheck("Input Mesh", s.m_drawMode == TitleMeshDRAWMODE_MESH, valid[TitleMeshDRAWMODE_MESH]) {
		s.m_drawMode = TitleMeshDRAWMODE_MESH
	}

	if s.gs.imguiCheck("Navmesh", s.m_drawMode == TitleMeshDRAWMODE_NAVMESH, valid[TitleMeshDRAWMODE_NAVMESH]) {
		s.m_drawMode = TitleMeshDRAWMODE_NAVMESH
	}

	if s.gs.imguiCheck("Navmesh Invis", s.m_drawMode == TitleMeshDRAWMODE_NAVMESH_INVIS, valid[TitleMeshDRAWMODE_NAVMESH_INVIS]) {
		s.m_drawMode = TitleMeshDRAWMODE_NAVMESH_INVIS
	}

	if s.gs.imguiCheck("Navmesh Trans", s.m_drawMode == TitleMeshDRAWMODE_NAVMESH_TRANS, valid[TitleMeshDRAWMODE_NAVMESH_TRANS]) {
		s.m_drawMode = TitleMeshDRAWMODE_NAVMESH_TRANS
	}

	if s.gs.imguiCheck("Navmesh BVTree", s.m_drawMode == TitleMeshDRAWMODE_NAVMESH_BVTREE, valid[TitleMeshDRAWMODE_NAVMESH_BVTREE]) {
		s.m_drawMode = TitleMeshDRAWMODE_NAVMESH_BVTREE
	}

	if s.gs.imguiCheck("Navmesh Nodes", s.m_drawMode == TitleMeshDRAWMODE_NAVMESH_NODES, valid[TitleMeshDRAWMODE_NAVMESH_NODES]) {
		s.m_drawMode = TitleMeshDRAWMODE_NAVMESH_NODES
	}

	if s.gs.imguiCheck("Navmesh Portals", s.m_drawMode == TitleMeshDRAWMODE_NAVMESH_PORTALS, valid[TitleMeshDRAWMODE_NAVMESH_PORTALS]) {
		s.m_drawMode = TitleMeshDRAWMODE_NAVMESH_PORTALS
	}

	if s.gs.imguiCheck("Voxels", s.m_drawMode == TitleMeshDRAWMODE_VOXELS, valid[TitleMeshDRAWMODE_VOXELS]) {
		s.m_drawMode = TitleMeshDRAWMODE_VOXELS
	}

	if s.gs.imguiCheck("Walkable Voxels", s.m_drawMode == TitleMeshDRAWMODE_VOXELS_WALKABLE, valid[TitleMeshDRAWMODE_VOXELS_WALKABLE]) {
		s.m_drawMode = TitleMeshDRAWMODE_VOXELS_WALKABLE
	}

	if s.gs.imguiCheck("Compact", s.m_drawMode == TitleMeshDRAWMODE_COMPACT, valid[TitleMeshDRAWMODE_COMPACT]) {
		s.m_drawMode = TitleMeshDRAWMODE_COMPACT
	}

	if s.gs.imguiCheck("Compact Distance", s.m_drawMode == TitleMeshDRAWMODE_COMPACT_DISTANCE, valid[TitleMeshDRAWMODE_COMPACT_DISTANCE]) {
		s.m_drawMode = TitleMeshDRAWMODE_COMPACT_DISTANCE
	}

	if s.gs.imguiCheck("Compact Regions", s.m_drawMode == TitleMeshDRAWMODE_COMPACT_REGIONS, valid[TitleMeshDRAWMODE_COMPACT_REGIONS]) {
		s.m_drawMode = TitleMeshDRAWMODE_COMPACT_REGIONS
	}

	if s.gs.imguiCheck("Region Connections", s.m_drawMode == TitleMeshDRAWMODE_REGION_CONNECTIONS, valid[TitleMeshDRAWMODE_REGION_CONNECTIONS]) {
		s.m_drawMode = TitleMeshDRAWMODE_REGION_CONNECTIONS
	}

	if s.gs.imguiCheck("Raw Contours", s.m_drawMode == TitleMeshDRAWMODE_RAW_CONTOURS, valid[TitleMeshDRAWMODE_RAW_CONTOURS]) {
		s.m_drawMode = TitleMeshDRAWMODE_RAW_CONTOURS
	}

	if s.gs.imguiCheck("Both Contours", s.m_drawMode == TitleMeshDRAWMODE_BOTH_CONTOURS, valid[TitleMeshDRAWMODE_BOTH_CONTOURS]) {
		s.m_drawMode = TitleMeshDRAWMODE_BOTH_CONTOURS
	}

	if s.gs.imguiCheck("Contours", s.m_drawMode == TitleMeshDRAWMODE_CONTOURS, valid[TitleMeshDRAWMODE_CONTOURS]) {
		s.m_drawMode = TitleMeshDRAWMODE_CONTOURS
	}

	if s.gs.imguiCheck("Poly Mesh", s.m_drawMode == TitleMeshDRAWMODE_POLYMESH, valid[TitleMeshDRAWMODE_POLYMESH]) {
		s.m_drawMode = TitleMeshDRAWMODE_POLYMESH
	}

	if s.gs.imguiCheck("Poly Mesh Detail", s.m_drawMode == TitleMeshDRAWMODE_POLYMESH_DETAIL, valid[TitleMeshDRAWMODE_POLYMESH_DETAIL]) {
		s.m_drawMode = TitleMeshDRAWMODE_POLYMESH_DETAIL
	}

	if unavail > 0 {
		s.gs.imguiValue("Tick 'Keep Itermediate Results'")
		s.gs.imguiValue("rebuild some tiles to see")
		s.gs.imguiValue("more debug mode options.")
	}
}
func (s *Sample_TileMesh) buildAllTiles() {
	if s.m_geom == nil {
		return
	}
	if s.m_navMesh == nil {
		return
	}

	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	var gw = 0
	var gh = 0
	recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
	ts := int(s.m_tileSize)
	tw := (gw + ts - 1) / ts
	th := gh + ts - 1/ts
	tcs := s.m_tileSize * s.m_cellSize

	// Start the build process.
	now := time.Now()
	for y := 0; y < th; y++ {
		for x := 0; x < tw; x++ {
			s.m_lastBuiltTileBmin[0] = bmin[0] + float64(x)*tcs
			s.m_lastBuiltTileBmin[1] = bmin[1]
			s.m_lastBuiltTileBmin[2] = bmin[2] + float64(y)*tcs

			s.m_lastBuiltTileBmax[0] = bmin[0] + float64(x+1)*tcs
			s.m_lastBuiltTileBmax[1] = bmax[1]
			s.m_lastBuiltTileBmax[2] = bmin[2] + float64(y+1)*tcs
			data := s.buildTileMesh(x, y, s.m_lastBuiltTileBmin[:], s.m_lastBuiltTileBmax[:])
			if data != nil {
				// Remove any previous data (navmesh owns and deletes the data).
				s.m_navMesh.RemoveTile(s.m_navMesh.GetTileRefAt(x, y, 0))
				// Let the navmesh own the data.
				s.m_navMesh.AddTile(data, detour.Dt_TILE_FREE_DATA, 0)
			}
		}
	}

	s.m_totalBuildTimeMs = time.Since(now)
}
func (s *Sample_TileMesh) removeAllTiles() {
	if s.m_geom == nil || s.m_navMesh == nil {
		return
	}

	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	gw := 0
	gh := 0
	recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
	ts := int(s.m_tileSize)
	tw := (gw + ts - 1) / ts
	th := (gh + ts - 1) / ts

	for y := 0; y < th; y++ {
		for x := 0; x < tw; x++ {
			s.m_navMesh.RemoveTile(s.m_navMesh.GetTileRefAt(x, y, 0))
		}
	}

}
func (s *Sample_TileMesh) buildTileMesh(tx, ty int, bmin, bmax []float64) (data *recast.NavMeshData) {
	if s.m_geom == nil || s.m_geom.getMesh() == nil || s.m_geom.getChunkyMesh() == nil {
		log.Printf("buildNavigation: Input rcMeshLoaderObj is not specified.")
		return
	}
	now := time.Now()
	s.m_tileMemUsage = 0
	s.m_tileBuildTime = 0

	s.cleanup()

	verts := s.m_geom.getMesh().getVerts()
	nverts := s.m_geom.getMesh().getVertCount()
	ntris := s.m_geom.getMesh().getTriCount()
	chunkyMesh := s.m_geom.getChunkyMesh()

	// Init build configuration from GUI
	s.m_cfg = &recast.RcConfig{}
	m_cfg := s.m_cfg
	m_cfg.Cs = s.m_cellSize
	m_cfg.Ch = s.m_cellHeight
	m_cfg.WalkableSlopeAngle = s.m_agentMaxSlope
	m_cfg.WalkableHeight = int(math.Ceil(s.m_agentHeight / m_cfg.Ch))
	m_cfg.WalkableClimb = int(math.Floor(s.m_agentMaxClimb / m_cfg.Ch))
	m_cfg.WalkableRadius = int(math.Ceil(s.m_agentRadius / m_cfg.Cs))
	m_cfg.MaxEdgeLen = int(s.m_edgeMaxLen / s.m_cellSize)
	m_cfg.MaxSimplificationError = s.m_edgeMaxError
	m_cfg.MinRegionArea = int(common.Sqr(s.m_regionMinSize))     // Note: area = size*size
	m_cfg.MergeRegionArea = int(common.Sqr(s.m_regionMergeSize)) // Note: area = size*size
	m_cfg.MaxVertsPerPoly = int(s.m_vertsPerPoly)
	m_cfg.TileSize = int(s.m_tileSize)
	m_cfg.BorderSize = m_cfg.WalkableRadius + 3 // Reserve enough padding.
	m_cfg.Width = m_cfg.TileSize + m_cfg.BorderSize*2
	m_cfg.Height = m_cfg.TileSize + m_cfg.BorderSize*2
	m_cfg.DetailSampleDist = s.m_cellSize * s.m_detailSampleDist
	if s.m_detailSampleDist < 0.9 {
		m_cfg.DetailSampleDist = 0
	}
	m_cfg.DetailSampleMaxError = s.m_cellHeight * s.m_detailSampleMaxError

	// Expand the heighfield bounding box by border size to find the extents of geometry we need to build this tile.
	//
	// This is done in order to make sure that the navmesh tiles connect correctly at the borders,
	// and the obstacles close to the border work correctly with the dilation process.
	// No polygons (or contours) will be created on the border area.
	//
	// IMPORTANT!
	//
	//   :''''''''':
	//   : +-----+ :
	//   : |     | :
	//   : |     |<--- tile to build
	//   : |     | :
	//   : +-----+ :<-- geometry needed
	//   :.........:
	//
	// You should use this bounding box to query your input geometry.
	//
	// For example if you build a navmesh for terrain, and want the navmesh tiles to match the terrain tile size
	// you will need to pass in data from neighbour terrain tiles too! In a simple case, just pass in all the 8 neighbours,
	// or use the bounding box below to only pass in a sliver of each of the 8 neighbours.
	copy(m_cfg.Bmin[:], bmin)
	copy(m_cfg.Bmax[:], bmax)
	m_cfg.Bmin[0] -= float64(m_cfg.BorderSize) * m_cfg.Cs
	m_cfg.Bmin[2] -= float64(m_cfg.BorderSize) * m_cfg.Cs
	m_cfg.Bmax[0] += float64(m_cfg.BorderSize) * m_cfg.Cs
	m_cfg.Bmax[2] += float64(m_cfg.BorderSize) * m_cfg.Cs

	// Reset build times gathering.
	// Start the build process.
	log.Printf("Building navigation:")
	log.Printf(" - %d x %d cells", m_cfg.Width, m_cfg.Height)
	log.Printf(" - %.1fK verts, %.1fK tris", nverts/1000.0, ntris/1000.0)

	// Allocate voxel heightfield where we rasterize our input data to.
	s.m_solid = recast.RcCreateHeightfield(m_cfg.Width, m_cfg.Height, m_cfg.Bmin[:], m_cfg.Bmax[:], m_cfg.Cs, m_cfg.Ch)
	if s.m_solid == nil {
		log.Printf("buildNavigation: Could not create solid heightfield.")
		return
	}

	// Allocate array that can hold triangle flags.
	// If you have multiple meshes you need to process, allocate
	// and array which can hold the max number of triangles you need to process.
	s.m_triareas = []int{}
	var tbmin [2]float64
	var tbmax [2]float64
	tbmin[0] = m_cfg.Bmin[0]
	tbmin[1] = m_cfg.Bmin[2]
	tbmax[0] = m_cfg.Bmax[0]
	tbmax[1] = m_cfg.Bmax[2]
	cid := make([]int, 512) // TODO: Make grow when returning too many items.
	ncid := rcGetChunksOverlappingRect(chunkyMesh, tbmin, tbmax, cid, 512)
	if ncid == 0 {
		return
	}

	s.m_tileTriCount = 0

	for i := 0; i < ncid; i++ {
		node := chunkyMesh.nodes[cid[i]]
		ctris := chunkyMesh.tris[node.i*3:]
		nctris := node.n

		s.m_tileTriCount += nctris
		s.m_triareas = make([]int, nctris)
		recast.RcMarkWalkableTriangles(m_cfg.WalkableSlopeAngle,
			verts, nverts, ctris, nctris, s.m_triareas)

		if !recast.RcRasterizeTriangles(verts, nverts, ctris, s.m_triareas, nctris, s.m_solid, m_cfg.WalkableClimb) {
			return
		}

	}

	if !s.m_keepInterResults {

		s.m_triareas = []int{}
	}

	// Once all geometry is rasterized, we do initial pass of filtering to
	// remove unwanted overhangs caused by the conservative rasterization
	// as well as filter spans where the character cannot possibly stand.
	if s.m_filterLowHangingObstacles {
		recast.RcFilterLowHangingWalkableObstacles(m_cfg.WalkableClimb, s.m_solid)
	}

	if s.m_filterLedgeSpans {
		recast.RcFilterLedgeSpans(m_cfg.WalkableHeight, m_cfg.WalkableClimb, s.m_solid)
	}

	if s.m_filterWalkableLowHeightSpans {
		recast.RcFilterWalkableLowHeightSpans(m_cfg.WalkableHeight, s.m_solid)
	}

	// Compact the heightfield so that it is faster to handle from now on.
	// This will result more cache coherent data as well as the neighbours
	// between walkable cells will be calculated.
	s.m_chf = &recast.RcCompactHeightfield{}
	if !recast.RcBuildCompactHeightfield(m_cfg.WalkableHeight, m_cfg.WalkableClimb, s.m_solid, s.m_chf) {
		log.Printf("buildNavigation: Could not build compact data.")
		return nil
	}

	if !s.m_keepInterResults {
		s.m_solid = nil
	}

	// Erode the walkable area by agent radius.
	if !recast.RcErodeWalkableArea(m_cfg.WalkableRadius, s.m_chf) {
		log.Printf("buildNavigation: Could not erode.")
		return
	}

	// (Optional) Mark areas.
	vols := s.m_geom.getConvexVolumes()
	for i := 0; i < s.m_geom.getConvexVolumeCount(); i++ {
		recast.RcMarkConvexPolyArea(vols[i].verts[:], vols[i].nverts, vols[i].hmin, vols[i].hmax, vols[i].area, s.m_chf)
	}

	// Partition the heightfield so that we can use simple algorithm later to triangulate the walkable areas.
	// There are 3 martitioning methods, each with some pros and cons:
	// 1) Watershed partitioning
	//   - the classic Recast partitioning
	//   - creates the nicest tessellation
	//   - usually slowest
	//   - partitions the heightfield into nice regions without holes or overlaps
	//   - the are some corner cases where this method creates produces holes and overlaps
	//      - holes may appear when a small obstacles is close to large open area (triangulation can handle this)
	//      - overlaps may occur if you have narrow spiral corridors (i.e stairs), this make triangulation to fail
	//   * generally the best choice if you precompute the nacmesh, use this if you have large open areas
	// 2) Monotone partioning
	//   - fastest
	//   - partitions the heightfield into regions without holes and overlaps (guaranteed)
	//   - creates long thin polygons, which sometimes causes paths with detours
	//   * use this if you want fast navmesh generation
	// 3) Layer partitoining
	//   - quite fast
	//   - partitions the heighfield into non-overlapping regions
	//   - relies on the triangulation code to cope with holes (thus slower than monotone partitioning)
	//   - produces better triangles than monotone partitioning
	//   - does not have the corner cases of watershed partitioning
	//   - can be slow and create a bit ugly tessellation (still better than monotone)
	//     if you have large open areas with small obstacles (not a problem if you use tiles)
	//   * good choice to use for tiled navmesh with medium and small sized tiles

	if s.m_partitionType == SAMPLE_PARTITION_WATERSHED {
		// Prepare for region partitioning, by calculating distance field along the walkable surface.
		if !recast.RcBuildDistanceField(s.m_chf) {
			log.Printf("buildNavigation: Could not build distance field.")
			return
		}

		// Partition the walkable surface into simple regions without holes.
		if !recast.RcBuildRegions(s.m_chf, m_cfg.BorderSize, m_cfg.MinRegionArea, m_cfg.MergeRegionArea) {
			log.Printf("buildNavigation: Could not build watershed regions.")
			return
		}
	} else if s.m_partitionType == SAMPLE_PARTITION_MONOTONE {
		// Partition the walkable surface into simple regions without holes.
		// Monotone partitioning does not need distancefield.
		if !recast.RcBuildRegionsMonotone(s.m_chf, m_cfg.BorderSize, m_cfg.MinRegionArea, m_cfg.MergeRegionArea) {
			log.Printf("buildNavigation: Could not build monotone regions.")
			return
		}
	} else // SAMPLE_PARTITION_LAYERS
	{
		// Partition the walkable surface into simple regions without holes.
		if !recast.RcBuildLayerRegions(s.m_chf, m_cfg.BorderSize, m_cfg.MinRegionArea) {
			log.Printf("buildNavigation: Could not build layer regions.")
			return
		}
	}

	// Create contours.
	s.m_cset = &recast.RcContourSet{}
	if !recast.RcBuildContours(s.m_chf, m_cfg.MaxSimplificationError, m_cfg.MaxEdgeLen, s.m_cset, 0) {
		log.Printf("buildNavigation: Could not create contours.")
		return
	}

	if s.m_cset.Nconts == 0 {
		return
	}

	// Build polygon navmesh from the contours.
	s.m_pmesh = &recast.RcPolyMesh{}
	if !recast.RcBuildPolyMesh(s.m_cset, m_cfg.MaxVertsPerPoly, s.m_pmesh) {
		log.Printf("buildNavigation: Could not triangulate contours.")
		return
	}

	// Build detail rcMeshLoaderObj.
	s.m_dmesh = &recast.RcPolyMeshDetail{}
	if !recast.RcBuildPolyMeshDetail(s.m_pmesh, s.m_chf, m_cfg.DetailSampleDist, m_cfg.DetailSampleMaxError,
		s.m_dmesh) {
		log.Printf("buildNavigation: Could build polymesh detail.")
		return
	}

	if !s.m_keepInterResults {
		s.m_chf = nil
		s.m_cset = nil
	}
	if m_cfg.MaxVertsPerPoly <= detour.Dt_VERTS_PER_POLYGON {
		if s.m_pmesh.Nverts >= 0xffff {
			// The vertex indices are ushorts, and cannot point to more than 0xffff vertices.
			log.Printf("Too many vertices per tile %d (max: %d).", s.m_pmesh.Nverts, 0xffff)
			return
		}

		// Update poly flags from areas.
		for i := 0; i < s.m_pmesh.Npolys; i++ {
			if s.m_pmesh.Areas[i] == recast.RC_WALKABLE_AREA {
				s.m_pmesh.Areas[i] = SAMPLE_POLYAREA_GROUND
			}

			if s.m_pmesh.Areas[i] == SAMPLE_POLYAREA_GROUND ||
				s.m_pmesh.Areas[i] == SAMPLE_POLYAREA_GRASS ||
				s.m_pmesh.Areas[i] == SAMPLE_POLYAREA_ROAD {
				s.m_pmesh.Flags[i] = SAMPLE_POLYFLAGS_WALK
			} else if s.m_pmesh.Areas[i] == SAMPLE_POLYAREA_WATER {
				s.m_pmesh.Flags[i] = SAMPLE_POLYFLAGS_SWIM
			} else if s.m_pmesh.Areas[i] == SAMPLE_POLYAREA_DOOR {
				s.m_pmesh.Flags[i] = SAMPLE_POLYFLAGS_WALK | SAMPLE_POLYFLAGS_DOOR
			}
		}

		params := &detour.DtNavMeshCreateParams{}
		params.Verts = s.m_pmesh.Verts
		params.VertCount = s.m_pmesh.Nverts
		params.Polys = s.m_pmesh.Polys
		params.PolyAreas = s.m_pmesh.Areas
		params.PolyFlags = s.m_pmesh.Flags
		params.PolyCount = s.m_pmesh.Npolys
		params.Nvp = s.m_pmesh.Nvp
		params.DetailMeshes = s.m_dmesh.Meshes
		params.DetailVerts = s.m_dmesh.Verts
		params.DetailVertsCount = s.m_dmesh.Nverts
		params.DetailTris = s.m_dmesh.Tris
		params.DetailTriCount = s.m_dmesh.Ntris
		params.OffMeshConVerts = s.m_geom.getOffMeshConnectionVerts()
		params.OffMeshConRad = s.m_geom.getOffMeshConnectionRads()
		params.OffMeshConDir = s.m_geom.getOffMeshConnectionDirs()
		params.OffMeshConAreas = s.m_geom.getOffMeshConnectionAreas()
		params.OffMeshConFlags = s.m_geom.getOffMeshConnectionFlags()
		params.OffMeshConUserID = s.m_geom.getOffMeshConnectionId()
		params.OffMeshConCount = s.m_geom.getOffMeshConnectionCount()
		params.WalkableHeight = s.m_agentHeight
		params.WalkableRadius = s.m_agentRadius
		params.WalkableClimb = s.m_agentMaxClimb
		params.TileX = tx
		params.TileY = ty
		params.TileLayer = 0
		copy(params.Bmin[:], s.m_pmesh.Bmin)
		copy(params.Bmax[:], s.m_pmesh.Bmax)
		params.Cs = m_cfg.Cs
		params.Ch = m_cfg.Ch
		params.BuildBvTree = true
		var ok bool
		data, ok = detour.DtCreateNavMeshData(params)
		if !ok {
			log.Printf("Could not build Detour navmesh.")
			return
		}
	}
	// Show performance stats.
	log.Printf(">> Polymesh: %d vertices  %d polygons", s.m_pmesh.Nverts, s.m_pmesh.Npolys)

	s.m_tileBuildTime = time.Since(now)
	return data
}
