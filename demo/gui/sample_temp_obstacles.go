package gui

import (
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"gonavamesh/common"
	"gonavamesh/debug_utils"
	"gonavamesh/recast"
	"log"
	"math"
	"time"
)

// This value specifies how many layers (or "floors") each navmesh tile is expected to have.
const EXPECTED_LAYERS_PER_TILE = 4

func drawTiles(dd debug_utils.DuDebugDraw, tc *recast.DtTileCache) {
	fcol := make([]int, 6)
	bmin := make([]float64, 3)
	bmax := make([]float64, 3)

	for i := 0; i < tc.GetTileCount(); i++ {
		tile := tc.GetTile(i)
		if tile.Header == nil {
			continue
		}

		tc.CalcTightTileBounds(tile.Header, bmin, bmax)

		col := debug_utils.DuIntToCol(i, 64)
		debug_utils.DuCalcBoxColors(fcol, col, col)
		debug_utils.DuDebugDrawBox(dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], fcol)
	}

	for i := 0; i < tc.GetTileCount(); i++ {
		tile := tc.GetTile(i)
		if tile.Header == nil {
			continue
		}

		tc.CalcTightTileBounds(tile.Header, bmin, bmax)

		col := debug_utils.DuIntToCol(i, 255)
		pad := tc.GetParams().Cs * 0.1
		debug_utils.DuDebugDrawBoxWire(dd, bmin[0]-pad, bmin[1]-pad, bmin[2]-pad,
			bmax[0]+pad, bmax[1]+pad, bmax[2]+pad, col, 2.0)
	}

}
func hitTestObstacle(tc *recast.DtTileCache, sp, sq []float64) recast.DtObstacleRef {
	tmin := math.MaxFloat64
	var obmin *recast.DtTileCacheObstacle
	for i := 0; i < tc.GetObstacleCount(); i++ {
		ob := tc.GetObstacle(i)
		if ob.State == recast.DT_OBSTACLE_EMPTY {
			continue
		}

		bmin := make([]float64, 3)
		bmax := make([]float64, 3)
		var t0, t1 float64
		tc.GetObstacleBounds(ob, bmin, bmax)

		if isectSegAABB(sp, sq, bmin, bmax, t0, t1) {
			if t0 < tmin {
				tmin = t0
				obmin = ob
			}
		}
	}
	return tc.GetObstacleRef(obmin)
}
func drawDetailOverlay(gs *guiState, tc *recast.DtTileCache, tx, ty int, proj, model []float64, view []int) {
	tiles := make([]recast.DtCompressedTileRef, MAX_LAYERS)
	ntiles := tc.GetTilesAt(tx, ty, tiles, MAX_LAYERS)
	if ntiles == 0 {
		return
	}
	for i := 0; i < ntiles; i++ {
		tile := tc.GetTileByRef(tiles[i])

		pos := make([]float64, 3)
		pos[0] = (tile.Header.Bmin[0] + tile.Header.Bmax[0]) / 2.0
		pos[1] = tile.Header.Bmin[1]
		pos[2] = (tile.Header.Bmin[2] + tile.Header.Bmax[2]) / 2.0

		res := common.GluProject([]float64{pos[0], pos[1], pos[2]}, model, proj, view)

		if len(res) != 0 {
			x, y := int(res[0]), int(res[1])
			text := fmt.Sprintf("(%d,%d)/%d", tile.Header.Tx, tile.Header.Ty, tile.Header.Tlayer)
			gs.imguiDrawText(x, y-25, IMGUI_ALIGN_CENTER, text, imguiRGBA(0, 0, 0, 220))
			//text = fmt.Sprintf("Compressed: %.1f kB", tile.DataSize/1024.0)
			gs.imguiDrawText(x, y-45, IMGUI_ALIGN_CENTER, text, imguiRGBA(0, 0, 0, 128))
			//text = fmt.Sprintf("Raw:%.1fkB", rawSize/1024.0)
			gs.imguiDrawText(x, y-65, IMGUI_ALIGN_CENTER, text, imguiRGBA(0, 0, 0, 128))
		}
	}
}

func drawObstacles(dd debug_utils.DuDebugDraw, tc *recast.DtTileCache) {
	// Draw obstacles
	for i := 0; i < tc.GetObstacleCount(); i++ {
		ob := tc.GetObstacle(i)
		if ob.State == recast.DT_OBSTACLE_EMPTY {
			continue
		}
		bmin := make([]float64, 3)
		bmax := make([]float64, 3)
		tc.GetObstacleBounds(ob, bmin, bmax)

		var col int
		if ob.State == recast.DT_OBSTACLE_PROCESSING {
			col = debug_utils.DuRGBA(255, 255, 0, 128)
		} else if ob.State == recast.DT_OBSTACLE_PROCESSED {
			col = debug_utils.DuRGBA(255, 192, 0, 192)
		} else if ob.State == recast.DT_OBSTACLE_REMOVING {
			col = debug_utils.DuRGBA(220, 0, 0, 128)
		}

		debug_utils.DuDebugDrawCylinder(dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], col)
		debug_utils.DuDebugDrawCylinderWire(dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], debug_utils.DuDarkenCol(col), 2)
	}
}

const (
	DRAWDETAIL_AREAS = iota
	DRAWDETAIL_REGIONS
	DRAWDETAIL_CONTOURS
	DRAWDETAIL_MESH
)

type TileCacheBuildContext struct {
	layer *recast.DtTileCacheLayer
	lcset *recast.DtTileCacheContourSet
	lmesh *recast.DtTileCachePolyMesh
}

func newTileCacheBuildContext() *TileCacheBuildContext {
	return &TileCacheBuildContext{
		layer: &recast.DtTileCacheLayer{},
		lcset: &recast.DtTileCacheContourSet{},
		lmesh: &recast.DtTileCachePolyMesh{},
	}
}
func (t *TileCacheBuildContext) purge() {
	t.layer = nil
	t.lcset = nil
	t.lmesh = nil
}
func drawDetail(dd debug_utils.DuDebugDraw, tc *recast.DtTileCache, tx, ty, tType int) {

	tiles := make([]recast.DtCompressedTileRef, MAX_LAYERS)
	ntiles := tc.GetTilesAt(tx, ty, tiles, MAX_LAYERS)
	params := tc.GetParams()
	for i := 0; i < ntiles; i++ {
		tile := tc.GetTileByRef(tiles[i])
		bc := newTileCacheBuildContext()
		walkableClimbVx := int(params.WalkableClimb / params.Ch)
		// Decompress tile layer data.
		data := tile.Data
		bc.layer = &recast.DtTileCacheLayer{
			Header:   data.Header,
			RegCount: data.RegCount,
			Heights:  data.Heights,
			Areas:    data.Areas,
			Cons:     data.Cons,
			Regs:     data.Regs,
		}
		if tType == DRAWDETAIL_AREAS {
			debug_utils.DuDebugDrawTileCacheLayerAreas(dd, bc.layer, params.Cs, params.Ch)
			continue
		}

		// Build navmesh
		status := recast.DtBuildTileCacheRegions(bc.layer, walkableClimbVx)
		if status.DtStatusFailed() {
			return
		}

		if tType == DRAWDETAIL_REGIONS {
			debug_utils.DuDebugDrawTileCacheLayerRegions(dd, bc.layer, params.Cs, params.Ch)
			continue
		}

		bc.lcset = &recast.DtTileCacheContourSet{}

		status = recast.DtBuildTileCacheContours(bc.layer, walkableClimbVx,
			params.MaxSimplificationError, bc.lcset)
		if status.DtStatusFailed() {
			return
		}

		if tType == DRAWDETAIL_CONTOURS {
			debug_utils.DuDebugDrawTileCacheContours(dd, bc.lcset, tile.Header.Bmin[:], params.Cs, params.Ch)
			continue
		}

		bc.lmesh = &recast.DtTileCachePolyMesh{}
		status = recast.DtBuildTileCachePolyMesh(bc.lcset, bc.lmesh)
		if status.DtStatusFailed() {
			return
		}

		if tType == DRAWDETAIL_MESH {
			debug_utils.DuDebugDrawTileCachePolyMesh(dd, bc.lmesh, tile.Header.Bmin[:], params.Cs, params.Ch)
			continue
		}

	}
}

type SampleTempObstacleDrawMode int

const (
	SampleTempObstacleDRAWMODE_NAVMESH SampleTempObstacleDrawMode = iota
	SampleTempObstacleDRAWMODE_NAVMESH_TRANS
	SampleTempObstacleDRAWMODE_NAVMESH_BVTREE
	SampleTempObstacleDRAWMODE_NAVMESH_NODES
	SampleTempObstacleDRAWMODE_NAVMESH_PORTALS
	SampleTempObstacleDRAWMODE_NAVMESH_INVIS
	SampleTempObstacleDRAWMODE_MESH
	SampleTempObstacleDRAWMODE_CACHE_BOUNDS
	SampleTempObstacleMAX_DRAWMODE
)

type SampleTempObstacles struct {
	*Sample
	m_keepInterResults bool
	m_tcomp            *FastLZCompressor
	m_tmproc           *MeshProcess

	m_tileCache *recast.DtTileCache

	m_cacheBuildTimeMs   time.Duration
	m_cacheLayerCount    int
	m_cacheBuildMemUsage int

	m_drawMode SampleTempObstacleDrawMode

	m_maxTiles        int
	m_maxPolysPerTile int
	m_tileSize        float64
}

func newSampleTempObstacles() *SampleTempObstacles {
	s := &SampleTempObstacles{}
	return s
}

func (s *SampleTempObstacles) handleSettings() {
	s.Sample.handleCommonSettings()

	if s.gs.imguiCheck("Keep Itermediate Results", s.m_keepInterResults) {
		s.m_keepInterResults = !s.m_keepInterResults
	}

	s.gs.imguiLabel("Tiling")
	s.gs.imguiSlider("TileSize", &s.m_tileSize, 16.0, 128.0, 8.0)

	gridSize := 1
	if s.m_geom != nil {
		bmin := s.m_geom.getNavMeshBoundsMin()
		bmax := s.m_geom.getNavMeshBoundsMax()
		gw := 0
		gh := 0
		recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
		ts := int(s.m_tileSize)
		tw := (gw + ts - 1) / ts
		th := (gh + ts - 1) / ts
		text := fmt.Sprintf("Tiles  %d x %d", tw, th)
		s.gs.imguiValue(text)

		// Max tiles and max polys affect how the tile IDs are caculated.
		// There are 22 bits available for identifying a tile and a polygon.
		tileBits := common.Min(common.Ilog2(common.NextPow2(tw*th*EXPECTED_LAYERS_PER_TILE)), 14)
		if tileBits > 14 {
			tileBits = 14
		}
		polyBits := 22 - tileBits
		s.m_maxTiles = 1 << tileBits
		s.m_maxPolysPerTile = 1 << polyBits
		text = fmt.Sprintf("Max Tiles  %d", s.m_maxTiles)
		s.gs.imguiValue(text)
		text = fmt.Sprintf("Max Polys  %d", s.m_maxPolysPerTile)
		s.gs.imguiValue(text)
		gridSize = tw * th
	} else {
		s.m_maxTiles = 0
		s.m_maxPolysPerTile = 0
	}

	s.gs.imguiSeparator()

	s.gs.imguiLabel("Tile Cache")

	//compressionRatio := float64(s.m_cacheCompressedSize) / float64(s.m_cacheRawSize+1)

	msg := fmt.Sprintf("Layers  %d", s.m_cacheLayerCount)
	s.gs.imguiValue(msg)
	msg = fmt.Sprintf("Layers (per tile)  %.1f", float64(s.m_cacheLayerCount)/float64(gridSize))
	s.gs.imguiValue(msg)

	//msg = fmt.Sprintf("Memory  %.1f kB / %.1f kB (%.1f%%)", float64(s.m_cacheCompressedSize)/1024.0, float64(s.m_cacheRawSize)/1024.0, compressionRatio*100.0)
	s.gs.imguiValue(msg)
	msg = fmt.Sprintf("Navmesh Build Time  %.1f ms", s.m_cacheBuildTimeMs)
	s.gs.imguiValue(msg)
	msg = fmt.Sprintf("Build Peak Mem Usage  %.1f kB", float64(s.m_cacheBuildMemUsage)/1024.0)
	s.gs.imguiValue(msg)

	s.gs.imguiSeparator()

	s.gs.imguiIndent()
	s.gs.imguiIndent()

	if s.gs.imguiButton("Save") {
		s.saveAll("all_tiles_tilecache.bin")
	}

	if s.gs.imguiButton("Load") {
		s.loadAll("all_tiles_tilecache.bin")
		s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048)
	}

	s.gs.imguiUnindent()
	s.gs.imguiUnindent()

	s.gs.imguiSeparator()
}
func (s *SampleTempObstacles) handleTools() {
	Type := s.m_tool.Type()
	if s.m_tool == nil {
		Type = TOOL_NONE
	}

	if s.gs.imguiCheck("Test Navmesh", Type == TOOL_NAVMESH_TESTER) {
		s.setTool(newNavMeshTesterTool(s.gs))
	}
	if s.gs.imguiCheck("Highlight Tile Cache", Type == TOOL_TILE_HIGHLIGHT) {
		s.setTool(newTempObstacleHilightTool(s.gs))
	}
	if s.gs.imguiCheck("Create Temp Obstacles", Type == TOOL_TEMP_OBSTACLE) {
		s.setTool(newTempObstacleCreateTool(s.gs))
	}
	if s.gs.imguiCheck("Create Off-Mesh Links", Type == TOOL_OFFMESH_CONNECTION) {
		s.setTool(newOffMeshConnectionTool(s.gs))
	}
	if s.gs.imguiCheck("Create Convex Volumes", Type == TOOL_CONVEX_VOLUME) {
		s.setTool(newConvexVolumeTool(s.gs))
	}
	if s.gs.imguiCheck("Create Crowds", Type == TOOL_CROWD) {
		s.setTool(newCrowdTool(s.gs))
	}

	s.gs.imguiSeparatorLine()

	s.gs.imguiIndent()

	if s.m_tool != nil {
		s.m_tool.handleMenu()
	}

	s.gs.imguiUnindent()
}
func (s *SampleTempObstacles) handleDebugMode() {
	// Check which modes are valid.
	valid := make([]bool, SampleTempObstacleMAX_DRAWMODE)
	for i := 0; i < int(SampleTempObstacleMAX_DRAWMODE); i++ {
		valid[i] = false
	}

	if s.m_geom != nil {
		valid[SampleTempObstacleDRAWMODE_NAVMESH] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_TRANS] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_BVTREE] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_NODES] = s.m_navQuery != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_PORTALS] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_INVIS] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_MESH] = true
		valid[SampleTempObstacleDRAWMODE_CACHE_BOUNDS] = true
	}

	unavail := 0
	for i := 0; i < int(SampleTempObstacleMAX_DRAWMODE); i++ {
		if !valid[i] {
			unavail++
		}
	}

	if unavail == int(SampleTempObstacleMAX_DRAWMODE) {
		return
	}

	s.gs.imguiLabel("Draw")
	if s.gs.imguiCheck("Input Mesh", s.m_drawMode == SampleTempObstacleDRAWMODE_MESH, valid[SampleTempObstacleDRAWMODE_MESH]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_MESH
	}

	if s.gs.imguiCheck("Navmesh", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH, valid[SampleTempObstacleDRAWMODE_NAVMESH]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH
	}

	if s.gs.imguiCheck("Navmesh Invis", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_INVIS, valid[SampleTempObstacleDRAWMODE_NAVMESH_INVIS]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_INVIS
	}

	if s.gs.imguiCheck("Navmesh Trans", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_TRANS, valid[SampleTempObstacleDRAWMODE_NAVMESH_TRANS]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_TRANS
	}

	if s.gs.imguiCheck("Navmesh BVTree", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_BVTREE, valid[SampleTempObstacleDRAWMODE_NAVMESH_BVTREE]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_BVTREE
	}

	if s.gs.imguiCheck("Navmesh Nodes", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_NODES, valid[SampleTempObstacleDRAWMODE_NAVMESH_NODES]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_NODES
	}

	if s.gs.imguiCheck("Navmesh Portals", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_PORTALS, valid[SampleTempObstacleDRAWMODE_NAVMESH_PORTALS]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_PORTALS
	}

	if s.gs.imguiCheck("Cache Bounds", s.m_drawMode == SampleTempObstacleDRAWMODE_CACHE_BOUNDS, valid[SampleTempObstacleDRAWMODE_CACHE_BOUNDS]) {
		s.m_drawMode = SampleTempObstacleDRAWMODE_CACHE_BOUNDS
	}

	if unavail > 0 {
		s.gs.imguiValue("Tick 'Keep Itermediate Results'")
		s.gs.imguiValue("rebuild some tiles to see")
		s.gs.imguiValue("more debug mode options.")
	}
}
func (s *SampleTempObstacles) handleRender() {
	if s.m_geom == nil || s.m_geom.getMesh() == nil {
		return
	}

	texScale := 1.0 / (s.m_cellSize * 10.0)

	// Draw mesh
	if s.m_drawMode != SampleTempObstacleDRAWMODE_NAVMESH_TRANS {
		// Draw mesh
		debug_utils.DuDebugDrawTriMeshSlope(s.m_dd, s.m_geom.getMesh().getVerts(), s.m_geom.getMesh().getVertCount(),
			s.m_geom.getMesh().getTris(), s.m_geom.getMesh().getNormals(), s.m_geom.getMesh().getTriCount(),
			s.m_agentMaxSlope, texScale)
		s.m_geom.drawOffMeshConnections(s.m_dd)
	}

	if s.m_tileCache != nil && s.m_drawMode == SampleTempObstacleDRAWMODE_CACHE_BOUNDS {
		drawTiles(s.m_dd, s.m_tileCache)
	}

	if s.m_tileCache != nil {
		drawObstacles(s.m_dd, s.m_tileCache)
	}

	gl.DepthMask(false)

	// Draw bounds
	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	debug_utils.DuDebugDrawBoxWire(s.m_dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], debug_utils.DuRGBA(255, 255, 255, 128), 1.0)

	// Tiling grid.
	gw := 0
	gh := 0
	recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
	tw := int((float64(gw) + s.m_tileSize - 1) / s.m_tileSize)
	th := int((float64(gh) + s.m_tileSize - 1) / s.m_tileSize)
	ss := float64(s.m_tileSize) * float64(s.m_cellSize)
	debug_utils.DuDebugDrawGridXZ(s.m_dd, bmin[0], bmin[1], bmin[2], tw, th, ss, debug_utils.DuRGBA(0, 0, 0, 64), 1.0)

	if s.m_navMesh != nil && s.m_navQuery != nil &&
		(s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH ||
			s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_TRANS ||
			s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_BVTREE ||
			s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_NODES ||
			s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_PORTALS ||
			s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_INVIS) {
		if s.m_drawMode != SampleTempObstacleDRAWMODE_NAVMESH_INVIS {
			debug_utils.DuDebugDrawNavMeshWithClosedList(s.m_dd, s.m_navMesh, s.m_navQuery, s.m_navMeshDrawFlags /*|DU_DRAWNAVMESH_COLOR_TILES*/)
		}
		if s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_BVTREE {
			debug_utils.DuDebugDrawNavMeshBVTree(s.m_dd, s.m_navMesh)
		}

		if s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_PORTALS {
			debug_utils.DuDebugDrawNavMeshPortals(s.m_dd, s.m_navMesh)
		}

		if s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_NODES {
			debug_utils.DuDebugDrawNavMeshNodes(s.m_dd, s.m_navQuery)
		}

		debug_utils.DuDebugDrawNavMeshPolysWithFlags(s.m_dd, s.m_navMesh, SAMPLE_POLYFLAGS_DISABLED, debug_utils.DuRGBA(0, 0, 0, 128))
	}

	gl.DepthMask(true)

	s.m_geom.drawConvexVolumes(s.m_dd)

	if s.m_tool != nil {
		s.m_tool.handleRender()
	}

	s.renderToolStates()

	gl.DepthMask(true)
}
func (s *SampleTempObstacles) handleRenderOverlay(proj, model []float64, view []int) {
	if s.m_tool != nil {
		s.m_tool.handleRenderOverlay(proj, model, view)
	}

	s.renderOverlayToolStates(proj, model, view)

	// Stats
	/*	imguiDrawRect(280,10,300,100,imguiRGBA(0,0,0,64));

		char text[64];
		int y = 110-30;

		snprintf(text,64,"Lean Data: %.1fkB", m_tileCache.getRawSize()/1024.0f);
		imguiDrawText(300, y, IMGUI_ALIGN_LEFT, text, imguiRGBA(255,255,255,255));
		y -= 20;

		snprintf(text,64,"Compressed: %.1fkB (%.1f%%)", m_tileCache.getCompressedSize()/1024.0f,
				 m_tileCache.getRawSize() > 0 ? 100.0f*(float)m_tileCache.getCompressedSize()/(float)m_tileCache.getRawSize() : 0);
		imguiDrawText(300, y, IMGUI_ALIGN_LEFT, text, imguiRGBA(255,255,255,255));
		y -= 20;

		if (m_rebuildTileCount > 0 && m_rebuildTime > 0.0f)
		{
			snprintf(text,64,"Changed obstacles, rebuild %d tiles: %.3f ms", m_rebuildTileCount, m_rebuildTime);
			imguiDrawText(300, y, IMGUI_ALIGN_LEFT, text, imguiRGBA(255,192,0,255));
			y -= 20;
		}
	*/
}
func (s *SampleTempObstacles) handleMeshChanged(geom *InputGeom) {
	s.Sample.handleMeshChanged(geom)
	s.m_tileCache = nil

	s.m_navMesh = nil

	if s.m_tool != nil {
		s.m_tool.reset()
		s.m_tool.init(s.Sample)
		s.m_tmproc = newMeshProcess(s.m_geom)
	}
	s.resetToolStates()
	s.initToolStates(s.Sample)
}
func (s *SampleTempObstacles) handleBuild() bool {

	if s.m_geom == nil || s.m_geom.getMesh() == nil {
		log.Printf("buildTiledNavigation: No vertices and triangles.")
		return false
	}

	s.m_tmproc = newMeshProcess(s.m_geom)

	// Init cache
	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	gw := 0
	gh := 0
	recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh)
	ts := int(s.m_tileSize)
	tw := (gw + ts - 1) / ts
	th := (gh + ts - 1) / ts

	// Generation params.
	var cfg recast.RcConfig
	cfg.Cs = s.m_cellSize
	cfg.Ch = s.m_cellHeight
	cfg.WalkableSlopeAngle = s.m_agentMaxSlope
	cfg.WalkableHeight = int(math.Ceil(s.m_agentHeight / cfg.Ch))
	cfg.WalkableClimb = int(math.Floor(s.m_agentMaxClimb / cfg.Ch))
	cfg.WalkableRadius = int(math.Ceil(s.m_agentRadius / cfg.Cs))
	cfg.MaxEdgeLen = int(s.m_edgeMaxLen / s.m_cellSize)
	cfg.MaxSimplificationError = s.m_edgeMaxError
	cfg.MinRegionArea = int(common.Sqr(s.m_regionMinSize))     // Note: area = size*size
	cfg.MergeRegionArea = int(common.Sqr(s.m_regionMergeSize)) // Note: area = size*size
	cfg.MaxVertsPerPoly = int(s.m_vertsPerPoly)
	cfg.TileSize = int(s.m_tileSize)
	cfg.BorderSize = cfg.WalkableRadius + 3 // Reserve enough padding.
	cfg.Width = cfg.TileSize + cfg.BorderSize*2
	cfg.Height = cfg.TileSize + cfg.BorderSize*2
	cfg.DetailSampleDist = s.m_cellSize * s.m_detailSampleDist
	if s.m_detailSampleDist < 0.9 {
		cfg.DetailSampleDist = 0
	}
	cfg.DetailSampleMaxError = s.m_cellHeight * s.m_detailSampleMaxError
	copy(cfg.Bmin[:], bmin)
	copy(cfg.Bmax[:], bmax)

	// Tile cache params.
	var tcparams recast.DtTileCacheParams
	copy(tcparams.Orig[:], bmin)
	tcparams.Cs = s.m_cellSize
	tcparams.Ch = s.m_cellHeight
	tcparams.Width = int(s.m_tileSize)
	tcparams.Height = int(s.m_tileSize)
	tcparams.WalkableHeight = s.m_agentHeight
	tcparams.WalkableRadius = s.m_agentRadius
	tcparams.WalkableClimb = s.m_agentMaxClimb
	tcparams.MaxSimplificationError = s.m_edgeMaxError
	tcparams.MaxTiles = tw * th * EXPECTED_LAYERS_PER_TILE
	tcparams.MaxObstacles = 128
	s.m_tileCache = &recast.DtTileCache{}
	if s.m_tileCache == nil {
		log.Printf("buildTiledNavigation: Could not allocate tile cache.")
		return false
	}
	status := s.m_tileCache.Init(&tcparams, s.m_tcomp, s.m_tmproc)
	if status.DtStatusFailed() {
		log.Printf("buildTiledNavigation: Could not init tile cache.")
		return false
	}

	var params recast.NavMeshParams
	copy(params.Orig[:], bmin)
	params.TileWidth = s.m_tileSize * s.m_cellSize
	params.TileHeight = s.m_tileSize * s.m_cellSize
	params.MaxTiles = s.m_maxTiles
	params.MaxPolys = s.m_maxPolysPerTile

	s.m_navMesh, status = recast.NewDtNavMeshWithParams(&params)
	if status.DtStatusFailed() {
		log.Printf("buildTiledNavigation: Could not init navmesh.")
		return false
	}

	s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048)
	// Preprocess tiles.
	s.m_cacheLayerCount = 0
	for y := 0; y < th; y++ {
		for x := 0; x < tw; x++ {
			tiles := s.rasterizeTileLayers(x, y, &cfg, MAX_LAYERS)

			for i := 0; i < len(tiles); i++ {
				tile := tiles[i]
				status = s.m_tileCache.AddTile(tile, recast.DT_COMPRESSEDTILE_FREE_DATA, nil)
				if status.DtStatusFailed() {
					continue
				}

				s.m_cacheLayerCount++
			}
		}
	}

	// Build initial meshes
	now := time.Now()
	for y := 0; y < th; y++ {
		for x := 0; x < tw; x++ {
			s.m_tileCache.BuildNavMeshTilesAt(x, y, s.m_navMesh)
		}
	}
	s.m_cacheBuildTimeMs = time.Since(now)
	if s.m_tool != nil {
		s.m_tool.init(s.Sample)
	}

	s.initToolStates(s.Sample)

	return true
}

func (s *SampleTempObstacles) handleUpdate(dt float64) {
	s.Sample.handleUpdate(dt)

	if s.m_navMesh == nil {
		return
	}

	if s.m_tileCache == nil {
		return
	}

	s.m_tileCache.Update(dt, s.m_navMesh)
}

func (s *SampleTempObstacles) getTilePos(pos []float64, tx, ty *int) {
	if s.m_geom == nil {
		return
	}

	bmin := s.m_geom.getNavMeshBoundsMin()

	ts := s.m_tileSize * s.m_cellSize
	*tx = int((pos[0] - bmin[0]) / ts)
	*ty = int((pos[2] - bmin[2]) / ts)
}

func (s *SampleTempObstacles) renderCachedTile(tx, ty, Type int) {
	if s.m_tileCache != nil {
		drawDetail(s.m_dd, s.m_tileCache, tx, ty, Type)
	}

}
func (s *SampleTempObstacles) renderCachedTileOverlay(tx, ty int, proj, model []float64, view []int) {
	if s.m_tileCache != nil {
		drawDetailOverlay(s.gs, s.m_tileCache, tx, ty, proj, model, view)
	}

}

func (s *SampleTempObstacles) addTempObstacle(pos []float64) {
	if s.m_tileCache == nil {
		return
	}

	p := make([]float64, 3)
	copy(p, pos)
	p[1] -= 0.5
	s.m_tileCache.AddObstacle(p, 1.0, 2.0, nil)
}
func (s *SampleTempObstacles) removeTempObstacle(sp, sq []float64) {
	if s.m_tileCache == nil {
		return
	}

	ref := hitTestObstacle(s.m_tileCache, sp, sq)
	s.m_tileCache.RemoveObstacle(ref)
}

func (s *SampleTempObstacles) clearAllTempObstacles() {
	if s.m_tileCache == nil {
		return
	}

	for i := 0; i < s.m_tileCache.GetObstacleCount(); i++ {
		ob := s.m_tileCache.GetObstacle(i)
		if ob.State == recast.DT_OBSTACLE_EMPTY {
			continue
		}
		s.m_tileCache.RemoveObstacle(s.m_tileCache.GetObstacleRef(ob))
	}
}

func (s *SampleTempObstacles) saveAll(p string) {}
func (s *SampleTempObstacles) loadAll(p string) {}

type RasterizationContext struct {
	solid    *recast.RcHeightfield
	triareas []int
	lset     *recast.RcHeightfieldLayerSet
	chf      *recast.RcCompactHeightfield
	tiles    []*recast.DetourTitleCacheLayerData
}

func newRasterizationContext() *RasterizationContext {
	return &RasterizationContext{
		tiles: make([]*recast.DetourTitleCacheLayerData, MAX_LAYERS),
	}
}

type FastLZCompressor struct {
}

func (FastLZCompressor) Compress(buffer []byte) ([]byte, error) {
	return buffer, nil
}
func (FastLZCompressor) Decompress(compressed []byte) ([]byte, error) {
	return compressed, nil
}
func (s *SampleTempObstacles) rasterizeTileLayers(tx, ty int, cfg *recast.RcConfig, maxTiles int) (datas []*recast.DetourTitleCacheLayerData) {
	if s.m_geom == nil || s.m_geom.getMesh() == nil || s.m_geom.getChunkyMesh() == nil {
		log.Printf("buildTile: Input mesh is not specified.")
		return nil
	}

	var comp FastLZCompressor
	rc := newRasterizationContext()

	verts := s.m_geom.getMesh().getVerts()
	nverts := s.m_geom.getMesh().getVertCount()
	chunkyMesh := s.m_geom.getChunkyMesh()

	// Tile bounds.
	tcs := float64(cfg.TileSize) * cfg.Cs

	var tcfg recast.RcConfig
	tcfg.Bmin[0] = cfg.Bmin[0] + float64(tx)*tcs
	tcfg.Bmin[1] = cfg.Bmin[1]
	tcfg.Bmin[2] = cfg.Bmin[2] + float64(ty)*tcs
	tcfg.Bmax[0] = cfg.Bmin[0] + float64(tx+1)*tcs
	tcfg.Bmax[1] = cfg.Bmax[1]
	tcfg.Bmax[2] = cfg.Bmin[2] + float64(ty+1)*tcs
	tcfg.Bmin[0] -= float64(tcfg.BorderSize) * tcfg.Cs
	tcfg.Bmin[2] -= float64(tcfg.BorderSize) * tcfg.Cs
	tcfg.Bmax[0] += float64(tcfg.BorderSize) * tcfg.Cs
	tcfg.Bmax[2] += float64(tcfg.BorderSize) * tcfg.Cs

	// Allocate voxel heightfield where we rasterize our input data to.
	rc.solid = recast.RcCreateHeightfield(tcfg.Width, tcfg.Height, tcfg.Bmin[:], tcfg.Bmax[:], tcfg.Cs, tcfg.Ch)
	// Allocate array that can hold triangle flags.
	// If you have multiple meshes you need to process, allocate
	// and array which can hold the max number of triangles you need to process.
	rc.triareas = make([]int, chunkyMesh.maxTrisPerChunk)
	var tbmin [2]float64
	var tbmax [2]float64
	tbmin[0] = tcfg.Bmin[0]
	tbmin[1] = tcfg.Bmin[2]
	tbmax[0] = tcfg.Bmax[0]
	tbmax[1] = tcfg.Bmax[2]
	cid := make([]int, 512) // TODO: Make grow when returning too many items.
	ncid := rcGetChunksOverlappingRect(chunkyMesh, tbmin, tbmax, cid, 512)
	if ncid == 0 {
		return // empty
	}

	for i := 0; i < ncid; i++ {
		node := chunkyMesh.nodes[cid[i]]
		tris := chunkyMesh.tris[node.i*3:]
		ntris := node.n
		recast.RcMarkWalkableTriangles(tcfg.WalkableSlopeAngle,
			verts, nverts, tris, ntris, rc.triareas)

		if !recast.RcRasterizeTriangles(verts, nverts, tris, rc.triareas, ntris, rc.solid, tcfg.WalkableClimb) {
			return
		}

	}

	// Once all geometry is rasterized, we do initial pass of filtering to
	// remove unwanted overhangs caused by the conservative rasterization
	// as well as filter spans where the character cannot possibly stand.
	if s.m_filterLowHangingObstacles {
		recast.RcFilterLowHangingWalkableObstacles(tcfg.WalkableClimb, rc.solid)
	}

	if s.m_filterLedgeSpans {
		recast.RcFilterLedgeSpans(tcfg.WalkableHeight, tcfg.WalkableClimb, rc.solid)
	}

	if s.m_filterWalkableLowHeightSpans {
		recast.RcFilterWalkableLowHeightSpans(tcfg.WalkableHeight, rc.solid)
	}

	rc.chf = &recast.RcCompactHeightfield{}
	if !recast.RcBuildCompactHeightfield(tcfg.WalkableHeight, tcfg.WalkableClimb, rc.solid, rc.chf) {
		log.Printf("buildNavigation: Could not build compact data.")
		return
	}

	// Erode the walkable area by agent radius.
	if !recast.RcErodeWalkableArea(tcfg.WalkableRadius, rc.chf) {
		log.Printf("buildNavigation: Could not erode.")
		return
	}

	// (Optional) Mark areas.
	vols := s.m_geom.getConvexVolumes()
	for i := 0; i < s.m_geom.getConvexVolumeCount(); i++ {
		recast.RcMarkConvexPolyArea(vols[i].verts, vols[i].nverts,
			vols[i].hmin, vols[i].hmax,
			vols[i].area, rc.chf)
	}

	rc.lset = &recast.RcHeightfieldLayerSet{}
	if !recast.RcBuildHeightfieldLayers(rc.chf, tcfg.BorderSize, tcfg.WalkableHeight, rc.lset) {
		log.Printf("buildNavigation: Could not build heighfield layers.")
		return
	}

	for i := 0; i < common.Min(rc.lset.Nlayers, MAX_LAYERS); i++ {
		layer := rc.lset.Layers[i]

		// Store header
		var header recast.DtTileCacheLayerHeader
		header.Magic = recast.DT_TILECACHE_MAGIC
		header.Version = recast.DT_TILECACHE_VERSION

		// Tile layer location in the navmesh.
		header.Tx = tx
		header.Ty = ty
		header.Tlayer = i
		header.Bmin = layer.Bmin
		header.Bmax = layer.Bmax

		// Tile info.
		header.Width = layer.Width
		header.Height = layer.Height
		header.Minx = layer.Minx
		header.Maxx = layer.Maxx
		header.Miny = layer.Miny
		header.Maxy = layer.Maxy
		header.Hmin = layer.Hmin
		header.Hmax = layer.Hmax
		rc.tiles = append(rc.tiles, recast.DtBuildTileCacheLayer(comp, &header, layer.Heights, layer.Areas, layer.Cons))
	}

	// Transfer ownsership of tile data from build context to the caller.

	for i := 0; i < common.Min(len(rc.tiles), maxTiles); i++ {

		datas = append(datas, rc.tiles[i])
		rc.tiles[i] = &recast.DetourTitleCacheLayerData{}
	}

	return datas
}

type MeshProcess struct {
	m_geom *InputGeom
}

func newMeshProcess(m_geom *InputGeom) *MeshProcess {
	return &MeshProcess{m_geom: m_geom}
}
func (p *MeshProcess) Process(params *recast.DtNavMeshCreateParams, polyAreas []int, polyFlags []int) {
	// Update poly flags from areas.
	for i := 0; i < params.PolyCount; i++ {
		if polyAreas[i] == recast.DT_TILECACHE_WALKABLE_AREA {
			polyAreas[i] = SAMPLE_POLYAREA_GROUND

			if polyAreas[i] == SAMPLE_POLYAREA_GROUND ||
				polyAreas[i] == SAMPLE_POLYAREA_GRASS ||
				polyAreas[i] == SAMPLE_POLYAREA_ROAD {
				polyFlags[i] = SAMPLE_POLYFLAGS_WALK
			} else if polyAreas[i] == SAMPLE_POLYAREA_WATER {
				polyFlags[i] = SAMPLE_POLYFLAGS_SWIM
			} else if polyAreas[i] == SAMPLE_POLYAREA_DOOR {
				polyFlags[i] = SAMPLE_POLYFLAGS_WALK | SAMPLE_POLYFLAGS_DOOR
			}
		}

		// Pass in off-mesh connections.
		if p.m_geom != nil {
			params.OffMeshConVerts = p.m_geom.getOffMeshConnectionVerts()
			params.OffMeshConRad = p.m_geom.getOffMeshConnectionRads()
			params.OffMeshConDir = p.m_geom.getOffMeshConnectionDirs()
			params.OffMeshConAreas = p.m_geom.getOffMeshConnectionAreas()
			params.OffMeshConFlags = p.m_geom.getOffMeshConnectionFlags()
			params.OffMeshConUserID = p.m_geom.getOffMeshConnectionId()
			params.OffMeshConCount = p.m_geom.getOffMeshConnectionCount()
		}
	}
}

const MAX_LAYERS = 32

type TempObstacleHilightTool struct {
	m_sample    *SampleTempObstacles
	m_hitPos    []float64
	m_hitPosSet bool
	m_drawType  int
	gs          *guiState
}

func newTempObstacleHilightTool(gs *guiState) *TempObstacleHilightTool {
	return &TempObstacleHilightTool{
		gs:         gs,
		m_hitPos:   make([]float64, 3),
		m_drawType: DRAWDETAIL_AREAS,
		m_sample:   newSampleTempObstacles(),
	}
}
func (s *TempObstacleHilightTool) Type() int { return TOOL_TILE_HIGHLIGHT }

func (s *TempObstacleHilightTool) init(sample *Sample) {
	s.m_sample.Sample = sample
}

func (s *TempObstacleHilightTool) reset() {}

func (s *TempObstacleHilightTool) handleMenu() {
	s.gs.imguiLabel("Highlight Tile Cache")
	s.gs.imguiValue("Click LMB to highlight a tile.")
	s.gs.imguiSeparator()
	if s.gs.imguiCheck("Draw Areas", s.m_drawType == DRAWDETAIL_AREAS) {
		s.m_drawType = DRAWDETAIL_AREAS
	}

	if s.gs.imguiCheck("Draw Regions", s.m_drawType == DRAWDETAIL_REGIONS) {
		s.m_drawType = DRAWDETAIL_REGIONS
	}

	if s.gs.imguiCheck("Draw Contours", s.m_drawType == DRAWDETAIL_CONTOURS) {
		s.m_drawType = DRAWDETAIL_CONTOURS
	}

	if s.gs.imguiCheck("Draw Mesh", s.m_drawType == DRAWDETAIL_MESH) {
		s.m_drawType = DRAWDETAIL_MESH
	}

}

func (s *TempObstacleHilightTool) handleClick(ss, p []float64, shift bool) {
	s.m_hitPosSet = true
	copy(s.m_hitPos[:], p)
}

func (s *TempObstacleHilightTool) handleToggle() {}

func (s *TempObstacleHilightTool) handleStep() {}

func (s *TempObstacleHilightTool) handleUpdate(dt float64) {}

func (s *TempObstacleHilightTool) handleRender() {
	if s.m_hitPosSet && s.m_sample != nil {
		ss := s.m_sample.getAgentRadius()
		gl.Color4ub(0, 0, 0, 128)
		gl.LineWidth(2.0)
		gl.vertext(gl.LINES)
		gl.Vertex3f(s.m_hitPos[0]-s.m_hitPos[1]+0.1, s.m_hitPos[2])
		glVertex3f(s.m_hitPos[0]+s.m_hitPos[1]+0.1, s.m_hitPos[2])
		glVertex3f(s.m_hitPos[0], s.m_hitPos[1]-ss+0.1, s.m_hitPos[2])
		glVertex3f(s.m_hitPos[0], s.m_hitPos[1]+ss+0.1, s.m_hitPos[2])
		glVertex3f(s.m_hitPos[0], s.m_hitPos[1]+0.1, s.m_hitPos[2]-ss)
		glVertex3f(s.m_hitPos[0], s.m_hitPos[1]+0.1, s.m_hitPos[2]+ss)
		gl.End()
		gl.LineWidth(1.0)

		tx := 0
		ty := 0
		s.m_sample.getTilePos(s.m_hitPos, &tx, &ty)
		s.m_sample.renderCachedTile(tx, ty, s.m_drawType)
	}
}

func (s *TempObstacleHilightTool) handleRenderOverlay(proj, model []float64, view []int) {
	if s.m_hitPosSet {
		if s.m_sample != nil {
			tx := 0
			ty := 0
			s.m_sample.getTilePos(s.m_hitPos, &tx, &ty)
			s.m_sample.renderCachedTileOverlay(tx, ty, proj, model, view)
		}
	}
}

type TempObstacleCreateTool struct {
	m_sample *SampleTempObstacles
	gs       *guiState
}

func newTempObstacleCreateTool(gs *guiState) *TempObstacleCreateTool {
	return &TempObstacleCreateTool{gs: gs, m_sample: newSampleTempObstacles()}
}
func (s *TempObstacleCreateTool) Type() int { return TOOL_TEMP_OBSTACLE }

func (s *TempObstacleCreateTool) init(sample *Sample) {
	s.m_sample.Sample = sample
}

func (s *TempObstacleCreateTool) reset() {}

func (s *TempObstacleCreateTool) handleMenu() {
	s.gs.imguiLabel("Create Temp Obstacles")

	if s.gs.imguiButton("Remove All") {
		s.m_sample.clearAllTempObstacles()
	}
	s.gs.imguiSeparator()
	s.gs.imguiValue("Click LMB to create an obstacle.")
	s.gs.imguiValue("Shift+LMB to remove an obstacle.")
}

func (s *TempObstacleCreateTool) handleClick(ss, p []float64, shift bool) {
	if s.m_sample != nil {
		if shift {
			s.m_sample.removeTempObstacle(ss, p)
		} else {
			s.m_sample.addTempObstacle(p)
		}

	}
}

func (s *TempObstacleCreateTool) handleToggle()                                         {}
func (s *TempObstacleCreateTool) handleStep()                                           {}
func (s *TempObstacleCreateTool) handleUpdate(dt float64)                               {}
func (s *TempObstacleCreateTool) handleRender()                                         {}
func (s *TempObstacleCreateTool) handleRenderOverlay(proj, model []float64, view []int) {}
