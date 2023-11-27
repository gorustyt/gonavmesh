package gui

import (
	"fmt"
	"gonavamesh/common"
	"gonavamesh/debug_utils"
	"gonavamesh/recast"
	"log"
	"math"
	"time"
)

// This value specifies how many layers (or "floors") each navmesh tile is expected to have.
 const  EXPECTED_LAYERS_PER_TILE = 4;


func  drawTiles(dd debug_utils.DuDebugDraw ,tc *recast.DtTileCache) {
 fcol:=make([]int,6)
bmin:=make([]float64,3)
bmax:=make([]float64,3)

for  i := 0; i < tc.GetTileCount(); i++{
 tile := tc.GetTile(i);
if (tile.Header==nil) {continue;}

tc.CalcTightTileBounds(tile.Header, bmin, bmax);

col := debug_utils.DuIntToCol(i,64);
	debug_utils.DuCalcBoxColors(fcol, col, col);
	debug_utils.DuDebugDrawBox(dd, bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2], fcol);
}

for  i := 0; i < tc.GetTileCount(); i++{
tile := tc.GetTile(i);
if (tile.Header==nil){ continue;}

tc.CalcTightTileBounds(tile.Header, bmin, bmax);

 col := debug_utils.DuIntToCol(i,255);
 pad := tc.GetParams().Cs * 0.1;
	debug_utils.DuDebugDrawBoxWire(dd, bmin[0]-pad,bmin[1]-pad,bmin[2]-pad,
bmax[0]+pad,bmax[1]+pad,bmax[2]+pad, col, 2.0);
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
func  drawDetailOverlay(gs *guiState,tc *recast. DtTileCache ,   tx,  ty int ,  proj, model []float64, view []int ){
tiles:=make([]recast.DtCompressedTileRef ,MAX_LAYERS);
 ntiles := tc.GetTilesAt(tx,ty,tiles,MAX_LAYERS);
if (ntiles==0){return;}


rawSize := calcLayerBufferSize(tc.GetParams().Width, tc.GetParams().Height);



for  i := 0; i < ntiles; i++{
 tile := tc.GetTileByRef(tiles[i]);

		pos:=make([]float64,3)
pos[0] = (tile.Header.Bmin[0]+tile.Header.Bmax[0])/2.0;
pos[1] = tile.Header.Bmin[1];
pos[2] = (tile.Header.Bmin[2]+tile.Header.Bmax[2])/2.0;

res:=common.GluProject([]float64{pos[0], pos[1], pos[2]}, model, proj, view,)

if (len(res)!=0){
	x,y:=int(res[0]),int(res[1])
	text:=fmt.Sprintf("(%d,%d)/%d", tile.Header.Tx,tile.Header.Ty,tile.Header.Tlayer);
	gs.imguiDrawText(x, y-25, IMGUI_ALIGN_CENTER, text, imguiRGBA(0,0,0,220));
	text=fmt.Sprintf("Compressed: %.1f kB", tile.DataSize/1024.0);
	gs.imguiDrawText(x, y-45, IMGUI_ALIGN_CENTER, text, imguiRGBA(0,0,0,128));
	text=fmt.Sprintf("Raw:%.1fkB", rawSize/1024.0);
	gs.imguiDrawText(x, y-65, IMGUI_ALIGN_CENTER, text, imguiRGBA(0,0,0,128));
}
}
}

func drawObstacles(dd debug_utils.DuDebugDraw ,tc *recast.DtTileCache ) {
// Draw obstacles
for  i := 0; i < tc.GetObstacleCount(); i++{
ob := tc.GetObstacle(i);
if (ob.State == recast.DT_OBSTACLE_EMPTY) {continue;}
bmin:=make([]float64,3)
bmax:=make([]float64,3)
tc.GetObstacleBounds(ob, bmin,bmax);

var  col int
if (ob.State == recast.DT_OBSTACLE_PROCESSING){
	col = debug_utils.DuRGBA(255,255,0,128);} else if (ob.State ==recast.DT_OBSTACLE_PROCESSED){
		col = debug_utils.DuRGBA(255,192,0,192);} else if (ob.State == recast.DT_OBSTACLE_REMOVING){
			col = debug_utils.DuRGBA(220,0,0,128);}


	debug_utils.DuDebugDrawCylinder(dd, bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2], col);
	debug_utils.DuDebugDrawCylinderWire(dd, bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2], debug_utils.DuDarkenCol(col), 2);
}
}

const(
	DRAWDETAIL_AREAS=iota
	DRAWDETAIL_REGIONS
	DRAWDETAIL_CONTOURS
	DRAWDETAIL_MESH
)
 type  TileCacheBuildContext struct{
	 layer *recast.DtTileCacheLayer ;
	 lcset *recast.DtTileCacheContourSet ;
	 lmesh *recast.DtTileCachePolyMesh ;
};

func newTileCacheBuildContext()*TileCacheBuildContext  {
	return &TileCacheBuildContext{
		layer: &recast.DtTileCacheLayer{},
		lcset: &recast.DtTileCacheContourSet{},
		lmesh:&recast.DtTileCachePolyMesh{},
	}
}
func (t*TileCacheBuildContext)purge() {
	t.layer = nil
	t.lcset = nil
	t.lmesh =nil
}
func drawDetail(dd debug_utils.DuDebugDraw , tc *recast.DtTileCache , tx,  ty, tType int ){

tiles:=make([]recast.DtCompressedTileRef ,MAX_LAYERS);
ntiles := tc.GetTilesAt(tx,ty,tiles,MAX_LAYERS);

 tcomp := tc.GetCompressor();
params := tc.GetParams();

for  i := 0; i < ntiles; i++{
tile := tc.GetTileByRef(tiles[i]);
bc :=newTileCacheBuildContext()
walkableClimbVx := int(params.WalkableClimb / params.Ch);
var status recast.DtStatus
// Decompress tile layer data.
	bc.layer,status  = recast.DtDecompressTileCacheLayer(tcomp, tile.Data, tile.DataSize, );
if (status.DtStatusFailed()){return;}

if (tType == DRAWDETAIL_AREAS){
	debug_utils.DuDebugDrawTileCacheLayerAreas(dd, bc.layer, params.Cs, params.Ch);
continue;
}

// Build navmesh
status = recast.DtBuildTileCacheRegions( bc.layer, walkableClimbVx);
if (status.DtStatusFailed()){return;}

if (tType == DRAWDETAIL_REGIONS){
	debug_utils.DuDebugDrawTileCacheLayerRegions(dd, bc.layer,  params.Cs, params.Ch);
continue;
}

bc.lcset = &recast.DtTileCacheContourSet{}


status = recast.DtBuildTileCacheContours( bc.layer, walkableClimbVx,
params.MaxSimplificationError, bc.lcset);
if (status.DtStatusFailed()){return;}

if (tType == DRAWDETAIL_CONTOURS) {
	debug_utils.DuDebugDrawTileCacheContours(dd, bc.lcset, tile.Header.Bmin[:],  params.Cs, params.Ch);
continue;
}

bc.lmesh = &recast.DtTileCachePolyMesh{}
status = recast.DtBuildTileCachePolyMesh( bc.lcset, bc.lmesh);
if (status.DtStatusFailed()){return;}


if (tType== DRAWDETAIL_MESH){
	debug_utils.DuDebugDrawTileCachePolyMesh(dd, bc.lmesh, tile.Header.Bmin[:], params.Cs, params.Ch);
continue;
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

	m_talloc *LinearAllocator
	m_tcomp  *FastLZCompressor
	m_tmproc *MeshProcess

	m_tileCache *recast.DtTileCache

	m_cacheBuildTimeMs   time.Duration
	m_cacheCompressedSize int
	m_cacheRawSize        int
	m_cacheLayerCount     int
	m_cacheBuildMemUsage  int

	m_drawMode SampleTempObstacleDrawMode

	m_maxTile         int
	m_maxPolysPerTile int
	m_tileSize        float64
}

func newSampleTempObstacles() *SampleTempObstacles {
	s := &SampleTempObstacles{}
	return s
}

func (s *SampleTempObstacles) handleSettings()                                       {
Sample::handleCommonSettings();

	if (imguiCheck("Keep Itermediate Results", m_keepInterResults))
		m_keepInterResults = !m_keepInterResults;

	imguiLabel("Tiling");
	imguiSlider("TileSize", &m_tileSize, 16.0f, 128.0f, 8.0f);

	int gridSize = 1;
	if (m_geom)
	{
		const float* bmin = m_geom->getNavMeshBoundsMin();
		const float* bmax = m_geom->getNavMeshBoundsMax();
		char text[64];
		int gw = 0, gh = 0;
		rcCalcGridSize(bmin, bmax, m_cellSize, &gw, &gh);
		const int ts = (int)m_tileSize;
		const int tw = (gw + ts-1) / ts;
		const int th = (gh + ts-1) / ts;
		snprintf(text, 64, "Tiles  %d x %d", tw, th);
		imguiValue(text);

		// Max tiles and max polys affect how the tile IDs are caculated.
		// There are 22 bits available for identifying a tile and a polygon.
		int tileBits = rcMin((int)dtIlog2(dtNextPow2(tw*th*EXPECTED_LAYERS_PER_TILE)), 14);
		if (tileBits > 14) tileBits = 14;
		int polyBits = 22 - tileBits;
		m_maxTiles = 1 << tileBits;
		m_maxPolysPerTile = 1 << polyBits;
		snprintf(text, 64, "Max Tiles  %d", m_maxTiles);
		imguiValue(text);
		snprintf(text, 64, "Max Polys  %d", m_maxPolysPerTile);
		imguiValue(text);
		gridSize = tw*th;
	}
	else
	{
		m_maxTiles = 0;
		m_maxPolysPerTile = 0;
	}

	imguiSeparator();

	imguiLabel("Tile Cache");
	char msg[64];

	const float compressionRatio = (float)m_cacheCompressedSize / (float)(m_cacheRawSize+1);

	snprintf(msg, 64, "Layers  %d", m_cacheLayerCount);
	imguiValue(msg);
	snprintf(msg, 64, "Layers (per tile)  %.1f", (float)m_cacheLayerCount/(float)gridSize);
	imguiValue(msg);

	snprintf(msg, 64, "Memory  %.1f kB / %.1f kB (%.1f%%)", m_cacheCompressedSize/1024.0f, m_cacheRawSize/1024.0f, compressionRatio*100.0f);
	imguiValue(msg);
	snprintf(msg, 64, "Navmesh Build Time  %.1f ms", m_cacheBuildTimeMs);
	imguiValue(msg);
	snprintf(msg, 64, "Build Peak Mem Usage  %.1f kB", m_cacheBuildMemUsage/1024.0f);
	imguiValue(msg);

	imguiSeparator();

	imguiIndent();
	imguiIndent();

	if (imguiButton("Save"))
	{
		saveAll("all_tiles_tilecache.bin");
	}

	if (imguiButton("Load"))
	{
		dtFreeNavMesh(m_navMesh);
		dtFreeTileCache(m_tileCache);
		loadAll("all_tiles_tilecache.bin");
		m_navQuery->init(m_navMesh, 2048);
	}

	imguiUnindent();
	imguiUnindent();

	imguiSeparator();
}
func (s *SampleTempObstacles) handleTools()                                          {
	 Type :=  s.m_tool.Type();
	 if s.m_tool==nil  {
		 Type=TOOL_NONE
	}

	if (s.gs.imguiCheck("Test Navmesh", Type == TOOL_NAVMESH_TESTER)){
		s.setTool(newNavMeshTesterTool(s.gs));
	}
	if (s.gs.imguiCheck("Highlight Tile Cache", Type == TOOL_TILE_HIGHLIGHT)){
		s.setTool(newTempObstacleHilightTool);
	}
	if (s.gs.imguiCheck("Create Temp Obstacles", Type == TOOL_TEMP_OBSTACLE)){
	s.setTool(newTempObstacleCreateTool);
	}
	if (s.gs.imguiCheck("Create Off-Mesh Links", Type == TOOL_OFFMESH_CONNECTION)){
		s.setTool(newOffMeshConnectionTool(s.gs));
	}
	if (s.gs.imguiCheck("Create Convex Volumes",Type == TOOL_CONVEX_VOLUME)) {
		s.setTool(newConvexVolumeTool(s.gs));
	}
	if (s.gs.imguiCheck("Create Crowds",Type == TOOL_CROWD)) {
		s.setTool(newCrowdTool(s.gs));
	}

	s.gs.imguiSeparatorLine();

	s.gs.imguiIndent();

	if (s.m_tool!=nil){
		s.m_tool.handleMenu();
	}


	s.gs.imguiUnindent();
}
func (s *SampleTempObstacles) handleDebugMode()                                      {
	// Check which modes are valid.
	valid:=make([]bool,MAX_DRAWMODE);
	for i := 0; i < MAX_DRAWMODE;i++{
		valid[i] = false;
	}


	if (s.m_geom!=nil) {
		valid[SampleTempObstacleDRAWMODE_NAVMESH] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_TRANS] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_BVTREE] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_NODES] = s.m_navQuery != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_PORTALS] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_NAVMESH_INVIS] = s.m_navMesh != nil
		valid[SampleTempObstacleDRAWMODE_MESH] = true;
		valid[SampleTempObstacleDRAWMODE_CACHE_BOUNDS] = true;
	}

	 unavail := 0;
	for  i := 0; i < int(SampleTempObstacleMAX_DRAWMODE); i++{
		if (!valid[i]){ unavail++;}
	}


	if (unavail == int(SampleTempObstacleMAX_DRAWMODE)){return;}


	s.gs.imguiLabel("Draw");
	if (s.gs.imguiCheck("Input Mesh", s.m_drawMode == SampleTempObstacleDRAWMODE_MESH, valid[SampleTempObstacleDRAWMODE_MESH])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_MESH;
	}

	if (s.gs.imguiCheck("Navmesh", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH, valid[SampleTempObstacleDRAWMODE_NAVMESH])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH;
	}

	if (s.gs.imguiCheck("Navmesh Invis", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_INVIS, valid[SampleTempObstacleDRAWMODE_NAVMESH_INVIS])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_INVIS;
	}

	if (s.gs.imguiCheck("Navmesh Trans", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_TRANS, valid[SampleTempObstacleDRAWMODE_NAVMESH_TRANS])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_TRANS;
	}

	if (s.gs.imguiCheck("Navmesh BVTree", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_BVTREE, valid[SampleTempObstacleDRAWMODE_NAVMESH_BVTREE])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_BVTREE;
	}

	if (s.gs.imguiCheck("Navmesh Nodes", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_NODES, valid[SampleTempObstacleDRAWMODE_NAVMESH_NODES])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_NAVMESH_NODES;
	}

	if (s.gs.imguiCheck("Navmesh Portals", s.m_drawMode == SampleTempObstacleDRAWMODE_NAVMESH_PORTALS, valid[SampleTempObstacleDRAWMODE_NAVMESH_PORTALS])){
		s.m_drawMode =SampleTempObstacleDRAWMODE_NAVMESH_PORTALS;
	}

	if (s.gs.imguiCheck("Cache Bounds", s.m_drawMode == SampleTempObstacleDRAWMODE_CACHE_BOUNDS, valid[SampleTempObstacleDRAWMODE_CACHE_BOUNDS])){
		s.m_drawMode = SampleTempObstacleDRAWMODE_CACHE_BOUNDS;
	}


	if (unavail>0) {
		s.gs.imguiValue("Tick 'Keep Itermediate Results'");
		s.gs.imguiValue("rebuild some tiles to see");
		s.gs.imguiValue("more debug mode options.");
	}
}
func (s *SampleTempObstacles) handleRender()                                         {
	if (!m_geom || !m_geom->getMesh())
	return;

	const float texScale = 1.0f / (m_cellSize * 10.0f);

	// Draw mesh
	if (m_drawMode != DRAWMODE_NAVMESH_TRANS)
	{
		// Draw mesh
		duDebugDrawTriMeshSlope(&m_dd, m_geom->getMesh()->getVerts(), m_geom->getMesh()->getVertCount(),
			m_geom->getMesh()->getTris(), m_geom->getMesh()->getNormals(), m_geom->getMesh()->getTriCount(),
			m_agentMaxSlope, texScale);
		m_geom->drawOffMeshConnections(&m_dd);
	}

	if (m_tileCache && m_drawMode == DRAWMODE_CACHE_BOUNDS)
		drawTiles(&m_dd, m_tileCache);

	if (m_tileCache)
		drawObstacles(&m_dd, m_tileCache);


	glDepthMask(GL_FALSE);

	// Draw bounds
	const float* bmin = m_geom->getNavMeshBoundsMin();
	const float* bmax = m_geom->getNavMeshBoundsMax();
	duDebugDrawBoxWire(&m_dd, bmin[0],bmin[1],bmin[2], bmax[0],bmax[1],bmax[2], duRGBA(255,255,255,128), 1.0f);

	// Tiling grid.
	int gw = 0, gh = 0;
	rcCalcGridSize(bmin, bmax, m_cellSize, &gw, &gh);
	const int tw = (gw + (int)m_tileSize-1) / (int)m_tileSize;
	const int th = (gh + (int)m_tileSize-1) / (int)m_tileSize;
	const float s = m_tileSize*m_cellSize;
	duDebugDrawGridXZ(&m_dd, bmin[0],bmin[1],bmin[2], tw,th, s, duRGBA(0,0,0,64), 1.0f);

	if (m_navMesh && m_navQuery &&
		(m_drawMode == DRAWMODE_NAVMESH ||
			m_drawMode == DRAWMODE_NAVMESH_TRANS ||
			m_drawMode == DRAWMODE_NAVMESH_BVTREE ||
			m_drawMode == DRAWMODE_NAVMESH_NODES ||
			m_drawMode == DRAWMODE_NAVMESH_PORTALS ||
			m_drawMode == DRAWMODE_NAVMESH_INVIS))
	{
		if (m_drawMode != DRAWMODE_NAVMESH_INVIS)
			duDebugDrawNavMeshWithClosedList(&m_dd, *m_navMesh, *m_navQuery, m_navMeshDrawFlags/*|DU_DRAWNAVMESH_COLOR_TILES*/);
		if (m_drawMode == DRAWMODE_NAVMESH_BVTREE)
			duDebugDrawNavMeshBVTree(&m_dd, *m_navMesh);
		if (m_drawMode == DRAWMODE_NAVMESH_PORTALS)
			duDebugDrawNavMeshPortals(&m_dd, *m_navMesh);
		if (m_drawMode == DRAWMODE_NAVMESH_NODES)
			duDebugDrawNavMeshNodes(&m_dd, *m_navQuery);
		duDebugDrawNavMeshPolysWithFlags(&m_dd, *m_navMesh, SAMPLE_POLYFLAGS_DISABLED, duRGBA(0,0,0,128));
	}


	glDepthMask(GL_TRUE);

	m_geom->drawConvexVolumes(&m_dd);

	if (m_tool)
		m_tool->handleRender();
	renderToolStates();

	glDepthMask(GL_TRUE);
}
func (s *SampleTempObstacles) handleRenderOverlay(proj, model []float64, view []int) {
	if (s.m_tool!=nil){s.m_tool.handleRenderOverlay(proj, model, view);}

	s.renderOverlayToolStates(proj, model, view);

	// Stats
	/*	imguiDrawRect(280,10,300,100,imguiRGBA(0,0,0,64));

		char text[64];
		int y = 110-30;

		snprintf(text,64,"Lean Data: %.1fkB", m_tileCache->getRawSize()/1024.0f);
		imguiDrawText(300, y, IMGUI_ALIGN_LEFT, text, imguiRGBA(255,255,255,255));
		y -= 20;

		snprintf(text,64,"Compressed: %.1fkB (%.1f%%)", m_tileCache->getCompressedSize()/1024.0f,
				 m_tileCache->getRawSize() > 0 ? 100.0f*(float)m_tileCache->getCompressedSize()/(float)m_tileCache->getRawSize() : 0);
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
func (s *SampleTempObstacles) handleMeshChanged(geom *InputGeom)                     {
s.Sample.handleMeshChanged(geom);
s.m_tileCache = nil


	s.m_navMesh = nil

	if (s.m_tool!=nil) {
		s.m_tool.reset();
		s.m_tool.init(s.Sample);
		s.m_tmproc=newMeshProcess(s.m_geom);
	}
	s.resetToolStates();
	s.initToolStates(s.Sample);
}
func (s *SampleTempObstacles) handleBuild() bool                                     {


	if (s.m_geom==nil || s.m_geom.getMesh()==nil) {
		log.Printf( "buildTiledNavigation: No vertices and triangles.");
		return false;
	}

	s.m_tmproc=newMeshProcess(s.m_geom)

	// Init cache
	 bmin := s.m_geom.getNavMeshBoundsMin();
	 bmax := s.m_geom.getNavMeshBoundsMax();
	gw := 0
	gh := 0;
	recast.RcCalcGridSize(bmin, bmax, s.m_cellSize, &gw, &gh);
	 ts := int(s.m_tileSize);
	 tw := (gw + ts-1) / ts;
	th := (gh + ts-1) / ts;

	// Generation params.
	var cfg recast.RcConfig ;
	cfg.Cs = s.m_cellSize;
	cfg.Ch = s.m_cellHeight;
	cfg.WalkableSlopeAngle = s.m_agentMaxSlope;
	cfg.WalkableHeight = int(math.Ceil(s.m_agentHeight / cfg.Ch));
	cfg.WalkableClimb = int(math.Floor(s.m_agentMaxClimb / cfg.Ch));
	cfg.WalkableRadius = int(math.Ceil(s.m_agentRadius / cfg.Cs));
	cfg.MaxEdgeLen = int(s.m_edgeMaxLen / s.m_cellSize);
	cfg.MaxSimplificationError = s.m_edgeMaxError;
	cfg.MinRegionArea = int(common.Sqr(s.m_regionMinSize));		// Note: area = size*size
	cfg.MergeRegionArea = int(common.Sqr(s.m_regionMergeSize));	// Note: area = size*size
	cfg.MaxVertsPerPoly = int(s.m_vertsPerPoly);
	cfg.TileSize = int(s.m_tileSize);
	cfg.BorderSize = cfg.WalkableRadius + 3; // Reserve enough padding.
	cfg.Width = cfg.TileSize + cfg.BorderSize*2;
	cfg.Height = cfg.TileSize + cfg.BorderSize*2;
	cfg.DetailSampleDist =  s.m_cellSize * s.m_detailSampleDist;
	if  s.m_detailSampleDist < 0.9 {
		cfg.DetailSampleDist = 0
	}
	cfg.DetailSampleMaxError = s.m_cellHeight * s.m_detailSampleMaxError;
	copy(cfg.Bmin[:], bmin);
	copy(cfg.Bmax[:], bmax);

	// Tile cache params.
	var tcparams recast.DtTileCacheParams ;
	copy(tcparams.Orig[:], bmin);
	tcparams.Cs = s.m_cellSize;
	tcparams.Ch = s.m_cellHeight;
	tcparams.Width =int(s.m_tileSize);
	tcparams.Height = int(s.m_tileSize);
	tcparams.WalkableHeight = s.m_agentHeight;
	tcparams.WalkableRadius = s.m_agentRadius;
	tcparams.WalkableClimb = s.m_agentMaxClimb;
	tcparams.MaxSimplificationError = s.m_edgeMaxError;
	tcparams.MaxTiles = tw*th*EXPECTED_LAYERS_PER_TILE;
	tcparams.MaxObstacles = 128;

	dtFreeTileCache(m_tileCache);

	m_tileCache = dtAllocTileCache();
	if (s.m_tileCache==nil) {
		log.Printf(  "buildTiledNavigation: Could not allocate tile cache.");
		return false;
	}
	status := s.m_tileCache.Init(&tcparams, m_talloc, m_tcomp, m_tmproc);
	if (status.DtStatusFailed()) {
		log.Printf( "buildTiledNavigation: Could not init tile cache.");
		return false;
	}

	var params recast.NavMeshParams ;
	copy(params.Orig[:], bmin);
	params.TileWidth = s.m_tileSize*s.m_cellSize;
	params.TileHeight = s.m_tileSize*s.m_cellSize;
	params.MaxTiles = s.m_maxTiles;
	params.MaxPolys = s.m_maxPolysPerTile;

	s.m_navMesh,status = recast.NewDtNavMeshWithParams(&params);
	if (status.DtStatusFailed()) {
		log.Printf( "buildTiledNavigation: Could not init navmesh.");
		return false;
	}

		s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048);
	// Preprocess tiles.
	s.m_cacheLayerCount = 0;
		s.m_cacheCompressedSize = 0;
		s.m_cacheRawSize = 0;

	for y := 0; y < th; y++{
	for  x := 0; x < tw; x++{
	 tiles:=make([]TileCacheData,MAX_LAYERS);
	memset(tiles, 0, sizeof(tiles));
	ntiles := rasterizeTileLayers(x, y, cfg, tiles, MAX_LAYERS);

	for i := 0; i < ntiles; i++{
	 tile := tiles[i];
	status = s.m_tileCache.AddTile(tile.data, tile.dataSize, recast.DT_COMPRESSEDTILE_FREE_DATA, 0);
	if (status.DtStatusFailed()) {
	tile.data = nil
	continue;
	}

					s.m_cacheLayerCount++;
					s.m_cacheCompressedSize += tile.dataSize;
	s.m_cacheRawSize += calcLayerBufferSize(tcparams.width, tcparams.height);
	}
	}
	}

	// Build initial meshes
	now:=time.Now()
	for  y := 0; y < th; y++{
			for  x := 0; x < tw; x++{	s.m_tileCache.BuildNavMeshTilesAt(x,y, s.m_navMesh);}
		}




	s.m_cacheBuildTimeMs =time.Since(now)
	s.m_cacheBuildMemUsage = static_cast<unsigned int>(m_talloc->high);


	nav := s.m_navMesh;
	 navmeshMemUsage := 0;
	for  i := 0; i < nav.GetMaxTiles(); i++{
	tile := nav.GetTile(i);
	if (tile.Header!=nil){
		navmeshMemUsage += tile.DataSize;
	}

	}
	log.Printf("navmeshMemUsage = %.1f kB", navmeshMemUsage/1024.0);


	if (s.m_tool!=nil){s.m_tool.init(s.Sample);}

		s.initToolStates(s.Sample);

	return true;
}

func (s *SampleTempObstacles) handleUpdate(dt  float64 )  {
	s.Sample.handleUpdate(dt);

	if (s.m_navMesh==nil){return;}

	if (s.m_tileCache==nil){return;}


	s.m_tileCache.Update(dt, s.m_navMesh);
}

func (s *SampleTempObstacles) getTilePos(pos []float64, tx, ty *int) {
	if (s.m_geom==nil) {return;}

	 bmin := s.m_geom.getNavMeshBoundsMin();

	 ts := s.m_tileSize*s.m_cellSize;
	*tx = int((pos[0] - bmin[0]) / ts);
	*ty = int((pos[2] - bmin[2]) / ts);
}

func (s *SampleTempObstacles) renderCachedTile(tx, ty, Type int) {
	if (s.m_tileCache!=nil){
		drawDetail(s.m_dd,s.m_tileCache,tx,ty,Type);
	}

}
func (s *SampleTempObstacles) renderCachedTileOverlay(tx, ty int, proj, model []float64, view []int) {
	if (s.m_tileCache!=nil){
		drawDetailOverlay(s.gs,s.m_tileCache, tx, ty, proj, model, view);}

}

func (s *SampleTempObstacles) addTempObstacle(pos []float64) {
	if (s.m_tileCache==nil){return;}

	 p:=make([]float64,3)
	copy(p, pos);
	p[1] -= 0.5;
	s.m_tileCache.AddObstacle(p, 1.0, 2.0, nil);
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
func (s *SampleTempObstacles) rasterizeTileLayers(tx, ty int, cfg *recast.RcConfig, tiles *recast.tileCacheData, maxTiles int) int {
	if (!m_geom || !m_geom->getMesh() || !m_geom->getChunkyMesh())
	{
		m_ctx->log(RC_LOG_ERROR, "buildTile: Input mesh is not specified.");
		return 0;
	}

	FastLZCompressor comp;
	RasterizationContext rc;

	const float* verts = m_geom->getMesh()->getVerts();
	const int nverts = m_geom->getMesh()->getVertCount();
	const rcChunkyTriMesh* chunkyMesh = m_geom->getChunkyMesh();

	// Tile bounds.
	const float tcs = cfg.tileSize * cfg.cs;

	rcConfig tcfg;
	memcpy(&tcfg, &cfg, sizeof(tcfg));

	tcfg.bmin[0] = cfg.bmin[0] + tx*tcs;
	tcfg.bmin[1] = cfg.bmin[1];
	tcfg.bmin[2] = cfg.bmin[2] + ty*tcs;
	tcfg.bmax[0] = cfg.bmin[0] + (tx+1)*tcs;
	tcfg.bmax[1] = cfg.bmax[1];
	tcfg.bmax[2] = cfg.bmin[2] + (ty+1)*tcs;
	tcfg.bmin[0] -= tcfg.borderSize*tcfg.cs;
	tcfg.bmin[2] -= tcfg.borderSize*tcfg.cs;
	tcfg.bmax[0] += tcfg.borderSize*tcfg.cs;
	tcfg.bmax[2] += tcfg.borderSize*tcfg.cs;

	// Allocate voxel heightfield where we rasterize our input data to.
	rc.solid = rcAllocHeightfield();
	if (!rc.solid)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'solid'.");
		return 0;
	}
	if (!rcCreateHeightfield(m_ctx, *rc.solid, tcfg.width, tcfg.height, tcfg.bmin, tcfg.bmax, tcfg.cs, tcfg.ch))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not create solid heightfield.");
		return 0;
	}

	// Allocate array that can hold triangle flags.
	// If you have multiple meshes you need to process, allocate
	// and array which can hold the max number of triangles you need to process.
	rc.triareas = new unsigned char[chunkyMesh->maxTrisPerChunk];
	if (!rc.triareas)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'm_triareas' (%d).", chunkyMesh->maxTrisPerChunk);
		return 0;
	}

	float tbmin[2], tbmax[2];
	tbmin[0] = tcfg.bmin[0];
	tbmin[1] = tcfg.bmin[2];
	tbmax[0] = tcfg.bmax[0];
	tbmax[1] = tcfg.bmax[2];
	int cid[512];// TODO: Make grow when returning too many items.
	const int ncid = rcGetChunksOverlappingRect(chunkyMesh, tbmin, tbmax, cid, 512);
	if (!ncid)
	{
		return 0; // empty
	}

	for (int i = 0; i < ncid; ++i)
	{
	const rcChunkyTriMeshNode& node = chunkyMesh->nodes[cid[i]];
	const int* tris = &chunkyMesh->tris[node.i*3];
	const int ntris = node.n;

	memset(rc.triareas, 0, ntris*sizeof(unsigned char));
	rcMarkWalkableTriangles(m_ctx, tcfg.walkableSlopeAngle,
	verts, nverts, tris, ntris, rc.triareas);

	if (!rcRasterizeTriangles(m_ctx, verts, nverts, tris, rc.triareas, ntris, *rc.solid, tcfg.walkableClimb))
	return 0;
	}

	// Once all geometry is rasterized, we do initial pass of filtering to
	// remove unwanted overhangs caused by the conservative rasterization
	// as well as filter spans where the character cannot possibly stand.
	if (m_filterLowHangingObstacles)
		rcFilterLowHangingWalkableObstacles(m_ctx, tcfg.walkableClimb, *rc.solid);
	if (m_filterLedgeSpans)
		rcFilterLedgeSpans(m_ctx, tcfg.walkableHeight, tcfg.walkableClimb, *rc.solid);
	if (m_filterWalkableLowHeightSpans)
		rcFilterWalkableLowHeightSpans(m_ctx, tcfg.walkableHeight, *rc.solid);


	rc.chf = rcAllocCompactHeightfield();
	if (!rc.chf)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'chf'.");
		return 0;
	}
	if (!rcBuildCompactHeightfield(m_ctx, tcfg.walkableHeight, tcfg.walkableClimb, *rc.solid, *rc.chf))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build compact data.");
		return 0;
	}

	// Erode the walkable area by agent radius.
	if (!rcErodeWalkableArea(m_ctx, tcfg.walkableRadius, *rc.chf))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not erode.");
		return 0;
	}

	// (Optional) Mark areas.
	const ConvexVolume* vols = m_geom->getConvexVolumes();
	for (int i  = 0; i < m_geom->getConvexVolumeCount(); ++i)
	{
	rcMarkConvexPolyArea(m_ctx, vols[i].verts, vols[i].nverts,
	vols[i].hmin, vols[i].hmax,
	(unsigned char)vols[i].area, *rc.chf);
	}

	rc.lset = rcAllocHeightfieldLayerSet();
	if (!rc.lset)
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Out of memory 'lset'.");
		return 0;
	}
	if (!rcBuildHeightfieldLayers(m_ctx, *rc.chf, tcfg.borderSize, tcfg.walkableHeight, *rc.lset))
	{
		m_ctx->log(RC_LOG_ERROR, "buildNavigation: Could not build heighfield layers.");
		return 0;
	}

	rc.ntiles = 0;
	for (int i = 0; i < rcMin(rc.lset->nlayers, MAX_LAYERS); ++i)
	{
	TileCacheData* tile = &rc.tiles[rc.ntiles++];
	const rcHeightfieldLayer* layer = &rc.lset->layers[i];

	// Store header
	dtTileCacheLayerHeader header;
	header.magic = DT_TILECACHE_MAGIC;
	header.version = DT_TILECACHE_VERSION;

	// Tile layer location in the navmesh.
	header.tx = tx;
	header.ty = ty;
	header.tlayer = i;
	dtVcopy(header.bmin, layer->bmin);
	dtVcopy(header.bmax, layer->bmax);

	// Tile info.
	header.width = (unsigned char)layer->width;
	header.height = (unsigned char)layer->height;
	header.minx = (unsigned char)layer->minx;
	header.maxx = (unsigned char)layer->maxx;
	header.miny = (unsigned char)layer->miny;
	header.maxy = (unsigned char)layer->maxy;
	header.hmin = (unsigned short)layer->hmin;
	header.hmax = (unsigned short)layer->hmax;

	dtStatus status = dtBuildTileCacheLayer(&comp, &header, layer->heights, layer->areas, layer->cons,
	&tile->data, &tile->dataSize);
	if (dtStatusFailed(status))
	{
	return 0;
	}
	}

	// Transfer ownsership of tile data from build context to the caller.
	int n = 0;
	for (int i = 0; i < rcMin(rc.ntiles, maxTiles); ++i)
	{
	tiles[n++] = rc.tiles[i];
	rc.tiles[i].data = 0;
	rc.tiles[i].dataSize = 0;
	}

	return n;
}
func (s *SampleTempObstacles) buildNavMeshTilesAt(tx, ty int, navmesh *recast.DtNavMesh) recast.DtStatus {
	const MAX_TILES = 32
	var tiles [MAX_TILES]dtCompressedTileRef
	ntiles := s.getTilesAt(tx, ty, tiles, MAX_TILES)

	for i := 0; i < ntiles; i++ {
		status := s.buildNavMeshTile(tiles[i], navmesh)
		if status.DtStatusFailed() {
			return status
		}

	}

	return DT_SUCCESS
}

type MeshProcess struct {
	m_geom *InputGeom ;
}

func newMeshProcess(m_geom *InputGeom)*MeshProcess {
	return &MeshProcess{m_geom:m_geom}
}
func (p*MeshProcess)Process(params *recast.DtNavMeshCreateParams ,
	 polyAreas []int ,  polyFlags[]int )  {
	// Update poly flags from areas.
	for i := 0; i < params.PolyCount; i++{
	if (polyAreas[i] == recast.DT_TILECACHE_WALKABLE_AREA){
		polyAreas[i] = SAMPLE_POLYAREA_GROUND;


	if (polyAreas[i] == SAMPLE_POLYAREA_GROUND ||
	polyAreas[i] == SAMPLE_POLYAREA_GRASS ||
	polyAreas[i] == SAMPLE_POLYAREA_ROAD) {
	polyFlags[i] = SAMPLE_POLYFLAGS_WALK;
	} else if (polyAreas[i] == SAMPLE_POLYAREA_WATER) {
	polyFlags[i] = SAMPLE_POLYFLAGS_SWIM;
	} else if (polyAreas[i] == SAMPLE_POLYAREA_DOOR) {
	polyFlags[i] = SAMPLE_POLYFLAGS_WALK | SAMPLE_POLYFLAGS_DOOR;
	}
	}

	// Pass in off-mesh connections.
	if (p.m_geom!=nil) {
		params.OffMeshConVerts = p.m_geom.getOffMeshConnectionVerts();
		params.OffMeshConRad = p.m_geom.getOffMeshConnectionRads();
		params.OffMeshConDir = p.m_geom.getOffMeshConnectionDirs();
		params.OffMeshConAreas = p.m_geom.getOffMeshConnectionAreas();
		params.OffMeshConFlags = p.m_geom.getOffMeshConnectionFlags();
		params.OffMeshConUserID = p.m_geom.getOffMeshConnectionId();
		params.OffMeshConCount = p.m_geom.getOffMeshConnectionCount();
	}
}
}

const  MAX_LAYERS = 32;

 type TileCacheData struct {
 data []byte
 dataSize int
};
