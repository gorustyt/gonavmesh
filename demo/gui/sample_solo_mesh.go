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

type SoloMeshDrawMode int

const (
	SOLOMESH_DRAWMODE_NAVMESH SoloMeshDrawMode = iota
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
	m_totalBuildTimeMs time.Duration

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
		s.setTool(newCrowdTool(s.gs))
	}

	s.gs.imguiSeparatorLine()

	s.gs.imguiIndent()

	if s.m_tool != nil {
		s.m_tool.handleMenu()
	}

	s.gs.imguiUnindent()

}

func (s *SampleSoloMesh) handleDebugMode() {
	// Check which modes are valid.
	valid := make([]bool, SOLOMESH_MAX_DRAWMODE)
	for i := 0; i < int(SOLOMESH_MAX_DRAWMODE); i++ {
		valid[i] = false
	}

	if s.m_geom != nil {
		valid[SOLOMESH_DRAWMODE_NAVMESH] = s.m_navMesh != nil
		valid[SOLOMESH_DRAWMODE_NAVMESH_TRANS] = s.m_navMesh != nil
		valid[SOLOMESH_DRAWMODE_NAVMESH_BVTREE] = s.m_navMesh != nil
		valid[SOLOMESH_DRAWMODE_NAVMESH_NODES] = s.m_navQuery != nil
		valid[SOLOMESH_DRAWMODE_NAVMESH_INVIS] = s.m_navMesh != nil
		valid[SOLOMESH_DRAWMODE_MESH] = true
		valid[SOLOMESH_DRAWMODE_VOXELS] = s.m_solid != nil
		valid[SOLOMESH_DRAWMODE_VOXELS_WALKABLE] = s.m_solid != nil
		valid[SOLOMESH_DRAWMODE_COMPACT] = s.m_chf != nil
		valid[SOLOMESH_DRAWMODE_COMPACT_DISTANCE] = s.m_chf != nil
		valid[SOLOMESH_DRAWMODE_COMPACT_REGIONS] = s.m_chf != nil
		valid[SOLOMESH_DRAWMODE_REGION_CONNECTIONS] = s.m_cset != nil
		valid[SOLOMESH_DRAWMODE_RAW_CONTOURS] = s.m_cset != nil
		valid[SOLOMESH_DRAWMODE_BOTH_CONTOURS] = s.m_cset != nil
		valid[SOLOMESH_DRAWMODE_CONTOURS] = s.m_cset != nil
		valid[SOLOMESH_DRAWMODE_POLYMESH] = s.m_pmesh != nil
		valid[SOLOMESH_DRAWMODE_POLYMESH_DETAIL] = s.m_dmesh != nil
	}

	unavail := 0
	for i := 0; i < int(SOLOMESH_MAX_DRAWMODE); i++ {
		if !valid[i] {
			unavail++
		}
	}

	if unavail == int(SOLOMESH_MAX_DRAWMODE) {
		return
	}

	s.gs.imguiLabel("Draw")
	if s.gs.imguiCheck("Input Mesh", s.m_drawMode == SOLOMESH_DRAWMODE_MESH, valid[SOLOMESH_DRAWMODE_MESH]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_MESH
	}

	if s.gs.imguiCheck("Navmesh", s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH, valid[SOLOMESH_DRAWMODE_NAVMESH]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_NAVMESH
	}

	if s.gs.imguiCheck("Navmesh Invis", s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_INVIS, valid[SOLOMESH_DRAWMODE_NAVMESH_INVIS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_NAVMESH_INVIS
	}

	if s.gs.imguiCheck("Navmesh Trans", s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_TRANS, valid[SOLOMESH_DRAWMODE_NAVMESH_TRANS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_NAVMESH_TRANS
	}

	if s.gs.imguiCheck("Navmesh BVTree", s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_BVTREE, valid[SOLOMESH_DRAWMODE_NAVMESH_BVTREE]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_NAVMESH_BVTREE
	}

	if s.gs.imguiCheck("Navmesh Nodes", s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_NODES, valid[SOLOMESH_DRAWMODE_NAVMESH_NODES]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_NAVMESH_NODES
	}

	if s.gs.imguiCheck("Voxels", s.m_drawMode == SOLOMESH_DRAWMODE_VOXELS, valid[SOLOMESH_DRAWMODE_VOXELS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_VOXELS
	}

	if s.gs.imguiCheck("Walkable Voxels", s.m_drawMode == SOLOMESH_DRAWMODE_VOXELS_WALKABLE, valid[SOLOMESH_DRAWMODE_VOXELS_WALKABLE]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_VOXELS_WALKABLE
	}

	if s.gs.imguiCheck("Compact", s.m_drawMode == SOLOMESH_DRAWMODE_COMPACT, valid[SOLOMESH_DRAWMODE_COMPACT]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_COMPACT
	}

	if s.gs.imguiCheck("Compact Distance", s.m_drawMode == SOLOMESH_DRAWMODE_COMPACT_DISTANCE, valid[SOLOMESH_DRAWMODE_COMPACT_DISTANCE]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_COMPACT_DISTANCE
	}

	if s.gs.imguiCheck("Compact Regions", s.m_drawMode == SOLOMESH_DRAWMODE_COMPACT_REGIONS, valid[SOLOMESH_DRAWMODE_COMPACT_REGIONS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_COMPACT_REGIONS
	}

	if s.gs.imguiCheck("Region Connections", s.m_drawMode == SOLOMESH_DRAWMODE_REGION_CONNECTIONS, valid[SOLOMESH_DRAWMODE_REGION_CONNECTIONS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_REGION_CONNECTIONS
	}

	if s.gs.imguiCheck("Raw Contours", s.m_drawMode == SOLOMESH_DRAWMODE_RAW_CONTOURS, valid[SOLOMESH_DRAWMODE_RAW_CONTOURS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_RAW_CONTOURS
	}

	if s.gs.imguiCheck("Both Contours", s.m_drawMode == SOLOMESH_DRAWMODE_BOTH_CONTOURS, valid[SOLOMESH_DRAWMODE_BOTH_CONTOURS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_BOTH_CONTOURS
	}

	if s.gs.imguiCheck("Contours", s.m_drawMode == SOLOMESH_DRAWMODE_CONTOURS, valid[SOLOMESH_DRAWMODE_CONTOURS]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_CONTOURS
	}

	if s.gs.imguiCheck("Poly Mesh", s.m_drawMode == SOLOMESH_DRAWMODE_POLYMESH, valid[SOLOMESH_DRAWMODE_POLYMESH]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_POLYMESH
	}

	if s.gs.imguiCheck("Poly Mesh Detail", s.m_drawMode == SOLOMESH_DRAWMODE_POLYMESH_DETAIL, valid[SOLOMESH_DRAWMODE_POLYMESH_DETAIL]) {
		s.m_drawMode = SOLOMESH_DRAWMODE_POLYMESH_DETAIL
	}

	if unavail > 0 {
		s.gs.imguiValue("Tick 'Keep Itermediate Results'")
		s.gs.imguiValue("to see more debug mode options.")
	}
}

func (s *SampleSoloMesh) handleRender() {
	if s.m_geom == nil || s.m_geom.getMesh() == nil {
		return

	}

	gl.Enable(GL_FOG)
	gl.DepthMask(true)

	texScale := 1.0 / (s.m_cellSize * 10.0)

	if s.m_drawMode != SOLOMESH_DRAWMODE_NAVMESH_TRANS {
		// Draw mesh
		debug_utils.DuDebugDrawTriMeshSlope(s.m_dd, s.m_geom.getMesh().getVerts(), s.m_geom.getMesh().getVertCount(),
			s.m_geom.getMesh().getTris(), s.m_geom.getMesh().getNormals(), s.m_geom.getMesh().getTriCount(),
			s.m_agentMaxSlope, texScale)
		s.m_geom.drawOffMeshConnections(s.m_dd)
	}

	glDisable(GL_FOG)
	gl.DepthMask(false)

	// Draw bounds
	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	debug_utils.DuDebugDrawBoxWire(s.m_dd, bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2], debug_utils.DuRGBA(255, 255, 255, 128), 1.0)
	s.m_dd.Begin(debug_utils.DU_DRAW_POINTS, 5.0)
	s.m_dd.Vertex1(bmin[0], bmin[1], bmin[2], debug_utils.DuRGBA(255, 255, 255, 128))
	s.m_dd.End()

	if s.m_navMesh != nil && s.m_navQuery != nil &&
		(s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH ||
			s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_TRANS ||
			s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_BVTREE ||
			s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_NODES ||
			s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_INVIS) {
		if s.m_drawMode != SOLOMESH_DRAWMODE_NAVMESH_INVIS {
			debug_utils.DuDebugDrawNavMeshWithClosedList(s.m_dd, s.m_navMesh, s.m_navQuery, s.m_navMeshDrawFlags)
		}

		if s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_BVTREE {
			debug_utils.DuDebugDrawNavMeshBVTree(s.m_dd, s.m_navMesh)
		}

		if s.m_drawMode == SOLOMESH_DRAWMODE_NAVMESH_NODES {
			debug_utils.DuDebugDrawNavMeshNodes(s.m_dd, s.m_navQuery)
		}

		debug_utils.DuDebugDrawNavMeshPolysWithFlags(s.m_dd, s.m_navMesh, SAMPLE_POLYFLAGS_DISABLED, debug_utils.DuRGBA(0, 0, 0, 128))
	}

	gl.DepthMask(true)

	if s.m_chf != nil && s.m_drawMode == SOLOMESH_DRAWMODE_COMPACT {
		debug_utils.DuDebugDrawCompactHeightfieldSolid(s.m_dd, s.m_chf)
	}

	if s.m_chf != nil && s.m_drawMode == SOLOMESH_DRAWMODE_COMPACT_DISTANCE {
		debug_utils.DuDebugDrawCompactHeightfieldDistance(s.m_dd, s.m_chf)
	}

	if s.m_chf != nil && s.m_drawMode == SOLOMESH_DRAWMODE_COMPACT_REGIONS {
		debug_utils.DuDebugDrawCompactHeightfieldRegions(s.m_dd, s.m_chf)
	}

	if s.m_solid != nil && s.m_drawMode == SOLOMESH_DRAWMODE_VOXELS {
		gl.Enable(gl.FOG)
		debug_utils.DuDebugDrawHeightfieldSolid(s.m_dd, s.m_solid)
		glDisable(GL_FOG)
	}
	if s.m_solid != nil && s.m_drawMode == SOLOMESH_DRAWMODE_VOXELS_WALKABLE {
		glEnable(GL_FOG)
		debug_utils.DuDebugDrawHeightfieldWalkable(s.m_dd, s.m_solid)
		glDisable(GL_FOG)
	}
	if s.m_cset != nil && s.m_drawMode == SOLOMESH_DRAWMODE_RAW_CONTOURS {
		gl.DepthMask(false)
		debug_utils.DuDebugDrawRawContours(s.m_dd, s.m_cset)
		gl.DepthMask(true)
	}
	if s.m_cset != nil && s.m_drawMode == SOLOMESH_DRAWMODE_BOTH_CONTOURS {
		gl.DepthMask(false)
		debug_utils.DuDebugDrawRawContours(s.m_dd, s.m_cset, 0.5)
		debug_utils.DuDebugDrawContours(s.m_dd, s.m_cset)
		gl.DepthMask(true)
	}
	if s.m_cset != nil && s.m_drawMode == SOLOMESH_DRAWMODE_CONTOURS {
		gl.DepthMask(false)
		debug_utils.DuDebugDrawContours(s.m_dd, s.m_cset)
		gl.DepthMask(true)
	}
	if s.m_chf != nil && s.m_cset != nil && s.m_drawMode == SOLOMESH_DRAWMODE_REGION_CONNECTIONS {
		debug_utils.DuDebugDrawCompactHeightfieldRegions(s.m_dd, s.m_chf)

		gl.DepthMask(false)
		debug_utils.DuDebugDrawRegionConnections(s.m_dd, s.m_cset)
		gl.DepthMask(true)
	}
	if s.m_pmesh != nil && s.m_drawMode == SOLOMESH_DRAWMODE_POLYMESH {
		gl.DepthMask(false)
		debug_utils.DuDebugDrawPolyMesh(s.m_dd, s.m_pmesh)
		gl.DepthMask(true)
	}
	if s.m_dmesh != nil && s.m_drawMode == SOLOMESH_DRAWMODE_POLYMESH_DETAIL {
		gl.DepthMask(false)
		debug_utils.DuDebugDrawPolyMeshDetail(s.m_dd, s.m_dmesh)
		gl.DepthMask(true)
	}

	s.m_geom.drawConvexVolumes(s.m_dd)

	if s.m_tool != nil {
	}
	s.m_tool.handleRender()

	s.renderToolStates()

	gl.DepthMask(true)
}

func (s *SampleSoloMesh) handleRenderOverlay(proj, model []float64, view []int) {
	if s.m_tool != nil {
		s.m_tool.handleRenderOverlay(proj, model, view)
	}

	s.renderOverlayToolStates(proj, model, view)
}
func (s *SampleSoloMesh) handleMeshChanged(geom *InputGeom) {
	s.Sample.handleMeshChanged(geom)
	s.m_navMesh = nil

	if s.m_tool != nil {
		s.m_tool.reset()
		s.m_tool.init(s.Sample)
	}
	s.resetToolStates()
	s.initToolStates(s.Sample)
}

func (s *SampleSoloMesh) handleBuild() bool {
	if s.m_geom == nil || s.m_geom.getMesh() == nil {
		log.Printf("buildNavigation: Input mesh is not specified.")
		return false
	}

	s.cleanup()

	bmin := s.m_geom.getNavMeshBoundsMin()
	bmax := s.m_geom.getNavMeshBoundsMax()
	verts := s.m_geom.getMesh().getVerts()
	nverts := s.m_geom.getMesh().getVertCount()
	tris := s.m_geom.getMesh().getTris()
	ntris := s.m_geom.getMesh().getTriCount()

	//
	// Step 1. Initialize build config.
	//
	s.m_cfg = &recast.RcConfig{}

	// Init build configuration from GUI

	s.m_cfg.Cs = s.m_cellSize
	s.m_cfg.Ch = s.m_cellHeight
	s.m_cfg.WalkableSlopeAngle = s.m_agentMaxSlope
	s.m_cfg.WalkableHeight = int(math.Ceil(s.m_agentHeight / s.m_cfg.Ch))
	s.m_cfg.WalkableClimb = int(math.Floor(s.m_agentMaxClimb / s.m_cfg.Ch))
	s.m_cfg.WalkableRadius = int(math.Ceil(s.m_agentRadius / s.m_cfg.Cs))
	s.m_cfg.MaxEdgeLen = int(s.m_edgeMaxLen / s.m_cellSize)
	s.m_cfg.MaxSimplificationError = s.m_edgeMaxError
	s.m_cfg.MinRegionArea = int(common.Sqr(s.m_regionMinSize))     // Note: area = size*size
	s.m_cfg.MergeRegionArea = int(common.Sqr(s.m_regionMergeSize)) // Note: area = size*size
	s.m_cfg.MaxVertsPerPoly = int(s.m_vertsPerPoly)
	s.m_cfg.DetailSampleDist = s.m_cellSize * s.m_detailSampleDist
	if s.m_detailSampleDist < 0.9 {
		s.m_cfg.DetailSampleDist = 0
	}
	s.m_cfg.DetailSampleMaxError = s.m_cellHeight * s.m_detailSampleMaxError
	now := time.Now()
	// Set the area where the navigation will be build.
	// Here the bounds of the input mesh are used, but the
	// area could be specified by an user defined box, etc.
	copy(s.m_cfg.Bmin[:], bmin)
	copy(s.m_cfg.Bmax[:], bmax)
	recast.RcCalcGridSize(s.m_cfg.Bmin[:], s.m_cfg.Bmax[:], s.m_cfg.Cs, &s.m_cfg.Width, &s.m_cfg.Height)

	log.Printf("Building navigation:")
	log.Printf(" - %d x %d cells", s.m_cfg.Width, s.m_cfg.Height)
	log.Printf(" - %.1fK verts, %.1fK tris", nverts/1000.0, ntris/1000.0)

	//
	// Step 2. Rasterize input polygon soup.
	//

	// Allocate voxel heightfield where we rasterize our input data to.
	s.m_solid = recast.RcCreateHeightfield(s.m_cfg.Width, s.m_cfg.Height, s.m_cfg.Bmin[:], s.m_cfg.Bmax[:], s.m_cfg.Cs, s.m_cfg.Ch)
	// Allocate array that can hold triangle area types.
	// If you have multiple meshes you need to process, allocate
	// and array which can hold the max number of triangles you need to process.
	s.m_triareas = make([]int, ntris)
	// Find triangles which are walkable based on their slope and rasterize them.
	// If your input data is multiple meshes, you can transform them here, calculate
	// the are type for each of the meshes and rasterize them.
	recast.RcMarkWalkableTriangles(s.m_cfg.WalkableSlopeAngle, verts, nverts, tris, ntris, s.m_triareas)
	if !recast.RcRasterizeTriangles(verts, nverts, tris, s.m_triareas, ntris, s.m_solid, s.m_cfg.WalkableClimb) {
		log.Printf("buildNavigation: Could not rasterize triangles.")
		return false
	}

	if !s.m_keepInterResults {
		s.m_triareas = nil
	}

	//
	// Step 3. Filter walkable surfaces.
	//

	// Once all geometry is rasterized, we do initial pass of filtering to
	// remove unwanted overhangs caused by the conservative rasterization
	// as well as filter spans where the character cannot possibly stand.
	if s.m_filterLowHangingObstacles {
		recast.RcFilterLowHangingWalkableObstacles(s.m_cfg.WalkableClimb, s.m_solid)
	}

	if s.m_filterLedgeSpans {
		recast.RcFilterLedgeSpans(s.m_cfg.WalkableHeight, s.m_cfg.WalkableClimb, s.m_solid)
	}

	if s.m_filterWalkableLowHeightSpans {
		recast.RcFilterWalkableLowHeightSpans(s.m_cfg.WalkableHeight, s.m_solid)
	}

	//
	// Step 4. Partition walkable surface to simple regions.
	//

	// Compact the heightfield so that it is faster to handle from now on.
	// This will result more cache coherent data as well as the neighbours
	// between walkable cells will be calculated.
	s.m_chf = &recast.RcCompactHeightfield{}

	if !recast.RcBuildCompactHeightfield(s.m_cfg.WalkableHeight, s.m_cfg.WalkableClimb, s.m_solid, s.m_chf) {
		log.Printf("buildNavigation: Could not build compact data.")
		return false
	}

	if !s.m_keepInterResults {

		s.m_solid = nil
	}

	// Erode the walkable area by agent radius.
	if !recast.RcErodeWalkableArea(s.m_cfg.WalkableRadius, s.m_chf) {
		log.Printf("buildNavigation: Could not erode.")
		return false
	}

	// (Optional) Mark areas.
	vols := s.m_geom.getConvexVolumes()
	for i := 0; i < s.m_geom.getConvexVolumeCount(); i++ {
		recast.RcMarkConvexPolyArea(vols[i].verts, vols[i].nverts, vols[i].hmin, vols[i].hmax,
			vols[i].area, s.m_chf)

	}

	// Partition the heightfield so that we can use simple algorithm later to triangulate the walkable areas.
	// There are 3 partitioning methods, each with some pros and cons:
	// 1) Watershed partitioning
	//   - the classic Recast partitioning
	//   - creates the nicest tessellation
	//   - usually slowest
	//   - partitions the heightfield into nice regions without holes or overlaps
	//   - the are some corner cases where this method creates produces holes and overlaps
	//      - holes may appear when a small obstacles is close to large open area (triangulation can handle this)
	//      - overlaps may occur if you have narrow spiral corridors (i.e stairs), this make triangulation to fail
	//   * generally the best choice if you precompute the navmesh, use this if you have large open areas
	// 2) Monotone partitioning
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
			return false
		}

		// Partition the walkable surface into simple regions without holes.
		if !recast.RcBuildRegions(s.m_chf, 0, s.m_cfg.MinRegionArea, s.m_cfg.MergeRegionArea) {
			log.Printf("buildNavigation: Could not build watershed regions.")
			return false
		}
	} else if s.m_partitionType == SAMPLE_PARTITION_MONOTONE {
		// Partition the walkable surface into simple regions without holes.
		// Monotone partitioning does not need distancefield.
		if !recast.RcBuildRegionsMonotone(s.m_chf, 0, s.m_cfg.MinRegionArea, s.m_cfg.MergeRegionArea) {
			log.Printf("buildNavigation: Could not build monotone regions.")
			return false
		}
	} else // SAMPLE_PARTITION_LAYERS
	{
		// Partition the walkable surface into simple regions without holes.
		if !recast.RcBuildLayerRegions(s.m_chf, 0, s.m_cfg.MinRegionArea) {
			log.Printf("buildNavigation: Could not build layer regions.")
			return false
		}
	}

	//
	// Step 5. Trace and simplify region contours.
	//

	// Create contours.
	s.m_cset = &recast.RcContourSet{}
	if !recast.RcBuildContours(s.m_chf, s.m_cfg.MaxSimplificationError, s.m_cfg.MaxEdgeLen, s.m_cset) {
		log.Printf("buildNavigation: Could not create contours.")
		return false
	}

	//
	// Step 6. Build polygons mesh from contours.
	//

	// Build polygon navmesh from the contours.
	s.m_pmesh = &recast.RcPolyMesh{}

	if !recast.RcBuildPolyMesh(s.m_cset, s.m_cfg.MaxVertsPerPoly, s.m_pmesh) {
		log.Printf("buildNavigation: Could not triangulate contours.")
		return false
	}

	//
	// Step 7. Create detail mesh which allows to access approximate height on each polygon.
	//

	s.m_dmesh = &recast.RcPolyMeshDetail{}
	if !recast.RcBuildPolyMeshDetail(s.m_pmesh, s.m_chf, s.m_cfg.DetailSampleDist, s.m_cfg.DetailSampleMaxError, s.m_dmesh) {
		log.Printf("buildNavigation: Could not build detail mesh.")
		return false
	}

	if !s.m_keepInterResults {
		s.m_chf = nil
		s.m_cset = nil
	}

	// At this point the navigation mesh data is ready, you can access it from m_pmesh.
	// See duDebugDrawPolyMesh or dtCreateNavMeshData as examples how to access the data.

	//
	// (Optional) Step 8. Create Detour data from Recast poly mesh.
	//

	// The GUI may allow more max points per polygon than Detour can handle.
	// Only build the detour navmesh if we do not exceed the limit.
	if s.m_cfg.MaxVertsPerPoly <= recast.DT_VERTS_PER_POLYGON {
		// Update poly flags from are as.
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

		var params recast.DtNavMeshCreateParams
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
		copy(params.Bmin[:], s.m_pmesh.Bmin)
		copy(params.Bmax[:], s.m_pmesh.Bmax)
		params.Cs = s.m_cfg.Cs
		params.Ch = s.m_cfg.Ch
		params.BuildBvTree = true

		navData, ok := recast.DtCreateNavMeshData(&params)
		if !ok {
			log.Printf("Could not build Detour navmesh.")
			return false
		}

		s.m_navMesh, _, _ = recast.NewDtNavMesh(navData, recast.DT_TILE_FREE_DATA)
		s.m_navQuery = recast.NewDtNavMeshQuery(s.m_navMesh, 2048)

	}
	// Show performance stats.
	s.m_totalBuildTimeMs = time.Since(now)
	log.Printf(">> Polymesh: %d vertices  %d polygons", s.m_pmesh.Nverts, s.m_pmesh.Npolys)

	if s.m_tool != nil {
		s.m_tool.init(s.Sample)
	}

	s.initToolStates(s.Sample)

	return true
}
