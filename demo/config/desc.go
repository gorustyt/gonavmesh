package config

type SamplePolyFlags int

const (
	SAMPLE_POLYFLAGS_WALK     = 0x01   // Ability to walk (ground, grass, road)
	SAMPLE_POLYFLAGS_SWIM     = 0x02   // Ability to swim (water).
	SAMPLE_POLYFLAGS_DOOR     = 0x04   // Ability to move through doors.
	SAMPLE_POLYFLAGS_JUMP     = 0x08   // Ability to jump.
	SAMPLE_POLYFLAGS_DISABLED = 0x10   // Disabled polygon
	SAMPLE_POLYFLAGS_ALL      = 0xffff // All abilities.

	DESC_SAMPLE_POLYFLAGS_WALK = "Walk"
	DESC_SAMPLE_POLYFLAGS_SWIM = "Swim"
	DESC_SAMPLE_POLYFLAGS_DOOR = "Door"
	DESC_SAMPLE_POLYFLAGS_JUMP = "Jump"
)

type ToolMode int

const (
	TOOLMODE_PATHFIND_FOLLOW = iota
	TOOLMODE_PATHFIND_STRAIGHT
	TOOLMODE_PATHFIND_SLICED
	TOOLMODE_RAYCAST
	TOOLMODE_DISTANCE_TO_WALL
	TOOLMODE_FIND_POLYS_IN_CIRCLE
	TOOLMODE_FIND_POLYS_IN_SHAPE
	TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD

	Desc_TOOLMODE_PATHFIND_FOLLOW          = "Pathfind Follow"
	Desc_TOOLMODE_PATHFIND_STRAIGHT        = "Pathfind Straight"
	Desc_TOOLMODE_PATHFIND_SLICED          = "Pathfind Sliced"
	Desc_TOOLMODE_RAYCAST                  = "Raycast"
	Desc_TOOLMODE_DISTANCE_TO_WALL         = "Distance to Wall"
	Desc_TOOLMODE_FIND_POLYS_IN_CIRCLE     = "Find Polys in Circle"
	Desc_TOOLMODE_FIND_POLYS_IN_SHAPE      = "Find Polys in Shape"
	Desc_TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD = "Find Local Neighbourhood"
)

const (
	DT_STRAIGHTPATH_NONE_CROSSINGS = "none"
	DT_STRAIGHTPATH_AREA_CROSSINGS = "Area"
	DT_STRAIGHTPATH_ALL_CROSSINGS  = "All"
)

type CrowdToolToolMode int

const (
	CrowdTool_TOOLMODE_CREATE CrowdToolToolMode = iota
	CrowdTool_TOOLMODE_MOVE_TARGET
	CrowdTool_TOOLMODE_SELECT
	CrowdTool_TOOLMODE_TOGGLE_POLYS

	DescCrowdTool_TOOLMODE_CREATE       = "Create Agents"
	DescCrowdTool_TOOLMODE_MOVE_TARGET  = "Move Target"
	DescCrowdTool_TOOLMODE_SELECT       = "Select Agent"
	DescCrowdTool_TOOLMODE_TOGGLE_POLYS = "Toggle Polys"
)

const (
	ExpandOptionsOptimizeVisibility = "Optimize Visibility"
	ExpandOptionsOptimizeTopology   = "Optimize Topology"
	ExpandOptionsAnticipateTurns    = "Anticipate Turns"
	ExpandOptionsObstacleAvoidance  = "Obstacle Avoidance"
	ExpandOptionsSeparation         = "Separation"

	ExpandSelectedDebugDrawShowCorners          = "Show Corners"
	ExpandSelectedDebugDrawShowCollisionSegs    = "Show Collision Segs"
	ExpandSelectedDebugDrawShowPath             = "Show Path"
	ExpandSelectedDebugDrawShowVO               = "Show VO"
	ExpandSelectedDebugDrawShowPathOptimization = "Show Path Optimization"
	ExpandSelectedDebugDrawShowNeighbours       = "Show Neighbours"

	ExpandDebugDrawShowLabels    = "Show Labels"
	ExpandDebugDrawShowProxGrid  = "Show Prox Grid"
	ExpandDebugDrawShowNodes     = "Show Nodes"
	ExpandDebugDrawShowPerfGraph = "Show Perf Graph"
	ExpandDebugDrawShowDetailAll = "Show Detail All"
)

type SamplePartitionType int

const (
	SAMPLE_PARTITION_WATERSHED = iota
	SAMPLE_PARTITION_MONOTONE
	SAMPLE_PARTITION_LAYERS

	DescSAMPLE_PARTITION_WATERSHED = "Watershed"
	DescSAMPLE_PARTITION_MONOTONE  = "Monotone"
	DescSAMPLE_PARTITION_LAYERS    = "Layers"
)

const (
	FilteringLowHangingObstacles    = "Low Hanging Obstacles"
	FilteringLedgeSpans             = "Ledge Spans"
	FilteringWalkableLowHeightSpans = "Walkable Low Height Spans"
)

const (
	DrawInputMesh         = "Input Mesh"
	DrawNavmesh           = "Navmesh"
	DrawNavmeshInvis      = "Navmesh Invis"
	DrawNavmeshTrans      = "Navmesh Trans"
	DrawNavmeshBVTree     = "Navmesh BVTree"
	DrawNavmeshNodes      = "Navmesh Nodes"
	DrawVoxels            = "Voxels"
	DrawWalkableVoxels    = "Walkable Voxels"
	DrawCompact           = "Compact"
	DrawCompactDistance   = "Compact Distance"
	DrawCompactRegions    = "Compact Regions"
	DrawRegionConnections = "Region Connections"
	DrawRawContours       = "Raw Contours"
	DrawBothContours      = "Both Contours"
	DrawContours          = "Contours"
	DrawPolyMesh          = "Poly Mesh"
	DrawPolyMeshDetail    = "Poly Mesh Detail"
	DrawPortals           = "Navmesh Portals"
	DrawCacheBounds       = "Cache Bounds"

	KeepItermediateResults = "Keep Itermediate Results"
	KeepBuildAllTiles      = "Build All Tiles"
)

const (
	SAMPLE_POLYAREA_GROUND = iota
	SAMPLE_POLYAREA_WATER
	SAMPLE_POLYAREA_ROAD
	SAMPLE_POLYAREA_DOOR
	SAMPLE_POLYAREA_GRASS
	SAMPLE_POLYAREA_JUMP

	Desc_SAMPLE_POLYAREA_GROUND = "Ground"
	Desc_SAMPLE_POLYAREA_WATER  = "Water"
	Desc_SAMPLE_POLYAREA_ROAD   = "Road"
	Desc_SAMPLE_POLYAREA_DOOR   = "Door"
	Desc_SAMPLE_POLYAREA_GRASS  = "Grass"
	Desc_SAMPLE_POLYAREA_JUMP   = "Jump"
)

const (
	OneWay        = "One Way"
	Bidirectional = "Bidirectional"

	ShowLog   = "show logs"
	ShowTools = "show tools"
)

const (
	Desc_TOOL_NAVMESH_TESTER      = "Test Navmesh"
	Desc_TOOL_NAVMESH_PRUNE       = "Prune Navmesh"
	Desc_TOOL_OFFMESH_CONNECTION  = "Create Off-Mesh Connections"
	Desc_TOOL_OFFMESH_Links       = "Create Off-Mesh Links"
	Desc_TOOL_CONVEX_VOLUME       = "Create Convex Volumes"
	Desc_TOOL_CROWD               = "Create Crowds"
	Desc_TOOL_CreateTiles         = "Create Tiles"
	Desc_TOOL_CreateTempObstacles = "Create Temp Obstacles"
	Desc_TOOL_HighlightTileCache  = "Highlight Tile Cache"
)

// / Tool types.
type SampleToolType int

const (
	TOOL_NONE SampleToolType = iota
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

const (
	SampleSoloMesh      = "Solo Mesh"
	SampleTileMesh      = "Tile Mesh"
	SampleTempObstacles = "Temp Obstacles"
)
const (
	HighLightTitleCacheDrawAreas = "Draw Areas"
	HighLightTitleDrawRegions    = "Draw Regions"
	HighLightTitleDrawContours   = "Draw Contours"
	HighLightTitleDrawMesh       = "Draw Mesh"
)
const (
	MeshObjExt = ".obj"
	MeshSetExt = ".gset"
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
