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
	Desc_OOLMODE_PATHFIND_STRAIGHT         = "Pathfind Straight"
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
