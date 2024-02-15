package config

import (
	"fmt"
	"github.com/gorustyt/fyne/v2/data/binding"
)

type Config struct {
	ToolsConfig *ToolsConfig
	PropsConfig *PropsConfig
}

func NewConfig() *Config {
	c := &Config{
		ToolsConfig: &ToolsConfig{},
		PropsConfig: &PropsConfig{
			VertLabelData:  binding.NewString(),
			BuildTimeLabel: binding.NewString(),
		},
	}
	c.PropsConfig.SetVertLabelData(0, 0)
	c.PropsConfig.SetBuildTimeLabelData(0)
	return c
}

type ToolsConfig struct {
	Uid string
	//TOOL_NAVMESH_TESTER
	PathFindToolMode string
	PathfindStraight string

	OnSetRandomStartClick         func()
	OnSetRandomEndClick           func()
	OnMakeRandomPointsAroundClick func()
	OnMakeRandomPointsClick       func()

	IncludeFlags []string
	ExcludeFlags []string
	//TOOL_NAVMESH_PRUNE

	//TOOL_OFFMESH_CONNECTION
	Bidir string
	//TOOL_CONVEX_VOLUME
	BoxHeight         float64
	BoxDescent        float64
	PolyOffset        float64
	AreaType          string
	OnClearShapeClick func()
	//TOOL_CROWD
	ToolModel               string
	ExpandOptions           []string
	ObstacleAvoidanceType   float64
	SeparationWeight        float64
	ExpandSelectedDebugDraw []string
	ExpandDebugDraw         []string
	//TOOL_CreateTiles
	OnCreateTilesCreateAllClick func()
	OnCreateTilesRemoveAllClick func()
	//CreateTempObstacles
	OneTempObstaclesRemoveAllClick func()
	//HighlightTileCache
	HighlightDrawType string
}

type PropsConfig struct {
	ShowLogOrShowTool string
	SampleType        string
	InputMeshPath     string

	VertLabelData binding.String

	RasterizationCellSize   float64
	RasterizationCellHeight float64

	AgentHeight   float64
	AgentRadius   float64
	AgentMaxClimb float64
	AgentMaxSlope float64

	RegionMinRegionSize    float64
	RegionMergedRegionSize float64

	Partitioning string
	Filtering    []string

	PolygonizationMaxEdgeLength float64
	PolygonizationMaxEdgeError  float64
	PolygonizationVertsPerPoly  float64

	DetailMeshSampleDistance       float64
	DetailMeshSampleMaxSampleError float64

	TileSize         float64
	KeepInterResults []string
	OnSaveClick      func()
	OnLoadClick      func()
	BuildTimeLabel   binding.String
	OnBuildClick     func()
	DrawMode         string
}

func (cfg *PropsConfig) SetVertLabelData(vertCount, triCount float64) {
	cfg.VertLabelData.Set(fmt.Sprintf("Verts: %.1fk  Tris: %.1fk", vertCount/1000.0, triCount/1000.0))
}

func (cfg *PropsConfig) SetBuildTimeLabelData(totalBuildTimeMs float64) {
	cfg.VertLabelData.Set(fmt.Sprintf("Build Time: %.1fms", totalBuildTimeMs))
}
