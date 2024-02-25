package config

import (
	"fmt"
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/gonavmesh/detour_crowd"
	"io/fs"
	"path/filepath"
)

type Config struct {
	ToolsConfig *ToolsConfig
	PropsConfig *PropsConfig
}

func (cfg *Config) Reset() {
	cfg.ToolsConfig.Reset()
	cfg.PropsConfig.Reset()
}
func NewConfig() *Config {
	c := &Config{
		ToolsConfig: &ToolsConfig{},
		PropsConfig: &PropsConfig{
			VertLabelData:                    binding.NewString(),
			BuildTimeLabel:                   binding.NewString(),
			TileSizeLabel:                    binding.NewString(),
			TileSizeMaxTitlesLabel:           binding.NewString(),
			TileSizeMaxPolysLabel:            binding.NewString(),
			TitleCacheLayersLabel:            binding.NewString(),
			TitleCacheLayerPerTileLabel:      binding.NewString(),
			TitleCacheMemoryLabel:            binding.NewString(),
			TitleCacheNavmeshBuildTimeLabel:  binding.NewString(),
			TitleCacheBuildPeakMemUsageLabel: binding.NewString(),
		},
	}
	c.Reset()
	return c
}

type ToolsConfig struct {
	Uid         string
	OnUidChange func()
	//TOOL_NAVMESH_TESTER
	PathFindToolMode string
	PathfindStraight string

	OnSetRandomStartClick         func()
	OnSetRandomEndClick           func()
	OnMakeRandomPointsAroundClick func()
	OnMakeRandomPointsClick       func()

	IncludeFlags  []string
	ExcludeFlags  []string
	OnFlagsChange func()
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
	OnToolModelChange       func()
	ExpandOptions           []string
	ExpandOptionsOnchange   func()
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

func (cfg *ToolsConfig) GetBidir() bool {
	if cfg.Bidir == OneWay {
		return false
	}
	if cfg.Bidir == Bidirectional {
		return true
	}
	return false
}
func (cfg *ToolsConfig) GetShowCorners() bool {
	return cfg.HasExpandOptionsChecked(ExpandSelectedDebugDrawShowCorners)
}
func (cfg *ToolsConfig) GetShowCollisionSegments() bool {
	return cfg.HasExpandOptionsChecked(ExpandSelectedDebugDrawShowCollisionSegs)
}
func (cfg *ToolsConfig) GetShowPath() bool {
	return cfg.HasExpandOptionsChecked(ExpandSelectedDebugDrawShowPath)
}
func (cfg *ToolsConfig) GetShowVO() bool {
	return cfg.HasExpandOptionsChecked(ExpandSelectedDebugDrawShowVO)
}
func (cfg *ToolsConfig) GetShowOpt() bool {
	return cfg.HasExpandOptionsChecked(ExpandSelectedDebugDrawShowPathOptimization)
}
func (cfg *ToolsConfig) GetShowNeis() bool {
	return cfg.HasExpandOptionsChecked(ExpandSelectedDebugDrawShowNeighbours)
}

func (cfg *ToolsConfig) GetShowLabels() bool {
	return cfg.HasExpandOptionsChecked(ExpandDebugDrawShowLabels)
}
func (cfg *ToolsConfig) GetShowGrid() bool {
	return cfg.HasExpandOptionsChecked(ExpandDebugDrawShowProxGrid)
}
func (cfg *ToolsConfig) GetShowNodes() bool {
	return cfg.HasExpandOptionsChecked(ExpandDebugDrawShowNodes)
}
func (cfg *ToolsConfig) GetShowPerfGraph() bool {
	return cfg.HasExpandOptionsChecked(ExpandDebugDrawShowPerfGraph)
}
func (cfg *ToolsConfig) GetShowDetailAll() bool {
	return cfg.HasExpandOptionsChecked(ExpandDebugDrawShowDetailAll)
}

func (cfg *ToolsConfig) GetAnticipateTurns() bool {
	return cfg.HasExpandOptionsChecked(ExpandOptionsAnticipateTurns)
}
func (cfg *ToolsConfig) GetOptimizeVis() bool {
	return cfg.HasExpandOptionsChecked(ExpandOptionsOptimizeVisibility)
}
func (cfg *ToolsConfig) GetOptimizeTopo() bool {
	return cfg.HasExpandOptionsChecked(ExpandOptionsOptimizeTopology)
}
func (cfg *ToolsConfig) GetObstacleAvoidance() bool {
	return cfg.HasExpandOptionsChecked(ExpandOptionsObstacleAvoidance)
}
func (cfg *ToolsConfig) GetSeparation() bool {
	return cfg.HasExpandOptionsChecked(ExpandOptionsSeparation)
}

func (cfg *ToolsConfig) HasExpandOptionsChecked(desc string) bool {

	for _, v := range cfg.ExpandDebugDraw {
		if v == desc {
			return true
		}
	}

	for _, v := range cfg.ExpandSelectedDebugDraw {
		if v == desc {
			return true
		}
	}

	for _, v := range cfg.ExpandOptions {
		if v == desc {
			return true
		}
	}
	return false
}

func (cfg *ToolsConfig) GetExpandOptionsUpdateFlags() int {
	updateFlags := 0
	for _, v := range cfg.ExpandOptions {
		switch v {
		case ExpandOptionsOptimizeVisibility:
			updateFlags |= detour_crowd.DT_CROWD_OPTIMIZE_VIS

		case ExpandOptionsOptimizeTopology:
			updateFlags |= detour_crowd.DT_CROWD_OPTIMIZE_TOPO

		case ExpandOptionsAnticipateTurns:
			updateFlags |= detour_crowd.DT_CROWD_ANTICIPATE_TURNS

		case ExpandOptionsObstacleAvoidance:
			updateFlags |= detour_crowd.DT_CROWD_OBSTACLE_AVOIDANCE
		case ExpandOptionsSeparation:
			updateFlags |= detour_crowd.DT_CROWD_SEPARATION

		}
	}
	return updateFlags
}

func (cfg *ToolsConfig) Reset() {
	cfg.AreaType = Desc_SAMPLE_POLYAREA_GRASS
	cfg.ExpandOptions = []string{
		ExpandOptionsOptimizeVisibility,
		ExpandOptionsOptimizeTopology,
		ExpandOptionsAnticipateTurns,
		ExpandOptionsObstacleAvoidance,
	}
	cfg.ToolModel = DescCrowdTool_TOOLMODE_CREATE
	cfg.Bidir = Bidirectional
	cfg.HighlightDrawType = HighLightTitleCacheDrawAreas
	cfg.ObstacleAvoidanceType = 3
	cfg.SeparationWeight = 2
	cfg.BoxHeight = 6
	cfg.BoxDescent = 1
	cfg.PathfindStraight = DT_STRAIGHTPATH_NONE_CROSSINGS
	cfg.IncludeFlags = []string{
		DESC_SAMPLE_POLYFLAGS_WALK,
		DESC_SAMPLE_POLYFLAGS_SWIM,
		DESC_SAMPLE_POLYFLAGS_DOOR,
		DESC_SAMPLE_POLYFLAGS_JUMP,
	}
}

type PropsConfig struct {
	ShowLogAndShowTool []string
	SampleType         string
	InputMeshLists     map[string]string
	InputMeshPath      string
	OnInputMesh        func()

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

	TileSize                         float64
	TileSizeLabel                    binding.String
	TileSizeMaxTitlesLabel           binding.String
	TileSizeMaxPolysLabel            binding.String
	TitleCacheLayersLabel            binding.String
	TitleCacheLayerPerTileLabel      binding.String
	TitleCacheMemoryLabel            binding.String
	TitleCacheNavmeshBuildTimeLabel  binding.String
	TitleCacheBuildPeakMemUsageLabel binding.String

	KeepInterResults []string
	OnSaveClick      func(writer fyne.URIWriteCloser)
	OnLoadClick      func(reader fyne.URIReadCloser)
	BuildTimeLabel   binding.String
	OnBuildClick     func()
	DrawMode         string
}

func (cfg *PropsConfig) Reset() {
	cfg.InputMeshLists = GetInputMeshList()
	cfg.Filtering = []string{
		FilteringLowHangingObstacles,
		FilteringLedgeSpans,
		FilteringWalkableLowHeightSpans,
	}
	cfg.Partitioning = DescSAMPLE_PARTITION_WATERSHED
	cfg.SampleType = SampleSoloMesh
	cfg.ShowLogAndShowTool = []string{ShowTools}
	cfg.RasterizationCellSize = 0.3
	cfg.RasterizationCellHeight = 0.2

	cfg.AgentHeight = 2
	cfg.AgentRadius = 0.6
	cfg.AgentMaxClimb = 0.9
	cfg.AgentMaxSlope = 45

	cfg.RegionMinRegionSize = 8
	cfg.RegionMergedRegionSize = 20

	cfg.PolygonizationMaxEdgeLength = 1
	cfg.PolygonizationMaxEdgeError = 1.3
	cfg.PolygonizationVertsPerPoly = 6.

	cfg.DetailMeshSampleDistance = 6
	cfg.DetailMeshSampleMaxSampleError = 1

	cfg.TileSize = 48
	cfg.SetVertLabelData(0, 0)
	cfg.SetBuildTimeLabelData(0)
	cfg.SetTileSizeLabel(0, 0)
	cfg.SetTileSizeMaxTitlesLabel(0)
	cfg.SetTileSizeMaxPolysLabel(0)
	cfg.SetTitleCacheLayersLabel(0)
	cfg.SetTitleCacheLayerPerTileLabel(0, 0)
	cfg.SetTitleCacheMemoryLabel(0, 0, 0)
	cfg.SetTitleCacheNavmeshBuildTimeLabel(0)
	cfg.SetTitleCacheBuildPeakMemUsageLabel(0)
}
func (cfg *PropsConfig) SetVertLabelData(vertCount, triCount float64) {
	cfg.VertLabelData.Set(fmt.Sprintf("Verts: %.1fk  Tris: %.1fk", vertCount/1000.0, triCount/1000.0))
}

func (cfg *PropsConfig) SetBuildTimeLabelData(totalBuildTimeMs float64) {
	cfg.VertLabelData.Set(fmt.Sprintf("Build Time: %.1fms", totalBuildTimeMs))
}

func (cfg *PropsConfig) SetTileSizeLabel(tw, th int) {
	cfg.TileSizeLabel.Set(fmt.Sprintf("Tiles  %d x %d", tw, th))
}
func (cfg *PropsConfig) SetTileSizeMaxTitlesLabel(maxTiles int) {
	cfg.TileSizeMaxTitlesLabel.Set(fmt.Sprintf("Max Tiles  %d", maxTiles))
}
func (cfg *PropsConfig) SetTileSizeMaxPolysLabel(maxPolysPerTile int) {
	cfg.TileSizeMaxPolysLabel.Set(fmt.Sprintf("Max Polys  %d", maxPolysPerTile))

}
func (cfg *PropsConfig) SetTitleCacheLayersLabel(cacheLayerCount int) {
	cfg.TitleCacheLayersLabel.Set(fmt.Sprintf("Layers  %d", cacheLayerCount))

}
func (cfg *PropsConfig) SetTitleCacheLayerPerTileLabel(cacheLayerCount float64, gridSize int) {
	cfg.TitleCacheLayerPerTileLabel.Set(fmt.Sprintf("Layers (per tile)  %.1f", cacheLayerCount/float64(gridSize)))

}
func (cfg *PropsConfig) SetTitleCacheMemoryLabel(cacheCompressedSize, cacheRawSize int, compressionRatio float64) {
	cfg.TitleCacheMemoryLabel.Set(fmt.Sprintf("Memory  %.1f kB / %.1f kB (%.1f%%)", float64(cacheCompressedSize)/1024.0, float64(cacheRawSize)/1024.0, compressionRatio*100.0))

}
func (cfg *PropsConfig) SetTitleCacheNavmeshBuildTimeLabel(cacheBuildTimeMs float64) {
	cfg.TitleCacheNavmeshBuildTimeLabel.Set(fmt.Sprintf("Navmesh Build Time  %.1f ms", cacheBuildTimeMs))

}
func (cfg *PropsConfig) SetTitleCacheBuildPeakMemUsageLabel(cacheBuildMemUsage int) {
	cfg.TitleCacheBuildPeakMemUsageLabel.Set(fmt.Sprintf("Build Peak Mem Usage  %.1f kB", float64(cacheBuildMemUsage)/1024.0))
}

func GetInputMeshList() map[string]string {
	res := map[string]string{}
	err := filepath.Walk("./demo/objs/Meshes", func(path string, info fs.FileInfo, err error) error {
		if err != nil {
			panic(err)
		}
		if filepath.Ext(path) == MeshObjExt || filepath.Ext(path) == MeshSetExt {
			res[info.Name()] = path
		}
		return nil
	})
	if err != nil {
		panic(err)
	}
	return res
}
