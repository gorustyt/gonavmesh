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
	c.PropsConfig.SetVertLabelData(0, 0)
	c.PropsConfig.SetBuildTimeLabelData(0)
	c.PropsConfig.SetTileSizeLabel(0, 0)
	c.PropsConfig.SetTileSizeMaxTitlesLabel(0)
	c.PropsConfig.SetTileSizeMaxPolysLabel(0)
	c.PropsConfig.SetTitleCacheLayersLabel(0)
	c.PropsConfig.SetTitleCacheLayerPerTileLabel(0, 0)
	c.PropsConfig.SetTitleCacheMemoryLabel(0, 0, 0)
	c.PropsConfig.SetTitleCacheNavmeshBuildTimeLabel(0)
	c.PropsConfig.SetTitleCacheBuildPeakMemUsageLabel(0)
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
