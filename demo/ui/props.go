package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/widget"
	"github.com/gorustyt/gonavmesh/demo/config"
)

type Props struct {
	c   *fyne.Container
	ctx *Context

	drawGroup        *widget.RadioGroup
	itermediateGroup *widget.CheckGroup
}

func NewProps(ctx *Context) *Props {
	p := &Props{ctx: ctx}
	ctx.AppendSampleChange(p)
	g := widget.NewCheckGroup([]string{
		config.ShowLog,
		config.ShowTools,
	}, func(s []string) {
		ctx.Config().PropsConfig.ShowLogAndShowTool = s
	})
	g.Selected = ctx.Config().PropsConfig.ShowLogAndShowTool
	p.c = container.NewVBox(
		g,
	)

	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetSample() {
		p.c.Add(v)
	}
	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetRasterization() {
		p.c.Add(v)
	}

	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetAgent() {
		p.c.Add(v)
	}
	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetPartitioning() {
		p.c.Add(v)
	}
	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetFiltering() {
		p.c.Add(v)
	}
	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetPolygonization() {
		p.c.Add(v)
	}

	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetDetailMesh() {
		p.c.Add(v)
	}

	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetItermediateResults() {
		p.c.Add(v)
	}
	p.c.Add(widget.NewSeparator())
	for _, v := range p.GetDraw() {
		p.c.Add(v)
	}
	return p
}

func (p *Props) SampleChange(sample string) {
	p.setDraw(sample)
	p.setItermediateGroup(sample)
}

func (p *Props) GetSample() (res []fyne.CanvasObject) {
	s := widget.NewSelect([]string{
		config.SampleSoloMesh,
		config.SampleTileMesh,
		config.SampleTempObstacles,
	}, func(s string) {
		p.ctx.Config().PropsConfig.SampleType = s
		p.ctx.OnSampleChange(s)
	})
	s.Selected = p.ctx.Config().PropsConfig.SampleType
	p.ctx.AppendAfterInit(func() {
		p.ctx.OnSampleChange(s.Selected)
		p.ctx.Config().PropsConfig.OnInputMesh()
	})
	var ss []string
	for k := range p.ctx.Config().PropsConfig.InputMeshLists {
		ss = append(ss, k)
	}
	return []fyne.CanvasObject{
		widget.NewLabel("Sample"),
		s,
		widget.NewLabel("Input Mesh"),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.VertLabelData),
		widget.NewSelect(ss, func(s string) {
			p.ctx.Config().PropsConfig.InputMeshPath = s
			p.ctx.Config().PropsConfig.OnInputMesh()
		}),
	}
}
func (p *Props) GetRenderObj() fyne.CanvasObject {
	s := container.NewVScroll(p.c)
	s.SetMinSize(fyne.NewSize(100, 780))
	return container.NewVBox(widget.NewLabel("Properties"), s)
}

func (p *Props) GetRasterization() (res []fyne.CanvasObject) {
	s1 := widget.NewSliderWithData(0.1, 1, binding.BindFloat(&p.ctx.Config().PropsConfig.RasterizationCellSize))
	s1.Step = 0.01
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSliderWithData(0.1, 1, binding.BindFloat(&p.ctx.Config().PropsConfig.RasterizationCellHeight))
	s2.Step = 0.01
	s2.OnChanged = func(f float64) {

	}
	return []fyne.CanvasObject{
		widget.NewLabel("Rasterization"),
		widget.NewLabel("Cell Size"),
		PackSlider(s1),
		widget.NewLabel("Cell Height"),
		PackSlider(s2),
	}
}
func (p *Props) GetAgent() (res []fyne.CanvasObject) {
	s1 := widget.NewSliderWithData(0.1, 5, binding.BindFloat(&p.ctx.Config().PropsConfig.AgentHeight))
	s1.Step = 0.1
	s1.OnChanged = func(f float64) {

	}

	s2 := widget.NewSliderWithData(0, 5, binding.BindFloat(&p.ctx.Config().PropsConfig.AgentRadius))
	s2.Step = 0.1
	s2.OnChanged = func(f float64) {

	}

	s3 := widget.NewSliderWithData(0.1, 5, binding.BindFloat(&p.ctx.Config().PropsConfig.AgentMaxClimb))
	s3.Step = 0.1
	s3.OnChanged = func(f float64) {

	}

	s4 := widget.NewSliderWithData(0.0, 90, binding.BindFloat(&p.ctx.Config().PropsConfig.AgentMaxSlope))
	s4.Step = 1
	s4.OnChanged = func(f float64) {

	}
	return []fyne.CanvasObject{
		widget.NewLabel("Height"),
		PackSlider(s1),
		widget.NewLabel("Radius"),
		PackSlider(s2),
		widget.NewLabel("Max Climb"),
		PackSlider(s3),
		widget.NewLabel("Max Slope"),
		PackSlider(s4),
	}
}

func (p *Props) GetRegion() (res []fyne.CanvasObject) {
	s1 := widget.NewSliderWithData(0, 150, binding.BindFloat(&p.ctx.Config().PropsConfig.RegionMinRegionSize))
	s1.Step = 1
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSliderWithData(0, 150, binding.BindFloat(&p.ctx.Config().PropsConfig.RegionMergedRegionSize))
	s2.Step = 1
	s2.OnChanged = func(f float64) {

	}
	return []fyne.CanvasObject{
		widget.NewLabel("Region"),
		widget.NewLabel("Min Region Size"),
		PackSlider(s1),
		widget.NewLabel("Merged Region Size"),
		PackSlider(s2),
	}
}

func (p *Props) GetPartitioning() (res []fyne.CanvasObject) {
	g := widget.NewRadioGroup([]string{
		config.DescSAMPLE_PARTITION_WATERSHED,
		config.DescSAMPLE_PARTITION_MONOTONE,
		config.DescSAMPLE_PARTITION_LAYERS,
	}, func(s string) {
		p.ctx.Config().PropsConfig.Partitioning = s
	})
	g.Selected = p.ctx.Config().PropsConfig.Partitioning
	return []fyne.CanvasObject{
		widget.NewLabel("Partitioning"),
		g,
	}
}

func (p *Props) GetFiltering() (res []fyne.CanvasObject) {
	g := widget.NewCheckGroup([]string{
		config.FilteringLowHangingObstacles,
		config.FilteringLedgeSpans,
		config.FilteringWalkableLowHeightSpans,
	}, func(s []string) {
		p.ctx.Config().PropsConfig.Filtering = s
	})
	g.Selected = p.ctx.Config().PropsConfig.Filtering
	return []fyne.CanvasObject{
		widget.NewLabel("Filtering"),
		g,
	}
}

func (p *Props) GetPolygonization() (res []fyne.CanvasObject) {
	s1 := widget.NewSliderWithData(0, 50, binding.BindFloat(
		&p.ctx.Config().PropsConfig.PolygonizationMaxEdgeLength))
	s1.Step = 1
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSliderWithData(0.1, 3, binding.BindFloat(
		&p.ctx.Config().PropsConfig.PolygonizationMaxEdgeError))
	s2.Step = 0.1
	s2.OnChanged = func(f float64) {

	}

	s3 := widget.NewSliderWithData(3, 12, binding.BindFloat(
		&p.ctx.Config().PropsConfig.PolygonizationVertsPerPoly))
	s3.Step = 1
	s3.OnChanged = func(f float64) {

	}
	return []fyne.CanvasObject{
		widget.NewLabel("Polygonization"),
		widget.NewLabel("Max Edge Length"),
		PackSlider(s1),
		widget.NewLabel("Max Edge Error"),
		PackSlider(s2),
		widget.NewLabel("Verts Per Poly"),
		PackSlider(s3),
	}
}

func (p *Props) GetDetailMesh() (res []fyne.CanvasObject) {
	s1 := widget.NewSliderWithData(0, 16, binding.BindFloat(
		&p.ctx.Config().PropsConfig.DetailMeshSampleDistance))
	s1.Step = 1
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSliderWithData(0, 16, binding.BindFloat(
		&p.ctx.Config().PropsConfig.DetailMeshSampleMaxSampleError))
	s2.Step = 1
	s2.OnChanged = func(f float64) {

	}
	return []fyne.CanvasObject{
		widget.NewLabel("Detail Mesh"),
		widget.NewLabel("Sample Distance"),
		PackSlider(s1),
		widget.NewLabel("Max Sample Error"),
		PackSlider(s2),
	}
}
func (p *Props) setItermediateGroup(sample string) {
	switch sample {
	case config.SampleSoloMesh:
		p.itermediateGroup.Selected = []string{config.KeepItermediateResults}
		p.itermediateGroup.Options = []string{config.KeepItermediateResults}
	case config.SampleTileMesh:
		p.itermediateGroup.Selected = []string{config.KeepBuildAllTiles}
		p.itermediateGroup.Options = []string{config.KeepItermediateResults, config.KeepBuildAllTiles}
	case config.SampleTempObstacles:
		p.itermediateGroup.Selected = []string{}
		p.itermediateGroup.Options = []string{config.KeepItermediateResults}
	}

	p.itermediateGroup.Refresh()
}

func (p *Props) GetTileCache() (res []fyne.CanvasObject) {
	vs := []*widget.Label{
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TitleCacheLayersLabel),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TitleCacheLayerPerTileLabel),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TitleCacheMemoryLabel),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TitleCacheNavmeshBuildTimeLabel),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TitleCacheBuildPeakMemUsageLabel),
	}

	res = append(res, widget.NewLabel("Tile Cache"))
	for _, v := range vs {
		v.Alignment = fyne.TextAlignTrailing
		v.TextStyle.Monospace = true
		v.TextStyle.Italic = true
		res = append(res, v)
	}
	res = append(res, widget.NewSeparator())
	return res
}
func (p *Props) GetTiling() (res []fyne.CanvasObject) {
	s := widget.NewSliderWithData(16.0, 1024.0, binding.BindFloat(&p.ctx.Config().PropsConfig.TileSize))
	s.Step = 16
	vs := []*widget.Label{
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TileSizeLabel),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TileSizeMaxTitlesLabel),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.TileSizeMaxPolysLabel),
	}
	res = []fyne.CanvasObject{
		widget.NewLabel("Tiling"),
		widget.NewLabel("TileSize"),
		PackSlider(s),
		widget.NewSeparator(),
	}
	for _, v := range vs {
		v.Alignment = fyne.TextAlignTrailing
		v.TextStyle.Monospace = true
		res = append(res, v)
	}
	res = append(res, widget.NewSeparator())
	return res
}

func (p *Props) GetItermediateResults() (res []fyne.CanvasObject) {
	p.itermediateGroup = widget.NewCheckGroup([]string{}, func(strings []string) {
		p.ctx.Config().PropsConfig.KeepInterResults = strings
	})
	res = append(res, p.itermediateGroup)
	res = append(res, widget.NewSeparator())
	for _, v := range p.GetTiling() {
		res = append(res, v)
		p.ctx.AppendShow(config.SampleTileMesh, v)
		p.ctx.AppendShow(config.SampleTempObstacles, v)
	}
	for _, v := range p.GetTileCache() {
		res = append(res, v)
		p.ctx.AppendShow(config.SampleTempObstacles, v)
	}
	return append(res, []fyne.CanvasObject{
		widget.NewSeparator(),
	}...)
}

func (p *Props) setDraw(sample string) {
	switch sample {
	case config.SampleSoloMesh:
		p.drawGroup.Options = []string{
			config.DrawInputMesh,
			config.DrawNavmesh,
			config.DrawNavmeshInvis,
			config.DrawNavmeshTrans,
			config.DrawNavmeshBVTree,
			config.DrawNavmeshNodes,
			config.DrawVoxels,
			config.DrawWalkableVoxels,
			config.DrawCompact,
			config.DrawCompactDistance,
			config.DrawCompactRegions,
			config.DrawRegionConnections,
			config.DrawRawContours,
			config.DrawBothContours,
			config.DrawContours,
			config.DrawPolyMesh,
			config.DrawPolyMeshDetail,
		}
	case config.SampleTileMesh:
		p.drawGroup.Options = []string{
			config.DrawInputMesh,
			config.DrawNavmesh,
			config.DrawNavmeshInvis,
			config.DrawNavmeshTrans,
			config.DrawNavmeshBVTree,
			config.DrawNavmeshNodes,
			config.DrawPortals,
			config.DrawVoxels,
			config.DrawWalkableVoxels,
			config.DrawCompact,
			config.DrawCompactDistance,
			config.DrawCompactRegions,
			config.DrawRegionConnections,
			config.DrawRawContours,
			config.DrawBothContours,
			config.DrawContours,
			config.DrawPolyMesh,
			config.DrawPolyMeshDetail,
		}
	case config.SampleTempObstacles:
		p.drawGroup.Options = []string{
			config.DrawInputMesh,
			config.DrawNavmesh,
			config.DrawNavmeshInvis,
			config.DrawNavmeshTrans,
			config.DrawNavmeshBVTree,
			config.DrawNavmeshNodes,
			config.DrawPortals,
			config.DrawCacheBounds,
		}
	}
	p.drawGroup.Selected = config.DrawNavmesh
	p.drawGroup.Refresh()
}
func (p *Props) GetDraw() (res []fyne.CanvasObject) {
	p.drawGroup = widget.NewRadioGroup([]string{}, func(s string) {
		p.ctx.Config().PropsConfig.DrawMode = s
	})
	return []fyne.CanvasObject{
		widget.NewLabel("Draw"),
		p.drawGroup,
		widget.NewLabel("Tick 'Keep Itermediate Results' "),
		widget.NewLabel("to see more debug mode options."),
	}
}
