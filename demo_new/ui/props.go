package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
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
	p.c = container.NewVBox(
		widget.NewLabel("prop"),
		widget.NewRadioGroup([]string{
			config.ShowLog,
			config.ShowTools,
		}, func(s string) {
			ctx.Config().PropsConfig.ShowLogOrShowTool = s
		}),
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
	s.Selected = config.SampleSoloMesh
	p.ctx.AppendAfterInit(func() {
		p.ctx.OnSampleChange(s.Selected)
	})
	return []fyne.CanvasObject{
		widget.NewLabel("Sample"),
		s,
		widget.NewLabel("Input Mesh"),
		widget.NewLabelWithData(p.ctx.Config().PropsConfig.VertLabelData),
		widget.NewSelect([]string{}, func(s string) {
			p.ctx.Config().PropsConfig.InputMeshPath = s
		}),
	}
}
func (p *Props) GetRenderObj() fyne.CanvasObject {
	s := container.NewVScroll(p.c)
	s.SetMinSize(fyne.NewSize(100, 900))
	return s
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
		s1,
		widget.NewLabel("Cell Height"),
		s2,
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
		s1,
		widget.NewLabel("Radius"),
		s2,
		widget.NewLabel("Max Climb"),
		s3,
		widget.NewLabel("Max Slope"),
		s4,
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
		s1,
		widget.NewLabel("Merged Region Size"),
		s2,
	}
}

func (p *Props) GetPartitioning() (res []fyne.CanvasObject) {
	return []fyne.CanvasObject{
		widget.NewLabel("Partitioning"),
		widget.NewRadioGroup([]string{
			config.DescSAMPLE_PARTITION_WATERSHED,
			config.DescSAMPLE_PARTITION_MONOTONE,
			config.DescSAMPLE_PARTITION_LAYERS,
		}, func(s string) {
			p.ctx.Config().PropsConfig.Partitioning = s
		}),
	}
}

func (p *Props) GetFiltering() (res []fyne.CanvasObject) {
	return []fyne.CanvasObject{
		widget.NewLabel("Filtering"),
		widget.NewCheckGroup([]string{
			config.FilteringLowHangingObstacles,
			config.FilteringLedgeSpans,
			config.FilteringWalkableLowHeightSpans,
		}, func(s []string) {
			p.ctx.Config().PropsConfig.Filtering = s
		}),
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
		s1,
		widget.NewLabel("Max Edge Error"),
		s2,
		widget.NewLabel("Verts Per Poly"),
		s3,
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
		s1,
		widget.NewLabel("Max Sample Error"),
		s2,
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

func (p *Props) GetItermediateResults() (res []fyne.CanvasObject) {
	b1 := widget.NewButton("Save", p.ctx.Config().PropsConfig.OnSaveClick)
	b1.Importance = widget.SuccessImportance
	b2 := widget.NewButton("Load", p.ctx.Config().PropsConfig.OnLoadClick)
	b3 := widget.NewButton("Build", p.ctx.Config().PropsConfig.OnBuildClick)
	b3.Importance = widget.HighImportance
	p.itermediateGroup = widget.NewCheckGroup([]string{}, func(strings []string) {
		p.ctx.Config().PropsConfig.KeepInterResults = strings
	})
	return []fyne.CanvasObject{
		p.itermediateGroup,
		b1,
		b2,
		widget.NewSeparator(),
		b3,
	}
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
