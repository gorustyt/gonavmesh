package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type Props struct {
	c   *fyne.Container
	cfg *config.Config
}

func NewProps(cfg *config.Config) *Props {
	p := &Props{cfg: cfg}

	p.c = container.NewVBox(
		widget.NewLabel("prop"),
		widget.NewRadioGroup([]string{
			config.ShowLog,
			config.ShowTools,
		}, func(s string) {

		}),
	)
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
func (p *Props) GetSample() (res []fyne.CanvasObject) {
	return []fyne.CanvasObject{
		widget.NewLabel("Sample"),
		widget.NewSelect([]string{
			config.SampleSoloMesh,
			config.SampleTileMesh,
			config.SampleTempObstacles,
		}, func(s string) {

		}),
		widget.NewLabel("Input Mesh"),
		widget.NewLabelWithData(p.cfg.PropsConfig.VertLabelData),
		widget.NewSelect([]string{}, func(s string) {

		}),
	}
}
func (p *Props) GetRenderObj() fyne.CanvasObject {
	s := container.NewVScroll(p.c)
	s.SetMinSize(fyne.NewSize(100, 900))
	return s
}

func (p *Props) GetRasterization() (res []fyne.CanvasObject) {
	s1 := widget.NewSlider(0.1, 1)
	s1.Step = 0.01
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSlider(0.1, 1)
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
	s1 := widget.NewSlider(0.1, 5)
	s1.Step = 0.1
	s1.OnChanged = func(f float64) {

	}

	s2 := widget.NewSlider(0, 5)
	s2.Step = 0.1
	s2.OnChanged = func(f float64) {

	}

	s3 := widget.NewSlider(0.1, 5)
	s3.Step = 0.1
	s3.OnChanged = func(f float64) {

	}

	s4 := widget.NewSlider(0.0, 90)
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
	s1 := widget.NewSlider(0, 150)
	s1.Step = 1
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSlider(0, 150)
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

		}),
	}
}

func (p *Props) GetPolygonization() (res []fyne.CanvasObject) {
	s1 := widget.NewSlider(0, 50)
	s1.Step = 1
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSlider(0.1, 3)
	s2.Step = 0.1
	s2.OnChanged = func(f float64) {

	}

	s3 := widget.NewSlider(3, 12)
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
	s1 := widget.NewSlider(0, 16)
	s1.Step = 1
	s1.OnChanged = func(f float64) {

	}
	s2 := widget.NewSlider(0, 16)
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
func (p *Props) GetItermediateResults() (res []fyne.CanvasObject) {
	b1 := widget.NewButton("Save", func() {

	})
	b2 := widget.NewButton("Load", func() {

	})
	b3 := widget.NewButton("Build", func() {

	})
	b3.Importance = widget.HighImportance
	return []fyne.CanvasObject{
		widget.NewRadioGroup([]string{
			config.KeepItermediateResults,
		}, func(strings string) {

		}),
		b1,
		b2,
		widget.NewSeparator(),
		b3,
	}
}
func (p *Props) GetDraw() (res []fyne.CanvasObject) {
	return []fyne.CanvasObject{
		widget.NewLabel("Draw"),
		widget.NewRadioGroup([]string{config.DrawInputMesh,
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
			config.DrawPolyMeshDetail}, func(s string) {

		}),
		widget.NewLabel("Tick 'Keep Itermediate Results' "),
		widget.NewLabel("to see more debug mode options."),
	}
}
