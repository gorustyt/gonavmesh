package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/widget"
)

const (
	SAMPLE_POLYAREA_GROUND = "Ground"
	SAMPLE_POLYAREA_WATER  = "Water"
	SAMPLE_POLYAREA_ROAD   = "Road"
	SAMPLE_POLYAREA_DOOR   = "Door"
	SAMPLE_POLYAREA_GRASS  = "Grass"
	SAMPLE_POLYAREA_JUMP   = "Jump"
)

type ToolsCreateConvexVolumes struct {
	c *fyne.Container
}

func NewToolsCreateConvexVolumes() *ToolsCreateConvexVolumes {
	button := widget.NewButton("Clear Shape", func() {

	})
	var f1 float64
	s1 := widget.NewSliderWithData(0.1, 20.0, binding.BindFloat(&f1))
	s1.Step = 0.1
	var f2 float64
	s2 := widget.NewSliderWithData(0.1, 20.0, binding.BindFloat(&f2))
	s2.Step = 0.1
	var f3 float64
	s3 := widget.NewSliderWithData(0.0, 20.0, binding.BindFloat(&f3))
	s3.Step = 0.1
	button.Importance = widget.DangerImportance
	group := widget.NewRadioGroup([]string{
		SAMPLE_POLYAREA_GROUND,
		SAMPLE_POLYAREA_WATER,
		SAMPLE_POLYAREA_ROAD,
		SAMPLE_POLYAREA_DOOR,
		SAMPLE_POLYAREA_GRASS,
		SAMPLE_POLYAREA_JUMP,
	}, func(s string) {

	})
	c := &ToolsCreateConvexVolumes{}
	c.c = container.NewVBox(
		widget.NewLabel("Shape Height"),
		s1,
		widget.NewLabel("Shape Descent"),
		s2,
		widget.NewLabel("Poly Offset"),
		s3,
		widget.NewSeparator(),
		widget.NewLabel("Area Type"),
		group,
		button,
	)
	return c
}
func (t *ToolsCreateConvexVolumes) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
