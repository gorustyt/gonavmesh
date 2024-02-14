package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type ToolsCreateConvexVolumes struct {
	c *fyne.Container
}

func NewToolsCreateConvexVolumes(cfg *config.Config) *ToolsCreateConvexVolumes {
	button := widget.NewButton("Clear Shape", cfg.ToolsConfig.OnClearShapeClick)
	s1 := widget.NewSliderWithData(0.1, 20.0, binding.BindFloat(&cfg.ToolsConfig.BoxHeight))
	s1.Step = 0.1
	s2 := widget.NewSliderWithData(0.1, 20.0, binding.BindFloat(&cfg.ToolsConfig.BoxDescent))
	s2.Step = 0.1
	s3 := widget.NewSliderWithData(0.0, 20.0, binding.BindFloat(&cfg.ToolsConfig.PolyOffset))
	s3.Step = 0.1
	button.Importance = widget.DangerImportance
	group := widget.NewRadioGroup([]string{
		config.SAMPLE_POLYAREA_GROUND,
		config.SAMPLE_POLYAREA_WATER,
		config.SAMPLE_POLYAREA_ROAD,
		config.SAMPLE_POLYAREA_DOOR,
		config.SAMPLE_POLYAREA_GRASS,
		config.SAMPLE_POLYAREA_JUMP,
	}, func(s string) {
		cfg.ToolsConfig.AreaType = s
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
