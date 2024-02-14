package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type ToolsCreateTitles struct {
	c *fyne.Container
}

func NewToolsCreateTitles(cfg *config.Config) *ToolsCreateTitles {
	t := &ToolsCreateTitles{}
	b1 := widget.NewButton("Create All", func() {

	})
	b1.Importance = widget.SuccessImportance
	b2 := widget.NewButton("Remove All", func() {

	})
	b2.Importance = widget.DangerImportance
	t.c = container.NewVBox(
		widget.NewLabel("Create Tiles"),
		b1, b2,
	)
	return t
}

func (t *ToolsCreateTitles) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
