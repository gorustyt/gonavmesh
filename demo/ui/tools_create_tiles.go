package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
)

type ToolsCreateTitles struct {
	c *fyne.Container
}

func NewToolsCreateTitles(ctx *Context) *ToolsCreateTitles {
	t := &ToolsCreateTitles{}
	b1 := widget.NewButton("Create All", ctx.Config().ToolsConfig.OnCreateTilesCreateAllClick)
	b1.Importance = widget.SuccessImportance
	b2 := widget.NewButton("Remove All", ctx.Config().ToolsConfig.OnCreateTilesRemoveAllClick)
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
