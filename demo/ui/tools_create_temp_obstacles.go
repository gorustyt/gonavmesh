package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
)

type ToolsCreateTempObstacles struct {
	c *fyne.Container
}

func NewToolsToolsCreateTempObstacles(ctx *Context) *ToolsHighlightTileCache {
	c := &ToolsHighlightTileCache{}
	b := widget.NewButton("Remove All", ctx.Config().ToolsConfig.OneTempObstaclesRemoveAllClick)
	c.c = container.NewVBox(
		widget.NewLabel("Create Temp Obstacles"),
		b,
		widget.NewLabel("Click LMB to create an obstacle."),
		widget.NewLabel("Shift+LMB to remove an obstacle."),
	)
	return c
}
func (t ToolsCreateTempObstacles) GetRenderObjs() []fyne.CanvasObject {
	return []fyne.CanvasObject{
		t.c,
	}
}
