package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/widget"
	"github.com/gorustyt/gonavmesh/demo/config"
)

type ToolOffMeshConnection struct {
	c fyne.CanvasObject
}

func NewToolOffMeshConnection(ctx *Context) *ToolOffMeshConnection {
	c := widget.NewRadioGroup([]string{
		config.OneWay,
		config.Bidirectional,
	}, func(s string) {
		ctx.Config().ToolsConfig.Bidir = s
	})
	c.Selected = config.Bidirectional
	return &ToolOffMeshConnection{
		c: c,
	}
}
func (t *ToolOffMeshConnection) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
