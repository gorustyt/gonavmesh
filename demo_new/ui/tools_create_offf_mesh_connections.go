package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type ToolOffMeshConnection struct {
	c fyne.CanvasObject
}

func NewToolOffMeshConnection(ctx *Context) *ToolOffMeshConnection {
	return &ToolOffMeshConnection{
		c: widget.NewRadioGroup([]string{
			config.OneWay,
			config.Bidirectional,
		}, func(s string) {
			ctx.Config().ToolsConfig.Bidir = s
		}),
	}
}
func (t *ToolOffMeshConnection) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
