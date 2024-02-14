package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/widget"
)

type ToolOffMeshConnection struct {
	c fyne.CanvasObject
}

func NewToolOffMeshConnection() *ToolOffMeshConnection {
	return &ToolOffMeshConnection{
		c: widget.NewRadioGroup([]string{
			"One Way",
			"Bidirectional",
		}, func(s string) {

		}),
	}
}
func (t *ToolOffMeshConnection) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
