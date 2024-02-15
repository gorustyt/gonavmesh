package ui

import (
	"github.com/gorustyt/fyne/v2"
)

type ToolsPruneNavmesh struct {
}

func NewToolsPruneNavmesh(ctx *Context) *ToolsPruneNavmesh {
	return &ToolsPruneNavmesh{}
}
func (t *ToolsPruneNavmesh) GetRenderObjs() (res []fyne.CanvasObject) {
	return
}
