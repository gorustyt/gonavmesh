package ui

import (
	"github.com/gorustyt/fyne/v2"
	"gonavamesh/demo_new/config"
)

type ToolsPruneNavmesh struct {
}

func NewToolsPruneNavmesh(cfg *config.Config) *ToolsPruneNavmesh {
	return &ToolsPruneNavmesh{}
}
func (t *ToolsPruneNavmesh) GetRenderObjs() (res []fyne.CanvasObject) {
	return
}
