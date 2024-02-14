package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"gonavamesh/demo_new/mesh"
)

func GetMenu() fyne.CanvasObject {
	c := mesh.NewContent()
	cfg := c.GetConfig()
	InitToolsMap(cfg)
	return container.NewBorder(nil, nil, NewTools(cfg).GetRenderObj(), NewProps(cfg).GetRenderObj(), c)
}
