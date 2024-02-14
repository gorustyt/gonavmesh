package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"gonavamesh/demo_new/config"
)

func GetMenu() fyne.CanvasObject {
	InitToolsMap()
	return container.NewBorder(nil, nil, NewTools().GetRenderObj(), NewProps().GetRenderObj(), config.NewContent())
}
