package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
)

func GetMenu() fyne.CanvasObject {
	InitToolsMap()
	return container.NewBorder(nil, nil, NewTools().GetRenderObj(), NewContent(), NewProps().GetRenderObj())
}
