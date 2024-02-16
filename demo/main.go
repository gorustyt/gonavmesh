package main

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/app"
	"github.com/gorustyt/fyne/v2/theme"
	"github.com/gorustyt/gonavmesh/demo/ui"
)

func main() {
	size := fyne.NewSize(1200, 900)
	a := app.NewWithID("recast")
	w := a.NewWindow("recast")
	a.Settings().SetTheme(theme.DarkTheme())
	ui.SetUi(a, w)
	w.Resize(size)
	w.ShowAndRun()

}
