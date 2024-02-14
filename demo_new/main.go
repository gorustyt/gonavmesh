package main

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/app"
	"github.com/gorustyt/fyne/v2/theme"
	"gonavamesh/demo_new/ui"
)

func main() {
	size := fyne.NewSize(1200, 900)
	a := app.NewWithID("recast")
	w := a.NewWindow("recast")
	a.Settings().SetTheme(theme.DarkTheme())
	w.SetContent(ui.GetMenu())
	w.Resize(size)
	w.ShowAndRun()

}
