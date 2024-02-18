package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/dialog"
	"github.com/gorustyt/fyne/v2/storage"
	"github.com/gorustyt/fyne/v2/theme"
	"github.com/gorustyt/fyne/v2/widget"
	"log"
	"net/url"
	"path/filepath"
)

func SetMainMenu(a fyne.App, w fyne.Window, ctx *Context) {
	openItem := fyne.NewMenuItem("Load", func() {
		fd := dialog.NewFileOpen(func(reader fyne.URIReadCloser, err error) {
			if err != nil {
				dialog.ShowError(err, w)
				return
			}
			if reader == nil {
				log.Println("Cancelled")
				return
			}
			ctx.Config().PropsConfig.OnLoadClick()
		}, w)
		fd.SetFilter(storage.NewExtensionFileFilter([]string{".obj"}))
		p, err := filepath.Abs("./demo/objs/Meshes")
		if err != nil {
			panic(err)
		}
		testData := storage.NewFileURI(p)
		dir, err := storage.ListerForURI(testData)
		if err != nil {
			panic(err)
		}
		fd.SetLocation(dir)
		fd.Show()
	})
	saveItem := fyne.NewMenuItem("Save", func() {
		fd := dialog.NewFileSave(func(writer fyne.URIWriteCloser, err error) {
			if err != nil {
				dialog.ShowError(err, w)
				return
			}
			if writer == nil {
				log.Println("Cancelled")
				return
			}

			ctx.Config().PropsConfig.OnSaveClick()
		}, w)
		fd.SetFilter(storage.NewExtensionFileFilter([]string{".obj"}))
		fd.SetFileName("xx.obj")
		fd.Show()
	})
	label := widget.NewLabel("setting theme")
	label.Alignment = fyne.TextAlignCenter
	themes := container.NewGridWithColumns(2,
		widget.NewButton("Dark", func() {
			fyne.CurrentApp().Settings().SetTheme(theme.DarkTheme())
		}),
		widget.NewButton("Light", func() {
			fyne.CurrentApp().Settings().SetTheme(theme.LightTheme())
		}),
	)
	themeItem := fyne.NewMenuItem("Theme", func() {
		w1 := a.NewWindow("Theme Settings")
		w1.SetContent(container.NewVBox(
			label,
			themes,
		))
		w1.Resize(fyne.NewSize(200, 200))
		w1.Show()
	})
	settingMenu := fyne.NewMenu("Settings", themeItem)

	buildItem := fyne.NewMenuItem("build", func() {
		ctx.Config().PropsConfig.OnBuildClick()
	})
	runMenu := fyne.NewMenu("Run", buildItem)
	helpMenu := fyne.NewMenu("Help",
		fyne.NewMenuItem("Documentation", func() {
			u, _ := url.Parse("https://github.com/gorustyt/gonavmesh")
			_ = a.OpenURL(u)
		}),
		fyne.NewMenuItemSeparator(),
	)
	// a quit item will be appended to our first (File) menu
	file := fyne.NewMenu("File", openItem, saveItem)
	main := fyne.NewMainMenu(
		file,
		runMenu,
		settingMenu,
		helpMenu,
	)
	w.SetMainMenu(main)
}
