package gui

import (
	g "github.com/AllenDang/giu"
	"github.com/go-gl/glfw/v3.3/glfw"
	"log"
	"math"
)

func init() {
	err := glfw.Init()
	if err != nil {
		panic(err)
	}
}

type Gui struct {
	title         string
	mainWindow    *g.MasterWindow
	gs            *guiState
	render        *imguiRender
	mesh          *mesh
	layout        *layout
	logger        *logger
	width, height int
	sample        *Sample
	InputGeom     *InputGeom
	font          *g.FontInfo
}

func NewGui(title string) *Gui {
	ui := &Gui{
		gs:        newGuiState(),
		title:     title,
		mesh:      mewMesh(),
		logger:    newLogger(),
		InputGeom: newInputGeom(),
		sample:    newSample(),
	}
	ui.width, ui.height = getWH()
	ui.layout = newLayout(ui)
	ui.render = newImguiRender(ui)
	ui.mainWindow = g.NewMasterWindow(ui.title, ui.width, ui.height, 0)
	ui.mainWindow.SetCloseCallback(ui.shutdown)
	ui.mainWindow.SetDropCallback(ui.fileDrop)
	ui.font = g.Context.FontAtlas.AddFont("DroidSans.ttf", 20)
	ui.mainWindow.SetBgColor(imguiRGBA(169, 169, 169, 5))
	return ui
}

func (ui *Gui) shutdown() bool {
	log.Printf("ui %v shuttdown.", ui.title)
	return true
}

func (ui *Gui) fileDrop(files []string) {
	log.Printf("find file %v ", files)
}

func (ui *Gui) Run() {
	ui.mainWindow.Run(func() {
		g.SingleWindow().Layout(g.Custom(func() {
			g.PushFont(ui.font)
		}), ui.layout, ui.render, g.Custom(func() {
			g.PopFont()
		}))
	})
}

func getWH() (width, height int) {
	sw := glfw.GetPrimaryMonitor().GetVideoMode().Width
	sh := glfw.GetPrimaryMonitor().GetVideoMode().Height
	aspect := 16.0 / 9.0
	width = int(math.Min(float64(sw), float64(sh)*aspect)) - 80
	height = sh - 80
	return width, height
}
