package gui

import (
	g "github.com/AllenDang/giu"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"github.com/go-gl/mathgl/mgl32"
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
	ui.mainWindow = g.NewMasterWindow(ui.title, ui.width, ui.height, g.MasterWindowFlagsNotResizable)
	ui.mainWindow.SetCloseCallback(ui.shutdown)
	ui.mainWindow.SetDropCallback(ui.fileDrop)
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
	cameraEulers := []float32{45, -45}
	cameraPos := []float64{0, 0, 0}
	camr := float32(1000.0)
	ui.mainWindow.Run(func() {
		gl.Viewport(0, 0, int32(ui.width), int32(ui.height))
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		gl.Enable(gl.BLEND)
		gl.BlendFunc(gl.SRC_ALPHA, gl.ONE_MINUS_SRC_ALPHA)
		gl.Disable(gl.TEXTURE_2D)
		gl.Enable(gl.DEPTH_TEST)

		// Compute the projection matrix.
		gl.MatrixLoadIdentityEXT(gl.PATH_PROJECTION_NV)
		mgl32.Perspective(50.0, float32(ui.width)/float32(ui.height), 1.0, camr)
		projectionMatrix := make([]float64, 16)
		gl.GetDoublev(gl.PATH_PROJECTION_NV, &projectionMatrix[0])
		// Compute the modelview matrix.
		gl.MatrixLoadIdentityEXT(gl.PATH_MODELVIEW_NV)
		gl.MatrixRotatefEXT(gl.PATH_MODELVIEW_NV, cameraEulers[0], 1, 0, 0)
		gl.MatrixRotatefEXT(gl.PATH_MODELVIEW_NV, cameraEulers[1], 0, 1, 0)
		gl.MatrixTranslatedEXT(gl.PATH_MODELVIEW_NV, -cameraPos[0], -cameraPos[1], -cameraPos[2])
		modelviewMatrix := make([]float64, 16)
		gl.GetDoublev(gl.PATH_MODELVIEW_NV, &modelviewMatrix[0])
		g.SingleWindow().Layout(ui.layout, ui.render)
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
