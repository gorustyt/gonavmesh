package mesh

import (
	"github.com/go-gl/mathgl/mgl32"
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/canvas3d/canvas3d_render"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/widget"
	"github.com/gorustyt/gonavmesh/demo/config"
	"log/slog"
)

type ContentRender struct {
	c *Content
}

func (c *ContentRender) Destroy() {

}

func (c *ContentRender) Layout(size fyne.Size) {
	c.c.Resize(size)
}

func (c *ContentRender) MinSize() fyne.Size {
	return c.c.MinSize()
}

func (c *ContentRender) Objects() []fyne.CanvasObject {
	return []fyne.CanvasObject{c.c.canvas3D.GetRenderObj()}
}

func (c *ContentRender) Refresh() {
	c.c.Refresh()
}

type Content struct {
	widget.BaseWidget
	labelC     fyne.CanvasObject
	label1Data binding.String
	label2Data binding.String
	canvas3D   canvas3d_render.ICanvas3d

	cfg    *config.Config
	sample ISample
	geom   *InputGeom

	state *ContentState
}

func (c *Content) CreateRenderer() fyne.WidgetRenderer {
	return &ContentRender{c: c}
}

func (c *Content) SetLabel1(s string) {
	c.label1Data.Set(s)
}

func (c *Content) SetLabel2(s string) {
	c.label2Data.Set(s)
}

func (c *Content) GetConfig() *config.Config {
	return c.cfg
}

func (c *Content) MinSize() fyne.Size {
	return c.state.GetSize()
}

func (c *Content) Move(position fyne.Position) {

}

func (c *Content) Position() fyne.Position {
	return fyne.Position{}
}

func (c *Content) Resize(size fyne.Size) {

}

func (c *Content) Size() fyne.Size {
	return fyne.Size{}
}

func (c *Content) Hide() {

}

func (c *Content) Visible() bool {
	return true
}

func (c *Content) Show() {

}

func (c *Content) Refresh() {
	c.canvas3D.Reset()
	c.sample.HandleRender()
	c.canvas3D.Refresh()
}

func (c *Content) InputMeshChange() {
	geom := newInputGeom()
	if !geom.load(c.cfg.PropsConfig.InputMeshLists[c.cfg.PropsConfig.InputMeshPath]) {
		slog.Info("geom load error")
		return
	}
	c.geom = geom
	if c.sample != nil {
		c.sample.handleMeshChanged(geom)
	}
	c.Refresh()
}

func (c *Content) ToolChange() {
	switch c.cfg.ToolsConfig.Uid {
	case config.Desc_TOOL_NAVMESH_TESTER:
		//c.sample.setTool(newNavMeshTesterTool(c))
	case config.Desc_TOOL_NAVMESH_PRUNE:
		//c.sample.setTool(newTempObstacleHilightTool(c))
	case config.Desc_TOOL_OFFMESH_CONNECTION:
		//c.sample.setTool(newOffMeshConnectionTool(c))
	case config.Desc_TOOL_OFFMESH_Links:
		//c.sample.setTool(newOffMeshConnectionTool(c))
	case config.Desc_TOOL_CONVEX_VOLUME:
		//c.sample.setTool(newConvexVolumeTool(c))
	case config.Desc_TOOL_CROWD:
		//c.sample.setTool(newCrowdTool(c))
	case config.Desc_TOOL_CreateTiles:
		//c.sample.setTool(newMeshTitleTool(c))
	case config.Desc_TOOL_CreateTempObstacles:
		//c.sample.setTool(newTempObstacleCreateTool(c))
	case config.Desc_TOOL_HighlightTileCache:

	}
}
func (c *Content) SampleChange(sample string) {
	switch sample {
	case config.SampleSoloMesh:
		c.sample = newSampleSoloMesh(c)
	case config.SampleTileMesh:
		//c.sample = newSampleTileMesh(c)
	case config.SampleTempObstacles:
		//c.sample = newSampleTempObstacles(c)
	}
}

func NewContent() *Content {
	cfg := config.NewConfig()
	c := &Content{
		state:      NewContentState(),
		label1Data: binding.NewString(),
		label2Data: binding.NewString(),
		canvas3D:   canvas3d_render.NewCanvas3d(),
		cfg:        cfg}
	cfg.PropsConfig.OnInputMesh = c.InputMeshChange
	cfg.ToolsConfig.OnUidChange = c.ToolChange
	c.ExtendBaseWidget(c)
	c.labelC = container.NewVBox(
		widget.NewLabelWithData(c.label1Data),
		widget.NewLabelWithData(c.label2Data),
	)
	return c
}

func (c *Content) refresh() {

}

type ContentState struct {
	camr     float32
	viewPort mgl32.Vec4
	model    mgl32.Mat4
	view     mgl32.Mat4
	project  mgl32.Mat4

	cameraEulers mgl32.Vec2
	cameraPos    mgl32.Vec3

	rayStart mgl32.Vec3
	rayEnd   mgl32.Vec3

	mousePos     mgl32.Vec2
	origMousePos mgl32.Vec2
}

func NewContentState() *ContentState {
	return &ContentState{
		camr:         1000,
		cameraEulers: mgl32.Vec2{45, -45},
		cameraPos:    mgl32.Vec3{},
		viewPort: mgl32.Vec4{
			0, 0, 600, 1080,
		},
	}
}

func (s *ContentState) GetSize() fyne.Size {
	return fyne.Size{s.viewPort[2], s.viewPort[3]}
}
func (s *ContentState) Update() {
	var err error
	size := s.GetSize()
	s.project = mgl32.Perspective(50, size.Width/size.Height, 1, s.camr)
	s.model = mgl32.Translate3D(
		-s.cameraPos[0],
		-s.cameraPos[1],
		-s.cameraPos[2],
	).Mul4(mgl32.HomogRotate3D(
		mgl32.DegToRad(s.cameraEulers[1]),
		mgl32.Vec3{0, 1, 0})).Mul4(
		mgl32.HomogRotate3D(
			mgl32.DegToRad(s.cameraEulers[0]),
			mgl32.Vec3{1, 0, 0}))
	s.rayStart, err = mgl32.UnProject(
		mgl32.Vec3{
			s.mousePos[0],
			s.mousePos[1],
			0,
		},
		s.model,
		s.project,
		int(s.viewPort[0]),
		int(s.viewPort[1]),
		int(s.viewPort[2]),
		int(s.viewPort[3]),
	)
	if err != nil {
		panic(err)
	}
	s.rayEnd, err = mgl32.UnProject(
		mgl32.Vec3{
			s.mousePos[0],
			s.mousePos[1],
			0,
		},
		s.model,
		s.project,
		int(s.viewPort[0]),
		int(s.viewPort[1]),
		int(s.viewPort[2]),
		int(s.viewPort[3]),
	)
	if err != nil {
		panic(err)
	}

}
