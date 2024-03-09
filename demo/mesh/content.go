package mesh

import (
	"github.com/go-gl/mathgl/mgl32"
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/canvas3d/canvas3d_render"
	"github.com/gorustyt/fyne/v2/canvas3d/context"
	"github.com/gorustyt/fyne/v2/canvas3d/context/enum"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/data/binding"
	"github.com/gorustyt/fyne/v2/driver/desktop"
	"github.com/gorustyt/fyne/v2/widget"
	"github.com/gorustyt/gonavmesh/common"
	"github.com/gorustyt/gonavmesh/debug_utils"
	"github.com/gorustyt/gonavmesh/demo/config"
	"log/slog"
	"math"
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
	w fyne.Window
	widget.BaseWidget
	labelC     fyne.CanvasObject
	label1Data binding.String
	label2Data binding.String
	canvas3D   canvas3d_render.ICanvas3d

	cfg    *config.Config
	sample ISample
	geom   *InputGeom

	state *ContentState
	ctx   context.Context
	dd    debug_utils.DuDebugDraw
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
	c.canvas3D.AppendDefaultRenderFunc(func(ctx context.Painter) {
		c.ctx = ctx.GetContext()
		if c.dd == nil {
			c.dd = NewDebugDrawGL(c.ctx)
		}
		c.sample.setDebugDraw(c.dd)
		c.sample.HandleRender()
	})
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
	if c.sample != nil {
		c.state.matchMeshBounds(geom)
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
	c.state.matchMeshBounds(c.geom)
}

func NewContent(w fyne.Window) *Content {
	cfg := config.NewConfig()
	c := &Content{
		w:          w,
		state:      NewContentState(),
		label1Data: binding.NewString(),
		label2Data: binding.NewString(),
		canvas3D:   canvas3d_render.NewCanvas3d(),
		cfg:        cfg}
	cfg.PropsConfig.OnInputMesh = c.InputMeshChange
	cfg.ToolsConfig.OnUidChange = c.ToolChange
	c.ExtendBaseWidget(c)
	c.setKeyEvent()
	c.labelC = container.NewVBox(
		widget.NewLabelWithData(c.label1Data),
		widget.NewLabelWithData(c.label2Data),
	)
	return c
}

func (c *Content) setKeyEvent() {
	if deskCanvas, ok := c.w.Canvas().(desktop.Canvas); ok {
		deskCanvas.SetOnKeyDown(c.state.OnKeyDown)
		deskCanvas.SetOnKeyUp(c.state.OnKeyUp)
	}
}

func (c *Content) updateMarker() {
	vec := c.state.projectMarkerPosition()
	// Draw marker circle
	c.ctx.ExtLineWidth(5.0)
	c.ctx.ExtColor4ub(240, 220, 0, 196)
	c.ctx.ExtBegin(enum.LineLoop)
	r := 25.0
	for i := float64(0); i < 20; i++ {
		a := i / 20.0 * math.Pi * 2
		fx := float64(vec[0]) + math.Cos(a)*r
		fy := float64(vec[1]) + math.Sin(a)*r
		c.ctx.ExtVertex2f(float32(fx), float32(fy))
	}
	c.ctx.ExtEnd()
	c.ctx.ExtLineWidth(1.0)
}

type ContentState struct {
	camr     float32
	viewPort mgl32.Vec4
	model    mgl32.Mat4
	view     mgl32.Mat4
	project  mgl32.Mat4

	cameraEulers   mgl32.Vec2
	cameraPos      mgl32.Vec3
	rayStart       mgl32.Vec3
	rayEnd         mgl32.Vec3
	mousePos       mgl32.Vec2
	origMousePos   mgl32.Vec2
	markerPosition mgl32.Vec3
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

func (s *ContentState) matchMeshBounds(geom *InputGeom) {
	if geom == nil {
		return
	}
	var bmin []float32
	var bmax []float32
	if geom != nil {
		bmin = geom.getNavMeshBoundsMin()
		bmax = geom.getNavMeshBoundsMax()
	}
	// Reset camera and fog to match the mesh bounds.
	if len(bmin) != 0 && len(bmax) != 0 {
		s.camr = float32(math.Sqrt(float64(common.Sqr(bmax[0]-bmin[0])+
			common.Sqr(bmax[1]-bmin[1])+
			common.Sqr(bmax[2]-bmin[2])))) / 2
		s.cameraPos[0] = (bmax[0]+bmin[0])/2 + s.camr
		s.cameraPos[1] = (bmax[1]+bmin[1])/2 + s.camr
		s.cameraPos[2] = (bmax[2]+bmin[2])/2 + s.camr
		s.camr *= 3
	}
	s.cameraEulers[0] = 45
	s.cameraEulers[1] = -45
}

func (s *ContentState) projectMarkerPosition() mgl32.Vec3 {
	return mgl32.Project(mgl32.Vec3{s.markerPosition[0], s.markerPosition[1], s.markerPosition[2]},
		s.model,
		s.project,
		int(s.viewPort[0]),
		int(s.viewPort[1]),
		int(s.viewPort[2]),
		int(s.viewPort[3]))
}

func (s *ContentState) OnKeyUp(ev *fyne.KeyEvent) {

}

func (s *ContentState) OnKeyDown(ev *fyne.KeyEvent) {

}

func (s *ContentState) update() {
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
