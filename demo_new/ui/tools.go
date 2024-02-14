package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type Tools struct {
	toolsCheckGroup *fyne.Container
	curContext      *fyne.Container
	cfg             *config.Config
}

var (
	toolsMap map[string]ToolsRender
)

func InitToolsMap(cfg *config.Config) {
	toolsMap = map[string]ToolsRender{
		config.TOOL_NAVMESH_TESTER:     NewToolsTestNavmesh(cfg),
		config.TOOL_NAVMESH_PRUNE:      NewToolsPruneNavmesh(cfg),
		config.TOOL_OFFMESH_CONNECTION: NewToolOffMeshConnection(cfg),
		config.TOOL_CONVEX_VOLUME:      NewToolsCreateConvexVolumes(cfg),
		config.TOOL_CROWD:              NewToolsCrowds(cfg),
		config.TOOL_CreateTiles:        NewToolsCreateTitles(cfg),
	}
}

type ToolsRender interface {
	GetRenderObjs() []fyne.CanvasObject
}

func NewTools(cfg *config.Config) *Tools {
	t := &Tools{curContext: container.NewVBox(), cfg: cfg}
	group := widget.NewRadioGroup([]string{
		config.TOOL_NAVMESH_TESTER,
		config.TOOL_NAVMESH_PRUNE,
		config.TOOL_CreateTiles,
		config.TOOL_OFFMESH_CONNECTION,
		config.TOOL_CONVEX_VOLUME,
		config.TOOL_CROWD}, func(c string) {
		v := toolsMap[c]
		if v == nil {
			return
		}
		cfg.ToolsConfig.Uid = c
		t.changeContext(v.GetRenderObjs()...)
	})
	group.Selected = config.TOOL_NAVMESH_TESTER
	t.changeContext(toolsMap[group.Selected].GetRenderObjs()...)
	t.toolsCheckGroup = container.NewVBox(
		widget.NewLabel("Tools"),
		group,
		widget.NewSeparator(),
	)

	return t

}

func (t *Tools) changeContext(cs ...fyne.CanvasObject) {
	t.curContext.RemoveAll()
	for _, v := range cs {
		t.curContext.Add(v)
	}
	if len(cs) > 0 {
		t.curContext.Refresh()
	}
}

func (t *Tools) GetRenderObj() fyne.CanvasObject {
	s := container.NewVScroll(t.curContext)
	s.SetMinSize(fyne.NewSize(200, 600))
	return container.NewVBox(t.toolsCheckGroup, s)
}
