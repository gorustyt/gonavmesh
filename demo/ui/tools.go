package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
	"github.com/gorustyt/gonavmesh/demo/config"
)

type Tools struct {
	toolsCheckGroup *fyne.Container
	curContext      *fyne.Container
	ctx             *Context

	toolGroup *widget.RadioGroup
}

var (
	toolsMap map[string]ToolsRender
)

func InitToolsMap(ctx *Context) {
	toolsMap = map[string]ToolsRender{
		config.TOOL_NAVMESH_TESTER:      NewToolsTestNavmesh(ctx),
		config.TOOL_NAVMESH_PRUNE:       NewToolsPruneNavmesh(ctx),
		config.TOOL_OFFMESH_CONNECTION:  NewToolOffMeshConnection(ctx),
		config.TOOL_OFFMESH_Links:       NewToolOffMeshConnection(ctx),
		config.TOOL_CONVEX_VOLUME:       NewToolsCreateConvexVolumes(ctx),
		config.TOOL_CROWD:               NewToolsCrowds(ctx),
		config.TOOL_CreateTiles:         NewToolsCreateTitles(ctx),
		config.TOOL_HighlightTileCache:  NewToolsHighlightTileCache(ctx),
		config.TOOL_CreateTempObstacles: NewToolsToolsCreateTempObstacles(ctx),
	}
}

type ToolsRender interface {
	GetRenderObjs() []fyne.CanvasObject
}

func NewTools(ctx *Context) *Tools {
	t := &Tools{curContext: container.NewVBox(), ctx: ctx}
	ctx.AppendSampleChange(t)
	t.toolGroup = widget.NewRadioGroup([]string{}, func(c string) {
		v := toolsMap[c]
		if v == nil {
			return
		}
		ctx.Config().ToolsConfig.Uid = c
		t.changeContext(v.GetRenderObjs()...)
	})

	t.toolsCheckGroup = container.NewVBox(
		widget.NewLabel("Tools"),
		t.toolGroup,
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
	s.SetMinSize(fyne.NewSize(100, 600))
	return container.NewVBox(t.toolsCheckGroup, s)
}

func (t *Tools) setToolsGroup(sample string) {
	switch sample {
	case config.SampleSoloMesh:
		t.toolGroup.Selected = config.TOOL_NAVMESH_TESTER
		t.toolGroup.Options = []string{
			config.TOOL_NAVMESH_TESTER,
			config.TOOL_NAVMESH_PRUNE,

			config.TOOL_OFFMESH_CONNECTION,
			config.TOOL_CONVEX_VOLUME,
			config.TOOL_CROWD,
		}
	case config.SampleTileMesh:
		t.toolGroup.Selected = config.TOOL_CreateTiles
		t.toolGroup.Options = []string{
			config.TOOL_NAVMESH_TESTER,
			config.TOOL_NAVMESH_PRUNE,
			config.TOOL_CreateTiles,
			config.TOOL_OFFMESH_Links,
			config.TOOL_CONVEX_VOLUME,
			config.TOOL_CROWD,
		}
	case config.SampleTempObstacles:
		t.toolGroup.Selected = config.TOOL_CreateTempObstacles
		t.toolGroup.Options = []string{
			config.TOOL_NAVMESH_TESTER,
			config.TOOL_HighlightTileCache,
			config.TOOL_CreateTempObstacles,
			config.TOOL_OFFMESH_Links,
			config.TOOL_CONVEX_VOLUME,
			config.TOOL_CROWD,
		}
	}
	t.toolGroup.Refresh()
	t.changeContext(toolsMap[t.toolGroup.Selected].GetRenderObjs()...)
}
func (t *Tools) SampleChange(sample string) {
	t.setToolsGroup(sample)
}
