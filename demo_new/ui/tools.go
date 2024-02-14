package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
)

const (
	TOOL_NAVMESH_TESTER     = "Test Navmesh"
	TOOL_NAVMESH_PRUNE      = "Prune Navmesh"
	TOOL_OFFMESH_CONNECTION = "Create Off-Mesh Connections"
	TOOL_CONVEX_VOLUME      = "Create Convex Volumes"
	TOOL_CROWD              = "Create Crowds"
)

type Tools struct {
	toolsCheckGroup *fyne.Container
	curContext      *fyne.Container
}

var (
	toolsMap map[string]ToolsRender
)

func InitToolsMap() {
	toolsMap = map[string]ToolsRender{
		TOOL_NAVMESH_TESTER:     NewToolsTestNavmesh(),
		TOOL_NAVMESH_PRUNE:      NewToolsPruneNavmesh(),
		TOOL_OFFMESH_CONNECTION: NewToolOffMeshConnection(),
		TOOL_CONVEX_VOLUME:      NewToolsCreateConvexVolumes(),
		TOOL_CROWD:              NewToolsCrowds(),
	}
}

type ToolsRender interface {
	GetRenderObjs() []fyne.CanvasObject
}

func NewTools() *Tools {
	t := &Tools{}
	group := widget.NewRadioGroup([]string{
		TOOL_NAVMESH_TESTER,
		TOOL_NAVMESH_PRUNE,
		TOOL_OFFMESH_CONNECTION,
		TOOL_CONVEX_VOLUME,
		TOOL_CROWD}, func(c string) {
		v := toolsMap[c]
		t.changeContext(v.GetRenderObjs()...)
	})
	t.toolsCheckGroup = container.NewVBox(
		widget.NewLabel("Tools"),
		group,
		widget.NewSeparator(),
	)
	t.toolsCheckGroup.Resize(fyne.Size{Width: 300, Height: 1080})
	t.curContext = container.NewVBox()

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
	return container.NewVBox(t.toolsCheckGroup, t.curContext)
}
