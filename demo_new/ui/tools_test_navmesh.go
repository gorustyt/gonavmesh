package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type ToolsTestNavmesh struct {
	c *fyne.Container
}

func NewToolsTestNavmesh() *ToolsTestNavmesh {
	t := &ToolsTestNavmesh{}
	checkGroup1 := widget.NewCheckGroup([]string{
		config.DESC_SAMPLE_POLYFLAGS_WALK,
		config.DESC_SAMPLE_POLYFLAGS_SWIM,
		config.DESC_SAMPLE_POLYFLAGS_DOOR,
		config.DESC_SAMPLE_POLYFLAGS_JUMP,
	}, func(strings []string) {

	})
	checkGroup2 := widget.NewCheckGroup([]string{
		config.DESC_SAMPLE_POLYFLAGS_WALK,
		config.DESC_SAMPLE_POLYFLAGS_SWIM,
		config.DESC_SAMPLE_POLYFLAGS_DOOR,
		config.DESC_SAMPLE_POLYFLAGS_JUMP,
	}, func(strings []string) {

	})
	b1 := widget.NewButton("Set Random Start", func() {

	})
	b2 := widget.NewButton("Set Random End", func() {

	})
	b3 := widget.NewButton("Make Random Points", func() {

	})
	b4 := widget.NewButton("Make Random Points Around", func() {

	})

	group1 := widget.NewRadioGroup([]string{
		config.DT_STRAIGHTPATH_NONE_CROSSINGS,
		config.DT_STRAIGHTPATH_AREA_CROSSINGS,
		config.DT_STRAIGHTPATH_ALL_CROSSINGS,
	}, func(v string) {

	})
	c1 := container.NewVBox(widget.NewLabel("Vertices at crossings"), group1, widget.NewSeparator())
	c1.Hide()
	group1.Selected = config.DT_STRAIGHTPATH_NONE_CROSSINGS
	group := widget.NewRadioGroup([]string{
		config.Desc_TOOLMODE_PATHFIND_FOLLOW,
		config.Desc_OOLMODE_PATHFIND_STRAIGHT,
		config.Desc_TOOLMODE_PATHFIND_SLICED,
		config.Desc_TOOLMODE_RAYCAST,
		config.Desc_TOOLMODE_DISTANCE_TO_WALL,
		config.Desc_TOOLMODE_FIND_POLYS_IN_CIRCLE,
		config.Desc_TOOLMODE_FIND_POLYS_IN_SHAPE,
		config.Desc_TOOLMODE_FIND_LOCAL_NEIGHBOURHOOD,
	}, func(s string) {
		if s == config.Desc_OOLMODE_PATHFIND_STRAIGHT {
			c1.Show()
		} else {
			c1.Hide()
		}
	})

	vbox := container.NewVBox(
		group,
		widget.NewSeparator(),
		c1,
		b1, b2, b3, b4,
		widget.NewSeparator(),
		widget.NewLabel("Include Flags"),
		checkGroup1,
		widget.NewLabel("Exclude Flags"),
		checkGroup2,
	)
	t.c = vbox
	return t
}
func (t *ToolsTestNavmesh) GetRenderObjs() (res []fyne.CanvasObject) {
	return append(res, t.c)
}
