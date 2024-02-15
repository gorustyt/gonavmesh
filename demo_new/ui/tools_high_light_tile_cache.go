package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
	"gonavamesh/demo_new/config"
)

type ToolsHighlightTileCache struct {
	c *fyne.Container
}

func (t ToolsHighlightTileCache) GetRenderObjs() []fyne.CanvasObject {
	return []fyne.CanvasObject{
		t.c,
	}
}

func NewToolsHighlightTileCache(ctx *Context) *ToolsHighlightTileCache {
	c := &ToolsHighlightTileCache{}
	b := widget.NewRadioGroup([]string{
		config.HighLightTitleCacheDrawAreas,
		config.HighLightTitleDrawRegions,
		config.HighLightTitleDrawContours,
		config.HighLightTitleDrawMesh,
	}, func(s string) {
		ctx.Config().ToolsConfig.HighlightDrawType = s
	})
	c.c = container.NewVBox(
		widget.NewLabel("Highlight Tile Cache"),
		widget.NewLabel("Click LMB to highlight a tile."),
		b,
	)
	return c
}
