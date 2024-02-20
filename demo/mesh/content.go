package mesh

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/gonavmesh/demo/config"
)

type Content struct {
	cfg    *config.Config
	sample *Sample
}

func (c *Content) GetConfig() *config.Config {
	return c.cfg
}

func (c *Content) MinSize() fyne.Size {
	return fyne.Size{
		Width:  600,
		Height: 1080,
	}
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

}

func NewContent() *Content {
	cfg := config.NewConfig()
	c := &Content{
		cfg: cfg}
	c.sample = NewSample(c)
	c.cfg.PropsConfig.OnLoadClick = c.sample.Load
	return c
}
