package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"gonavamesh/demo_new/config"
	"gonavamesh/demo_new/mesh"
)

func GetMenu() fyne.CanvasObject {
	c := mesh.NewContent()
	ctx := NewContext(c.GetConfig())
	InitToolsMap(ctx)
	root := container.NewBorder(nil, nil, NewTools(ctx).GetRenderObj(), NewProps(ctx).GetRenderObj(), c)
	ctx.AfterInit()
	return root
}

type SampleChange interface {
	SampleChange(Sample string)
}

type Context struct {
	root fyne.CanvasObject

	Sample  string
	Changes []fyne.CanvasObject
	Show    map[string][]fyne.CanvasObject
	cfg     *config.Config

	SampleChanges []SampleChange
	afterInit     []func()
}

func NewContext(c *config.Config) *Context {
	return &Context{
		Show: map[string][]fyne.CanvasObject{},
		cfg:  c,
	}
}
func (s *Context) AppendAfterInit(ss ...func()) {
	s.afterInit = append(s.afterInit, ss...)
}
func (s *Context) AppendSampleChange(ss ...SampleChange) {
	s.SampleChanges = append(s.SampleChanges, ss...)
}

func (s *Context) OnSampleChange(sample string) {
	s.Sample = sample
	for _, v := range s.SampleChanges {
		v.SampleChange(sample)
	}
	for _, v := range s.Changes {
		v.Hide()
	}
	for _, v := range s.Show[sample] {
		v.Show()
	}
}

func (s *Context) Config() *config.Config {
	return s.cfg
}

func (s *Context) Refresh() {
	s.root.Refresh()
}

func (s *Context) AfterInit() {
	for _, v := range s.afterInit {
		v()
	}
}
