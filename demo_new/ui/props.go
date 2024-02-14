package ui

import (
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
)

type Props struct {
	c *fyne.Container
}

func NewProps() *Props {
	p := &Props{}
	p.c = container.NewVBox(
		widget.NewLabel("prop"),
		widget.NewRadioGroup([]string{
			"show logs",
			"show tools",
		}, func(s string) {

		}),
	)
	return p
}
func (p *Props) GetRenderObj() fyne.CanvasObject {
	return container.NewScroll(p.c)
}
