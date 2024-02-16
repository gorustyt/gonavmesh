package ui

import (
	"fmt"
	"github.com/gorustyt/fyne/v2"
	"github.com/gorustyt/fyne/v2/container"
	"github.com/gorustyt/fyne/v2/widget"
)

func PackSlider(slider *widget.Slider, changes ...func(f float64)) fyne.CanvasObject {
	l := widget.NewLabel(fmt.Sprintf("%.1f", slider.Value))
	slider.OnChanged = func(f float64) {
		l.SetText(fmt.Sprintf("%.1f", f))
		for _, v := range changes {
			v(f)
		}
	}
	return container.NewBorder(nil, nil, nil, l, slider)
}
