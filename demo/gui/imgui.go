package gui

import (
	"fmt"
	g "github.com/AllenDang/giu"
	"github.com/AllenDang/imgui-go"
	"image"
	"image/color"
)

func onClickMe() {
	fmt.Println("Hello world!")
}

func onImSoCute() {
	fmt.Println("Im sooooooo cute!!")
}

func loop() {
	h := ""
	c := g.GetCanvas()
	c.AddText(image.Point{X: 11, Y: 12}, color.Black, "helllofdfd")
	g.SingleWindow().Layout(
		g.Label("Hello world from giu"),
		g.Row(
			g.Button("Click Me").OnClick(onClickMe),
			g.Button("I'm so cute").OnClick(onImSoCute),
			g.InputText(&h).Callback(func(data imgui.InputTextCallbackData) int32 {
				if data.EventKey() == int(g.KeyEnter) {
					fmt.Println(h)
				}
				fmt.Println("=====fdfdf")
				return 0
			})),
		g.Menu("hello world").Layout(g.MenuItem("hello")),
		c,
	)

}

func InitImgui() {
	wnd := g.NewMasterWindow("Hello world", 400, 200, g.MasterWindowFlagsNotResizable)
	wnd.Run(loop)
}
