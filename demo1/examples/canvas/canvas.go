package main

import (
	"image"
	"image/color"

	g "github.com/AllenDang/giu"
)

var texture *g.Texture

func loop() {
	g.SingleWindow().Layout(
		g.Label("Canvas demo"),
		g.Custom(func() {
			canvas := g.GetCanvas()
			pos := g.GetCursorScreenPos()
			col := color.RGBA{200, 75, 75, 255}
			canvas.AddLine(pos, pos.Add(image.Pt(100, 100)), col, 1)
			canvas.AddRect(pos.Add(image.Pt(110, 0)), pos.Add(image.Pt(200, 100)), col, 5, g.DrawFlagsRoundCornersAll, 1)
			canvas.AddRectFilled(pos.Add(image.Pt(220, 0)), pos.Add(image.Pt(320, 100)), col, 0, 0)

			pos0 := pos.Add(image.Pt(0, 110))
			cp0 := pos.Add(image.Pt(80, 110))
			cp1 := pos.Add(image.Pt(50, 210))
			pos1 := pos.Add(image.Pt(120, 210))
			canvas.AddBezierCubic(pos0, cp0, cp1, pos1, col, 1, 0)

			p1 := pos.Add(image.Pt(160, 110))
			p2 := pos.Add(image.Pt(120, 210))
			p3 := pos.Add(image.Pt(210, 210))
			p4 := pos.Add(image.Pt(210, 150))
			// canvas.AddTriangle(p1, p2, p3, col, 2)
			canvas.AddQuad(p1, p2, p3, p4, col, 1)

			p1 = p1.Add(image.Pt(120, 60))
			canvas.AddCircleFilled(p1, 50, col)

			p1 = pos.Add(image.Pt(10, 400))
			p2 = pos.Add(image.Pt(50, 440))
			p3 = pos.Add(image.Pt(200, 500))
			canvas.PathLineTo(p1)
			canvas.PathLineTo(p2)
			canvas.PathBezierCubicCurveTo(p2.Add(image.Pt(40, 0)), p3.Add(image.Pt(-50, 0)), p3, 0)
			canvas.PathStroke(col, 0, 1)

			if texture != nil {
				canvas.AddImage(texture, image.Pt(pos.X+350, pos.Y+25), image.Pt(pos.X+500, pos.Y+125))
				p1 = image.Pt(pos.X+350, pos.Y+200)
				p2 = image.Pt(pos.X+350, pos.Y+400)
				p3 = image.Pt(pos.X+550, pos.Y+300)
				p4 = image.Pt(pos.X+500, pos.Y+200)
				canvas.AddImageV(texture, p1, p2, p3, p4,
					image.Pt(0, 0),
					image.Pt(1, 0),
					image.Pt(1, 1),
					image.Pt(0, 1),
					color.RGBA{255, 255, 255, 100},
				)
			}
		}),
	)
}

func main() {
	wnd := g.NewMasterWindow("Canvas", 600, 600, g.MasterWindowFlagsNotResizable)

	img, _ := g.LoadImage("gopher.png")
	g.EnqueueNewTextureFromRgba(img, func(tex *g.Texture) {
		texture = tex
	})

	wnd.Run(loop)
}
