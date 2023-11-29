package gui

import (
	"github.com/AllenDang/giu"
	"github.com/go-gl/gl/v4.1-core/gl"
	"image/color"
	"math"
)

const (
	CIRCLE_VERTS     = 8 * 4
	TEMP_COORD_COUNT = 100
)

// 负责渲染
type imguiRender struct {
	circleVerts []float64
	tempCoords  []float64
	tempNormals []float64
	ftex        uint32
	ui *Gui
}

func newImguiRender(ui *Gui) *imguiRender {
	render := &imguiRender{
		circleVerts: make([]float64, CIRCLE_VERTS*2),
		tempCoords:  make([]float64, TEMP_COORD_COUNT*2),
		tempNormals: make([]float64, TEMP_COORD_COUNT*2),
		ui:          ui,
	}
	render.init("")
	return render
}
func (render *imguiRender) init(fontPath string) {
	if fontPath == "" {
		return
	}
	for i := 0; i < CIRCLE_VERTS; i++ {
		a := float64(i) / CIRCLE_VERTS * math.Pi * 2
		render.circleVerts[i*2+0] = math.Cos(a)
		render.circleVerts[i*2+1] = math.Sin(a)
	}
	return
}

func (render *imguiRender) render(cmds ...*imguiGfxCmd) {
	s := 1.0 / 8.0
	gl.Disable(gl.SCISSOR_TEST)
	for _, v := range cmds {
		render.doRender(v, s)
	}
	gl.Enable(gl.SCISSOR_TEST)
}

func (render *imguiRender) doRender(cmd *imguiGfxCmd, s float64) {

	switch cmd.Type {
	case IMGUI_GFXCMD_RECT:
		if cmd.rect.r == 0 {
			render.drawRect(float64(cmd.rect.x)*s+0.5, float64(cmd.rect.y)*s+0.5,
				float64(cmd.rect.w)*s-1, float64(cmd.rect.h)*s-1,
				1.0, cmd.col)
		} else {
			render.drawRoundedRect(float64(cmd.rect.x)*s+0.5, float64(cmd.rect.y)*s+0.5,
				float64(cmd.rect.w)*s-1, float64(cmd.rect.h)*s-1,
				float64(cmd.rect.r)*s, 1.0, cmd.col)
		}
	case IMGUI_GFXCMD_LINE:
		render.drawLine(float64(cmd.line.x0)*s, float64(cmd.line.y0)*s, float64(cmd.line.x1)*s, float64(cmd.line.y1)*s, float64(cmd.line.r)*s, 1.0, cmd.col)
	case IMGUI_GFXCMD_TRIANGLE:
		if cmd.flags == 1 {
			verts := []float64{
				float64(cmd.rect.x)*s + 0.5, float64(cmd.rect.y)*s + 0.5,
				float64(cmd.rect.x)*s + 0.5 + float64(cmd.rect.w)*s - 1, float64(cmd.rect.y)*s + 0.5 + float64(cmd.rect.h)*s/2 - 0.5,
				float64(cmd.rect.x)*s + 0.5, float64(cmd.rect.y)*s + 0.5 + float64(cmd.rect.h)*s - 1,
			}
			render.drawPolygon(verts, 3, 1.0, cmd.col)
		}
		if cmd.flags == 2 {
			verts := []float64{
				float64(cmd.rect.x)*s + 0.5, float64(cmd.rect.y)*s + 0.5 + float64(cmd.rect.h)*s - 1,
				float64(cmd.rect.x)*s + 0.5 + float64(cmd.rect.w)*s/2 - 0.5, float64(cmd.rect.y)*s + 0.5,
				float64(cmd.rect.x)*s + 0.5 + float64(cmd.rect.w)*s - 1, float64(cmd.rect.y)*s + 0.5 + float64(cmd.rect.h)*s - 1,
			}
			render.drawPolygon(verts, 3, 1.0, cmd.col)
		}
	case IMGUI_GFXCMD_TEXT:
		render.drawText(float64(cmd.text.x), float64(cmd.text.y), cmd.text.text, cmd.text.align, cmd.col)
	case IMGUI_GFXCMD_SCISSOR:
		if cmd.flags > 0 {
			gl.Enable(gl.SCISSOR_TEST)
			gl.Scissor(int32(cmd.rect.x), int32(cmd.rect.y), int32(cmd.rect.w), int32(cmd.rect.h))
		} else {
			gl.Disable(gl.SCISSOR_TEST)
		}
	}
}

func (render *imguiRender) drawPolygon(coords []float64, numCoords int, r float64, col color.RGBA) {
	if numCoords > TEMP_COORD_COUNT {
		numCoords = TEMP_COORD_COUNT
	}
	i := 0
	j := numCoords - 1
	for i < numCoords {
		v0 := coords[j*2:]
		v1 := coords[i*2:]
		dx := v1[0] - v0[0]
		dy := v1[1] - v0[1]
		d := math.Sqrt(dx*dx + dy*dy)
		if d > 0 {
			d = 1.0 / d
			dx *= d
			dy *= d
		}
		render.tempNormals[j*2+0] = dy
		render.tempNormals[j*2+1] = -dx
		j = i
		i++
	}
	i = 0
	j = numCoords - 1
	for i < numCoords {
		dlx0 := render.tempNormals[j*2+0]
		dly0 := render.tempNormals[j*2+1]
		dlx1 := render.tempNormals[i*2+0]
		dly1 := render.tempNormals[i*2+1]
		dmx := (dlx0 + dlx1) * 0.5
		dmy := (dly0 + dly1) * 0.5
		dmr2 := dmx*dmx + dmy*dmy
		if dmr2 > 0.000001 {
			scale := 1.0 / dmr2
			if scale > 10.0 {
				scale = 10.0
			}
			dmx *= scale
			dmy *= scale
		}
		render.tempCoords[i*2+0] = coords[i*2+0] + dmx*r
		render.tempCoords[i*2+1] = coords[i*2+1] + dmy*r
		j = i
		i++
	}
	glBegin(GL_TRIANGLES);
	glColor4ubv((GLubyte*)&col);
	i = 0
	j = numCoords - 1
	for i < numCoords {
		glVertex2fv(&coords[i*2]);
		glVertex2fv(&coords[j*2]);
		gl.Color4ubv((GLubyte*)&colTrans);
		glVertex2fv(&g_tempCoords[j*2]);

		glVertex2fv(&g_tempCoords[j*2]);
		glVertex2fv(&g_tempCoords[i*2]);

		glColor4ubv((GLubyte*)&col);
		glVertex2fv(&coords[i*2]);
		//render.pen.AddTriangle(imguiPoint(coords[i*2:], render.ui.height), imguiPoint(coords[j*2:], render.ui.height), imguiPoint(render.tempCoords[j*2:], render.ui.height), col, float32(r))
		//render.pen.AddTriangle(imguiPoint(render.tempCoords[j*2:], render.ui.height), imguiPoint(render.tempCoords[i*2:], render.ui.height), imguiPoint(coords[i*2:], render.ui.height), col, float32(r))
		j = i
		i++
	}
	for i = 2; i < numCoords; i++ {
		glVertex2fv(&coords[0]);
		glVertex2fv(&coords[(i-1)*2]);
		glVertex2fv(&coords[i*2]);
		//render.pen.AddTriangle(imguiPoint(coords, render.ui.height), imguiPoint(coords[(i-1)*2:], render.ui.height), imguiPoint(coords[i*2:], render.ui.height), col, float32(r))
	}
}

func (render *imguiRender) drawRect(x, y, w float64, h float64, fth float64, col color.RGBA) {
	verts :=
		[]float64{
			x + 0.5, y + 0.5,
			x + w - 0.5, y + 0.5,
			x + w - 0.5, y + h - 0.5,
			x + 0.5, y + h - 0.5,
		}
	render.drawPolygon(verts, 4, fth, col)
}

func (render *imguiRender) drawRoundedRect(x, y, w, h, r, fth float64, col color.RGBA) {
	n := CIRCLE_VERTS / 4
	verts := make([]float64, (n+1)*4*2)
	cverts := render.circleVerts
	v := 0

	for i := 0; i <= n; i++ {
		verts[v] = x + w - r + cverts[i*2]*r
		v++
		verts[v] = y + h - r + cverts[i*2+1]*r
		v++
	}

	for i := n; i <= n*2; i++ {
		verts[v] = x + r + cverts[i*2]*r
		v++
		verts[v] = y + h - r + cverts[i*2+1]*r
		v++
	}

	for i := n * 2; i <= n*3; i++ {
		verts[v] = x + r + cverts[i*2]*r
		v++
		verts[v] = y + r + cverts[i*2+1]*r
		v++
	}

	for i := n * 3; i < n*4; i++ {
		verts[v] = x + w - r + cverts[i*2]*r
		v++
		verts[v] = y + r + cverts[i*2+1]*r
		v++
	}
	verts[v] = x + w - r + cverts[0]*r
	v++
	verts[v] = y + r + cverts[1]*r
	v++

	render.drawPolygon(verts, (n+1)*4, fth, col)
}

func (render *imguiRender) drawLine(x0, y0, x1, y1, r, fth float64, col color.RGBA) {
	dx := x1 - x0
	dy := y1 - y0
	d := math.Sqrt(dx*dx + dy*dy)
	if d > 0.0001 {
		d = 1.0 / d
		dx *= d
		dy *= d
	}
	nx := dy
	ny := -dx
	verts := make([]float64, 4*2)
	r -= fth
	r *= 0.5
	if r < 0.01 {
		r = 0.01
	}
	dx *= r
	dy *= r
	nx *= r
	ny *= r

	verts[0] = x0 - dx - nx
	verts[1] = y0 - dy - ny

	verts[2] = x0 - dx + nx
	verts[3] = y0 - dy + ny

	verts[4] = x1 + dx + nx
	verts[5] = y1 + dy + ny

	verts[6] = x1 + dx - nx
	verts[7] = y1 + dy - ny

	render.drawPolygon(verts, 4, fth, col)
}
func  getBakedQuad(stbtt_bakedchar *chardata,  pw,  ph,  char_index int , xpos, ypos []float64, q *stbtt_aligned_quad ) {
stbtt_bakedchar *b = chardata + char_index;
 round_x := math.Floor(*xpos + b.xoff);
round_y := math.Floor(*ypos - b.yoff);

q.x0 = (float)round_x;
q.y0 = (float)round_y;
q.x1 = (float)round_x + b.x1 - b.x0;
q.y1 = (float)round_y - b.y1 + b.y0;

q.s0 = b.x0 / (float)pw;
q.t0 = b.y0 / (float)pw;
q.s1 = b.x1 / (float)ph;
q.t1 = b.y1 / (float)ph;

*xpos += b.xadvance;
}

var  g_tabStops = [4]float64{150, 210, 270, 330};

func (render *imguiRender) drawText(x, y float64, text string, align giu.AlignmentType, col color.RGBA) {
	if (!g_ftex) {return;}
	if (!text) {return;}

	if (align == IMGUI_ALIGN_CENTER){
		x -= getTextLength(g_cdata, text)/2;
	} else if (align == IMGUI_ALIGN_RIGHT){
		x -= getTextLength(g_cdata, text);
	}
	glColor4ub(col&0xff, (col>>8)&0xff, (col>>16)&0xff, (col>>24)&0xff);
	glEnable(GL_TEXTURE_2D);
	// assume orthographic projection with units = screen pixels, origin at top left
	glBindTexture(GL_TEXTURE_2D, g_ftex);
	glBegin(GL_TRIANGLES);
	 ox := x;

	for  (*text){
		int c = (unsigned char)*text;
		if (c == '\t') {
			for  i := 0; i < 4; i++{
			if (x < g_tabStops[i]+ox) {
			x = g_tabStops[i]+ox;
			break;
			}
			}
		} else if (c >= 32 && c < 128) {
		stbtt_aligned_quad q;
		getBakedQuad(g_cdata, 512,512, c-32, &x,&y,&q);

		glTexCoord2f(q.s0, q.t0);
		glVertex2f(q.x0, q.y0);
		glTexCoord2f(q.s1, q.t1);
		glVertex2f(q.x1, q.y1);
		glTexCoord2f(q.s1, q.t0);
		glVertex2f(q.x1, q.y0);

		glTexCoord2f(q.s0, q.t0);
		glVertex2f(q.x0, q.y0);
		glTexCoord2f(q.s0, q.t1);
		glVertex2f(q.x0, q.y1);
		glTexCoord2f(q.s1, q.t1);
		glVertex2f(q.x1, q.y1);
	}
		text++
	}
	glEnd();
	glDisable(GL_TEXTURE_2D);
}

func (render *imguiRender) Update() {
	render.render(render.ui.gs.getRenderCmd()...)

}
