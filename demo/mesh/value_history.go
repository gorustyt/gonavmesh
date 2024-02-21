package mesh

import (
	"fmt"
	"image/color"
)

const (
	MAX_HISTORY = 256
)

type GraphParams struct {
	x, y, w, h, pad int
	vmin, vmax      float64
	ndiv            int
	units           [16]int
	gs              *guiState
}

func newGraphParams(gs *guiState) *GraphParams {
	return &GraphParams{gs: gs}
}

func (g *GraphParams) setRect(ix, iy, iw, ih, ipad int) {
	g.x = ix
	g.y = iy
	g.w = iw
	g.h = ih
	g.pad = ipad
}

func (g *GraphParams) setValueRange(ivmin, ivmax float64, indiv int, iunits []int) {
	g.vmin = ivmin
	g.vmax = ivmax
	g.ndiv = indiv
	copy(g.units[:], iunits)
}

type ValueHistory struct {
	m_samples  []float64
	m_hsamples int
}

func newValueHistory() *ValueHistory {
	return &ValueHistory{m_samples: make([]float64, MAX_HISTORY)}
}
func (h *ValueHistory) addSample(val float64) {
	h.m_hsamples = (h.m_hsamples + MAX_HISTORY - 1) % MAX_HISTORY
	h.m_samples[h.m_hsamples] = val
}

func (h *ValueHistory) getSampleCount() int {
	return MAX_HISTORY
}

func (h *ValueHistory) getSample(i int) float64 {
	return h.m_samples[(h.m_hsamples+i)%MAX_HISTORY]
}

func (h *ValueHistory) getSampleMin() float64 {
	val := h.m_samples[0]
	for i := 1; i < MAX_HISTORY; i++ {
		if h.m_samples[i] < val {
			val = h.m_samples[i]
		}
	}

	return val
}

func (h *ValueHistory) getSampleMax() float64 {
	val := h.m_samples[0]
	for i := 1; i < MAX_HISTORY; i++ {
		if h.m_samples[i] > val {
			val = h.m_samples[i]
		}
	}

	return val
}

func (h *ValueHistory) getAverage() float64 {
	val := 0.0
	for i := 0; i < MAX_HISTORY; i++ {
		val += h.m_samples[i]
	}

	return val / MAX_HISTORY
}

func (g *GraphParams) drawGraphBackground() {
	// BG
	g.gs.imguiDrawRoundedRect(float64(g.x), float64(g.y), float64(g.w), float64(g.h), float64(g.pad), imguiRGBA(64, 64, 64, 128))

	sy := (float64(g.h) - float64(g.pad)*2) / float64(g.vmax-g.vmin)
	oy := float64(g.y) + float64(g.pad) - g.vmin*sy
	// Divider Lines
	for i := 0; i <= g.ndiv; i++ {
		u := float64(i) / float64(g.ndiv)
		v := g.vmin + (g.vmax-g.vmin)*u
		text := fmt.Sprintf("%.2f %s", v, g.units)
		fy := oy + v*sy
		g.gs.imguiDrawText(g.x+g.w-g.pad, int(fy-4), IMGUI_ALIGN_RIGHT, text, imguiRGBA(0, 0, 0, 255))
		g.gs.imguiDrawLine(float64(g.x)+float64(g.pad), fy, float64(g.x)+float64(g.w)-float64(g.pad)-50, fy, 1.0, imguiRGBA(0, 0, 0, 64))
	}
}

func (g *GraphParams) drawGraph(graph *ValueHistory,
	idx int, label string, col color.RGBA) {
	sx := float64(g.w-g.pad*2) / float64(graph.getSampleCount())
	sy := float64(g.h-g.pad*2) / float64(g.vmax-g.vmin)
	ox := float64(g.x) + float64(g.pad)
	oy := float64(g.y) + float64(g.pad) - g.vmin*sy

	// Values
	px := 0.0
	py := 0.0
	for i := 0; i < graph.getSampleCount()-1; i++ {
		x := ox + float64(i)*sx
		y := oy + graph.getSample(i)*sy
		if i > 0 {
			g.gs.imguiDrawLine(px, py, x, y, 2.0, col)
		}

		px = x
		py = y
	}

	// Label
	size := 15
	spacing := 10
	ix := g.x + g.w + 5
	iy := g.y + g.h - (idx+1)*(size+spacing)

	g.gs.imguiDrawRoundedRect(float64(ix), float64(iy), float64(size), float64(size), 2.0, col)

	text := fmt.Sprintf("%.2f %v", graph.getAverage(), g.units)
	g.gs.imguiDrawText(ix+size+5, iy+3, IMGUI_ALIGN_LEFT, label, imguiRGBA(255, 255, 255, 192))
	g.gs.imguiDrawText(ix+size+150, iy+3, IMGUI_ALIGN_RIGHT, text, imguiRGBA(255, 255, 255, 128))
}
