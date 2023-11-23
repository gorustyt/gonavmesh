package gui

import (
	"container/list"
	"fmt"
	"github.com/AllenDang/giu"
	"image"
	"image/color"
	"math"
)

// Pull render interface.
type imguiGfxCmdType int

const (
	IMGUI_GFXCMD_RECT = iota
	IMGUI_GFXCMD_TRIANGLE
	IMGUI_GFXCMD_LINE
	IMGUI_GFXCMD_TEXT
	IMGUI_GFXCMD_SCISSOR
)

func imguiPoint(points []float64, height int) image.Point {
	return image.Point{
		X: int(points[0]),
		Y: int(float64(height) - points[1]),
	}
}

func imguiRGBA(r, g, b, a uint8) color.RGBA {
	return color.RGBA{
		R: r,
		G: g,
		B: b,
		A: a,
	}
}

type imguiGfxRect struct {
	x, y, w, h, r int
}

type imguiGfxText struct {
	x, y  int
	align giu.AlignmentType
	text  string
}

type imguiGfxLine struct {
	x0, y0, x1, y1, r int
}

type imguiGfxCmd struct {
	Type  imguiGfxCmdType
	flags int
	pad   [2]int
	col   color.RGBA
	line  *imguiGfxLine
	rect  *imguiGfxRect
	text  *imguiGfxText
}

func newImguiGfxCmd() *imguiGfxCmd {
	return &imguiGfxCmd{
		line: &imguiGfxLine{},
		rect: &imguiGfxRect{},
		text: &imguiGfxText{},
	}
}

const GFXCMD_QUEUE_SIZE = 5000

func (gs *guiState) addGfxCmdScissor(x, y, w, h int) {
	if gs.queue.Len() >= GFXCMD_QUEUE_SIZE {
		return
	}
	cmd := newImguiGfxCmd()
	cmd.Type = IMGUI_GFXCMD_SCISSOR
	cmd.flags = 1 // on/off flag.
	if x < 0 {
		cmd.flags = 0
	}
	cmd.col = color.RGBA{}
	cmd.rect.x = x
	cmd.rect.y = y
	cmd.rect.w = w
	cmd.rect.h = h
	gs.queue.PushBack(cmd)
}

func (gs *guiState) addGfxCmdRect(x, y, w, h float64, color color.RGBA) {
	if gs.queue.Len() >= GFXCMD_QUEUE_SIZE {
		return
	}
	cmd := newImguiGfxCmd()
	cmd.Type = IMGUI_GFXCMD_RECT
	cmd.flags = 0
	cmd.col = color
	cmd.rect.x = int(x * 8.0)
	cmd.rect.y = int(y * 8.0)
	cmd.rect.w = int(w * 8.0)
	cmd.rect.h = int(h * 8.0)
	cmd.rect.r = 0
	gs.queue.PushBack(cmd)
}

func (gs *guiState) addGfxCmdLine(x0, y0, x1, y1, r float64, color color.RGBA) {
	if gs.queue.Len() >= GFXCMD_QUEUE_SIZE {
		return
	}
	cmd := newImguiGfxCmd()
	cmd.Type = IMGUI_GFXCMD_LINE
	cmd.flags = 0
	cmd.col = color
	cmd.line.x0 = int(x0 * 8.0)
	cmd.line.y0 = int(y0 * 8.0)
	cmd.line.x1 = int(x1 * 8.0)
	cmd.line.y1 = int(y1 * 8.0)
	cmd.line.r = int(r * 8.0)
	gs.queue.PushBack(cmd)
}

func (gs *guiState) addGfxCmdRoundedRect(x, y, w, h, r float64, color color.RGBA) {
	if gs.queue.Len() >= GFXCMD_QUEUE_SIZE {
		return
	}
	cmd := newImguiGfxCmd()
	cmd.Type = IMGUI_GFXCMD_RECT
	cmd.flags = 0
	cmd.col = color
	cmd.rect.x = int(x * 8.0)
	cmd.rect.y = int(y * 8.0)
	cmd.rect.w = int(w * 8.0)
	cmd.rect.h = int(h * 8.0)
	cmd.rect.r = int(r * 8.0)
	gs.queue.PushBack(cmd)
}

func (gs *guiState) addGfxCmdTriangle(x, y, w, h, flags int, color color.RGBA) {
	if gs.queue.Len() >= GFXCMD_QUEUE_SIZE {
		return
	}
	cmd := newImguiGfxCmd()
	cmd.Type = IMGUI_GFXCMD_TRIANGLE
	cmd.flags = flags
	cmd.col = color
	cmd.rect.x = x * 8.0
	cmd.rect.y = y * 8.0
	cmd.rect.w = w * 8.0
	cmd.rect.h = h * 8.0
	gs.queue.PushBack(cmd)
}

func (gs *guiState) addGfxCmdText(x, y int, align giu.AlignmentType, text string, color color.RGBA) {
	if gs.queue.Len() >= GFXCMD_QUEUE_SIZE {
		return
	}
	cmd := newImguiGfxCmd()
	cmd.Type = IMGUI_GFXCMD_TEXT
	cmd.flags = 0
	cmd.col = color
	cmd.text.x = x
	cmd.text.y = y
	cmd.text.align = align
	cmd.text.text = text
	gs.queue.PushBack(cmd)
}

type guiState struct {
	left                      bool
	leftPressed, leftReleased bool
	mx, my                    int
	scroll                    int
	active                    int
	hot                       int
	hotToBe                   int
	wentActive                bool
	dragX, dragY              int
	dragOrig                  float64
	widgetX, widgetY, widgetW int
	insideCurrentScroll       bool
	areaId                    int
	widgetId                  int
	isHot                     bool
	isActive                  bool
	queue                     list.List
}

func newGuiState() *guiState {
	return &guiState{
		widgetW: 100,
		mx:      -1,
		my:      -1,
	}
}

func (gs *guiState) anyActive() bool {
	return gs.active != 0
}

func (gs *guiState) IsActive(id int) bool {
	return gs.active == id
}

func (gs *guiState) IsHot(id int) bool {
	return gs.hot == id
}

func (gs *guiState) inRect(x, y, w, h int, checkScrolls ...bool) bool {
	checkScroll := true
	if len(checkScrolls) > 0 {
		checkScroll = checkScrolls[0]
	}
	return (!checkScroll || gs.insideCurrentScroll) && gs.mx >= x && gs.mx <= x+w && gs.my >= y && gs.my <= y+h
}

func (gs *guiState) clearInput() {
	gs.leftPressed = false
	gs.leftReleased = false
	gs.scroll = 0
}

func (gs *guiState) clearActive() {
	gs.active = 0
	// mark all UI for this frame as processed
	gs.clearInput()
}

func (gs *guiState) setActive(id int) {
	gs.active = id
	gs.wentActive = true
}

func (gs *guiState) setHot(id int) {
	gs.hotToBe = id
}

func (gs *guiState) buttonLogic(id int, over bool) bool {
	res := false
	// process down
	if !gs.anyActive() {
		if over {
			gs.setHot(id)
		}

		if gs.IsHot(id) && gs.leftPressed {
			gs.setActive(id)
		}

	}

	// if button is active, then react on left up
	if gs.IsActive(id) {
		gs.isActive = true
		if over {
			gs.setHot(id)
		}

		if gs.leftReleased {
			if gs.IsHot(id) {
				res = true
			}

			gs.clearActive()
		}
	}

	if gs.IsHot(id) {
		gs.isHot = true
	}

	return res
}

func (gs *guiState) updateInput(mx, my int, left bool, scroll int) {
	gs.mx = mx
	gs.my = my
	gs.leftPressed = !gs.left && left
	gs.leftReleased = gs.left && !left
	gs.left = left
	gs.scroll = scroll
}

func (gs *guiState) beginFrame(mx, my int, left bool, scroll int) {
	//开始渲染
	gs.updateInput(gs.mx, my, left, scroll)

	gs.hot = gs.hotToBe
	gs.hotToBe = 0

	gs.wentActive = false
	gs.isActive = false
	gs.isHot = false

	gs.widgetX = 0
	gs.widgetY = 0
	gs.widgetW = 0

	gs.areaId = 1
	gs.widgetId = 1
	gs.queue = list.List{}
}
func (gs *guiState) getRenderCmd() (res []*imguiGfxCmd) {
	e := gs.queue.Front()
	for e != nil {
		res = append(res, e.Value.(*imguiGfxCmd))
		next := e.Next()
		gs.queue.Remove(e)
		e = next
	}
	return res
}
func (gs *guiState) endFrame() {
	gs.clearInput()
}

// ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
const (
	BUTTON_HEIGHT       = 20
	SLIDER_HEIGHT       = 20
	SLIDER_MARKER_WIDTH = 10
	CHECK_SIZE          = 8
	DEFAULT_SPACING     = 4
	TEXT_HEIGHT         = 8
	SCROLL_AREA_PADDING = 6
	INDENT_SIZE         = 16
	AREA_HEADER         = 28
)

var (
	g_scrollTop        = 0
	g_scrollBottom     = 0
	g_scrollRight      = 0
	g_scrollAreaTop    = 0
	g_scrollVal        *int
	g_focusTop         = 0
	g_focusBottom      = 0
	g_scrollId         = 0
	g_insideScrollArea = false
)

func (gs *guiState) imguiIndent() {
	gs.widgetX += INDENT_SIZE
	gs.widgetW -= INDENT_SIZE
}

func (gs *guiState) imguiUnindent() {
	gs.widgetX -= INDENT_SIZE
	gs.widgetW += INDENT_SIZE
}

func (gs *guiState) imguiSeparator() {
	gs.widgetY -= DEFAULT_SPACING * 3
}

func (gs *guiState) imguiSeparatorLine() {
	x := gs.widgetX
	y := gs.widgetY - DEFAULT_SPACING*2
	w := gs.widgetW
	h := 1
	gs.widgetY -= DEFAULT_SPACING * 4
	gs.addGfxCmdRect(float64(x), float64(y), float64(w), float64(h), color.RGBA{R: 255, G: 255, B: 255, A: 32})
}

func (gs *guiState) imguiDrawText(x, y int, align giu.AlignmentType, text string, color color.RGBA) {
	gs.addGfxCmdText(x, y, align, text, color)
}

func (gs *guiState) imguiDrawLine(x0, y0, x1, y1, r float64, color color.RGBA) {
	gs.addGfxCmdLine(x0, y0, x1, y1, r, color)
}

func (gs *guiState) imguiDrawRect(x, y, w, h float64, color color.RGBA) {
	gs.addGfxCmdRect(x, y, w, h, color)
}

func (gs *guiState) imguiDrawRoundedRect(x, y, w, h, r float64, color color.RGBA) {
	gs.addGfxCmdRoundedRect(x, y, w, h, r, color)
}

func (gs *guiState) imguiLabel(text string) {
	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	gs.widgetY -= BUTTON_HEIGHT
	gs.addGfxCmdText(x, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(255, 255, 255, 255))
}

func (gs *guiState) imguiValue(text string) {
	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	w := gs.widgetW
	gs.widgetY -= BUTTON_HEIGHT
	gs.addGfxCmdText(x+w-BUTTON_HEIGHT/2, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignRight, text, imguiRGBA(255, 255, 255, 200))
}
func (gs *guiState) imguiSlider(text string, val *float64, vmin, vmax, vinc float64, enabled bool) bool {
	gs.widgetId++
	id := (gs.areaId << 16) | gs.widgetId

	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	w := gs.widgetW
	h := SLIDER_HEIGHT
	gs.widgetY -= SLIDER_HEIGHT + DEFAULT_SPACING

	gs.addGfxCmdRoundedRect(float64(x), float64(y), float64(w), float64(h), 4.0, imguiRGBA(0, 0, 0, 128))

	rangef := w - SLIDER_MARKER_WIDTH

	u := (*val - vmin) / (vmax - vmin)
	if u < 0 {
		u = 0
	}
	if u > 1 {
		u = 1
	}
	m := int(u * float64(rangef))

	over := enabled && gs.inRect(x+m, y, SLIDER_MARKER_WIDTH, SLIDER_HEIGHT)
	res := gs.buttonLogic(id, over)
	valChanged := false

	if gs.IsActive(id) {
		if gs.wentActive {
			gs.dragX = gs.mx
			gs.dragOrig = u
		}
		if gs.dragX != gs.mx {
			u = gs.dragOrig + float64(gs.mx-gs.dragX)/float64(rangef)
			if u < 0 {
				u = 0
			}
			if u > 1 {
				u = 1
			}
			*val = vmin + u*(vmax-vmin)
			*val = math.Floor(*val/vinc+0.5) * vinc // Snap to vinc
			m = (int)(u * float64(rangef))
			valChanged = true
		}
	}

	if gs.IsActive(id) {
		gs.addGfxCmdRoundedRect(float64(x+m), float64(y), SLIDER_MARKER_WIDTH, SLIDER_HEIGHT, 4.0, imguiRGBA(255, 255, 255, 255))
	} else {
		c := imguiRGBA(255, 255, 255, 64)
		if gs.IsHot(id) {
			c = imguiRGBA(255, 196, 0, 128)
		}
		gs.addGfxCmdRoundedRect(float64(x+m), float64(y), SLIDER_MARKER_WIDTH, SLIDER_HEIGHT, 4.0, c)
	}

	// TODO: fix this, take a look at 'nicenum'.
	digits := (int)(math.Ceil(math.Log10(vinc)))
	if digits >= 0 {
		digits = 0
	} else {
		digits = -digits
	}
	format := fmt.Sprintf("%%.%df", digits)
	msg := fmt.Sprintf(format, *val)
	if enabled {
		c := imguiRGBA(255, 255, 255, 200)
		if gs.IsHot(id) {
			c = imguiRGBA(255, 196, 0, 255)
		}
		gs.addGfxCmdText(x+SLIDER_HEIGHT/2, y+SLIDER_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, c)
		c = imguiRGBA(255, 255, 255, 200)
		if gs.IsHot(id) {
			c = imguiRGBA(255, 196, 0, 255)
		}
		gs.addGfxCmdText(x+w-SLIDER_HEIGHT/2, y+SLIDER_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignRight, msg, c)
	} else {
		gs.addGfxCmdText(x+SLIDER_HEIGHT/2, y+SLIDER_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(128, 128, 128, 200))
		gs.addGfxCmdText(x+w-SLIDER_HEIGHT/2, y+SLIDER_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignRight, msg, imguiRGBA(128, 128, 128, 200))
	}

	return res || valChanged
}

func (gs *guiState) imguiBeginScrollArea(name string, x, y, w, h int, scroll *int) bool {
	gs.areaId++
	gs.widgetId = 0
	g_scrollId = (gs.areaId << 16) | gs.widgetId

	gs.widgetX = x + SCROLL_AREA_PADDING
	gs.widgetY = y + h - AREA_HEADER + (*scroll)
	gs.widgetW = w - SCROLL_AREA_PADDING*4
	g_scrollTop = y - AREA_HEADER + h
	g_scrollBottom = y + SCROLL_AREA_PADDING
	g_scrollRight = x + w - SCROLL_AREA_PADDING*3
	g_scrollVal = scroll

	g_scrollAreaTop = gs.widgetY

	g_focusTop = y - AREA_HEADER
	g_focusBottom = y - AREA_HEADER + h

	g_insideScrollArea = gs.inRect(x, y, w, h, false)
	gs.insideCurrentScroll = g_insideScrollArea

	gs.addGfxCmdRoundedRect(float64(x), float64(y), float64(w), float64(h), 6, imguiRGBA(0, 0, 0, 192))

	gs.addGfxCmdText(x+AREA_HEADER/2, y+h-AREA_HEADER/2-TEXT_HEIGHT/2, giu.AlignLeft, name, imguiRGBA(255, 255, 255, 128))

	gs.addGfxCmdScissor(x+SCROLL_AREA_PADDING, y+SCROLL_AREA_PADDING, w-SCROLL_AREA_PADDING*4, h-AREA_HEADER-SCROLL_AREA_PADDING)

	return g_insideScrollArea
}

func (gs *guiState) imguiEndScrollArea() {
	// Disable scissoring.
	gs.addGfxCmdScissor(-1, -1, -1, -1)

	// Draw scroll bar
	x := g_scrollRight + SCROLL_AREA_PADDING/2
	y := g_scrollBottom
	w := SCROLL_AREA_PADDING * 2
	h := g_scrollTop - g_scrollBottom

	stop := g_scrollAreaTop
	sbot := gs.widgetY
	sh := stop - sbot // The scrollable area height.

	barHeight := float64(h) / float64(sh)

	if barHeight < 1 {
		barY := float64(y-sbot) / float64(sh)
		if barY < 0 {
			barY = 0
		}
		if barY > 1 {
			barY = 1
		}

		// Handle scroll bar logic.
		hid := g_scrollId
		hx := x
		hy := y + int(barY*float64(h))
		hw := w
		hh := int(barHeight * float64(h))

		rangef := h - (hh - 1)
		over := gs.inRect(hx, hy, hw, hh)
		gs.buttonLogic(hid, over)
		if gs.IsActive(hid) {
			u := float64((hy - y)) / float64(rangef)
			if gs.wentActive {
				gs.dragY = gs.my
				gs.dragOrig = u
			}
			if gs.dragY != gs.my {
				u = gs.dragOrig + float64(gs.my-gs.dragY)/float64(rangef)
				if u < 0 {
					u = 0
				}
				if u > 1 {
					u = 1
				}
				*g_scrollVal = int((1 - u) * float64(sh-h))
			}
		}

		// BG
		gs.addGfxCmdRoundedRect(float64(x), float64(y), float64(w), float64(h), float64(w)/2-1, imguiRGBA(0, 0, 0, 196))
		// Bar
		if gs.IsActive(hid) {
			gs.addGfxCmdRoundedRect(float64(hx), float64(hy), float64(hw), float64(hh), float64(w)/2-1, imguiRGBA(255, 196, 0, 196))
		} else {
			c := imguiRGBA(255, 255, 255, 64)
			if gs.IsHot(hid) {
				c = imguiRGBA(255, 196, 0, 96)
			}
			gs.addGfxCmdRoundedRect(float64(hx), float64(hy), float64(hw), float64(hh), float64(w)/2-1, c)
		}

		// Handle mouse scrolling.
		if g_insideScrollArea { // && !anyActive())
			if gs.scroll != 0 {
				*g_scrollVal += 20 * gs.scroll
				if *g_scrollVal < 0 {
					*g_scrollVal = 0
				}
				if *g_scrollVal > (sh - h) {
					*g_scrollVal = (sh - h)
				}
			}
		}
	}
	gs.insideCurrentScroll = false
}

func (gs *guiState) imguiButton(text string, enableds ...bool) bool {
	enabled := true
	if len(enableds) > 0 {
		enabled = enableds[0]
	}
	gs.widgetId++
	id := (gs.areaId << 16) | gs.widgetId

	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	w := gs.widgetW
	h := BUTTON_HEIGHT
	gs.widgetY -= BUTTON_HEIGHT + DEFAULT_SPACING

	over := enabled && gs.inRect(x, y, w, h)
	res := gs.buttonLogic(id, over)
	c := uint8(96)
	if gs.IsActive(id) {
		c = 196
	}
	gs.addGfxCmdRoundedRect(float64(x), float64(y), float64(w), float64(h), float64(BUTTON_HEIGHT)/float64(2)-1, imguiRGBA(128, 128, 128, c))
	if enabled {
		c := imguiRGBA(255, 255, 255, 200)
		if gs.IsHot(id) {
			c = imguiRGBA(255, 196, 0, 255)
		}
		gs.addGfxCmdText(x+BUTTON_HEIGHT/2, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, c)

	} else {
		gs.addGfxCmdText(x+BUTTON_HEIGHT/2, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(128, 128, 128, 200))
	}

	return res
}

func (gs *guiState) imguiItem(text string, enabled bool) bool {
	gs.widgetId++
	id := (gs.areaId << 16) | gs.widgetId

	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	w := gs.widgetW
	h := BUTTON_HEIGHT
	gs.widgetY -= BUTTON_HEIGHT + DEFAULT_SPACING

	over := enabled && gs.inRect(x, y, w, h)
	res := gs.buttonLogic(id, over)

	if gs.IsHot(id) {
		c := uint8(96)
		if gs.IsActive(id) {
			c = 196
		}
		gs.addGfxCmdRoundedRect(float64(x), float64(y), float64(w), float64(h), 2.0, imguiRGBA(255, 196, 0, c))
	}

	if enabled {
		gs.addGfxCmdText(x+BUTTON_HEIGHT/2, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(255, 255, 255, 200))
	} else {
		gs.addGfxCmdText(x+BUTTON_HEIGHT/2, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(128, 128, 128, 200))
	}

	return res
}

func (gs *guiState) imguiCheck(text string, checked bool, enableds ...bool) bool {
	enabled := true
	if len(enableds) > 0 {
		enabled = enableds[0]
	}
	gs.widgetId++
	id := (gs.areaId << 16) | gs.widgetId

	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	w := gs.widgetW
	h := BUTTON_HEIGHT
	gs.widgetY -= BUTTON_HEIGHT + DEFAULT_SPACING

	over := enabled && gs.inRect(x, y, w, h)
	res := gs.buttonLogic(id, over)

	cx := x + BUTTON_HEIGHT/2 - CHECK_SIZE/2
	cy := y + BUTTON_HEIGHT/2 - CHECK_SIZE/2
	c := uint8(96)
	if gs.IsActive(id) {
		c = 196
	}
	gs.addGfxCmdRoundedRect(float64(cx)-3, float64(cy)-3, CHECK_SIZE+6, CHECK_SIZE+6, 4, imguiRGBA(128, 128, 128, c))
	if checked {
		if enabled {
			c := uint8(200)
			if gs.IsActive(id) {
				c = 255
			}
			gs.addGfxCmdRoundedRect(float64(cx), float64(cy), CHECK_SIZE, CHECK_SIZE, CHECK_SIZE/2-1, imguiRGBA(255, 255, 255, c))
		} else {
			gs.addGfxCmdRoundedRect(float64(cx), float64(cy), CHECK_SIZE, CHECK_SIZE, CHECK_SIZE/2-1, imguiRGBA(128, 128, 128, 200))
		}
	}

	if enabled {
		c := imguiRGBA(255, 255, 255, 200)
		if gs.IsHot(id) {
			c = imguiRGBA(255, 196, 0, 255)
		}
		gs.addGfxCmdText(x+BUTTON_HEIGHT, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, c)
	} else {
		gs.addGfxCmdText(x+BUTTON_HEIGHT, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(128, 128, 128, 200))
	}

	return res
}

func (gs *guiState) imguiCollapse(text string, subtext string, checked, enabled bool) bool {
	gs.widgetId++
	id := (gs.areaId << 16) | gs.widgetId

	x := gs.widgetX
	y := gs.widgetY - BUTTON_HEIGHT
	w := gs.widgetW
	h := BUTTON_HEIGHT
	gs.widgetY -= BUTTON_HEIGHT // + DEFAULT_SPACING;

	cx := x + BUTTON_HEIGHT/2 - CHECK_SIZE/2
	cy := y + BUTTON_HEIGHT/2 - CHECK_SIZE/2

	over := enabled && gs.inRect(x, y, w, h)
	res := gs.buttonLogic(id, over)

	if checked {
		c := uint8(200)
		if gs.IsActive(id) {
			c = 255
		}
		gs.addGfxCmdTriangle(cx, cy, CHECK_SIZE, CHECK_SIZE, 2, imguiRGBA(255, 255, 255, c))
	} else {
		c := uint8(200)
		if gs.IsActive(id) {
			c = 255
		}
		gs.addGfxCmdTriangle(cx, cy, CHECK_SIZE, CHECK_SIZE, 1, imguiRGBA(255, 255, 255, c))
	}

	if enabled {
		c := imguiRGBA(255, 255, 255, 200)
		if gs.IsHot(id) {
			c = imguiRGBA(255, 196, 0, 255)
		}
		gs.addGfxCmdText(x+BUTTON_HEIGHT, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, c)
	} else {
		gs.addGfxCmdText(x+BUTTON_HEIGHT, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignLeft, text, imguiRGBA(128, 128, 128, 200))
	}

	if subtext != "" {
		gs.addGfxCmdText(x+w-BUTTON_HEIGHT/2, y+BUTTON_HEIGHT/2-TEXT_HEIGHT/2, giu.AlignRight, subtext, imguiRGBA(255, 255, 255, 128))
	}

	return res
}
