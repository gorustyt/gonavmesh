package gui

import (
	"github.com/AllenDang/giu"
)

// ui界面的管理和控制
type layout struct {
	ui            *Gui
	showLog       bool
	showMenu      bool
	showTools     bool
	showLevels    bool
	showSample    bool
	showTestCases bool

	// Window scroll positions.
	propScroll  int
	logScroll   int
	toolsScroll int

	mouseOverMenu bool
}

func newLayout(ui *Gui) *layout {
	return &layout{ui: ui, showTools: true, showMenu: true}
}

func (l *layout) init() {

}
func (l *layout) Build() {
	var mouseScroll int
	//var processHitTest bool
	//var processHitTestShift bool
	point := giu.GetMousePos()
	l.ui.gs.beginFrame(point.X, point.Y, giu.IsMouseClicked(giu.MouseButtonLeft), mouseScroll)
	l.ShowMenu()
	l.ShowLog()
	l.ShowTools()
	l.ui.gs.endFrame()
}

func (l *layout) ShowMenu() {
	// Help text.
	if !l.showMenu {
		return
	}
	msg := "W/S/A/D: Move  RMB: Rotate"
	l.ui.gs.imguiDrawText(280, l.ui.height-20, giu.AlignLeft, msg, imguiRGBA(255, 255, 255, 128))
	if l.ui.gs.imguiBeginScrollArea("Properties", l.ui.width-250-10, 10, 250, l.ui.height-20, &l.propScroll) {
		l.mouseOverMenu = true
	}

	if l.ui.gs.imguiCheck("Show Log", l.showLog, true) {
		l.showLog = !l.showLog
	}

	if l.ui.gs.imguiCheck("Show Tools", l.showTools, true) {
		l.showTools = !l.showTools
	}

}

func (l *layout) ShowLog() {
	if !(l.showLog && l.showMenu) {
		return
	}
	if l.ui.gs.imguiBeginScrollArea("Log", 250+20, 10, l.ui.width-300-250, 200, &l.logScroll) {
		l.mouseOverMenu = true
	}

	for i := 0; i < l.ui.logger.getLogCount(); i++ {
		l.ui.gs.imguiLabel(l.ui.logger.getLogText(i))
	}

	l.ui.gs.imguiEndScrollArea()
}

func (l *layout) ShowTools() {
	// Left column tools menu
	if !(!l.showTestCases && l.showTools && l.showMenu) { // && geom && sample)
		return
	}
	if l.ui.gs.imguiBeginScrollArea("Tools", 10, 10, 250, l.ui.height-20, &l.toolsScroll) {
		l.mouseOverMenu = true
	}

	if l.ui.sample != nil {
		l.ui.sample.handleTools()
	}

	l.ui.gs.imguiEndScrollArea()
}
