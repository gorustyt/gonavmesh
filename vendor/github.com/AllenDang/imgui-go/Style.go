package imgui

// #include "StyleWrapper.h"
import "C"

// StyleVarID identifies a style variable in the UI style.
type StyleVarID int

const (
	StyleVarAlpha StyleVarID = iota
	StyleVarDisabledAlpha
	StyleVarWindowPadding
	StyleVarWindowRounding
	StyleVarWindowBorderSize
	StyleVarWindowMinSize
	StyleVarWindowTitleAlign
	StyleVarChildRounding
	StyleVarChildBorderSize
	StyleVarPopupRounding
	StyleVarPopupBorderSize
	StyleVarFramePadding
	StyleVarFrameRounding
	StyleVarFrameBorderSize
	StyleVarItemSpacing
	StyleVarItemInnerSpacing
	StyleVarIndentSpacing
	StyleVarCellPadding
	StyleVarScrollbarSize
	StyleVarScrollbarRounding
	StyleVarGrabMinSize
	StyleVarGrabRounding
	StyleVarTabRounding
	StyleVarButtonTextAlign
	StyleVarSelectableTextAlign
)

// StyleColorID identifies a color in the UI style.
type StyleColorID int

// StyleColor identifier
const (
	StyleColorText StyleColorID = iota
	StyleColorTextDisabled
	StyleColorWindowBg // Background of normal windows
	StyleColorChildBg  // Background of child windows
	StyleColorPopupBg  // Background of popups, menus, tooltips windows
	StyleColorBorder
	StyleColorBorderShadow
	StyleColorFrameBg // Background of checkbox, radio button, plot, slider, text input
	StyleColorFrameBgHovered
	StyleColorFrameBgActive
	StyleColorTitleBg
	StyleColorTitleBgActive
	StyleColorTitleBgCollapsed
	StyleColorMenuBarBg
	StyleColorScrollbarBg
	StyleColorScrollbarGrab
	StyleColorScrollbarGrabHovered
	StyleColorScrollbarGrabActive
	StyleColorCheckMark
	StyleColorSliderGrab
	StyleColorSliderGrabActive
	StyleColorButton
	StyleColorButtonHovered
	StyleColorButtonActive
	StyleColorHeader // Header* colors are used for CollapsingHeader, TreeNode, Selectable, MenuItem
	StyleColorHeaderHovered
	StyleColorHeaderActive
	StyleColorSeparator
	StyleColorSeparatorHovered
	StyleColorSeparatorActive
	StyleColorResizeGrip
	StyleColorResizeGripHovered
	StyleColorResizeGripActive
	StyleColorTab
	StyleColorTabHovered
	StyleColorTabActive
	StyleColorTabUnfocused
	StyleColorTabUnfocusedActive
	StyleColorPlotLines
	StyleColorPlotLinesHovered
	StyleColorPlotHistogram
	StyleColorPlotHistogramHovered
	StyleColorTableHeaderBg     // Table header background
	StyleColorTableBorderStrong // Table outer and header borders (prefer using Alpha=1.0 here)
	StyleColorTableBorderLight  // Table inner borders (prefer using Alpha=1.0 here)
	StyleColorTableRowBg        // Table row background (even rows)
	StyleColorTableRowBgAlt     // Table row background (odd rows)
	StyleColorTextSelectedBg
	StyleColorDragDropTarget
	StyleColorNavHighlight          // Gamepad/keyboard: current highlighted item
	StyleColorNavWindowingHighlight // Highlight window when using CTRL+TAB
	StyleColorNavWindowingDimBg     // Darken/colorize entire screen behind the CTRL+TAB window list, when active
	StyleColorModalWindowDimBg      // Darken/colorize entire screen behind a modal window, when one is active
)

// Style describes the overall graphical representation of the user interface.
type Style uintptr

func (style Style) handle() C.IggGuiStyle {
	return C.IggGuiStyle(style)
}

// ItemSpacing is the horizontal and vertical spacing between widgets/lines.
func (style Style) ItemSpacing() Vec2 {
	var value Vec2
	valueArg, valueFin := value.wrapped()
	C.iggStyleGetItemSpacing(style.handle(), valueArg)
	valueFin()
	return value
}

// ItemInnerSpacing is the horizontal and vertical spacing between elements of
// a composed widget (e.g. a slider and its label).
func (style Style) ItemInnerSpacing() Vec2 {
	var value Vec2
	valueArg, valueFin := value.wrapped()
	C.iggStyleGetItemInnerSpacing(style.handle(), valueArg)
	valueFin()
	return value
}

func (style Style) WindowPadding() Vec2 {
	var value Vec2
	valueArg, valueFin := value.wrapped()
	C.iggStyleGetWindowPadding(style.handle(), valueArg)
	valueFin()
	return value
}

func (style Style) FramePadding() Vec2 {
	var value Vec2
	valueArg, valueFin := value.wrapped()
	C.iggStyleGetFramePadding(style.handle(), valueArg)
	valueFin()
	return value
}

// SetColor sets a color value of the UI style.
func (style Style) SetColor(id StyleColorID, value Vec4) {
	valueArg, _ := value.wrapped()
	C.iggStyleSetColor(style.handle(), C.int(id), valueArg)
}

func (style Style) GetColor(id StyleColorID) Vec4 {
	var col Vec4
	colArg, colFin := col.wrapped()
	C.iggStyleGetColor(style.handle(), C.int(id), colArg)
	colFin()
	return col
}

// ScaleAllSizes applies a scaling factor to all sizes.
// To scale your entire UI (e.g. if you want your app to use High DPI or generally be DPI aware) you may use this helper function.
// Scaling the fonts is done separately and is up to you.
//
// Important: This operation is lossy because all sizes are rounded to integer.
// If you need to change your scale multiples, call this over a freshly initialized style rather than scaling multiple times.
func (style Style) ScaleAllSizes(scale float32) {
	C.iggStyleScaleAllSizes(style.handle(), C.float(scale))
}
