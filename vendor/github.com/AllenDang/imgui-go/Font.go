package imgui

// #include "FontWrapper.h"
import "C"

// Font describes one loaded font in an atlas.
type Font uintptr

// DefaultFont can be used to refer to the default font of the current font atlas without
// having the actual font reference.
const DefaultFont Font = 0

func (font Font) handle() C.IggFont {
	return C.IggFont(font)
}

func (font Font) FindGlyph(c rune) C.IggFontGlyph {
	return C.iggFontFindGlyph(font.handle(), C.uint(c))
}
