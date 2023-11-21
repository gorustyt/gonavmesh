#include "imguiWrappedHeader.h"
#include "FontWrapper.h"

IggFontGlyph iggFontFindGlyph(IggFont handle, unsigned int c)
{
  ImFont *font = reinterpret_cast<ImFont *>(handle);
  return static_cast<IggFontGlyph>(const_cast<ImFontGlyph *>(font->FindGlyph(c)));
}
