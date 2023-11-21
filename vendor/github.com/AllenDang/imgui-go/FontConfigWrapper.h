#pragma once

#include "imguiWrapperTypes.h"

#ifdef __cplusplus
extern "C"
{
#endif

extern IggFontConfig iggNewFontConfig();
extern void iggFontConfigDelete(IggFontConfig handle);

extern void iggFontConfigSetSize(IggFontConfig handle, float sizePixels);
extern void iggFontConfigSetOversampleH(IggFontConfig handle, int value);
extern void iggFontConfigSetOversampleV(IggFontConfig handle, int value);
extern void iggFontConfigSetPixelSnapH(IggFontConfig handle, IggBool value);
extern void iggFontConfigSetGlyphMinAdvanceX(IggFontConfig handle, float value);
extern void iggFontConfigSetGlyphMaxAdvanceX(IggFontConfig handle, float value);
extern void iggFontConfigSetMergeMode(IggFontConfig handle, IggBool value);
extern int iggFontConfigGetFontDataOwnedByAtlas(IggFontConfig handle);
extern void iggFontConfigSetRasterizerMultiply(IggFontConfig handle, float value);
extern unsigned int iggFontConfigGetFontBuilderFlags(IggFontConfig handle);
extern void iggFontConfigSetFontBuilderFlags(IggFontConfig handle, unsigned int flags);

#ifdef __cplusplus
}
#endif
