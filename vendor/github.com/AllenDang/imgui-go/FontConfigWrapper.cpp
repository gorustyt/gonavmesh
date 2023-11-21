#include "imguiWrappedHeader.h"
#include "FontConfigWrapper.h"

IggFontConfig iggNewFontConfig()
{
   ImFontConfig *fontConfig = new ImFontConfig();
   return static_cast<IggFontConfig>(fontConfig);
}

void iggFontConfigDelete(IggFontConfig handle)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   delete fontConfig;
}

void iggFontConfigSetSize(IggFontConfig handle, float sizePixels)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->SizePixels = sizePixels;
}

void iggFontConfigSetOversampleH(IggFontConfig handle, int value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->OversampleH = value;
}

void iggFontConfigSetOversampleV(IggFontConfig handle, int value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->OversampleV = value;
}

void iggFontConfigSetPixelSnapH(IggFontConfig handle, IggBool value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->PixelSnapH = value;
}

void iggFontConfigSetGlyphMinAdvanceX(IggFontConfig handle, float value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->GlyphMinAdvanceX = value;
}

void iggFontConfigSetGlyphMaxAdvanceX(IggFontConfig handle, float value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->GlyphMaxAdvanceX = value;
}

void iggFontConfigSetMergeMode(IggFontConfig handle, IggBool value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->MergeMode = value;
}

int iggFontConfigGetFontDataOwnedByAtlas(IggFontConfig handle)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   return fontConfig->FontDataOwnedByAtlas;
}

void iggFontConfigSetRasterizerMultiply(IggFontConfig handle, float value)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->RasterizerMultiply = value;
}

unsigned int iggFontConfigGetFontBuilderFlags(IggFontConfig handle)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   return fontConfig->FontBuilderFlags;
}

void iggFontConfigSetFontBuilderFlags(IggFontConfig handle, unsigned int flags)
{
   ImFontConfig *fontConfig = reinterpret_cast<ImFontConfig *>(handle);
   fontConfig->FontBuilderFlags = flags;
}
