#pragma once

#include "imguiWrapperTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

extern unsigned int iggDrawCommandIdxOffset(IggDrawCmd handle);
extern void iggDrawCommandGetVertexOffset(IggDrawCmd handle, unsigned int *count);
extern void iggDrawCommandGetElementCount(IggDrawCmd handle, unsigned int *count);
extern void iggDrawCommandGetClipRect(IggDrawCmd handle, IggVec4 *rect);
extern void iggDrawCommandGetTextureID(IggDrawCmd handle, IggTextureID *id);
extern IggBool iggDrawCommandHasUserCallback(IggDrawCmd handle);
extern void iggDrawCommandCallUserCallback(IggDrawCmd handle, IggDrawList listHandle);

#ifdef __cplusplus
}
#endif
