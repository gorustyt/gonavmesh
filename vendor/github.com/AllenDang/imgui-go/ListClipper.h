#pragma once

#include "imguiWrapperTypes.h"

#ifdef __cplusplus
extern "C" {
#endif

extern IggListClipper iggNewListClipper();

extern int iggListClipperDisplayStart(IggListClipper handle);
extern int iggListClipperDisplayEnd(IggListClipper handle);

extern void iggListClipperDelete(IggListClipper handle);
extern IggBool iggListClipperStep(IggListClipper handle);
extern void iggListClipperBegin(IggListClipper handle, int items_count, float items_height);
extern void iggListClipperEnd(IggListClipper handle);

#ifdef __cplusplus
}
#endif
