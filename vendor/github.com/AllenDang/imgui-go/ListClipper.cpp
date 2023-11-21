#include "ListClipper.h"
#include "WrapperConverter.h"
#include "imguiWrappedHeader.h"

IggListClipper iggNewListClipper() {
  ImGuiListClipper *clipper = new ImGuiListClipper();
  return static_cast<IggListClipper>(clipper);
}

int iggListClipperDisplayStart(IggListClipper handle) {
  ImGuiListClipper *clipper = reinterpret_cast<ImGuiListClipper *>(handle);
  return clipper->DisplayStart;
}

int iggListClipperDisplayEnd(IggListClipper handle) {
  ImGuiListClipper *clipper = reinterpret_cast<ImGuiListClipper *>(handle);
  return clipper->DisplayEnd;
}

void iggListClipperDelete(IggListClipper handle) {
  ImGuiListClipper *clipper = reinterpret_cast<ImGuiListClipper *>(handle);
  delete clipper;
}

IggBool iggListClipperStep(IggListClipper handle) {
  ImGuiListClipper *clipper = reinterpret_cast<ImGuiListClipper *>(handle);
  return clipper->Step();
}

void iggListClipperBegin(IggListClipper handle, int items_count, float items_height) {
  ImGuiListClipper *clipper = reinterpret_cast<ImGuiListClipper *>(handle);
  clipper->Begin(items_count, items_height);
}

void iggListClipperEnd(IggListClipper handle) {
  ImGuiListClipper *clipper = reinterpret_cast<ImGuiListClipper *>(handle);
  clipper->End();
}
