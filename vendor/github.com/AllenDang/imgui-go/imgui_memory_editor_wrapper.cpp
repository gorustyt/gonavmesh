#include "imgui.h"
#include "imgui_memory_editor.h"
#include "imgui_memory_editor_wrapper.h"

IggMemoryEditor IggNewMemoryEditor()
{
  MemoryEditor *mem_editor = new MemoryEditor();
  return static_cast<IggMemoryEditor>(mem_editor);
}

void IggMemoryEditorDrawContents(IggMemoryEditor handle, void *mem_data, size_t mem_size)
{
  MemoryEditor *mem_editor = reinterpret_cast<MemoryEditor *>(handle);
  mem_editor->DrawContents(mem_data, mem_size);
}
