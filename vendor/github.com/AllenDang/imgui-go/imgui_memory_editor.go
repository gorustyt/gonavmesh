package imgui

// #include "imgui_memory_editor_wrapper.h"
import "C"
import "unsafe"

type MemoryEditor uintptr

func NewMemoryEditor() MemoryEditor {
	handle := C.IggNewMemoryEditor()
	return MemoryEditor(handle)
}

func (me MemoryEditor) handle() C.IggMemoryEditor {
	return C.IggMemoryEditor(me)
}

func (me MemoryEditor) DrawContents(data []uint8) {
	C.IggMemoryEditorDrawContents(me.handle(), unsafe.Pointer(&data[0]), C.size_t(len(data)))
}
