package imgui

type GocusedFlags int

const (
	FocusedFlagsNone             = 0
	FocusedFlagsChildWindows     = 1 << 0 // Return true if any children of the window is focused
	FocusedFlagsRootWindow       = 1 << 1 // Test from root window (top most parent of the current hierarchy)
	FocusedFlagsAnyWindow        = 1 << 2 // Return true if any window is focused. Important: If you are trying to tell how to dispatch your low-level inputs, do NOT use this. Use 'io.WantCaptureMouse' instead! Please read the FAQ!
	FocusedFlagsNoPopupHierarchy = 1 << 3 // Do not consider popup hierarchy (do not treat popup emitter as parent of popup) (when used with ChildWindows or RootWindow)
	//FocusedFlagsDockHierarchy               = 1 << 4   // Consider docking hierarchy (treat dockspace host as parent of docked window) (when used with ChildWindows or RootWindow)
	FocusedFlagsRootAndChildWindows = FocusedFlagsRootWindow | FocusedFlagsChildWindows
)
