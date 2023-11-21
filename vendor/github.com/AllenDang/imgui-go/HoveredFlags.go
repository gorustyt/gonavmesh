package imgui

const (
	HoveredFlagsNone             = 0      // Return true if directly over the item/window, not obstructed by another window, not obstructed by an active popup or modal blocking inputs under them.
	HoveredFlagsChildWindows     = 1 << 0 // IsWindowHovered() only: Return true if any children of the window is hovered
	HoveredFlagsRootWindow       = 1 << 1 // IsWindowHovered() only: Test from root window (top most parent of the current hierarchy)
	HoveredFlagsAnyWindow        = 1 << 2 // IsWindowHovered() only: Return true if any window is hovered
	HoveredFlagsNoPopupHierarchy = 1 << 3 // IsWindowHovered() only: Do not consider popup hierarchy (do not treat popup emitter as parent of popup) (when used with ChildWindows or RootWindow)
	//HoveredFlagsDockHierarchy               = 1 << 4   // IsWindowHovered() only: Consider docking hierarchy (treat dockspace host as parent of docked window) (when used with ChildWindows or RootWindow)
	HoveredFlagsAllowWhenBlockedByPopup = 1 << 5 // Return true even if a popup window is normally blocking access to this item/window
	//HoveredFlagsAllowWhenBlockedByModal     = 1 << 6   // Return true even if a modal popup window is normally blocking access to this item/window. FIXME-TODO: Unavailable yet.
	HoveredFlagsAllowWhenBlockedByActiveItem = 1 << 7 // Return true even if an active item is blocking access to this item/window. Useful for Drag and Drop patterns.
	HoveredFlagsAllowWhenOverlapped          = 1 << 8 // IsItemHovered() only: Return true even if the position is obstructed or overlapped by another window
	HoveredFlagsAllowWhenDisabled            = 1 << 9 // IsItemHovered() only: Return true even if the item is disabled
	HoveredFlagsRectOnly                     = HoveredFlagsAllowWhenBlockedByPopup | HoveredFlagsAllowWhenBlockedByActiveItem | HoveredFlagsAllowWhenOverlapped
	HoveredFlagsRootAndChildWindows          = HoveredFlagsRootWindow | HoveredFlagsChildWindows
)
