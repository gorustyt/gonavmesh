package imgui

const (
	// ComboFlagsNone default = 0
	ComboFlagsNone = 0
	// ComboFlagsPopupAlignLeft aligns the popup toward the left by default.
	ComboFlagsPopupAlignLeft = 1 << 0
	// ComboFlagsHeightSmall has max ~4 items visible.
	// Tip: If you want your combo popup to be a specific size you can use SetNextWindowSizeConstraints() prior to calling BeginCombo().
	ComboFlagsHeightSmall = 1 << 1
	// ComboFlagsHeightRegular has max ~8 items visible (default).
	ComboFlagsHeightRegular = 1 << 2
	// ComboFlagsHeightLarge has max ~20 items visible.
	ComboFlagsHeightLarge = 1 << 3
	// ComboFlagsHeightLargest has as many fitting items as possible.
	ComboFlagsHeightLargest = 1 << 4
	// ComboFlagsNoArrowButton displays on the preview box without the square arrow button.
	ComboFlagsNoArrowButton = 1 << 5
	// ComboFlagsNoPreview displays only a square arrow button.
	ComboFlagsNoPreview = 1 << 6

	ComboFlagsHeightMask = ComboFlagsHeightSmall | ComboFlagsHeightRegular | ComboFlagsHeightLarge | ComboFlagsHeightLargest
)
