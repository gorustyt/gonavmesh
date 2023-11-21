package imgui

import (
	"image"

	"github.com/go-gl/glfw/v3.3/glfw"
)

type Platform interface {
	// ShouldStop is regularly called as the abort condition for the program loop.
	ShouldStop() bool

	// SetShouldStop sets whether window should be closed
	SetShouldStop(bool)

	// ProcessEvents is called once per render loop to dispatch any pending events.
	ProcessEvents()

	// DisplaySize returns the dimension of the display.
	DisplaySize() [2]float32

	// FramebufferSize returns the dimension of the framebuffer.
	FramebufferSize() [2]float32

	// NewFrame marks the begin of a render pass. It must update the imgui IO state according to user input (mouse, keyboard, ...)
	NewFrame()

	// PostRender marks the completion of one render pass. Typically this causes the display buffer to be swapped.
	PostRender()

	// Dispose
	Dispose()

	// Set size change callback
	SetSizeChangeCallback(func(width, height int))

	// Set pos change callback
	SetPosChangeCallback(func(x, y int))

	// Set drop callback
	SetDropCallback(func(names []string))

	// Set input callback
	SetInputCallback(func(key glfw.Key, mods glfw.ModifierKey, action glfw.Action))

	// Set close callback, returned value will be used to close or cancel the window
	SetCloseCallback(func() bool)

	// Force Update
	Update()

	// Get content from system clipboard
	GetClipboard() string

	// Set content to system clipboard
	SetClipboard(content string)

	// Get the event pulling ticks per second
	GetTPS() int

	// Set the event pulling ticks per second
	SetTPS(tps int)

	// Set icon to master window
	SetIcon(icons []image.Image)

	// SetSizeLimits sets the size limits of the client area of the specified window.
	SetSizeLimits(minw, minh, maxw, maxh int)

	// SetTitle sets the title of platform window.
	SetTitle(title string)

	// Get window position
	GetPos() (x, y int)

	// Get DPI scale factor
	GetContentScale() float32

	// Check whehter window is minimized
	IsMinimized() bool

	// Check whether window is visible
	IsVisible() bool
}
