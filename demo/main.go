package main

import (
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"gonavamesh/demo/gui"
)

func init() {
}

func main() {
	gui.InitGui()
	defer glfw.Terminate()
	window := gui.CreateWindow()
	for !window.ShouldClose() {
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		gui.Run(window)
		// Do OpenGL stuff.
		window.SwapBuffers()
		glfw.PollEvents()
	}
}
