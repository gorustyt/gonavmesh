package main

import (
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"gonavamesh/demo/lib/canvas"
	"gonavamesh/demo/lib/glfont"
	"log"
	"runtime"
)

const windowWidth = 800
const windowHeight = 600

func init() {
	// GLFW event handling must run on the main OS thread
	runtime.LockOSThread()
}
func main() {
	if err := glfw.Init(); err != nil {
		log.Fatalln("failed to initialize glfw:", err)
	}
	defer glfw.Terminate()

	glfw.WindowHint(glfw.Resizable, glfw.False)
	glfw.WindowHint(glfw.ContextVersionMajor, 4)
	glfw.WindowHint(glfw.ContextVersionMinor, 1)
	glfw.WindowHint(glfw.OpenGLProfile, glfw.OpenGLCoreProfile)
	glfw.WindowHint(glfw.OpenGLForwardCompatible, glfw.True)
	window, err := glfw.CreateWindow(windowWidth, windowHeight, "Cube", nil, nil)
	if err != nil {
		panic(err)
	}
	window.MakeContextCurrent()

	// Initialize Glow
	if err := gl.Init(); err != nil {
		panic(err)
	}

	version := gl.GoStr(gl.GetString(gl.VERSION))
	fmt.Println("OpenGL version", version)
	font, err := glfont.LoadFont("./demo/bin/DroidSans.ttf", 50, windowWidth, windowHeight)
	if err != nil {
		panic(err)
	}
	t := canvas.NewTriangle()
	for !window.ShouldClose() {
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		t.Draw([]float32{
			-0.5, -0.5, 0.0, // Left
			0.5, -0.5, 0.0, // Right
			0.0, 0.5, 0.0, // Top
		})
		font.Printf(100, 100, 1.0, "helllo wolrd ")
		// Maintenance
		window.SwapBuffers()
		glfw.PollEvents()
	}
}
