package gui

import (
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"log"
	"math"
)

func InitGui() {
	err := glfw.Init()
	if err != nil {
		panic(err)
	}
	glfw.WindowHint(glfw.Resizable, glfw.False)
	glfw.WindowHint(glfw.ContextVersionMajor, 4)
	glfw.WindowHint(glfw.ContextVersionMinor, 1)
	glfw.WindowHint(glfw.OpenGLProfile, glfw.OpenGLCoreProfile)
	glfw.WindowHint(glfw.OpenGLForwardCompatible, glfw.True)
}

func CreateWindow() (window *glfw.Window) {

	var err error
	w, h := getWH()
	window, err = glfw.CreateWindow(w, h, "Testing", nil, nil)
	if err != nil {
		panic(err)
	}
	window.MakeContextCurrent()
	registerEvent(window)
	// Initialize Glow
	if err = gl.Init(); err != nil {
		panic(err)
	}
	version := gl.GoStr(gl.GetString(gl.VERSION))
	log.Println("OpenGL version", version)
	gl.Enable(gl.CULL_FACE)
	gl.DepthFunc(gl.LEQUAL)
	InitImgui()
	return window
}

func getWH() (width, height int) {
	sw := glfw.GetPrimaryMonitor().GetVideoMode().Width
	sh := glfw.GetPrimaryMonitor().GetVideoMode().Height
	aspect := 16.0 / 9.0
	width = int(math.Min(float64(sw), float64(sh)*aspect)) - 80
	height = sh - 80
	return width, height
}

func registerEvent(window *glfw.Window) {
	//文件拖拽事件
	window.SetDropCallback(func(w *glfw.Window, names []string) {
		fmt.Println("文件拖拽", names)
	})
	//鼠标事件
	window.SetMouseButtonCallback(func(w *glfw.Window, button glfw.MouseButton, action glfw.Action, mods glfw.ModifierKey) {
		switch action {
		case glfw.Press:
		case glfw.Release:
		}
	})
	//鼠标滚轮或者触摸板，鼠标滚轮只有yoff，表示垂直滚动了多少，触摸板有xoff和yoff。
	window.SetScrollCallback(func(w *glfw.Window, xoff float64, yoff float64) {
		log.Printf("滚动了======%v\n", yoff)
	})
	//键盘事件
	window.SetKeyCallback(func(w *glfw.Window, key glfw.Key, scancode int, action glfw.Action, mods glfw.ModifierKey) {
		switch key {
		case glfw.Key9:
			fmt.Println("================")
		}
	})
	//关闭事件
	window.SetCloseCallback(func(w *glfw.Window) {
		log.Printf("shuttdown ...")
	})
}
