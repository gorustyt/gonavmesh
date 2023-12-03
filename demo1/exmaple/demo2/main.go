package main

import (
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"github.com/go-gl/mathgl/mgl32"
	"gonavamesh/demo/lib/canvas"
	"image"
	"image/draw"
	_ "image/jpeg"
	_ "image/png"
	"log"
	"math"
	"os"
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

	//m := mgl32.Scale3D(0.5, 0.5, 0.5)
	axis := mgl32.Vec3{1.0, 0.0, 0}
	model := mgl32.HomogRotate3D(mgl32.DegToRad(-55), axis)
	view := mgl32.Translate3D(0, 0, -3)
	projection := mgl32.Perspective(mgl32.DegToRad(45), float32(windowWidth)/float32(windowHeight), 0.1, 100)
	window.MakeContextCurrent()
	// Initialize Glow
	if err := gl.Init(); err != nil {
		panic(err)
	}
	vers := []float32{
		//     ---- 位置 ----       ---- 颜色 ----     - 纹理坐标 -
		0.5, 0.5, 0.0, 1.0, 0.0, 0.0, 1, 1, // 右上
		0.5, -0.5, 0.0, 0.0, 1.0, 0.0, 1, 0, // 右下
		-0.5, -0.5, 0.0, 0.0, 0.0, 1.0, 0, 0, // 左下
		-0.5, 0.5, 0.0, 1.0, 1.0, 0.0, 0, 1, // 左上
	}
	indices := []uint32{
		0, 1, 3, // first triangle
		1, 2, 3, // second triangle
	}

	t := makeTexture("./demo1/assets/awesomeface.png", gl.TEXTURE0)
	t1 := makeTexture("./demo1/assets/container.jpg", gl.TEXTURE1)
	vao := canvas.MakeVao(vers)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 8*4, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointer(1, 3, gl.FLOAT, false, 8*4, gl.PtrOffset(3*4))
	gl.EnableVertexAttribArray(1)
	gl.VertexAttribPointer(2, 2, gl.FLOAT, false, 8*4, gl.PtrOffset(6*4))
	gl.EnableVertexAttribArray(2)
	prg, err := canvas.NewProgram(vertexShader, fragmentShader)
	if err != nil {
		panic(err)
	}
	gl.UseProgram(prg)
	gl.Uniform1i(gl.GetUniformLocation(prg, gl.Str("ourTexture1\x00")), 0)
	gl.Uniform1i(gl.GetUniformLocation(prg, gl.Str("ourTexture2\x00")), 1)
	gl.Enable(gl.DEPTH_TEST)
	gl.DepthFunc(gl.LESS)
	tt := float32(0.2)
	for !window.ShouldClose() {
		gl.ClearColor(0.2, 0.3, 0.3, 1.0)
		gl.UseProgram(prg)
		if window.GetKey(glfw.KeySpace) == glfw.Press {

		} else if window.GetKey(glfw.KeyUp) == glfw.Press {
			tt = float32(math.Min(1.0, float64(tt)+0.001))
		} else if window.GetKey(glfw.KeyDown) == glfw.Press {
			tt = float32(math.Max(0.1, float64(tt)-0.001))
		}
		gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("model\x00")),
			1, false, &model[0])
		gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("view\x00")),
			1, false, &view[0])
		gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("projection\x00")),
			1, false, &projection[0])
		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("t\x00")), tt)
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		gl.ActiveTexture(t)
		gl.BindTexture(gl.TEXTURE_2D, t)
		gl.ActiveTexture(t1)
		gl.BindTexture(gl.TEXTURE_2D, t1)
		gl.BindVertexArray(vao)
		gl.DrawElements(gl.TRIANGLES, 6, gl.UNSIGNED_INT, gl.Ptr(indices))
		// Maintenance
		window.SwapBuffers()
		glfw.PollEvents()
	}
}

var vertexShader = `
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 color;
layout (location = 2) in vec2 texCoord;

out vec3 ourColor;
out vec2 TexCoord;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main()
{
	gl_Position = projection*view*model*vec4(position, 1.0f);
	ourColor = color;

	TexCoord = vec2(texCoord.x, 1.0 - texCoord.y);
}
` + "\x00"
var fragmentShader = `
#version 330 core
in vec3 ourColor;
in vec2 TexCoord;

out vec4 color;

uniform sampler2D ourTexture1;
uniform sampler2D ourTexture2;
uniform float t;
void main()
{
	color = mix(texture(ourTexture1, TexCoord), texture(ourTexture2, TexCoord), t);
}
` + "\x00"

func makeTexture(p string, index uint32) (te uint32) {
	f, err := os.Open(p)
	if err != nil {
		panic(err)
	}
	img, _, err := image.Decode(f)
	if err != nil {
		panic(err)
	}

	rgba := image.NewRGBA(img.Bounds())
	if rgba.Stride != rgba.Rect.Size().X*4 {
		panic("unsupported stride")
	}

	draw.Draw(rgba, rgba.Bounds(), img, image.Point{0, 0}, draw.Src)
	defer f.Close()
	gl.GenTextures(1, &te)
	gl.ActiveTexture(index)
	gl.BindTexture(gl.TEXTURE_2D, te)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST)
	gl.TexImage2D(gl.TEXTURE_2D, 0, gl.RGBA, int32(rgba.Rect.Size().X), int32(rgba.Rect.Size().Y), 0, gl.RGBA, gl.UNSIGNED_BYTE, gl.Ptr(rgba.Pix))
	gl.GenerateMipmap(gl.TEXTURE_2D)
	return te
}
