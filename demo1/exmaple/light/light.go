package main

import (
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"github.com/go-gl/mathgl/mgl32"
	"gonavamesh/demo/lib/canvas"
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
	prg, err := canvas.NewProgram(vertexShader, fragmentShader)
	if err != nil {
		panic(err)
	}
	lightPrg, err := canvas.NewProgram(vertexShader, lightFragmentShader)
	if err != nil {
		panic(err)
	}
	vs := []float32{
		-0.5, -0.5, -0.5,
		0.5, -0.5, -0.5,
		0.5, 0.5, -0.5,
		0.5, 0.5, -0.5,
		-0.5, 0.5, -0.5,
		-0.5, -0.5, -0.5,

		-0.5, -0.5, 0.5,
		0.5, -0.5, 0.5,
		0.5, 0.5, 0.5,
		0.5, 0.5, 0.5,
		-0.5, 0.5, 0.5,
		-0.5, -0.5, 0.5,

		-0.5, 0.5, 0.5,
		-0.5, 0.5, -0.5,
		-0.5, -0.5, -0.5,
		-0.5, -0.5, -0.5,
		-0.5, -0.5, 0.5,
		-0.5, 0.5, 0.5,

		0.5, 0.5, 0.5,
		0.5, 0.5, -0.5,
		0.5, -0.5, -0.5,
		0.5, -0.5, -0.5,
		0.5, -0.5, 0.5,
		0.5, 0.5, 0.5,

		-0.5, -0.5, -0.5,
		0.5, -0.5, -0.5,
		0.5, -0.5, 0.5,
		0.5, -0.5, 0.5,
		-0.5, -0.5, 0.5,
		-0.5, -0.5, -0.5,

		-0.5, 0.5, -0.5,
		0.5, 0.5, -0.5,
		0.5, 0.5, 0.5,
		0.5, 0.5, 0.5,
		-0.5, 0.5, 0.5,
		-0.5, 0.5, -0.5,
	}
	lightColor := mgl32.Vec3{1.0, 1.0, 1.0}
	lightPos := mgl32.Vec3{1.2, 1.0, 2.0}
	objectColor := mgl32.Vec3{1.0, 0.5, 0.31}
	model := mgl32.Scale3D(0.2, 0.2, 0.2)
	model := model.Mul4(mgl32.Translate3D(lightPos[0], lightPos[1], lightPos[2]))
	view :=
		gl.Uniform1fv(gl.GetUniformLocation(prg, gl.Str("objectColor\x00")), 1, &objectColor[0])
	gl.Uniform1fv(gl.GetUniformLocation(prg, gl.Str("lightColor\x00")), 1, &lightColor[0])
	lightVao := MakeVaoLight()
	objVao := canvas.MakeVao(vs)
	gl.Enable(gl.DEPTH_TEST)
	gl.DepthFunc(gl.LESS)
	for !window.ShouldClose() {
		gl.ClearColor(0.2, 0.3, 0.3, 1.0)
		gl.UseProgram(prg)
		// Maintenance
		window.SwapBuffers()
		glfw.PollEvents()
	}
}

var vertexShader = `
#version 330 core
layout (location = 0) in vec3 position;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main()
{
	gl_Position = projection*view*model*vec4(position, 1.0f);

}
` + "\x00"
var fragmentShader = `
#version 330 core
out vec4 FragColor;

uniform vec3 objectColor;
uniform vec3 lightColor;

void main()
{
    FragColor = vec4(lightColor * objectColor, 1.0);
}
` + "\x00"

var lightFragmentShader = `
#version 330 core
out vec4 FragColor;

void main()
{
    FragColor = vec4(1.0); // 将向量的四个分量全部设置为1.0
}
` + "\x00"

func MakeVaoLight() uint32 {
	var vbo uint32

	// 在显卡中开辟一块空间，创建顶点缓存对象，个数为1，变量vbo会被赋予一个ID值。
	gl.GenBuffers(1, &vbo)

	// 将 vbo 赋值给 gl.ARRAY_BUFFER，要知道这个对象会被赋予不同的vbo，因此其值是变化的
	// 可选类型：GL_ARRAY_BUFFER, GL_ELEMENT_ARRAY_BUFFER, GL_PIXEL_PACK_BUFFER, GL_PIXEL_UNPACK_BUFFER
	gl.BindBuffer(gl.ARRAY_BUFFER, vbo)

	// 将内存中的数据传递到显卡中的gl.ARRAY_BUFFER对象上，其实是把数据传递到绑定在其上面的vbo对象上。
	// 4*len(points) 代表总的字节数，因为是32位的
	//gl.BufferData(gl.ARRAY_BUFFER, 3*len(points), gl.Ptr(points), gl.STATIC_DRAW)

	var vao uint32
	// 创建顶点数组对象，个数为1，变量vao会被赋予一个ID值。
	gl.GenVertexArrays(1, &vao)

	// 后面的两个函数都是要操作具体的vao的，因此需要先将vao绑定到opengl上。
	// 解绑：gl.BindVertexArray(0)，opengl中很多的解绑操作都是传入0
	gl.BindVertexArray(vao)
	return vao
}
