package gui

import (
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"image"
	"image/draw"
	"log"
	"math"
	"os"
	"strings"
)

type UpdateWindow interface {
	Update()
}

type Gui struct {
	title         string
	gs            *guiState
	render        *imguiRender
	layout        *layout
	logger        *logger
	width, height int
	window        *glfw.Window
	updates       []UpdateWindow
}

func NewGui(title string) *Gui {
	ui := &Gui{
		gs:     newGuiState(),
		title:  title,
		logger: newLogger(),
	}
	ui.layout = newLayout(ui)
	ui.render = newImguiRender(ui)
	ui.CreateWindow()
	ui.updates = append(ui.updates, ui.layout)
	return ui
}

func (ui *Gui) shutdown(w *glfw.Window) {
	log.Printf("ui %v shuttdown.", ui.title)
}

func (ui *Gui) Run() {
	for !ui.window.ShouldClose() {
		for _, v := range ui.updates {
			v.Update()
		}
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		// Maintenance
		ui.window.SwapBuffers()
		glfw.PollEvents()
	}
}

func (ui *Gui) initGui() {
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

func (ui *Gui) CreateWindow() (window *glfw.Window) {
	ui.initGui()
	var err error
	ui.width, ui.height = ui.getWH()
	window, err = glfw.CreateWindow(ui.width, ui.height, "Testing", nil, nil)
	if err != nil {
		panic(err)
	}
	window.MakeContextCurrent()
	// Initialize Glow
	if err = gl.Init(); err != nil {
		panic(err)
	}
	version := gl.GoStr(gl.GetString(gl.VERSION))
	log.Println("OpenGL version", version)
	gl.Enable(gl.CULL_FACE)
	gl.DepthFunc(gl.LEQUAL)
	//注册关闭事件
	window.SetCloseCallback(ui.shutdown)
	return window
}

func (ui *Gui) getWH() (width, height int) {
	sw := glfw.GetPrimaryMonitor().GetVideoMode().Width
	sh := glfw.GetPrimaryMonitor().GetVideoMode().Height
	aspect := 16.0 / 9.0
	width = int(math.Min(float64(sw), float64(sh)*aspect)) - 80
	height = sh - 80
	return width, height
}

func Draw(vao uint32, window *glfw.Window, prog uint32, triangles []int) {
	gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

	// 使用这个程序
	gl.UseProgram(prog)

	// 绑定VAO，可能会感到奇怪，明明在 makeVao 中已经调用 gl.BindVertexArray(vao) 进行了绑定，为什么这里还要再绑定一次呢？
	// 因为绑定的操作是为了后续的操作服务的，并且有可能在中途又绑定了别的VAO，所以最好是在每次调用跟VAO有关的函数之前绑定一次。
	gl.BindVertexArray(vao)

	// 绘制的类型mode：
	//   1、gl.TRIANGLES：每三个顶点之间绘制三角形，之间不连接。
	//   2、gl.TRIANGLE_FAN：以V0,V1,V2;V0,V2,V3;V0,V3,V4，……的形式绘制三角形。
	//   3、gl.TRIANGLE_STRIP：以V0,V1,V2;V1,V2,V3;V2,V3,V4……的形式绘制三角形。
	// first：一般从第一个顶点开始
	// count：除以3为顶点个数
	gl.DrawArrays(gl.TRIANGLES, 0, int32(len(triangles)/3))

	glfw.PollEvents()

	window.SwapBuffers()
}

func makeVao(points []float32) uint32 {
	var vbo uint32

	// 在显卡中开辟一块空间，创建顶点缓存对象，个数为1，变量vbo会被赋予一个ID值。
	gl.GenBuffers(1, &vbo)

	// 将 vbo 赋值给 gl.ARRAY_BUFFER，要知道这个对象会被赋予不同的vbo，因此其值是变化的
	// 可选类型：GL_ARRAY_BUFFER, GL_ELEMENT_ARRAY_BUFFER, GL_PIXEL_PACK_BUFFER, GL_PIXEL_UNPACK_BUFFER
	gl.BindBuffer(gl.ARRAY_BUFFER, vbo)

	// 将内存中的数据传递到显卡中的gl.ARRAY_BUFFER对象上，其实是把数据传递到绑定在其上面的vbo对象上。
	// 4*len(points) 代表总的字节数，因为是32位的
	gl.BufferData(gl.ARRAY_BUFFER, 4*len(points), gl.Ptr(points), gl.STATIC_DRAW)

	var vao uint32
	// 创建顶点数组对象，个数为1，变量vao会被赋予一个ID值。
	gl.GenVertexArrays(1, &vao)

	// 后面的两个函数都是要操作具体的vao的，因此需要先将vao绑定到opengl上。
	// 解绑：gl.BindVertexArray(0)，opengl中很多的解绑操作都是传入0
	gl.BindVertexArray(vao)

	// 使vao去引用到gl.ARRAY_BUFFER上面的vbo，这一步完成之后vao就建立了对特定vbo的引用，后面即使gl.ARRAY_BUFFER 的值发生了变化也不影响vao的使用
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 0, nil)
	// 设置 vertex attribute 的状态enabled，默认是disabled，后面会有具体解释
	gl.EnableVertexAttribArray(0)

	return vao
}

func MakeVaoWithEbo(points []float32, indexs []uint32) uint32 {
	var vbo uint32
	gl.GenBuffers(1, &vbo)
	gl.BindBuffer(gl.ARRAY_BUFFER, vbo)
	gl.BufferData(gl.ARRAY_BUFFER, 4*len(points), gl.Ptr(points), gl.STATIC_DRAW)

	var ebo uint32
	gl.GenBuffers(1, &ebo)
	gl.BindBuffer(gl.ELEMENT_ARRAY_BUFFER, ebo)
	gl.BufferData(gl.ELEMENT_ARRAY_BUFFER, 4*len(indexs), gl.Ptr(indexs), gl.STATIC_DRAW)

	var vao uint32
	gl.GenVertexArrays(1, &vao)
	gl.BindVertexArray(vao)

	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 0, nil)

	return vao
}

func MakeTexture(filepath string) uint32 {
	var texture uint32
	gl.GenTextures(1, &texture)
	gl.BindTexture(gl.TEXTURE_2D, texture)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.REPEAT)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.REPEAT)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR)
	imgFile2, _ := os.Open(filepath)
	defer imgFile2.Close()
	img2, _, _ := image.Decode(imgFile2)
	rgba2 := image.NewRGBA(img2.Bounds())
	draw.Draw(rgba2, rgba2.Bounds(), img2, image.Point{0, 0}, draw.Src)
	gl.TexImage2D(gl.TEXTURE_2D, 0, gl.RGBA, int32(rgba2.Rect.Size().X), int32(rgba2.Rect.Size().Y), 0, gl.RGBA, gl.UNSIGNED_BYTE, gl.Ptr(rgba2.Pix))
	gl.GenerateMipmap(gl.TEXTURE_2D)

	return texture
}

func compileShader(source string, shaderType uint32) (uint32, error) {
	shader := gl.CreateShader(shaderType)

	csources, free := gl.Strs(source)
	gl.ShaderSource(shader, 1, csources, nil)
	free()
	gl.CompileShader(shader)

	var status int32
	gl.GetShaderiv(shader, gl.COMPILE_STATUS, &status)
	if status == gl.FALSE {
		var logLength int32
		gl.GetShaderiv(shader, gl.INFO_LOG_LENGTH, &logLength)

		log := strings.Repeat("\x00", int(logLength+1))
		gl.GetShaderInfoLog(shader, logLength, nil, gl.Str(log))

		return 0, fmt.Errorf("failed to compile %v: %v", source, log)
	}

	return shader, nil
}

func newProgram(vertexShaderSource, fragmentShaderSource string) (uint32, error) {
	vertexShader, err := compileShader(vertexShaderSource, gl.VERTEX_SHADER)
	if err != nil {
		return 0, err
	}

	fragmentShader, err := compileShader(fragmentShaderSource, gl.FRAGMENT_SHADER)
	if err != nil {
		return 0, err
	}

	program := gl.CreateProgram()

	gl.AttachShader(program, vertexShader)
	gl.AttachShader(program, fragmentShader)
	gl.LinkProgram(program)

	var status int32
	gl.GetProgramiv(program, gl.LINK_STATUS, &status)
	if status == gl.FALSE {
		var logLength int32
		gl.GetProgramiv(program, gl.INFO_LOG_LENGTH, &logLength)

		log := strings.Repeat("\x00", int(logLength+1))
		gl.GetProgramInfoLog(program, logLength, nil, gl.Str(log))

		return 0, fmt.Errorf("failed to link program: %v", log)
	}

	gl.DeleteShader(vertexShader)
	gl.DeleteShader(fragmentShader)

	return program, nil
}

var vertexShader = `
#version 330

uniform mat4 projection;
uniform mat4 camera;
uniform mat4 model;

in vec3 vert;
in vec2 vertTexCoord;

out vec2 fragTexCoord;

void main() {
    fragTexCoord = vertTexCoord;
    gl_Position = projection * camera * model * vec4(vert, 1);
}
` + "\x00"

var fragmentShader = `
#version 330

uniform sampler2D tex;

in vec2 fragTexCoord;

out vec4 outputColor;

void main() {
    outputColor = texture(tex, fragTexCoord);
}
` + "\x00"
