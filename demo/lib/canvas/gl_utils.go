package canvas

import (
	"fmt"
	"github.com/go-gl/gl/v4.1-core/gl"
	"github.com/go-gl/glfw/v3.3/glfw"
	"image"
	"image/draw"
	"os"
	"strings"
)

func NewProgram(vertexShaderSource, fragmentShaderSource string) (uint32, error) {
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

func NewTexture(file string) (uint32, error) {
	imgFile, err := os.Open(file)
	if err != nil {
		return 0, fmt.Errorf("texture %q not found on disk: %v", file, err)
	}
	img, _, err := image.Decode(imgFile)
	if err != nil {
		return 0, err
	}

	rgba := image.NewRGBA(img.Bounds())
	if rgba.Stride != rgba.Rect.Size().X*4 {
		return 0, fmt.Errorf("unsupported stride")
	}
	draw.Draw(rgba, rgba.Bounds(), img, image.Point{0, 0}, draw.Src)

	var texture uint32
	gl.GenTextures(1, &texture)
	gl.ActiveTexture(gl.TEXTURE0)
	gl.BindTexture(gl.TEXTURE_2D, texture)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.LINEAR)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE)
	gl.TexParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE)
	gl.TexImage2D(
		gl.TEXTURE_2D,
		0,
		gl.RGBA,
		int32(rgba.Rect.Size().X),
		int32(rgba.Rect.Size().Y),
		0,
		gl.RGBA,
		gl.UNSIGNED_BYTE,
		gl.Ptr(rgba.Pix))

	return texture, nil
}

func MakeVao(points []float32) uint32 {
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
	return vao
}

func Draw(vao uint32, window *glfw.Window, prog uint32, points []float32) {
	gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)

	// 使用这个程序
	gl.UseProgram(prog)

	// 绑定VAO，可鞥会感到奇怪，明明在 makeVao 中已经调用 gl.BindVertexArray(vao) 进行了绑定，为什么这里还要再绑定一次呢？
	// 因为绑定的操作是为了后续的操作服务的，并且有可能在中途又绑定了别的VAO，所以最好是在每次调用跟VAO有关的函数之前绑定一次。
	gl.BindVertexArray(vao)

	// 绘制的类型mode：
	//   1、gl.TRIANGLES：每三个顶点之间绘制三角形，之间不连接。
	//   2、gl.TRIANGLE_FAN：以V0,V1,V2;V0,V2,V3;V0,V3,V4，……的形式绘制三角形。
	//   3、gl.TRIANGLE_STRIP：以V0,V1,V2;V1,V2,V3;V2,V3,V4……的形式绘制三角形。
	// first：一般从第一个顶点开始
	// count：除以3为顶点个数
	gl.DrawArrays(gl.TRIANGLES, 0, int32(len(points)/3))
}
