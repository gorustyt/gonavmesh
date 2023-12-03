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
	//axis := mgl32.Vec3{0.5, 1, 0}
	//view := mgl32.Translate3D(0, 0, -3)
	window.MakeContextCurrent()
	// Initialize Glow
	if err := gl.Init(); err != nil {
		panic(err)
	}
	vers := []float32{
		-0.5, -0.5, -0.5, 0.0, 0.0,
		0.5, -0.5, -0.5, 1.0, 0.0,
		0.5, 0.5, -0.5, 1.0, 1.0,
		0.5, 0.5, -0.5, 1.0, 1.0,
		-0.5, 0.5, -0.5, 0.0, 1.0,
		-0.5, -0.5, -0.5, 0.0, 0.0,

		-0.5, -0.5, 0.5, 0.0, 0.0,
		0.5, -0.5, 0.5, 1.0, 0.0,
		0.5, 0.5, 0.5, 1.0, 1.0,
		0.5, 0.5, 0.5, 1.0, 1.0,
		-0.5, 0.5, 0.5, 0.0, 1.0,
		-0.5, -0.5, 0.5, 0.0, 0.0,

		-0.5, 0.5, 0.5, 1.0, 0.0,
		-0.5, 0.5, -0.5, 1.0, 1.0,
		-0.5, -0.5, -0.5, 0.0, 1.0,
		-0.5, -0.5, -0.5, 0.0, 1.0,
		-0.5, -0.5, 0.5, 0.0, 0.0,
		-0.5, 0.5, 0.5, 1.0, 0.0,

		0.5, 0.5, 0.5, 1.0, 0.0,
		0.5, 0.5, -0.5, 1.0, 1.0,
		0.5, -0.5, -0.5, 0.0, 1.0,
		0.5, -0.5, -0.5, 0.0, 1.0,
		0.5, -0.5, 0.5, 0.0, 0.0,
		0.5, 0.5, 0.5, 1.0, 0.0,

		-0.5, -0.5, -0.5, 0.0, 1.0,
		0.5, -0.5, -0.5, 1.0, 1.0,
		0.5, -0.5, 0.5, 1.0, 0.0,
		0.5, -0.5, 0.5, 1.0, 0.0,
		-0.5, -0.5, 0.5, 0.0, 0.0,
		-0.5, -0.5, -0.5, 0.0, 1.0,

		-0.5, 0.5, -0.5, 0.0, 1.0,
		0.5, 0.5, -0.5, 1.0, 1.0,
		0.5, 0.5, 0.5, 1.0, 0.0,
		0.5, 0.5, 0.5, 1.0, 0.0,
		-0.5, 0.5, 0.5, 0.0, 0.0,
		-0.5, 0.5, -0.5, 0.0, 1.0,
	}
	cubePositions := []mgl32.Mat4{
		mgl32.Translate3D(0.0, 0.0, -0.0),
		mgl32.Translate3D(2.0, 5.0, -15.0),
		mgl32.Translate3D(-1.5, -2.2, -2.5),
		mgl32.Translate3D(-3.8, -2.0, -12.3),
		mgl32.Translate3D(2.4, -0.4, -3.5),
		mgl32.Translate3D(-1.7, 3.0, -7.5),
		mgl32.Translate3D(1.3, -2.0, -2.5),
		mgl32.Translate3D(1.5, 2.0, -2.5),
		mgl32.Translate3D(1.5, 0.2, -1.5),
		mgl32.Translate3D(-1.3, 1.0, -1.5),
	}
	t := makeTexture("./demo1/assets/awesomeface.png", gl.TEXTURE0)
	t1 := makeTexture("./demo1/assets/container.jpg", gl.TEXTURE1)
	vao := canvas.MakeVao(vers)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 5*4, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointer(1, 2, gl.FLOAT, false, 5*4, gl.PtrOffset(3*4))
	gl.EnableVertexAttribArray(1)
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
	fov := float32(45.0)
	lastFrame := glfw.GetTime()
	cameraFront := mgl32.Vec3{0, 0, -1}
	cameraPos := mgl32.Vec3{0, 0, 3}
	cameraUp := mgl32.Vec3{0, 1, 0}
	window.SetScrollCallback(func(w *glfw.Window, xoff float64, yoffset float64) {
		if fov >= 1.0 && fov <= 45.0 {
			fov -= float32(yoffset)
		}
		if fov <= 1.0 {
			fov = 1.0
		}
		if fov >= 45.0 {
			fov = 45.0
		}

	})
	firstMouse := true
	lastX := float32(windowWidth / 2)
	lastY := float32(windowHeight / 2)
	yaw := float32(-90.0) // yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
	pitch := float32(0.0)
	window.SetCursorPosCallback(func(w *glfw.Window, xp float64, yp float64) {

		xpos := float32(xp)
		ypos := float32(yp)
		if firstMouse {
			lastX = xpos
			lastY = ypos
			firstMouse = false
		}
		xoffset := xpos - lastX
		yoffset := lastY - ypos
		lastX = xpos
		lastY = ypos

		sensitivity := float32(0.05)
		xoffset *= sensitivity
		yoffset *= sensitivity

		yaw += xoffset
		pitch += yoffset

		if pitch > 89.0 {
			pitch = 89.0
		}

		if pitch < -89.0 {
			pitch = -89.0
		}
		if window.GetMouseButton(glfw.MouseButtonLeft) != glfw.Press &&
			window.GetMouseButton(glfw.MouseButtonRight) != glfw.Press {
			return
		}
		front := mgl32.Vec3{
			float32(math.Cos(float64(mgl32.DegToRad(yaw))) * math.Cos(float64(mgl32.DegToRad(pitch)))),
			float32(math.Sin(float64(mgl32.DegToRad(pitch)))),
			float32(math.Sin(float64(mgl32.DegToRad(yaw))) * math.Cos(float64(mgl32.DegToRad(pitch)))),
		}
		cameraFront = front.Normalize()
	})
	for !window.ShouldClose() {
		gl.ClearColor(0.2, 0.3, 0.3, 1.0)
		gl.UseProgram(prg)
		speed := float32(glfw.GetTime()-lastFrame) * 2.5
		lastFrame = glfw.GetTime()
		if window.GetKey(glfw.KeySpace) == glfw.Press {

		} else if window.GetKey(glfw.KeyUp) == glfw.Press {
			fov += 00.01
			tt = float32(math.Min(1.0, float64(tt)+0.001))
		} else if window.GetKey(glfw.KeyDown) == glfw.Press {
			fov -= 00.01
			tt = float32(math.Max(0.1, float64(tt)-0.001))
		} else if window.GetKey(glfw.KeyW) == glfw.Press {
			cameraPos = cameraPos.Add(cameraFront.Mul(speed))
		} else if window.GetKey(glfw.KeyA) == glfw.Press {
			cameraPos = cameraPos.Sub(cameraFront.Cross(cameraUp).Normalize().Mul(speed))
		} else if window.GetKey(glfw.KeyS) == glfw.Press {
			cameraPos = cameraPos.Sub(cameraFront.Mul(speed))
		} else if window.GetKey(glfw.KeyD) == glfw.Press {
			cameraPos = cameraPos.Add(cameraFront.Cross(cameraUp).Normalize().Mul(speed))
		}

		//radius := 6.0
		//cx := float32(math.Sin(glfw.GetTime()) * radius)
		//	cz := float32(math.Cos(glfw.GetTime()) * radius)

		projection := mgl32.Perspective(mgl32.DegToRad(fov), float32(windowWidth)/float32(windowHeight), 0.1, 100)
		gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("projection\x00")),
			1, false, &projection[0])
		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("t\x00")), tt)
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		gl.ActiveTexture(t)
		gl.BindTexture(gl.TEXTURE_2D, t)
		gl.ActiveTexture(t1)
		gl.BindTexture(gl.TEXTURE_2D, t1)
		gl.BindVertexArray(vao)
		for i := range cubePositions {
			model := cubePositions[i]
			if i%3 == 0 {
				model = model.Mul4(mgl32.HomogRotate3D(mgl32.DegToRad(float32(glfw.GetTime())*float32(i+1)), mgl32.Vec3{1.0, 0.3, 0.5}))
			}
			view := mgl32.LookAtV(cameraPos, cameraPos.Add(cameraFront), cameraUp)
			gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("view\x00")),
				1, false, &view[0])
			gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("model\x00")),
				1, false, &model[0])
			gl.DrawArrays(gl.TRIANGLES, 0, 36)
		}

		// Maintenance
		window.SwapBuffers()
		glfw.PollEvents()
	}
}

var vertexShader = `
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec2 texCoord;

out vec2 TexCoord;
uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
void main()
{
	gl_Position = projection*view*model*vec4(position, 1.0f);
	TexCoord = vec2(texCoord.x, 1.0 - texCoord.y);
}
` + "\x00"
var fragmentShader = `
#version 330 core
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
