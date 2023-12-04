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
	lightVs := []float32{
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
	objVs := []float32{
		// positions          // normals           // texture coords
		-0.5, -0.5, -0.5, 0.0, 0.0, -1.0, 0.0, 0.0,
		0.5, -0.5, -0.5, 0.0, 0.0, -1.0, 1.0, 0.0,
		0.5, 0.5, -0.5, 0.0, 0.0, -1.0, 1.0, 1.0,
		0.5, 0.5, -0.5, 0.0, 0.0, -1.0, 1.0, 1.0,
		-0.5, 0.5, -0.5, 0.0, 0.0, -1.0, 0.0, 1.0,
		-0.5, -0.5, -0.5, 0.0, 0.0, -1.0, 0.0, 0.0,

		-0.5, -0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0,
		0.5, -0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 0.0,
		0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 1.0,
		0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 1.0, 1.0,
		-0.5, 0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 1.0,
		-0.5, -0.5, 0.5, 0.0, 0.0, 1.0, 0.0, 0.0,

		-0.5, 0.5, 0.5, -1.0, 0.0, 0.0, 1.0, 0.0,
		-0.5, 0.5, -0.5, -1.0, 0.0, 0.0, 1.0, 1.0,
		-0.5, -0.5, -0.5, -1.0, 0.0, 0.0, 0.0, 1.0,
		-0.5, -0.5, -0.5, -1.0, 0.0, 0.0, 0.0, 1.0,
		-0.5, -0.5, 0.5, -1.0, 0.0, 0.0, 0.0, 0.0,
		-0.5, 0.5, 0.5, -1.0, 0.0, 0.0, 1.0, 0.0,

		0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0,
		0.5, 0.5, -0.5, 1.0, 0.0, 0.0, 1.0, 1.0,
		0.5, -0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 1.0,
		0.5, -0.5, -0.5, 1.0, 0.0, 0.0, 0.0, 1.0,
		0.5, -0.5, 0.5, 1.0, 0.0, 0.0, 0.0, 0.0,
		0.5, 0.5, 0.5, 1.0, 0.0, 0.0, 1.0, 0.0,

		-0.5, -0.5, -0.5, 0.0, -1.0, 0.0, 0.0, 1.0,
		0.5, -0.5, -0.5, 0.0, -1.0, 0.0, 1.0, 1.0,
		0.5, -0.5, 0.5, 0.0, -1.0, 0.0, 1.0, 0.0,
		0.5, -0.5, 0.5, 0.0, -1.0, 0.0, 1.0, 0.0,
		-0.5, -0.5, 0.5, 0.0, -1.0, 0.0, 0.0, 0.0,
		-0.5, -0.5, -0.5, 0.0, -1.0, 0.0, 0.0, 1.0,

		-0.5, 0.5, -0.5, 0.0, 1.0, 0.0, 0.0, 1.0,
		0.5, 0.5, -0.5, 0.0, 1.0, 0.0, 1.0, 1.0,
		0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 1.0, 0.0,
		0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 1.0, 0.0,
		-0.5, 0.5, 0.5, 0.0, 1.0, 0.0, 0.0, 0.0,
		-0.5, 0.5, -0.5, 0.0, 1.0, 0.0, 0.0, 1.0,
	}
	//lightColor := mgl32.Vec3{1.0, 1.0, 1.0}
	origin := mgl32.Vec3{3.2, 3.0, 2.0}
	lightPos := origin
	//objectColor := mgl32.Vec3{1.0, 0.5, 0.31}
	lightVao := canvas.MakeVao(lightVs)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 3*4, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(0)
	objVao := canvas.MakeVao(objVs)
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 8*4, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(0)
	gl.VertexAttribPointer(1, 3, gl.FLOAT, false, 8*4, gl.PtrOffset(3*4))
	gl.EnableVertexAttribArray(1)
	gl.VertexAttribPointer(2, 2, gl.FLOAT, false, 8*4, gl.PtrOffset(6*4))
	gl.EnableVertexAttribArray(2)
	firstMouse := true
	lastX := float32(windowWidth / 2)
	lastY := float32(windowHeight / 2)
	yaw := float32(-90.0) // yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
	pitch := float32(0.0)
	tt := float32(0.2)
	fov := float32(45.0)
	lastFrame := glfw.GetTime()
	cameraFront := mgl32.Vec3{0, 0, -1}
	cameraPos := mgl32.Vec3{0, 0, 3}
	cameraUp := mgl32.Vec3{0, 1, 0}
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
	diffuseMap := makeTexture("demo1/assets/container2.png", 0)
	diffuseMap1 := makeTexture("demo1/assets/container2_specular.png", 1)
	gl.Enable(gl.DEPTH_TEST)
	gl.DepthFunc(gl.LESS)
	cubePositions := []mgl32.Vec3{
		mgl32.Vec3{0.0, 0.0, 0.0},
		mgl32.Vec3{2.0, 5.0, -15.0},
		mgl32.Vec3{-1.5, -2.2, -2.5},
		mgl32.Vec3{-3.8, -2.0, -12.3},
		mgl32.Vec3{2.4, -0.4, -3.5},
		mgl32.Vec3{-1.7, 3.0, -7.5},
		mgl32.Vec3{1.3, -2.0, -2.5},
		mgl32.Vec3{1.5, 2.0, -2.5},
		mgl32.Vec3{1.5, 0.2, -1.5},
		mgl32.Vec3{-1.3, 1.0, -1.5},
	}
	for !window.ShouldClose() {
		gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
		gl.ClearColor(0.2, 0.3, 0.3, 1.0)

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
		view := mgl32.LookAtV(cameraPos, cameraPos.Add(cameraFront), cameraUp)
		projection := mgl32.Perspective(mgl32.DegToRad(45), float32(windowWidth)/float32(windowHeight), 0.1, 100)
		gl.UseProgram(prg)
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("lightColor\x00")), 1.0, 1.0, 1.0)
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("objectColor\x00")), 1.0, 0.5, 0.31)

		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("material.ambient\x00")), 1.0, 0.5, 0.31)
		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("material.diffuse\x00")), 0)
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("material.specular\x00")), 0.5, 0.5, 0.5)
		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("material.shininess\x00")), 32)

		diffuseColor := mgl32.Vec3{0.5, 0.5, 0.5} // 降低影响
		ambientColor := mgl32.Vec3{0.2, 0.2, 0.2} // 很低的影响
		//gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("light.direction\x00")), -0.2, -1, -0.3)
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("light.ambient\x00")), ambientColor[0], ambientColor[1], ambientColor[2])
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("light.diffuse\x00")), diffuseColor[0], diffuseColor[1], diffuseColor[2])
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("light.specular\x00")), 1.0, 1.0, 1.0)
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("light.position\x00")), 1.1, 1, 2)

		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("light.constant\x00")), 1)
		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("light.linear\x00")), 0.09)
		gl.Uniform1f(gl.GetUniformLocation(prg, gl.Str("light.quadratic\x00")), 0.032)

		//gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("lightPos\x00")), lightPos[0], lightPos[1], lightPos[2])
		gl.Uniform3f(gl.GetUniformLocation(prg, gl.Str("viewPos\x00")), cameraPos[0], cameraPos[1], cameraPos[2])
		gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("projection\x00")),
			1, false, &projection[0])
		//modle := mgl32.Translate3D(0, 0, 0.0).Mul4(mgl32.HomogRotate3D(float32(45), mgl32.Vec3{1.0, 0.3, 0.5}.Normalize()))
		//model := mgl32.Scale3D(0.8, 0.8, 0.8).Mul4(mgl32.HomogRotate3D(mgl32.DegToRad(45), mgl32.Vec3{0, 1, 1}))

		gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("view\x00")),
			1, false, &view[0])
		gl.ActiveTexture(gl.TEXTURE0)
		gl.BindTexture(gl.TEXTURE_2D, diffuseMap)
		gl.ActiveTexture(gl.TEXTURE1)
		gl.BindTexture(gl.TEXTURE_2D, diffuseMap1)
		gl.BindVertexArray(objVao)

		for i, v := range cubePositions {
			angle := 20 * i
			model := mgl32.Translate3D(v[0], v[1], v[2]).Mul4(mgl32.HomogRotate3D(mgl32.DegToRad(float32(angle)), mgl32.Vec3{1, 0.3, 0.5}))
			gl.UniformMatrix4fv(gl.GetUniformLocation(prg, gl.Str("model\x00")),
				1, false, &model[0])
			gl.DrawArrays(gl.TRIANGLES, 0, 36)
		}
		gl.BindVertexArray(0)

		gl.UseProgram(lightPrg)
		gl.UniformMatrix4fv(gl.GetUniformLocation(lightPrg, gl.Str("projection\x00")),
			1, false, &projection[0])
		lightModel := mgl32.Scale3D(0.2, 0.2, 0.2)
		lightModel = lightModel.Mul4(mgl32.Translate3D(lightPos[0], lightPos[1], lightPos[2]))
		gl.UniformMatrix4fv(gl.GetUniformLocation(lightPrg, gl.Str("model\x00")),
			1, false, &lightModel[0])
		gl.UniformMatrix4fv(gl.GetUniformLocation(lightPrg, gl.Str("view\x00")),
			1, false, &view[0])
		gl.BindVertexArray(lightVao)
		gl.DrawArrays(gl.TRIANGLES, 0, 36)
		gl.BindVertexArray(0)
		// Maintenance
		window.SwapBuffers()
		glfw.PollEvents()
	}
}

var vertexShader = `
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 innormal;
layout (location = 2) in vec2 aTexCoords;
out vec2 TexCoords;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
out vec3 norm;
out vec3 FragPos;  
void main()
{
	gl_Position = projection*view*model*vec4(position, 1.0f);
	norm=mat3(transpose(inverse(model))) * normalize(innormal);
	FragPos = vec3(model * vec4(position, 1.0));
	TexCoords = aTexCoords;
}
` + "\x00"
var fragmentShader = `
#version 330 core
in vec3 norm;
in vec3 FragPos;

out vec4 FragColor;

struct Material {
    sampler2D diffuse;
   sampler2D specular;
    float     shininess;
}; 
in vec2 TexCoords;
uniform Material material;
struct Light {
    vec3 position; // 使用定向光就不再需要了
    //vec3 direction;
    vec3 ambient;
    vec3 diffuse;
    vec3 specular;

 float constant;
    float linear;
    float quadratic;
};

uniform Light light;
uniform vec3 objectColor;
uniform vec3 lightColor;
uniform vec3 lightPos;
uniform vec3 viewPos;
void main()
{
	float distance    = length(light.position - FragPos);
float attenuation = 1.0 / (light.constant + light.linear * distance + 
                light.quadratic * (distance * distance));
    //vec3 lightDir = normalize(-light.direction);
	vec3 lightDir =normalize(light.position - FragPos);
    float diff = max(dot(norm, lightDir), 0.0);

vec3 viewDir = normalize(viewPos - FragPos);
    vec3 reflectDir = reflect(-lightDir, norm);  
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), material.shininess);
vec3 ambient  = light.ambient  * vec3(texture(material.diffuse, TexCoords));
vec3 diffuse  = light.diffuse  * diff * vec3(texture(material.diffuse, TexCoords));  
vec3 specular = light.specular * spec * vec3(texture(material.specular, TexCoords));
ambient  *= attenuation; 
diffuse  *= attenuation;
specular *= attenuation;
FragColor = vec4(ambient + diffuse + specular, 1.0);
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
