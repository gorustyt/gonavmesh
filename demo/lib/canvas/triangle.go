package canvas

import "github.com/go-gl/gl/v4.1-core/gl"

const (
	GL_FLOAT32_SIZE = 4
)

type Triangle struct {
	program  uint32
	vao, vbo uint32
}

func NewTriangle() *Triangle {
	t := &Triangle{}
	err := t.initGL()
	if err != nil {
		panic(err)
	}
	return t
}

func (ht *Triangle) createBuffers(vertices []float32) (uint32, uint32) {
	var vao, vbo uint32
	gl.GenVertexArrays(1, &vao)
	gl.GenBuffers(1, &vbo)

	gl.BindVertexArray(vao)

	gl.BindBuffer(gl.ARRAY_BUFFER, vbo)
	gl.BufferData(gl.ARRAY_BUFFER, len(vertices)*GL_FLOAT32_SIZE, gl.Ptr(vertices), gl.STATIC_DRAW)

	//vertAttrib := uint32(gl.GetAttribLocation(ht.program, gl.Str("position\x00")))
	// here we can skip computing the vertAttrib value and use 0 since our shader declares layout = 0 for
	// the uniform
	gl.VertexAttribPointer(0, 3, gl.FLOAT, false, 3*GL_FLOAT32_SIZE, gl.PtrOffset(0))
	gl.EnableVertexAttribArray(0)

	gl.BindVertexArray(0)
	return vao, vbo
}

func (ht *Triangle) initGL() error {
	var vertexShader = `
	#version 330 core
	layout (location = 0) in vec3 position;
	void main() {
	  gl_Position = vec4(position.x, position.y, position.z, 1.0);
	}` + "\x00"

	var fragShader = `
	#version 330 core
	out vec4 color;
	void main() {
	  color = vec4(1.0f, 0.0f, 0.0f, 0.0f);
	}` + "\x00"

	var err error
	ht.program, err = NewProgram(vertexShader, fragShader)
	if err != nil {
		return err
	}
	return nil
}

func (ht *Triangle) Draw(points []float32) {
	ht.vao = MakeVao(points)
	//gl.Clear(gl.COLOR_BUFFER_BIT | gl.DEPTH_BUFFER_BIT)
	//gl.ClearColor(ht.Color32.R, ht.Color32.G, ht.Color32.B, ht.Color32.A)
	// Draw our first triangle
	gl.UseProgram(ht.program)
	gl.BindVertexArray(ht.vao)
	gl.DrawArrays(gl.TRIANGLES, 0, 3)
	gl.BindVertexArray(0)
}

func (ht *Triangle) Close() {
	gl.DeleteVertexArrays(1, &ht.vao)
	//gl.DeleteBuffers(1, &ht.vbo)
	gl.DeleteProgram(ht.program)
}
