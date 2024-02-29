package mesh

var (
	objVertexShader = `
#version 330 core
layout (location = 0) in vec3 position;
layout (location = 1) in vec3 color;
layout (location = 2) in vec2 texCoord;

out vec2 TexCoord;
out vec3 Vcolor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

void main()
{
    gl_Position = projection * view * model * vec4(position, 1.0f);
    TexCoord = vec2(texCoord.x, 1.0 - texCoord.y);
	Vcolor=color;
}`

	objFragShader = `#version 330 core
in vec2 TexCoord;
in vec3 Vcolor;

out vec4 FragColor;

uniform sampler2D texture1;
void main()
{
	FragColor = texture(texture1, TexCoord);
	FragColor*=Vcolor;
}`
)
