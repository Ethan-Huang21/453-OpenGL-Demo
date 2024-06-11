#version 330 core
layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 color;

out vec3 fragColor;

uniform mat4 viewMatrix;

void main() {
	fragColor = color;
	gl_Position = viewMatrix * vec4(pos, 1.0);
}