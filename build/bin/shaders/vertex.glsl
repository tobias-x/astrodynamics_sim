#version 150 core

in vec3 aPos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out float lightIntensity;

void main() {
    gl_Position = projection * view * model * vec4(aPos, 1.0);
    vec3 worldPos = (model * vec4(aPos, 1.0)).xyz;
    vec3 normal = normalize(aPos);
    vec3 dirToCenter = normalize(-worldPos);
    lightIntensity = max(dot(normal, dirToCenter), 0.15);
}
