#version 150 core

in float lightIntensity;

out vec4 FragColor;

uniform vec4 objectColor;
uniform bool isGrid;
uniform bool GLOW;

void main() {
    if (isGrid) {
        FragColor = objectColor;
    } else if (GLOW) {
        FragColor = vec4(objectColor.rgb * 100000.0, objectColor.a);
    } else {
        float fade = smoothstep(0.0, 10.0, lightIntensity * 10.0);
        FragColor = vec4(objectColor.rgb * fade, objectColor.a);
    }
}
