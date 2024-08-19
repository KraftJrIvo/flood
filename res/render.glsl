#version 430
in vec2 fragTexCoord;
out vec4 finalColor;

layout(std430, binding = 1) readonly buffer colorsLayout
{
    float colorsBuffer[];
};

uniform vec2 RESOLUTION;

void main()
{
    ivec2 coords = ivec2(fragTexCoord * RESOLUTION);
    float r = colorsBuffer[coords.y * uvec2(RESOLUTION).x * 3 + coords.x * 3 + 0];
    float g = colorsBuffer[coords.y * uvec2(RESOLUTION).x * 3 + coords.x * 3 + 1];
    float b = colorsBuffer[coords.y * uvec2(RESOLUTION).x * 3 + coords.x * 3 + 2];
    finalColor = vec4(r, g, b, 1.0);
}