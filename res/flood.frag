#version 430

in vec2 fragTexCoord;
in vec4 fragColor;

out vec4 outColor;

uniform layout(binding=3, rgba8) image2D colors_0;
uniform layout(binding=4, rgba32f) image2D masVelFix_0;
uniform layout(binding=5, rgba8) image2D colors_1;
uniform layout(binding=6, rgba32f) image2D masVelFix_1;

uniform vec2 RESOLUTION;
uniform int IMG_IDX;

void main() {
    ivec2 xy = ivec2(int(fragTexCoord.x * RESOLUTION.x), int(fragTexCoord.y * RESOLUTION.y));
    outColor = vec4(imageLoad(colors_0, xy).rgb, 1.0);
    imageStore(colors_0, xy, vec4(1.0, 0., 0., 1.0));
}