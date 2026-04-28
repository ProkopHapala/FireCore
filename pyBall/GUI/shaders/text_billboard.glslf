#version 330 core

in vec2 v_texCoord;
uniform sampler2D fontAtlas;
uniform vec4 textColor;        // New uniform for text color

out vec4 FragColor;

void main() {
    vec4 tex = texture(fontAtlas, v_texCoord);
    FragColor = vec4(textColor.rgb, tex.a+0.0);
    // FragColor = vec4(textColor.rgb, tex.a+0.1); // debug UV mapping
    // FragColor = vec4(1.,0.,0.,1.);              // debug billboard placement
}
