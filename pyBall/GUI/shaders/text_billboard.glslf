#version 330 core

in vec2 v_texCoord;
uniform sampler2D fontAtlas;
uniform vec4 textColor; // New uniform for text color

out vec4 FragColor;

void main() {
    // DEBUG: output the raw texture color so we can see if the atlas is bound
    vec4 tex = texture(fontAtlas, v_texCoord);
    //FragColor = tex;
    FragColor = vec4(textColor.rgb*tex.a, tex.a);
    //FragColor = vec4(1.,0.,0.,1.);
}
