#version 330 core

in vec2 v_texCoord;
uniform sampler2D fontAtlas;
uniform vec4 textColor; // New uniform for text color

out vec4 FragColor;

void main() {
    // Sample the font atlas texture
    float alpha = texture(fontAtlas, v_texCoord).r; // Assuming a single-channel (red) texture for the font atlas
    //if (alpha < 0.1) discard;
    // Use the alpha channel for transparency, and apply the uniform text color
    //FragColor = vec4(textColor.rgb, alpha * textColor.a);
    FragColor = vec4(1.,0.,1.,1.);
}
