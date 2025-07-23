#version 330 core
out vec4 FragColor;

in vec3 fNormal_world;
in vec3 fpos_world;
in vec4 fColor;

uniform vec3 lightPos;   // Passed from BaseGLWidget
uniform vec3 viewPos;    // Passed from BaseGLWidget
uniform vec3 lightColor; // Passed from BaseGLWidget

void main() {
    // Ambient
    float ambientStrength = 0.3;
    vec3 ambient = ambientStrength * lightColor;

    // Diffuse
    vec3 norm = normalize(fNormal_world);
    vec3 lightDir = normalize(lightPos - fpos_world);
    float diff = max(dot(norm, lightDir), 0.0);
    vec3 diffuse = diff * lightColor;

    // Specular
    float specularStrength = 0.5;
    vec3 viewDir = normalize(viewPos - fpos_world);
    vec3 reflectDir = reflect(-lightDir, norm);
    float spec = pow(max(dot(viewDir, reflectDir), 0.0), 32);
    vec3 specular = specularStrength * spec * lightColor;

    vec3 result = (ambient + diffuse + specular) * fColor.rgb;
    FragColor = vec4(result, fColor.a);
}