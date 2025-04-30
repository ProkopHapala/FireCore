
#ifndef SHADER_H
#define SHADER_H

#include "GLES.h"
#include "Mat3.h"
#include "Vec2.h"
#include "Vec3.h"
#include "Mat4.h"
#include "quaternion.h"
#include <string>

#define SHADER_ATTRIB_POSITION 0
#define SHADER_ATTRIB_NORMAL 1
#define SHADER_ATTRIB_COLOR 2
#define SHADER_ATTRIB_UV 3

template<unsigned int attrib_flags>
class Shader {
private:
    GLuint programId        = 0;
    GLuint vertexShaderId   = 0;
    GLuint fragmentShaderId = 0;

    GLuint uMVPloc = -1;
    GLuint uColorloc = -1;

    char* __vertexShaderSource   = nullptr;
    char* __fragmentShaderSource = nullptr;

    struct uniform{
        enum {f1, f2, f3, f4, i1, i2, i3, i4, ui1, ui2, ui3, ui4, m2, m3, m4} type;
        union data{
            GLfloat f1;
            Vec2T<GLfloat> f2;
            Vec3T<GLfloat> f3;
            Vec4T<GLfloat> f4;

            GLint i1;
            Vec2T<GLint> i2;
            Vec3T<GLint> i3;
            Vec4T<GLint> i4;

            GLuint ui1;
            Vec2T<GLuint> ui2;
            Vec3T<GLuint> ui3;
            Vec4T<GLuint> ui4;

            //Mat2T<GLfloat> m2;
            Mat3T<GLfloat> m3;
            Mat4T<GLfloat> m4;
        }data;

        bool operator==(const uniform& other) const {
            if (type != other.type) return false;
            switch (type) {
                case f1: return data.f1 == other.data.f1;
                case f2: return data.f2 == other.data.f2;
                case f3: return data.f3 == other.data.f3;
                case f4: return data.f4 == other.data.f4;
                case i1: return data.i1 == other.data.i1;
                case i2: return data.i2 == other.data.i2;
                case i3: return data.i3 == other.data.i3;
                case i4: return data.i4 == other.data.i4;
                case ui1: return data.ui1 == other.data.ui1;
                case ui2: return data.ui2 == other.data.ui2;
                case ui3: return data.ui3 == other.data.ui3;
                case ui4: return data.ui4 == other.data.ui4;
                case m3: return data.m3 == other.data.m3;
                case m4: return data.m4 == other.data.m4;
                default: return false;
            }
        };
    };

    unsigned int update_uniforms = 0;
    std::vector<uniform> uniforms;

    inline void setUniform(const char* name, uniform value){
        ensure_handle();
        GLuint loc = getUniformLocation(name);
        if (loc == -1) return;

        if (loc >= uniforms.size()){ uniforms.resize(loc + 1); }
        
        update_uniforms |= (1 << loc) * (uniforms[loc] == value); // TODO: check that loc is not > 32
        uniforms[loc] = value;
    }

public:
    Shader(const char* vertexShaderSource, const char* fragmentShaderSource);
    ~Shader();
    
    inline GLuint getProgramId() const { return programId; }

    void ensure_handle();
    void use();
    
    // Uniform setters
    GLint getUniformLocation(const char* name) {
        GLint loc = glGetUniformLocation(programId, name);
        if (loc == -1) {
            printf("Uniform '%s' not found\n", name);
            return -1;
        }
        return loc;
    }

    inline void setuMVPMatrix(Mat4f mat){ ensure_handle(); glUniformMatrix4fv(uMVPloc, 1, GL_FALSE, mat.array); }
    inline void setuColor(Vec3f color){ ensure_handle(); glUniform3f(uColorloc, color.x, color.y, color.z); }


    inline void setUniform1f(const char* name, GLfloat value)       { setUniform(name, (uniform){.type=uniform::f1, .data={.f1=value}}); }
    inline void setUniform2f(const char* name, Vec2T<GLfloat> value){ setUniform(name, (uniform){.type=uniform::f2, .data={.f2=value}}); }
    inline void setUniform3f(const char* name, Vec3T<GLfloat> value){ setUniform(name, (uniform){.type=uniform::f3, .data={.f3=value}}); }
    inline void setUniform4f(const char* name, Vec4T<GLfloat> value){ setUniform(name, (uniform){.type=uniform::f4, .data={.f4=value}}); }

    inline void setUniform1i(const char* name, GLint value)         { setUniform(name, (uniform){.type=uniform::i1, .data={.i1=value}}); }
    inline void setUniform2i(const char* name, Vec2T<GLint> value)  { setUniform(name, (uniform){.type=uniform::i2, .data={.i2=value}}); }
    inline void setUniform3i(const char* name, Vec3T<GLint> value)  { setUniform(name, (uniform){.type=uniform::i3, .data={.i3=value}}); }
    inline void setUniform4i(const char* name, Vec4T<GLint> value)  { setUniform(name, (uniform){.type=uniform::i4, .data={.i4=value}}); }

    inline void setUniform1ui(const char* name, GLuint value)       { setUniform(name, (uniform){.type=uniform::ui1, .data={.ui1=value}}); }
    inline void setUniform2ui(const char* name, Vec2T<GLuint> value){ setUniform(name, (uniform){.type=uniform::ui2, .data={.ui2=value}}); }
    inline void setUniform3ui(const char* name, Vec3T<GLuint> value){ setUniform(name, (uniform){.type=uniform::ui3, .data={.ui3=value}}); }
    inline void setUniform4ui(const char* name, Vec4T<GLuint> value){ setUniform(name, (uniform){.type=uniform::ui4, .data={.ui4=value}}); }

    inline void setUniform3m(const char* name, Mat3T<GLfloat> value){ setUniform(name, (uniform){.type=uniform::m3, .data={.m3=value}}); }
    inline void setUniform4m(const char* name, Mat4T<GLfloat> value){ setUniform(name, (uniform){.type=uniform::m4, .data={.m4=value}}); }

private:
    GLuint compileShader(GLenum shaderType, const char* source);
    GLuint linkProgram(GLuint vertexShader, GLuint fragmentShader);
};







// default shaders
constexpr const std::string buildDefaultVertexShaderSource(unsigned int attrib_flags){
    std::string source = std::string("");
    source += "#version 100\n";

    // uniforms
    source += "uniform mat4 uMVPMatrix;\n";

    // attributes
    source += "attribute vec4 vPosition;\n";
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "attribute vec3 vNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "attribute vec3 vColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "attribute vec2 vUV    ;\n";

    // varyings
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "varying vec3 fNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "varying vec3 fColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "varying vec2 fUV    ;\n";

    // void main()
    source += "void main() {\n";
    source += "gl_Position = uMVPMatrix * vPosition;\n";
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "fNormal = vNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "fColor  = vColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "fUV     = vUV    ;\n";
    source += "}";

    return source;
}

constexpr const std::string buildDefaultFragmentShaderSource(unsigned int attrib_flags){
    std::string source = std::string("");
    source += "#version 100\n";

    // uniforms
    source += "uniform mediump vec3 uColor;\n";
    if (attrib_flags & GLMESH_FLAG_TEX) source += "uniform sampler2D uTexture;\n";

    // varyings
    if (attrib_flags & GLMESH_FLAG_NORMAL) source += "varying mediump vec3 fNormal;\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "varying mediump vec3 fColor ;\n";
    if (attrib_flags & GLMESH_FLAG_UV    ) source += "varying mediump vec2 fUV    ;\n";

    // void main()
    source += "void main() {\n";
    source += "gl_FragColor = vec4(uColor, 1.0);\n";
    if (attrib_flags & GLMESH_FLAG_TEX && attrib_flags & GLMESH_FLAG_UV) source += "gl_FragColor = gl_FragColor*texture2D(uTexture, fUV);\n";
    if (attrib_flags & GLMESH_FLAG_COLOR ) source += "gl_FragColor = gl_FragColor*vec4(fColor, 1.0);\n";
    if (attrib_flags & GLMESH_FLAG_NORMAL){
        source += "mediump float light = dot(fNormal, vec3(1.0, -1.0, 1.0));\n";
        source += "light = (light+1.0)/2.0;\n"; // normalised to range (0; 1)
        source += "light = 0.3 + light*0.6;\n"; // to range (0.3, 1.1)
        //source += "gl_FragColor = gl_FragColor*vec4(light, light, light, 1.0);\n"; // TODO: remove or implement lighting
    }
    source += "}";

    return source;
}

template <unsigned int attrib_flags>
Shader<attrib_flags>* defaultShader = new Shader<attrib_flags>(buildDefaultVertexShaderSource(attrib_flags).c_str(), buildDefaultFragmentShaderSource(attrib_flags).c_str());

#endif // SHADER_H
