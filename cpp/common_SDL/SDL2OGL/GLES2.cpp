#include "GLES2.h"

CameraT<float>* GLES2::active_camera = nullptr;
GLuint GLES2::currentGL_ARRAY_BUFFER = 0;
Vec2i GLES2::screen_size = {1820, 980};
