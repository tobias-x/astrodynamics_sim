#include "opengl_utils.h"

#include <iostream>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

GLFWwindow* StartGLU() {
    if (!glfwInit()) {
        std::cerr << "GLFW init failure\n";
        return nullptr;
    }
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
#endif

    GLFWwindow* w = glfwCreateWindow(800, 600, "GeodesicSim", nullptr, nullptr);
    if (!w) {
        std::cerr << "GLFW window creation failure\n";
        glfwTerminate();
        return nullptr;
    }
    glfwMakeContextCurrent(w);

    glewExperimental = true;
    if (glewInit() != GLEW_OK) {
        std::cerr << "GLEW init failure\n";
        glfwTerminate();
        return nullptr;
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    return w;
}

GLuint CreateShaderProgram(const char* vsSrc, const char* fsSrc) {
    auto compile = [&](GLenum type, const char* src) {
        GLuint s = glCreateShader(type);
        glShaderSource(s, 1, &src, nullptr);
        glCompileShader(s);
        GLint ok; glGetShaderiv(s, GL_COMPILE_STATUS, &ok);
        if (!ok) {
            char buf[512];
            glGetShaderInfoLog(s, 512, nullptr, buf);
            std::cerr << "Shader compile error:\n" << buf << "\n";
        }
        return s;
    };
    GLuint vs = compile(GL_VERTEX_SHADER, vsSrc);
    GLuint fs = compile(GL_FRAGMENT_SHADER, fsSrc);

    GLuint prog = glCreateProgram();
    glAttachShader(prog, vs);
    glAttachShader(prog, fs);
    glLinkProgram(prog);
    GLint ok; glGetProgramiv(prog, GL_LINK_STATUS, &ok);
    if (!ok) {
        char buf[512];
        glGetProgramInfoLog(prog, 512, nullptr, buf);
        std::cerr << "Shader link error:\n" << buf << "\n";
    }

    glDeleteShader(vs);
    glDeleteShader(fs);
    return prog;
}

void UpdateCam(GLuint shaderProgram, const glm::vec3& cameraPos) {
    glm::mat4 view = glm::lookAt(cameraPos, glm::vec3(0.0f), glm::vec3(0.0f, 1.0f, 0.0f));
    GLint loc = glGetUniformLocation(shaderProgram, "view");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(view));
}

void keyCallback(GLFWwindow* w, int key, int, int action, int) {
    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
        glfwSetWindowShouldClose(w, GLFW_TRUE);
}

void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const float* verts, size_t count) {
    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glBindVertexArray(VAO);
      glBindBuffer(GL_ARRAY_BUFFER, VBO);
      glBufferData(GL_ARRAY_BUFFER, count * sizeof(float), verts, GL_STATIC_DRAW);
      glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), nullptr);
      glEnableVertexAttribArray(0);
    glBindVertexArray(0);
}

glm::vec3 sphericalToCartesian(float r, float theta, float phi) {
    return {
        r * std::sin(theta) * std::cos(phi),
        r * std::cos(theta),
        r * std::sin(theta) * std::sin(phi)
    };
}
