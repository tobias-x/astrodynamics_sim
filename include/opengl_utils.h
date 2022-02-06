#pragma once
#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <glm/glm.hpp>
#include <vector>
#include "types.h"

GLFWwindow* StartGLU();
GLuint CreateShaderProgram(const char* vsSrc, const char* fsSrc);
void CreateVBOVAO(GLuint& VAO, GLuint& VBO, const float* verts, size_t count);
glm::vec3 sphericalToCartesian(float r, float theta, float phi);
void UpdateCam(GLuint shaderProgram, const glm::vec3& cameraPos);
void keyCallback(GLFWwindow* w, int key, int scancode, int action, int mods);