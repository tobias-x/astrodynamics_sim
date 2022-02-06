#pragma once
#include <vector>
#include "file_io.h"
#include <GL/glew.h>
#include <glm/glm.hpp>
class Object {
public:
    Object(const Body& b);
    ~Object();
    void updateVertices();
    void setPosition(const glm::vec3& p);
    void draw(GLuint shaderProgram) const;
    void setColor(const glm::vec4& c);
    GLuint getVAO() const { return VAO; }
    GLsizei getVertexCount() const { return vertexCount; }


private:
    GLuint VAO, VBO;
    size_t vertexCount;
    glm::vec4 color;
    float     radius;
    glm::vec3 position, velocity;
    float density;

    std::vector<float> buildMesh() const;
};
