#include "renderer.h"
#include "opengl_utils.h"
#include <glm/gtc/constants.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>

Object::Object(const Body& b)
  : position(b.x, b.y, b.z),
    velocity(b.vx, b.vy, b.vz),
    density(5515.0f),
    color(1, 1, 1, 1)
{
    radius = 0.2;
    auto mesh = buildMesh();
    vertexCount = mesh.size() / 3; // <<< 🔥 Corrected
    CreateVBOVAO(VAO, VBO, mesh.data(), mesh.size());
}

void Object::setColor(const glm::vec4& c) {
    color = c;
}

Object::~Object(){
    glDeleteBuffers(1, &VBO);
    glDeleteVertexArrays(1, &VAO);
}

std::vector<float> Object::buildMesh() const {
    std::vector<float> verts;
    int stacks = 10;
    int sectors = 10;

    for (int i = 0; i < stacks; ++i) {
        float theta1 = (float)i / stacks * glm::pi<float>();
        float theta2 = (float)(i + 1) / stacks * glm::pi<float>();

        for (int j = 0; j < sectors; ++j) {
            float phi1 = (float)j / sectors * 2.0f * glm::pi<float>();
            float phi2 = (float)(j + 1) / sectors * 2.0f * glm::pi<float>();

            glm::vec3 v1 = sphericalToCartesian(radius, theta1, phi1);
            glm::vec3 v2 = sphericalToCartesian(radius, theta1, phi2);
            glm::vec3 v3 = sphericalToCartesian(radius, theta2, phi1);
            glm::vec3 v4 = sphericalToCartesian(radius, theta2, phi2);

            if (i != 0) {
                verts.insert(verts.end(), {
                    v1.x, v1.y, v1.z,
                    v2.x, v2.y, v2.z,
                    v3.x, v3.y, v3.z
                });
            }

            if (i != (stacks - 1)) {
                verts.insert(verts.end(), {
                    v2.x, v2.y, v2.z,
                    v4.x, v4.y, v4.z,
                    v3.x, v3.y, v3.z
                });
            }
        }
    }

    return verts;
}


void Object::updateVertices(){
    auto mesh = buildMesh();
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, mesh.size() * sizeof(float), mesh.data(), GL_STATIC_DRAW);
}

void Object::setPosition(const glm::vec3& p){
    position = p;
}

void Object::draw(GLuint shader) const {
    glUseProgram(shader); // 🔥 Make sure the shader is bound (redundant but safe)

    glm::mat4 model = glm::translate(glm::mat4(1.0f), position);
    GLint loc = glGetUniformLocation(shader, "model");
    glUniformMatrix4fv(loc, 1, GL_FALSE, glm::value_ptr(model));
    glUniform4fv(glGetUniformLocation(shader, "objectColor"), 1, glm::value_ptr(color));
    glUniform1i(glGetUniformLocation(shader, "isGrid"), 0);

    glBindVertexArray(VAO); // 🔥 Bind the object's VAO
    glDrawArrays(GL_TRIANGLES, 0, vertexCount);
    glBindVertexArray(0);   // 🔥 Unbind VAO after drawing to avoid messing up next draw
}

