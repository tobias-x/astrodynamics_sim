#include "file_io.h"
#include "simulation_engine.h"
#include "renderer.h"
#include "opengl_utils.h"
#include "types.h"

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <iostream>
#include <thread>
#include <fstream>
#include <sstream>
#include <vector>
#include <limits>

std::vector<Object> rendererObjects;
glm::vec3 cameraPos = glm::vec3(0.0f, 0.0f, 1.0f);
float yaw = 90.0f;
float pitch = -20.0f;
float radius = 20.0f;

static std::string LoadShaderSource(const std::string& filepath) {
    std::ifstream file(filepath);
    if (!file.is_open()) {
        throw std::runtime_error("Failed to open shader file: " + filepath);
    }
    std::stringstream buffer;
    buffer << file.rdbuf();
    return buffer.str();
}

static bool dragging = false;
static double lastX = 400, lastY = 300;

static void cursorPositionCallback(GLFWwindow* window, double xpos, double ypos) {
    if (!dragging) return;
    float dx = float(xpos - lastX);
    float dy = float(ypos - lastY);
    yaw += dx * 0.2f;
    pitch -= dy * 0.2f;
    pitch = glm::clamp(pitch, -89.0f, 89.0f);
    lastX = xpos; lastY = ypos;
}

static void mouseButtonCallback(GLFWwindow* window, int button, int action, int) {
    if (button == GLFW_MOUSE_BUTTON_LEFT) {
        if (action == GLFW_PRESS) {
            dragging = true;
            glfwGetCursorPos(window, &lastX, &lastY);
        } else {
            dragging = false;
        }
    }
}

static void scrollCallback(GLFWwindow* window, double, double yoffset) {
    radius -= float(yoffset);
    radius = glm::clamp(radius, 2.0f, 100.0f);
}

std::vector<GridBody> createInitialGridBodies() {
    std::vector<GridBody> gridBodies;
    float half = 10.0f;
    int divisions = 25;
    float step = (half * 2.0f) / divisions;
    float scale = 1e11f;

    for (int i = 0; i <= divisions; ++i) {
        float z = -half + i * step;
        for (int j = 0; j <= divisions; ++j) {
            float x = -half + j * step;
            gridBodies.push_back(GridBody{
                x * scale, 0.0f, z * scale,
                x * scale, 0.0f, z * scale,
                i, j  // record (i, j) indices
            });
        }
    }

    return gridBodies;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Usage: geodesic_gl <config.csv> <dt>\n";
        return 1;
    }
    double dt = std::stod(argv[2]);
    auto bodies = parseConfig(argv[1]);
    if (bodies.empty()) {
        std::cerr << "No bodies\n";
        return 1;
    }

    auto gridBodies = createInitialGridBodies();
    SimulationEngine engine(std::move(bodies), std::move(gridBodies), dt);
    std::thread simThread(&SimulationEngine::run, &engine);

    GLFWwindow* window = StartGLU();
    if (!window) return 2;

    auto vsSrc = LoadShaderSource("shaders/vertex_debug.glsl");
    auto fsSrc = LoadShaderSource("shaders/fragment_debug.glsl");
    GLuint shader = CreateShaderProgram(vsSrc.c_str(), fsSrc.c_str());
    glUseProgram(shader);

    glfwSetKeyCallback(window, keyCallback);
    glfwSetCursorPosCallback(window, cursorPositionCallback);
    glfwSetScrollCallback(window, scrollCallback);
    glfwSetMouseButtonCallback(window, mouseButtonCallback);

    GLint projLoc = glGetUniformLocation(shader, "projection");
    glm::mat4 proj = glm::perspective(glm::radians(45.0f), 800.0f / 600.0f, 0.1f, 100.0f);
    glUniformMatrix4fv(projLoc, 1, GL_FALSE, glm::value_ptr(proj));

    const float scale = 1e11f;

    // Setup grid VAO/VBO once
    GLuint gridVAO = 0, gridVBO = 0;
    glGenVertexArrays(1, &gridVAO);
    glGenBuffers(1, &gridVBO);

    glBindVertexArray(gridVAO);
    glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float) * 3 * 10000, nullptr, GL_DYNAMIC_DRAW);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);
    glBindVertexArray(0);

    auto initialState = engine.latest();
    rendererObjects.reserve(initialState.size());
    for (auto& b : initialState) {
        rendererObjects.emplace_back(b);
    }

    while (!glfwWindowShouldClose(window)) {
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        float ry = glm::radians(yaw), rp = glm::radians(pitch);
        cameraPos.x = radius * cos(rp) * sin(ry);
        cameraPos.y = radius * sin(rp);
        cameraPos.z = radius * cos(rp) * cos(ry);
        UpdateCam(shader, cameraPos);

        auto gridState = engine.latestGrid();

        // Upload grid lines (connect neighbors)
        std::vector<float> gridVerts;
        gridVerts.reserve(gridState.size() * 6); // two vertices per line

        for (const auto& gb : gridState) {
            for (const auto& neighbor : gridState) {
                if (neighbor.i == gb.i && neighbor.j == gb.j + 1) { // right neighbor
                    gridVerts.push_back(gb.x / scale);
                    gridVerts.push_back(gb.y / scale);
                    gridVerts.push_back(gb.z / scale);
                    gridVerts.push_back(neighbor.x / scale);
                    gridVerts.push_back(neighbor.y / scale);
                    gridVerts.push_back(neighbor.z / scale);
                }
                if (neighbor.i == gb.i + 1 && neighbor.j == gb.j) { // bottom neighbor
                    gridVerts.push_back(gb.x / scale);
                    gridVerts.push_back(gb.y / scale);
                    gridVerts.push_back(gb.z / scale);
                    gridVerts.push_back(neighbor.x / scale);
                    gridVerts.push_back(neighbor.y / scale);
                    gridVerts.push_back(neighbor.z / scale);
                }
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER, gridVBO);
        glBufferSubData(GL_ARRAY_BUFFER, 0, gridVerts.size() * sizeof(float), gridVerts.data());

        glm::mat4 identity = glm::mat4(1.0f);
        glUniformMatrix4fv(glGetUniformLocation(shader, "model"), 1, GL_FALSE, glm::value_ptr(identity));

        // Draw grid
        glUseProgram(shader);
        glUniform4f(glGetUniformLocation(shader, "objectColor"), 1.0f, 1.0f, 1.0f, 1.0f);
        glUniform1i(glGetUniformLocation(shader, "isGrid"), 1);
        glUniform1i(glGetUniformLocation(shader, "GLOW"), 0);

        glLineWidth(1.0f);
        glBindVertexArray(gridVAO);
        glDrawArrays(GL_LINES, 0, gridVerts.size() / 3);
        glBindVertexArray(0);

        // Draw planets
        auto state = engine.latest();

        // for (const auto& s : state) {
        //     std::cout << "Rendered body: " << s.name
        //               << " position: (" << s.x / scale << ", " << s.y / scale << ", " << s.z / scale << ")\n";
        // }

        for (size_t i = 0; i < state.size(); ++i) {
            glm::vec3 pos = {
                state[i].x / scale,
                state[i].y / scale,
                state[i].z / scale
            };
            
            // Set color
            if (state[i].name == "Sun") {
                rendererObjects[i].setColor({1.0f, 1.0f, 0.0f, 1.0f}); // Yellow (Red + Green)
            }
            else if (state[i].name == "Earth") {
                rendererObjects[i].setColor({0.0f, 0.0f, 1.0f, 1.0f}); // Blue
            }
            else {
                rendererObjects[i].setColor({1.0f, 0.0f, 0.0f, 1.0f}); // Red
            }
        
            // Set position
            rendererObjects[i].setPosition(pos);
        
            // Log position
            // std::cout << "Planet " << state[i].name << " position: ("
            //           << state[i].x / scale << ", " << state[i].y / scale << ", " << state[i].z / scale << ")\n";
        
            // Draw
            rendererObjects[i].draw(shader);
        }        

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    engine.stop();
    simThread.join();
    glfwTerminate();
    return 0;
}
