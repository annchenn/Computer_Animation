#include "graphic.h"

#include <algorithm>

#include "Eigen/Geometry"
#include "glad/glad.h"

#include "../simulation/particle.h"
#include "../util/helper.h"

namespace {
constexpr int g_SphereSectors = 50;
constexpr int g_SphereStacks = 50;
constexpr int g_SphereRadius = 1;
}  // namespace

namespace gfx {

std::vector<GLfloat> Plane::vertices;
std::vector<GLfloat> Elevator::vertices; 

std::vector<GLfloat> SkyBox::vertices = {-1.0f, -1.0f, -1.0f, 1.0f, -1.0f, -1.0f, 1.0f, 1.0f, -1.0f, -1.0f, 1.0f, -1.0f,
                                         -1.0f, -1.0f, 1.0f,  1.0f, -1.0f, 1.0f,  1.0f, 1.0f, 1.0f,  -1.0f, 1.0f, 1.0f};
std::vector<GLuint> SkyBox::indices = {0, 1, 3, 3, 1, 2, 1, 5, 2, 2, 5, 6, 5, 4, 6, 6, 4, 7,
                                       4, 0, 7, 7, 0, 3, 3, 2, 7, 7, 2, 6, 4, 5, 0, 0, 5, 1};

Rigidbody::Rigidbody() {
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);

    glBindVertexArray(vao);
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    // Vertex
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), nullptr);
    // normal
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void*>(3 * sizeof(GLfloat)));
    // Texture
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(GLfloat), reinterpret_cast<void*>(6 * sizeof(GLfloat)));
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

Rigidbody::~Rigidbody() {
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
}

void Rigidbody::setModelMatrix(const Eigen::Matrix4f& _modelMatrix) {
    modelMatrix = _modelMatrix;
    inverseTransposeModel = modelMatrix.block(0, 0, 3, 3).inverse().transpose();
}

void Rigidbody::setTexture(const std::shared_ptr<Texture>& _texture) { texture = _texture; }

void Rigidbody::setTexture(const Eigen::Vector3f& color) { baseColor = color; }

void Rigidbody::initialize() {
    generateVertices();
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertexPointer->size() * sizeof(GLfloat), vertexPointer->data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

Plane::Plane() {
    if (vertexPointer == nullptr) {
        initialize();
    }
}

void Plane::generateVertices() {
    vertices = {-0.5f, 0.0f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f, 1.0f, -0.5f, 0.0f, 0.5f,  0.0f, 1.0f, 0.0f, 1.0f, 1.0f,
                0.5f,  0.0f, 0.5f,  0.0f, 1.0f, 0.0f, 1.0f, 0.0f, 0.5f,  0.0f, -0.5f, 0.0f, 1.0f, 0.0f, 0.0f, 0.0f};
    vertexPointer = &vertices;
}

void Plane::render(Program* shaderProgram) {
    if (texture) {
        shaderProgram->setUniform("useTexture", 1);
        shaderProgram->setUniform("diffuseTexture", texture->getIndex());
    } else {

        shaderProgram->setUniform("useTexture", 0);
        shaderProgram->setUniform("baseColor", baseColor);
    }
    shaderProgram->setUniform("obj", 2);
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}

SoftJelly::SoftJelly(simulation::Jelly* _jelly) {
    jelly = _jelly;

    glGenVertexArrays(6, jellyVAO);
    glGenVertexArrays(1, &structVAO);
    glGenVertexArrays(1, &shearVAO);
    glGenVertexArrays(1, &bendingVAO);
    glGenVertexArrays(1, &particleVAO);

    glGenBuffers(1, &vertexVBO);
    glGenBuffers(1, &textureCoordVBO);
    glGenBuffers(6, jellyVBO);
    glGenBuffers(6, normalVBO);

    glGenBuffers(1, &structEBO);
    glGenBuffers(1, &bendingEBO);
    glGenBuffers(1, &shearEBO);
    glGenBuffers(2, jellyEBO);

    allocateBuffers();
    calculateTextureCoords();
    calculateIndices();
    initializeBuffers();
    update();
}

SoftJelly::~SoftJelly() {
    glDeleteVertexArrays(6, jellyVAO);
    glDeleteVertexArrays(1, &structVAO);
    glDeleteVertexArrays(1, &shearVAO);
    glDeleteVertexArrays(1, &bendingVAO);
    glDeleteVertexArrays(1, &particleVAO);

    glDeleteBuffers(1, &vertexVBO);
    glDeleteBuffers(1, &textureCoordVBO);
    glDeleteBuffers(6, jellyVBO);
    glDeleteBuffers(6, normalVBO);

    glDeleteBuffers(1, &structEBO);
    glDeleteBuffers(1, &bendingEBO);
    glDeleteBuffers(1, &shearEBO);
    glDeleteBuffers(2, jellyEBO);
}

void SoftJelly::allocateBuffers() {
    int edgeNum = jelly->getNumAtEdge();
    int totalSize = edgeNum * edgeNum * sizeof(GLfloat);
    int vertexSize = 3 * jelly->getParticleNum() * sizeof(GLfloat);
    int indexSize = 6 * edgeNum * edgeNum * sizeof(GLuint);
    // Temporary
    faceVertices = std::vector<GLfloat>(3 * edgeNum * edgeNum);
    normalBuffer = std::vector<float>(3 * edgeNum * edgeNum);
    fullVertices = std::vector<float>(3 * jelly->getParticleNum());
    // All vertices
    glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
    glBufferData(GL_ARRAY_BUFFER, vertexSize, nullptr, GL_STREAM_DRAW);
    // Texture
    glBindBuffer(GL_ARRAY_BUFFER, textureCoordVBO);
    glBufferData(GL_ARRAY_BUFFER, 2 * totalSize, nullptr, GL_STATIC_DRAW);

    for (int face = 0; face < 6; ++face) {
        glBindBuffer(GL_ARRAY_BUFFER, jellyVBO[face]);
        glBufferData(GL_ARRAY_BUFFER, 3 * totalSize, nullptr, GL_STREAM_DRAW);
        // Normal
        glBindBuffer(GL_ARRAY_BUFFER, normalVBO[face]);
        glBufferData(GL_ARRAY_BUFFER, 3 * totalSize, nullptr, GL_STREAM_DRAW);
    }
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    // Front and back of cube
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, jellyEBO[0]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexSize, nullptr, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, jellyEBO[1]);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indexSize, nullptr, GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void SoftJelly::calculateIndices() {
    int edgeNum = jelly->getNumAtEdge();
    int dividor = edgeNum - 1;
    {
        std::vector<GLuint> temp(6 * (edgeNum - 1) * (edgeNum - 1));
        int currentPos = -1;
        for (int i = 0; i < edgeNum - 1; ++i) {
            for (int j = 0; j < edgeNum - 1; ++j) {
                unsigned int k1 = i * edgeNum + j;
                unsigned int k2 = (i + 1) * edgeNum + j;
                temp[++currentPos] = k2 + 1;
                temp[++currentPos] = k2;
                temp[++currentPos] = k1;
                temp[++currentPos] = k1 + 1;
                temp[++currentPos] = k2 + 1;
                temp[++currentPos] = k1;
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, jellyEBO[0]);
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, temp.size() * sizeof(GLuint), temp.data());
        currentPos = -1;
        for (int i = 0; i < edgeNum - 1; ++i) {
            for (int j = 0; j < edgeNum - 1; ++j) {
                unsigned int k1 = i * edgeNum + j;
                unsigned int k2 = (i + 1) * edgeNum + j;
                temp[++currentPos] = k1;
                temp[++currentPos] = k2;
                temp[++currentPos] = k2 + 1;
                temp[++currentPos] = k1;
                temp[++currentPos] = k2 + 1;
                temp[++currentPos] = k1 + 1;
            }
        }
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, jellyEBO[1]);
        glBufferSubData(GL_ELEMENT_ARRAY_BUFFER, 0, temp.size() * sizeof(GLuint), temp.data());
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
    }
    for (int i = 0; i < jelly->getSpringNum(); ++i) {
        std::vector<GLuint>* currentIndices = &structIndices;
        switch (jelly->getSpring(i).getType()) {
            case simulation::Spring::SpringType::STRUCT:
                currentIndices = &structIndices;
                break;
            case simulation::Spring::SpringType::SHEAR:
                currentIndices = &shearIndices;
                break;
            case simulation::Spring::SpringType::BENDING:
                currentIndices = &bendingIndices;
                break;
        }
        currentIndices->emplace_back(jelly->getSpring(i).getSpringStartID());
        currentIndices->emplace_back(jelly->getSpring(i).getSpringEndID());
    }
    // Spring indices
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, structEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, structIndices.size() * sizeof(GLuint), structIndices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shearEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, shearIndices.size() * sizeof(GLuint), shearIndices.data(), GL_STATIC_DRAW);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bendingEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, bendingIndices.size() * sizeof(GLuint), bendingIndices.data(),
                 GL_STATIC_DRAW);
}

void SoftJelly::calculateTextureCoords() {
    int edgeNum = jelly->getNumAtEdge();
    int dividor = edgeNum - 1;
    std::vector<GLfloat> temp(2 * edgeNum * edgeNum);
    int currentPos = -1;
    for (int i = 0; i < edgeNum; ++i) {
        for (int j = 0; j < edgeNum; ++j) {
            temp[++currentPos] = static_cast<GLfloat>(i) / dividor;
            temp[++currentPos] = static_cast<GLfloat>(j) / dividor;
        }
    }
    glBindBuffer(GL_ARRAY_BUFFER, textureCoordVBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, temp.size() * sizeof(GLfloat), temp.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void SoftJelly::initializeBuffers() {
    for (int face = 0; face < 6; ++face) {
        glBindVertexArray(jellyVAO[face]);
        // Vertex
        glBindBuffer(GL_ARRAY_BUFFER, jellyVBO[face]);
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
        // Normal
        glBindBuffer(GL_ARRAY_BUFFER, normalVBO[face]);
        glEnableVertexAttribArray(1);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
        // Texture
        glBindBuffer(GL_ARRAY_BUFFER, textureCoordVBO);
        glEnableVertexAttribArray(2);
        glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GLfloat), nullptr);
        // Index
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, jellyEBO[face & 1]);
    }

    // Spring - Struct
    glBindVertexArray(structVAO);
    glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, structEBO);
    // Spring - Shear
    glBindVertexArray(shearVAO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, shearEBO);
    // Spring - Bending
    glBindVertexArray(bendingVAO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, bendingEBO);
    // Particles
    glBindVertexArray(particleVAO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
    // Unbind all
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

void SoftJelly::setTextures(const std::array<std::shared_ptr<Texture>, 6>& _textures) { textures = _textures; }

void SoftJelly::update() {
    // All vertex
    int nParticle = jelly->getParticleNum();
    auto particles = jelly->getParticlePointer();
    int currentPos = -1;
    for (int i = 0; i < nParticle; ++i) {
        const auto& pos = (*particles)[i].getPosition();
        fullVertices[++currentPos] = pos[0];
        fullVertices[++currentPos] = pos[1];
        fullVertices[++currentPos] = pos[2];
    }
    glBindBuffer(GL_ARRAY_BUFFER, vertexVBO);
    glBufferSubData(GL_ARRAY_BUFFER, 0, fullVertices.size() * sizeof(GLfloat), fullVertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    for (int face = 1; face <= 6; ++face) {
        updateSingleFace(face);
    }
}

void SoftJelly::updateSingleFace(int face) {
    int edgeNum = jelly->getNumAtEdge();
    // Vertex per face
    int currentPos = -1;
    for (int i = 0; i < edgeNum; ++i) {
        for (int j = 0; j < edgeNum; ++j) {
            int pos = 3 * jelly->getPointMap(face, i, j);
            faceVertices[++currentPos] = fullVertices[pos++];
            faceVertices[++currentPos] = fullVertices[pos++];
            faceVertices[++currentPos] = fullVertices[pos];
        }
    }

    glBindBuffer(GL_ARRAY_BUFFER, jellyVBO[face - 1]);
    glBufferSubData(GL_ARRAY_BUFFER, 0, faceVertices.size() * sizeof(GLfloat), faceVertices.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    // Reset to 0
    auto temp = std::vector<std::vector<Eigen::Vector3f>>(
        edgeNum, std::vector<Eigen::Vector3f>(edgeNum, Eigen::Vector3f::Zero()));
    // Use to flip normal
    float dFaceFactor = (face & 1) ? -1.0f : 1.0f;
    // Accumulate normals
    for (int i = 0; i < edgeNum - 1; ++i) {
        for (int j = 0; j < edgeNum - 1; ++j) {
            // V1   V4
            // V2   V3
            Eigen::Vector3f V1 = jelly->getParticle(jelly->getPointMap(face, i, j)).getPosition();
            Eigen::Vector3f V2 = jelly->getParticle(jelly->getPointMap(face, i + 1, j)).getPosition();
            Eigen::Vector3f V3 = jelly->getParticle(jelly->getPointMap(face, i + 1, j + 1)).getPosition();
            Eigen::Vector3f V4 = jelly->getParticle(jelly->getPointMap(face, i, j + 1)).getPosition();

            Eigen::Vector3f N = (V3 - V2).cross(V1 - V2);
            N *= dFaceFactor;
            N.normalize();
            temp[i][j] += N;
            temp[i + 1][j] += N;
            temp[i + 1][j + 1] += N;

            N = (V1 - V4).cross(V3 - V4);
            N *= dFaceFactor;
            N.normalize();
            temp[i][j] += N;
            temp[i + 1][j + 1] += N;
            temp[i][j + 1] += N;
        }
    }
    // Update normal buffer
    currentPos = -1;
    for (int i = 0; i < edgeNum; ++i) {
        for (int j = 0; j < edgeNum; ++j) {
            temp[i][j].normalize();
            normalBuffer[++currentPos] = temp[i][j][0];
            normalBuffer[++currentPos] = temp[i][j][1];
            normalBuffer[++currentPos] = temp[i][j][2];
        }
    }
    glBindBuffer(GL_ARRAY_BUFFER, normalVBO[face - 1]);
    glBufferSubData(GL_ARRAY_BUFFER, 0, normalBuffer.size() * sizeof(GLfloat), normalBuffer.data());
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void SoftJelly::renderJelly(Program* shaderProgram) {
    shaderProgram->setUniform("obj", 3);
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);
    shaderProgram->setUniform("useTexture", 1);
    int dividor = jelly->getNumAtEdge() - 1;
    for (int face = 0; face < 6; ++face) {
        glBindVertexArray(jellyVAO[face]);
        shaderProgram->setUniform("diffuseTexture", textures[face]->getIndex());
        glDrawElements(GL_TRIANGLES, 6 * dividor * dividor, GL_UNSIGNED_INT, 0);
    }
    glBindVertexArray(0);
}

void SoftJelly::renderLines(Program* shaderProgram, simulation::Spring::SpringType springType) {
    using SpringType = simulation::Spring::SpringType;
    shaderProgram->setUniform("useTexture", 0);
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);
    std::size_t currentEBOSize = 0;
    switch (springType) {
        case SpringType::STRUCT:
            shaderProgram->setUniform("baseColor", Eigen::Vector3f(0.0f, 1.0f, 1.0f));
            glBindVertexArray(structVAO);
            currentEBOSize = structIndices.size();
            break;
        case SpringType::SHEAR:
            shaderProgram->setUniform("baseColor", Eigen::Vector3f(1.0f, 1.0f, 0.0f));
            glBindVertexArray(shearVAO);
            currentEBOSize = shearIndices.size();
            break;
        case SpringType::BENDING:
            shaderProgram->setUniform("baseColor", Eigen::Vector3f(0.0f, 1.0f, 0.0f));
            glBindVertexArray(bendingVAO);
            currentEBOSize = bendingIndices.size();
    }
    glDrawElements(GL_LINES, static_cast<GLsizei>(currentEBOSize), GL_UNSIGNED_INT, 0);
    glBindVertexArray(0);
}

void SoftJelly::renderPoints(Program* shaderProgram) {
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);
    shaderProgram->setUniform("useTexture", 0);
    shaderProgram->setUniform("baseColor", Eigen::Vector3f(1.0f, 0.0f, 0.0f));
    glBindVertexArray(particleVAO);
    glDrawArrays(GL_POINTS, 0, jelly->getParticleNum());
    glBindVertexArray(0);
}

SkyBox::SkyBox() {
    glGenVertexArrays(1, &vao);
    glGenBuffers(1, &vbo);
    glGenBuffers(1, &ebo);
    glBindVertexArray(vao);
    // Vertex
    glBindBuffer(GL_ARRAY_BUFFER, vbo);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(GLfloat), vertices.data(), GL_STATIC_DRAW);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * sizeof(GLfloat), nullptr);
    // Index
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(GLuint), indices.data(), GL_STATIC_DRAW);
    // Unbind all
    glBindVertexArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
}

SkyBox::~SkyBox() {
    glDeleteVertexArrays(1, &vao);
    glDeleteBuffers(1, &vbo);
    glDeleteBuffers(1, &ebo);
}

void SkyBox::setTexture(const std::shared_ptr<CubeTexture>& _texture) { texture = _texture; }

void SkyBox::render(Program* shaderProgram) {
    glDepthFunc(GL_LEQUAL);  // Put on the background
    glBindVertexArray(vao);
    shaderProgram->setUniform("skybox", texture->getIndex());
    glDrawElements(GL_TRIANGLES, static_cast<GLsizei>(indices.size()), GL_UNSIGNED_INT, nullptr);
    glBindVertexArray(0);
    glDepthFunc(GL_LESS);
}

Elevator::Elevator() {
    if (vertexPointer == nullptr) {
        initialize();
    }
}

void Elevator::generateVertices() {
    vertices = {-0.5f, 0.0f, -0.5f, 
                0.0f, 1.0f, 0.0f, 
                0.0f, 1.0f, 
                -0.5f, 0.0f, 0.5f, 
                0.0f, 1.0f, 0.0f, 
                1.0f, 1.0f,
                0.5f,  0.0f, 0.5f,  
                0.0f, 1.0f, 0.0f, 
                1.0f, 0.0f, 
                0.5f,  0.0f, -0.5f, 
                0.0f, 1.0f, 0.0f, 
                0.0f, 0.0f};
    vertexPointer = &vertices;
}

void Elevator::render(Program* shaderProgram) {
    if (texture) {
        shaderProgram->setUniform("useTexture", 1);
        shaderProgram->setUniform("diffuseTexture", texture->getIndex());
    } else {
        shaderProgram->setUniform("useTexture", 0);
        shaderProgram->setUniform("baseColor", baseColor);
    }
    shaderProgram->setUniform("obj", 2);
    shaderProgram->setUniform("model", modelMatrix);
    shaderProgram->setUniform("invtransmodel", inverseTransposeModel);
    glBindVertexArray(vao);
    glDrawArrays(GL_TRIANGLE_FAN, 0, 4);
}

}  // namespace gfx
