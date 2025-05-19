#include "terrain.h"

#include <stdexcept>

#include "../util/helper.h"
#include <iostream>

namespace simulation {
// Factory
std::unique_ptr<Terrain> TerrainFactory::CreateTerrain(TerrainType type) {
    switch (type) {
        case simulation::TerrainType::Plane:
            return std::make_unique<PlaneTerrain>();
        case simulation::TerrainType::Elevator:
            return std::make_unique<ElevatorTerrain>();

        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}
// Terrain

Eigen::Matrix4f Terrain::getModelMatrix() { return modelMatrix; }

float Terrain::getMass() const { return mass; }

Eigen::Vector3f Terrain::getPosition() const { return position; }

Eigen::Vector3f Terrain::getVelocity() const { return velocity; }

Eigen::Vector3f Terrain::getAcceleration() const { return force / mass; }

Eigen::Vector3f Terrain::getForce() const { return force; }

void Terrain::setMass(const float _mass) { mass = _mass; }

void Terrain::setPosition(const Eigen::Vector3f& _position) { 
    modelMatrix = util::translate(_position - position) * modelMatrix;
    position = _position;
}

void Terrain::setVelocity(const Eigen::Vector3f& _velocity) { velocity = _velocity; }

void Terrain::setAcceleration(const Eigen::Vector3f& _acceleration) { force = _acceleration * mass; }

void Terrain::setForce(const Eigen::Vector3f& _force) { force = _force; }

void Terrain::addPosition(const Eigen::Vector3f& _position) { 
    position += _position;
    modelMatrix = util::translate(_position) * modelMatrix;
}

void Terrain::addVelocity(const Eigen::Vector3f& _velocity) { velocity += _velocity; }

void Terrain::addAcceleration(const Eigen::Vector3f& _acceleration) { force += _acceleration * mass; }

void Terrain::addForce(const Eigen::Vector3f& _force) { force += _force; }

// Note:
// You should update each particles' velocity (base on the equation in
// slide) and force (contact force : resist + friction) in handleCollision function

// PlaneTerrain //

PlaneTerrain::PlaneTerrain() { reset(); }

void PlaneTerrain::reset() { 
    modelMatrix = util::translate(0.0f, position[1], 0.0f) * util::rotateDegree(0, 0, -20) * util::scale(30, 1, 30);
}


TerrainType PlaneTerrain::getType() { return TerrainType::Plane; }

void PlaneTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.5f;
    // TODO#3-1: Handle collision when a particle collide with the plane terrain.
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    //   3. The plane spans 30x30 units in the XZ plane and is rotated -20 degrees around the Z-axis.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);

        Eigen::Vector3f x = particle.getPosition();
        Eigen::Vector3f v = particle.getVelocity();

        float dis = normal.dot(x);

        if (dis < eEPSILON) {
        //collision occurs
            float v_along_nor = v.dot(normal);
            if (v_along_nor < 0) {
                
                Eigen::Vector3f vn = v_along_nor * normal;
                Eigen::Vector3f vt = v - vn;
                Eigen::Vector3f reflect_v = vt - coefResist * vn;
                particle.setVelocity(reflect_v);

                if (dis < 0) {
                    particle.addPosition(-dis * normal);
                }
            }

            if (std::abs(dis) < eEPSILON) {
                Eigen::Vector3f force = particle.getForce();
                float force_nor = force.dot(normal);

                //force pushes particle into the plane
                if (force_nor < 0) {
                    Eigen::Vector3f contactForce = -force_nor * normal;
                    ////fc = -(N¡Pf)N
                    Eigen::Vector3f vt = v - v.dot(normal) * normal;
                    if (vt.norm() > 1e-6) {
                        Eigen::Vector3f friction = -coefFriction * (-force_nor) * vt.normalized();
                        // f = -kf(-N¡Pf)vt
                        particle.addForce(contactForce + friction);
                    }
                    else {
                        particle.addForce(contactForce);
                    }
                }
            }
        }
        
    }
}

ElevatorTerrain::ElevatorTerrain() { 
    reset();
}

void ElevatorTerrain::reset() {
    modelMatrix = util::translate(0.0f, 1.0f, 0.0f) * util::rotateDegree(0, 0, 0) * util::scale(5, 1, 5);
    position = Eigen::Vector3f(0.0f, 1.0f, 0.0f);
    velocity = Eigen::Vector3f(0.0f, 0.0f, 0.0f);
}

TerrainType ElevatorTerrain::getType() { return TerrainType::Elevator; }

void ElevatorTerrain::handleCollision(const float delta_T, Jelly& jelly) {
    constexpr float eEPSILON = 0.01f;
    constexpr float coefResist = 0.8f;
    constexpr float coefFriction = 0.5f;

    // TODO#3-2: Implement the collision handling between the jelly and the elevator
    //   If collision happens:
    //      1. Directly update particles' velocity
    //      2. Apply contact force to particles when needed
    // Note:
    //   1. There are `jelly.getParticleNum()` particles.
    //   2. See TODOs in `Jelly::computeInternalForce` and functions in `particle.h` if you don't know how to access
    //   data.
    //   3. The elevator plane spans 5x5 units in the XZ plane.
    // Hint:
    //   1. Review "particles.pptx" from p.14 - p.19
    //   1. Use a.norm() to get length of a.
    //   2. Use a.normalize() to normalize a inplace.
    //          a.normalized() will create a new vector.
    //   3. Use a.dot(b) to get dot product of a and b.

    for (int i = 0; i < jelly.getParticleNum(); i++) {
        Particle& particle = jelly.getParticle(i);

        Eigen::Vector3f x = particle.getPosition();
        Eigen::Vector3f v = particle.getVelocity();

        float halfLen = 2.5f;
        //half width of elevator
        if (std::abs(x[0]) <= halfLen && std::abs(x[2]) <= halfLen) {
        //check whether the particle is in the range of elevator
            
            float dis = (x - position).dot(normal);
            // N¡P(x-p) < £`
            if (dis < eEPSILON) {
                Eigen::Vector3f relative_v = v - velocity;
                //relative velocity of particle and evaletor

                float nor_v = relative_v.dot(normal);

                if (nor_v < 0) {
                    //particle is moving toward the elevator
                    Eigen::Vector3f new_v = relative_v - (1 + coefResist) *nor_v * normal;

                    new_v += velocity;//add  velocity of elevator to get the absolute velocity
                    particle.setVelocity(new_v);

                    if (dis < 0) {
                        Eigen::Vector3f correction = -dis * normal;
                        particle.addPosition(correction);
                    }
                }

                if (std::abs(dis) < eEPSILON) {
                    //apply contact force

                    Eigen::Vector3f force = particle.getForce();

                    float force_nor = force.dot(normal);

                    //a force push particle toward the surface (N¡Pf < 0)
                    if (force_nor < 0) {
                        Eigen::Vector3f contactForce = -force_nor * normal;//fc = -(N¡Pf)N

                        Eigen::Vector3f normalComponent = nor_v * normal;
                        Eigen::Vector3f tangentComponent = relative_v - normalComponent;

                        Eigen::Vector3f friction;
                        if (tangentComponent.norm() > 1e-6) {
                            friction = -coefFriction * (-force_nor) * tangentComponent.normalized();
                            // f = -kf(-N¡Pf)vt
                        }
                        else {
                            friction = Eigen::Vector3f::Zero();
                        }

                        particle.addForce(contactForce + friction);
                    }
                }
            }
        }

    }
}

}  // namespace simulation
