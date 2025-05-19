#include "integrator.h"
#include <iostream>

#include <vector>
namespace simulation {
// Factory
std::unique_ptr<Integrator> IntegratorFactory::CreateIntegrator(IntegratorType type) {
    switch (type) {
        case simulation::IntegratorType::ExplicitEuler:
            return std::make_unique<ExplicitEulerIntegrator>();
        case simulation::IntegratorType::ImplicitEuler:
            return std::make_unique<ImplicitEulerIntegrator>();
        case simulation::IntegratorType::MidpointEuler:
            return std::make_unique<MidpointEulerIntegrator>();
        case simulation::IntegratorType::RungeKuttaFourth:
            return std::make_unique<RungeKuttaFourthIntegrator>();
        default:
            throw std::invalid_argument("TerrainFactory::CreateTerrain : invalid TerrainType");
            break;
    }
    return nullptr;
}

//
// ExplicitEulerIntegrator
//

IntegratorType ExplicitEulerIntegrator::getType() { return IntegratorType::ExplicitEuler; }

void ExplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-1: Integrate position and velocity
    //   1. Integrate position using current velocity.
    //   2. Integrate velocity using current acceleration.
    //   3. Clear force
    // Note:
    //   1. You should do this first because it is very simple. Then you can check whether your collision is correct or
    //   not.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Review "ODE_basics.pptx" from p.15 - p.16
    Jelly* jelly = particleSystem.getJellyPointer(0);
    float dt = particleSystem.deltaTime;
    
    for (int i = 0; i < jelly->getParticleNum(); i++) {
        Particle& particle = jelly->getParticle(i);

        Eigen::Vector3f x = particle.getPosition();
        Eigen::Vector3f v = particle.getVelocity();
        //x(t+dt) = x(t) + v(t) * dt
        x += v * dt;
        particle.setPosition(x);

        Eigen::Vector3f f = particle.getForce();
        float m = particle.getMass();
        // v(t+dt) = v(t) + a(t) * dt
        v += (f / m) * dt;
        particle.setVelocity(v);

        //clear force
        particle.setForce(Eigen::Vector3f::Zero());
    }

    // if it's the elevator scene, update elevator's position
    if (particleSystem.sceneIdx == 1) {
        Eigen::Vector3f x = particleSystem.elevatorTerrain->getPosition();
        Eigen::Vector3f v = particleSystem.elevatorTerrain->getVelocity();
        Eigen::Vector3f f = particleSystem.elevatorTerrain->getForce();
        float m = particleSystem.elevatorTerrain->getMass();

        // v(t+dt) = v(t) + a(t) * dt
        v += (f / m) * dt;
        particleSystem.elevatorTerrain->setVelocity(v);

        // x(t+dt) = x(t) + v(t) * dt
        x += v * dt;
        particleSystem.elevatorTerrain->setPosition(x);

        // 清除力
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    }
}

//
// ImplicitEulerIntegrator
//

IntegratorType ImplicitEulerIntegrator::getType() { return IntegratorType::ImplicitEuler; }

void ImplicitEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-2: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 and Vn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_implicit.pptx" from p.18 - p.19

    Jelly* jelly = particleSystem.getJellyPointer(0);
    float dt = particleSystem.deltaTime;

    std::vector<Eigen::Vector3f> x;
    std::vector<Eigen::Vector3f> v;
    std::vector<Eigen::Vector3f> f;
    int n = jelly->getParticleNum();

    x.reserve(n);
    v.reserve(n);
    f.reserve(n);

    // store original states of particles
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        x.push_back(p.getPosition());
        v.push_back(p.getVelocity());
        f.push_back(p.getForce());
    }

    // perform explicit euler method first
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        Eigen::Vector3f xi = x[i];
        Eigen::Vector3f vi = v[i];
        Eigen::Vector3f fi = f[i];
        float m = p.getMass();

        xi += vi * dt;
        vi += (fi / m) * dt;

        p.setPosition(xi);
        p.setVelocity(vi);
        p.setForce(Eigen::Vector3f::Zero());
    }

    int elevator = particleSystem.elevatorCounter;
    Eigen::Vector3f ex, ev, ef, ex1, ev1, ef1;
    if (particleSystem.sceneIdx == 1) {
        ex = particleSystem.elevatorTerrain->getPosition();
        ev = particleSystem.elevatorTerrain->getVelocity();
        ef = particleSystem.elevatorTerrain->getForce();
        float em = particleSystem.elevatorTerrain->getMass();

        // explicit euler method
        ex1 = ex + ev * dt;
        ev1 = ev + (ef / em) * dt;
        particleSystem.elevatorTerrain->setPosition(ex1);
        particleSystem.elevatorTerrain->setVelocity(ev1);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    }

    // jelly force at t=t+delta_t
    if (particleSystem.sceneIdx == 1) particleSystem.computeElevatorForce();
    particleSystem.computeJellyForce(*jelly);

    // perform implicit method
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        Eigen::Vector3f f1 = p.getForce();
        float m = p.getMass();

        Eigen::Vector3f v1 = v[i] + (f1 / m) * dt;  // vi+1 = vi + ai+1 * t
        Eigen::Vector3f x1 = x[i] + v1 * dt;        // xi+1 = xi + vi+1 * t

        p.setPosition(x1);
        p.setVelocity(v1);
        p.setForce(Eigen::Vector3f::Zero());
    }

    if (particleSystem.sceneIdx == 1) {
        float em = particleSystem.elevatorTerrain->getMass();
        Eigen::Vector3f ef1 = particleSystem.elevatorTerrain->getForce();
        Eigen::Vector3f ev1 = ev + (ef1 / em) * dt;
        Eigen::Vector3f ex1 = ex + ev1 * dt;
        particleSystem.elevatorTerrain->setPosition(ex1);
        particleSystem.elevatorTerrain->setVelocity(ev1);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
        particleSystem.elevatorCounter = elevator;
    }
}

//
// MidpointEulerIntegrator
//

IntegratorType MidpointEulerIntegrator::getType() { return IntegratorType::MidpointEuler; }

void MidpointEulerIntegrator::integrate(MassSpringSystem& particleSystem) {
    // TODO#4-3: Integrate position and velocity
    //   1. Backup original particles' data.
    //   2. Integrate position and velocity using explicit euler to get Xn+1.
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. Review "ODE_basics.pptx" from p.18 - p.19
    Jelly* jelly = particleSystem.getJellyPointer(0);
    float dt = particleSystem.deltaTime;

    std::vector<Eigen::Vector3f> x;
    std::vector<Eigen::Vector3f> v;
    std::vector<Eigen::Vector3f> f;
    int n = jelly->getParticleNum();

    x.reserve(n);
    v.reserve(n);
    f.reserve(n);

    // store original states of particles
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        x.push_back(p.getPosition());
        v.push_back(p.getVelocity());
        f.push_back(p.getForce());
    }

    float halft = dt / 2.0f;
    //perform explicit euler method in t + delta t/2
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        Eigen::Vector3f xi = x[i];
        Eigen::Vector3f vi = v[i];
        Eigen::Vector3f fi = f[i];
        float m = p.getMass();

        xi += vi * halft;
        vi += (fi / m) * halft;

        p.setPosition(xi);
        p.setVelocity(vi);
        p.setForce(Eigen::Vector3f::Zero());
    }

    int elevator = particleSystem.elevatorCounter;
    Eigen::Vector3f ex, ev, ef, ex1, ev1, ef1;
    if (particleSystem.sceneIdx == 1) {
        ex = particleSystem.elevatorTerrain->getPosition();
        ev = particleSystem.elevatorTerrain->getVelocity();
        ef = particleSystem.elevatorTerrain->getForce();
        float em = particleSystem.elevatorTerrain->getMass();

        // explicit euler method
        ex1 = ex + ev * halft;
        ev1 = ev + (ef / em) * halft;
        particleSystem.elevatorTerrain->setPosition(ex1);
        particleSystem.elevatorTerrain->setVelocity(ev1);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
    }
    // jelly force at t=t+dt/2
    if (particleSystem.sceneIdx == 1) particleSystem.computeElevatorForce();
    particleSystem.computeJellyForce(*jelly);

    // perform midpoint method
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        Eigen::Vector3f f_mid = p.getForce();
        Eigen::Vector3f v_mid = p.getVelocity();
        float m = p.getMass();

        Eigen::Vector3f x1 = x[i] + v_mid * dt;
        Eigen::Vector3f v1 = v[i] + (f_mid / m) * dt;  
        //x(t+dt) = x(t) + v(t+dt/2) * dt
       // v(t + dt) = v(t) + a(t + dt / 2) * dt
        
        p.setPosition(x1);
        p.setVelocity(v1);
        p.setForce(Eigen::Vector3f::Zero());
    }

    if (particleSystem.sceneIdx == 1) {
        float em = particleSystem.elevatorTerrain->getMass();
        Eigen::Vector3f ef1 = particleSystem.elevatorTerrain->getForce();
        ex1 = ex + ev1 * dt;
        ev1 = ev + (ef1 / em) * dt;
        particleSystem.elevatorTerrain->setPosition(ex1);
        particleSystem.elevatorTerrain->setVelocity(ev1);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
        particleSystem.elevatorCounter = elevator;
    }
}

//
// RungeKuttaFourthIntegrator
//

IntegratorType RungeKuttaFourthIntegrator::getType() { return IntegratorType::RungeKuttaFourth; }

void RungeKuttaFourthIntegrator::integrate(MassSpringSystem& particleSystem) {
    struct StateStep {
        Eigen::Vector3f deltaVel;
        Eigen::Vector3f deltaPos;
    };
    // TODO#4-4: Integrate velocity and acceleration
    //   1. Backup original particles' data.
    //   2. Compute k1, k2, k3, k4
    //   3. Compute refined Xn+1 using (1.) and (2.).
    // Note:
    //   1. Use `MassSpringSystem::computeJellyForce`and `MassSpringSystem::computeElevatorForce`
    //      with modified position and velocity to get Xn+1.
    //   2. See functions in `particle.h` if you don't know how to access data.
    //   3. You can access data and function in particleSystem.elevatorTerrain
    //   4. Remember to reset particleSystem.elevatorCounter back to the original state
    //   5. StateStep struct is just a hint, you can use whatever you want.
    //   6. Review "ODE_basics.pptx" from p.21

    Jelly* jelly = particleSystem.getJellyPointer(0);
    float dt = particleSystem.deltaTime;

    std::vector<Eigen::Vector3f> x;
    std::vector<Eigen::Vector3f> v;
    int n = jelly->getParticleNum();

    x.reserve(n);
    v.reserve(n);

    // store original states of particles
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        x.push_back(p.getPosition());
        v.push_back(p.getVelocity());
        p.setForce(Eigen::Vector3f::Zero());
    }

    int eleCounter = particleSystem.elevatorCounter;
    Eigen::Vector3f ex = particleSystem.elevatorTerrain->getPosition();
    Eigen::Vector3f ev = particleSystem.elevatorTerrain->getVelocity();

    // k value for each particle
    std::vector<std::vector<StateStep> > k(4, std::vector<StateStep>(n));

    for (int i = 0; i < n; i++) {
        jelly->getParticle(i).setForce(Eigen::Vector3f::Zero());
        // reset the force to aviod force accumulation
    }

    if (particleSystem.sceneIdx == 1) {
        particleSystem.computeElevatorForce();
    }
    particleSystem.computeAllForce();

    // k1 = hf(x₀, t₀)
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        k[0][i].deltaVel = p.getForce() / p.getMass() * dt;
        // k1.v = h * a = h * F/m
        k[0][i].deltaPos = p.getVelocity() * dt;
        // k1.x = h * v
    }

    // k2 = hf(x₀ + k₁/2, t₀ + h/2)
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        p.setPosition(x[i] + k[0][i].deltaPos * 0.5f);
        p.setVelocity(v[i] + k[0][i].deltaVel * 0.5f);
        p.setForce(Eigen::Vector3f::Zero());
    }

    if (particleSystem.sceneIdx == 1) {
        particleSystem.computeElevatorForce();
    }
    particleSystem.computeAllForce();

    // k2 = hf(x₀ + k₁/2, t₀ + h/2)
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        k[1][i].deltaVel = p.getForce() / p.getMass() * dt;
        k[1][i].deltaPos = p.getVelocity() * dt;
    }

    // k3 = hf(x₀ + k₂/2, t₀ + h/2)
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        p.setPosition(x[i] + k[1][i].deltaPos * 0.5);
        p.setVelocity(v[i] + k[1][i].deltaVel * 0.5);
        p.setForce(Eigen::Vector3f::Zero());
    }

    if (particleSystem.sceneIdx == 1) {
        particleSystem.computeElevatorForce();
    }
    particleSystem.computeAllForce();

    // k3 = hf(x₀ + k₂/2, t₀ + h/2)
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        k[2][i].deltaVel = p.getForce() / p.getMass() * dt;
        k[2][i].deltaPos = p.getVelocity() * dt;
    }

    // k4 = hf(x₀ + k₃, t₀ + h)
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        p.setPosition(x[i] + k[2][i].deltaPos);
        p.setVelocity(v[i] + k[2][i].deltaVel);
        p.setForce(Eigen::Vector3f::Zero());
    }
    if (particleSystem.sceneIdx == 1) {
        particleSystem.computeElevatorForce();
    }
    particleSystem.computeAllForce();

    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);
        k[3][i].deltaVel = p.getForce() / p.getMass() * dt;
        k[3][i].deltaPos = p.getVelocity() * dt;
    }

    //(k1 + 2*k2 + 2*k3 + k4)/6
    for (int i = 0; i < n; i++) {
        Particle& p = jelly->getParticle(i);

        Eigen::Vector3f dx =
            (k[0][i].deltaPos + 2.0f * k[1][i].deltaPos + 2.0f * k[2][i].deltaPos + k[3][i].deltaPos) / 6.0f;
        Eigen::Vector3f dv =
            (k[0][i].deltaVel + 2.0f * k[1][i].deltaVel + 2.0f * k[2][i].deltaVel + k[3][i].deltaVel) / 6.0f;

        p.setPosition(x[i] + dx);
        p.setVelocity(v[i] + dv);
        p.setForce(Eigen::Vector3f::Zero());
    }

    //apply RK4 for elevator
    if (particleSystem.sceneIdx == 1) {
        particleSystem.elevatorCounter = eleCounter;
        particleSystem.elevatorTerrain->setPosition(ex);
        particleSystem.elevatorTerrain->setVelocity(ev);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
        particleSystem.computeElevatorForce();

        // Store k values for elevator
        std::vector<Eigen::Vector3f> exk(4);
        std::vector<Eigen::Vector3f> evk(4);

        // k1 
        Eigen::Vector3f ef = particleSystem.elevatorTerrain->getForce();
        float em = particleSystem.elevatorTerrain->getMass();
        evk[0] = ef / em * dt;
        exk[0] = ev * dt;

        // k2
        particleSystem.elevatorCounter = eleCounter;
        particleSystem.elevatorTerrain->setPosition(ex + exk[0] * 0.5f);
        particleSystem.elevatorTerrain->setVelocity(ev + evk[0] * 0.5f);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
        particleSystem.computeElevatorForce();
        ef = particleSystem.elevatorTerrain->getForce();
        evk[1] = ef / em * dt;
        exk[1] = (ev + evk[0] * 0.5f) * dt;

        // k3
        particleSystem.elevatorCounter = eleCounter;
        particleSystem.elevatorTerrain->setPosition(ex + exk[1]*0.5f);
        particleSystem.elevatorTerrain->setVelocity(ev + evk[1]*0.5f);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
        particleSystem.computeElevatorForce();
        ef = particleSystem.elevatorTerrain->getForce();
        evk[2] = ef / em* dt;
        exk[2] = (ev + evk[1] * 0.5f) * dt;

        // k4 already calculated in fourth stage
        particleSystem.elevatorCounter = eleCounter;
        particleSystem.elevatorTerrain->setPosition(ex + exk[2]);
        particleSystem.elevatorTerrain->setVelocity(ev + evk[2]);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());
        particleSystem.computeElevatorForce();
        ef = particleSystem.elevatorTerrain->getForce();
        evk[3] = ef / em * dt;
        exk[3] = (ev + evk[2]) * dt;

        // Final update
        Eigen::Vector3f dx =
            (exk[0] + 2.0f * exk[1] + 2.0f * exk[2] + exk[3]) / 6.0f;

        Eigen::Vector3f dv =
            (evk[0] + 2.0f * evk[1] + 2.0f * evk[2] + evk[3]) / 6.0f;

        particleSystem.elevatorTerrain->setPosition(ex + dx);
        particleSystem.elevatorTerrain->setVelocity(ev + dv);
        particleSystem.elevatorTerrain->setForce(Eigen::Vector3f::Zero());

        particleSystem.elevatorCounter = eleCounter;
    }
}
}  // namespace simulation
