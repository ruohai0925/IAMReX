/*
 * SPDX-FileCopyrightText: 2023 - 2026 Yadong Zeng<zdsjtu@gmail.com> & ZhuXu Li<1246206018@qq.com>
 *
 * SPDX-License-Identifier: BSD-3-Clause
 */

#ifndef COLLISION_H
#define COLLISION_H

#include <AMReX_REAL.H>
#include <AMReX_RealVect.H>
#include <array>
#include <list>

//collision object type
enum CollisionType{
    COLLISION_PARTICLE = 0,
    COLLISION_PLANE,
};

//particle
typedef struct{
    CollisionType type{COLLISION_PARTICLE};
    amrex::RealVect location{0., 0., 0.};
    amrex::RealVect velocity{0., 0., 0.};
    amrex::RealVect preForece{0., 0., 0.};
    amrex::RealVect fc{0., 0., 0.};
    amrex::Real radius{0.};
    amrex::Real rho{0.};
}CollisionParticle;

// mesh cell to store particle
typedef struct{
    bool hasParticle{false};
    std::list<CollisionParticle*> collectParticle;
}CollisionCell;

// collision paire to calculate
typedef struct{
    CollisionParticle* solid1;
    CollisionParticle* solid2;
}CollisionPair;

/**
 * @defgroup collision Particle collision manager for DEM system
 *
 * Particle collision class
 */

/// @{
class ParticleCollision{
public:
    using Real = amrex::Real;
    using RealVect = amrex::RealVect;
    /**
     * Set the collision model to be used for particle interactions.
     * @param model Integer identifier for the collision model
     */
    void takeModel(int model);

    /**
     * Apply the Discrete Element Method (DEM) collision model.
     * This function implements the DEM approach for particle-particle collisions.
     */
    void DKTModel();

    explicit ParticleCollision() = default;

    /**
     * Generate background mesh for collision detection.
     * @param gm_lo Lower bounds of the computational domain
     * @param gm_hi Upper bounds of the computational domain
     * @param size Particle size parameter
     * @param l Grid spacing for collision mesh
     */
    void SetGeometry(RealVect gm_lo, RealVect gm_hi, Real size, Real l);

    /**
     * Insert a particle into the collision detection system.
     * @param location Position of the particle
     * @param velocity Velocity of the particle
     * @param radius Radius of the particle
     * @param rho Density of the particle
     */
    void InsertParticle(RealVect location, RealVect velocity, Real radius, Real rho);

    /**
     * Insert a collision pair into the collision collector.
     * @param p1 Pointer to first particle in collision pair
     * @param p2 Pointer to second particle in collision pair
     */
    void InsertCollision(CollisionParticle* p1, CollisionParticle* p2);

    /**
     * Generate all possible collision pairs between particles.
     * This function identifies which particles are in close proximity
     * and may collide during the current timestep.
     */
    void GenerateCollisionPairs();

    /**
     * Resolve all collision pairs using the selected collision model.
     * This function applies collision forces and updates particle states
     * based on the detected collisions.
     */
    void ResolveCollisionPairs();
    //all the Geometry split mesh
    std::vector<CollisionCell> *Cells;
    std::list<CollisionParticle> Particles;
    //mesh x,y,z cell size, loop for i,j,k
    int Nx, Ny, Nz;
    Real cell_size, mesh_size;
    //lo,hi
    RealVect lo, hi;
    //store collision type
    std::list<CollisionPair> CollisionCollector;

};

///@}

#endif
