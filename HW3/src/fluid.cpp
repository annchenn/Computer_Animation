#include "fluid.h"
#include <iostream>
void Fluid::integrate()
{
    Eigen::Vector2f gravityVec(0, -gravity * dt);

    for (int i = 0; i < num_particles; ++i)
    {
        particle_vel[i] += gravityVec;
        particle_pos[i] += dt * particle_vel[i];
    }
}

void Fluid::handleCollisions(bool lmbDown)
{
    Eigen::Vector2f minVal(cell_dim + particle_radius, cell_dim + particle_radius),
        maxVal(viewport_w - cell_dim - particle_radius, viewport_h - cell_dim - particle_radius);
    
    float min_obstacle_d2 = (obstacle_r + particle_radius) * (obstacle_r + particle_radius);
    Eigen::Vector2f obstacle_pos(obstacle_x, obstacle_y), obstacle_vel(obstacle_vx, obstacle_vy);

    for (int i = 0; i < num_particles; ++i)
    {
        Eigen::Vector2f& pos = particle_pos[i];
        Eigen::Vector2f& vel = particle_vel[i];

        // collision with obstacle
        if (lmbDown)
        {
            float d2_to_obstacle = (pos - obstacle_pos).squaredNorm();
            if (d2_to_obstacle < min_obstacle_d2)
                vel = obstacle_vel;
        }

        // collision with edges
        Eigen::Vector2f clamped = pos.cwiseMax(minVal).cwiseMin(maxVal);

        if (clamped.x() != pos.x()) vel.x() = 0.0f;
        if (clamped.y() != pos.y()) vel.y() = 0.0f;

        pos = clamped;
    }
}

void Fluid::particleRelaxation()
{
    /*TODO*/
    /*
    This function handles inter-particle collisions, here's what you'll have to implement:


    1.	Load all particles indices into relaxation cells

    2.	Implement the iterative solver:

        At each iteration, go through all relaxation cells.
        For each particle in the cell, compare its position to every other particle within a 3x3 region.
        If two particles are colliding, i.e. |d| < 2 * r, shift them away from each other so the two particles are 2 * r units apart.
        [DEBUGGING HINT: Two particles may completely overlap with each other!]

        You can directly apply the displacement vectors right after their calculation.


    Some of the variables you might need, but not limited to:

    Relaxation Cells:
        relaxation_cell_col				// Number of relaxation cell columns
        relaxation_cell_rows			// Number of relaxation cell rows
        relaxation_cell_particle_ids	// Index of particles in each cell
        relaxation_cell_dim				// Dimension of each relaxation cell (same width and height)


    Particle Information:
        particle_pos					// Position (x, y) of each particle
        particle_radius	
        num_particles
    */

    // TODO: Assign particles to cells

    //clear the relaxation cell list
    for (int i = 0; i < relaxation_cell_rows; i++) {
        for (int j = 0; j < relaxation_cell_cols; j++) {
            relaxation_cell_particle_ids[i][j].clear();
        }
    }

    //assign particle to cells
    for (int i = 0; i < num_particles; i++) {
        int row = static_cast<int>(particle_pos[i].y() / relaxation_cell_dim);
        int col = static_cast<int>(particle_pos[i].x() / relaxation_cell_dim);

        row = std::max(0, std::min(row, relaxation_cell_rows - 1));
        col = std::max(0, std::min(col, relaxation_cell_cols - 1));
        relaxation_cell_particle_ids[row][col].push_back(i);
    }
    

    //perform relaxation
    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int i = 0; i < relaxation_cell_rows; ++i) for (int j = 0; j < relaxation_cell_cols; ++j)
        {
            // TODO: Perform particle relaxation
            for (size_t p = 0; p < relaxation_cell_particle_ids[i][j].size(); p++) {
                int pid = relaxation_cell_particle_ids[i][j][p];
                Eigen::Vector2f pos_p = particle_pos[pid];

                for (int ni = std::max(0, i - 1); ni <= std::min(relaxation_cell_rows - 1, i + 1); ni++) {
                    for (int nj = std::max(0, j - 1); nj <= std::min(relaxation_cell_cols - 1, j + 1);nj++) {
                        for (size_t q = 0; q < relaxation_cell_particle_ids[ni][nj].size(); q++) {
                            int qid = relaxation_cell_particle_ids[ni][nj][q];
                            if (pid == qid) continue;
                            Eigen::Vector2f pos_q = particle_pos[qid];
                            Eigen::Vector2f d = pos_p - pos_q;

                            float d2 = d.squaredNorm();
                            float min_d = 2.0f * particle_radius;
                            if (d2 < min_d * min_d&&d2>0.0f) {
                                float dist = std::sqrt(d2);
                                //(1/2) * (2r - |d|) * (d/|d|)
                                Eigen::RowVector2f corr = 0.5f * (min_d - dist) * (d / dist);

                                particle_pos[pid] += corr;
                                particle_pos[qid] -= corr;
                            }
                        }
                    }
                }
            }
        }
    }
}

void Fluid::transferVelocities(bool to_cell)
{
    /*TODO*/
    /*
    This function has a parameter to_cell that determines whether the particle velocities should be transfered to the cells 
    or in the other direction.

    1.	Some setup you need to do when to_cell is true:

        Reset all cell_velocities.
        Assign each cell their proper cell type. A cell should be assigned the FLUID type if it contains any particles and AIR if not.
        SOLID cells are present on initialization and should not be modified.

    2.	In main loop:

        For each of x and y components, calculate the 4 cells that should contribute to the biliear interpolation of each particle.
        From this step, we can obtain 4 weights, w_1 ~ w_4. These weights will be used in both scenarios (to_cell = true or false).

        When transfering particle velocities to the cells, add w_n * v to each of the 4 cell velocity components. 

        We should also accumulate w_n of each sampling point in a separate buffer to normalize the cell velocities after 
        iterating through all the particles.

        When transfering cell velocities back to the particles, before doing the usual bilinear interpolation, we need to test if 
        the velocities are valid.
        
        Let's take a look at an example, where we calculate the x velocity component
    
                |		|
        -------------------------
                |		|			vx_1: Vx(i, j)
            F	•->	A	•->	A		vx_2: Vx(i + 1, j)
           Vx_4	|   	| Vx_3		vx_3: Vx(i + 1, j + 1)
        -------------------------	vx_4: Vx(i, j + 1)
                |	  ⦾	|			
            F	•->	F <-•	A		F: Fluid Cells
           Vx_1	|		| Vx_2		A: Air Cells
        -------------------------
                |		|

        In this case, only vx_3 is invalid since the value it carries is derived from 2 AIR cells.
        Note, since the Vy is sampled at the center of the horizontal edge, we should test the upper and lower cells in that case.

        After this we can finally obtain the interpolated velocity. In the provided example, since Vx_3 is invalid, the interpolated velocity
        will be:

                    w_1 * vx_1 + w_2 * vx_2 + w_4 * vx_4
         Vx_pic =	------------------------------------
                              w_1 + w_2 + w_4

        We're almost done!
        The particle velocity was calculated using the PIC method. To obtain FLIP velocity use this formula:

                            w_1 * (vx_1 - prev_v1) + w_2 * (vx_2 - prev_v2) + w_4 * (vx_4 - prev_v4)
         Vx_flip =	Vx_p +	-----------------------------------------------------------------------
                                                       w_1 + w_2 + w_4

        We can finally set the particle velocity using a blend of PIC and FLIP velocities:

         Vx_p = (1 - flip_ratio) * Vx_pic + flip_ratio * Vx_flip

        Don't forget to do this in the Y direction as well!

    3.	After the loop, if we're transfering cell velocities to the particle, we'll have to normalize each velocity component using 
        the stored normalizing values mentioned before.

        Backup cell velocities in prev_cell_velocities for FLIP calculation.


    Some of the variables you might need, but not limited to:

    flip_ratio					# used during pic/flip interpolation

    Particle Information:
        particle_pos
        particle_vel
        num_particles

    MAC Cells:
        cell_rows
        cell_cols
        cell_dim
        cell_types				# Type of each cell. {CellType::FLUID, CellType::AIR, CellType::SOLID}
        cell_velocities			# Vx (sampled at the left edge of the cell)
                                # Vy (sampled at the bottom edge of the cell)

        prev_cell_velocities

    A 2d vector to store the normalizing terms for cell velicities when to_cell is true

    */


    if (to_cell)
    {
        // TODO: Reset velocities and update cell types
        for (int i = 0; i < cell_rows; i++) for (int j = 0; j < cell_cols; j++) {
            cell_velocities[i][j] = Eigen::Vector2f::Zero();
            //default to be air
            if (cell_types[i][j] != CellType::SOLID)
                cell_types[i][j] = CellType::AIR;
        }
    }
    //array to store w_accumulation
    std::vector<std::vector<float>> weight_sum_x(cell_rows, std::vector<float>(cell_cols+1, 0.0f));
    std::vector<std::vector<float>> weight_sum_y(cell_rows, std::vector<float>(cell_cols+1, 0.0f));

    
    // If you find iterating over x and y components unintuitive, feel free to change the structure of the code
    for (int component = 0; component < 2; ++component)
    {
        // Might need some setup here

        for (int i = 0; i < num_particles; ++i)
        {
            // TODO: Bilinear interpolation on staggered grid
            Eigen::Vector2f pos = particle_pos[i];

            float gx, gy;
            if (component == 0) { // Vx component (located at (j*cell_dim, (i+0.5)*cell_dim))
                gx = pos.x() / cell_dim;
                gy = (pos.y() - 0.5f * cell_dim) / cell_dim;
            }
            else { // Vy component (located at ((j+0.5)*cell_dim, i*cell_dim))
                gx = (pos.x() - 0.5f * cell_dim) / cell_dim;
                gy = pos.y() / cell_dim;
            }

            int j_lower = static_cast<int>(gx);
            int i_lower = static_cast<int>(gy);

            float dx = gx - j_lower; // Fractional distance in x
            float dy = gy - i_lower; // Fractional distance in y

            float w1 = (1.0f - dx) * (1.0f - dy); 
            float w2 = dx * (1.0f - dy);
            float w3 = dx * dy;
            float w4 = (1.0f - dx) * dy;

            float weights[4] = { w1, w2, w3, w4 };

            Eigen::Vector2i v_indices[4] = {
                {i_lower, j_lower},         // Bottom-left grid node for this component's stencil
                {i_lower, j_lower + 1},     // Bottom-right
                {i_lower + 1, j_lower + 1}, // Top-right
                {i_lower + 1, j_lower}
            };
            if (to_cell)//粒子到網格
            {
                float v_p = (component == 0) ? particle_vel[i].x() : particle_vel[i].y();
                int p_cell_row = static_cast<int>(pos.y() / cell_dim);
                int p_cell_col = static_cast<int>(pos.x() / cell_dim);
                p_cell_row = std::max(0, std::min(p_cell_row, cell_rows - 1));
                p_cell_col = std::max(0, std::min(p_cell_col, cell_cols - 1));
                
                if (cell_types[p_cell_row][p_cell_col] != CellType::SOLID) {
                    cell_types[p_cell_row][p_cell_col]=CellType::FLUID;
                }
                for (int k = 0; k < 4; ++k) {
                    int r_idx = v_indices[k].x();
                    int c_idx = v_indices[k].y();

                    bool valid_access = false;
                    if (component == 0) { // Vx: r in [0, cell_rows-1], c in [0, cell_cols]
                        valid_access = (r_idx >= 0 && r_idx < cell_rows && c_idx >= 0 && c_idx <= cell_cols);
                    }
                    else { // Vy: r in [0, cell_rows], c in [0, cell_cols-1]
                        valid_access = (r_idx >= 0 && r_idx <= cell_rows && c_idx >= 0 && c_idx < cell_cols);
                    }

                    if (valid_access) {
                        cell_velocities[r_idx][c_idx][component] += weights[k] *v_p;
                        if (component == 0) weight_sum_x[r_idx][c_idx] += weights[k];
                        else weight_sum_y[r_idx][c_idx] += weights[k];
                    }
                }
            }
            else
            {
                // TODO: Transfer valid cell velocities back to the particles using a mixture of PIC and FLIP
                //網格到粒子
                float v_p = (component == 0) ? particle_vel[i].x() : particle_vel[i].y();
                float pic_v_sum = 0.0f;
                float flip_delta_v_sum = 0.0f;
                float total_valid_w = 0.0f;

                for (int k = 0; k < 4; k++) {
                    int r_idx = v_indices[k].x();
                    int c_idx = v_indices[k].y();

                    bool valid_access = false;
                    if (component == 0) { valid_access = (r_idx >= 0 && r_idx < cell_rows && c_idx >= 0 && c_idx <= cell_cols); }
                    else { valid_access = (r_idx >= 0 && r_idx <= cell_rows && c_idx >= 0 && c_idx < cell_cols); }

                    if (!valid_access) continue;

                    // Velocity sample validity (based on adjacent cell types)
                    bool is_v_sample_valid = false;
                    if (component == 0) { // Vx at (r_idx, c_idx) separates cell (r_idx, c_idx-1) and (r_idx, c_idx)
                        // Valid if c_idx is not a solid boundary itself and at least one adjacent cell is FLUID
                        if (c_idx > 0 && c_idx < cell_cols) { // Interior Vx
                            is_v_sample_valid = (cell_types[r_idx][c_idx - 1] == CellType::FLUID || cell_types[r_idx][c_idx] == CellType::FLUID);
                        }
                        else if (c_idx == 0 && cell_types[r_idx][0] != CellType::SOLID) { // Left boundary Vx (if not solid itself)
                            is_v_sample_valid = (cell_types[r_idx][0] == CellType::FLUID);
                        }
                        else if (c_idx == cell_cols && cell_types[r_idx][cell_cols - 1] != CellType::SOLID) { // Right boundary Vx
                            is_v_sample_valid = (cell_types[r_idx][cell_cols - 1] == CellType::FLUID);
                        }
                    }
                    else { // Vy at (r_idx, c_idx) separates cell (r_idx-1, c_idx) and (r_idx, c_idx)
                        if (r_idx > 0 && r_idx < cell_rows) { // Interior Vy
                            is_v_sample_valid = (cell_types[r_idx - 1][c_idx] == CellType::FLUID || cell_types[r_idx][c_idx] == CellType::FLUID);
                        }
                        else if (r_idx == 0 && cell_types[0][c_idx] != CellType::SOLID) { // Bottom boundary Vy
                            is_v_sample_valid = (cell_types[0][c_idx] == CellType::FLUID);
                        }
                        else if (r_idx == cell_rows && cell_types[cell_rows - 1][c_idx] != CellType::SOLID) { // Top boundary Vy
                            is_v_sample_valid = (cell_types[cell_rows - 1][c_idx] == CellType::FLUID);
                        }
                    }

                    if (is_v_sample_valid) {
                        pic_v_sum += weights[k] * cell_velocities[r_idx][c_idx][component];
                        flip_delta_v_sum += weights[k] * (cell_velocities[r_idx][c_idx][component] - prev_cell_velocities[r_idx][c_idx][component]);
                        total_valid_w += weights[k];
                    }
                }
                if (total_valid_w > 0.0f) {
                    float pic_v = pic_v_sum / total_valid_w;
                    float flip_v_change = flip_delta_v_sum / total_valid_w;
                    float flip_v = v_p + flip_v_change;

                    float final_v = (1.0f - flip_ratio) * pic_v + flip_ratio * flip_v;
                    if (component == 0) particle_vel[i].x() = final_v;
                    else particle_vel[i].y() = final_v;
                }
            }
        }

        if (to_cell)
        {
            // TODO: Normalize cell velocities and store a backup in prev_velocities.
            for (int i = 0; i < cell_rows; i++) for (int j = 0; j < cell_cols; j++) {
                float w = (component == 0) ? weight_sum_x[i][j] : weight_sum_y[i][j];
                if (w > 0.0f) {
                    cell_velocities[i][j][component] /= w;
                }


                prev_cell_velocities[i][j][component] = cell_velocities[i][j][component];
            }
        }
    }
}

void Fluid::updateDensity()
{
    /*TODO*/
    /*
    Here we will update the cell densities, which will be used to determine if a cell is overly compressed in the next step.
    The density of each cell is sampled at the center, and it uses the same concept as the transferVelocities function. 

    1.	Perform bilinear interpolation on each particle, each particle will contribute 1 * w_n to each of the sampling points.

    2.	Set the the resting density of water cells
        
        particle rest density = sum of all water cell densities / number of water cells

        This only needs to be done once, so simply test if particle_rest_density is 0 (the initial value).
        The particles are initialized to be tightly packed, so if you changed the layout of the particles, this step will have incorrect behavior.


    Some of the variables you might need, but not limited to:

    particle_rest_density		# updated only once at the first pass (test if particle_rest_density == 0)

    MAC cells:
        cell_densities
        cell_types
        cell_row
        cell_cols
        cell_dim

    Particle Information:
        num_particles
        particle_pos
    */

    for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
        cell_densities[i][j] = 0.0f;

    
    for (int i = 0; i < num_particles; ++i)
    {
        // TODO: Perform bilinear interpolation

        //the base cell of the particle
        float gx = (particle_pos[i].x() - 0.5f * cell_dim) / cell_dim;
        float gy = (particle_pos[i].y() - 0.5f * cell_dim) / cell_dim;

        int cellX = static_cast<int>(gx); // Column index of bottom-left cell center
        int cellY = static_cast<int>(gy); // Row index of bottom-left cell center

        float dx = gx - cellX;
        float dy = gy - cellY;

        float w1 = (1 - dx) * (1 - dy);
        float w2 = dx * (1 - dy); 
        float w3 = dx * dy;
        float w4 = (1 - dx) * dy;

        //add desity to surrounding cells
        if (cellY >= 0 && cellY < cell_rows && cellX >= 0 && cellX < cell_cols)
            cell_densities[cellY][cellX] += w1;

        if (cellY >= 0 && cellY < cell_rows && cellX + 1 >= 0 && cellX + 1 < cell_cols)
            cell_densities[cellY][cellX + 1] += w2;

        if (cellY + 1 >= 0 && cellY + 1 < cell_rows && cellX + 1 >= 0 && cellX + 1 < cell_cols)
            cell_densities[cellY + 1][cellX + 1] += w3;

        if (cellY + 1 >= 0 && cellY + 1 < cell_rows && cellX >= 0 && cellX < cell_cols)
            cell_densities[cellY + 1][cellX] += w4;

    }

    if (particle_rest_density == 0.0) 
    {
        // TODO: Calculate resting particle densities in fluid cells.
        float total_density = 0.0f;
        int numFluidCells = 0;

        for (int i = 0; i < cell_rows; i++) for (int j = 0; j < cell_cols; j++) {
            if (cell_types[i][j] == CellType::FLUID) {
                total_density += cell_densities[i][j];
                numFluidCells += 1;
            }
        }

        if (numFluidCells) {
            particle_rest_density = total_density / numFluidCells;
        }
        else {
            particle_rest_density = 1.0f;
        }
    }

}

void Fluid::solveIncompressibility()
{
    /*TODO*/
    /*
    We will implement an iterative solver to enforce incompressibility on the cells.
    Here's what you will need to do:

    For each iteration, iterate over all fluid cells and perform the following:

        Sum up the number of potential flow directions (number of non-solid neighbor cells).

        Calculate the divergence of a cell using cell velocities.

        Calculate how much the cell density differs from the rest density and if the cell is overly compressed,
        we treat it as having more inflow by adding a penalty of stiffness_coefficient * compression to the divergence 
        to promote the expansion of fluids within the cell.

        The divergence needs to be corrected using the available flow directions, so the correction magnitude will be:
            divergence / normalizing_term * over_relaxation
        
        Apply the correction to cell velocities to try and solve for divergence = 0.

    Some of the variables you might need, but not limited to:

    particle_rest_density		
    density_correction			# Boolean that controls whether density correction is applied
    stiffness_coefficient		# Constant you can set with sliders
    over_relaxation				# Constant you can set with sliders

    MAC Cells:
        cell_rows
        cell_cols
        cell_types
        cell_velocities
        cell_densities
        
    */

    // TODO: Store velocities
    std::vector<std::vector<Eigen::Vector2f> > origin_v = cell_velocities;

    for (int iter = 0; iter < iterations; ++iter)
    {
        for (int i = 1; i < cell_rows - 1; ++i) for (int j = 1; j < cell_cols - 1; ++j)
        {
            // TODO: Calculate divergence of fluid cells
            if (cell_types[i][j] != CellType::FLUID) continue;

            int flow_dirs = 0;
            float div = 0.0f;
            if (cell_types[i - 1][j] != CellType::SOLID) {//bottom
                div -= cell_velocities[i][j].y();
                flow_dirs++; 
            }
            if (cell_types[i + 1][j] != CellType::SOLID) {//top
                div += cell_velocities[i + 1][j].y();
                flow_dirs++; 
            }
            if (cell_types[i][j - 1] != CellType::SOLID) {//left
                flow_dirs++;
                div -= cell_velocities[i][j].x();
            }
            if (cell_types[i][j + 1] != CellType::SOLID) {//right
                div += cell_velocities[i][j + 1].x();
                flow_dirs++;
            }
            if (flow_dirs == 0) continue;

            // TODO: Add bias to ouflow if density_correction is true
            if (density_correction && (cell_densities[i][j] > particle_rest_density)) {
                div -= stiffness_coefficient * (cell_densities[i][j] - particle_rest_density);
            }
            float corr = div / flow_dirs * over_relaxation;
            // TODO: Correct the cell velocitiesif (cell_types[i-1][j] != CellType::SOLID)
            if (cell_types[i][j-1] != CellType::SOLID)
                cell_velocities[i][j].x() += corr;    // Left velocity

            if (cell_types[i][j+1] != CellType::SOLID)
                cell_velocities[i][j+1].x() -= corr;  // Right velocity

            if (cell_types[i-1][j] != CellType::SOLID)
                cell_velocities[i][j].y() += corr;    // Bottom velocity

            if (cell_types[i+1][j] != CellType::SOLID)
                cell_velocities[i+1][j].y() -= corr;  // Top velocity
        }
    }
}

Fluid::Fluid(int num_particles, float radius, float obstacle_radius, float flip_ratio, float cell_dim, float relaxation_cell_dim, int iterations, float viewport_w, 
    float viewport_h, float dt, float gravity, float stiffness, bool density_correction, float over_relaxation)
    : num_particles(num_particles), particle_radius(radius), cell_dim(cell_dim), relaxation_cell_dim(relaxation_cell_dim), viewport_w(viewport_w), viewport_h(viewport_h), 
    dt(dt), gravity(gravity), obstacle_x(viewport_w), obstacle_y(viewport_h), obstacle_vx(0.0f), obstacle_vy(0.0f), obstacle_r(obstacle_radius), iterations(iterations), 
    flip_ratio(flip_ratio), stiffness_coefficient(stiffness), density_correction(density_correction), over_relaxation(over_relaxation)
{
    cell_rows = static_cast<int>(viewport_h / cell_dim);
    cell_cols = static_cast<int>(viewport_w / cell_dim);

    relaxation_cell_rows = static_cast<int>(viewport_h / relaxation_cell_dim);
    relaxation_cell_cols = static_cast<int>(viewport_w / relaxation_cell_dim);

    relaxation_cell_particle_ids = std::vector<std::vector<std::vector<int>>>(relaxation_cell_rows, std::vector<std::vector<int>>(relaxation_cell_cols));
    cell_types = std::vector<std::vector<CellType>>(cell_rows, std::vector<CellType>(cell_cols, CellType::AIR));
    cell_velocities = std::vector<std::vector<Eigen::Vector2f>>(cell_rows, std::vector<Eigen::Vector2f>(cell_cols, Eigen::Vector2f::Zero()));
    prev_cell_velocities = std::vector<std::vector<Eigen::Vector2f>>(cell_rows, std::vector<Eigen::Vector2f>(cell_cols, Eigen::Vector2f::Zero()));
    cell_colors = std::vector<Eigen::Vector3f>(cell_rows * cell_cols, Eigen::Vector3f::Zero());
    cell_densities = std::vector<std::vector<float>>(cell_rows, std::vector<float>(cell_cols, 0.0f));
    cell_centers_rendering.reserve(cell_rows * cell_cols);

    particle_pos = std::vector<Eigen::Vector2f>(num_particles, Eigen::Vector2f::Zero());
    particle_vel = std::vector<Eigen::Vector2f>(num_particles, Eigen::Vector2f::Zero());
    particle_colors = std::vector<Eigen::Vector3f>(num_particles, Eigen::Vector3f::Ones());

    for (int i = 1; i < cell_rows; ++i)
    {
        cell_types[i][0] = CellType::SOLID;
        cell_types[i][cell_cols - 1] = CellType::SOLID;
    }

    for (int j = 0; j < cell_cols; ++j)
    {
        cell_types[0][j] = CellType::SOLID;
        cell_types[cell_rows - 1][j] = CellType::SOLID;
    }

    // You can modify the way particle positions are initialized,
    // but the layout of the particles should be the same since it is used for rest density calculation.
    float x = cell_dim * 10;
    float y = cell_dim * 10;
    float x_max = 0.5f * viewport_w;

    float dx = 2 * radius;
    float dy = sqrt(3.0) / 2.0 * dx;

    bool stagger = false;

    for (int i = 0; i < num_particles; i++) 
    {
        particle_pos[i] = Eigen::Vector2f(x, y);
        x += dx;
        if (x > x_max)
        {
            y += dy;
            x = stagger ? cell_dim * 10 + radius : cell_dim * 10;
            stagger = !stagger;
        }
    }
}

void Fluid::update(int render_option, bool lmbDown)
{

    // physics
    integrate();
    particleRelaxation();
    handleCollisions(lmbDown);
    transferVelocities(true);
    updateDensity();
    solveIncompressibility();
    transferVelocities(false);

    
    // graphics
    switch (render_option)
    {
    case 0: // cells only
        updateCellColors();
        updateCellColorBuffers();
        renderCells();
        break;
    case 1: // both
        updateCellColors();
        updateCellColorBuffers();
        renderCells();
    case 2: // particles only
        updateParticleColors();
        updateParticleBuffers();
        renderParticles();
        break;
    default:
        return;
    }
    
    renderObstacle();
}

void Fluid::setObstacle(float mouse_x, float mouse_y, float mouse_vx, float mouse_vy)
{
    obstacle_x = mouse_x;
    obstacle_y = viewport_h - mouse_y;
    obstacle_vx = mouse_vx;
    obstacle_vy = -mouse_vy;
}

void Fluid::updateParticleColors()
{/*
    for (int i = 0; i < num_particles; ++i)
    {
        float vel = particle_vel[i].norm();
        float mult = std::min(vel, 50.0f) / 50.0f;
        particle_colors[i] = Eigen::Vector3f(mult, mult, 1.0f);
    }
*/
    for (int i = 0; i < num_particles; ++i) {
        float vel = particle_vel[i].norm();
        float height = particle_pos[i].y() / viewport_h; // Normalized height

        // Blend between deep blue and light blue based on height
        Eigen::Vector3f waterColor = Eigen::Vector3f(0.0f, 0.3f + 0.7f * height, 0.8f + 0.2f * height);

        // Add white foam based on velocity
        float foamFactor = std::min(vel / 50.0f, 1.0f);
        particle_colors[i] = waterColor * (1.0f - foamFactor) + Eigen::Vector3f::Ones() * foamFactor;
    }
}

void Fluid::updateCellColors()
{
    float max_density_estimate = cell_dim / particle_radius;
    for (int i = 0; i < cell_rows; ++i) for (int j = 0; j < cell_cols; ++j)
    {
        int ind = i * cell_cols + j;
        if (cell_types[i][j] == CellType::SOLID)
            cell_colors[ind] = Eigen::Vector3f(0.5f, 0.5f, 0.5f);
        else if (cell_types[i][j] == CellType::FLUID)
        {
            float mult = std::min(cell_densities[i][j], max_density_estimate) / max_density_estimate;
            
            cell_colors[ind] = Eigen::Vector3f(mult, mult, 1.0f);
        }
        else
            cell_colors[ind] = Eigen::Vector3f::Zero();
    }
}

void Fluid::setupRendering(GLuint pointShaders, GLuint obstacleShaders)
{
    // particles
    point_shaders = pointShaders;
    obstacle_shaders = obstacleShaders;

    glGenVertexArrays(1, &particles_vao);
    glBindVertexArray(particles_vao);

    glGenBuffers(1, &particles_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, particles_vbo);
    glBufferData(GL_ARRAY_BUFFER, num_particles * sizeof(Eigen::Vector2f),
        particle_pos.data(), GL_DYNAMIC_DRAW);

    GLint particleAttrPositionLoc = glGetAttribLocation(point_shaders, "attrPosition");
    if (particleAttrPositionLoc != -1) {
        glVertexAttribPointer(particleAttrPositionLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
        glEnableVertexAttribArray(particleAttrPositionLoc);
    }

    glGenBuffers(1, &particle_colors_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, particle_colors_vbo);
    glBufferData(GL_ARRAY_BUFFER, num_particles * sizeof(Eigen::Vector3f),
        particle_colors.data(), GL_DYNAMIC_DRAW);

    GLint particleAttrColorLoc = glGetAttribLocation(point_shaders, "attrColor");
    if (particleAttrColorLoc != -1) {
        glVertexAttribPointer(particleAttrColorLoc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
        glEnableVertexAttribArray(particleAttrColorLoc);
    }

    glBindVertexArray(0);

    // cells
    int num_cells = cell_rows * cell_cols;
    glGenVertexArrays(1, &cells_vao);
    glBindVertexArray(cells_vao);

    glGenBuffers(1, &cells_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, cells_vbo);

    for (int i = 0; i < cell_rows; ++i)
        for (int j = 0; j < cell_cols; ++j)
            cell_centers_rendering[i * cell_cols + j] = Eigen::Vector2f(j + 0.5, i + 0.5) * cell_dim;

    glBufferData(GL_ARRAY_BUFFER, num_cells * sizeof(Eigen::Vector2f),
        cell_centers_rendering.data(), GL_DYNAMIC_DRAW);

    GLint cellAttrPositionLoc = glGetAttribLocation(point_shaders, "attrPosition");

    glVertexAttribPointer(cellAttrPositionLoc, 2, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(cellAttrPositionLoc);

    glGenBuffers(1, &cell_colors_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, cell_colors_vbo);

    glBufferData(GL_ARRAY_BUFFER, num_cells * sizeof(Eigen::Vector3f),
        cell_colors.data(), GL_DYNAMIC_DRAW);

    GLint cellAttrColorLoc = glGetAttribLocation(point_shaders, "attrColor");

    glVertexAttribPointer(cellAttrColorLoc, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
    glEnableVertexAttribArray(cellAttrColorLoc);

    glBindVertexArray(0);

    // obstacle dummy vao
    glGenVertexArrays(1, &obstacle_vao);
}

void Fluid::updateParticleBuffers()
{
    glBindBuffer(GL_ARRAY_BUFFER, particles_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, num_particles * sizeof(Eigen::Vector2f), particle_pos.data());

    glBindBuffer(GL_ARRAY_BUFFER, particle_colors_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, num_particles * sizeof(Eigen::Vector3f), particle_colors.data());
}

void Fluid::renderParticles()
{
    glUseProgram(point_shaders);

    glUniform2f(glGetUniformLocation(point_shaders, "domainSize"), viewport_w, viewport_h);
    glUniform1f(glGetUniformLocation(point_shaders, "pointSize"), 2.0f * particle_radius);
    glUniform1f(glGetUniformLocation(point_shaders, "drawDisk"), 1.0f);

    glBindVertexArray(particles_vao);
    glDrawArrays(GL_POINTS, 0, num_particles);

    glBindVertexArray(0);
}

void Fluid::updateCellColorBuffers()
{
    glBindBuffer(GL_ARRAY_BUFFER, cell_colors_vbo);
    glBufferSubData(GL_ARRAY_BUFFER, 0, cell_rows * cell_cols * sizeof(Eigen::Vector3f), cell_colors.data());

}
void Fluid::renderCells()
{
    glUseProgram(point_shaders);

    glUniform2f(glGetUniformLocation(point_shaders, "domainSize"), viewport_w, viewport_h);
    glUniform1f(glGetUniformLocation(point_shaders, "pointSize"), cell_dim);
    glUniform1f(glGetUniformLocation(point_shaders, "drawDisk"), 0.0f);

    glBindVertexArray(cells_vao);
    glDrawArrays(GL_POINTS, 0, cell_rows * cell_cols);

    glBindVertexArray(0);
}

void Fluid::renderObstacle()
{
    glUseProgram(obstacle_shaders);

    glBindVertexArray(obstacle_vao);

    glUniform2f(glGetUniformLocation(obstacle_shaders, "domainSize"), viewport_w, viewport_h);
    glUniform2f(glGetUniformLocation(obstacle_shaders, "attrPosition"), obstacle_x, obstacle_y);
    glUniform1f(glGetUniformLocation(obstacle_shaders, "pointSize"), 2 * obstacle_r);
    
    glDrawArrays(GL_POINTS, 0, 1);
    glBindVertexArray(0);
}



