#include "simulation/kinematics.h"

#include <iostream>
#include "Eigen/Dense"
#include "acclaim/bone.h"
#include "util/helper.h"

namespace kinematics {

void dfs(acclaim::Bone* current, acclaim::Posture& posture) {
    // Bone: point to the current bone
    // parent end: the point parents bone end
    // posture:　all angles for all bones
    if (!current) {
        return;
    }

    if (current->idx == 0) {
        // root bone
        current->start_position = posture.bone_translations[current->idx];

        Eigen::Quaterniond local_rotation = util::rotateDegreeZYX(posture.bone_rotations[current->idx]);
        current->rotation = local_rotation * current->rot_parent_current;
        current->end_position = current->start_position + current->rotation * (current->dir * current->length);
    } else {
        current->start_position = current->parent->end_position;

        Eigen::Quaternion local_rotation = util::rotateDegreeZYX(posture.bone_rotations[current->idx]);
        current->rotation = current->parent->rotation * current->rot_parent_current * local_rotation;
        current->end_position = current->start_position + current->rotation * (current->dir * current->length);
    }

    if (current->child) {
        dfs(current->child, posture);
    }
    if (current->sibling) {
        dfs(current->sibling, posture);
    }
}

void forwardSolver(const acclaim::Posture& posture, acclaim::Bone* bone) {
    // TODO#1: Forward Kinematic
    // Hint:
    // - Traverse the skeleton tree from root to leaves.
    // - Compute each bone's global rotation and global position.
    // - Use local rotation (from posture) and bone hierarchy (parent rotation, offset, etc).
    // - Remember to update both bone->start_position and bone->end_position.
    // - Use bone->rotation to store global rotation (after combining parent, local, etc).
    acclaim::Posture posture1 = posture;
    dfs(bone, posture1);
}

Eigen::VectorXd pseudoInverseLinearSolver(const Eigen::Matrix4Xd& Jacobian, const Eigen::Vector4d& target) {
    Eigen::VectorXd deltatheta;
    // TODO#2: Inverse linear solver (find x which min(| jacobian * x - target |))
    // Hint:
    //   1. Linear algebra - least squares solution
    //   2. https://en.wikipedia.org/wiki/Moore%E2%80%93Penrose_inverse#Construction
    // Note:
    //   1. SVD or other pseudo-inverse method is useful
    //   2. Some of them have some limitation, if you use that method you should check it.
    unsigned int compute_option = Eigen::ComputeThinU + Eigen::ComputeThinV;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(Jacobian, compute_option);

    double epsilon = 1e-6;
    deltatheta = svd.setThreshold(epsilon).solve(target);
    return deltatheta;
}

/**
 * @brief Perform inverse kinematics (IK)
 *
 * @param target_pos The position where `end_bone` will move to.
 * @param obs_pos The position where the obstacle is at
 * @param obsActive Whether the obstacle is active or not
 * @param start_bone This bone is the last bone you can move while doing IK
 * @param end_bone This bone will try to reach `target_pos`
 * @param posture The original AMC motion's reference, you need to modify this
 *
 * @return True if IK is stable (HW2 bonus)
 */
bool inverseJacobianIKSolver(const Eigen::Vector4d& target_pos, const Eigen::Vector4d& obs_pos, bool obsActive,
                             acclaim::Bone* start_bone, acclaim::Bone* end_bone, acclaim::Posture& posture) {
    constexpr int max_iteration = 1000;
    constexpr double epsilon = 1E-3;
    constexpr double step = 0.1;
    constexpr double obsAvoidThreshold = 1.01;  // if bone is within 1 unit from obstacle

    // Since bone stores in bones[i] that i == bone->idx, we can use bone - bone->idx to find bones[0] which is the
    // root.
    acclaim::Bone* root_bone = start_bone - start_bone->idx;

    // TODO#3:
    // Perform inverse kinematics (IK)
    // HINTs will tell you what should do in that area.
    // Of course you can ignore it (Any code below this line) and write your own code.
    acclaim::Posture original_posture(posture);

    size_t bone_num = 0;
    std::vector<acclaim::Bone*> boneList;
    acclaim::Bone* current = end_bone;
    // TODO#3-1:
    // Calculate number of bones need to move to perform IK, store in `bone_num`
    // (a.k.a. how may bones from end_bone to its parent than to start_bone (include both side))
    // Store the bones need to move to perform IK into boneList
    // Hint:
    //   1. Traverse from end_bone to start_bone is easier than start to end (since there is only 1 parent)
    //   2. If start bone is not reachable from end. Go to root first.
    // Note:
    //   1. Both start and end should be in the list
    boneList.push_back(current);
    // traverse from end bone to start bone
    while (current != start_bone && current != nullptr) {
        current = current->parent;
        if (current != nullptr) {
            boneList.push_back(current);
            if (current == start_bone) break;
        }
    }

    // if reach the nullptr without finding the start bone, go to the root
    if (current == nullptr || (current != start_bone && boneList.back() != start_bone)) {
        boneList.clear();

        // collect bones from end_bone to root
        current = end_bone;
        boneList.push_back(current);

        while (current->parent != nullptr) {
            current = current->parent;
            boneList.push_back(current);
        }

        // collect bones from start_bone to root
        std::vector<acclaim::Bone*> startToRoot;
        current = start_bone;
        startToRoot.push_back(current);

        while (current->parent != nullptr) {
            current = current->parent;
            startToRoot.push_back(current);
        }

        // Find the common ancestor
        int commonIndex = -1;
        for (size_t i = 0; i < boneList.size(); i++) {
            for (size_t j = 0; j < startToRoot.size(); j++) {
                if (boneList[i] == startToRoot[j]) {
                    commonIndex = i;
                    break;
                }
            }
            if (commonIndex != -1) break;
        }
        if (commonIndex != -1) {
            std::vector<acclaim::Bone*> newList;
            // Add bones from end_bone to common ancestor
            for (int i = 0; i <= commonIndex; i++) {
                newList.push_back(boneList[i]);
            }

            int start_to_root_idx = -1;
            for (size_t i = 0; i < startToRoot.size(); i++) {
                if (startToRoot[i] == boneList[commonIndex]) {
                    start_to_root_idx = i;
                    break;
                }
            }
            // Add bones from common ancestor to start_bone (in reverse order)
            if (start_to_root_idx != -1) {
                for (int i = start_to_root_idx - 1; i >= 0; i--) {
                    newList.push_back(startToRoot[i]);
                }
            }
            boneList = newList;
        }
    }

    bone_num = boneList.size();

    int total_dofs = 0;
    for (int i = 0; i < bone_num; i++) {
        acclaim::Bone* current = boneList[i];
        if (current->dofrx) total_dofs++;
        if (current->dofry) total_dofs++;
        if (current->dofrz) total_dofs++;
    }
    Eigen::Matrix4Xd Jacobian(4, total_dofs);
    Jacobian.setZero();
    for (int iter = 0; iter < max_iteration; ++iter) {
        forwardSolver(posture, root_bone);
        Eigen::Vector4d desiredVector = target_pos - end_bone->end_position;

        if (desiredVector.norm() < epsilon) {
            break;
        }
        // TODO#3-2 (compute jacobian)
        //   1. Compute arm vectors
        //   2. Compute jacobian columns, store in `Jacobian`
        // Hint:
        //   1. You should not put rotation in jacobian if it doesn't have that DoF.
        //   2. jacobian.col(/* some column index */) = /* jacobian column */
        Jacobian.setZero();
        int dof_idx = 0;
        for (int i = 0; i < bone_num; i++) {
            acclaim::Bone* current = boneList[i];
            // skip the bone that doesn't have rotation DOF
            if (!current->dofrx && !current->dofry && !current->dofrz) continue;

            Eigen::Vector3d rotate_axis;
            Eigen::Vector3d arm = (end_bone->end_position - current->start_position).head<3>();

            // x
            if (current->dofrx) {
                rotate_axis = current->rotation.rotation() * Eigen::Vector3d::UnitX();
                // ∂θi/∂p = ri​×(p−qi​)
                Eigen::Vector3d jacobian_col = rotate_axis.cross(arm);

                // Store in Jacobian, converting to 4D by padding with 0
                Jacobian.col(dof_idx) << jacobian_col, 0.0;
                dof_idx++;
            }

            // y
            if (current->dofry) {
                rotate_axis = current->rotation.rotation() * Eigen::Vector3d::UnitY();
                // ∂θi/∂p​ = ri​×(p−qi​)
                Eigen::Vector3d jacobian_col = rotate_axis.cross(arm);

                // Store in Jacobian, converting to 4D by padding with 0
                Jacobian.col(dof_idx) << jacobian_col, 0.0;
                dof_idx++;
            }

            // z
            if (current->dofrz) {
                rotate_axis = current->rotation.rotation() * Eigen::Vector3d::UnitZ();
                // ∂θi/∂p​ = ri​×(p−qi​)
                Eigen::Vector3d jacobian_col = rotate_axis.cross(arm);

                // Store in Jacobian, converting to 4D by padding with 0
                Jacobian.col(dof_idx) << jacobian_col, 0.0;
                dof_idx++;
            }
        }

        // TODO#3-3 (obstacle avoidance)
        //  1. Iterate through all bones in `boneList`.
        //  2. Compute the center of each bone (average of start and end positions).
        //  3. Calculate the vector from obstacle center to bone center.
        //  4. If distance is below threshold, compute repulsive vector.
        //  5. Add this repulsive vector to `desiredVector`.
        // Hint:
        // - Use a constant threshold distance to determine proximity.
        // - The repulsive vector should point away from the obstacle.
        // - Use `.head<3>().norm()` to compute 3D distance from a 4D vector.
        // - Normalize the repulsive vector and scale it based on how close it is.
        if (obsActive) {
            double half = 0.5;
            for (int i = 0; i < bone_num; i++) {
                acclaim::Bone* bone = boneList[i];
                Eigen::Vector4d mid_point = (bone->start_position +bone->start_position)/2.0;
                Eigen::Vector4d vector = mid_point - obs_pos;  // vector frobm obstacble to bone mid point
                Eigen::Vector3d dist = vector.head<3>().cwiseAbs();//x, y, z distance to the obstacle center
                bool inside_x = dist.x() <= half;
                bool inside_y = dist.y() <= half;
                bool inside_z = dist.z() <= half;

                if (inside_x && inside_y && inside_z) {
                    double d = dist.norm();
                    if (d < obsAvoidThreshold) {
                        Eigen::Vector4d inv_dir = vector.normalized();
                        Eigen::Vector4d repulse = inv_dir * (obsAvoidThreshold - d)*20.0;
                        desiredVector += repulse;
                    }
                }

                
            }
        }

        Eigen::VectorXd deltatheta = step * pseudoInverseLinearSolver(Jacobian, desiredVector);
        // TODO#3-4 (update rotation)
        //   Update `posture.bone_rotation` (in euler angle / degrees) using deltaTheta
        // Hint:
        //   1. You can ignore rotation limit of the bone.
        // Bonus:
        //   1. You cannot ignore rotation limit of the bone.
        int delta_idx = 0;
        for (int i = 0; i < bone_num; i++) {
            acclaim::Bone* current = boneList[i];
            int dof = 0;
            // x
            if (current->dofrx) {
                double new_rotate = posture.bone_rotations[current->idx][0] + deltatheta[delta_idx];
                if (new_rotate < current->rxmin) {
                    new_rotate = current->rxmin;
                } else if (new_rotate > current->rxmax) {
                    new_rotate = current->rxmax;
                }
                posture.bone_rotations[current->idx][0] = new_rotate;
                delta_idx++;
                dof++;
            }
            // y
            if (current->dofry) {
                double new_rotate = posture.bone_rotations[current->idx][1] + deltatheta[delta_idx];
                if (new_rotate < current->rymin) {
                    new_rotate = current->rymin;
                } else if (new_rotate > current->rymax) {
                    new_rotate = current->rymax;
                }
                posture.bone_rotations[current->idx][1] = new_rotate;
                delta_idx++;
                dof++;
            }
            // z
            if (current->dofrz) {
                double new_rotate = posture.bone_rotations[current->idx][2] + deltatheta[delta_idx];
                if (new_rotate < current->rzmin) {
                    new_rotate = current->rzmin;
                } else if (new_rotate > current->rzmax) {
                    new_rotate = current->rzmax;
                }
                posture.bone_rotations[current->idx][2] = new_rotate;
                delta_idx++;
                dof++;
            }
        }
    }
    // TODO#3-5
    // Return whether IK is stable
    // i.e. whether the ball is reachable
    // Hint:
    //      1. comment out the line here and return whether the IK is stable or not
    //      2. if the ball is reachable,  swinging its hand in air
    bool isStable = false;
    Eigen::Vector4d final_dist = target_pos - end_bone->end_position;
    double final_dist_norm = final_dist.norm();
    if (final_dist_norm < epsilon) {
        isStable = true;
    }

    return isStable;
}

}  // namespace kinematics
