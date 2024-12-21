#ifndef WOS_LAPLACIAN_BARY_H
#define WOS_LAPLACIAN_BARY_H

#include <Eigen/Core>

// Solve ∆u = 0 over space at given points P subject to B on the given boundary
// mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
//   B  #V by 1 list of Dirichlet boundary conditions
//   P  #P by 3 list of query positions
// Outputs:
//   U  #P by 1 list of values at query positions
void WoS_Laplacian_Bary(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	const Eigen::VectorXd& B,
	const Eigen::MatrixXd& P,
	std::vector<Eigen::VectorXd>& U,
	std::map<int, std::vector<Eigen::MatrixXd>>& sampled_yk_d,
    std::map<int, std::vector<Eigen::MatrixXd>>& scalars_wk_d,
    std::map<int, std::vector<Eigen::MatrixXd>>& basis_fn_d);

#endif // !WOS_LAPLACIAN_BARY_H