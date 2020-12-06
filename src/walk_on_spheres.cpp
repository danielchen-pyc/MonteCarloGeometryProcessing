#include <walk_on_spheres.h>
#include <pts_on_spheres.h>
#include <interpolate.h>

#include <igl/AABB.h>
#include <igl/slice_mask.h>

#include <iostream>

void walk_on_spheres(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
	const Eigen::MatrixXd& P, Eigen::VectorXd& U) 
{
	igl::AABB<Eigen::MatrixXd, 3> aabb;
	aabb.init(V, F);

	const float tol = 1e-4;
	const int max_itr = 5;
	const int n_walks = 64;

	U = Eigen::VectorXd::Zero(P.rows());
	
	// start the random walk
	for (int i = 0; i < n_walks; i++) {
		Eigen::MatrixXd X = P;
		Eigen::VectorXd R = std::numeric_limits<double>::max() * Eigen::VectorXd::Ones(P.rows());

		int itr = 0;
		Eigen::VectorXd sqrD;
		Eigen::VectorXi I;
		Eigen::MatrixXd C;
		while (itr < max_itr) {
			// find closest point on the boundary and change radius
			aabb.squared_distance(V, F, X, sqrD, I, C);
			R = R.cwiseMin(sqrD.cwiseSqrt());
			pts_on_spheres(R, X);
			itr++;
		}

		for (int j = 0; j < I.rows(); j++) {
			// interpolate between the three vertices and find value
			Eigen::Vector3d phi;
			interpolate(C.row(j), V.row(F(I(j), 0)), V.row(F(I(j), 1)), V.row(F(I(j), 2)), phi);
			U(j) += phi(0) * B(F(I(j),0)) + phi(1) * B(F(I(j),1)) + phi(2) * B(F(I(j),2));
		}
	}
	U = U / n_walks;
}