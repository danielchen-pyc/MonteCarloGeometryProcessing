#include <WoS_Laplacian.h>
#include <sample_on_spheres.h>
#include <interpolate.h>

#include <igl/AABB.h>

void WoS_Laplacian_Bary(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& B,
	const Eigen::MatrixXd& P, Eigen::VectorXd& U) 
{
	igl::AABB<Eigen::MatrixXd, 3> aabb;
	aabb.init(V, F);

	const int max_itr = 5;
	const int n_walks = 20;

    Eigen::MatrixXd query_points;
    Eigen::MatrixXd closest_matrix;
    Eigen::VectorXd radii;
    Eigen::VectorXd distances_squared;
    Eigen::VectorXi closest_indices;
    U = Eigen::VectorXd::Zero(P.rows());

    const auto& handle_done_pts = [&](const Eigen::VectorXi &I_) {
        Eigen::Vector3d phi;
        for (int j = 0; j < I_.rows(); j++) {
            // interpolate between the three vertices and find value
            int f = I_(j);
            interpolate(closest_matrix.row(j), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), phi);
            U(j) += phi(0) * B(F(f,0)) + phi(1) * B(F(f,1)) + phi(2) * B(F(f,2));
        }
    };

	// start the random walks
	for (int i = 0; i < n_walks; i++) {
		query_points = P;
        int itr = 0;
		while (itr < max_itr) {
			aabb.squared_distance(V, F, query_points, distances_squared, closest_indices, closest_matrix);
            radii = distances_squared.cwiseSqrt();
			sample_on_spheres(radii, query_points);
			itr++;
		}
		handle_done_pts(closest_indices);
	}
	U = U / n_walks;
}