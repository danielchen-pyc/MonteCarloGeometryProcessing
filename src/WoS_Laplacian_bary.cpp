#include <WoS_Laplacian.h>
#include <sample_on_spheres.h>
#include <interpolate.h>
#include <iostream>
#include <map>


#include <igl/AABB.h>

void WoS_Laplacian_Bary(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXd& Boundary,
	const Eigen::MatrixXd& P, std::vector<Eigen::VectorXd>& U, std::map<int, std::vector<Eigen::MatrixXd>>& sampled_yk_d,
    std::map<int, std::vector<Eigen::MatrixXd>>& scalars_wk_d,
    std::map<int, std::vector<Eigen::MatrixXd>>& basis_fn_d)
{
	igl::AABB<Eigen::MatrixXd, 3> aabb;
	aabb.init(V, F);

	const int max_itr = 20;
	const int n_walks = 50;

    Eigen::Vector3d min_corner = V.colwise().minCoeff();
    Eigen::Vector3d max_corner = V.colwise().maxCoeff();
    double cage_threshold = 0.001 * (max_corner - min_corner).norm(); // 0.1% of diagonal
    // double cage_threshold = 0.01;
    std::cout << "Cage Threshold: " << cage_threshold << std::endl;

    Eigen::MatrixXd query_points;
    Eigen::MatrixXd closest_matrix;
    Eigen::VectorXd radii;
    Eigen::VectorXd distances_squared;
    Eigen::VectorXi closest_indices;
    Eigen::VectorXd U_m = Eigen::VectorXd::Zero(P.rows());

    int d_plus_1 = V.cols()+1;
    Eigen::MatrixXd basis_fn_unfiltered = Eigen::MatrixXd::Zero(V.rows(), d_plus_1);
    
    // Vector to store converged points
    Eigen::MatrixXd converged_points_matrix;

    // Helper to handle converged points
    const auto& handle_done_pts = [&](const Eigen::VectorXi& idx) {
        Eigen::VectorXd phi(d_plus_1);
        phi.setZero();
        for (int j = 0; j < idx.rows(); j++) {
            // interpolate between the three vertices and find value
            int f = idx(j);

            interpolate(closest_matrix.row(j), V.row(F(f, 0)), V.row(F(f, 1)), V.row(F(f, 2)), phi);
            U_m(j) += phi(0) * Boundary(F(f,0)) + phi(1) * Boundary(F(f,1)) + phi(2) * Boundary(F(f,2));
            basis_fn_unfiltered.row(j) = phi;
        }
        U.push_back(U_m);
    };

    // need to keep sampling until we have m = 50 converged sampled yk
    // use a dictionary to keep track of how many converged yk each query point x has

    std::map<int, int> d;
    int count = 0;
    int curr_itr = 0;

	// start the random walks
	while (count < P.rows()) {
		query_points = P;
        int itr = 0;
		while (itr < max_itr) {
			aabb.squared_distance(V, F, query_points, distances_squared, closest_indices, closest_matrix);
            radii = distances_squared.cwiseSqrt();
			sample_on_spheres(radii, query_points);
			itr++;
		}

        handle_done_pts(closest_indices);

        // Check convergence
        Eigen::VectorXi converged = (radii.array() < cage_threshold).cast<int>();

        // int curr_row = 0;
        for (int i = 0; i < converged.rows(); i++) {
            if (converged[i]) {
                if (d.count(i) > 0) {
                    if (d[i] <= n_walks && d[i] > 0) {
                        d[i]++;
                        sampled_yk_d[i].push_back(closest_matrix.row(i).transpose());
                        Eigen::VectorXd s = Eigen::VectorXd::Zero(converged.rows());
                        s(i) = 1;
                        scalars_wk_d[i].push_back(s);
                        basis_fn_d[i].push_back(basis_fn_unfiltered.row(i).transpose());
                    } else if (d[i] == -1) continue;
                    else {
                        // For this query point at index i, 
                        // it reached n_walks sampled yk points
                        d[i] = -1;
                        count++;
                        // std::cout << "Curr_itr: " << curr_itr << ", Count: " << count << std::endl;
                        continue;
                    }
                } else d[i] = 1;
            }
        }

        curr_itr++;
	}

    for (int i = 0; i < U.size(); i++) {
        U[i] /= n_walks;
    }

    // std::cout << sampled_yk_d.size() << ", " << scalars_wk_d.size() << ", " << basis_fn_d.size() << ", " << U.size() << std::endl;
    for (int pt = 0; pt < sampled_yk_d.size(); pt++) {
        // std::cout << pt << ", " << sampled_yk_d[pt].size() << ", " << scalars_wk_d[pt].size() << ", " << basis_fn_d[pt].size() << std::endl;
        assert(sampled_yk_d[pt].size() == 50);
        assert(scalars_wk_d[pt].size() == 50);
        assert(basis_fn_d[pt].size() == 50);
    }

    std::cout << "sampled_yk_d is a " << sampled_yk_d.size() << "x" << sampled_yk_d[0].size() << "x" << sampled_yk_d[0][0].rows() << "x" << sampled_yk_d[0][0].cols() << " matrix." << std::endl;
    std::cout << "scalars_wk_d is a " << scalars_wk_d.size() << "x" << scalars_wk_d[0].size() << "x" << scalars_wk_d[0][0].rows() << "x" << scalars_wk_d[0][0].cols() << " matrix." << std::endl;
    std::cout << "basis_fn_d is a " << basis_fn_d.size() << "x" << basis_fn_d[0].size() << "x" << basis_fn_d[0][0].rows() << "x" << basis_fn_d[0][0].cols() << " matrix." << std::endl;
}