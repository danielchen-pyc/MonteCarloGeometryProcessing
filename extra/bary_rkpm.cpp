#define GL_SILENCE_DEPRECATION

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <cmath>
#include <glad/glad.h>
// #include <GLFW/glfw3.h>

#include <timer.h>
#include <ctime>
#include <boundary_setup.h>
#include <random_sampling.h>
#include <WoS_Laplacian.h>
#include <WoS_Poisson.h>
#include <WoS_biharmonic.h>

// #include <igl/signed_distance.h>
// #include <igl/offset_surface.h>
// #include <igl/copyleft/offset_surface.h>

#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>


struct Point3D {
    float x, y, z;
};

// Kernel function (e.g., Gaussian kernel)
double kernelFunction(const Eigen::VectorXd &x, const Eigen::VectorXd &y, double sigma) {
    double norm = (x - y).squaredNorm();
    return exp(-norm / (2 * sigma * sigma));
}

// Compute the weighted matrix M(x) and vector m(x)
void computeKernelMatrixAndVector(
    const std::vector<Eigen::VectorXd> &samplePoints, // Cage sample points
    const std::vector<Eigen::VectorXd> &vertexBasis, // Basis functions
    const Eigen::VectorXd &queryPoint,              // Query point x
    const std::vector<double> &weights,             // Weights for each sample point
    Eigen::MatrixXd &M,                             // Output matrix M(x)
    Eigen::VectorXd &m,                             // Output vector m(x)
    double sigma                                    // Kernel bandwidth parameter
) {
    size_t d = queryPoint.size();
    M = Eigen::MatrixXd::Zero(d + 1, d + 1);
    m = Eigen::VectorXd::Zero(d + 1);

    for (size_t i = 0; i < samplePoints.size(); ++i) {
        double kappa = kernelFunction(queryPoint, samplePoints[i], sigma);
        Eigen::VectorXd y = vertexBasis[i];
        Eigen::VectorXd weightedY = kappa * weights[i] * y;

        // Update M(x)
        M.block(0, 0, d, d) += kappa * weights[i] * y * y.transpose();
        M(d, d) += kappa * weights[i];

        // Update m(x)
        m.head(d) += weightedY;
        m(d) += kappa * weights[i];
    }
}

// Compute u_v(x) from M(x) and m(x)
Eigen::VectorXd computeUV(const Eigen::MatrixXd &M, const Eigen::VectorXd &m) {
    return M.inverse() * m;
}

// std::pair<Eigen::MatrixXd, Eigen::MatrixXi> getBoundingMesh(Eigen::MatrixXd V, Eigen::MatrixXd F, int s, double isolevel) {
//     Eigen::MatrixXd SV;
//     Eigen::MatrixXi SF;
//     Eigen::MatrixXd GV;
//     Eigen::VectorXd S;
//     Eigen::Vector3i side;

//     isolevel = isolevel < 0 ? 0 : isolevel;
//     s = s < 0 ? 0 : s;
//     //igl::SignedDistanceType signed_distance_type = igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_DEFAULT;
//     igl::SignedDistanceType signed_distance_type = igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_FAST_WINDING_NUMBER;
//     std::cout << "Start Offset Surface - isoLevel: " << isolevel << ", s: " << s << std::endl;
//     clock_t startTime = clock();
//     igl::copyleft::offset_surface(V, F, isolevel, s, signed_distance_type, SV, SF, GV, side, S);
//     std::cout << "End offset surface: " << clock() - startTime << std::endl;
//     return std::make_pair(SV, SF);
// }


int main(int argc, char* argv[]) {
    Eigen::MatrixXd V;
	Eigen::MatrixXi F;
    size_t dimension = 3;    // Dimension of the problem (e.g., 3D)

	igl::opengl::glfw::Viewer viewer;
    const int xid = viewer.selected_data_index;
    viewer.data().point_size = 5;
	viewer.append_mesh();
	const int yid = viewer.selected_data_index;
    viewer.data().point_size = 5;

    std::cout << R"(
[space] Solves the PDE with current sampled points
H,h     Hide boundary points
S,s     Show boundary points
+       Double the number of sample points
-       Halve the number of sample points
R,r     Resets the sampled points
G,g     Changes the pre-listed boundary conditions
M,m     Changes the WoS solver (Laplacian, Poisson, biharmonic)
C       Double the constant c for Screened Poisson
c       Halve the constant c for Screened Poisson
X,x     Switch between regular Poisson and Screened Poisson
    )";
    std::cout << "\n";

	igl::read_triangle_mesh((argc > 1 ? argv[1] : "data/bunny-offset.off"), V, F);

    std::cout << "V is a " << V.rows() << "x" << V.cols() << " matrix." << std::endl;
    std::cout << "F is a " << F.rows() << "x" << F.cols() << " matrix." << std::endl;

    std::cout << V.row(0) << std::endl;

    // Get Bounding Mesh
    // const std::pair<Eigen::MatrixXd, Eigen::MatrixXi> p = getBoundingMesh(V, F, 100, 0.02);
    // Eigen::MatrixXd SV = p.first;
    // Eigen::MatrixXi SF = p.second;


    // Monte Carlo
    const Eigen::Vector3d point_source = Eigen::Vector3d(0.5, 0.5, 0.5);
    std::vector<std::function<double(const Eigen::Vector3d&)>> source_terms;
    source_terms.emplace_back([&](const Eigen::Vector3d &x) -> double {
        double r2 = (x - point_source).squaredNorm();
        return 1e4 * std::pow(exp(1.0), -r2);
    });
    std::vector<std::string> source_strings;
    source_strings.emplace_back("f(y) = Dirac delta");

    // Configurable variables
    int n_samples = 1048576;
    int boundary_id = 0;
    int solver_id = 0;
    double screened_c = 0;
    bool show_boundary = true;
    bool use_pt_source = true;
    igl::ColorMapType cmap_type = igl::COLOR_MAP_TYPE_MAGMA;

    // sets up the built-in PDE and boundary conditions
    std::vector<std::pair<std::string, std::vector<f_pairs>>> solver_funcs;
    boundary_setup(point_source, solver_funcs);

    double solve_time;
    Eigen::MatrixXd P;
    Eigen::VectorXd U;
    Eigen::VectorXd B = Eigen::VectorXd::Zero(V.rows());
    Eigen::VectorXd Bh = Eigen::VectorXd::Zero(V.rows()); // for Biharmonic only
    Eigen::MatrixXd BCM;

    // compute and show the current boundary
    const auto &set_boundary = [&]() {
        auto &pair = solver_funcs[solver_id].second[boundary_id];
        const auto &g = pair.second;

        if(solver_funcs[solver_id].first == "biharmonic") {
            auto &pair2 = solver_funcs[solver_id].second[boundary_id+1];
            const auto &h = pair2.second;
            for (int i = 0; i < V.rows(); i++) {
                B(i) = g(V.row(i).transpose());
                Bh(i) = h(V.row(i).transpose());
            }
            std::cout << "Boundary condition:\n\t" << pair.first << "\n\t" << pair2.first << "\n";
        } else {
            for (int i = 0; i < V.rows(); i++)
                B(i) = g(V.row(i).transpose());
            std::cout << "Boundary condition: " << pair.first << "\n";
        }

        igl::colormap(cmap_type, B, B.minCoeff(), B.maxCoeff(), BCM);
        viewer.data_list[yid].set_points(V, BCM);
    };

    // shows the boundary
    const auto& show_bndry_pts = [&]()
    {
        if (show_boundary)
            viewer.data_list[yid].set_points(V, BCM);
        else
            viewer.data_list[yid].clear_points();
    };

    // resample and set the data points
	const auto& set_points = [&]()
	{
		random_sampling(V, F, n_samples, P);
        const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
        viewer.data_list[xid].set_points(P, (1. - (1. - orange.array()) * .8));

        igl::colormap(cmap_type, B, B.minCoeff(), B.maxCoeff(), BCM);
        if(show_boundary)
            viewer.data_list[yid].set_points(V, BCM);
	};

    const auto& solve = [&]()
	{
	    auto &pair = solver_funcs[solver_id];
	    if(pair.first == "Laplacian") {
            WoS_Laplacian(V, F, B, P, U);
	    } else if(pair.first == "Poisson") {
            WoS_Poisson(V, F, B, source_terms[0], screened_c, use_pt_source, point_source, P, U);
	    } else if(pair.first == "biharmonic") {
            WoS_biharmonic(V, F, B, Bh, use_pt_source, point_source, P, U);
	    }

        double min_scale = std::min(U.minCoeff(), B.minCoeff());
        double max_scale = std::max(U.maxCoeff(), B.maxCoeff());

		Eigen::MatrixXd CM;
		igl::colormap(cmap_type, U, min_scale, max_scale, CM);
		viewer.data_list[xid].set_points(P, CM);

        igl::colormap(cmap_type, B, min_scale, max_scale, BCM);
        if (show_boundary)
            viewer.data_list[yid].set_points(V, BCM);
	};

    // updates
    const auto update = [&]() { show_bndry_pts(); };

	viewer.callback_key_pressed = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
	{
		switch (key)
		{
		case ' ':
		    std::cout << "Solving...\n";
		    solve_time = timer(solve);
			std::cout << "Solve time: " << solve_time << " sec\n";
			break;
		case '+':
			n_samples *= 2;
			set_points();
			std::cout << "# Sample points: " << P.rows() << "\n";
			break;
		case '-':
			n_samples /= 2;
			set_points();
			std::cout << "# Sample points: " << P.rows() << "\n";
			break;
        case 'S':
		case 's':
			show_boundary = true;
            std::cout << "Showing Boundry." << "\n";
			break;
		case 'H':
		case 'h':
			show_boundary = false;
            std::cout << "Hiding Boundry." << "\n";
			break;
        case 'G':
        case 'g':
            boundary_id = (boundary_id + 1) % solver_funcs[solver_id].second.size();
            set_boundary();
            break;
        case 'M':
        case 'm':
            solver_id = (solver_id + 1) % solver_funcs.size();
            std::cout << "\nChanged to solver: " << solver_funcs[solver_id].first << "\n";
            if(solver_funcs[solver_id].first != "Laplacian" && use_pt_source)
                std::cout << "Point source: ON at (" << point_source.transpose() << ")\n";
            if(solver_funcs[solver_id].first == "Poisson" && screened_c > 0)
                std::cout << "Screened Poisson c: " << screened_c << "\n";
            boundary_id = 0;
            set_boundary();
		case 'R':
		case 'r':
			set_points();
			break;
        case 'Q':
        case 'q':
            use_pt_source = !use_pt_source;
            std::cout << (use_pt_source ? "Enabled " : "Disabled ") << "point source at (" << point_source.transpose() << ")\n";
            break;
        case 'C':
            if(solver_funcs[solver_id].first == "Poisson") {
                screened_c *= 2;
                std::cout << "Screened Poisson c: " << screened_c << "\n";
            }
            break;
        case 'c':
            if(solver_funcs[solver_id].first == "Poisson") {
                screened_c /= 2;
                std::cout << "Screened Poisson c: " << screened_c << "\n";
            }
            break;
        case 'X':
        case 'x':
            std::cout << ((screened_c > 0) ? "Turned off " : "Turned on ") << "Screened Poisson\n";
            screened_c = (screened_c > 0) ? 0 : 1.;
            if(screened_c > 0.)
                std::cout << "Screened Poisson c: " << screened_c << "\n";
            break;
		default:
			return false;
		}
		update();
		return true;
	};


    // RKPM
    Eigen::VectorXd queryPoint = Eigen::VectorXd::Random(dimension);

    // Parameters
    double sigma = 0.1; // Kernel bandwidth

    // Compute M(x) and m(x)
    // Eigen::MatrixXd M;
    // Eigen::VectorXd m;
    // computeKernelMatrixAndVector(samplePoints, vertexBasis, queryPoint, weights, M, m, sigma);

    // // Compute u_v(x)
    // Eigen::VectorXd uv = computeUV(M, m);

    // TODO: Figure out above

    viewer.data().set_mesh(V, F);
    // viewer.data().add_points(SV, Eigen::RowVector3d(1,0,0));
    std::cout << "Changed to solver: " << solver_funcs[solver_id].first << "\n";

    set_boundary();
	set_points();
    std::cout << "# Sample points: " << P.rows() << "\n";

    viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}



// // Main function
// int main() {
//     // Example setup
//     size_t numSamples = 100; // Number of cage sample points
//     size_t dimension = 2;    // Dimension of the problem (e.g., 2D)

//     // Generate some sample points (example)
//     vector<Eigen::VectorXd> samplePoints(numSamples);
//     vector<Eigen::VectorXd> vertexBasis(numSamples);
//     vector<double> weights(numSamples, 1.0); // Uniform weights

//     for (size_t i = 0; i < numSamples; ++i) {
//         samplePoints[i] = Eigen::VectorXd::Random(dimension);
//         vertexBasis[i] = Eigen::VectorXd::Random(dimension);
//     }

//     // Query point
//     Eigen::VectorXd queryPoint = Eigen::VectorXd::Random(dimension);

//     // Parameters
//     double sigma = 0.1; // Kernel bandwidth

//     // Compute M(x) and m(x)
//     Eigen::MatrixXd M;
//     Eigen::VectorXd m;
//     computeKernelMatrixAndVector(samplePoints, vertexBasis, queryPoint, weights, M, m, sigma);

//     // Compute u_v(x)
//     Eigen::VectorXd uv = computeUV(M, m);

//     // Output the result
//     cout << "Computed u_v(x): \n" << uv << endl;

//     return 0;
// }