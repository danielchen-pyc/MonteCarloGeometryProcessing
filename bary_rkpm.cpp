#define GL_SILENCE_DEPRECATION

#include <iostream>
#include <vector>
#include <map>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <igl/AABB.h>
#include <cmath>
#include <glad/glad.h>
// #include <GLFW/glfw3.h>

#include <timer.h>
#include <ctime>
#include <boundary_setup.h>
#include <random_sampling.h>
#include <WoS_Laplacian_bary.h>
#include <WoS_Poisson_bary.h>
#include <WoS_biharmonic_bary.h>

// #include <igl/signed_distance.h>
// #include <igl/offset_surface.h>
// #include <igl/copyleft/offset_surface.h>

#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/colormap.h>


// Helper function to compute barycentric coordinates
void interpolate_phi(const Eigen::Vector3d& X, const Eigen::Vector3d& X0, const Eigen::Vector3d& X1, 
    const Eigen::Vector3d& X2, Eigen::Vector4d& phi) {
    // Triangle edges
    Eigen::Vector3d v0 = X1 - X0;
    Eigen::Vector3d v1 = X2 - X0;
    Eigen::Vector3d v2 = X - X0;

    // Areas for barycentric coordinates
    double d00 = v0.dot(v0);
    double d01 = v0.dot(v1);
    double d11 = v1.dot(v1);
    double d20 = v2.dot(v0);
    double d21 = v2.dot(v1);

    double denom = d00 * d11 - d01 * d01;

    // Compute barycentric coordinates
    phi(1) = (d11 * d20 - d01 * d21) / denom; // Weight for vertex X1
    phi(2) = (d00 * d21 - d01 * d20) / denom; // Weight for vertex X2
    phi(0) = 1.0 - phi(1) - phi(2);           // Weight for vertex X0
}


// Compute u_v(x) from M(x) and m(x)
Eigen::VectorXd computeUV(const Eigen::MatrixXd &M, const Eigen::VectorXd &m) {
    return M.inverse() * m;
}


// Compute the weighted matrix M(x) and vector m(x)
void computeMmv(
    const std::vector<Eigen::MatrixXd> &samplePoints,            // Cage sample points
    const std::vector<Eigen::MatrixXd> &vertexBasis,             // Basis functions
    const Eigen::MatrixXd &queryPoint,                           // Query point x
    const std::vector<Eigen::MatrixXd> &weights,                 // Weights for each sample point
    std::vector<Eigen::MatrixXd> &bary_coordinates
) {
    // std::cout << "38" << std::endl;
    size_t d = queryPoint.row(0).size();
    size_t n_v = vertexBasis[0].cols();

    // std::cout << samplePoints.rows() << "x" << samplePoints.cols() << " " << samplePoints.row(0) << std::endl;
    // std::cout << vertexBasis.rows() << "x" << vertexBasis.cols() << std::endl;
    // std::cout << weights.rows() << "x" << weights.cols() << std::endl;

    int weight = 1;



    for (int xi = 0; xi < samplePoints.size(); xi++) {
        bary_coordinates.push_back(Eigen::MatrixXd::Zero(n_v, d));
        for (size_t v = 0; v < vertexBasis[xi].cols(); v++) {
            Eigen::MatrixXd M_v = Eigen::MatrixXd::Zero(d+1, d+1);
            Eigen::VectorXd m_v = Eigen::VectorXd::Zero(d+1);
            bool test = false;
            for (size_t k = 0; k < vertexBasis[xi].rows(); k++) {
                Eigen::VectorXd y_0 = samplePoints[xi].row(k);
                double phi = vertexBasis[xi](k, v);
    
                Eigen::VectorXd y(d+1);
                y << y_0, 1;
    
                // Update M(x)
                M_v += y * y.transpose();
                
                // Update m(x)
                if (phi != 0.0) {
                    std::cout << phi << std::endl;
                    test = true;
                }
                m_v += phi * y;
            }
            // std::cout << "M: " << M_v << std::endl;
    
            size_t d = 3;
            Eigen::VectorXd x(d+1);
            Eigen::VectorXd x_0 = queryPoint.row(xi);
            // std::cout << "x: " << x_0 << std::endl;
            x << x_0, 1;
    
            Eigen::VectorXd bary_cord_0 = x.cwiseProduct(computeUV(M_v, m_v));
            Eigen::VectorXd bary_cord(d);
            bary_cord << bary_cord_0(0), bary_cord_0(1), bary_cord_0(2);

            if (test) {
                std::cout << "m: " << m_v(0) << " " << m_v(1) << " " << m_v(2) << " " << m_v(3) << " " << std::endl;
                std::cout << "bary " << v << ": " << bary_cord(0) << " " << bary_cord(1) << " " << bary_cord(2) << std::endl;
            }
            bary_coordinates[xi].row(v) = bary_cord;
        }
        std::cout << "Bary step: " << xi << std::endl;
    }

    std::cout << "Bary Completed! " << std::endl;
}

void reorganize_walks_by_point(
    const std::vector<Eigen::VectorXd>& U,
    const std::map<int, std::vector<Eigen::MatrixXd>>& sampled_yk_d,
    const std::map<int, std::vector<Eigen::MatrixXd>>& scalars_wk_d,
    // const std::map<int, std::vector<Eigen::MatrixXd>>& basis_fn_d,
    std::vector<Eigen::VectorXd>& U_by_point,
    std::vector<Eigen::MatrixXd>& sampled_yk_by_point,
    std::vector<Eigen::MatrixXd>& scalars_by_point,
    std::vector<Eigen::MatrixXd>& basis_fn_by_point
) {
    int n_walks = sampled_yk_d.begin()->second.size(); // should be 50
    int n_points = sampled_yk_d.size();
    int dimension = sampled_yk_d.begin()->second[0].rows();

    std::cout << n_walks << ", " << n_points << ", " << dimension << std::endl;

    // Initialize storage for each point
    for (int p = 0; p < n_points; ++p) {
        // U_by_point[p] = Eigen::VectorXd::Zero(n_walks);
        sampled_yk_by_point.push_back(Eigen::MatrixXd::Zero(n_walks, dimension)); // Rows = walks, rows = dimensions
        scalars_by_point.push_back(Eigen::MatrixXd::Zero(n_walks, n_points));
        basis_fn_by_point.push_back(Eigen::MatrixXd::Zero(n_walks, dimension+1));
    }

    // // Iterate over all walks
    for (auto const& x: sampled_yk_d) {
        int p = x.first;
        std::vector<Eigen::MatrixXd> value = x.second;
        for (int walk = 0; walk < value.size(); walk++)
            sampled_yk_by_point[p].row(walk) = value[walk].transpose();
    }

    for (auto const& x: scalars_wk_d) {
        int p = x.first;
        std::vector<Eigen::MatrixXd> value = x.second;
        for (int walk = 0; walk < value.size(); walk++)
            scalars_by_point[p].row(walk) = value[walk].transpose();
    }

    std::cout << "sampled_yk_by_point is a " << sampled_yk_by_point.size() << "x" << sampled_yk_by_point[0].rows() << "x" << sampled_yk_by_point[0].cols() << std::endl;
    std::cout << "scalars_by_point is a " << scalars_by_point.size() << "x" << scalars_by_point[0].rows() << "x" << scalars_by_point[0].cols() << std::endl;
}


// Function to compute vertex-based basis functions
void computeVertexBasedBasisFunctions(
    const Eigen::MatrixXd& SV,               // Cage vertices (n x 3)
    const Eigen::MatrixXi& SF,               // Cage faces (m x 3)
    const Eigen::MatrixXd& sampled_yk,       // Sample points y_k (p x 3)
    Eigen::MatrixXd& basis_fn                // Output basis functions (p x n)
) {
    // Initialize the AABB tree for fast spatial queries
    igl::AABB<Eigen::MatrixXd, 3> aabb;
    aabb.init(SV, SF);

    // Dimensions
    int p = sampled_yk.rows(); // Number of sample points
    int n = SV.rows();         // Number of cage vertices

    // Initialize the output matrix
    basis_fn = Eigen::MatrixXd::Zero(p, n);

    // Storage for closest face and barycentric coordinates
    Eigen::MatrixXd closest_points;
    Eigen::VectorXi closest_faces;
    Eigen::VectorXd distances_squared;

    // Find the closest triangle for each sample point
    aabb.squared_distance(SV, SF, sampled_yk, distances_squared, closest_faces, closest_points);

    // Iterate over each sample point
    for (int k = 0; k < p; ++k) {
        int f = closest_faces(k); // Closest face index
        Eigen::Vector3i face = SF.row(f); // Vertices of the closest face

        // Vertices of the triangle containing y_k
        Eigen::Vector3d v0 = SV.row(face(0));
        Eigen::Vector3d v1 = SV.row(face(1));
        Eigen::Vector3d v2 = SV.row(face(2));

        // Compute barycentric coordinates for y_k
        Eigen::Vector3d y_k = sampled_yk.row(k);
        Eigen::Vector4d phi; // Barycentric coordinates
        phi.setZero();

        // Interpolation helper function
        interpolate_phi(y_k, v0, v1, v2, phi);

        // Assign barycentric coordinates to the corresponding vertices
        basis_fn(k, face(0)) = phi(0); // Vertex 0 of the face
        basis_fn(k, face(1)) = phi(1); // Vertex 1 of the face
        basis_fn(k, face(2)) = phi(2); // Vertex 2 of the face
    }
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

    Eigen::MatrixXd SV;
    Eigen::MatrixXi SF;

	igl::read_triangle_mesh("../data/teapot.obj", V, F);

    // Get Bounding Mesh
	igl::read_triangle_mesh("../data/teapot-offset.obj", SV, SF);

    // Merge The meshes
    Eigen::MatrixXd mergedV(V.rows() + SV.rows(), V.cols());
    mergedV << V, SV;
    Eigen::MatrixXi mergedF(F.rows() + SF.rows(), F.cols());
    mergedF << F, (SF.array() + V.rows());

    std::cout << "V is a " << V.rows() << "x" << V.cols() << " matrix." << std::endl;
    std::cout << "F is a " << F.rows() << "x" << F.cols() << " matrix." << std::endl;

    std::cout << "SV is a " << SV.rows() << "x" << SV.cols() << " matrix." << std::endl;
    std::cout << "SF is a " << SF.rows() << "x" << SF.cols() << " matrix." << std::endl;

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
    int n_samples =2048;
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

    // Initialize Monte Carlo Variables
    Eigen::MatrixXd points;
    std::vector<Eigen::VectorXd> U;
    Eigen::VectorXd boundary = Eigen::VectorXd::Zero(SV.rows());
    Eigen::VectorXd Bh = Eigen::VectorXd::Zero(SV.rows()); // for Biharmonic only
    Eigen::MatrixXd BCM;

    // Initialize RKPM Variables
    std::map<int, std::vector<Eigen::MatrixXd>> sampled_yk_d;
    std::map<int, std::vector<Eigen::MatrixXd>> scalars_wk_d;
    std::map<int, std::vector<Eigen::MatrixXd>> basis_fn_d;


    // Output
    std::vector<Eigen::MatrixXd> bary_coordinates;
    

    // Reorganized
    std::vector<Eigen::VectorXd> U_by_point;
    std::vector<Eigen::MatrixXd> sampled_yk_by_point;
    std::vector<Eigen::MatrixXd> scalars_by_point;
    std::vector<Eigen::MatrixXd> basis_fn_by_point;


    // compute and show the current boundary
    const auto &set_boundary = [&]() {
        auto &pair = solver_funcs[solver_id].second[boundary_id];
        const auto &g = pair.second;

        if(solver_funcs[solver_id].first == "biharmonic") {
            auto &pair2 = solver_funcs[solver_id].second[boundary_id+1];
            const auto &h = pair2.second;
            for (int i = 0; i < SV.rows(); i++) {
                boundary(i) = g(SV.row(i).transpose());
                Bh(i) = h(V.row(i).transpose());
            }
            std::cout << "Boundary condition:\n\t" << pair.first << "\n\t" << pair2.first << "\n";
        } else {
            for (int i = 0; i < SV.rows(); i++)
                boundary(i) = g(SV.row(i).transpose());
            std::cout << "Boundary condition: " << pair.first << "\n";
        }

        igl::colormap(cmap_type, boundary, boundary.minCoeff(), boundary.maxCoeff(), BCM);
        viewer.data_list[yid].set_points(SV, BCM);
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
		random_sampling(V, F, n_samples, points);
        const Eigen::RowVector3d orange(1.0, 0.7, 0.2);
        viewer.data_list[xid].set_points(points, (1. - (1. - orange.array()) * .8));

        igl::colormap(cmap_type, boundary, boundary.minCoeff(), boundary.maxCoeff(), BCM);
        if(show_boundary)
            viewer.data_list[yid].set_points(V, BCM);
	};

    const auto& solve = [&]()
	{
	    auto &pair = solver_funcs[solver_id];
	    if (pair.first == "Laplacian") {
            WoS_Laplacian_Bary(SV, SF, boundary, points, U, sampled_yk_d, scalars_wk_d, basis_fn_d);
        }
	    // } else if (pair.first == "Poisson") {
        //     WoS_Poisson_Bary(SV, SF, boundary, source_terms[0], screened_c, use_pt_source, point_source, points, U);
	    // } else if (pair.first == "biharmonic") {
        //     WoS_biharmonic_Bary(SV, SF, boundary, Bh, use_pt_source, point_source, points, U);
	    // }

        std::cout << "WoS Completed!" << std::endl;

        reorganize_walks_by_point(U, sampled_yk_d, scalars_wk_d, 
            U_by_point, sampled_yk_by_point, scalars_by_point, basis_fn_by_point);
    
        // Compute the vertex-based basis functions
        std::cout << "Computing Basis Functions" << std::endl;
        for (int i = 0; i < sampled_yk_by_point.size(); i++) {
            computeVertexBasedBasisFunctions(SV, SF, sampled_yk_by_point[i], basis_fn_by_point[i]);
            std::cout << "Basis Function: " << i << std::endl;
        }
        std::cout << "Completed: Basis Functions" << std::endl;

        // RKPM Method
        computeMmv(sampled_yk_by_point, basis_fn_by_point, points, scalars_by_point, bary_coordinates);


        double min_scale = std::min(U[0].minCoeff(), boundary.minCoeff());
        double max_scale = std::max(U[0].maxCoeff(), boundary.maxCoeff());


        for (int i = 0; i < bary_coordinates.size(); i++) {
            Eigen::MatrixXd CM;
            Eigen::MatrixXd points = bary_coordinates[i];
            const Eigen::RowVector3d orange(1.0, 0.7, i / bary_coordinates.size());

            // igl::colormap(cmap_type, points, min_scale, max_scale, CM);
            viewer.data_list[xid].add_points(points, (1. - (1. - orange.array()) * .8));
            // viewer.data_list[xid].add_points(sampled_yk_by_point[i], (1. - (1. - orange.array()) * .8));
        }

        std::cout << "Added " << bary_coordinates.size() * bary_coordinates[0].rows() << " points" << std::endl;
	
        igl::colormap(cmap_type, boundary, min_scale, max_scale, BCM);
        if (show_boundary)
            viewer.data_list[yid].set_points(SV, BCM);
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
			std::cout << "# Sample points: " << points.rows() << "\n";
			break;
		case '-':
			n_samples /= 2;
			set_points();
			std::cout << "# Sample points: " << points.rows() << "\n";
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
    // Parameters
    double sigma = 0.1; // Kernel bandwidth

    viewer.data().set_mesh(SV, SF);
    // viewer.data().add_points(SV, Eigen::RowVector3d(1,0,0));
    std::cout << "Changed to solver: " << solver_funcs[solver_id].first << "\n";

    set_boundary();
	set_points();
    std::cout << "# Sample points: " << points.rows() << "\n";

    viewer.data().show_lines = false;
	viewer.launch();

	return EXIT_SUCCESS;
}