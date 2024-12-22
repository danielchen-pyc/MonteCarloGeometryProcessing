#include <igl/read_triangle_mesh.h>
#include <igl/write_triangle_mesh.h>
#include <igl/signed_distance.h>
#include <igl/offset_surface.h>

int main(int argc, char *argv[])
{
    Eigen::MatrixXd V,SV;
    Eigen::MatrixXi F,SF;
    Eigen::MatrixXd GV;
    Eigen::VectorXd S;
    Eigen::RowVector3i side;

    double isolevel = 0.02;
    int s = 100;
    igl::SignedDistanceType signed_distance_type = igl::SignedDistanceType::SIGNED_DISTANCE_TYPE_DEFAULT;
   
    igl::read_triangle_mesh("data/teapot.obj", V,F);

    // https://github.com/libigl/libigl-legacy/blob/master/include/igl/copyleft/offset_surface.cpp
    igl::offset_surface(V, F, isolevel, s, signed_distance_type, SV, SF, GV, side, S);
    igl::write_triangle_mesh("data/teapot-offset.obj", SV, SF);

    return 0;
}