// CRYSTAL.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "file.h"

#include "unit_cell.h"
#include "niggli.h"
#include <fstream>

//#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
//#include <CGAL/Periodic_3_Delaunay_triangulation_traits_3.h>
//#include <CGAL/Periodic_3_Delaunay_triangulation_3.h>
//#include <CGAL/periodic_3_triangulation_3_io.h>
//#include <fstream>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel       K;
//typedef CGAL::Periodic_3_Delaunay_triangulation_traits_3<K>       Gt;
//typedef CGAL::Periodic_3_Delaunay_triangulation_3<Gt>             P3DT3;
//typedef P3DT3::Point             Point;
//typedef P3DT3::Iso_cuboid        Iso_cuboid;
//typedef P3DT3::Vertex_handle     Vertex_handle;
//typedef P3DT3::Cell_handle       Cell_handle;
//typedef P3DT3::Locate_type       Locate_type;

void cgal(math::dvec3 point) {
	//std::list<Point> L;
	//L.push_back(Point(point.x, point.y, point.z));

	//Iso_cuboid domain(-9, -9, -9, 8, 8, 8);

	//P3DT3 T(L.begin(), L.end(), domain);
}
// check stability for niggli - perturb -> niggli change but voronoi doesn't
//voronoic polyhedra
// inkscape.org - drawing

int main()
{
	int numReductionFailures = 0;
	int numAngleFailures = 0;

	std::ofstream outfile("global_reduction_results");

	for (int i = 1; i <= 5688; i++) {
		{
			std::string file = "molGeom/T2_" + std::to_string(i) + "_num_molGeom.cif";
			std::ifstream infile{ file };

			infile.open(file);
			
			//std::string file_contents{ std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>() };
			std::stringstream buffer;
			buffer << infile.rdbuf();
			std::string file_contents = buffer.str();

			infile.close();

			const auto v = io::split(file_contents);

			double c1x, c1y, c1z;

			c1x = 0.99;
			c1y = 0.05;
			c1z = 0.135;

			unit_cell uc;
			uc.init(v);

			niggli n;
			if (n.reduce(uc, 0.0000001) == false)
				numReductionFailures++;
			else
				outfile << uc.a << " " << uc.b << " " << uc.c << " " << uc.alpha << " " << uc.beta << " " << uc.gamma << std::endl;

			if (uc.alpha <= 60 || uc.beta <= 60 || uc.gamma <= 60)
				numAngleFailures++;
		}
	}

	outfile << std::endl;
	outfile << "Reduction Failures: " << numReductionFailures << std::endl;
	outfile << "Angle Failures (<=60): " << numAngleFailures;

	return 0;
}