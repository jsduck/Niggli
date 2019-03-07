//--- Standard Includes
#include <fstream>

//--- Niggli Includes
#include "file.h"
#include "cifio.h"
#include "niggli.h"


int main(int argc, char* argv[])
{
	int numReductionFailures = 0;
	int numAngleFailures = 0;

	std::ofstream outfile("global_reduction_results");

	CIFHandler ch;

	std::string file = "molGeom/T2_1_num_molGeom.cif";
	std::ifstream infile{ file };

	infile.open(file);

	std::string file_contents{ std::istreambuf_iterator<char>(infile), std::istreambuf_iterator<char>() };

	infile.close();

	auto cf = ch.copyCIFrom(cif::read_file(file));
	auto gb = cf->blocks[0];

	if (cf->blocks.size() != 1)
		gb = cf->blocks[1];

	std::vector<double> vals;

	vals.push_back(io::stod(*gb.find_value("_cell_length_a")));
	vals.push_back(io::stod(*gb.find_value("_cell_length_b")));
	vals.push_back(io::stod(*gb.find_value("_cell_length_c")));
	vals.push_back(io::stod(*gb.find_value("_cell_angle_alpha")));
	vals.push_back(io::stod(*gb.find_value("_cell_angle_beta")));
	vals.push_back(io::stod(*gb.find_value("_cell_angle_gamma")));

	if (!niggli::reduce(vals, 0.0000001, 100))
		numReductionFailures++;
	else {
		for (auto &val : vals)
			outfile << val << " ";
		outfile << std::endl;
	}
		
	if (vals[3] < 60 || vals[4] < 60 || vals[5] < 60)
		numAngleFailures++;

	outfile << std::endl;
	outfile << "Reduction Failures: " << numReductionFailures << std::endl;
	outfile << "Angle Failures (<60): " << numAngleFailures;

	return 0;
}
