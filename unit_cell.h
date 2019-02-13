#pragma once

#include <vector>
#include <string>

#include "math.h"
#include "file.h"

class unit_cell
{
public:
	double a, b, c;
	double alpha, beta, gamma;
	double volume;

	void init(std::vector<std::string> file) {
		if(!file.empty()) {
			load_file(file);
		}
		else {
			a = b = c = alpha = beta = gamma = volume = 0;
		}
	}

private:
	void load_file(std::vector<std::string> file) {
		std::vector<std::string> list;
		list.emplace_back("data_");
		list.emplace_back("_cell_length_a");
		list.emplace_back("_cell_length_b");
		list.emplace_back("_cell_length_c");
		list.emplace_back("_cell_angle_alpha");
		list.emplace_back("_cell_angle_beta");
		list.emplace_back("_cell_angle_gamma");
		list.emplace_back("_cell_volume");

		for (const auto& l : file) {
			auto token = io::tokenize(l, " ()\r\t");

			if (!token.empty()) {
				if (token[0] == list[1])
					a = io::stod(token[1]);
				else if (token[0] == list[2])
					b = io::stod(token[1]);
				else if (token[0] == list[3])
					c = io::stod(token[1]);
				else if (token[0] == list[4])
					alpha = io::stod(token[1]);
				else if (token[0] == list[5])
					beta = io::stod(token[1]);
				else if (token[0] == list[6])
					gamma = io::stod(token[1]);
				else if (token[0] == list[7])
					volume = io::stod(token[1]);
			}
		}
	}
};
