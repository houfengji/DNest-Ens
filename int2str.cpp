#include <iostream>
#include <sstream>
#include <string>
#include "int2str.h"

void int2str( const unsigned long i,
              const size_t width,   // how many digits to output
              std::string & str) {
	if (width >= 20) {
		std::cerr << "A width of " << width << " is too large!" << std::endl;
	}
	std::stringstream tmp;
	tmp << i;
	str = tmp.str();
	while (str.size() < width) {
		str = "0" + str;
	}
}

std::string int2str(const unsigned long i,
                    const size_t width){   // how many digits to output
	if (width >= 20) {
		std::cerr << "A width of " << width << " is too large!" << std::endl;
	}
	std::stringstream tmp;
	std::string str;
	tmp << i;
	str = tmp.str();
	while (str.size() < width) {
		str = "0" + str;
	}
	return str;
}

std::string int2str(const unsigned long i){
	std::stringstream tmp;
	std::string str;
	tmp << i;
	str = tmp.str();
	return str;
}
