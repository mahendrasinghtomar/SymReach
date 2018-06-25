#include <iostream>
#include <vector>
#include <boost/numeric/interval.hpp>

using namespace boost::numeric;
typedef boost::numeric::interval<double> intervalD;
typedef std::vector<intervalD> vec_interval;

extern "C" void func_from_file(const vec_interval& x, vec_interval& Hessi){
	Hessi[107] = -8.0000000000000004e-01;
	Hessi[113] = -8.0000000000000004e-01;
	Hessi[164] = -1.3000000000000000e+00;
	Hessi[170] = -1.3000000000000000e+00;
	Hessi[221] = -1.0;
	Hessi[227] = -1.0;
	Hessi[307] = -1.5000000000000000e+00;
	Hessi[337] = -1.5000000000000000e+00;
}