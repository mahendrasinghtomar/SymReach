#include <iostream>
#include <vector>
#include <boost/numeric/interval.hpp>

using namespace boost::numeric;
typedef boost::numeric::interval<double> intervalD;
typedef std::vector<intervalD> vec_interval;

extern "C" void func_from_file(const vec_interval& x, vec_interval& Hessi){
	Hessi[4] = -2.0*x[1];
	Hessi[5] = -2.0*x[0];
	Hessi[6] = -2.0*x[0];
}