/*
 * ReachableSet2.h
 *
 *  Created on: Jun 9, 2018
 *      Author: MahendraSinghTomar
 */

#ifndef REACHABLESET2_H_
#define REACHABLESET2_H_

#include <iostream>
#include <ostream>
#include <fstream>
#include <string>
#include <numeric>
#include <algorithm>
#include <Eigen/Dense>
#include <unsupported/Eigen/MatrixFunctions>
#include "UsingGnuplot.h"
#include <boost/numeric/interval.hpp>
#include <ginac/ginac.h>
#include "TicToc.hh"

// namespace Eigen {
  // namespace internal {
    // template<typename X, typename S, typename P>
    // struct is_convertible<X,boost::numeric::interval<S,P> > ;

    // template<typename S, typename P1, typename P2>
    // struct is_convertible<boost::numeric::interval<S,P1>,boost::numeric::interval<S,P2> > ;
  // }
// }
namespace vnodelp{
	typedef boost::numeric::interval<double> interval;

	double sup(interval& I);

	double inf(interval& I);

	double mag(interval& I);

	typedef Eigen::Matrix<interval,Eigen::Dynamic,Eigen::Dynamic> iMatrix;
	typedef Eigen::MatrixXd pMatrix;
	typedef Eigen::VectorXd pVector;

	double midpoint(interval& I);

	void midpoint(pMatrix& M, iMatrix& iM);

	double rad(interval& I);

	void rad(pVector& M, iMatrix iM);

}

// #include "systemFunction.cpp"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXld;
typedef Eigen::Matrix<vnodelp::interval,Eigen::Dynamic,Eigen::Dynamic> MatrixXint;

void ginac_function_to_file(GiNaC::ex *e,const int& dim);

namespace mstom {

	Eigen::MatrixXd unitNorm(Eigen::MatrixXd& Min, int m);
		// returns unitNorm columnwise(m=1) or rowwise(m=2)

	Eigen::MatrixXd infNorm(Eigen::MatrixXd& Min, int m);
		// returns infinity norm columnwise(m=1) or  rowwise(m=2)

    class zonotope ;

    // product of matrix with Zonotope
    zonotope operator * (const Eigen::MatrixXd& M, const zonotope& Z);

    // vector<interval> to zonotope
    mstom::zonotope vecIntToZono(std::vector<vnodelp::interval> Kprime);

    class intervalMatrix{
    public:
        MatrixXld lb;
        MatrixXld ub;
        intervalMatrix(){};
        zonotope operator * ( const zonotope& Z) const;

        intervalMatrix operator * (double a) const;

        intervalMatrix operator + (const MatrixXld& M) const;

        intervalMatrix operator + (const intervalMatrix& Mi) const ;
    };
	
    zonotope convexHull(const zonotope& Z1, const Eigen::MatrixXd& eAr);

    zonotope convexHull(const zonotope& Z1, const zonotope& Z2);

    zonotope convexHull(std::vector<zonotope>& stora);

    double factorial(const double& n);

    double compute_epsilon(const Eigen::MatrixXd& A, const double& r, int& p);

    double p_adjust_Er_bound(Eigen::MatrixXd& A, double& r, int& p, double& epsilone);

    void matrix_product(double* M1, double* M2, double* Mresult, unsigned int m, unsigned int n, unsigned int q);

    void sum_matrix(double M1[], double M2[], unsigned int m, unsigned int n);

    void matrix_exponential(const MatrixXld& A, const double r, const int& p, intervalMatrix& Er, std::vector<MatrixXld>& Apower);

    intervalMatrix compute_F(const int& p, const double& r, const MatrixXld& A, const intervalMatrix& Er, const std::vector<MatrixXld>& Apower);

    intervalMatrix compute_F_tilde(const int& p, const double& r, const MatrixXld& A, const intervalMatrix& Er, const std::vector<MatrixXld>& Apower, int isOriginContained);

    intervalMatrix compute_Data_interm(const intervalMatrix& Er, const double& r, const int& p, const MatrixXld& A, const std::vector<MatrixXld>& Apower);

    intervalMatrix IntervalHull(const zonotope& Z);

    zonotope project(const zonotope& Z, const int& a, const int& b);

	std::vector<zonotope> project(const std::vector<zonotope>& Zv, const int& a, const int& b);

	template <typename T>
    std::vector<size_t> sort_indexes(const T &v) ;

    zonotope deletezeros(const zonotope& Z);

    void vertices(const zonotope& Z, Eigen::MatrixXd& p2);

	std::vector<std::pair<double, double>> vertices_pair(const zonotope& Z);

	Eigen::MatrixXd verticesH(const zonotope& Z);
        // vertices for H-representation; same as vertices(), difference is only in the return type

    void H_rep(const zonotope& Z, Eigen::VectorXd& M);
        // H representation for 2D zonotope

    bool isOriginInZonotope(const zonotope& Z);
        // only for 2D zonotopes

    void plot(const zonotope& Z, const int& a1, const int& a2);

	void plot(const std::vector<zonotope>& Zv, const int& a1, const int& a2);
        // a1, a2 : dimensions to plot

    void plotfilled(const std::vector<zonotope>& Zv, const int& a1, const int& a2);
        // a1, a2 : dimensions to plot

	std::vector<double> project(const std::vector<mstom::zonotope>& Zv, const int& a);
        // project on to the dimension a(begins from 1); returns the end points of the line segment

	std::vector<std::pair<double, double>> vertices(const std::vector<double>& ve, const double& tau, const double& k);
		// for plot w.r.t. time. Returns vertices of a rectangle of time width tau
		//k = the time instant

	void plotfilled(const std::vector<std::vector<mstom::zonotope>>& Ztp, const int& a1, const double& tau);
        // a1 : dimension to plot w.r.t. time

    void plot(const std::vector<zonotope>& Zv, const int& a1, const int& a2, bool tb);
        // a1, a2 : dimensions to plot

    void plot(const std::vector<double>& L);

    void plotstore(std::vector<zonotope>& PlotStorage, const zonotope& Z);

    void plotstore(std::vector<zonotope>& PlotStorage, const std::vector<zonotope>& Zv);

    void printVector(const std::vector<double>& v);

	void reduce(zonotope& Z, const int& morder);
		// reduces Z to order = morder, if it is greater than that

	void reduce(std::vector<zonotope>& Zv, const int& morder);

	void wfile(const zonotope& Z, const std::string& str1, const int& flag);
		// Eigen::IOFormat HeavyFmt(Eigen::FullPrecision, 0, " ", "\n", "", " ", "[", "]");
		//Eigen::IOFormat LightFmt(Eigen::StreamPrecision, 0, " ", "\n", "", " ", "[", "]");
		//flag: 1(write), 2(append)

	void wfile(const std::vector<std::vector<mstom::zonotope>>& Zti);
		// to plot with MATLAB

	void wfile_gnuplot(std::vector<std::vector<mstom::zonotope>>& Zti);

	void wfile_time(std::vector<std::vector<mstom::zonotope>>& Zti, int a1, double tau);
		// to plot one dimension vs time with MATLAB

} // end namespace mstom

//#############################################################################
// Derivative Hessian

void computeJacobian(double A[], const Eigen::VectorXd& x_bar, Eigen::VectorXd uin);


template<typename Fu>
void computeJacobian_Lu_array(vnodelp::interval xin[], Fu u[], Eigen::MatrixXd& L,const int dim);
    // using arrays
    // jacobian for L(u) for growth bound

template<typename Fu>
void computeJacobian_Lu_array2(vnodelp::interval xin[], Fu u[], const int dim, int jin, double* LuStore);
    // without eigen::matrix
    //  using arrays
    // jacobian for L(u) for growth bound

template<typename T2>
void compute_J_abs_max(const mstom::intervalMatrix& iM, Eigen::MatrixXd J_abs_max[], T2 u);

// void compute_H(const mstom::intervalMatrix& iM, std::vector<vnodelp::iMatrix>& H, Eigen::VectorXd uin);

Eigen::MatrixXd pMatrix_to_MatrixXd(vnodelp::pMatrix pM);

Eigen::VectorXd pV_to_V(vnodelp::pVector pV);

mstom::zonotope compute_quad(mstom::zonotope Z, std::vector<Eigen::MatrixXd> H_mid);

Eigen::VectorXd maxAbs(const mstom::intervalMatrix& IH);

// Eigen::VectorXd compute_L_Hat1(mstom::zonotope Rtotal1, Eigen::VectorXd x_bar, int state_dim, Eigen::VectorXd uin);

Eigen::VectorXd compute_L_Hat3(std::vector<vnodelp::interval> Kprime, std::vector<vnodelp::interval> uin);
    // global L_hat computation

mstom::zonotope compute_L_Hat2(mstom::zonotope Rtotal1, Eigen::VectorXd x_bar, int state_dim, Eigen::VectorXd u);
    // 2nd L_hat computation method (less interval arithmatic)

// mstom::zonotope compute_Rerr_bar(int state_dim, mstom::intervalMatrix& Data_interm, mstom::zonotope& Rhomt, Eigen::VectorXd x_bar,
		// Eigen::VectorXd f_bar, Eigen::VectorXd u, Eigen::VectorXd& L_hat, int LinErrorMethod, mstom::intervalMatrix& F_tilde,
		// Eigen::VectorXd& L_max, int& nr, double& perfInd);
	// nr tells if split needed
	// updates: nr, perfInd, Rhomt

mstom::zonotope compute_L_hatB(int state_dim, Eigen::VectorXd& x_bar, mstom::zonotope& Z0, mstom::zonotope& exprAX0, double r, Eigen::VectorXd& fAx_bar,double Datab, double Datac, double Datad, int LinErrorMethod, Eigen::VectorXd& L_hat, Eigen::VectorXd u);
    // Guernic Girard

//----------------------------------------------------------------------------------------
//########################################################################################
// Reachable set



void splitz(mstom::zonotope& Z0in, mstom::zonotope& Z01, mstom::zonotope& Z02, int mIndex);
	// splitted sets returned in Z01 and Z02; split along dimension mIndex by first taking interval hull

void splitz2(const mstom::zonotope& Z0in, mstom::zonotope& Z01, mstom::zonotope& Z02,
				int mIndex);

template<class state_type>
Eigen::VectorXd computeM(double tau, state_type lower_left, state_type upper_right, state_type inp_lower_left, state_type inp_upper_right);


int computeKprime(double tau, vnodelp::interval* xin,  vnodelp::interval u[], int dim, vnodelp::interval* Kbprime, int KprimeLimit, vnodelp::interval* x_initial);

template<class Lugbx, class Lugbu>
void computeLu(Lugbx xin, Lugbu uin, Lugbx rin, double tau_in, int dim, int dimInput, double* Lu);

template<class Lugbx, class Lugbu, class F4>
Eigen::MatrixXd LuOverSS(Lugbx& lower_left, Lugbx& upper_right, Lugbu& uin, int dim, int dimInput);

template<class Lugbx, class Lugbu>
Eigen::MatrixXd LuOverSS_array(Lugbx& lower_left, Lugbx& upper_right, Lugbu& uin, int& dim, int& dimInput);

template<class Lugbx, class Lugbu>
void LuOverSS_array2(Lugbx& lower_left, Lugbx& upper_right, Lugbu& uin, int& dim, int& dimInput, int jin, double* LuStore);
    // without eigen::matrix
    // using arrays

// template<class CL>
// double one_iteration(mstom::zonotope Z0, Eigen::VectorXd u, int state_dim, double r, int& p, Eigen::VectorXd L_max,
		// std::vector<mstom::zonotope>& stora, std::vector<mstom::zonotope>& Zti_stora, int& count1, int LinErrorMethod, CL& L_hat_storage, const Eigen::VectorXd& ss_eta,
		// int recur, int morder);

// template<class CL>
// std::vector<mstom::zonotope> one_iteration_s(std::vector<mstom::zonotope> Z0, Eigen::VectorXd u,
				// int state_dim, double r, int& p, Eigen::VectorXd L_max, int& count1,
				// int LinErrorMethod, CL& L_hat_storage, const Eigen::VectorXd& ss_eta,
				// int morder, std::vector<mstom::zonotope>& Zti);

// int ReachableSet(const int dim, const int dimInput, double tau, double rr[], double x[], double uu[], int no_of_steps, int LinErrorMethod, double l_bar, int morder, int taylorTerms, std::vector<std::vector<mstom::zonotope>>& Zti, mstom::zonotope& Z0);



#endif /* REACHABLESET2_H_ */
