/* 					line 1654 system(g++,...)
 * ReachableSet2.cpp
 *
 *  Created on: Apr 15, 2018
 *      Author: MahendraSinghTomar
 */

#ifndef REACHABLESET2CPP_H_
#define REACHABLESET2CPP_H_

#include <dlfcn.h>
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
#include "ReachableSet2.h"

namespace Eigen {
  namespace internal {
    template<typename X, typename S, typename P>
    struct is_convertible<X,boost::numeric::interval<S,P> > {
      enum { value = is_convertible<X,S>::value };
    };

    template<typename S, typename P1, typename P2>
    struct is_convertible<boost::numeric::interval<S,P1>,boost::numeric::interval<S,P2> > {
      enum { value = true };
    };
  }
}

template<typename Tx, typename Tu, typename Txx >
Txx funcLj_system(Tx& x, Tu& u, Txx& xx);

extern std::vector<GiNaC::symbol> xs;    // x_symbol
extern std::vector<GiNaC::symbol> us;
// extern std::vector<ex> f;
// extern int PERFINDS;
int PERFINDS = 1;

namespace SymReach{
	
	int state_dim;
	int input_dim;
	
	typedef boost::numeric::interval<double> interval;

	double sup(interval& I){
		return I.upper();
	}

	double inf(interval& I){
		return I.lower();
	}

	double mag(interval& I){
		return boost::numeric::norm(I);
	}

	typedef Eigen::Matrix<interval,Eigen::Dynamic,Eigen::Dynamic> iMatrix;
	typedef Eigen::MatrixXd pMatrix;
	typedef Eigen::VectorXd pVector;

	double midpoint(interval& I){
		return boost::numeric::median(I);
	}

	void midpoint(pMatrix& M, iMatrix& iM){
		unsigned int r = iM.rows();
		unsigned int c = iM.cols();
		for(unsigned int i=0;i<r;i++)
			for(unsigned int j=0;j<c;j++)
				M(i,j) = midpoint(iM(i,j));
	}

	double rad(interval& I){
		return boost::numeric::width(I);
	}

	void rad(pVector& M, iMatrix iM){
		unsigned int r = iM.rows();
		for(unsigned int i=0;i<r;i++)
			M(i) = rad(iM(0,i));
	}
}


// typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> MatrixXld;
// typedef Eigen::Matrix<SymReach::interval,Eigen::Dynamic,Eigen::Dynamic> iMatrix;

Eigen::IOFormat HeavyFmt2(Eigen::FullPrecision, 0, " ", "\n", "", " ", "[", "]");
Eigen::IOFormat LightFmt(Eigen::StreamPrecision, 0, " ", "\n", "", " ", "[", "]");

// std::basic_stringstream<char>& operator << (std::basic_stringstream<char>& coutte, SymReach::interval intd){
// coutte << "[" << intd.lower() << ", " << intd.upper() << "]" << "\n";
// return coutte;
// }

std::ostream& operator <<(std::ostream& os, SymReach::interval& intv){
	os << "[" << intv.lower() << ", " << intv.upper() << "]";
	return os;
}


namespace SymReach {

	Eigen::MatrixXd unitNorm(Eigen::MatrixXd& Min, int m){
		// returns unitNorm columnwise(m=1) or rowwise(m=2)
		Eigen::MatrixXd M;
		if(m==1)
			M = Min.cwiseAbs().colwise().sum();
		else
			M = Min.cwiseAbs().rowwise().sum();
		return M;
	}

	Eigen::MatrixXd infNorm(Eigen::MatrixXd& Min, int m){
		// returns infinity norm columnwise(m=1) or  rowwise(m=2)
		Eigen::MatrixXd M;
		if(m==1)
			M = Min.cwiseAbs().colwise().maxCoeff();
		else
			M = Min.cwiseAbs().rowwise().maxCoeff();
		return M;
	}

    class zonotope {
    public:
        //int n;
        Eigen::MatrixXd generators; //Each column represents a generator
        Eigen::VectorXd centre;
        zonotope (){};
        zonotope (const Eigen::VectorXd& c, const Eigen::VectorXd& eta){
            centre = c;
            unsigned int dim = eta.rows();
            generators = Eigen::MatrixXd::Zero(dim,dim);
            for (unsigned int i=0; i<dim; i++){
                generators(i,i) = eta(i)/2;
            }
        };

        template<class state_type>
        zonotope (state_type ct, state_type etat, int dim){
            Eigen::VectorXd c(dim), eta(dim);
            for(unsigned int i=0;i<dim;i++)
            {
                c(i) = ct[i];
                eta(i) = etat[i];
            }
            centre = c;
            generators = Eigen::MatrixXd::Zero(dim,dim);
            for (unsigned int i=0; i<dim; i++){
                generators(i,i) = eta(i)/2;
            }
        };


        //Minkowskii sum
        zonotope operator + (const zonotope& Zb){
            zonotope Ztemp;
            Ztemp.centre = centre + Zb.centre;
            Ztemp.generators = Eigen::MatrixXd::Zero(generators.rows(), generators.cols() + Zb.generators.cols());
            Ztemp.generators.block(0,0,generators.rows(),generators.cols()) = generators;
            Ztemp.generators.block(0,generators.cols(),Zb.generators.rows(),Zb.generators.cols()) = Zb.generators;
            return Ztemp;
        }

        zonotope operator + (Eigen::VectorXd& V){
            // zonotope + vector
            zonotope Ztemp;
            Ztemp.centre = centre + V;
            Ztemp.generators = generators;
            return Ztemp;
        }

        zonotope operator - (Eigen::VectorXd& V){
            // zonotope - vector
            zonotope Ztemp;
            Ztemp.centre = centre - V;
            Ztemp.generators = generators;
            return Ztemp;
        }

        zonotope operator - (const zonotope& Zb){
            zonotope Ztemp;
            Ztemp.centre = centre - Zb.centre;
            Ztemp.generators = Eigen::MatrixXd::Zero(generators.rows(), generators.cols() + Zb.generators.cols());
            Ztemp.generators.block(0,0,generators.rows(),generators.cols()) = generators;
            Ztemp.generators.block(0,generators.cols(),Zb.generators.rows(),Zb.generators.cols()) = (-1 * Zb.generators);
            return Ztemp;
        }

        zonotope operator * (double a){
            // zonotope * scalar
            zonotope Ztemp;
            Ztemp.centre = centre * a;
            Ztemp.generators = generators * a;
            return Ztemp;
        }

        void display(){
            unsigned int r = generators.rows();
            unsigned int gc = generators.cols();
            Eigen::MatrixXd M(r,1+gc);
            M.block(0,0,r,1) = centre;
            M.block(0,1,r,gc) = generators;
            std::cout<< "no. of gen = " << gc << "\n"<< M.format(HeavyFmt2) << std::endl;        }
    };

    // product of matrix with Zonotope
    zonotope operator * (const Eigen::MatrixXd& M, const zonotope& Z){
        zonotope Ztemp;
        Ztemp.centre = M * Z.centre;
        Ztemp.generators = M * Z.generators;
        return Ztemp;
    }

    // // vector<interval> to zonotope
    // SymReach::zonotope vecIntToZono(std::vector<SymReach::interval> Kprime){
        // unsigned int dim = Kprime.size();
        // Eigen::VectorXd c(dim), lb(dim), ub(dim);
        // for(unsigned int i=0;i<dim;i++)
        // {
            // lb(i) = SymReach::inf(Kprime[i]);
            // ub(i) = SymReach::sup(Kprime[i]);
        // }
        // c = (lb + ub) * 0.5;
        // SymReach:: zonotope Z(c, ub-lb);
        // return Z;
    // }

    // class intervalMatrix{
    // public:
        // MatrixXld lb;
        // MatrixXld ub;
        // intervalMatrix(){};
        // zonotope operator * ( const zonotope& Z){    //interval matrix map for zonotope
            // unsigned int row =  Z.generators.rows();
            // unsigned int col =  Z.generators.cols();
            // zonotope Ztemp;
            // MatrixXld M_tilde = (lb + ub)/2;
            // MatrixXld M_hat = ub - M_tilde;
            // Ztemp.centre = M_tilde * Z.centre;
            // Ztemp.generators = Eigen::MatrixXd::Zero(row, col+row);
            // Ztemp.generators.block(0,0,row,col) = M_tilde * Z.generators;
            // MatrixXld v = MatrixXld::Zero(row,row);
            // for (unsigned int i=0; i<row; i++){
                // v(i,i) = M_hat.row(i) * (Z.centre.cwiseAbs() + Z.generators.cwiseAbs().rowwise().sum() ) ;
            // }
            // Ztemp.generators.block(0,col,row,row) = v;
            // return Ztemp;
        // }

        // intervalMatrix operator * (double a){
            // intervalMatrix Mitemp;
            // MatrixXld temp1 = lb*a;
            // MatrixXld temp2 = ub*a;
            // Mitemp.lb = temp1.array().min(temp2.array()).matrix();  //converted to array for min, then reconverted to matrix
            // Mitemp.ub = temp1.array().max(temp2.array()).matrix();
            // return Mitemp;
        // }

        // intervalMatrix operator + (const MatrixXld& M){
            // intervalMatrix Mitemp;
            // Mitemp.lb = lb + M;
            // Mitemp.ub = ub + M;
            // return Mitemp;
        // }

        // intervalMatrix operator + (const intervalMatrix& Mi){
            // intervalMatrix Mitemp;
            // Mitemp.lb = lb + Mi.lb;
            // Mitemp.ub = ub + Mi.ub;
            // return Mitemp;
        // }

    // };
	
	// class intervalMatrix{
    // public:
        // MatrixXld lb;
        // MatrixXld ub;
        // intervalMatrix(){};
        // zonotope operator * ( const zonotope& Z);

        // intervalMatrix operator * (double a);

        // intervalMatrix operator + (const MatrixXld& M);

        // intervalMatrix operator + (const intervalMatrix& Mi);
    // };
	
	// intervalMatrix::intervalMatrix(){};

	zonotope intervalMatrix::operator * ( const zonotope& Z) const{    //interval matrix map for zonotope
            unsigned int row =  Z.generators.rows();
            unsigned int col =  Z.generators.cols();
            zonotope Ztemp;
            MatrixXld M_tilde = (lb + ub)/2;
            MatrixXld M_hat = ub - M_tilde;
            Ztemp.centre = M_tilde * Z.centre;
            Ztemp.generators = Eigen::MatrixXd::Zero(row, col+row);
            Ztemp.generators.block(0,0,row,col) = M_tilde * Z.generators;
            MatrixXld v = MatrixXld::Zero(row,row);
            for (unsigned int i=0; i<row; i++){
                v(i,i) = M_hat.row(i) * (Z.centre.cwiseAbs() + Z.generators.cwiseAbs().rowwise().sum() ) ;
            }
            Ztemp.generators.block(0,col,row,row) = v;
            return Ztemp;
        }

    intervalMatrix intervalMatrix::operator * (double a) const {
            intervalMatrix Mitemp;
            MatrixXld temp1 = lb*a;
            MatrixXld temp2 = ub*a;
            Mitemp.lb = temp1.array().min(temp2.array()).matrix();  //converted to array for min, then reconverted to matrix
            Mitemp.ub = temp1.array().max(temp2.array()).matrix();
            return Mitemp;
        }

    intervalMatrix intervalMatrix::operator + (const MatrixXld& M) const {
            intervalMatrix Mitemp;
            Mitemp.lb = lb + M;
            Mitemp.ub = ub + M;
            return Mitemp;
        }

    intervalMatrix intervalMatrix::operator + (const intervalMatrix& Mi) const {
            intervalMatrix Mitemp;
            Mitemp.lb = lb + Mi.lb;
            Mitemp.ub = ub + Mi.ub;
            return Mitemp;
        }

    zonotope convexHull(const zonotope& Z1, const Eigen::MatrixXd& eAr){
        // requires initial zonotope and exponential matirx (state-transition matrix)
        unsigned int r1 =  Z1.generators.rows();
        unsigned int c1 =  Z1.generators.cols();
        zonotope Ztemp;
        Ztemp.centre = (Z1.centre + eAr * Z1.centre)/2;
        Ztemp.generators = Eigen::MatrixXd::Zero(r1, 1+2*c1);
        Ztemp.generators.block(0,0,r1,c1) = (Z1.generators + eAr * Z1.generators)/2;
        Ztemp.generators.block(0,c1,r1,1) = (Z1.centre - eAr * Z1.centre)/2;
        Ztemp.generators.block(0,1+c1,r1,c1) = (Z1.generators - eAr * Z1.generators)/2;
        return Ztemp;
    }

    zonotope convexHull(const zonotope& Z1, const zonotope& Z2){
        // when zonotopes may have different number of generators
        zonotope Ztemp;
        unsigned int r = (int)Z1.generators.rows();
        unsigned int c1 = (int)Z1.generators.cols();
        unsigned int c2 = (int)Z2.generators.cols();
        unsigned int c = (c1 > c2 ? c1 : c2);
        Eigen::MatrixXd M = Eigen::MatrixXd::Zero(r,c);
        Ztemp.generators = Eigen::MatrixXd::Zero(r,1+2*c);
        if(c1<c2)
        {
            M.block(0,0,r,c1) = Z1.generators;
            Ztemp.generators.block(0,0,r,c) = 0.5 * (M + Z2.generators);
            Ztemp.generators.block(0,c,r,1) = 0.5 * (Z1.centre - Z2.centre);
            Ztemp.generators.block(0,c+1,r,c) = 0.5 * (M - Z2.generators);
        }
        else
        {
            M.block(0,0,r,c2) = Z2.generators;
            Ztemp.generators.block(0,0,r,c) = 0.5 * (M + Z1.generators);
            Ztemp.generators.block(0,c,r,1) = 0.5 * (Z2.centre - Z1.centre);
            Ztemp.generators.block(0,c+1,r,c) = 0.5 * (M - Z1.generators);
        }
        Ztemp.centre = 0.5*(Z1.centre + Z2.centre);
        return Ztemp;
    }

    zonotope convexHull(std::vector<zonotope>& stora){
        zonotope Z = stora[0];
        for(unsigned int i=1;i<stora.size();i++)
        {
            Z = convexHull(Z, stora[i]);
        }

        return Z;
    }

    double factorial(const double& n){
        return (n==1 || n==0) ? 1 : factorial(n-1) * n;
    }

    double compute_epsilon(const Eigen::MatrixXd& A, const double& r, int& p){
        double norm_rA = (r*A).cwiseAbs().rowwise().sum().maxCoeff();
        double epsilone = norm_rA/(p+2);
        while (epsilone >= 1){
            p += 1;
            epsilone = norm_rA/(p+2);
         }
        return epsilone;
    }

    double p_adjust_Er_bound(Eigen::MatrixXd& A, double& r, int& p, double& epsilone){
        // adjusts p, epsilone
        Eigen::MatrixXd rA = r * A;
        double norm_rA = rA.cwiseAbs().rowwise().sum().maxCoeff();
        double temp = std::pow(norm_rA,p+1) / factorial(p+1);
        double bound = temp / (1-epsilone);
        while(bound > std::pow(10,-12))
        {
            p++;
            epsilone = norm_rA/(p+2);
            temp = temp * norm_rA / (p+1);
            bound = temp / (1-epsilone);
        }
        return bound;
    }

    void matrix_product(double* M1, double* M2, double* Mresult, unsigned int m, unsigned int n, unsigned int q){
        // size M1 = mxn; size M2 = nxq; M1xM2
        double temp;
        for(unsigned int i=0;i<m;i++)
            for(unsigned int j=0;j<q;j++)
            {
                temp = 0;
                for(unsigned int k=0;k<n;k++)
                {
                    temp += M1[i*n+k] * M2[k*q+j];
                }
                Mresult[i*q+j] = temp;
            }
    }

    void sum_matrix(double M1[], double M2[], unsigned int m, unsigned int n){
        // result stored in M1; size = m x n
        for(unsigned int i=0;i<m;i++)
            for(unsigned int j=0;j<n;j++)
                M1[i*n+j] = M1[i*n+j] + M2[i*n+j];
    }

    void matrix_exponential(const MatrixXld& A, const double r, const int& p, intervalMatrix& Er, std::vector<MatrixXld>& Apower){
        // Apower and Er updated
        // unsigned int state_dim = A.rows();
		std::vector<MatrixXld> Apower_abs(p+1);
		MatrixXld M = MatrixXld::Identity(state_dim,state_dim);
		Apower_abs[0] = A.cwiseAbs();
		Apower[0] = A;
		//int fac = 1;
		for(int i=0;i<p;i++)
		{
			Apower[i+1] = Apower[i] * A;
			Apower_abs[i+1] = Apower_abs[i] * Apower_abs[0];
			M = M + Apower_abs[i] * std::pow(r,i+1) / factorial(i+1);
			//fac = fac*(i+2);
		}
		MatrixXld W = (Apower_abs[0]*r).exp() - M;
		W = W.cwiseAbs();
		Er.lb = -W;
		Er.ub = W;

		//std::cout<<"M:\n" << M << std::endl;
		//std::cout <<"exp(Aabs*r):\n" << (Apower_abs[0]*r).exp() << std::endl;
        return ;
    }

    intervalMatrix compute_F(const int& p, const double& r, const MatrixXld& A, const intervalMatrix& Er, const std::vector<MatrixXld>& Apower){
        // unsigned int state_dim = A.rows();
        intervalMatrix Asum;
        Asum.ub = MatrixXld::Zero(state_dim, state_dim);
        Asum.lb = MatrixXld::Zero(state_dim, state_dim);
		//int fac = 2;
		for (double i=2; i<=p; i++)
		{
             double data = (std::pow(i,-i/(i-1)) - std::pow(i,-1/(i-1))) * std::pow(r,i) / factorial(i);
			 //fac = fac * (i+1);
			 for(unsigned int j=0;j<state_dim;j++)
				 for(unsigned int k=0;k<state_dim;k++)
					 if(Apower[i-1](j,k)<0)
						 Asum.ub(j,k) = Asum.ub(j,k) + data * Apower[i-1](j,k);
					 else
						 Asum.lb(j,k) = Asum.lb(j,k) + data * Apower[i-1](j,k);
		}
		Asum = Er + Asum;
        // Eigen::MatrixXd temp(state_dim,state_dim);
        // for (int i=2; i<=p; i++){
              // double data = (pow(i,-i/(i-1)) - pow(i,-1/(i-1)));
            // for(int i2=0;i2<state_dim;i2++)
                // for(int j=0;j<state_dim;j++)
                    // temp(i2,j) = data * Ar_powers_fac[i*state_dim*state_dim + i2*state_dim + j];

            // for(int i1=0; i1<A.rows(); i1++){
                // for(int i2=0; i2<A.rows(); i2++){
                    // if (temp(i1,i2) < 0)
                        // Ftemp.lb(i1,i2) += temp(i1,i2);
                    // else
                        // Ftemp.ub(i1,i2) += temp(i1,i2);
                // }
            // }
        // }
        // Ftemp = Ftemp + Er;
        return Asum;
    }

    intervalMatrix compute_F_tilde(const int& p, const double& r, const MatrixXld& A, const intervalMatrix& Er, const std::vector<MatrixXld>& Apower, int isOriginContained=0){
        // unsigned int state_dim = A.rows();
        intervalMatrix Asum;
        Asum.ub = MatrixXld::Zero(state_dim,state_dim);
        Asum.lb = MatrixXld::Zero(state_dim,state_dim);
		for(double i=2;i<= p+1;i++)
		{
			double data = (std::pow(i,-i/(i-1)) - std::pow(i,-1/(i-1))) * std::pow(r,i)/factorial(i);
			for(unsigned int j=0;j<state_dim;j++)
				for(unsigned int k=0;k<state_dim;k++)
					if(Apower[i-1-1](j,k)<0)
						Asum.ub(j,k) += data * Apower[i-1-1](j,k);
					else
						Asum.lb(j,k) += data * Apower[i-1-1](j,k);
		}
		Asum = Er * (r / A.cwiseAbs().rowwise().sum().maxCoeff()) + Asum;

        // if(!isOriginContained)
        // {
            // Eigen::MatrixXd temp(state_dim,state_dim);
            // for (int i=2; i<=p; i++)
			// {
                   // double data = (pow(i,-i/(i-1)) - pow(i,-1/(i-1)));
                // for(int i2=0;i2<state_dim;i2++)
                    // for(int j=0;j<state_dim;j++)
                        // temp(i2,j) = data * Ar_powers_fac[(i-1)*state_dim*state_dim + i2*state_dim + j] * r / i;

                // for(int i1=0; i1<A.rows(); i1++){
                    // for(int i2=0; i2<A.rows(); i2++){
                        // if (temp(i1,i2) < 0)
                            // Ftemp.lb(i1,i2) += temp(i1,i2);
                        // else
                            // Ftemp.ub(i1,i2) += temp(i1,i2);
                    // }
                // }
            // }
            // Eigen::VectorXd temp2 = A.cwiseAbs().rowwise().sum();
            // Ftemp = Ftemp + Er * pow(temp2.maxCoeff(),-1);
        // }
         return Asum;
    }

    intervalMatrix compute_Data_interm(const intervalMatrix& Er, const double& r, const int& p, const MatrixXld& A, const std::vector<MatrixXld>& Apower){
        // unsigned int state_dim = A.rows();
		MatrixXld Asum = r * MatrixXld::Identity(state_dim,state_dim);
		//int fac = 2;
		for(int i=1;i<=p;i++)
		{
			Asum = Asum + Apower[i-1] * std::pow(r,i+1) / factorial(i+1);
			//fac = fac * (i+2);
		}
		intervalMatrix eAtInt = (Er * r) + Asum;
		return eAtInt;
    }

    intervalMatrix IntervalHull(const zonotope& Z)
    {
        intervalMatrix iM;
        iM.lb = Z.centre - Z.generators.cwiseAbs().rowwise().sum();
        iM.ub = Z.centre + Z.generators.cwiseAbs().rowwise().sum();
        return iM;
    }

    zonotope project(const zonotope& Z, const int& a, const int& b){
        // project on to the dimensions a and b
        zonotope Zp;
        Eigen::Vector2d c;
        Eigen::MatrixXd M(2,Z.generators.cols());
        c(0) = Z.centre(a-1);
        c(1) = Z.centre(b-1);
        M.row(0) = Z.generators.row(a-1);
        M.row(1) = Z.generators.row(b-1);
        Zp.centre = c;
        Zp.generators = M;
        return Zp;
    }

	std::vector<zonotope> project(const std::vector<zonotope>& Zv, const int& a, const int& b){
		unsigned int num = Zv.size();
		std::vector<zonotope> Z(num);
		for(unsigned int i=0;i<num;i++)
			Z[i] = project(Zv[i], a, b);
		return Z;
	}

	template <typename T>
    std::vector<size_t> sort_indexes(const T &v) {
// https://stackoverflow.com/questions/1577475/c-sorting-and-keeping-track-of-indexes?rq=1
        // sorts std::vector and Eigen::VectorXd
		// initialize original index locations
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);

        // sort indexes based on comparing values in v
        std::sort(idx.begin(), idx.end(),
             [&v](size_t i1, size_t i2) {return v[i1] < v[i2];});

        return idx;
    }

    zonotope deletezeros(const zonotope& Z){
        // delete zero generators
        zonotope Zd;
        Eigen::MatrixXd vtemp = Z.generators.cwiseAbs().colwise().sum();
        int ncol = (vtemp.array() == 0).select(vtemp, 1).sum();
        Eigen::MatrixXd M(Z.generators.rows(),ncol);
        int index = 0;
        for(int i=0;i<vtemp.cols();i++)
        {
            if(vtemp(i) != 0)
            {
                M.col(index) = Z.generators.col(i);
                index++;
            }
        }
        Zd.centre = Z.centre;
        Zd.generators = M;
        return Zd;
    }

    void vertices(const zonotope& Z, Eigen::MatrixXd& p2){
        // std::vector<std::pair<double, double>> v;
        const Eigen::VectorXd& c = Z.centre;
        const Eigen::MatrixXd& g = Z.generators;
        int n = g.cols();   // number of generators
        double xmax = g.row(0).cwiseAbs().sum();
        double ymax = g.row(1).cwiseAbs().sum();
        Eigen::MatrixXd gup = g;
        for(int i=0;i<g.cols();i++)
        {
            if(g(1,i)<0)    // if 2nd column of g is negative, then reverse the vector to make it point up
                gup.col(i) = -1 * g.col(i);
        }
        std::vector<double> angles;
        for(int i =0;i<n;i++)
        {
            double angletemp = atan2(gup(1,i), gup(0,i));
            if(angletemp<0)
                angletemp += 2 * M_PI;
            angles.push_back(angletemp);
        }
        std::vector<size_t> sortIndexes = sort_indexes(angles);
        Eigen::MatrixXd p = Eigen::MatrixXd::Zero(2, n+1);
        for(int i=0;i<n;i++)
        {
            p.col(i+1) = p.col(i) + 2 * gup.col(sortIndexes[i]);
        }
        p.row(0) = p.row(0) + Eigen::MatrixXd::Constant(1,p.cols(), (xmax - p.row(0).maxCoeff()));
        p.row(1) = p.row(1) - Eigen::MatrixXd::Constant(1,p.cols(), ymax);
        //Eigen::MatrixXd p2(2, (p.cols()*2));
        p2.row(0).head(p.cols()) = p.row(0);
        p2.row(0).tail(p.cols()) = Eigen::MatrixXd::Constant(1,p.cols(), (p.row(0).tail(1) + p.row(0).head(1))(0)) - p.row(0);
        p2.row(1).head(p.cols()) = p.row(1);
        p2.row(1).tail(p.cols()) = Eigen::MatrixXd::Constant(1,p.cols(), (p.row(1).tail(1)+p.row(1).head(1))(0)) - p.row(1);
        p2.row(0) += Eigen::MatrixXd::Constant(1,p2.cols(), c(0));
        p2.row(1) += Eigen::MatrixXd::Constant(1,p2.cols(), c(1));
        return ;
    }

	std::vector<std::pair<double, double>> vertices_pair(const zonotope& Z){
		int n = Z.generators.cols();   // number of generators
		int p2cols = (n+1)*2;
		Eigen::MatrixXd p2(2, p2cols);
		vertices(Z, p2);
		std::vector<std::pair<double, double>> v(p2cols);
		for(int i=0;i<p2cols;i++)
            v[i] = (std::make_pair(p2(0,i),p2(1,i)));
		return v;
	}

	Eigen::MatrixXd verticesH(const zonotope& Z){
        // vertices for H-representation; same as vertices(), difference is only in the return type
        const Eigen::VectorXd& c = Z.centre;
        const Eigen::MatrixXd& g = Z.generators;
        int n = g.cols();   // number of generators
        double xmax = g.row(0).cwiseAbs().sum();
        double ymax = g.row(1).cwiseAbs().sum();
        Eigen::MatrixXd gup = g;
        for(int i=0;i<g.cols();i++)
        {
            if(g(1,i)<0)    // if 2nd column of g is negative, then reverse the vector to make it point up
                gup.col(i) = -1 * g.col(i);
        }
        std::vector<double> angles;
        for(int i =0;i<n;i++)
        {
            double angletemp = atan2(gup(1,i), gup(0,i));
            if(angletemp<0)
                angletemp += 2 * M_PI;
            angles.push_back(angletemp);
        }
        std::vector<size_t> sortIndexes = sort_indexes(angles);
        Eigen::MatrixXd p = Eigen::MatrixXd::Zero(2, n+1);
        for(int i=0;i<n;i++)
        {
            p.col(i+1) = p.col(i) + 2 * gup.col(sortIndexes[i]);
        }
        p.row(0) = p.row(0) + Eigen::MatrixXd::Constant(1,p.cols(), (xmax - p.row(0).maxCoeff()));
        p.row(1) = p.row(1) - Eigen::MatrixXd::Constant(1,p.cols(), ymax);
        Eigen::MatrixXd p2(2, (p.cols()*2));
        p2.row(0).head(p.cols()) = p.row(0);
        p2.row(0).tail(p.cols()) = Eigen::MatrixXd::Constant(1,p.cols(), (p.row(0).tail(1) + p.row(0).head(1))(0)) - p.row(0);
        p2.row(1).head(p.cols()) = p.row(1);
        p2.row(1).tail(p.cols()) = Eigen::MatrixXd::Constant(1,p.cols(), (p.row(1).tail(1)+p.row(1).head(1))(0)) - p.row(1);
        p2.row(0) += Eigen::MatrixXd::Constant(1,p2.cols(), c(0));
        p2.row(1) += Eigen::MatrixXd::Constant(1,p2.cols(), c(1));
        // vertex in column
        return p2;
    }

    void H_rep(const zonotope& Z, Eigen::VectorXd& M){
        // H representation for 2D zonotope
        if(Z.centre.rows() > 2)
            std::cout << "Please use a 2D zonotope\n";
        Eigen::MatrixXd v = verticesH(Z);
        int n = v.cols();   // number of vertices
        Eigen::MatrixXd temp(v.rows(),v.cols());
        temp.block(0,0,v.rows(),v.cols()-1) = v.block(0,1,v.rows(),v.cols()-1);
        temp.col(v.cols()-1) = v.col(0);
        temp = temp - v;    // column as edge vector
        Eigen::Matrix2d temp2;
        temp2 << 0, 1, -1, 0;
        temp = temp2 * temp;    // perpendicular vectors
        Eigen::VectorXd Ma(n);
        for(int i=0;i<n;i++)
        {
            double dd = temp.col(i).transpose() * v.col(i);
            int j = ((i+n/2)< n)? (i+n/2) : (i+n/2-n);
            if(temp.col(i).transpose() * v.col(j) < dd)
            {
                  Ma(i) = dd;
            }
            else
            {
                  Ma(i) = -dd;
            }
        }
        M = Ma;
    }

    bool isOriginInZonotope(const zonotope& Z){
        // only for 2D zonotopes
        Eigen::VectorXd M;
        H_rep(Z,M); // H representation; last row of M stores d; h.x < d inside the Z
        // for origin only sign of d need to be checked
        bool chk = (M.minCoeff()>=0);
        return chk;
    }

    void plot(const zonotope& Z, const int& a1, const int& a2){
        Gnuplot gp;
        //gp << "set terminal lua\n";
        std::vector<std::pair<double, double>> v = vertices_pair(project(Z, a1, a2));
        //gp << "set output 'my_graph_1.png'\n";
        //gp << "set xrange [-2:2]\nset yrange [-2:2]\n";
        gp << "plot '-' with lines title 'cubic'\n";
        gp.send1d(v);
    }

	void plot(const std::vector<zonotope>& Zv, const int& a1, const int& a2){
        // a1, a2 : dimensions to plot
        Gnuplot gp;
        gp << "set grid\n";
        //gp << "set output 'my_graph_1.png'\n";
        gp << "unset key\n";
		int numbe = Zv.size();
        for(int i=0;i<numbe;i++)
        {
            if(i==0)
                gp << "plot '-' with lines lc rgb \"blue\"";
            // else if(i==numbe-1)
				// gp << "'-' with filledcurves closed lc rgb \"white\"";
			else
                gp << "'-' with lines lc rgb \"blue\"";
            if(i==numbe-1)
                // gp << "\n";
			gp << std::endl;
            else
                gp << ", ";
        }
        for(int i=0;i<numbe;i++)
        {
            std::vector<std::pair<double, double>> v = vertices_pair(project(Zv[i], a1, a2));
            gp.send1d(v);
        }
    }

    void plotfilled(const std::vector<zonotope>& Zv, const int& a1, const int& a2){
        // a1, a2 : dimensions to plot
        Gnuplot gp;
        gp << "set grid\n";
        //gp << "set output 'my_graph_1.png'\n";
        gp << "unset key\n";

		gp <<"set style fill noborder\n";
		gp << "set style function filledcurves closed\n";
		gp << "set style fs solid\n";
		int numbe = Zv.size();
        for(int i=0;i<numbe;i++)
        {
            if(i==0)
                gp << "plot '-' with filledcurves closed lc rgb \"grey\"";
            else if(i==numbe-1)
				gp << "'-' with filledcurves closed lc rgb \"white\"";
			else
                gp << "'-' with filledcurves closed lc rgb \"grey\"";
            if(i==numbe-1)
                gp << "\n";
            else
                gp << ", ";
        }
        for(int i=0;i<numbe;i++)
        {
            std::vector<std::pair<double, double>> v = vertices_pair(project(Zv[i], a1, a2));
            gp.send1d(v);
        }
    }

	std::vector<double> project(const std::vector<SymReach::zonotope>& Zv, const int& a){
        // project on to the dimension a(begins from 1); returns the end points of the line segment
        //int dim = Zv[0].centre.rows();
		std::vector<double> v(2), vs(2);
		double c;
		for(unsigned int i=0;i<Zv.size();i++)
		{
			Eigen::VectorXd M(Zv[i].generators.cols());
			M = Zv[i].generators.row(a-1).transpose();
			//M.row(1) = Z.generators.row(b-1);
			c = Zv[i].centre(a-1);
			vs[0] = c + M.cwiseAbs().sum();	// ub
			vs[1] = c - M.cwiseAbs().sum();	// lb
			if(i==0)
			{
				v[0] = vs[0];
				v[1] = vs[1];
			}
			else
			{
				v[0] = (vs[0] > v[0]) ? (vs[0]) : (v[0]);
				v[1] = (vs[1] < v[1]) ? (vs[1]) : (v[1]);
			}
		}
        return v;
    }

	std::vector<std::pair<double, double>> vertices(const std::vector<double>& ve, const double& tau, const double& k){
		// for plot w.r.t. time. Returns vertices of a rectangle of time width tau
		//k = the time instant
		// double v1 = Z.centre + Z.generators.cwiseAbs().sum();
		// double v2 = Z.centre - Z.generators.cwiseAbs().sum();
		std::vector<std::pair<double, double>> v(4);
		v[0] = std::make_pair(k*tau, ve[0]);
		v[1] = std::make_pair((k+1)*tau, ve[0]);
		v[2] = std::make_pair((k+1)*tau, ve[1]);
		v[3] = std::make_pair(k*tau, ve[1]);
		return v;
	}

	void plotfilled(const std::vector<std::vector<SymReach::zonotope>>& Ztp, const int& a1, const double& tau){
        // a1 : dimension to plot w.r.t. time
        Gnuplot gp;
        gp << "set grid\n";
        //gp << "set output 'my_graph_1.png'\n";
        gp << "unset key\n";

		gp <<"set style fill noborder\n";
		gp << "set style function filledcurves closed\n";
		gp << "set style fs solid\n";
		int steps = Ztp.size();
        for(int i=0;i<steps;i++)
        {
            if(i==0)
                gp << "plot '-' with filledcurves closed lc rgb \"grey\"";
            // else if(i==steps-1)
				// gp << "'-' with filledcurves closed lc rgb \"white\"";
			else
                gp << "'-' with filledcurves closed lc rgb \"grey\"";
            if(i==steps-1)
                gp << "\n";
            else
                gp << ", ";
        }
        for(double i=0;i<steps;i++)
        {
            std::vector<std::pair<double, double>> v = vertices(project(Ztp[i], a1), tau, i);
            gp.send1d(v);
        }
    }

    void plot(const std::vector<zonotope>& Zv, const int& a1, const int& a2, bool tb){
        // a1, a2 : dimensions to plot
        Gnuplot gp;
        gp << "set grid\n";

        //gp << "set output 'my_graph_1.png'\n";
        for(unsigned int i=0;i<Zv.size();i++)
        {
            if(i==0)
                gp << "plot '-' with lines title " << ((tb)? "'True'":"'False'");
            else
                gp << "'-' with lines";
            if(i==Zv.size()-1)
                gp << "\n";
            else
                gp << ", ";
        }
        for(unsigned int i=0;i<Zv.size();i++)
        {
            std::vector<std::pair<double, double>> v = vertices_pair(project(Zv[i], a1, a2));
            gp.send1d(v);
        }
    }

    void plot(const std::vector<double>& L){
        Gnuplot gp;
        gp << "set output 'my_graph.png'\n";
        gp << "plot '-' with points\n";
        gp.send1d(L);
    }

    void plotstore(std::vector<zonotope>& PlotStorage, const zonotope& Z){
        PlotStorage.push_back(Z);
    }

    void plotstore(std::vector<zonotope>& PlotStorage, const std::vector<zonotope>& Zv){
        for(unsigned int i=0;i<Zv.size();i++)
        {
            PlotStorage.push_back(Zv[i]);
        }
    }

    void printVector(const std::vector<double>& v){
        std::cout<< "The vector is: \n";
        for(unsigned int i=0;i<v.size();i++)
            std::cout << v[i] << ", ";
        std::cout << std::endl;
    }

// std::vector<zonotope> PlotStorage2;

	void reduce(zonotope& Z, const int& morder){
		// reduces Z to order = morder, if it is greater than that
		int dim = Z.centre.rows();
		if(Z.generators.cols() <= (morder * dim))
			return;
		Eigen::MatrixXd Gred(dim, dim * morder);
		Eigen::VectorXd h = (unitNorm(Z.generators,1) - infNorm(Z.generators,1)).transpose();
		std::vector<size_t> idx = sort_indexes(h);
		int nUnreduced = dim * (morder-1);
		int nReduced = Z.generators.cols() - nUnreduced;
		Eigen::VectorXd d_ihr(dim);
		d_ihr = Z.generators.col(idx[0]).cwiseAbs();
		for(int i=1;i<nReduced;i++)
			d_ihr += Z.generators.col(idx[i]).cwiseAbs();
		for(unsigned int i=nReduced; i<idx.size(); i++)
			Gred.col(i-nReduced) = Z.generators.col(idx[i]);	//
		Gred.block(0,nUnreduced,dim,dim) = d_ihr.asDiagonal();
		Z.generators = Gred;
		return;
	}

	void reduce(std::vector<zonotope>& Zv, const int& morder){
		for(unsigned int i=0;i<Zv.size();i++)
		{
			reduce(Zv[i], morder);
		}
	}

	void wfile(const zonotope& Z, const std::string& str1, const int& flag){
		// Eigen::IOFormat HeavyFmt2(Eigen::FullPrecision, 0, " ", "\n", "", " ", "[", "]");
		//Eigen::IOFormat LightFmt(Eigen::StreamPrecision, 0, " ", "\n", "", " ", "[", "]");
		//flag: 1(write), 2(append)
		std::ofstream myfile;
		if(flag==1)
			myfile.open("datastorage.txt");
		else
			myfile.open("datastorage.txt",std::ios::app);
		int r = Z.generators.rows();
        int gc = Z.generators.cols();
        Eigen::MatrixXd M(r,1+gc);
        M.block(0,0,r,1) = Z.centre;
        M.block(0,1,r,gc) = Z.generators;
		myfile << str1 << "=";
        myfile << "zonotope("<< M.format(LightFmt) << ");" << std::endl;
		// myfile << "no_of_gen = " << gc << ";" << std::endl;
		myfile.close();
		return;
	}

	void wfile(const std::vector<std::vector<SymReach::zonotope>>& Zti){
		// to plot with MATLAB
		int flag;	// 1(write), 2(append)
		for(unsigned int i=0;i<Zti.size();i++)
		{
			for(unsigned int j=0;j<Zti[i].size();j++)
			{
				if(i==0 && j==0)
					flag = 1;
				else
					flag = 2;
				std::string stri = "Zti{" + std::to_string(i+1) + "}{" + std::to_string(j+1) + "}";
				wfile(Zti[i][j], stri, flag);
			}
		}
	}

	void wfile_gnuplot(std::vector<std::vector<SymReach::zonotope>>& Zti){
		Eigen::IOFormat GnuplotFmt(Eigen::StreamPrecision, 0, " ", "\n", "", "", "", "\n");
		// precision,flags,coeffSeparator,rowSeparator,rowPrefix,rowSuffix,matPrefix,matSuffix
		std::ofstream myfile;
		myfile.open("datastorage_gnuplot.txt");
		myfile << "# x y\n";
		for(unsigned int i=0;i<Zti.size();i++)
		{
			for(unsigned int j=0;j<Zti[i].size();j++)
			{
				unsigned int n = Zti[i][j].generators.cols();   // number of generators
				unsigned int p2cols = (n+1)*2;
				Eigen::MatrixXd p2(2, p2cols);
				vertices(Zti[i][j], p2);
				myfile << p2.transpose().format(GnuplotFmt) << "\n";
			}
		}
		myfile.close();
	}

	void wfile_time(std::vector<std::vector<SymReach::zonotope>>& Zti, int a1, double tau){
		// to plot one dimension vs time (with MATLAB)
		int flag;	// 1(write), 2(append)
		for(double i=0;i<Zti.size();i++)
		{
			if(i==0)
				flag = 1;
			else
				flag = 2;
			std::string stri = "Zti1D{" + std::to_string(i+1) + "}" ;
			std::vector<double> v = project(Zti[i], a1);
			//vertices(project(Zti[i], a1), tau, i);
			Eigen::Vector2d c;
			c(0) = (i+1)*tau - 0.5*tau;
			c(1) = 0.5*(v[0]+v[1]);
			Eigen::Vector2d g;
			g(0) = 0.5*tau;
			g(1) = 0.5*(v[0]-v[1]);
			zonotope Z;
			Z.centre = c;
			Z.generators = g.asDiagonal();
			wfile(Z, stri, flag);
		}

	}

} // end namespace SymReach

//#############################################################################
// Derivative Hessian

// std::vector<SymReach::zonotope> PlotStorage;
int ZorDeltaZ;  // 1 if Z, 0 if deltaZ

namespace SymReach{
    
Eigen::MatrixXd pMatrix_to_MatrixXd(SymReach::pMatrix pM)
{
    int n = pM.size();
    Eigen::MatrixXd M(n,n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            M(i,j) = pM(i,j);
        }
    }
    return M;
}

Eigen::VectorXd pV_to_V(SymReach::pVector pV)
{
    int n = pV.size();
    Eigen::VectorXd V(n);
    for(int i=0;i<n;i++)
        V(i) = pV[i];
    return V;
}

SymReach::zonotope compute_quad(SymReach::zonotope Z, std::vector<Eigen::MatrixXd> H_mid)
{
    int row, gcol; // gcol = number of generators
    row = Z.centre.rows();
    gcol = Z.generators.cols();
    SymReach::zonotope Zq;
    Eigen::MatrixXd Zmat(row, 1+gcol);
    Zmat.block(0,0,row,1) = Z.centre;
    Zmat.block(0,1,row,gcol) = Z.generators;
    Eigen::MatrixXd quadMat, centre(row,1), gen(row,int((gcol*gcol + 3*gcol)*0.5));
    for(unsigned int i =0;i<H_mid.size();i++)
    {
        quadMat = Zmat.transpose() * H_mid[i] * Zmat; // (gcol+1) x (gcol+1)
        centre(i,0) = quadMat(0,0) + 0.5 * quadMat.diagonal().tail(gcol).sum();
        gen.row(i).head(gcol) = 0.5 * quadMat.diagonal().tail(gcol);
        int count =0;
        for(int j=0;j<gcol+1;j++)
        {
            gen.row(i).segment(gcol+count,gcol-j) = quadMat.row(j).segment(j+1,gcol-j)  + quadMat.col(j).segment(j+1,gcol-j).transpose();
            count += gcol - j;
        }
    }
    Zq.centre = centre;
    Zq.generators = gen;
    return SymReach::deletezeros(Zq);
}

Eigen::VectorXd maxAbs(const SymReach::intervalMatrix& IH){
    Eigen::MatrixXd Mtemp2(IH.lb.rows(),2);
    Mtemp2.col(0) = IH.lb;
    Mtemp2.col(1) = IH.ub;
    return Mtemp2.cwiseAbs().rowwise().maxCoeff();
}

void ginac_function_to_file(GiNaC::ex *e, const int& dimHess, const int& dim){
	std::ofstream myfile;
	myfile.open("func_from_file.cpp");
	myfile << "#include <iostream>\n"
	<< "#include <vector>\n"
	<< "#include <boost/numeric/interval.hpp>\n\n" 
	<< "using namespace boost::numeric;\n"
	<< "using namespace interval_lib;\n"
	<< "typedef boost::numeric::interval<double, policies<save_state<rounded_transc_std<double> >, checking_base<double> > > intervalD;\n"
	<< "typedef std::vector<intervalD> vec_interval;\n\n"
	<< "extern \"C\" void func_from_file(const vec_interval& x, const vec_interval& u, vec_interval& Hessi){\n";
	int i2;
	for(int i=0; i<dim; i++)
		for(int j=0;j<dimHess;j++)
			for(int k=0;k<dimHess;k++)
				{
					i2 = i*dimHess*dimHess + j*dimHess + k;
					if(!e[i2].is_zero())
						myfile << "\tHessi[" << i2 << "] = " << GiNaC::csrc_double << e[i2] << ";\n";
				}
	myfile << "}";
	myfile.close();
	return;
	}

typedef void (*func_from_file_t)(const std::vector<SymReach::interval>& , const std::vector<SymReach::interval>& , std::vector<SymReach::interval>& );

func_from_file_t func_from_file;

void compute_H(const SymReach::intervalMatrix& iM, std::vector<SymReach::iMatrix>& H, const SymReach::intervalMatrix& iMu, const int& dimInput)
{
    int nm = iM.lb.rows(); // state_dimension
	int dimHess = nm + dimInput;
	std::vector<SymReach::interval> x(dimHess), u(dimInput);
    for(int j=0;j<nm;j++)
    {
        x[j] = SymReach::interval(iM.lb(j,0), iM.ub(j,0));
    }
	for(int j=0;j<dimInput;j++)
		u[j] = SymReach::interval(iMu.lb(j,0), iMu.ub(j,0));
	std::vector<SymReach::interval> Hessian(nm*dimHess*dimHess, 0);
	func_from_file(x, u, Hessian);
    for(int i=0;i<nm;i++)
    {
        for(int j=0;j<dimHess;j++)
            for(int k=0;k<dimHess;k++)
            {
                H[i](j,k) = Hessian[i*dimHess*dimHess + j*dimHess + k];
            }
	}
}


Eigen::VectorXd compute_L_Hat1(const SymReach::zonotope& Rtotal1, const Eigen::VectorXd& x_bar, const int& state_dim, const int& input_dim, const Eigen::VectorXd& uin, const SymReach::zonotope& delta_u){
	SymReach::intervalMatrix totalIntu;
	int dimHess = input_dim + state_dim;
	SymReach::intervalMatrix RIH = SymReach::IntervalHull(Rtotal1);
	
	// std::cout << "x_bar:\n" << x_bar<<"\n";
	// std::cout <<"RIH.lb, ub\n"<< RIH.lb<<"\n"<<RIH.ub<<"\n";
	
	SymReach::intervalMatrix totalInt = RIH + x_bar;
    Eigen::VectorXd L_hat(state_dim);
    {
        Eigen::VectorXd Gamma(dimHess);
        if (ZorDeltaZ)  //1 if Z, 0 if deltaZ;
		{
            Gamma.head(state_dim) = maxAbs(RIH);
			//std::cout<<"gamma:\n" << Gamma << std::endl;
			//Gamma = (Rtotal1.centre - x_bar).cwiseAbs() + Rtotal1.generators.cwiseAbs().rowwise().sum();
        }
		else
            Gamma = (Rtotal1.centre).cwiseAbs() + Rtotal1.generators.cwiseAbs().rowwise().sum();
		if(input_dim != 0)
		{
			SymReach::intervalMatrix duIH = SymReach::IntervalHull(delta_u);
			totalIntu = duIH + uin;
			Gamma.tail(input_dim) = duIH.ub ; //maxAbs(duIH);

		}
       // Eigen::MatrixXd J_abs_max[state_dim];
        //compute_J_abs_max(totalInt, J_abs_max, uin);
		std::vector<SymReach::iMatrix> H(state_dim, SymReach::iMatrix::Zero(dimHess,dimHess));
		iMatrix Htemp(dimHess, dimHess);
		compute_H(totalInt, H, totalIntu, input_dim);

		iMatrix error(state_dim,1), Gammaint(dimHess,1);
		for(int i=0;i<dimHess;i++)
			Gammaint(i,0) = SymReach::interval(Gamma(i));
        for(int i=0; i<state_dim;i++)
        {
			// error(i,0) = 0.5 * (Gammaint.transpose()* H[i] * Gammaint);
				
			iMatrix ertemp;
			SymReach::interval sctemp(0,0);
			ertemp = (H[i] * Gammaint);
			iMatrix ertemptemp = Gammaint.transpose();
			// ertemp = (ertemptemp * ertemp);
			for(int i=0;i<ertemp.rows();i++)
				sctemp += ertemptemp(0,i) * ertemp(i,0);
			sctemp = 0.5 * sctemp;
			error(i,0) = sctemp;
			//std::cout <<"error: " << ertemp(0,0) << std::endl;
		}
		for(int i=0;i<state_dim;i++)
			L_hat(i) = SymReach::mag(error(i,0));
    }
    return L_hat;
}

SymReach::zonotope compute_L_Hat2(SymReach::zonotope Rtotal1, Eigen::VectorXd x_bar, int state_dim, Eigen::VectorXd u){
    // 2nd L_hat computation method (less interval arithmatic)
    SymReach::intervalMatrix RIH = SymReach::IntervalHull(Rtotal1);
    Eigen::VectorXd L_hat(state_dim);
         std::vector<SymReach::iMatrix> H;
        // compute_H(RIH, H, u);
        int dim = H[1].size();
        std::vector<Eigen::MatrixXd> H_mid, H_rad;
        Eigen::MatrixXd Mtemp(dim,dim);
        SymReach::pMatrix pMtemp(dim,dim);
        // SymReach::sizeM(pMtemp,dim);
        SymReach::pVector pV(dim);
        // SymReach::sizeV(pV,dim);
        for(unsigned int i=0;i<H.size();i++)
        {
            SymReach::midpoint(pMtemp, H[i]);
            H_mid.push_back(pMatrix_to_MatrixXd(pMtemp));
            for(int ii=0;ii<dim;ii++)
            {
                SymReach::rad(pV, H[i].row(ii));
                Mtemp.block(ii,0,1,dim) = pV_to_V(pV).transpose();
            }
            H_rad.push_back(Mtemp);
        }
        SymReach::zonotope Zq = compute_quad(Rtotal1-x_bar, H_mid);
        SymReach::zonotope error_mid;
        error_mid.centre = 0.5 * Zq.centre;
        error_mid.generators = 0.5 * Zq.generators;
        Eigen::VectorXd dz_abs = maxAbs(SymReach::IntervalHull(Rtotal1-x_bar));
        Eigen::VectorXd error_rad(H.size());
        for(unsigned int i=0;i<H.size();i++)
            error_rad(i) = 0.5 * dz_abs.transpose() * H_rad[i] * dz_abs;
        SymReach::zonotope error  = error_mid + SymReach::zonotope(Eigen::VectorXd::Zero(dim), 2*error_rad);
    return error;
}

SymReach::zonotope compute_Rerr_bar(int& state_dim, int& input_dim, SymReach::intervalMatrix& Data_interm, SymReach::zonotope& Rhomt, Eigen::VectorXd x_bar,
		Eigen::VectorXd f_bar, Eigen::VectorXd u, Eigen::VectorXd& L_hat, int LinErrorMethod, SymReach::intervalMatrix& F_tilde,
		Eigen::VectorXd& L_max, int& nr, double& perfInd, SymReach::zonotope& delta_u){
	// nr tells if split needed
	// updates: nr, perfInd, Rhomt
	// SymReach::zonotope F_tilde_f_bar = F_tilde * SymReach::zonotope(f_bar, Eigen::VectorXd::Zero(state_dim));
	SymReach::zonotope Rhom;
    Eigen::VectorXd appliedError;

//    if(L_hat_previous.rows() == 0)
       appliedError  = Eigen::MatrixXd::Constant(state_dim,1,0);
//    else
//        appliedError = L_hat_previous;

    SymReach::zonotope Verror, RerrAE, Rtotal1, ErrorZonotope;
    SymReach::intervalMatrix IMtemp;
    Eigen::VectorXd trueError;
    Verror.centre = Eigen::VectorXd::Zero(state_dim);
    double perfIndCurr = 2;
    while(perfIndCurr > 1)
    {
        Rhom = Rhomt;
        //if((f_bar-appliedError).maxCoeff() > 0 || (f_bar+appliedError).minCoeff() < 0)
            // Rhom = Rhom + F_tilde_f_bar;    // when f_bar + [-L,L] does not contain origin
        Verror.generators = Eigen::MatrixXd(appliedError.asDiagonal());
        RerrAE = Data_interm * Verror;
        Rtotal1 = (RerrAE + Rhom);
        if(LinErrorMethod == 1)
            trueError = compute_L_Hat1(Rtotal1, x_bar, state_dim, input_dim, u, delta_u);
        else
        {
        	 ErrorZonotope = compute_L_Hat2(Rtotal1, x_bar, state_dim,u);// maxAbs(IH(error zonotope))
        	 IMtemp = SymReach::IntervalHull(ErrorZonotope);
        	 trueError = maxAbs(IMtemp);
        }
           // trueError = compute_L_Hat2(Rtotal1, x_bar, state_dim,u);
        perfIndCurr = (trueError.cwiseProduct(appliedError.cwiseInverse())).maxCoeff(); // max(trueError./appliedError)
        perfInd = (trueError.cwiseProduct(L_max.cwiseInverse())).maxCoeff(); // max(trueError./L_max)
        if(perfInd > 1)
        {
        	nr = -2;
        	break;
        }
        // appliedError = 1.02 * trueError;
		appliedError = 1.02 * trueError;
    }
    L_hat = trueError;
   // L_hat_previous = trueError;
    Verror.generators = Eigen::MatrixXd(trueError.asDiagonal());
    SymReach::zonotope Rerror = Data_interm * Verror;

    Rhom = Rhomt;
        if(LinErrorMethod == 1)
        {
            //if((f_bar-trueError).maxCoeff() > 0 || (f_bar+trueError).minCoeff() < 0)
            //    Rhom = Rhom + F_tilde_f_bar;    // when f_bar + [-L,L] does not contain origin
            Rhomt = (Rerror + Rhom) + SymReach::zonotope(x_bar, Eigen::VectorXd::Zero(state_dim));    // R[0,r] returned in Rhomt
        }
        else
        {
            SymReach::zonotope Ztemp = ErrorZonotope + f_bar;
            SymReach::zonotope F_tilde_u_tilde = F_tilde * SymReach::zonotope(Ztemp.centre, Eigen::VectorXd::Zero(state_dim));
            Ztemp.centre = Eigen::VectorXd::Zero(state_dim);
    //        if((f_bar+IMtemp.lb).maxCoeff() > 0 || (f_bar+IMtemp.ub).minCoeff() < 0)
                Rhom = Rhom + F_tilde_u_tilde;    // when f_bar + [-L,L] does not contain origin
            Rerror = Data_interm * Ztemp;
            Rhomt = Rhom + Rerror;
            Rerror = Data_interm * ErrorZonotope;
        }
        return Rerror;
}

//----------------------------------------------------------------------------------------
//########################################################################################
// Reachable set

void splitz(SymReach::zonotope& Z0in, SymReach::zonotope& Z01, SymReach::zonotope& Z02, int mIndex)
{
	// splitted sets returned in Z01 and Z02; split along dimension mIndex by first taking interval hull
	SymReach::intervalMatrix Z0ih = IntervalHull(Z0in);
	SymReach::zonotope Z0;	// Z0 = hyperinterval
	Z0.centre = 0.5*(Z0ih.lb + Z0ih.ub);
	Z0.generators = (0.5*(Z0ih.ub - Z0ih.lb)).asDiagonal();

    Z01.centre = Z0.centre - 0.5 * Z0.generators.col(mIndex);
    Z02.centre = Z0.centre + 0.5 * Z0.generators.col(mIndex);
    Z01.generators = Z0.generators;
    Z02.generators = Z0.generators;
    Z01.generators.col(mIndex) = 0.5 * Z0.generators.col(mIndex);
    Z02.generators.col(mIndex) = 0.5 * Z0.generators.col(mIndex);
}

void splitz2(const SymReach::zonotope& Z0in, SymReach::zonotope& Z01, SymReach::zonotope& Z02,
				int mIndex){
	Z01.centre = Z0in.centre - 0.5 * Z0in.generators.col(mIndex);
	Z01.generators = Z0in.generators;
	Z01.generators.col(mIndex) = 0.5 * Z0in.generators.col(mIndex);
	Z02.centre = Z0in.centre + 0.5 * Z0in.generators.col(mIndex);
	Z02.generators = Z0in.generators;
	Z02.generators.col(mIndex) = 0.5 * Z0in.generators.col(mIndex);
	}



double one_iteration(SymReach::zonotope Z0, Eigen::VectorXd u, double r, int& p, Eigen::VectorXd L_max,
		std::vector<SymReach::zonotope>& stora, std::vector<SymReach::zonotope>& Zti_stora, int& count1, int LinErrorMethod, const Eigen::VectorXd& ss_eta,
		int recur, int morder, GiNaC::ex *JacobianB, GiNaC::ex *HessianB, double ru[],  GiNaC::ex *JacobianBu)
{
    // TicToc timehat;
    count1++;
    Eigen::VectorXd c = Z0.centre;
    // std::vector<double> cv(state_dim) ;  // to pass c as c++ vector to func
    // SymReach::VXdToV(cv, c);
    Eigen::VectorXd x_bar(state_dim);
    x_bar = c + 0.5 * r * funcLj_system(c, u, x_bar);
	
	// SymReach::interval

    MatrixXld A(state_dim,state_dim), B(state_dim, input_dim);

    // double A_array[state_dim*state_dim];

    // computeJacobian(A_array, x_bar, u);
    // for(int i=0;i<state_dim;i++)
        // for(int j=0;j<state_dim;j++)
            // A(i,j) = A_array[i*state_dim+j];

	// computing Jacobian at x_bar
	GiNaC::exmap mb;
	for(int i=0;i<state_dim;i++)
		mb[xs[i]] = x_bar(i);
	for(int i=0;i<input_dim;i++)
		mb[us[i]] = u(i);
	for(int i=0;i<state_dim;i++)
	{
		for(int j=0;j<state_dim;j++)
			A(i,j) = GiNaC::ex_to<GiNaC::numeric>(JacobianB[i*state_dim+j].subs(mb)).to_double();	
		for(int j=0;j<input_dim;j++)
			B(i,j) = GiNaC::ex_to<GiNaC::numeric>(JacobianBu[i*input_dim+j].subs(mb)).to_double();
	}
	Eigen::VectorXd ru_vec(input_dim);
	for(int i=0;i<input_dim;i++)
		ru_vec(i) = ru[i];
	SymReach::zonotope delta_u;
    Eigen::VectorXd f_bar(state_dim);
    funcLj_system(x_bar, u, f_bar);
	
	
	// std::cout <<"A:\n" << A<<"\n";
	// std::cout<<"f_bar:\n"<<f_bar<<"\n";
	
	
    SymReach::zonotope Rtotal_tp;
    Eigen::VectorXd L_hat;

	// std::cout << "x_bar = ";
	// for(int i=0;i<state_dim;i++)
		// std::cout << x_bar(i) << ", ";
	// std::cout << std::endl;

	//std::cout << "Jacobian, A: \n" << A.format(HeavyFmt2)  << std::endl;

    SymReach::intervalMatrix Er;   //E(r)
    SymReach::intervalMatrix Data_interm;  //(Aˆ-1)*(exp(Ar) - I)
    //double epsilone = SymReach::compute_epsilon(A,r,p);
    // int isOriginContained = 0;
    //    Ar_powers_fac; // pow(A*r,i)/factorial(i)

    //double bound = SymReach::p_adjust_Er_bound(A,r,p,epsilone);
    //double A_powers_arr[state_dim*state_dim*(p+1)];
	std::vector<MatrixXld> Apower(p+1);

    SymReach::matrix_exponential(A, r, p, Er, Apower);  //updates Er, Apower
    SymReach::intervalMatrix F = SymReach::compute_F(p, r, A, Er, Apower);
    SymReach::intervalMatrix F_tilde = SymReach::compute_F_tilde(p,r,A,Er, Apower);
    SymReach::zonotope F_tilde_f_bar = F_tilde * SymReach::zonotope(f_bar, Eigen::VectorXd::Zero(state_dim));

    // Data_interm = (Aˆ-1)*(exp(Ar) - I)
    Data_interm = compute_Data_interm(Er, r, p, A, Apower);

    SymReach::zonotope Z0delta = Z0 + SymReach::zonotope(-x_bar,Eigen::VectorXd::Zero(A.rows()));
	SymReach::zonotope Bdelta_uPf_bar;
	if(input_dim == 0)
		Bdelta_uPf_bar = SymReach::zonotope(f_bar, Eigen::VectorXd::Zero(state_dim));
	else
	{
		delta_u = SymReach::zonotope(Eigen::VectorXd::Zero(input_dim), 2*ru_vec);
		Bdelta_uPf_bar = ((B*delta_u) + f_bar);	// B * delta_u + f_bar
	}
	SymReach::zonotope Rtrans = Data_interm * Bdelta_uPf_bar;
    SymReach::zonotope Rhom_tp = (A*r).exp() * Z0delta + Rtrans ;

    SymReach::zonotope Rtotal;
	SymReach::zonotope Rhom = SymReach::convexHull(Z0delta, Rhom_tp) + F * Z0delta + F_tilde * Bdelta_uPf_bar;
		
	// std::cout<<"Rhom:\n";
	// Rhom.display();

	SymReach::reduce(Rhom, morder);
	SymReach::reduce(Rhom_tp, morder);

    ZorDeltaZ = 1;  //1 if Z, 0 if deltaZ
    SymReach::zonotope RV;

    int nr = -1;	// nr = -1 (empty), -2 (split needed)
	double perfInd;

    RV = compute_Rerr_bar(state_dim, input_dim, Data_interm, Rhom, x_bar, f_bar, u, L_hat,LinErrorMethod, F_tilde, L_max, nr, perfInd, delta_u);  // Rerr_bar; L_hat, Rhom updated in the call

    Rtotal_tp = Rhom_tp + RV  + x_bar;
    Rtotal = Rhom;

    if(nr == -1)	// no split needed
    {
		SymReach::reduce(Rtotal_tp, morder);
		SymReach::reduce(Rtotal, morder);
        stora.push_back(Rtotal_tp);
		Zti_stora.push_back(Rtotal);
		// SymReach::plotstore(PlotStorage, Rtotal);
		// std::cout << "stora.size()=" << stora.size() << "\n";
		// std::cout << "Zti_stora.size()=" << Zti_stora.size() << "\n";
    }
	
	// if(recur == 1)
	// {
		// std::cout << "L_hat = " ;
		// for(int i=0;i<state_dim;i++)
			// std::cout << L_hat(i) << ", ";
		// std::cout << std::endl;
	// }
	
    SymReach::zonotope Z01, Z02;
    if(nr == -2 && recur == 1) // split needed. (recur == 0) => (return perfInd, no further splits)
    {
		std::cout << "###########  Split  #############\n" ;

		// int number_split = state_dim;
		int number_split = Z0.generators.cols();

		std::vector<std::vector<SymReach::zonotope>> split_stora1(number_split);
		std::vector<std::vector<SymReach::zonotope>> split_Zti1(number_split);
		std::vector<std::vector<SymReach::zonotope>> split_Zti2(number_split);
		//if(PERFINDS==2)
			std::vector<std::vector<SymReach::zonotope>> split_stora2(number_split);

		std::vector<double> perfInd_split(number_split);
		std::vector<double> perftemp2(number_split);
		std::vector<SymReach::zonotope> Zsplit1(number_split);	// only one set from each split
		std::vector<SymReach::zonotope> Zsplit2(number_split);	// only one set from each split
		for(int i=0;i<number_split;i++)
		{
			// splitz(Z0, Zsplit1[i], Zsplit2[i], i);	// split along state_dim
			splitz2(Z0, Zsplit1[i], Zsplit2[i], i); // split along ith generator
			double perftemp;
			perftemp = one_iteration(Zsplit1[i], u, r, p, L_max, split_stora1[i], split_Zti1[i],
									count1, LinErrorMethod, ss_eta, 0, morder, JacobianB, HessianB, ru, JacobianBu);
			if(PERFINDS==2)
			{
				perftemp2[i] = one_iteration(Zsplit2[i], u, r, p, L_max, split_stora2[i], split_Zti2[i],
											count1, LinErrorMethod, ss_eta, 0, morder, JacobianB, HessianB, ru, JacobianBu);
				perftemp = perftemp * perftemp2[i];
			}
			perfInd_split[i] = perftemp;
		}
		auto r_iterat = std::min_element(perfInd_split.begin(), perfInd_split.end());
		int sel_index = r_iterat - perfInd_split.begin();

		// std::cout << " psel= " << perfInd_split[sel_index]
			// <<", p1= " <<perfInd_split[sel_index]/perftemp2[sel_index]<<", p2= "<<perftemp2[sel_index]<<std::endl;

        Z01 = Zsplit1[sel_index];
		Z02 = Zsplit2[sel_index];
		//plot(Z01,1,2);
		//plot(Z02,1,2);
		
		// std::cout << std::endl << "perftemp2: " ;
		// for(int i=0; i<number_split; i++)
			// std::cout << perftemp2[i] << ", ";
		// std::cout << "---------" << std::endl;
		// std::cin.get();
		
		
		

		if((PERFINDS==1 and perfInd_split[sel_index] < 1) or (PERFINDS==2 && perfInd_split[sel_index]/perftemp2[sel_index] < 1 && perftemp2[sel_index] < 1))
		{
			if(stora.size()==0)
			{
				stora = split_stora1[sel_index];
				Zti_stora = split_Zti1[sel_index];
				if(PERFINDS==2)
				{
					stora.insert(stora.end(), split_stora2[sel_index].begin(), split_stora2[sel_index].end());
					Zti_stora.insert(Zti_stora.end(), split_Zti2[sel_index].begin(), split_Zti2[sel_index].end());
				}

			}
			else
			{
				stora.insert(stora.end(), split_stora1[sel_index].begin(), split_stora1[sel_index].end());
				Zti_stora.insert(Zti_stora.end(), split_Zti1[sel_index].begin(), split_Zti1[sel_index].end());
				if(PERFINDS==2)
				{
					stora.insert(stora.end(), split_stora2[sel_index].begin(), split_stora2[sel_index].end());
					Zti_stora.insert(Zti_stora.end(), split_Zti2[sel_index].begin(), split_Zti2[sel_index].end());
				}

			}
			if(PERFINDS==1)
			{
			std::vector<SymReach::zonotope> split_sto2, split_Ztio2;
			one_iteration(Z02, u, r, p, L_max, split_sto2, split_Ztio2, count1, LinErrorMethod, ss_eta, 1, morder, JacobianB, HessianB, ru, JacobianBu);
			stora.insert(stora.end(), split_sto2.begin(), split_sto2.end());
			Zti_stora.insert(Zti_stora.end(), split_Ztio2.begin(), split_Ztio2.end());
			}
		}
		else
		{
			// to edit: split as we already know that the split is needed.
			one_iteration(Z01, u, r, p, L_max, stora, Zti_stora, count1, LinErrorMethod, ss_eta, 1, morder, JacobianB, HessianB, ru, JacobianBu);
			one_iteration(Z02, u, r, p, L_max, stora, Zti_stora, count1, LinErrorMethod, ss_eta, 1, morder, JacobianB, HessianB, ru, JacobianBu);
		}
    }

    return perfInd;
}

std::vector<SymReach::zonotope> one_iteration_s(std::vector<SymReach::zonotope> Z0, Eigen::VectorXd u,
				double r, int& p, Eigen::VectorXd L_max, int& count1,
				int LinErrorMethod, const Eigen::VectorXd& ss_eta,
				int morder, std::vector<SymReach::zonotope>& Zti, GiNaC::ex *JacobianB, GiNaC::ex *HessianB, double ru[],  GiNaC::ex *JacobianBu)
{
	std::vector<SymReach::zonotope> Zoutput, Ztemp;
	for(unsigned int i=0;i<Z0.size();i++)
	{
		int recur = 1;	// 1 (keep on iterating for reachable set), 0 (return perfInd, no further splits)
		std::vector<SymReach::zonotope> stora;
		std::vector<SymReach::zonotope> Zti_stora;
		one_iteration(Z0[i],u, r,p, L_max, stora, Zti_stora, count1, LinErrorMethod, ss_eta, recur, morder, JacobianB, HessianB, ru, JacobianBu);	// reachable sets returned in stora
		if(i==0)
		{
			Zoutput = stora;
			Zti = Zti_stora;
		}
		else
		{
			Zoutput.insert(Zoutput.end(),stora.begin(),stora.end());
			Zti.insert(Zti.end(),Zti_stora.begin(),Zti_stora.end());
		}
	}
	// std::cout << "Zoutput.size = " << Zoutput.size() << std::endl;
	return Zoutput;
}

int ReachableSet(const int dim, const int dimInput, double tau, double rr[], double x[], double ru[], double uu[], int no_of_steps, int LinErrorMethod, double l_bar, int morder, int taylorTerms, 
std::vector<std::vector<SymReach::zonotope>>& Zti, SymReach::zonotope& Z0)
{
	state_dim = dim;
	input_dim = dimInput;
	
	std::vector<GiNaC::ex> f(dim);
	funcLj_system(xs, us, f);
	GiNaC::ex *JacobianB = new GiNaC::ex[dim*dim];
	for(int i=0;i<dim;i++)
		for(int j=0;j<dim;j++)
		{
			JacobianB[i*dim+j] = f[i].diff(xs[j]);
		}
	GiNaC::ex *JacobianBu;
	if(dimInput != 0)
	{
		JacobianBu = new GiNaC::ex[dim*dimInput];	// Bu 	(n x m)
		for(int i=0;i<dim;i++)
			for(int j=0;j<dimInput;j++)
				JacobianBu[i*dimInput+j] = f[i].diff(us[j]);
	}
	const int dimHess(dim+dimInput);	// Hessian dimension z = [x;u]
	
	GiNaC::ex *JacobianHtemp; // 
	if(dimInput == 0)
		JacobianHtemp = JacobianB;
	else
	{
		for(int i=0;i<dimInput;i++)
			xs.push_back(us[i]);
		// Now xs = zs
		JacobianHtemp = new GiNaC::ex[dim*dimHess];
		for(int i=0;i<dim;i++)
			for(int j=0;j<dim;j++)
				JacobianHtemp[i*dimHess+j] = f[i].diff(xs[j]);
	}
	GiNaC::ex *HessianB = new GiNaC::ex[dim*dimHess*dimHess];
	for(int i=0;i<dim;i++)
		for(int j=0;j<dimHess;j++)
			for(int k=0;k<dimHess;k++)
				HessianB[i*dimHess*dimHess + j*dimHess + k] = JacobianHtemp[i*dimHess+j].diff(xs[k]);
	ginac_function_to_file(HessianB, dimHess, dim);
  // system("clang++ -shared func_from_file.cpp -o func_from_file.so");
    system("g++ -shared func_from_file.cpp -o func_from_file.so");
    
    using std::cout;
    using std::cerr;
    void* handle = dlopen("./func_from_file.so",RTLD_LAZY);
    if(!handle){
        cerr << "Cannot open library: " << dlerror() << '\n';
        return 1;
    }
    // typedef void (*func_from_file_t)(const std::vector<SymReach::interval>& , std::vector<SymReach::interval>& );
    dlerror();
    func_from_file = (func_from_file_t) dlsym(handle, "func_from_file");
    const char *dlsym_error = dlerror();
    if(dlsym_error){
        cerr << "Cannot load symbol 'hello': " << dlsym_error << '\n';
        dlclose(handle);
        return 1;
    }


    int p = taylorTerms; // matrix exponential terminated after p terms
    // appliedError = 1.02 * trueError for both Rerr_bar and L_hatB
    // appliedError begins with _0_ for Rerr_bar

    // int state_dim = dim;
    // int input_dim = dimInput;
    double r = tau;    //sampling period = r

    Eigen::VectorXd L_max(state_dim);
    for(int i=0;i<state_dim;i++)
    	L_max(i) = l_bar;
    Eigen::VectorXd c(state_dim);  // c = centre of the cell
    Eigen::VectorXd ss_eta(state_dim); // state space eta
    Eigen::VectorXd u(input_dim);  //the selected value of input
    for(int i=0;i<state_dim;i++)
    {
        c(i) = x[i];
        ss_eta(i) = 2 * rr[i];
    }
    for(int i=0;i<input_dim;i++)
        u(i) = uu[i];

    Z0 = SymReach::zonotope(c,ss_eta);

    std::vector<SymReach::zonotope> Zn;
    Zn.push_back(Z0);
	// int no_of_steps = finaltime/tau;
	// no_of_steps = (no_of_steps != finaltime/tau) ? ((finaltime/tau)+1) : (finaltime/tau);
	// std::vector<std::vector<SymReach::zonotope>> Zti(no_of_steps);
    std::vector<std::vector<SymReach::zonotope>> Ztp(no_of_steps);
   // Ztp[0] = Zn;
    // std::cout << "r, finaltime, no of steps = " << r << ","
    		// << finaltime << "," << no_of_steps << std::endl;

				// std::cout << 4 << "\n";

			
    for(int i=0;i<no_of_steps;i++)
    {
        int count1 = 0;
        Zn = one_iteration_s(Zn, u, r, p, L_max, count1,
							LinErrorMethod, ss_eta, morder, Zti[i], JacobianB, HessianB, ru, JacobianBu);
			Ztp[i] = Zn;

		std::cout << "step = " << i+1 << ", time = " << (i+1)*r << std::endl;
     }

    std::cout << "Zn.size() = " << Zn.size() << std::endl;
    //for(int i=0;i<Ztp.size();i++)   
    //	plot(Ztp[i],1,2);
//    char ask2 = 'y';
//    while(ask2 == 'y')
//    {
//    	int inda;
//    	std::cout << "\n index of Ztp to plot = ";
//    	std::cin >> inda;
//    	plot(Ztp[inda],1,2);
//    	std::cout << "\n want to plot more?(y/n) = ";
//    	std::cin >> ask2;
//    }


	/* char ask3 = 'y';
	while(ask3 == 'y')
	{
		int ste;
		std::cout << "\nenter the step to plot = ";
		std::cin >> ste;
		plot(Ztp[ste],1,2);
		std::cout << "\nWant to plot more(y/n) = ";
		std::cin >> ask3;
	} */
    dlclose(handle);
	delete[] JacobianB;
	if(dimInput != 0)
	{
		delete[] JacobianBu;
		delete[] JacobianHtemp;
	}
	delete[] HessianB;
    return 0;
}
    
} // end namespace SymReach

#endif  /* REACHABLESET2CPP_H_ */

