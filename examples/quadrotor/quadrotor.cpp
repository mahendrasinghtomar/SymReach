/*
 * main.cpp
 *
 *  Created on: May 17, 2018
 *      Author: MahendraSinghTomar
 */

#include "ReachableSet2.h"
#include "ReachableSet2cpp.h"
#include <iostream>
#include <array>
#include <cmath>  

using namespace GiNaC; 

const int state_dim = 12;
const int input_dim = 3;	

struct QuadrotorParam{
	double g = 9.81;
	double R = 0.1;
	double l = 0.5;
	double Mrotor = 0.1;
	double M = 1;
	double m = M+4*Mrotor;
	double Jx = 2.0*M*R*R*std::pow(5,-1) + 2*l*l*Mrotor;
	double Jy = Jx;
	double Jz = 2.0*M*R*R*std::pow(5,-1) + 4*l*l*Mrotor;
};

QuadrotorParam gpar;

template<typename Tx, typename Tu>
Tx funcLj_system(Tx x, Tu u, Tx xx){
	
	auto F = gpar.m*gpar.g - 10*(x[2]-u[0]) + 3*x[5];  
	auto tau_phi = -(x[6]-u[1]) - x[9]; 
	auto tau_theta = -(x[7]-u[2]) - x[10];
	double tau_si = 0; 
	
	xx[0] = cos(x[7])*cos(x[8])*x[3] + (sin(x[6])*sin(x[7])*sin(x[8])-cos(x[6])*sin(x[8]))*x[4] + (cos(x[6])*sin(x[7])*cos(x[8])+sin(x[6])*sin(x[8]))*x[5];
	xx[1] = cos(x[7])*sin(x[8])*x[3] + (sin(x[6])*sin(x[7])*sin(x[8])+cos(x[6])*cos(x[8]))*x[4] + (cos(x[6])*sin(x[7])*sin(x[8])-sin(x[6])*cos(x[8]))*x[5];
	xx[2] = sin(x[7])*x[3] - sin(x[6])*cos(x[7])*x[4] - cos(x[6])*cos(x[7])*x[5];
	xx[3] = x[11]*x[4] - x[10] * x[5] - gpar.g*sin(x[7]);
	xx[4] = x[9]*x[5] - x[11]*x[3] + gpar.g*cos(x[7])*sin(x[6]);
	xx[5] = x[10]*x[3] - x[9]*x[4] + gpar.g*cos(x[7])*cos(x[6]) - F*std::pow(gpar.m,-1);
	xx[6] = x[9]+sin(x[6])*tan(x[7])*x[10] + cos(x[6])*tan(x[7])*x[11];
	xx[7] = cos(x[6])*x[10] - sin(x[6])*x[11];
	xx[8] = sin(x[6])*std::pow(cos(x[7]),-1)*x[10] + cos(x[6])*std::pow(cos(x[7]),-1)*x[11];
	xx[9] = (gpar.Jy-gpar.Jz)*std::pow(gpar.Jx,-1)*x[10]*x[11] + tau_phi*std::pow(gpar.Jx,-1);
	xx[10] = (gpar.Jz-gpar.Jx)*std::pow(gpar.Jy,-1)*x[9]*x[11] + tau_theta*std::pow(gpar.Jy,-1);
	xx[11] = (gpar.Jx-gpar.Jy)*std::pow(gpar.Jz,-1)*x[9]*x[10] + tau_si*std::pow(gpar.Jz,-1);
	return xx;
}

symbol x0("x[0]"), x1("x[1]"), x2("x[2]"), x3("x[3]"), x4("x[4]"), x5("x[5]"), x6("x[6]"), x7("x[7]"), x8("x[8]"), x9("x[9]"), x10("x[10]"), x11("x[11]");
std::vector<symbol> xs{x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11};

symbol u0("u[0]"), u1("u[1]"), u2("u[2]");
std::vector<symbol> us{u0, u1, u2};

	auto F = gpar.m*gpar.g - 10*(xs[2]-us[0]) + 3*xs[5];
	auto tau_phi = -(xs[6]-us[1]) - xs[9];
	auto tau_theta = -(xs[7]-us[2]) - xs[10];
	double tau_si = 0; 

std::vector<ex> f{cos(xs[7])*cos(xs[8])*xs[3] + (sin(xs[6])*sin(xs[7])*sin(xs[8])-cos(xs[6])*sin(xs[8]))*xs[4] + (cos(xs[6])*sin(xs[7])*cos(xs[8])+sin(xs[6])*sin(xs[8]))*xs[5],
	cos(xs[7])*sin(xs[8])*xs[3] + (sin(xs[6])*sin(xs[7])*sin(xs[8])+cos(xs[6])*cos(xs[8]))*xs[4] + (cos(xs[6])*sin(xs[7])*sin(xs[8])-sin(xs[6])*cos(xs[8]))*xs[5],
	sin(xs[7])*xs[3] - sin(xs[6])*cos(xs[7])*xs[4] - cos(xs[6])*cos(xs[7])*xs[5],
	xs[11]*xs[4] - xs[10] * xs[5] - gpar.g*sin(xs[7]),
	xs[9]*xs[5] - xs[11]*xs[3] + gpar.g*cos(xs[7])*sin(xs[6]),
	xs[10]*xs[3] - xs[9]*xs[4] + gpar.g*cos(xs[7])*cos(xs[6]) - F*pow(gpar.m,-1),	
	xs[9]+sin(xs[6])*tan(xs[7])*xs[10] + cos(xs[6])*tan(xs[7])*xs[11],
	cos(xs[6])*xs[10] - sin(xs[6])*xs[11],
	sin(xs[6])*pow(cos(xs[7]),-1)*xs[10] + cos(xs[6])*pow(cos(xs[7]),-1)*xs[11],
	(gpar.Jy-gpar.Jz)*pow(gpar.Jx,-1)*xs[10]*xs[11] + tau_phi*pow(gpar.Jx,-1),
	(gpar.Jz-gpar.Jx)*pow(gpar.Jy,-1)*xs[9]*xs[11] + tau_theta*pow(gpar.Jy,-1),
(gpar.Jx-gpar.Jy)*pow(gpar.Jz,-1)*xs[9]*xs[10] + tau_si*pow(gpar.Jz,-1)};


int main() {
	char ask = 'y';
	double tau, finaltime;
	double l_max;
	int morder, taylorTerms;
	int ask_chk = 0;
	while(ask == 'y')
	{
    double x[state_dim], r[state_dim], u[input_dim], ru[input_dim];

	{
		if(ask_chk==0)
		{
		tau = 0.1; 
		finaltime = 5;
		morder = 5;	// zonotope order above which order reduction
		l_max = 0.05;	
		taylorTerms = 4;
		}
			x[0] = 0;	//inertial(north) position
			x[1] = 0; //inertial(east) position
			x[2] = 0; // altitude
			x[3] = 0; //longitudinal velocity
			x[4] = 0; // lateral velocity
			x[5] = 0; // vertical velocity
			x[6] = 0; // roll angle
			x[7] = 0; // pitch angle
			x[8] = 0; // yaw angle
			x[9] = 0; // roll rate
			x[10] = 0; // pitch rate
			x[11] = 0; // yaw rate
			r[0] = 0.4; // lateral veloc
			r[1] = 0.4;
			r[2] = 0.4;
			r[3] = 0.4;
			r[4] = 0.4;
			r[5] = 0.4;
			r[6] = 0;
			r[7] = 0;
			r[8] = 0;
			r[9] = 0;
			r[10] = 0;
			r[11] = 0;
			
			u[0] = 1;	// desired value of height (input to the controller)
			u[1] = 0;	// of roll
			u[2] = 0;	// of pitch
			ru[0] = 0;
			ru[1] = 0;
			ru[2] = 0;
	}
		int no_of_steps = finaltime/tau;
		no_of_steps = (no_of_steps != finaltime/tau) ? ((finaltime/tau)+1) : (finaltime/tau);
		std::vector<std::vector<mstom::zonotope>> Zti(no_of_steps);
		mstom::zonotope Z0;
		
		TicToc reachtime; 
		reachtime.tic();
		ReachableSet(state_dim, input_dim, tau, r, x, ru, u, no_of_steps, 1, l_max, morder, taylorTerms, Zti, Z0);
		reachtime.toc();

		// Storing projection of Zti along dimension 3. Dimension count starts from 1 (not 0).
		int dimplot = 3;
		std::cout << "Writing File" << std::endl;
		wfile_time(Zti, dimplot, tau); //
		std::cout << "Plotting" << std::endl;
		plotfilled(Zti, dimplot, tau); // Plots the dimension 3 vs time

		std::cout << "Do you want to try another tau and finaltime?(y,n):";
		std::cin >> ask;
		if(ask!='y')
			return 0;
		char ask2 = 'n';
		std::cout << "Repeat with previous parameters?(y,n):";
		std::cin >> ask2;
		if(ask2 != 'y')
		{
			std::cout << "enter tau = ";
			std::cin >> tau;
			std::cout << "enter finaltime =  ";
			std::cin >> finaltime;
			std::cout << "enter morder = ";
			std::cin >> morder;
			std::cout << "enter l_max =  ";
			std::cin >> l_max;
		}

		ask_chk += 1;
	}
    return 0;
}

