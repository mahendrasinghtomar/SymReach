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

const int state_dim = 7;   
const int input_dim = 0;
template<typename Tx, typename Tu, typename Txx>
Txx funcLj_system(Tx& x, Tu& u, Txx& xx){
    xx[0] = 1.4*x[2] - 0.9*x[0];
    xx[1] = 2.5*x[4] - 1.5*x[1];
	xx[2] = 0.6*x[6] - 0.8*x[1]*x[2];
	xx[3] = 2 - 1.3*x[2]*x[3];
	xx[4] = 0.7*x[0] - x[3]*x[4];
	xx[5] = 0.3*x[0] - 3.1*x[5];
	xx[6] = 1.8*x[5] - 1.5*x[1]*x[6];
    return xx;
}
GiNaC::symbol x0("x[0]"), x1("x[1]"), x2("x[2]"), x3("x[3]"), x4("x[4]"), x5("x[5]"), x6("x[6]");
std::vector<GiNaC::symbol> xs{x0, x1, x2, x3, x4, x5, x6};	// x_symbol

std::vector<GiNaC::symbol> us;


int main() {
	char ask = 'y';
	double tau, finaltime;
	double l_max;
	int morder, taylorTerms;
	int ask_chk = 0;
	while(ask == 'y')
	{
    double x[state_dim], r[state_dim], u[input_dim], ru[input_dim];

	double W;
	{
		if(ask_chk==0)
		{
		tau = 0.1; 
		finaltime = 20;
		morder = 40;	// zonotope order above which order reduction
		l_max = 0.05;	
		taylorTerms = 4;
		W = 0.01;	// initial box radius
		}
			x[0] = 1.2;
			x[1] = 1.05;
			x[2] = 1.5;
			x[3] = 2.4;
			x[4] = 1;
			x[5] = 0.1;
			x[6] = 0.45;
			for(int i=0;i<state_dim;i++)
				r[i] = W;
	}
		int no_of_steps = finaltime/tau;
		no_of_steps = (no_of_steps != finaltime/tau) ? ((finaltime/tau)+1) : (finaltime/tau);
		std::vector<std::vector<SymReach::zonotope>> Zti(no_of_steps);
		SymReach::zonotope Z0;

		TicToc reachtime; 
		reachtime.tic();
        SymReach::ReachableSet(state_dim, input_dim, tau, r, x, ru, u, no_of_steps, 1, l_max, morder, taylorTerms, Zti, Z0);
		reachtime.toc();

		// Storing projection of Zti along dimension 4. Dimension count starts from 1 (not 0).
		std::cout << "Writing File" << std::endl;
		SymReach::wfile_time(Zti, 4, tau); //
		std::cout << "Plotting" << std::endl;
		SymReach::plotfilled(Zti, 4, tau); // Plots the dimension 4 vs time

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
			{
				std::cout << "enter W = ";
				std::cin >> W;
			}
		}

		ask_chk += 1;
	}
    return 0;
}

