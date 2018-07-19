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

const int state_dim = 2;
const int input_dim = 0;
template<typename Tx, typename Tu>
Tx funcLj_system(Tx x, Tu u, Tx xx){
    xx[0] = x[1];
    xx[1] = x[1] - x[0] - x[0]*x[0]*x[1];
    return xx;
}
symbol x0("x[0]"), x1("x[1]");
std::vector<symbol> xs{x0, x1};	// x_symbol
std::vector<symbol> us;
std::vector<ex> f{xs[1],
					(xs[1] - xs[0] - xs[0]*xs[0]*xs[1])};

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
		tau = 0.01;
		finaltime = 7; 
		morder = 20;	// zonotope order above which order reduction
		l_max = 0.01;
		taylorTerms = 4;
		}

			x[0] = (1.55+1.25)*0.5;
			x[1] = (2.35+2.45)*0.5;
			r[0] = (1.55 - 1.25)*0.5;
			r[1] = (2.45 - 2.35)*0.5;
	}
		int no_of_steps = finaltime/tau;
		no_of_steps = (no_of_steps != finaltime/tau) ? ((finaltime/tau)+1) : (finaltime/tau);
		std::vector<std::vector<mstom::zonotope>> Zti(no_of_steps);
		mstom::zonotope Z0;

		TicToc reachtime;
		reachtime.tic();
		ReachableSet(state_dim, input_dim, tau, r, x, ru, u, no_of_steps, 1, l_max, morder, taylorTerms, Zti, Z0);
		reachtime.toc();

		//<-----------Plot and file write
		int plotorder = 3;
		std::vector<mstom::zonotope> PlotStorage;
		mstom::reduce(Zti[0], plotorder); // reduces zonotope order
		PlotStorage = Zti[0];
		for(unsigned int i=1;i<Zti.size();i++)
		{
				Zti[i] = project(Zti[i], 1, 2);	// projection along dimension 1 and 2
				// now Zti has 2D zonotopes
				mstom::reduce(Zti[i], plotorder);
				PlotStorage.insert(PlotStorage.end(),Zti[i].begin(), Zti[i].end());
		}
		mstom::plotstore(PlotStorage, project(Z0, 1, 2));
		std::cout << "Writing File" << std::endl;
		wfile(Zti);
		std::cout << "Plotting" << std::endl;
		plotfilled(PlotStorage, 1, 2); // plots x1 vs x2 
		//------------<

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

