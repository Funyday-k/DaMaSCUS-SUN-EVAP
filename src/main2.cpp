#include <chrono>
#include <cmath>
#include <cstring>	 // for strlen
#include <iostream>
#include <mpi.h>

#include "libphysica/Natural_Units.hpp"
#include "libphysica/Special_Functions.hpp"
#include "libphysica/Utilities.hpp"

#include "Data_Generation.hpp"
#include "Parameter_Scan.hpp"
#include "Reflection_Spectrum.hpp"
#include "Solar_Model.hpp"
#include "version.hpp"

#include "Simulation_Trajectory.hpp"
#include <algorithm>

#include "libphysica/Special_Functions.hpp"
#include "libphysica/Statistics.hpp"

#include "obscura/Astronomy.hpp"

using namespace DaMaSCUS_SUN;
using namespace libphysica::natural_units;

double angular_momentum = 0.0;
double time_step = 10.0;
double v_radial =100.0;
double radius = 100000000;
double phi = 0.0;
double G = 6.67e-15;
double delta = 1.0;
int maximum_time_steps = 1000000;
int time_steps = 0;
double dr_dt(double v){return v;}
double dv_dt(double r, double mass){return angular_momentum * angular_momentum / r / r / r - G * mass / r / r;} // G_Newton
double dphi_dt(double r){return angular_momentum / r / r;}

double m_dm = 1.1;
double m_H = 1.0;
double mu_H = m_dm / m_H;
double mu_Hp = (mu_H + 1) / 2;
double mu_Hm = (mu_H - 1) / 2;
double T_s = 1.0e-6;
double n_H = 1.0;
double sigma = 1.0;
double u_H = 1.0;

int main(int argc, char* argv[])
{
	MPI_Init(NULL, NULL);
	int mpi_processes, mpi_rank;
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_processes);
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);

	// Initial terminal output
	auto time_start	  = std::chrono::system_clock::now();
	auto time_start_t = std::chrono::system_clock::to_time_t(time_start);
	auto* ctime_start = ctime(&time_start_t);
	if(ctime_start[std::strlen(ctime_start) - 1] == '\n')
		ctime_start[std::strlen(ctime_start) - 1] = '\0';
	if(mpi_rank == 0)
		std::cout << "[Started on " << ctime_start << "]" << std::endl
				  << PROJECT_NAME << "-" << PROJECT_VERSION << "\tgit:" << GIT_BRANCH << "/" << GIT_COMMIT_HASH << std::endl
				  << DAMASCUS_SUN_LOGO
				  << std::endl
				  << "MPI processes:\t" << mpi_processes << std::endl;

	// Configuration parameters
	Configuration cfg(argv[1], mpi_rank);
	Solar_Model SSM;
	cfg.Print_Summary(mpi_rank);
	MPI_Barrier(MPI_COMM_WORLD);
	////////////////////////////////////////////////////////////////////////
	void Runge_Kutta_45_Step(double mass);


	Solar_Model solar_model;


	/* while(time_steps < maximum_time_steps)
	{
		time_steps++;
		double r_before = radius;
		Runge_Kutta_45_Step(solar_model.Mass(r_before));
		std::cout << "mass" << solar_model.Mass(r_before) << "\n";
		r_before = radius;
	}
*/
	double func1(double w, double v);
	double func2(double w, double v);
	double integral(double(*f)(double,double), double(*g)(double,double), double min, double max);

	integral(func1,func2,0.0, 10.0);
	return 0;
}

void Runge_Kutta_45_Step(double mass)
{
	// RK coefficients:
	double k_r[6];
	double k_v[6];
	double k_p[6];
    double atime = 0.0;

	
	k_r[0] = time_step * dr_dt(v_radial);
	k_v[0] = time_step * dv_dt(radius, mass);
	
	k_p[0] = time_step * dphi_dt(radius);

	k_r[1] = time_step * dr_dt(v_radial + k_v[0] / 4.0);
	k_v[1] = time_step * dv_dt(radius + k_r[0] / 4.0, mass);
	std::cout << "mass" << mass << "\n";
	// k_p[1]=	dt*dphi_dt(radius+k_r[0]/4.0,J);

	k_r[2] = time_step * dr_dt(v_radial + 3.0 / 32.0 * k_v[0] + 9.0 / 32.0 * k_v[1]);
	k_v[2] = time_step * dv_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1], mass);
	k_p[2] = time_step * dphi_dt(radius + 3.0 / 32.0 * k_r[0] + 9.0 / 32.0 * k_r[1]);

	k_r[3] = time_step * dr_dt(v_radial + 1932.0 / 2197.0 * k_v[0] - 7200.0 / 2197.0 * k_v[1] + 7296.0 / 2197.0 * k_v[2]);
	k_v[3] = time_step * dv_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2], mass);
	k_p[3] = time_step * dphi_dt(radius + 1932.0 / 2197.0 * k_r[0] - 7200.0 / 2197.0 * k_r[1] + 7296.0 / 2197.0 * k_r[2]);

	k_r[4] = time_step * dr_dt(v_radial + 439.0 / 216.0 * k_v[0] - 8.0 * k_v[1] + 3680.0 / 513.0 * k_v[2] - 845.0 / 4104.0 * k_v[3]);
	k_v[4] = time_step * dv_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3], mass);
	k_p[4] = time_step * dphi_dt(radius + 439.0 / 216.0 * k_r[0] - 8.0 * k_r[1] + 3680.0 / 513.0 * k_r[2] - 845.0 / 4104.0 * k_r[3]);

	k_r[5] = time_step * dr_dt(v_radial - 8.0 / 27.0 * k_v[0] + 2.0 * k_v[1] - 3544.0 / 2565.0 * k_v[2] + 1859.0 / 4104.0 * k_v[3] - 11.0 / 40.0 * k_v[4]);
	k_v[5] = time_step * dv_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4], mass);
	k_p[5] = time_step * dphi_dt(radius - 8.0 / 27.0 * k_r[0] + 2.0 * k_r[1] - 3544.0 / 2565.0 * k_r[2] + 1859.0 / 4104.0 * k_r[3] - 11.0 / 40.0 * k_r[4]);

	// New values with Runge Kutta 4 and Runge Kutta 5
	double radius_4	  = radius + 25.0 / 216.0 * k_r[0] + 1408.0 / 2565.0 * k_r[2] + 2197.0 / 4101.0 * k_r[3] - 1.0 / 5.0 * k_r[4];
	double v_radial_4 = v_radial + 25.0 / 216.0 * k_v[0] + 1408.0 / 2565.0 * k_v[2] + 2197.0 / 4101.0 * k_v[3] - 1.0 / 5.0 * k_v[4];
	double phi_4	  = phi + 25.0 / 216.0 * k_p[0] + 1408.0 / 2565.0 * k_p[2] + 2197.0 / 4101.0 * k_p[3] - 1.0 / 5.0 * k_p[4];

	double radius_5	  = radius + 16.0 / 135.0 * k_r[0] + 6656.0 / 12825.0 * k_r[2] + 28561.0 / 56430.0 * k_r[3] - 9.0 / 50.0 * k_r[4] + 2.0 / 55.0 * k_r[5];
	double v_radial_5 = v_radial + 16.0 / 135.0 * k_v[0] + 6656.0 / 12825.0 * k_v[2] + 28561.0 / 56430.0 * k_v[3] - 9.0 / 50.0 * k_v[4] + 2.0 / 55.0 * k_v[5];
	double phi_5	  = phi + 16.0 / 135.0 * k_p[0] + 6656.0 / 12825.0 * k_p[2] + 28561.0 / 56430.0 * k_p[3] - 9.0 / 50.0 * k_p[4] + 2.0 / 55.0 * k_p[5];

	std::cout << "k_v," << k_v[0] << k_v[1] << k_v[2] << k_v[3] << k_v[4] << "\n";
	// Error and adapting the time step
	//std::vector<double> errors = {fabs(radius_5 - radius_4), fabs(v_radial_5 - v_radial_4), fabs(phi_5 - phi_4)};
	//std::vector<double> deltas;
	//for(int i = 0; i < 3; i++)
	//	deltas.push_back(0.84 * pow(error_tolerances[i] / errors[i], 1.0 / 4.0));
	//double delta		 = *std::min_element(std::begin(deltas), std::end(deltas));
	double time_step_new = delta * time_step;

	
	atime	  = atime + time_step;
	radius	  = radius_4;
	v_radial  = v_radial_4;
	phi		  = phi_4;
	time_step = time_step_new;

    printf("%f,%f\n",radius,v_radial);
	
}

double func1(double w, double v){
	return pow(mu_Hp,2)/mu_H * v/w * n_H * sigma * (erf((mu_Hp*v+mu_Hm*w)/u_H)-erf((mu_Hp*v-mu_Hm*w)/u_H)+(erf((mu_Hm*v+mu_Hp*w)/u_H)-erf((mu_Hm*v-mu_Hp*w)/u_H)) * pow(M_E,mu_H*(pow(w,2)-pow(v,2))/u_H));
}
 
double func2(double w, double v){
	return 1.0/sqrt(pow(M_PI,2)) * pow(m_dm / (2 * T_s),3/2) * pow(M_E, -m_dm*pow(w,2) / (2*T_s));
}

double func3(double w, double v){
	return w;
}

double func4(double w, double v){
	return v;
}
 
double integral(double(*f)(double,double), double(*g)(double,double), double min, double max){
	double result = 0;
	const int N = 1000;
	double delta = (max - min) / N;
	for(double j = 0; j < max ; j +=delta)
	{
		result = 0;
		for (double i = min+delta; i < max; i+=delta)
		{
			result += f(i,j)*g(i,j)*4*M_PI* pow(i,2)* delta;
			//std::cout << "f(i,j)" << f(i,j) << ", g(i,j)" << g(i,j)<< "\n";
			//std::cout << result << "\n";

		}
		std::cout << j << "," << result << "\n";
		
	}
	
	return 0;
}



