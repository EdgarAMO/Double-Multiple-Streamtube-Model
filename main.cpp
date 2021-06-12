
#include <iostream>
#include <cmath>
#include <string>
#include <iomanip>
#include <sstream>
#include <string>
#include <fstream>
#include <vector>
#include "shape.h"

const int ROWS = 117; 	// number of rows in the coefficients files
const int COLS = 11;	// number of cols in the coefficients files
const int NZ = 21;		// number of steps
const int NT = 35;		// number of tubes
const int RS = 0.01;	// residual

float integrate(float q[NZ][NT], int NZ, int NT, float DZ, float DT);

void parse(std::string cname, float C[ROWS][COLS]); 

float get_coeff(float aa, float re, float C[ROWS][COLS], float AA[ROWS], float RN[COLS]);

int main() {

	float u1[NZ][NT];	// upwind induction factor
	float u2[NZ][NT];	// downwind induction factor
	float w1[NZ][NT];	// upwind dimensionless relative velocity
	float w2[NZ][NT];	// downwind dimensionless relative velocity
	float a1[NZ][NT];	// upwind angle of attack
	float a2[NZ][NT]; 	// downwind angle of attack
	float r1[NZ][NT];	// upwind local Re
	float r2[NZ][NT];	// downwind local Re
	float l1[NZ][NT]; 	// upwind lift coefficient
	float l2[NZ][NT];	// downwind lift coefficient
	float d1[NZ][NT];	// upwind drag coefficient
	float d2[NZ][NT];	// downwind drag coefficient
	float n1[NZ][NT]; 	// upwind normal coefficient
	float n2[NZ][NT];	// downwind normal coefficient
	float t1[NZ][NT];	// upwind tangential coefficient
	float t2[NZ][NT];	// downwind tangential coefficient
	float q1[NZ][NT];   	// upwind local torque
	float q2[NZ][NT];	// downwind local torque

	float CL_TABLE[ROWS][COLS];
	float CD_TABLE[ROWS][COLS];
	float AA[ROWS];
	float RN[COLS];

	/* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + */

	std::cout << "..... Double Multiple-Stream Tube (DoMuST)" << "\n";
	std::cout << "\n";
	std::cout << "............. Choose blade shape: ";
	std::cout << "\n";
	std::cout << "                     Parabolic [0]";
	std::cout << "\n";
	std::cout << "                      Straight [1]";
	std::cout << "\n";
	std::cout << ".................... Type number: ";
	int GEO; std::cin >> GEO;
	std::cout << "\n";
	std::cout << "............... Turbine's radius: ";
	float R; std::cin >> R;
	std::cout << "\n";
	std::cout << ".......... Turbine's half-length: ";
	float H; std::cin >> H;
	std::cout << "\n";
	std::cout << ".......... Turbine's lower point: ";
	float ZL; std::cin >> ZL;
	std::cout << "\n";
	std::cout << "..... Turbine's number of blades: ";
	int N; std::cin >> N;
	std::cout << "\n";
	std::cout << "............... Equatorial chord: ";
	float C1; std::cin >> C1;
	std::cout << "\n";
	std::cout << "...................... Tip chord: ";
	float C2; std::cin >> C2;
	std::cout << "\n";
	std::cout << "............ Choose NACA profile: ";
	std::cout << "\n";
	std::cout << "                    NACA0012 [12]";
	std::cout << "\n";
	std::cout << "                    NACA0015 [15]";
	std::cout << "\n";
	std::cout << "                    NACA0018 [18]";
	std::cout << "\n";
	std::cout << "                    NACA0021 [21]";
	std::cout << "\n";
	std::cout << "                    DU06W200 [20]";
	std::cout << "\n";
	std::cout << ".................... Type number: ";
	int NACA; std::cin >> NACA;
	std::cout << "\n";
	std::cout << "........ Initial tip-speed ratio: ";
	float XI; std::cin >> XI;
	std::cout << "\n";
	std::cout << ".......... Final tip-speed ratio: ";
	float XF; std::cin >> XF;
	std::cout << "\n";
	std::cout << "......... Revolutions-per-minute: ";
	float OM; std::cin >> OM;
	std::cout << "\n";
	std::cout << "............ Wind shear exponent: ";
	float ALPHA; std::cin >> ALPHA;
	std::cout << "\n";

	/* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + */

	// Pi constant:
	float PI = 4.0 * atan(1.0);

	// Global Re:
	float RE = (PI / 30.0) * (R * OM * C1) / 1.5e-5;

	// Turbine's aspect ratio:
	float B = R / H;

	// Equatorial height:
	float ZEQ = ZL + H;

	// Height step size
	float DZ = 2.0 / NZ;					

	// Sector step size
	float DT = PI / NT;						

	// Swept area:
	float S = 0;

	if (GEO == 0) {S = (2.0 / 3.0) * 4 * H * R;}
	if (GEO == 1) {S = 4 * H * R;}

	// Name of the file:
	std::string clname = " ";
	std::string cdname = " ";

	// Parse tables of coefficients
	switch (NACA) {
		case 12:
			clname = "naca0012cl.csv";
			cdname = "naca0012cd.csv";
			break;
		case 15:
			clname = "naca0015cl.csv";
			cdname = "naca0015cd.csv";
			break;
		case 18:
			clname = "naca0018cl.csv";
			cdname = "naca0018cd.csv";
			break;
		case 21:
			clname = "naca0021cl.csv";
			cdname = "naca0021cd.csv";
			break;
		case 20:
			clname = "du06w200cl.csv";
			cdname = "du06w200cd.csv";
			break;
		default:
			std::cout << "\nWrong choice!\n";
	}

	parse(clname, CL_TABLE);
	parse(cdname, CD_TABLE);

	// Fill AA[]
    std::ifstream aaFile;
    aaFile.open("nacaaa.csv");

    int n = 0;
    std::string line;

    while(aaFile.good()){
        std::getline(aaFile, line, '\n');
        AA[n] = atof(line.c_str());
        ++n;
    }

    aaFile.close();

    // Fill RN[]
    std::ifstream reFile;
    reFile.open("nacare.csv");

    n = 0;

    while(reFile.good()){
        std::getline(reFile, line, '\n');
        RN[n] = atof(line.c_str());
        ++n;
    }

    reFile.close();

	/* + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + + */

	std::cout << "* * * * * * * * * * * * * * * * * * * * *\n";
	std::cout << "* * XEQ * * * CP1 * * * CP2 * * * CPT * *\n"; 
	std::cout << "* * * * * * * * * * * * * * * * * * * * *\n";

	// Loop through all tip-speed ratios:
	for (int XEQ = XI; XEQ <= XF; XEQ+=1) {

		// Free stream velocity:
		float vinf = (PI / 30.0) * (R * OM * C1) / XEQ;

		// Loop through all vertical steps:
		for (int i = 0; i < NZ; i++) {

			// Dimensionless height (-0.9, 0.9)
			float zeta = -0.9 + 1.8 * i / (NZ - 1);	

			// Local height
			float zi = (ZL + 0.1 * H) + i * 1.8 * H / (NZ - 1);

			// Height parameter used in the tip-loss factor:
			float symH = H - fabs(zi - ZL - H);

			// Free stream ratio
			float vz = std::pow((zi / ZEQ), ALPHA);

			// Local free stream velocity:
			float vinfi = vinf * vz; 

			// Dimensionless radius
			float eta = get_eta(zeta, GEO);	

			// Blade's angle over the equator
			float del = get_del(zeta, B, GEO);   	

			// Local chord
			float c = 0.0;

			if (zeta > 0.00) { c = zeta * (C2 - C1) + C1; }
			if (zeta <= 0.0) { c = zeta * (C1 - C2) + C1; }

			// Loop through upwind tubes:
			for (int j = 0; j < NT; j++) {

				// Theta in radians
				float t;

				t = (-PI / 2.0 + DT / 2.0) + j * (PI - DT) / (NT - 1);

				// Initialize induction factor:
				float a = 0;

				int iter = 0;

				while (1) {

					u1[i][j] = 1.0 - a;

					// Tip losses:
					float v = (2 * u1[i][j] - 1) * vinfi;
					float s = (PI * v) / (N * OM * PI / 30.0);
					float arg1 = std::exp(-PI * symH / s);
					float term1 = std::acos(arg1);
					float arg2 = std::exp(-PI * H / s);
					float term2 = std::acos(arg2);
					float f = term1 / term2; 

					// Blade's aspect ratio:
					float AR = get_ar(B, H, C1, GEO);
					float a0 = 1.8 * PI * (1 + 0.8 * NACA / 100.0);

					// Local tip-speed ratio:
					float x = (XEQ * eta) / (u1[i][j] * vz);

					// Dimensionless relative local velocity:
					float ls = std::pow(x - sin(t), 2);
					float rs = std::pow(cos(t), 2) * std::pow(cos(del), 2);
					w1[i][j] = std::sqrt(ls + f * f * rs);

					// Angle of attack:
					a1[i][j] = asin(f * cos(t) * cos(del) / w1[i][j]);

					// Local Re:
					r1[i][j] = (RE * eta / x) * w1[i][j];

					// Get lift and drag:
					l1[i][j] = get_coeff(a1[i][j], r1[i][j], CL_TABLE, AA, RN);
					d1[i][j] = get_coeff(a1[i][j], r1[i][j], CD_TABLE, AA, RN);

					// Finite aspect coefficients:
					float cl = l1[i][j] / (1 + a0 / (AR * PI));
					float cd = d1[i][j] + std::pow(cl, 2) / (AR * PI);

					// Effective angle of attack:
					a1[i][j] -= cl / (PI * AR);

					// Normal and tangential:
					n1[i][j] = cl * cos(a1[i][j]) + cd* sin(a1[i][j]);
					t1[i][j] = cl * sin(a1[i][j]) - cd* cos(a1[i][j]);

					// New induction factors:
					float anew = 0;

					// Auxiliary coefficient:
					float coeff = (N * c) / (8.0 * PI * R * eta);

					// Upwind function:
					ls = (n1[i][j] * cos(t)) / fabs(cos(t));
					rs = (t1[i][j] * sin(t)) / (fabs(cos(t)) * cos(del));

					float fu;

					fu = (ls + rs) * w1[i][j] * w1[i][j];

					if (a <= 0.3) {

						float term = std::pow(a, 2);
						anew = coeff * fu * std::pow(1.0 - a, 2) + term; 
					}

					else { 

						float term = std::pow(a, 2) * (5.0 - 3.0 * a) / 4.0;
						anew = coeff * fu * std::pow(1.0 - a, 2) + term; 
					}

					iter++;

					if (iter > 50) { 

						u1[i][i] = -1.0 * u1[i][j];
						anew = 1.0 - u1[i][j];
						break; 

					}

					if (anew > 1.0) {

						u1[i][j] = 0.0;
						anew = 1.0 - u1[i][j];
						break;

					}

					if (abs(anew - a) / a < RS) { break; }

					a = anew;
					
				}

				// Integrand's coefficient:
				float coeff = (N * c * H) / (2.0 * PI * S);
			
				// W / Vinf:
				float winf = std::pow(w1[i][j] * u1[i][j] * vz, 2);

				// eta / cos(del):
				float quot = eta / cos(del);

				// local torque:
				q1[i][j] = coeff * t1[i][j] * winf * quot;

			} // Upwind loop ends...



			/* * * * * * * * * * Turbine's half delimiter * * * * * * * * * */


			// Loop through downwind tubes:
			for (int j = 0; j < NT; j++) {

				// Theta in radians
				float t;

				t = (PI / 2.0 + DT / 2.0) + j * (PI - DT) / (NT - 1);

				// Initialize induction factor according to upstream tube:
				float a = 1.0 - u1[i][NT - 1 - j];

				int iter = 0;

				while (1) {

					u2[i][j] = 1.0 - a;

					// Tip losses:
					float uu = u1[i][NT - 1 - j];
					float ud = u2[i][j];
					float v = (2 * uu - 1) * (2 * ud - 1) * vinfi;

					if (v < 0.) { v = 1.0; }

					float s = (PI * v) / (N * OM * PI / 30.0);
					float arg1 = std::exp(-PI * symH / s);
					float term1 = std::acos(arg1);
					float arg2 = std::exp(-PI * H / s);
					float term2 = std::acos(arg2);
					float f = term1 / term2;

					// Blade's aspect ratio:
					float AR = get_ar(B, H, C1, GEO);
					float a0 = 1.8 * PI * (1 + 0.8 * NACA / 100.0);

					// Local tip-speed ratio:
					float den = u2[i][j] * (2 * u1[i][NT -1 - j] - 1.0);
					float x = (XEQ * eta) / (den * vz);

					// Dimensionless relative local velocity:
					float ls = std::pow(x - sin(t), 2);
					float rs = std::pow(cos(t), 2) * std::pow(cos(del), 2);
					w2[i][j] = std::sqrt(ls + f * f * rs);	

					// Angle of attack:
					a2[i][j] = asin(f * cos(t) * cos(del) / w2[i][j]);

					// Local Re:
					r2[i][j] = (RE * eta / x) * w2[i][j];

					// Get lift and drag:
					l2[i][j] = get_coeff(a2[i][j], r2[i][j], CL_TABLE, AA, RN);
					d2[i][j] = get_coeff(a2[i][j], r2[i][j], CD_TABLE, AA, RN);

					// Finite aspect coefficients:
					float cl = l2[i][j] / (1 + a0 / (AR * PI));
					float cd = d2[i][j] + std::pow(cl, 2) / (AR * PI);

					// Effective angle of attack:
					a2[i][j] -= cl / (PI * AR);

					// Normal and tangential
					n2[i][j] = cl * cos(a2[i][j]) + cd* sin(a2[i][j]);
					t2[i][j] = cl * sin(a2[i][j]) - cd* cos(a2[i][j]);

					// New induction factors:
					float anew = 0;

					// Auxiliary coefficient:
					float coeff = (N * c) / (8.0 * PI * R * eta);

					// Downwind function:
					ls = (n2[i][j] * cos(t)) / fabs(cos(t));
					rs = (t2[i][j] * sin(t)) / (fabs(cos(t)) * cos(del));

					float fd;

					fd = (ls + rs) * w2[i][j] * w2[i][j];

					if (a <= 0.3) {

						float term = std::pow(a, 2);
						anew = coeff * fd * std::pow(1.0 - a, 2) + term; 
					}

					else { 

						float term = std::pow(a, 2) * (5.0 - 3.0 * a) / 4.0;
						anew = coeff * fd * std::pow(1.0 - a, 2) + term; 
					}

					iter++;

					if (iter > 50) { 

						u2[i][i] = -1.0 * u2[i][j];
						anew = 1.0 - u2[i][j];
						break; 

					}

					if (anew > 1.0) {

						u2[i][j] = 0.0;
						anew = 1.0 - u2[i][j];
						break;

					}

					if (abs(anew - a) / a < RS) { break; }

					a = anew;
					
				}

				// save u2 as V'/Vinf instead of V'/Ve:
				u2[i][j] = u2[i][j] * (2 * u1[i][NT -1 - j] - 1.0);

				// Integrand's coefficient:
				float coeff = (N * c * H) / (2.0 * PI * S);
			
				// W' / Vinf:
				float winf = std::pow(w2[i][j] * u2[i][j] * vz, 2);

				// eta / cos(del):
				float quot = eta / cos(del);

				// local torque:
				q2[i][j] = coeff * t2[i][j] * winf * quot;

			} // Downwind loop ends...

		} // Step loop ends...

		float cq1 = integrate(q1, NZ, NT, DZ, DT);
		float cq2 = integrate(q2, NZ, NT, DZ, DT);
		float cp1 = XEQ * cq1;
		float cp2 = XEQ * cq2;
		float cpt = cp1 + cp2;

		std::cout << "  " << std::setw(7) << std::setprecision(3) << XEQ;
		std::cout << "   " << std::setw(7) << std::setprecision(3) << cp1;
		std::cout << "   " << std::setw(7) << std::setprecision(3) << cp2;
		std::cout << "   " << std::setw(7) << std::setprecision(3) << cpt;
		std::cout << "\n";

	} // Tip-speed ratio loop ends...

	return 0;

}


float integrate(float q[NZ][NT], int NZ, int NT, float DZ, float DT) {

	float sums[NT];

	for (int j = 0; j < NT; j++) {

		float count = 0.0;

		for (int i = 1; i <= (NZ - 1) / 2; i++) {

			int I = (2 * i) - 2;
			int J = (2 * i) - 1;
			int K = 2* i;

			count += q[I][j] + 4 * q[J][j] + q[K][j];

		}

		sums[j] = (DZ / 3.0) * count;

	}

	float count = 0.0;

	for (int j = 1; j <= (NT - 1) / 2; j++) {

		int I = (2 * j) - 2;
		int J = (2 * j) - 1;
		int K = 2* j;

		count += sums[I] + 4 * sums[J] + sums[K];
	}

	count *= (DT / 3.0); 

	return count;
}

void parse(std::string cname, float C[ROWS][COLS]) {
    
    // cname: ...... coefficients file name
    // C[][]: ...... coefficients arrays

    std::ifstream cFile;
    cFile.open(cname.c_str());

    std::string line, field, element;

    std::vector< std::vector<std::string> > array;
    std::vector< std::string > v;

    while ( std::getline(cFile, line) )

    {
        v.clear();
        std::stringstream ss(line);

        while ( std::getline(ss, field, ',') )
        {
            v.push_back(field);
        }

        array.push_back(v);
    }

    for (std::size_t i = 0; i < array.size(); ++i)
    {
        for (std::size_t j = 0; j < array[i].size(); ++j)
        {
            element = array[i][j];
            C[i][j] = atof(element.c_str());
        }
    }

    cFile.close();

}

float get_coeff(float aa, float re, float C[ROWS][COLS], 
									float AA[ROWS], 
									float RN[COLS]) {

    float TODEG = 180 / (4.0 * atan(1.0));

    // Convert radians to degrees
    aa *= TODEG;

    // Normalize Re
    re /= 1e6;

    int i0 = 0;
    int i1 = 0;
    int j0 = 0;
    int j1 = 0;

    // AA and RE locii
    float aLoc = 0.;
    float rLoc = 0.;

    // coefficient at aa
    float cRej0 = 0.;
    float cRej1 = 0.;

    // Find i0 and i1
    for (int k = 0; k < ROWS; k++){

        if (aa >= AA[k] && aa<= AA[k + 1]) {
            i0 = k;
            i1 = k + 1;
            aLoc = (aa - AA[k]) / (AA[k + 1] - AA[k]);
            break;
        }
    }

    // Find j0 and j1
    for (int k = 0; k < COLS; k++){

        if (re >= RN[k] && re<= RN[k + 1]) {
            j0 = k;
            j1 = k + 1;
            rLoc = (re - RN[k]) / (RN[k + 1] - RN[k]);
            break;
        }
    }

    // Coefficients at aa
    cRej0 = C[i0][j0] + aLoc * (C[i1][j0] - C[i0][j0]);
    cRej1 = C[i0][j1] + aLoc * (C[i1][j1] - C[i0][j1]);

    // Coefficient at re
    return cRej0 + rLoc * (cRej1 - cRej0);

}

	



	
