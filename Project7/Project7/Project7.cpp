// Project7.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "fdm.h"
#include <iostream>
#include <Eigen>

using namespace std;
using namespace Eigen;

int main()
{
	//third para 0 indicates american, 1 euro
	//fourth para 0 indicates put, 1 call
	/*cout << "EFDM Euro put" << endl;
	cout << "k = 1: " << EFDM(10, 1, 1, 0) << endl;
	cout << "k = 3: " << EFDM(10, 3, 1, 0) << endl;
	cout << "k = 4: " << EFDM(10, 4, 1, 0) << endl;

	cout << "IFDM Euro put" << endl;
	cout << "k = 1: " << IFDM(10, 1, 1, 0) << endl;
	cout << "k = 3: " << IFDM(10, 3, 1, 0) << endl;
	cout << "k = 4: " << IFDM(10, 4, 1, 0) << endl;*/

	/*cout << "CNFDM" << endl;
	cout << "k = 1: " << CNFDM(10, 1) << endl;
	cout << "k = 3: " << CNFDM(10, 3) << endl;
	cout << "k = 4: " << CNFDM(10, 4) << endl;*/

	cout << "EFDM American put" << endl;
	cout << "ds = 0.25: " << EFDM_S(10, 0.25, 1, 0) << endl;
	cout << "ds = 1: " << EFDM_S(10, 1, 1, 0) << endl;
	cout << "ds = 1.25: " << EFDM_S(10, 1.25, 1, 0) << endl;
	cout << "EFDM American call" << endl;
	cout << "ds = 0.25: " << EFDM_S(10, 0.25, 1, 1) << endl;
	/*cout << "IFDM American put" << endl;
	cout << "k = 1: " << IFDM(10, 1, 0, 0) << endl;
	cout << "k = 3: " << IFDM(10, 3, 0, 0) << endl;
	cout << "k = 4: " << IFDM(10, 4, 0, 0) << endl;*/
}
