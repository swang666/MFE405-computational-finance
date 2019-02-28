// Project7.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

#include "fdm.h"
#include <iostream>
#include <Eigen>
#include <stdlib.h> 
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <string>

using namespace std;
using namespace Eigen;

const static IOFormat CSVFormat(StreamPrecision, DontAlignCols, ", ", "\n");

int main()
{
	ofstream myfile;
	double* prices = new double[13];
	for (int i = 0; i < 13; ++i) {
		prices[i] = 16 - i;
	}
	//third para 0 indicates american, 1 euro
	//fourth para 0 indicates put, 1 call
	cout << "EFDM Euro put" << endl;
	MatrixXd p1a1 = EFDM(10, 1, 1, 0);
	MatrixXd p1a2 = EFDM(10, 3, 1, 0);
	MatrixXd p1a3 = EFDM(10, 4, 1, 0);
	/*myfile.open("p1a1.csv", fstream::out);
	myfile << p1a1.format(CSVFormat);
	myfile.close();
	
	myfile.open("p1a2.csv", fstream::out);
	myfile << p1a2.format(CSVFormat);
	myfile.close();
	
	myfile.open("p1a3.csv", fstream::out);
	myfile << p1a3.format(CSVFormat);
	myfile.close();*/
	cout << "k = 1: " << p1a1((p1a1.rows()-1)/2, 1) << endl;
	cout << "k = 3: " << p1a2((p1a2.rows() - 1) / 2, 1) << endl;
	cout << "k = 4: " << p1a3((p1a3.rows() - 1) / 2, 1) << endl << endl;
	
	cout << "IFDM Euro put" << endl;
	MatrixXd p1b1 = IFDM(10, 1, 1, 0);
	MatrixXd p1b2 = IFDM(10, 3, 1, 0);
	MatrixXd p1b3 = IFDM(10, 4, 1, 0);
	/*myfile.open("p1b1.csv", fstream::out);
	myfile << p1b1.format(CSVFormat);
	myfile.close();
	
	myfile.open("p1b2.csv", fstream::out);
	myfile << p1b2.format(CSVFormat);
	myfile.close();
	
	myfile.open("p1b3.csv", fstream::out);
	myfile << p1b3.format(CSVFormat);
	myfile.close();*/
	cout << "k = 1: " << p1b1((p1b1.rows() - 1) / 2, 1) << endl;
	cout << "k = 3: " << p1b2((p1b2.rows() - 1) / 2, 1) << endl;
	cout << "k = 4: " << p1b3((p1b3.rows() - 1) / 2, 1) << endl << endl;

	cout << "CNFDM Euro put" << endl;
	MatrixXd p1c1 = CNFDM(10, 1);
	MatrixXd p1c2 = CNFDM(10, 3);
	MatrixXd p1c3 = CNFDM(10, 4);
	//myfile.open("p1c1.csv", fstream::out);
	//myfile << p1c1.format(CSVFormat);
	//myfile.close();
	//
	//myfile.open("p1c2.csv", fstream::out);
	//myfile << p1c2.format(CSVFormat);
	//myfile.close();
	//
	//myfile.open("p1c3.csv", fstream::out);
	//myfile << p1c3.format(CSVFormat);
	//myfile.close();
	
	cout << "k = 1: " << p1c1((p1c1.rows() - 1) / 2, 1) << endl;
	cout << "k = 3: " << p1c2((p1c2.rows() - 1) / 2, 1) << endl;
	cout << "k = 4: " << p1c3((p1c3.rows() - 1) / 2, 1) << endl << endl;

	cout << "GFD American put" << endl;
	MatrixXd p2a1 = GFD(10, 0.25, 1, 0);
	MatrixXd p2b1 = GFD(10, 0.25, 0, 0);
	MatrixXd p2c1 = GFD(10, 0.25, 0.5, 0);
	/*myfile.open("p2a1.csv", fstream::out);
	myfile << p2a1.format(CSVFormat);
	myfile.close();
	myfile.open("p2b1.csv", fstream::out);
	myfile << p2b1.format(CSVFormat);
	myfile.close();
	myfile.open("p2c1.csv", fstream::out);
	myfile << p2c1.format(CSVFormat);
	myfile.close();*/
	cout << "explicit: " << p2a1((p2a1.rows() - 1) / 2, 1) << endl;
	cout << "implicit: " << p2b1((p2b1.rows() - 1) / 2, 1) << endl;
	cout << "crank-nicolson: " << p2c1((p2c1.rows() - 1) / 2, 1) << endl << endl;

	cout << "GFD American call" << endl;
	MatrixXd p2a2 = GFD(10, 0.25, 1, 1);
	MatrixXd p2b2 = GFD(10, 0.25, 0, 1);
	MatrixXd p2c2 = GFD(10, 0.25, 0.5, 1);
	/*myfile.open("p2a2.csv", fstream::out);
	myfile << p2a2.format(CSVFormat);
	myfile.close();
	myfile.open("p2b2.csv", fstream::out);
	myfile << p2b2.format(CSVFormat);
	myfile.close();
	myfile.open("p2c2.csv", fstream::out);
	myfile << p2c2.format(CSVFormat);
	myfile.close();*/
	
	cout << "explicit: " << p2a2((p2a2.rows() - 1) / 2, 1) << endl;
	cout << "implicit: " << p2b2((p2b2.rows() - 1) / 2, 1) << endl;
	cout << "crank-nicolson: " << p2c2((p2c2.rows() - 1) / 2, 1) << endl << endl;

}
