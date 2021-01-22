/*********************************************************************************

			   Main View for Point Cloud Remeshing and Mesh 
			                Reconstruction

						Updating in 2020/09/09

						   By Dr. Chenlei Lv

			The functions includes:
			1. The pipeline of point cloud-based shape reconstruction;
			


*********************************************************************************/
#pragma once
#include "View.h"
#include <stdio.h>
#include <stdlib.h>
#include "trackball.h"
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include "Method_AIVS_RNPro.hpp"
using namespace std;

int main()
{
	cout << "start!" << endl;

	//input parameters
	//support form: obj, xyz, ply, txt, off
	string pathObj = "woodMan";	
	string pathFile = "data//shrec//" + pathObj + ".obj"; 
	int pointNum = 10000;//simplification Number
	int classificationIndex = 1;//1: isotropic; 2: Curvature sensitive; 3: Edge. 
	double multi = 10;//default multi = 2 ~ 10
	int optIter = 5;//default optIter = 5
	
	//init AIVS
	AIVS_RNPro ar;
	ar.AIVS_RNPro_init(pathFile);		
	clock_t t1;
	clock_t t2;
	t1 = clock();	
	ar.AIVS_RNPro_Remesh(pathObj, pointNum, classificationIndex, multi, optIter);
	t2 = clock();
	std::cout << "Running time:" << (t2 - t1) / 1000.0 << "s" << endl;	
  	cout << "Finish!" << endl;

}


















