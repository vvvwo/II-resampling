/*********************************************************************************
             
			 Mesh Optimization to achieve isotropic remeshing

						Updating in 2020/08/28

						   By Dr. Chenlei Lv

			The functions includes:
			Cropping mesh
			[Botsch 2004] A Remeshing Approach to Multiresolution Modeling.

*********************************************************************************/
#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <math.h> 
#include <algorithm>
#include <pcl/io/pcd_io.h>  //File input$output
#include <pcl/octree/octree_search.h>  //octree define
#include <pcl/point_types.h>  //point type
#include <pcl/kdtree/kdtree_flann.h>
using namespace std;

class MeshOptimization {

private:

	int originalNumber;
	string fileName;
	vector<vector<double>> points;//point cloud
	vector<vector<int>> faceInfor;//trangulars
	vector<vector<int>> pointNeighbor;//neighbor faces
	double L_ave; //the ave edge length of all borders
	int optIter;

public:

	void MeshOptimization_init(string fileN, vector<vector<double>> pointsi, vector<vector<int>> faceInfori) {
		originalNumber = pointsi.size();
		fileName = fileN;
		cout << "Mesh optimization start!" << endl;
		points = pointsi;
		faceInfor = faceInfori;//the index of first point is 1		
		//init neighbor structure, correct point index to first point index is 0
		pointNeighbor.resize(points.size());
		for (int i = 0; i < faceInfor.size(); i++) {
			//faceInfor[i][0]--;
			//faceInfor[i][1]--;
			//faceInfor[i][2]--;
			int b1 = faceInfor[i][0];
			int b2 = faceInfor[i][1];
			int b3 = faceInfor[i][2];			
			pointNeighbor[b1].push_back(b2);
			pointNeighbor[b1].push_back(b3);
			pointNeighbor[b2].push_back(b3);
			pointNeighbor[b2].push_back(b1);
			pointNeighbor[b3].push_back(b1);
			pointNeighbor[b3].push_back(b2);		
		}
		//init ave length of trangular border
		MeshOptimization_AveEdge_init();
		cout << "L_ave:" << L_ave << endl;
	}
	
	void MeshOptimization_Start(int iter, bool judge) {

		optIter = iter;

		cout << "Mesh Split_Collapse start!" << endl;
		while (iter) {
			cout << "iter:" << iter << endl;
			//processing split, collapse, Flip and Tangent smoothing
			cout << "Split start:";
			MeshOptimization_Split();			
			
			cout << "Collapse start:";
			MeshOptimization_Collapse(judge);
			
			cout << "Flip start:";
			MeshOptimization_Flip();			
			
			cout << "Tangent smoothing start:";
			MeshOptimization_TangentSmoothing();
			iter--;		
		}
		//cout << "tangent moving:" << endl;
		//MeshOptimization_TangentSmoothing();
		MeshOptimization_UpdateStructure_Face();		
		//cout << "Mesh Save." << endl;
		//MeshOptimization_SaveOBJ_Face();
		cout << endl;

	}

	vector<vector<double>> MeshOptimization_GetPoints() {
		return points;
	}

	vector<vector<int>> MeshOptimization_GetFaces() {
		return faceInfor;
	}

private:	

	int MeshOptimization_Split_Collapse_Old() {//Split-Collapse edges based on ave

		vector<bool> pointActive(points.size(), true);//state of point, if ture, can be Split-Collapse, alse can not		
		vector<int> point_Delete;//store the index of delete points.

		//store all edges satisfy the condition of Split-Collapse
		vector<vector<int>> Split_list;
		vector<vector<int>> Collapse_List;

		//cout << "Achieve all edges for Split-Collapse" << endl;
		//achieve all edges for Split-Collapse
		for (int i = 0; i < pointNeighbor.size(); i++) {
			vector<int> b2_list;
			int b1 = i;
			if (!pointActive[b1]) {//if b1 has been labeled by other points for Split-Collapse, then dont process it
				continue;
			}
			for (int j = 0; j < pointNeighbor[i].size(); j++) {
				int b2 = pointNeighbor[i][j];
				if (!pointActive[b2]) {//if b2 has been labeled by other points for Split-Collapse, then dont process it
					continue;
				}
				if (!MeshOptimization_Exist(b2, b2_list)) {
					b2_list.push_back(b2);
					double length12 = MeshOptimization_EdgeLength(b1, b2);
					if (length12 >= 1.33 * L_ave) {
						vector<int> edgeij;
						edgeij.push_back(b1);
						edgeij.push_back(b2);
						Split_list.push_back(edgeij);

						//label the points, which are waitiing for Split-Collapse, to false, 
						pointActive[b1] = false;
						pointActive[b2] = false;
						vector<int> pointNeighbor_b1 = pointNeighbor[b1];
						vector<int> pointNeighbor_b2 = pointNeighbor[b2];
						for (int k = 0; k < pointNeighbor_b1.size(); k++) {
							int b3 = pointNeighbor_b1[k];
							pointActive[b3] = false;
						}

						for (int k = 0; k < pointNeighbor_b2.size(); k++) {
							int b3 = pointNeighbor_b2[k];
							pointActive[b3] = false;
						}
					}
					if (length12 < 0.8 * L_ave) {
						vector<int> edgeij;
						edgeij.push_back(b1);
						edgeij.push_back(b2);
						Collapse_List.push_back(edgeij);

						//label the points, which are waitiing for Split-Collapse, to false, 
						pointActive[b1] = false;
						pointActive[b2] = false;
						vector<int> pointNeighbor_b1 = pointNeighbor[b1];
						vector<int> pointNeighbor_b2 = pointNeighbor[b2];
						for (int k = 0; k < pointNeighbor_b1.size(); k++) {
							int b3 = pointNeighbor_b1[k];
							pointActive[b3] = false;
						}

						for (int k = 0; k < pointNeighbor_b2.size(); k++) {
							int b3 = pointNeighbor_b2[k];
							pointActive[b3] = false;
						}

					}
				}
			}
		}

		int numberProcess = Split_list.size() < Collapse_List.size() ?
			Split_list.size() : Collapse_List.size();//store the min value of two lists Split_list and Collapse_list

		cout << "Pair number:"<< numberProcess << endl;		

		//cout << "Pair processing for Split-Collapse" << endl;
		//Pair processing for Split-Collapse
		int pointsEndIndex = points.size();//store last index of points list
		for (int i = 0; i < numberProcess; i++) {

			//cout << "Pair index:" << i << ",";
			//************************************Split****************************************
			int bs1 = Split_list[i][0];
			int bs2 = Split_list[i][1];
			vector<double> middlePoint(3);
			middlePoint[0] = (points[bs1][0] + points[bs2][0]) / 2;
			middlePoint[1] = (points[bs1][1] + points[bs2][1]) / 2;
			middlePoint[2] = (points[bs1][2] + points[bs2][2]) / 2;
			vector<int> pointNeighbor_b1 = pointNeighbor[bs1];

			vector<int> pointNeighbor_middle;

			//1. Split: add middle point into points list
			points.push_back(middlePoint);

			//2. Split: update trangular structure
			//vector<vector<int>> bs3_list_T;
			vector<int> bs3_list;
			bs3_list.push_back(bs1);
			bs3_list.push_back(bs2);
			//2.1 store related points
			for (int i = 0; i < pointNeighbor_b1.size()/2; i++) {
				int b1i1 = pointNeighbor_b1[2 * i];
				int b1i2 = pointNeighbor_b1[2 * i + 1];
				if (b1i1 == bs2) {
					//vector<int> bs3_list_i;
					//bs3_list_i.push_back(bs1);
					//bs3_list_i.push_back(b1i1);
					//bs3_list_i.push_back(b1i2);
					//bs3_list_T.push_back(bs3_list_i);
					bs3_list.push_back(b1i2);									
				}
				else if(b1i2 == bs2){
					//vector<int> bs3_list_i;
					//bs3_list_i.push_back(bs1);
					//bs3_list_i.push_back(b1i1);
					//bs3_list_i.push_back(b1i2);
					//bs3_list_T.push_back(bs3_list_i);
					bs3_list.push_back(b1i1);				
				}
				else {
					continue;				
				}
			}

			//2.2 update structure
			//update 4 related points
			vector<vector<int>> T_n;
 			for (int i = 0; i < bs3_list.size(); i++) {
				int b_index = bs3_list[i];
				if (b_index == bs1 || b_index == bs2) {
					for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
						int b_index_j1 = pointNeighbor[b_index][2 * j];
						int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
						if (b_index_j1 == bs1 || b_index_j1 == bs2) {
							pointNeighbor[b_index][2 * j] = pointsEndIndex;
							pointNeighbor_middle.push_back(b_index_j2);
							pointNeighbor_middle.push_back(b_index);

						}
						if (b_index_j2 == bs1 || b_index_j2 == bs2) {
							pointNeighbor[b_index][2 * j + 1] = pointsEndIndex;
							pointNeighbor_middle.push_back(b_index);
							pointNeighbor_middle.push_back(b_index_j1);
						}
					}					
				}
				else {
					for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
						int b_index_j1 = pointNeighbor[b_index][2 * j];
						int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
						if ((b_index_j1 == bs1 || b_index_j1 == bs2)&&(b_index_j2 == bs1 || b_index_j2 == bs2)) {
							pointNeighbor[b_index][2 * j] = pointsEndIndex;
							pointNeighbor[b_index].push_back(b_index_j1);
							pointNeighbor[b_index].push_back(pointsEndIndex);							
							pointNeighbor_middle.push_back(b_index_j2);
							pointNeighbor_middle.push_back(b_index);
							pointNeighbor_middle.push_back(b_index);
							pointNeighbor_middle.push_back(b_index_j1);
						}						
					}				
				}						
			}

			//2.3 update neighbor structure
			pointNeighbor.push_back(pointNeighbor_middle);
			//Update last point index
			pointsEndIndex++;//update last index of points list

			//************************************Collapse*************************************
			int bc1 = Collapse_List[i][0];
			int bc2 = Collapse_List[i][1];
			vector<double> middlePointc(3);
			middlePointc[0] = (points[bc1][0] + points[bc2][0]) / 2;
			middlePointc[1] = (points[bc1][1] + points[bc2][1]) / 2;
			middlePointc[2] = (points[bc1][2] + points[bc2][2]) / 2;
			vector<int> pointNeighbor_bc1 = pointNeighbor[bc1];
			point_Delete.push_back(bc1);
			point_Delete.push_back(bc2);

			vector<int> pointNeighbor_middlec;

			//1. Collapse: add middle point into points list
			points.push_back(middlePointc);

			vector<int> pc_list;
			
			//2. Collapse: update trangular structure
			//2.1 add the related points in the list
			for (int i = 0; i < pointNeighbor[bc1].size(); i++) {
				int p_i = pointNeighbor[bc1][i];
				if (p_i == bc1 || p_i == bc2) {
					continue;				
				}
				if (!MeshOptimization_Exist(p_i, pc_list)) {
					pc_list.push_back(p_i);
				}
			}

			for (int i = 0; i < pointNeighbor[bc2].size(); i++) {
				int p_i = pointNeighbor[bc2][i];
				if (p_i == bc1 || p_i == bc2) {
					continue;
				}
				if (!MeshOptimization_Exist(p_i, pc_list)) {
					pc_list.push_back(p_i);
				}
			}

			for (int i = 0; i < pc_list.size();i++) {
				int p_index = pc_list[i];
				bool rejudge = false;
				for (int j = 0; j < pointNeighbor[p_index].size() / 2; j++) {
					int pc1 = pointNeighbor[p_index][2 * j];
					int pc2 = pointNeighbor[p_index][2 * j + 1];
					if ((pc1 == bc1 || pc1 == bc2) && (pc2 == bc1 || pc2 == bc2)) {
						pointNeighbor[p_index][2 * j] = -1;
						pointNeighbor[p_index][2 * j + 1] = -1;
						rejudge = true;
					}
					else if (pc1 == bc1 || pc1 == bc2) {
						pointNeighbor[p_index][2 * j] = pointsEndIndex;
						pointNeighbor_middlec.push_back(pc2);
						pointNeighbor_middlec.push_back(p_index);
					}
					else if (pc2 == bc1 || pc2 == bc2) {
						pointNeighbor[p_index][2 * j + 1] = pointsEndIndex;
						pointNeighbor_middlec.push_back(p_index);
						pointNeighbor_middlec.push_back(pc1);					
					}
					else {
						continue;					
					}									
				}
				if (rejudge) {
					vector<int> pointNeighbor_new;
					for (int j = 0; j < pointNeighbor[p_index].size(); j++) {
						if (pointNeighbor[p_index][j] != -1) {
							pointNeighbor_new.push_back(pointNeighbor[p_index][j]);						
						}
					}
					pointNeighbor[p_index].clear();
					pointNeighbor[p_index] = pointNeighbor_new;
				}
			}

			pointNeighbor.push_back(pointNeighbor_middlec);
			//Update 
			pointsEndIndex++;//update last index of points list
		}
		//cout << endl;
		//cout << "Update structure." << endl;	

		MeshOptimization_UpdateStructure(point_Delete);
		return numberProcess;
	}

	//Split and Collapse
	int MeshOptimization_Split_Collapse() {//Split-Collapse edges based on ave

		vector<bool> pointActive(points.size(), true);//state of point, if ture, can be Split-Collapse, alse can not		
		vector<int> point_Delete;//store the index of delete points.

		//store all edges satisfy the condition of Split-Collapse
		vector<vector<int>> Split_list;
		vector<vector<int>> Collapse_List;

		//cout << "Achieve all edges for Split-Collapse" << endl;
		//achieve all edges for Split-Collapse
		for (int i = 0; i < pointNeighbor.size(); i++) {
			vector<int> b2_list;
			int b1 = i;			
			for (int j = 0; j < pointNeighbor[i].size(); j++) {
				int b2 = pointNeighbor[i][j];
				if (!MeshOptimization_Exist(b2, b2_list)&& b1 < b2) {
					b2_list.push_back(b2);
					double length12 = MeshOptimization_EdgeLength(b1, b2);
					if (length12 >= 1.33 * L_ave) {
						vector<int> edgeij;
						edgeij.push_back(b1);
						edgeij.push_back(b2);
						Split_list.push_back(edgeij);										
					}
					if (length12 < 0.8 * L_ave) {
						//Test if there have invaild condition
						int valence = MeshOptimization_Edge_Valence(b1, b2);						
						if (valence <= 2) {
							vector<int> edgeij;
							edgeij.push_back(b1);
							edgeij.push_back(b2);
							Collapse_List.push_back(edgeij);						
						}						
					}
				}
			}
		}

		//store the processing index of two lists
		int Split_list_index = 0;
		int Collapse_List_index = 0;
		int pointsEndIndex = points.size();//the index for last point index, by update point list
		int numberProcess = 0;
		while (1) {					
			//pay attention the pointActive list update
			//Split test:			
			if (Split_list_index >= Split_list.size()) {
				break;
			}
			int bs1 = -1;
			int bs2 = -1;
			for (int i = Split_list_index; i < Split_list.size(); i++) {
				int bs1_i = Split_list[i][0];
				int bs2_i = Split_list[i][1];
				if (pointActive[bs1_i] && pointActive[bs2_i]) {
					bs1 = bs1_i;
					bs2 = bs2_i;
					Split_list_index = i + 1;					
					break;
				}			
			}
			if (bs1 == -1) {
				break;
			}

			//Collapse test
			if (Collapse_List_index >= Collapse_List.size()) {
				break;
			}
			int bc1 = -1;
			int bc2 = -1;
			for (int i = Collapse_List_index; i < Collapse_List.size(); i++) {
				int bc1_i = Collapse_List[i][0];
				int bc2_i = Collapse_List[i][1];
				int valenceE = MeshOptimization_Edge_Valence(bc1_i, bc2_i);
				if (bc1_i == bs1 || bc1_i == bs2 || bc2_i == bs1 || bc2_i == bs2) {
					continue;				
				}
				if (pointActive[bc1_i] && pointActive[bc2_i]&& valenceE <= 2) {
					bc1 = bc1_i;
					bc2 = bc2_i;
					Collapse_List_index = i + 1;
					break;
				}
			}
			if (bc1 == -1) {
				break;
			}

			//Split processing
			//Split: 1 achieve middle point
			vector<double> middlePoint(3);
			middlePoint[0] = (points[bs1][0] + points[bs2][0]) / 2;
			middlePoint[1] = (points[bs1][1] + points[bs2][1]) / 2;
			middlePoint[2] = (points[bs1][2] + points[bs2][2]) / 2;
			points.push_back(middlePoint);
			pointActive.push_back(true);

			vector<int> pointNeighbor_b1 = pointNeighbor[bs1];
			vector<int> pointNeighbor_middle;//store the structure of middle point

			//Split: 2 store related points
			vector<int> bs3_list;
			bs3_list.push_back(bs1);
			bs3_list.push_back(bs2);			
			for (int i = 0; i < pointNeighbor_b1.size() / 2; i++) {
				int b1i1 = pointNeighbor_b1[2 * i];
				int b1i2 = pointNeighbor_b1[2 * i + 1];
				if (b1i1 == bs2) {					
					bs3_list.push_back(b1i2);
				}
				else if (b1i2 == bs2) {					
					bs3_list.push_back(b1i1);
				}
				else {
					continue;
				}
			}
			
			//Split: 3 update the structure			
			vector<vector<int>> T_n;
			for (int i = 0; i < bs3_list.size(); i++) {
				int b_index = bs3_list[i];				
				if (b_index == bs1 || b_index == bs2) {
					for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
						int b_index_j1 = pointNeighbor[b_index][2 * j];
						int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
						if (b_index_j1 == bs1 || b_index_j1 == bs2) {
							pointNeighbor[b_index][2 * j] = pointsEndIndex;
							pointNeighbor_middle.push_back(b_index_j2);
							pointNeighbor_middle.push_back(b_index);
							//if (!MeshOptimization_Check_Structure(pointNeighbor[b_index])) {
								//cout << "Repeat!" << endl;							
							//};
							//if (!MeshOptimization_Check_Structure(pointNeighbor[b_index])) {
								//cout << "Repeat!" << endl;
							//};
						}
						if (b_index_j2 == bs1 || b_index_j2 == bs2) {
							pointNeighbor[b_index][2 * j + 1] = pointsEndIndex;
							pointNeighbor_middle.push_back(b_index);
							pointNeighbor_middle.push_back(b_index_j1);							
							//if (!MeshOptimization_Check_Structure(pointNeighbor[b_index])) {
								//cout << "Repeat!" << endl;
							//};
							//if (!MeshOptimization_Check_Structure(pointNeighbor[b_index])) {
								//cout << "Repeat!" << endl;
							//};
						}
					}
				}
				else {
					for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
						int b_index_j1 = pointNeighbor[b_index][2 * j];
						int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
						if ((b_index_j1 == bs1 || b_index_j1 == bs2) && (b_index_j2 == bs1 || b_index_j2 == bs2)) {
							pointNeighbor[b_index][2 * j] = pointsEndIndex;
							pointNeighbor[b_index].push_back(b_index_j1);
							pointNeighbor[b_index].push_back(pointsEndIndex);
							pointNeighbor_middle.push_back(b_index_j2);
							pointNeighbor_middle.push_back(b_index);
							pointNeighbor_middle.push_back(b_index);
							pointNeighbor_middle.push_back(b_index_j1);
						}
					}
				}
			}

			//Split: 4 update neighbor structure
			pointNeighbor.push_back(pointNeighbor_middle);
			pointsEndIndex++;

            //Collapse processing;
			//Collapse: 1 init
			point_Delete.push_back(bc1);
			point_Delete.push_back(bc2);
			pointActive[bc1] = false;
			pointActive[bc2] = false;
			vector<int> pointNeighbor_bc1 = pointNeighbor[bc1];
			vector<int> pointNeighbor_middlec;//store the structure of new point
			vector<double> middlePointc(3);
			middlePointc[0] = (points[bc1][0] + points[bc2][0]) / 2;
			middlePointc[1] = (points[bc1][1] + points[bc2][1]) / 2;
			middlePointc[2] = (points[bc1][2] + points[bc2][2]) / 2;
			points.push_back(middlePointc);
			pointActive.push_back(true);
			
			//Collapse: 2 define the related point
			//add the related points in the list
			vector<int> pc_list;
			//bc1 related point
			for (int i = 0; i < pointNeighbor[bc1].size(); i++) {
				int p_i = pointNeighbor[bc1][i];
				if (p_i == bc1 || p_i == bc2) {
					continue;
				}
				if (!MeshOptimization_Exist(p_i, pc_list)&& pointActive[p_i]) {
					pc_list.push_back(p_i);
				}
			}
			//bc2 related point
			for (int i = 0; i < pointNeighbor[bc2].size(); i++) {
				int p_i = pointNeighbor[bc2][i];
				if (p_i == bc1 || p_i == bc2) {
					continue;
				}
				if (!MeshOptimization_Exist(p_i, pc_list) && pointActive[p_i]) {
					pc_list.push_back(p_i);
				}
			}

			//Collapse: 3 update the structure
			for (int i = 0; i < pc_list.size(); i++) {
				int p_index = pc_list[i];				
				bool rejudge = false;
				
				for (int j = 0; j < pointNeighbor[p_index].size() / 2; j++) {
					int pc1 = pointNeighbor[p_index][2 * j];
					int pc2 = pointNeighbor[p_index][2 * j + 1];
					
					if ((pc1 == bc1 || pc1 == bc2) && (pc2 == bc1 || pc2 == bc2)) {						
						pointNeighbor[p_index][2 * j] = -1;
						pointNeighbor[p_index][2 * j + 1] = -1;						
						rejudge = true;
					}
					else if (pc1 == bc1 || pc1 == bc2) {
						pointNeighbor[p_index][2 * j] = pointsEndIndex;
						pointNeighbor_middlec.push_back(pc2);
						pointNeighbor_middlec.push_back(p_index);
						//if (!MeshOptimization_Check_Structure(pointNeighbor_middlec)) {
							//cout << "Repeat!" << endl;
						//};
						//if (!MeshOptimization_Check_Structure(pointNeighbor[p_index])) {
							//cout << "Repeat!" << endl;
						//};
					}
					else if (pc2 == bc1 || pc2 == bc2) {
						pointNeighbor[p_index][2 * j + 1] = pointsEndIndex;
						pointNeighbor_middlec.push_back(p_index);
						pointNeighbor_middlec.push_back(pc1);
						//if (!MeshOptimization_Check_Structure(pointNeighbor_middlec)) {
							//cout << "Repeat!" << endl;
						//};
						//if (!MeshOptimization_Check_Structure(pointNeighbor[p_index])) {
							//cout << "Repeat!" << endl;
						//};
					}
					else {
						continue;
					}
				}
				if (rejudge) {
					vector<int> pointNeighbor_new;
					for (int j = 0; j < pointNeighbor[p_index].size(); j++) {
						if (pointNeighbor[p_index][j] != -1) {
							pointNeighbor_new.push_back(pointNeighbor[p_index][j]);
						}
					}
					pointNeighbor[p_index].clear();
					pointNeighbor[p_index] = pointNeighbor_new;
					//if (!MeshOptimization_Check_Structure(pointNeighbor_new)) {
						//cout << "Repeat!" << endl;
					//};
				}
			}

			pointNeighbor[bc1].clear();
			pointNeighbor[bc1].push_back(-1);
			pointNeighbor[bc1].push_back(-1);

			pointNeighbor[bc2].clear();
			pointNeighbor[bc2].push_back(-1);
			pointNeighbor[bc2].push_back(-1);

			pointNeighbor.push_back(pointNeighbor_middlec);
			pointsEndIndex++;

			//store the pair step
			numberProcess++;
		
		}//the while end		

		//Remove Repeat from neighbor
		for (int i = originalNumber; i < pointNeighbor.size(); i++) {
			vector<int> v = pointNeighbor[i];
			vector<int> v_RemoveRepeat = MeshOptimization_Structure_RemoveRepeat(v);
			pointNeighbor[i].clear();
			pointNeighbor[i] = v_RemoveRepeat;
		}			

		//MeshOptimization_Flip(point_Delete);
		MeshOptimization_UpdateStructure(point_Delete);
		return numberProcess;

	}

	//Split, return the sum of point number
	void MeshOptimization_Split() {

		while (1) {

			int Split_list_index = 0;
			int pointsEndIndex = points.size();//the index for last point index, by update point list
			int numberProcess = 0;

			vector<vector<int>> Split_list;
			for (int i = 0; i < pointNeighbor.size(); i++) {
				vector<int> b2_list;
				int b1 = i;
				for (int j = 0; j < pointNeighbor[i].size(); j++) {
					int b2 = pointNeighbor[i][j];
					if (!MeshOptimization_Exist(b2, b2_list) && b1 < b2) {
						b2_list.push_back(b2);
						double length12 = MeshOptimization_EdgeLength(b1, b2);
						if (length12 >= 1.33 * L_ave) {
							vector<int> edgeij;
							edgeij.push_back(b1);
							edgeij.push_back(b2);
							Split_list.push_back(edgeij);
						}
					}
				}
			}

			cout << Split_list.size();
			if (Split_list.size() == 0) {
				cout << ".";
				break;			
			}
			else {
				cout << ",";			
			}

			//Split processing
			for (int i = 0; i < Split_list.size(); i++) {
				int bs1 = Split_list[i][0];
				int bs2 = Split_list[i][1];
				//Split: 1 achieve middle point
				vector<double> middlePoint(3);
				middlePoint[0] = (points[bs1][0] + points[bs2][0]) / 2;
				middlePoint[1] = (points[bs1][1] + points[bs2][1]) / 2;
				middlePoint[2] = (points[bs1][2] + points[bs2][2]) / 2;
				points.push_back(middlePoint);

				vector<int> pointNeighbor_b1 = pointNeighbor[bs1];
				vector<int> pointNeighbor_middle;//store the structure of middle point

				//Split: 2 store related points
				vector<int> bs3_list;
				bs3_list.push_back(bs1);
				bs3_list.push_back(bs2);
				for (int i = 0; i < pointNeighbor_b1.size() / 2; i++) {
					int b1i1 = pointNeighbor_b1[2 * i];
					int b1i2 = pointNeighbor_b1[2 * i + 1];
					if (b1i1 == bs2) {
						bs3_list.push_back(b1i2);
					}
					else if (b1i2 == bs2) {
						bs3_list.push_back(b1i1);
					}
					else {
						continue;
					}
				}

				//Split: 3 update the structure			
				vector<vector<int>> T_n;
				for (int i = 0; i < bs3_list.size(); i++) {
					int b_index = bs3_list[i];
					if (b_index == bs1 || b_index == bs2) {
						for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
							int b_index_j1 = pointNeighbor[b_index][2 * j];
							int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
							if (b_index_j1 == bs1 || b_index_j1 == bs2) {
								pointNeighbor[b_index][2 * j] = pointsEndIndex;
								pointNeighbor_middle.push_back(b_index_j2);
								pointNeighbor_middle.push_back(b_index);
							}
							if (b_index_j2 == bs1 || b_index_j2 == bs2) {
								pointNeighbor[b_index][2 * j + 1] = pointsEndIndex;
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index_j1);
							}
						}
					}
					else {
						for (int j = 0; j < pointNeighbor[b_index].size() / 2; j++) {
							int b_index_j1 = pointNeighbor[b_index][2 * j];
							int b_index_j2 = pointNeighbor[b_index][2 * j + 1];
							if ((b_index_j1 == bs1 || b_index_j1 == bs2) && (b_index_j2 == bs1 || b_index_j2 == bs2)) {
								pointNeighbor[b_index][2 * j] = pointsEndIndex;
								pointNeighbor[b_index].push_back(b_index_j1);
								pointNeighbor[b_index].push_back(pointsEndIndex);
								pointNeighbor_middle.push_back(b_index_j2);
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index);
								pointNeighbor_middle.push_back(b_index_j1);
							}
						}
					}
				}
				//Split: 4 update neighbor structure
				pointNeighbor.push_back(pointNeighbor_middle);
				pointsEndIndex++;
			}

			//Remove Repeat from neighbor
			for (int i = originalNumber; i < pointNeighbor.size(); i++) {
				vector<int> v = pointNeighbor[i];
				vector<int> v_RemoveRepeat = MeshOptimization_Structure_RemoveRepeat(v);
				pointNeighbor[i].clear();
				pointNeighbor[i] = v_RemoveRepeat;
			}		
		}
		cout << endl;
	}

	//Collapse, return the point number which has been deleted from the original point cloud
	void MeshOptimization_Collapse(bool judge) {
		
		int CollapseNum = points.size() - originalNumber;

		while (1) {	

			if (CollapseNum == 0 && judge) {
				break;//control the certain point number			
			}

			int Collapse_Number = 0;

			vector<bool> pointActive(points.size(), true);//state of point, if ture, can be Split-Collapse, alse can not		
			vector<int> point_Delete;//store the index of delete points.
			int pointsEndIndex = points.size();

			//store all edges satisfy the condition of Split-Collapse		
			vector<vector<int>> Collapse_List;

			//cout << "Achieve all edges for Collapse" << endl;
			//achieve all edges for Collapse
			for (int i = 0; i < pointNeighbor.size(); i++) {				

				vector<int> b2_list;
				int b1 = i;
				for (int j = 0; j < pointNeighbor[i].size(); j++) {
					int b2 = pointNeighbor[i][j];
					if (!MeshOptimization_Exist(b2, b2_list) && b1 < b2) {
						b2_list.push_back(b2);
						double length12 = MeshOptimization_EdgeLength(b1, b2);
						if (length12 < 0.8 * L_ave) {
							//Test if there have invaild condition
							int valence = MeshOptimization_Edge_Valence(b1, b2);
							if (valence <= 2) {
								vector<int> edgeij;
								edgeij.push_back(b1);
								edgeij.push_back(b2);
								Collapse_List.push_back(edgeij);
							}
						}
					}
				}
			}			

			for (int i = 0; i < Collapse_List.size(); i++) {
				//*****************************Collapse test**********************************			
				int bc1 = -1;
				int bc2 = -1;
				int bc1_i = Collapse_List[i][0];
				int bc2_i = Collapse_List[i][1];
				int valenceE = MeshOptimization_Edge_Valence(bc1_i, bc2_i);
				if (pointActive[bc1_i] && pointActive[bc2_i] && valenceE <= 2) {
					bc1 = bc1_i;
					bc2 = bc2_i;
				}
				else {
					continue;
				}			

				vector<int> pointNeighbor_bc1 = pointNeighbor[bc1];
				vector<int> pointNeighbor_middlec;//store the structure of new point
				vector<double> middlePointc(3);
				middlePointc[0] = (points[bc1][0] + points[bc2][0]) / 2;
				middlePointc[1] = (points[bc1][1] + points[bc2][1]) / 2;
				middlePointc[2] = (points[bc1][2] + points[bc2][2]) / 2;

				//Collapse test
				//add the related points in the list
				vector<int> pc_list;
				//bc1 related point
				for (int i = 0; i < pointNeighbor[bc1].size(); i++) {
					int p_i = pointNeighbor[bc1][i];
					if (p_i == bc1 || p_i == bc2) {
						continue;
					}
					if (!MeshOptimization_Exist(p_i, pc_list) && pointActive[p_i]) {
						pc_list.push_back(p_i);
					}
				}
				//bc2 related point
				for (int i = 0; i < pointNeighbor[bc2].size(); i++) {
					int p_i = pointNeighbor[bc2][i];
					if (p_i == bc1 || p_i == bc2) {
						continue;
					}
					if (!MeshOptimization_Exist(p_i, pc_list) && pointActive[p_i]) {
						pc_list.push_back(p_i);
					}
				}

				//check the long edge border length should not longer than 4/3
				bool Collapse_Invaild = false;
				for (int i = 0; i < pc_list.size(); i++) {
					vector<double> pN2 = points[pc_list[i]];
					double radius12 = sqrt((middlePointc[0] - pN2[0]) * (middlePointc[0] - pN2[0]) +
						(middlePointc[1] - pN2[1]) * (middlePointc[1] - pN2[1]) +
						(middlePointc[2] - pN2[2]) * (middlePointc[2] - pN2[2]));
					if (radius12 > 1.33 * L_ave) {
						Collapse_Invaild = true;
						break;
					}
				}
				if (Collapse_Invaild) {
					continue;
				}				

				Collapse_Number++;
				//Collapse processing;
				//Collapse: 1 init
				point_Delete.push_back(bc1);
				point_Delete.push_back(bc2);
				pointActive[bc1] = false;
				pointActive[bc2] = false;				
				points.push_back(middlePointc);
				pointActive.push_back(true);

				//Collapse: 2 update the structure
				for (int i = 0; i < pc_list.size(); i++) {
					int p_index = pc_list[i];
					bool rejudge = false;

					for (int j = 0; j < pointNeighbor[p_index].size() / 2; j++) {
						int pc1 = pointNeighbor[p_index][2 * j];
						int pc2 = pointNeighbor[p_index][2 * j + 1];

						if ((pc1 == bc1 || pc1 == bc2) && (pc2 == bc1 || pc2 == bc2)) {
							pointNeighbor[p_index][2 * j] = -1;
							pointNeighbor[p_index][2 * j + 1] = -1;
							rejudge = true;
						}
						else if (pc1 == bc1 || pc1 == bc2) {
							pointNeighbor[p_index][2 * j] = pointsEndIndex;
							pointNeighbor_middlec.push_back(pc2);
							pointNeighbor_middlec.push_back(p_index);
						}
						else if (pc2 == bc1 || pc2 == bc2) {
							pointNeighbor[p_index][2 * j + 1] = pointsEndIndex;
							pointNeighbor_middlec.push_back(p_index);
							pointNeighbor_middlec.push_back(pc1);
						}
						else {
							continue;
						}
					}
					if (rejudge) {
						vector<int> pointNeighbor_new;
						for (int j = 0; j < pointNeighbor[p_index].size(); j++) {
							if (pointNeighbor[p_index][j] != -1) {
								pointNeighbor_new.push_back(pointNeighbor[p_index][j]);
							}
						}
						pointNeighbor[p_index].clear();
						pointNeighbor[p_index] = pointNeighbor_new;
					}
				}
				pointNeighbor[bc1].clear();
				pointNeighbor[bc1].push_back(-1);
				pointNeighbor[bc1].push_back(-1);

				pointNeighbor[bc2].clear();
				pointNeighbor[bc2].push_back(-1);
				pointNeighbor[bc2].push_back(-1);

				pointNeighbor.push_back(pointNeighbor_middlec);
				pointsEndIndex++;
				CollapseNum--;
				if (CollapseNum == 0) {
					break;//control the certain point number			
				}
			}

			cout << Collapse_Number;
			if (Collapse_Number <= 0) {
				cout << ".";
				break;			
			}
			else {
				cout << ",";			
			}
			//Remove Repeat from neighbor
			for (int i = originalNumber; i < pointNeighbor.size(); i++) {
				vector<int> v = pointNeighbor[i];
				vector<int> v_RemoveRepeat = MeshOptimization_Structure_RemoveRepeat(v);
				pointNeighbor[i].clear();
				pointNeighbor[i] = v_RemoveRepeat;
			}
			MeshOptimization_UpdateStructure(point_Delete);		
		}
		cout << endl;
		
    }
	
	//Flip
	void MeshOptimization_Flip() {
		while (1) {
			int FilpNum = 0;
			for (int i = 0; i < points.size(); i++) {
				//Flip test					
				int b1 = i;
				//cout << b1;			
				//Flip start
				vector<int> pNb1;//add points into list
				for (int j = 0; j < pointNeighbor[b1].size(); j++) {
					int b2_j = pointNeighbor[b1][j];
					if (!MeshOptimization_Exist(b2_j, pNb1)) {
						pNb1.push_back(b2_j);
					}
				}

				for (int j = 0; j < pNb1.size(); j++) {
					int b2 = pNb1[j];
					if (b2 < b1) {
						continue;
					}
					//cout << b2;
					vector<int> b3List;//store the other two points for flip
					for (int k = 0; k < pointNeighbor[b1].size() / 2; k++) {
						int b11 = pointNeighbor[b1][2 * k];
						int b12 = pointNeighbor[b1][2 * k + 1];
						if (b11 == b2 && !MeshOptimization_Exist(b12, b3List)) {
							b3List.push_back(b12);
						}
						if (b12 == b2 && !MeshOptimization_Exist(b11, b3List)) {
							b3List.push_back(b11);
						}
					}
					if (b3List.size() == 2) {
						//Flip processing
						//compute 4 point valence
						int b3 = b3List[0];
						int b4 = b3List[1];
						//cout <<"Flip:"<< b1 << "," << b2 << "," << b3 << "," << b4 << endl;
						if (MeshOptimization_If_Flip(b1, b2, b3, b4)) {
							//Flip update
							//cout << b1 << "," << b2 << "," << b3 << "," << b4 << endl;
							MeshOptimization_Flip_Processing(b1, b2, b3, b4);
							FilpNum++;
						}
					}
					else {
						continue;
					}
				}
			}
			cout << FilpNum;
			if (FilpNum == 0) {
				cout << ".";
				break;			
			}
			else {
				cout << ",";
			}
		}
		cout << endl;		
	}

	//Flip related edge
	void MeshOptimization_Flip_old(vector<int> point_Delete) {

		for (int i = originalNumber; i < points.size(); i++) {
			//Flip test
			//cout << i << endl;			
			int b1 = i;	
			if (MeshOptimization_Exist(b1, point_Delete)) {
				continue;			
			}
			//Flip start
			vector<int> pNb1;//add points into list
			for (int j = 0; j < pointNeighbor[b1].size(); j++) {
				int b2_j = pointNeighbor[b1][j];
				if (!MeshOptimization_Exist(b2_j, pNb1)) {
					pNb1.push_back(b2_j);
				}
			}

			for (int j = 0; j < pNb1.size(); j++) {
				int b2 = pNb1[j];
				if (MeshOptimization_Exist(b2, point_Delete)) {
					continue;
				}
				vector<int> b3List;//store the other two points for flip
				for (int k = 0; k < pointNeighbor[b1].size() / 2; k++) {
					int b11 = pointNeighbor[b1][2 * k];
					int b12 = pointNeighbor[b1][2 * k + 1];	
					if (b11 == b2 && !MeshOptimization_Exist(b12, b3List)) {
						b3List.push_back(b12);
					}
					if (b12 == b2 && !MeshOptimization_Exist(b11, b3List)) {
						b3List.push_back(b11);					
					}
				}
				if (b3List.size() == 2) {
					//Flip processing
					//compute 4 point valence
					int b3 = b3List[0];
					int b4 = b3List[1];	
					if (MeshOptimization_Exist(b3, point_Delete)|| MeshOptimization_Exist(b4, point_Delete)) {
						continue;
					}
					if (MeshOptimization_If_Flip(b1, b2, b3, b4)) {
						//Flip update
						cout << b1 << "," << b2 << "," << b3 << "," << b4 << endl;						
						MeshOptimization_Flip_Processing(b1, b2, b3, b4);
					}				
				}
				else {
					continue;				
				}
			}
		}	
	}

	//Update structure
	void MeshOptimization_UpdateStructure(vector<int> pointDelete) {
		
		vector<int> pointsIndex(points.size());
		vector<bool> pointjudge(points.size(), true);
		for (int i = 0; i < pointsIndex.size(); i++) {
			pointsIndex[i] = i;		
		}		
		
		for (int i = 0; i < pointDelete.size(); i++) {
			int p_index = pointDelete[i];
			pointjudge[p_index] = false;
			pointsIndex[p_index] = -1;
			for (int j = p_index + 1; j < pointsIndex.size(); j++) {
				pointsIndex[j]--;			
			}	
		}
		
		//update point list
		vector<vector<double>> pointNew;
		vector<vector<int>> pointNeighborNew;
		vector<vector<int>> faceInforNew;
		for (int i = 0; i < pointjudge.size(); i++) {
			if (pointjudge[i]) {
				pointNew.push_back(points[i]);	
				vector<int> pointNeighborNew_i = pointNeighbor[i];
				for (int j = 0; j < pointNeighborNew_i.size(); j++) {
					pointNeighborNew_i[j] = pointsIndex[pointNeighborNew_i[j]];
					//int index = pointNeighbor[i][j];
					//pointNeighborNew_i.push_back(pointsIndex[index]);
				}
				pointNeighborNew.push_back(pointNeighborNew_i);
			}		
		}
		points.clear();
		points = pointNew;
		pointNeighbor.clear();
		pointNeighbor = pointNeighborNew;
	}

	//Finally output the 
	void MeshOptimization_UpdateStructure_Face() {

		vector<vector<int>> pointNeighborNew = pointNeighbor;
		vector<vector<int>> faceInforNew;
		//update face list
		for (int i = 0; i < pointNeighborNew.size(); i++) {
			int b1 = i;
			vector<int> b1_N = pointNeighborNew[b1];
			for (int j = 0; j < pointNeighborNew[i].size() / 2; j++) {
				int b2 = pointNeighborNew[i][2 * j];
				int b3 = pointNeighborNew[i][2 * j + 1];
				if (b2 == -1 || b3 == -1) {
					continue;
				}
				else {
					//f1 << "f " << i + 1 << " " << pointNeighborNew[i][2 * j] + 1 << " " << pointNeighborNew[i][2 * j + 1] + 1 << endl;
					vector<int> faceInforNew_i;
					faceInforNew_i.push_back(i);
					faceInforNew_i.push_back(pointNeighborNew[i][2 * j]);
					faceInforNew_i.push_back(pointNeighborNew[i][2 * j + 1]);
					faceInforNew.push_back(faceInforNew_i);
					for (int k = 0; k < pointNeighborNew[b2].size() / 2; k++) {
						int b2n2 = pointNeighborNew[b2][2 * k];
						int b2n3 = pointNeighborNew[b2][2 * k + 1];
						if ((b2n2 == b1 || b2n2 == b3) && (b2n3 == b1 || b2n3 == b3)) {
							pointNeighborNew[b2][2 * k] = -1;
							pointNeighborNew[b2][2 * k + 1] = -1;
						}
					}
					for (int k = 0; k < pointNeighborNew[b3].size() / 2; k++) {
						int b3n2 = pointNeighborNew[b3][2 * k];
						int b3n3 = pointNeighborNew[b3][2 * k + 1];
						if ((b3n2 == b1 || b3n2 == b2) && (b3n3 == b1 || b3n3 == b2)) {
							pointNeighborNew[b3][2 * k] = -1;
							pointNeighborNew[b3][2 * k + 1] = -1;
						}
					}
				}
			}
		}
		faceInfor.clear();
		faceInfor = faceInforNew;	
	}

	//*******************************************Basic Function*****************************************************
	//save mesh
	void MeshOptimization_SaveOBJ() {

		string fileName = "Remesh\\optimization\\" + fileName + ".obj";

		ofstream f1(fileName);
		//ofstream f1(fileName, ios::app);//write in the end

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}

		for (int i = 0; i < pointNeighbor.size(); i++) {
			int b1 = i;
			vector<int> b1_N = pointNeighbor[b1];
			for (int j = 0; j < pointNeighbor[i].size()/2; j++) {
				int b2 = pointNeighbor[i][2 * j];
				int b3 = pointNeighbor[i][2 * j + 1];
				if (b2 == -1 || b3 == -1) {
					continue;
				}
				else {
					f1 << "f " << i + 1 << " " << pointNeighbor[i][2 * j] + 1 << " " << pointNeighbor[i][2 * j + 1] + 1 << endl;
					for (int k = 0; k < pointNeighbor[b2].size()/2; k++) {
						int b2n2 = pointNeighbor[b2][2 * k];
						int b2n3 = pointNeighbor[b2][2 * k + 1];
						if ((b2n2 == b1 || b2n2 == b3) && (b2n3 == b1 || b2n3 == b3)) {
							pointNeighbor[b2][2 * k] = -1;
							pointNeighbor[b2][2 * k + 1] = -1;
						}					
					}
					for (int k = 0; k < pointNeighbor[b3].size() / 2; k++) {
						int b3n2 = pointNeighbor[b3][2 * k];
						int b3n3 = pointNeighbor[b3][2 * k + 1];
						if ((b3n2 == b1 || b3n2 == b2) && (b3n3 == b1 || b3n3 == b2)) {
							pointNeighbor[b3][2 * k] = -1;
							pointNeighbor[b3][2 * k + 1] = -1;
						}
					}
				}				
			}			
		}

		f1.close();


	}
	
	void MeshOptimization_SaveOBJ_Face() {

		string fin = "Remesh\\optimization\\" + fileName+"_" + to_string(optIter) + ".obj";

		ofstream f1(fin);

		//f1 << points.size() << " " << facet.size() << " " << 0 << endl;

		for (int i = 0; i < points.size(); i++) {
			f1 << "v " << points[i][0] << " " << points[i][1] << " " << points[i][2] << " " << endl;
		}

		for (int i = 0; i < faceInfor.size(); i++) {
			f1 << "f " << faceInfor[i][0]+1 << " " << faceInfor[i][1]+1 << " " << faceInfor[i][2]+1 << endl;
		}

		f1.close();


	}
	//compute ave edge length of a mesh
	void MeshOptimization_AveEdge_init() {//computer the ave edge l for trangular edge estimation

		//select 1000 trangulars for estimation
		
		double l_sum = 0;
		for (int i = 0; i < faceInfor.size(); i++) {			
			int b1 = faceInfor[i][0];
			int b2 = faceInfor[i][1];
			int b3 = faceInfor[i][2];
			double l1 = MeshOptimization_EdgeLength(b1, b2);
			double l2 = MeshOptimization_EdgeLength(b2, b3);
			double l3 = MeshOptimization_EdgeLength(b3, b1);
			double l_ave_i = (l1 + l2 + l3) / 3;
			l_sum = l_sum + l_ave_i;
		}
		L_ave = l_sum / faceInfor.size();
	}

	//return the length of an edge (b1, b2)
	double MeshOptimization_EdgeLength(int b1, int b2) {

		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];

		double lengthE = sqrt((p1[0] - p2[0]) * (p1[0] - p2[0]) +
			(p1[1] - p2[1]) * (p1[1] - p2[1]) + (p1[2] - p2[2]) * (p1[2] - p2[2]));

		return lengthE;

	}

	//judge the point is in a list
	bool MeshOptimization_Exist(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return true;
			}
		}
		return false;
	}

	int MeshOptimization_Exist_Index(int b, vector<int> v) {

		for (int i = 0; i < v.size(); i++) {
			if (b == v[i]) {
				return i;
			}
		}
		return -1;
	}

	//return the Valence value of a edge (p1, p2). 1: edge; 2: normal; 3: invalid
	int MeshOptimization_Edge_Valence(int p1, int p2) {
		int valence = 0;
		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {
			if (pointNeighbor[p1][2 * i] == p2 || pointNeighbor[p1][2 * i + 1] == p2) {
				valence++;
			}
		}
		return valence;
	}

	//return the Valence value of a point
	int MeshOptimization_Point_Valence(int p1) {

		
		vector<int> b_num;
		for (int i = 0; i < pointNeighbor[p1].size() / 2; i++) {
			int p2 = pointNeighbor[p1][2 * i];
			int p3 = pointNeighbor[p1][2 * i + 1];
			if (!MeshOptimization_Exist(p2, b_num)) {
				b_num.push_back(p2);			
			}
			if (!MeshOptimization_Exist(p3, b_num)) {
				b_num.push_back(p3);
			}
		}
		return b_num.size();

	}

	//return if a point is an edge point
	bool MeshOptimization_Edge_Point(int p1) {

		for (int i = 0; i < pointNeighbor[p1].size(); i++) {
			int p2 = pointNeighbor[p1][i];
			int valenceE = MeshOptimization_Edge_Valence(p1, p2);
			if (valenceE == 1) {
				return true;			
			}								
		}
		return false;
	}

	//Flip judge and processing
	bool MeshOptimization_If_Flip(int b1, int b2, int b3, int b4) {

		double disb34 = MeshOptimization_EdgeLength(b3, b4);
		if (disb34 > 1.33 * L_ave || disb34 < 0.8 * L_ave) {
			return false;		
		}
		bool b1e = MeshOptimization_Edge_Point(b1);
		bool b2e = MeshOptimization_Edge_Point(b2);
		bool b3e = MeshOptimization_Edge_Point(b3);
		bool b4e = MeshOptimization_Edge_Point(b4);
		int b1ValenceTarget = 6;
		if (b1e) {
			b1ValenceTarget = 4;		
		}
		int b2ValenceTarget = 6;
		if (b2e) {
			b2ValenceTarget = 4;
		}
		int b3ValenceTarget = 6;
		if (b3e) {
			b3ValenceTarget = 4;
		}
		int b4ValenceTarget = 6;
		if (b4e) {
			b4ValenceTarget = 4;
		}
		int b1ValenceBefore = MeshOptimization_Point_Valence(b1);
		int b2ValenceBefore = MeshOptimization_Point_Valence(b2);
		int b3ValenceBefore = MeshOptimization_Point_Valence(b3);
		int b4ValenceBefore = MeshOptimization_Point_Valence(b4);
		int b1ValenceAfter = b1ValenceBefore - 1;
		int b2ValenceAfter = b2ValenceBefore - 1;
		int b3ValenceAfter = b3ValenceBefore + 1;
		int b4ValenceAfter = b4ValenceBefore + 1;

		int BeforeValence = 0;
		BeforeValence = abs(b1ValenceBefore - b1ValenceTarget) +
			abs(b2ValenceBefore - b2ValenceTarget) +
			abs(b3ValenceBefore - b3ValenceTarget) +
			abs(b4ValenceBefore - b4ValenceTarget);

		int AfterValence = 0;
		AfterValence = abs(b1ValenceAfter - b1ValenceTarget) +
			abs(b2ValenceAfter - b2ValenceTarget) +
			abs(b3ValenceAfter - b3ValenceTarget) +
			abs(b4ValenceAfter - b4ValenceTarget);
		//object Valence	

		if (AfterValence < BeforeValence) {
			return true;		
		}
		else {
			return false;		
		}	
	}

	void MeshOptimization_Flip_Processing(int b1, int b2, int b3, int b4) {

		vector<int> b1N = pointNeighbor[b1];
		vector<int> b2N = pointNeighbor[b2];
		vector<int> b3N = pointNeighbor[b3];
		vector<int> b4N = pointNeighbor[b4];

		//update b1 structure
		vector<int> b1N_New;
		for (int i = 0; i < b1N.size() / 2; i++) {
			int b1N2 = b1N[2 * i];
			int b1N3 = b1N[2 * i + 1];
			if (b1N2 == b2 && (b1N3 == b3 || b1N3 == b4)) {
				b1N[2 * i] = -1;
				b1N[2 * i + 1] = -1;
				if (b1N3 == b3) {
					b1N_New.push_back(b4);
					b1N_New.push_back(b3);				
				}
				else {
					b1N_New.push_back(b3);
					b1N_New.push_back(b4);				
				}
			}
			if (b1N3 == b2 && (b1N2 == b3 || b1N2 == b4)) {
				b1N[2 * i] = -1;
				b1N[2 * i + 1] = -1;
				if (b1N2 == b3) {
					b1N_New.push_back(b3);
					b1N_New.push_back(b4);
				}
				else {
					b1N_New.push_back(b4);
					b1N_New.push_back(b3);
				}			
			}		
		}
		if (b1N_New.size() == 0) {
			cout << "Hello!";		
		}
		vector<int> b1N_final;
		for (int i = 0; i < b1N.size() / 2; i++) {
			int b1N2 = b1N[2 * i];
			int b1N3 = b1N[2 * i + 1];
			if (b1N2 != -1 && b1N3 != -1) {
				b1N_final.push_back(b1N2);
				b1N_final.push_back(b1N3);			
			}		
		}
		b1N_final.push_back(b1N_New[0]);
		b1N_final.push_back(b1N_New[1]);

		//update b2 structure
		vector<int> b2N_New;
		for (int i = 0; i < b2N.size() / 2; i++) {
			int b2N2 = b2N[2 * i];
			int b2N3 = b2N[2 * i + 1];
			if (b2N2 == b1 && (b2N3 == b3 || b2N3 == b4)) {
				b2N[2 * i] = -1;
				b2N[2 * i + 1] = -1;
				if (b2N3 == b3) {
					b2N_New.push_back(b4);
					b2N_New.push_back(b3);
				}
				else {
					b2N_New.push_back(b3);
					b2N_New.push_back(b4);
				}
			}
			if (b2N3 == b1 && (b2N2 == b3 || b2N2 == b4)) {
				b2N[2 * i] = -1;
				b2N[2 * i + 1] = -1;
				if (b2N2 == b3) {
					b2N_New.push_back(b3);
					b2N_New.push_back(b4);
				}
				else {
					b2N_New.push_back(b4);
					b2N_New.push_back(b3);
				}
			}
		}
		if (b2N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b2N_final;
		for (int i = 0; i < b2N.size() / 2; i++) {
			int b2N2 = b2N[2 * i];
			int b2N3 = b2N[2 * i + 1];
			if (b2N2 != -1 && b2N3 != -1) {
				b2N_final.push_back(b2N2);
				b2N_final.push_back(b2N3);
			}
		}
		b2N_final.push_back(b2N_New[0]);
		b2N_final.push_back(b2N_New[1]);

		//update b3 structure	
		vector<int> b3N_New;
		for (int i = 0; i < b3N.size() / 2; i++) {
			int b3N2 = b3N[2 * i];
			int b3N3 = b3N[2 * i + 1];
			if ((b3N2 == b1 || b3N2 == b2) && (b3N3 == b1 || b3N3 == b2)) {
				b3N[2 * i] = -1;
				b3N[2 * i + 1] = -1;
				b3N_New.push_back(b3N2);
				b3N_New.push_back(b4);
				b3N_New.push_back(b4);
				b3N_New.push_back(b3N3);
			}			
		}
		if (b3N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b3N_final;
		for (int i = 0; i < b3N.size() / 2; i++) {
			int b3N2 = b3N[2 * i];
			int b3N3 = b3N[2 * i + 1];
			if (b3N2 != -1 && b3N3 != -1) {
				b3N_final.push_back(b3N2);
				b3N_final.push_back(b3N3);
			}
		}
		b3N_final.push_back(b3N_New[0]);
		b3N_final.push_back(b3N_New[1]);
		b3N_final.push_back(b3N_New[2]);
		b3N_final.push_back(b3N_New[3]);

		//update b4 structure	
		vector<int> b4N_New;
		for (int i = 0; i < b4N.size() / 2; i++) {
			int b4N2 = b4N[2 * i];
			int b4N3 = b4N[2 * i + 1];
			if ((b4N2 == b1 || b4N2 == b2) && (b4N3 == b1 || b4N3 == b2)) {
				b4N[2 * i] = -1;
				b4N[2 * i + 1] = -1;
				b4N_New.push_back(b4N2);
				b4N_New.push_back(b3);
				b4N_New.push_back(b3);
				b4N_New.push_back(b4N3);
			}
		}
		if (b4N_New.size() == 0) {
			cout << "Hello!";
		}
		vector<int> b4N_final;
		for (int i = 0; i < b4N.size() / 2; i++) {
			int b4N2 = b4N[2 * i];
			int b4N3 = b4N[2 * i + 1];
			if (b4N2 != -1 && b4N3 != -1) {
				b4N_final.push_back(b4N2);
				b4N_final.push_back(b4N3);
			}
		}
		b4N_final.push_back(b4N_New[0]);
		b4N_final.push_back(b4N_New[1]);
		b4N_final.push_back(b4N_New[2]);
		b4N_final.push_back(b4N_New[3]);

		pointNeighbor[b1].clear();
		pointNeighbor[b2].clear();
		pointNeighbor[b3].clear();
		pointNeighbor[b4].clear();

		pointNeighbor[b1] = b1N_final;
		pointNeighbor[b2] = b2N_final;
		pointNeighbor[b3] = b3N_final;
		pointNeighbor[b4] = b4N_final;

	}

	//Relocate vertices on the surface by tangential smoothing
	void MeshOptimization_TangentSmoothing() {
		
		double landa = 0.6;
		vector<vector<double>> pointNew = points;

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			vector<double> pi = points[i];
			int p1 = i;
			if (MeshOptimization_Edge_Point(p1)) {
				pointNew[i] = pi;
				continue;			
			}
			vector<int> pointsN = pointNeighbor[i];
			vector<int> pointsN_withoutzRepeat;//neighbor points without repeat point
			vector<double> pointN_Area;//correspoinding area of neighbor point
			double AreaSum = 0;
			for (int j = 0; j < pointsN.size() / 2; j++) {
				int p2 = pointsN[2 * j];
				int p3 = pointsN[2 * j + 1];	
				double areai = MeshOptimization_TrangularArea(p1, p2, p3);				
				//double areai2 = areai / (double)2;
				int indexp2 = MeshOptimization_Exist_Index(p2, pointsN_withoutzRepeat);
				int indexp3 = MeshOptimization_Exist_Index(p3, pointsN_withoutzRepeat);
				if (indexp2 == -1) {
					pointsN_withoutzRepeat.push_back(p2);
					pointN_Area.push_back(areai);
				}else{
					pointN_Area[indexp2] = (pointN_Area[indexp2] + areai)/2;
				}
				if (indexp3 == -1) {
					pointsN_withoutzRepeat.push_back(p3);
					pointN_Area.push_back(areai);
				}
				else {
					pointN_Area[indexp3] = (pointN_Area[indexp3] + areai)/2;
				}
			}

			for (int j = 0; j < pointN_Area.size(); j++) {

				AreaSum = AreaSum + pointN_Area[j];
			
			}
		
			vector<double> gi(3, 0);
			for (int j = 0; j < pointsN_withoutzRepeat.size(); j++) {
				gi[0] = gi[0] + (pointN_Area[j] / AreaSum) * points[pointsN_withoutzRepeat[j]][0];
				gi[1] = gi[1] + (pointN_Area[j] / AreaSum) * points[pointsN_withoutzRepeat[j]][1];
				gi[2] = gi[2] + (pointN_Area[j] / AreaSum) * points[pointsN_withoutzRepeat[j]][2];
			}

			//achieve normal of p1
			vector<double> ni = MeshOptimization_Normal_Point(p1);
			//achieve move vertor v_pg for gi and pi
			vector<double> v_pg(3);
			v_pg[0] = gi[0] - pi[0];
			v_pg[1] = gi[1] - pi[1];
			v_pg[2] = gi[2] - pi[2];
			double mapN_dis = v_pg[0] * ni[0] + v_pg[1] * ni[1] + v_pg[2] * ni[2];//v_pg map into n distance.
			vector<double> pi_new(3, 0);
			pi_new[0] = pi[0] + landa * (v_pg[0] - mapN_dis * ni[0]);
			pi_new[1] = pi[1] + landa * (v_pg[1] - mapN_dis * ni[1]);
			pi_new[2] = pi[2] + landa * (v_pg[2] - mapN_dis * ni[2]);
			if (!(pi_new[0]<9999 && pi_new[0]>-9999)) {
				continue;			
			}
			pointNew[i] = pi_new;
		}

#pragma omp parallel for
		for (int i = 0; i < points.size(); i++) {
			points[i][0] = pointNew[i][0];
			points[i][1] = pointNew[i][1];
			points[i][2] = pointNew[i][2];		
		}
	}

	//achieve the trangualr area
	double MeshOptimization_TrangularArea(int b1, int b2, int b3) {
		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];
		vector<double> p3 = points[b3];

		vector<double> v1(3);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];

		vector<double> v2(3);
		v2[0] = p3[0] - p2[0];
		v2[1] = p3[1] - p2[1];
		v2[2] = p3[2] - p2[2];

		vector<double> v3(3);
		v3[0] = p1[0] - p3[0];
		v3[1] = p1[1] - p3[1];
		v3[2] = p1[2] - p3[2];

		double v1d = sqrt(v1[0] * v1[0] + v1[1] * v1[1] + v1[2] * v1[2]);
		double v2d = sqrt(v2[0] * v2[0] + v2[1] * v2[1] + v2[2] * v2[2]);
		double v3d = sqrt(v3[0] * v3[0] + v3[1] * v3[1] + v3[2] * v3[2]);

		double p = (v1d + v2d + v3d) / 2;
		double S = sqrt(p * (p - v1d) * (p - v2d) * (p - v3d));

		return S;
	}

	//normal return
	vector<double> MeshOptimization_Normal_Point(int b1) {

		vector<int> p1N = pointNeighbor[b1];
		vector<double> p1N_Weight;//weights of the normal list
		vector<vector<double>> p1N_Normal;//normal list
		double Area_sum = 0;
		for (int i = 0; i < p1N.size() / 2; i++) {
			int b2 = p1N[2 * i];
			int b3 = p1N[2 * i + 1];
			vector<double> n_i = MeshOptimization_Normal_Face(b1, b2, b3);
			double area_i = MeshOptimization_TrangularArea(b1, b2, b3);
			Area_sum = Area_sum + area_i;
			p1N_Weight.push_back(area_i);
			p1N_Normal.push_back(n_i);
		}
		vector<double> n(3, 0);
		for (int i = 0; i < p1N_Weight.size(); i++) {
			n[0] = n[0] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][0];
			n[1] = n[1] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][1];
			n[2] = n[2] + (p1N_Weight[i] / Area_sum) * p1N_Normal[i][2];
		}
		return n;
	}

	vector<double> MeshOptimization_Normal_Face(int b1, int b2, int b3) {

		vector<double> p1 = points[b1];
		vector<double> p2 = points[b2];
		vector<double> p3 = points[b3];

		vector<double> v1(3, 0);
		v1[0] = p2[0] - p1[0];
		v1[1] = p2[1] - p1[1];
		v1[2] = p2[2] - p1[2];

		vector<double> v2(3, 0);
		v2[0] = p3[0] - p2[0];
		v2[1] = p3[1] - p2[1];
		v2[2] = p3[2] - p2[2];

		vector<double> n(3, 0);
		n[0] = v2[2] * v1[1] - v2[1] * v1[2];
		n[1] = -v2[2] * v1[0] + v2[0] * v1[2];
		n[2] = v2[1] * v1[0] - v2[0] * v1[1];

		vector<double> n_Unit = MeshOptimization_Unit_Normal(n);
		return n_Unit;

	}

	vector<double> MeshOptimization_Unit_Normal(vector<double> n) {

		double length = sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]);
		n[0] = n[0] / length;
		n[1] = n[1] / length;
		n[2] = n[2] / length;
		return n;

	}

	//Check structure and remove the repeat item from the structure
	vector<int> MeshOptimization_Structure_RemoveRepeat(vector<int> v) {

		vector<int> result;
		for (int i = 0; i < v.size() / 2; i++) {
			int b1 = v[2 * i];
			int b2 = v[2 * i + 1];
			if (b1 == -1 || b2 == -1) {
				continue;			
			}
			for (int j = i + 1; j < v.size() / 2; j++) {
				int b11 = v[2 * j];
				int b12 = v[2 * j + 1];
				if (b11 == -1 || b12 == -1) {
					continue;
				}
				if ((b1 == b11 || b1 == b12) && (b2 == b11 || b2 == b12)) {
					v[2 * j] = -1;
					v[2 * j + 1] = -1;
				}
			}		
		}

		for (int i = 0; i < v.size() / 2; i++) {
			int b1 = v[2 * i];
			int b2 = v[2 * i + 1];
			if (b1 == -1 || b2 == -1) {
				continue;
			}
			else {
				result.push_back(b1);
				result.push_back(b2);			
			}		
		}
		return result;
	}

};


