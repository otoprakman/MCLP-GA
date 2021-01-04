#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <cmath>
#include <time.h>       
#include <assert.h>				
#include <fstream>

using namespace std;

float min = 1000;
float max = 0;
const int num_node =708;
const int num_chromosome = 10;
const int num_parents = 4;
const int range = 3;
const int num_facility = 5;
const float	x_size = 30;
const float	y_size = 30;
const int max_demand = 200;
const float mutation_rate = 0.1;
const float replacement_prob = 0.8;
const int num_iter = 10;
const float xover_rate = 0.7;
const int loc_src_const = 3;

struct node {
	double x_cord;
	double y_cord;
	int demand;
	//int index;
	//bool flaged = false;
	int satisfied_demand;
};

node Array_node[num_node];
float distance_node[num_node][num_node];

struct chromosome {
	int soln[num_node];
	int satisfied[num_node];
	int fitness;
	float prob;
	int index;
	bool selected;
};

struct population {
	chromosome chros[num_chromosome];
};

struct parent {
	chromosome parents[num_parents];
};

parent mating_pool;

int rnd_demand_func(int max, int min)
{
	return rand() % (max - min + 1) + min;
}
float rnd_cord_func(int max, int min) {
	return ((float)rand() / RAND_MAX)*(max - min) + min;
}

float rnd_num(int max, int min) {
	return ((float)rand() / RAND_MAX)*(max - min) + min;
}

void rand_node_gen(void) {

	for (int i = 0; i < num_node; i++) {
		Array_node[i].x_cord = rnd_cord_func(x_size, 0); //30x30 area size
		Array_node[i].y_cord = rnd_cord_func(y_size, 0);
		Array_node[i].demand = rnd_demand_func(max_demand, 1);
	}
}

void dist_calc(void) {
	ifstream get_coordinates_x;
	ifstream get_coordinates_y;
	get_coordinates_y.open("sjc708_coord_Y.txt");
	get_coordinates_x.open("sjc708_coord_X.txt");
	for (int i = 0; i < num_node; i++)
	{
		get_coordinates_x >> Array_node[i].x_cord;
		get_coordinates_y >> Array_node[i].y_cord;
	}
	
	for (int i = 0; i < num_node; i++) { //Calculate distances of each nodes from each other
		for (int j = 0; j < num_node; j++) {
			float x = Array_node[i].x_cord - Array_node[j].x_cord;
			float y = Array_node[i].y_cord - Array_node[j].y_cord;
			distance_node[i][j] = sqrt(pow(x, 2) + pow(y, 2));
		}
	}
}

int main() {
	//rand_node_gen();
	dist_calc();

	/*ofstream out("so1800_demand.txt");
	
	for (int i = 0; i < num_node; i++)
	{
		out << Array_node[i].demand << endl;
	}
	*/
	ofstream out2("sjc708_distance.txt");

	for (int i = 0; i < num_node; i++)
	{
		for (int j = 0; j < num_node; j++)
		{
			out2 << distance_node[i][j] << endl;
		}
	}

	/*ofstream out3("Coordinates2.txt");
	for (int i = 0; i < num_node; i++)
	{
		for (int j = 0; j < 2; j++)
		{
			out3 << Array_node[i].x_cord << "\t" << Array_node[i].y_cord << endl;
		}
	}*/
	return 0;
}