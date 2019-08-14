#include <time.h>
#include <ctime>
#include <sstream>
#include <math.h>
#include <stdio.h>
#include <ilcplex/ilocplex.h>



//define global constants
#define infinity 10000000000
#define MODLUS 2147483647
#define MULT1      24112
#define MULT2      26143


bool sirma;
ILOSTLBEGIN
IloEnv env;


int N = 5;  // number of nodes 
int P = 3; // number of facilities that can be opened


IloExprArray expression(env);
IloExpr expr(env);
IloExpr expr2(env);
IloRangeArray scenarioconstsets(env);

IloExpr demand(env);
IloNumArray prob(env, 1);
ofstream output("debug_output.txt");
const int num_node = 100;

struct node
{
	double x_cord;
	double y_cord;
};

int num_node;
node Array_node[num_node];
float distance_node[num_node][num_node];


void dist_calc(void)
{
	for (int i = 0; i < N; i++) { //Calculate distances of each nodes from each other
		for (int j = 0; j < N; j++)
		{
			float x_cord = Array_node[i].x_cord - Array_node[j].x_cord;
			float y_cord = Array_node[i].y_cord - Array_node[j].y_cord;
			distance_node[i][j] = sqrt(pow(x_cord, 2) + pow(y_cord, 2));
		}
	}
}
ofstream distance_node("distances.txt");

string getvarname_x(int j)
{
	char ch[40];

	sprintf_s(ch, "x%d", j);
	return(ch);
}
string getvarname_z(int i)
{
	char ch[40];

	sprintf_s(ch, "z%d", i);
	return(ch);
}

////////////////////////////////////


void const1(IloModel model, IloIntVarArray x)
{
	int j;
	for (j = 0; j < N; j++)
	{
		expr = expr + x[j] - P;
		scenarioconstsets.add(expr == 0);
		expr.clear();
	}
}

void const2(IloModel model, IloIntVarArray z, IloIntVarArray x)
{
	int i, j;
	for (i = 0; i < N; i++)
	{
		expr = expr + z[i];
		for (j = 0; j < N; j++)
		{
			expr2 += x[j];
		}
		expr =expr - expr2;
		scenarioconstsets.add(expr == 0);
		expr.clear();
		expr2.clear();
	}
}

void objfun(IloModel model, const IloIntArray h, IloIntVarArray z)
{

	for (int i = 0; i < N; i++)
	{
		expr = expr + h[i] * z[i];
	}

	expr.clear();
}

void solution_values(IloCplex cplex, IloIntVarArray x, IloIntArray x_var, IloIntVarArray z, IloIntArray z_var)

{
	env.out() << "Solution status = " << cplex.getStatus() << endl;
	env.out() << "Solution value = " << cplex.getObjValue() << endl;
	int i, j;
	for (i = 0; i < N; i++)
	{
		x_var[i] = cplex.getIntValue(x[i]);
		output << "x" << i << "=" << " " << x_var[i] << endl;
		z_var[i] = cplex.getIntValue(z[i]);
		output << "z" << i << "=" << " " << z_var[i] << endl;

	}
}
void model(const IloIntArray h, IloIntArray x, IloIntArray z)
{
		int i, j;
		string var_name;
		IloTimer elapsed_time(env);
		IloNum objfun; //objective value for DEM

		IloIntArray x_var(env, N);
		IloIntArray z_var(env, N);
		elapsed_time.start();
		IloIntVarArray x(env, N, 0, 1);
		IloIntVarArray z(env, N, 0, 1);

		for (i = 0; i < N; i++)
		{
			var_name = getvarname_x(i);
			x[i].setName(var_name.c_str());
			var_name = getvarname_z(i);
			z[i].setName(var_name.c_str());
		}

		IloModel model(env);
		{
			//constraints
			const1(model, x);
			const2(model, z, x);
			objfun(model, h, z);
			expr.clear();
		}


		IloObjective objective = IloMaximize(env, demand);
		model.add(objective);
		IloCplex cplex(model);
		//		cplex.setOut(env.getNullStream()); // turn off the cplex screen outputs
		cplex.exportModel("Facilitymax.lp");
		//	cplex.setParam(cplex.EpGap, 0.0008); // for testing purposes
		//      cplex.setParam(cplex.TiLim, 3600); 
		cplex.solve();
		//		cout << "elapsed time: "<< elapsed_time.getTime() << endl;
		elapsed_time.reset();
		solution_values(cplex, x, x_var, z, z_var);
		scenarioconstsets.end();
		demand.end();
		cplex.end();
		objective.end();
		model.end();
}


int main(int, char**)
{
	int i, j, ind;

	//////////////////////////////////////////////////////////////////////////////
	// now define all parameters necessary for the model
	IloIntArray h(env, 1); // demand amounts
	IloIntArray d(env, 1);  // distance between nodes

	h.setSize(N);
	d.setSize(N*N);
	////////////////////////////////////////////////////////

	ind = 0;
	ifstream in;
	ifstream in2;
	in.open("demand_amonuts.txt");

	for (i = 0; i < N; i++)
	{
		in >> h[i];
		//		cout << h[i]<<endl;
		ind++;
	}

in.close();
ind = 0;

in2.open("distances.txt");
for (i = 0; i < N; i++)
{
	for (j = 0; j < N; j++)
	{
		in2 >> d[i, j];
		ind++;
	}
}
in2.close();


model(h);
env.end();
return 0;
}