#include <iostream>
#include <iomanip>
#include <cstring>
#include <cmath>
#include <fstream>

using namespace std;

class vec
{
public:
	double x[3];
	vec & operator=(const vec& B)
	{
	    x[0] = B.x[0];
	    x[1] = B.x[1];
	    x[2] = B.x[2];
	    return *this;
	}
	vec & operator=(double* B)
	{
	    x[0] = B[0];
	    x[1] = B[1];
	    x[2] = B[2];
	    return *this;
	}
	double operator*(const vec& B)
	{
		double res;
		res = x[0]*B.x[0] + x[1]*B.x[1] + x[2]*B.x[2];
		return res;
	}
	vec operator+(const vec& B)
	{
		vec res;
		res.x[0] = x[0] + B.x[0];
		res.x[1] = x[1] + B.x[1];
		res.x[2] = x[2] + B.x[2];
		return res;
	}
	vec operator-(const vec& B)
	{
		vec res;
		res.x[0] = x[0] - B.x[0];
		res.x[1] = x[1] - B.x[1];
		res.x[2] = x[2] - B.x[2];
		return res;
	}
	vec operator*(double B)
	{
		vec res;
		res.x[0] = x[0]*B;
		res.x[1] = x[1]*B;
		res.x[2] = x[2]*B;
		return res;
	}
	vec operator/(double B)
	{
		vec res;
		res.x[0] = x[0]/B;
		res.x[1] = x[1]/B;
		res.x[2] = x[2]/B;
		return res;
	}
	double norm()
	{
		double res;
		res = x[0]*x[0] + x[1]*x[1] + x[2]*x[2];
		return sqrt(res);
	}
	void print()
	{
		cout<<x[0]<<'\t'<<x[1]<<'\t'<<x[2]<<endl;
	}
	void print(ofstream& out)
	{
		out<<setw(14)<<setprecision(9)<<fixed<<x[0]<<setw(14)<<setprecision(9)<<fixed<<x[1]<<setw(14)<<setprecision(9)<<fixed<<x[2]<<endl;
	}
};

int main ()
{
	string flag_cell="space primitive vectors";
	string flag_coord="reduced coordinates";
	string flag_force="cartesian forces";
	string flag_stress="Cartesian components of stress tensor";

	string flag_nat="mqgrid", check_nat="natom";
	string flag_nel="ntime", check_nel="ntypat";
	string flag_atom="tolmxf", check_atom="typat";
	string flag_el="opening atomic psp file";
	const double b2a = 0.529177249;

	ofstream out_cell, out_coord;

	string tmp;
	int nat, nel;
	int *at_l;
	double *el_l;	//used for mass
	vec * coord, *fcoord;
	vec * force;
	vec cell[3];
	vec cell_new[3], cell_norm_new, cell_unit_new[3];
	vec stress[3];
	double tot_ene;
	const double lambda_f = 40.0;
	const double lambda_t = 50.0;

	out_cell.open("cell.dat");
	out_coord.open("coord.dat");

	//find number of atoms
	getline(cin,tmp);
	while(tmp.find(flag_nat) == string::npos)
		getline(cin,tmp);
	cin>>tmp;
	if (tmp != check_nat)
	{
		cout<<"Error, could not find natom"<<endl;
		exit(0);
	}
	else
	{
		cin>>tmp>>nat;
	}
	//find number of elements
	while(tmp.find(flag_nel) == string::npos)
		getline(cin,tmp);
	cin>>tmp;
	if (tmp != check_nel)
	{
		cout<<"Error, could not find ntypat"<<endl;
		exit(0);
	}
	else
	{
		cin>>nel;
	}
	//define array of atoms and elements;
	at_l = new int[nat];
	el_l = new double[nel];
	coord = new vec[nat];
	fcoord = new vec[nat];
	force = new vec[nat];
	//find atom type
	while(tmp.find(flag_atom) == string::npos)
		getline(cin,tmp);
	cin>>tmp;
	if (tmp != check_atom)
	{
		cout<<"Error, could not find typat"<<endl;
		exit(0);
	}
	else
	{
		for (int t1=0; t1<nat; t1++)
		{
			cin>>at_l[t1];
			at_l[t1]--;
		}
	}
	// find atom position and cell paramter
	while(tmp.find(flag_cell) == string::npos)
		getline(cin,tmp);
	for (int t1=0; t1<3; t1++)
	{
		cin>>tmp>>cell[t1].x[0]>>cell[t1].x[1]>>cell[t1].x[2];
		getline(cin,tmp);
	}
	// find element symbol
	for (int t1=0; t1<nel; t1++)
	{
		while(tmp.find(flag_el) == string::npos)
			getline(cin,tmp);
		getline(cin,tmp);
		cin>>tmp>>el_l[t1];
		getline(cin,tmp);
	}
	// find coordinate in fractional unit
	while (tmp.find(flag_coord) == string::npos)
		getline(cin,tmp);
	for(int t1=0; t1<nat; t1++)
	{
		cin>>fcoord[t1].x[0]>>fcoord[t1].x[1]>>fcoord[t1].x[2];
//		coord[t1] = cell[0]*fcoord[t1].x[0] + cell[1]*fcoord[t1].x[1] + cell[2]*fcoord[t1].x[2];
	}
	// find force
	while(tmp.find(flag_force) == string::npos)
		getline(cin,tmp);
	for(int t1=0; t1<nat; t1++)
	{
		cin>>tmp>>force[t1].x[0]>>force[t1].x[1]>>force[t1].x[2];
	}
	// find stress
	while(tmp.find(flag_stress) == string::npos)
		getline(cin,tmp);
	cin>>tmp>>tmp>>stress[0].x[0]>>tmp>>tmp>>stress[2].x[1];
	cin>>tmp>>tmp>>stress[1].x[1]>>tmp>>tmp>>stress[2].x[0];
	cin>>tmp>>tmp>>stress[2].x[2]>>tmp>>tmp>>stress[1].x[0];
	stress[0].x[1] = stress[1].x[0];
	stress[0].x[2] = stress[2].x[0];
	stress[1].x[2] = stress[2].x[1];
	// calculate new cell parameter
	for (int t1=0; t1<3; t1++)
	{
		cell_new[t1] = cell[t1] + stress[t1] * lambda_t;
		cell_norm_new.x[t1] = cell_new[t1].norm();
		cell_unit_new[t1] = cell_new[t1] / cell_norm_new.x[t1];
	}
	// calculate cartesian coordinates
	for (int t1=0; t1<nat; t1++)
		coord[t1] = cell_new[0]*fcoord[t1].x[0] + cell_new[1]*fcoord[t1].x[1] + cell_new[2]*fcoord[t1].x[2];
	// calculate new position
	for(int t1=0; t1<nat; t1++)
		coord[t1] = coord[t1] + force[t1]*lambda_f;

	// print out cell and coord
	out_cell<<"acell"<<endl;
	cell_norm_new.print(out_cell);
//	out_cell<<setw(14)<<setprecision(9)<<fixed<<cell_norm_new.x[0]<<setw(14)<<setprecision(9)<<fixed<<cell_norm_new.x[1]<<setw(14)<<setprecision(9)<<fixed<<cell_norm_new.x[2];
	out_cell<<"rprim"<<endl;
	for (int t1=0; t1<3; t1++)
		cell_unit_new[t1].print(out_cell);
	//
	out_coord<<"xcart"<<endl;
	for(int t1=0; t1<nat; t1++)
		coord[t1].print(out_coord);

	out_cell.close();
	out_coord.close();
	return 0;
}
