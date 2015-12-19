#include <stdio.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>

#include <iostream>
#include <string>
#include <vector>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/LU> 
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>

#define ToDeg 180/M_PI
#define ToRad M_PI/180

#define BODY 1
#define RLEG_J0 2
#define RLEG_J1 3
#define RLEG_J2 4
#define RLEG_J3 5
#define RLEG_J4 6
#define RLEG_J5 7

#define LLEG_J0 8
#define LLEG_J1 9
#define LLEG_J2 10
#define LLEG_J3 11
#define LLEG_J4 12
#define LLEG_J5 13

using namespace Eigen;
using namespace std;


class bipedLINK {
public:
  
  VectorXd p;			//Position in World Coordinates
  MatrixXd R;			//Attitude in World Coordinates
  VectorXd v;			//Linear Velocity in World Coordinates
  VectorXd w;			//Angular Velocity in World Coordinates
  
  VectorXd a;			//Joint Axis Vector(Relative to Parent)
  VectorXd b;			//Joint Relative Position(Relative to Parent)

  VectorXd c;			//Center of Mass(Link Local)
  MatrixXd I;			//Moment of Inertia(Linl Local)
  
  const char *name;			//Name of the link
  int m;					//Mass
  int sister;				//Sister ID
  int child;				//Child ID
  int mother;				//Mother ID
  
  double q;					//Joint Angle
  double dq;				//Joint Velocity
  double ddq;				//Joint Acceleration
  
  //VectorXd vertex(3);			//Shape(Vertex Information, Link Local)
  //VectorXd face(3);			//Shape(Vertex Information (Point Connections)
  
  bipedLINK(): p(3), R(3,3), v(3), w(3), a(3), b(3), c(3), I(3,3) {}
};

bipedLINK uLINK[14];
bipedLINK Rfoot;
bipedLINK Lfoot;

void PrintLinkName(int j) {	
	if(j!=0){
		cout << "j=" << j << " " << uLINK[j].name << "\n" << uLINK[j].p << "\n" << endl;
		PrintLinkName(uLINK[j].child);
		PrintLinkName(uLINK[j].sister);
		
		}

}

int TotalMass(int j) {
	if(j==0)return 0;
	else return uLINK[j].m + TotalMass(uLINK[j].sister) + TotalMass(uLINK[j].child);
	}
	

MatrixXd Hat(const VectorXd& x)
{
  MatrixXd xhat(3,3);
  xhat << 0,-x(2),x(1),
	      x(2),0,-x(0),
		  -x(1),x(0),0;
  
  return xhat;
}

MatrixXd Rodrigues(const VectorXd& x, double teta)
{
  MatrixXd E = MatrixXd::Identity(3,3);
  return E + Hat(x)*sin(teta) + Hat(x)*Hat(x)*(1-cos(teta));
}


void ForwardKinematics(int j){
	if(j==0)return;
	if(j!=1){
		int i = uLINK[j].mother;
		uLINK[j].p = uLINK[i].R * uLINK[j].b + uLINK[i].p;
		uLINK[j].R = uLINK[i].R * Rodrigues(uLINK[j].a, uLINK[j].q);
		}
	ForwardKinematics(uLINK[j].sister);
    ForwardKinematics(uLINK[j].child);		
	}
	
vector<int> FindRoute(int to){
	
  int j;
  int length = 1;

  vector<int> idx;
  idx.resize(1);
  idx[0]=to;
  
  j = uLINK[to].mother;

  while(j != 1){
  length++;
  idx.resize(length);
  idx[length-1]=j;
  j = uLINK[j].mother;
  }
  reverse(idx.begin(), idx.end()); 	
  return idx;
}



void SetJointAngles(const vector<int> &idx, const VectorXd& q){

// se pasa idx y el vector q
for (unsigned n = 0;n < idx.size();n++){
    int j = idx[n];
    uLINK[j].q = q[n];
}
ForwardKinematics(1);
}

//Jacobian matrix of current configration in World frame
MatrixXd CalcJacobian(const vector<int> &idx){
	Vector3d a;
	Vector3d target;
	int jsize = idx.size();
	
	target = uLINK[idx[jsize-1]].p;
	MatrixXd J = MatrixXd::Zero(6,jsize);
		
	for(int n=0; n<jsize; n++){
		int j = idx[n];
		int mom = uLINK[j].mother;
		a = uLINK[mom].R * uLINK[j].a; //joint axis in world frame		
		J.col(n) << a.cross(target - uLINK[j].p), a;
		}
		
	return J;
}

VectorXd rot2omega(const MatrixXd& R){
// T.Sugihara "Solvability-unconcerned Inverse Kinemacics based on 
// Levenberg-Marquardt method with Robust Damping," Humanoids 2009
	Vector3d w;
	Vector3d el;
	
	el << R(2,1)-R(1,2), R(0,2)-R(2,0), R(1,0)-R(0,1);
	double norm_el = el.norm();

	if (norm_el > 2.2204e-16){
		w = atan2(norm_el, R.trace()-1)/norm_el * el;
	}
	else if (R(0,0)>0 && R(1,1)>0 && R(2,2)>0)
		w << 0, 0, 0;
	else{
		w << R(0,0)+1, R(1,1)+1, R(2,2)+1;
		w = M_PI/2*w;
	}
	return w;
}

VectorXd CalcVWerr(const bipedLINK& Cref, const bipedLINK& Cnow){
	Vector3d perr;
	Vector3d werr;
	MatrixXd Rerr;
	VectorXd err(6);

	perr = Cref.p - Cnow.p;
	Rerr = Cnow.R.transpose() * Cref.R;
	werr = Cnow.R * rot2omega(Rerr);
	
	err << perr, werr;

	return err;
}

void MoveJoints(const vector<int> &idx, const VectorXd& dq){

	for(unsigned n=0; n<idx.size(); n++){
    int j = idx[n];
    uLINK[j].q = uLINK[j].q + dq(n);
	}
}

double InverseKinematics(int to, const bipedLINK& Target){
	
	VectorXd err;
	VectorXd dq;
	vector<int> idx;
	MatrixXd J;

	double lambda = 0.9;
	
	idx = FindRoute(to);
	ForwardKinematics(1);
	err = CalcVWerr(Target, uLINK[to]);
	
	for(int n = 1;n <= 10; n++){
	if (err.norm() < 1e-6) break;	
	J  = CalcJacobian(idx);
	dq = lambda * (J.colPivHouseholderQr().solve(err));
	cout << "dq:" << dq << endl; 
	MoveJoints(idx, dq);
	ForwardKinematics(1);
	err = CalcVWerr(Target, uLINK[to]);
	}

	
	return err.norm();
}

int main() {
	
  //VectorXd a(3);
  
  //a << 0, 1, 0;
  //double q = M_PI/3;
  
  //cout << "R: " <<  Rodrigues(a,q) << endl;
   
  uLINK[1].name = "BODY";
  uLINK[1].m = 4;
  uLINK[1].sister = 0;
  uLINK[1].child = 2;
  uLINK[1].mother = 0;
  uLINK[1].p << 0, 0, 0.6;
  uLINK[1].R = MatrixXd::Identity(3,3);
  uLINK[1].v << 0, 0, 0;
  uLINK[1].w << 0, 0, 0;
  uLINK[1].q = 0;
  uLINK[1].dq = 0;
  uLINK[1].ddq = 0;
  uLINK[1].a << 0, 0, 0;
  uLINK[1].b << 0, 0, 0;
  uLINK[1].c << 0, 0, 0;
  uLINK[1].I = MatrixXd::Identity(3,3);
  
  uLINK[2].name = "RLEG_J0";
  uLINK[2].m = 1;
  uLINK[2].sister = 8;
  uLINK[2].child = 3;
  uLINK[2].mother = 1;
  uLINK[2].p << 0, -0.1, 0.6;
  uLINK[2].R = MatrixXd::Identity(3,3);
  uLINK[2].v << 0, 0, 0;
  uLINK[2].w << 0, 0, 0;
  uLINK[2].q = 0;
  uLINK[2].dq = 0;
  uLINK[2].ddq = 0;
  uLINK[2].a << 0, 0, 1;
  uLINK[2].b << 0, -0.1, 0;
  uLINK[2].c << 0, 0, 0;
  uLINK[2].I = MatrixXd::Identity(3,3);
  
  uLINK[3].name = "RLEG_J1";
  uLINK[3].m = 1;
  uLINK[3].sister = 0;
  uLINK[3].child = 4;
  uLINK[3].mother = 2;
  uLINK[3].p << 0, -0.1, 0.6;
  uLINK[3].R = MatrixXd::Identity(3,3);
  uLINK[3].v << 0, 0, 0;
  uLINK[3].w << 0, 0, 0;
  uLINK[3].q = 0;
  uLINK[3].dq = 0;
  uLINK[3].ddq = 0;
  uLINK[3].a << 1, 0, 0;
  uLINK[3].b << 0, 0, 0;
  uLINK[3].c << 0, 0, 0;
  uLINK[3].I = MatrixXd::Identity(3,3);
  
  uLINK[4].name = "RLEG_J2";
  uLINK[4].m = 1;
  uLINK[4].sister = 0;
  uLINK[4].child = 5;
  uLINK[4].mother = 3;
  uLINK[4].p << 0, -0.1, 0.6;
  uLINK[4].R = MatrixXd::Identity(3,3);
  uLINK[4].v << 0, 0, 0;
  uLINK[4].w << 0, 0, 0;
  uLINK[4].q = 0;
  uLINK[4].dq = 0;
  uLINK[4].ddq = 0;
  uLINK[4].a << 0, 1, 0;
  uLINK[4].b << 0, 0, 0;
  uLINK[4].c << 0, 0, 0;
  uLINK[4].I = MatrixXd::Identity(3,3);
  
  uLINK[5].name = "RLEG_J3";
  uLINK[5].m = 1;
  uLINK[5].sister = 0;
  uLINK[5].child = 6;
  uLINK[5].mother = 4;
  uLINK[5].p << 0, -0.1, 0.3;
  uLINK[5].R = MatrixXd::Identity(3,3);
  uLINK[5].v << 0, 0, 0;
  uLINK[5].w << 0, 0, 0;
  uLINK[5].q = 0;
  uLINK[5].dq = 0;
  uLINK[5].ddq = 0;
  uLINK[5].a << 0, 1, 0;
  uLINK[5].b << 0, 0, -0.3;
  uLINK[5].c << 0, 0, 0;
  uLINK[5].I = MatrixXd::Identity(3,3);
  
  uLINK[6].name = "RLEG_J4";
  uLINK[6].m = 1;
  uLINK[6].sister = 0;
  uLINK[6].child = 7;
  uLINK[6].mother = 5;
  uLINK[6].p << 0, -0.1, 0;
  uLINK[6].R = MatrixXd::Identity(3,3);
  uLINK[6].v << 0, 0, 0;
  uLINK[6].w << 0, 0, 0;
  uLINK[6].q = 0;
  uLINK[6].dq = 0;
  uLINK[6].ddq = 0;
  uLINK[6].a << 0, 1, 0;
  uLINK[6].b << 0, 0, -0.3;
  uLINK[6].c << 0, 0, 0;
  uLINK[6].I = MatrixXd::Identity(3,3);
  
  uLINK[7].name = "RLEG_J5";
  uLINK[7].m = 1;
  uLINK[7].sister = 0;
  uLINK[7].child = 0;
  uLINK[7].mother = 6;
  uLINK[7].p << 0, -0.1, 0;
  uLINK[7].R = MatrixXd::Identity(3,3);
  uLINK[7].v << 0, 0, 0;
  uLINK[7].w << 0, 0, 0;
  uLINK[7].q = 0;
  uLINK[7].dq = 0;
  uLINK[7].ddq = 0;
  uLINK[7].a << 1, 0, 0;
  uLINK[7].b << 0, 0, 0;
  uLINK[7].c << 0, 0, 0;
  uLINK[7].I = MatrixXd::Identity(3,3);
  
  uLINK[8].name = "LLEG_J0";
  uLINK[8].m = 1;
  uLINK[8].sister = 0;
  uLINK[8].child = 9;
  uLINK[8].mother = 1;
  uLINK[8].p << 0, 0.1, 0.6;
  uLINK[8].R = MatrixXd::Identity(3,3);
  uLINK[8].v << 0, 0, 0;
  uLINK[8].w << 0, 0, 0;
  uLINK[8].q = 0;
  uLINK[8].dq = 0;
  uLINK[8].ddq = 0;
  uLINK[8].a << 0, 0, 1;
  uLINK[8].b << 0, 0.1, 0;
  uLINK[8].c << 0, 0, 0;
  uLINK[8].I = MatrixXd::Identity(3,3);
  
  uLINK[9].name = "LLEG_J1";
  uLINK[9].m = 1;
  uLINK[9].sister = 0;
  uLINK[9].child = 10;
  uLINK[9].mother = 8;
  uLINK[9].p << 0, 0.1, 0.6; 
  uLINK[9].R = MatrixXd::Identity(3,3);
  uLINK[9].v << 0, 0, 0;
  uLINK[9].w << 0, 0, 0;
  uLINK[9].q = 0;
  uLINK[9].dq = 0;
  uLINK[9].ddq = 0;
  uLINK[9].a << 1, 0, 0;
  uLINK[9].b << 0, 0, 0;
  uLINK[9].c << 0, 0, 0;
  uLINK[9].I = MatrixXd::Identity(3,3);
  
  uLINK[10].name = "LLEG_J2";
  uLINK[10].m = 1;
  uLINK[10].sister = 0;
  uLINK[10].child = 11;
  uLINK[10].mother = 9;
  uLINK[10].p << 0, 0.1, 0.6;
  uLINK[10].R = MatrixXd::Identity(3,3);
  uLINK[10].v << 0, 0, 0;
  uLINK[10].w << 0, 0, 0;
  uLINK[10].q = 0;
  uLINK[10].dq = 0;
  uLINK[10].ddq = 0;
  uLINK[10].a << 0, 1, 0;
  uLINK[10].b << 0, 0, 0;
  uLINK[10].c << 0, 0, 0;
  uLINK[10].I = MatrixXd::Identity(3,3);
  
  uLINK[11].name = "LLEG_J3";
  uLINK[11].m = 1;
  uLINK[11].sister = 0;
  uLINK[11].child = 12;
  uLINK[11].mother = 10;
  uLINK[11].p << 0, 0.1, 0.3;
  uLINK[11].R = MatrixXd::Identity(3,3);
  uLINK[11].v << 0, 0, 0;
  uLINK[11].w << 0, 0, 0;
  uLINK[11].q = 0;
  uLINK[11].dq = 0;
  uLINK[11].ddq = 0;
  uLINK[11].a << 0, 1, 0;
  uLINK[11].b << 0, 0, -0.3;
  uLINK[11].c << 0, 0, 0;
  uLINK[11].I = MatrixXd::Identity(3,3);
  
  uLINK[12].name = "LLEG_J4";
  uLINK[12].m = 1;
  uLINK[12].sister = 0;
  uLINK[12].child = 13;
  uLINK[12].mother = 11;
  uLINK[12].p << 0, 0.1, 0;
  uLINK[12].R = MatrixXd::Identity(3,3);
  uLINK[12].v << 0, 0, 0;
  uLINK[12].w << 0, 0, 0;
  uLINK[12].q = 0;
  uLINK[12].dq = 0;
  uLINK[12].ddq = 0;
  uLINK[12].a << 0, 1, 0;
  uLINK[12].b << 0, 0, -0.3;
  uLINK[12].c << 0, 0, 0;
  uLINK[12].I = MatrixXd::Identity(3,3);
  
  uLINK[13].name = "LLEG_J5";
  uLINK[13].m = 1;
  uLINK[13].sister = 0;
  uLINK[13].child = 0; 
  uLINK[13].mother = 12;
  uLINK[13].p << 0, 0.1, 0; 
  uLINK[13].R = MatrixXd::Identity(3,3);
  uLINK[13].v << 0, 0, 0;
  uLINK[13].w << 0, 0, 0;
  uLINK[13].q = 0;
  uLINK[13].dq = 0;
  uLINK[13].ddq = 0;
  uLINK[13].a << 1, 0, 0;
  uLINK[13].b << 0, 0, 0;
  uLINK[13].c << 0, 0, 0;
  uLINK[13].I = MatrixXd::Identity(3,3);
 
  cout << "name " << "sister " << "child " << "\n"; 
 
  ofstream outLinkParameters("LinkParameters.dat");
  ofstream outPosition("Position.dat");
  
  ForwardKinematics(1);
 
  for(int i=1; i < 14; i++) {
    cout << uLINK[i].name << " " << uLINK[i].mother << " " << uLINK[i].sister << 
	 " " << uLINK[i].child << "\n" << uLINK[i].p.transpose() << "\n" << uLINK[i].R <<"\n";
    outLinkParameters << uLINK[i].name << "," << uLINK[i].mother << "," << uLINK[i].sister << 
	 "," << uLINK[i].child << "," << uLINK[i].p.transpose() <<
	 "," << uLINK[i].v.transpose() << "," << uLINK[i].w.transpose() << "," << uLINK[i].q <<
	 "," << uLINK[i].dq << "," << uLINK[i].ddq << "," << uLINK[i].a.transpose() <<
	 "," << uLINK[i].b.transpose() << "," << uLINK[i].c.transpose() << "\n";
	 outPosition << uLINK[i].p(0) << " " << uLINK[i].p(1) << " " <<uLINK[i].p(2) << "\n";
  }
   
  //cout << "grandchild:" << uLINK[uLINK[uLINK[1].child].child].name << endl;
  
  cout << "Total mass:" << TotalMass(1) << "\n" << endl;
  
  uLINK[1].p << 0.0, 0.0, 0.65;
  ForwardKinematics(1);


  //set non singular posture
  uLINK[RLEG_J2].q = -5.0*ToRad;
  uLINK[RLEG_J3].q = 10.0*ToRad;
  uLINK[RLEG_J4].q = -5.0*ToRad;

  uLINK[LLEG_J2].q = -5.0*ToRad;
  uLINK[LLEG_J3].q = 10.0*ToRad;
  uLINK[LLEG_J4].q = -5.0*ToRad;

  //uLINK[BODY].p << 0.0, 0.0, 0.7;
  //uLINK[BODY].R = eye(3);

  uLINK[BODY].p << 0.0, 0.0, 0.5;
  //uLINK(BODY).R = eye(3);
  
  double rerr_norm;
  double lerr_norm;
    
  Rfoot.p << 0.0900, -0.1538, 0.0214;
  Rfoot.R << 0.7487,   -0.4117,    0.5195, 
             0.3269,    0.9111,    0.2509, 
			-0.5767,   -0.0180,    0.8168;

  rerr_norm = InverseKinematics(RLEG_J5, Rfoot);
  cout << "rerr_norm:" << rerr_norm << endl;
    

  Lfoot.p << -0.0044, 0.0519, 0.0321;
  Lfoot.R << 0.8820,   -0.4550,    0.1227,
			 0.4354,    0.8864,    0.1573,
			-0.1803,   -0.0853,    0.9799;
			
  lerr_norm = InverseKinematics(LLEG_J5, Lfoot);
  cout << "lerr_norm:" << lerr_norm << endl;


			

/*			
  vector<int> idx;
  idx = FindRoute(14);
  cout << "idx contains:";
  for (unsigned i=0;i<idx.size();i++)
  cout << ' ' << idx[i];
  cout << '\n';
  //SetJointAngles(idx,qf);
  //MatrixXd J(6,6);
  //J = CalcJacobian(idx);
  
  //cout << "\nJ: " << J << "\n";
   */

return 0;
}

