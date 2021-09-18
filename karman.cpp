#include <stdio.h>
#include <iostream>
#include "pico/stdlib.h"
#include <Eigen/Dense>

using Eigen::MatrixXd;
//大文字予測、小文字推定

int ekf(MatrixXd Xk, MatrixXd Pk, float ax, float ay, float az, float mx, float my, float mz, float p, float q, float r){
//z
	z(0,1) = ax;
	z(1,1) = ay;
	z(2,1) = az;
	z(3,1) = mx;
	z(4,1) = my;
	z(5,1) = mz;



//G
	G(0,0) = 0;
	G(0,1) = 0;
	G(0,2) = 0;
	G(1,0) = 0;
	G(1,1) = 0;
	G(1,2) = 0;
	G(2,0) = 0;
	G(2,1) = 0;
	G(2,2) = 0;
	G(3,0) = 0;
	G(3,1) = 0;
	G(3,2) = 0;
	G(4,0) = 1;
	G(4,1) = 0;
	G(4,2) = 0;
	G(5,0) = 0;
	G(5,1) = 1;
	G(5,2) = 0;
	G(6,0) = 0;
	G(6,1) = 0;
	G(6,2) = 1;

//Q
	Q(0,0) = shigumax;
	Q(0,1) = 0;
	Q(0,2) = 0;
	Q(1,0) = 0;
	Q(1,1) = shigumay;
	Q(1,2) = 0;
	Q(2,0) = 0;
	Q(2,1) = 0;
	Q(2,2) = shigumaz;


//R
	R(0,0) = shiguma_ax;
	R(1,0) = 0;
	R(2,0) = 0;
	R(3,0) = 0;
	R(4,0) = 0;
	R(5,0) = 0;
	R(0,1) = 0;
	R(1,1) = shiguma_ay;
	R(2,1) = 0;
	R(3,1) = 0;
	R(4,1) = 0;
	R(5,1) = 0;
	R(0,2) = 0;
	R(1,2) = 0;
	R(2,2) = shiguma_az;
	R(3,2) = 0;
	R(4,2) = 0;
	R(5,2) = 0;
	R(1,3) = 0;
	R(2,3) = 0;
	R(3,3) = shiguma_mx;
	R(4,3) = 0;
	R(5,3) = 0;
	R(1,4) = 0;
	R(2,4) = 0;
	R(3,4) = 0;
	R(4,4) = shiguma_my;
	R(5,4) = 0;
	R(1,5) = 0;
	R(2,5) = 0;
	R(3,5) = 0;
	R(4,5) = 0;
	R(5,5) = shiguma_mz;


	return 0;
}

int karuman5gyou(){

	//	for(short k =0; k<7; k++){
	Xk=Xk+Fk*h;
	Pk=Fk*pk*Fk.transpose() + G * Q * G.transpose();
	HkT= Hk.transpose();
	Kk=Pk * Hk * (Hk * Pk * HkT + R).inverse();
	Xk=Xk+Kk*(z-h);
	Pk=Pk-Kk*Hk*Pk;



	//}
}

int fk(qzero, qone, qtwo, qthree, derutap, derutaq, derutar){
//fk
	fk(0,0) = 0.5*(-p*qone-q*qtwo-r*qthree),
	fk(1,0) = 0.5*(p*qzero+r*qtwo-q*qthree),
	fk(2,0) = 0.5*(q*qzero-r*qone+p*qthree),
	fk(3,0) = 0.5*(r*qzero+q*qone-pqtwo),
	fk(4,0) = wx+(-bx*derutap),
	fk(5,0) = wy+(-by*derutaq),
	fk(6,0) = wz+(-bz*derutar);
	return 0;
}

int h(){

//h
	h(6,0) = 2*(qone*qthree+qzero*qtwo)*g,
	h(6,1) = 2*(qtwo*qthree+qzero*qone)*g,
	h(6,2) = (qzero*qzero-qone*qone-qtwo*qtwo+qthree*qthree)*g,
	h(6,3) = mx*(qzero*qzero+qone*qone-qtwo*qtwo-qthree*qthree)+mz*2*(qone*qthree+qzero*qtwo),
	h(6,4) = 2*(qone*qtwo-qzero*qthree)*mx+2*(qtwo*qthree+qzero*qone)*mz,
	h(6,5) = 2*(qone*qthree+qzero*qtwo)*mx+(qzero*qzero-qone*qone-qtwo*qtwo+qthree*qthree)*mz;
	return 0;
}

int F(){

//Fk
	Fk(0,0) = 0;
	Fk(1,0) = -0.5*p;
	Fk(2,0) = -0.5*q;
	Fk(3,0) = -0.5*r;
	Fk(4,0) = 0;
	Fk(5,0) = 0;
	Fk(6,0) = 0;
	Fk(0,1) = 0.5*p;
	Fk(1,1) = 0;
	Fk(2,1) = 0.5*r;
	Fk(3,1) = -0.5*q;
	Fk(4,1) = 0;
	Fk(5,1) = 0;
	Fk(6,1) = 0;
	Fk(0,2) = -0.5*q;
	Fk(1,2) = 0.5*r;
	Fk(2,2) = 0;
	Fk(3,2) = -0.5*p;
	Fk(4,2) = 0;
	Fk(5,2) = 0;
	Fk(6,2) = 0;
	Fk(0,3) = -0.5*r;
	Fk(1,3) = -0.5*q;
	Fk(2,3) = 0.5*p;
	Fk(3,3) = 0;
	Fk(4,3) = 0;
	Fk(5,3) = 0;
	Fk(6,3) = 0;
	Fk(1,4) = 0;
	Fk(2,4) = 0;
	Fk(3,4) = 0;
	Fk(4,4) = wx-bx;
	Fk(5,4) = wy;
	Fk(6,4) = wz;
	Fk(1,5) = 0;
	Fk(2,5) = 0;
	Fk(3,5) = 0;
	Fk(4,5) = wx;
	Fk(5,5) = wy-by;
	Fk(6,5) = wz;
	Fk(0,6) = 0;
	fk(1,6) = 0;
	fk(2,6) = 0;
	fk(3,6) = 0;
	fk(4,6) = wx;
	fk(5,6) = wy;
	fk(6,6) = wz-bz;
	return 0;
}

int Hk(){

//Hk
	Hk(0,0) = 2*qtwo;
	Hk(1,0) = 2*qone*g;
	Hk(2,0) = 2*qzero*g;
	Hk(3,0) = 2*mx*qzero+2*qtwo*mz;
	Hk(4,0) = -2*qthree*mx+2*qone*mz;
	Hk(5,0) = 2*qtwo*mx+2*qzero*mz;
	Hk(0,1) = 2*qthree*g; 
	Hk(1,1) = 2*qzero*g;
	Hk(2,1) = -2*qone*g;
	Hk(3,1) = mx*2*qone+2*mz*qthree;
	Hk(4,1) = 2*qtwo*mx+2*qzero*mz;
	Hk(5,1) = 2*qthree*mx-2*qone*mz;
	Hk(0,2) = 2*qzero*g;
	Hk(1,2) = 2*qthree*g;
	Hk(2,2) = -2*qtwo*g;
	Hk(3,2) = -2*qtwo*mx+2*mz*qzero;
	Hk(4,2) = 2*qone*mx-2*qthree*mz;
	Hk(5,2) = 2*qzero*mx-2*qtwo*mz;
	Hk(0,3) = 2*qone*g;
	Hk(1,3) = 2*qtwo*g;
	Hk(2,3) = 2*qthree*g;
	Hk(3,3) = -2*qthree*mx+2*mz*qone;
	Hk(4,3) = -2*qzero*mx+2*qthree*mz;
	Hk(5,3) = 2*qone*mx+2*qthree*mz;
	Hk(1,4) = 0;
	Hk(2,4) = 0;
	Hk(3,4) = 0;
	Hk(4,4) = 0;
	Hk(5,4) = 0;
	Hk(1,5) = 0;
	Hk(2,5) = 0;
	Hk(3,5) = 0;
	Hk(4,5) = 0;
	Hk(5,5) = 0;
	return 0;
}

int main(){

	float qone,qtwo,qthree,qzero,derutap,derutaq,derutar;
	float g=9.80665,shigumax,shigumay,shigumaz,shiguma_ax,shiguma_ay,shiguma_az,shiguma_mx,shiguma_my,shiguma_mz;
	float mx,mz,my,wx,wy,wz,p,q,r,ax,ay,az,bx,by,bz;
	double x[7] = {qzero,qone,qtwo,qthree,derutap,derutaq,derutar};
	MatrixXd Xk(7,1);
	MatrixXd Pk(7,7);
	MatrixXd Kk(7,6);
	MatrixXd xk(7,1);
	MatrixXd pk(6,1);
	MatrixXd HkT(7,6);
	MatrixXd z(6,1);
	MatrixXd h(6,1);
	MatrixXd fk(7,1);
	MatrixXd Fk(7,7);
	MatrixXd Hk(6,6);
	MatrixXd G(7,3);
	MatrixXd Q(3,3);
	MatrixXd R(6,6)	;

	// while(1){
	ekf();
	fk();
	h();
	Fk();
	Hk();



	//}
}
