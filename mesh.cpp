#include <iostream>
#include <fstream>

#include "mesh.hpp"

using namespace std;

class GammaFunctorTwoSided {
  double dF_min_des;
  double n;
  GammaFunctorTwoSided(double init_val,int nx) : dF_min_des(init_val), n(nx) {} 
  double operator()(double gam) {
    return gam/(n*(tanh(gam/2)-tanh(-gam/2)))*pow(1/cosh(-gam/2),2)-dF_min_des; 
  }
};

class GammaFunctorOneSided {
  double dF_min_des;
  double n;
  GammaFunctorOneSided(double init_val,int nx) : dF_min_des(init_val), n(nx) {}
  double operator()(double gam) {
    return tanh(gam/2.*(1./n - 1.))/tanh(gam/2.) + 1. - dF_min_des;
  }

};

void Mesh::genXmesh(string type) {

  //Generate Uniform Mesh Points
  if (type == "uniform") {
    //staggered mesh points
    for (int i=0; i<N+1; i++) {
      x_vect_sX[i] = double(i)/double(N);
      dxdz_sX[i] = 1./N;
      dzdx_sX[i] = 1.0/dxdz_sX[i]; 
    }
    //collocated mesh points
    for (int i=0; i<N; i++) {
      x_vect_cX[i] = (double(i)+0.5)/double(N);
      dxdz_cX[i] = 1./N;
      dzdx_cX[i] = 1.0/dxdz_cX[i];
    }
    dx_min = x_vect_sX[1] - x_vect_sX[0];
    dx_max = x_vect_sX[1] - x_vect_sX[0];
  }

  //Generate One Sided Xmesh
  if (type == "exponential") {
    //other mesh constants
    double x_end = (dx_max/beta)*log((dx_min/dx_max)*(exp(beta*N)-1)+1);
    double dx_max_des   = dx_max/x_end;
    double dx_min_des   = dx_min/x_end;

    //staggered mesh points
    double delta_cur;
    double dxMax = 0.0;
    for( int i = 0; i < N+1; i++){
      double x = double(i)/double(2.*N);
      x_vect_sX[i] = (dx_max_des/beta)*log((dx_min_des/dx_max_des)*(exp(beta*(x*2*N))-1)+1);
      dxdz_sX[i] = dx_min_des*exp(beta*(x*2*N))/((dx_min_des/dx_max_des)*(exp(beta*(x*2*N))-1)+1);
      dzdx_sX[i] = 1.0/dxdz_sX[i];
      if(i > 0){
            delta_cur = x_vect_sX[i] - x_vect_sX[i-1];
            if(delta_cur > dxMax) dxMax = delta_cur;
        }
   }
   //collocated mesh points
   for( int i = 0; i < N; i++){
        double x = (double(i)+0.5)/double(2.*N);
        x_vect_cX[i] = (dx_max_des/beta)*log((dx_min_des/dx_max_des)*(exp(beta*(x*2*N))-1)+1);
        dxdz_cX[i] = dx_min_des*exp(beta*(x*2*N))/((dx_min_des/dx_max_des)*(exp(beta*(x*2*N))-1)+1);
        dzdx_cX[i] = 1.0/dxdz_cX[i];
    }
    dx_min = x_vect_sX[1] - x_vect_sX[0];
   

  } 
}

void Mesh::genYmesh(string type) {

  //Generate Uniform Mesh Points
  if (type == "uniform") {
    //staggered mesh points
    for ( int j = 0; j < M+1; j++){
      y_vect_sY[j] = double(j)/double(M);
      dydn_sY[j] = 1./M;
      dndy_sY[j] = 1.0/dydn_sY[j];
    }
    //collocated mesh points
    for ( int j = 0; j < M; j++){ 
      y_vect_cY[j] = (double(j)+0.5)/double(M);
      dydn_cY[j] = 1./M;
      dndy_cY[j] = 1.0/dydn_cY[j];
    }
    dy_min = y_vect_sY[1] - y_vect_sY[0];
  }

  //Generate Two Sided Mesh


}

void Mesh::printXmesh(string type) {
  if (type == "x_vect_cX") {
    cout << "Printing x_vect_cX: " << endl;
    for (int i=0; i<N; i++) {
      cout << x_vect_cX[i] << " " << endl;
    }
    cout << endl;
  }
  if (type == "x_vect_sX") {
    cout << "Printing x_vect_sX: " << endl;
    for (int i=0; i<N+1; i++) {
      cout << x_vect_sX[i] << " " << endl;
    }
    cout << endl;
  }
}

void Mesh::printYmesh(string type) {
  if (type == "y_vect_cY") {
    cout << "Printing y_vect_cY: " << endl;
    for (int j=0; j<M; j++) {
      cout << y_vect_cY[j] << " " << endl;
    }
    cout << endl;
  }
  if (type == "y_vect_sY") {
    cout << "Printing y_vect_sY: " << endl;
    for (int j=0; j<M+1; j++) {
      cout << y_vect_sY[j] << " " << endl;
    }
    cout << endl;
  }
}
