#include <iostream>
#include <fstream>
#include <math.h>

#include "mesh.hpp"
#include "boost/math/tools/roots.hpp"

using namespace std;

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
  if (type == "twoSided") {
    //Calculate gamma
    GammaFunctorTwoSided f(dy_min, M);
    typedef std::pair<double, double> Result;
    boost::uintmax_t max_iter=5000;
    boost::math::tools::eps_tolerance<double> tol(50);
    Result gamma_bounds = boost::math::tools::toms748_solve(f, 0.0000001, 100000.0, tol, max_iter);
    double gamma = (gamma_bounds.first + gamma_bounds.second)/2.0;

    //Other Mesh Constants
    double a = double(1.0)/(tanh(gamma/double(2.0))-tanh(double(-1.0)*gamma/double(2.0)));
    double b = double(0.5);

    // Generate Staggered Mesh Points
    double delta_cur;
    double dyMax = 0.0;
    for( int j = 0; j < M+1; j++){
      double y = double(j)/double(M);
      y_vect_sY[j] = a*tanh(gamma*(y-.5)) + b;
      dydn_sY[j] = a*pow(1/cosh(gamma*(y-0.5)),2)*gamma/M;
      dndy_sY[j] = 1.0/dydn_sY[j];
      if(j > 0){
        delta_cur = y_vect_sY[j] - y_vect_sY[j-1];
        if(delta_cur > dyMax) dyMax = delta_cur;
      }
  }
  //Generate Collocated Mesh Points
  for( int j = 0; j < M; j++){
        double y = (double(j)+0.5)/double(M);
        y_vect_cY[j] = a*tanh(gamma*(y-.5))+b;
        dydn_cY[j] = a*pow(1/cosh(gamma*(y-0.5)),2)*gamma/M;
        dndy_cY[j] = 1.0/dydn_cY[j];
    }

    dy_min = y_vect_sY[1] - y_vect_sY[0];   
  }
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
  if (type == "dxdz_cX") {
    cout << "Printing dxdz_cX: " << endl;
    for (int i=0; i<N; i++) {
      cout << dxdz_cX[i] << " " << endl;
    }
  }
  if (type == "dxdz_sX") {
    cout << "Printing dxdz_sX: " << endl;
    for (int i=0; i<N+1; i++) {
      cout << dxdz_sX[i] << " " << endl;
    }
  }
  if (type == "dzdx_cX") {
    cout << "Printing dzdx_cX: " << endl;
    for (int i=0; i<N; i++) {
      cout << dzdx_cX[i] << " " << endl;
    }
  }
  if (type == "dzdx_sX") {
    cout << "Printing dzdx_sX: " << endl;
    for (int i=0; i<N+1; i++) {
      cout << dzdx_sX[i] << " " << endl;
    }
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
  if (type == "dydn_cY") {
    cout << "Printing dydn_cY: " << endl;
    for (int j=0; j<M; j++) {
      cout << dydn_cY[j] << " " << endl;
    }
  }
  if (type == "dydn_sY") {
    cout << "Printing dydn_sY: " << endl;
    for (int j=0; j<M+1; j++) {
      cout << dydn_sY[j] << " " << endl;
    }
  }
}
