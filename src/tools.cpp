#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {

  VectorXd rmse(4);
  if(estimations.size()<0){
	rmse << -1,-1,-1,-1;
    return rmse;};
  if(estimations.size()!=ground_truth.size()){
	rmse << -1,-1,-1,-1;
    return rmse;};
  
  //accumulate squared residuals
  for(int i=0; i < estimations.size(); ++i){
	for(int j=0; j<4;++j){
	  rmse[j] += pow(estimations[i][j]-ground_truth[i][j],2);
    }
  }

  //calculate the mean
  for(int j=0; j<4;++j){
	  rmse[j] /= (1. * estimations.size());
    }

  //calculate the squared root
  for(int j=0; j<4;++j){
	  rmse[j] = pow(rmse[j],0.5);
    }
  
  //return the result
  return rmse;
}