//Project: Reduced_Rank_Envelope_Regression
//liyixi
//2018.11.15
//Reduced_rank_envelope




#include <stdio.h>
#include "ReadMatrix615.h"
#include <RcppEigen.h>

// [[Rcpp:depends(RcppEigen)]]

using namespace std;
using namespace boost;
using namespace Eigen;

MatrixXd sqrtSingularvalues(MatrixXd D){
	MatrixXd sqrtD = MatrixXd::Zero(D.rows(),D.rows()); 
	ArrayXXd sqrtArrayD = D.array().sqrt();
  	int i=0;
  	while(i<D.rows()){
  		sqrtD(i,i) = sqrtArrayD(i,0);
  		i++;
  	}
  	return sqrtD;
}

class ReduceRankEnvelope{
public:
	MatrixXd X;
	MatrixXd Y;
	MatrixXd Gamma;
	MatrixXd Gamma0;
	MatrixXd Sx_sqrt;
	MatrixXd Syg_sqrt;
	MatrixXd Syg_sqrt2;
	MatrixXd Cancov_truncated;
	int u, d;
	MatrixXd Beta;
	MatrixXd FitValues;
	MatrixXd Residuals;
	MatrixXd ResidualVariance;
	ReduceRankEnvelope(MatrixXd dX, MatrixXd dY, int du, int dd, MatrixXd dGamma, MatrixXd dGamma0){
		X = dX;
		Y = dY;
		u = du;
		d = dd;
		Gamma = dGamma;
		Gamma0 = dGamma0;
	}

	void initialize(){
	cout << "---------Reduced-rank Envelope Regression--------" << endl;
	int p = X.rows(); //X: p*N
	int r = Y.rows(); //Y: r*N
	int N = Y.cols();

	//estimate envelope given u
	MatrixXd Y_gamma(u,N); 
	if(u==r){ 
		Y_gamma = Y;
	}
	else{
		Y_gamma = Gamma.transpose() * Y;
	}
	int r_gamma = Y_gamma.rows(); // Y_gamma: u * N

	//estimate the reduced rank envelope estimator
	//variance and covariance
	MatrixXd Sx = (X * X.transpose())/ double(N);  //p*p
	MatrixXd Syg = (Y_gamma*Y_gamma.transpose())/ double(N); //u*u 
	MatrixXd Sxyg = (X * Y_gamma.transpose()) /double(N); //p*u

	//canonical covariance
	
	JacobiSVD<MatrixXd> svdx(Sx.inverse(),ComputeThinU | ComputeThinV);
	JacobiSVD<MatrixXd> svdy(Syg.inverse(),ComputeThinU | ComputeThinV);
	JacobiSVD<MatrixXd> svdyp(Syg,ComputeThinU | ComputeThinV);
	
	Sx_sqrt = svdx.matrixU() * sqrtSingularvalues(svdx.singularValues()) * svdx.matrixU().inverse();
	Syg_sqrt = svdy.matrixU() * sqrtSingularvalues(svdy.singularValues()) * svdy.matrixU().inverse();
	Syg_sqrt2 = svdyp.matrixU() * sqrtSingularvalues(svdyp.singularValues())* svdyp.matrixU().inverse();
	MatrixXd Cancov = Syg_sqrt * Sxyg.transpose() * Sx_sqrt; // u*u * u*p * p*p = u*p
	
	//truncate canonical covariance
	JacobiSVD<MatrixXd> svd(Cancov,ComputeThinU | ComputeThinV);
	MatrixXd truncated = MatrixXd::Zero(u,u);
  	int i=0;
  	while(i<d){
  		truncated(i,i) = svd.singularValues()(i,0);
  		i++;
  	}
  		
  	Cancov_truncated = svd.matrixU() * truncated * svd.matrixV().transpose(); // u*p, with rank = d
  	}

  	//estimate beta
  	MatrixXd beta(){
  		cout << "------------Beta Estimation-----------" << endl;
  		Beta = Gamma * Syg_sqrt2 * Cancov_truncated * Sx_sqrt; //r*u * u*u * u*p * p*p = r*p
  		return Beta;
  	}

  	//estimate fit values and residuals
  	MatrixXd fitValues(){
  		cout << "-------------Fitted Values------------" << endl;
  		FitValues = Beta * X;
  		return FitValues;
  	}

  	MatrixXd residualValues(){
  		cout << "---------------Residuals--------------" << endl;
  		Residuals = Y - FitValues;
  		return Residuals;
  	}

  	MatrixXd residualVariance(){
  		cout << "---------Variance of Residuals---------" << endl;
  		MatrixXd Sy = Y * Y.transpose() / double(Y.cols()); //r*r
  		if(u== Y.rows()){
  			ResidualVariance = Gamma * Syg_sqrt2 * (MatrixXd::Identity(u,u) - Cancov_truncated * Cancov_truncated.transpose()) * Syg_sqrt2 * Gamma.transpose();
  		}
  		else{
  			ResidualVariance = Gamma * Syg_sqrt2 * (MatrixXd::Identity(u,u) - Cancov_truncated * Cancov_truncated.transpose()) * Syg_sqrt2 * Gamma.transpose() + Gamma0 * Gamma0.transpose() *  Sy * Gamma0 * Gamma0.transpose();
  		}
  		return ResidualVariance;
  	}
  	
};

extern "C"{
void R_reduced_rank_envelope(MatrixXd* R_X,MatrixXd* R_Y,MatrixXd* R_Gamma,MatrixXd* R_Gamma0){
 
	MatrixXd X = *R_X;
	MatrixXd Y = *R_Y;
	MatrixXd Gamma = *R_Gamma;
	MatrixXd Gamma0 = *R_Gamma0;

//parameters(TODO: use BIC to choose u,d)
	int u,d;
	u = 30;
	d = 10; 

//reduced rank envelope
ReduceRankEnvelope rre(X, Y, u, d, Gamma, Gamma0); //input(TODO: use BIC to choose u,d)
rre.initialize();
MatrixXd Beta_hat = rre.beta(); 
MatrixXd FitValues = rre.fitValues();
MatrixXd Residuals = rre.residualValues();
MatrixXd ResidualVariance = rre.residualVariance(); 

}
}