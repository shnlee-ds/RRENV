//Project: Reduced_Rank_Envelope_Regression
//liyixi
//2018.11.15
//Reduced_rank_envelope

#include <stdio.h>
#include "ReadMatrix615.h"

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
	ReduceRankEnvelope(MatrixXd dX, MatrixXd dY, int du, int dd, MatrixXd dGamma, MatrixXd dGamma0){
		X = dX;
		Y = dY;
		u = du;
		d = dd;
		Gamma = dGamma;
		Gamma0 = dGamma0;
	}

	void initialize(){
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
	//MatrixXd Syg = ((Y_gamma.rowwise() - Y_gamma.colwise().mean()) * (Y_gamma.rowwise() - Y_gamma.colwise().mean()).transpose())/ double(N); //u*u 
	MatrixXd Syg = (Y_gamma*Y_gamma.transpose())/ double(N); //u*u 
	//MatrixXd Sxyg = (X * (Y_gamma.rowwise() - Y_gamma.colwise().mean()).transpose()) /double(N); //p*u
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
  		MatrixXd Beta = Gamma * Syg_sqrt2 * Cancov_truncated * Sx_sqrt; //r*u * u*u * u*p * p*p = r*p
  		return Beta;
  	}

  	
  	MatrixXd betaVariance(){
  		MatrixXd Sy = Y * Y.transpose() / double(Y.cols()); //r*r
  		MatrixXd BetaVariance = Gamma * Syg_sqrt2 * (MatrixXd::Identity(u,u) - Cancov_truncated * Cancov_truncated.transpose()) * Syg_sqrt2 * Gamma.transpose() + Gamma0 * Gamma0.transpose() *  Sy * Gamma0 * Gamma0.transpose();
  		return BetaVariance;
  	}
  	
};


int main(int argc, char *argv[]){
 
 //read in simulated data
 MatrixXd X = readFromFile<double>("X.txt"); //X from data_simulation.cpp
 MatrixXd Y = readFromFile<double>("Y.txt"); //Y from data_simulation.cpp
 MatrixXd Beta = readFromFile<double>("Betas.txt"); // Beta from data_simulation.cpp
 MatrixXd Beta_R = readFromFile<double>("Beta_R.csv");

//parameters
	int u,d;
	u = atoi(argv[1]);
	d = atoi(argv[2]); 

//read in estimated envelope from R
MatrixXd Gamma = readFromFile<double>("Gamma.csv"); //Gamma from red.rank.env.R
MatrixXd Gamma0 = readFromFile<double>("Gamma0.csv"); //Gamma0 from red.rank.env.R

//reduced rank envelope
ReduceRankEnvelope rre(X, Y, u, d, Gamma, Gamma0);
rre.initialize();
MatrixXd Beta_hat = rre.beta();
//cout << Beta_hat << endl;
MatrixXd Beta_variance_hat = rre.betaVariance();
//cout << Beta_variance_hat << endl;


cout<< (Beta_hat - Beta).norm() << endl;
cout<< (Beta_R - Beta).norm() << endl;

ofstream Beta_cpp;
   Beta_cpp.open ("Beta_cpp.txt");
       Beta_cpp << Beta_hat;
        Beta_cpp.close();

ofstream Beta_variance_cpp;
   Beta_variance_cpp.open ("Beta_variance_cpp.txt");
        Beta_variance_cpp << Beta_hat;
        Beta_variance_cpp.close();

}