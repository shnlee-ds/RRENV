//Project: Reduced_Rank_Envelope_Regression
//liyixi
//2018.11.15
//Reduced_rank_envelope

#include <stdio.h>
#include "ReadMatrix615.h"

using namespace std;
using namespace boost;
using namespace Eigen;


MatrixXd reduceRankEnvelope(MatrixXd X, MatrixXd Y, int u, int d, MatrixXd Gamma){

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
	
	MatrixXd sqrtD = MatrixXd::Zero(p,p);
  	ArrayXXd D = svdx.singularValues().array().sqrt();
  	int i=0;
  	while(i<p){
  		sqrtD(i,i) = D(i,0);
  		i++;
  	}	
	MatrixXd Sx_sqrt = svdx.matrixU() * sqrtD * svdx.matrixU().inverse();
	
	sqrtD = MatrixXd::Zero(u,u);
  	D = svdy.singularValues().array().sqrt();
  	i=0;
  	while(i<u){
  		sqrtD(i,i) = D(i,0);
  		i++;
  	}
	MatrixXd Syg_sqrt = svdy.matrixU() * sqrtD * svdy.matrixU().inverse();
	
  	D = svdyp.singularValues().array().sqrt();
  	i=0;
  	while(i<u){
  		sqrtD(i,i) = D(i,0);
  		i++;
  	}
	MatrixXd Syg_sqrt2 = svdyp.matrixU() * sqrtD * svdyp.matrixU().inverse();

		
	MatrixXd Cancov = Syg_sqrt * Sxyg.transpose() * Sx_sqrt; // u*u * u*p * p*p = u*p
	
	//truncate canonical covariance
	JacobiSVD<MatrixXd> svd(Cancov,ComputeThinU | ComputeThinV);
	MatrixXd truncated = MatrixXd::Zero(u,u);
  	i=0;
  	while(i<d){
  		truncated(i,i) = svd.singularValues()(i,0);
  		i++;
  	}	
  	MatrixXd Cancov_truncated = svd.matrixU() * truncated * svd.matrixV().transpose(); // u*p, with rank = d
  	

  	//estimate beta
  	MatrixXd Beta = Gamma * Syg_sqrt2 * Cancov_truncated * Sx_sqrt; //r*u * u*u * u*p * p*p = r*p
  	return Beta;
}


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

//reduced rank envelope
MatrixXd Beta_hat = reduceRankEnvelope(X,Y,u,d,Gamma);


cout<< (Beta_hat - Beta).norm() << endl;
cout<< (Beta_R - Beta).norm() << endl;

ofstream Bfile;
   Bfile.open ("Beta_cpp.txt");
        Bfile << Beta_hat;
        Bfile.close();

}