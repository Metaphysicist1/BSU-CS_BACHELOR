#include <iostream>
#include<math.h>
#include<iomanip>
#define EPSILON 0.01
using namespace std;




// Initialization
void Init(double ** matrix,int n,double * vector,double * start_values){

	for(int i=0;i<n;i++){
		
		matrix[i]= new double[n];
		for(int j=0;j<n;j++){
			matrix[i][j]=0;
		}
		vector[i]=0;
		start_values[i]=0;
	}
}

void InputValues(double **matrix,int n,double * B){
	cout<<"Enter Matrix Values"<<endl;
	for(int i=0;i<n;i++){
		
		for(int j=0;j<n;j++){
			cout<<"A["<<i<<"]["<<j<<"]:";
			cin>>matrix[i][j];
			
		}
		cout<<"B["<<i<<"]:";
			cin>>B[i];
	}
}
// Output values
void PrintMatrix(double ** A,int n,double * B){
	for(int i=0;i<n;i++){
		for(int j=0;j<n;j++){
			
			
			cout<<A[i][j]<<'\t';
		}
		cout<<B[i]<<'\t'<<endl;
		
	}
}

//extra function for next function if elements on left are more than elements on right
double CheckRow(int row,int n,double ** a){
	double sum=0;	
	for(int j=0;j<n;j++){
			
			if(row!=j){
			sum+=a[row][j];
			}
		}
	return sum;
}

// checking if processing jacob method give result
bool isMethodApllied(int n,double ** A){
	
	for(int i=0;i<n;i++){
		
		for(int j=0;j<n;j++){
			
			if(A[i][i]>CheckRow(i,n,A)){
				break;
			}
			
			else{
				cout<<"Jacob method is not Applicable here";
				return false;
				}
		}		
	}
	cout<<"Jacob method is Applicable here";
	return true;
}

void Print_roots(double * A,int n){
	for(int i=0;i<n;i++){
		
		cout<<A[i]<<'\t';
	
	}
	cout<<endl;
}

// main function for computing jacobs method
void ComputeJacob(double **a,int maxIter,double * start_values,int n,double * B){
	double error=0;
	int flag=0;
	int i,j,k,iteration; //variables
	double sum=0; //summing elements at right side from equal sign (without B element)
	double new_values[n]; // for saving new computed values after each iteration
	for(iteration=0 ;iteration<maxIter;iteration++){ //main loop for iterations
		
		for(i=0;i<n;i++){ // loop for computing values at each iteration
			
			sum=0; // initializing sum after finish 	  of computation after each row
			
			for(j=0;j<n;j++){ // computing values at each row
				if(i!=j) //without diagonal values
				sum+=a[i][j]*start_values[j]; //summing
			}
			new_values[i]=(B[i]-sum)/a[i][i]; //new value through Jacob's method
		}
		
		for(k=0;k<n;k++){  //loop for finding the most precisied values 
			error = fabs(start_values[k]-new_values[k]);
			if(error<EPSILON){ //comparing with epsilon value which is defined at start
				continue; // if values are less in that case we skipping step
							// if all values are less than epsilon we stopping computing and getting values and iteration number
			}
			
			else{ // if not we are continuing computing at next iteration
				flag=1; // temporary varable
				break;
				
			}
			
		}
		for(k=0;k<n;k++){
			start_values[k]=new_values[k];
		}	
		cout<<"Iteration("<<iteration<<")"<<'\t'<<"x1="<<new_values[0]<<'\t';
		cout<<"x2="<<new_values[1]<<'\t'<<"x3="<<new_values[2]<<'\t'<<"Error: "<<error;
		if(flag==0){
		cout<<'\t'<<"less than EPSILON";
		}
		cout<<endl;
		flag =0; //init
	}
}

//void ComputeSeidel(double **a,int maxIter,double * start_values,int n,double * B){
//	
//	int flag=0;
//	int i,j,k,iteration; //variables
//	double sum=0; //summing elements at right side from equal sign (without B element)
//	double new_values[n]; // for saving new computed values after each iteration
//	
//	for(iteration=1;iteration<maxIter;iteration++){ //main loop for iterations
//		new_values=start_values;
//		for(i=0;i<n;i++){ // loop for computing values at each iteration
//			
//			sum=0; // initializing sum after finish of computation after each row
//			
//			for(j=0;j<n;j++){ // computing values at each row
//				if(i!=j) //without diagonal values
//				sum+=a[i][j]*new_values[j]; //summing
//				new_values[j]=(B[i]-sum)/a[i][i]; //new value through Jacob's method
//			}
//			
//		}
//		
//		for(k=0;k<n;k++){  //loop for finding the most precisied values 
//			
//			if(fabs(start_values[k]-new_values[k])<EPSILON){ //comparing with epsilon value which is defined at start
//				continue; // if values are less in that case we skipping step
//							// if all values are less than epsilon we stopping computing and getting values and iteration number
//			}
//			
//			else{ // if not we are continuing computing at next iteration
//				flag=1; // temporary varable
//				break;
//				
//			}
//			
//		}
//		if(flag==0){
//			cout<<iteration<<new_values;	//results
//			return ;
//		}
//		flag =0; //init
//		
//		for(k=0;k<n;k++){
//			start_values[k]=new_values[k];
//		}
//		printf("Iteration: ",iteration,'\t',new_values);
//	}
//		
//	
//}

int main(void){
	int n=0;
	int maxIter=11;
	
	cout<<"------------------- START ------------------------"<<endl;
	cout<<"Enter square Matrix size: ";
	cin>>n;
	cout<<"Matrix size: "<<n<<endl;
	
	double ** A;
	double * B;
	double * start_values;
	
	A = new double *[n];
	B = new double [n];
	start_values = new double [n];
	
	
	cout<<"------------------- Print and Init Matrix ------------------------"<<endl;
	
	Init(A,n,B,start_values);
	PrintMatrix(A,n,B);
	
	cout<<"------------------- Input Values of Matrix and Print them ------------------------"<<endl;
	
	InputValues(A,n,B);
	PrintMatrix(A,n,B);

	cout<<"------------------- Check applicability of matrix ------------------------"<<endl;
	
//	if(isMethodApllied(n,A)){
	cout<<"Before computation"<<endl;
	Print_roots(start_values,n);
	ComputeJacob(A,maxIter,start_values,n,B);
	cout<<"After computation"<<endl;
	Print_roots(start_values,n);
//	}
	return 0;
}
