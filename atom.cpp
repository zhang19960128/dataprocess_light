#include "atom.h"
#include <iostream>
#include <cmath>
#include <fstream>
#include <list>
//a=origin,b=end
double* distance(atom* a,atom* b,double* p){
	double* dist=new double[3];
	double temp;
	for(size_t i=0;i<3;i++){
		temp=b->position[i]-a->position[i];
		temp=(temp/p[i]-round(temp/p[i]))*p[i];
		dist[i]=temp;
	}
	return dist;
}
//compute the distance from a to b
double far(atom* a,atom* b,double* p){
	double* temp;
	double sum=0;
	temp=distance(a,b,p);
	for(size_t i=0;i<3;i++){
		sum=sum+temp[i]*temp[i];
	}
	return sqrt(sum);
}
int* changeindex(int index,int cell){
	int* re=new int[3];
	re[2]=floor(index/(cell*cell));
	index=index-re[2]*cell*cell;
	re[1]=floor(index/cell);
	re[0]=index-cell*re[1];
	return re;
}
void sort(double* input,int dim){
	double good;
	for(size_t i=0;i<dim-1;i++){
		for(size_t j=0;j<dim-i-1;j++){
			if(input[j]>input[j+1]){
			good=input[j];
			input[j]=input[j+1];
			input[j+1]=good;
			}
			else continue;
		}
	}
}
int changeback(int x,int y, int z,int cell){
	return (x+cell)%cell+(y+cell)%cell*cell+(z+cell)%cell*cell*cell;
}
int* neighbor_o_forB(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[6];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0],index_3D[1],index_3D[2]+1,cell);
	nei[2]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[3]=changeback(index_3D[0],index_3D[1]+1,index_3D[2],cell)+cell*cell*cell;
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	nei[5]=changeback(index_3D[0]+1,index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	return nei;
}
int* neighbor_o_forA(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[12];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0]-1,index_3D[1],index_3D[2],cell);
	nei[2]=changeback(index_3D[0],index_3D[1]-1,index_3D[2],cell);
	nei[3]=changeback(index_3D[0]-1,index_3D[1]-1,index_3D[2],cell);
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[5]=changeback(index_3D[0]-1,index_3D[1],index_3D[2],cell)+cell*cell*cell;
	nei[6]=changeback(index_3D[0]-1,index_3D[1],index_3D[2]-1,cell)+cell*cell*cell;
	nei[7]=changeback(index_3D[0],index_3D[1],index_3D[2]-1,cell)+cell*cell*cell;
	nei[8]=changeback(index_3D[0],index_3D[1],index_3D[2],cell)+2*cell*cell*cell;
	nei[9]=changeback(index_3D[0],index_3D[1]-1,index_3D[2],cell)+2*cell*cell*cell;
	nei[10]=changeback(index_3D[0],index_3D[1],index_3D[2]-1,cell)+2*cell*cell*cell;
	nei[11]=changeback(index_3D[0],index_3D[1]-1,index_3D[2]-1,cell)+2*cell*cell*cell;
	return nei;
}
int* neighbor_A_forB(int index,int cell){
	int* index_3D=changeindex(index,cell);
	int* nei=new int[8];
	nei[0]=changeback(index_3D[0],index_3D[1],index_3D[2],cell);
	nei[1]=changeback(index_3D[0]+1,index_3D[1],index_3D[2],cell);
	nei[2]=changeback(index_3D[0],index_3D[1]+1,index_3D[2],cell);
	nei[3]=changeback(index_3D[0]+1,index_3D[1]+1,index_3D[2],cell);
	nei[4]=changeback(index_3D[0],index_3D[1],index_3D[2]+1,cell);
	nei[5]=changeback(index_3D[0]+1,index_3D[1],index_3D[2]+1,cell);
	nei[6]=changeback(index_3D[0],index_3D[1]+1,index_3D[2]+1,cell);
	nei[7]=changeback(index_3D[0]+1,index_3D[1]+1,index_3D[2]+1,cell);
	return nei;
}
void sum_together(double* sum,double* add,int len){
	for(size_t i=0;i<len;i++){
		sum[i]=sum[i]+add[i];
	}
}
double average(std::list<double> &input){
	double sum=0.0;
	for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
		sum=sum+*a;
	}
	return sum/input.size();
}
double variance(std::list<double> &input){
	double sum=0.0;
	double ave=average(input);
	for(std::list<double>::iterator a=input.begin();a!=input.end();a++){
		sum=sum+(*a-ave)*(*a-ave);
	}
	return sum/input.size();
}
double* polar_average(atom *A,atom *B,atom *oxygen,double *p,int cell,std::vector<std::list<double> >& polar_x,std::vector<std::list<double> >& polar_y,std::vector<std::list<double> >& polar_z){
	std::list<double> px;
	std::list<double> py;
	std::list<double> pz;
	double volume=1.0;
	for(size_t k=0;k<3;k++){
		volume=volume*(p[k]/cell);
	}
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
			for(size_t k=0;k<3;k++){
				dist[k]=dist[k]*((oxygen+neighbor[j])->charge)/2.0;
			}
			sum_together(sum,dist,3);
		}
		neighbor=neighbor_A_forB(i,cell);
		for(size_t j=0;j<8;j++){
			dist=distance(B+i,A+neighbor[j],p);
			for(size_t k=0;k<3;k++){
				//now this guy turn into polar.
				dist[k]=dist[k]*((A+neighbor[j])->charge)/8.0;
			}
			sum_together(sum,dist,3);
		}
		polar_x[i].push_back(sum[0]/volume*16);
		polar_y[i].push_back(sum[1]/volume*16);
		polar_z[i].push_back(sum[2]/volume*16);
		px.push_back(sum[0]/volume*16);//16 is aim at converting the units from e to C
		py.push_back(sum[1]/volume*16);//16 is aim at converting the units from e to C
		pz.push_back(sum[2]/volume*16);//16 is aim at converting the units from e to C
	}
	double* pall=new double[3];
	pall[0]=average(px);
	pall[1]=average(py);
	pall[2]=average(pz);
	return pall;
}
double* displace_average_A(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
double* displace_average_Ca(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='c'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
		}
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
double* displace_average_Ba(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='b'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t t=0;t<12;t++){
			dist=distance(A+i,neighbor[t]+oxygen,p);
			sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/12.0);
		dy.push_back(sum[1]/12.0);
		dz.push_back(sum[2]/12.0);
		}
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
//compute the angle a---b----c
double tiltangle(atom* a,atom* b,atom* c,double* p){
	double ab=far(a,b,p);
	double bc=far(b,c,p);
	double ac=far(a,c,p);
	double theta;
	theta=(ab*ab+bc*bc-ac*ac)/2/ab/bc;
	if(theta>-1 && theta<1){
	   theta=(180-acos(theta)/(3.141592653)*180.0);
	}
	else{
		 theta=0;
	}
	return fabs(theta);
}
//giving the scalar average version of this code
double displace_average_Ca_scalar(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	double all=0;
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='c'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		all=0.0;
		for(size_t j=0;j<3;j++){
			all=all+sum[j]/12.0*sum[j]/12.0;
		}
		dall.push_back(sqrt(all));
		}
	}
	return average(dall);
}
double displace_average_Ba_scalar(atom* A,atom* oxygen,double *p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	double all=0;
	for(size_t i=0;i<cell*cell*cell;i++){
		if(A[i].type=='b'){
		neighbor=neighbor_o_forA(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<12;j++){
			dist=distance(A+i,neighbor[j]+oxygen,p);
			sum_together(sum,dist,3);
		}
		all=0.0;
		for(size_t j=0;j<3;j++){
			all=all+sum[j]/12.0*sum[j]/12.0;
		}
		dall.push_back(sqrt(all));
		}
	}
	return average(dall);
}
double displace_average_B_scalar(atom* B,atom* oxygen,double* p,int cell){
	std::list<double> dall;
	int* neighbor;
	double* dist;
	double all;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
		}
		all=0.0;
		for(size_t k=0;k<3;k++){
			all=all+sum[k]/6.0*sum[k]/6.0;
		}
		dall.push_back(sqrt(all));
	}
	return average(dall);
}
double* displace_average_B(atom* B,atom* oxygen,double* p,int cell){
	std::list<double> dx;
	std::list<double> dy;
	std::list<double> dz;
	int* neighbor;
	double* dist;
	double* sum=new double[3];
	for(size_t i=0;i<cell*cell*cell;i++){
		neighbor=neighbor_o_forB(i,cell);
		for(size_t k=0;k<3;k++){
			sum[k]=0.0;
		}
		for(size_t j=0;j<6;j++){
			dist=distance(B+i,neighbor[j]+oxygen,p);
		  sum_together(sum,dist,3);
		}
		dx.push_back(sum[0]/6.0);
		dy.push_back(sum[1]/6.0);
		dz.push_back(sum[2]/6.0);
	}
	double* dm=new double[3];
	dm[0]=average(dx);
	dm[1]=average(dy);
	dm[2]=average(dz);
	return dm;
}
