/*

  RAKOTONIRINA TOKINANTENAINA MATHIEU RAZOKINY
  L3 MSIA
  
  
*/

#include <iostream>
#include <vector>
#include <fstream>
#include <cmath>
#include <iomanip>
using namespace std;
void ResolutionEliminationGauss(vector<vector<float>> A,vector<float> B, int dim);
void afficheSytemeEquation(vector<vector<float>> A,vector<float> B, int dim);
int main()
{
	cout << "give me a bottle of rum!" << endl;
	 vector<vector<float>> A;
	vector<float> temp;
	vector<float> B;
	int dim = 0;
	float number;
	ifstream file("data.txt");
	if(!file.is_open()){
		cout << "Le fichier est introuvable " << endl;
	}else{
		if(file >> number){
		    dim = number;
			cout << number << endl;
			for(int i(0); i<dim ; i++){
				
				for(int j(0); j<dim ; j++){
					file >> number;
					temp.push_back(number);
	
				}
				A.push_back(temp);
				temp.clear();    
			}
			for(int i(0); i<dim ; i++){
				file >> number;
				B.push_back(number);
				
			}
	    }
		file.close();
    }  
	 
	cout << "----------------------------------------VOTRE SYSTEME D'EQUATION------------------------" << endl << endl;
	afficheSytemeEquation(A,B,dim);
	cout << "----------------------------------------------------------------------------------------" << endl << endl;
	ResolutionEliminationGauss(A,B,dim);


	

	return 0;
}

void ResolutionEliminationGauss(vector<vector<float>> A,vector<float> B, int dim){
// epsilon
	float p = pow(10,-6);
// Déclaration des pointeurs et variable que nous allons utiliser dans la résolution	
	vector<float> *temp;
	vector<float> *pk;
	vector<float> *piv=nullptr;
	float *pb = nullptr;

// tp est une variable temporaire pour faire une permutation de même que pb
	float solution[dim] = {0};
	vector<float> tp({});
	float temporaire = 0;
	float max = 0;
		int pos = 0;

		for(int k(0); k<dim-1 ; k++){
			pb = &B[k];
			pos = k;
			pk = &A[k];
			piv = &A[k];
			max = (*(pk))[k];
	// Nous allons chercher le  plus grand pivot
			for(int l=k+1; l<dim ; l++){
				temp = &A[l];
				if(max<fabs((*temp)[k])){
					max = fabs((*temp)[k]);
					piv = temp;
	// on prendra la position du plus grand pivot pour pouvoir faire une permutation dans le second membre
					pos = l;
					
				}
				
			}
			if(fabs(max)<p){
				cout << "La matrice est numériquement non inversible" << endl;
				cout << "Nous avons une solution vide" << endl;
				return exit(1);
			}else{
				// permutation des pivot dans le second membre
			temporaire = *pb;
			*pb = B[pos];
			B[pos] = temporaire;
	// permutation des pivots dans la matrice
			tp = *pk;
			*pk = *piv;
			*piv = tp;
	//Maintenant nous allons faire l'elimination de gauss
			for(int i=k+1;i<dim;i++){
	// implément pk+dim-i pour eviter toutes erreurs de segmentation
				float alpha = (*(pk+i-k))[k]/(*(pk))[k];
				for(int j=k+1;j<dim;j++){
					(*(pk+i-k))[j] -=  alpha * (*(pk))[j];
					
				}
				(*(pk+i-k))[k] = 0;
				B[i] -= alpha * B[k]; 
			}
		}
	}

// Ainsi a la fin de l'elimination nous obtenons le système d'équation 
	cout << "--------------------------------------MATRICE TRIANGULARISER------------------------" << endl << endl;
	afficheSytemeEquation(A,B,dim);
	cout << "-------------------------------------------------------------------------------------" << endl << endl;

// Maintenant nous allons résoudre le systeme
	for(int i=dim-1;i>=0;i--){
		float s=0;
// implementant l'algorithme pour trouver la solution
		for(int j=i+1;j<dim;j++){
			s += A[i][j] * solution[j];
		}		
 		solution[i]= ((B[i]-s)/A[i][i]);
 	}
// Nous allons afficher la solution
	cout << "Nous avons comme resultat : (";
	for(int i=0; i<dim-1 ; i++)
		cout << solution[i] << ", ";
	cout << solution[dim-1] << ")";

}
// Methode pour affiche le systeme d'equation
void afficheSytemeEquation(vector<vector<float>> A,vector<float> B, int dim){
	
	for(int i=0; i<dim ; i++){
		for(int j=0 ; j < dim ; j++){
		
			cout << setw(10)<< A[i][j];
		}
		cout << setw(10)<<"= \t" << B[i] <<setw(10) << endl << endl;
	}
		
}