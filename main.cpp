/*
*	Projet Info 1: Filtre 
*	Auteurs: ROMIGUIERES William & VILLON Alexandra
*	Date: 05/12/17 - 16h25
*
*
**/

#include <iostream.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h> 
#include <string.h> 
#include <complex>


#include "structures.h" //contient les structures

#define M_PI 3.14159265358979323846 //fixe la valeur de PI pour les calculs trigonométriques

#define F_MIN 1e6
#define F_MAX 1e9
#define F_STEP 1e6

using namespace std;


/* Prototype des fonctions */

void dataInput(TchebychevInput *dataFilter);
void fillMatrix(TchebychevInput *inputData, TchebychevCoef *coeffFiltre, int indice, int frequence, float tabChaineI[2][2]);
void calculProduitMatriceChaine(TchebychevInput *inputData, TchebychevCoef *coeffFiltre);
void printMatrix22(float ** mat);

/*Fonction principale */

void main(){
	
	/*Initialisation des structures contenant les paramètres du filtre */
	TchebychevInput dataFilter = {0.0,0.0,0.0,0};
	TchebychevCoef coefficients_filtre;
	
	/* Initialisation de la matrice finale*/
	int i;
	float** matTt;
	matTt = new float* [2];
		for(i = 0;i<2;i++){
			matTt[i] = new float[2];
		}


	/* Saisie des données du filtre par l'utilisateur*/
	
	dataInput(&dataFilter);

	/* Définition des tableaux  de tailles n contenant les paramètres ak, gk, et bk,  avec l'allocation dynamique*/
	
	long double* ak = (long double*)malloc((dataFilter.n)*sizeof(long double));
	float* gk = (float*)malloc((dataFilter.n)*sizeof(float));
	float* bk = (float*)malloc((dataFilter.n)*sizeof(float));
	
	/* Copie de tableaux dans notre structure contenant les coefficients du filtre*/
	coefficients_filtre.ak = ak;
	coefficients_filtre.gk = gk;
	coefficients_filtre.bk = bk;


	/*  Calcul des coefficients ak  */

	int k,p,m;

	for(k = 1;k<=dataFilter.n;k++){

		coefficients_filtre.ak[k-1] = sin(((2*k-1)*M_PI)/(2*dataFilter.n));
		//printf("\nak%d = %lf\n",k,coefficients_filtre.ak[k-1]); //debug
	}
	
	/* Calcul du coefficient Beta*/
	
	coefficients_filtre.beta = -1*log(tanh((dataFilter.Lar) / 17.37));

	
	/* Calcul du coefficient Gamma*/

	coefficients_filtre.gamma = sinh(coefficients_filtre.beta / (2*dataFilter.n));

	/* Calcul des coefficients bk */

	for(p = 1;p<=dataFilter.n;p++){

		coefficients_filtre.bk[p-1] =pow(coefficients_filtre.gamma,2)+pow(sin(p*M_PI/dataFilter.n),2);
		printf("\nbk%d = %lf\n",p,coefficients_filtre.bk[p-1]);
	}
	
	/* Calcul des coefficients k */
	coefficients_filtre.gk[0]=2*(coefficients_filtre.ak[0])/coefficients_filtre.gamma;
		
	for(m = 2;m<=dataFilter.n;m++){

		coefficients_filtre.gk[m-1] =(4*coefficients_filtre.ak[m-2]*coefficients_filtre.ak[m-1])/(coefficients_filtre.bk[m-2]*coefficients_filtre.gk[m-2]);   
		printf("\ngk%d = %lf\n",m,coefficients_filtre.gk[m-1]);
	}

	
	/*Calcul final:
	*
	* - Remplissage de chaque matrice
	* - Produit des matrices
	* - Calcul des paramètres S11 et S12
	* - Exportation dans un fichier
	*/

	calculProduitMatriceChaine(&dataFilter,&coefficients_filtre);


}


/*Cette fonction permet à l'utilisateur de saisir les données du filtres qui seront stockés dans notre structure*/
void dataInput(TchebychevInput *dataFilter){

	/* Variables temporaires */

	float fc_buff;
	double Lar_buff;
	float R1_buff;
	int n_buff;
	
	/* Saisie utilisateur */

	printf("Frequence de coupure fc?\n");
	scanf("%f",&fc_buff);

	printf("\nOrdre n du filtre ?\n");
	scanf("%d",&n_buff);

	printf("\nOndulation Lar en dB?\n");
	scanf("%lf",&Lar_buff);

	printf("\nValeur de R1?\n");
	scanf("%f",&R1_buff);
	
	/* Copie des données */
	dataFilter->fc = fc_buff;
	dataFilter->n = n_buff;
	dataFilter->Lar = Lar_buff;
	dataFilter->R1 = R1_buff;

	printf("%f",Lar_buff);

}

/* Cette fonction va retourner et remplir la ième matrice de la matrice chaîne en fonction de la fréquence donnée */

void fillMatrix(TchebychevInput *inputData, TchebychevCoef *coeffFiltre, int indice, int frequence, float tabChaineI[2][2]){	
	
	float ZliImPart;
	float YciImpart;


	if(indice%2 == 0){ //indice i pair
	
		ZliImPart = (coeffFiltre->gk[indice-1] * inputData->R1 * frequence) / inputData->fc;

		tabChaineI[0][0] = 1;
		tabChaineI[1][1] = 1;
		tabChaineI[1][0] = 0;
		tabChaineI[0][1] = ZliImPart;

		/*
		[1		Zl_indice]
		[0		        1]	
		
		*/

	
	}else{ //indice i impair
		YciImpart = (coeffFiltre->gk[indice-1] * frequence) / (inputData->fc* inputData->R1);
	
		tabChaineI[0][0] = 1;
		tabChaineI[1][1] = 1;
		tabChaineI[0][1] = 0;
		tabChaineI[1][0] = YciImpart;	

		/*
		[1		        0]
		[Yc_indice		1]	
		
		*/


	}

	//tabChaineI n'est pas retourné, mais modifié et donc utilisable par la suite
	

}

/* Cette fonction effectue le produit des n matrices ce qui donne la matrice chaîne pour une fréquence donnée.*/

void calculProduitMatriceChaine(TchebychevInput *inputData, TchebychevCoef *coeffFiltre){

	float MatChaineGlobal[2][2] = {{0,0},{0,0}};

	float ProdMatAdmi[2][2];
	float ProdMatImpe[2][2];

	float buffMatAdmi[2][2];
	float buffMatImpe[2][2];

	long double S11_module;
	long double S12_module;

	float A,B,C,D; //paramètres pour une lecture plus aisée

	FILE *fichierOutput = NULL; //fichier de sortie
	
	
	int i,j,k,f;

	/* Produit des N matrices */
	
	for(f = F_MIN;f<=F_MAX;f=f+F_STEP){

			fillMatrix(inputData, coeffFiltre, 1, f, ProdMatAdmi); //première des matrices admittances
			fillMatrix(inputData, coeffFiltre, 2, f,ProdMatImpe); //première des matrices impedances
			
				

			for(i = 3;i<=inputData->n;i = i+2){
				//produit des matrices admittances

				fillMatrix(inputData, coeffFiltre, i, f, buffMatAdmi);

				ProdMatAdmi[1][0] = buffMatAdmi[1][0] + ProdMatAdmi[1][0];

			}



			for(j = 4;j<=inputData->n;j = j+2){
				//produit des matrices impédances

				 fillMatrix(inputData, coeffFiltre, j, f, buffMatImpe);

				ProdMatImpe[0][1] = buffMatImpe[0][1] + ProdMatImpe[0][1];
			}


			
			//produit matrice admittance et impédance
			
			MatChaineGlobal[0][0] = 1; //[A]
			MatChaineGlobal[0][1] = ProdMatImpe[0][1]; //complexe - imaginaire pur!! [B]
			MatChaineGlobal[1][0] = ProdMatAdmi[1][0]; // complexe - imaginaire pur!! [C]
			MatChaineGlobal[1][1] = 1-ProdMatImpe[0][1]*ProdMatAdmi[1][0]; //réel pur [D]

			A = MatChaineGlobal[0][0];
			B = MatChaineGlobal[0][1];
			C = MatChaineGlobal[1][0];
			D = MatChaineGlobal[1][1];



			/*Calcul des paramètres S11 et S12 en fonction de la fréquence
			*
			*
			* Notez qu'ici nous travaillons directement avec le module pour pouvoir s'affranchir des complexes
			*/

			
			S11_module = sqrt(pow((A-D)*inputData->R1,2) + pow((B-C*pow(inputData->R1,2)),2)) / sqrt(pow((A+D)*inputData->R1,2) + pow((B-C*pow(inputData->R1,2)),2));
			S12_module = (2*(A*D+B*C)) / sqrt(pow((A+D)*inputData->R1,2) + pow((B-C*pow(inputData->R1,2)),2));

			S11_module = 20*log10(S11_module);
			S12_module = 20*log10(S12_module);

			/*Enregistrement dans un fichier en vue de Matlab*/
			
			/*création du fichier contenant S11*/
			fichierOutput=fopen("S11_out.m","a");
			fprintf(fichierOutput,"%f",S11_module);
			fprintf(fichierOutput,"%c",' ');
			fclose(fichierOutput);

			/*création du fichier contenant S12*/
			fichierOutput=fopen("S12_out.m","a");
			fprintf(fichierOutput,"%f",S12_module);
			fprintf(fichierOutput,"%c",' ');
			fclose(fichierOutput);

			/*création du fichier contenant les frequences*/
			fichierOutput=fopen("freq_out.m","a");
			fprintf(fichierOutput,"%d",f);
			fprintf(fichierOutput,"%c",' ');
			fclose(fichierOutput);

			
			printf("\n\nCalcul termine!\nA la frequence f = %d, S11 = %f et S12 = %f\n",f,S11_module, S12_module);
	}
		
}
