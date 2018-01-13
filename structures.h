/* D�finition des structures */


typedef struct TchebychevInput{ // cette structure contiendra les param�tres d'entr�e du filtre saisis par l'utilisateur
		float fc;
		double Lar;
		float R1;
		int n;
} TchebychevInput;

typedef struct TchebychevCoef{ //cette structure contiendra les param�tres calcul�s du filtre
		long double* ak;
		float* gk;
		float* bk;
		long double beta;
		long double gamma;
		float* MatChaine;

} TchebychevCoef;

/*typedef struct complexMatrix {

	float complex m[2][2];
	
}*/