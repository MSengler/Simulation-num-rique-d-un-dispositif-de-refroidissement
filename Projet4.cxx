#include <iostream>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <string>


// Donnée
// Donnée
const double kappa = 164;
double h_c = 200;
double Phi_p = 125000;
double T_e = 20;
const double C_p = 940;
const double rho = 2700;

double L_x = 0.08;
double L_y = 0.004;
double L_z = 0.05;
double S = L_y*L_z;
double p = 2 * (L_y + L_z);

double tfinal = 300;
int N = 600;
int M = 100;
bool stationary = 0;
int Mx = 10; int My = 5; int Mz = 3;

void decomposition_LU(double* M, double* L, double* U, int taille){
    for (int i = 0; i < (taille+1)*(taille+1); i++){
        L[i] = 0;
        U[i] = 0;}

    L[0] = M[0];
    U[1*(taille+1)] = M[1*(taille+1)]/L[0];
    U[0] = 1;
    L[1] = M[1];
    for (int i = 1; i < taille; i++){
        U[i*(taille+1)+i] = 1;
        L[i*(taille+1)+i+1] = M[i*(taille+1)+i+1];
        L[i*(taille+1)+i] = M[i*(taille+1)+i] - L[i*(taille+1)+i+1] * U[i*(taille+1)+i-1];
        U[(i+1)*(taille+1)+i] = M[(i+1)*(taille+1)+i] / L[i*(taille+1)+i]; }
    U[taille*(taille+1)+taille] = 1;
    L[(taille-1)*(taille+1)+(taille-1)] = M[(taille-1)*(taille+1)+(taille-1)] - L[(taille-2)*(taille+1)+(taille-1)] * U[(taille-1)*(taille+1)+(taille-2)];
    U[(taille)*(taille+1)+(taille-1)] = M[(taille)*(taille+1)+(taille-1)] / L[(taille-1)*(taille+1)+(taille-1)] ;
    L[(taille)*(taille+1)+(taille)] = M[(taille)*(taille+1)+(taille)] - L[(taille-1)*(taille+1)+(taille)] * U[(taille)*(taille+1)+(taille-1)];
}   

// Resolution du systeme
double* resoudre (double* L, double* U, double* F, int taille){
    double* X = new double[taille+1];
    double* Y = new double[taille+1];
    
    Y[0] = F[0] / L[0];
    for (int i = 1; i < taille+1; i++){
        Y[i] = (F[i] - L[(i-1)*(taille+1)+i]*Y[i-1]) / L[i*(taille+1)+i];}
    X[taille] = Y[taille];
    for (int i = taille-1; i > -1; i--){
        X[i] = Y[i] - U[(i+1)*(taille+1)+i]*X[i+1];}
    delete [] Y;
    return X;
}
// Opération sur les matrices
void printMatrix(double* M, int ligne, int colonne) {
    for (int i=0; i<ligne; ++i){
        for (int j=0; j<colonne; ++j)
            std::cout <<std::setw(12)<< M[j*ligne+i] << " ";
        std::cout << std::endl;}}


double* matmul (double* A, double* B, int ligA, int colA, int colB){
    double* C = new double[ligA*colB];
    for (int i=0; i<ligA; i++){
        for (int j=0; j<colB; j++){
            C[j*ligA+i] = 0;
            for (int k=0; k<colA; k++){
                C[j*ligA+i] += A[k*ligA+i]*B[j*colA +k];}}}
    return C;}


// Modèle stationnaire

class Stationnaire {
    public:
    // Constructeur
    Stationnaire (int taille): Mtaille(taille) {
        M = new double[(taille+1)*(taille+1)];
        L = new double[(taille+1)*(taille+1)];
        U = new double[(taille+1)*(taille+1)];
        F = new double[taille+1];
        set_M();
        set_F();
        decomposition_LU(M,L,U,Mtaille);
        T = resoudre(L,U,F,Mtaille);
        }

    // Destructeur
    ~Stationnaire(){delete [] L; delete [] U; delete [] M; delete [] F;};

    // Assesseur, Mutatateur
    
    int taille() const {return Mtaille;} 
    void matrice() const {printMatrix(M,Mtaille+1,Mtaille+1);}
    void decompL() const {printMatrix(L,Mtaille+1,Mtaille+1);}
    void decompU() const {printMatrix(U,Mtaille+1,Mtaille+1);}
    void vecteurF() const {printMatrix(F,Mtaille+1,1);}
    void vecteurT() const {printMatrix(T,Mtaille+1,1);}
    void LU() const {printMatrix(matmul(L,U,Mtaille+1,Mtaille+1,Mtaille+1),Mtaille+1,Mtaille+1);}
    void MT(){ printMatrix(matmul(M,T,Mtaille+1,Mtaille+1,1),Mtaille+1,1);}


    // Initialisation de M, L, U et F
    void set_M(){
        const double h = L_x / Mtaille;
        for (int i = 0; i < (Mtaille+1)*(Mtaille+1); i++){
            M[i] = 0;}
        for (int i = 0; i < Mtaille; i++){
            M[i*(Mtaille+1)+i+1] = -kappa / pow(h,2); //a_i
            M[(i+1)*(Mtaille+1)+i] = -kappa / pow(h,2); //c_i
            M[i*(Mtaille+1)+i] = 2*kappa / pow(h,2) + h_c*p/S; } //b_i
        M[0] = kappa / h;  //b_0
        M[(Mtaille+1)*(Mtaille+1)-1] = kappa / h; //b_M
        M[(Mtaille-1)*(Mtaille+1)+Mtaille] = - kappa / h; //a_M
        M[Mtaille+1] = - kappa / h; //c_0
    }

    
    void set_F(){
        for (int i = 1; i < Mtaille; i++){
            F[i] =  h_c*p/S*T_e;}
        F[0] = Phi_p;
        F[Mtaille] = 0;
    }
    

    // Solution exacte
    static double Texact(double x){
        double a = h_c*p / (S*kappa);
        return T_e + Phi_p/kappa * cosh(std::sqrt(a)*L_x)/(sqrt(a)*sinh(sqrt(a)*L_x)) * cosh(sqrt(a)*(L_x-x))/cosh(sqrt(a)*L_x);
    }

    // Erreur
    double erreur() {
        double erreur = 0;
        for (int i = 0; i < Mtaille; i++){
            erreur += T[i]-Stationnaire::Texact(i*L_x/(Mtaille+1));}
        return erreur/Mtaille;}

    // format csv
    void save_csv(std::string format = "csv") {
        std::string nomfichier = "resultat/Stationnaire.";
        nomfichier += format;
        std :: ofstream ofile(nomfichier, std :: ios :: out);
        if (ofile){
            ofile << "x"<<" "<<"sol_app"<<" "<<"sol_exact"<<" \n";
             for (int i = 0; i < Mtaille+1; i++){
                double a = i*L_x/Mtaille;
                double b = T[i];
                double c = Stationnaire::Texact(i*L_x/Mtaille);
                ofile<<std :: fixed<< a<<" "<<b<<" "<<c<<" \n";}
            ofile.close ();}
    }

    // Opérations :
    friend void decomposition_LU(double* M, double* L, double* U, int taille);
    friend double* resoudre (double* L, double* U, double* F, int taille);
    friend void printMatrix(double* M, int ligne, int colonne);
    friend double* matmul (double* A, double* B, int ligA, int colA, int colB);
   
    protected:
    double* M;
    double* L;
    double* U;
    double* F;
    double* T;
    int Mtaille;
};



// Modèle instationnaire
class NonStationnaire {
        public:
    // Constructeur
    NonStationnaire (int taille, int temps = N, bool phiconst = 1): Mtaille(taille), Mtemps(temps),  phi_const(phiconst)  {
        M = new double[(taille+1)*(taille+1)];
        L = new double[(taille+1)*(taille+1)];
        U = new double[(taille+1)*(taille+1)];
        F = new double[taille+1];
        Ti = new double[taille+1];
        Ttps = new double [3*temps];

        set_Ti();
        set_M();
        set_F();
        decomposition_LU(M,L,U,Mtaille);
        iteration();
    }

    // Destructeur
    ~NonStationnaire() {delete [] L; delete [] U; delete [] M; delete [] F; delete[] Ti; delete [] Ttps;};

    // Assesseur 
    int taille() const {return Mtaille;} 
    void matrice() const {printMatrix(M,Mtaille+1,Mtaille+1);}
    void decompL() const {printMatrix(L,Mtaille+1,Mtaille+1);}
    void decompU() const {printMatrix(U,Mtaille+1,Mtaille+1);}
    void vecteurF() const {printMatrix(F,Mtaille+1,1);}
    void vecteurT() const {printMatrix(Ti,Mtaille+1,1);}
    void vecteurTtps() const {printMatrix(Ttps,Mtemps,3);}
    void LU() const {printMatrix(matmul(L,U,Mtaille+1,Mtaille+1,Mtaille+1),Mtaille+1,Mtaille+1);}

    double & operator[](int k) { return Ti[k]; }
    double const& operator()(int k) const { return Ti[k]; }

    // inilisation de M, F, Ti
    void set_Ti(){
        for (int i=0; i<Mtaille+1; i++){
            Ti[i] = T_e;}
    }

    void set_M(){
        const double h = L_x / Mtaille;
        const double deltaT = tfinal/N;
        for (int i = 0; i < (Mtaille+1)*(Mtaille+1); i++){
            M[i] = 0;}
        for (int i = 0; i < Mtaille; i++){
            M[i*(Mtaille+1)+i+1] = -kappa / pow(h,2); //a_i
            M[(i+1)*(Mtaille+1)+i] = -kappa / pow(h,2); //c_i
            M[i*(Mtaille+1)+i] = 2*kappa / pow(h,2) + h_c*p/S + rho*C_p/deltaT; } //b_i
        M[0] = kappa / h;   //b_0
        M[(Mtaille+1)*(Mtaille+1)-1] = kappa / h; //b_M
        M[(Mtaille-1)*(Mtaille+1)+Mtaille] = - kappa / h; //a_M
        M[Mtaille+1] = - kappa / h;  //c_0
    }

    
    
    void set_F(){
        const double deltaT = tfinal/N;
        for (int i = 1; i < Mtaille; i++){
            F[i] =  h_c*p/S*T_e + rho*C_p/deltaT*Ti[i];}
        F[0] = Phi_p;
        F[Mtaille] = 0;
    }

    void set_F(bool m){
        if (m){F[0] = Phi_p;}
        else {F[0] = 0;}
    }

    //calcule de la solution après N itérations
    void iteration(){
        bool k = 1;
        Ttps[0] = Ti[0]; Ttps[Mtemps] = Ti[(Mtaille+1)/2]; Ttps[2*Mtemps] = Ti[Mtaille];
        for (int n=1; n<Mtemps; n++){
            Ti = resoudre(L,U,F,Mtaille);
            Ttps[n] = Ti[0];
            Ttps[Mtemps + n] = Ti[Mtaille/2];
            Ttps[2*Mtemps + n] = Ti[Mtaille];
            set_F();
            if (!phi_const){
                set_F(k);
                if (n % 30 == 0){
                    if (k){k=0;}
                    else {k=1;}}}
            } 
    }

    // Calcul d'une itération
    void uneiteration(int n, bool flux){
        set_F();        
        if (!flux){
            if ((n /30)%2 == 0){
                set_F(1); }
            else {set_F(0);}}
         
        Ti = resoudre(L,U,F,Mtaille);
        }


        // Erreur
    double erreur() {
        double erreur = 0;
        for (int i = 0; i < Mtaille; i++){
            erreur += std::abs(Ti[i]-Stationnaire::Texact(i*L_x/(Mtaille+1)));}
        return erreur/Mtaille;}


        // format csv
    void save_csv(std::string format = "csv") {
        std::string nomfichier = "resultat/Nonstationnaire";
        nomfichier += std::to_string(Mtemps); 
        nomfichier += "."; nomfichier += format;
        std :: ofstream ofile(nomfichier,std :: ios :: out);
        if (ofile){
            ofile << "x"<<" "<<"sol_app"<<" \n";
             for (int i = 0; i < Mtaille+1; i++){
                double a = i*L_x/Mtaille;
                double b = Ti[i];
                ofile<<std :: fixed<< a<<" "<<b<<" \n";}
            ofile.close ();}
    }

        void save_csv_temps(std::string format = "csv") {
        std::string nomfichier = "resultat/Nonstationnairetemps";
        nomfichier += std::to_string(Mtemps); 
        nomfichier += "."; nomfichier += format;
        std :: ofstream ofile(nomfichier,std :: ios :: out);
        if (ofile){
            ofile << "t"<<" "<<"x_0"<<" "<<"x_M/2"<<" "<<"x_M"<< "\n";
             for (int i = 0; i < Mtemps; i++){
                int t = i;
                double x_0 = Ttps[i];
                double x_M2 = Ttps[Mtemps + i];
                double x_M = Ttps[2*Mtemps + i];
                ofile<<std :: fixed<< t<<" "<<x_0<<" "<<x_M2<<" "<<x_M<<" \n";}
            ofile.close ();}
    }

     // Opérations :
    friend void decomposition_LU(double* M, double* L, double* U, int taille);
    friend double* resoudre (double* L, double* U, double* F, int taille);
    friend void printMatrix(double* M, int ligne, int colonne);
    friend double* matmul (double* A, double* B, int ligA, int colA, int colB);


    protected:
    double* M;
    double* L;
    double* U;
    double* F;
    double* Ti;
    double* Ttps;
    int Mtaille;
    int Mtemps;
    bool phi_const;
};

class Points {
    public:
    Points(int Mx, int My, int Mz): Mx(Mx), My(My), Mz(Mz), temperature(new double*[Mx+1]){
        for (int i=0; i<Mx+1; i++ ){temperature[i] = new double[(My+1)*(Mz+1)];}};

    ~Points(){
        for (int i=0; i<Mx; i++){delete [] temperature[i];}
        delete [] temperature;}

    double & operator()(int i, int j, int k) { return temperature[i][k*(My+1) +j]; }
    double const& operator()(int i, int j, int k) const { return temperature[i][k*(My+1) +j]; }
    
    void vecteurTemp(int i) const {printMatrix(temperature[i],My+1,Mz+1);}

    int Mx; int My; int Mz; double** temperature;
};

class visualisation3D   {
	public:
	visualisation3D(int M, int temps, int Mx, int My, int Mz, bool phiconst = true): syst(new NonStationnaire(M,1,phiconst)), temps(temps), Mx(Mx), My(My), Mz(Mz), T_3D(new Points(Mx,My,Mz)), phiconst(phiconst){};
	
    ~visualisation3D(){delete T_3D;};
    void taille3D() const {std::cout << "Mx : " << Mx << " My : " << My << " Mz : " <<Mz;}

    void calculTijk(){
        int k = 0;
        for (int i=0; i < Mx+1; i++ ){
            k = (int) std::floor(i * ((*syst).taille() + 1) / (Mx + 1));
            double a = ((*syst)[k+1] - (*syst)[k]) / (L_x/((*syst).taille() + 1));
            double b = ((k+1)*L_x/((*syst).taille() + 1)*(*syst)[k] - k*L_x/((*syst).taille() + 1)*(*syst)[k+1]) /(L_x/((*syst).taille() + 1));
            for (int j=0; j<My+1; j++){for (k=0; k<Mz+1; k++){(*T_3D)(i,j,k) = a*i*L_x/(Mx+1) + b;}}
            }
        }

    void vecteurT3D(int i) const {printMatrix(T_3D->temperature[i],My+1,Mz+1);}

    void save_vtk(){
        for (int i=0; i<temps; i++){
            syst->uneiteration(i+1, phiconst);
            this->calculTijk();
            std::string nomfichier = "resultat/solution3D";
            nomfichier += std::to_string(i+1);
            nomfichier += ".vtk";
            std :: ofstream ofile(nomfichier, std :: ios :: out);
            if (ofile){
                ofile << "# vtk DataFile Version 2.0\n";
                ofile << "vtk output\n";
                ofile << "ASCII\n";
                ofile << "DATASET STRUCTURED_GRID\n";
                ofile << "DIMENSIONS " << Mx << " " << My << " " << Mz << "\n";
                ofile << "POINTS "<<Mx*My*Mz <<" float\n";
                for (int k=0; k<Mz; k++){for(int j=0; j<My; j++){for (int i=0; i<Mx; i++){
                    ofile << i << " " << j << " " << k << "\n";}}}
                //solution numérique
                ofile << "POINT_DATA "<< Mx*My*Mz << "\n";
                ofile << "FIELD FieldData 1\n";
                ofile << "sol1 1 " << Mx*My*Mz << " float\n";
                for (int k=0; k<Mz; k++){for(int j=0; j<My; j++){for (int i=0; i<Mx; i++){
                    ofile << (*T_3D)(i,j,k) <<"\n";}}}
            }
            ofile.close();
    }}

	protected:
    NonStationnaire* syst;
	Points* T_3D;
    int M;
    int temps;
    int Mx;
    int My;
    int Mz;
    bool phiconst;
	};



int main(int argc , char* argv []) {

    // Fichier de configuration
    std::string nomfichier;
    if (argc>1){
        nomfichier = argv[1];
        std::ifstream configFile(nomfichier);
        if (configFile){
            std::string monstring; int monint; int monint1; int monint2; double mondouble; bool monbool;
            configFile >> monstring >> monint >> monstring >> monint1 >> monstring >> monint2;
            L_x = monint; L_y = monint1; L_z = monint2;
            configFile >> monstring >> monint;
            M = monint;
            configFile >> monstring >> mondouble;
            Phi_p = mondouble;
            configFile >> monstring >> mondouble;
            h_c = mondouble;
            configFile >> monstring >> monint;
            T_e = monint;
            configFile >> monstring >> monbool;
            stationary = monbool;
            configFile >> monstring >> monint;
            tfinal = monint;
            configFile >> monstring >> monint;
            N = monint;
            configFile >> monstring >> monint >> monstring >> monint1 >> monstring >> monint2;
            Mx = monint; My = monint1; Mz = monint2;
            double S = L_y*L_z;
            double p = 2 * (L_y + L_z);
        }
        configFile.close();
    }
    
    /*
    Stationnaire syststat(M);
    syststat.save_csv();

    NonStationnaire systnonstat(M,N,true);
    systnonstat.save_csv();
    systnonstat.save_csv_temps();

    visualisation3D syst3D(M,N,Mx,My,Mz,true);
    syst3D.save_vtk()*/
    
    return 0;
}