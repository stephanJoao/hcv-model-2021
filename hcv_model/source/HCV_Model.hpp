#ifndef HCV_Model_H
#define HCV_Model_H

#include <iostream>
#include <string>

#define AGE 500
#define buffer 2

class HCV_Model
{
    // Variables
    private:

        double T;     // Target cells
        double V;     // Virus
        double A_alt; // ALT
        
        double I[buffer][AGE];  // Infected cells

        double Rp[buffer][AGE]; // RNA intracelular positivo
        double Rn[buffer][AGE]; // RNA intracelular negativo
        double Rt[buffer][AGE]; // RNA positivo traduzido


        int simCase;
        int days;
        int points;
        double deltaT;
        double deltaA;
        int iterPerDay;
        double tol;
        int vardelta;
        int varrho;

        double  N; // Initial values

        // Target cells
        double s;     // Target cells production rate
        double d;     // Target cells natural decay
        double beta;  // Infection rate (of target cells)

        // Virus (is coupled to the rest of the model with the integral in its equation)
        double c; // Virus elimination rate (by the imune system)
        
        // Infected cells
        double delta; // Infected cells decay (should depend on age)
        
        // ALT
        double s_alt; //ALT liberado no sangue por fatores externos 
        double c_alt; //ALT removido da circulação
        double alpha_alt; //ALT liberado pelas células infectadas
        
        // Intracelular positive-strand RNA
        double alpha; // Positive-strand RNA replication rate (should depend on age)
        double rho;   // Positive-strand RNA exportation rate (should depend on age)
        double mu_p; // Natural decay of positive and negative RNA
        
        // Intracelular negative-strand RNA
        double r; // Negative-strand RNA replication rate //TODO(tem uma dependência?)
        double Rmax; // Maximum number of negative RNA

        // Translated positive-strand RNA
        double sigma; // Replication complexes production rate (by translated positive RNA)
        double mu_t; // Natural decay of translated positive RNA


        double k; // coeficiente da fun��o exponencial de atraso na exporta��o de RNA positivo
        double tau; // Delay time for positive RNA exportation
        double n; //atraso de delta
        double theta; // taxa de disponibilidade para tradu��o

        // Therapy 
        double epsilon_s;     // Therapy effectiveness in reducing the exportation of positive RNA
        double epsilon_alpha; // Therapy effectiveness in reducing the replication of positive RNA
        double epsilon_r;     // Therapy effectiveness in reducing the replication of negative RNA

        double kappa_t; // Factor that increases positive RNA (ready for translation) degradation
        double kappa_c; // Factor that increases positive and negative RNA (translation complex) degradation
        
        int saveFiles;
        char* dir;
        FILE* dataInfected;
        FILE* dataVirus;
        FILE* dataTarget;
        FILE* dataRNA_Positivo;
        FILE* dataRNA_Negativo;
        FILE* dataRNA_Traduzido;
        FILE* datadelta;
        FILE* datarho;

        std::string Header();
        std::string Footer(long int t);
        int checkFile(FILE* theFile);
        
    // Methods
    private:

        void initialize();
        double calcIntegral(double vec1[][AGE], double vec2[][AGE], double vec3[][AGE]);
        double calcIntegral(double vec1[][AGE]);
        double calcIntegral(double a, double b, double vec1[][AGE], double vec2[][AGE], double delta, double rho, double deltaA);
        void update(double vec[][AGE]);

    public:
        
        HCV_Model();
        int solve();

};

#endif
