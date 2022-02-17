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

        double Rp[buffer][AGE]; // Positive-strand intracellular RNA 
        double Rn[buffer][AGE]; // Negative-strand intracellular RNA 
        double Rt[buffer][AGE]; // Positive-strand translated RNA 


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
        
        // Intracellular positive-strand RNA
        double alpha; // Positive-strand RNA replication rate (should depend on age)
        double rho;   // Positive-strand RNA exportation rate (is delayed according to tau variable)
        double mu_p;  // Natural decay of positive and negative RNA
        
        // Intracellular negative-strand RNA
        double r;    // Negative-strand RNA replication rate
        double Rmax; // Maximum number of negative RNA

        // Translated positive-strand RNA
        double theta; // Translation rate
        double sigma; // Replication complexes production rate (by translated positive RNA)
        double mu_t;  // Natural decay of translated positive RNA

        // Exportation delay
        double tau; // Delay time for positive RNA exportation
        double k; // Exponential function coefficient in exportation delay
        
        //TODO tentaram atrasar delta, verificar se ele varia no .cpp
        double n; //atraso de delta

        // Therapy 
        double epsilon_s;     // Therapy effectiveness in reducing the exportation of positive RNA
        double epsilon_alpha; // Therapy effectiveness in reducing the replication of positive RNA
        double epsilon_r;     // Therapy effectiveness in reducing the replication of negative RNA

        double kappa_t; // Factor that increases positive RNA (ready for translation) degradation
        double kappa_c; // Factor that increases positive and negative RNA (translation complex) degradation

        std::string parameters_file;
        
    // Methods
    private:

        void initialize();
        
        double calcIntegral(double vec1[][AGE], double vec2[][AGE], double vec3[][AGE]);
        double calcIntegral(double vec1[][AGE]);
        double calcIntegral(double a, double b, double vec1[][AGE], double vec2[][AGE], double delta, double rho, double deltaA);
        
        void update(double vec[][AGE]);

    public:
        
        HCV_Model(std::string parameters_file);
        void solve();

};

#endif
