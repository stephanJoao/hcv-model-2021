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

        double I[buffer][AGE];  // Infected cells

        double Rp[buffer][AGE]; // RNA intracelular positivo
        double Rn[buffer][AGE]; // RNA intracelular negativo
        double Rt[buffer][AGE]; // RNA positivo traduzido

        double T;     // Target cells
        double V;     // Virus
        double A_alt; // ALT

        int simCase;
        int days;
        int points;
        double deltaT;
        double deltaA;
        int iterPerDay;
        double tol;
        int vardelta;
        int varrho;

        double  N; //valores iniciais

        double s; // taxa constante de producao de celulas alvo
        double d; // decaimento natural de celulas alvo
        double beta; // taxa de infecao
        double delta; // decaimento de celulas infectadas
        double rho; // taxa de exporta��o de RNA positivo
        double c; // taxa de eliminacaoo do virus pelo SI

        double alpha; // taxa de replica��o de RNA positivo

        double r; // taxa de replica��o de RNA negativo / complexo de replica��o (Rn)

        double k; // coeficiente da fun��o exponencial de atraso na exporta��o de RNA positivo
        double tau; // tempo de atraso para a exporta��o de RNA positivo
        double n; //atraso de delta
        double Rmax; // n�mero m�ximo de RNA negativo / complexo de replica��o (Rn)
        double sigma; // taxa de forma��o de complexos de replica��o
        double mu_t; // decaimento natural de Rt
        double theta; // taxa de disponibilidade para tradu��o
        double mu_c; // decaimento natural de Rc e Rn

        double epsilon_s;     // Efetividade da terapia em diminuir ou bloquear a exportação de RNA positivo
        double epsilon_alpha; // Efetividade da terapia em diminuir ou bloquear a replicação de RNA positivo
        double epsilon_r;     // Efetividade da terapia em diminuir ou bloquear a replica��o de RNA negativo

        double kappa_t; // fator para aumentar a degrada��o de RNA positivo dispon�vel para tradu��o
        double kappa_c; // fator para aumentar a degrada��o de RNA positivo e negativo no complexo de replica��o

        double s_alt; //ALT liberado no sangue por fatores externos 
        double c_alt; //ALT removido da circulação
        double alpha_alt; //ALT liberado pelas células infectadas

        
        int saveFiles;
        char *dir;
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
        double calcIntegral1(double vec1[][AGE]);
        double calcIntegral2(double a, double b, double vec1[][AGE], double vec2[][AGE], double delta, double rho, double deltaA);
        void update(double vec[][AGE]);

    public:
        
        HCV_Model();
        int solve();

};

#endif
