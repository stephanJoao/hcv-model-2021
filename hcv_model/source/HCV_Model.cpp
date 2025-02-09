#include <fstream>

#include "HCV_Model.hpp"
#include "GaussLegendreQuadrature.hpp"

/**
* Constructor
*/
HCV_Model::HCV_Model(std::string input, std::string output)
{
    this->input = input;
    this->output = output;
}

/**
* Set conditions and parameter values
*/
void HCV_Model::initialize()
{
    std::cout << "Reading from \"" << input << "\"" << std::endl;
    
    std::ifstream parameters_file;
    parameters_file.open(input.c_str());
    
    std::string aux_string;

    getline(parameters_file, aux_string, ',');
    double V0 = atof(aux_string.c_str());

    getline(parameters_file, aux_string, ',');
    epsilon_r = atof(aux_string.c_str());
    
    getline(parameters_file, aux_string, ',');
    epsilon_alpha = atof(aux_string.c_str());
    
    getline(parameters_file, aux_string, ',');
    epsilon_s = atof(aux_string.c_str());

    getline(parameters_file, aux_string, ',');
    alpha = atof(aux_string.c_str());

    getline(parameters_file, aux_string, ',');
    r = atof(aux_string.c_str());
        
    getline(parameters_file, aux_string, ',');
    delta = atof(aux_string.c_str());
    
    getline(parameters_file, aux_string, ',');
    mu_p = atof(aux_string.c_str());

    getline(parameters_file, aux_string, ',');
    rho = atof(aux_string.c_str());
    
    getline(parameters_file, aux_string, ',');
    theta = atof(aux_string.c_str());
    
    getline(parameters_file, aux_string, ',');
    sigma = atof(aux_string.c_str());

    getline(parameters_file, aux_string, ',');
    c = atof(aux_string.c_str());

    parameters_file.close();

    s_alt= 1.0;
    c_alt = 1.0;
    alpha_alt = 2.0;

    /**
    * Number of days simulated
    */
    days = 72;
    /**
    * Number of points saved
    */
    points = 10;
    /**
    * Simulation
    */
    deltaT = 0.01; //pow(10,-4);
    /**
    * Age
    */
    deltaA = 0.1;
    /**
    * number of iterations per day
    */
    iterPerDay = 100;

    tol = pow(10,-4);

    /**
     * Vary delta 1, fix delta 0
     */
    vardelta = 0;
    /**
     * Vary rho 1, fix delta 0
     */
    varrho = 1;

    /**
    * Baseline parameters
    */
    d     = 0.010;
    s     = 130000;
    // delta = 0.58;
    beta  = 5*pow(10,-8);
    // c     = 22.30;
    // rho   = 8.180;
    // alpha = 30.0;
    Rmax  = 50.0;
    // r     = 1.49; 
    tau   = 0.50;
    n     = 1.00;
    k     = 0.80;
    mu_t  = 0.89; 
    // mu_p  = 2.55; 
    // sigma = 1.30;
    // theta = 1.20;
    /**
    * therapy parameters
    */
    // epsilon_alpha = 0.928;
    // epsilon_r     = 0.47; 
    // epsilon_s     = 0.998;
    kappa_t       = 1.00;
    kappa_c       = 1.00;
    
    if ((sigma + rho + mu_p - (sigma * theta) / (theta + rho + mu_t) <= 0) || 
            (alpha*r - (sigma + rho + mu_p - (sigma * theta) / (theta + rho + mu_t)) * mu_p <= 0))
    {
        sigma = 1.30;
        theta = 1.20;
        epsilon_alpha = 0.928;
        epsilon_r =0.47;
        // cout << "caiu nas condiçoes de estabilidade!!!!" << endl;
    } else {
        // cout << "Não entrou nas condiçoes de estabilidade!!!!" << endl;
    }
    
    /**
    * Initial Conditions
    */
    for(int j = 0; j < buffer; j++){
        for(int i = 0; i < AGE; i++){
            I[j][i] = 0.0;
            Rt[j][i] = 0.0;
            Rp[j][i] = 0.0;
            Rn[j][i] = 0.0;
        }
    }
    Rt[0][0] = 1.00;

    double rho1;
    if(varrho) 
        rho1 = 0.00;
    else 
        rho1 = rho;

    double soma = 0;
    soma += Rp[0][0] + Rt[0][0];

    for(int a = 1; a < AGE; a++){
        if(varrho){
            if((deltaA*(double) a)<tau){
                rho1=0.0;
            }else{
                rho1=(1-exp(-k*((deltaA*(double)a)-tau)))*rho;
            }
        }else rho1 = rho;
        Rn[0][a] = (r*Rp[0][a-1] - r*Rp[0][a-1]*(Rn[0][a-1]/Rmax) -
                mu_p*Rn[0][a-1])*deltaA + Rn[0][a-1];

        Rp[0][a] = (alpha*Rn[0][a-1] + sigma*Rt[0][a-1] -
                theta*Rp[0][a-1]- rho1*Rp[0][a-1] -
                mu_p*Rp[0][a-1])*deltaA + Rp[0][a-1];

        Rt[0][a] = (theta*Rp[0][a-1] - sigma*Rt[0][a-1] -
                rho1*Rt[0][a-1] - mu_t*Rt[0][a-1])*deltaA + Rt[0][a-1];

        //if (a<250) soma+=(Rp[0][a] + Rt[0][a]);

        //printf("a = %d Rt = %.4lf Rp = %.4lf Rn = %.4lf \n", a, Rt[0][a], Rp[0][a], Rn[0][a]);
    }
    
    N = calcIntegral(0.0, AGE, Rp, Rt, delta, rho, deltaA); 
    T = c / (beta * N);                                       
    V = V0; //fixar baseado na outra DE
    I[0][0] = beta * T * V;
    A_alt = (s+I[0][0])/c_alt; //MATHEUS linha 233 do docx

    double delta1;
    if (vardelta) delta1 = 0.01;
    else delta1 = delta;
    for(int a=1; a<AGE;a++){
        //delta variavel
        if(vardelta) delta1=(1-exp(-k*(a*deltaA)))*delta;
        else delta1 = delta;
        I[0][a] = (beta*V*T*exp(-delta1*deltaA*(double)a))*deltaA + I[0][a-1];
       // printf("a = %d I0 = %lf \n", a, I[0][a]);
    }
}

/**
 * Updates the current results to position 1
 */
void HCV_Model::update(double vec[][AGE])
{
    for(int a = 0; a < AGE; a++) {
        vec[0][a] = vec[1][a];
	}
}

double HCV_Model::calcIntegral(double vec1[][AGE], double vec2[][AGE],double vec3[][AGE])
{
    int a;
    double sum = 0.0;
    for(a = 0; a < AGE; a++){
        sum += (vec1[0][a] + vec2[0][a]) * vec3[0][a];
    }
    return sum / (2.0 * AGE);
}

double HCV_Model::calcIntegral(double vec1[][AGE])
{ 
    //MATHEUS integral do ALT
    int a;
    double sum = 0.0;
    for(a = 0; a < AGE; a++){
        sum += vec1[0][a];
    }
    return sum / (2.0*AGE);
}

double HCV_Model::calcIntegral(double a, double b, double vec1[][AGE], double vec2[][AGE], double delta, double rho, double deltaA)
{
    Rosetta::GaussLegendreQuadrature<5> gl5;
    std::setprecision(6);
    //gl5.print_roots_and_weights(std::cout);
    return gl5.integrate(a, b, vec1, vec2, Rosetta::RosettaExp, delta, rho, deltaA);
}

/**
 * Solve model equations
 */
void HCV_Model::solve()
{
    long int t = 0;
    long int a = 0;

    // Set initial conditions
    initialize();

    // Opens output file
    std::ofstream output_file;
    output_file.open(output.c_str());
    output_file << "time,viral_load" << std::endl;

    std::cout << "Writing in \"" << output << "\"" << std::endl;

    /**
    * Begin time loop
    */
    do {

        //TODO Dá warning de "not used" / pode tirar eu acho
        // int value = ((int)iterPerDay * days) / points;
        float time_save = (float)t / (float)iterPerDay;
        
        if(V < 0){
            V = -1 * V;
        }
        
        // Saving data
        //TODO ALT na próxima coluna ou novo arquivo/ verificar no plotDados como que plota os dados
        output_file << time_save << "," << V << std::endl;
        
        /**
        * ODEs
        */
        //TODO verificar se tá certo os parenteses com a Bárbara p.74 da tese CERTO
        T = (s - d * T - beta * V * T) 
                * deltaT + T;
        V = ((1 - epsilon_s) * rho * calcIntegral(I,Rp,Rt) - c * V) 
                * deltaT + V;
        A_alt = (alpha_alt * delta * calcIntegral(I) - c_alt * A_alt) 
                * deltaT + A_alt; // MATHEUS linha 100. delta da linha 70
        
        //printf("t = %ld V = %lf T = %lf \n", t,V,T);
        
        /**
        * Begin age loop
        */
        double rho_func;
        double delta_func;
        if(varrho)
            rho_func = 0.0;
        else
            rho_func = rho;
        
        a = 0;

        Rp[1][a]  = 0.0;
        Rn[1][a]  = 0.0;
        Rt[1][a]  = 1.0;
        I[1][a]   = beta * V *T;

        if(vardelta)
            delta_func = 0.01;
        else
            delta_func = delta;

        for(a = 1; a < AGE-1; a++){
            if(varrho){
                if(((double)a * deltaA) < tau){
                    rho_func = 0;
                } else {
                    rho_func = (1 - exp(-k * (((double)a * deltaA) - tau))) * rho;
                }
            } else {
                rho_func = rho;
            }
            
            if(vardelta)
                delta_func = delta * (1 - exp(-k * (a * deltaA)));
            else
                delta_func = delta;
            
        
            I[1][a]  = (-delta_func * I[0][a] 
                        - (I[0][a] - I[0][a - 1]) / (deltaA)) //TODO Diferenças finitas
                        * deltaT + I[0][a];

		    Rn[1][a] = ((1 - epsilon_r) * r * Rp[0][a] - (1 - epsilon_r) * r * Rp[0][a] * (Rn[0][a] / Rmax) - kappa_c * mu_p * Rn[0][a] 
                        - (Rn[0][a] - Rn[0][a - 1]) / (deltaA)) 
                        * deltaT + Rn[0][a];

		    Rp[1][a] = ((1 - epsilon_alpha) * alpha * Rn[0][a] + sigma * Rt[0][a] - theta * Rp[0][a] - (1 - epsilon_s)
                        * rho_func * Rp[0][a] - kappa_c * mu_p * Rp[0][a] 
                        - (Rp[0][a] - Rp[0][a - 1]) / (deltaA)) 
                        * deltaT + Rp[0][a];

            Rt[1][a] = (theta * Rp[0][a] - sigma * Rt[0][a] - (1 - epsilon_s) * rho_func * Rt[0][a] - kappa_t 
                        * mu_t * Rt[0][a] - ((Rt[0][a] - Rt[0][a - 1]) / (deltaA))) 
                        * deltaT + Rt[0][a];
	    }

	    update(I);
	    update(Rn);
	    update(Rp);
	    update(Rt);

        t++;

    } 
    while(t < (iterPerDay * days));

    output_file.close();

    return;
 }