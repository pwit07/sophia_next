#include <iostream>
#include <fstream>
#include <string>

#include "sophia_interface.h"

const double pi = 3.1415926;

const double eta_0 = 0.313;

// speed of light value [cm/s]
const double c = 2.99792458e10; 

// m_p c^2 value [GeV]
const double mpc2 = 0.938272088;

// (m_p c^2)^2 value [GeV^2]
const double mp2c4 = 0.880354511;

double eta(double Ep, double eps);

double eps_prime(double eta, double theta);

double crossection(double eps_prime);

// class computing histogram
class sophia_histogram{
    private:
        bool check_make_his = false;
        int eventCounter = 0, nbin, partID;
        double eta, theta, Ep, log10xmin, log10xmax, dlog10x;
        std::vector <int> eventID;
        std::vector <int> PDG;
        std::vector <double> EGeV;
        std::vector <int> H;
    public:
        sophia_histogram(int npartID, double neta, double ntheta, double nEp, int nnbin, double nxmin, double nxmax);
        sophia_histogram(int npartID, int nnbin, double nsmin, double nsmax);
        ~sophia_histogram(){};
        
        int find_index(double x);
        double x_i(int index);
        
        void add(const sophiaevent_output &seo);
        void print_data();
        void make_hist();
        void make_hist_s();
        void print_hist();

        double value(double x);
        double integrand(double x);

        void save_integrand(std::string path);
};

std::string partName(int partID);

// double integrand(sophia_histogram &HIST, int partID, double x, double eta, double theta, int N = 10000, int nbin = 10, double Ep = 1., double xmin = 1.e-6);

double compute_integrand(int partID, double x, double eta, double theta, int N = 10000, int nbin = 10, double Ep = 1., double xmin = 1.e-6);

double compute_Phi(int partID, double x, double eta, int N = 10000, int nbin = 10, double Ep = 1., int Ntheta = 2*18, double xmin = 1.e-6);

void Phi2File(std::string path, int partID, double x_a, double x_b, double eta, int N){
    std::ofstream outfile;
    outfile.open(path.c_str());

    // std::string id = partName(partID);
    // outfile << "#Particle type: " << id << "\n";
    
    outfile << "#x\t" << "Phi [cm^3/s]\n";
    
    double log10xmin = std::log10(x_a), log10xmax = std::log10(x_b), dlog10x = (log10xmax - log10xmin )/double(N);
    double log10x = log10xmin, x = std::pow(10.,log10x);

    while (log10x <= log10xmax)
    {
        outfile << x << "\t" << x*compute_Phi(partID,x,eta,1e4,100,1e2,18) << "\n";
        log10x+=dlog10x;
        x = std::pow(10.,log10x);
    }
    
    outfile.close();
    return;
}

// void integrand2File(std::string path, int partID, double x_a, double x_b, double eta, double theta, int N);

// void compute_s_tabs(int partID, int N, double s, int nbin=20, double smin=1.e-3, double smax=1.e3);

void compute_tabs(std::string path, int partID, int N = 10000, int nbin = 10, double Ep = 1., double xmin = 1.e-6);

// If you want to see warnings please uncomment line 283 in sophia_interface.cpp file

int main() {
    std::cout.precision(10);

    // double x = 1.e-1;
    // double eta = 1.5*eta_0;

    // std::cout<<"Phi = "<<x*compute_Phi(22,x,eta,10000,10,1e4,18)<<"\n";

    // Phi2File("./src/kelner_aharonian_2008/fig2_values/gamma_1.5eta0_sophia.txt",22,1.e-4,1.,1.5*eta_0,20);
    // Phi2File("./src/kelner_aharonian_2008/fig2_values/gamma_30eta0_sophia.txt",22,1.e-4,1.,30.*eta_0,30);

    // Phi2File("./src/kelner_aharonian_2008/fig3_values/positron_1.5eta0_sophia.txt",-11,1.e-4,1.,1.5*eta_0,50);
    // Phi2File("./src/kelner_aharonian_2008/fig3_values/positron_30eta0_sophia.txt",-11,1.e-4,1.,30.*eta_0,50); ///!!!

    // Phi2File("./src/kelner_aharonian_2008/fig4_values/muon_antineutrino_1.5eta0_sophia.txt",-14,1.e-4,1.,1.5*eta_0,50);
    // Phi2File("./src/kelner_aharonian_2008/fig4_values/muon_antineutrino_30eta0_sophia.txt",-14,1.e-4,1.,30.*eta_0,50);

    // Phi2File("./src/kelner_aharonian_2008/fig5_values/muon_neutrino_1.5eta0_sophia.txt",14,1.e-4,1.,1.5*eta_0,50);
    // Phi2File("./src/kelner_aharonian_2008/fig5_values/muon_neutrino_30eta0_sophia.txt",14,1.e-4,1.,30.*eta_0,50); // !!!

    // Phi2File("./src/kelner_aharonian_2008/fig_values/electron_3eta0_sophia.txt",11,1.e-4,1.,3*eta_0,50);
    //  Phi2File("./src/kelner_aharonian_2008/fig_values/electron_10eta0_sophia.txt",11,1.e-4,1.,10*eta_0,50);
    // Phi2File("./src/kelner_aharonian_2008/fig_values/electron_30eta0_sophia.txt",11,1.e-4,1.,30*eta_0,50);
    // Phi2File("./src/kelner_aharonian_2008/fig_values/electron_100eta0_sophia.txt",11,1.e-4,1.,100*eta_0,50);

    // compute_s_tabs(22,1e0,1e1);


    // integrand2File("./src/integrand_values/gamma_30eta0_sophia.txt",22,1.e-4,1.,40.*eta_0,50.,30);
    // integrand2File("./src/integrand_values/gamma_30eta0_sophia_.txt",22,1.e-4,1.,40.*eta_0,50.,30);

    compute_tabs("./src/integrand_values/",22,1e6,30,1e2,1e-6);
    compute_tabs("./src/integrand_values/",11,1e6,30,1e2,1e-6);
    compute_tabs("./src/integrand_values/",-11,1e6,30,1e2,1e-6);
    compute_tabs("./src/integrand_values/",12,1e6,30,1e2,1e-6);
    compute_tabs("./src/integrand_values/",-12,1e6,30,1e2,1e-6);
    compute_tabs("./src/integrand_values/",14,1e6,30,1e2,1e-6);
    compute_tabs("./src/integrand_values/",-14,1e6,30,1e2,1e-6);
    
    return 0;
}

double eta(double Ep, double eps){
// ****** INPUT ******************************************
// Ep = energy of incident proton (in lab frame) [GeV]
// eps = energy of incident photon (in lab frame) [Gev]
// ****** OUTPUT *****************************************
// eta = normalized value of Ep*eps []
// *******************************************************
    return 4.*Ep*eps/mp2c4;
}

double eps_prime(double eta, double theta){
// ****** INPUT ******************************************
// eta = eta(Ep,eps) []
// theta = angle between incident proton and photon (in lab frame) [degrees]
// ****** OUTPUT *****************************************
// eps_prime = energy of photon in PRF (proton rest frame) [GeV]
// *******************************************************
    return 0.25*eta*mpc2*(1-std::cos(theta* pi / 180.));
}

double crossection(double eps_prime){
// ****** INPUT ******************************************
// eps_prime = eps_prime(eta, theta) energy of photon in PRF (proton rest frame) [GeV]
// ****** OUTPUT *****************************************
// total corssection of photomeson interaction [cm^2]
// *******************************************************
    static sophia_interface SI; 
                                          // micro barn to cm^2
    return SI.crossection(eps_prime,3,13)*1.e-30; 
}

sophia_histogram::sophia_histogram(int npartID, double neta, double ntheta, double nEp, int nnbin, double nxmin, double nxmax){
    partID = npartID;

    eta = neta;
    Ep = nEp;
    theta = ntheta;

    nbin = nnbin;

    log10xmin = std::log10(nxmin);
    log10xmax = std::log10(nxmax);
    dlog10x = (log10xmax - log10xmin)/nbin;

    for(int i=0;i<nbin;i++)
    {
        H.push_back(0);
    }
}

sophia_histogram::sophia_histogram(int npartID, int nnbin, double nsmin, double nsmax){
    partID = npartID;

    nbin = nnbin;

    log10xmin = std::log10(nsmin);
    log10xmax = std::log10(nsmax);
    dlog10x = (log10xmax - log10xmin)/nbin;

    for(int i=0;i<nbin;i++)
    {
        H.push_back(0);
    }
}

void sophia_histogram::add(const sophiaevent_output &seo){
// ****** OUTPUT *****************************************
// add sophiaevent_output seo data to this class
// *******************************************************
    for (int i = 0; i < seo.Nout; ++i) {
        eventID.push_back(eventCounter+1);
        PDG.push_back(ID_sophia_to_PDG(seo.outPartID[i]));
        EGeV.push_back(seo.outPartP[3][i]);    
    }
    check_make_his = false;
    eventCounter++;
}

void sophia_histogram::print_data(){
// ****** OUTPUT *****************************************
// printed data from the all sophiaevent_output seo passed to this class
// *******************************************************
    for(int i=0;i<eventID.size();i++)
    {
        std::cout << eventID[i] << " " << PDG[i]<< " " << EGeV[i] << "\n";
    }
}

void sophia_histogram::make_hist(){
// ****** OUTPUT *****************************************
// table H[i] filled with a number of particles corresponding to the value x_i(i) 
// *******************************************************
    static bool check_next_use = false;

    if(check_make_his == false)
    {
        if(check_next_use==false)
        {
            check_next_use = true;
        }
        else
        {
            for(int i=0;i<nbin;i++)
            {
                H[i] = 0;
            }
        }

        for(int i=0;i<eventID.size();i++)
        {
            if(PDG[i] == partID)
            {
                H[find_index(EGeV[i]/Ep)]++;
            }
        }
        check_make_his = true;
    }
}

void sophia_histogram::make_hist_s(){
// ****** OUTPUT *****************************************
// table H[i] filled with a number of particles corresponding to the value x_i(i) 
// *******************************************************
    static bool check_next_use = false;

    if(check_make_his == false)
    {
        if(check_next_use==false)
        {
            check_next_use = true;
        }
        else
        {
            for(int i=0;i<nbin;i++)
            {
                H[i] = 0;
            }
        }

        for(int i=0;i<eventID.size();i++)
        {
            if(PDG[i] == partID)
            {
                H[find_index(EGeV[i])]++;
            }
        }
        check_make_his = true;
    }
}

int sophia_histogram::find_index(double x){
// ****** INPUT ******************************************
// EGeV = energy of produced particle (in lab frame) [GeV]
// ****** OUTPUT *****************************************
// index = index corresponding to the energy EGeV of produced particle []
// *******************************************************
    int index = static_cast<int>((std::log10(x)-log10xmin)/dlog10x);
        
    if(index >= nbin){index = nbin-1;}
    if(index < 0){index = 0;}

    return index;
}

double sophia_histogram::x_i(int index){
// ****** INPUT ******************************************
// index = index corresponding to the energy EGeV of produced particle []
// ****** OUTPUT *****************************************
// EGeV = energy of produced particle (in lab frame) [GeV]
// *******************************************************
    return std::pow(10.,log10xmin+double(index)*dlog10x);
}

void sophia_histogram::print_hist(){
// ****** OUTPUT *****************************************
// printed histogram in table
// *******************************************************
    for(int i=0;i<nbin;i++)
    {
        std::cout<<x_i(i)<<" "<<H[i]<<"\n";
    }
}

double sophia_histogram::value(double x){
// ****** INPUT ******************************************
// x = EGeV / Ep = energy EGeV of produced particle divided by energy of proton (both in lab frame) [GeV/GeV]
// ****** OUTPUT *****************************************
// histogram []
// *******************************************************
// dN/dlog(x) = dN/dlog10(x) / log(10) = E * dN/dE
// dN/dx = (dN/dlog(x))/x = (dN/dlog10(x))/(x*log(10))
// *******************************************************
    return (((double(H[find_index(x)])/double(eventCounter))/dlog10x)/x)/std::log(10);
}

double sophia_histogram::integrand(double x){
// ****** INPUT ******************************************
// partID = sophia ID of produced particles: 
//  1000010010 -> proton
// -1000010010 -> anti-proton
//  1000000010 -> neutron
// -1000000010 -> anti-neutron
//          22 -> photon
//          11 -> electron
//         -11 -> positron
//          12 -> nu_e
//         -12 -> anti-nu_e
//          14 -> nu_mu
//         -14 -> anti-nu_mu
//         111 -> pi0
//         211 -> pi+
//        -211 -> pi-
// 
// ****** OUTPUT *****************************************
// integrand for computation Phi function (integral over the angles) [cm^3/(s rad)]
// *******************************************************
    return c*(1.-std::cos(theta*pi/180.))*crossection(eps_prime(eta,theta))*(this->value(x));
}

void sophia_histogram::save_integrand(std::string path){
    static double x;

    std::ofstream outfile;
    outfile.open(path.c_str());

    outfile << "#x\t" << "f [cm^3/s]\n";
    
    for(int i=0;i<nbin;i++)
    {
        x = x_i(i);
        outfile<<x<<" "<<this->integrand(x)<<"\n";
    }

    outfile.close();
}

std::string partName(int partID){
    switch (partID)
    {
        case 1000010010:
            return "proton";
            break;
        case 1000000010:
            return "neutron";
            break;
        case -1000010010:
            return "antiproton";
            break;
        case -1000000010:
            return "antineutron";
            break;
        case 22:
            return "photon";
            break;
        case 11:
            return "electron";
            break;
        case -11:
            return "positron";
            break;
        case 12:
            return "electron_neutrino";
            break;
        case -12:
            return "electron_antineutrino";
            break;
        case 14:
            return "muon_neutrino";
            break;
        case -14:
            return "muon_antineutrino";
            break;
        case 111:
            return "pi0";
            break;
        case 211:
            return "pi+";
            break;
        case -211:
            return "pi-";
            break;
        default:
            throw std::runtime_error("Unknown particle ID!");
    }
}

// double integrand(sophia_histogram &HIST, int partID, double x, double eta, double theta, int N, int nbin, double Ep, double xmin){
// // ****** INPUT ******************************************
// // partID = sophia ID of produced particles: 
// //  1000010010 -> proton
// // -1000010010 -> anti-proton
// //  1000000010 -> neutron
// // -1000000010 -> anti-neutron
// //          22 -> photon
// //          11 -> electron
// //         -11 -> positron
// //          12 -> nu_e
// //         -12 -> anti-nu_e
// //          14 -> nu_mu
// //         -14 -> anti-nu_mu
// //         111 -> pi0
// //         211 -> pi+
// //        -211 -> pi-
// // 
// // ****** OUTPUT *****************************************
// // integrand for computation Phi function (integral over the angles) [cm^3/(s rad)]
// // *******************************************************
//     return c*(1.-std::cos(theta*pi/180.))*crossection(eps_prime(eta,theta))*HIST.value(x);
// }

double compute_integrand(int partID, double x, double eta, double theta, int N, int nbin, double Ep, double xmin){
// ****** OUTPUT *****************************************
// integrand(theta) value
// *******************************************************    
    static sophia_interface SI;
    static sophiaevent_output seo;

    double eps = (mp2c4*eta/(4.*Ep));

    sophia_histogram HIST(partID,eta,theta,Ep,nbin,xmin,1.);

    for(int i=0;i<N;i++){
        seo = SI.sophiaevent_mod(Ep, eps, theta, false);
        HIST.add(seo);
    }

    HIST.make_hist();
    // HIST.print_hist();

    // return integrand(HIST,partID,x,eta,theta,N,nbin,Ep);
    return HIST.integrand(x);
}

double compute_Phi(int partID, double x, double eta, int N, int nbin, double Ep, int Ntheta, double xmin){
// ****** OUTPUT *****************************************
// value of Phi function from K&A 2008 [cm^3/s]
// *******************************************************    
    static sophia_interface SI;
    static sophiaevent_output seo;

    double eps = (mp2c4*eta/(4.*Ep));
    std::vector <double> Y;

    Y.push_back(0.);

    double theta_i = 0., dtheta = 180./double(Ntheta);

    for(int i=0;i<Ntheta-1;i++)
    {
        theta_i+=dtheta;

        Y.push_back(compute_integrand(partID,x,eta,theta_i,N,nbin,Ep,xmin)*std::sin(theta_i*pi/180.));
    }

    Y.push_back(0.);

    double F = 0.;
    
    for(int i=0;i<Y.size()-1;i++)
    {
        F+=(Y[i+1]+Y[i]);
    }
    
    return (F*(pi/180.)*dtheta*0.5)/(2.);
}

// void compute_s_tabs(int partID, int N, double s, int nbin, double smin, double smax)
// {
//     static sophia_interface SI;
//     static sophiaevent_output seo;

//     sophia_histogram HIST(partID,nbin,smin,smax);

//     for(int i=0;i<N;i++){
//         seo = SI.sophiaevent_mod2(s, false);
//         HIST.add(seo);
//     }
//     // HIST.print_data();
//     HIST.make_hist_s();
//     HIST.print_hist();
    
//     return;
// }

// void integrand2File(std::string path, int partID, double x_a, double x_b, double eta, double theta, int N){
//     std::ofstream outfile;
//     outfile.open(path.c_str());

//     std::string id = partName(partID);
//     // outfile << "#Particle type: " << id << "\n";
    
//     std::cout<<"=================================================================\n";
//     std::cout<<"id = "<<id<<" | x_a = "<<x_a<<" | x_b = "<<x_b<<" | "<<eta<<"  | theta = "<<theta<<" deg\n";
//     std::cout<<"=================================================================\n";
    
//     std::cout << "#x\t" << "f [cm^3/s]\n";
//     outfile << "#x\t" << "f [cm^3/s]\n";
    
//     double log10xmin = std::log10(x_a), log10xmax = std::log10(x_b), dlog10x = (log10xmax - log10xmin )/double(N);
//     double log10x = log10xmin, x = std::pow(10.,log10x);
//     double f;

//     while (log10x <= log10xmax)
//     {
//         f = compute_integrand(partID,x,eta,theta,1e4,100,1e2,x_a);
        
//         outfile << x << "\t" << f << "\n";
//         std::cout << x << "\t" << f << "\n";

//         log10x+=dlog10x;
//         x = std::pow(10.,log10x);
//     }
    
//     outfile.close();
//     std::cout<<"=================================================================\n";
    
//     return;
// }

void compute_tabs(std::string path, int partID, int N, int nbin, double Ep, double xmin){
    static std::string _eta[22] = {"1.1", "1.2", "1.3", "1.4", "1.5", "1.6", "1.7", "1.8", "1.9", "2.0", "3.0", "4.0", 
                                   "5.0", "6.0", "7.0", "8.0", "9.0", "10.0", "20.0", "30.0", "40.0", "100.0"};
    static sophia_interface SI;
    static sophiaevent_output seo;

    double eps;
    std::string full_path;

    std::string id = partName(partID);
    
    for(int i=0;i<22;i++)
    {
        // std::cout<<_eta[i]<<"\n";
        for(int theta=0;theta<=180;theta+=10)
        {
            full_path = path+id+"/"+id+"_"+_eta[i]+"eta0_"+std::to_string(theta)+"theta.txt";
            std::cout<<full_path<<"\n";

            eps = (mp2c4*std::stod(_eta[i])/(4.*Ep));

            sophia_histogram HIST(partID,std::stod(_eta[i])*eta_0,theta,Ep,nbin,xmin,1.);

            for(int i=0;i<N;i++){
                seo = SI.sophiaevent_mod(Ep, eps, theta, false);
                HIST.add(seo);
            }

            HIST.make_hist();
            HIST.save_integrand(full_path);
        }
    }

}