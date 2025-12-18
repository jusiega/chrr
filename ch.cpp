//CH_CPP
//kompilacja z debugiem
//icx /march=native /mtune=nativ /Ofast /fp:fast=2 /Qopenmp /Qopenmp-simd /Wall /Wextra /Wpedantic /Wshadow /Wconversion /F0x2137000  ch.cpp -o ch.exe
//do debug dodać /fsanitize=address
//albo /link /STACK:0x2137000 na koniec do zwiększania stacka
#if defined(_WIN32) || defined(_WIN64)
    #define  _CRT_SECURE_NO_WARNINGS
    #define NOMINMAX
    #include <Windows.h>
#endif

#include <iostream>
#include <complex>
#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <fstream>
#include <vector>
#include <string>
#include <ctime>
#include <filesystem>
#include <algorithm>
#include <omp.h>
#include <cfloat>
#include <math.h>



#define  tttt 32//ile wątków wariacie

#define  pot "ch.chi.potenszjals2.h"
#define  inp "ch.chi.h"
#include inp


using namespace std::complex_literals;



struct characteristic
{
    double E_c, dE_c, x1, p1, x9, p9, dt_c;// E_next, E_previous;// procesem myślowym odbytym o piątej w nocy nie ma sensu rozważać E_next
    double s1, s9;
};

double f0(double x, double p)
{
    bool ctrl=false;
    double ligma=4.;//ligma balsss
    if(p0+ligma*sqrt(1./2./β)<pmax || p0-ligma*sqrt(1./2./β)>-pmax) ctrl=true;
    if(false==ctrl) std::cout<<"AAAAA ŹLEEEEEE PĘD UCIEKŁ POZA SKALĘ\n", exit(EXIT_FAILURE);
    
    double a=2.*λ/(1.-ωr*ωr);
    double b=2.*β/(1.-ωr*ωr);
    double ω=sqrt(a)*sqrt(b)*ωr;
    //for (int i=0;i<npointsx;i++)//i odpowiada za położenia (zawsze i wszędzie)
    {
      //  for (int j=0;j<npointsp;j++)//j odpowiada za pedy
        {
            {
                //fxp[i*npointsp+j]=
                return sqrt(a*b-ω*ω)/(2.*π)*exp(-a/2.*(x-x0)*(x-x0)-b/2.*(p-p0)*(p-p0)-ω*(x-x0)*(p-p0));
            }
        }
        //double obw=1.-exp(-a*pow(x[i]-xmin,p))-exp(-a*pow(x[i]-xmax,p));
        //std::cout<<obw<<"\n";
        //f[i]*=obw;//obwiednia gasząca zbędna
    }
}

double H(double x, double p)
{
    return p*p/2./m+potencjał(x,tpp);
}


int main(int argc, char* argv[])
{
    #if defined(_WIN32) || defined(_WIN64)
    SetConsoleOutputCP(65001);
    #endif
    
    std::ifstream drff("DATA_ROOT_FOLDER", std::ios::in);
    std::string DATA_ROOT_FOLDER;
    std::getline(drff, DATA_ROOT_FOLDER);
    //std::cout<<DATA_ROOT_FOLDER<<"\n";//dziala
    
    std::string folder="data";
    if(argc==2) folder=argv[1];
    folder=DATA_ROOT_FOLDER+"/"+folder;
    std::cout<<"output catalog: "<<folder<<"/\n";
    std::filesystem::create_directory(folder);
    std::filesystem::copy(inp,folder+"/"+inp,std::filesystem::copy_options::overwrite_existing);  
    std::filesystem::copy(pot,folder+"/"+pot,std::filesystem::copy_options::overwrite_existing);  
    //stworzyliśmy katalog wyjściowy i skopiowaliśmy tam input
    //problem jest tylko taki że kopiuje on plik wejściowy z wtedy gdy program się kompilował
    //a nie z wtedy co się wykonywał
    //wystarczy se przekomplilować przed każdym odpaleniem i nie ma problemu    
    
    //inne pliki tam potrzebne: E0.txt, info.txt. 
    //fxpl.txt można sb darować bo NIE będzie się to zmieniało w cale i podobnie na razie split.txt
    
    
    

    //skanowanie brzegu
    #include "ch.brz.cpp"
    //dostaje z tego const double E_L, E_P, E_max_którą_mogę_odpowiedzialnie_symulować
    //właściwie tylko ta ostatnia jest ważna
    
        
    
    
    ///////rozkład energii-początkowa
    //const int spoints2=16000;
    const double eeps=(E_max_którą_mogę_odpowiedzialnie_symulować-E_min_const)/double(beans-1);
    static double PE[beans]{};//p-stwo
    static double DE[beans]{};//instrybuanta
    //wyrażenie zaś na indeks energii wewnątrz beans [emin: emax] -> [0:beans-1] jest
    //indeks(H)=int((H-E_min+0.5*eeps)/eeps)
    
    //double Eeeeeeeeiiiii[beans]{};
    //for (int k=0;k<beans;k++)
    //{
    //    Eeeeeeeeiiiii[k]=E_min+double(k)*eeps;
    //}
    
    //double dx2=(xmax-xmin)/2./double(spoints2);
    //double dp2=pmax/double(spoints2);
    //double norm=0;//a se sprawdze przy okazji
    //a nie muszę skoro histogram energii się sumuje do idealnej perfekcyjnej jedynki
    //to całkowanie możnaby trochę ładniej na brzzegu
    const double wprobability=prop_to_denisity_energy_spectrum_weight/(prop_to_denisity_energy_spectrum_weight+system_energy_spectrum_weight+equal_energy_steps_weight);
    const double wuniform=system_energy_spectrum_weight/(prop_to_denisity_energy_spectrum_weight+system_energy_spectrum_weight+equal_energy_steps_weight);
    const double wequidistant=equal_energy_steps_weight/(prop_to_denisity_energy_spectrum_weight+system_energy_spectrum_weight+equal_energy_steps_weight);
    
    double pokryciecałki=1;
    for (int i=0;i<spoints_x;i++)
    {
        for (int j=0;j<spoints_p;j++)
        {
            x=double(i)*dx+xmin;
            p=double(j)*dp;
            if (i==0 ||i==spoints_x-1 || j==spoints_p-1) {/*std::cout << "x p "<<x<<" "<<p<<"\n";*/ pokryciecałki=0.5;}
            else if ((i==0 && j==spoints_p-1) || (i==spoints_x-1 && j==spoints_p-1)) pokryciecałki=0.25;
            else pokryciecałki=1.;
            //przedziały hist lewostronnie domknięte
            PE[int((H(x,p)-E_min_const+0.5*eeps)/eeps)]+=pokryciecałki*(
            wprobability*dx*dp*f0(x,p)/eeps
            +wuniform*dx*dp/(xmax-xmin)/(2.*pmax)/eeps);
            //double eeeeee_można_łatwo_utknąć_w_korku=H(x,p);
            //for (int k=0;k<beans;k++)
            //{
            //    //if(eeeeee_można_łatwo_utknąć_w_korku<=Eeeeeeeeiiiii[k]+eeps/2.&&eeeeee_można_łatwo_utknąć_w_korku>Eeeeeeeeiiiii[k]-eeps/2.) PE[k]+=dx*dp*f0(x,p), break;
            //}
        }
        for (int j=1;j<spoints_p;j++)///od jeden żeby nie iść dwa razy do tej samej rzeki
        {
            x=double(i)*dx+xmin;
            p=-double(j)*dp;
            if (i==0 ||i==spoints_x-1 || j==spoints_p-1) {/*std::cout << "x p "<<x<<" "<<p<<"\n";*/ pokryciecałki=0.5;}
            else if ((i==0 && j==spoints_p-1) || (i==spoints_x-1 && j==spoints_p-1)) pokryciecałki=0.25;
            else pokryciecałki=1.;
            PE[int((H(x,p)-E_min_const+0.5*eeps)/eeps)]+=pokryciecałki*(
            wprobability*dx*dp*f0(x,p)/eeps
            +wuniform*dx*dp/(xmax-xmin)/(2.*pmax)/eeps);
        }
    }
    
    for (int k=0;k<beans;k++)
    {
        PE[k]+=wequidistant/(E_max_którą_mogę_odpowiedzialnie_symulować-E_min_const);
    }
    
    
    
    double soom=0;
    for (int k=0;k<beans;k++)
    {
        double cx=1;
        if (k==0||k==beans-1) cx=0.5;
        soom+=eeps*PE[k]*cx;//bez de bo to histogram pefekt jedynka
    }
    std::cout<<"norma rozkładu energii: "<<soom<<"\n";
    for (int k=0;k<beans;k++)
    {
        //dwie pętle żeby nei nałożyć błędów numerycznych
        for (int q=0;q<=k;q++)
        {
            double cx=1;
            if (q==0||q==k) cx=0.5;
            DE[k]+=eeps*PE[q]*cx;
        }
    }
    
    //printuj rozkład prawdopodobieństwa energiiii
    {
        FILE* p_od_E;
        p_od_E=fopen((folder+"/pe.dat").c_str(),"w");
        
        for (int k=0;k<beans;k++)
        {
            fprintf(p_od_E,"%d %e %e %e\n", k, E_min_const+double(k)*eeps, PE[k], DE[k]);
        }
    }
    //to do: kroki energetyczne respektujące rozkład energetyczny warunku począstkowego 
    //na razie to zlewam
    //robimy równe kroki energetyczne
    //moment to bez znaczenia
    
    //i tak musze odwrócić instrybuantę
    //pacze na tablice DE[beans]
    const double dDE=(DE[beans-1]-DE[0])/double(n_energy_levels);
    //std::cout<<"dDE "<<dDE<<"\n"; cośtam jest
    static double fractionalenergyindices[n_energy_levels];//nie biore minimalnej energii bo to będzie dążyć do absurdalnie krótkiej charakterystyki
    //może jeszcze zmienię zdanie ntt
    int base=0;
    int energycounter=0;
    for (int i=0;i<=n_energy_levels-1;i++)
    {
        double D=DE[0]+dDE*double(i+1);
        if (i==n_energy_levels-1) D=DE[beans-1];//ewentualnie D-=1E-15;//żeby nie wydarzyło się to co sie ma nie wydarzyć
        bornholm:;
        if (base<beans-1) 
        {
            if (DE[base+1]>=D && DE[base]<=D)
            {
                //std::cin.get();
                //std::cout<<"eeeee";//nigdy tu nie wchodzi
                double DDL=(D-DE[base])/(DE[base+1]-DE[base]);
                //std::cout<<DE[base+1]<<" "<<DE[base]<<" "<<D<<" "<<DDL<<"\n";
                fractionalenergyindices[i]=(double(base)+DDL)/double(beans);
                energycounter++;
            }
            else
            {
                base++;
                goto bornholm;
            }
        }
        else
        {
            std::cout<<"to miało się nidgy nie wydarzyć";
            exit(EXIT_FAILURE);
        }
    }
    std::cout<<"odnalezniono poziomów energetycznych: "<<energycounter<<" od "<<fractionalenergyindices[0]*eeps*double(beans)+E_min_const<<" do "<<fractionalenergyindices[n_energy_levels-1]*eeps*double(beans)+E_min_const<<"\n";
    static double energylevels[n_energy_levels];
    for (int i=0;i<=n_energy_levels-1;i++)
    {
        energylevels[i]=fractionalenergyindices[i]*eeps*double(beans)+E_min_const;//-0.5*eeps;
        //std::cout<<"frakcjonalny indeks "<<fractionalenergyindices[i]<<"\n";
        //std::cout<<"energia "<<energylevels[i]<<"\n";
        //std::cin.get();
    }
    //printuj poziomy energetyczne
    {
        FILE* enerdżylevels_s;
        enerdżylevels_s=fopen((folder+"/energylevs.dat").c_str(),"w");
        
        //for (int n=0;n<n_characteristics;n++)
        for (int n=0;n<=n_energy_levels-1;n++)
        {
            //do gnuplota with vectors nohead idzie tak:
            //4 columns:  x  y  xdelta  ydelta
            //fprintf(enerdżylevels_s,"%d %e %d %d\n", OUR_system[n].s1, OUR_system[n].E_c, OUR_system[n].s9-OUR_system[n].s1, 0);
            fprintf(enerdżylevels_s,"%d %e %d %d\n", 0, energylevels[n], spoints_p*2+spoints_x, 0);
        }
    }
    
    //const int maksymalnamożliwaliczbacharakterystyk=n_energy_levels*ileminimówwariacie;
    static characteristic OUR_system[4*n_energy_levels];//do 4 minimów obsłuży
    //characteristic* __restrict R_system; //czy to jest szybciej? może zobaczymy jak mi się będzie chciało
    
    {//sprawdzam odszukanie minimów i maksimów
        
        std::cout<<"znaleziono minimów: "<<ileminimówwariacie<<", w (x p s E): ";
        for (int n=0;n<ileminimówwariacie;n++)
        {
            std::cout<<",   ("<<XX[indeksy_minimów[n]]<<", "<<PP[indeksy_minimów[n]]<<", "<<indeksy_minimów[n]<<", "<<energiabrzegowa[indeksy_minimów[n]]<<")";
        }
        std::cout<<".\n";
        
        std::cout<<"znaleziono maksimów: "<<ilemaksimówwariacie<<", w (x p s E): ";
        for (int n=0;n<ilemaksimówwariacie;n++)
        {
            std::cout<<",   ("<<XX[indeksy_maksimów[n]]<<", "<<PP[indeksy_maksimów[n]]<<", "<<indeksy_maksimów[n]<<", "<<energiabrzegowa[indeksy_maksimów[n]]<<"), ";
        }
        std::cout<<".\n";
        if (ileminimówwariacie+1!=ilemaksimówwariacie)
        {
            std::cout<<"jakiś dziwny ten układ żeś zadał"<<"\n";
            exit(EXIT_FAILURE);
        }
        
    }
    
    //odszukuję początki i końce charakterystyk
    int characteristic_count=0;
    int klold, kpold;
    bool edging;
    double* __restrict__ elevels=energylevels;
    double* __restrict__ eb=energiabrzegowa;
    for (int i=0;i<ileminimówwariacie;i++)
    {
        edging=false;
        const int s0=indeksy_minimów[i];
        const int L0=indeksy_maksimów[i];//jako że tak se szedłem zbierając minima to w tym samym indeksie bd minimum i maksimum po jego lewicy
        klold=s0, kpold=s0;
        int j=s0;
        while (eb[j]<E_max_którą_mogę_odpowiedzialnie_symulować && j>L0) j--;//j wskazuje NA maksimum albo na pierwszy indeks który przekracza E_max_blablabla
        int eli=0;
        while (elevels[eli]<eb[s0]) eli++;//eli teraz wskazuje na pierwszy indeks poz. energii powyżej energii itego minimum
        shadowmoses:;
        OUR_system[characteristic_count].E_c=elevels[eli];
        if (eli!=0) OUR_system[characteristic_count].dE_c=elevels[eli]-elevels[eli-1]; else OUR_system[characteristic_count].dE_c=elevels[eli]-E_min_const;
        int k=klold;
        while (eb[k]<elevels[eli]&&k>j) k--;
        double fr=(elevels[eli]-eb[k+1])/(eb[k]-eb[k+1]);
        OUR_system[characteristic_count].x1=XX[k+1]-dx*fr*(XX[k+1]-XX[k]);
        OUR_system[characteristic_count].p1=PP[k+1]-dp*fr*(PP[k+1]-PP[k]);
        OUR_system[characteristic_count].s1=(double(k+1)-fr);
        klold=k+1;//safety net
        if (k==j) edging=true;
        k=kpold;
        while (eb[k]<elevels[eli]&&k<2*spoints_p+spoints_x-1) k++;
        fr=(elevels[eli]-eb[k-1])/(eb[k]-eb[k-1]);
        OUR_system[characteristic_count].x9=XX[k-1]+dx*fr*(XX[k+1]-XX[k]);
        OUR_system[characteristic_count].p9=PP[k-1]+dp*fr*(PP[k+1]-PP[k]);
        OUR_system[characteristic_count].s9=(double(k-1)+fr);
        kpold=k-1;
        characteristic_count++;
        if (!edging) {eli++; goto shadowmoses;}
    }
    //jest prawie dobrze tylko ten ostatni jest spierdolony. pierdole to na razie
    
    //printuj poziomy energetyczne z granicami
    {
        FILE* enerdżylevels_s;
        enerdżylevels_s=fopen((folder+"/energylevsbounds.dat").c_str(),"w");
        
        //for (int n=0;n<=n_energy_levels-1;n++)
        for (int n=0;n<characteristic_count;n++)
        {
            //do gnuplota with vectors nohead idzie tak:
            //4 columns:  x  y  xdelta  ydelta
            double x, y, xdelta, ydelta;
            x=OUR_system[n].s1;
            y=OUR_system[n].E_c;
            xdelta=OUR_system[n].s9-x;
            ydelta=0;
            
            fprintf(enerdżylevels_s,"%e %e %e %e\n", x, y, xdelta, ydelta);
        }
    }
}
















