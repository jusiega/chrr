#ifndef CH_CHI_H
#define CH_CHI_H
#include "ch.chi.constants.h"

#include "ch.chi.potenszjals2.h"

#define potencjał woodssaxon
//TABLICA PARAMETRÓW POTENCJAŁU
constexpr double tpp[7]={1E-7,00,10,3};//rozmiar 7 oczywiście na szczęście

/////////VARIABLES////////////////////
constexpr double m=1.;
constexpr double λ=0.0002;//"spread" in position
//(DWA ALFA!!!!!!!!)
constexpr double β=1./(λ*hbar*hbar);//"spread" in momentum
//"spread"=1/(2σ^2)
constexpr double ωr=0.;//rotacja gausjana
constexpr double x0=0.;//lokalizacja w x
constexpr double p0=0.05;//inaczej p0/h czy coś

constexpr double xmin=-300;
constexpr double xmax=+300;
constexpr double pmax=0.2;//pmin: mirrored

constexpr double tmax=1E+4;



/////////////////////////////energy spectrum layout coefficients
//the following variables allow us to control how will the characteristics be laid out in the energy spectrum
constexpr double equal_energy_steps_weight=1;//uniform dE in the whole system; uniform energy intervals between characteristics
constexpr double prop_to_denisity_energy_spectrum_weight=1;//characteristics inferred from the energy distribution of the initial condition in our system (respecting the entire hamiltonian)
constexpr double system_energy_spectrum_weight=0;//similar to above, but as if the initial condition was laid out uniformly on the phase space



////////////////////////precision parameters
constexpr int spoints_p=5000;//mniej więcej points w jeden direction
constexpr int spoints_x=10000;//mniej więcej points w jeden direction for initial boundary reconassaince
constexpr int beans=30000;//musimy histogram. enerdży histogram poinc. determines the amount of characteristics but not exactly: determines the lower bound on their number which may rise in a system with many minimae
constexpr int n_energy_levels=1000;//raczej lepiej mniej niż jest beans. nie mylić z ilością charakterystyk która zależna bd także od tego ile jest minimów

constexpr int    timeres_internal=16000;//wewnętrzna rozdzielczość czasowa układu. Dlaczego nie taka sama jak rozdzielczość zewnętrzna? Bo nie dla wszystkich charakterystyk dt będzie takie samo, może się wahać +/-0.5dt0. Chyba. W ten sposób po prostu lepiej dopasujemy stan każdej charakterystyki do czasu ewolucji w którym ma być robiony rysunek
constexpr double dt0=tmax/double(timeres_internal);


/////////////////////////plotting parameters
constexpr int    timeres_1D_ext=4000;//ile punktów od początku do końca symulki dla rysunków f(t)
constexpr int    timeres_2D_ext=400;//ile punktów czasowych dla map f(a, t)
constexpr int    full_map_count=6;//ile map układu wariacie
//constexpr int    przydałoby się coś ograniczającego rozdzielczość zapewne ale 


/////int    tresmap=640;//bdb pomysł
/////int    tres1d=6000;//bdb pomysł, szypko i ładnie
/////const int    stepnode=1;
/////const int    stepnodemap=1;
///////a mury runą, runą runą
/////const double mury=0.0;//to miała być zmienna mówiąca o tym gdzie jest bariera jakby transmisji, tzn. od którego 
//////////////////////////////////////
#endif
