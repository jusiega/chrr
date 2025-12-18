//CH_BRZ_CPP    
    ////////////////////////////////////////////////////
    //skanowanie brzegu
    //////////////////////////////////////////////////
    //tu dzielimy układ jako preludium do skanowania brzegu
    //points between max and min of the system for integration
    static double energiabrzegowa[spoints_x+2*spoints_p];//lewa krawędź, zerowa (podwójnej długości), prawa krawędź
    const double dx=(xmax-xmin)/double(spoints_x-1);
    const double dp=pmax/double(spoints_p-1);//double(spoints-2);
    
    //double x=(xmax+xmin)/2.;
    //NIE idziemy po górnej krawędzi
    double x=xmin;
    double p=pmax;
    
    std::cout<<std::setprecision(10)<<"x max reconstructed "<<xmin+dx*double(spoints_x-1)<<"\n";
    std::cout<<std::setprecision(10)<<"p max reconstructed "<<dp*double(spoints_p-1)<<"\n";
    std::cout<<std::setprecision(10)<<"minus p max reconstructed "<<-dp*double(spoints_p-1)<<"\n";
    //potrzebna mi jest tablica gdzie jest każdy ten kurwa x oraz p
    //lbo może nie potrzebna
    //nie no przyda się
    static double XX[spoints_x+2*spoints_p]{};//wypełniam to w brz idąc po brz
    static double PP[spoints_x+2*spoints_p]{};//wypełniam to w brz idąc po brz
    
    double E_min=DBL_MAX;//(std::numeric_limits<double>::max)();
    double E_max=-DBL_MAX;//(std::numeric_limits<double>::lowest)();
    
    
    constexpr int nminmax=500;//więcej chyba nie będzie
    static int indeksy_minimów[nminmax]{};
    static int indeksy_maksimów[nminmax]{};
    
    FILE* filee_od_es;
    filee_od_es=fopen((folder+"/es.dat").c_str(),"w");
    

    //counter!!!!!!
    int cou=0;
    
    //verbose
    bool vvverbose=false;
    
    //potrzebne mi są jeszcze minima i maksima
    
    int indeks_tablicy_z_indeksami_minimów=0;
    for (int i=0;i<nminmax;i++) indeksy_minimów[i]=-1;
    //std::fill(indeksy_minimów, indeksy_minimów+nminmax, std::numeric_limits<double>::quiet_NaN());
    //ja pierdole ja brałem double quiet_NaN
    
    int indeks_tablicy_z_indeksami_maksimów=0;
    for (int i=0;i<nminmax;i++) indeksy_maksimów[i]=-1;
    
    
    //na pewno na brzegu s jest maksimum i nie ma ani max ani min aż do lewego rogu
    indeksy_maksimów[indeks_tablicy_z_indeksami_maksimów]=cou, indeks_tablicy_z_indeksami_maksimów++;
    
    //for (int i=0;i<spoints;i++)
    //{
    //    double E=H(x,p);
    //    if (E<E_min) E_min=E;
    //    if (E>E_max) E_max=E;
    //    fprintf(filee_od_es,"%d %e\n",cou,E);
    //    x-=dx;
    //    cou++;
    //}
    
    const double E_L=H(x,p);//lewa max energia
    for (int i=0;i<=spoints_p-2;i++)//lewy brzeg
    {
        p=pmax-(dp*double(i));
        XX[cou]=x;
        PP[cou]=p; 
        if (vvverbose) std::cout<<"x p "<<x<<" "<<p<<"\n";
        double E=H(x,p);
        energiabrzegowa[cou]=E;
        if (E<E_min) E_min=E;
        if (E>E_max) E_max=E;
        fprintf(filee_od_es,"%d %e\n",cou,E);
        cou++;
    }
    
    {
    //lewy rug. jestem w nim ale nie zapisałem energii
        p=0;
        XX[cou]=x;
        PP[cou]=p; 
        if (vvverbose) std::cout<<"x p "<<x<<" "<<p<<"\n";
        
        double E=H(x,p);
        energiabrzegowa[cou]=E;
        if (E<E_min) E_min=E;
        if (E>E_max) E_max=E;
        fprintf(filee_od_es,"%d %e\n",cou,E);
        cou++;
    }
    
    for (int i=1;i<=spoints_x-2;i++)//duł
    {
        x=xmin+dx*double(i);
        XX[cou]=x;
        PP[cou]=p; 
        //std::cout<<cou<<" "<<XX[cou]<<"\n";
        if (vvverbose) std::cout<<"x p "<<x<<" "<<p<<"\n";
        
        double E=H(x,p);        
        energiabrzegowa[cou]=E;
        //tu już potrzebuję zacząsć sprawdzić czy nie ma minimum
        //sprawdzam tak trochę w tył bo potrzebuję 3 pkt
        //trzeba się też zabezpieczyć na wypadek gdyby ekstremum wypadało dokładnie pomiędzy punktami skanowania
        //najprościej zagęścić siatke
        if(E>energiabrzegowa[cou-1]&&energiabrzegowa[cou-1]<energiabrzegowa[cou-2])
        {
            indeksy_minimów[indeks_tablicy_z_indeksami_minimów]=cou-1, indeks_tablicy_z_indeksami_minimów++;
        }
        else if(E<energiabrzegowa[cou-1]&&energiabrzegowa[cou-1]>energiabrzegowa[cou-2])
        {
            indeksy_maksimów[indeks_tablicy_z_indeksami_maksimów]=cou-1, indeks_tablicy_z_indeksami_maksimów++;
        }
        if (E<E_min) E_min=E;
        if (E>E_max) E_max=E;
        fprintf(filee_od_es,"%d %e\n",cou,E);
        cou++;
    }
    {
        //prawy rug
        x=xmax;
        XX[cou]=x;
        PP[cou]=p; 
        if (vvverbose) std::cout<<"x p "<<x<<" "<<p<<"\n";
        
        double E=H(x,p);
        energiabrzegowa[cou]=E;
        //INNY troszkę warunek ekstremum. moze być tylko minimum i będzie jeśli przedostatni był większy
        if(E<energiabrzegowa[cou-1])
        {
            indeksy_minimów[indeks_tablicy_z_indeksami_minimów]=cou-1, indeks_tablicy_z_indeksami_minimów++;
        }
        if (E<E_min) E_min=E;
        if (E>E_max) E_max=E;
        fprintf(filee_od_es,"%d %e\n",cou,E);
    }
    
    for (int i=1;i<=spoints_p-1;i++)//pierwsza iteracja tu to w zasadzie dotyczny jeszcze przemieszczenia w iksie
    {//prawy bok
        p=dp*double(i);
        cou++;
        XX[cou]=x;
        PP[cou]=p; 
        if (vvverbose) std::cout<<"x p "<<x<<" "<<p<<"\n";
        
        double E=H(x,p);
        energiabrzegowa[cou]=E;
        //na pewno nie ma ekstremów tu
        if (E<E_min) E_min=E;
        if (E>E_max) E_max=E;
        fprintf(filee_od_es,"%d %e\n",cou,E);
    }
    //aż do maks pędu ma się rozumieć. kolejne maksymum to
    indeksy_maksimów[indeks_tablicy_z_indeksami_maksimów]=cou;//, indeks_tablicy_z_indeksami_maksimów++;
    
    std::cout<<std::endl;
    //std::cout<<XX[cou]<<std::endl;
    //std::cout<<PP[cou]<<std::endl;
    
    //czy dobrze przeszedłem?
    //std::cout<<"x p "<<XX[1]<<" "<<PP[1]<<"\n";
    //std::cout<<"cou "<<cou<<"\n";
    //std::cout<<"ssss"<<XX[cou-500]<<" "<<PP[cou-500]<<"\n";
    //printf("%e %e\n",XX[1],PP[1]);
    
    const double E_P=H(x,p);//prawa max energia
    const double E_max_którą_mogę_odpowiedzialnie_symulować=std::min({E_L,E_P});
    const double E_min_const=E_min;
    //nie chcemy w układzie mieć energii niż zarówno lewa max energia ani prawa max energia
    //bo jak dojdzie do brzegu lewego czy prawego to wyjdzie poza ten obszar który zadawałem
    //jak chce na większym obszarze to se trzeba po prostu powiększyć granice (głównie pmax) w pliku wejściowym i tyle
    //tak to oznacza że może nie cały obszar prostokąta rysowanego być pokryty charakterystykami
    //ale jeśli norma warunku początkowego na siatce charakterystyk jest odtwarzana do jedności to co mnie to obchodzi? nic
    
    
    //teraz sobie wypiszmy te maksima i minima bo się martwię że się pierdolnąłem
    int ilemaksimówwariacie=0;
    {
        FILE* maksyma;
        maksyma=fopen((folder+"/max.dat").c_str(),"w");
        int i=0;
        while (indeksy_maksimów[i]!=-1) 
        {
            //std::cout<<indeksy_maksimów[i]<<"\n";
            //std::cin.get();
            fprintf(maksyma,"%d %e\n", indeksy_maksimów[i], energiabrzegowa[indeksy_maksimów[i]]);
            i++;
            ilemaksimówwariacie++;
        }
    }
    int ileminimówwariacie=0;
    {
        FILE* minima;
        minima=fopen((folder+"/min.dat").c_str(),"w");
        int i=0;
        while (indeksy_minimów[i]!=-1) 
        {
            fprintf(minima,"%d %e\n", indeksy_minimów[i], energiabrzegowa[indeksy_minimów[i]]);
            i++;
            ileminimówwariacie++;
        }
    }
    
    //for (int i=0;i<spoints;i++)
    //{
    //    double E=H(x,p);
    //    if (E<E_min) E_min=E;
    //    if (E>E_max) E_max=E;
    //    fprintf(filee_od_es,"%d %e\n",cou,E);
    //    x-=dx;
    //    cou++;
    //}
    std::cout<<"brzeg podzielono na punktów w ilości: "<<cou+1<<"\n";    
    std::cout<<"energia minimalna: "<<E_min<<"\n";
    std::cout<<"energia maksymalna: "<<E_max<<"\n";
    std::cout<<"energia lewego brzegu a max pędu: "<<E_L<<"\n";
    std::cout<<"energia prawego brzegu a max pędu: "<<E_P<<"\n";
    std::cout<<"energia maksymalna odpowiedzialna: "<<E_max_którą_mogę_odpowiedzialnie_symulować<<"\n";
    
    
    /////////////////////////////////////////////////////////////////////////////
    //end skanowanie brzegu
    /////////////////////////////////////////////////////////////////////////////
