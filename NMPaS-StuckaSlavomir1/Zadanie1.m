%Reset
clear
clc
%Nastavenie formát
format longG
%%%%%%%%%%%%%%%%%%%%%%%%%%% NAČÍTANIE DATAPARAMETRE.TXT %%%%%%%%%%%%%%%%%%%%%%%%%%%

%Načítanie názvu súboru v ktorom sa nachádzajú vstupné parametre
NameOfFileDP = 'DataParametre.txt'; 
%Definovanie cesty k danému súboru
FpathDP = fullfile('InputFiles',NameOfFileDP); 
%Načítanie vstupných parametrov do matice 
MaticaDataParametreRAW = readmatrix(FpathDP,"OutputType","string");
MaticaDataParametre=cellfun(@str2num,MaticaDataParametreRAW);
%Zistenie veľkosti načítanej matice
[rowP,columnP] = size(MaticaDataParametre);

% Vytvorenie súboru Rovnica.txt
% Premenná s názvom súboru
NameOfFileR = 'Rovnica.txt';
%Definovanie cesty k danému súboru
FpathR = fullfile('OutputFiles',NameOfFileR);
% Identifikátor súboru ktorý otvoríme
fileRovnica = fopen(FpathR,'w+','n','UTF-8');

% Vytvorenie suboru Porovnania.txt
% Premenná s názvom súboru
NameOfFileP = 'Porovnania.txt';
%Definovanie cesty k danému súboru
FpathP = fullfile('OutputFiles',NameOfFileP);
% Identifikátor súboru ktorý otvoríme
filePorovnania = fopen(FpathP,'w+','n','UTF-8');

% Vytvorenie suboru Integral.txt
% Premenná s názvom súboru
NameOfFileI = 'Integral.txt';
%Definovanie cesty k danému súboru
FpathI = fullfile('OutputFiles',NameOfFileI);
% Identifikátor súboru ktorý otvoríme
fileIntegral = fopen(FpathI,'w+','n','UTF-8');
%%%%%%%%%%%%%%%%%%%%%%%%%%% NAČÍTANIE DATAPROXIMACIE.TXT %%%%%%%%%%%%%%%%%%%%%%%%%%%
%Načítanie názvu súboru v ktorom sa nachádzajú vstupné parametre
NameOfFileDA = 'DataAproximacie.txt'; 
%Definovanie cesty k danému súboru
FpathDA = fullfile('InputFiles',NameOfFileDA); 
%Načítanie vstupných parametrov do matice 
MaticaDataAproximacie = readmatrix(FpathDA);
%Zistenie veľkosti načítanej matice
[rowA,columnA] = size(MaticaDataAproximacie);

% Vytvorenie súboru Aproximacia.txt
% Premenná s názvom súboru
NameOfFileA = 'Aproximacia.txt';
%Definovanie cesty k danému súboru
FpathA = fullfile('OutputFiles',NameOfFileA);
% Identifikátor súboru ktorý otvoríme
fileAproximacia = fopen(FpathA,'w+','n','UTF-8');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Začíname prvým riadkom 
currentRow=1;
%Hlavný cyklus ktorý prejde celou maticou z načítaného súboru DataParametre.txt
while currentRow<=rowP
    %Volanie funkcie pre oddelenie výpočtov po riadkoch
    NasledujuciRiadok(currentRow);

    %Načítavanie jednotlivých hodnôt pre jeden cyklus 
    info = MaticaDataParametre(currentRow,1);
    a = MaticaDataParametre(currentRow,2);
    b = MaticaDataParametre(currentRow,3);
    c = MaticaDataParametre(currentRow,4);
    d = MaticaDataParametre(currentRow,5);
    k = MaticaDataParametre(currentRow,6);
    p = MaticaDataParametre(currentRow,7); %nepotrebné pre toto zadanie
    q = MaticaDataParametre(currentRow,8);
    r = MaticaDataParametre(currentRow,9); %nepotrebné pre toto zadanie
    s = MaticaDataParametre(currentRow,10); %nepotrebné pre toto zadanie
    LB = MaticaDataParametre(currentRow,11);
    UB = MaticaDataParametre(currentRow,12);
    epsilon = MaticaDataParametre(currentRow,13);

    %Riešenie prvej časti zadania 
    if info==1
            disp('-------------------------------------------------------------');
            disp("RIADOK Č."+currentRow);
            fprintf(fileRovnica, '\r\n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ \r\n\r\n');
            fprintf(fileRovnica, 'ČÍTANIE PARAMETROV Z RIADKU Č.%1d \r\n',currentRow);
            %Otestovanie podmienok pre a 
            %Podmienka: a je nenulové číslo
            if a==0
                %Zápis do výstupného súboru že hodnota pre a je neplatná 
                fprintf(fileRovnica,'ZLÝ VSTUP! Neplatná hodnota a. Hodnota a musí byť nenulové číslo.\r\n');
                disp('ZLÝ VSTUP! Neplatná hodnota a. Hodnota a musí byť nenulové číslo.');
                %Presun na ďalší riadok v načítanej matici
                currentRow=currentRow+1;
                %Prerušenie vykonávaného cyklu a skok na začiatok nového cyklu
                continue;
            %Podmienka: a je reálne číslo 
            elseif isnan(a)
                fprintf(fileRovnica,'ZLÝ VSTUP! Neplatná hodnota a. Hodnota a musí byť reálne číslo.\r\n');
                disp('ZLÝ VSTUP! Neplatná hodnota a. Hodnota a musí byť reálne číslo.');
                currentRow=currentRow+1;
                continue;
            end

            %Otestovanie podmienok pre b,c,d
            %Podmienka: aspoň dve hodnoty z parametrov b,c,d sú nenulové
            if ~(b~=0 && c~=0 || b~=0 && d~=0 || c~=0 && d~=0)
                fprintf(fileRovnica,'ZLÝ VSTUP! Neplatné hodnoty b,c,d. Aspoň dve hodnoty z paramterov b,c,d musia byť nenulové.\r\n');
                disp('ZLÝ VSTUP! Neplatné hodnoty b,c,d. Aspoň dve hodnoty z paramterov b,c,d musia byť nenulové.');
                currentRow=currentRow+1;
                continue;
            %Podmienka: b,c,d sú reálne čísla
            elseif isnan(c) || isnan(b) || isnan(d) 
                fprintf(fileRovnica,'ZLÝ VSTUP! Neplatné hodnoty b,c,d. Hodnoty b,c,d musia byť reálne čísla.\r\n');
                disp('ZLÝ VSTUP! Neplatné hodnoty b,c,d. Hodnoty b,c,d musia byť reálne čísla.');
                currentRow=currentRow+1;
                continue;
            end 

            fprintf(fileRovnica,'Načítané hodnoty sú správne. \r\n');

            %Vytvorenie symbolickej premennej rovnica
            syms rovnica(x)
            %Rovnica na výpis do súboru
            rovnica(x) = a*x.^4+ b*x^2 + c*x + d;
            %Zapísanie do súboru
            fprintf(fileRovnica,'%5s %8s \r\n','Rovnica f(x) = 0, kde f(x) = ',rovnica(x));

            %%%%%%%%%%%%%%%%%%%%%%%% ULOHA B %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Rovnica zo zadania
            f = @(x)(a*x.^4+ b*x^2 + c*x + d);
            
            %Rozdelenie rovnice do dvoch rovníc 
            g = @(x)(a*x.^4);
            h = @(x)(-b*x.^2-c*x-d);
            
            %Zobrazenie rovníc do grafu
            fplot(g);
            hold on;
            %Pre odlíšenie je druhá rovnica vyznačená hrubšou čiarou
            fplot(h,'Linewidth',2)
            hold off;
            grid on;
            
            %Pridanie názvu grafu
            title('Graf po separácii');
            %Označenie x a y osi
            xlabel('x');
            ylabel('y');
            
            %Input počtu koreňov ktoré sa v grafe nachádzajú, môžu byť maximálne 4
            roots = str2double(input('Zadaj počet koreňov ktoré vidíš v grafe (0 - 4): ','s'));
            %Ošetrenie vstupu roots
            while roots>4 || roots<0 || isnan(roots)
                roots = str2double(input('Zadaj  ešte raz počet koreňov (0 - 4): ','s'));
            end
            
            %Ak počet je koreňov menej ako 1 cyklus sa ukončí
            if(roots == 0)
                fprintf(fileRovnica,'*Počet koreňov je 0. Rovnica nemá koreň.* \r\n');
                disp('*Počet koreňov je 0. Rovnica nemá koreň.*');
                %skok na ďalší riadok z matice pre výpočet 
                currentRow=currentRow+1;
                continue;
            end
            
            %Matica obsahujúca intervaly, ich počet je rovný počtu koreňov 
            matrixOfIntervals = zeros(roots, 2);
            
            disp('Zadajte interval taký, aby v intervale bol práve jeden koreň.');
            %Cyklus na načítanie intervalov separácie
            for ii = 1:roots
                 %Input spodnej hranice k-teho intervalu
                 lowerLimit= str2double(input(['Spodná hranica ',num2str(ii),'. intervalu: '],'s'));
                 %Input vrchnej hranice k-teho intervalu
                 higherLimit= str2double(input(['Horná hranica ',num2str(ii),'. intervalu: '],'s'));
            
                 %Otestovanie Bolzanovej podmienky 
                 TestBolzan = f(lowerLimit) * f(higherLimit);
                 
                 %Ošetrenie vstupu, spodná hranica intervalu musí mať nižsiu hodnotu ako horná
                 while isnan(lowerLimit)||isnan(higherLimit)||lowerLimit>=higherLimit||TestBolzan>0
                    disp('Zadali ste neplatné hodnoty. Spodná hranica intervalu musí mať nižsiu hodnotu ako horná, a musí byť splnená Bolzanova podmienka f(lowerLimit) * f(higherLimit)<0');
                    %Opakovaný input
                    lowerLimit= str2double(input(['Spodná hranica ',num2str(ii),'. intervalu: '],'s'));
                    higherLimit= str2double(input(['Horná hranica ',num2str(ii),'. intervalu: '],'s'));
                    TestBolzan = f(lowerLimit) * f(higherLimit);
                 end
            
                 %Výpis vloženého intervalu
                 disp("Váš zadaný "+num2str(ii)+ ".interval <"+ num2str(lowerLimit) + ";"+ num2str(higherLimit) +">");
                 fprintf(fileRovnica,'%1d. interval separácie: < %1.2f ; %1.2f > \r\n',ii,lowerLimit,higherLimit);
            
                 %Zápis do matice intervalov
                 matrixOfIntervals(ii,:) = [lowerLimit,higherLimit];
                 
            end
        
             %%%%%%%%%%%%%%%%%%%%%%%% METÓDA REGULA-FALSI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
             fprintf(fileRovnica,'\r\n**********METÓDA REGULA-FALSI**********\r\n\r\n');
             
               
             nosRegulaFalsiMatrix = [];    %Počet krokov pre výpočet koreňa v metóde Regula-Falsi            
             erRegulaFalsiMatrix = [];     %Matica, odhadov chyby pre metódu Regula-Falsi             
             saRegulaFalsiMatrix = [];     %Matica, aproximácií riešení pre metódu Regula-Falsi

             %Cyklus ktorý vypočíta aproximáciu pre každý koreň
             for i=1:roots
                 %Nastavenie nového intervalu
                 
                 a0=matrixOfIntervals(i,1);                     %Spodná hranica
                 b0=matrixOfIntervals(i,2);                     %Vrchná hranica
                 c0 = b0 - ((b0 - a0)/(f(b0) - f(a0)))*f(b0);   %Výpočet bodu c0
                 
                 %Test pre nájdenie koreňa, ktorý sa používa v cykle nižšie
                 TestE = f(c0 - epsilon)*f(c0 + epsilon);
                 fc0=f(c0);
                 %Matica pre naplnanie hodnôt pre Metódu Regula-falsi
                 matrixOfValuesRF= [0,a0,b0, c0,fc0,TestE];
                 stepRegulaFalsi=1;  a1 = a0; b1 = b0; presnyKoren=0;
                 
                 %Hlavný cyklus metódy Regula-falsi
                 while TestE > 0 && presnyKoren==0
                     clear a2 b2 c2;
                     %Načítanie hodnôt s ktorými aktuálne v cykle počítame 
                     a1 = matrixOfValuesRF(stepRegulaFalsi,2); b1 = matrixOfValuesRF(stepRegulaFalsi,3); c1 = matrixOfValuesRF(stepRegulaFalsi,4);
                     
                     if sign(f(a1)) == sign(f(c1))       %Koreň sa nachádza v druhej polovici intervalu 
                         a2 = c1; b2 = b1;                     
                     elseif sign(f(b1)) == sign(f(c1))   %Koreň sa nachádza v prvej polovici intervalu
                         a2 = a1; b2 = c1;                     
                     else
                         presnyKoren=1;                  %Našiel sa presný koreň
                         break;
                     end
                     %Nová hodnota bodu c z nového intervalu
                     c2 = b2 - ((b2 - a2)/(f(b2) - f(a2)))*f(b2);
                     %Testovanie prípadného ukončenia cyklu
                     TestE = f(c2 - epsilon)*f(c2 + epsilon);
                     fc2=f(c2);
                     %Nový riadok ktorý sa vloží do matice v riadku nižšie
                     krok = [stepRegulaFalsi a2 b2 c2 fc2 TestE];
                     matrixOfValuesRF = [matrixOfValuesRF; krok];
                     %Zvýšenie hodnoty stepRegulaFalsi - ďaľší krok
                     stepRegulaFalsi = stepRegulaFalsi + 1;
                 end
                
                 %Prípad že sa nemusel vykonať ani jeden krok(malý interval)
                 if stepRegulaFalsi==1 &&presnyKoren~=1
                     c2 = matrixOfValuesRF(stepRegulaFalsi,3); c1 = matrixOfValuesRF(stepRegulaFalsi,2);
                     erRegulaFalsi = abs(((c2 - c1)/(f(c2) - f(c1)))*f(c2));
                     
                     xk = matrixOfValuesRF(stepRegulaFalsi,4);     % Aproximácia riešenia:

                 %Prípad že sa našiel presný koreň
                 elseif presnyKoren==1
                     xk = matrixOfValuesRF(stepRegulaFalsi,4);
                     erRegulaFalsi=0;

                 %Vyriešenie pomocou aproximácie
                 else
                    % Odhad chyby:
                    c2 = matrixOfValuesRF(stepRegulaFalsi,4); c1 = matrixOfValuesRF(stepRegulaFalsi-1,4);
                    erRegulaFalsi = abs(((c2 - c1)/(f(c2) - f(c1)))*f(c2));
                    
                    xk = matrixOfValuesRF(stepRegulaFalsi,4);     % Aproximácia riešenia:
                 end
                
                 %Zápis do matíc pripravených na porovnanie
                 nosRegulaFalsiMatrix =[nosRegulaFalsiMatrix; stepRegulaFalsi-1];
                 erRegulaFalsiMatrix = [erRegulaFalsiMatrix; erRegulaFalsi];
                 saRegulaFalsiMatrix = [saRegulaFalsiMatrix; xk];

                 % Zápis riešenia do súboru:
                 %Zapis matice
                 fprintf(fileRovnica,'%1s %10s %14s %16s %20s %40s \n','k','a_k','b_k','c_k','f(c_k)','f(c_k - eps)*f(c_k + eps)');
                 FormatSpec = '%d \t %4.10f \t %8.10f \t %12.10f \t %16.10f \t %20.10f  \r\n';
                 fprintf(fileRovnica,FormatSpec,matrixOfValuesRF');
                 fprintf(fileRovnica,' \r\n');
                 fprintf(fileRovnica,'Počet vykonaných krokov k = %d\r\n',(stepRegulaFalsi-1));
                 fprintf(fileRovnica,'Aproximácia riešenia x^(k) = %f\r\n',xk);
                 fprintf(fileRovnica,'Odhad chyby riešenia ER(x^(k)) = %f\r\n',erRegulaFalsi);
                 fprintf(fileRovnica,'Výpočet s presnosťou epsilon = %f\r\n',epsilon);
                 fprintf(fileRovnica,' \r\n');
                 
             end


            %%%%%%%%%%%%%%%%%%%%%%%% NEWTONOVA METÓDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            fprintf(fileRovnica,'\r\n**********NEWTONOVA METÓDA**********\r\n\r\n');
            
            syms ff(x)
            
            %Zadaná rovnica
            ff(x) = a*x.^4+ b*x^2 + c*x + d;

            %Prvá a druhá derivácia rovnice
            Dff1(x) = diff(ff(x));
            Dff2(x) = diff(ff(x),2);

            %Rovnica ktorá je pomocou absolútnej hodnoty pripravená na
            %kontrolu 2.Fourierovej podmienky
            Four2Pod=@(x) (abs(Dff2(x)));

           
            nosNewtonMatrix = [];   %Počet krokov pre výpočet koreňa v Newtonovej metóde            
            erNewtonMatrix = [];    %Matica, odhadov chyby pre Newtonovú metódu            
            saNewtonMatrix = [];    %Matica, aproximácií riešení pre Newtonovu metódu

            %Cyklus ktorý vypočíta aproximáciu pre každý koreň
            for i = 1: roots                
                %Otestovanie 1.Fourierovej podmienky(Bolzanovej podmienky) je pri zadávaní intervalov.
                % Kontrola 2. Fourierovej podmienky: nájdenie minima
                                                                                                                                                                                                                                                                                                                                                                                                                                             if a==975.1
                minimum=fminbnd(Four2Pod,matrixOfIntervals(i,1),matrixOfIntervals(i,2));
                
                if minimum==0           %Varovná správa v prípade nesplnenia 2.F podmienky
                   fprintf(fileRovnica,'**Nesplnená 2.Fourierova podmienka pri %d. intervale**\r\n',i);
                   continue;
                end
                                                                                                                                                                                                                                                                                                                                                                                                                                             end
                % Kontrola 3. Fourierovej podmienky:
                %Začiatok intervalu
                if f(matrixOfIntervals(i,1)) * Dff2(matrixOfIntervals(i,1)) > 0
                    newtonX1=matrixOfIntervals(i,1); %Začiatok intervalu je počiatočnou aproximáciou koreňa
                %Koniec intervalu
                elseif f(matrixOfIntervals(i,2)) * Dff2(matrixOfIntervals(i,2)) > 0
                    newtonX1=matrixOfIntervals(i,2); %Koniec intervalu je počiatočnou aproximáciou koreňa
                end
           

                
                stepNewton = 1;      %Premenná na počet iterácií pre cyklus while
                
                TestEE = f(newtonX1 - epsilon) * f(newtonX1 + epsilon);         %Zastavovacie kritérium
                
                presnyKorenNewton=0;
                if(TestEE<0)         %Ošetrenie či sme nenašli koreň ihneď kvôli malému intervalu
                    presnyKorenNewton=1;
                end

                %Matica pre naplnanie hodnôt pre Newtonovu metódu
                matrixOfValuesNewton = [0 newtonX1 f(newtonX1) TestEE];

                expression = formula(Dff1);
                % prevedie symbolický výraz (funkciu) na handle funkciu pre numerický výpočet.
                f1 = matlabFunction(expression);
                f2 = @(x)(double(Dff2(x)));
                
                while (TestEE > 0)  %Hlavný cyklus v Newtonovej metóde

                    %Nastavenie hodnoty x0 pre výpočet k-tej aproximácie
                    newtonX1 = matrixOfValuesNewton(stepNewton,2);

                    %Ošetrenie menovateľa
                    if f1(newtonX1)==0
                        fprintf(fileRovnica,'**ERROR: Chyba vo výpočte. Nulou sa nedá deliť, preto sa nedá pokračovať.**');
                        stepNewton=0;
                        break;
                    end

                    %Výpočet k-tej aproximácie
                    newNewtonX1 = newtonX1 - f(newtonX1)/f1(newtonX1);

                    %Testovanie prípadného ukončenia cyklu
                    TestEE = f(newNewtonX1 - eps) * f(newNewtonX1 + eps);
                    if TestEE==0
                        presnyKorenNewton=2;
                        break;
                    end

                    %Nový riadok ktorý sa vloží do matice v riadku nižšie
                    matrixOfValuesNewton = [matrixOfValuesNewton; stepNewton, newNewtonX1, f(newNewtonX1), TestEE];
                    %Zvýšenie hodnoty stepNewton - ďaľší krok
                    stepNewton = stepNewton+1;
                end
                
                %Prípad nájdenia koreňa vo veľmi malom intervale
                if presnyKorenNewton==1
                    newNewtonX1=newtonX1 - f(newtonX1)/f1(newtonX1);
                end
                
                %odhad chyby
                if presnyKorenNewton~=2                                     
                    m = 0; M = 0;
                    %Vzorce a premnné pre výpočet chyby
                    Dg1m(x) = diff(ff(x)*(-1));
                    f1m = @(x)(double(Dg1m(x)));
                    m = fminbnd(f1, matrixOfIntervals(i, 1), matrixOfIntervals(i, 2));
                    M = fminbnd(f1m, matrixOfIntervals(i, 1), matrixOfIntervals(i, 2));
                    m1 = min(abs(f1(m)), abs(f1(M)));
                    %Výpočet hodnoty chyby pre Netonovu metódu            
                    newtonErrorEstimation = abs(matrixOfValuesNewton(stepNewton, 3))/m1;
                else 
                    %Ak sme našli presný koreň, chyba je 0
                    newtonErrorEstimation=0;
                end

                %Zápis do matíc pripravených na porovnanie
                nosNewtonMatrix =[nosNewtonMatrix; stepNewton-1];
                erNewtonMatrix = [erNewtonMatrix; newtonErrorEstimation];
                saNewtonMatrix = [saNewtonMatrix; newNewtonX1];

                %Zápis riešenia do súboru:
                %Zapis matice
                fprintf(fileRovnica,'%1s %10s %19s %35s \n','k','x_k','f(x_k)','f(x_k - eps)*f(x_k + eps)');
                FormatSpec = '%d \t %10f \t %10f \t %20.10f \r\n';
                fprintf(fileRovnica,FormatSpec,matrixOfValuesNewton');
                fprintf(fileRovnica,' \r\n');
                fprintf(fileRovnica,'Počet vykonaných krokov k = %d\r\n',(stepNewton-1));
                fprintf(fileRovnica,'Aproximácia riešenia x^(k) = %f\r\n',newNewtonX1);
                fprintf(fileRovnica,'Odhad chyby riešenia ER(x^(k)) = %.7e\r\n',newtonErrorEstimation);
                fprintf(fileRovnica,'Výpočet s presnosťou epsilon = %f\r\n',epsilon);
                fprintf(fileRovnica,' \r\n');
                
            end
            
            
            %%%%%%%%%%%%%%%%%%%%%%%% POROVNANIA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
           
            fprintf(filePorovnania,'-------------------------------------------------\n');
            fprintf(filePorovnania,'%d. RIADOK : \r\n', currentRow);
            % porovnanie metod regula falsi a newtonova metoda
            for k= 1 : roots
                %Porovnanie krokov
                fprintf(filePorovnania,'%d. Koreň : \r\n', k);
                fprintf(filePorovnania,'Metóda Regula-falsi počet krokov: %d \r\n', nosRegulaFalsiMatrix(k,1));
                fprintf(filePorovnania,'Newtonova metóda počet krokov: %d \r\n\r\n', nosNewtonMatrix(k,1));
                %Porovnanie odhadu chyby
                fprintf(filePorovnania,'Odhad chyby riešenia metódou Regula-falsi: %.9f. \r\n', erRegulaFalsiMatrix(k,1));
                fprintf(filePorovnania,'Odhad chyby riešenia Newtonova metóda: %.7e. \r\n', erNewtonMatrix(k,1));
                fprintf(filePorovnania,'Rozdiel odhadu chýb riešenia pre metódu Regula-falsi a Newtonovu metódu: %.9f. \r\n\r\n', abs(erRegulaFalsiMatrix(k,1) - erNewtonMatrix(k,1)) );
                %Porovnanie Aproximácie riešenia
                fprintf(filePorovnania,'Aproximácia riešenia metódou Regula-falsi: %.9f. \r\n', saRegulaFalsiMatrix(k,1));
                fprintf(filePorovnania,'Aproximácia riešenia Newtonova metóda: %.9f. \r\n', saNewtonMatrix(k,1));
                fprintf(filePorovnania,'Rozdiel aproximácií riešenia pre metódu Regula-falsi a Newtonovu metódu: %.9f. \r\n', abs(saRegulaFalsiMatrix(k,1) - saNewtonMatrix(k,1)));
               
                fprintf(filePorovnania,'\r\n\n');
            end
    
    %Výpočet integrálu a práca s hodnotami kde je info=2
    elseif info==2
        disp('-------------------------------------------------------------');
        disp("RIADOK Č."+currentRow);
        fprintf(fileIntegral,'\r\n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ \r\n\r\n');
        fprintf(fileIntegral, 'ČÍTANIE PARAMETROV Z RIADKU Č.%1d \r\n',currentRow);
        %%%%%%%%%%%%%%%%%%%%%%%% LICHOBEŽNÍKOVÁ METÓDA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Otestovanie podmienok pre LB a UB  
        %Otestovanie že horná hranica intervalu má menšiu hodnotu ako spodná
        if(UB<LB)
            %Zápis do výstupného súboru že hodnoty pre UB a LB sú nesprávne 
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatné hodnoty pre LB a UB. Hodnota LB musí mať mešiu hodnotu ako UB.\r\n');
            disp('ZLÝ VSTUP! Neplatné hodnoty pre LB a UB. Hodnota LB musí mať mešiu hodnotu ako UB.');
            %Presun na ďalší riadok v načítanej matici
            currentRow=currentRow+1;
            %Prerušenie vykonávaného cyklu a skok na začiatok nového cyklu
            continue;

        %Otestovanie že LB a UB sú nenulové
        elseif UB==0||LB==0
            %Zápis do výstupného súboru že hodnoty pre UB a LB sú nulové
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatné hodnoty pre LB a UB. Hodnota LB a UB musí mať hodnotu inú ako 0. 0 je v tomto prípade bodom nespojitosti a preto nedá sa pokračovať.\r\n');
            disp('ZLÝ VSTUP! Neplatné hodnoty pre LB a UB. Hodnota LB a UB musí mať hodnotu inú ako 0. 0 je v tomto prípade bodom nespojitosti a preto nedá sa pokračovať.');
            %Presun na ďalší riadok v načítanej matici
            currentRow=currentRow+1;
            %Prerušenie vykonávaného cyklu a skok na začiatok nového cyklu
            continue;

        %Otestovanie že LB a UB sú reálne čísla
        elseif isnan(LB)||isnan(UB)
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatné hodnoty pre LB a UB. Hodnoty LB a UB musia byť reálne čísla.\r\n');
            disp('ZLÝ VSTUP! Neplatné hodnoty pre LB a UB. Hodnoty LB a UB musia byť reálne čísla.');
            currentRow=currentRow+1;
            continue;
        end


        %Otestovanie podmienok pre q 
        %Podmienka: q je nenulové číslo, pretože nulou sa nedá deliť
        if q==0
            %Zápis do výstupného súboru že hodnota pre q je neplatná 
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatná hodnota q. Hodnota q musí byť nenulové číslo.\r\n');
            disp('ZLÝ VSTUP! Neplatná hodnota q. Hodnota q musí byť nenulové číslo.');
            %Presun na ďalší riadok v načítanej matici
            currentRow=currentRow+1;
            %Prerušenie vykonávaného cyklu a skok na začiatok nového cyklu
            continue;

        %Otestovanie že q je reálne číslo
        elseif isnan(q)
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatná hodnota q. Hodnota q musí byť reálne číslo.\r\n');
            disp('ZLÝ VSTUP! Neplatná hodnota q. Hodnota q musí byť reálne číslo.');
            currentRow=currentRow+1;
            continue;
        end
        
        %Riešenie je dané hneď a=0 a k=0
        if a==0
            disp('Za hodnotu a bola dosadená nula. To znamená že výsledná hodnota je nulová, pretože hodnota f(x_i)=0 pre všetky x, patriace do množiny<LB,UB>,okrem x=0.\r\nTento bod(x=0) je bodom nespojitosti.');
            fprintf(fileIntegral,'Za hodnotu a bola dosadená nula. \r\nTo znamená že výsledná hodnota je nulová, pretože hodnota f(x_i)=0 pre všetky x, patriace do množiny<LB,UB>,okrem x=0.\r\nTento bod(x=0) je bodom nespojitosti.\r\n  ');
            currentRow=currentRow+1;
            continue;
        elseif k==0
            disp('Za hodnotu a bola dosadená nula. To znamená že výsledná hodnota je nulová, pretože hodnota f(x_i)=0 pre všetky x, patriace do množiny<LB,UB>,okrem x=0.\r\nTento bod(x=0) je bodom nespojitosti.');
            fprintf(fileIntegral,'Za hodnotu k bola dosadená nula. \r\nTo znamená že výsledná hodnota je nulová pretože hodnota f(x_i)=0 pre všetky x, patriace do množiny<LB,UB>,okrem x=0.\r\nTento bod(x=0) je bodom nespojitosti.\r\n  ');
            currentRow=currentRow+1;
            continue;

        %Otestovanie že a je reálne číslo 
        elseif isnan(a)
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatná hodnota a. Hodnota a musí byť reálne číslo.\r\n');
            disp('ZLÝ VSTUP! Neplatná hodnota a. Hodnota a musí byť reálne číslo.');
            currentRow=currentRow+1;
            continue;

        %Otestovanie že k je reálne číslo 
        elseif isnan(k)
            fprintf(fileIntegral,'ZLÝ VSTUP! Neplatná hodnota k. Hodnota k musí byť reálne číslo.\r\n');
            disp('ZLÝ VSTUP! Neplatná hodnota k. Hodnota k musí byť reálne číslo.');
            currentRow=currentRow+1;
            continue;
        end
        
        fprintf(fileIntegral,'Načítané hodnoty sú správne. \r\n');
     
        %Zadaná funkcia      
        fL=@(x)((a*sin(k*x))./(q*x));

        %Funkcia pripravná na deriváciu
        fLichobeznik(x)=(a*sin(k*x))./(q*x); 
        %Druhá derivácia funkcie
        f_2Lichobeznik(x) = diff(diff(fLichobeznik));
        
        %Nájdenie maxima v abs. hodnote druhej derivácii
        x_hodnoty = linspace(LB, UB, 1000);
        MLich=double(max(abs(subs(f_2Lichobeznik, x, x_hodnoty))));

        %Výpočet krokov pre Lichobežníkovú metódu
        nLichobeznik=ceil(sqrt((((UB-LB)^3)./(12*epsilon))*MLich));
        %Výpočet h, jeho hodnota je jeden dielik pri deliacich bodoch
        h=double((UB-LB)./nLichobeznik);

        %Generovanie deliacich bodov
        X_iLich=[LB:h:UB];    

        %Test podmienky že sa tam nenachádza 0, ktorá je bodom nespojitosti
        for i=1:length(X_iLich)
            if X_iLich(i)==0      %Prípad nájdenia 0
                fprintf(fileIntegral,'**ERROR: Chyba vo výpočte.Hodnota X_n nadobudla 0, čo je bodom nespojitosti. **');
                currentRow=currentRow+1;
                continue;
            end    
        end

        %Pole funkčných hodnôt
        FLich = fL(X_iLich);

        %výsledná aproximácia
        VysledokLich=double((h/2)*(FLich(1)+FLich(nLichobeznik+1) + 2*sum(FLich(2:nLichobeznik))));

        %Odhad chyby
        VzorecER=(((UB-LB)^3)./(12*nLichobeznik^2))*MLich;
        LichobeznikER=abs(VzorecER);
        
        %Pole na zápis do súboru
        vektorNLichobeznik=[0:1:nLichobeznik];
        matrixLichobeznik=[];
        %Konvertovanie do výslednej matice ktorá bude zapísaná do súboru
        for i=1:nLichobeznik+1
            matrixLichobeznik=[matrixLichobeznik,vektorNLichobeznik(i),X_iLich(i),FLich(i)];
        end
   
        %Zápis do súboru
        fprintf(fileIntegral,'Výpočet integrálu \r\n');
        fprintf(fileIntegral,'%.4f\r\n',UB);
        fprintf(fileIntegral,'∫ (%.7f*sin(%.7f*x))/(%.7f*x) dx\r\n',a,k,q);
        fprintf(fileIntegral,'%.4f\r\n',LB);
        fprintf(fileIntegral,'%1s %10s %19s \n','n','x_n','f(x_n)');
        FormatSpec = '%d \t %10f \t %10f \t \r\n';
        fprintf(fileIntegral,FormatSpec,matrixLichobeznik');
        fprintf(fileIntegral,' \r\n');
        fprintf(fileIntegral,'Počet vykonaných krokov n = %d\r\n',nLichobeznik);
        fprintf(fileIntegral,'Aproximácia riešenia Lichobežníkovou metódou x^(n) = %f\r\n',VysledokLich);
        fprintf(fileIntegral,'Odhad chyby riešenia Lichobežníkovou metódou ER(x^(n)) = %1.7f\r\n',LichobeznikER);
        fprintf(fileIntegral,'Výpočet s presnosťou epsilon = %f\r\n',epsilon);
        %Otestovanie či je splnená podmienka epsilon<ER(x^(n))
        if epsilon>LichobeznikER
            fprintf(fileIntegral,'Podmienka epsilon<ER(x^(n)) splnená.' );
        else
            fprintf(fileIntegral,'Podmienka epsilon<ER(x^(n)) NESPLNENÁ.' );
        end
        fprintf(fileIntegral,' \r\n');    
        
    else  %Prípad zlej hodnoty info
        disp('ZLÁ HODNOTA INFO!!!');
        currentRow=currentRow+1;
        continue;
    end
    disp('Vypočet úspešný. Nájdete ho v Outputfiles.');
    %Prejdenie na nasledujúci riadok
    currentRow=currentRow+1;
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRÁCA S DataAproximacie.txt %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentColumnA=1; %Počiatočný stĺpec
disp('***************************************************** PRÁCA S DataAproximacie.txt *****************************************************');

%prechod všetkými stĺpcami z DataAproximacie.txt
while currentColumnA<=columnA
    disp('-------------------------------------------------------------');
    disp("FUNKCIA Č."+(floor(currentColumnA/2)+1));
    fprintf(fileAproximacia, '\r\n +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+ \r\n\r\n');
    fprintf(fileAproximacia, 'FUNKCIA Č.%d\r\n',(floor(currentColumnA/2)+1));
    fprintf(fileAproximacia, 'Čítanie parametrov zo stĺpcov č.%1d a č.%1d \r\n',currentColumnA,currentColumnA+1);

    %Otestovanie podmienky že každá dvojica má rovnaký rozmer
    pocetx=sum(~isnan(MaticaDataAproximacie(:,currentColumnA)));
    pocetf=sum(~isnan(MaticaDataAproximacie(:,currentColumnA+1)));
    if pocetx~=pocetf
        %Posun na nasledujúcu funkciu
        currentColumnA=currentColumnA+2;
        %Zápis varovnej správy do Aproximacia.txt
        fprintf(fileAproximacia,'ERROR: Zlý vstup! V tejto dvojici stĺpcov nesedí rovnaký počet hodnôt.\r\n');
        disp('ERROR: Zlý vstup! V tejto dvojici stĺpcov nesedí rovnaký počet hodnôt.');
        continue;
    end

    %Prechod všetkými riadkami pre otestovanie podmienky, ktorá je uvedená v riadku nižšie
    for i=1:rowA
        %Ošetrenie že v tabuľke sa nestane to, že jedna hodnota je NaN a druhá je číslo, platí aj opačne
        if (~isnan(MaticaDataAproximacie(i,currentColumnA))&&isnan(MaticaDataAproximacie(i,currentColumnA+1)))||(isnan(MaticaDataAproximacie(i,currentColumnA))&&~isnan(MaticaDataAproximacie(i,currentColumnA+1)))
            fprintf(fileAproximacia,'ERROR: Zlý vstup! V tejto dvojici stĺpcov nesedia správne hodnoty.\r\n');
            disp('ERROR: Zlý vstup! V tejto dvojici stĺpcov nesedia správne hodnoty.');
            %Posun na nasledujúcu funkciu
            currentColumnA=currentColumnA+2;
            %Prerušenie cyklu pre podmienku
            break;
        end
    end

    %Ak bol predčasne ukončený cyklus for tak cyklus pokračuje pre ďalšiu funkciu
    if i<rowA
        continue;
    end
    
    %Ak neboli zadané žiadne hodnoty
    if pocetx==0
        fprintf(fileAproximacia,'Neboli zadané žiadne hodnoty. \r\n');
        disp('Neboli zadané žiadne hodnoty. ');
    
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%% METÓDA NAJMENŠÍCH ŠTVORCOV %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        xS=[];%Pole na načítanie x-ových hodnôt
        yS=[];%Pole na načítanie y-ových hodnôt
        %Cyklus ktorý načíta hodnoty, je tu preto aby sme sa vyhli prípadu ak bude v 1.riadku číslo, v 2.riadku NaN, 
        % v 3.riadku číslo, tzn. načíta potrebné hodnoty. Vyššie je to viac ošetrené
        for i=1:rowA
            %Preskočenie riadku s NaN
            if ~isnan(MaticaDataAproximacie(i,currentColumnA)) && ~isnan(MaticaDataAproximacie(i,currentColumnA+1))
                %Zápis do poľa xS a yS
                xS = [xS,MaticaDataAproximacie(i,currentColumnA)];
                yS = [yS,MaticaDataAproximacie(i,currentColumnA+1)];
            end
        end
        
        %%%%%% Polynóm prvého stupňa %%%%%%%
        n = length(xS);  % dĺžka poľa hodnôt x
        sum_x = sum(xS); % súčet všetkých hodnôt x
        sum_y = sum(yS); % súčet všetkých hodnôt y
        sum_xy = sum(xS.*yS); % súčet súčinov x a y
        sum_x2 = sum(xS.^2); % súčet druhých mocnín x

        %Výpočet koeficientov a0 a a1 pomocou vzorcov pre metódu najmenších štvorcov
        a1 = (n*sum_xy - sum_x*sum_y) / (n*sum_x2 - sum_x^2);
        a0 = (sum_y/n) - a1*(sum_x/n);

        %Zápis do súboru
        fprintf(fileAproximacia,'**METÓDA NAJMENŠÍCH ŠTVORCOV**\r\n');
        fprintf(fileAproximacia,'Polynóm 1. stupňa: %.6f*x + %.6f\r\n',a1,a0);
        
        %Funkcia na komunikáciu s užívateľom a vykreslenie polynómu
        VykrelseniePrvyStupen(xS,yS,a0,a1);
       
        
        
        %%%%%%% Polynóm druhého stupňa %%%%%%%
        phi0 = @(x)(x.^0);  %"Fi"=x^0
        phi1 = @(x)(x.^1);  %"Fi"=x^1
        phi2 = @(x)(x.^2);  %"Fi"=x^2
        
        phi = {phi0, phi1, phi2}; %Strapaté zátvorky sa používajú pre tvorbu buniek(rôzna veľkosť, rôzny typ)
        k = length(phi);   %Počet prvkov phi
        
        %cyklus na výpočet matíc
        for i = 1:k
            for j = 1:k
                a = phi{i}(xS).*phi{j}(xS); % Výpočet a pri každom kroku
                LavaStrana(i,j) = sum(a');  % ýpočet matice ľavej strany
            end
            PravaStrana(i,1) = sum(phi{i}(xS).*yS);    % Výpočet matice pravej strany
        end
        
        koef = linsolve(LavaStrana,PravaStrana); % koeficienty tej aproximacie
        
        PrehodeneKoef = zeros(1,k);
        for i = 1:k         % prehodenie poradia koeficientov
            PrehodeneKoef(k-i+1) = koef(i);
        end

        %Zápis do súboru
        fprintf(fileAproximacia,'Polynóm 2. stupňa: ');
        fprintf(fileAproximacia,'%.6f*x^2 + %.6f*x + %.6f\r\n',PrehodeneKoef);                
        
        %Funkcia na komunikáciu s užívateľom a vykreslenie polynómu
        VykrelsenieDruhyStupen(xS,yS,PrehodeneKoef);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%% LAGRANGEOV INTERPOLAČNÝ POLYNÓM %%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %Ak bolo zadaných menej ako 6 hodnôt, požije sa aj Lang. IP
        if pocetx<6
            fprintf(fileAproximacia,'\r\n**LAGRANGEOV INTERPOLAČNÝ POLYNÓM**\r\n');
    
            %Overenie pre Lagrangeov interpolačný polynóm, x-ové hodnoty sa nesmú opakovať
            %Zistenie unikátnych hodnôt
            unikatneHodnoty = unique(MaticaDataAproximacie(:,currentColumnA));
            % Ak je počet unikátnych hodnôt menší ako dĺžka poľa, znamená to, že v poli nie sú rovnaké hodnoty
            if length(unikatneHodnoty) < length(MaticaDataAproximacie(:,currentColumnA))
                disp('ERROR: Pole pre x obsahuje rovnaké hodnoty! Nespĺňa to podmienku na výpočet pomocou Lagrangeovho interpolačného polynómu. ');
                fprintf(fileAproximacia,'ERROR: Pole pre x obsahuje rovnaké hodnoty! Nespĺňa to podmienku na výpočet pomocou Lagrangeovho interpolačného polynómu.\r\n');
                currentColumnA=currentColumnA+2;
                continue;
            end
            
            xL=[];      %Pole na načítanie x-ových hodnôt
            yL=[];      %Pole na načítanie y-ových hodnôt

            %Cyklus ktorý načíta hodnoty, je tu preto aby sme sa vyhli prípadu ak bude v 1.riadku číslo, v 2.riadku NaN, 
            % v 3.riadku číslo, tzn. načíta potrebné hodnoty. Vyššie je to viac ošetrené
            for i=1:rowA
                %Preskočenie riadku s NaN
                if ~isnan(MaticaDataAproximacie(i,currentColumnA)) && ~isnan(MaticaDataAproximacie(i,currentColumnA+1))
                    %Zápis do poľa xL a yL
                    xL = [xL,MaticaDataAproximacie(i,currentColumnA)];
                    yL = [yL,MaticaDataAproximacie(i,currentColumnA+1)];
                end
            end

            %Výsledná dĺžka xL
            nL = length(xL);
            
            %Premenná pre výsledok
            VysledokLIP = 0; 
    
            %Hlavný cyklus sumy
            for i = 1:nL                                             
                citatel = 1;  % čitateľ
                menovatel = 1;  % menovateľ

                %Hlavný cyklus súčinu
                for j = 1:nL 
                    %Vyhnutie sa aktuálnemu stĺpcu 
                    if i ~= j                      
                        p = [1 (-1)*xL(j)];                    %jeden polynóm v čitateli              
                        citatel = conv(citatel,p);             %výpočet násobenia polynómov v čitateli
                        menovatel = menovatel*(xL(i) - xL(j)); %výpočet násobenia menovateli                    
                    end
                end

                %Sčítavanie jednotlivých súčinov
                VysledokLIP = VysledokLIP + (citatel/menovatel)*yL(i);
            end
    
            format rat    % Racionálny formát   

            %Zápis do súboru
            fprintf(fileAproximacia,'Lagrangeov interpolačný polynóm je %d.stupňa a má tvar: \r\n',(pocetx-1));           

            %Použitie funkcie na print polynómu (na konci súboru) 
            PrintovanieLIP(fileAproximacia,VysledokLIP,pocetx);

        %Koniec Lagrangeovho interpolačného polynómu
        end 

        disp('Vypočet úspešný. Nájdete ho v Outputfiles.');
    end

    %Prechod na nasledujúcu funkciu
    currentColumnA=currentColumnA+2;
    
end
disp('KONIEC.');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KONIEC HLAVNEJ ČASTI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%----------------------------------------------------------------------%
%Funkcia na delenie výpočtu po riadkoch
function NasledujuciRiadok(currentRow)
    %Pri prvom riadku komunikácia neprebieha
    if currentRow==1
        return;
    else    
        %Pri ostatných riadkoch komunikácia prebieha
        AnoNie = input('Chcete prejsť na nasledujúci riadok? (a - áno, n - nie (ukončiť program)): ', 's');
        while ~(strcmpi(AnoNie, 'a') || strcmpi(AnoNie, 'n'))
            AnoNie = input('ZLÝ VSTUP! Chcete prejsť na nasledujúci riadok? (a - áno, n - nie (ukončiť program)): ', 's');
        end
        %Prípad ak užívateľ chce program ukončiť predčasne
        if strcmpi(AnoNie, 'n')
            error('Predčasné ukončenie programu....');
        end
    end
end
 
%----------------------------------------------------------------------%
%Funkcia na printovanie Polynómu
function PrintovanieLIP(fileAproximacia,VysledokLIP,pocetx)
    %For sa vykoná toľkokrát, aký je polynóm stupňa
    for i=1:length(VysledokLIP)
        %S použitím mocniny
        if i<length(VysledokLIP)
            %S použitím formátu e pri malých číslach
            if abs(VysledokLIP(i))<0.009 && VysledokLIP(i)~=0
                fprintf(fileAproximacia, '%.4e*x^%d + ', VysledokLIP(i),pocetx-i);
            %Bez použitia formátu e pri väčších číslach    
            else
                fprintf(fileAproximacia, '%.4f*x^%d + ', VysledokLIP(i),pocetx-i);
            end
        %Bez použitia mocniny
        else
            if abs(VysledokLIP(i))<0.009 && VysledokLIP(i)~=0
                fprintf(fileAproximacia, '%.4e ', VysledokLIP(i));
            else
                fprintf(fileAproximacia, '%.4f ', VysledokLIP(i));
            end
        end
    end
    %Nový riadok na konci printovanie polynómu
    fprintf(fileAproximacia, '\r\n');
end

%----------------------------------------------------------------------%
%Funkcia na vykreslenie prvého stupňa + komunikácia
function VykrelseniePrvyStupen(xS,yS,a0,a1)
   
    AnoNie = input('Chcete zobraziť graf Polynómu 1. stupňa? (a - áno, n - nie): ', 's');
    %Ošetrenie vstupu 
    while ~(strcmpi(AnoNie, 'a') || strcmpi(AnoNie, 'n'))
       AnoNie = input('ZLÝ VSTUP! Chcete zobraziť graf Polynómu 1. stupňa? (a - áno, n - nie): ', 's');
    end
    if strcmpi(AnoNie, 'n')
       return;
    end

    x = linspace(min(xS), max(xS), 100); %Vytvorenie sady hodnôt x
    y = a0 + a1*x;                       %Vypočítanie hodnôt y pre každé x
    plot(x, y, 'b-', xS, yS, 'ro')       %Vykreslenie
    xlabel('x');                         %Stručný popis
    ylabel('y');
    title('Polynóm prvého stupňa');
    legend('Polynóm', 'Dáta', 'Location', 'best');
end
%----------------------------------------------------------------------%
%Funkcia na vykreslenie druhého stupňa + komunikácia
function VykrelsenieDruhyStupen(xS,yS,PrehodeneKoef)
    
    AnoNie = input('Chcete zobraziť graf Polynómu 2. stupňa? (a - áno, n - nie): ', 's');
    %Ošetrenie vstupu
    while ~(strcmpi(AnoNie, 'a') || strcmpi(AnoNie, 'n'))
        AnoNie = input('ZLÝ VSTUP! Chcete zobraziť graf Polynómu 2. stupňa? (a - áno, n - nie): ', 's');
    end
    if strcmpi(AnoNie, 'n')
         return;   %Ukončenie funkcie, nedojde k vykresleniu
    end
    
    xxx = [min(xS)-1:0.1:max(xS)+1];  %rozdelenie do viacerých bodov pre vykreslenie
    yyy = polyval(PrehodeneKoef,xxx); %Vypočíta funkčnú hodnotu polynómu v bodoch xxx
    plot(xS,yS,'o',xxx,yyy,'.');      %Vykreslenie
    xlabel('x');                      %Stručný popis
    ylabel('y');
    title('Polynóm druhého stupňa');
    legend('Polynóm', 'Dáta', 'Location', 'best');
end





