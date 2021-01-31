function  [Bolting_point] = phenology (hour,T,sunrise,sunset,geno)

% Fixed parameters
%_________________

global p

Tb = p(1);
Tvmin = p(2);
Tvmax = p(3);
Vsat = p(4);
v = p(5);
sigma = p(6);
m = p(7);
Dld = p(8);
CSDL = p(9);
CLDL = p(10);

if  geno == 1
    
    Fb = p(11);
    Dsd = p(12);
    Threshold = p(13);
    Night = p(14);
    Phot_a = p(73);
    Phot_b = p(74);
    Phot_c = p(75);
    Phot_n = p(76);
    
else

    Fb = p(15);
    Dsd = p(16);
    Threshold = p(17);
    Night = p(18);
    Phot_a = p(77);
    Phot_b = p(78);
    Phot_c = p(79);
    Phot_n = p(80);

end

%Calculation begins
%__________________


    n=1;
    Th(1)=0;
    

    while	Th (n) <  Threshold
        
            i=n;
        
            %Calculating thermal component
            %_____________________________
    
            if      sunrise(i)>=hour(i) || sunset(i)<=hour(i)-1
                    fraction_light(i)=0;
            elseif  sunrise(i)<=hour(i)-1 && sunset(i)>hour(i)
                    fraction_light(i)=1;
            elseif  sunrise(i)>=hour(i)-1
                    fraction_light(i)=hour(i)-sunrise(i);
            else
                    fraction_light(i)=sunset(i)-hour(i)+1;
            end
            
           
    
            if      fraction_light(i)==0
                    Thrm(i)=Night*max(T(i)-Tb,0);
            elseif  fraction_light(i)==1    
                    Thrm(i)=max(T(i)-Tb,0);
            else 
                    Thrm(i)=max(0,(T(i)-Tb)*fraction_light(i)) + Night*max(0,(T(i)-Tb)*(1-fraction_light(i)));
            end
            
               
      
    
            %Calculating photoperiod component
            %_________________________________
            
            dl(i) = sunset(i) - sunrise(i);
            
                       
            if  i == 1
                        
                [FTarea24,yo] = link(dl(1),sunrise(1)); %initialise
                [FTarea1,yo] = sublink(hour(i),dl(i),sunrise(i),yo);
                FTarea24 = [FTarea24(2:end) FTarea1];
                dailyFTarea(i) = sum(FTarea24);                                       
            else
                [FTarea1,yo] = sublink(hour(i),dl(i),sunrise(i),yo);
                FTarea24 = [FTarea24(2:end) FTarea1];
                dailyFTarea(i) = sum(FTarea24);
            end
                             
                       
                              
            Phot(i) = Phot_a + Phot_b*Phot_c^Phot_n/(Phot_c^Phot_n+dailyFTarea(i)^Phot_n);
                            
                      
                    
    
            %Calculating vernalization effectiveness
            %_______________________________________

    
            if	T(i)>=Tvmin && T(i)<=Tvmax
	
                Ve(i)=exp(v)*(T(i)-Tvmin)^m*(Tvmax-T(i))^sigma*1;
            else	
                Ve(i)=0;
            end
    
             

            %Calculating cumulative vernalization hours
            %__________________________________________

     
            Vh(i)=sum(Ve(1:i));
    
      
    

            %Calculating vernalization fraction
            %__________________________________

            if	Vh(i)<=Vsat
	
                Vern(i) = Fb + Vh(i)*(1 - Fb)/Vsat;
            else
                Vern(i)= 1;	
            end
    
    
     
    

            %Calculating modified photothermal unit
            %______________________________________

            mptu(i) = Vern(i)*Phot(i)*Thrm(i);
      
      
      

      
            %Calculating switching threshold
            %_______________________________

            Th(n+1)=sum(mptu(1:n));
            n=n+1;  
    
    end


    %Output:
    %-------

    DtB = (n-1)/24;
    Bolting_point = n-1;

            


