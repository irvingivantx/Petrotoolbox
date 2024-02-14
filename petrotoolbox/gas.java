
package petrotoolbox;


import java.io.FileNotFoundException;
import java.util.HashMap;

public class gas {
    public Datareader properties= new Datareader();
    public HashMap<String,Double> Reservoirproperties;
   
   
    
 
    public gas()throws FileNotFoundException{
        //properties.AssignVariables();
    }
    public double pseudopressure() throws FileNotFoundException{
    double Mpp=0;
    double Pold=0;
    double Xold=0;
    double Xnew=0;
    
    double Pstep = properties.Pressurepsi / 20;
    System.out.println(properties.Pressurepsi+"Printng this pressure");
    double Pnew=0;
    
    int n=1;
   
    while (n<21){
      
        Pnew = Pold + Pstep;
        properties.Pressurepsi=Pnew;
        Xnew = 2 * Pnew / zfactor()/viscositygas_cp();
        Mpp = Mpp + (Xold + Xnew) / 2 * Pstep;
        Pold = Pnew;
        Xold = Xnew;
        n++;
        
        }
        //print (jMpp, "the Mpp")
        
    return Mpp;
    }
    public double Pseudocriticalpressure() throws FileNotFoundException{
        
        double sg=properties.gasgravity;
        double co2=properties.carbondioxdedpercent/100;
        double h2s=properties.hydrogensulfitepercent/100;
        double x=co2+h2s;
        double y=h2s;
        double eps=120*(Math.pow(x, .9) - Math.pow(x, 1.6)) + 15 * (Math.pow(y, .5) - Math.pow(y, 4));
        double ptca = 169.2 + 349.5 * sg - 74 * Math.pow(sg, 2);
        double ptc = ptca - eps;
        double ppca = 756.8 - 131 * sg - 3.6 * Math.pow(sg,2);
        double pseudocriticalpressure_psi = (ppca * ptc) / (ptca + h2s * (1 - h2s) * eps);
        
        return pseudocriticalpressure_psi;
    }
    public double Pseudocriticaltemperature()throws FileNotFoundException{
         
        double sg=properties.gasgravity;
        double co2=properties.carbondioxdedpercent/100;
        double h2s=properties.hydrogensulfitepercent/100;
        double x=co2+h2s;
        double y=h2s;
        double eps=120*(Math.pow(x, .9) - Math.pow(x, 1.6)) + 15 * (Math.pow(y, .5) - Math.pow(y, 4));
        double tc=169.2+349.5*sg-74*Math.pow(sg, 2);
        double criticaltemperature=tc-eps;
        return criticaltemperature;
    }
    public double zfactor()throws FileNotFoundException{
        
        double Tpc=326+315.7*(properties.gasgravity-.5)-240*properties.nitrogenpercent-83.3*properties.carbondioxdedpercent+133.3*properties.hydrogensulfitepercent;
        double Ppc=678-50*(properties.gasgravity-.5)-206.7*properties.nitrogenpercent+440*properties.carbondioxdedpercent+606.7*properties.hydrogensulfitepercent;
     
        double Tpr=(properties.Temperaturef+460)/Tpc;
        double tr=1/Tpr;
        double Ppr=properties.Pressurepsi/Ppc;
        double A=.06125*tr*Math.exp(-1.2*Math.pow((1-tr),2));
        double B=tr*(14.76-9.76*tr+4.58*Math.pow(tr, 2));
        double C=tr*(90.7-242.2*tr+42.4*Math.pow(tr, 2));
        double D=2.18+2.82*tr;
        double zero;
        double zeroderivative;
        double z;
        double Y1;
        double Y=5;
        double power=D-1;
        int i=0;
        //for (int i=0;i<100;i++){
        while(i<10000){
        zero=(Y+Math.pow(Y, 2)+Math.pow(Y, 3)-Math.pow(Y, 4))/Math.pow((1-Y),3)-A*Ppr-B*Math.pow(Y, 2)+C*Math.pow(Y, D);
        zeroderivative=(1+4*Y+4*Math.pow(Y, 2)-4*Math.pow(Y, 3)+Math.pow(Y, 4))/(Math.pow(1-Y, 4))-2*B*Y+C*Y*Math.pow(Y, D-1);
        Y1=Y-zero/(zeroderivative);
        Y=Y1;
       
        i++;
        }
        z=A*Ppr/Y;
       
        return z;
        
        }
    public double viscositygas_cp() throws FileNotFoundException{
   
        double zfactor=zfactor();
        double a = (9.379 + 0.01607 * properties.molecularweight) * Math.pow((properties.Temperaturef + 459.67),1.5)/(209.2 + 19.26 * properties.molecularweight + (properties.Temperaturef + 459.67));
        double b = 3.448 + 986.4 / (properties.Temperaturef + 459.67) + 0.01009 * properties.molecularweight;
        double c = 2.447 - 0.2224 * b;
        double rho = properties.Pressurepsi*properties.molecularweight/zfactor/669.8/(properties.Temperaturef + 459.67);
        //rho=1
        double viscosity = (a*Math.exp(b * Math.pow(rho, c))) / 10000;
        
        return viscosity;
    }   
    public double densitygas_lbmft3() throws FileNotFoundException{
        double z_=zfactor();
        double densitygas=2.7*properties.gasdensity*properties.Pressurepsi/(z_*(properties.Temperaturef+460));
        return densitygas;
    }
    public double Bg_RBperSCF() throws FileNotFoundException{
        double zfactor=zfactor();
        double Bg = (zfactor * (properties.Temperaturef + 459.67) * properties.Pressurestandardpsi) / (properties.Pressurepsi * (properties.Temperaturestandardf + 459.67) * 5.6145833333); //RB/scf
        return Bg;
    }
    public double Cgcompressibilitygas() throws FileNotFoundException{
    double TR = (properties.Temperaturef + 459.67) / Pseudocriticaltemperature();
    double PR = properties.Pressurepsi/ Pseudocriticalpressure();
    double a = 0.064225133;
    double b = 0.53530771 * TR - 0.61232032;
    double c = 0.31506237 * TR - 1.0467099 - 0.57832729 / Math.pow(TR, 2);
    double d = TR;
    double e = 0.68157001 / Math.pow(TR,2);
    double f = 0.68446549;
    double g = 0.27 * PR;
    double rho = 0.27 * PR / TR ;//#Initial guess;
    double rhoold = rho;
    double Cg;
    int i=0;
    while (i<1000){
        double frho = a * Math.pow(rho, 6) + b * Math.pow(rho, 3)+ c * Math.pow(rho, 2) + d * rho + e * Math.pow(rho, 3) * (1 + f * Math.pow(rho, 2)) * Math.exp(-f *Math.pow(rho, 2)) - g;
        double dfrho = 6 * a * Math.pow(rho, 5) + 3 * b * Math.pow(rho, 2) + 2 * c * rho + d + e * Math.pow(rho, 2) * (3 + f * Math.pow(rho,2) * (3 - 2 * f * Math.pow(rho, 2))) * Math.exp(-f * Math.pow(rho, 2));
        rho = rho - frho / dfrho;
        double test = Math.abs((rho - rhoold) / rho);
        if (test < 1e-05){
            break;
                    }
        rhoold = rho;
        }
        double zm = 0.27 * PR / rho / TR;
        double der = 1 / rho / TR * (5 * a * Math.pow(rho, 5) + 2 * b * Math.pow(rho, 2)  + c * rho + 2 * e * Math.pow(rho, 2)  * (1 + f * Math.pow(rho, 2)  - Math.pow(f, 2)  * Math.pow(rho, 4) ) * Math.exp(-f * Math.pow(rho, 2) ));
        double cr = 1 / PR / (1 + rho / zm * der);
        Cg = zm;
        Cg = cr / Pseudocriticalpressure();
        i++;
    
    return Cg;
    }
    public double Ctcompressibility() throws FileNotFoundException{
        double sw= properties.Watersaturation/100;
        double Cg=Cgcompressibilitygas();
        double ct=Cg*(1-sw);
        return ct;
    }
    public double Gasbottomholepressurepsi() throws FileNotFoundException{
    properties.AssignVariables();
    double pwf=0;
    double Qg = properties.gasratescfday; //Gas Rate in Mcf/day
    double SurfTempR = properties.Surfacetemperaturef + 459.67;
    double BHTempR = properties.Temperaturef + 459.67;
    double BHTempf=properties.Temperaturef;
    double L = properties.Flowlengthft;
    double h = properties.TrueVerticaldepthft;
    double d = properties.TubingIDin;
    double D1 = properties.TubingODin;
    double D2 = properties.CasingIDin;
    double Gw = properties.gasgravity;//'specific gravity of gas at the wellhead
    double WGR = properties.wgr_stbmmcf / 1000000; //'converts water-gas-ratio from STB/MMcf to STB/scf
    double DPL;
    double DM;
    double DE5;
    double OMEGAWH;
    double XIBH=0;
    double NRE;
    double RR;
    double f;
    
    double cheetos = 0;
    double ALPHA=0;
    double Viscg = 0;
    double Z =0; //zfactor();
    //#Z=1
    double PTZ = 0;
    double XIWH =0;
    double TMPF;
    double TMPR;
    double OMEGAMP;
    double XIMP=0;
    double PMPN;
    double DP;
   //Calculate the first estimate of midpoint pressure
    double PMP = properties.Surfacepressurepsi+(ALPHA/(XIWH+XIWH));
        if (properties.Flowtype == 0){ //FlowType = 0 for flow up the tubing or casing, = 1for flow up the tubing/casing annulus
            DPL = d;
            DM = d;
            DE5 = Math.pow(d, 5);

            //Initialize equivalent flow diameters for flow up the tubing/casing annulus
                }
        else{
            DPL = D2 + D1;
            DM = D2 - D1;
            DE5 = (Math.pow(DM, 3))* (Math.pow(DPL, 2));

            }
    //Calculate ALPHA
    ALPHA = (Gw / 53.34 + WGR * 86.27) * L;
    //Calculate OMEGA at wellhead conditions
    properties.Pressurepsi=properties.Surfacepressurepsi;
    properties.Temperaturef=properties.Surfacetemperaturef;
    Viscg=viscositygas_cp();
    
    //Viscg = viscositygas_cP(3000, 200, 5,2,2)
    if (Qg == 0)
        {
        OMEGAWH = 0;
        }
    else{
        NRE = 20.011 * Gw * Qg / (Viscg * DPL);
        RR = 0.0023 / DM;
        f = 1 / Math.pow(((2.28 - 4 * Math.log10(RR + 21.25 / (Math.pow(NRE, .9))))),2);
        OMEGAWH = 0.0026665 * f * Math.pow(Qg, 2) / DE5;
        }
    //Evaluate the integral I at wellhead conditions
    Z=zfactor();
    //#Z=1
    PTZ = properties.Surfacepressurepsi / (SurfTempR * Z);
    XIWH = (199.3 * WGR * Math.pow(PTZ, 2) + PTZ) / (OMEGAWH + (Math.pow(PTZ, 2)) * h / L);
   //Calculate the first estimate of midpoint pressure
    PMP = properties.Surfacepressurepsi+(ALPHA/(XIWH+XIWH));
   
    //Evaluate the intergral I at the midpoint conditions
    int i=1;
 
    while (i<3){
        TMPF = (properties.Surfacetemperaturef + BHTempf) / 2;
        TMPR = TMPF + 459.67;
        properties.Pressurepsi=PMP;

        properties.Temperaturef=TMPF;

        Z =zfactor();
        //#Z=1
        Viscg =viscositygas_cp();
            if (Qg == 0){
                OMEGAMP = 0;
                        }
            else{
                NRE = 20.011 * Gw * Qg / (Viscg * DPL);
                RR = 0.0023 / DM;
                f = 1 / Math.pow(((2.28 - 4 * Math.log10(RR + 21.25 / (Math.pow(NRE, .9))))),2);
                OMEGAMP = 0.0026665 * f * Math.pow(Qg, 2)/ DE5;
                }
        PTZ = PMP / (TMPR * Z);
        XIMP = (199.3 * WGR * Math.pow(PTZ, 2) + PTZ) / (OMEGAMP + (Math.pow(PTZ, 2)) * h / L);
        //Recalculate the midpoint pressure until convergence by iteration;
        PMPN = ALPHA / (XIMP + XIWH) + properties.Surfacepressurepsi;
        DP = PMP - PMPN;
        PMP = PMPN;
            if (Math.abs(DP) < .00001)
                {
                    break;
                }
            
         i++;     
        }
    //Calculate a first estimate of bottomhole pressure
    pwf=PMP+(ALPHA/(XIMP+XIMP));
    //Evaluate the integral I at bottomhole conditions
    double j=1;
    double OMEGABH;
    double PBHN;
    while (j<3){
            properties.Pressurepsi=pwf;
            properties.Temperaturef=BHTempf;
            //Z = zfactor(pwf, BHTempF, GasGravity, H2S_Percent,CO2_Percent)
            //Viscg = viscositygas_cP(pwf, BHTempF, Gw,H2S_Percent,CO2_Percent)
            Z=zfactor();
            Viscg=viscositygas_cp();

            if (Qg == 0){
                OMEGABH = 0;
                  }
            else{
                NRE = 20.011 * Gw * Qg / (Viscg * DPL);
                RR = 0.0023 / DM;
                f = 1 / Math.pow(((2.28 - 4 * Math.log10(RR + 21.25 / (Math.pow(NRE, .9))))),2);
                OMEGABH = 0.0026665 * f * Math.pow(Qg, 2) / DE5;
                        }
            PTZ = pwf / (BHTempR * Z);
            XIBH = (199.3 * WGR * Math.pow(PTZ, 2)+ PTZ) / (OMEGABH + (Math.pow(PTZ, 2)) * h / L);
            //Recalculate the bottomhole pressure until convergence by iteration
            PBHN = ALPHA / (XIBH + XIMP) + PMP;
            DP = pwf - PBHN;
            pwf = PBHN;
            //System.out.print(PTZ+"Pwf");
            if (Math.abs(DP) < .00005){
                break;
            }
            j++;
        }
    //#Apply Simpsons rule to obtain a more accurte estimate of PWF
    pwf = properties.Surfacepressurepsi + (6 * ALPHA) / (XIWH + 4 * XIMP + XIBH);
    double bottomholepressure = pwf;
    //#print(jBHPRESS, "psia")
    
    
    return bottomholepressure;
    }
    /*public double Tubingheadpressure()throws FileNotFoundException{
        double NewFTP = 5500;
        double OldFTPhigh = 13000;
        double OldFTPlow = 50;
        double NewPwf=0;
        int i=1;
        while(i<50){
           
            properties.Pressurepsi=NewPwf;
            NewPwf =Gasbottomholepressurepsi();
            double check =  Math.abs((NewPwf - flowingbottomholepressure) / NewPwf);
            if (check < 0.0001){
                break;
            }
            if (NewPwf > flowingbottomholepressure){
                OldFTPhigh = NewFTP;
                        }
            if (NewPwf < flowingbottomholepressure){
                OldFTPlow = NewFTP;
            }
            NewFTP = (OldFTPlow + OldFTPhigh) / 2;
            i++;
           
        }
    
    return NewFTP;
    }
    */
    
    public double GasflowratePSS_Mcfd() throws FileNotFoundException{
        double Pbar = properties.Pressurepsi; //#'average reservor pressure, psia
        double TresR = properties.Temperaturef + 459.67; //#'reservoir temperature, R
        double sg = properties.gasgravity;
        double pwf = properties.flowingbottomholepressure;// #'flowing bottomhole pressure, psia
        double Pay = properties.Reservoirheight_ft;// #'Net Pay, feet
        double Porosity_ = properties.porosity/100; //#'porosity as a fraction
        double sw_ = properties.Watersaturation / 100;// #'water saturation as a fraction
        double rw = properties.Wellboreradius_ft;// #'wellbore radius, feet
        double Perm_md=properties.Permeability_md;
        double mpbar =0; 
        mpbar=pseudopressure(); //Mpp(Pbar, ResTemp_F, sg, H2S_Percent, CO2_Percent) #'pseudo-pressure of the average reservoir pressure
        
        
        double mpwf =pseudopressure(); //Mpp(pwf, ResTemp_F, sg,H2S_Percent, CO2_Percent) #'pseudo-pressure of flowing bottomhole pressure
        
        double ugpwf =0.014468928229831857;//viscositygas_cp(); //viscositygas_cP(pwf, ResTemp_F, sg, H2S_Percent, CO2_Percent) #'viscosity of gas at flowing bottomhole pressure
        double re = Math.pow((properties.Area_ac * 43560 / Math.PI),0.5); //#'radius of the reservoir's drainage area assuming a circular area in feet
        double ks =0.19310202717008523; //skinperm_md(Perm_md, WellRadius_ft, skin, DamageRadius_ft) #'permeability of the damaged zone in md
        double Dee = (6e-05 * sg * (Math.pow(ks, -.1)) * Pay) / (ugpwf * rw * Math.pow(Pay, 2));
        double d = (1424 * TresR * Dee) / (Perm_md * Pay);
        double c = (1424 * TresR / (Perm_md * Pay)) * (Math.log(re / rw) - 0.75 + properties.Skin);
        double jQgPSS_Mcfd = (-c + Math.pow((Math.pow(c, 2)+ 4 * d * (mpbar - mpwf)),0.5))/ (2 * d);
        
        return jQgPSS_Mcfd;
    }
    }
    
