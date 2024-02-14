
package petrotoolbox;

public class test {
    public double number1;
    public test(){
        
    }
    public void first(){
        System.out.println(number1);
    }
    public void second(){
    double pwf=0;
    double Qg; //Gas Rate in Mcf/day
    double SurfTempR ;
    double BHTempR;
    double BHTempf;
    double L ;
    double h ;
    double d ;
    double D1;
    double D2 ;
    double Gw ;//'specific gravity of gas at the wellhead
    double WGR; //'converts water-gas-ratio from STB/MMcf to STB/scf
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
    double XIMP;
    double PMPN;
    double DP;
    int i=1;
    while (i<200){
        double temp=pwf;
        pwf=temp+1;
        
        System.out.println(pwf);
        i++;
    }
    //pwf=4+2;
    System.out.println(pwf);
    }
}
