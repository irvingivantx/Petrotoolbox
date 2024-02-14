
package petrotoolbox;

import java.io.FileNotFoundException;
import java.util.HashMap;

public class Reservoir {
    public Datareader properties= new Datareader();
    public gas gas1=new gas();
  
    public Reservoir()throws FileNotFoundException
    {
        //properties.AssignVariables();
    }  
    
    public double Wellboreapparentradius_ft(){
        double rw=properties.Wellboreradius_ft;
        double apparent=rw*Math.exp(-properties.Skin);
        return apparent;
    }
    public double Permeabilitypoint()throws FileNotFoundException{
    /*
    'Calculates formation permeabililty from a single rate from a flow test
    'From SPE 012847 "Estimating Formation Permeability From Single-Point Flow Data"
    '1984 by Lee, Holditch, and McVay
    '''*/

    double Pi = properties.Pressurepsi; //#'initial reservor pressure, psia
    double TresR = properties.Temperaturef + 459.67; //#'reservoir temperature, R
    double sg = properties.gasgravity;
    double Qg = properties.gasratescfday; //#'Gas rate, Mcf/d
    double Timehrs = properties.flowtime ; //#'length of flow time, hours
    // #'flowing bottomhole pressure, psia
    double pwf=properties.flowingbottomholepressure;
    double Pay = properties.Reservoirheight_ft;// #'Net Pay, feet
    double Porosity = properties.porosity / 100;// #'porosity as a fraction
    double sw = properties.Watersaturation / 100;// #'water saturation as a fraction
    double rw = properties.Wellboreradius_ft;// #'wellbore radius, feet
    double mpi = gas1.pseudopressure();// #'pseudo-pressure of initial reservoir pressure
    properties.Pressurepsi=pwf;
    System.out.println(pwf);
    double mpwf = gas1.pseudopressure();// #'pseudo-pressure of flowing bottomhole pressure
    properties.Pressurepsi=Pi;
    double ugi = gas1.viscositygas_cp();// #'viscosity of gas at inital reservoir pressure
    double Cti = gas1.Ctcompressibility();//#'total compressibiility at initial reservoir pressure
    double Kold = 1;// #'Initial guess for permeability
    int i=1;
    double Knew=0;
    while (i<2){
        double rd = Math.pow(((Kold * Timehrs) / (377 * Porosity * ugi * Cti)),0.5);
        double ONE = 1422 * Qg * TresR;
        double TWO = Pay * (mpi - mpwf);
        double THREE = Math.log(rd / rw)/Math.log(2.71828182846) - 0.75 + properties.Skin;
        Knew = (ONE / TWO) * THREE;
        double Kdiff = Math.abs((Knew - Kold) / Knew);
        
        if (Kdiff < 1e-05){
            break;
        }
        Kold = Knew;
        i++;
         }
    double permeability = Knew;
    
    return permeability;
    }
}
