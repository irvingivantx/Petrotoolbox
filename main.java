
package petrotoolbox;

import java.io.FileNotFoundException;
import java.util.HashMap;


public class main {

  
    public static void main(String[] args) throws FileNotFoundException{
      
     
      Reservoir reservoir=new Reservoir();
      //oil fluid1= new oil();
      gas gas1=new gas();
      //test yu= new test();
      BoreholeCollapse well= new BoreholeCollapse();
      //Reservoir properties
      double SurfacePressure=100;
      double SurfaceTemperature=80;
      double Pressure=6500;
      double Temperature=250;
      double Flowingbottomholepressure=500;
      double densityoilstandard=49;
      double api=35;
      double oilgravity=.65;
      double densitygas=.7;
      double gasgravity=.7;
      double oilrate=500;
      double gasrate=0;
      double gor=600;
      double wgr=0;
      double molecularweight=20.278999999999996;
      double Bubblepressure=2745;
      String saturated="Yes";
      //Tubing Properties
      double flowlenght_ft=10000;
      double TVD_ft=10000;
      double TubingIDin=3;
      double TubingODin=4;
      double CasingIDin=5;
      double flowtype=0;
      //Reservoir Properties
      double Reservoirheight_ft=50;
      double Porosity=15;
      double Sw=30;
      double rw=.3;
      double Permeability=0.08;
      double Area=40;
      double skin=-3;
    
      //yu.second();
      //gas1.compileallvariables(SurfacePressure,SurfaceTemperature,Pressure, Temperature,  densityoilstandard, oilrate, gasrate, gor,  gasgravity, api,oilgravity,Bubblepressure,saturated,molecularweight,flowlenght_ft,TVD_ft,TubingIDin,TubingODin,CasingIDin,wgr,flowtype,Flowingbottomholepressure,Reservoirheight_ft,Porosity,Sw,rw,Permeability,Area,skin);
      //fluid1.densityoil();
      //fluid1.viscosityoil();
      //fluid1.solutiongasoilratio();
      //gas1.zfactor(0,0,0,0,0,0);
      //gas1.Cgcompressibilitygas();
      //gas1.Gasbottomholepressurepsi();
      //gas1.Tubingheadpressure();
      //gas1.GasflowratePSS_Mcfd();6
      //fluid1.densityoil();
      //gas1.pseudopressure();
      //gas1.GasflowratePSS_Mcfd();
      //reservoir.Permeabilitypoint();
      //well.calculate(3,1);
      //well.calculate(3,2);
      well.startBoreholecollapse(30,30,2);
      
      }
    
}
