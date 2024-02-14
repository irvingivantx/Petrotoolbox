
package petrotoolbox;

import java.io.FileNotFoundException;
import java.util.HashMap;


public class oil {
    public Datareader properties= new Datareader();
    public HashMap<String,Double> Reservoirproperties;
    public double Pressurepsi,Temperaturef,Bubblepressurepsi;
    public double densityoil_standardconditions,api,oilgravity,gasdensity,gasgravity;
    public double oilratestbday,gasratescfday,GORscfbbl;
    public String saturated;
    
    public oil() throws FileNotFoundException{
        properties.AssignVariables();
    }
    public void compileallvariables(double Pressure,double Temperature, double densityoil_stc,double oilrate,double gasrate,double GOR, double gasgravity_,double api_,double oilgravity_,double Bubblepressurepsi_, String saturated_){
        Pressurepsi=Pressure;
        Temperaturef=Temperature;
        Bubblepressurepsi=Bubblepressurepsi_;
        densityoil_standardconditions=densityoil_stc;
        gasgravity=gasgravity_;
        oilratestbday=oilrate;
        gasratescfday=gasrate;
        GORscfbbl=GOR;
        api=api_;
        oilgravity=oilgravity_;
        saturated=saturated_;
                    
        if (densityoil_standardconditions==0){
            if(oilgravity==0){
            oilgravity=141.5/(api+131.5);
            densityoil_standardconditions=oilgravity*62.4;
            if(api==0){
            densityoil_standardconditions=oilgravity*62.4;
                }
            }
            else{
                    densityoil_standardconditions=oilgravity*62.4;
                }
        }
        if (oilgravity==0){
            if (api==0){
                oilgravity=densityoil_standardconditions/62.4;
            }
            if(densityoil_standardconditions==0){
                oilgravity=141.5/(api+131.5);
                
            }}
        }
    public double densityoil()throws FileNotFoundException{
        /*
        Reservoirproperties=properties.readreservoirproperties();
        
        GORscfbbl=Reservoirproperties.get("gor");
        oilgravity=Reservoirproperties.get("oilgravity");
        Temperaturef=Reservoirproperties.get("Temperature");
        */
                double densityoil=(62.4*properties.gasgravity+.0136*properties.GORscfbbl*properties.gasgravity)/(.972+.000147*Math.pow(properties.GORscfbbl*(Math.sqrt(properties.gasgravity/properties.oilgravity))+1.25*properties.Temperaturef,1.175));
                System.out.println(properties.gasgravity);
                System.out.println(properties.GORscfbbl);
                System.out.println(properties.oilgravity);
                System.out.println(properties.Temperaturef);
                return densityoil;
        
    }
    public void Bo()throws FileNotFoundException{
        Reservoirproperties=properties.readreservoirproperties();
        double bo=.9759+.00012*Math.pow(Reservoirproperties.get("gor")*Math.sqrt(Reservoirproperties.get("gasgravity")/Reservoirproperties.get("oilgravity"))+1.25*Reservoirproperties.get("Temperature"),1.2);
        System.out.print(bo);
    }
    public double viscosityoil()throws FileNotFoundException{
        Reservoirproperties=properties.readreservoirproperties();
        double t=Reservoirproperties.get("Temperature");
        
        double A=Math.pow(10,(.43+8.33/Reservoirproperties.get("api")));
        double yod=(.32+1.8*Math.pow(10,7)/Math.pow(api,4.53))*Math.pow(360/(t+200),A);
        double e_=3.74*Math.pow(10, -3)*Reservoirproperties.get("gor");
        double d=1.1*Math.pow(10,-3)*Reservoirproperties.get("gor");
        double c=8.62*Math.pow(10,-5)*Reservoirproperties.get("gor");
        double b=(.68/Math.pow(10, c))+(.25/Math.pow(10, d))+(.062/Math.pow(10, e_));
        double a=Reservoirproperties.get("gor")*(2.2*Math.pow(10,-7)*Reservoirproperties.get("gor")-7.4*Math.pow(10, -4));
        double yob=Math.pow(10, a)*Math.pow(yod, b);
        double yo=yob+.001*(Reservoirproperties.get("Pressurepsi")-Reservoirproperties.get("Bubblepressure"))*(.024*Math.pow(yob, 1.6)+.38*Math.pow(yob,.56));
        if (saturated=="Yes"){
            return yob;
        }
        else{
            return yo;
        }
        
    }
    public double solutiongasoilratio()throws FileNotFoundException{
        Reservoirproperties=properties.readreservoirproperties();
        double gor=Reservoirproperties.get("gasgravity")*Math.pow(Reservoirproperties.get("Pressurepsi")*Math.pow(10,(.0125*Reservoirproperties.get("api")))/(18*Math.pow(10, (.00091*Reservoirproperties.get("Temperature")))),1.2048);
        return gor;
    }
}
