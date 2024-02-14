
package petrotoolbox;
import java.io.*;
import java.util.*;

public class Datareader {
    
        public double Surfacepressurepsi;
        public double Surfacetemperaturef;
        public double Pressurepsi;
        public double Temperaturef;
        public double Bubblepressurepsi;
        public double flowingbottomholepressure;
        public double densityoil_standardconditions;
        public double api;
        public double oilgravity;
        public double gasdensity;
        public double gasgravity;
        public double oilratestbday;
        public double gasratescfday;
        public double GORscfbbl;
        public double molecularweight;
        public double nitrogenpercent;
        public double carbondioxdedpercent;
        public double hydrogensulfitepercent; 
        public double wgr_stbmmcf;
        public double Pressurestandardpsi=14.7;
        public double Temperaturestandardf=60;
        public double Flowlengthft;
        public double TrueVerticaldepthft;
        public double TubingIDin;
        public double TubingODin;
        public double CasingIDin;
        public double Flowtype;
        public double Reservoirheight_ft;
        public double Wellboreradius_ft;
        public double porosity;
        public double Watersaturation;
        public double Permeability_md;
        public double Area_ac;
        public double Skin;
        public String saturated;
        public double flowtime;
        
    public Datareader(){
        
    }
    
    public HashMap<String, Double> readreservoirproperties() throws FileNotFoundException{
       
        Scanner scanner= new Scanner(new FileReader("Untitled.txt"));
        HashMap<String,Double> map= new HashMap<>();
        
        while (scanner.hasNextLine()){
            String[] columns=scanner.nextLine().split(" ");
            map.put(columns[0],Double.parseDouble(columns[1]));
        }
        return map;
        
    }
    public void AssignVariables() throws FileNotFoundException{
        HashMap<String,Double> properties= new HashMap<>();
        properties=readreservoirproperties();
        
        Surfacepressurepsi=properties.get("SurfacePressure");
        Surfacetemperaturef=properties.get("SurfaceTemperature");
        Pressurepsi=properties.get("Pressure");
        Temperaturef=properties.get("Temperature");
        Bubblepressurepsi=properties.get("Bubblepressure");
        flowingbottomholepressure=properties.get("Flowingbottomholepressure");
        densityoil_standardconditions=properties.get("densityoilstandard");
        api=properties.get("api");
        oilgravity=properties.get("oilgravity");
        gasdensity=properties.get("densitygas");
        gasgravity=properties.get("gasgravity");
        oilratestbday=properties.get("oilrate_stbday");
        gasratescfday=properties.get("gasrate_scfday");
        GORscfbbl=properties.get("GORscfbbl");
        molecularweight=properties.get("molecularweight");
        nitrogenpercent=properties.get("Nitrogenpercent");
        carbondioxdedpercent=properties.get("Carbondioxidepercent");
        hydrogensulfitepercent=properties.get("Hydrogensulfidepercent"); 
        wgr_stbmmcf=properties.get("wgr_stbmmcf");
        Pressurestandardpsi=properties.get("Pressurestandardpsi");
        Temperaturestandardf=properties.get("Temperaturestandardf");
        Flowlengthft=properties.get("flowlenght_ft");
        TrueVerticaldepthft=properties.get("TVD_ft");
        TubingIDin=properties.get("TubingIDin");
        TubingODin=properties.get("TubingODin");
        CasingIDin=properties.get("CasingIDin");
        Flowtype=properties.get("flowtype");
        Reservoirheight_ft=properties.get("Reservoirheight_ft");
        Wellboreradius_ft=properties.get("Wellboreradius_ft");
        porosity=properties.get("Porosity");
        Watersaturation=properties.get("WaterSaturation_percent");
        Permeability_md=properties.get("Permeability_md");
        Area_ac=properties.get("Area_ac");
        Skin=properties.get("skin");
        flowtime=properties.get("flowtime_hrs");
        
    }
    
}
