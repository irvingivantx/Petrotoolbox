/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package petrotoolbox;

import java.text.DecimalFormat;
import java.util.ArrayList;
import javax.swing.table.DefaultTableModel;

/**
 *
 * @author irving
 */
public class BoreholeCollapse {
    public double s1,s2, s3, alpha, beta, gamma;

    public double shmaxazimuth=119.4582013;
    public double shmax =33.1714;
    public double sv=31.3796;
    public double shmin=29.4087;
    public double wellboreinclination=0.49864228654592857;
    public double wellboreazimuth=243.27830905805945;
    public double StressRegime;
    public double pp=22.4209;
    public double internalfriction=26.61230750766;
    public double ucs=6.0419;
    public double inc2=.5;
    private double PI = 4 * Math.atan(1);
    public double ss[][] = new double[3][3];
    public double rs[][] = new double[3][3];
    public double rst[][] = new double[3][3];
    public double rb[][] = new double[3][3];
    public double rbt[][] = new double[3][3];
    public double s[][] = new double[3][3];
    public boolean cont =true;
    public double rr,tt,zz,tz,tmax,tmin,w;
    public double poissonratio=.25;
    public double incp;
    public int failed[][] = new int[721][2];
    public final int _resetArr[][] = new int[721][2];
    public int FNodes[] = new int[8];
    public double mudWeight, BKOAngle, BKODir, TensANgle, TensDir, firstBKO;
    private Sigma sigmas = new Sigma();
//    Sigma[] sigmas = new Sigma[721];

    public double md;
   
    public int changes = 5;
    public int ch = 0;
    public double previousMW = Double.NaN;
    public double mwIncrimentFactor = 1;
    public double divFactorForIncrDecrFactor = 2;
    public double inc = 1;
   
    public double angleTolerance = 0.5;
//    public double tol = 0.0005;
    public double tol = 0.005;
    public int f1 = 0;
    public double mu;
    public double muT;
    public double fmuTi;
   
    public double F;
    public int maxVal = 0;
    public int count = 0;
    public double angleInRad;
    public static int MC2D_MODE = 1;
    public static int LADE_MODE = 2;
    public static int MC3D_MODE = 3;

    public int calculationMode = MC2D_MODE;

//    double tmp[][] = new double[3][3];
//Lade Setings
    public double lade_S1;
    public double lade_h;
    public double lade_tail;
    public double lade_I1;
    public double lade_I3;
    public double ladeFailVal;

    //MC3D Settings
    public double mu_i_sin;
    public double mohr_S1;
    public double mohr_3d_I1;
    public double mohr_3d_I2;
    public double mohr_3d_P;
    public double mohr_3d_J;
    public double mohr_3d_theta;
    public double mohr_3d_Func_val_theta;

    public ArrayList<Double> gradientMWList = new ArrayList<Double>();

    public boolean isNewCalculation = true;//Flag for new calculation 
//    public boolean isPolarPlotCalMode = false;//Flag for new calculation 
    public int maxSteps =40;
    public String measuredDepthUnit;
    public double tvd;
  
    public DecimalFormat twoDecimalFormat = new DecimalFormat("##.##");
    public DefaultTableModel polatResultsModel;
    public boolean updateDisplay = false;
    public double polarChartAngleIncriment;

    public int currentAzimuth = 0;
    public int currentInclination = 0;
 
    public BoreholeCollapse(){
        
    }
    public void startBoreholecollapse(double Horizontal, double Vertical, int failure_criteria){
        double Angle=0;
        Angle= Horizontal+(90-wellboreinclination)*(Vertical-Horizontal)/90;
        calculate(Angle, failure_criteria);
        
    }
    public void calculate(double BKO,int failurecriteria) {
        
        calculationMode=failurecriteria;
        //System.out.println(mudWeight);
//        int err, iu;
        int steps = 0;
        getRegime();
        //calculate Stresses matrices
        _initMatrices();
       
        //reset variavles
       _reAssignVariables();

        int n = (int) (360 / inc2) + 1;
        int iTmax = 0;
        int iTmin = 0;
        //System.out.println(n);
        int mwSign = 1;
        while (cont) {
            //steps++;
            _resetInnerVariables(n);
            iTmax = 0;
            iTmin = 0;
            for (int i = 0; i < n; i++) {
                angleInRad = i * inc2;
                angleInRad = angleInRad * Math.PI / 180;
                //Compute Stresses
                _computeStresses(i, angleInRad);

                if (sigmas.tmax > tmax) {
                    tmax = sigmas.tmax;
                    iTmax = i;
                }
                if (sigmas.tmin < tmin) {
                    tmin = sigmas.tmin;
                    iTmin = i;
                }

                _getPrincipalStresses(i);
                //Calculate Failed
                _calculatedFailed(i);
                //Detect changes
                _detectChanges(i);
                
            }
            if (maxVal > 8) {
                mudWeight = Double.NaN;
                return;
            }
            if (isNewCalculation) {
                _computeAngles2(n, iTmax, iTmin, inc2, maxVal);
            } else {
                _computeAngles(n, iTmax, iTmin, inc2, maxVal);
            }
            //calculate mudweight     
            _calculateMudWeight(BKO, mwSign);
            if (steps > maxSteps) {
                cont = false;
            }
            
        }
        System.out.println("Mudweight: "+mudWeight);
        System.out.println("Breakout Direction : "+BKODir);
        System.out.println("Breakout Angle: "+BKOAngle);
    }
    public int getRegime(){ 
        if (sv >= shmax && sv >= shmin) {
            s1 = sv;
            s2 = shmax;
            s3 = shmin;
            alpha = shmaxazimuth - 90;
            beta = 90;
            gamma = 0;
            alpha = Math.toRadians(alpha);
            beta = Math.toRadians(beta);
            gamma = Math.toRadians(gamma);
//            System.out.println(s1);
//            System.out.println(s2);
//            System.out.println(s3);
//            System.out.println(alpha);
//            System.out.println(beta);
//            System.out.println(gamma);
//            alpha = alpha * PI / 180;
//            beta = beta * PI / 180;
//            gamma = gamma * PI / 180;
            return 1;
//        } else if (shmax > shmin && shmin > sv) {
        } else if (shmax >= shmin && shmin >= sv) {
            s1 = shmax;
            s2 = shmin;
            s3 = sv;
            alpha = 360 - shmaxazimuth;
            beta = 0;
            gamma = 0;
            alpha = Math.toRadians(alpha);
            beta = Math.toRadians(beta);
            gamma = Math.toRadians(gamma);
            return 2;
//        } else if (shmax > sv && sv > shmin) {
        } else if (shmax > sv && sv >= shmin) {
            s1 = shmax;
            s2 = shmin;
            s3 = sv;
            alpha = 360 - shmaxazimuth;
            beta = 0;
            gamma = 90;
            alpha = Math.toRadians(alpha);
            beta = Math.toRadians(beta);
            gamma = Math.toRadians(gamma);
            return 3;
        } else {
            return 7;  // 7 is error
    }
}
    public void _initMatrices(){
        _iniS();
        _sRotate();

        for (int i = 0; i < 3; i++) {
            s[i][i] = s[i][i] - pp;
        }
        
       
    }
    public void _iniS(){
        ss[0][0] = s1;
        ss[1][1] = s2;
        ss[2][2] = s3;
      
    }
   private void _sRotate() {
        rs[0][0] = Math.cos(alpha) * Math.cos(beta);
        rs[0][1] = Math.sin(alpha) * Math.cos(beta);
        rs[0][2] = -Math.sin(beta);
        rs[1][0] = Math.cos(alpha) * Math.sin(beta) * Math.sin(gamma) - Math.sin(alpha) * Math.cos(gamma);
        rs[1][1] = Math.sin(alpha) * Math.sin(beta) * Math.sin(gamma) + Math.cos(alpha) * Math.cos(gamma);
        rs[1][2] = Math.cos(beta) * Math.sin(gamma);
        rs[2][0] = Math.cos(alpha) * Math.sin(beta) * Math.cos(gamma) + Math.sin(alpha) * Math.sin(gamma);
        rs[2][1] = Math.sin(alpha) * Math.sin(beta) * Math.cos(gamma) - Math.cos(alpha) * Math.sin(gamma);
        rs[2][2] = Math.cos(beta) * Math.cos(gamma);

        rst = getTranspose(rs);

        double delta1 = wellboreazimuth * PI / 180;
        double fi1 = wellboreinclination * PI / 180;

        rb[0][0] = -Math.cos(delta1) * Math.cos(fi1);
        rb[0][1] = -Math.sin(delta1) * Math.cos(fi1);
        rb[0][2] = Math.sin(fi1);

        rb[1][0] = Math.sin(delta1);
        rb[1][1] = -Math.cos(delta1);
        rb[1][2] = 0;

        rb[2][0] = Math.cos(delta1) * Math.sin(fi1);
        rb[2][1] = Math.sin(delta1) * Math.sin(fi1);
        rb[2][2] = Math.cos(fi1);

        rbt = getTranspose(rb);

        s = getMulti(rs, rbt);
        s = getMulti(ss, s);
        s = getMulti(rst, s);
        s = getMulti(rb, s);
    }
    public double[][] getTranspose(double[][] matrix3X3) {
        double tmp[][] = new double[3][3];
//        for (double[] failed1 : tmp) {
//            failed1[0] = Double.NaN;
//            failed1[1] = Double.NaN;
//            failed1[2] = Double.NaN;
//        }
        for (int c = 0; c < 3; c++) {
            for (int d = 0; d < 3; d++) {
                tmp[d][c] = matrix3X3[c][d];
            }
        }
        return tmp;
    }
    public double[][] getMulti(double[][] matrix3X3_a, double[][] matrix3X3_b) {
        double tmp[][] = new double[3][3];
//        for (double[] failed1 : tmp) {
//            failed1[0] = Double.NaN;
//            failed1[1] = Double.NaN;
//            failed1[2] = Double.NaN;
//        }
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                tmp[i][j] = 0;
                for (int k = 0; k < 3; k++) {
                    tmp[i][j] = tmp[i][j] + (matrix3X3_a[i][k] * matrix3X3_b[k][j]);
                }
            }
        }
        return tmp;
    }
   private void _computeStresses(int i, double m) {

        sigmas.rr = incp;
        sigmas.tt = s[0][0] + s[1][1] - 2 * (s[0][0] - s[1][1]) * Math.cos(2 * m) - 4 * s[0][1] * Math.sin(2 * m) - incp;
        sigmas.zz = s[2][2] - 2 * poissonratio * (s[0][0] - s[1][1]) * Math.cos(2 * m) - 4 * poissonratio * s[0][1] * Math.sin(2 * m);
        sigmas.tz = 2 * (s[1][2] * Math.cos(m) - s[0][2] * Math.sin(m));

        sigmas.tmax = (sigmas.zz + sigmas.tt + Math.sqrt(Math.pow(sigmas.zz - sigmas.tt, 2) + 4 * Math.pow(sigmas.tz, 2))) / 2;
        sigmas.tmin = (sigmas.zz + sigmas.tt - Math.sqrt(Math.pow(sigmas.zz - sigmas.tt, 2) + 4 * Math.pow(sigmas.tz, 2))) / 2;
        sigmas.w = 2 * sigmas.tz / (sigmas.zz - sigmas.tt);
        //New Changes
        sigmas.w = Math.atan(sigmas.w) / 2;
    }
private void _getPrincipalStresses(int i) {

        sigmas.s1 = Math.max(sigmas.tmax, sigmas.tmin);
        sigmas.s1 = Math.max(sigmas.s1, sigmas.rr);

        sigmas.s3 = Math.min(sigmas.tmax, sigmas.tmin);
        sigmas.s3 = Math.min(sigmas.s3, sigmas.rr);

        sigmas.s2 = sigmas.tmax + sigmas.tmin + sigmas.rr - (sigmas.s1 + sigmas.s3);
     
    }
private void _calculatedFailed(int i) {
        //calculate MC2D failed
        if (calculationMode == MC2D_MODE) {
            F = sigmas.s1 - (fmuTi * sigmas.s3 + ucs);
        } else if (calculationMode == LADE_MODE) {
            lade_S1 = (ucs * (Math.sqrt(Math.pow(muT, 2) + 1) - muT)) / (2 * muT);
            lade_h = (4 * Math.pow(muT, 2) * (9 * Math.sqrt(Math.pow(muT, 2) + 1) - 7 * muT)) / (Math.sqrt(Math.pow(muT, 2) + 1) - muT);
            lade_tail = 27 + lade_h;
            lade_I1 = sigmas.s1 + sigmas.s2 + sigmas.s3 + (3 * lade_S1);
//            lade_I1 = sigmas.s1 + lade_S1 + sigmas.s2 + lade_S1 + sigmas.s3 + lade_S1;
            lade_I3 = (sigmas.s1 + lade_S1) * (sigmas.s2 + lade_S1) * (sigmas.s3 + lade_S1);
            ladeFailVal = Math.pow(lade_I1, 3) / lade_I3;
            F = ladeFailVal - lade_tail;
        } else if (calculationMode == MC3D_MODE) {
            mu_i_sin = Math.sin(Math.toRadians(mu));
            mohr_S1 = (ucs * (Math.sqrt(Math.pow(muT, 2) + 1) - muT)) / (2 * muT);
            mohr_3d_I1 = sigmas.s1 + sigmas.s2 + sigmas.s3;
            mohr_3d_I2 = -(sigmas.s1 * sigmas.s2 + sigmas.s2 * sigmas.s3 + sigmas.s1 * sigmas.s3);
//            double I3 = sigmas.s1 * sigmas.s2 * sigmas.s3;
            mohr_3d_P = mohr_3d_I1 / 3;
            mohr_3d_J = Math.pow((Math.sqrt(2 * (Math.pow(mohr_3d_I1, 2)) + 6 * mohr_3d_I2) / 3), 2);
            mohr_3d_J = Math.sqrt(mohr_3d_J * 3 / 2);
            mohr_3d_theta = Math.atan((2 * sigmas.s2 - sigmas.s1 - sigmas.s3) / (Math.sqrt(3) * (sigmas.s1 - sigmas.s3)));
            mohr_3d_Func_val_theta = mu_i_sin / (Math.cos(mohr_3d_theta) + Math.sin(mohr_3d_theta) * mu_i_sin / Math.sqrt(3));
            F = (mohr_3d_J - (mohr_3d_P + mohr_S1) * mohr_3d_Func_val_theta);

        }
        //System.out.println(F);
    }
   private void _detectChanges(int i) {
   
        if (F > 0) {
            failed[i][0] = 1;
          
            if (i > 0) {
                 //System.out.println(failed[i][0]+" "+failed[i-1][0]);
                if (failed[i][0] != failed[i - 1][0]) {
                    f1++;
                    failed[i][1] = f1;
                    if (maxVal < f1) {
                        maxVal = f1;
                     
                    }
                }

            }
        } else {
            failed[i][0] = 0;
            
            if (i > 0) {
                if (failed[i][0] != failed[i - 1][0]) {
                    f1++;
                    failed[i - 1][1] = f1;
                    if (maxVal < f1) {
                        
                        maxVal = f1;
                    }
                }

            }
        }

    }
       private void _computeAngles(int n, int iTmax, int iTmin, double inc2, int m) {
        int itMin2;
        int itMax2;
        if (m > 0) {
            for (int i = 0; i < FNodes.length; i++) {
                FNodes[i] = 0;
            }
            for (int i = 0; i < n; i++) {
                if (failed[i][1] != 0) {
                    FNodes[failed[i][1] - 1] = i;  // Our index starts with 0
                }
            }

            itMin2 = iTmin + (int) n / 2;
            if (itMin2 > n) {
//                itMin2 -= n;
                itMin2 = iTmin - (int) n / 2;  //Enric changes
            }
            itMax2 = iTmax + (int) n / 2;
            if (itMax2 > n) {
                itMax2 = itMax2 - (int) n / 2;
                //itMax2=itMax2-n;
            }

            // Detect case
            if (failed[0][0] == 1) {
                if (m == 4) {
                    int ierr = isInside(m, FNodes, 1, 2, iTmin, itMin2);
                    if (ierr == 1) {
                        BKOAngle = 0;
                        BKODir = Double.NaN;
                        TensANgle = computeAngles(1, 2, 3, 0, FNodes, n, inc2);
                        TensDir = iTmin * inc2;
                    } else {
                        BKOAngle = computeAngles(1, 2, 3, 0, FNodes, n, inc2);
                        BKODir = iTmax * inc2;
                        TensANgle = 0;
                        TensDir = Double.NaN;

                    }
                } else {
                    int ierr = isInside(m, FNodes, 1, 2, iTmin, itMin2);
                    if (ierr == 1) {
                        BKOAngle = computeAngles(3, 4, 7, 0, FNodes, n, inc2);
                        BKODir = iTmax * inc2;
                        TensANgle = computeAngles(1, 2, 5, 6, FNodes, n, inc2);
                        TensDir = iTmin * inc2;
                    } else {
                        BKOAngle = computeAngles(1, 2, 5, 6, FNodes, n, inc2);
                        BKODir = iTmax * inc2;
                        TensANgle = computeAngles(3, 4, 7, 0, FNodes, n, inc2);
                        TensDir = iTmin * inc2;

                    }
                }
            } else {
                if (m == 4) {
                    int ierr = isInside(m, FNodes, 0, 1, iTmin, itMin2);
                    if (ierr == 1) {
                        BKOAngle = 0;
                        BKODir = Double.NaN;// iTmin
                        TensANgle = computeAngles(0, 1, 2, 3, FNodes, n, inc2);
                        TensDir = iTmin * inc2;
                    } else {
                        BKOAngle = computeAngles(0, 1, 2, 3, FNodes, n, inc2);
                        BKODir = iTmax * inc2;
                        TensANgle = 0;
                        TensDir = Double.NaN;// itmin 
                    }
                } else {
                    int ierr = isInside(m, FNodes, 0, 1, iTmin, itMin2);
                    if (ierr == 1) {
                        BKOAngle = computeAngles(2, 3, 6, 7, FNodes, n, inc2);
                        BKODir = iTmax * inc2;
                        TensANgle = computeAngles(0, 1, 4, 5, FNodes, n, inc2);
                        TensDir = iTmin * inc2;
                    } else {
                        BKOAngle = computeAngles(0, 1, 4, 5, FNodes, n, inc2);
                        BKODir = iTmax * inc2;
                        TensANgle = computeAngles(2, 3, 6, 7, FNodes, n, inc2);
                        TensDir = iTmin * inc2;

                    }
                }
            }

        } else {
            if (failed[0][0] == 0) {
                BKOAngle = 0;
            } else {
                BKOAngle = 180;
            }
        }

    }

    private void _computeAngles2(int n, int iTmax, int iTmin, double inc2, int m) {
        int itMin2;
        int itMax2;
        int ncase = 2;
        int iDmax;
        int iDmax2;
        int ierr;
        int ierr2;
        int a, b;

        if (m >= 4) {
            FNodes = new int[m];
            for (int i = 0; i < FNodes.length; i++) {
                FNodes[i] = 0;
            }
            for (int i = 0; i < n; i++) {
                if (failed[i][1] != 0) {
                    FNodes[failed[i][1] - 1] = i;  // Our index starts with 0
                }
            }

            itMin2 = iTmin + (int) n / 2;
            if (itMin2 > n) {
                itMin2 = iTmin - (int) n / 2;  //Enric changes
            }

            itMax2 = iTmax + (int) n / 2;
            if (itMax2 > n) {
                itMax2 = iTmax - (int) n / 2;//itMax2=itMax2-n;
            }

            if (ncase == 1) {
                iDmax = iTmax;
                iDmax2 = itMax2;
            } else {
                iDmax = iTmin;
                iDmax2 = itMin2;

            }

            // Detect case
            if (failed[0][0] == 1) {
                /*
                 !m=8 then failed areas are 12,34,56,78
                 !m=4 then failed areas are 12,34
                 */
                if (m == 4) {
                    // !Test when Tensile area and breakout area are joined
                    ierr = isInside(m, FNodes, 1, 2, iTmax, itMax2);
                    ierr2 = isInside(m, FNodes, 1, 2, iTmin, itMin2);
                    if (ierr == 1 && ierr2 == 1) {
                        /*
                         !case D
                         !Get the two distances to BKO direction
                        
                         */
                        a = getValueInside(iTmax, itMax2, FNodes[1], FNodes[2]);
                        b = getValueInside(iTmin, itMin2, FNodes[1], FNodes[2]);
                        //     !Compute angles
                        if (a > b) { //!case D.1.1
                            BKOAngle = inc2 * 2 * FNodes[2] - a;
                            BKODir = iTmax * inc2;
                            TensANgle = inc2 * 2 * b - FNodes[1];
                            TensDir = iTmin * inc2;
                        } else {
                            BKOAngle = inc2 * 2 * a - FNodes[1];
                            BKODir = iTmax * inc2;
                            TensANgle = inc2 * 2 * FNodes[2] - b;
                            TensDir = iTmin * inc2;

                        }

                    } else {
                        //!Test if iTmin or iTmin2 are in the interval of 1-2 nodes
                        ierr = isInside(m, FNodes, 1, 2, iDmax, iDmax2);
                        if (ierr == 1) {
                            //!The areas are tensile areas, then return angles
                            if (ncase == 1) {
                                BKOAngle = computeAngles(1, 2, 3, 0, FNodes, n, inc2);
                                BKODir = iTmax * inc2;
                                TensANgle = 0;
                                TensDir = iTmin * inc2;

                            } else {
                                BKOAngle = 0;
                                BKODir = iTmax * inc2;
                                TensANgle = computeAngles(1, 2, 3, 0, FNodes, n, inc2);
                                TensDir = iTmin * inc2;

                            }
                        } else {
                            if (ncase == 1) {
                                //  !The areas are tensile areas, then return angles
                                BKOAngle = computeAngles(1, 2, 3, 0, FNodes, n, inc2);
                                BKODir = iTmax * inc2;
                                TensANgle = 0;
                                TensDir = iTmax * inc2;

                            } else {
                                // !The areas are tensile areas, then return angles
                                BKOAngle = computeAngles(1, 2, 3, 0, FNodes, n, inc2);
                                BKODir = iTmax * inc2;
                                TensANgle = 0;
                                TensDir = iTmax * inc2;

                            }
                        }
                    }
                } else {
                    /*
                     !m=8
                     !Test if iTmin or iTmin2 are in the interval of 1-2 nodes
                     */
                    ierr = isInside(m, FNodes, 1, 2, iTmin, itMin2);
                    if (ierr == 1) {
                        if (ncase == 1) {
                            // !The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(1, 2, 5, 6, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(3, 4, 7, 0, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        } else {
                            //!The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(3, 4, 7, 0, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(1, 2, 5, 6, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        }

                    } else {
                        if (ncase == 1) {
                            // !The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(3, 4, 7, 0, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(1, 2, 5, 6, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        } else {
                            //!The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(1, 2, 5, 6, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(3, 4, 7, 0, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        }
                    }
                }
            } else {
                /*
                 !m=4 then failed areas are 23,41
                 !m=8 then failed areas are 23,45,67,81
                 */
                if (m == 4) {
                    // !Test when Tensile area and breakout area are joined
                    ierr = isInside(m, FNodes, 0, 1, iTmax, itMax2);
                    ierr2 = isInside(m, FNodes, 0, 1, iTmin, itMin2);

                    if (ierr == 1 && ierr2 == 1) {
                        /*
                         !case D
                         !Get the two distances to BKO direction                        
                         */
                        a = getValueInside(iTmax, itMax2, FNodes[0], FNodes[1]);
                        b = getValueInside(iTmin, itMin2, FNodes[0], FNodes[1]);
                        //     !Compute angles
                        if (a > b) {
                            BKOAngle = inc2 * 2 * (FNodes[1] - a);
                            BKODir = iTmax * inc2;// iTmin
                            TensANgle = inc2 * 2 * (b - FNodes[0]);
                            TensDir = iTmin * inc2;
                        } else {
                            BKOAngle = inc2 * 2 * (a - FNodes[0]);
                            BKODir = iTmax * inc2;// iTmin
                            TensANgle = inc2 * 2 * (FNodes[1] - b);
                            TensDir = iTmin * inc2;
                        }
                    } else {
                        //!Test if iTmin or iTmin2 are in the interval of 1-2 nodes
                        ierr = isInside(m, FNodes, 0, 1, iDmax, iDmax2);
                        if (ierr == 1) {
                            if (ncase == 1) {
                                // !The areas are tensile areas, then return angles
                                BKOAngle = computeAngles(0, 1, 2, 3, FNodes, n, inc2);
                                BKODir = iTmax * inc2;
                                TensANgle = 0;
                                TensDir = iTmin * inc2;

                            } else {
                                // !The areas are tensile areas, then return angles
                                BKOAngle = 0;
                                BKODir = iTmax * inc2;
                                TensANgle = computeAngles(0, 1, 2, 3, FNodes, n, inc2);
                                TensDir = iTmin * inc2;

                            }
                        } else {
                            if (ncase == 1) {
                                // !The areas are tensile areas, then return angles
                                BKOAngle = 0;
                                BKODir = iTmax * inc2;
                                TensANgle = computeAngles(0, 1, 2, 3, FNodes, n, inc2);
                                TensDir = iTmin * inc2;

                            } else {
                                // !The areas are tensile areas, then return angles
                                BKOAngle = computeAngles(0, 1, 2, 3, FNodes, n, inc2);
                                BKODir = iTmax * inc2;
                                TensANgle = 0;
                                TensDir = iTmin * inc2;

                            }
                        }

                    }
                } else {
                    /*
                     !m=8
                     !Test if iTmin or iTmin2 are in the interval of 1-2 nodes
                     */
                    ierr = isInside(m, FNodes, 0, 1, iDmax, iDmax2);
                    if (ierr == 1) {
                        if (ncase == 1) {
                            // !The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(0, 1, 4, 5, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(2, 3, 6, 7, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        } else {
                            // !The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(2, 3, 6, 7, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(0, 1, 4, 5, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        }

                    } else {
                        if (ncase == 1) {
                            // !The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(2, 3, 6, 7, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(0, 1, 4, 5, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        } else {
                            // !The areas are tensile areas, then return angles
                            BKOAngle = computeAngles(0, 1, 4, 5, FNodes, n, inc2);
                            BKODir = iTmax * inc2;
                            TensANgle = computeAngles(2, 3, 6, 7, FNodes, n, inc2);
                            TensDir = iTmin * inc2;

                        }

                    }
                }
            }

        } else {
            if (failed[0][0] == 0) {
                BKOAngle = 0;
                BKODir = iTmax * inc2;
                TensANgle = 0;
                TensDir = iTmin * inc2;
            } else {
                BKOAngle = 180;
                BKODir = iTmax * inc2;
                TensANgle = 0;
                TensDir = iTmin * inc2;
            }
        }

    }

    private double computeAngles(int node1, int node2, int node3, int node4, int[] fNodes, int n, double inc) {
        double computeAngle = Double.NaN;
        double angle1 = Double.NaN;
        double angle2 = Double.NaN;

        angle1 = (fNodes[node2] - fNodes[node1]) * inc;

        if (fNodes[node4] < fNodes[node3]) {
            angle2 = (fNodes[node4] + (n - fNodes[node3])) * inc;
        } else {
            angle2 = (fNodes[node4] - fNodes[node3]) * inc;
        }

        // compute test angle
        if (Math.abs(angle1 - angle2) > inc) {
            computeAngle = angle1;
        } else {
            computeAngle = (angle1 + angle2) / 2;
        }
        //System.out.println(computeAngle);
        return computeAngle;
    }

   public void _calculateMudWeight(double BKO, int mwSign) {
//        if (Math.abs(BKO - firstBKO) < angleTolerance * 2) {
//            cont = false;
//        } else if (firstBKO > BKO) {
        if (Math.abs(BKO - BKOAngle) < angleTolerance * 2) {
            cont = false;
        } else if (BKOAngle > BKO) {

            if (mwSign == -1) {
                inc = inc / divFactorForIncrDecrFactor;
                mwSign = 1;
                ch++;
            }
            mudWeight = mudWeight + inc;
        } else {

            if (mwSign == 1) {
                inc = inc / divFactorForIncrDecrFactor;
                mwSign = -1;
                ch++;
            }
            mudWeight = mudWeight - inc;
        }
        if (Math.abs(previousMW - mudWeight) < tol) {
            cont = false;

        } else if (ch > changes) {
            cont = false;
        }
        previousMW = mudWeight;
        count++;
        
    }
     public void _resetInnerVariables(int n) {
        incp = mudWeight - pp;
        f1 = 0;
        failed = new int[n][2];
        for (int[] failed1 : failed) {
            failed1[0] = 0;
            failed1[1] = 0;
        }
    }
     public void _reAssignVariables() {
        cont = true;
        changes = 5;
        ch = 0;
        previousMW = Double.NaN;
        inc = mwIncrimentFactor;
//        inc2 = 0.5;
//        tol = 0.0005;
        f1 = 0;
        mudWeight = pp;
        mu = internalfriction;
        muT = Math.tan(Math.toRadians(mu));
        fmuTi = Math.pow(Math.sqrt(Math.pow(muT, 2) + 1) + muT, 2);
        incp = Double.NaN;
        F = Double.NaN;
        maxVal = 0;
        count = 0;
    }
     public int isInside(int nodes, int[] fNodes, int node1, int node2, int dir1, int dir2) {

        if (((fNodes[node1] <= dir1) && (fNodes[node2] >= dir1)) || ((fNodes[node1] <= dir2) && (fNodes[node2] >= dir2))) {
//        if (((fNodes[node1] >= dir1) && (fNodes[node2] <= dir1)) || ((fNodes[node1] >= dir2) && (fNodes[node2] <= dir2))) {
            return 1;
        } else {
            return 0;
        }
    }
     public int getValueInside(int iTmax, int itMax2, int low, int height) {
        int getValueInside;
        if ((iTmax >= low) && (iTmax <= height)) {
            getValueInside = iTmax;
        } else if ((itMax2 >= low) && (itMax2 <= height)) {
            getValueInside = itMax2;
        } else {
            getValueInside = -1;

        }
        return getValueInside;
    }
     class Sigma{
         double tt;
        double tz;
        double rr;
        double zz;
        double tmax;
        double tmin;
        double w;
        double s1;
        double s2;
        double s3;
     }
}

