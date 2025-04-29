/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package CEC2022;

import java.io.BufferedReader;  
import java.io.FileReader;  
import java.io.IOException;  
import java.util.ArrayList;
import java.util.Arrays; 
import java.util.List;

/**
 *
 * @author farzaneh
 */
public class cec22fit {
    int func_num;
    public cec22fit(int function){
        this.func_num = function;
    }
    
    private static double[] sr_func(double[] x, int nx, double[] os, double[][] Mr, double sh_rate){
        double[] xshift = shiftfunc(x, nx, os);
        for (int i = 0; i < nx; i++) {  
                    xshift[i] = xshift[i] * sh_rate;  
        } 
        double[] rotatedX = rotatefunc(xshift, nx, Mr); 
        return rotatedX;
    }
    
    private static double[] shiftfunc(double[] x, int nx, double[] Os) {  
        double[] xshift = new double[nx];  
        for (int i = 0; i < nx; i++) {  
            xshift[i] = x[i] - Os[i];  
        }  
        return xshift;  
    }  
    
    private static double[] rotatefunc(double[] x, int nx, double[][] Mr) {  
        double[] xrot = new double[nx];  
        for (int i = 0; i < nx; i++) {  
            for (int j = 0; j < nx; j++) {  
                xrot[i] += x[j] * Mr[i][j];  
            }  
        }  
        return xrot;  
    }  
    private static double zakharov_func(double[] x, int nx, double[] Os, double[][] Mr) {     //F1
            double[] z = sr_func(x, nx, Os, Mr, 1.0);  
            double sum1 = 0.0;  
            double sum2 = 0.0; 
            for (int i = 0; i < nx; i++) {  
                sum1 += Math.pow(z[i], 2);  
                sum2 += 0.5 * z[i];  
            }  
            return sum1 + Math.pow(sum2, 2) + Math.pow(sum2, 4);   
        }      
    
    private static double rosenbrock_func(double[] x, int nx, double[] Os, double[][] Mr) {     //F2
        double f = 0.0;  
        double[] z = sr_func(x, nx, Os, Mr, 2.048 / 100.0); 
        for(int j=0; j<z.length; j++){
            z[j] += 1;
        }
        for (int i = 0; i < nx - 1; i++) {  
            double tmp1 = z[i] * z[i] - z[i + 1];  
            double tmp2 = z[i+1] - 1.0;  
            f += (100.0 * Math.pow(tmp1, 2)) + Math.pow(tmp2, 2);  
        }  
        return f;  
    }
   // Shifted and full Rotated Expanded Schaffer’s F7     
   public static double schaffer_F7_func(double[] x, int nx, double[] Os, double[][] Mr) {    // F3            Multi-modal
        double[] z = sr_func(x, nx, Os, Mr, 0.5/100.0);    
        double sum = 0.0;
            for (int i = 0; i < nx-1; i++) {
                double numerator = Math.pow(Math.sin(Math.sqrt(z[i] * z[i] + z[i + 1] * z[i + 1])), 2)   - 0.5;
                double denominator = Math.pow(1 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]),2);
                sum += 0.5 + (numerator / denominator) ;

            }
            int d = nx-1 ;
            double nume1 = Math.pow(Math.sin(Math.sqrt(z[d] * z[d] + z[0] * z[0])), 2)   - 0.5;
            double den1 = Math.pow(1 + 0.001 * (z[d] * z[d] + z[0] * z[0]), 2) ;
            double sum1 = 0.5 + (nume1 / den1) ;  // g(x[d], x[1])
            sum += sum1; 
          return    sum;  
           
    }  
    //   Shifted and Rotated Non-Continuous Rastrigin’s Function      F4
    public static double rastrigin_func(double[] x, int nx, double[] Os, double[][] Mr) {  //Multi-modal
        double fit = 0.0;  
        double[] z = sr_func(x, nx, Os, Mr, 5.12 / 100.0);  
        for (int i = 0; i < nx; i++) {
            fit += Math.pow(z[i] , 2) - (10 * Math.cos(2 * Math.PI * z[i])) + 10;
        }
        return fit;  
    } 
    //   Shifted and Rotated Levy Function
    public static double levy_func(double[] x, int nx, double[] Os, double[][] Mr) {  //  F5 
        double[] z = sr_func(x, nx, Os, Mr, 5.12/100.0);
        double[] w = new double[nx];   
        for (int i = 0; i < nx; i++) {  
            w[i] = (z[i] - 1) / 4.0 + 1.0 ;  
        }  
        double term1 = Math.pow((Math.sin(Math.PI * w[0])), 2);  
        double term3 = Math.pow((w[nx-1] - 1), 2) * (1 + Math.pow((Math.sin(2 * Math.PI * w[nx - 1])), 2));  
        double Sum = 0.0;  

        for (int i = 0; i < nx - 1; i++) {  
            double newv = Math.pow((w[i] - 1), 2) * (1 + 10 * Math.pow((Math.sin(Math.PI * w[i] - 1)), 2));  
            Sum = Sum + newv;  
        }  
        return term1 + Sum + term3;  
          
    }  
    public static double Hybrid_Function_1(double[] x, int nx, double[] Os, double[][] Mr, double[] S) {  //f6
        int cf_num = 3;  
        double[] fit = new double[3];  
        int[] G = new int[3];  
        int[] G_nx = new int[3];  
        double[] Gp = {0.4, 0.4, 0.2};
        double tmp = 0.0;  
        for (int i = 0; i < cf_num - 1; i++) {  
            G_nx[i] = (int) Math.ceil(Gp[i] * nx);  
            tmp += G_nx[i];  
        }  
        G_nx[cf_num - 1] = nx - (int) tmp;
        G[0] = 0;  
        for (int i = 1; i < cf_num; i++) {  
            G[i] = G[i - 1] + G_nx[i - 1];  
        }  
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        int[] sInt = new int[S.length];  
        for (int i = 0; i < S.length; i++) {  
            sInt[i] = (int) S[i] - 1;  
        }  
        double [] y = new double[nx];  
        for (int i = 0; i < nx; i++) {  
            y[i] = z[sInt[i]];  
        }
        int i = 0;  
        fit[0] = bent_cigar_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 1;  
        fit[i] = hgbat_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 2;  
        fit[i] = rastrigin_func(Arrays.copyOfRange(y, G[i], nx), G_nx[i], Os, Mr);  
        double f = 0.0;  
        for (i = 0; i < cf_num; i++) {  
            f += fit[i];  
        }  
        return f;  
    }
    public static double happycat_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double alpha = 1.0 / 8.0;  
        double[] z = sr_func(x, nx, Os, Mr, 5.0 / 100.0);  
        double r2 = 0.0;  
        double sum_z = 0.0;
        for (int i = 0; i < nx; i++) {  
            z[i] = z[i] - 1.0;  
            r2 += z[i] * z[i];  
            sum_z += z[i];  
        }
        double f = Math.pow(Math.abs(r2 - nx), 2 * alpha) + (0.5 * r2 + sum_z) / nx + 0.5;  
        return f;  
    }
    //  Expanded Rosenbrock’s plus Griewangk’s Function
    public static double grie_rosen_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double sol1 = 0.0;  double sol2 = 0.0; 
        double[] z = sr_func(x, nx, Os, Mr, 1.0);
        for(int j=0; j<z.length; j++){
            z[j] += 1;
        }
        for (int i = 0; i < nx - 1; i++) {  
            double tmp1 = Math.pow(z[i], 2) - z[i + 1];  
            double tmp2 = z[i+1] - 1.0;  
            sol1 += (100.0 * Math.pow(tmp1, 2)) + Math.pow(tmp2, 2);
            sol2 += (sol1*sol1)/4000.0-Math.cos(sol1) +1 ;
        }  
        double tmp1 = Math.pow(z[nx-1], 2) - z[0];  
        double tmp2 = z[nx - 1] - 1.0;  
        double temp = (100.0 * Math.pow(tmp1, 2)) + Math.pow(tmp2, 2);
        sol2 += (temp * temp) / 4000.0 - Math.cos(temp) + 1.0; 
        return sol2;  
    }  
//  F8
    public static double Hybrid_Function_3(double[] x, int nx, double[] Os, double[][] Mr, double[] S) {  
        int cf_num = 5;  
        double[] fit = new double[cf_num];  
        double[] G = new double[cf_num];  
        double[] G_nx = new double[cf_num];  
        double[] Gp = {0.3, 0.2, 0.2, 0.1, 0.2}; 
        double tmp = 0;  
        for (int i = 0; i < cf_num - 1; i++) {  
            G_nx[i] = Math.ceil(Gp[i] * nx);  
            tmp += G_nx[i];  
        }  
        G_nx[cf_num - 1] = nx - tmp;  
        int[] G_nx_int = new int[cf_num];  
        for (int i = 0; i < cf_num; i++) {  
            G_nx_int[i] = (int) G_nx[i]; // Explicit cast to int  
        } 
        G[0] = 0;  
        for (int i = 1; i < cf_num; i++) {  
            G[i] = G[i - 1] + G_nx_int[i - 1];  
        }
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        double []   y = new double[nx];  
        // Convert S to int[] and adjust indices (subtract 1)  
        int[] S_int = new int[nx];  
        for (int i = 0; i < nx; i++) {  
            S_int[i] = (int) S[i] - 1; // Subtract 1 to adjust indices  
        }  
        for (int i = 0; i < nx; i++) {  
            y[i] = z[S_int[i]];  
        } 
        // Extract subarrays using Arrays.copyOfRange and call the functions  
        int i = 0;  
        fit[i] = katsuura_func(Arrays.copyOfRange(y, (int)G[i], (int)G[i + 1]), G_nx_int[i], Os, Mr);  
        i = 1;  
        fit[i] = happycat_func(Arrays.copyOfRange(y, (int)G[i], (int)G[i + 1]), G_nx_int[i], Os, Mr);  
        i = 2;  
        fit[i] = grie_rosen_func(Arrays.copyOfRange(y, (int)G[i], (int)G[i + 1]), G_nx_int[i], Os, Mr);  
        i = 3;  
        fit[i] = schwefel_func(Arrays.copyOfRange(y, (int)G[i], (int)G[i + 1]), G_nx_int[i], Os, Mr);  
        i = 4;  
        fit[i] = ackley_func(Arrays.copyOfRange(y, (int)G[i], nx), G_nx_int[i], Os, Mr);  

        double f = 0.0;  
        for (i = 0; i < cf_num; i++) {  
            f += fit[i];  
        }  
        return f;  
    }  

    public static double katsuura_func(double[] x, int nx, double[] Os, double[][] Mr) {   //f9
        double f = 1.0; 
        double term1 = 10/Math.pow( nx, 2);
        double term2 = 10/Math.pow( nx, 1.2);
        double[] z = sr_func(x, nx, Os, Mr, 1.0);
        for (int i = 0; i < nx; i++) {  
            double temp = 0.0;  
            for (int j = 0; j < 32; j++) {    
                double tmp1 = Math.pow(2.0, j);  
                double tmp2 = tmp1 * z[i]; 
                double tmp3 = Math.abs(tmp2 - Math.round(tmp2 ));
                temp +=  tmp3/ tmp1;  
            }  
            double tmp4 = 1+(i*temp);
            double tmp5 = Math.pow(tmp4, term2);
            f *= tmp5;
        }  
        double tmp1 = term1 * f - term1; 
        
        return tmp1;  
    } 
    public static double ackley_func(double[] x, int nx, double[] Os, double[][] Mr) {   //f13
        double sum1 = 0.0;  
        double sum2 = 0.0;
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        for (int i = 0; i < nx; i++) {  
            sum1 += Math.pow(z[i], 2) ;  
            sum2 += Math.cos(2.0 * Math.PI * z[i]);  
        } 
        double sum3 = -0.2 * Math.sqrt(sum1 / nx);  
        double sum4 = -20 * Math.exp(sum3);
        double sum5 = Math.exp(sum2);
        double fit = sum4 - sum5 +20 + Math.E; 
        return fit;  
    }  
    // F7
    public static double Hybrid_Function_2(double[] x, int nx, double[] Os, double[][] Mr, double[] S) {  
        int cf_num = 6;  
        double[] fit = new double[6];  
        int[] G = new int[6];  
        int[] G_nx = new int[6];  
        double[] Gp = {0.1, 0.2, 0.2, 0.2, 0.1, 0.2};
        double tmp = 0.0;  
        for (int i = 0; i < cf_num - 1; i++) {  
            G_nx[i] = (int) Math.ceil(Gp[i] * nx);  
            tmp += G_nx[i];  
        }  
        G_nx[cf_num - 1] = nx - (int) tmp;  

        G[0] = 0;  
        for (int i = 1; i < cf_num; i++) {  
            G[i] = G[i - 1] + G_nx[i - 1];  
        }  
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        int[] sInt = new int[S.length];  
        for (int i = 0; i < S.length; i++) {  
            sInt[i] = (int) S[i] - 1;  
        }  
        double [] y = new double[nx];  
        for (int i = 0; i < nx; i++) {  
            y[i] = z[sInt[i]];  
        }  
        int i = 0;  
        fit[i] = hgbat_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 1;  
        fit[i] = katsuura_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 2;  
        fit[i] = ackley_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 3;  
        fit[i] = rastrigin_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 4;  
        fit[i] = schwefel_func(Arrays.copyOfRange(y, G[i], G[i + 1]), G_nx[i], Os, Mr);  
        i = 5;  
        fit[i] = schaffer_F7_func(Arrays.copyOfRange(y, G[i], nx), G_nx[i], Os, Mr);  

        double f = 0.0;  
        for (i = 0; i < cf_num; i++) {  
            f += fit[i];  
        }  
        return f;  
    }  

    
    //High Conditioned Elliptic Function
    public static double ellips_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double fit = 0.0;  
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        for (int i = 0; i < nx; i++) { 
            double term1 = Math.pow(10.0, 6.0);
            double term2 = (i-1)/(nx-1);
            double term3 = Math.pow(z[i], 2);
            fit += Math.pow(term1, term2) * term3;
        }  
        return fit;  
    }
    //  Bent Cigar Function
    public static double bent_cigar_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        double term1 = z[0] * z[0]; 
        double term2 = Math.pow(10.0, 6.0);
        double f=0;
        for (int i = 1; i < nx; i++) {  
            f += z[i] * z[i];  
        }  
        return term1 + term2 * f;  
    }  
    public static double discus_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        double f = Math.pow(10.0, 6.0) * z[0] * z[0];  
        for (int i = 1; i < nx; i++) {  
            f += z[i] * z[i];  
        }  
        return f;  
    } 
    
    
    private static double[] subArray(double[] array, int startIndex, int endIndex) {  
    int length = endIndex - startIndex;  
    double[] result = new double[length];  
    System.arraycopy(array, startIndex, result, 0, length);  
    return result;  
    }  
    private static double[][] subMatrix(double[][] matrix, int startRow, int endRow, int startCol, int endCol) {  
        int numRows = endRow - startRow;  
        int numCols = endCol - startCol;  
        double[][] result = new double[numRows][numCols];  
        for (int i = 0; i < numRows; i++) {  
            System.arraycopy(matrix[startRow + i], startCol, result[i], 0, numCols);  
        }  
        return result;  
    } 
    

public static double cf_cal(double[] x, int nx, double[] Os, double[] delta, double[] bias, double[] fit, int cf_num) {  
    double w_max = 0;  
    double w_sum = 0;  
    double[] w = new double[cf_num];  
    double INF = Double.MAX_VALUE;
    for (int i = 0; i < cf_num; i++) {  
        fit[i] += bias[i];  
        w[i] = 0;  
        for (int j = 0; j < nx; j++) {  
            w[i] += Math.pow(x[j] - Os[i * nx + j], 2.0);  
        }  
        if (w[i] != 0) {
            
            w[i] = Math.pow(1.0 / w[i], 0.5) * Math.exp(-w[i] / (2.0 * nx * Math.pow(delta[i], 2.0)));  
        } else {  
            w[i] = INF;  
        }  
        if (w[i] > w_max) {  
            w_max = w[i];  
        }  
    }  

    for (int i = 0; i < cf_num; i++) {  
        w_sum = w_sum + w[i];  
    }  
    if (w_max == 0) {  
        for (int i = 0; i < cf_num; i++) {  
            w[i] = 1;  
        }  
        w_sum = cf_num;  
    }  

    double f = 0.0;  
    for (int i = 0; i < cf_num; i++) {  
        f = f + (w[i] / w_sum) * fit[i];  
    }    
    return f;  
}  
  //Rotated Rosenbrock’s Function f2+High Conditioned Elliptic Function f8+Rotated Bent Cigar Function f6+
// Rotated Discus Function f14+High Conditioned Elliptic Function f8
    public static double Composition_Function_1(double[] x, int nx, double[] Os, double[][] Mr) { //F9
       int cf_num=5;  
       double[] fit = new double[cf_num];  
       double[] delta = {10, 20, 30, 40, 50};
       double [] bias = {0, 200, 300, 100, 400};
       int i = 0;  
        fit[i] = rosenbrock_func(x, nx, subArray(Os, i * nx, (i + 1) * nx),  subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        fit[i] =  fit[i] *1;  
        i = 1;  
        fit[i] = ellips_func(x, nx, subArray(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        fit[i] =  fit[i] * 1e-6;  
        i = 2;  
        fit[i] = bent_cigar_func(x, nx, subArray(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        fit[i] = fit[i] * 1e-6;    
        i = 3;  
        fit[i] = discus_func(x, nx, subArray(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        fit[i] = fit[i] * 1e-6;   
        i = 4;  
        fit[i] = ellips_func(x, nx, subArray(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        fit[i] = fit[i] * 1e-6;   
        
        double f = cf_cal(x, nx, Os, delta, bias, fit, cf_num);  
        return f;   
    }
    
    public static double schwefel_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double fit = 0.0;
        double gz =0 ;
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  

        for (int i = 0; i < nx; i++) {  
            z[i] += 420.9687462275036 ; 
            if (z[i] > 500) {  
                double term1 = Math.IEEEremainder(z[i], 500);
                double term3 = (500.0 - term1) * Math.sin(Math.sqrt(Math.abs(term1)));  
                double tmp = Math.pow((z[i] - 500.0),2)/(10000*nx) ;  
                 gz = term3 - tmp;  
                
            } else if (z[i] < -500) {  
                double term1 = Math.IEEEremainder(Math.abs(z[i]), 500);
                double term3 = (term1- 500.0) *Math.sin(Math.abs(term1));
                double tmp = Math.pow((z[i] + 500.0),2)/(10000*nx) ; 
                gz = term3 - tmp;
            } else if (Math.abs(z[i]) <= 500) {  
                double term1 = z[i] * Math.sin(Math.pow(Math.abs(z[i]), 0.5));  
                gz = term1;
            } 
            fit =418.9829 * nx - gz ; 
        }  
         
        return fit;  
    } 
  
    //  HGBat Function
    public static double hgbat_func(double[] x, int nx, double[] Os, double[][] Mr) {   
        double[] z = sr_func(x, nx, Os, Mr, 5.0 / 100.0); 
        double term1 =0;
        double term2 =0;
        for(int i=0; i< nx; i++){
            term1 += Math.pow(z[i], 2);
            term2 += z[i];
        }
        double cal1 = Math.pow(term1, 2)- Math.pow(term2, 2);
        double cal2 = Math.pow(Math.abs(cal1), 0.5) ;
         double cal3 = ((term1 * 0.5)+ term2) / nx;
        
        double f = cal2 + cal3 + 0.5; 
        return f;  
    } 
    //F10
    public static double Composition_Function_2(double[] x, int nx, double[] Os, double[][] Mr) {  
        int cf_num = 3;  
        double[] fit = new double[3];  
        double[] delta = {20, 10, 10};  
        double[] bias = {0, 200, 100};  

        int i = 0;  
        fit[i] = schwefel_func(x, nx, Arrays.copyOfRange(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        i = 1;  
        fit[i] = rastrigin_func(x, nx, Arrays.copyOfRange(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, i * nx + nx, 0, nx));  
        i = 2;  
        fit[i] = hgbat_func(x, nx, Arrays.copyOfRange(Os, i * nx, (i + 1) * nx), subMatrix(Mr, i * nx, (i + 1) * nx, 0, nx));  
        double f = cf_cal(x, nx, Os, delta, bias, fit, cf_num);  
        return f;  
    }  
    public static double escaffer6_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double[] z = sr_func(x, nx, Os, Mr, 1.0);  
        double f = 0.0;
        for (int i = 0; i < nx - 1; i++) {  
            double temp1 = Math.sin(Math.sqrt(z[i] * z[i] + z[i + 1] * z[i + 1]));  
            temp1 = temp1 * temp1;  
            double temp2 = 1.0 + 0.001 * (z[i] * z[i] + z[i + 1] * z[i + 1]);  
            f += 0.5 + (temp1 - 0.5) / (temp2 * temp2);  
        }
        double temp1 = Math.sin(Math.sqrt(z[nx - 1] * z[nx - 1] + z[0] * z[0]));  
        temp1 = temp1 * temp1;  
        double temp2 = 1.0 + 0.001 * (z[nx - 1] * z[nx - 1] + z[0] * z[0]);  
        f += 0.5 + (temp1 - 0.5) / (temp2 * temp2); 
        return f;  
    }
    public static double griewank_func(double[] x, int nx, double[] Os, double[][] Mr) {  
        double s = 0.0;  
        double p = 1.0; 
        double[] z = sr_func(x, nx, Os, Mr, 600.0 / 100.0);  
        for (int i = 0; i < nx; i++) {  
            s += z[i] * z[i];  
            p *= Math.cos(z[i] / Math.sqrt(1.0 + i));  
        }  

        double f = 1.0 + s / 4000.0 - p;  
        return f;  
    }  
//F11
    public static double Composition_Function_3(double[] x, int nx, double[] Os, double[][] Mr) {  
        int cf_num = 5;  
        double[] fit = new double[5];  
        double[] delta = {20, 20, 30, 30, 20};  
        double[] bias = {0, 200, 300, 400, 200};    
        fit[0] = escaffer6_func(x, nx,Arrays.copyOfRange(Os, 0 * nx, (0 + 1) * nx),subMatrix(Mr, 0 * nx, (0 + 1) * nx, 0, nx));  
        fit[0] = 10000 * fit[0] / 2e+7;   
        fit[1] = schwefel_func(x, nx,Arrays.copyOfRange(Os, 1 * nx, (1 + 1) * nx), subMatrix(Mr, 1 * nx, (1 + 1) * nx, 0, nx));   
        fit[2] = griewank_func(x, nx, Arrays.copyOfRange(Os, 2 * nx, (2 + 1) * nx),subMatrix(Mr, 2 * nx, (2 + 1) * nx, 0, nx));  
        fit[2] = 1000 * fit[2] / 100;   
        fit[3] = rosenbrock_func(x, nx, Arrays.copyOfRange(Os, 3 * nx, (3 + 1) * nx), subMatrix(Mr, 3 * nx, (3 + 1) * nx, 0, nx));  
        fit[4] = rastrigin_func(x, nx,Arrays.copyOfRange(Os, 4 * nx, (4 + 1) * nx), subMatrix(Mr, 4 * nx, (4 + 1) * nx, 0, nx));  
        fit[4] = 10000 * fit[4] / 1e+3;
        double f = cf_cal(x, nx, Os, delta, bias, fit, cf_num);  
        return f;  
    }
    //F12
    public static double Composition_Function_4(double[] x, int nx, double[] Os, double[][] Mr) {  
    int cf_num = 6;  
    double[] fit = new double[6];  
    double[] delta = {10, 20, 30, 40, 50, 60};  
    double[] bias = {0, 300, 500, 100, 400, 200};  
    fit[0] = hgbat_func(x, nx, Arrays.copyOfRange(Os, 0 * nx, (0 + 1) * nx),subMatrix(Mr, 0 * nx, (0 + 1) * nx, 0, nx));  
    fit[0] = 10000 * fit[0] / 1000;
    fit[1] = rastrigin_func(x, nx,Arrays.copyOfRange(Os, 1 * nx, (1 + 1) * nx), subMatrix(Mr, 1 * nx, (1 + 1) * nx, 0, nx));  
    fit[1] = 10000 * fit[1] / 1e+3;  
    fit[2] = schwefel_func(x, nx, Arrays.copyOfRange(Os, 2 * nx, (2 + 1) * nx), subMatrix(Mr, 2 * nx, (2 + 1) * nx, 0, nx));  
    fit[2] = 10000 * fit[2] / 4e+3;  
    fit[3] = bent_cigar_func(x, nx, Arrays.copyOfRange(Os, 3 * nx, (3 + 1) * nx), subMatrix(Mr, 3 * nx, (3 + 1) * nx, 0, nx));  
    fit[3] = 10000 * fit[3] / 1e+30;  
    fit[4] = ellips_func(x, nx,Arrays.copyOfRange(Os, 4 * nx, (4 + 1) * nx), subMatrix(Mr, 4 * nx, (4 + 1) * nx, 0, nx));  
    fit[4] = 10000 * fit[4] / 1e+10;  
    fit[5] = escaffer6_func(x, nx, Arrays.copyOfRange(Os, 5 * nx, (5 + 1) * nx), subMatrix(Mr, 5 * nx, (5 + 1) * nx, 0, nx));  
    fit[5] = 10000 * fit[5] / 2e+7; 
    double f = cf_cal(x, nx, Os, delta, bias, fit, cf_num);  
    return f;  
}  

    public double calculate_fit(double[] x ){
        int nx= x.length;
        double ff=0.0; 
        int cf_num = 10;
        //double[][]  M= new double[nx][nx]; 
        double[] OShift = new double[nx];
        double[] OShift1 = new double[(cf_num - 1) * nx];
        double[] SS = new double[nx]; 
        
        List<double[]> rows = new ArrayList<>();  
        double[][] M= null;

        // load matrix Rotate
        String fileNam = String.format("C:\\Users\\farzaneh\\Documents\\NetBeansProjects\\Hurestic\\src\\CEC2022\\input_data\\M_%d_D%d.txt", func_num, nx);  
        try (BufferedReader br = new BufferedReader(new FileReader(fileNam))) {  
            //List<double[]> rows = new ArrayList<>();  
            String line;  
            while ((line = br.readLine()) != null) {  
                String[] values = line.trim().split("\\s+");  
                double[] row = new double[values.length];  
                for (int i = 0; i < values.length; i++) {  
                    row[i] = Double.parseDouble(values[i]);  
                }  
                rows.add(row);  
            } 
            M= new double[rows.size()][nx];
            for (int i = 0; i < rows.size(); i++) {  
                    M[i] = rows.get(i);  
                }
            

        } catch (IOException e) {  
            System.err.println("\n Error: Cannot open " + fileNam + " for reading \n");  
             
        } 
        
        // بارگذاری OShift  
        String fileName = String.format("C:\\Users\\farzaneh\\Documents\\NetBeansProjects\\Hurestic\\src\\CEC2022\\input_data\\shift_data_%d.txt", func_num);  
         try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {  
            List<Double> oShiftTempList = new ArrayList<>();  
            String line; 
            while ((line = br.readLine()) != null) {  
                if (line.trim().isEmpty()) {  
                    continue;  // گذر از خط‌های خالی  
                } 
                // Split the line by spaces and parse each value  
                String[] values = line.trim().split("\\s+");  // Use \\s+ to handle multiple spaces  
                for (String value : values) {  
                    try {  
                        oShiftTempList.add(Double.valueOf(value));  
                    } catch (NumberFormatException e) {  
                        System.err.println("*** Error 1 parsing double: " + value + " - " + e.getMessage()); 
                    }  
                }  
            }  
            double[] OShift_temp = oShiftTempList.stream().mapToDouble(Double::doubleValue).toArray();  

            if (func_num < 9) {                  
                System.arraycopy(OShift_temp, 0, OShift, 0, nx);  
            } else {                
                for (int i = 0; i < cf_num - 1; i++) {  
                    for (int j = 0; j < nx; j++) {  
                        OShift1[i * nx + j] = OShift_temp[i * nx + j];  
                    }  
                }  
            }  
        } catch (IOException e) {  
            System.out.println("Error reading file: " + e.getMessage());  
        }  
          
        
        if (func_num >= 6 && func_num <= 8) {  
            fileName = String.format("C:\\Users\\farzaneh\\Documents\\NetBeansProjects\\Hurestic\\src\\CEC2022\\input_data\\shuffle_data_%d_D%d.txt", func_num, nx);  
            try (BufferedReader br = new BufferedReader(new FileReader(fileName))) {  

                List<Double> ssList = new ArrayList<>();  
                String line;  
                while ((line = br.readLine()) != null) {  
                     // Split the line by spaces and parse each value  
                    String[] values = line.trim().split("\\s+");  // Use \\s+ to handle multiple spaces  
                    for (String value : values) {  
                        try {  
                            ssList.add(Double.valueOf(value));  
                        } catch (NumberFormatException e) {  
                            System.err.println("Error 3 parsing double: " + value);  
                            // Handle the parsing error as needed  
                        }  
                    }  
                }  
                SS = ssList.stream().mapToDouble(Double::doubleValue).toArray();  
                } catch (IOException e) {  
            System.out.println("Error reading file: " + e.getMessage());  
        }  
    }  
     
 //--------------------------------------------------------------------
    if (M != null) {  // چک کنید که آرایه مقداردهی شده است       
         
            switch (func_num) {  
                case 1 -> {
                    ff = zakharov_func(x, nx, OShift, M)+ 300.0;
                }  
                case 2 -> {
                     ff = rosenbrock_func(x, nx, OShift, M)+ 400.0; 
                 }      
                case 3 -> { 
                      ff = schaffer_F7_func(x, nx, OShift, M)+ 600.0; 
                  }
                case 4 ->{
                      ff = rastrigin_func(x,nx, OShift, M)+ 800.0; 
                  }
                case 5 ->{
                      ff = levy_func(x, nx, OShift, M)+ 900.0; 
                  }  
                case 6 ->{
                      ff = Hybrid_Function_1(x, nx, OShift, M, SS)+ 1800.0;
                  } 
                case 7 ->{
                      ff = Hybrid_Function_2(x, nx, OShift, M, SS)+ 2000.0;
                  }  
                case 8 ->{
                      ff = Hybrid_Function_3(x, nx, OShift, M, SS)+ 2200.0; 
                  }  
                case 9 ->{
                      ff = Composition_Function_1(x, nx, OShift1, M)+ 2300.0;  
                  } 
                case 10 ->{
                      ff = Composition_Function_2(x, nx, OShift1, M)+ 2400.0;  
                  } 
                case 11 ->{
                      ff = Composition_Function_3(x, nx, OShift1, M)+ 2600.0; 
                  } 
                case 12 ->{
                      ff = Composition_Function_4(x, nx, OShift1, M)+ 2700.0;
                  }  
            }  
         
    }else{
         System.out.println("Array M is uninitialized and is Null.");  
    }
    return ff;
    }
}
