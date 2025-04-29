/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package CEC2022;


import java.util.Random;
import static java.util.concurrent.ThreadLocalRandom.current;

/**
 *
 * @author farzaneh
 */
public class TLBO {
    private final double lb, ub;
    private final int pop, dim, iter;
    private double[][] x_stud;
    private double[] cost_stud;
    private double[] mean;
    private double[] x_teacher;
    private double fit_teacher;
    Random rand = new Random();
    double[] best_fit_iter_TLB;
    int typefunc;
    
    public  TLBO(int population, double lb, double ub, int dimention, int maxiter, double[][] x_init, int typfunc){
        this.pop = population;
        this.iter = maxiter;
        this.dim = dimention;
        this.lb = lb;
        this.ub = ub;
        this.x_stud = x_init.clone();
        this.cost_stud = new double[pop];        
        this.x_teacher = new double[dim];
        this.fit_teacher= Double.MAX_VALUE;
        this.mean = new double[dim];
        this.typefunc = typfunc;
        this.best_fit_iter_TLB = new double[iter];
    }
    
  
    private double[] check_bound(double[] s){
        for(int j=0; j< dim ; j++){
                if( s[j]< lb ) {
                    s[j]= lb;
                } else
                    if(s[j]> ub){
                       s[j]= ub; 
                    }                       
        }
        return s;
    }
  
    private double[] calculateMeanSolution() {  
        double[] mean = new double[dim];  
        for (int j = 0; j < dim; j++) {  
            mean[j] = 0;  
        }  
        for (int i = 0; i < pop; i++) {  
            for (int j = 0; j < dim; j++) {  
                mean[j] += x_stud[i][j]; // جمع کردن مقادیر  
            }  
        }  
        for (int j = 0; j < dim; j++) {  
            mean[j] /= pop;  
        }  
        return mean;  
    }  
     

    
    public double[] execute_TLBO(){
        cec22fit fit= new cec22fit(typefunc);        
        
        for(int k=0; k<pop ; k++){
           cost_stud[k]= fit.calculate_fit(x_stud[k]);
           if( cost_stud[k] < fit_teacher){
             x_teacher = x_stud[k];
             fit_teacher= cost_stud[k];
            }
        }
        for( int t=0; t< iter; t++){
           
        //    Teacher phase
            for(int i=0; i<pop ; i++){
                
               mean= calculateMeanSolution();
                // Teaching Factor
                double Tf=  current().nextInt(1, 3);
                double[] new_sol= new double[dim];
                // Teaching (moving towards teacher)
                for(int j=0; j<dim ; j++){
                    new_sol[j]= x_stud[i][j]+ Math.random()*(x_teacher[j]-(Tf*mean[j]));
                }
                new_sol = check_bound(new_sol);
                double fit_new_sol= fit.calculate_fit(new_sol);
                if(fit_new_sol < cost_stud[i]){
                    x_stud[i] = new_sol.clone();
                    cost_stud[i]= fit_new_sol;
                }
            }
        //    Learner phase
            for(int i=0; i<pop ; i++){
                int rnd1 = rand.nextInt(pop);
                int rnd2 = rand.nextInt(pop);
                while(rnd1==rnd2){
                    rnd2=rand.nextInt(pop);
                }
                double[] xrand1 = x_stud[rnd1].clone();
                double[] xrand2 = x_stud[rnd2].clone();
                double fit_rnd1=fit.calculate_fit(xrand1);
                double fit_rnd2=fit.calculate_fit(xrand2);
                double[] new_sol= new double[dim];
                if(fit_rnd1 < fit_rnd2){
                     for(int j=0; j<dim ; j++){
                         double r= rand.nextDouble();
                         new_sol[j]= x_stud[i][j]+ (r*(xrand1[j]-xrand2[j]));
                     }                    
                }
                else{
                    for(int j=0; j<dim ; j++){
                        double r= rand.nextDouble();
                         new_sol[j]= x_stud[i][j]+ (r*(xrand2[j]-xrand1[j]));
                     } 
                }
                new_sol = check_bound(new_sol);
                double fit_new_sol= fit.calculate_fit(new_sol);
                if(fit_new_sol < cost_stud[i]){
                    x_stud[i] = new_sol.clone();
                    cost_stud[i]= fit_new_sol;
                }

            }
            for(int k=0; k<pop ; k++){                   
                   if( cost_stud[k] < fit_teacher){
                     x_teacher= x_stud[k];
                     fit_teacher= cost_stud[k];
                    }
                }
        best_fit_iter_TLB[t]=fit_teacher;    System.out.print(fit_teacher + "  ");
       }
       
       System.out.println("end _NNA");
     return best_fit_iter_TLB;
    }
    
    
}
