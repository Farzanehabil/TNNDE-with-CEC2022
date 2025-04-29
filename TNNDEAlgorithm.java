/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package CEC2022;


import java.util.Random;
import java.util.concurrent.ThreadLocalRandom;
import static java.util.concurrent.ThreadLocalRandom.current;

/**
 *
 * @author farzaneh
 */
public class TNNDEAlgorithm {
    private int pop,iter,dim;
    private double lb=0, ub;
    private double[][] x;
    private double[] cost;        
    private double[][] worst;
    private double[][] best;
    private double[] cost_worst;
    private double[] cost_best;
    private double[][] weight;
    private double[] mean;
    private double beta = 1.0;
    Random rand = new Random();
    private double[] xglobal;
    private double globalfit;
    private double[] best_fit_iter;
    int typefunc;
        
        
    public TNNDEAlgorithm(int population, int dimention, double lb, double ub, int iter,double[][] x_init, int typfunc){
        this.typefunc = typfunc;
        this.pop = population;
        this.iter = iter;
        this.dim = dimention;
        this.ub = ub; 
        this.lb = lb;
        this.x = new double[pop][dim];
        this.cost = new double[pop];   
        this.best = new double[pop/2][dim];
        this.worst =new double[pop/2][dim];  
        this.cost_worst= new double[pop/2]; 
        this.cost_best= new double[pop/2];
        this.weight= new double[pop/2][pop/2];
        this.mean = new double[dim];
        this.xglobal= new double[dim];
        this.globalfit = Double.MAX_VALUE;
        this.best_fit_iter = new double[iter];
    }
    
    private void initPopulation() {          
	for (int i = 0; i < pop; i++) {
	    for (int j = 0; j < dim; j++) {
		x[i][j]= lb + (ub - lb)*rand.nextDouble();   
	    }
       }        
    }
    
    private void init_weight(){        
        for(int i=0; i<pop/2 ; i++){
            for(int j=0 ; j< pop/2 ; j++){
                weight[i][j]=  ThreadLocalRandom.current().nextDouble(0.0, 1.0);  
            }
        }
             weight= sumation_weight(weight); 

    }    
    private double[][] sumation_weight(double[][] w){
         // constraint of Summation each column = 1
         double Npop = pop/2;
        for(int j=0; j< Npop; j++){
            double sum=0;
            for(int i=0; i< Npop; i++){
                sum += w[i][j];
            }
            if (sum > 0) { 
                for (int i = 0; i < Npop; i++) {  
                    w[i][j] /= sum;  
                }  
            } else { 
                for (int i = 0; i < Npop; i++) {  
                    w[i][j] = Math.random()+0.001;  
                }  
                sum = 0;  
                for (int i = 0; i < Npop; i++) {  
                    sum += w[i][j];  
                }  
                for (int i = 0; i < Npop; i++) {  
                    w[i][j] /= sum;  
                }     
            }    
        }
        return w;
    
    }
    

    private void sort_divid_array(){
        double temp;
        double[] tem;
        for(int i=0; i<pop; i++){
            for(int j=0; j<pop-i-1;j++){
                if(cost[j] < cost[j+1]){
                      tem = x[j].clone();
                      x[j] = x[j+1].clone();
                      x[j+1]= tem.clone();
                      temp = cost[j];
                      cost[j]=cost[j+1];
                      cost[j+1]= temp;                                    
                }
            }
        }
        //Divid Array best And Worst            
        for(int i=0; i<pop/2; i++){
            worst[i]= x[i].clone();  cost_worst[i]=cost[i];    
            best[i]= x[i+pop/2].clone();  cost_best[i]=cost[i+pop/2];
       }            
    } 
  
    
    private int min_index(double[] fit ){
         int index=0;
        double minfit= Double.MAX_VALUE;
        for(int i=0; i< pop/2; i++){
            if(fit[i]< minfit){
                index=i;
                minfit= fit[i];
            }
        }
        return index;
    }
    private int maxcost(){
         int index=0;
        double maxfit= Double.MIN_VALUE;
        for(int i=0; i< pop/2; i++){
            if(cost_worst[i]> maxfit){
                index=i;
                maxfit= cost_worst[i];
            }
        }
        return index;
    }
    
 
    private double[] check_bound(double[] s){
        for(int j=0; j< dim ; j++){
            if( s[j]< lb )   s[j]= lb;
            if(s[j]> ub)     s[j]= ub;
        }
        return s;
    }
    
    
    
    
    private double[] crossover(double[] targetVector, double[] referenceVector) {  
        double cr = 0.7; 
        double[] child = new double[dim]; 
        int s = rand.nextInt(dim); 
        for (int i = 0; i < dim; i++) {  
            double r = rand.nextDouble(); 
            if (r <= cr || i == s) {  
                child[i] = targetVector[i]; 
            } else {  
                child[i] = referenceVector[i]; 
            }  
        }  
        return child; 
    }  
  
    
    private double[] calculate_mean(){
        double[] average = new double[dim];
        for(int j=0; j<dim; j++){
            double sum=0;
            for(int i=0; i<pop/2; i++){
                sum+=best[i][j];
            }
            average[j]= sum/(pop/2);
        }
        return average;
    }
    
    
    
public double[] execute_TLNND(){ 
        initPopulation();
        init_weight();
        cec22fit fit= new cec22fit(typefunc);
         for(int i=0; i<pop ; i++){
            cost[i]= fit.calculate_fit(x[i]);
            if(cost[i] < globalfit){
                xglobal = x[i].clone();
                globalfit= cost[i];                
            }
        }
        int index;
        int Npop = pop/2 ;
        double[] x_teacher;
        double target;
        double[] x_target;
        double[] w_target= new double[Npop];
    // iteration program    
   for(int t=0; t< iter; t++){                  
        sort_divid_array();
   //process population with NNA for worst half ----------------------------------------
            target = globalfit;
            x_target = xglobal.clone();
            if(t==0){
                index= rand.nextInt(Npop);
                for(int k=0; k<Npop ; k++){
                    w_target[k]= weight[k][index];
                }
            }
            double[][] x_new = new double[Npop][dim]; 
            for (int i = 0; i < Npop; i++) {  
                for (int j = 0; j < dim; j++) { 
                     x_new[i][j] = 0;
                    for (int k = 0; k < Npop; k++) {  
                         x_new[i][j] += weight[k][i] * worst[k][j];  //???????????
                    }  
                } check_bound(x_new[i]);
            } 
            for(int i=0; i<Npop; i++){
                for(int j=0; j<dim; j++){
                    worst[i][j]=( x_new[i][j] + worst[i][j])* rand.nextDouble();
                }check_bound(worst[i]);
            }
            double[][] test= new double[Npop][Npop];
            for(int i=0; i< Npop; i++){   
               for(int j=0; j<Npop; j++){ 
                  double value =  weight[j][i]+ ( (w_target[j]- weight[j][i])* 2 * rand.nextDouble()); 
                    if (value < 0){
                        test[j][i]= Math.abs(value);
                    }else 
                        test[j][i]= value;
                 }
            }
        weight = test.clone();
        weight= sumation_weight(weight);
        for(int i=0; i<Npop; i++){                
            if(Math.random() < beta){
                int Nb = (int) Math.round(dim*beta);
                for(int j=0; j<Nb; j++ ){
                    worst[i][rand.nextInt(dim)]= lb + (ub - lb)* rand.nextDouble();
                }
                int Nw = (int) Math.round((Npop)*beta);
                for(int j=0; j< Nw; j++){
                    weight[j][rand.nextInt(Npop)]= rand.nextDouble()+ 0.001;
                }
            }else{                    
                 for(int j=0; j<dim; j++){
                        worst[i][j]= worst[i][j]+(2*Math.random()*(x_target[j]- worst[i][j]));
                    }                    
                }
                check_bound(worst[i]);
            }
            weight= sumation_weight(weight);
            beta= 0.99 * beta;
            if(beta<0.01) {
            beta = rand.nextDouble()* 0.05;
            }
            for(int i=0; i< Npop; i++){
                cost_worst[i]= fit.calculate_fit(worst[i]);
            }
            index = min_index(cost_worst);
            
            if(cost_worst[index] < target){
                target = cost_worst[index];
                x_target = worst[index].clone();                         
                for(int k=0; k<Npop ; k++){
                    w_target[k]= weight[k][index];
                }             
            }else{
                 index = maxcost();
                 worst[index] = x_target.clone();
                 cost_worst[index] = target;
                 weight[index] = w_target;
             }
            
           if(target < globalfit){
                xglobal = x_target.clone();
                globalfit = target;                
            } 
    //process population with tlbo for half best
        x_teacher=xglobal;
        double fit_teacher= globalfit;
        //    Teacher phase
        for(int i=0; i<pop/2; i++){
            mean= calculate_mean();
            // Teaching Factor
            double Tf=  current().nextInt(1, 3);
            double[] new_sol= new double[dim];
            double fit_new_sol;
            // Teaching (moving towards teacher)
            for(int j=0; j<dim ; j++){
                new_sol[j]= best[i][j]+ Math.random()*(x_teacher[j]-(Tf*mean[j]));
            }
            new_sol = check_bound(new_sol);
            fit_new_sol= fit.calculate_fit(best[i]);
            if(fit_new_sol < cost_best[i]){
                best[i] = new_sol.clone();
                cost_best[i]= fit_new_sol;
            }           
        }
         //    Learner phase
            for(int i=0; i<pop/2 ; i++){
                int rnd1 = rand.nextInt(pop/2);
                int rnd2 = rand.nextInt(pop/2);
                while(rnd1==rnd2){
                    rnd2=rand.nextInt(pop/2);
                }
                double[] xrand1 = best[rnd1].clone();
                double[] xrand2 = best[rnd2].clone();
                double fit_rnd1 = fit.calculate_fit(xrand1);
                double fit_rnd2 = fit.calculate_fit(xrand2);
                double[] new_sol= new double[dim];
                if(fit_rnd1 < fit_rnd2){
                     for(int j=0; j<dim ; j++){
                         double r= rand.nextDouble();
                         new_sol[j]= best[i][j]+ (r*(xrand1[j]-xrand2[j]));
                     }                    
                }
                else{
                    for(int j=0; j<dim ; j++){
                        double r= rand.nextDouble();
                         new_sol[j]= best[i][j]+ (r*(xrand2[j]-xrand1[j]));
                     } 
                }
                new_sol = check_bound(new_sol);
                double fit_new_sol = fit.calculate_fit(new_sol);
                if(fit_new_sol < cost_best[i]){
                    best[i] = new_sol.clone();
                    cost_best[i]= fit_new_sol;
                }
            }
            for(int i=0; i<pop/2; i++){
               if(cost_best[i] < fit_teacher){
                    x_teacher= best[i].clone();
                    fit_teacher= cost_best[i];
                }
            }
            if(fit_teacher < globalfit ){
                globalfit = fit_teacher;
                xglobal = x_teacher.clone();
            }
     // End TLBO           
           
    //===========================================================================
        for(int i=0; i< Npop; i++){
            if(cost_worst[i]< cost[i]){
                x[i]= worst[i].clone();
                cost[i]= cost_worst[i];
            }else {
                worst[i]= crossover(worst[i], xglobal);     check_bound(worst[i]);
                cost_worst[i]= fit.calculate_fit(worst[i]);
                if(cost_worst[i]< cost[i]){
                    x[i]= worst[i];
                    cost[i]= cost_worst[i];
                }
            }    
            if(cost_best[i]< cost[i+Npop]){
                x[i+Npop] = best[i];
                cost[i+Npop]= cost_best[i];
            }else{
                best[i]= crossover(best[i], xglobal);     check_bound(best[i]);
                cost_best[i]= fit.calculate_fit(best[i]);
                if(cost_best[i]< cost[i]){
                    x[i]= best[i];
                    cost[i]= cost_best[i];
                }
            }
        }
     best_fit_iter[t] = globalfit;           System.out.print(best_fit_iter[t] + "  ");
      }
        System.out.println();
    return best_fit_iter ;
   }
}
