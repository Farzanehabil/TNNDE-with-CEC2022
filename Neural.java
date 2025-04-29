/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package CEC2022;




import java.util.Random;


/**
 *
 * @author farzaneh
 */
public class Neural {
        private final int pop,iter,dim;
        private final double lb, ub;
        private double[][] x_patern;
        private double[] cost;        
        private double[][] weight;
        private double target= Double.MAX_VALUE;
        private double[] x_target;   
        private double[] w_target;        
        private double[][] x_new; 
        double beta = 1.0;
        Random rand = new Random();
        private double[] best_fit_iter_NNA;
        int typefunc;
        
    public Neural(int population, int dimention, double lb, double ub, int iter,double[][] xinit, int typfunc){
    this.typefunc = typfunc;
    this.best_fit_iter_NNA = new double[iter];
    this.pop = population;
    this.iter = iter;
    this.dim = dimention;
    this.lb = lb;
    this.ub = ub;
    this. x_patern = xinit.clone();
    this.cost = new double[pop];        
    this.weight= new double[pop][pop];
    this.x_target= new double [dim];       
    this.w_target= new double[pop];
    this.x_new= new double[pop][dim];


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
     private void init_weight(){        
        for(int i=0; i<pop ; i++){
            for(int j=0 ; j< pop ; j++){
                weight[i][j]= rand.nextDouble()*0.5;
            }
        }
        weight=sumation_weight(weight);       
    }    
    private double[][] sumation_weight(double[][] w){
        for(int j=0; j< pop; j++){
            double sum=0;
            for(int i=0; i< pop; i++){
                sum += w[i][j];
            }
            for(int i=0; i< pop; i++){
                w[i][j]= w[i][j]/sum;
            }
        }
        return w;
    } 
    private int min_index(){
         int index=0;
        double minfit= Double.MAX_VALUE;
        for(int i=0; i< pop; i++){
            if(cost[i]< minfit){
                index=i;
                minfit= cost[i];
            }
        }
        return index;
    }
    private int maxcost(){
         int index=0;
        double maxfit= Double.MIN_VALUE;
        for(int i=0; i< pop; i++){
            if(cost[i]> maxfit){
                index=i;
                maxfit= cost[i];
            }
        }
        return index;
    }
    
    private void calculate_x_new(){
        for(int i=0; i<pop ;i++){
                    for(int j=0 ; j< dim ; j++){//                      
                         double sum=0;
                        for(int k=0; k<pop; k++){
                            sum+= (weight[k][i]*  x_patern[k][j]) ;                           
                        } 
                         x_new[i][j] = sum;  
                    }
                check_bound(x_new[i]);    
            }
    }
    
    private void update_x_pattern(){
        
        for(int i=0; i<pop; i++){
            for(int j=0; j<dim; j++){
                     x_patern[i][j]= x_new[i][j] + x_patern[i][j];
            }
                  check_bound(x_patern[i]);
        }
        
    }
    private void udate_weight_matrix(){
        for(int j=0; j< pop; j++){
                for(int i=0; i<pop; i++){ 
                    double r1= rand.nextDouble();
                    weight[i][j]=Math.abs(weight[i][j]+(2*  r1 *(w_target[i]- weight[i][j])));  
                }
            }   
            weight= sumation_weight(weight);
    }
    
    
    
    private double bias_reduction(double beta){
       
        beta= 0.99 * beta;
        if(beta<0.01) {
            beta=0.05;
        }
       return beta;
    }
   

    public double[] execute_neural(){ 
        cec22fit fit= new cec22fit(typefunc);
        
         init_weight();
          
        for(int i=0; i<pop ; i++){
            cost[i]= fit.calculate_fit(x_patern[i]);
        }
         
        int index = min_index();
        target = cost[index];
        x_target =  x_patern[index].clone();
        
        for(int k=0; k<pop ; k++){
            w_target[k]= weight[k][index];
        }
          
   for(int t=0; t< iter; t++){     
          calculate_x_new();  
          update_x_pattern(); 
          udate_weight_matrix();
           
        //---------------------------------------------------
        for(int i=0; i<pop; i++){
                double r= Math.random();
                if(r < beta){
                    int Nb = (int) Math.round(dim*beta);
                    for(int j=0; j<Nb; j++ ){
                         x_patern[i][rand.nextInt(dim)]= lb + (ub - lb)* Math.random();
                    }
                    
                    int Nw = (int) Math.round(pop*beta);
                    for(int j=0; j< Nw; j++){
                        weight[j][rand.nextInt(pop)]= rand.nextDouble();
                    }
                    weight= sumation_weight(weight);
                }                       
                else{                    
                    for(int j=0; j<dim; j++){
                         x_patern[i][j]=  x_patern[i][j]+(2*Math.random()*(x_target[j]-  x_patern[i][j]));
                    }                    
                }
                check_bound( x_patern[i]);
        }
    //---------------------------------------------------------
             beta = bias_reduction(beta);
             for(int i=0; i< pop; i++){
                cost[i]= fit.calculate_fit( x_patern[i]);
            }
             
             index = min_index();
             if(cost[index] < target){
                  target = cost[index];
                  x_target =  x_patern[index].clone();                         
                  for(int k=0; k<pop ; k++){
                        w_target[k]= weight[k][index];
                   }             
            }
             
             else{
                int index2= maxcost();
                 x_patern[index2]= x_target.clone();
                cost[index2]=target;                
                for(int k=0; k<pop ; k++){
                    weight[k][index2]= w_target[k];
                }
            }
            best_fit_iter_NNA[t]= target;    System.out.print(target + " ");
        }
         
      System.out.println("end _NNA");
         
       return best_fit_iter_NNA;
    }
    

}
