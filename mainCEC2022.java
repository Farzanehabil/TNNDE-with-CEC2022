/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Main.java to edit this template
 */
package CEC2022;


import java.awt.Color;

import org.jfree.data.xy.XYSeries;
import java.awt.BasicStroke;
import java.awt.Font;
import org.jfree.data.xy.XYSeriesCollection;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.chart.axis.ValueAxis;

import org.jfree.chart.ChartFactory;
import org.jfree.chart.ChartFrame;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.plot.PlotOrientation;
import org.jfree.chart.plot.XYPlot;





/**
 *
 * @author farzaneh
 */
public class mainCEC2022 {
    
    /**
     * @param args the command line arguments
     */
        
    public static void main(String[] args) {
        
        int numberOfRuns = 1;
        int type = 9; 
        
        double lb_plot = 2000;
        double ub_plot = 4000;
        
        final int maxiter = 1000;
        final int pop = 100;
        final int dim = 10;                  // Dim: 10 and 20
        final double lb = -100;
        final double ub = 100;
        double[][] Results_NNA   = new double[numberOfRuns][maxiter];  
        double[][] Results_TLBO  = new double[numberOfRuns][maxiter];
        double[][] Results_TNNDE = new double[numberOfRuns][maxiter];
        
        
        
        double[][] x= new double[pop][dim];
        for (int i = 0; i < pop; i++) {
	    for (int j = 0; j < dim; j++) {
		x[i][j] = lb + (ub - lb) * Math.random();
	    }
	}
   
        for (int i = 0; i < numberOfRuns; i++) {  
            System.out.println("Run CEC2022 : " + (i + 1));
            Neural nna= new Neural(pop, dim ,lb, ub,  maxiter, x , type);
            double[] neural= nna.execute_neural();
            Results_NNA[i] = neural;
            TLBO tlbo= new TLBO(pop ,lb, ub, dim, maxiter, x, type);
            double[] teach = tlbo.execute_TLBO();
            Results_TLBO[i] = teach;
            TNNDEAlgorithm tnnde= new  TNNDEAlgorithm(pop ,dim, lb, ub, maxiter, x , type);
            double[] Tnnde = tnnde.execute_TNNDE();
            Results_TNNDE[i]= Tnnde;
    //--------------------------------
        }
        double[] MeanEvolution_NNA = calculateMeanEvolution(Results_NNA, maxiter,numberOfRuns);  
        double[] MeanEvolution_TLBO = calculateMeanEvolution(Results_TLBO, maxiter,numberOfRuns);
        double[] MeanEvolution_TNNDE = calculateMeanEvolution(Results_TNNDE, maxiter,numberOfRuns);
        
        double mean_NNA = calculateMean(MeanEvolution_NNA);
        System.out.println("Mean Objective Value Neural: " + mean_NNA);  
        double SD_NNA = calculateStandardDeviation(MeanEvolution_NNA, mean_NNA);  
        System.out.println("Standard Deviation Neural: " + SD_NNA); 
    
        double mean_TLBO = calculateMean(MeanEvolution_TLBO);
        System.out.println("Mean Objective Value TLBO: " + mean_TLBO);  
        double SD_TLBO = calculateStandardDeviation(MeanEvolution_TLBO, mean_TLBO);  
        System.out.println("Standard Deviation TLBO: " + SD_TLBO); 
        
        double mean_TNNDE = calculateMean(MeanEvolution_TNNDE);
        System.out.println("Mean Objective Value TNNDE: " + mean_TNNDE);  
        double SD_TNNDE = calculateStandardDeviation(MeanEvolution_TNNDE, mean_TNNDE);  
        System.out.println("Standard Deviation TNNDE: " + SD_TNNDE);

 
        XYSeriesCollection  dataset=new XYSeriesCollection();
         XYSeries series0= new XYSeries("NNA");
         XYSeries series1= new XYSeries("TLBO");         
         XYSeries series2= new XYSeries("TNNDE");
        
       for(int k=0; k< MeanEvolution_NNA.length; k++){
             series0.add(k+1,MeanEvolution_NNA[k]);
             series1.add(k+1,MeanEvolution_TLBO[k]);
             series2.add(k+1,MeanEvolution_TNNDE[k]);
        }       
        dataset.addSeries(series0);
        dataset.addSeries(series1);
        dataset.addSeries(series2);
        
     JFreeChart chart= ChartFactory.createXYLineChart("F" + type , "iteration", "Objective Function Value" ,
                dataset,PlotOrientation.VERTICAL, true,true, false );
        
        XYPlot plot= (XYPlot) chart.getPlot();
        
        ValueAxis yAxis = plot.getRangeAxis();
        
        
        yAxis.setRange(lb_plot, ub_plot); 
        
        Font font= new Font("Adobe Arabic" , Font.BOLD,20);
        Font num_font= new Font("Adobe Arabic" , Font.BOLD, 20);
        plot.getDomainAxis().setLabelFont(font);
        plot.getRangeAxis().setLabelFont(font);
        plot.getDomainAxis().setTickLabelFont(num_font);
        plot.getRangeAxis().setTickLabelFont(num_font);
        LegendTitle legend= chart.getLegend();
        legend.setItemFont(font);
        
        XYLineAndShapeRenderer render= (XYLineAndShapeRenderer) plot.getRenderer();
        plot.setBackgroundPaint(Color.WHITE);
        render.setSeriesPaint(0, Color.GREEN);    
        render.setSeriesPaint(1, Color.red);         
        render.setSeriesPaint(2, Color.BLUE);
        
        for(int i=0; i<3; i++){
            plot.getRenderer().setSeriesStroke(i, new BasicStroke(2.0f));
        }
        ChartFrame frame= new ChartFrame("Line plot" , chart);
        frame.setSize(800,900);
        frame.setVisible(true);
        
    }
    public static double[] calculateMeanEvolution(double[][] data, int maxIterations, int num_run) {  
        double[] meanEvolution = new double[maxIterations];  
        for (int i = 0; i < maxIterations; i++) {  
            double sum = 0; 
            for (int j=0; j< num_run; j++) {
                sum += data[j][i];
            }  
            meanEvolution[i] = sum / num_run;  
        }  
        return meanEvolution;  
    }  
    public static double calculateMean(double[] data) {  
        double sum = 0;  
        for (int i=0; i<data.length ; i++) {  
            sum += data[i];  
        }  
        return sum / data.length;  
    } 
    public static double calculateStandardDeviation(double[] data, double mean) {  
        double sumOfSquaredDifferences = 0;  
        for (double value : data) {  
            sumOfSquaredDifferences += Math.pow(value - mean, 2);  
        }  
        double variance = sumOfSquaredDifferences / data.length;  
        return Math.sqrt(variance);  
    }  
    
}
