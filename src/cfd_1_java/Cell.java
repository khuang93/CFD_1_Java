/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cfd_1_java;

/**
 *
 * @author Kailin Huang
 */
public class Cell {
    int x;
    int index;
    double[] Ts;
    double[] Ts_star;
    double[] Tf;
    double[] Tf_star;
    
    Cell(){
        
    }
    Cell(int numCells){
        this.Ts=new double[numCells];
        this.Tf=new double[numCells];
        this.Ts_star=new double[numCells];
        this.Tf_star=new double[numCells];
    }
    
}
