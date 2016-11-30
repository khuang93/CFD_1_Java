/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cfd_1_java;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

/**
 *
 * @author Kailin Huang
 */
public class CFD_1_Java {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        // TODO code application logic here
        double height, diameter, initTemp;
        final int numberCells;
        Scanner sc = new Scanner(System.in);
        height = sc.nextDouble();
        diameter = sc.nextDouble();
        numberCells = sc.nextInt();
        initTemp = sc.nextDouble();

        double deltaX = height / numberCells;

        File file1 = new File("plot_data_pr1_JAVA.csv");

        double t_charging = 2500;
        double t_idle_C_Dis = 2000;
        double t_discharging = 500;
        double t_idle_Dis_C = 1000;
        ArrayList<Double> t_states = new ArrayList<>();
        t_states.add(t_charging);
        t_states.add(t_idle_C_Dis);
        t_states.add(t_discharging);
        t_states.add(t_idle_Dis_C);
        int numberCycles = 1;
        int timeStepsPerCycle = 1;
        for (double t : t_states) {
            timeStepsPerCycle += t;
        }
        timeStepsPerCycle -= 1;
        int total_t_per_cycle = timeStepsPerCycle * numberCycles;

        double delta_t = total_t_per_cycle / timeStepsPerCycle;

        double alpha_f = 0.1;
        double alpha_s = 0.1;
        double uf = 0;
        int total_timesteps = timeStepsPerCycle * numberCycles;

        class CellNode {

            Cell thisCell;
            CellNode next;
            CellNode prev;
            
            CellNode(int numCells){
                this.thisCell=new Cell(numCells);
            }

        }

        CellNode cellT0 = new CellNode(numberCells);
        CellNode cellCurrentT;
        CellNode cellLastT=new CellNode(numberCells);
        cellT0.prev = null;
        cellCurrentT = cellT0;
        cellCurrentT.thisCell.Tf = new double[numberCells];
            cellCurrentT.thisCell.Ts = new double[numberCells];
        for (int i = 0; i < total_timesteps; i++) {
            cellCurrentT.thisCell.Tf = new double[numberCells];
            cellCurrentT.thisCell.Ts = new double[numberCells];
            for (int xi = 0; xi < numberCells; xi++) {
                cellCurrentT.thisCell.Tf[xi] = initTemp;
                cellCurrentT.thisCell.Ts[xi] = initTemp;
                System.out.print("i"+i+"xi"+xi+"Tf"+cellCurrentT.thisCell.Tf[xi]+"\n");
                
            }
            cellCurrentT.next=new CellNode(numberCells);
            cellCurrentT.next.prev=cellCurrentT;
            cellCurrentT=cellCurrentT.next;
            
            if(i==total_timesteps-1){
                cellLastT=cellCurrentT;
            }

        }
        cellLastT.next=null;
        cellCurrentT=cellT0.next;
        int ii = 0;
        while(cellCurrentT.next!=null &&cellCurrentT.prev!=null){
            System.out.println(ii+"numberCells"+numberCells);
            for(int n = 1; n<numberCells-1;n++){
                cellCurrentT.thisCell.Tf[n+1]= cellCurrentT.thisCell.Tf[n];
                
                System.out.print("n"+n+"Tf"+cellCurrentT.thisCell.Tf[n]+"\n");
                
                
            }
            cellCurrentT=cellCurrentT.next;
            ii++;
        }
        
        

    }

}
