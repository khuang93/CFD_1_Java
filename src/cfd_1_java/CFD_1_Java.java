/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package cfd_1_java;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Scanner;
import java.math.*;
import java.util.logging.Level;
import java.util.logging.Logger;

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
        double alpha_f = 2E-7;
        double alpha_s = 9E-7;
        double uf = 0.1;
        double Tf_in = 773;
        double hv=1000;
        double epsilon = 0.4;
        double rho_f=1835.6;
        double rho_s = 2600;
        double Cp_f=1511.8;
        double Cs=900;
        double ks = 2;
        double kf=0.52;
        double m_f_dot = 0.1;
        double mu_f=2.63;
        
        double hv_f = hv/(epsilon*rho_f*Cp_f);
        double hv_s = hv/((1-epsilon)*rho_s*Cs);

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
        int numberCycles = 10;
        int timeStepsPerCycle = 1;
        for (double t : t_states) {
            timeStepsPerCycle += t;
        }
        timeStepsPerCycle -= 1;
        int total_t_per_cycle = timeStepsPerCycle * numberCycles;

        double delta_t = total_t_per_cycle / timeStepsPerCycle;

        int total_timesteps = timeStepsPerCycle * numberCycles;

        class CellNode {

            Cell thisCell;
            CellNode next;
            CellNode prev;

            CellNode(int numCells) {
                this.thisCell = new Cell(numCells);
            }

        }

        int initNumberCells = numberCells + 1;

        CellNode cellT0 = new CellNode(initNumberCells);
        CellNode cellCurrentT;
        CellNode cellLastT = new CellNode(initNumberCells);
        cellT0.prev = null;
        cellCurrentT = cellT0;
        cellCurrentT.thisCell.Tf = new double[initNumberCells];
        cellCurrentT.thisCell.Ts = new double[initNumberCells];

        int num = 1;
        double k = 2 * Math.PI * num / height;

        for (int i = 0; i < total_timesteps + 1; i++) {
            cellCurrentT.thisCell.Tf = new double[initNumberCells];
            cellCurrentT.thisCell.Ts = new double[initNumberCells];
            cellCurrentT.thisCell.Tf[0] = Tf_in;
            cellCurrentT.thisCell.Ts[0] = initTemp;
            for (int xi = 1; xi < numberCells; xi++) {
                cellCurrentT.thisCell.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
                cellCurrentT.thisCell.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
                cellCurrentT.thisCell.Tf_star[xi] = 0;//Math.cos(k * xi * deltaX);
                cellCurrentT.thisCell.Ts_star[xi] = 0;//Math.sin(k * xi * deltaX);
                //  System.out.print("i" + i + "xi" + xi + "Tf" + cellCurrentT.thisCell.Tf[xi] + "\n");

            }
            cellCurrentT.next = new CellNode(initNumberCells);
            cellCurrentT.next.prev = cellCurrentT;
            cellCurrentT = cellCurrentT.next;

            if (i == total_timesteps - 1) {
                cellLastT = cellCurrentT;
            }

        }

        cellLastT.next = null;
        cellCurrentT = cellT0.next;

        double[][] matrixM = new double[2][2];
        matrixM[0][0] = 1 + hv_f * delta_t;
        matrixM[1][0] = 0 - hv_s * delta_t;
        matrixM[0][1] = -hv_f * delta_t;
        matrixM[1][1] = 1 + hv_s * delta_t;

        double OneDivDetM = 1 / (matrixM[0][0] * matrixM[1][1] - matrixM[0][1] * matrixM[1][0]);

        double[][] matrixM_inv = new double[2][2];
        matrixM_inv[0][0] = OneDivDetM * matrixM[1][1];
        matrixM_inv[1][0] = 0 - OneDivDetM * matrixM[1][0];
        matrixM_inv[0][1] = 0 - OneDivDetM * matrixM[0][1];
        matrixM_inv[1][1] = OneDivDetM * matrixM[0][0];

        int ii = 0;

        while (cellCurrentT.next != null && cellCurrentT.prev != null) {
            if (ii % 500 == 0) {
                try {
                    PrintStream fout = new PrintStream(new File("plot_data_pr3_t" + ii + "java.csv"));
                    fout.println("x, Ts, Tf");
                    for (int n = 0; n < numberCells; n++) {
                        fout.println(n * deltaX + "," + cellCurrentT.prev.thisCell.Ts[n] + "," + cellCurrentT.prev.thisCell.Tf[n]);
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

            for (int n = 0; n < numberCells; n++) {
                //backup
                /* cellCurrentT.thisCell.Tf[n + 1] = cellCurrentT.thisCell.Tf[n]
                        - uf * delta_t / deltaX * (cellCurrentT.thisCell.Tf[n] - cellCurrentT.prev.thisCell.Tf[n])
                        + alpha_f * delta_t / (deltaX * deltaX)
                        * (cellCurrentT.next.thisCell.Tf[n] - 2 * cellCurrentT.thisCell.Tf[n] + cellCurrentT.prev.thisCell.Tf[n]);

                //System.out.print("n" + n + "Tf" + cellCurrentT.thisCell.Tf[n] + "\n");
                cellCurrentT.thisCell.Ts[n + 1] = cellCurrentT.thisCell.Ts[n] + alpha_s * delta_t / (deltaX * deltaX)
                        * (cellCurrentT.next.thisCell.Ts[n] - 2 * cellCurrentT.thisCell.Ts[n] + cellCurrentT.prev.thisCell.Ts[n]);
                 */
                cellCurrentT.thisCell.Tf_star[n] = cellCurrentT.thisCell.Tf[n] - uf * delta_t / deltaX * (cellCurrentT.thisCell.Tf[n] - cellCurrentT.prev.thisCell.Tf[n])
                        + alpha_f * delta_t / (deltaX * deltaX)
                        * (cellCurrentT.next.thisCell.Tf[n] - 2 * cellCurrentT.thisCell.Tf[n] + cellCurrentT.prev.thisCell.Tf[n]);

                //System.out.print("n" + n + "Tf" + cellCurrentT.thisCell.Tf[n] + "\n");
                cellCurrentT.thisCell.Ts_star[n] = cellCurrentT.thisCell.Ts[n] + alpha_s * delta_t / (deltaX * deltaX)
                        * (cellCurrentT.next.thisCell.Ts[n] - 2 * cellCurrentT.thisCell.Ts[n] + cellCurrentT.prev.thisCell.Ts[n]);

                //now use inv matrix to get Tfn+1 and tsn+1
                cellCurrentT.thisCell.Tf[n + 1] = matrixM_inv[0][0] * cellCurrentT.thisCell.Tf_star[n] + matrixM_inv[0][1] * cellCurrentT.thisCell.Ts_star[n];

                cellCurrentT.thisCell.Ts[n + 1] = matrixM_inv[1][0] * cellCurrentT.thisCell.Tf_star[n] + matrixM_inv[1][1] * cellCurrentT.thisCell.Ts_star[n];

            }

            cellCurrentT = cellCurrentT.next;
            ii++;
        }
        try {
            PrintStream fout = new PrintStream(new File("plot_data_pr3_final_java.csv"));
            fout.println("x, Ts, Tf");
            for (int n = 0; n < numberCells; n++) {
                fout.println(n * deltaX + "," + cellCurrentT.prev.thisCell.Ts[n] + "," + cellCurrentT.prev.thisCell.Tf[n]);
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

}
