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

    static final int CHARGING = 0;
    static final int IDLE_C_DIS = 1;
    static final int DISCHARGING = 2;
    static final int IDLE_DIS_C = 3;

    static double delta_t;

    static double Re;
    static double Pr;
    static double Nu_fs;
    static double h_fs;
    static double h;

    static double hv;

    static double hv_f;
    static double hv_s;

    static double alpha_f = 2E-7;
    static double alpha_s = 9E-7;
    static double uf = 0;
    static double uf_charging = 0.1;
    static double uf_discharging = -0.1;
    static double Tf_in = 873;

    static double epsilon = 0.4;
    static double rho_f = 1835.6;
    static double rho_s = 2600;
    static double Cp_f = 1511.8;
    static double Cs = 900;
    static double ks = 2;
    static double kf = 0.52;
    static double m_f_dot = 0.1;
    static double mu_f = 2.63;
    static double ds = 0.03;

    static double[][] matrixM = new double[2][2];
    static double[][] matrixM_inv = new double[2][2];

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

        updateParameters();

        File file1 = new File("plot_data_pr1_JAVA.csv");

        Double t_charging = 2500.;
        Double t_idle_C_Dis = 2000.;
        Double t_discharging = 500.;
        Double t_idle_Dis_C = 1000.;
        ArrayList<Double> t_states = new ArrayList<>();
        t_states.add(t_charging);
        t_states.add(t_idle_C_Dis);
        t_states.add(t_discharging);
        t_states.add(t_idle_Dis_C);

        int numberCycles = 1;

        int timeStepsPerCycle = 1;
        int total_t_per_cycle = 1;
        for (double t : t_states) {
            total_t_per_cycle += t;
            timeStepsPerCycle = (int) (timeStepsPerCycle + (t * 10));
        }
        timeStepsPerCycle -= 1;
        total_t_per_cycle -= 1;

        delta_t = total_t_per_cycle * 1.0 / timeStepsPerCycle;

        int total_timesteps = timeStepsPerCycle * numberCycles;

        int initNumberCells = numberCells + 2;

        //creation of the linked list. each node represents a timestep.
        TimeStepNode T0 = new TimeStepNode(initNumberCells);
        TimeStepNode currentTimeStepNode;
        TimeStepNode TLast = new TimeStepNode(initNumberCells);
        T0.prev = null;
        currentTimeStepNode = T0;
        currentTimeStepNode.thisTS.Tf = new double[initNumberCells];
        currentTimeStepNode.thisTS.Ts = new double[initNumberCells];

        int num = 1;
        double k = 2 * Math.PI * num / height;

        //initialization of Tf and Ts
        for (int ti = 0; ti < total_timesteps + 1; ti++) {
            currentTimeStepNode.thisTS.Tf = new double[initNumberCells];
            currentTimeStepNode.thisTS.Ts = new double[initNumberCells];
            currentTimeStepNode.thisTS.Tf[0] = Tf_in;
            currentTimeStepNode.thisTS.Ts[0] = initTemp;
            for (int xi = 1; xi < initNumberCells; xi++) {
                currentTimeStepNode.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
                currentTimeStepNode.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
//                currentTimeStepNode.thisTS.Tf_star[xi] = 0;//Math.cos(k * xi * deltaX);
//                currentTimeStepNode.thisTS.Ts_star[xi] = 0;//Math.sin(k * xi * deltaX);
                //  System.out.print("i" + i + "xi" + xi + "Tf" + cellCurrentTime.thisTS.Tf[xi] + "\n");

            }
            currentTimeStepNode.next = new TimeStepNode(initNumberCells);
            currentTimeStepNode.next.prev = currentTimeStepNode;
            currentTimeStepNode = currentTimeStepNode.next;

            if (ti == total_timesteps - 1) {
                TLast = currentTimeStepNode;
            }

        }

        TLast.next = null;
        currentTimeStepNode = T0.next;

//        double[][] matrixM = new double[2][2];
//        double[][] matrixM_inv = new double[2][2];
        int currentTimeStep = 0;

        while (currentTimeStepNode.next != null && currentTimeStepNode.prev != null) {
            int status = getCurrentState(t_states, currentTimeStep * delta_t);
//                foutT.print(currentTimeStep * delta_t + "," + status + "\n");

            if (currentTimeStep%500 == 0) {
                try {
                    PrintStream fout = new PrintStream(new File("plot_data_pr3_t_" + currentTimeStep + "_java.csv"));
                    fout.println("x, Ts, Tf");
                    for (int n = 0; n < initNumberCells; n++) {
                        fout.println(n * deltaX + "," + currentTimeStepNode.prev.thisTS.Ts[n] + "," + currentTimeStepNode.prev.thisTS.Tf[n]);
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

            if (status == CHARGING) {
                uf = uf_charging;
                updateParameters();
                for (int i = 1; i < numberCells + 1; i++) {
                    //backup
                    /* cellCurrentTime.thisTS.Tf[n + 1] = cellCurrentTime.thisTS.Tf[n]
                        - uf * delta_t / deltaX * (cellCurrentTime.thisTS.Tf[n] - cellCurrentTime.prev.thisTS.Tf[n])
                        + alpha_f * delta_t / (deltaX * deltaX)
                        * (cellCurrentTime.next.thisTS.Tf[n] - 2 * cellCurrentTime.thisTS.Tf[n] + cellCurrentTime.prev.thisTS.Tf[n]);

                //System.out.print("n" + n + "Tf" + cellCurrentTime.thisTS.Tf[n] + "\n");
                cellCurrentTime.thisTS.Ts[n + 1] = cellCurrentTime.thisTS.Ts[n] + alpha_s * delta_t / (deltaX * deltaX)
                        * (cellCurrentTime.next.thisTS.Ts[n] - 2 * cellCurrentTime.thisTS.Ts[n] + cellCurrentTime.prev.thisTS.Ts[n]);
                     */
                    currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i - 1])
                            + alpha_f * delta_t / (deltaX * deltaX)
                            * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                    //System.out.print("n" + n + "Tf" + cellCurrentTime.thisTS.Tf[n] + "\n");
                    currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                            * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);

                    //now use inv matrix to get Tfn+1 and tsn+1
                    currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * currentTimeStepNode.thisTS.Ts_star[i];

                    currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];

                    // System.out.println("t" + currentTimeStep + "n" + i + "Tfs" + currentTimeStepNode.thisTS.Tf_star[i] + "Tss" + currentTimeStepNode.thisTS.Ts_star[i]);
                }

            }

//            if (status == DISCHARGING) {
//                uf = uf_discharging;
//                updateParameters();
//
//                for (int n = numberCells - 1; n >= numberCells; n--) {
//
//                    
//
//                }
//
//            }
            currentTimeStepNode = currentTimeStepNode.next;
            currentTimeStep++;
        }

//        } catch (FileNotFoundException ex) {
//            Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
//        }
        try {
            PrintStream fout = new PrintStream(new File("plot_data_pr3_final_java.csv"));
            fout.println("x, Ts, Tf");
            for (int n = 0; n < numberCells; n++) {
                fout.println(n * deltaX + "," + currentTimeStepNode.prev.thisTS.Ts[n] + "," + currentTimeStepNode.prev.thisTS.Tf[n]);
            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    private static int getCurrentState(ArrayList<Double> stateCycle, double t) {

        int numOfStates = stateCycle.size();
        double totalTimePerCycle = 0;
        int returnState = CHARGING;
        for (int i = 0; i < numOfStates; i++) {
            totalTimePerCycle += stateCycle.get(i);
        }
        double t_lastCycle = t % totalTimePerCycle;
        for (int i = 0; i < numOfStates && t_lastCycle > 0; i++) {
            t_lastCycle -= stateCycle.get(i);
            returnState = i;
        }

//         System.out.println("getCurrentState" + returnState);
        return returnState;
    }

    private static void updateParameters() {
        Re = epsilon * rho_f * ds * uf / mu_f;
        Pr = mu_f * Cp_f / kf;
        Nu_fs = 0.255 / epsilon * Math.pow(Pr, 1 / 3) * Math.pow(Re, 2 / 3);
        h_fs = Nu_fs * kf / ds;
        h = 1 / (1 / h_fs + ds / (10 * ks));
        hv = 6 * (1 - epsilon) * h / ds;

        hv_f = hv / (epsilon * rho_f * Cp_f);
        hv_s = hv / ((1 - epsilon) * rho_s * Cs);

        matrixM[0][0] = 1 + hv_f * delta_t;
        matrixM[1][0] = -hv_s * delta_t;
        matrixM[0][1] = -hv_f * delta_t;
        matrixM[1][1] = 1 + hv_s * delta_t;

        double OneDivDetM = 1 / (matrixM[0][0] * matrixM[1][1] - matrixM[0][1] * matrixM[1][0]);

        matrixM_inv[0][0] = OneDivDetM * matrixM[1][1];
        matrixM_inv[1][0] = -OneDivDetM * matrixM[1][0];
        matrixM_inv[0][1] = -OneDivDetM * matrixM[0][1];
        matrixM_inv[1][1] = OneDivDetM * matrixM[0][0];
    }
}
