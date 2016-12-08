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

    static Double t_charging = 15000.;
    static Double t_idle_C_Dis = 0.;//2000.;
    static Double t_discharging = 0.;//500.;
    static Double t_idle_Dis_C = 0.;// 1000.;

    static int numberCycles = 1;

    static final int TS_per_sec = 1000;

    static double delta_t;

    static double Re;
    static double Pr;
    static double Nu_fs;
    static double h_fs;
    static double h;

    static double hv;

    static double hv_f;
    static double hv_s;

    static double alpha_f = 0;
    static double alpha_s = 0;
    static double uf = 0;
    static double uf_charging = 0.1;
    static double uf_discharging = -0.1;
    static double Tf_in = 873;

    static double height, diameter, initTemp;

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
        // double height, diameter, initTemp;
        final int numberCells;
        Scanner sc = new Scanner(System.in);
        height = sc.nextDouble();
        diameter = sc.nextDouble();
        numberCells = sc.nextInt();
        initTemp = sc.nextDouble();

        double deltaX = height * 1.0f / numberCells;

        updateParameters(CHARGING);

//        File file1 = new File("plot_data_pr1_JAVA.csv");
        ArrayList<Double> t_states = new ArrayList<>();
        t_states.add(t_charging);
        t_states.add(t_idle_C_Dis);
        t_states.add(t_discharging);
        t_states.add(t_idle_Dis_C);

        int timeStepsPerCycle = 1;
        int total_t_per_cycle = 1;
        for (double t : t_states) {
            total_t_per_cycle += t;
            timeStepsPerCycle = (int) (timeStepsPerCycle + (t * TS_per_sec));
        }
        timeStepsPerCycle -= 1;
        total_t_per_cycle -= 1;

        delta_t = 1.0 / TS_per_sec;//total_t_per_cycle * 1.0 / timeStepsPerCycle;

        int total_timesteps = timeStepsPerCycle * numberCycles;

        int initNumberCells = numberCells + 2;

        //creation of the linked list. each node represents a timestep.
        TimeStepNode T0 = new TimeStepNode(initNumberCells);
        TimeStepNode currentTimeStepNode;
        TimeStepNode TLast = new TimeStepNode(initNumberCells);
        T0.prev = null;

        currentTimeStepNode = T0;//.next;
        currentTimeStepNode.thisTS.Tf[0] = Tf_in;
        currentTimeStepNode.thisTS.Ts[0] = initTemp;
        for (int xi = 1; xi < initNumberCells; xi++) {
            currentTimeStepNode.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
            currentTimeStepNode.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
        }

//        double[][] matrixM = new double[2][2];
//        double[][] matrixM_inv = new double[2][2];
        int currentTimeStep = 0;

        while (currentTimeStep < total_timesteps + 1) {
            int status = getCurrentState(t_states, currentTimeStep * delta_t);
//                foutT.print(currentTimeStep * delta_t + "," + status + "\n");
            currentTimeStepNode.next = new TimeStepNode(initNumberCells);
            currentTimeStepNode.next.next = null;
            currentTimeStepNode.next.prev = currentTimeStepNode;

            currentTimeStepNode.next.thisTS.Tf[0] = Tf_in;
//            currentTimeStepNode.next.thisTS.Tf_star[0] = -initTemp;
//            currentTimeStepNode.next.thisTS.Ts_star[0] = -initTemp;
            for (int xi = 1; xi < initNumberCells; xi++) {
                currentTimeStepNode.next.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
                currentTimeStepNode.next.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
//                currentTimeStepNode.next.thisTS.Tf_star[xi] = -initTemp;
//                currentTimeStepNode.next.thisTS.Ts_star[xi] = -initTemp;
            }

            //different states
            if (status == CHARGING) { //charging
//                uf = uf_charging;

                updateParameters(CHARGING);
                for (int i = 0; i < numberCells + 1; i++) {
                    if (i == 0) {
                        currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i])
                                + alpha_f * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i]);
                        //System.out.print("n" + n + "Tf" + cellCurrentTime.thisTS.Tf[n] + "\n");
                        currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i]);
                        currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];
                    } else {
                        currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i - 1])
                                + alpha_f * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                        //System.out.print("n" + n + "Tf" + cellCurrentTime.thisTS.Tf[n] + "\n");
                        currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);

                        //now use inv matrix to get Tfn+1 and tsn+1
                        currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * 1.0 * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * 1.0 * currentTimeStepNode.thisTS.Ts_star[i];

                        currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];

//                        System.out.println("t" + currentTimeStep + "n" + i + "Tfs " + currentTimeStepNode.thisTS.Tf_star[i] + "Tss " + currentTimeStepNode.thisTS.Ts_star[i]);
//                        System.out.println("t" + (currentTimeStep + 1) + "n" + i + "TfNext " + currentTimeStepNode.next.thisTS.Tf[i] + "TsNext " + currentTimeStepNode.next.thisTS.Ts[i]);
//                        System.out.println("t" + currentTimeStep + "n" + i + "M00: " + matrixM_inv[0][0] + "M01: " + matrixM_inv[0][1]);
                    }
                }

            } else if (status == DISCHARGING) { //discharging
                //uf = uf_discharging;
                updateParameters(DISCHARGING);

                for (int i = numberCells; i > 0; i--) {
                    currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i + 1] - currentTimeStepNode.thisTS.Tf[i])
                            + alpha_f * delta_t / (deltaX * deltaX)
                            * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                    //System.out.print("n" + n + "Tf" + cellCurrentTime.thisTS.Tf[n] + "\n");
                    currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                            * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);

                    //now use inv matrix to get Tfn+1 and tsn+1
                    currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * currentTimeStepNode.thisTS.Ts_star[i];

                    currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];

                }

            } else { //idle
                currentTimeStepNode.next.thisTS = currentTimeStepNode.thisTS;
            }
            if (currentTimeStep % (500 * TS_per_sec) == 0 || currentTimeStep < 10) {
                System.out.println("current_t = " + currentTimeStep * delta_t);
            }

            if (currentTimeStep % (100 * TS_per_sec) == 0 || currentTimeStep < 10) {
                try {
                    PrintStream fout = new PrintStream(new File("plot_data_pr3_t_" + currentTimeStep * 1.0 / TS_per_sec + "_java.csv"));
                    fout.println("x, Tf, Ts , Tf* , Ts*");
                    for (int n = 0; n < initNumberCells; n++) {
                        fout.println(n * deltaX + "," + currentTimeStepNode.thisTS.Tf[n] + "," + currentTimeStepNode.thisTS.Ts[n] + "," + currentTimeStepNode.thisTS.Tf_star[n] + "," + currentTimeStepNode.thisTS.Ts_star[n]);
                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

            currentTimeStepNode = currentTimeStepNode.next;
            //System.out.println("NextTf1 " + currentTimeStepNode.thisTS.Tf[1] + "\n");
            currentTimeStepNode.prev = null;
            currentTimeStep++;
        }

//        } catch (FileNotFoundException ex) {
//            Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
//        }
        try {
            PrintStream fout = new PrintStream(new File("plot_data_pr3_final_java.csv"));
            fout.println("x, Tf, Ts");
            for (int n = 0; n < numberCells; n++) {
                fout.println(n * deltaX + "," + currentTimeStepNode.thisTS.Tf[n] + "," + currentTimeStepNode.thisTS.Ts[n]);
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

    private static void updateParameters(int status) {
        double area = Math.PI * diameter * diameter / 4.;
        uf_charging = 2.5*m_f_dot / (area * rho_f);
        uf_discharging = -uf_charging;
        if (status == CHARGING) {
            uf = uf_charging;
        } else if (status == DISCHARGING) {
            uf = uf_discharging;

        }
        Re = epsilon * rho_f * ds * uf / mu_f;
        Pr = mu_f * Cp_f / kf;
        Nu_fs = (0.255 / epsilon) * Math.pow(Pr, 0.3333333333) * Math.pow(Re, 0.6666666667);
        h_fs = Nu_fs * kf / ds;
        h = 1 / (1 / h_fs + ds / (10 * ks));
        hv = 6. * (1. - epsilon) * h / ds;

        hv_f = hv / (epsilon * rho_f * Cp_f);
        hv_s = hv / ((1. - epsilon) * rho_s * Cs);

        alpha_f = kf / (epsilon * rho_f * Cp_f);

        alpha_s = ks / ((1. - epsilon) * rho_s * Cs);

        alpha_f = 0;
        alpha_s = 0;

        matrixM[0][0] = 1 + hv_f * delta_t;
        matrixM[1][0] = -hv_s * delta_t;
        matrixM[0][1] = -hv_f * delta_t;
        matrixM[1][1] = 1 + hv_s * delta_t;

        // System.out.println("Matrix M = [" + matrixM[0][0] + " " + matrixM[0][1] + ";" + matrixM[1][0] + " " + matrixM[1][1] + "\n");
        double OneDivDetM = 1 / (matrixM[0][0] * matrixM[1][1] - matrixM[0][1] * matrixM[1][0]);

        matrixM_inv[0][0] = OneDivDetM * matrixM[1][1];
        matrixM_inv[1][0] = -OneDivDetM * matrixM[1][0];
        matrixM_inv[0][1] = -OneDivDetM * matrixM[0][1];
        matrixM_inv[1][1] = OneDivDetM * matrixM[0][0];
        //  System.out.println("Matrix M-1 = [" + matrixM_inv[0][0] + " " + matrixM_inv[0][1] + ";" + matrixM_inv[1][0] + " " + matrixM_inv[1][1] + "\n");
    }
}
