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

    static Double t_charging = 3600.;
    static Double t_idle_C_Dis = 2000.;
    static Double t_discharging = 3600.;
    static Double t_idle_Dis_C = 2000.;

    static int t_charging_h = 0;
    static int t_discharging_h = 0;

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
    static double T_ref = 288.15; //Reference Temp for calculation of the exergy

    //declaration of the exerfy variables
    static double eta = 0;
    static double ex_d_out = 0;
    static double ex_d_in = 0;
    static double ex_c_out = 0;
    static double ex_c_in = 0;
    static double Q_charged = 0;
    static double Q_discharged = 0;
    static double capacity_factor = 0;
    static double Q_max = 0;

    static double height, diameter, initTemp;

    static double epsilon = 0.4;
    static double rho_f = 1835.6;
    static double rho_s = 2600;
    static double Cp_f = 1511.8;
    static double Cs = 900;
    static double ks = 2;
    static double kf = 0.52;
    static double m_f_dot = 10;
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
        numberCycles = sc.nextInt();
        t_charging_h = sc.nextInt();
        t_discharging_h = t_charging_h;

        t_charging = t_charging_h * 3600.;
        t_discharging = t_discharging_h * 3600.;

        sc.close();

        double deltaX = height * 1.0f / numberCells;

        updateParameters(CHARGING);
        Q_max = (epsilon * rho_f * Cp_f + (1 - epsilon) * rho_s * Cs) * Math.PI / 4 * diameter * diameter * height * (Tf_in - initTemp);

//        File file1 = new File("plot_data_pr1_JAVA.csv");
        int total_t_per_cycle = 24 * 3600;
        int timeStepsPerCycle = total_t_per_cycle * TS_per_sec;

        t_idle_C_Dis = total_t_per_cycle / 2 - t_charging;
        t_idle_Dis_C = t_idle_C_Dis;
        ArrayList<Double> t_states = new ArrayList<>();
        t_states.add(t_charging);
        t_states.add(t_idle_C_Dis);
        t_states.add(t_discharging);
        t_states.add(t_idle_Dis_C);

//        for (double t : t_states) {
//            total_t_per_cycle += t;
//            timeStepsPerCycle = (int) (timeStepsPerCycle + (t * TS_per_sec));
//        }
//        timeStepsPerCycle -= 1;
//        total_t_per_cycle -= 1;
        delta_t = 1.0 / TS_per_sec;//total_t_per_cycle * 1.0 / timeStepsPerCycle;

        int total_timesteps = timeStepsPerCycle * numberCycles;

        int initNumberCells = numberCells + 2;

        //creation of the linked list. each node represents a timestep.
        TimeStepNode T_Node_0 = new TimeStepNode(initNumberCells);
        TimeStepNode currentTimeStepNode;
        TimeStepNode TLast = new TimeStepNode(initNumberCells);
        T_Node_0.prev = null;

        currentTimeStepNode = T_Node_0;//.next;
        currentTimeStepNode.thisTS.Tf[0] = Tf_in;
        currentTimeStepNode.thisTS.Ts[0] = initTemp;
        for (int xi = 1; xi < initNumberCells; xi++) {
            currentTimeStepNode.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
            currentTimeStepNode.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
        }

//        double[][] matrixM = new double[2][2];
//        double[][] matrixM_inv = new double[2][2];
        int currentTimeStep = 0;
        int currentCycleNumer = 0;

        while (currentTimeStep < total_timesteps + 1) {
            int status = getCurrentState(t_states, currentTimeStep * delta_t);
//                foutT.print(currentTimeStep * delta_t + "," + status + "\n");
            currentTimeStepNode.next = new TimeStepNode(initNumberCells);
            currentTimeStepNode.next.next = null;
            currentTimeStepNode.next.prev = currentTimeStepNode;

//            currentTimeStepNode.next.thisTS.Tf_star[0] = -initTemp;
//            currentTimeStepNode.next.thisTS.Ts_star[0] = -initTemp;
            for (int xi = 1; xi < initNumberCells; xi++) {
                currentTimeStepNode.next.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
                currentTimeStepNode.next.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
                currentTimeStepNode.next.thisTS.Tf_star[xi] = initTemp;
                currentTimeStepNode.next.thisTS.Ts_star[xi] = initTemp;
            }

            //different states
            if (status == CHARGING) { //charging
                double Tf_integral = 0;
                double Ts_integral = 0;
                currentTimeStepNode.next.thisTS.Tf[0] = Tf_in;
                updateParameters(CHARGING);
                for (int i = 0; i < numberCells + 1; i++) {
                    if (i == 0) {
                        currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i])
                                + alpha_f * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i]);

                        currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i]);
                        currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];
                    } else {
                        currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i - 1])
                                + alpha_f * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                        currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);

                        //now use inv matrix to get Tfn+1 and tsn+1
                        currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * 1.0 * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * 1.0 * currentTimeStepNode.thisTS.Ts_star[i];

                        currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];
                    }
                    //calc the Q
                    Tf_integral += epsilon * rho_f * Cp_f * (currentTimeStepNode.thisTS.Tf[i] - initTemp)*deltaX;
                    Ts_integral += (1 - epsilon) * rho_s * Cs * (currentTimeStepNode.thisTS.Ts[i] - initTemp)*deltaX;
                    
                }
                Q_charged = Math.PI / 4 * diameter * diameter * (Tf_integral + Ts_integral);
            } else if (status == DISCHARGING) { //discharging
                double Tf_integral = 0;
                double Ts_integral = 0;
                updateParameters(DISCHARGING);
                //Boundary Condition at x = n, i = nu
                currentTimeStepNode.next.thisTS.Tf[0] = initTemp;
                currentTimeStepNode.thisTS.Tf[numberCells + 1] = initTemp;

                //for the BC at i = 0
                for (int i = 0; i < numberCells + 1; i++) {

                    //FTFS
                    if (i == 0) {
                        currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i + 1] - currentTimeStepNode.thisTS.Tf[i])
                                + alpha_f * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i]);

                        currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX) * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i]);
                    } else {
                        currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i + 1] - currentTimeStepNode.thisTS.Tf[i])
                                + alpha_f * delta_t / (deltaX * deltaX)
                                * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                        currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX) * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);
                    }

                    //now use inv matrix to get Tfn+1 and tsn+1
                    currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * 1.0 * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * 1.0 * currentTimeStepNode.thisTS.Ts_star[i];

                    currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];

                    //calc the Q
                    Tf_integral += epsilon * rho_f * Cp_f * (currentTimeStepNode.thisTS.Tf[i] - initTemp)*deltaX;
                    Ts_integral += (1 - epsilon) * rho_s * Cs * (currentTimeStepNode.thisTS.Ts[i] - initTemp)*deltaX;
                    
                }
Q_discharged = Math.PI / 4 * diameter * diameter * (Tf_integral + Ts_integral);
            } else { //idle
                currentTimeStepNode.next.thisTS = currentTimeStepNode.thisTS;
            }

            //Pr6 Exergy stuffs
            //discharge : out is at x=0 left end,in is at right end
            //charge: in is at x=0 left end, out is right end
            double T_left = currentTimeStepNode.thisTS.Tf[0];
            double T_right = currentTimeStepNode.thisTS.Tf[numberCells + 1];
            if (status == DISCHARGING) {
                ex_d_out = ex_d_out + delta_t * (m_f_dot * Cp_f * (T_left - T_ref - T_ref * Math.log(T_left / T_ref)));
                ex_d_in = ex_d_in + delta_t * (m_f_dot * Cp_f * (T_right - T_ref - T_ref * Math.log(T_right / T_ref)));
            } else if (status == CHARGING) {
                ex_c_out = ex_c_out + delta_t * (m_f_dot * Cp_f * (T_right - T_ref - T_ref * Math.log(T_right / T_ref)));
                ex_c_in = ex_c_in + delta_t * (m_f_dot * Cp_f * (T_left - T_ref - T_ref * Math.log(T_left / T_ref)));
            }
            currentCycleNumer = currentTimeStep / timeStepsPerCycle + 1;

            if (currentTimeStep % timeStepsPerCycle == (timeStepsPerCycle - 1)) {
                eta = (ex_d_out - ex_d_in) / (ex_c_in - ex_c_out);
                ex_d_out = 0;
                ex_d_in = 0;
                ex_c_out = 0;
                ex_c_in = 0;
            }
            double t_lastCycle = (currentTimeStep * TS_per_sec) % total_t_per_cycle;

            capacity_factor = (Q_charged - Q_discharged) / Q_max;

            if ((currentTimeStep + 1) % (1200 * TS_per_sec) == 0) {
                System.out.println("Cycle = " + currentCycleNumer);
                System.out.println("current_t = " + (currentTimeStep +1)* delta_t);
//                System.out.println("Matrix M-1 = [" + matrixM_inv[0][0] + " " + matrixM_inv[0][1] + ";" + matrixM_inv[1][0] + " " + matrixM_inv[1][1] + "\n");
            }

            if ((currentTimeStep + 1) % (1200 * TS_per_sec) == 0 || currentTimeStep < 10 || (4500 * TS_per_sec < currentTimeStep && currentTimeStep < (4500 * TS_per_sec + 20))) {
                try {
                    PrintStream fout = new PrintStream(new File("plot_data_pr3_t_" + (currentTimeStep * 1.0+1) / TS_per_sec + "_java.csv"));
                    fout.println("cycleNumber, x, Tf, Ts , Tf* , Ts*, sigma, d, u_f, Re, Pr, Nu_fs, hv, hv_f, hv_s, alpha_f, alpha_s,ex_d_out,ex_d_in,ex_c_out,ex_c_in, eta, Q_charged, Q_discharged,Q_max,capacity_factor");

                    for (int n = 0; n < initNumberCells; n++) {
                        fout.print(currentCycleNumer + ",");
                        fout.print(n * deltaX + "," + currentTimeStepNode.thisTS.Tf[n] + "," + currentTimeStepNode.thisTS.Ts[n] + "," + currentTimeStepNode.thisTS.Tf_star[n] + "," + currentTimeStepNode.thisTS.Ts_star[n] + ",");

                        fout.print(uf * delta_t / deltaX + ",");
                        fout.print(alpha_f * delta_t / (deltaX * deltaX) + ",");
                        fout.print(uf + ",");
                        fout.print(Re + ",");
                        fout.print(Pr + ",");
                        fout.print(Nu_fs + ",");
                        fout.print(hv + ",");
                        fout.print(hv_f + ",");
                        fout.print(hv_s + ",");
                        fout.print(alpha_f + ",");
                        fout.print(alpha_s + ",");
                        fout.print(ex_d_out + ",");
                        fout.print(ex_d_in + ",");
                        fout.print(ex_c_out + ",");
                        fout.print(ex_c_in + ",");
                        fout.print(eta + ",");
                        fout.print(Q_charged + ",");
                        fout.print(Q_discharged + ",");
                        fout.print(Q_max + ",");
                        fout.print(capacity_factor + "\n");

                    }
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
                }

            }

            currentTimeStepNode = currentTimeStepNode.next;
            //System.out.println("NextTf1 " + currentTimeStepNode.thisTS.Tf[1] + "\n");
            currentTimeStepNode.prev.prev = null;
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
        double area = epsilon * Math.PI * diameter * diameter / 4.;
        uf_charging = m_f_dot / (area * rho_f);
        uf_discharging = -uf_charging;
        uf = uf_charging;
        if (status == CHARGING) {
            uf = uf_charging;
        } else if (status == DISCHARGING) {
            uf = uf_discharging;
        }
        Re = epsilon * rho_f * ds * Math.abs(uf) / mu_f;
        Pr = mu_f * Cp_f / kf;
        Nu_fs = (0.255 / epsilon) * Math.pow(Pr, 0.3333333333) * Math.pow(Re, 0.6666666667);
        h_fs = Nu_fs * kf / ds;
        h = 1 / (1 / h_fs + ds / (10 * ks));
        hv = 6. * (1. - epsilon) * h / ds;

        hv_f = hv / (epsilon * rho_f * Cp_f);
        hv_s = hv / ((1. - epsilon) * rho_s * Cs);

        alpha_f = kf / (epsilon * rho_f * Cp_f);

        alpha_s = ks / ((1. - epsilon) * rho_s * Cs);

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
