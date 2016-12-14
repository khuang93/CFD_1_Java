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
    static final int writeOutFreq = 3600;

    static Double t_charging = 3600.;
    static Double t_idle_C_Dis = 2000.;
    static Double t_discharging = 3600.;
    static Double t_idle_Dis_C = 2000.;

    static int t_charging_h = 0;
    static int t_discharging_h = 0;

    static int numberCycles = 1;

    static final int TS_per_sec = 100;

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
    static double deltaT = 0;

    //declaration of the Pr6 variables
    static ArrayList<Double> eta = new ArrayList<Double>();
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
        //read in height, diameter, numberCells, initTemo, numberCycles and t_charging_h
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
        //calculation of Q_max for the capacity factor calculation
        Q_max = (epsilon * rho_f * Cp_f + (1 - epsilon) * rho_s * Cs) * Math.PI / 4 * diameter * diameter * height * (Tf_in - initTemp);

        int total_t_per_cycle = 24 * 3600;
        int timeStepsPerCycle = total_t_per_cycle * TS_per_sec;

        t_idle_C_Dis = total_t_per_cycle / 2 - t_charging;
        t_idle_Dis_C = t_idle_C_Dis;
        ArrayList<Double> t_states = new ArrayList<>();
        t_states.add(t_charging);
        t_states.add(t_idle_C_Dis);
        t_states.add(t_discharging);
        t_states.add(t_idle_Dis_C);

        delta_t = 1.0 / TS_per_sec;
        int total_timesteps = timeStepsPerCycle * numberCycles;

        int initNumberCells = numberCells + 1;

        //creation of the linked list. each node represents a timestep.
        TimeStepNode T_Node_0 = new TimeStepNode(initNumberCells);
        TimeStepNode currentTimeStepNode;
        TimeStepNode TLast = new TimeStepNode(initNumberCells);
        T_Node_0.prev = null;

        currentTimeStepNode = T_Node_0;
        currentTimeStepNode.thisTS.Tf[0] = Tf_in; //for the start is always charging
        currentTimeStepNode.thisTS.Ts[0] = initTemp;
        for (int xi = 1; xi < initNumberCells; xi++) {
            currentTimeStepNode.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
            currentTimeStepNode.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
        }

        int currentTimeStep = 0;
        double currentTime = currentTimeStep * delta_t;
        int currentCycleNumber = 0;
        double Tf_rightEnd_before_charge = initTemp;
        try {
            while (currentTimeStep < total_timesteps + 1) {
                int status = getCurrentState(t_states, currentTimeStep * delta_t);
                currentTimeStepNode.next = new TimeStepNode(initNumberCells);
                currentTimeStepNode.next.next = null;
                currentTimeStepNode.next.prev = currentTimeStepNode;

                for (int xi = 1; xi < initNumberCells; xi++) {
                    currentTimeStepNode.next.thisTS.Tf[xi] = initTemp;//Math.cos(k * xi * deltaX);
                    currentTimeStepNode.next.thisTS.Ts[xi] = initTemp;//Math.sin(k * xi * deltaX);
                    currentTimeStepNode.next.thisTS.Tf_star[xi] = initTemp;
                    currentTimeStepNode.next.thisTS.Ts_star[xi] = initTemp;
                }

                //different states
                switch (status) {
                    case CHARGING: {
                        //charging
                        //declaration of variables needed to calculate Q
                        double Tf_integral = 0;
                        double Ts_integral = 0;
                        //BC: Tf at x=0 always Tf_in for Charging
                        currentTimeStepNode.thisTS.Tf[0] = Tf_in;
                        updateParameters(CHARGING);
                        for (int i = 0; i <= numberCells; i++) {
                            if (i == 0) {
                                //left boundary
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i])
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i]);
                                currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];
                            } else if (i == numberCells) {
                                //right boundary
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i - 1])
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Ts[i] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);
                            } else {
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i - 1])
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);
                            }
                            //now use inv matrix to get Tfn+1 and tsn+1
                            currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * 1.0 * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * 1.0 * currentTimeStepNode.thisTS.Ts_star[i];

                            currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];
                            //calc the Q
                            Tf_integral += epsilon * rho_f * Cp_f * (currentTimeStepNode.thisTS.Tf[i] - initTemp) * deltaX;
                            Ts_integral += (1 - epsilon) * rho_s * Cs * (currentTimeStepNode.thisTS.Ts[i] - initTemp) * deltaX;

                        }
                        //Q(t) for charging
                        Q_charged = Math.PI / 4 * diameter * diameter * (Tf_integral + Ts_integral);
                        break;
                    }
                    case DISCHARGING: {
                        //discharging
                        double Tf_integral = 0;
                        double Ts_integral = 0;
                        updateParameters(DISCHARGING);
                        //Boundary Condition at x = n, i = nu
                        currentTimeStepNode.thisTS.Tf[numberCells] = initTemp;
                        for (int i = 0; i <= numberCells; i++) {

                            //FTFS
                            if (i == 0) {
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i + 1] - currentTimeStepNode.thisTS.Tf[i])
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX) * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i]);
                            } else if (i == numberCells) {
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] - uf * delta_t / deltaX * (currentTimeStepNode.thisTS.Tf[i] - currentTimeStepNode.thisTS.Tf[i])
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX) * (currentTimeStepNode.thisTS.Ts[i] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);
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
                            Tf_integral += epsilon * rho_f * Cp_f * (currentTimeStepNode.thisTS.Tf[i] - initTemp) * deltaX;
                            Ts_integral += (1 - epsilon) * rho_s * Cs * (currentTimeStepNode.thisTS.Ts[i] - initTemp) * deltaX;

                        }
                        //calculates Q(t) for discharge.
                        Q_discharged = Math.PI / 4 * diameter * diameter * (Tf_integral + Ts_integral);
                        break;
                    }
                    default:
                        //idle
                        updateParameters(IDLE_C_DIS);
                        //iterate over cells
                        //not the advection term is completely deleted
                        for (int i = 0; i <= numberCells; i++) {
                            if (i == 0) {
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i] + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i]);
                                currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];
                            } else if (i == numberCells) {
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i]
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Ts[i] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);

                            } else {
                                currentTimeStepNode.thisTS.Tf_star[i] = currentTimeStepNode.thisTS.Tf[i]
                                        + alpha_f * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Tf[i + 1] - 2 * currentTimeStepNode.thisTS.Tf[i] + currentTimeStepNode.thisTS.Tf[i - 1]);

                                currentTimeStepNode.thisTS.Ts_star[i] = currentTimeStepNode.thisTS.Ts[i] + alpha_s * delta_t / (deltaX * deltaX)
                                        * (currentTimeStepNode.thisTS.Ts[i + 1] - 2 * currentTimeStepNode.thisTS.Ts[i] + currentTimeStepNode.thisTS.Ts[i - 1]);

                            }
                            //now use inv matrix to get Tfn+1 and Tsn+1
                            currentTimeStepNode.next.thisTS.Tf[i] = matrixM_inv[0][0] * 1.0 * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[0][1] * 1.0 * currentTimeStepNode.thisTS.Ts_star[i];

                            currentTimeStepNode.next.thisTS.Ts[i] = matrixM_inv[1][0] * currentTimeStepNode.thisTS.Tf_star[i] + matrixM_inv[1][1] * currentTimeStepNode.thisTS.Ts_star[i];

                        }
                        //this will be updated each idle period and will be left at the value of the last idle time step before used to calculate deltaT
                        Tf_rightEnd_before_charge = currentTimeStepNode.thisTS.Tf[numberCells];
                        break;
                }

                //Pr6 Exergy efficiency
                //discharge : out is at x=0 left end,in is at right end
                //charge: in is at x=0 left end, out is right end
                double T_left = currentTimeStepNode.thisTS.Tf[0];
                double T_right = currentTimeStepNode.thisTS.Tf[numberCells];

                //integration over time for the exergy flux
                if (status == DISCHARGING) {
                    ex_d_out = ex_d_out + delta_t * (m_f_dot * Cp_f * (T_left - T_ref - T_ref * Math.log(T_left / T_ref)));
                    ex_d_in = ex_d_in + delta_t * (m_f_dot * Cp_f * (T_right - T_ref - T_ref * Math.log(T_right / T_ref)));

                } else if (status == CHARGING) {
                    ex_c_out = ex_c_out + delta_t * (m_f_dot * Cp_f * (T_right - T_ref - T_ref * Math.log(T_right / T_ref)));
                    ex_c_in = ex_c_in + delta_t * (m_f_dot * Cp_f * (T_left - T_ref - T_ref * Math.log(T_left / T_ref)));
                    //calculates the temperature increase. when the charging phase ends the variable will have the value of deltaT at the end of charging
                    //initially i thought this is the difference to initTemp. However it makes more sense to calculate the delta of T to the beginning of this charge cycle
                    deltaT = currentTimeStepNode.thisTS.Tf[numberCells] - Tf_rightEnd_before_charge;//initTemp;
                }

                currentCycleNumber = currentTimeStep / timeStepsPerCycle + 1;

                //will be executed after each cycle, calculates capacity factor and eta. writes out in File, then resets the variables to 0 if needed
                if (currentTimeStep % timeStepsPerCycle == (timeStepsPerCycle - 1)) {
                    eta.add((ex_d_out - ex_d_in) / (ex_c_in - ex_c_out));
                    ex_d_out = 0;
                    ex_d_in = 0;
                    ex_c_out = 0;
                    ex_c_in = 0;
                    capacity_factor = (Q_charged - Q_discharged) / Q_max;
                    System.out.println("delta_T = " + deltaT + " eta = " + eta.get(currentCycleNumber - 1) + " c_factor = " + capacity_factor);
                    PrintStream fout = new PrintStream(new File("deltaT_eta_cFactor_" + currentCycleNumber + ".csv"));
                    fout.println("delta_T,eta, c_factor");
                    fout.println(deltaT + "," + eta.get(currentCycleNumber - 1) + " , " + capacity_factor);
                    deltaT = 0;//reset deltaT
                    //Stopping Criteria for convergence of eta
                    if (currentCycleNumber > 2) {
                        if (eta.get(currentCycleNumber - 1) - eta.get(currentCycleNumber - 2) < 1E-4) {
                            break;
                        }
                    }
                }
                double t_lastCycle = (currentTimeStep * TS_per_sec) % total_t_per_cycle;

                //console output to show the progress of the current execution
                if ((currentTimeStep + 1) % (writeOutFreq * TS_per_sec) == 0) {
                    System.out.println("Cycle = " + currentCycleNumber);
                    System.out.println("current_t = " + (currentTimeStep + 1) * delta_t);
                }

                //Write out file with all the variables every 3600s
                if (((currentTimeStep + 1) % (writeOutFreq * TS_per_sec) == 0 || currentTimeStep == 0) || (currentTimeStep + 1) % timeStepsPerCycle == 0) {

                    PrintStream fout = new PrintStream(new File("data_t_" + (currentTimeStep * 1.0 + 1) / (TS_per_sec * 3600) + "h_tCharge_" + t_charging_h + ".csv"));
                    fout.println("cycleNumber, x, Tf, Ts , Tf* , Ts*,status, sigma, d, u_f, Re, Pr, Nu_fs,h_fs,h, hv, hv_f, hv_s, alpha_f, alpha_s,ex_d_out,ex_d_in,ex_c_out,ex_c_in, eta, Q_charged, Q_discharged,Q_max,capacity_factor");

                    for (int n = 0; n < initNumberCells; n++) {
                        fout.print(currentCycleNumber + ",");
                        fout.print(n * deltaX + "," + currentTimeStepNode.thisTS.Tf[n] + "," + currentTimeStepNode.thisTS.Ts[n] + "," + currentTimeStepNode.thisTS.Tf_star[n] + "," + currentTimeStepNode.thisTS.Ts_star[n] + "," + status + ",");

                        fout.print(uf * delta_t / deltaX + ",");
                        fout.print(alpha_f * delta_t / (deltaX * deltaX) + ",");
                        fout.print(uf + ",");
                        fout.print(Re + ",");
                        fout.print(Pr + ",");
                        fout.print(Nu_fs + ",");
                        fout.print(h_fs + ",");
                        fout.print(h + ",");
                        fout.print(hv + ",");
                        fout.print(hv_f + ",");
                        fout.print(hv_s + ",");
                        fout.print(alpha_f + ",");
                        fout.print(alpha_s + ",");
                        fout.print(ex_d_out + ",");
                        fout.print(ex_d_in + ",");
                        fout.print(ex_c_out + ",");
                        fout.print(ex_c_in + ",");
                        fout.print((eta.size() == 0 ? 0 : eta.get(eta.size() - 1)) + ",");
                        fout.print(Q_charged + ",");
                        fout.print(Q_discharged + ",");
                        fout.print(Q_max + ",");
                        fout.print(capacity_factor + "\n");

                    }
                }

                //increment of the time step
                currentTimeStepNode = currentTimeStepNode.next;
                currentTimeStepNode.prev.prev = null;
                currentTimeStep++;

            }

            //final output file
            PrintStream fout = new PrintStream(new File("plot_data_pr3_final_java.csv"));
            fout.println("cycleNumber, x, Tf, Ts , Tf* , Ts*,status, sigma, d, u_f, Re, Pr, Nu_fs, hv, hv_f, hv_s, alpha_f, alpha_s,ex_d_out,ex_d_in,ex_c_out,ex_c_in, eta, Q_charged, Q_discharged,Q_max,capacity_factor");

            for (int n = 0; n < initNumberCells; n++) {
                fout.print(currentCycleNumber + ",");
                fout.print(n * deltaX + "," + currentTimeStepNode.thisTS.Tf[n] + "," + currentTimeStepNode.thisTS.Ts[n] + "," + currentTimeStepNode.thisTS.Tf_star[n] + "," + currentTimeStepNode.thisTS.Ts_star[n] + "," + "final" + ",");

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
                fout.print((eta.size() == 0 ? 0 : eta.get(eta.size() - 1)) + ",");
                fout.print(Q_charged + ",");
                fout.print(Q_discharged + ",");
                fout.print(Q_max + ",");
                fout.print(capacity_factor + "\n");

            }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(CFD_1_Java.class.getName()).log(Level.SEVERE, null, ex);
        }

    }

    //returns int representing the current state beased on the current time t
    private static int getCurrentState(ArrayList<Double> stateCycle, double t) {

        int numOfStates = stateCycle.size();
        //hard coded total time per cycle for the design study, which equalls to 24h
        double totalTimePerCycle = 86400;
        int returnState = CHARGING;
        //if not hard coded, then must use the following:
//        for (int i = 0; i < numOfStates; i++) {
//            totalTimePerCycle += stateCycle.get(i);
//        }
        double t_lastCycle = t % totalTimePerCycle;
        for (int i = 0; i < numOfStates && t_lastCycle > 0; i++) {
            t_lastCycle -= stateCycle.get(i);
            returnState = i;
        }

        return returnState;
    }

    //updates all the parameters based on the current status
    private static void updateParameters(int status) {
        double area = epsilon * Math.PI * diameter * diameter / 4.;
        uf_charging = m_f_dot / (area * rho_f);
        uf_discharging = -uf_charging;
        uf = uf_charging;
        switch (status) {
            case CHARGING:
                uf = uf_charging;
                break;
            case DISCHARGING:
                uf = uf_discharging;
                break;
            default:
                uf = 0.;
                break;
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
