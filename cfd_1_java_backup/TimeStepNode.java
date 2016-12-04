package cfd_1_java;


import cfd_1_java.TimeStep;

/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/**
 *
 * @author Kailin
 */
public class TimeStepNode {
     TimeStep thisTS;
            TimeStepNode next;
            TimeStepNode prev;

            TimeStepNode(int numCells) {
                this.thisTS = new TimeStep(numCells);
            }
    
}
