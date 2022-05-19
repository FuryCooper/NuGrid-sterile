#!/bin/sh

mpirun -np 16 ./Gadget2 ./rui/Rui_twotypes/Rui_cmbonly-m020-xi000-lcdmgas-twotypes.param49
mpirun -np 16 ./Gadget2 ./rui/Rui_twotypes/Rui_cmbonly-m040-xi000-lcdmgas-twotypes.param49
mpirun -np 16 ./Gadget2 ./rui/Rui_twotypes/Rui_cmbonly-m060-xi000-lcdmgas-twotypes.param49
mpirun -np 16 ./Gadget2 ./rui/Rui_twotypes/Rui_cmbonly-m080-xi000-lcdmgas-twotypes.param49
mpirun -np 16 ./Gadget2 ./rui/Rui_twotypes/Rui_cmbonly-m100-xi000-lcdmgas-twotypes.param49

