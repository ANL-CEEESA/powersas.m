function mpc = kundur
%KUNDUR
%   PSS(R)E 32 RAW created by rawd32  THU, JAN 30 2020  10:51
%   MODIFIED KUNDUR'S TWO-AREA TEST SYSTEM, DISTRIBUTED WITH ANDES
%   SEE THE BOOK "POWER SYSTEM STABILITY AND CONTROL" FOR ORIGINAL DATA
%
%   Converted by MATPOWER 7.1 using PSSE2MPC on 16-Mar-2023
%   from 'kundur.raw' using PSS/E rev 32 format.
%
%   WARNINGS:
%       Skipped 1 line of zone data.
%       Skipped 1 line of owner data.
%       Using default voltage magnitude limits: VMIN = 0.9 p.u., VMAX = 1.1 p.u.
%
%   See CASEFORMAT for details on the MATPOWER case file format.

%% MATPOWER Case Format : Version 2
mpc.version = '2';

%%-----  Power Flow Data  -----%%
%% system MVA base
mpc.baseMVA = 100;

%% bus data
%	bus_i	type	Pd	Qd	Gs	Bs	area	Vm	Va	baseKV	zone	Vmax	Vmin
mpc.bus = [
	1	3	0	0	0	0	1	1	32.6732	20	1	1.1	0.9;
	2	2	0	0	0	0	1	1	21.6548	20	1	1.1	0.9;
	3	2	0	0	0	0	2	1	11.2148	20	1	1.1	0.9;
	4	2	0	0	0	0	2	1	21.6398	20	1	1.1	0.9;
	5	1	0	0	0	0	1	0.98337	27.6488	230	1	1.1	0.9;
	6	1	0	0	0	0	1	0.96908	16.8176	230	1	1.1	0.9;
	7	1	1159	-73.5	0	0	1	0.95621	8.1662	230	1	1.1	0.9;
	8	1	1575	-89.9	0	0	2	0.954	-2.1295	230	1	1.1	0.9;
	9	1	0	0	0	0	2	0.96856	6.3774	230	1	1.1	0.9;
	10	1	0	0	0	0	2	0.98377	16.8036	230	1	1.1	0.9;
];

%% generator data
%	bus	Pg	Qg	Qmax	Qmin	Vg	mBase	status	Pmax	Pmin	Pc1	Pc2	Qc1min	Qc1max	Qc2min	Qc2max	ramp_agc	ramp_10	ramp_30	ramp_q	apf
mpc.gen = [
	1	745.861	143.612	600	0	1	900	1	900	0	0	0	0	0	0	0	0	0	0	0	0;
	2	700	300	600	-600	1	900	1	900	0	0	0	0	0	0	0	0	0	0	0	0;
	3	700	550	600	-600	1	900	1	900	0	0	0	0	0	0	0	0	0	0	0	0;
	4	700	-100	600	-600	1	900	1	900	0	0	0	0	0	0	0	0	0	0	0	0;
];

%% branch data
%	fbus	tbus	r	x	b	rateA	rateB	rateC	ratio	angle	status	angmin	angmax
mpc.branch = [
	5	6	0.005	0.05	0.075	0	0	0	0	0	1	-360	360;
	5	6	0.00501	0.05001	0.075	0	0	0	0	0	1	-360	360;
	6	7	0.002	0.02	0.03	0	0	0	0	0	1	-360	360;
	6	7	0.00201	0.02001	0.03	0	0	0	0	0	1	-360	360;
	7	8	0.02201	0.22001	0.33	0	0	0	0	0	1	-360	360;
	7	8	0.02202	0.22002	0.33	0	0	0	0	0	1	-360	360;
	7	8	0.022	0.22	0.33	0	0	0	0	0	1	-360	360;
	8	9	0.002	0.02	0.03	0	0	0	0	0	1	-360	360;
	8	9	0.00201	0.02001	0.03	0	0	0	0	0	1	-360	360;
	9	10	0.005	0.05	0.075	0	0	0	0	0	1	-360	360;
	9	10	0.00501	0.05001	0.075	0	0	0	0	0	1	-360	360;
	1	5	0.001	0.012	0	0	0	0	1	0	1	-360	360;
	2	6	0.001	0.012	0	0	0	0	1	0	1	-360	360;
	3	9	0.001	0.012	0	0	0	0	1	0	1	-360	360;
	4	10	0.001	0.012	0	0	0	0	1	0	1	-360	360;
];

%% bus names
mpc.bus_name = {
	'1           ';
	'2           ';
	'12          ';
	'11          ';
	'101         ';
	'102         ';
	'3           ';
	'13          ';
	'112         ';
	'111         ';
};
