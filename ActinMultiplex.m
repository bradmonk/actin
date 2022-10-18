function [varargout] = ActinMultiplex(AMX,AMS,varargin)
%% ActinMultiplex2
%==========================================================================
% This function simulates actin dynamics using DMC / MCMC methods...
% https://en.wikipedia.org/wiki/Gillespie_algorithm
% https://en.wikipedia.org/wiki/Dynamic_Monte_Carlo_method
% https://en.wikipedia.org/wiki/Kinetic_Monte_Carlo
% 
% time is measured in Monte Carlo steps (MCS)
% http://www.roentzsch.org/SurfDiff/index.html
%==========================================================================
clc, close all;



if (~isdeployed)
    cd(fileparts(which(mfilename))); 
end


HOME = fileparts(which(mfilename));
cd(HOME);
PROJECT_PATH = genpath(HOME); 
addpath(PROJECT_PATH);






%==========================================================================
%% Actin Simulation Tools
%==========================================================================

if exist('AMX','var') && nargin > 0
    fprintf('Starting Actin Simulation \n');
else
    fprintf('Launching Actin GUI: ActinMultiplexGUI.m \n');

    ActinMultiplexGUI

    varargout={};
  return
end






%% ACTIN MULTIPLEX TIP DATA
%==========================================================================

ActSteps = AMX{4};

BTs = [];
AFMx = [];
Nsteps = ActSteps;
dt = AMX{47};

SaveTipsAfter = AMX{19};
SaveTipsRate = AMX{20};
SvTpMx = AMX{25};

doActCounts = AMX{35};
doDelOrig = AMX{39};
doDelOrigT = AMX{40};

%-- LIVE PLOTS --
doLive3DActinPlot = AMS{1};
doLiveHullPlot = AMS{2};
LiveTipMod = AMS{3};
LiveHullMod = AMS{4};

LivePlotMod = AMS{6};
TriHullMod  = AMS{7};
doTM = 0;
doFinalPlots  = AMS{12};



runTest = AMS{5};
if runTest
    Nsteps = 30000;
    LivePlotMod = 1000;
    SaveTipsAfter = 26000;
    SaveTipsRate = 10;
    doLiveHullPlot = 0;
    doLive3DActinPlot = 0;
end





% Angles Setup
%--------------------------------------------------------------------------
% unitsratio('rad', 'deg')
d2r = 1/(180/pi);
Ov = [2 29 57 85 112 140 168 195 223 251 278 306 334];
Ov = [Ov -15 -40 -65 -100 -125 -150 -175 -205 -230 -260 -295 -320 -350];



% Spine Dimensions
%--------------------------------------------------------------------------
SPYnXY = AMX{10}/2;    % Spy neck XY radius
SPYhXY = AMX{13}/2;    % Spy head XY radius
SPYhZN = AMX{11};      % Spy head north
SPYhZS = AMX{12};      % Spy head south


PSDproxy = AMX{15};
PSDproxyXY = AMX{27};
inPSD = SPYhZN - PSDproxy;

SPYH = [SPYhXY SPYhXY];
AcMx = zeros(SPYhXY*2/5,SPYhXY*2/5);

dims = [SPYnXY SPYhZN SPYhZS SPYhXY SPYhXY PSDproxy inPSD PSDproxyXY];





% INITIALIZE AND TAG STARTING FILAMENTS
%--------------------------------------------------------------------------

NStFils = AMX{42};

Actin = zeros(NStFils,20);  


%--------------------------------------------------------------------------
%                               Actin(Nfil,19)
% N  Xa  Xo  Xt  Ya  Yo  Yt  Za  Zo  Zt  MomID  ID  Fkd Born Died Lif MaxL MeanL Null Lgth
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
%--------------------------------------------------------------------------




% Starting Length Loc & Angles
StartMonos = AMX{41};
fXYo = AMX{45};
fZo  = AMX{43};
fXYa = AMX{46};
fZa  = AMX{44};

% Branching Angles
TPi = d2r*AMX{9};
PPi = 0;


% MAKES STARTING FILAMENTS
Actin = MakeStartFils(Actin,NStFils,StartMonos,d2r,fZa,fXYo,SPYhZN,SPYhZS,SPYhXY);



TagN = numel(Actin(:,1));
Actin(:,12) = 1:TagN;
TagN = TagN+1;





%% Actin Morphology
%--------------------------------------------------------------------------
% NOTES
%{
From: Andre Kamkin, Irina Kiseleva
Springer Science & Business Media, Nov 18, 2010 - Biochemistry - 395 pages
Chapter Authors: Luo and Robinson
2.2 Microstructures and Deformations of the Actin Cyctoskeleton
- Gactin is 5 nm in diameter
- Factin filaments are 8 nm wide with a left-handed helical morphology
- 13 actin monomers per pseudo-repeat
- 1 pseudo-repeat length of 37 nm 
- Alternatively the actin filament can be considered to have a right-handed 
helical structure with two strands slowly twisting around each other. 
Each actin monomer is rotated 166 degrees with respect to its nearest neighbors 
across the strand (Holmes 1990). Within the strand, subdomains 2 and 4 contact 
subdomains 1 and 3 in the next monomer in the strand, and each monomer reaches 
across to the other strand through a hydrophobic plug that links the two 
strands together. 

37 nm / 13 p = 2.85 nm/p
13 p / 37 nm = 0.35 p/nm

From: Lodish, Principles of Molecular Biology
Adapted from C. E. Schutt et al. 1993, Nature 365:810
- There are 2 strands with 14 units per strand, 
- so 28 monomers in each 360 degree turn,
- with a length of 72 nm per 360 degree turn

72 nm / 28 p = 2.57 nm/p
28 p / 72 nm = 0.39 p/nm

Factin particle size in filaments:
0.369 p/nm   (Factin particles per nm of filament)
2.71  nm/p    (nm of filament per Factin particle)

%}
%--------------------------------------------------------------------------


p_per_nm = 0.369;    % (Factin particles per nm of filament)
nm_per_p = 2.71;     % (nm of filament per Factin particle)
Actin(:,20) = Actin(:,1) .* nm_per_p;   % Store filament length



% TRIG: branch XYZ tip coordinates
Actin(:,4) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
Actin(:,7) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
Actin(:,10) = Actin(:,1) .* nm_per_p .* cos(Actin(:,8)) + Actin(:,9);


ACTs = Actin;
oActin = Actin;
DiedACTs = Actin(1,:);







%% Conversion Factors & VCP Function
%--------------------------------------------------------------------------
% Conversion Factors Notes
%{

pNuM = VCP(vol,uM,pN)

% [pNuM] = VCP(vol,uM,pN);
% INPUTS
% vol: volume in um^3
% uM: concentration in uM
% pN: particle count
% 
% enter zero for the unknown value
% if pN is unknown enter... pN = VCP(.1,10,0)
% if uM is unknown enter... uM = VCP(.1,0,6e5)
 
%}
%--------------------------------------------------------------------------



MOL = 6.022e23;     % 1 mol Avagadro's number
mol = 6e23;         % 1 mol rounded






%% SPINE MORPHOLOGY
%--------------------------------------------------------------------------
% Spine Morphology Notes
%{
%---------------------------------------------------
% A typical spine volume is 0.1 um^3 or 1e-16 L
% volume of cylinder
% V = pi * r^2 * h
%---------------------------------------------------


% MATH - Spine Volume

This spine has a neck volume equivalent to a cylendar of dimensions: 
50d x 200

And a head volume equivalent to a 3D disk with dimensions:
100d x 100


% Units
mol = 6e23;        % in N
mM = 1e-3;
uM = 1e-6;
upM = 1e3;
dnM = 1e-3;

% 1 cm^3 = 1 mL
% SI units conversion from m^3 to L are:
% 1 cm^3 = 1 mL
% To convert cubic volume to mL you may first need to convert
% by an order of magnitude; here's a reminder of cubic volume conversions:
% 1 m^3 = 1 cm^3 * .01^3;        % 1e-6
% 1 m^3 = 1 mm^3 * .001^3;        % 1e-9
% 1 cm^3 = 1 mm^3 * .1^3;        % 1e-3
% 1 cm^3 = 1 um^3 * .0001^3;    % 1e-12
% Thus, if dendritic spines have an average volume of 0.1 um^3
% that would be equivalent to X uL
%
% 0.1 um^3 * (1 cm^3 / 1e12 um^3) * (1 mL / 1 cm^3) * (1000 uL / 1 mL)
% 0.1*(1/1e12)*(1/1)*(1000/1)
% 0.1 um^3 = .1e-12 cm^3 = .1e-9 uL
%
% and 0.1 um^3 equivalent to X L
% 
% 0.1*(1/1e12)*(1/1)*(1/1000)
% 1e-16


% Spine volume (0.1 um^3) in L and uL
SpyV = 1e-16;    % in L
SpyVu = .1e-9;    % in uL

% Actin Polymerization (12 p/µM*s)
% Act_Na = 1e5;
% BeKanT = 1000;

% Actin Depolymerization (2 p/s)    
BeKdnT = 200;
DePSum = 0;

% Actin Poly Math
% Cytosolic concentration of actin in cells ranges from .1 to .5 mM
% Given a spine volume of 'SpyV' we can find how many actin monomers
% are in an average spine:
%
% .1 mM (mmol/L) * (1 mol / 1000 mmol) * SpyV = 1e-17 mol
% 1e-17 mol * (6e23 units / 1 mol) = 6e3 monomer units

molAct = AMX{16};
Act_N = molAct * (1/1000) * SpyV * 6e23; % 6e3 monomer units

% we can check our math starting with a set number of actin monomers
% and calculate the spine molarity (6e3 monomer units as an example):
% 
% 6e3 units/SpyV * (1 mol / 6e23 units) * (1000 mmol / 1 mol)
% 6e3/SpyV*(1/6e23) 

Gactin_uM = Act_N/SpyV * (1/6e23) * (1000/1);    % 1.6e-10 

% Gactin_uM = Act_N / SpyVu / mol;    % 1.6e-10 
% Act_N = Act_N / SpyVu / mol;
BeKa = 12 * Gactin_uM * dt;
BeKd = 2 * dt;

volume of cylinder
V = pi * r^2 * h

volume of a sphere
V = (4*pi*r^3)/3



%}
%--------------------------------------------------------------------------


Vneck = pi * SPYnXY^2 * SPYhZS;
Vhead = pi * SPYhXY^2 * (SPYhZN-SPYhZS);
SpyV = (Vneck+Vhead) * 1e-24;    % nm^3 to L conversion (nm: 1e-24; um: 1e-15)
SpyV = 1e-16;


VuMOL = SpyV * MOL / 1e6;

uMOLV = 1e6 / MOL / SpyV;

% uM = p / SpyV * (1/MOL) * 1e6;
% Pon = Ka / 1e6 * MOL * SpyV;





%% ACTIN CONCENTRATIONS
%--------------------------------------------------------------------------

uM_Act = AMX{16};                           % Actin uM
GActinN0 = uM_Act / 1e6 * MOL * SpyV;       % t0 N Gactin monomer units (6e4)
Gactin_uM = GActinN0 / SpyV *(1/MOL)*1e6;      % Check uM_Act == Gactin_uM

FActinN = sum(Actin(:,1));                  % current number FActins
GActinN = ceil(GActinN0) - FActinN;         % current number GActins

Nfi = numel(Actin(:,1));                    % current number of branches






%% ACTIN POLYMERIZATION & DEPOLYMERIZATION RATES
%--------------------------------------------------------------------------
% ACTIN NOTES
%{
Pollard:    
Actin+	ATP      ADPP       ADP	    | MEAN
BeKa: 	11.6     3.4        2.9	    | 6.0 p/µM*s
BeKd:   1.4      0.2		5.4	    | 2.3 p/s
PeKa:    	1.3      0.11       0.09	| 0.5 p/µM*s
PeKd: 	0.8      0.02       0.25	| 0.4 p/s
          BeKaATP  BeKaADPP   BeKaADP
          BeKdATP  BeKdADPP   BeKdADP
          PeKaATP  PeKaADPP   PeKaADP
          PeKdATP  PeKdADPP   PeKdADP

%}
%--------------------------------------------------------------------------

BeKaATP = AMX{60};  BeKaADPP = AMX{61};  BeKaADP = AMX{62};
BeKdATP = AMX{63};  BeKdADPP = AMX{64};  BeKdADP = AMX{65};
PeKaATP = AMX{66};  PeKaADPP = AMX{67};  PeKaADP = AMX{68};
PeKdATP = AMX{69};  PeKdADPP = AMX{70};  PeKdADP = AMX{71};

BeKaATD = mean([BeKaATP  mean([BeKaADPP  BeKaADP])]);
BeKdATD = mean([BeKdATP  mean([BeKdADPP  BeKdADP])]);
PeKaATD = mean([PeKaATP  mean([PeKaADPP  PeKaADP])]);
PeKdATD = mean([PeKdATP  mean([PeKdADPP  PeKdADP])]);

KaATD = roundn(mean([BeKaATD PeKaATD]),-2);
KdATD = roundn(mean([BeKdATD PeKdATD]),-2);


doKez = AMX{72};
KaBE = AMX{23};                 % Barbed End On-Rate Scalar
KdBE = AMX{24};                 % Barbed End Off-Rate Scalar
KaPE = AMX{58};                 % Pointed End On-Rate Scalar
KdPE = AMX{59};                 % Pointed End Off-Rate Scalar
KaTIP = (KaBE + KaPE) / 2;      % Tip Mean Ka On-Rate Scalar
KdTIP = (KdBE + KdPE) / 2;      % Tip Mean Kd Off-Rate Scalar

BeKa = KaBE * Gactin_uM * dt;   % Ka Barbed End ON-Rate
BeKd = KdBE * dt;               % Kd Barbed End OFF-Rate
PeKa = KaPE * Gactin_uM * dt;   % Ka Pointed End ON-Rate
PeKd = KdPE * dt;               % Kd Pointed End OFF-Rate




if doKez
	TKa = KaTIP;
	TKd = KdTIP;
else
	TKa = KaATD;
	TKd = KdATD;
end

fKa = TKa * Gactin_uM * dt;    % Ka Fil ON-Rate
fKd = TKd * dt;                % Kd Fil OFF-Rate





%% THYMOSIN ACTIN-SEQUESTERING VARIABLES
%--------------------------------------------------------------------------
% THYMOSIN REACTION RATE NOTES
%{
The Law of Mass Action (*LMA) describes the rate at which chemicals collide
and interact to form different chemical combinations. When two different
chemicals can collide to form a dimer product, and the dimer can dissociate
reversably back into the individual component reactants, is described as:

        Ka>
T + A <---->  TA
       <Kd  

Where
    T : thymosin                (thymosin monomers)
    A : actin                   (Gactin monomers)
    TA: thymosin-actin          (thymosin-actin dimers)
    Ka: forward rate constant
    Kd: backward rate constant


The rate equations describing the CHANGE in molecular concentration per dt are 
[https://en.wikipedia.org/wiki/Rate_equation]...

TA/dt = Ka[T][A] - Kd[TA]      (LMA forward reaction: TA accumulation)
A/dt  = Ka[T][A] - Kd[TA]      (LMA reverse reaction: A accumulation)
T/dt  = Ka[T][A] - Kd[TA]      (LMA reverse reaction: T accumulation)

For many reactions the rate is given by a power law such as:

r = k * [A]^x * [B]^y

where [A] and [B] express the concentration of the species A and B, respectively (usually in moles
per liter (molarity, M)). The exponents x and y are the partial reaction orders and must be
determined experimentally; they are often not equal to the stoichiometric coefficients. The constant
k is the rate coefficient or rate constant of the reaction. For elementary reactions, which consist
of a single step, the order equals the molecularity as predicted by collision theory. For example, a
bimolecular elementary reaction A + B products will be second order overall and first order in
each reactant, with rate equation r = k * [A] * [B]. For multistep reactions, the
order of each step equals the molecularity, but this is not generally true for the overall rate.

% DIFFUSION RATE ESTIMATION
% k = 1.38064852e-23; % Boltzmann Constant
% T = 310;            % approximately 98 Degrees Fahrenheit
% nu = 3;             % centipoise viscosity
% r = 3e-9;           % particle radius in meters (1 Angstrom == 1e-10 meters)
% 
% DR = k*T / (6*pi*nu*r) * (10^6^2); % microns squared per second

%}
%--------------------------------------------------------------------------


D_Thym_Basal    = AMX{94};
D_Thym_T2	    = AMX{96};
D_Thym_T3	    = AMX{97};
D_Thym_T4	    = AMX{98};
D_Thym_T5	    = AMX{99};

Pct_TA_Basal	= AMX{95};
Pct_TA_T2	    = AMX{100};
Pct_TA_T3	    = AMX{101};
Pct_TA_T4	    = AMX{102};
Pct_TA_T5	    = AMX{103};

ThymTimeStart   = AMX{104};


Ka_Thymosin_T2	= AMX{86};
Ka_Thymosin_T3	= AMX{87};
Ka_Thymosin_T4	= AMX{88};
Ka_Thymosin_T5	= AMX{89};

Kd_Thymosin_T2	= AMX{90};
Kd_Thymosin_T3	= AMX{91};
Kd_Thymosin_T4	= AMX{92};
Kd_Thymosin_T5	= AMX{93};


uM0_Thymosin  = AMX{73};      % Thymosin concentration (uM) total
N0_Thymosin   = round(uM0_Thymosin / 1e6 * MOL * SpyV);    % Thymosin count in spine (N)
Thymosin_uM   = N0_Thymosin / SpyV *(1/MOL)*1e6;    % Thymosin concentration in spine (uM)



THYM_ACT_uM = AMX{74};                               % [TA] Concentration (uM)
THYM_uM     = Thymosin_uM;                           % [T] Concentration (uM)
THYM_N      = round(THYM_uM / 1e6 * MOL * SpyV);     % Thymosin monomers (N)
THYM_ACT_N  = round(THYM_ACT_uM / 1e6 * MOL * SpyV); % Thymosin-Actin dimers (N)



Ka_Thymosin     = AMX{75};      % Thymosin Ka (basal)
Kd_Thymosin     = AMX{76};      % Thymosin Kd (basal)

THYM_N0         = THYM_N;
THYM_ACT_N0     = THYM_ACT_N;
Ka_Thymosin0    = Ka_Thymosin;
Kd_Thymosin0    = Kd_Thymosin;


tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dt; % Thymosin+Actin Association
tKd = Kd_Thymosin * THYM_ACT_uM * dt;         % Thymosin+Actin Dissociation






THYM_D      = D_Thym_Basal;  % Proportion of THYMOSIN replaced every second
Pct_TA_Rep  = Pct_TA_Basal;  % Proportion of Replacement THYMOSIN bound to actin;

ntt = 1;










%% ARP2/3 FILAMENT BRANCHING VARIABLES
%--------------------------------------------------------------------------
% ARP NOTES
%{
According to Smith & Gelles 2012:
Arp helps to nucleate new daughter filament branches 
at a rate of 2.5 to 9.7 df/(mM_Arp * s * ummf)
2.5 with no WASp and 9.7 with 300 nM WASp

At some intermediate WASp level, branching may proceed at
5 daughter filaments (df)
per mM of Arp
per second 
per micrometer of mother filament (ummf)

Example calculations with 9 uM Arp and 21600 nm of mf
2.5 df/(Arp_mM*s*ummf) * ((9/1000) * dt * (21600/1000))
9.7 df/(Arp_mM*s*ummf) * ((9/1000) * dt * (21600/1000))
5.0 df/(Arp_mM*s*ummf) * ((9/1000) * dt * (21600/1000))

5 * ((9/1000) * dt * (21600/1000))
5 * ((Arp_uM/1000) * dt * (nmmf/1000))



Friederich reports an Arp2/3 association rate of:
5.4e-4 um^-3 S^-1 (crossref from Carlsson 2004) but use a value of
1e-5 b/µM*s as a baseline parameter in their simulation software ActSimChem
ASRT:	1e-5 b/µM*s

These 3 values give branching rates that differ orders of magnitude:
5 * (9/1000) * 1	% .045
1e-5 * 9 * 1	    % .00009
5.4e-4 * (9^3)        % 0.39

9000
But the Gelles and Carlsson value are fairly close when the umf is 9
5 * (9/1000) * 9	% 0.40
5.4e-4 * (9^3)        % 0.39

which would be the case if the Factin concentration was 51 uM
(9000/2.7) / 1e-16 *(1/6e23)*1e6

5.4e-4 * Arp_uM^3 * dt

%}
%--------------------------------------------------------------------------




Arp_Sc = AMX{17} * dt;                    % Arp empirical branching rate scalar
Arp_uM = AMX{33};                        % Arp uM
nmmf = Actin(:,20);                        % Length of filaments (nm)

FArpN = NStFils;                        % #of F-Arp
GArpN = Arp_uM / 1e6 * 6e23 * SpyV;        % #of G-Arp
GArpN0 = ceil(GArpN);                    % #of G-Arp (starting)
Arp_uM = GArpN / SpyV *(1/6e23)*1e6;    % Check uM_Act == Gactin_uM

ArpAdd = AMX{29};                        % Add X units to new branches
ARPmax = AMX{22};                        % Maximum Arp branches

ArpBR  = Arp_Sc * ((Arp_uM/1000) .* (nmmf/1000));    % Arp branching rate (per fil)


% MATH - Arp Branching
%--------------------------------------------------------------------------
%{

Tnmmf = sum(Actin(:,20));                % Length of filaments (nm) Combined
ArpBRT = Arp_Sc * ((Arp_uM/1000) .* (Tnmmf/1000));    % Arp rate scalar    (total)

ArpN = 1e3;
ArpOn = 5;
ArpOff = 1;
Arp_uM = ArpN / SpyVu / mol;    % 1.6 - 16 uM
Arp_PR = ArpOn * Arp_uM * dt;
Arp_DR = ArpOff * dt;
ArpR = AMX{17}* dt;    % Arp activity rate
ArpR0 = ArpR;            % Arp activity rate
ArpScalar = AMX{21};    % Arp filament length scalar
%}







%% COFILIN & GENERAL DEPOLYMERIZATION VARIABLES
%--------------------------------------------------------------------------
% fdKd = fKd*10;          % Depoly rate when Fil is Fkd
% 
% CofR = AMX{18}* dt;     % cofilin activity rate
% CofS = AMX{21};         % cofilin delete Nunits
% CofN = AMX{30};         % cofilin delete Nunits
% delOr =AMX{14}* dt;     % delete from origin rate

CofSMax = AMX{26};
TheoMaxFact = ceil(SPYhZN / nm_per_p);
ScFil = ceil(TheoMaxFact/CofSMax);

uM0_Cofilin     = AMX{82};      % Cofilin concentration (uM) total
N0_Cofilin = uM0_Cofilin / 1e6 * MOL * SpyV;    % Cofilin count in spine (N)
Cofilin_uM = N0_Cofilin / SpyV *(1/MOL)*1e6;    % Cofilin concentration in spine (uM)


Ka_Cofilin     = AMX{83};      % Cofilin Ka (basal)
Ka_Cofilin_LTP = AMX{84};      % Cofilin Ka (LTP)
Ka_Cofilin_LTD = AMX{85};      % Cofilin Ka (LTD)

Factin_uM = FActinN / SpyV *(1/MOL)*1e6;    % FActin uM
KaCof = (dt * Ka_Cofilin * Factin_uM * Cofilin_uM); % Cofilin+Actin Association










% ADVANCED PARAMETERS (SAVE)
%--------------------------------------------------------------------------
%{

    	Abbreviations
-------------------------------------
(+)end: Barbed end of filament
(-)end: Pointed end of filament
mono: 	monomer, free protein monomer not associated with a filament
proto: 	protomer, polymerized protein or protein associated with a filament 
Pi: 	Phosphate
ATG: 	G-actin + ATP (same as ATM)
ADG: 	G-actin + ADP (same as ADM)
ATM: 	G-actin + ATP
ADM: 	G-actin + ADP 
ATF: 	F-actin proto + ATP
APF: 	F-actin proto + ADP-Pi
ADF: 	F-actin proto + ADP
FTB:     (+)ends terminating in ATP-actin
FPB:     (+)ends terminating in ADP-Pi-actin
FDB:     (+)ends terminating in ADP-actin
FTP:     (-)ends terminating in ATP-actin
FPP:     (-)ends terminating in ADP-Pi-actin
FDP:     (-)ends terminating in ADP-actin
CBM:     (+)end capping mono
CBF:     (+)end capper proto
CPM,     (-)end capping mono
CPF:     (-)end capper proto
FOM: 	Formin mono
FOF: 	Formin proto at (+)end
ARM: 	Arp in mono
ARF: 	Arp proto
FRP: 	Arp proto at (-)end (complex has bound actin at growth end)
FRB: 	Arp proto at (+)end (complex has no bound actin at growth end)



SNUC:	Spontanious Nucleation
FNUC:	Formin nucleation
CBNU:	Nucleation by barbed cap
CPNU:	Nucleation by pointed cap
ASTB:	Poly on rate of ATG at (+)end
ASDB:	Poly on rate of ADG at (+)end
ASTP:	Poly on rate of ATG at (-)end
ASDP:	Poly on rate of ATG at (-)end
DITB:	DPoly off rate of ATF at (+)end
DIPB:	DPoly off rate of APF at (+)end
DIDB:	DPoly off rate of ADF at (+)end
DITP:	DPoly off rate of ATF at (-)end
DIPP:	DPoly off rate of APF at (-)end
DIDP:	DPoly off rate of ADF at (-)end
TTOP:	ATP hydrolysis (ATF -> APF)
PTOD:	P release (APF -> ADF)
DTOT:	Recharge of G-actin by ATF
ASRT:	Arp association rate
DIRP:	Arp dissociation rate
FRGM:	Spontanious fragmentation
ANNL:	Annealing
ASFB:	Association of formin to (+)end
DIFB:	Dissocaition of formin
FASB:	Formin-aided filament growth
ASCB:	Capping of (+)end
DICB:	Uncapping of (-)end
ASCP:	Capping of (+)end
DICP:	Uncapping of (-)end


    	Default Rates
-------------------------------------
SNUC:	1e-8	p/(uM^2*s^1)
FNUC:	7e-5	p/(uM^3*s^1)
CBNU:	1e-5	p/(uM^3*s^1)
CPNU:	1e-5	p/(uM^3*s^1)
ASTB: 	11.5	p/µM*s
ASDB: 	3.8		p/µM*s
ASTP: 	1.3		p/µM*s
ASDP: 	0.16	p/µM*s
DITB:	0.90	p/s
DIPB:	0.90	p/s
DIDB:	1.50	p/s
DITP:	0.19	p/s
DIPP:	0.19	p/s
DIDP:	0.26	p/s
FRGM:	1.8e-8	p/s
ANNL:	1e-8	p/µM*s
TTOP:	0.30	p/s
PTOD:	2.6e-3	p/s
DTOT:	1e-2	p/s
ASRT:	1.0e-5	p/µM*s
DIRP:	1.0e-3	p/s
ASCB:	3.0		p/µM*s
DICB:	4e-4	p/s
ASCP:	1.0		p/µM*s
DICP:	1e-2	p/s
ASFB:	3.0		p/µM*s
DIFB:	1e-4	p/s
FASB:	110		p/µM*s



    	Rate Calculation Notes
-------------------------------------
exponents of µ^-1 == 1/µ
exponents of µ^-2 == 1/µ/µ
exponents of µ^-2 == 1/µ/µ/µ
[note these are sequential divisions, not 1/(µ/µ)]
Thus if something has a rate of [0.1 µM^-3 s^-1]
to calculate the current rate do [p=.1	uM=uM	s=dt]

	p/(uM^-3)/(s^-1)   or   p*(uM^3)*(s^1)
-------------------------------------



    	Default Rates
-------------------------------------
NUCLEATIONS
SNUC:	1e-8	p/(uM^2*s^1)
FNUC:	7e-5	p/(uM^3*s^1)
CBNU:	1e-5	p/(uM^3*s^1)
CPNU:	1e-5	p/(uM^3*s^1)

ACTIN ASSOCIATIONS
ASTB: 	11.5	p/µM*s
ASDB: 	3.8		p/µM*s
ASTP: 	1.3		p/µM*s
ASDP: 	0.16	p/µM*s

ACTIN DISSOCIATIONS
DITB:	0.90	p/s
DIPB:	0.90	p/s
DIDB:	1.50	p/s
DITP:	0.19	p/s
DIPP:	0.19	p/s
DIDP:	0.26	p/s

FRAGMENTATION/ANNEALING
FRGM:	1.8e-8	p/s
ANNL:	1e-8	p/µM*s

ATP/ADP
TTOP:	0.30	p/s
PTOD:	2.6e-3	p/s
DTOT:	1e-2	p/s

ARP2/3
ASRT:	1.0e-5	p/µM*s
DIRP:	1.0e-3	p/s

CAPPING
ASCB:	3.0		p/µM*s
DICB:	4e-4	p/s
ASCP:	1.0		p/µM*s
DICP:	1e-2	p/s

FORMIN
ASFB:	3.0		p/µM*s
DIFB:	1e-4	p/s
FASB:	110		p/µM*s








Abbreviations Full Description: 
F-actin, Filamentous actin; 
nSRF model, non-structurally-resolved filament model; 
SRF model, structurally-resolved filament model; 
Pi, Phosphate; 
ATM: Globular actin (monomeric form) with incorporated ATP; 
ADM: Globular actin (monomeric form) with incorporated ADP; 
ATF: Filamentous actin protomer (F-actin) with incorporated ATP; 
APF: Filamentous actin protomer (F-actin) with incorporated ADP-Pi; 
ADF: Filamentous actin protomer (F-actin) with incorporated ADP; 
FTB: Barbed ends of filaments, terminating by ATP-actin; 
FPB: Barbed ends of filaments, terminating by ADP-Pi-actin; 
FDB: Barbed ends of filaments, terminating by ADP-actin; 
FTP: Pointed ends of filaments, terminating by ATP-actin; 
FPP: Pointed ends of filaments, terminating by ADP-Pi-actin; 
FDP: Pointed ends of filaments, terminating by ADP-actin; 
CBM: Barbed-end capping protein (capper) in monomer (free) form; 
CBF: Barbed-end capper bound to filament; CPM, Pointed-end capper in monomer form; 
CPF: Pointed-end capper bound to filament; 
FOM: Formin in monomer (free) form; 
FOF: Formin, bound to filament barbed ends; 
ARM: Arp2/3 in monomer (free) form; 
ARF: Arp2/3 associated with filament; 
FRP: Arp2/3 associated with filament pointed end (fil terminates at ARP2/3); 
FRB: Arp2/3 associated with filament with no bound actins (fil terminates at ARP2/3); 

%}
%--------------------------------------------------------------------------



%% GET TOTAL PROTEIN COUNTS IN CLOSED SYSTEM


TOTAL_ACTIN = GActinN + FActinN + THYM_ACT_N;

TOTAL_THYMOSIN = THYM_N + THYM_ACT_N;

TOTAL_ARP = GArpN + FArpN;

TOTAL_COFILIN = N0_Cofilin;



TOTAL_ACTIN0 = TOTAL_ACTIN;

TOTAL_THYMOSIN0 = TOTAL_THYMOSIN;

TOTAL_ARP0 = TOTAL_ARP;

TOTAL_COFILIN0 = TOTAL_COFILIN;


%% LTP related variables

doActT  = AMX{48};
doArpT  = AMX{49};
doThymT = AMX{79};

T2start = AMX{50};
T3start = AMX{51};
T4start = AMX{52};
T5start = AMX{53};

T2pct   = AMX{54};
T3pct   = AMX{55};
T4pct   = AMX{56};
T5pct   = AMX{57};

Tbaseline       = AMX{77};
TbaselinePct    = AMX{78};


T0_TOTAL_ACTIN = GActinN + FActinN + THYM_ACT_N;
T2_TOTAL_ACTIN = T2pct * GActinN + FActinN + THYM_ACT_N;
T3_TOTAL_ACTIN = T3pct * GActinN + FActinN + THYM_ACT_N;
T4_TOTAL_ACTIN = T4pct * GActinN + FActinN + THYM_ACT_N;
T5_TOTAL_ACTIN = T5pct * GActinN + FActinN + THYM_ACT_N;



GArpT2 = GArpN0 * T2pct;
GArpT3 = GArpN0 * T3pct;
GArpT4 = GArpN0 * T4pct;
GArpT5 = GArpN0 * T5pct;








%% PREALLOCATE COUNTERS
nT = 1:Nsteps;
NumFils = nT;    % Number of branch filaments
NumFAct = nT;    % Number FActins
NumGAct = nT;    % Number GActins
NumCOFACT = nT;    % Number of Depoly events
NumARPACT = nT;    % Number of Poly events
NumdelFi = nT;
muFlength = nT;

ArpBRsum = nT;
NumAct_uM = nT;
NumfKa = nT;
NumFnowFtot = nT;
NumFArp = nT;    % Number of Arp in Filaments
NumGArp = nT;    % Number of Free Arp

NumHullV = zeros(3,round(Nsteps/LivePlotMod));

Num_nmmf = nT;
NumFded = nT;

NumGFActinN = nT;    % Number Unbound Actins
NumTHYM = nT;       % Number Unbound Thymosin
NumTHYMACT = nT;    % Number Thymosin+GActins
NumKaTHYM = nT;     % Ka Thymosin
NumKdTHYM = nT;     % Kd Thymosin
%----------------------------------------





%%      ANIMATED REAL-TIME FIGURE SETUP
close all

    doLiveLinePlot = AMS{10}; 
    LiveLinePlotMod = AMS{11};

% GET SCREEN SIZE AND CALCULATE DIMS FOR FIGURE
ssz = get(groot,'ScreenSize');
pos = [30 100 ssz(3)-100 ssz(4)-200];

% CREATE FIGURE 'f1'
fh1 = figure(1); 
    set(fh1,'OuterPosition',pos,'Color',[1 1 1]);

% CREATE AXES 'hax1' FOR MAIN ACTIN FILAMENT PLOT
hax1 = axes('Position',[.02 .08 .38 .50],'Color','none','NextPlot','replacechildren');
    rot0 = [5 0]; rota = rot0; azel = [-32 12]; vLim = [-32 12];
    AxLim = [-dims(4) dims(4) -dims(5) dims(5) 0 dims(2)].*1.2;
    axis(AxLim); view(vLim); axis vis3d; grid off;
    hold on

% CREATE AXES 'hax2' FOR MEMBRANE HULL PLOT
% hax2 = axes('Position',[.02 .70 .20 .29],'Color','none','NextPlot','replacechildren');
%     axis(AxLim); view(vLim); axis vis3d; grid off;
%     hold on
    
% CREATE AXES 'hax6' FOR FILAMENT COUNT LINE PLOT
hax2 = axes('Position',[.05 .70 .15 .25],'Color','none','NextPlot','replacechildren');
    axis([0 Nsteps*dt/60 0 100]); 
    xlabel('time (min)'); ylabel('actin particles');
    title('Mean Filament Length');
    hold on    

% CREATE AXES 'hax3' FOR ACTIN LEVELS LINE PLOT
hax3 = axes('Position',[.43 .58 .25 .35],'Color','none','NextPlot','replacechildren');
    axis([0 Nsteps*dt/60 0 TOTAL_ACTIN]);
    xlabel('time (min)'); ylabel('actin particles');
    title('G-actin vs F-actin');
    hold on

% CREATE AXES 'hax4' FOR ARP2/3 LEVELS LINE PLOT
hax4 = axes('Position',[.73 .58 .25 .35],'Color','none','NextPlot','replacechildren');
    axis([0 Nsteps*dt/60 0 TOTAL_ARP]);
    xlabel('time (min)'); ylabel('arp_{2/3} particles');
    title('Mobile (blue) vs Filamentous (red) Arp_{2/3}');
    hold on

% CREATE AXES 'hax5' FOR THYMOSIN LEVELS LINE PLOT
hax5 = axes('Position',[.43 .10 .25 .35],'Color','none','NextPlot','replacechildren');
    axis([0 Nsteps*dt/60 0 TOTAL_THYMOSIN]);
    xlabel('time (min)'); ylabel('thymosin particles');
    title('Free (red) vs Bound (blue) Thymosin');
    hold on

% CREATE AXES 'hax6' FOR FILAMENT COUNT LINE PLOT
hax6 = axes('Position',[.73 .10 .25 .35],'Color','none','NextPlot','replacechildren');
    axis([0 Nsteps*dt/60 0 TOTAL_ARP]); 
    xlabel('time (min)'); ylabel('filament count');
    title('Number of Actin Filaments');
    hold on

% CREATE AXES 'hax7' FOR TEXT OR BAR PLOT
hax7 = axes('Position',[.23 .70 .15 .28],'Color','none','XTick', [],'YTick',[]);
    axis([0 4 0 100]);
    hold on

    
figure(fh1)
delete(hax7.Children);
hax7lims = axis;
h7y = hax7lims(4);

spf1  = sprintf(' Tot Actin:  % 9.6g  ',TOTAL_ACTIN);
spf2  = sprintf(' GActin:     % 9.6g  ',GActinN);
spf3  = sprintf(' FActin:     % 9.6g  ',FActinN);
spf4  = sprintf(' Thymosin:   % 9.6g  ',THYM_N);
spf5  = sprintf(' ThymAct:    % 9.6g  ',THYM_ACT_N);
spf6  = sprintf(' GArpN:      % 9.6g  ',GArpN);
spf7  = sprintf(' FArpN:      % 9.6g  ',FArpN);
spf8  = sprintf(' Ka Actin:   % 9.6g  ',fKa);
spf9  = sprintf(' Ka Thymo:   % 9.6g  ',tKa);
spf10 = sprintf(' Kd Thymo:   % 9.6g  ',tKd);
text(.1,(h7y*.92),spf1 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.84),spf2 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.76),spf3 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.68),spf4 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.60),spf5 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.52),spf6 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.44),spf7 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.36),spf8 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.28),spf9 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
text(.1,(h7y*.20),spf10,'FontSize',14,'FontName','FixedWidth','Parent',hax7);


%% PRE-PLOT FILAMENTS



ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,~] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,~] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);

P3x = [Actin(:,3) Actin(:,4)]';
P3y = [Actin(:,6) Actin(:,7)]';
P3z = [Actin(:,9) Actin(:,10)]';
P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;


axes(hax1)

ph_spyTip = scatter3(hax1,SPYTips(:,1)', SPYTips(:,2)', SPYTips(:,3)',20,'ob');
ph_psdTip = scatter3(hax1,PSDTips(:,1)', PSDTips(:,2)', PSDTips(:,3)',20,'or');
ph_fil = line('XData',P3x(:),'YData',P3y(:),'ZData',P3z(:),'Parent',hax1);

set(ph_spyTip,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph_psdTip,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph_fil,'LineStyle','-','Color',[.5 .5 .5],'LineWidth',1.5);
view(azel+rota)





%% PRE-LOOP SETTINGS FOR MAIN LOOP

mempressure = 10;

axes(hax1) % make axes 'hax1' top-level plot

%==========================================================================
%%                    	MAIN OUTER LOOP
for nT = 1:Nsteps
%--------------------------------------------------------------------------


	if doActT
      if nT==T2start; TOTAL_ACTIN = T2_TOTAL_ACTIN; end
      if nT==T3start; TOTAL_ACTIN = T3_TOTAL_ACTIN; end
      if nT==T4start; TOTAL_ACTIN = T4_TOTAL_ACTIN; end
      if nT==T5start; TOTAL_ACTIN = T5_TOTAL_ACTIN; end
	end
	if doArpT
      if nT==T2start; GArpN0 = GArpT2; end
      if nT==T3start; GArpN0 = GArpT3; end
      if nT==T4start; GArpN0 = GArpT4; end
      if nT==T5start; GArpN0 = GArpT5; end
	end
	if doThymT      
      if nT==T2start; Kd_Thymosin = Kd_Thymosin_T2; Ka_Thymosin = Ka_Thymosin_T2; end
      if nT==T3start; Kd_Thymosin = Kd_Thymosin_T3; Ka_Thymosin = Ka_Thymosin_T3; end
      if nT==T4start; Kd_Thymosin = Kd_Thymosin_T4; Ka_Thymosin = Ka_Thymosin_T4; end
      if nT==T5start; Kd_Thymosin = Kd_Thymosin_T5; Ka_Thymosin = Ka_Thymosin_T5; end
	end
    
	NFact = numel(Actin(:,1));    % Current #of Filaments 
      %  fN - used in inner loop below as "for fN=1:NFact..."
      %  Nfi - used in outer loop below as "Nfi = numel(Actin(:,1));"
      %  note: do not re-compute NFact anywhere else but here
      %  use Nfi in final counters since filaments are being added/deleted
    

	ACTdepoly = 0; COFdepoly = 0; ARPpoly = 0; ACTpoly = 0;
    
    
    %----------------------------------------------
	Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];

    % radial distance to spine shaft membrane
	XYtipLoc = sqrt(Tip_xyz(:,1).^2 + Tip_xyz(:,2).^2);

	ZtipInHead = Tip_xyz(:,3) >= SPYhZS;
	ZtipInNeck = Tip_xyz(:,3) < SPYhZS;

	XYneckOut = XYtipLoc > SPYnXY;        % Logical array of filaments beyond SPYnXY
	XYheadOut = XYtipLoc > SPYhXY;        % Logical array of filaments beyond SPYhXY
	ZtopOut = Tip_xyz(:,3) > SPYhZN;    % Logical array of filaments above SPYhZN
	ZbotOut = Tip_xyz(:,3) < 0;            % Logical array of filaments below zero
    

	TipOut =((XYneckOut & ZtipInNeck)+(XYheadOut & ZtipInHead)+ZtopOut+ZbotOut)>0;
	TipOK = ~TipOut;
    %----------------------------------------------
    

	LngXYneckOut = (XYtipLoc - SPYnXY) .* ZtipInNeck;        % how far are filaments beyond SPYnXY
	LngXYheadOut = (XYtipLoc - SPYhXY) .* ZtipInHead;        % how far are filaments beyond SPYhXY
	LngZtopOut = (Tip_xyz(:,3) - SPYhZN) .* ZtipInHead;     % how far are filaments above SPYhZN
    
    LngXYneckOut(LngXYneckOut<0) = 0;
    LngXYheadOut(LngXYheadOut<0) = 0;
    LngZtopOut(LngZtopOut<0) = 0;
    
    LngOut = LngXYneckOut + LngXYheadOut + LngZtopOut;
    
    LO = LngOut ./ 100 .* mempressure;    
    

    
    Actin(:,20) = Actin(:,1) .* nm_per_p;
    

    % assure no negative actin values
	Actin(:,1) = Actin(:,1) .* (Actin(:,1)>0);
	Actin(:,1) = Actin(:,1) .* (Actin(:,10) > -20);

    
    % MCMC requires a random filament at each decision step
	rN1 = randi(NFact,1,NFact); % actin polymerization
    rN2 = randi(NFact,1,NFact); % arp branching
    rN3 = randi(NFact,1,NFact); % actin depolymerization
    rN4 = randi(NFact,1,NFact); % cofilin severing




    %======================================================================
    %                    	MAIN INNER LOOP
	for fN = 1:NFact
    %----------------------------------------------------------------------
    
        rv=rand(9,1); % Generate a few random vaules from uniform{0:1}


        
        % POLYMERIZATION
        %------------------------------------------------------------------
        aN = rN1(fN); % Get random filament
        if nT < 10000
            

            if  (fKa > rv(1)) && TipOK(aN) && (GActinN>1)
                if fKa > NFact
                    Actin(aN,1) = Actin(aN,1) + 2;
                    ACTpoly = ACTpoly + 2;
                else
                    Actin(aN,1) = Actin(aN,1) + 1;
                    ACTpoly = ACTpoly + 1;
                end
            end
        
        
        else
        
        
            if  (fKa-LO(aN) > rv(1)) && ~ZbotOut(aN) && (GActinN>1)
                if fKa > NFact
                    Actin(aN,1) = Actin(aN,1) + 2;
                    ACTpoly = ACTpoly + 2;
                else
                    Actin(aN,1) = Actin(aN,1) + 1;
                    ACTpoly = ACTpoly + 1;
                end
            end
        
        
        end
        
        



        
        % BRANCHING
        %------------------------------------------------------------------
        aN = rN2(fN);
     	if  (ArpBR(aN) > rv(2)) && (GActinN>ArpAdd)

        	fNf = NFact+1;

            
        	Actin(fNf,1) = ArpAdd;        % create branch: add N(ArpAdd) subunits to new branch
        	Actin(fNf,11)=Actin(aN,12);    % tag branch with MomID
        	Actin(fNf,12) = TagN;        % tag branch with ID
        	Actin(fNf,14) = nT;            % tag branch with Born time (nT)
        	TagN = TagN+1;


        	Nmono = Actin(aN,1);        % N monomers in mother filament
        	Rmm = ceil(Nmono * rand);    % Random mother monomer (== vector_length) along segment



            % Rotational point of current branch
        	Ct_x = (Rmm) .* nm_per_p * sin(Actin(aN,8)) * cos(Actin(aN,2)) + Actin(aN,3);
        	Ct_y = (Rmm) .* nm_per_p * sin(Actin(aN,8)) * sin(Actin(aN,2)) + Actin(aN,6);
        	Ct_z = (Rmm) .* nm_per_p * cos(Actin(aN,8)) + Actin(aN,9);



            % Origin of current branch
        	Co_x = Actin(aN,3);    % X origin (old branch)
        	Co_y = Actin(aN,6);    % Y origin (current branch)
        	Co_z = Actin(aN,9);    % Z origin (current branch)
            
            
        	Po = [Co_x; Co_y; Co_z];
        	Pt = [Ct_x; Ct_y; Ct_z];
        	Pv = Pt - Po;
            
        	tL = sqrt(sum((Pv).^2));        % Length of vector PoPt (aka Pv)
        	Pu = (Pv) ./ tL;                % Unit vector of Pv

        	tTheta = acos(Pu(3))  +TPi;        % angle theta
        	tPhi = atan2(Pu(2),Pu(1)) +PPi;    % angle phi

        	x = sin(tTheta) * cos(tPhi);
        	y = sin(tTheta) * sin(tPhi);
        	z = cos(tTheta);
        	Pr = [x;y;z]+Pt;
            
            
            % Random rotational angle (theta)
        	Otta = Ov(randi(26,1));        
            

            dv = [Pt(1);Pt(2);Pt(3)] - [Po(1);Po(2);Po(3)];

            u=dv(1);  
            v=dv(2);
            w=dv(3);



            % Rotation Matrix Function
            Pn = RotateVec(Pr(1),Pr(2),Pr(3),Pt(1),Pt(2),Pt(3),u,v,w,Otta);


            

            % Lengths & Angles
        	Po2 = Pt;
        	Pt2 = Pn;
        	Pv2 = Pt2-Po2;
        	tL2 = sqrt(sum((Pv2).^2));        % Length of vector PoPt (aka Pv)
        	Pu2 = (Pv2) ./ tL2;                % Unit vector of Pv
        	tTheta2 = acos(Pu2(3));            % angle theta	+TPi;
        	tPhi2 = atan2(Pu2(2),Pu2(1));    % angle phi	    +PPi;
            
            
            % New branch Angles
            %-------------------
        	Actin(fNf,2) = tPhi2;
        	Actin(fNf,8) = tTheta2;
            
            % New branch Origin
        	Actin(fNf,3) = Po2(1);
        	Actin(fNf,6) = Po2(2);
        	Actin(fNf,9) = Po2(3);
            
            % New branch Tip
        	Actin(fNf,4) = (ArpAdd) .* nm_per_p * sin(tTheta2) * cos(tPhi2) + Po2(1);
        	Actin(fNf,7) = (ArpAdd) .* nm_per_p * sin(tTheta2) * sin(tPhi2) + Po2(2);
        	Actin(fNf,10) = (ArpAdd) .* nm_per_p * cos(tTheta2) + Po2(3);
            
        	ARPpoly = ARPpoly + ArpAdd;
            
    	end

        
        
        
        
        
        
        % ACTIN PASSIVE DEPOLYMERIZATION
        %------------------------------------------------------------------
        aN = rN3(fN);
    	if (fKd > rv(3))

        	Actin(aN,1) = Actin(aN,1)-1;
        	ACTdepoly = ACTdepoly + 1;
    	end
        


        
        
        % COFILIN ASSISTED DEPOLYMERIZATION
        %------------------------------------------------------------------
        aN = rN4(fN);
    	if  KaCof > rv(4) && Actin(aN,1) > 10
            
            ActBoundCof = randi(Actin(aN,1));    % cofilin binds random Factin along filament

            ActDelCof = Actin(aN,1)-ActBoundCof; % depoly all Factin beyond cofilin site

            Actin(aN,1) = ActBoundCof;           % remaining Factin count == cofilin site
            
        	COFdepoly = COFdepoly + ActDelCof;   % keep tally of cofilin depolymerized Factin

            GActinN = GActinN + ActDelCof;       % add depoly actin to Gactin pool (temporary)

    	end



        

        % THYMOSIN+ACTIN ASSOCIATION
        %------------------------------------------------------------------
        if nT > ThymTimeStart
        if (tKa > rv(8)) && (GActinN>1) && (THYM_N>1)

            THYM_N = THYM_N - 1;
            THYM_ACT_N = THYM_ACT_N + 1;
            GActinN = GActinN - 1;

        end

        % THYMOSIN+ACTIN DISSOCIATION
        %------------------------------------------------------------------
        if (tKd > rv(9))  && (THYM_ACT_N>1)

            THYM_N = THYM_N + 1;
            THYM_ACT_N = THYM_ACT_N - 1;
            GActinN = GActinN + 1;

        end
        end

        
    


        % ADJUST RATE VALUES
        %------------------------------------------------------------------
        if ~mod(fN,10)

            % uM = p / SpyV * (1/MOL) * 1e6;

            % ACTIN VALUES
            FActinN = sum(Actin(:,1));                    % Current #of FActins

            TA_plus_Factin_N = FActinN + THYM_ACT_N;    % Current #of Thymosin:Actin + Factin

            GActinN = TOTAL_ACTIN - TA_plus_Factin_N;   % Current #of GActins

            Factin_uM = FActinN / SpyV * (1/MOL) * 1e6; % FActin uM
            Gactin_uM = GActinN / SpyV * (1/MOL) * 1e6; % Gactin uM

            fKa = TKa * Gactin_uM * dt;                    % actin polymerization rate


            % ARP VALUES
            FArpN = numel(Actin(:,1));              % #of F-Arp
            GArpN = GArpN0 - FArpN;                 % #of G-Arp
            Arp_uM = GArpN / SpyV * (1/MOL) * 1e6;  % Arp uM
            nmmf = Actin(:,20);                                   % Length of each filament (nm)
            Fact_Branch_uM = Actin(:,1) ./ SpyV * (1/MOL) * 1e6;  % conc. of each filament (uM)
            ArpBR = Arp_Sc * (Arp_uM/1000 .* Fact_Branch_uM);


            % THYMOSIN VALUES
            THYM_uM     = THYM_N     / SpyV * (1/MOL) * 1e6;     % Thymosin uM
            THYM_ACT_uM = THYM_ACT_N / SpyV * (1/MOL) * 1e6;     % Thymosin+Actin uM

            tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dt/NFact;  % uM associated per second
            tKd = Kd_Thymosin * THYM_ACT_uM * dt/NFact;          % uM dissociated per second

            TA_Non  = round(tKa / 1e6 * MOL * SpyV);
            TA_Noff = round(tKd / 1e6 * MOL * SpyV);


            % COFILIN VALUES
            KaCof = (dt/NFact * Ka_Cofilin * Factin_uM * Cofilin_uM);  % Cofilin+Actin Association

        end








	end % MAIN INNER LOOP
    %----------------------------------------------------------------------
    
    



    % ADJUST RATE VALUES
    %----------------------------------------------------------------------    
        
    % ACTIN VALUES
	FActinN = sum(Actin(:,1));                   % Current #of FActins

    TA_plus_Factin_N = FActinN + THYM_ACT_N;   % Current #of Thymosin/actin + Factin

    GActinN = TOTAL_ACTIN - TA_plus_Factin_N;  % Current #of GActins

    Factin_uM = FActinN / SpyV * (1/MOL) * 1e6;               % FActin uM
	Gactin_uM = GActinN / SpyV * (1/MOL) * 1e6;               % Gactin uM
	fKa = TKa * Gactin_uM * dt;                   % actin polymerization rate
    

    % ARP VALUES
    %----------------------------------------------------------------------
    FArpN = numel(Actin(:,1));              % #of F-Arp
    GArpN = GArpN0 - FArpN;                 % #of G-Arp
    Arp_uM = GArpN / SpyV * (1/MOL) * 1e6;                 % Arp uM
    nmmf = Actin(:,20);                     % Length of each filament (nm)
    % Arp Branching Rate is dependent on uM Arp and lengths of existing filaments
    Fact_Branch_uM = Actin(:,1) ./ SpyV * (1/MOL) * 1e6;   % conc. of each filament (uM)
    ArpBR = Arp_Sc * (Arp_uM/1000 .* Fact_Branch_uM);
    


    
       
    % COFILIN VALUES
    %----------------------------------------------------------------------
    KaCof = (dt/NFact * Ka_Cofilin * Factin_uM * Cofilin_uM);  % Cofilin+Actin Association

 
    
    
    % THYMOSIN VALUES
    %----------------------------------------------------------------------
    THYM_uM     = THYM_N     / SpyV * (1/MOL) * 1e6;    % Thymosin uM
    THYM_ACT_uM = THYM_ACT_N / SpyV * (1/MOL) * 1e6;    % Thymosin+Actin uM
    
    tKa = Ka_Thymosin * THYM_uM * Gactin_uM * dt;     % Thymosin+Actin association rate
    tKd = Kd_Thymosin * THYM_ACT_uM * dt;             % Thymosin+Actin dissociation rate
    
    
    
    
    % THYMOSIN DIFFUSION DURING LTP
    %----------------------------------------------------------------------
    if doThymT
        if nT==T2start
            TAeq = THYM_ACT_N;
            Teq = THYM_N;
            PctTeq = Teq / (TAeq + Teq);
            PctTAeq = TAeq / (TAeq + Teq);
        end
      if nT>T2start && nT<T3start % && ntt > 1/dt
          
        % Nt = No * exp(-Lam*dt); % [1]
        % Nt = No * exp(-dt/Tao); % [2]
        % 
        % Lam     : exponential decay constant
        % Tao     : exponential mean element lifetime
        % Nt      : quantity remaining at time t
        % No      : quantity at start
        % dt      : elapsed time from start


        NrepT = round(THYM_N * THYM_D);
        NrepTA = round(THYM_ACT_N * THYM_D);
        NrepTAT = NrepT + NrepTA;

        THYM_N = THYM_N + -NrepT + round(PctTeq*NrepTAT);
        THYM_ACT_N = THYM_ACT_N + -NrepTA + round(PctTAeq*NrepTAT);
        
        TOTAL_ACTIN = GActinN + FActinN + THYM_ACT_N;
        hax3.YLim = [0 TOTAL_ACTIN+100];
        
        ntt = 0;
                
      end
      ntt = ntt+(ntt/1);
    end
    
    

    

    % DELETE ORIGINAL FILAMENTS
    %----------------------------------------------------------------------
	if doDelOrig && (nT == doDelOrigT)
        Actin(1:NStFils,:) = []; 
	end


    
    % PROCESS STORE AND REMOVE DEL BRANCHES
    %----------------------------------------------------------------------
    
    % Filament Lifetime
	FLif = Actin(:,16) + 1;
	Actin(:,16) = FLif;
    
    % Mean Length
	Actin(:,18) = (Actin(:,18).*(FLif-1) + Actin(:,1))./FLif;
    
    %Longest Ever Length
	MaxL = Actin(:,1) > Actin(:,17);
	Actin(MaxL,17) = Actin(MaxL,1); 
    
    
    % Delete Filaments With No Actin
	delFi = (Actin(:,1)<1);
	Actin(delFi,15) = nT;    % Death Timestamp
        
    Actin(delFi,:) = [];
    
	Nfi = numel(Actin(:,1));    % Current #of Filaments ("NFact" used above loop)
    
    
    
    
    % COMPUTE XYZ TIP LOCATION
    %----------------------------------------------------------------------
    % MATH - branch XYZ tip coordinates
	Actin(:,4) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* cos(Actin(:,2)) + Actin(:,3);
	Actin(:,7) = Actin(:,1) .* nm_per_p .* sin(Actin(:,8)) .* sin(Actin(:,2)) + Actin(:,6);
	Actin(:,10) = Actin(:,1) .* nm_per_p .* cos(Actin(:,8)) + Actin(:,9);
    
    
    

    
    % SAVE TipMatrix
    %----------------------------------------------------------------------
	if SvTpMx && (nT >SaveTipsAfter) && (mod(nT,SaveTipsRate) == 0)

    	ActMx = TipMatrix(fh1,doTM,nT,Actin,dims,AcMx,SPYH,rota,azel);
            BTs{numel(BTs)+1} = ActMx;
            AFMx{numel(AFMx)+1} = Actin;

	end
    
    
    
    

    % Counters
    %----------------------------------------------------------------------
	if doActCounts
        
    	NumFils(nT)     = Nfi;                    % Number of branch filaments
    	NumFAct(nT)     = FActinN;                % Number FActins
    	NumGAct(nT)     = GActinN;                % Number GActins
    	NumCOFACT(nT)   = COFdepoly+ACTdepoly;    % Number of Depoly events
    	NumARPACT(nT)   = ACTpoly+ARPpoly;        % Number of Poly events
    	NumdelFi(nT)    = sum(delFi);
        
        muFlength(nT)   = mean(Actin(:,1));
        
    	ArpBRsum(nT)    = sum(ArpBR);
    	NumAct_uM(nT)   = Gactin_uM;
    	NumfKa(nT)      = fKa;
    	NumFnowFtot(nT) = Nfi / Actin(end,12);
    	NumFArp(nT)     = FArpN;
    	NumGArp(nT)     = GArpN;
    	Num_nmmf(nT)    = sum(nmmf);
    	NumFded(nT)     = numel(DiedACTs(:,1));

        NumGFActinN(nT) = GActinN+FActinN;      % Number Unbound Actins
    	NumTHYM(nT)     = THYM_N;               % Number Unbound Thymosin
    	NumTHYMACT(nT)  = THYM_ACT_N;           % Number Thymosin+GActins
        NumKaTHYM(nT)   = tKa;                  % Ka Thymosin
    	NumKdTHYM(nT)   = tKd;                  % Kd Thymosin
	end

    
    
    
    

    % LIVE PLOTS
    %----------------------------------------------------------------------
    if doLive3DActinPlot && ~mod(nT,LiveTipMod)

        Live3DActinPlot(fh1,hax1,hax2,azel,AxLim,rota,nT,...
              Actin,inPSD,SPYhZN,SPYhZS,SPYnXY,SPYhXY,doLiveHullPlot,...
                ph_spyTip, ph_psdTip, ph_fil);
    	rota = rota + rot0;

    end
    

    if doLiveLinePlot && ~mod(nT,LiveLinePlotMod)

        LiveLinePlot(fh1,hax3,hax4,hax5,hax6,hax7,nT,dt,Nsteps,...
            NumGAct,NumFAct,NumGArp,NumFArp,NumTHYM,NumTHYMACT,...
            NumfKa,ArpBRsum,Num_nmmf,NumFils, hax2, muFlength)

    end
    
    
    if ~mod(nT,LivePlotMod)

        delete(hax7.Children);
        spf0  = sprintf(' Step:  % 3.6g  Min:  % 2.6g  ',[nT nT*dt/60]);
        spf1  = sprintf(' Tot Actin:  % 9.6g  ',TOTAL_ACTIN);
        spf2  = sprintf(' GActin:     % 9.6g  ',GActinN);
        spf3  = sprintf(' FActin:     % 9.6g  ',FActinN);
        spf4  = sprintf(' Thymosin:   % 9.6g  ',THYM_N);
        spf5  = sprintf(' ThymAct:    % 9.6g  ',THYM_ACT_N);
        spf6  = sprintf(' GArpN:      % 9.6g  ',GArpN);
        spf7  = sprintf(' FArpN:      % 9.6g  ',FArpN);
        spf8  = sprintf(' Ka Actin:   % 9.5f  ',fKa);
        spf9  = sprintf(' Ka Thymo:   % 9.3f  ',tKa);
        spf10 = sprintf(' Kd Thymo:   % 9.3f  ',tKd);
        text(.1,(h7y*.92),spf1 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.84),spf2 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.76),spf3 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.68),spf4 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.60),spf5 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.52),spf6 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.44),spf7 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.36),spf8 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.28),spf9 ,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        text(.1,(h7y*.20),spf10,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
        
        text(.1,(h7y*.05),spf0,'FontSize',14,'FontName','FixedWidth','Parent',hax7);
    
        spfd1 = sprintf('  nT: %1.5g  ',nT);
        spfd2 = sprintf('  GActinN: %1.5g  ',GActinN);
        spfd3 = sprintf('  FActinN: %1.5g  ',FActinN);
        spfd4 = sprintf('  NFils: %1.4g  ',NFact);
        spfd5 = sprintf('  GArpN: %1.4g  ',GArpN);
        spfd6 = sprintf('  FArpN: %1.4g  ',FArpN);
        spfd7 = sprintf('  ThymN: %1.0f  ',THYM_N);
        spfd8 = sprintf('  ThymActN: %1.0f  ',THYM_ACT_N);
    
        disp([spfd1 spfd2 spfd3 spfd4 spfd5 spfd6 spfd7 spfd8]);
        pause(.001)
    end
    %--------------------------------------------------
    
    
    



end
%%                    	MAIN OUTER LOOP
%==========================================================================





% Lifetime of remaining filaments
Actin(:,15) = (nT + (nT-Actin(:,14)));      % Artificial Death Time
Actin(:,16) = (Actin(:,15) - Actin(:,14));  % Lifetime (nTdied - nTborn)







%%    	Counter Crunch and Plot
%==========================================================================
if doActCounts


    doNdel    = AMS{8};
    nDelSteps = AMS{9};

    if doNdel
            % Ndel = SaveTipsAfter; % Deletes the first 'Ndel' steps
            Ndel = nDelSteps;       % Deletes the first 'Ndel' steps
            nT0=nT; Ns0=Nsteps;     % Save the original Nt value
            nT = nT-Ndel;           % Reduce Nt by number of deleted steps
            Nsteps = nT;

            NumFils(1:Ndel) = [];
            NumFAct(1:Ndel) = [];
            NumGAct(1:Ndel) = [];
            NumCOFACT(1:Ndel) = [];
            NumARPACT(1:Ndel) = [];
            NumdelFi(1:Ndel) = [];

            ArpBRsum(1:Ndel) = [];
            NumAct_uM(1:Ndel) = [];
            NumfKa(1:Ndel) = [];
            NumFnowFtot(1:Ndel) = [];
            NumFArp(1:Ndel) = [];
            NumGArp(1:Ndel) = [];
            Num_nmmf(1:Ndel) = [];
            NumFded(1:Ndel) = [];

            NumGFActinN(1:Ndel) = [];
            NumTHYM(1:Ndel) = [];
            NumTHYMACT(1:Ndel) = [];
    end


    % PLOT OUTPUT DATA
    if doFinalPlots



    end

end




%--------------------------------------------------------------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);





AxLP = {fh1,nT,Actin,inPSD,rota,azel,dims};
Ax = {0,nT,Actin,dims,AcMx,SPYH,rota,azel};





%% AFMx
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% N  Xa  Xo  Xt  Ya  Yo  Yt  Za  Zo  Zt  MomID  ID  Fkd Born Died Lif MaxL MeanL Null Lgth
%--------------------------------------------------------------------------
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
% 1  2   3   4   5   6   7   8   9   10  11     12  13  14   15   16  17   18    19   20
%--------------------------------------------------------------------------
% DISPLAY WHAT IS CONTAINED IN EACH CELL OF AFMx
% The first row in this preview table is the column number
% those values aren't actually part of the Nx20 matrix

Tnams={'N';'Xa';'Xo';'Xt';'Ya';'Yo';'Yt';'Za';'Zo';'Zt';...
       'MomID';'ID';'Fkd';'Born';'Died';'Lif';'MaxL';'MeanL';'Null';'Lgth'};
l=[1:20; Actin(50:60,:)]; g=1:8;
T=table(l(g,1),l(g,2),l(g,3),l(g,4),l(g,5),l(g,6),l(g,7),l(g,8),l(g,9),l(g,10),l(g,11),l(g,12),...
l(g,13),l(g,14),l(g,15),l(g,16),l(g,17),l(g,18),l(g,19),l(g,20),'VariableNames',Tnams(1:20));
disp(' '); disp(T); disp(' ');




%% GATHER AND FORMAT FILAMENT DATA
%--------------------------------------------------------------------------

PSD_Zdim = dims(7);
SPYtips = {};
PSDtips = {};
ActXYZ = {};

for nn = 1:numel(AFMx)

    Ain = AFMx{nn};

    ActXYZ{nn} = [Ain(:,3) Ain(:,4) Ain(:,6) Ain(:,7) Ain(:,9) Ain(:,10)];
    ActinTips = [ActXYZ{nn}(:,2) ActXYZ{nn}(:,4) ActXYZ{nn}(:,6)];

    [Zrow,Zcol] = find(ActinTips(:,3) < PSD_Zdim);

    SPYtps = ActinTips(Zrow,:);
    SPYtips{nn} = SPYtps;

    [Zrow,Zcol] = find(ActinTips(:,3) >= PSD_Zdim);

    PSDtps = ActinTips(Zrow,:);
    PSDtips{nn} = PSDtps;

end








%% SAVE QUANTITATIVE DATA TO DISK




ActinData.AMX = AMX;                    % Original Input Parameters A
ActinData.AMS = AMS;                    % Original Input Parameters B
ActinData.Nsteps = Nsteps;              % Simulated Steps
ActinData.dt = dt;                      % Real time per step
ActinData.ACTs = ACTs;                  % Full Actin matrix
ActinData.Actin = Actin;                % Last-frame Actin matrix

ActinData.NumFils = NumFils;            % Number of branch filaments
ActinData.NumFAct = NumFAct;            % Number FActins
ActinData.NumGAct = NumGAct;            % Number GActins
ActinData.NumCOFACT = NumCOFACT;        % Number of Depoly events
ActinData.NumARPACT = NumARPACT;        % Number of Poly events
ActinData.NumdelFi = NumdelFi;

ActinData.ArpBRsum = ArpBRsum;          % Arp Branching rate
ActinData.NumAct_uM = NumAct_uM;        % Actin concentration uM
ActinData.NumfKa = NumfKa;              % fKa
ActinData.NumFnowFtot = NumFnowFtot;    % Filaments now
ActinData.NumFArp = NumFArp;            % Number of fArp
ActinData.NumGArp = NumGArp;            % Number of gArp
ActinData.Num_nmmf = Num_nmmf;          % nanometers of filament
ActinData.NumFded = NumFded;            % Number of filaments died
ActinData.NumGFActinN = NumGFActinN;    % Number Unbound Actins
ActinData.NumTHYM = NumTHYM;            % Number Unbound Thymosin
ActinData.NumTHYMACT = NumTHYMACT;      % Number Thymosin+GActins
ActinData.NumKaTHYM = NumKaTHYM;        % Ka Thymosin
ActinData.NumKdTHYM = NumKdTHYM;        % Kd Thymosin
ActinData.muFlength = muFlength;        % Mean Filament Length




%% SAVE MAT FILE TO DISK

% ActinMx.Ax = Ax;
% ActinMx.AFMx = AFMx;
ActinMx.ActXYZ = ActXYZ;
ActinMx.SPYtips = SPYtips;
ActinMx.PSDtips = PSDtips;


save('ActinMx.mat', 'ActinMx')
save('ActinData.mat', 'ActinData')









        


%% RETURN VALUES AS VARARGOUT

varargout = {BTs,AxLP,Ax,AFMx,ActinMx};

end






%##########################################################################
%% FUNCTIONS
%##########################################################################




% RotateVec()
%==========================================================================
function [Pn] = RotateVec(x,y,z,a,b,c,u,v,w,t)


Pn = ... 
	([(((a*(v^2+w^2)-u*(b*v+c*w-u*x-v*y-w*z))*(1-cos(t))+(u^2+v^2+w^2)*x*cos(t)+...
		sqrt((u^2+v^2+w^2))*(-c*v+b*w-w*y+v*z)*sin(t))/(u^2+v^2+w^2));
	
	(((b*(u^2+w^2)-v*(a*u+c*w-u*x-v*y-w*z))*(1-cos(t))+(u^2+v^2+w^2)*y*cos(t)+...
		sqrt((u^2+v^2+w^2))*(c*u-a*w+w*x-u*z)*sin(t))/(u^2+v^2+w^2));
	
	(((c*(u^2+v^2)-w*(a*u+b*v-u*x-v*y-w*z))*(1-cos(t))+(u^2+v^2+w^2)*z*cos(t)+...
		sqrt((u^2+v^2+w^2))*(-b*u+a*v-v*x+u*y)*sin(t))/(u^2+v^2+w^2))]);

end





% Live3DActinPlot()
%==========================================================================
function Live3DActinPlot(fh1,hax1,hax2,azel,AxLim,rota,nT,...
	Actin,inPSD,SPYhZN,SPYhZS,SPYnXY,SPYhXY,doLiveHullPlot,...
	ph_spyTip, ph_psdTip, ph_fil)



    ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
    [Zrow1,~] = find(ActinTips(:,3) > inPSD);
    PSDTips = ActinTips(Zrow1,:);
    [Zrow2,~] = find(ActinTips(:,3) < inPSD);
    SPYTips = ActinTips(Zrow2,:);

    P3x = [Actin(:,3) Actin(:,4)]';
    P3y = [Actin(:,6) Actin(:,7)]';
    P3z = [Actin(:,9) Actin(:,10)]';
    P3x(3,:)=NaN; P3y(3,:)=NaN; P3z(3,:)=NaN;

    

set(ph_spyTip,'XData',SPYTips(:,1)','YData',SPYTips(:,2)','ZData',SPYTips(:,3)');
set(ph_psdTip,'XData',PSDTips(:,1)','YData',PSDTips(:,2)','ZData',PSDTips(:,3)');
set(ph_fil,'XData',P3x(:),'YData',P3y(:),'ZData',P3z(:));
    view(azel+rota)



% SPINE MEMBRANE PLOT (USING: delaunayTriangulation & freeBoundary)
if doLiveHullPlot

    Ori_xyz = [Actin(:,3) Actin(:,6) Actin(:,9)];
    Tip_xyz = [Actin(:,4) Actin(:,7) Actin(:,10)];
    ATPxyz = [Ori_xyz; Tip_xyz];


    ATPxyz(end+1,:) = 0;
    szATP = size(ATPxyz,1);
    ATPxyz(szATP+1,:) = [SPYnXY/2 SPYnXY/2 0];
    ATPxyz(szATP+2,:) = [-SPYnXY/2 SPYnXY/2 0];
    ATPxyz(szATP+3,:) = [SPYnXY/2 -SPYnXY/2 0];
    ATPxyz(szATP+4,:) = [-SPYnXY/2 -SPYnXY/2 0];

    
    Uxyz = unique(ATPxyz,'rows');

    TinHead = Uxyz(:,3) >= (SPYhZS);
    TxyzHead = Uxyz(TinHead,:);

    TinNeck = Uxyz(:,3) < (SPYhZS-5);
    TxyzNeck = Uxyz(TinNeck,:);

    TinChin = (Uxyz(:,3) >= (SPYhZS-120)) & (Uxyz(:,3) <= (SPYhZS+100));
    TxyzChin = Uxyz(TinChin,:);

    TriSZh = size(TxyzHead,1);
    TriSZn = size(TxyzNeck,1);
    TriSZc = size(TxyzChin,1);
    TriSZall = TriSZn+TriSZc+TriSZh;


    delete(hax2.Children);
    %axes(hax2);


    if (TriSZall>=9) && (TriSZn<24) || (TriSZc<24) || (TriSZh<24)

        UTriTip = delaunayTriangulation(Uxyz);
        [Utri,Upoints] = freeBoundary(UTriTip);

      trimesh(Utri,Upoints(:,1),Upoints(:,2),Upoints(:,3), ...
            'FaceColor',[.2 1 .2],'FaceAlpha', 0.35,'Parent', hax2);
            hold on
    
    else

        TriHead = delaunayTriangulation(TxyzHead);
        [FBtHead,FBpHead] = freeBoundary(TriHead);

      trimesh(FBtHead,FBpHead(:,1),FBpHead(:,2),FBpHead(:,3), ...
            'FaceColor',[.2 1 .2],'FaceAlpha', 0.35,'Parent', hax2);
            hold on


        TriNeck = delaunayTriangulation(TxyzNeck);
        [FBtNeck,FBpNeck] = freeBoundary(TriNeck);

      trimesh(FBtNeck,FBpNeck(:,1),FBpNeck(:,2),FBpNeck(:,3), ...
            'FaceColor',[.2 1 .2],'FaceAlpha', 0.35,'Parent', hax2);
            hold on


        TriChin = delaunayTriangulation(TxyzChin);
        [FBtChin,FBpChin] = freeBoundary(TriChin);

      trimesh(FBtChin,FBpChin(:,1),FBpChin(:,2),FBpChin(:,3), ...
            'FaceColor',[.2 1 .2],'FaceAlpha', 0.3,'Parent', hax2);
            hold on

    end

    % scatter3(hax2,Uxyz(:,1),Uxyz(:,2),Uxyz(:,3),25,'fill');
    view(azel+rota)

end



pause(.001);
end






% LiveLinePlot()
%==========================================================================
function LiveLinePlot(fh1,hax3,hax4,hax5,hax6,hax7,nT,dt,Nsteps,...
            NumGAct,NumFAct,NumGArp,NumFArp,NumTHYM,NumTHYMACT,...
            NumfKa,ArpBRsum,Num_nmmf,NumFils, hax2, muFlength)


if nT == 20000; 
    
    hax2.XLim = [20 Nsteps*dt/60];
    hax3.XLim = [20 Nsteps*dt/60];
    hax4.XLim = [20 Nsteps*dt/60];
    hax5.XLim = [20 Nsteps*dt/60];
    hax6.XLim = [20 Nsteps*dt/60];
    % hax7.XLim = [10 Nsteps*dt/60];


end;

%%

tt = round(1/dt*60);
xax = 1:numel(NumFAct(1:tt:nT));
delete([hax3.Children,hax4.Children,hax5.Children]);
delete(hax6.Children); delete(hax7.Children); delete(hax2.Children);


% PLOT ACTIN LEVELS
line(xax,NumFAct(1:tt:nT),'Color',[1 0 0],'Parent',hax3,'LineWidth',3);
	hold on
line(xax,NumGAct(1:tt:nT),'Color',[0 0 1],'Parent',hax3,'LineWidth',3);


% PLOT ARP LEVELS
line(xax,NumFArp(1:tt:nT),'Color',[1 0 0],'Parent',hax4,'LineWidth',3);
	hold on
line(xax,NumGArp(1:tt:nT),'Color',[0 0 1],'Parent',hax4,'LineWidth',3);


% PLOT THYMOSIN LEVELS

line(xax,NumTHYM(1:tt:nT),'Color',[1 0 0],'Parent',hax5,'LineWidth',3);
	hold on
line(xax,NumTHYMACT(1:tt:nT),'Color',[0 0 1],'Parent',hax5,'LineWidth',3);


% PLOT FILAMENT NUMBER & LENGTH

line(xax,NumFils(1:tt:nT),'Color',[1 0 0],'Parent',hax6,'LineWidth',3);




% PLOT FILAMENT LENGTH

line(xax,muFlength(1:tt:nT),'Color',[1 0 0],'Parent',hax2,'LineWidth',3);



%%

% alternative way to live plot:
%{
    delete(hax3.Children);
    axes(hax3);

plot(NumFAct(1:tt:nT)','r');
	hold on
plot(NumGAct(1:tt:nT)','b');
    %set(hax3,'XLim', [0 Nsteps])

%}


pause(.001);
% if nT > 20000; keyboard; end
end





% MakeStartFils()
%==========================================================================
function Actin = MakeStartFils(Actin,NStFils,StartMonos,d2r,fZa,fXYo,SPYhZN,SPYhZS,SPYhXY)



% FIL-1
%--
Actin(1,1) = StartMonos;
Actin(1,2) = 1;                % Xa
Actin(1,5) = 1;                % Ya
Actin(1,8) = 91;            % Za
Actin(1,3) = 1;                % Xo
Actin(1,6) = 1;                % Yo
Actin(1,9) = 1;                % Zo
%--

if NStFils > 1
% FIL-2
%--
Actin(2,1) = StartMonos;
Actin(2,2) = -45;            % Xa
Actin(2,5) = 45;            % Ya
Actin(2,8) = d2r*-fZa;        % Za
Actin(2,3) = fXYo;            % Xo
Actin(2,6) = fXYo;            % Yo
Actin(2,9) = 1;                % Zo
%--
end

if NStFils > 2
% FIL-3
%--
Actin(3,1) = StartMonos;
Actin(3,2) = 45;            % Xa
Actin(3,5) = -45;            % Ya
Actin(3,8) = d2r*fZa;        % Za
Actin(3,3) = fXYo;            % Xo
Actin(3,6) = -fXYo;            % Yo
Actin(3,9) = 1;                % Zo
%--
end

if NStFils > 3
% FIL-4
%--
Actin(4,1) = StartMonos;
Actin(4,2) = 45;            % Xa
Actin(4,5) = 0;                % Ya
Actin(4,8) = d2r*-fZa;        % Za
Actin(4,3) = -fXYo;            % Xo
Actin(4,6) = fXYo;            % Yo
Actin(4,9) = 1;                % Zo
%--
end

if NStFils > 4
% FIL-5
%--
Actin(5,1) = StartMonos;
Actin(5,2) = 0;                % Xa
Actin(5,5) = 45;            % Ya
Actin(5,8) = d2r*fZa;        % Za
Actin(5,3) = -fXYo;            % Xo
Actin(5,6) = -fXYo;            % Yo
Actin(5,9) = 1;                % Zo
%--
end


if NStFils > 5
% FIL-6
%--
Actin(6,1)=round(SPYhZN-SPYhZS-2);
Actin(6,2) = 225*d2r;        % Xa
Actin(6,5) = 0*d2r;            % Ya
Actin(6,8) = 45*d2r;        % Za
Actin(6,3) = (SPYhXY/2)-10;    % Xo
Actin(6,6) = (SPYhXY/2)-10;    % Yo
Actin(6,9) = SPYhZS+1;        % Zo
end

if NStFils > 6
% FIL-7
%--
Actin(7,1)=round(SPYhZN-SPYhZS-2);
Actin(7,2) = 45*d2r;        % Xa
Actin(7,5) = 0*d2r;            % Ya
Actin(7,8) = 45*d2r;        % Za
Actin(7,3) = 10-(SPYhXY/2);    % Xo
Actin(7,6) = 10-(SPYhXY/2);    % Yo
Actin(7,9) = SPYhZS+1;        % Zo
end

if NStFils > 7
% % FIL-8
% %--
Actin(8,1)=round(SPYhZN-SPYhZS-2);
Actin(8,2) = 315*d2r;        % Xa
Actin(8,5) = 0*d2r;            % Ya
Actin(8,8) = 45*d2r;        % Za
Actin(8,3) = 10-(SPYhXY/2);    % Xo
Actin(8,6) = (SPYhXY/2)-10;    % Yo
Actin(8,9) = SPYhZS+1;        % Zo
end

if NStFils > 8
% FIL-9
%--
Actin(9,1)=round(SPYhZN-SPYhZS-2);
Actin(9,2) = 135*d2r;        % Xa
Actin(9,5) = 0*d2r;            % Ya
Actin(9,8) = 45*d2r;        % Za
Actin(9,3) = (SPYhXY/2)-10;    % Xo
Actin(9,6) = 10-(SPYhXY/2);    % Yo
Actin(9,9) = SPYhZS+1;        % Zo
end


end





% TipMatrix()
%==========================================================================
function [varargout] = TipMatrix(Fh,doTM,nT,Actin,dims,AcMx,SPYH,varargin)


%{
% keyboard

scsz = get(0,'ScreenSize');
basepos = scsz./[2.5e-3 2.5e-3 1.5 2];
baseazel = [-32 12];
baserot=[0 0];

if nargin >= 2
	rot=varargin{1};
	azel=varargin{2};
elseif nargin == 1 
	rot=varargin{1};
	azel = baseazel;
else
	rot=baserot;
	azel = baseazel;
end






%--------------------
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];
[Zrow1,Zcol1] = find(ActinTips(:,3) > inPSD);
PSDTips = ActinTips(Zrow1,:);
[Zrow2,Zcol2] = find(ActinTips(:,3) < inPSD);
SPYTips = ActinTips(Zrow2,:);
%--------------------
figure(Fh)
subplot('Position',[.08 .15 .45 .70]), 
ph11c = plot3([Actin(:,3) Actin(:,4)]', [Actin(:,6) Actin(:,7)]', [Actin(:,9) Actin(:,10)]');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel)
grid off
set(gca,'Color',[1,1,1])
hold on;
ph11a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph11b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
set(gcf,'Color',[1,1,1])
xlabel('X');ylabel('Y');zlabel('Z');
view(azel+rot)
grid off
set(gca,'Color',[1,1,1])
set(ph11a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph11b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
set(ph11c,'LineStyle','-','Color',[.7 .7 .7],'LineWidth',.1);
hold off;
%--------------------
figure(Fh)
subplot('Position',[.6 .55 .38 .38]), 
ph12a = scatter3([SPYTips(:,1)]', [SPYTips(:,2)]', [SPYTips(:,3)]',7,'ob');
hold on;
ph12b = scatter3([PSDTips(:,1)]', [PSDTips(:,2)]', [PSDTips(:,3)]',7,'or');
axis([-500 500 -500 500 0 1000])
view([0 90])
grid off
set(gca,'Color',[1,1,1])
set(ph12a,'Marker','o','MarkerEdgeColor',[.1 .1 .9],'MarkerFaceColor',[.1 .1 .9]);
set(ph12b,'Marker','o','MarkerEdgeColor',[.9 .2 .2],'MarkerFaceColor',[.9 .2 .2]);
%--------------------
set(gca,'XTickLabel', sprintf('%.1f|',nT),'FontSize',10)
hold off;
%--------------------
%}


%==================================================%
% Spine Dimensions
HX = dims(4);
HY = dims(5);
SPYz = dims(7);
SPYxy = dims(8);
%-----------------------------------%
ActinTips = [Actin(:,4) Actin(:,7) Actin(:,10)];


[Zin,~] = find(ActinTips(:,3) > SPYz);
XYsq = sqrt(ActinTips(:,1).^2 + ActinTips(:,2).^2);
[XYin,~] = find(XYsq < (HX-SPYxy));
XYZin = intersect(XYin, Zin);
% XYt = intersect(Xrow1, Yrow1);
%[tf, loc] = ismember(Xrow1, Yrow1)


PSDTips = ActinTips(XYZin,:);    % get just the tips, just to see how they feel
PSDXY = round([PSDTips(:,1) PSDTips(:,2)]);
PSDY = PSDXY(:,2)+HY;    % Add Y radius width
PSDX = PSDXY(:,1)+HX;    % Add X radius width
SPYHX = SPYH(1)*2;        % I'm specifying a radius in SPYH(), must double
SPYHY = SPYH(2)*2;        % I'm specifying a radius in SPYH(), must double
PSDX(PSDX>SPYHX) = SPYHX-1;    % Everything here that is somehow greater SPYHX place at edge of SPYHX
PSDY(PSDY>SPYHY) = SPYHY-1;    % Everything here that is somehow greater SPYHY place at edge of SPYHY
PSDX(PSDX<1) = 1;    % Everything here that is somehow negative place at 1
PSDY(PSDY<1) = 1;    % Everything here that is somehow negative place at 1
%==================================================%

PSDX = round(PSDX/5);
PSDY = round(PSDY/5);

for mxp = 1:numel(PSDXY(:,1))
AcMx(PSDY(mxp), PSDX(mxp)) = 1;
end

% SAVE THIS (COULD BE USEFUL)
% ActMask=ones(4);
% ActMx = convn(AcMx,ActMask,'same');
% ActMx = (ActMx>0).*1.0;

ActMx = AcMx;
%===================================%
%            	FIGURE
%-----------------------------------%
if doTM
    clrmap = [1 1 1; .9 .2 .2];
    figure(Fh)
    subplot('Position',[.7 .1 .28 .38]), 
    imagesc(ActMx)
    colormap(clrmap);
    hold on
    subplot('Position',[.7 .1 .28 .38]), 
    scatter(PSDX,PSDY, 'r')
    hold off
end
%-----------------------------------%


%==================================================%
% if nT >20000; keyboard; end
varargout = {ActMx};
end




