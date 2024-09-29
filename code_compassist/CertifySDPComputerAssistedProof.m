%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: QETLAB and all .mat files like data_kslots_k_Aplusk_B_general
%Last update: 28/09/2024

%Function which constructs a computer assisted proof to certify that the quantum switch cannot be simulated by combs of a given scenario
%Below in this code you'll find some auxiliary function

%Input: scenario, a string with the scenario to be analysed (Available options : AA, AB, AAA, ABA, AAB, BAA, ABAB, AABB, ABBA, AAAB, AABA, ABAA, BAAA)
%Input (Optional): DecimalPrecision, set the number of decimal cases which will be used. The default value is 9
%Input (Optional): SaveOutputVariables, set 1 to save main variables, 2 to save all variables, and 0 to save nothing. The default value is 0
%Input (Optional): SymbolicCholesky, set 1 to force all Cholesky decompositions to be done symbolically (instead of arbitrary-precision arithmetic)
%Output: [p, pOK, pHP, isEverythingOK], % p is the value which comes from the SDP
% pOK is a symbolic variable which is certified by a computer-assisted proof
% pHP is the High-Precision value, where the SDP is certified using floating-point arithmetic
% isEverythingOK is a Boolean variable which states if everything worked correctly

function [p, pOK, pHP, isEverythingOK] = CertifySDPComputerAssistedProof(scenario,DecimalPrecision,SaveOutputVariables,SymbolicCholesky)
TotalTimeThisFunction = tic;
scenario = num2str(scenario); % This is just to prevent bugs resulting from the difference between "AB" and 'AB'
% First we check if the variables DecimalPrecision and SaveOutputVariables were initialised, if not set the default values
if exist('DecimalPrecision')
else
    DecimalPrecision = 9;
end
if exist('SaveOutputVariables')
else
    SaveOutputVariables = 0;
end
if exist('SymbolicCholesky')
else
    SymbolicCholesky=0;
end
DecimalPrecision = DecimalPrecision
SaveOutputVariables = SaveOutputVariables
SymbolicCholesky = SymbolicCholesky
disp('We start by running the High Precision code, because it is fast and interesting')
[pHP] = CertifySDPHighPrecision(scenario,DecimalPrecision);
'Great, the High-Precision code has finished!! We now run the computer assisted code';
%%%%%%%%%%%% Conventions for Li, L, and t
% In this code, we don't use Rij and J^A_i, and J^B_j. We group (i,j) as a single variable i, that is, i = (i,j)
% Also, the channels are given by A_i and B_i
% Li is defined as Li := S*(A_i ⊗ B_i)
% L is defined as L := \sum_i transpose(A_i{^⊗k1} ⊗ B_i{^⊗k2}) ⊗ R_i
% t is defined as t := \sum_i \tr(R_i L_i)
L=0;
switch scenario
    case 'AB'
        load('inputdata_2slots_1plus1_general') %Load variables
        Ri = Ri_AB;  G = G_AB;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),B(:,:,i))),Ri(:,:,i));
        end
    case 'AA'
        load('inputdata_2slots_2plus0_general') %Load variables
        Ri = Ri_AA;  G = G_AA;    	A = C; B = C;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    case 'AAB'
        load('inputdata_3slots_2plus1_general') %Load variables
        Ri = Ri_AAB;  G = G_AAB;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),A(:,:,i),B(:,:,i))),Ri(:,:,i));
        end
    case 'ABA'
        load('inputdata_3slots_2plus1_general') %Load variables
        Ri = Ri_ABA;  G = G_ABA;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),B(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    case 'BAA'
        load('inputdata_3slots_2plus1_general') %Load variables
        Ri = Ri_BAA;  G = G_BAA;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),B(:,:,i),B(:,:,i))),Ri(:,:,i));
        end
    case 'AAA'
        load('inputdata_3slots_3plus0_general') %Load variables
        Ri = Ri_AAA;  G = G_AAA;    	A = C; B = C;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),A(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    case 'AABB'
        load('inputdata_4slots_2plus2_general') %Load variables
        Ri = Ri_AABB;  G = G_AABB;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),A(:,:,i),B(:,:,i),B(:,:,i))),Ri(:,:,i));
        end
    case 'ABAB'
        load('inputdata_4slots_2plus2_general') %Load variables
        Ri = Ri_ABAB;  G = G_ABAB;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),B(:,:,i),A(:,:,i),B(:,:,i))),Ri(:,:,i));
        end
    case 'ABBA'
        load('inputdata_4slots_2plus2_general') %Load variables
        Ri = Ri_ABBA;  G = G_ABBA;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),B(:,:,i),B(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    case 'AAAB'
        load('inputdata_4slots_3plus1_general') %Load variables
        Ri = Ri_AAAB;  G = G_AAAB;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),A(:,:,i),A(:,:,i),B(:,:,i))),Ri(:,:,i));
        end
    case 'AABA'
        load('inputdata_4slots_3plus1_general') %Load variables
        Ri = Ri_AABA;  G = G_AABA;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),A(:,:,i),B(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    case 'ABAA'
        load('inputdata_4slots_3plus1_general') %Load variables
        Ri = Ri_ABAA;  G = G_ABAA;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(A(:,:,i),B(:,:,i),A(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    case 'BAAA'
        load('inputdata_4slots_3plus1_general') %Load variables
        Ri = Ri_BAAA;  G = G_BAAA;    	A = C1; B = C2;
        for i=1:size(A,3)
            L = L + kron(transpose(Tensor(B(:,:,i),A(:,:,i),A(:,:,i),A(:,:,i))),Ri(:,:,i));
        end
    otherwise
        error('You chose a scenario which is not covered by this code')
end
%%%%%%%%%%%% Create the reduced qubit Switch
if exist('SWITCH')
    SWITCH_loaded_from_input = SWITCH;
else
    'For some reason, the SWITCH was not defined in the loaded input variables.'
    'No worries, we will defined it!!'
end
d=2; %Set dimension as d=2
alpha = 1/sqrt(2); beta  = sqrt(1-alpha^2); % control in qubit |+>
target = zeros(d,1); target(1,1) = 1; % target in |0> in dimension d
ketId = [1 0 0 1]'; % Define |Id>> = |00>+|11>
% Choi of identity channel in dimension d
SWITCH = alpha*Tensor(target,ketId,ketId,[1;0]) + beta*PermuteSystems(Tensor(target,ketId,ketId,[0;1]),[3 4 1 2 5 6],[d d d d d 2]);
%SWITCH in Ai Ao Bi Bo Ft Fc,  %dim(Ai)=dim(Ao)=dim(Bi)=dim(Bo)=dim(Ft)=d, dim(Fc)=2
SWITCH = PartialTrace(SWITCH*SWITCH',2,[d^4 d 2]); %This is the reduced switch
if exist('SWITCH_loaded_from_input')
    if norm(SWITCH_loaded_from_input-SWITCH)==0
        disp("The Switch constructed in this function coincides with the Switch from input_variables" + newline + "All good regarding the SWITCH!  =)");
    else
        error('The Switch constructed in this function is different from the Switch from input_variables, there is a bug somewhere...')
    end
end
k = size(scenario,2); %evaluate the number of calls k
p = trace(G)/2/d^k; % This is the p which comes from the SDP output
pInitial = p
%%%%%%%%%%%% Create Li := S*(A_i ⊗ B_i)
Li=zeros(2,2,size(A,3));
for i=1:size(A,3)
    Li(:,:,i) = PartialTrace(SWITCH*transpose(Tensor(A(:,:,i),B(:,:,i),eye(2))),1,[(d*d)^2 2]);
end
%%%%%%%%%%%% MAKE EVERYTHING SYMBOLIC  (GS is Gamma Symbolic)
Id2sym = sym(eye(2));
SWITCHS=sym(chopMTQdetail(SWITCH,9)); %Switch Symbolic
AS=sym(chopMTQdetail(A,DecimalPrecision)); %Channel A Symbolic
for i=1:size(AS,3)
    if Id2sym==PartialTrace2SysMTQ(AS(:,:,i),2,[2 2])
    else
        error('ONE OF THE CHANNELS AS IS NOT PROPERLY NORMALISED!!!')
    end
end
disp('Great, all channels AS are properly normalised')
BS=sym(chopMTQdetail(B,DecimalPrecision)); %Channel B Symbolic
for i=1:size(AS,3)
    if Id2sym==PartialTrace2SysMTQ(BS(:,:,i),2,[2 2])
    else
        error('ONE OF THE CHANNELS BS IS NOT PROPERLY NORMALISED!!!')
    end
end
disp('Great, all channels BS are properly normalised')
GS = sym(chopMTQdetail(G,DecimalPrecision)); %Make G symbolic
GS = (GS + GS')/2; %Make GS Self-Adjoint
tic
'We now create LiS, the symbolic version of Li := S*(A_i{^⊗k1} ⊗ B_i{^⊗k2})'
LiS = sym(zeros(size(Li)));
for i=1:size(Li,3)
    LiS(:,:,i) = PartialTrace2SysMTQ(SWITCHS*transpose(Tensor(AS(:,:,i),BS(:,:,i),Id2sym)),1,[(d*d)^2 2]);
end
timeLiS_min = toc/60
'We now project GS into a valid subspace'
DIM = [d*ones(1,2*k) 2];
GSValid = ProjSeqSuperChannelNoFuture4Sym(GS,DIM);
clear GS;  %Here, we clear the variable GS, to save ram memory
LS=sym(0); %Initialise LS as a symbolic variable
tic
'We now create RiS and LS. L is defined as L := \sum_i transpose(A_i{^⊗k1} ⊗ B_i{^⊗k2}) ⊗ R_i'
!date
RiS = sym(chopMTQdetail(Ri,DecimalPrecision));  %Initialise RiS as the symbolic truncated version of Ri
switch scenario
    case 'AB'
        for i=1:size(A,3)
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),BS(:,:,i))),RiS(:,:,i));
        end
    case 'AA'
        for i=1:size(A,3)
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'AAB'
        for i=1:size(A,3)
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),AS(:,:,i),BS(:,:,i))),RiS(:,:,i));
        end
    case 'ABA'
        for i=1:size(A,3)
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),BS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'BAA'
        for i=1:size(A,3)
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(BS(:,:,i),AS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'AAA'
        for i=1:size(A,3)
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),AS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'ABAB'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),BS(:,:,i),AS(:,:,i),BS(:,:,i))),RiS(:,:,i));
        end
    case 'ABBA'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),BS(:,:,i),BS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'AABB'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),AS(:,:,i),BS(:,:,i),BS(:,:,i))),RiS(:,:,i));
        end
    case 'AAAB'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),AS(:,:,i),AS(:,:,i),BS(:,:,i))),RiS(:,:,i));
        end
    case 'AABA'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),AS(:,:,i),BS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'ABAA'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(AS(:,:,i),BS(:,:,i),AS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    case 'BAAA'
        tic
        for i=1:size(A,3)
            if i==1
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            end
            if mod(i,1000)==0
                disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
                time_for_1000_iteractionsHours = toc/60/60
                tic
            end
            RiS(:,:,i) = (RiS(:,:,i) + RiS(:,:,i)')/2;
            LS = LS + kron(transpose(Tensor(BS(:,:,i),AS(:,:,i),AS(:,:,i),AS(:,:,i))),RiS(:,:,i));
        end
    otherwise
        error('You chose a scenario which is not covered by this code')
end
time4Ris_and_LS_min = toc/60
time4Ris_and_LS_hours = toc/60/60
'We now create tS, LSOK, and RiSOK. t is defined as \sum_i \tr(R_i L_i)'
tS=sym(0);
for i=1:size(A,3)
    tS=tS+HS(RiS(:,:,i),LiS(:,:,i));  %HS(A,B) = trace(A'B)
end
if abs(double(tS)-1)>10^(-4)
    tSvalue = vpa(tS)
    error("The value of tS is fair away from one.... this is strange... please check your code")
else
    tSvalue = vpa(tS)
    'The value t := \sum_i \tr(R_i L_i) is close to one, this is a good sign! '
    'The initial outcome of your SDP is close to respect the equality constraints =)'
end
RiSOK = sym(zeros(size(Ri))); %RiSOK will lead to a valid Ri such that t := \sum_i \tr(R_i L_i) = 1
for i=1:size(A,3)
    RiSOK(:,:,i) = RiS(:,:,i)/tS; %We just define all RiS by tS, this will ensure that  t := \sum_i \tr(R_i L_i) = 1
end
clear RiS %Clear RiS to save RAM
LSOK = LS/tS;  %Make LS valid %L is defined as L := \sum_i transpose(A_i{^⊗k1} ⊗ B_i{^⊗k2}) ⊗ R_i
clear LS %Clear LS to save RAM
%%%%%%%%%%%% Double check that tsOK is indeed one
isEverythingOK = 1; %This variable will be 1 if all tests pass, and 0 if there is a problem
tic
tSOK=sym(0);
for i=1:size(A,3)
    tSOK=tSOK+HS(RiSOK(:,:,i),LiS(:,:,i)); %HS(A,B) = trace(A'B), Li := S*(A_i ⊗ B_i)
end
if tSOK
    disp('Great, tSOK:=\sum_i \tr(RiSOK LiS) is equals to one!!')
else
    isEverythingOK = 0;
    tSOK = tSOK
    error(['Failure... tSOK:=\sum_i \tr(RiSOK LiS) is NOT equals to one!!', newline, 'Its value is something like ', num2str(double(tSOK))]);
end
timetSOK_min = toc/60

%%%%%%%%%%%% Chose eta and make the final GSOK
minEigGSValidDouble = min(eig(double(GSValid)))
minEigGSVALID_minus_LS_Double = min(eig(double(GSValid) - double(LSOK)))

lambda = abs(min(minEigGSValidDouble,minEigGSVALID_minus_LS_Double));
% etaDouble=ceil(lambda*10^(DecimalPrecision+1))/10^(DecimalPrecision+1)
etaDouble=(ceil(lambda*10^(DecimalPrecision))/10^(DecimalPrecision) + 10^(-7)) %Here we add 10(-7) to help us to ^ensure positivity 
eta = sym(etaDouble);
% eta = sym(0)
GSOK = GSValid + eta*sym(eye(size(G)));
clear GSValid %Clear GSValid to save RAM, we now keep GSOK
minEigGSOKDouble = min(eig(double(GSOK)))
minEigGSOK_minus_LS_Double = min(eig(double(GSOK) - double(LSOK)))
pOK = trace(GSOK)/2/d^k
pOKVPAprecision7 = vpa(pOK,7)

tic
if SymbolicCholesky
    Gprime = GSOK;
    disp(["We will now certify that GSOK is PSD via Cholesky decomposition of GSOK."]);
!date
    [IsPSDOK_GSOK,GSSqrt] = IsPSDMTQ(GSOK);  %Perform the Cholesky decomposition with symbolic variables. Fast for k<4
else
G_VPA = vpa(GSOK); %by default, vpa has d=32 decimal digits of precision
Gprime =  (G_VPA+G_VPA')/2 - vpa(10^(vpa(-29))*ones(size(G_VPA)));
% [IsPSDOK_GSOK,GSSqrt] = IsPSDMTQ(GSOK);  %With this line, we perform the Cholesky decomposition with symbolic variables
% [IsPSDOK_GSOK,GprimeSqrt] = IsPSDMTQ(Gprime); %With this line, we perform the Cholesky decomposition with vpa 
disp(["We will now certify that GSOK is PSD via Cholesky decomposition of ((G_VPA+G_VPA')/2 - vpa(10^(vpa(-29))*ones(size(G_VPA)))."]);
!date
[IsPSDOK_GSOK] = IsPD_Rump(Gprime); %With this line, we perform the Rump Cholesky decomposition with vpa
end

if IsPSDOK_GSOK
    disp(['Great, GSOK (gamma) is a symbolic PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGSOKDouble)]);
else
    disp(['Failure... we could not certify that GSOK is a symbolic PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGSOKDouble)])
    isEverythingOK=0;
end
timeCheckGSOK_PSD_min = toc/60

M = GSOK-LSOK;
if SymbolicCholesky
    Mprime = M;
    !date
    disp(["We will now certify that GSOK-LSOK is PSD via Cholesky decomposition of M=GSOK-LSOK"]);
    [IsPSDOK_GSOKminusLS,MSqrt] = IsPSDMTQ(M); %Perform the Cholesky decomposition with symbolic variables. Feasible for k<4 (for k=3 this step may take 1 hour)
else
Mvpa = vpa(M); %by default, vpa has d=32 decimal digits of precision
Mprime = (Mvpa+Mvpa')/2 - vpa(10^(vpa(-29))*ones(size(Mvpa)));
disp(["We will now certify that GSOK-LSOK is PSD via Cholesky decomposition of (Mvpa+Mvpa')/2 - vpa(10^(vpa(-29))*ones(size(Mvpa)))."]);
!date
tic
[IsPSDOK_GSOKminusLS] = IsPD_Rump(Mprime); %With this line, we perform the Rump Cholesky decomposition with vpa
end
if IsPSDOK_GSOKminusLS
    disp(['Great, M = GSOK-LSOK is a symbolic PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGSOK_minus_LS_Double)]);
else
    isEverythingOK=0;
    disp(['Failure... we could not certify that M = GSOK-LSOK is a symbolic PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGSOK_minus_LS_Double)]);
end
timeCheckGSOKminusLS_PSD_min = toc/60

if SaveOutputVariables
    disp(['We will now save important varaibles, this may take a very long time. (The current scenario is: ',scenario,')'])
    !date
    % We will save the date in a format like outputdata_3slots_2plus1_generalABA
    tic
    kA=count(scenario, 'A');    	kB=count(scenario, 'B');
    kMax = max(kA,kB);          	kMin=k-kMax;
    saveFileName = ["outputdata_"+num2str(k)+"slots_"+num2str(kMax)+"plus"+num2str(kMin)+"_general"+scenario+"_DecimalPrecision_"+num2str(DecimalPrecision)];
%     save(saveFileName,'GSOK','RiSOK','p','pOK','pHP','AS','BS','GprimeSqrt','MprimeSqrt')
save(saveFileName,'GSOK','RiSOK','LSOK','p','pOK','pHP','AS','BS','Gprime','Mprime')
    Time2saveImportantVariables_min = toc/60
end

if isEverythingOK
    disp(['Everything is OK!!!  We could perform a computer-assisted proof, this is great!!!!', newline, '=)', newline, 'a HIGH-PRECISION pOK for the scenario ',scenario,' is something like ', num2str(pHP)]);
    disp(['pOK (computer-assisted proof) for the scenario ',scenario,' is something like ', num2str(double(pOK)), newline, 'You may want to compare it with the initial p, which is ', num2str(p)]);
else
    error('DANGER, Gamma is not positive or Gamma - LSA is not positive...')
end

pOKVPAprecision7 = vpa(pOK,7)

if SaveOutputVariables>1
    tic
    disp('We will now ALL used varaibles, this may take a long time')
    !date
    save("SwitchSimulationFullData_"+scenario+"_DecimalPrecision_"+num2str(DecimalPrecision))
    disp('The FullOutputData data was saved! You have all info!!')
    Time2saveALLVariables_min = toc/60
end

TotalTimeComputerAssistedCodeMinutes = toc(TotalTimeThisFunction)/60
TotalTimeComputerAssistedCodeHours = toc(TotalTimeThisFunction)/60/60
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%  AUXILIARY FUNCTIONS %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out=chopMTQdetail(in,DecimalDigits)
%Function that truncates the input in with some particular number of DecimalDigits
%Requires: Nothing
%Input: in, a matrix or number
%Input: DecimalDigits, the number of DecimalDigits in the truncation
%Output: out, the truncated input in

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2024
out=round(in*10^DecimalDigits)/10^DecimalDigits;
end

function [out,R] = IsPSDMTQ(M)
%Function that verifies if a symbolic matrix is positive semidefinite using Cholesky decomposition
%Requires: Nothing
%Input: a complex symbolic matrix M
%Output: out=1 if PSD, out=0 if not PSD, and a Matrix R such that R'*R==M
%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2024
try R=chol(M);
    %     	if R'*R==M
    %         	out=1;
    %     	else
    %         	out=0;
    %     	end
    out=1;
catch
    out=0;
    R=0;
end
end

function C=ProjSeqSuperChannelNoFuture4Sym(X, DIM)
%Function that projects a matrix into the linear space spanned by Choi operators of sequential superchannels
%Requires: PartialTraceSys2MTQ
%Input: X, input matrix X\in L(H_in_1 ⊗ H_out_1 ⊗ H_in_2 ⊗ H_out_2 ⊗ ... ⊗  H_in_N ⊗ H_out_N)
%Input: DIM, vector with all the local dimensions ordered as [input1 output1 input2 output2 ... inputN outputN]
%Output: C, matrix X after being projected onto the linear space spanned by Choi operators of sequential quantum channels without future

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2024

Nspaces=max(size(DIM));

C=0;
for i=Nspaces:-1:2
    d1 = prod(DIM(1:i-1));
    d2 = prod(DIM(i:Nspaces));
    if mod(i,2)==0
        C = C - kron(PartialTrace2SysMTQ(X,2,[d1 d2]),eye(d2))/d2;
    else
        C = C + kron(PartialTrace2SysMTQ(X,2,[d1 d2]),eye(d2))/d2;
    end
end

C = C + trace(X)*eye(prod(DIM))/prod(DIM);
end

function [out] = PartialTrace2SysMTQ(M,sys,DIM)
%Function that performs the partial trace on system 2 of a matrix
%Requires: Nothing
%Input: a matrix M, its dimensions DIM = [d1 d2]
%Output: out=1 if PSD, out=0 if not PSD, and a Matrix R such that R'*R==M

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2024
d1 = DIM(1);
d2 = DIM(2);
out=0;
if sys == 2
    Id1 = eye(d1);
    for i=1:d2
        keti=zeros(d2,1);
        keti(i)=1;
        Ki=kron(Id1,keti');
        out = out + Ki*M*Ki';
    end
elseif sys==1
    Id2 = eye(d2);
    for i=1:d1
        keti=zeros(d1,1);
        keti(i)=1;
        Ki=kron(keti',Id2);
        out = out + Ki*M*Ki';
    end
else
    error('This partial function trace only works for bipartite systems')
end
end

function out=HS(A,B)
%This function implements the Hilbert Schmidt inner product of two matrices living in the same matrix space
%Requires: Nothing
%Input:	A,B, matrices which can be multiplied
%Output: out=<<A|B>>=trace(A'*B), %This is useful to reduce complexity and reduce numerical imprecision

%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Last update: 21/08/2024
out=A(:)'*B(:);
end

function [res,R] = IsPD_Rump(A)  %We also output R, which respects R'*R = A~
% This function is extracted from the paper "Verification of Positive Definiteness" from Siegfried M. Rump.
% https://doi.org/10.1007/s10543-006-0056-1
% I have commented out the line "p = symamd(A); A = A(p,p);", so that the function also works for vpa (variable precision arithmetic)

%ISSPD logical function: Matrix A is positive definite
%
%Given real symmetric or Hermitian complex matrix A,
%
% res 1 Matrix A is proved to positive definite
% 0 positive definiteness could not be verified
%
% constants
n = size(A,1); Eps = 2^(-53); Eta = 2^(-1074);
% diagonal check
if any( diag(A)<=0 )
    res = 0; return
end
% scaling

d = 2.^(-ceil(0.5*log2(diag(A))));
maxdiagA = max(d); mindiagA = min(d);
if ( maxdiagA/mindiagA>sqrt(n) ) && ~( ( maxdiagA>1e100 ) ||    ( mindiagA<1e-100 ) )
    % apply van der Sluis scaling
    D = spdiags( d ,0,n,n );
    % D_ii are powers of 2
    A = D*A*D;
    % 0.25 <= abs(A_ii) < 1
    maxdiagA = 1;
end

% Minimum degree sorting
%p = symamd(A); A = A(p,p);  %This line is commented out to ensure this code works for VPA as well
[i,j] = find(A); index = find(diff(j));
t = [0 ; (2:n)'-i(index+1)];
% max #elts left of (or above) diag(A)
if any ( t>3e15 )
    % dimension check, make sure alpha<1
    res = 0; return
end

% hull of A
if n<67108861
    % alpha_k/(1-alpha_k) < (k+1)*Eps
    alpha = (t+3)*Eps;
    % exact
else
    alpha = (t+2)*Eps; alpha = ( alpha./(1-alpha) )/(1-2*Eps);
    alpha = ( alpha./(1-alpha) )/(1-3*Eps);
end
d = sqrt(alpha.*diag(A))/(1-3*Eps);

% Upper bound for norm(dA) and shift
c = ( (3*n)*(2*n+maxdiagA)/(1-3*Eps) ) * Eta;
c = ( (d'*d)/(1-(n+2)*Eps) + c )/(1-2*Eps);
% bound II)
% floating-point Cholesky
A = A - c*speye(n); A = A-diag(diag(A)*(Eps*(1+2*Eps)));
[R,p] = chol(A);
% floating-point Cholesky
res = ( p==0 );
% p=0 <=> successful completion
end
