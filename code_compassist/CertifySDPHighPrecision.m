%Author: Marco Túlio Quintino, https://github.com/mtcq, mtcq.mm@gmail.com
%Requires: QETLAB and all .mat files like data_kslots_k_Aplusk_B_general 
%Last update: 21/08/2024

%Function which constructs a computer assisted proof to certify that the quantum switch cannot be simulated by combs of a given scenario
%Below in this code you'll find some auxiliary function

%Input: scenario, a string with the scenario to be analysed (Available options : AA, AB, AAA, ABA, AAB, BAA, ABAB, AABB,ABBA, AAAB, AABA, ABAA, BAAA)
%Input (Optional): DecimalPrecision, set the number of decimal cases which will be used. The default value is 9
%Input (Optional): SaveOutputVariables, set 1 to save main variables, 2 to save all variables, and 0 to save nothing. The default value is 0
%Output: [p, pHP, isEverythingOK], % p is the value which comes from the SDP
% pHP is the High-Precision value, where the SDP is certified using floating-point arithmetic
% isEverythingOK is a Boolean variable which states if everything worked correctly

function [pHP, p, isEverythingOK] = CertifySDPHighPrecision(scenario,DecimalPrecision,SaveOutputVariables)
TotalTimeThisFunction = tic;
scenario = num2str(scenario); % This is just to prevent bugs resulting from the difference between "AB" and 'AB'
if exist('DecimalPrecision')
else
    DecimalPrecision = 9;
end
if exist('SaveOutputVariables')
else
    SaveOutputVariables = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Conventions for Li, L, and t
% Li is defined as Li := S*(A_i ⊗ B_i)
% L is defined as L := \sum_i transpose(A_i{^⊗k1} ⊗ B_i{^⊗k2}) ⊗ R_i
% t is defined as \sum_i \tr(R_i L_i)
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
        %              	load('inputdata_4slots_2plus2_general') %Load variables
%         load('inputdata_4slots_2plus2_general_old.mat')
        load('inputdata_4slots_2plus2_general_ABBA_10^-5.mat')
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
        disp("The Switch constructed in this function coincides with the Switch from input_variables" + newline + "All good regarding the SWITCH! =)");
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
%%%%%%%%%%%% Evaluate t=\sum_i \tr(R_i L_i)
t=0;
for i=1:size(A,3)
    t=t+HS(Ri(:,:,i),Li(:,:,i));	%HS(A,B) = trace(A'B)
end
%%%%%%%%%%%% Basic tests to see that comes out from the SDP
t=t
minEigG=min(eig(G))
minEigGK=min(eig(G-L))
DIM = [d*ones(1,2*k) 2];
normG_minus_projG = norm(G-ProjSeqSuperChannelNoFuture4Sym(G,DIM))
p = trace(G)/2/d^k
%%%%%%%%%%%% Truncate Gamma and Ri
G_trunc = chopMTQdetail(G,DecimalPrecision);
Ri_trunc = chopMTQdetail(Ri,DecimalPrecision);
%%%%%%%%%%%% Make Gamma Self-Adjoint  (GSA stands for "Gamma Seft-Adjoint")
GSA = (G_trunc + G_trunc')/2;
%%%%%%%%%%%% Make Ri Self-Adjoint and Create L_SA := \sum_i transpose(A_i{^⊗k1} ⊗ B_i{^⊗k2}) ⊗ R_iSA
%%%%%%%%%%%% R_iSA stands for Ri (truncated) Self-Adjoint)
RiSA = zeros(size(Ri));
LSA = 0;
switch scenario
    case 'AB'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),B(:,:,i))),RiSA(:,:,i));
        end
    case 'AA'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'ABA'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),B(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'ABB'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),B(:,:,i),B(:,:,i))),RiSA(:,:,i));
        end
    case 'BAA'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(B(:,:,i),A(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'AAB'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),A(:,:,i),B(:,:,i))),RiSA(:,:,i));
        end
    case 'AAA'
        for i=1:size(A,3)
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),A(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'ABAB'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),B(:,:,i),A(:,:,i),B(:,:,i))),RiSA(:,:,i));
        end
    case 'ABBA'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),B(:,:,i),B(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'AABB'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),A(:,:,i),B(:,:,i),B(:,:,i))),RiSA(:,:,i));
        end
    case 'AAAB'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),A(:,:,i),A(:,:,i),B(:,:,i))),RiSA(:,:,i));
        end
    case 'AABA'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),A(:,:,i),B(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'ABAA'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(A(:,:,i),B(:,:,i),A(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    case 'BAAA'
        for i=1:size(A,3)
            %             if i==1
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            %             if mod(i,1000)==0
            %                 disp(['The variable i ranges from 1 to ',num2str(size(A,3)),'. The current value of i is ', num2str(i)])
            %             end
            RiSA(:,:,i) = (Ri_trunc(:,:,i) + Ri_trunc(:,:,i)')/2;
            LSA = LSA + kron(transpose(Tensor(B(:,:,i),A(:,:,i),A(:,:,i),A(:,:,i))),RiSA(:,:,i));
        end
    otherwise
        error('You chose a scenario which is not covered by this code')
end
%%%%%%%%%%%% Evaluate tSA  (the value of t with Self-Adjoint truncated variables)
tOK=0;
for i=1:size(A,3)
    tOK=tOK+HS(RiSA(:,:,i),Li(:,:,i));  %HS(A,B) = trace(A'B)
end

if abs(tOK-1)>10^(-4)
    tOKvalue = vpa(tOK)
    error("The value of tS is fair away from one.... this is strange... please check your code")
else
    tOKvalue = vpa(tOK)
    'The value t := \sum_i \tr(R_i L_i) is close to one, this is a good sign! '
    'The initial outcome of your SDP is close to respect the equality constraints =)'
end
%%%%%%%%%%%% Create RiOK (just divide RiSA by tOK)
RiOK = RiSA;  %This is done just to initialise RiSOK with the correct dimension (everything will be overwritten anyway)
for i=1:size(A,3)
    RiOK(:,:,i) = RiSA(:,:,i)/tOK;
end
%%%%%%%%%%%% Create LOK (just divide LSA by tOK) L := \sum_i transpose(A_i{^⊗k1} ⊗ B_i{^⊗k2}) ⊗ R_i
LOK = LSA/tOK;
%%%%%%%%%%%% Double check that final t is indeed 1, t1 is the "t that should be one"
t1 = 0;
for i=1:size(A,3)
    t1=t1+HS(RiOK(:,:,i),Li(:,:,i)); %HS(A,B) = trace(A'B)
end
%%%%%%%%%%%% Project Gamma into the valid space
GValid = ProjSeqSuperChannelNoFuture4Sym(GSA,DIM);
%%%%%%%%%%%% Chose eta and make the final GOK
minEigGValidDouble = min(real(eig(GValid)))
minEigGSVALID_minus_LS= min(real(eig(double(GValid) - LOK)))
lambda = abs(min(minEigGValidDouble,minEigGSVALID_minus_LS));
eta=ceil(lambda*10^DecimalPrecision)/10^DecimalPrecision

GOK = GValid + eta*(eye(size(G))); %Add eta white noise into Gvalid to create GOK
minEigGOK = min(eig(double(GOK))) %Check the min eig of GOK
minEigGOK_minus_LOK = min(eig(GOK - LOK))  %Check the min eig of GOK-LOK
pHP = trace(GOK)/2/d^k %Value of p with high-precision

isEverythingOK = 1;
if abs(t1-1)<10^(-DecimalPrecision)
    disp('Great, t1:=\sum_i \tr(RiSOK Li) is close to one!!')
    t1=vpa(t1)
else
    disp(['Failure... t1:=\sum_i \tr(RiSOK Li) is NOT equals to one!!', newline, 'Its value is something like', num2str(double(t1))]);
    t1=vpa(t1)
    isEverythingOK=0;
end
% [IsPSDOK_GOK,GOKsqrt] = IsPSDSym(GOK);
if real(minEigGOK)>0 && abs(imag(minEigGOK))<10^(-10)
    disp(['Great, GOK (gamma) is a very close to a PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGOK)]);
else
    disp(['Failure... GOK (gamma) is NOT close to a PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGOK)])
    isEverythingOK=0;
end
% [IsPSDOK_GOKminusLOK,GOKminusLOKsqrt] = IsPSDSym(GOK-LOK);
if real(minEigGOK_minus_LOK)>0 && abs(imag(minEigGOK_minus_LOK))<10^(-10)
    disp(['Great, GOK-LOK is close to a PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGOK_minus_LOK)]);
else
    disp(['Failure... GOK-LSA is is NOT close to a PSD matrix!!', newline, 'Its minimal eigenvalue is something like ', num2str(minEigGOK_minus_LOK)]);
    isEverythingOK=0;
end

% if SaveOutputVariables
% 	disp(['We will now save important varaibles, this may take a very long time. (The current scenario is: ',scenario,')'])
% 	!date
% 	% We will save the date in a format like outputdata_3slots_2plus1_generalABA
% 	tic
% 	kA=count(scenario, 'A');    	kB=count(scenario, 'B');
% 	kMax = max(kA,kB);          	kMin=k-kMax;
% 	saveFileName = ["outputdata_"+num2str(k)+"slots_"+num2str(kMax)+"plus"+num2str(kMin)+"_general"+scenario+"_DecimalPrecision_"+num2str(DecimalPrecision)];
% 	save(saveFileName,'GOK','RiSOK','pHP','AS','BS','GOKsqrt','GOKminusLSsqrt','TotalTimeMaincode')
% 	Time2saveImportantVariables = toc
% end

if SaveOutputVariables
    % 	tic
    save("HighPrecisionSwitchSimulationFullData_"+scenario+"_DecimalPrecision_"+num2str(DecimalPrecision))
    disp('The HighPrecisionFullOutputData data was saved! You have all info!!')
    % 	Time2saveALLVariables = toc
end

if isEverythingOK
    disp(['Everything is OK!!!', newline, '=)', newline, 'a HIGH-PRECISION pHP for the scenario ',scenario,' is something like ', num2str(pHP), newline, 'You may want to compare it with the initial p (from the SDP), which is ', num2str(p)]);
else
    error('DANGER, Gamma is not positive or Gamma - LSA is not positive...')
end
TotalTimeHPfunctionSeconds = toc(TotalTimeThisFunction)
end % End of main function

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