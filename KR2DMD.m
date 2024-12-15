clear all, close all, clc
%%Dynamic mode decomposition
%A method of correlating dynamic systems (systems with many variables and
%no governing equation) through creating the best possible psuedolinear 
%dynamical system, sort of like a hybrid between an optimization problem
%and linear regression.

%%Load the data and convert the data to the right format
%Load the excel spreadsheet with all the data from the weather station 
[sFile, sPath] = uigetfile ('*xlsx', 'Select Data File');
sSelectedfile = fullfile (sPath,sFile);
%Use the function readtable to turn the selected file into a table
mTable = readtable(sSelectedfile);

%Turn the table into a matrix
mTable1=mTable{:,:};
%Transpose the matrix so that the columns vary by time 
mTable1=mTable1.';
%Replace NaN values with 0, if applicable and valid 
mTable1(isnan(mTable1)) = 0;
%Find the size of the rows and columns of the table using the size function
iRows = size(mTable1, 1);
iCols = size(mTable1, 2);

%Create a time vector outlining the time's represented by each column. In
%this example the unit of time is days so we are creating a time vector
%where we start at 0 days, end at 91 days, and points for each column and
%zero (iCols+1)
t=linspace(0,91,iCols+1);%13 values from 0 to 91 in days 

%% Create X and Xprime matrixes for a singular value decomposition test
%Singular value decomposition tests generalize the eigendecomposition of a
%matrix into two main veoctors (U and V*(V)) and a matrix (S). 
%Eigendecomposition is a method of matrix vectorization that produces
%eigenvectors and eigenvalues 
%Calculate the singular value decomposition of the matrix, use 'econ' to
%run an economy-size decomposition of the matrix 
[U,S,V] = svd(mTable1,'econ');
%% plot svd 
%plot diagonal data
%Figure 1 plots the sigma matrix log scale on y with index as x. The elbows
%show harmonic eigan values pairs and then with that you can decide how 
%many modes to keep and which modes to truncate at.
%Check these values with plot two that shows individual points indicating
%possible modes.
figure (1)
semilogy(diag(S),'LineWidth',2)
figure(2)
plot(diag(S)/(sum(diag(S))),'ro')

%Truncate modes at the value seen on figures one and two that represent the
%most amount of data without including unnecessary modes
r = 6;  

%%Here we plot the real U and V* vector values from 1 to r in an effort to
%%illustrate the lack of fixed oscillation necessary for creating a
%%psuedolinear equation necessary for modeling this data
figure(3)
subplot(2,1,1),plot(real(U(:,(1:r))))
subplot(2,1,2),plot(real(V(:,(1:r))))

%Create an X matrix that goes from the first time column to the second to
%last time column and an Xprime matrix that goes from the second time
%column to the last column
X = mTable1(:,(1:end-1));
Xprime = mTable1(:,(2:end));

%%run singular value decomposition test on just the X matrix
[u,s,v]=svd(X,'econ');  
%Truncate the vectors and matrices given from the SVD calculation at the r
%value
uR=u(:,1:r);
sR=s(1:r,1:r);
vR=v(:,1:r);
%%  Compute DMD (Phi are eigenvectors)
%Atilde is the best fit linear model that tells me how these 6 principal
%orthogonal decomposition modes are evolving in time 
aTilde = uR'*Xprime*vR/sR;
%compute eigan decomposition: to get eigan vectors (W) and eigan
%values(D)
[W,D]=eig(aTilde);
%reconstruct full DMD eiganvectors from W
Phi=Xprime*vR/sR*W; %DMD modes
%set dt to the change in time between n and n+1
dt=7;
%Calculate lambda and omega values needed for final DMD calculation
lambda=diag(D);
omega=log(lambda)/dt;

%Plot the real phi values to show the 2D version of the function we will
%use encapsulating the DMD modes 
figure(4)
plot(real(Phi),'LineWidth',[2])

%create a vector with the first row of our initial matrix
X1=mTable1(:,1);
%create a b vector highlighting the b values for each DMD mode
b=Phi\X1;

%Create a new time vector with iCol values from 0 to 91 to predict values
%inside the matrix to determine if it is working
t1=linspace(0,91,iCols);

%create a zero vector of time_dynamics1 from r to the length of the t1
%vector
time_dynamics1=zeros(r,length(t1));
%create a for loop that starts at one and advances by values in the length
%of t1 using the basic equation for DMD for each value in the column
%vectors
for iter=1:length(t1)
    time_dynamics1(:,iter)=(b.*exp(omega*t1(iter)));
end

%Reform the matrix of values by multiplying the time_dynamics value by the
%Phi matrix.
x_dmd1=Phi*time_dynamics1;

%Create a new time vector with iCol values from 0 to 147 to predict values
%inside the matrix to determine the predicted values
t2=linspace(0,147,iCols);
%apply the same process above to predict values
time_dynamics2=zeros(r,length(t2));
for iter=1:length(t2)
    time_dynamics2(:,iter)=(b.*exp(omega*t2(iter)));
end

x_dmd2=Phi*time_dynamics2;

%Display visualization of what DMD looks like through time in both cases.
%Look for similarities and differences
figure(5)
subplot(2,1,1),surfl(real(x_dmd1));shading interp,colormap("gray")
subplot(2,1,2),surfl(real(x_dmd2));shading interp,colormap("gray")

%%Load excel spreadsheet with actual values
[sFile, sPath] = uigetfile ('*xlsx', 'Select Data File');
sSelectedfile = fullfile (sPath,sFile);
%Use the function readtable to turn it into a table
mTableExp = readtable(sSelectedfile);
%Turn this table into a matrix and transpose it 
mExp=mTableExp{:,:};
mExp=mExp.';
mExp(isnan(mExp)) = 0;

%Create a vector with only the real values of x_dmd2 from the last column
%(our target column)
predictSep=real(x_dmd2(:,13));

%create a matrix with the expected and predicted vectors and add column
%names
actExpSepValues = [mExp(:), predictSep(:)];
colNames={'Actual Values', 'DMD expected values'};
%turn this matrix into a table
actExpSepTable = array2table(actExpSepValues,'VariableNames',colNames);

%%Do a correlation test to determine if there is valid correlation in this
%%data
%Valid correlation occures when the rho value is close to 1, a positive
%value shows a positive correlation and a negative rho value indicates a
%negative correlation
%We also tested for the p-value to determine if our data is statistically
%significant. if our pval is less than 0.05 we can assert that it
%statistically significant.
[rho,pval] = corr(mExp,predictSep);

%Plot this data to visualize correlations 
figure (6)
plot(actExpSepValues(:,1),actExpSepValues(:,2),'o')

