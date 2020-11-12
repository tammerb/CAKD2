%AA_Spatial_Orientation_Transformation_Matrix_and_Inverse for p
%{
Author: Edward J. Haug, Blake Hannah, Victor Honein, Michael Horton,
Abhay Negi, Trevor Vidano
Institution: University of California, Davis
Location: Davis, CA 95616
Date: October 23, 2020
%}
%% How to use this code:
%{
This script is a user command line interface for computing the spatial
orientation transformation matrix: A and it's Euler paramterization: p.
This script prompts the user for the proper inputs to pass into
convertRotationMatrix(). 

When prompted for input mode, mode 1 requires the user to provide A 
(orthogonal transformation matrix); mode 2 requires the user to provide rO,
rP, and rQ; mode 3 requires the user to provide uBar, and chi.

When prompted for output mode, mode 1 prints the A matrix to the command
window, mode 2 prints the p matrix to the command.
%}

%% User input section:

inputMode = -1;
inputModePrompt = 'Enter input mode (1, 2, or 3): ';
while not(ismember(inputMode,[1,2,3]))
    inputMode = input(inputModePrompt);
end

% Prompt user based on inputMode.
if inputMode==1         % Prompt user to provide an A matrix.
    A = input('Enter your A matrix: ');
    [A,p] = convertRotationMatrix('A',A);
end

if inputMode==2         % Prompt user for point definition (Section 2.5.5.)
    firstHalfPrompt = 'Enter vector from global origin to';
    rO = input(strcat(firstHalfPrompt,' body-fixed origin (rO): '));
    rP = input(strcat(firstHalfPrompt," a point on x' axis (rP): ")); 
    rQ = input(strcat(firstHalfPrompt," a point on x'-y' plane (rQ): "));
    [A,p] = convertRotationMatrix('rO',rO,'rP',rP,'rQ',rQ);
end

if inputMode==3         % Prompt user for definition from Eulers Theorem of 
                        % Section 2.5.3.
    ubar = input("Enter unit vector u about which to rotate x-y-z (uBar): ");             
    chi = input("Enter amount of rotation in radians (chi): ");
    [A,p] = convertRotationMatrix('ubar',ubar,'chi',chi);
end

outputMode = -1;
outputModePrompt = 'Enter output mode (1 or 2): ';
while not(ismember(outputMode,[1,2]))
    outputMode = input(outputModePrompt);
end
%% Print output based on user input.
if outputMode==1  % Present A (orthogonal transformation matrix).
    A
end

if outputMode==2  % Present p (Euler parameterization of A).
    p
end