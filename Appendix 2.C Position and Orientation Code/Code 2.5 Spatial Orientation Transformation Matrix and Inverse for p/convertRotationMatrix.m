function [A,p] = convertRotationMatrix(varargin)
% Pass in arguments like this: convert2EulerVector('A',eye(3))
parser = inputParser;
parser.CaseSensitive = true;
validMatrix = @(x) isnumeric(x) && isequal(size(x),[3,3]);
validVector = @(x) isnumeric(x) && isequal(size(x),[3,1]);
validScalar = @(x) isnumeric(x) && isscalar(x);
validEulerVector = @(x) isnumeric(x) && isequal(size(x),[4,1]);
addParameter(parser,'A',false,validMatrix);
addParameter(parser,'rO',false,validVector);
addParameter(parser,'rP',false,validVector);
addParameter(parser,'rQ',false,validVector);
addParameter(parser,'ubar',false,validVector);
addParameter(parser,'chi',false,validScalar);
addParameter(parser,'p',false,validEulerVector);
parse(parser,varargin{:});

% Unpack arguments.
A = parser.Results.A;
rO = parser.Results.rO;
rP = parser.Results.rP;
rQ = parser.Results.rQ;
ubar = parser.Results.ubar;
chi = parser.Results.chi;
p = parser.Results.p;

% Enforce mutually exclusive arguments: {A},{rO,rP,rQ},{ubar,chi},{p}
mode = -1;
if validMatrix(A)
    assert(islogical([rO,rP,rQ,ubar,chi,p]), 'too many arguments given.')
elseif any([validVector(rO),validVector(rP),validVector(rQ)])
    assert(all([validVector(rO),validVector(rP),validVector(rQ)]),...
        'not enough arguments provided. Define rO, rP, and rQ')
    assert(islogical([A,ubar,chi,p]), 'too many arguments given.')
    mode = 2;
elseif validVector(ubar) || validScalar(chi)
    assert(all([validVector(ubar),isscalar(chi)]), ...
        'not enough arguments provided. Define ubar and chi.')
    assert(islogical([A,rO,rP,rQ,p]), 'too many arguments given.')
    mode = 3;
elseif validEulerVector(p)
    assert(islogical([A,rO,rP,rQ,ubar,chi]), 'too many arguments given.')
    mode = 4;
elseif islogical([A,rO,rP,rQ,ubar,chi,p])
    errorMsg = ['No arguments given. Please define one of the following'...
        ' sets: {A},{rO,rP,rQ},{ubar,chi},{p}.'];
    error(errorMsg);
end

% Create A matrix if undefined.
if A==0
    % Logic to begin mode 2: Produce A give rQ, rP, rO
    if mode==2
        
        % Enforce rO, rP, and rQ definitions. (Section 2.5.5)
        assert(any(abs(rQ - rO) > ones(3,1)*10^-6),'rQ must not equal rO');
        assert(any(abs(rP - rO) > ones(3,1)*10^-6),'rP must not equal rO');
        assert(any(abs(rP - rQ) > ones(3,1)*10^-6),'rP must not equal rQ');
        assert(any(mod((rP-rO),(rQ-rO)) < ones(3,1)*10^-6),...
            "Point Q must not be on x' axis");
        
        %Calculate A.
        f=(1/norm(rP-rO)*(rP-rO));
        h=(1/norm(atil(f)*(rQ-rO)))*atil(f)*(rQ-rO);
        g=-atil(f)*h;

        A=[(f),(g),(h)];
        
    % Logic to initiate mode 3: Produce A given ubar and chi
    elseif mode==3
    
        u=(1/norm(ubar))*ubar;
        e0=cos(chi/2);
        e=sin(chi/2)*u;
        
        A=(e0^2-e'*e)*eye(3)+2*(e*e')+2*e0*atil(e);
    end
end

%% Given orthogonal A, compute Euler parameter vector p and vice versa

if mode==4                  % First checks if p was entered by the user
    A=eulerToA(p);          % Return tranfromation matix A if true
elseif issymmetric(A)       % Computes Euler parameter when Transpose A = A
    e0=1;
    e=zeros(3,1);
    p=[e0;e];
else                        % Computes Euler parameter from A, Eq(2.5.30)
    trA=trace(A);
    e0=0.5*sqrt(trA+1);
    e=(1/(2*sqrt(trA+1)))*[A(3,2)-A(2,3);A(1,3)-A(3,1);A(2,1)-A(1,2)];
    p=[e0;e];
end
end