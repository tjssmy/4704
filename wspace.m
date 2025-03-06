function w = wspace(t,nt); %t=time matrix with values at every point, nt = size of t

if (nargin<2) %if nt is not provided
    nt = length(t); %nt = length of t
    dt = t(2) - t(1); %calculates the step time
    t = t(nt) - t(1) + dt;
end

if (nargin == 2) %if nt is provided
    dt = t/nt; %calculates the step time
end

w = 2*pi*(0:nt-1)'/t; %computes a vector of angular frequencies for the time vector

%figure; plot(w);

kv = find(w >= pi/dt); 
w(kv) = w(kv) - 2*pi/dt;

%figure; plot(kv);
%figure; plot(w);
