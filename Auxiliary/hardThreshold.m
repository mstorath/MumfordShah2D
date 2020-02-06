function y = plHardThreshold( x, tau )
%plHardThreshold the classical hard-threshold function
y = x .* (abs(x) >= tau);
end

