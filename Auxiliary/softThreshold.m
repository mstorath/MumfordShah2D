function y = plSoftThreshold( x, tau )
%plSoftThreshold the classical soft-threshold function

y = max( abs(x) - tau, 0) .* sign(x);
end

