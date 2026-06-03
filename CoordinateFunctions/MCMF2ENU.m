function [enu, alt] = MCMF2ENU(X, landingLatDeg, landingLonDeg, isPosition, radBody)

lat0 = deg2rad(landingLatDeg);
lon0 = deg2rad(landingLonDeg);

[E0, N0, U0] = enuBasis(lat0,lon0);
R = [E0, N0, U0]';

if isPosition
    Xorig = X;
    r0 = radBody * U0; % position of landing site in MCMF
    X = X - r0;
    if nargout > 1
        alt = vecnorm(Xorig, 2, 1) - radBody;
    end
end

enu = R * X;
end
