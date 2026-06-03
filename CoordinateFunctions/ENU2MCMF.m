function mcmf = ENU2MCMF(enu, anchorLatDeg, anchorLonDeg, isPosition, radBody)

    lat0 = deg2rad(anchorLatDeg);
    lon0 = deg2rad(anchorLonDeg);

    [E0, N0, U0] = enuBasis(lat0, lon0);
    R = [E0, N0, U0]';

    if isPosition
        r0 = radBody * U0;
        mcmf = R.' * enu + r0;
    else
        mcmf = R.' * enu;
    end
end
