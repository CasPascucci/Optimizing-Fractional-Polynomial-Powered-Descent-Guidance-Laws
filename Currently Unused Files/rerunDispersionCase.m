function outputSingle = rerunDispersionCase(idx, PDINom, planetaryParams, targetState, vehicleNom, optimizationParams, beta, seeds)
    % 3-sigma ranges
    r_disp = 200; % m
    lon_disp = 0.25; % deg
    lat_disp = 0.25; % deg
    v_disp = 3; % m/s
    fpa_disp = 0.1; % deg
    azmth_disp = 0.2; % deg
    mass_disp = 0.0066; % fraction of nominal mass

    dalt = seeds.alt_seeds(idx) * (r_disp / 3);
    dlon = seeds.lon_seeds(idx) * (lon_disp / 3);
    dlat = seeds.lat_seeds(idx) * (lat_disp / 3);
    dv = seeds.v_seeds(idx) * (v_disp / 3);
    dfpa = seeds.fpa_seeds(idx) * (fpa_disp / 3);
    dazmth = seeds.azmth_seeds(idx) * (azmth_disp / 3);
    mass_mult = seeds.mass_seeds(idx) * (mass_disp / 3);

    PDI = PDINom;
    PDI.altitude_km = PDINom.altitude / 1000 + dalt/1000;
    PDI.lonInitDeg = PDINom.lonInitDeg + dlon;
    PDI.latInitDeg = PDINom.latInitDeg + dlat;
    PDI.inertialVelocity = PDINom.inertialVelocity + dv;
    PDI.flightPathAngleDeg = PDINom.flightPathAngleDeg + dfpa;
    PDI.azimuth = PDINom.azimuth + dazmth;

    vehicle = vehicleNom;
    vehicle.massInit = vehicleNom.massInit * (1+ mass_mult);
    vehicle.dryMass = vehicle.massInit - 8248;

    [gammaOpt, gamma2Opt, krOpt, tgoOpt, ~, ~, optFuel, simFuel, ~, finalPosSim, ~, ~] = getParams(PDI, planetaryParams, targetState, vehicle, optimizationParams, beta, true, true, false, true);
   
    outputSingle = struct;
    outputSingle.gamma = gammaOpt;
    outputSingle.gamma2 = gamma2Opt;
    outputSingle.kr = krOpt;
    outputSingle.tgo = tgoOpt;
    outputSingle.optFuel = optFuel;
    outputSingle.simFuel = simFuel;
    outputSingle.finalError = finalPosSim;
end