function [rfVirtual, vfVirtual, afVirtual, tgoVirtual] = computeBeyondTerminationTargeting(r, v, gamma, kr, rfStar, vfStar, afStar, delta_t, tgo_true, rMoonND)
    r = r(:); v = v(:);
    rfStar = rfStar(:); vfStar = vfStar(:); afStar = afStar(:);

    if gamma < 0, error('gamma must be >= 0'); end
    if kr < 2*(gamma + 2), error('kr must be >= 2*(gamma+2)'); end

    tgoVirtual = tgo_true + delta_t;
    tgo = tgoVirtual;

    if delta_t < 1e-15
        tgoVirtual = tgo_true;
        rfVirtual = rfStar;
        afVirtual = afStar;
        vfVirtual = vfStar;
        return
    end
    gamma1 = gamma;
    gamma2 = kr/(gamma + 2) - 2;

    if abs(gamma1 - gamma2) < 1e-15
        % kr = (gamma+2)^2 case: BTT should not be done per paper
        rfVirtual = rfStar; vfVirtual = vfStar; afVirtual = afStar;
        tgoVirtual = tgo_true;
        return;
    end

    % Local gravity from current position (matches Fortran reference)
    g = -(rMoonND^2) * r / (norm(r)^3);

    % Basis functions evaluated at tgo = tgo_true + delta_t
    phi1 = tgo^gamma1;
    phi2 = tgo^gamma2;
    phi1_bar = -(1/(gamma1 + 1)) * tgo^(gamma1 + 1);
    phi2_bar = -(1/(gamma2 + 1)) * tgo^(gamma2 + 1);
    phi1_hat = tgo^(gamma1 + 2) / ((gamma1 + 1)*(gamma1 + 2));
    phi2_hat = tgo^(gamma2 + 2) / ((gamma2 + 1)*(gamma2 + 2));
    delta = phi1_hat * phi2_bar - phi2_hat * phi1_bar;

    k1r = -phi2_bar / delta;
    k1v = (phi2_bar * tgo + phi2_hat) / delta;
    k1a = -(0.5 * tgo * phi2_bar + phi2_hat) * tgo / delta;

    k2r =  phi1_bar / delta;
    k2v = -(phi1_bar * tgo + phi1_hat) / delta;
    k2a =  (0.5 * tgo * phi1_bar + phi1_hat) * tgo / delta;

    d1 = ((r - 0.5*g*tgo^2).*phi2_bar - (v + g*tgo).*phi2_hat)/delta;
    d2 = -((r - 0.5*g*tgo^2).*phi1_bar - (v + g*tgo).*phi1_hat)/delta;

    % Basis functions re-evaluated at delta_t
    phi1 = delta_t^gamma1;
    phi2 = delta_t^gamma2;
    if gamma2 <= 1e-15 && delta_t < 1e-15
        phi2 = 0;
    end
    phi1_bar = -(1/(gamma1 + 1)) * delta_t^(gamma1 + 1);
    phi2_bar = -(1/(gamma2 + 1)) * delta_t^(gamma2 + 1);
    phi1_hat = (1/((gamma1 + 1)*(gamma1 + 2))) * delta_t^(gamma1 + 2);
    phi2_hat = (1/((gamma2 + 1)*(gamma2 + 2))) * delta_t^(gamma2 + 2);

    m11 = k1r * phi1_hat + k2r * phi2_hat + 1;
    m12 = k1v * phi1_hat + k2v * phi2_hat - delta_t;
    m13 = k1a * phi1_hat + k2a * phi2_hat + 0.5 * delta_t^2;

    m21 = k1r * phi1_bar + k2r * phi2_bar;
    m22 = k1v * phi1_bar + k2v * phi2_bar + 1;
    m23 = k1a * phi1_bar + k2a * phi2_bar - delta_t;

    b1 = phi1_hat * d1 + phi2_hat * d2 + 0.5 * g * delta_t^2;
    b2 = phi1_bar * d1 + phi2_bar * d2 - g * delta_t;

    I3 = eye(3);

    % Third row: acceleration at delta_t from virtual endpoint equals afStar
    % aT(delta_t) = afV + c1*delta_t^gamma1 + c2*delta_t^gamma2
    % where c1 = k1r*rfV + k1v*vfV + k1a*afV + d1 (and similarly c2)
    m31 = k1r * phi1 + k2r * phi2;
    m32 = k1v * phi1 + k2v * phi2;
    m33 = k1a * phi1 + k2a * phi2 + 1;

    b3 = d1 * phi1 + d2 * phi2;

    % 9x9 system: solve for rfVirtual, vfVirtual, afVirtual simultaneously
    MM = [m11*I3, m12*I3, m13*I3;
          m21*I3, m22*I3, m23*I3;
          m31*I3, m32*I3, m33*I3];
    bb = [rfStar - b1;
          vfStar - b2;
          afStar - b3];

    sol = MM \ bb;
    rfVirtual = sol(1:3);
    vfVirtual = sol(4:6);
    afVirtual = sol(7:9);
