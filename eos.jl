using Polynomials
using Polynomials: derivative


function roots_abc(p::Polynomial)
    """Roots of a quadratic polynomial p using the abc formula.
    
    Warning: This simple method is not guaranteed to be numerically stable.
    """
    c, b, a = p.coeffs
    # https://en.wikipedia.org/wiki/Quadratic_equation#Avoiding_loss_of_significance
    Δ = b^2 - 4a*c  # TODO the discriminant is prone to cancellation errors.
    if b ≥ 0
        x1 = -2c / (b + sqrt(Δ))
    else
        x1 = -2c / (b - sqrt(Δ))
    end
    x2 = (-b + sqrt(Δ)) / (2a)
    return [x1, x2]
end

function eos_birch_murnaghan_fit(volumes, energies; root_method=Polynomials.roots)
    # Based on https://github.com/aiidateam/acwf-verification-scripts/blob/main/3-analyze/eos_utils/eosfit_31_adapted.py
    p = fit(volumes.^(-2.0/3.0), energies, 3)
    deriv1 = derivative(p)
    deriv2 = derivative(deriv1)
    deriv3 = derivative(deriv2)

    # Find extrema by first derivative condition
    p_extrema = root_method(deriv1)  

    # Select local minimum by a second derivative test
    idx = findfirst(x -> x > 0 && deriv2(x) > 0, p_extrema)
    isnothing(idx) && error("BM fit failed: No minimum could be found")
    x = p_extrema[idx]
    volume0 = x^(-3.0/2.0)
    
    e0 = p(x)

    derivV2 = 4. / 9. * x^5. * deriv2(x)
    derivV3 = (-20. /  9. * x^(13. / 2.) * deriv2(x)
               - 8. / 27. * x^(15. / 2.) * deriv3(x))
    bulk_modulus0 = derivV2 / x^(3. / 2.)
    bulk_deriv0 = -1 - x^(-3. / 2.) * derivV3 / derivV2

    (; volume0, e0, bulk_modulus0, bulk_deriv0)
end

function eos_birch_murnaghan_energy(v; volume0, e0, bulk_modulus0, bulk_deriv0)
    r = (volume0 / v)^(2.0/3.0)
    (e0 + 9.0/16.0 * bulk_modulus0 * volume0 * (
            (r-1.)^3 * bulk_deriv0 + 
            (r-1.)^2 * (6.0 - 4.0 * r)))
end

function nu_metric(eos1, eos2; prefactor=100.0, weight_volume=1.0,
                   weight_bulk_modulus=1/20, weight_bulk_deriv=1/400)
    prefactor * sqrt(
        weight_volume * (eos1.volume0 - eos2.volume0)^2
        + weight_bulk_modulus * (eos1.bulk_modulus0 - eos2.bulk_modulus0)^2
        + weight_bulk_deriv * (eos1.bulk_deriv0 - eos2.bulk_deriv0)^2
    )
end
