# -*- coding: utf-8 -*-
"""
Created on Wed Jul  6 14:48:43 2022
Updated on Mon Aug 19 15:21:29 2024

@author: abosserelle, Jules
"""

import numpy as np
import math

"""
This script provides an analytical solution for groundwater response to sea-level rise (SLR)
in a sloping coastal aquifer. The calculations are based on the work of Morgan & Werner (2016),
who commented on Chesnaux (2015). The example is from the Pioneer Valley, Australia.

References:
- Morgan, L. K. and A. D. Werner (2016). "Comment on “Closed-form analytical solutions for
  assessing the consequences of sea-level rise on groundwater resources in sloping coastal
  aquifers”: paper published in Hydrogeology Journal (2015) 23:1399–1413, by R. Chesnaux."
  Hydrogeology journal 24(5): 1325-1328.
- Werner & Gallagher (2006). "Characterisation of sea-water intrusion in the Pioneer Valley,
  Australia using hydrochemistry and three-dim."
"""

# --- Aquifer Parameters ---
K = 100         # Hydraulic conductivity (m/d)
Z0 = 25         # Aquifer thickness, depth to the aquifer base from mean sea level (m)
RHO_F = 1000    # Density of freshwater (kg/m^3)
RHO_S = 1025    # Density of saltwater (kg/m^3)
W = 0.11/365.25 # Net recharge (m/d)

# --- Observation Well Data ---
HB = 2          # Hydraulic head at observation well (m above MSL)
XB = 2000       # Distance of observation well from the coast (m)

# --- Sea Level Rise Scenario ---
DELTA_Z = 1     # Sea Level Rise (m)

def calculate_delt():
    """Calculates the density difference ratio."""
    return (RHO_S - RHO_F) / RHO_F

def calculate_coastal_flow(hb, xb, z0, k, w, delt):
    """
    Calculates the lateral groundwater flow at the coast (q0).

    Args:
        hb (float): Hydraulic head at observation well (m).
        xb (float): Distance of observation well from coast (m).
        z0 (float): Aquifer thickness (m).
        k (float): Hydraulic conductivity (m/d).
        w (float): Net recharge (m/d).
        delt (float): Density difference ratio.

    Returns:
        float: Discharge to the sea (m^2/d).
    """
    qb = ((((hb + z0)**2 - ((1 + delt) * z0**2)) * k) - (w * xb**2)) / (2 * xb)
    q0 = qb + (w * xb)
    return q0

def calculate_groundwater_divide(q0, w):
    """
    Calculates the distance from the coast to the groundwater divide.

    Args:
        q0 (float): Discharge to the sea (m^2/d).
        w (float): Net recharge (m/d).

    Returns:
        float: Distance to groundwater divide (m).
    """
    return q0 / w

def calculate_mixed_convection_ratio(k, delt, z0, w, xn):
    """
    Calculates the mixed convection ratio (M), a measure of SWI vulnerability.

    Args:
        k (float): Hydraulic conductivity (m/d).
        delt (float): Density difference ratio.
        z0 (float): Aquifer thickness (m).
        w (float): Net recharge (m/d).
        xn (float): Distance to groundwater divide (m).

    Returns:
        float: Mixed convection ratio.
    """
    return k * delt * (1 + delt) * z0**2 / (w * xn**2)

def calculate_interface_toe_position(xn, m):
    """
    Calculates the position of the saltwater interface toe.

    Args:
        xn (float): Distance to groundwater divide (m).
        m (float): Mixed convection ratio.

    Returns:
        float: Position of the interface toe from the coast (m).
    """
    if 1 - m < 0:
        return xn # Toe is at the divide if M >= 1
    return xn * (1 - np.sqrt(1 - m))

def calculate_head_at_x(x, q0, w, k, z0, delt):
    """
    Calculates the hydraulic head at a distance x from the coast.

    Args:
        x (float): Distance from the coast (m).
        q0 (float): Discharge to the sea (m^2/d).
        w (float): Net recharge (m/d).
        k (float): Hydraulic conductivity (m/d).
        z0 (float): Aquifer thickness (m).
        delt (float): Density difference ratio.

    Returns:
        float: Hydraulic head (m).
    """
    term1 = (2 * q0 * x - w * x**2) / k
    term2 = (1 + delt) * z0**2
    return np.sqrt(term1 + term2) - z0

def analyze_flux_controlled(q0, w, k, z0, delt, delta_z):
    """
    Analyzes the impact of SLR under flux-controlled conditions.

    Args:
        q0 (float): Pre-SLR discharge to the sea (m^2/d).
        w (float): Net recharge (m/d).
        k (float): Hydraulic conductivity (m/d).
        z0 (float): Pre-SLR aquifer thickness (m).
        delt (float): Density difference ratio.
        delta_z (float): Sea level rise (m).

    Returns:
        dict: A dictionary containing post-SLR results.
    """
    z0_post_slr = z0 + delta_z
    xn = calculate_groundwater_divide(q0, w) # Unchanged

    m_post_slr = calculate_mixed_convection_ratio(k, delt, z0_post_slr, w, xn)
    xt_post_slr = calculate_interface_toe_position(xn, m_post_slr)

    return {
        "xn_post_slr": xn,
        "m_post_slr": m_post_slr,
        "xt_post_slr": xt_post_slr
    }

def analyze_head_controlled(hb, xb, w, k, z0, delt, delta_z):
    """
    Analyzes the impact of SLR under head-controlled conditions.

    Args:
        hb (float): Hydraulic head at observation well (m).
        xb (float): Distance of observation well from coast (m).
        w (float): Net recharge (m/d).
        k (float): Hydraulic conductivity (m/d).
        z0 (float): Pre-SLR aquifer thickness (m).
        delt (float): Density difference ratio.
        delta_z (float): Sea level rise (m).

    Returns:
        dict: A dictionary containing post-SLR results.
    """
    z0_post_slr = z0 + delta_z

    # Recalculate coastal flow q0 with the new z0
    qb_post_slr = ((((hb + z0)**2 - ((1 + delt) * z0_post_slr**2)) * k) - (w * xb**2)) / (2 * xb)
    q0_post_slr = qb_post_slr + (w * xb)

    xn_post_slr = calculate_groundwater_divide(q0_post_slr, w)
    m_post_slr = calculate_mixed_convection_ratio(k, delt, z0_post_slr, w, xn_post_slr)
    xt_post_slr = calculate_interface_toe_position(xn_post_slr, m_post_slr)

    return {
        "q0_post_slr": q0_post_slr,
        "xn_post_slr": xn_post_slr,
        "m_post_slr": m_post_slr,
        "xt_post_slr": xt_post_slr
    }


if __name__ == '__main__':
    # --- Pre-SLR Analysis ---
    delt = calculate_delt()
    q0_pre_slr = calculate_coastal_flow(HB, XB, Z0, K, W, delt)
    xn_pre_slr = calculate_groundwater_divide(q0_pre_slr, W)
    m_pre_slr = calculate_mixed_convection_ratio(K, delt, Z0, W, xn_pre_slr)
    xt_pre_slr = calculate_interface_toe_position(xn_pre_slr, m_pre_slr)

    print("--- Pre-SLR Conditions ---")
    print(f"Discharge to the sea (q0): {q0_pre_slr:.4f} m^2/d")
    print(f"Groundwater divide (xn): {xn_pre_slr:.2f} m")
    print(f"Mixed convection ratio (M): {m_pre_slr:.4f}")
    print(f"Interface toe position (xt): {xt_pre_slr:.2f} m")
    print(f"Ratio xt/xn: {xt_pre_slr/xn_pre_slr:.4f}")

    # --- Post-SLR Analysis: Flux-Controlled ---
    flux_results = analyze_flux_controlled(q0_pre_slr, W, K, Z0, delt, DELTA_Z)
    print("\n--- Post-SLR: Flux-Controlled ---")
    print(f"Groundwater divide (xn): {flux_results['xn_post_slr']:.2f} m (Unchanged)")
    print(f"Mixed convection ratio (M): {flux_results['m_post_slr']:.4f}")
    print(f"Interface toe position (xt): {flux_results['xt_post_slr']:.2f} m")


    # --- Post-SLR Analysis: Head-Controlled ---
    head_results = analyze_head_controlled(HB, XB, W, K, Z0, delt, DELTA_Z)
    print("\n--- Post-SLR: Head-Controlled ---")
    print(f"Discharge to the sea (q0): {head_results['q0_post_slr']:.4f} m^2/d")
    print(f"Groundwater divide (xn): {head_results['xn_post_slr']:.2f} m")
    print(f"Mixed convection ratio (M): {head_results['m_post_slr']:.4f}")
    print(f"Interface toe position (xt): {head_results['xt_post_slr']:.2f} m")

    # --- Head Change Calculation (Example at x=500m) ---
    x_point = 500
    h_pre_slr = calculate_head_at_x(x_point, q0_pre_slr, W, K, Z0, delt)

    # Flux-controlled head change
    h_post_slr_flux = calculate_head_at_x(x_point, q0_pre_slr, W, K, Z0 + DELTA_Z, delt)
    delta_h_flux = (h_post_slr_flux - h_pre_slr) + DELTA_Z
    print(f"\nHead change at x={x_point}m (Flux-Controlled): {delta_h_flux:.4f} m")

    # Head-controlled head change
    h_post_slr_head = calculate_head_at_x(x_point, head_results['q0_post_slr'], W, K, Z0 + DELTA_Z, delt)
    delta_h_head = (h_post_slr_head - h_pre_slr) + DELTA_Z
    print(f"Head change at x={x_point}m (Head-Controlled): {delta_h_head:.4f} m")
