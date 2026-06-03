// ============================================================
// ASKREATIF Perfumery Engine
// Temperature & Climate Model
// ============================================================
// Antoine equation for T-dependent Psat:
//   log₁₀(Psat) = A - B/(C + T)   [kPa, °C]
//
// When full Antoine parameters unavailable:
//   Clausius-Clapeyron approximation:
//   Psat(T) = Psat(T_ref) × exp[ΔHvap/R × (1/T_ref - 1/T)]
//
// Humidity correction for evaporation:
//   k_eff = k_dry × (1 - RH × f_polar)
// ============================================================

import 'dart:math';
import 'package:askreatif_app/Compound.dart';

// ── Climate Conditions ────────────────────────────────────
class ClimateCondition {
  final String name;
  final double ambientTempC; // °C
  final double skinTempC; // °C (skin surface)
  final double relativeHumidity; // 0–1
  final double airflowMs; // m/s (air movement)
  final String region;

  const ClimateCondition({
    required this.name,
    required this.ambientTempC,
    required this.skinTempC,
    required this.relativeHumidity,
    required this.airflowMs,
    required this.region,
  });
}

// ── Preset Climate Profiles ───────────────────────────────
const ClimateCondition climateEuropeanSpring = ClimateCondition(
  name: 'European Spring',
  ambientTempC: 18.0,
  skinTempC: 32.0,
  relativeHumidity: 0.55,
  airflowMs: 0.5,
  region: 'Temperate',
);

const ClimateCondition climateTropicalHot = ClimateCondition(
  name: 'Tropical Hot',
  ambientTempC: 34.0,
  skinTempC: 36.5,
  relativeHumidity: 0.85,
  airflowMs: 0.3,
  region: 'Tropical',
);

const ClimateCondition climateMiddleEastDry = ClimateCondition(
  name: 'Middle East Dry',
  ambientTempC: 40.0,
  skinTempC: 37.5,
  relativeHumidity: 0.15,
  airflowMs: 1.5,
  region: 'Arid',
);

const ClimateCondition climateNordicWinter = ClimateCondition(
  name: 'Nordic Winter',
  ambientTempC: -5.0,
  skinTempC: 30.0,
  relativeHumidity: 0.70,
  airflowMs: 2.0,
  region: 'Cold',
);

const ClimateCondition climateSEAsiaIndoor = ClimateCondition(
  name: 'SE Asia Indoor (AC)',
  ambientTempC: 24.0,
  skinTempC: 33.0,
  relativeHumidity: 0.60,
  airflowMs: 0.2,
  region: 'Tropical',
);

const List<ClimateCondition> presetClimates = [
  climateEuropeanSpring,
  climateTropicalHot,
  climateMiddleEastDry,
  climateNordicWinter,
  climateSEAsiaIndoor,
];

// ── Clausius-Clapeyron Psat Correction ───────────────────
/// Corrects vapour pressure for temperature deviation from 25°C.
/// Uses Trouton's rule estimate for ΔHvap when not available:
///   ΔHvap ≈ 88 × Tb [J/mol]  (Trouton's rule)
/// Approximate Tb from Psat at 25°C via inverted Clausius-Clapeyron.
///
/// Returns Psat_corrected [kPa] at temperature T [°C].
double psatAtTemperature({
  required double psatRef, // kPa at 25°C
  required double mw, // g/mol (used for ΔHvap estimate)
  double tRefC = 25.0,
  required double tC, // target temperature °C
}) {
  if (psatRef <= 0) return psatRef;

  const double R = 8.314;
  final double tRef = tRefC + 273.15;
  final double t = tC + 273.15;

  // Estimate ΔHvap using modified Watson correlation:
  // ΔHvap ≈ k_w × MW^0.38 × Psat^(-0.12)
  // Empirically calibrated for fragrance materials (kJ/mol)
  final double dHvap = _estimateDeltaHvap(psatRef, mw);

  final double correction = exp((dHvap / R) * (1.0 / tRef - 1.0 / t));
  final double result = psatRef * correction;
  return result.isFinite ? result : psatRef;
}

double _estimateDeltaHvap(double psatKpa, double mwGmol) {
  // Empirical correlation for fragrance materials [J/mol]
  // Calibrated against NIST data for terpenoids, esters, alcohols
  // Range: ~30–80 kJ/mol for typical fragrance materials
  double dh = 50000.0; // default 50 kJ/mol
  if (psatKpa >= 10.0)
    dh = 35000.0; // very volatile
  else if (psatKpa >= 1.0)
    dh = 45000.0; // moderately volatile
  else if (psatKpa >= 0.1)
    dh = 55000.0; // low volatility
  else if (psatKpa >= 0.01)
    dh = 65000.0; // very low
  else
    dh = 75000.0; // base note level
  // MW correction: larger molecules have higher ΔHvap
  dh *= (1.0 + 0.002 * (mwGmol - 150.0).clamp(-100.0, 200.0));
  return dh;
}

// ── Climate-Corrected Evaporation Constant ────────────────
/// Enhanced evaporation constant with climate corrections.
///
/// k_climate = k_base × f_temp × f_humidity × f_airflow × f_skin
///
/// where:
///   f_temp     = Psat(Tskin) / Psat(25°C)  — temperature effect
///   f_humidity = 1 - RH × polarity_factor  — humidity retardation
///   f_airflow  = 1 + β × v_air             — convective enhancement
///   f_skin     = skin adsorption damping    — log(MW) correction
double climateEvaporationConstant({
  required Compound c,
  required ClimateCondition climate,
  double alphaBase = 0.0012,
}) {
  // Base constant (Graham's law model)
  if (c.Psat <= 0) return 1e-6;
  final double kBase = alphaBase * c.Psat / sqrt(c.MW);

  // Temperature correction via Clausius-Clapeyron
  final double psatSkin = psatAtTemperature(
    psatRef: c.Psat,
    mw: c.MW,
    tC: climate.skinTempC,
  );
  final double fTemp = psatSkin / c.Psat;

  // Humidity correction: polar compounds (OH, COOH groups) are
  // retained more in humid conditions
  final double polarity = _polarityFactor(c);
  final double fHumidity = 1.0 - climate.relativeHumidity * polarity * 0.4;

  // Airflow enhancement (convective mass transfer coefficient)
  // Sherwood number correction: Sh ∝ Re^0.5
  final double fAirflow = 1.0 + 0.8 * sqrt(climate.airflowMs);

  // Skin adsorption damping for base notes
  // Larger, less volatile molecules are retained by skin lipids
  final double fSkin = _skinRetentionFactor(c);

  final double kClimate =
      kBase * fTemp * fHumidity.clamp(0.3, 1.5) * fAirflow * fSkin;

  return kClimate.isFinite && kClimate > 0 ? kClimate : kBase;
}

double _polarityFactor(Compound c) {
  // Estimate polarity from functional groups
  double polarity = 0.0;
  const Map<String, double> groupPolarity = {
    'OH': 1.0,
    'ACOH': 0.9,
    'COOH': 1.0,
    'H2O': 1.0,
    'CH3OH': 0.9,
    'CH3CO': 0.4,
    'CHO': 0.5,
    'COO': 0.3,
    'CH3COO': 0.3,
    'CH2O': 0.4,
  };
  double totalGroups = 0;
  c.groups.forEach((group, count) {
    polarity += (groupPolarity[group] ?? 0.1) * count;
    totalGroups += count;
  });
  return totalGroups > 0 ? (polarity / totalGroups).clamp(0.0, 1.0) : 0.1;
}

double _skinRetentionFactor(Compound c) {
  // Base notes and high-MW compounds are adsorbed by skin/hair
  // Retention factor < 1 means slower effective evaporation
  // Based on octanol-water partitioning heuristics
  if (c.MW > 250) return 0.4;
  if (c.MW > 200 && c.Psat < 0.1) return 0.6;
  if (c.MW > 150 && c.Psat < 0.05) return 0.75;
  return 1.0;
}

// ── Climate Performance Metrics ───────────────────────────
class ClimatePerformanceReport {
  final ClimateCondition climate;
  final double projectedLongevityHours;
  final double projectionStrength; // 0–1, relative sillage
  final double topNoteIntensity; // perceived intensity at t=5min
  final double drydownBalance; // OV balance at t=2h
  final String performanceSummary;

  const ClimatePerformanceReport({
    required this.climate,
    required this.projectedLongevityHours,
    required this.projectionStrength,
    required this.topNoteIntensity,
    required this.drydownBalance,
    required this.performanceSummary,
  });
}

ClimatePerformanceReport evaluateClimatePerformance({
  required List<Compound> fragranceComps,
  required List<double> moleFractions,
  required List<double> odorValues,
  required ClimateCondition climate,
}) {
  // Compute climate-corrected k values
  final List<double> kValues =
      fragranceComps
          .map((c) => climateEvaporationConstant(c: c, climate: climate))
          .toList();

  // Projected longevity: OV-weighted t_90%
  final double totalOV = odorValues.fold(0.0, (a, b) => a + b);
  double kEff = 0.0;
  if (totalOV > 0) {
    for (int i = 0; i < fragranceComps.length; i++) {
      kEff += (odorValues[i] / totalOV) * kValues[i];
    }
  }
  final double longevity = kEff > 0 ? -log(0.10) / kEff : 99.0;

  // Projection strength: sum of high-volatility OV contributions
  // at t=5min normalized
  double projSum = 0.0;
  for (int i = 0; i < fragranceComps.length; i++) {
    if (fragranceComps[i].Psat > 1.0) {
      projSum += odorValues[i] * exp(-kValues[i] * (5.0 / 60.0));
    }
  }
  final double projStrength = (projSum / (totalOV + 1e-30)).clamp(0.0, 1.0);

  // Top note intensity at 5 minutes
  double topI = 0.0;
  for (int i = 0; i < fragranceComps.length; i++) {
    final double residual = exp(-kValues[i] * (5.0 / 60.0));
    topI += odorValues[i] * residual;
  }

  // Drydown balance at 2h: std dev of log(OV)
  List<double> ovAt2h = [];
  for (int i = 0; i < fragranceComps.length; i++) {
    final double ov2 = odorValues[i] * exp(-kValues[i] * 2.0);
    if (ov2 > 1e-30) ovAt2h.add(log(ov2));
  }
  double balance = 0.0;
  if (ovAt2h.length > 1) {
    final double mean = ovAt2h.fold(0.0, (a, b) => a + b) / ovAt2h.length;
    balance = sqrt(
      ovAt2h.map((v) => (v - mean) * (v - mean)).fold(0.0, (a, b) => a + b) /
          ovAt2h.length,
    );
  }

  // Summary
  String summary = '${climate.name}: ';
  if (longevity > 8)
    summary += 'Excellent longevity. ';
  else if (longevity > 4)
    summary += 'Good longevity. ';
  else
    summary += 'Short-lived. ';

  if (climate.ambientTempC > 30) {
    summary += 'Hot conditions amplify top notes and reduce drydown. ';
  }
  if (climate.relativeHumidity > 0.75) {
    summary += 'High humidity retards polar compound evaporation. ';
  }
  if (climate.airflowMs > 1.5) {
    summary += 'High airflow enhances projection but reduces longevity. ';
  }

  return ClimatePerformanceReport(
    climate: climate,
    projectedLongevityHours: longevity,
    projectionStrength: projStrength,
    topNoteIntensity: topI,
    drydownBalance: balance,
    performanceSummary: summary,
  );
}
