#pragma once
#include "pch.h"
#include "sqlite3.h"
#include <functional>
#include <array>
#include <string>
#include <cmath>
#include <complex>


#if defined(_WIN32) || defined(_WIN64)
#ifdef WATERSTEAMEQUATIONOFSTATE_EXPORTS
#define WATERSTEAMEQUATIONOFSTATE_API __declspec(dllexport)
#else
#define WATERSTEAMEQUATIONOFSTATE_API __declspec(dllimport)
#endif
#else
#define WATERSTEAMEQUATIONOFSTATE_API
#endif

struct IceSublimationCoefficient
{
	int Index;
	double AiBase;
	int AiExponent;
	double BiBase;
	int BiExponent;
};

struct IceComplexCoefficient
{
	int Index;
	double RealBase;
	int RealExponent;
	double ImaginaryBase;
	int ImaginaryExponent;
};

struct IceRealCoefficient
{
	int Index;
	double RealBase;
	int RealExponent;
};

struct Region13Coefficient
{
	int Index;
	int Ii;
	int Ji;
	double NiBase;
	int NiExponent;
};

struct Region25IdealCoefficient
{
	int Index;
	int Ji;
	double NiBase;
	int NiExponent;
};

struct Region25ResidualCoefficient
{
	int Index;
	int Ii;
	int Ji;
	double NiBase;
	int NiExponent;
};

struct Region4B23Coefficient
{
	int Index;
	double NiBase;
	int NiExponent;
};

enum Region
{
	ICE,
	SUBCOOLED_WATER,
	SUPERCRITICAL_WATER_STEAM,
	SUPERHEATED_STEAM,
	SATURATION,
	HIGH_TEMPERATURE_STEAM,
	INVALID
};

class WATERSTEAMEQUATIONOFSTATE_API WaterSteamEquationOfState {
public:
	explicit WaterSteamEquationOfState(const std::string& databasePath);

	static constexpr double NORMAL_PRESSURE = 0.101325;
	static constexpr double SPECIFIC_GAS_CONSTANT = 0.461526;
	static constexpr double CRITICAL_TEMPERATURE = 647.096;
	static constexpr double CRITICAL_DENSITY = 322.0;
	static constexpr double CRITICAL_PRESSURE = 22.064;
	static constexpr double TRIPLE_POINT_SPECIFIC_ENTHALPY = 333.55;
	static constexpr double TRIPLE_POINT_TEMPERATURE = 273.16;
	static constexpr double TRIPLE_POINT_PRESSURE = 0.000611657;
	static constexpr double TRIPLE_POINT_PRESSURE_PASCALS = 611.657;
	static constexpr double UPPER_PRESSURE_LIMIT = 100.0;
	static constexpr double LOWER_TEMPERATURE_LIMIT = 273.15;
	static constexpr double LOWER_PRESSURE_LIMIT = 0.000611657;
	static constexpr double REGION123_TEMPERATURE_BOUNDARY = 623.15;
	static constexpr double REGION123_PRESSURE_BOUNDARY = 100.0;

	static constexpr double ICE_ZERO_POINT_ENTROPY_ABSOLUTE = 189.13;
	static constexpr double ICE_ZERO_POINT_ENTROPY_IAPWS95 = -3327.337565;
	static constexpr double ICE_EXPERIMENTAL_REDUCED_NORMAL_PRESSURE = NORMAL_PRESSURE / TRIPLE_POINT_PRESSURE;

	static constexpr double REGION1_REDUCING_TEMPERATURE = 1386.0;
	static constexpr double REGION1_REDUCING_PRESSURE = 16.53;

	static constexpr double REGION2_UPPER_TEMPERATURE_LIMIT = 1073.15;
	static constexpr double REGION2_REDUCING_TEMPERATURE = 540.0;
	static constexpr double REGION2_REDUCING_PRESSURE = 1.0;

	static constexpr double REGION23_BOUNDARY_REDUCING_TEMPERATURE = 1;
	static constexpr double REGION23_BOUNDARY_REDUCING_PRESSURE = 1;
	static constexpr double REGION23_BOUNDARY_LOWER_TEMPERATURE_LIMIT = 623.15;
	static constexpr double REGION23_BOUNDARY_LOWER_PRESSURE_LIMIT = 16.5292;
	static constexpr double REGION23_BOUNDARY_UPPER_TEMPERATURE_LIMIT = 863.15;
	static constexpr double REGION23_BOUNDARY_UPPER_PRESSURE_LIMIT = 100.0;

	static constexpr double REGION4_UPPER_TEMPERATURE_LIMIT = 647.096;

	static constexpr double REGION5_REDUCING_PRESSURE = 1.0;
	static constexpr double REGION5_REDUCING_TEMPERATURE = 1000.0;
	static constexpr double REGION5_UPPER_PRESSURE_LIMIT = 50.0;
	static constexpr double HIGH_TEMPERATURE_STEAM_LOWER_TEMPERATURE_LIMIT = 1073.15;
	static constexpr double HIGH_TEMPERATURE_STEAM_UPPER_TEMPERATURE_LIMIT = 2273.15;
	static constexpr double HIGH_TEMPERATURE_STEAM_UPPER_PRESSURE_LIMIT = 50.0;

	static constexpr int NUMBER_OF_ICE_SUBLIMATION_COEFFICIENTS = 3;
	static constexpr int NUMBER_OF_ICE_REAL_GIBBS_COEFFICIENTS = 5;
	static constexpr int NUMBER_OF_ICE_COMPLEX_FUNCTION_COEFFICIENTS = 3;
	static constexpr int NUMBER_OF_ICE_COMPLEX_CONSTANTS_T = 2;
	static constexpr int NUMBER_OF_ICE_COMPLEX_CONSTANTS_R = 1;
	static constexpr int NUMBER_OF_REGION1_COEFFICIENTS = 34;
	static constexpr int NUMBER_OF_REGION2_IDEAL_COEFFICIENTS = 9;
	static constexpr int NUMBER_OF_REGION2_RESIDUAL_COEFFICIENTS = 43;
	static constexpr int NUMBER_OF_REGION23_COEFFICIENTS = 5;
	static constexpr int NUMBER_OF_REGION3_COEFFICIENTS = 40;
	static constexpr int NUMBER_OF_REGION4_COEFFICIENTS = 10;
	static constexpr int NUMBER_OF_REGION5_COEFFICIENTS = 6;

	std::array<IceSublimationCoefficient, NUMBER_OF_ICE_SUBLIMATION_COEFFICIENTS> IceSublimationCoefficients;
	std::array<IceComplexCoefficient, NUMBER_OF_ICE_COMPLEX_FUNCTION_COEFFICIENTS> IceComplexFunctionCoefficients;
	std::array<IceRealCoefficient, NUMBER_OF_ICE_REAL_GIBBS_COEFFICIENTS> IceRealResidualGibbsCoefficients;
	std::array<IceComplexCoefficient, NUMBER_OF_ICE_COMPLEX_CONSTANTS_T> IceComplexConstantsT;
	std::array<IceComplexCoefficient, NUMBER_OF_ICE_COMPLEX_CONSTANTS_R> IceComplexConstantsR;

	[[nodiscard]] Region DetermineRegion(double temperature, double pressure) const;

	template <typename IceFunc, typename RegionFunc1, typename RegionFunc2, typename RegionFunc3, typename RegionFunc5>
	void CalculatePropertyArray(
		const double* temperatures,
		const double* pressures,
		double* results,
		const size_t length,
		IceFunc iceFunc,
		RegionFunc1 region1Func,
		RegionFunc2 region2Func,
		RegionFunc3 region3Func,
		RegionFunc5 region5Func) const
	{
		for (size_t i = 0; i < length; i++) {
			const double temperature = temperatures[i];

			switch (const double pressure = pressures[i]; DetermineRegion(temperature, pressure)) {
			case ICE:
				results[i] = (this->*iceFunc)(temperature, pressure);
				break;
			case SUBCOOLED_WATER:
				results[i] = (this->*region1Func)(temperature, pressure);
				break;
			case SUPERCRITICAL_WATER_STEAM:
				results[i] = (this->*region2Func)(temperature, pressure);
				break;
			case SUPERHEATED_STEAM: {
				results[i] = (this->*region3Func)(temperature, pressure);
				break;
			}
			case HIGH_TEMPERATURE_STEAM:
				results[i] = (this->*region5Func)(temperature, pressure);
				break;
			case SATURATION:
			case INVALID:
				results[i] = std::numeric_limits<double>::quiet_NaN();
				break;
			}
		}
	}

	void CalculateDensityArray(
		const double* temperatures,
		const double* pressures,
		double* densities,
		size_t length) const;

	void CalculateSpecificVolumeArray(
		const double* temperatures,
		const double* pressures,
		double* specificVolumes,
		size_t length) const;

	void CalculateSpecificInternalEnergyArray(
		const double* temperatures,
		const double* pressures,
		double* specificInternalEnergies,
		size_t length) const;

	void CalculateSpecificEnthalpyArray(
		const double* temperatures,
		const double* pressures,
		double* enthalpies,
		size_t length) const;

	void CalculateSpecificEntropyArray(
		const double* temperatures,
		const double* pressures,
		double* entropies,
		size_t length) const;

	void CalculateSpecificIsobaricHeatCapacityArray(
		const double* temperatures,
		const double* pressures,
		double* specificHeatCapacities,
		size_t length) const;

	void CalculateSpecificIsochoricHeatCapacityArray(
		const double* temperatures,
		const double* pressures,
		double* specificHeatCapacities,
		size_t length) const;

	void CalculateSpeedOfSoundArray(
		const double* temperatures,
		const double* pressures,
		double* speedsOfSound,
		size_t length) const;
	double CalculateIceSpeedOfSound(double temperature, double pressure) const;

#pragma region ICE REGION METHOD SIGNATURES
	double CalculateIceSublimationPressure(double temperature) const;
	static double CalculateIceReducedPressure(double pressure);
	static double CalculateIceReducedTemperature(double temperature);
	double CalculateIceResidualGibbsEnergy(double pressure) const;
	double GetIceResidualGibbsCoefficient(int index) const;
	double CalculateFirstDerivativeIceResidualGibbsEnergyPressure(double pressure) const;
	double CalculateSecondDerivativeIceResidualGibbsEnergyPressure(double pressure) const;
	auto CalculateIceResidualFunctionCoefficient(double pressure) const -> std::complex<double>;
	std::complex<double> CalculateFirstDerivativeIceResidualFunctionCoefficientPressure(double pressure) const;
	std::complex<double> CalculateSecondDerivativeIceResidualFunctionCoefficientPressure() const;
	double CalculateIceSpecificGibbsFreeEnergy(double temperature, double pressure) const;
	double CalculateIceSpecificVolume(const double temperature, const double pressure) const;
	double CalculateIceDensity(const double temperature, const double pressure) const;
	double CalculateIceSpecificEntropy(double temperature, double pressure) const;
	double CalculateIceSpecificIsobaricHeatCapacity(double temperature, double pressure) const;
	double CalculateIceSpecificIsochoricHeatCapacity(double temperature, double pressure) const;
	double CalculateIceSpecificEnthalpy(double temperature, double pressure) const;
	double CalculateIceSpecificInternalEnergy(double temperature, double pressure) const;
	double CalculateIceSpecificHelmholtzEnergy(double temperature, double pressure) const;
	double CalculateIceCubicExpansionCoefficient(double temperature, double pressure) const;
	double CalculateIcePressureCoefficient(double temperature, double pressure) const;
	double CalculateIceIsothermalCompressibility(double temperature, double pressure) const;
	double CalculateIceIsentropicCompressibility(double temperature, double pressure) const;
	std::complex<double> GetIceComplexConstantT(int index) const;
	std::complex<double> GetIceComplexFunctionConstantR(int index) const;
	std::complex<double> GetIceComplexCoefficientR(int index) const;
	double CalculateIceFirstDerivativeGibbsFreeEnergyTemperature(double temperature, double pressure) const;
	double CalculateIceFirstDerivativeGibbsFreeEnergyPressure(double temperature, double pressure) const;
	double CalculateIceSecondDerivativeGibbsFreeEnergyPressure(double temperature, double pressure) const;
	double CalculateIceSecondDerivativeGibbsFreeEnergyTemperature(double temperature, double pressure) const;
	double CalculateIceSecondMixedDerivativeGibbsFreeEnergyTemperaturePressure(const double temperature, const double pressure) const;
#pragma endregion

#pragma region REGION 1 METHOD SIGNATURES
	[[nodiscard]] double CalculateRegion1SpecificVolume(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpecificInternalEnergy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpecificGibbsFreeEnergy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpecificEntropy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpecificIsobaricHeatCapacity(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpecificIsochoricHeatCapacity(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpecificEnthalpy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion1SpeedOfSound(double temperature, double pressure) const;
	static double CalculateRegion1ReducedTemperature(double temperature);
	static double CalculateRegion1ReducedPressure(double pressure);
	static double CalculateRegion1ReciprocalReducedTemperature(double temperature);
	static double CalculateRegion1ReciprocalReducedPressure(double pressure);
	[[nodiscard]] double CalculateDimensionlessGibbsFreeEnergyRegion1(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondMixedDerivativeDimensionlessGibbsFreeEnergyRegion1(double temperature, double pressure) const;
	auto CalculateRegion1Density(double temperature, double pressure) const -> double;
#pragma endregion

#pragma region REGION 2 METHOD SIGNATURES
	auto CalculateRegion2Density(double temperature, double pressure) const -> double;
	[[nodiscard]] double CalculateRegion2SpecificVolume(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpecificInternalEnergy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpecificGibbsFreeEnergy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpecificEntropy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpecificEnthalpy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpecificIsochoricHeatCapacity(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpecificIsobaricHeatCapacity(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion2SpeedOfSound(double temperature, double pressure) const;
	static double CalculateRegion2ReducedPressure(double pressure);
	static double CalculateRegion2ReducedTemperature(double temperature);
	static double CalculateRegion2ReciprocalReducedTemperature(double temperature);
	static double CalculateRegion2ReciprocalReducedPressure(double pressure);
	[[nodiscard]] double CalculateDimensionlessGibbsFreeEnergyRegion2(double temperature, double pressure) const;
	[[nodiscard]] double CalculateDimensionlessIdealGibbsFreeEnergyRegion2(double temperature, double pressure) const;
	static double CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion2(const double pressure);
	static double CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion2(const double pressure);
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(const double temperature) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(const double temperature) const;
	[[nodiscard]] double CalculateDimensionlessResidualGibbsFreeEnergyRegion2(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion2(double temperature, double pressure) const;
#pragma endregion

#pragma region REGION BOUNDARY 2_3 METHOD SIGNATURES
	static double CalculateRegion23ReducedPressure(double pressure);
	static double CalculateRegion23ReducedTemperature(double temperature);
	static double CalculateRegion23ReciprocalReducedTemperature(double temperature);
	static double CalculateRegion23ReciprocalReducedPressure(double pressure);
	[[nodiscard]] double CalculateRegion23BoundaryTemperature(double pressure) const;
	[[nodiscard]] double CalculateRegion23BoundaryPressure(double temperature) const;
	auto CalculateRegion3SpecificVolume(double temperature, double pressure) const -> double;
	auto CalculateRegion3Density(double temperature, double pressure) const -> double;
#pragma endregion

#pragma region REGION 3 METHOD SIGNATURES
	static double CalculateRegion3ReducedTemperature(double temperature);
	static double CalculateRegion3ReciprocalReducedTemperature(double temperature);
	static double CalculateRegion3ReducedDensity(double density);
	[[nodiscard]] double CalculateRegion3SpecificHelmholtzFreeEnergy(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3DimensionlessHelmholtzEnergy(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyTemperature(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyTemperature(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyDensity(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyDensity(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3SecondMixedDerivativeDimensionlessHelmholtzEnergyDensityTemperature(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3Pressure(double temperature, double density) const;
	[[nodiscard]] double CalculateRegion3SpecificInternalEnergy(double temperature, double density) const;
	auto CalculateRegion3SpecificInternalEnergyViaPressure(double temperature, double pressure) const -> double;
	[[nodiscard]] double CalculateRegion3SpecificEntropy(double temperature, double density) const;
	auto CalculateRegion3SpecificEntropyViaPressure(double temperature, double pressure) const -> double;
	[[nodiscard]] double CalculateRegion3SpecificEnthalpy(double temperature, double density) const;
	auto CalculateRegion3SpecificEnthalpyViaPressure(double temperature, double pressure) const -> double;
	[[nodiscard]] double CalculateRegion3SpecificIsochoricHeatCapacity(double temperature, double density) const;
	auto CalculateRegion3SpecificIsochoricHeatCapacityViaPressure(double temperature, double pressure) const -> double;
	[[nodiscard]] double CalculateRegion3SpecificIsobaricHeatCapacity(double temperature, double density) const;
	auto CalculateRegion3SpecificIsobaricHeatCapacityViaPressure(double temperature, double pressure) const -> double;
	[[nodiscard]] double CalculateRegion3SpeedOfSound(double temperature, double density) const;
	auto CalculateRegion3SpeedOfSoundViaPressure(double temperature, double pressure) const -> double;
#pragma endregion

#pragma region REGION 4 METHOD SIGNATURES
	[[nodiscard]] double CalculateRegion4SaturationPressure(double temperature) const;
	[[nodiscard]] double CalculateRegion4SaturationTemperature(double pressure) const;
#pragma endregion

#pragma region REGION 5 METHOD SIGNATURES
	static double CalculateRegion5ReciprocalReducedTemperature(double temperature);
	static double CalculateRegion5ReducedPressure(double pressure);
	double CalculateRegion5Density(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificVolume(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificInternalEnergy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificEntropy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificEnthalpy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificIsochoricHeatCapacity(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificIsobaricHeatCapacity(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpeedOfSound(double temperature, double pressure) const;
	[[nodiscard]] double CalculateRegion5SpecificGibbsFreeEnergy(double temperature, double pressure) const;
	[[nodiscard]] double CalculateDimensionlessGibbsFreeEnergyRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateDimensionlessIdealGibbsFreeEnergyRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(double temperature) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(double temperature) const;
	[[nodiscard]] double CalculateDimensionlessResidualGibbsFreeEnergyRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(double temperature, double pressure) const;
	[[nodiscard]] double CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion5(double temperature, double pressure) const;
	static double CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion5(double pressure);
	static double CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion5(double pressure);
#pragma endregion

private:
	std::string DatabasePath;
	[[nodiscard]] std::string GetDatabasePath() const;
	void SetDatabasePath(const std::string& path);

	bool GetCoefficients();
	template<typename CoefficientType, size_t N>
	void LoadCoefficients(sqlite3* db, const std::string& tableName, std::array<CoefficientType, N>& coefficients, std::function<CoefficientType(sqlite3_stmt*)> rowMapper)
	{
		sqlite3_stmt* stmt;
		const std::string sql = "SELECT * FROM " + tableName;

		int rc = sqlite3_prepare_v2(db, sql.c_str(), -1, &stmt, nullptr);
		if (rc != SQLITE_OK) {
			throw std::runtime_error("Failed to prepare SQL statement for table " + tableName + ": " + sqlite3_errmsg(db));
		}

		while ((rc = sqlite3_step(stmt)) == SQLITE_ROW) {
			CoefficientType coefficient = rowMapper(stmt);
			const int index = coefficient.Index;
			if (index <= 0 || index > N) {
				throw std::runtime_error("Index out of bounds in table: " + tableName);
			}
			coefficients[index - 1] = coefficient;
		}

		if (rc != SQLITE_DONE) {
			throw std::runtime_error("Error reading from table: " + tableName + " - " + sqlite3_errmsg(db));
		}
	}

	void LoadRegion1Coefficients(sqlite3* db);
	void LoadRegion2IdealCoefficients(sqlite3* db);
	void LoadRegion2ResidualCoefficients(sqlite3* db);
	void LoadRegionB23Coefficients(sqlite3* db);
	void LoadRegion3Coefficients(sqlite3* db);
	void LoadRegion4Coefficients(sqlite3* db);
	void LoadRegion5IdealCoefficients(sqlite3* db);
	void LoadRegion5ResidualCoefficients(sqlite3* db);
	void LoadIceRealGibbsCoefficients(sqlite3* db);
	void LoadIceComplexFunctionCoefficients(sqlite3* db);
	void LoadIceComplexConstantsT(sqlite3* db);
	void LoadIceComplexConstantsR(sqlite3* db);
	void LoadIceSublimationCoefficients(sqlite3* db);

	static constexpr auto ICE_SUBLIMATION_COEFFICIENTS_TABLE_NAME = "ice_sublimation_coefficients";
	static constexpr auto ICE_COMPLEX_CONSTANT_R_TABLE_NAME = "ice_complex_constant_r";
	static constexpr auto ICE_COMPLEX_CONSTANTS_T_TABLE_NAME = "ice_complex_constants_t";
	static constexpr auto ICE_COMPLEX_FUNCTION_COEFFICIENTS = "ice_complex_function_coefficients";
	static constexpr auto ICE_REAL_GIBBS_CONSTANTS = "ice_real_gibbs_constants";
	static constexpr auto REGION_1_TABLE_NAME = "water_steam_coefficients_region1";
	static constexpr auto REGION_2_IDEAL_TABLE_NAME = "water_steam_ideal_coefficients_region2";
	static constexpr auto REGION_2_RESIDUAL_TABLE_NAME = "water_steam_residual_coefficients_region2";
	static constexpr auto REGION_B23_TABLE_NAME = "water_steam_b23_coefficients";
	static constexpr auto REGION_3_TABLE_NAME = "water_steam_coefficients_region3";
	static constexpr auto REGION_4_TABLE_NAME = "water_steam_coefficients_region4";
	static constexpr auto REGION_5_IDEAL_TABLE_NAME = "water_steam_ideal_coefficients_region5";
	static constexpr auto REGION_5_RESIDUAL_TABLE_NAME = "water_steam_residual_coefficients_region5";

	std::array<Region13Coefficient, NUMBER_OF_REGION1_COEFFICIENTS> Region1Coefficients;
	std::array<Region25IdealCoefficient, NUMBER_OF_REGION2_IDEAL_COEFFICIENTS> Region2IdealCoefficients;
	std::array<Region25ResidualCoefficient, NUMBER_OF_REGION2_RESIDUAL_COEFFICIENTS> Region2ResidualCoefficients;
	std::array<Region4B23Coefficient, NUMBER_OF_REGION23_COEFFICIENTS> Region23Coefficients;
	std::array<Region13Coefficient, NUMBER_OF_REGION3_COEFFICIENTS> Region3Coefficients;
	std::array<Region4B23Coefficient, NUMBER_OF_REGION4_COEFFICIENTS> Region4Coefficients;
	std::array<Region25IdealCoefficient, NUMBER_OF_REGION5_COEFFICIENTS> Region5IdealCoefficients;
	std::array<Region25ResidualCoefficient, NUMBER_OF_REGION5_COEFFICIENTS> Region5ResidualCoefficients;
};

template <typename Func>
int CalculatePropertyArray(
	const WaterSteamEquationOfState* instance,
	const double* temperatures,
	const double* pressures,
	double* results,
	const size_t length,
	Func func)
{
	if (!instance)
	{
		return -1; // Invalid instance
	}

	try
	{
		(instance->*func)(temperatures, pressures, results, length);
		return 0; // Success
	}
	catch (...)
	{
		return -1; // Error
	}
}

extern "C" {
	WATERSTEAMEQUATIONOFSTATE_API WaterSteamEquationOfState* CreateWaterSteamEquationOfState(const char* databasePath);
	WATERSTEAMEQUATIONOFSTATE_API void DestroyWaterSteamEquationOfState(const WaterSteamEquationOfState* instance);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateDensityArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* densities,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpecificVolumeArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* specificVolumes,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpecificInternalEnergyArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* specificInternalEnergies,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpecificEnthalpyArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* enthalpies,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpecificEntropyArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* entropies,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpecificIsobaricHeatCapacityArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* specificHeatCapacities,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpecificIsochoricHeatCapacityArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* specificHeatCapacities,
		size_t length);
	WATERSTEAMEQUATIONOFSTATE_API int CalculateSpeedOfSoundArray(
		const WaterSteamEquationOfState* instance,
		const double* temperatures,
		const double* pressures,
		double* speedsOfSound,
		size_t length);
}