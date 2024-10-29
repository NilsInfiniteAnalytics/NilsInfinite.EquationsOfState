#include "pch.h"
#include "WaterSteamEquationOfState.h"
extern "C" {
#include "sqlite3.h"
}

WaterSteamEquationOfState::WaterSteamEquationOfState(const std::string& databasePath) :
	Region1Coefficients{},
	Region2IdealCoefficients{},
	Region2ResidualCoefficients{},
	Region23Coefficients{},
	Region4Coefficients{},
	Region5IdealCoefficients{},
	Region5ResidualCoefficients{}
{
	if (databasePath.empty())
	{
		throw std::invalid_argument("Database path cannot be empty.");
	}
	SetDatabasePath(databasePath);
	const bool success = GetCoefficients();
	if (!success) {
		throw std::runtime_error("Failed to load coefficients from the database.");
	}
}

#pragma region Setup Methods
std::string WaterSteamEquationOfState::GetDatabasePath() const {
	return DatabasePath;
}

void WaterSteamEquationOfState::SetDatabasePath(const std::string& path) {
	DatabasePath = path;
}

bool WaterSteamEquationOfState::GetCoefficients()
{
	sqlite3* db = nullptr;
	const int rc = sqlite3_open(DatabasePath.c_str(), &db);
	if (rc != SQLITE_OK) {
		throw std::runtime_error("Can't open database: " + std::string(sqlite3_errmsg(db)));
	}
	try {
		LoadRegion1Coefficients(db);
		LoadRegion2IdealCoefficients(db);
		LoadRegion2ResidualCoefficients(db);
		LoadRegionB23Coefficients(db);
		LoadRegion3Coefficients(db);
		LoadRegion4Coefficients(db);
		LoadRegion5IdealCoefficients(db);
		LoadRegion5ResidualCoefficients(db);
		return true;
	}
	catch (...) {
		sqlite3_close_v2(db);
		return false;
	}
}

void WaterSteamEquationOfState::LoadRegion1Coefficients(sqlite3* db)
{
	this->LoadCoefficients<Region13Coefficient, NUMBER_OF_REGION1_COEFFICIENTS>(
		db,
		REGION_1_TABLE_NAME,
		Region1Coefficients,
		[](sqlite3_stmt* stmt) -> Region13Coefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const int ii = sqlite3_column_int(stmt, 1);
			const int ji = sqlite3_column_int(stmt, 2);
			const double niBase = sqlite3_column_double(stmt, 3);
			const int niExponent = sqlite3_column_int(stmt, 4);
			return { index, ii, ji, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegion2IdealCoefficients(sqlite3* db)
{
	this->LoadCoefficients<Region25IdealCoefficient, NUMBER_OF_REGION2_IDEAL_COEFFICIENTS>(
		db,
		REGION_2_IDEAL_TABLE_NAME,
		Region2IdealCoefficients,
		[](sqlite3_stmt* stmt) -> Region25IdealCoefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const int ji = sqlite3_column_int(stmt, 1);
			const double niBase = sqlite3_column_double(stmt, 2);
			const int niExponent = sqlite3_column_int(stmt, 3);
			return Region25IdealCoefficient{ index, ji, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegion2ResidualCoefficients(sqlite3* db)
{
	this->LoadCoefficients<Region25ResidualCoefficient, NUMBER_OF_REGION2_RESIDUAL_COEFFICIENTS>(
		db,
		REGION_2_RESIDUAL_TABLE_NAME,
		Region2ResidualCoefficients,
		[](sqlite3_stmt* stmt) -> Region25ResidualCoefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const int ii = sqlite3_column_int(stmt, 1);
			const int ji = sqlite3_column_int(stmt, 2);
			const double niBase = sqlite3_column_double(stmt, 3);
			const int niExponent = sqlite3_column_int(stmt, 4);
			return Region25ResidualCoefficient{ index, ii, ji, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegionB23Coefficients(sqlite3* db)
{
	this->LoadCoefficients<Region4B23Coefficient, NUMBER_OF_REGION23_COEFFICIENTS>(
		db,
		REGION_B23_TABLE_NAME,
		Region23Coefficients,
		[](sqlite3_stmt* stmt) -> Region4B23Coefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const double niBase = sqlite3_column_double(stmt, 1);
			const int niExponent = sqlite3_column_int(stmt, 2);
			return Region4B23Coefficient{ index, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegion3Coefficients(sqlite3* db)
{
	this->LoadCoefficients<Region13Coefficient, NUMBER_OF_REGION3_COEFFICIENTS>(
		db,
		REGION_3_TABLE_NAME,
		Region3Coefficients,
		[](sqlite3_stmt* stmt) -> Region13Coefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const int ii = sqlite3_column_int(stmt, 1);
			const int ji = sqlite3_column_int(stmt, 2);
			const double niBase = sqlite3_column_double(stmt, 3);
			const int niExponent = sqlite3_column_int(stmt, 4);
			return Region13Coefficient{ index, ii, ji, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegion4Coefficients(sqlite3* db)
{
	this->LoadCoefficients<Region4B23Coefficient, NUMBER_OF_REGION4_COEFFICIENTS>(
		db,
		REGION_4_TABLE_NAME,
		Region4Coefficients,
		[](sqlite3_stmt* stmt) -> Region4B23Coefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const double niBase = sqlite3_column_double(stmt, 1);
			const int niExponent = sqlite3_column_int(stmt, 2);
			return Region4B23Coefficient{ index, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegion5IdealCoefficients(sqlite3* db)
{
	this->LoadCoefficients<Region25IdealCoefficient, NUMBER_OF_REGION5_COEFFICIENTS>(
		db,
		REGION_5_IDEAL_TABLE_NAME,
		Region5IdealCoefficients,
		[](sqlite3_stmt* stmt) -> Region25IdealCoefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const int ji = sqlite3_column_int(stmt, 1);
			const double niBase = sqlite3_column_double(stmt, 2);
			const int niExponent = sqlite3_column_int(stmt, 3);
			return Region25IdealCoefficient{ index, ji, niBase, niExponent };
		});
}

void WaterSteamEquationOfState::LoadRegion5ResidualCoefficients(sqlite3* db)
{
	this->LoadCoefficients<Region25ResidualCoefficient, NUMBER_OF_REGION5_COEFFICIENTS>(
		db,
		REGION_5_RESIDUAL_TABLE_NAME,
		Region5ResidualCoefficients,
		[](sqlite3_stmt* stmt) -> Region25ResidualCoefficient {
			const int index = sqlite3_column_int(stmt, 0);
			const int ii = sqlite3_column_int(stmt, 1);
			const int ji = sqlite3_column_int(stmt, 2);
			const double niBase = sqlite3_column_double(stmt, 3);
			const int niExponent = sqlite3_column_int(stmt, 4);
			return Region25ResidualCoefficient{ index, ii, ji, niBase, niExponent };
		});
}
#pragma endregion

Region WaterSteamEquationOfState::DetermineRegion(const double temperature, const double pressure) const
{
	if (temperature > HIGH_TEMPERATURE_STEAM_UPPER_TEMPERATURE_LIMIT)
	{
		return INVALID;
	}

	if (pressure > REGION123_PRESSURE_BOUNDARY)
	{
		return INVALID;
	}

	if (temperature <= REGION123_TEMPERATURE_BOUNDARY)
	{
		const double saturationPressure = CalculateRegion4SaturationPressure(temperature);
		if (pressure > saturationPressure)
		{
			return SUBCOOLED_WATER;
		}
		if (abs(pressure - saturationPressure) < 1e-5)
		{
			return SATURATION;
		}
		return SUPERCRITICAL_WATER_STEAM;
	}

	if (temperature <= REGION2_UPPER_TEMPERATURE_LIMIT)
	{
		if (temperature <= REGION4_UPPER_TEMPERATURE_LIMIT)
		{
			const double saturationPressure = CalculateRegion4SaturationPressure(temperature);
			if (pressure > saturationPressure)
			{
				return SUPERHEATED_STEAM;
			}
			if (abs(pressure - saturationPressure) < 1e-5)
			{
				return SATURATION;
			}
			return SUPERCRITICAL_WATER_STEAM;
		}
		const double boundaryPressure = CalculateRegion23BoundaryPressure(temperature);
		if (pressure > boundaryPressure)
		{
			return SUPERHEATED_STEAM;
		}
		return SUPERCRITICAL_WATER_STEAM;
	}

	if (temperature <= HIGH_TEMPERATURE_STEAM_LOWER_TEMPERATURE_LIMIT)
	{
		if (pressure > HIGH_TEMPERATURE_STEAM_UPPER_PRESSURE_LIMIT)
		{
			throw std::invalid_argument("Pressure is above the upper limit for high temperature steam.");
		}
		return HIGH_TEMPERATURE_STEAM;
	}
	return INVALID;
}

#pragma region Region 1 Equations

double WaterSteamEquationOfState::CalculateRegion1ReciprocalReducedPressure(const double pressure)
{
	if (pressure <= 0.0) {
		throw std::invalid_argument("Pressure must be greater than zero.");
	}
	return REGION1_REDUCING_PRESSURE / pressure;
}

double WaterSteamEquationOfState::CalculateRegion1ReciprocalReducedTemperature(const double temperature)
{
	if (temperature <= 0.0) {
		throw std::invalid_argument("Temperature must be greater than zero.");
	}
	return REGION1_REDUCING_TEMPERATURE / temperature;
}

double WaterSteamEquationOfState::CalculateRegion1ReducedPressure(const double pressure)
{
	return pressure / REGION1_REDUCING_PRESSURE;
}

double WaterSteamEquationOfState::CalculateRegion1ReducedTemperature(const double temperature)
{
	return temperature / REGION1_REDUCING_TEMPERATURE;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificVolume(const double temperature, const double pressure) const
{
	const double pressureKiloPascals = pressure * 1000;
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double firstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(temperature, pressure);
	const double specificVolume = SPECIFIC_GAS_CONSTANT * temperature * reducedPressure * firstDerivativeGibbsPressure / pressureKiloPascals;
	return specificVolume;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificInternalEnergy(const double temperature, const double pressure) const
{
	const double firstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(temperature, pressure);
	const double firstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(temperature, pressure);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double specificInternalEnergy = SPECIFIC_GAS_CONSTANT * temperature * (reciprocalReducedTemperature * firstDerivativeGibbsTemperature - reducedPressure * firstDerivativeGibbsPressure);
	return specificInternalEnergy;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificGibbsFreeEnergy(const double temperature, const double pressure) const
{
	const double gibbsFreeEnergy = CalculateDimensionlessGibbsFreeEnergyRegion1(temperature, pressure);
	const double specificGibbsFreeEnergy = SPECIFIC_GAS_CONSTANT * temperature * gibbsFreeEnergy;
	return specificGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificEntropy(const double temperature, const double pressure) const
{
	const double dimensionlessGibbsEnergy = CalculateDimensionlessGibbsFreeEnergyRegion1(temperature, pressure);
	const double firstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(temperature, pressure);
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double specificEntropy = SPECIFIC_GAS_CONSTANT * (reciprocalReducedTemperature * firstDerivativeGibbsTemperature - dimensionlessGibbsEnergy);
	return specificEntropy;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificIsobaricHeatCapacity(const double temperature, const double pressure) const
{
	const double secondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(temperature, pressure);
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double specificIsobaricHeatCapacity = -1.0 * reciprocalReducedTemperature * reciprocalReducedTemperature * SPECIFIC_GAS_CONSTANT * secondDerivativeGibbsTemperature;
	return specificIsobaricHeatCapacity;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificIsochoricHeatCapacity(const double temperature, const double pressure) const
{
	const double firstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(temperature, pressure);
	const double secondDerivativeGibbsPressure = CalculateSecondDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(temperature, pressure);
	const double secondMixedDerivativeGibbs = CalculateSecondMixedDerivativeDimensionlessGibbsFreeEnergyRegion1(temperature, pressure);
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double specificIsobaricHeatCapacity = CalculateRegion1SpecificIsobaricHeatCapacity(temperature, pressure);
	const double specificIsochoricHeatCapacity = specificIsobaricHeatCapacity
		+ SPECIFIC_GAS_CONSTANT
		* (pow(firstDerivativeGibbsPressure - reciprocalReducedTemperature * secondMixedDerivativeGibbs, 2.0) / secondDerivativeGibbsPressure);
	return specificIsochoricHeatCapacity;
}

double WaterSteamEquationOfState::CalculateRegion1SpecificEnthalpy(const double temperature, const double pressure) const
{
	const double firstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(temperature, pressure);
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double specificEnthalpy = SPECIFIC_GAS_CONSTANT * temperature * reciprocalReducedTemperature * firstDerivativeGibbsTemperature;
	return specificEnthalpy;
}

double WaterSteamEquationOfState::CalculateRegion1SpeedOfSound(const double temperature, const double pressure) const
{
	const double firstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(temperature, pressure);
	const double secondDerivativeGibbsPressure = CalculateSecondDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(temperature, pressure);
	const double secondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(temperature, pressure);
	const double secondMixedDerivativeGibbs = CalculateSecondMixedDerivativeDimensionlessGibbsFreeEnergyRegion1(temperature, pressure);
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double numerator = pow(firstDerivativeGibbsPressure, 2.0) * pow(reciprocalReducedTemperature, 2.0) * secondDerivativeGibbsTemperature;
	const double denTerm1 = pow(firstDerivativeGibbsPressure - reciprocalReducedTemperature * secondMixedDerivativeGibbs, 2.0);
	const double denTerm2 = secondDerivativeGibbsPressure * pow(reciprocalReducedTemperature, 2.0) * secondDerivativeGibbsTemperature;
	const double denominator = denTerm1 - denTerm2;
	const double speedOfSound = sqrt(SPECIFIC_GAS_CONSTANT * temperature * (numerator / denominator) * 1000);
	return speedOfSound;
}

double WaterSteamEquationOfState::CalculateDimensionlessGibbsFreeEnergyRegion1(const double temperature, const double pressure) const
{
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double iTerm = 7.1 - reducedPressure;
	const double jTerm = reciprocalReducedTemperature - 1.222;
	double dimensionlessGibbsFreeEnergy = 0.0;
	for (const auto& coefficient : Region1Coefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(iTerm, coefficient.Ii * 1.0)
			* pow(jTerm, coefficient.Ji * 1.0);
		dimensionlessGibbsFreeEnergy += term;
	}
	return dimensionlessGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(const double temperature, const double pressure) const
{
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double iTerm = 7.1 - reducedPressure;
	const double jTerm = reciprocalReducedTemperature - 1.222;
	double dimensionlessGibbsFreeEnergy = 0.0;
	for (const auto& coefficient : Region1Coefficients) {
		const double iTermExp = coefficient.Ii * 1.0 - 1.0;
		const double term = -1.0
			* coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* (coefficient.Ii * 1.0)
			* pow(iTerm, iTermExp)
			* pow(jTerm, coefficient.Ji * 1.0);
		dimensionlessGibbsFreeEnergy += term;
	}
	return dimensionlessGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessGibbsFreeEnergyPressureRegion1(const double temperature, const double pressure) const
{
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double iTerm = 7.1 - reducedPressure;
	const double jTerm = reciprocalReducedTemperature - 1.222;
	double dimensionlessGibbsFreeEnergy = 0.0;
	for (const auto& coefficient : Region1Coefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* (coefficient.Ii * 1.0)
			* (coefficient.Ii * 1.0 - 1)
			* pow(iTerm, coefficient.Ii * 1.0 - 2.0)
			* pow(jTerm, coefficient.Ji * 1.0);
		dimensionlessGibbsFreeEnergy += term;
	}
	return dimensionlessGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(const double temperature, const double pressure) const
{
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double iTerm = 7.1 - reducedPressure;
	const double jTerm = reciprocalReducedTemperature - 1.222;
	double dimensionlessGibbsFreeEnergy = 0.0;
	for (const auto& coefficient : Region1Coefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(iTerm, coefficient.Ii * 1.0)
			* (coefficient.Ji * 1.0)
			* pow(jTerm, coefficient.Ji * 1.0 - 1.0);
		dimensionlessGibbsFreeEnergy += term;
	}
	return dimensionlessGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessGibbsFreeEnergyTemperatureRegion1(const double temperature, const double pressure) const
{
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double iTerm = 7.1 - reducedPressure;
	const double jTerm = reciprocalReducedTemperature - 1.222;
	double dimensionlessGibbsFreeEnergy = 0.0;
	for (const auto& coefficient : Region1Coefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(iTerm, coefficient.Ii * 1.0)
			* (coefficient.Ji * 1.0)
			* (coefficient.Ji * 1.0 - 1.0)
			* pow(jTerm, coefficient.Ji * 1.0 - 2.0);
		dimensionlessGibbsFreeEnergy += term;
	}
	return dimensionlessGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateSecondMixedDerivativeDimensionlessGibbsFreeEnergyRegion1(const double temperature, const double pressure) const
{
	const double reciprocalReducedTemperature = CalculateRegion1ReciprocalReducedTemperature(temperature);
	const double reducedPressure = CalculateRegion1ReducedPressure(pressure);
	const double iTerm = 7.1 - reducedPressure;
	const double jTerm = reciprocalReducedTemperature - 1.222;
	double dimensionlessGibbsFreeEnergy = 0.0;
	for (const auto& coefficient : Region1Coefficients) {
		const double term = -1.0
			* coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* (coefficient.Ii * 1.0)
			* pow(iTerm, coefficient.Ii * 1.0 - 1.0)
			* (coefficient.Ji * 1.0)
			* pow(jTerm, coefficient.Ji * 1.0 - 1.0);
		dimensionlessGibbsFreeEnergy += term;
	}
	return dimensionlessGibbsFreeEnergy;
}

#pragma endregion

#pragma region Region 2 Equations

double WaterSteamEquationOfState::CalculateRegion2SpecificVolume(const double temperature, const double pressure) const
{
	const double pressureKiloPascals = pressure * 1000;
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	const double idealFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion2(pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(temperature, pressure);
	const double dimensionlessPressure = region2ReducedPressure * (idealFirstDerivativeGibbsPressure + residualFirstDerivativeGibbsPressure);
	const double specificVolume = SPECIFIC_GAS_CONSTANT * temperature * dimensionlessPressure / pressureKiloPascals;
	return specificVolume;
}

double WaterSteamEquationOfState::CalculateRegion2SpecificInternalEnergy(const double temperature, const double pressure) const
{
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double idealFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion2(pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(temperature, pressure);
	const double idealFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(temperature);
	const double residualFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(temperature, pressure);
	const double dimensionlessEnergy = region2ReciprocalReducedTemperature * (idealFirstDerivativeGibbsTemperature + residualFirstDerivativeGibbsTemperature)
		- region2ReducedPressure * (idealFirstDerivativeGibbsPressure + residualFirstDerivativeGibbsPressure);
	const double specificInternalEnergy = SPECIFIC_GAS_CONSTANT * temperature * dimensionlessEnergy;
	return specificInternalEnergy;
}

double WaterSteamEquationOfState::CalculateRegion2SpecificGibbsFreeEnergy(const double temperature, const double pressure) const
{
	return CalculateDimensionlessGibbsFreeEnergyRegion2(temperature, pressure) * SPECIFIC_GAS_CONSTANT * temperature;
}

double WaterSteamEquationOfState::CalculateRegion2SpecificEntropy(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double idealFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(temperature);
	const double residualFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(temperature, pressure);
	const double dimensionlessGibbsEnergy = CalculateDimensionlessGibbsFreeEnergyRegion2(temperature, pressure);
	const double dimensionlessEntropy = region2ReciprocalReducedTemperature * (idealFirstDerivativeGibbsTemperature + residualFirstDerivativeGibbsTemperature)
		- dimensionlessGibbsEnergy;
	const double specificEntropy = SPECIFIC_GAS_CONSTANT * dimensionlessEntropy;
	return specificEntropy;
}

double WaterSteamEquationOfState::CalculateRegion2SpecificEnthalpy(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double idealFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(temperature);
	const double residualFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(temperature, pressure);
	const double dimensionlessEnthalpy = region2ReciprocalReducedTemperature * (idealFirstDerivativeGibbsTemperature + residualFirstDerivativeGibbsTemperature);
	const double specificEnthalpy = SPECIFIC_GAS_CONSTANT * temperature * dimensionlessEnthalpy;
	return specificEnthalpy;
}

double WaterSteamEquationOfState::CalculateRegion2SpecificIsochoricHeatCapacity(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	const double specificIsobaricHeatCapacity = CalculateRegion2SpecificIsobaricHeatCapacity(temperature, pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(temperature, pressure);
	const double residualSecondDerivativeGibbsPressure = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(temperature, pressure);
	const double residualSecondMixedDerivativeGibbsPressureTemperature = CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion2(temperature, pressure);
	const double numerator = 1 + region2ReducedPressure * residualFirstDerivativeGibbsPressure - region2ReciprocalReducedTemperature * region2ReducedPressure * residualSecondMixedDerivativeGibbsPressureTemperature;
	const double denominator = 1 - pow(region2ReducedPressure, 2.0) * residualSecondDerivativeGibbsPressure;
	const double specificIsochoricHeatCapacity = specificIsobaricHeatCapacity - SPECIFIC_GAS_CONSTANT * (pow(numerator, 2.0) / denominator);
	return specificIsochoricHeatCapacity;
}

double WaterSteamEquationOfState::CalculateRegion2SpecificIsobaricHeatCapacity(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double idealSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(temperature);
	const double residualSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(temperature, pressure);
	const double specificIsobaricHeatCapacity = -1.0
		* region2ReciprocalReducedTemperature
		* region2ReciprocalReducedTemperature
		* SPECIFIC_GAS_CONSTANT
		* (idealSecondDerivativeGibbsTemperature + residualSecondDerivativeGibbsTemperature);
	return specificIsobaricHeatCapacity;
}

double WaterSteamEquationOfState::CalculateRegion2SpeedOfSound(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(temperature, pressure);
	const double residualSecondDerivativeGibbsPressure = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(temperature, pressure);
	const double idealSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(temperature);
	const double residualSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(temperature, pressure);
	const double residualSecondMixedDerivativeGibbsPressureTemperature = CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion2(temperature, pressure);
	const double pressureTermSquared = pow(region2ReducedPressure, 2.0);
	const double residualFirstTermSquared = pow(residualFirstDerivativeGibbsPressure, 2.0);
	const double reciprocalReducedTemperatureSquared = pow(region2ReciprocalReducedTemperature, 2.0);

	// Numerator calculation: (1 + 2 * pi * gr_p + pi^2 * (gr_p)^2)
	const double numerator = 1.0 + 2.0 * region2ReducedPressure * residualFirstDerivativeGibbsPressure + pressureTermSquared * residualFirstTermSquared;

	// Denominator term 1: (1 - pi^2 * gr_pp)
	const double denominatorTerm1 = 1.0 - pressureTermSquared * residualSecondDerivativeGibbsPressure;

	// Denominator term 2: (1 + pi * gr_p - tau * pi * gr_pt)^2 / (tau^2 * (g0_tt + gr_tt))
	const double denominatorTerm2Numerator = 1.0
		+ region2ReducedPressure * residualFirstDerivativeGibbsPressure
		- region2ReciprocalReducedTemperature * region2ReducedPressure * residualSecondMixedDerivativeGibbsPressureTemperature;
	const double denominatorTerm2Denominator = reciprocalReducedTemperatureSquared * (idealSecondDerivativeGibbsTemperature + residualSecondDerivativeGibbsTemperature);
	if (fabs(denominatorTerm2Denominator) < 1e-12) {
		throw std::runtime_error("Denominator term too small in speed of sound calculation.");
	}
	const double denominatorTerm2 = pow(denominatorTerm2Numerator, 2.0) / denominatorTerm2Denominator;
	const double denominator = denominatorTerm1 + denominatorTerm2;

	if (denominator <= 0) {
		throw std::runtime_error("Non-positive denominator in speed of sound calculation.");
	}
	const double speedOfSound = sqrt(SPECIFIC_GAS_CONSTANT * temperature * (numerator / denominator) * 1000.0);
	return speedOfSound;
}

double WaterSteamEquationOfState::CalculateRegion2ReducedPressure(const double pressure)
{
	return pressure / REGION2_REDUCING_PRESSURE;
}

double WaterSteamEquationOfState::CalculateRegion2ReciprocalReducedTemperature(const double temperature)
{
	return REGION2_REDUCING_TEMPERATURE / temperature;
}

double WaterSteamEquationOfState::CalculateRegion2ReciprocalReducedPressure(const double pressure)
{
	return REGION2_REDUCING_PRESSURE / pressure;
}

double WaterSteamEquationOfState::CalculateRegion2ReducedTemperature(const double temperature)
{
	return temperature / REGION2_REDUCING_TEMPERATURE;
}

double WaterSteamEquationOfState::CalculateDimensionlessGibbsFreeEnergyRegion2(const double temperature, const double pressure) const
{
	const double idealGibbsFreeEnergy = CalculateDimensionlessIdealGibbsFreeEnergyRegion2(temperature, pressure);
	const double residualGibbsFreeEnergy = CalculateDimensionlessResidualGibbsFreeEnergyRegion2(temperature, pressure);
	return idealGibbsFreeEnergy + residualGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateDimensionlessIdealGibbsFreeEnergyRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	const double logReducedPressure = log(region2ReducedPressure);
	double idealGibbsEnergy = logReducedPressure;
	for (const auto& coefficient : Region2IdealCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region2ReciprocalReducedTemperature, coefficient.Ji * 1.0);
		idealGibbsEnergy += term;
	}
	return idealGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion2(const double pressure)
{
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	return 1.0 / region2ReducedPressure;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion2(const double pressure)
{
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	return -1.0 * (1.0 / pow(region2ReducedPressure, 2.0));
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(const double temperature) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	double idealGibbsEnergy = 0.0;
	for (const auto& coefficient : Region2IdealCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* coefficient.Ji * 1.0
			* pow(region2ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 1.0);
		idealGibbsEnergy += term;
	}
	return idealGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion2(const double temperature) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	double idealGibbsEnergy = 0.0;
	for (const auto& coefficient : Region2IdealCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* coefficient.Ji * 1.0
			* (coefficient.Ji * 1.0 - 1.0)
			* pow(region2ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 2.0);
		idealGibbsEnergy += term;
	}
	return idealGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateDimensionlessResidualGibbsFreeEnergyRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	const double temperatureTerm = region2ReciprocalReducedTemperature - 0.5;
	for (const auto& coefficient : Region2ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region2ReducedPressure, coefficient.Ii * 1.0)
			* pow(temperatureTerm, coefficient.Ji * 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	const double temperatureTerm = region2ReciprocalReducedTemperature - 0.5;
	for (const auto& coefficient : Region2ResidualCoefficients) {
		const double pressureTerm = pow(region2ReducedPressure, coefficient.Ii - 1.0);
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent) * coefficient.Ii * pressureTerm * pow(temperatureTerm, coefficient.Ji);

		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	const double temperatureTerm = region2ReciprocalReducedTemperature - 0.5;
	for (const auto& coefficient : Region2ResidualCoefficients) {
		const double pressureTerm = pow(region2ReducedPressure, coefficient.Ii - 2.0);
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent) * coefficient.Ii * (coefficient.Ii - 1.0) * pressureTerm * pow(temperatureTerm, coefficient.Ji);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	const double temperatureTerm = region2ReciprocalReducedTemperature - 0.5;
	for (const auto& coefficient : Region2ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region2ReducedPressure, coefficient.Ii * 1.0)
			* (coefficient.Ji * 1.0)
			* pow(temperatureTerm, coefficient.Ji * 1.0 - 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	const double temperatureTerm = region2ReciprocalReducedTemperature - 0.5;
	for (const auto& coefficient : Region2ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region2ReducedPressure, coefficient.Ii * 1.0)
			* (coefficient.Ji * 1.0)
			* (coefficient.Ji * 1.0 - 1.0)
			* pow(temperatureTerm, coefficient.Ji * 1.0 - 2.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion2(const double temperature, const double pressure) const
{
	const double region2ReciprocalReducedTemperature = CalculateRegion2ReciprocalReducedTemperature(temperature);
	const double region2ReducedPressure = CalculateRegion2ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	const double temperatureTerm = region2ReciprocalReducedTemperature - 0.5;
	for (const auto& coefficient : Region2ResidualCoefficients) {
		const double pressureTerm = pow(region2ReducedPressure, coefficient.Ii - 1.0);
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* coefficient.Ii * pressureTerm
			* (coefficient.Ji * 1.0)
			* pow(temperatureTerm, coefficient.Ji * 1.0 - 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

#pragma endregion

#pragma region Region 23 Auxiliary Boundary Equations

double WaterSteamEquationOfState::CalculateRegion23ReducedPressure(const double pressure)
{
	return pressure / REGION23_BOUNDARY_REDUCING_PRESSURE;
}

double WaterSteamEquationOfState::CalculateRegion23ReducedTemperature(const double temperature)
{
	return temperature / REGION23_BOUNDARY_REDUCING_TEMPERATURE;
}

double WaterSteamEquationOfState::CalculateRegion23ReciprocalReducedTemperature(const double temperature)
{
	return REGION23_BOUNDARY_REDUCING_TEMPERATURE / temperature;
}

double WaterSteamEquationOfState::CalculateRegion23ReciprocalReducedPressure(const double pressure)
{
	return REGION23_BOUNDARY_REDUCING_PRESSURE / pressure;
}

double WaterSteamEquationOfState::CalculateRegion23BoundaryTemperature(const double pressure) const
{
	const double reducedPressure = CalculateRegion23ReducedPressure(pressure);
	const double n3 = Region23Coefficients[2].NiBase * pow(10, Region23Coefficients[2].NiExponent * 1.0);
	const double n4 = Region23Coefficients[3].NiBase * pow(10, Region23Coefficients[3].NiExponent * 1.0);
	const double n5 = Region23Coefficients[4].NiBase * pow(10, Region23Coefficients[4].NiExponent * 1.0);

	const double squaredTerm = (reducedPressure - n5) / n3;
	const double reducedTemperature = n4 + sqrt(squaredTerm);

	return reducedTemperature * REGION23_BOUNDARY_REDUCING_TEMPERATURE;
}

double WaterSteamEquationOfState::CalculateRegion23BoundaryPressure(const double temperature) const
{
	const double reducedTemperature = CalculateRegion23ReducedTemperature(temperature);
	const double n1 = Region23Coefficients[0].NiBase * pow(10, Region23Coefficients[0].NiExponent * 1.0);
	const double n2 = Region23Coefficients[1].NiBase * pow(10, Region23Coefficients[1].NiExponent * 1.0);
	const double n3 = Region23Coefficients[2].NiBase * pow(10, Region23Coefficients[2].NiExponent * 1.0);

	const double reducedPressure = n1 + n2 * reducedTemperature + n3 * pow(reducedTemperature, 2.0);

	return reducedPressure * REGION23_BOUNDARY_REDUCING_PRESSURE;
}

#pragma endregion

#pragma region Region 3 Equations

double WaterSteamEquationOfState::CalculateRegion3ReducedTemperature(const double temperature)
{
	return temperature / CRITICAL_TEMPERATURE;
}

double WaterSteamEquationOfState::CalculateRegion3ReciprocalReducedTemperature(const double temperature)
{
	return CRITICAL_TEMPERATURE / temperature;
}

double WaterSteamEquationOfState::CalculateRegion3ReducedDensity(const double density)
{
	return density / CRITICAL_DENSITY;
}

double WaterSteamEquationOfState::CalculateRegion3SpecificHelmholtzFreeEnergy(const double temperature, const double density) const
{
	const double dimensionlessHelmholtzEnergy = CalculateRegion3DimensionlessHelmholtzEnergy(temperature, density);
	return dimensionlessHelmholtzEnergy * temperature * SPECIFIC_GAS_CONSTANT;
}

double WaterSteamEquationOfState::CalculateRegion3DimensionlessHelmholtzEnergy(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessHelmholtzEnergy = Region3Coefficients[0].NiBase * pow(10.0, Region3Coefficients[0].NiExponent * 1.0) * log(reducedDensity);
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* pow(reducedDensity, Region3Coefficients[i].Ii * 1.0)
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji * 1.0);
		dimensionlessHelmholtzEnergy += term;
	}
	return dimensionlessHelmholtzEnergy;
}

double WaterSteamEquationOfState::CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyTemperature(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessHelmholtzEnergyDerivative = 0.0;
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* pow(reducedDensity, Region3Coefficients[i].Ii * 1.0)
			* Region3Coefficients[i].Ji * 1.0
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji * 1.0 - 1.0);
		dimensionlessHelmholtzEnergyDerivative += term;
	}
	return dimensionlessHelmholtzEnergyDerivative;

}

double WaterSteamEquationOfState::CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyTemperature(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessHelmholtzEnergyDerivative = 0.0;
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* pow(reducedDensity, Region3Coefficients[i].Ii * 1.0)
			* Region3Coefficients[i].Ji * 1.0
			* (Region3Coefficients[i].Ji * 1.0 - 1.0)
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji * 1.0 - 2.0);
		dimensionlessHelmholtzEnergyDerivative += term;
	}
	return dimensionlessHelmholtzEnergyDerivative;

}

double WaterSteamEquationOfState::CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyDensity(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessHelmholtzEnergyDerivative = Region3Coefficients[0].NiBase * pow(10.0, Region3Coefficients[0].NiExponent * 1.0) / reducedDensity;
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* Region3Coefficients[i].Ii * 1.0
			* pow(reducedDensity, Region3Coefficients[i].Ii * 1.0 - 1.0)
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji * 1.0);
		dimensionlessHelmholtzEnergyDerivative += term;
	}
	return dimensionlessHelmholtzEnergyDerivative;
}

double WaterSteamEquationOfState::CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyDensity(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessHelmholtzEnergyDerivative = -1.0 * Region3Coefficients[0].NiBase * pow(10.0, Region3Coefficients[0].NiExponent * 1.0) / pow(reducedDensity, 2.0);
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* Region3Coefficients[i].Ii * 1.0
			* (Region3Coefficients[i].Ii * 1.0 - 1.0)
			* pow(reducedDensity, Region3Coefficients[i].Ii * 1.0 - 2.0)
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji * 1.0);
		dimensionlessHelmholtzEnergyDerivative += term;
	}
	return dimensionlessHelmholtzEnergyDerivative;
}

double WaterSteamEquationOfState::CalculateRegion3SecondMixedDerivativeDimensionlessHelmholtzEnergyDensityTemperature(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessHelmholtzEnergyDerivative = 0.0;
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* Region3Coefficients[i].Ii * 1.0
			* pow(reducedDensity, Region3Coefficients[i].Ii * 1.0 - 1.0)
			* Region3Coefficients[i].Ji * 1.0
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji * 1.0 - 1.0);
		dimensionlessHelmholtzEnergyDerivative += term;
	}
	return dimensionlessHelmholtzEnergyDerivative;
}

double WaterSteamEquationOfState::CalculateRegion3Pressure(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	double dimensionlessPressure = Region3Coefficients[0].NiBase * pow(10.0, Region3Coefficients[0].NiExponent * 1.0) / reducedDensity;
	for (int i = 1; i < 40; i++)
	{
		const double term = Region3Coefficients[i].NiBase * pow(10.0, Region3Coefficients[i].NiExponent * 1.0)
			* Region3Coefficients[i].Ii * pow(reducedDensity, Region3Coefficients[i].Ii - 1.0)
			* pow(reciprocalReducedTemperature, Region3Coefficients[i].Ji);
		dimensionlessPressure += term;
	}
	return reducedDensity * dimensionlessPressure * SPECIFIC_GAS_CONSTANT * temperature * density * 10e-4;
}

double WaterSteamEquationOfState::CalculateRegion3SpecificInternalEnergy(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double helmholtzFirstDerivativeTemperature = CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyTemperature(temperature, density);
	return helmholtzFirstDerivativeTemperature * temperature * SPECIFIC_GAS_CONSTANT * reciprocalReducedTemperature;
}

double WaterSteamEquationOfState::CalculateRegion3SpecificEntropy(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double dimensionlessHelmholtzEnergy = CalculateRegion3DimensionlessHelmholtzEnergy(temperature, density);
	const double helmholtzFirstDerivativeTemperature = CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyTemperature(temperature, density);
	return (reciprocalReducedTemperature * helmholtzFirstDerivativeTemperature - dimensionlessHelmholtzEnergy) * SPECIFIC_GAS_CONSTANT;
}

double WaterSteamEquationOfState::CalculateRegion3SpecificEnthalpy(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	const double firstDerivativeHelmholtzTemperature = CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyTemperature(temperature, density);
	const double firstDerivativeHelmholtzDensity = CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyDensity(temperature, density);
	const double term1 = reciprocalReducedTemperature * firstDerivativeHelmholtzTemperature;
	const double term2 = reducedDensity * firstDerivativeHelmholtzDensity;
	return SPECIFIC_GAS_CONSTANT * temperature * (term1 + term2);
}

double WaterSteamEquationOfState::CalculateRegion3SpecificIsochoricHeatCapacity(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double secondDerivativeHelmholtzTemperature = CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyTemperature(temperature, density);
	return -1.0 * SPECIFIC_GAS_CONSTANT * pow(reciprocalReducedTemperature, 2.0) * secondDerivativeHelmholtzTemperature;
}

double WaterSteamEquationOfState::CalculateRegion3SpecificIsobaricHeatCapacity(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);
	const double firstDerivativeHelmholtzDensity = CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyDensity(temperature, density);
	const double secondDerivativeHelmholtzDensity = CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyDensity(temperature, density);
	const double mixedDerivativeHelmholtzDensityTemperature = CalculateRegion3SecondMixedDerivativeDimensionlessHelmholtzEnergyDensityTemperature(temperature, density);

	const double numerator = reducedDensity * firstDerivativeHelmholtzDensity - reducedDensity * reciprocalReducedTemperature * mixedDerivativeHelmholtzDensityTemperature;
	const double numeratorSquared = pow(numerator, 2.0);
	const double denominator = 2.0 * reducedDensity * firstDerivativeHelmholtzDensity + pow(reducedDensity, 2.0) * secondDerivativeHelmholtzDensity;
	const double specificIsochoricHeatCapacity = CalculateRegion3SpecificIsochoricHeatCapacity(temperature, density);
	return specificIsochoricHeatCapacity + SPECIFIC_GAS_CONSTANT * (numeratorSquared / denominator);
}

double WaterSteamEquationOfState::CalculateRegion3SpeedOfSound(const double temperature, const double density) const
{
	const double reciprocalReducedTemperature = CalculateRegion3ReciprocalReducedTemperature(temperature);
	const double reducedDensity = CalculateRegion3ReducedDensity(density);

	const double firstDerivativeHelmholtzDensity = CalculateRegion3FirstDerivativeDimensionlessHelmholtzEnergyDensity(temperature, density);
	const double secondDerivativeHelmholtzDensity = CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyDensity(temperature, density);
	const double mixedDerivativeHelmholtzDensityTemperature = CalculateRegion3SecondMixedDerivativeDimensionlessHelmholtzEnergyDensityTemperature(temperature, density);
	const double secondDerivativeHelmholtzTemperature = CalculateRegion3SecondDerivativeDimensionlessHelmholtzEnergyTemperature(temperature, density);

	const double numerator = reducedDensity * firstDerivativeHelmholtzDensity - reducedDensity * reciprocalReducedTemperature * mixedDerivativeHelmholtzDensityTemperature;
	const double numeratorSquared = pow(numerator, 2.0);
	const double denominator = pow(reciprocalReducedTemperature, 2.0) * secondDerivativeHelmholtzTemperature;
	const double rootedDimensionlessSpeed = 2.0 * reducedDensity * firstDerivativeHelmholtzDensity
		+ pow(reducedDensity, 2.0) * secondDerivativeHelmholtzDensity
		- numeratorSquared / denominator;
	const double speedOfSoundSquared = SPECIFIC_GAS_CONSTANT * temperature * rootedDimensionlessSpeed * 1000;
	return sqrt(speedOfSoundSquared);

}

#pragma endregion

#pragma region Region 4 Equations

double WaterSteamEquationOfState::CalculateRegion4SaturationPressure(const double temperature) const
{
	const double temperatureRatio = temperature / 1.0;
	const double theta = temperatureRatio + Region4Coefficients[8].NiBase * pow(10.0, Region4Coefficients[8].NiExponent * 1.0)
		/ (temperatureRatio - Region4Coefficients[9].NiBase * pow(10.0, Region4Coefficients[9].NiExponent * 1.0));

	const double rootCoefficientA = pow(theta, 2.0)
		+ theta * (Region4Coefficients[0].NiBase * pow(10.0, Region4Coefficients[0].NiExponent * 1.0))
		+ Region4Coefficients[1].NiBase * pow(10.0, Region4Coefficients[1].NiExponent * 1.0);

	const double rootCoefficientB = Region4Coefficients[2].NiBase * pow(10.0, Region4Coefficients[2].NiExponent * 1.0) * pow(theta, 2.0)
		+ theta * (Region4Coefficients[3].NiBase * pow(10.0, Region4Coefficients[3].NiExponent * 1.0))
		+ Region4Coefficients[4].NiBase * pow(10.0, Region4Coefficients[4].NiExponent * 1.0);

	const double rootCoefficientC = Region4Coefficients[5].NiBase * pow(10.0, Region4Coefficients[5].NiExponent * 1.0) * pow(theta, 2.0)
		+ theta * (Region4Coefficients[6].NiBase * pow(10.0, Region4Coefficients[6].NiExponent * 1.0))
		+ Region4Coefficients[7].NiBase * pow(10.0, Region4Coefficients[7].NiExponent * 1.0);

	const double discriminant = pow(rootCoefficientB, 2.0) - 4.0 * rootCoefficientA * rootCoefficientC;
	const double root = 2.0 * rootCoefficientC / (-rootCoefficientB + sqrt(discriminant));
	const double saturationPressure = 1.0 * pow(root, 4.0);
	return saturationPressure;
}

double WaterSteamEquationOfState::CalculateRegion4SaturationTemperature(const double pressure) const
{
	const double pressureRatio = pressure / 1.0;
	const double beta = pow(pressureRatio, 0.25);

	const double equationE = pow(beta, 2.0)
		+ Region4Coefficients[2].NiBase * pow(10.0, Region4Coefficients[2].NiExponent * 1.0) * beta
		+ Region4Coefficients[5].NiBase * pow(10.0, Region4Coefficients[5].NiExponent * 1.0);
	const double equationF = Region4Coefficients[0].NiBase * pow(10.0, Region4Coefficients[0].NiExponent * 1.0) * pow(beta, 2.0)
		+ Region4Coefficients[3].NiBase * pow(10.0, Region4Coefficients[3].NiExponent * 1.0) * beta
		+ Region4Coefficients[6].NiBase * pow(10.0, Region4Coefficients[6].NiExponent * 1.0);
	const double equationG = Region4Coefficients[1].NiBase * pow(10.0, Region4Coefficients[1].NiExponent * 1.0) * pow(beta, 2.0)
		+ Region4Coefficients[4].NiBase * pow(10.0, Region4Coefficients[4].NiExponent * 1.0) * beta
		+ Region4Coefficients[7].NiBase * pow(10.0, Region4Coefficients[7].NiExponent * 1.0);

	const double equationD = 2.0 * equationG / (-equationF - sqrt(pow(equationF, 2.0) - 4.0 * equationE * equationG));
	const double term = Region4Coefficients[9].NiExponent * pow(10.0, Region4Coefficients[9].NiExponent) + equationD;
	const double saturationTemperature = 0.5 * (term - sqrt(pow(term, 2.0) - 4 * (Region4Coefficients[8].NiExponent * pow(10.0, Region4Coefficients[8].NiExponent) + Region4Coefficients[9].NiExponent * pow(10.0, Region4Coefficients[9].NiExponent) * equationD)));
	return saturationTemperature;
}

#pragma endregion

#pragma region Region 5 Equations

double WaterSteamEquationOfState::CalculateRegion5SpecificVolume(const double temperature, const double pressure) const
{
	const double pressureKiloPascals = pressure * 1000;
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	const double idealFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion5(pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(temperature, pressure);
	const double dimensionlessPressure = region5ReducedPressure * (idealFirstDerivativeGibbsPressure + residualFirstDerivativeGibbsPressure);
	const double specificVolume = SPECIFIC_GAS_CONSTANT * temperature * dimensionlessPressure / pressureKiloPascals;
	return specificVolume;
}

double WaterSteamEquationOfState::CalculateRegion5SpecificInternalEnergy(const double temperature, const double pressure) const
{
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);

	const double idealFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion5(pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(temperature, pressure);

	const double idealFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(temperature);
	const double residualFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(temperature, pressure);

	const double dimensionlessEnergy = region5ReciprocalReducedTemperature * (idealFirstDerivativeGibbsTemperature + residualFirstDerivativeGibbsTemperature)
		- region5ReducedPressure * (idealFirstDerivativeGibbsPressure + residualFirstDerivativeGibbsPressure);
	const double specificInternalEnergy = SPECIFIC_GAS_CONSTANT * temperature * dimensionlessEnergy;
	return specificInternalEnergy;
}

double WaterSteamEquationOfState::CalculateRegion5SpecificEntropy(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);

	const double idealFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(temperature);
	const double residualFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(temperature, pressure);
	const double dimensionlessGibbsEnergy = CalculateDimensionlessGibbsFreeEnergyRegion5(temperature, pressure);
	const double dimensionlessEntropy = region5ReciprocalReducedTemperature * (idealFirstDerivativeGibbsTemperature + residualFirstDerivativeGibbsTemperature)
		- dimensionlessGibbsEnergy;
	const double specificEntropy = SPECIFIC_GAS_CONSTANT * dimensionlessEntropy;
	return specificEntropy;
}

double WaterSteamEquationOfState::CalculateRegion5SpecificEnthalpy(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double idealFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(temperature);
	const double residualFirstDerivativeGibbsTemperature = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(temperature, pressure);
	const double dimensionlessEnthalpy = region5ReciprocalReducedTemperature * (idealFirstDerivativeGibbsTemperature + residualFirstDerivativeGibbsTemperature);
	const double specificEnthalpy = SPECIFIC_GAS_CONSTANT * temperature * dimensionlessEnthalpy;
	return specificEnthalpy;
}

double WaterSteamEquationOfState::CalculateRegion5SpecificIsochoricHeatCapacity(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	const double specificIsobaricHeatCapacity = CalculateRegion5SpecificIsobaricHeatCapacity(temperature, pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(temperature, pressure);
	const double residualSecondDerivativeGibbsPressure = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(temperature, pressure);
	const double residualSecondMixedDerivativeGibbsPressureTemperature = CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion5(temperature, pressure);
	const double numerator = 1 + region5ReducedPressure * residualFirstDerivativeGibbsPressure - region5ReciprocalReducedTemperature * region5ReducedPressure * residualSecondMixedDerivativeGibbsPressureTemperature;
	const double denominator = 1 - pow(region5ReducedPressure, 2.0) * residualSecondDerivativeGibbsPressure;
	const double specificIsochoricHeatCapacity = specificIsobaricHeatCapacity - SPECIFIC_GAS_CONSTANT * (pow(numerator, 2.0) / denominator);
	return specificIsochoricHeatCapacity;
}

double WaterSteamEquationOfState::CalculateRegion5SpecificIsobaricHeatCapacity(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double idealSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(temperature);
	const double residualSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(temperature, pressure);
	const double specificIsobaricHeatCapacity = -1.0
		* region5ReciprocalReducedTemperature
		* region5ReciprocalReducedTemperature
		* SPECIFIC_GAS_CONSTANT
		* (idealSecondDerivativeGibbsTemperature + residualSecondDerivativeGibbsTemperature);
	return specificIsobaricHeatCapacity;
}

double WaterSteamEquationOfState::CalculateRegion5SpeedOfSound(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	const double residualFirstDerivativeGibbsPressure = CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(temperature, pressure);
	const double residualSecondDerivativeGibbsPressure = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(temperature, pressure);
	const double idealSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(temperature);
	const double residualSecondDerivativeGibbsTemperature = CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(temperature, pressure);
	const double residualSecondMixedDerivativeGibbsPressureTemperature = CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion5(temperature, pressure);
	const double pressureTermSquared = pow(region5ReducedPressure, 2.0);
	const double residualFirstTermSquared = pow(residualFirstDerivativeGibbsPressure, 2.0);
	const double reciprocalReducedTemperatureSquared = pow(region5ReciprocalReducedTemperature, 2.0);

	// Numerator calculation: (1 + 2 * pi * gr_p + pi^2 * (gr_p)^2)
	const double numerator = 1.0 + 2.0 * region5ReducedPressure * residualFirstDerivativeGibbsPressure + pressureTermSquared * residualFirstTermSquared;

	// Denominator term 1: (1 - pi^2 * gr_pp)
	const double denominatorTerm1 = 1.0 - pressureTermSquared * residualSecondDerivativeGibbsPressure;

	// Denominator term 2: (1 + pi * gr_p - tau * pi * gr_pt)^2 / (tau^2 * (g0_tt + gr_tt))
	const double denominatorTerm2Numerator = 1.0
		+ region5ReducedPressure * residualFirstDerivativeGibbsPressure
		- region5ReciprocalReducedTemperature * region5ReducedPressure * residualSecondMixedDerivativeGibbsPressureTemperature;
	const double denominatorTerm2Denominator = reciprocalReducedTemperatureSquared * (idealSecondDerivativeGibbsTemperature + residualSecondDerivativeGibbsTemperature);
	if (fabs(denominatorTerm2Denominator) < 1e-12) {
		throw std::runtime_error("Denominator term too small in speed of sound calculation.");
	}
	const double denominatorTerm2 = pow(denominatorTerm2Numerator, 2.0) / denominatorTerm2Denominator;
	const double denominator = denominatorTerm1 + denominatorTerm2;

	if (denominator <= 0) {
		throw std::runtime_error("Non-positive denominator in speed of sound calculation.");
	}
	const double speedOfSound = sqrt(SPECIFIC_GAS_CONSTANT * temperature * (numerator / denominator) * 1000.0);
	return speedOfSound;
}

double WaterSteamEquationOfState::CalculateRegion5SpecificGibbsFreeEnergy(const double temperature, const double pressure) const
{
	return CalculateDimensionlessGibbsFreeEnergyRegion5(temperature, pressure) * SPECIFIC_GAS_CONSTANT * temperature;
}

double WaterSteamEquationOfState::CalculateRegion5ReciprocalReducedTemperature(const double temperature)
{
	return REGION5_REDUCING_TEMPERATURE / temperature;
}

double WaterSteamEquationOfState::CalculateRegion5ReducedPressure(const double pressure)
{
	return pressure / REGION5_REDUCING_PRESSURE;
}

double WaterSteamEquationOfState::CalculateDimensionlessGibbsFreeEnergyRegion5(const double temperature, const double pressure) const
{
	const double idealGibbsFreeEnergy = CalculateDimensionlessIdealGibbsFreeEnergyRegion5(temperature, pressure);
	const double residualGibbsFreeEnergy = CalculateDimensionlessResidualGibbsFreeEnergyRegion5(temperature, pressure);
	return idealGibbsFreeEnergy + residualGibbsFreeEnergy;
}

double WaterSteamEquationOfState::CalculateDimensionlessIdealGibbsFreeEnergyRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	const double logReducedPressure = log(region5ReducedPressure);
	double idealGibbsEnergy = logReducedPressure;
	for (const auto& coefficient : Region5IdealCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0);
		idealGibbsEnergy += term;
	}
	return idealGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion5(const double pressure)
{
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	return 1.0 / region5ReducedPressure;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyPressureRegion5(const double pressure)
{
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	return -1.0 * (1.0 / pow(region5ReducedPressure, 2.0));
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(const double temperature) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	double idealGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5IdealCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* coefficient.Ji * pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 1.0);
		idealGibbsEnergy += term;
	}
	return idealGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessIdealGibbsFreeEnergyTemperatureRegion5(const double temperature) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	double idealGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5IdealCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* coefficient.Ji * (coefficient.Ji - 1.0) * pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 2.0);
		idealGibbsEnergy += term;
	}
	return idealGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateDimensionlessResidualGibbsFreeEnergyRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region5ReducedPressure, coefficient.Ii * 1.0)
			* pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent)
			* coefficient.Ii * pow(region5ReducedPressure, coefficient.Ii - 1.0)
			* pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyPressureRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent)
			* coefficient.Ii * (coefficient.Ii - 1.0) * pow(region5ReducedPressure, coefficient.Ii - 2.0)
			* pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateFirstDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region5ReducedPressure, coefficient.Ii * 1.0)
			* coefficient.Ji * pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondDerivativeDimensionlessResidualGibbsFreeEnergyTemperatureRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* pow(region5ReducedPressure, coefficient.Ii * 1.0)
			* (coefficient.Ji * 1.0) * (coefficient.Ji * 1.0 - 1.0) * pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 2.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

double WaterSteamEquationOfState::CalculateSecondMixedDerivativeDimensionlessResidualGibbsFreeEnergyPressureTemperatureRegion5(const double temperature, const double pressure) const
{
	const double region5ReciprocalReducedTemperature = CalculateRegion5ReciprocalReducedTemperature(temperature);
	const double region5ReducedPressure = CalculateRegion5ReducedPressure(pressure);
	double residualGibbsEnergy = 0.0;
	for (const auto& coefficient : Region5ResidualCoefficients) {
		const double term = coefficient.NiBase * pow(10.0, coefficient.NiExponent * 1.0)
			* coefficient.Ii * pow(region5ReducedPressure, coefficient.Ii - 1.0)
			* coefficient.Ji * pow(region5ReciprocalReducedTemperature, coefficient.Ji * 1.0 - 1.0);
		residualGibbsEnergy += term;
	}
	return residualGibbsEnergy;
}

#pragma endregion