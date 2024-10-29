#include "pch.h"
#include "CppUnitTest.h"
#include "WaterSteamEquationOfState.h"

using namespace Microsoft::VisualStudio::CppUnitTestFramework;


namespace WaterSteamEquationOfStateTests
{
	struct Region125TestCase
	{
		double Temperature;
		double Pressure;
		double ExpectedSpecificVolume;
		double ExpectedSpecificInternalEnergy;
		double ExpectedSpecificEntropy;
		double ExpectedSpecificEnthalpy;
		double ExpectedSpecificIsobaricHeatCapacity;
		double ExpectedSpeedOfSound;
	};

	struct Region3TestCase
	{
		double Temperature;
		double Density;
		double ExpectedPressure;
		double ExpectedSpecificInternalEnergy;
		double ExpectedSpecificEntropy;
		double ExpectedSpecificEnthalpy;
		double ExpectedSpecificIsobaricHeatCapacity;
		double ExpectedSpeedOfSound;
	};

	struct Region4PressureTestCase
	{
		double Temperature;
		double ExpectedSaturationPressure;
	};

	struct Region4TemperatureTestCase
	{
		double Pressure;
		double ExpectedSaturationTemperature;
	};

	TEST_CLASS(WaterSteamEquationOfStateTestClass)
	{
	public:
		static WaterSteamEquationOfState* WaterEquationOfState;
		static constexpr double TOLERANCE = 1e-5;

		TEST_CLASS_INITIALIZE(ClassInitialize)
		{
			WaterEquationOfState = nullptr;
			const std::string databasePath = "C:\\Research\\Databases\\thermodynamic_properties.db";
			WaterEquationOfState = new WaterSteamEquationOfState(databasePath);
		}

		TEST_CLASS_CLEANUP(ClassCleanup)
		{
			delete WaterEquationOfState;
			WaterEquationOfState = nullptr;
		}

		TEST_METHOD(Test_Region1_VerificationCases)
		{
			Region125TestCase testCases[] = {
				// IF97 Table 5 Case 1
				{300.0, 3.0, 0.00100215168, 112.324818, 0.392294792, 115.331273 , 4.17301218, 1507.73921},
				// IF97 Table 5 Case 2
				{300.0, 80.0, 0.000971180894, 106.448356, 0.368563852, 184.142828, 4.01008987, 1634.69054},
				// IF97 Table 5 Case 3
				{500.0, 3.0, 0.00120241800, 971.934985, 2.58041912, 975.542239, 4.65580682, 1240.71337}
			};
			for (const auto& testCase : testCases)
			{
				VerifyRegion1Properties(testCase);
			}
		}

		TEST_METHOD(Test_Region2_VerificationCases)
		{
			Region125TestCase testCases[] = {
				// IF97 Table 15 Case 1
				{300.0, 0.0035, 39.4913866, 2411.69160, 8.52238967, 2549.91145, 1.91300162, 427.920172},
				// IF97 Table 15 Case 2
				{700.0, 0.0035, 92.3015898, 3012.62819, 10.1749996, 3335.68375, 2.08141274, 644.289068},
				// IF97 Table 15 Case 3
				{700.0, 30.0, 0.00542946619, 2468.61076, 5.17540298, 2631.49474, 10.3505092, 480.386523}
			};
			for (const auto& testCase : testCases)
			{
				VerifyRegion2Properties(testCase);
			}
		}

		TEST_METHOD(Test_Region2_3_BoundaryEquations)
		{
			constexpr double expectedTemperature = 623.15;
			constexpr double expectedPressure = 16.5292;

			const double actualTemperature = WaterEquationOfState->CalculateRegion23BoundaryTemperature(expectedPressure);
			const double actualPressure = WaterEquationOfState->CalculateRegion23BoundaryPressure(expectedTemperature);

			Assert::IsTrue(AreDoublesEqual(expectedTemperature, actualTemperature));
			Assert::IsTrue(AreDoublesEqual(expectedPressure, actualPressure));
		}

		TEST_METHOD(Test_Region3_VerificationCases)
		{
			Region3TestCase testCases[] = {
				// IF97 Table 33 Case 1
				{650.0, 500.0, 25.5837018, 1812.26279, 4.05427273, 1863.43019, 13.8935717, 502.005554},
				// IF97 Table 33 Case 2
				{650.0, 200.0, 22.2930643, 2263.65868, 4.85438792, 2375.12401, 44.6579342, 383.444594},
				// IF97 Table 33 Case 3
				{750.0, 500.0, 78.3095639, 2102.06932, 4.46971906, 2258.68845, 6.34165359, 760.696041}
			};
			for (const auto& testCase : testCases)
			{
				VerifyRegion3Properties(testCase);
			}
		}

		TEST_METHOD(Test_Region4_PressureVerificationCases)
		{
			Region4PressureTestCase testCases[] = {
				// IF97 Table 35 Case 1
				{300.0, 0.00353658941},
				// IF97 Table 35 Case 2
				{500.0, 2.63889776},
				// IF97 Table 35 Case 3
				{600.0, 12.3443146}
			};
			for (const auto& testCase : testCases)
			{
				VerifyRegion4PressureProperties(testCase);
			}
		}

		TEST_METHOD(Test_Region4_TemperatureVerificationCases)
		{
			Region4TemperatureTestCase testCases[] = {
				// IF97 Table 36 Case 1
				{0.1, 372.755919},
				// IF97 Table 36 Case 2
				{1.0, 453.035632},
				// IF97 Table 36 Case 3
				{10.0, 584.149488}
			};
			for (const auto& testCase : testCases)
			{
				VerifyRegion4TemperatureProperties(testCase);
			}
		}

		TEST_METHOD(Test_Region5_VerificationCases)
		{
			Region125TestCase testCases[] = {
				// IF97 Table 41 Case 1
				{1500.0, 0.5, 1.38455090, 4527.49310, 9.65408875, 5219.76855, 2.61609445, 917.068690},
				// IF97 Table 41 Case 2
				{1500.0, 30.0, 0.0230761299, 4474.95124, 7.72970133, 5167.23514, 2.72724317, 928.548002},
				// IF97 Table 41 Case 3
				{2000.0, 30.0, 0.0311385219, 5637.07038, 8.53640523, 6571.22604, 2.88569882, 1067.36948}
			};
			for (const auto& testCase : testCases)
			{
				VerifyRegion5Properties(testCase);
			}
		}
	private:
		static bool AreDoublesEqual(const double expected, const double actual)
		{
			const double difference = fabs(expected - actual);
			const double maxAbsValue = std::max(fabs(expected), fabs(actual));
			return (difference <= TOLERANCE * maxAbsValue);
		}

		static void VerifyRegion1Properties(const Region125TestCase& testCase)
		{
			const double specificVolume = WaterEquationOfState->CalculateRegion1SpecificVolume(testCase.Temperature, testCase.Pressure);
			const double specificInternalEnergy = WaterEquationOfState->CalculateRegion1SpecificInternalEnergy(testCase.Temperature, testCase.Pressure);
			const double specificEntropy = WaterEquationOfState->CalculateRegion1SpecificEntropy(testCase.Temperature, testCase.Pressure);
			const double specificEnthalpy = WaterEquationOfState->CalculateRegion1SpecificEnthalpy(testCase.Temperature, testCase.Pressure);
			const double specificIsobaricHeatCapacity = WaterEquationOfState->CalculateRegion1SpecificIsobaricHeatCapacity(testCase.Temperature, testCase.Pressure);
			const double speedOfSound = WaterEquationOfState->CalculateRegion1SpeedOfSound(testCase.Temperature, testCase.Pressure);

			Assert::AreEqual(testCase.ExpectedSpecificVolume, specificVolume, TOLERANCE, L"Specific volume mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificInternalEnergy, specificInternalEnergy, TOLERANCE, L"Specific internal energy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEntropy, specificEntropy, TOLERANCE, L"Specific entropy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEnthalpy, specificEnthalpy, TOLERANCE, L"Specific enthalpy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificIsobaricHeatCapacity, specificIsobaricHeatCapacity, TOLERANCE, L"Specific isobaric heat capacity mismatch.");
			Assert::AreEqual(testCase.ExpectedSpeedOfSound, speedOfSound, TOLERANCE, L"Speed of sound mismatch.");
		}

		static void VerifyRegion2Properties(const Region125TestCase& testCase)
		{
			const double specificVolume = WaterEquationOfState->CalculateRegion2SpecificVolume(testCase.Temperature, testCase.Pressure);
			const double specificInternalEnergy = WaterEquationOfState->CalculateRegion2SpecificInternalEnergy(testCase.Temperature, testCase.Pressure);
			const double specificEntropy = WaterEquationOfState->CalculateRegion2SpecificEntropy(testCase.Temperature, testCase.Pressure);
			const double specificEnthalpy = WaterEquationOfState->CalculateRegion2SpecificEnthalpy(testCase.Temperature, testCase.Pressure);
			const double specificIsobaricHeatCapacity = WaterEquationOfState->CalculateRegion2SpecificIsobaricHeatCapacity(testCase.Temperature, testCase.Pressure);
			const double speedOfSound = WaterEquationOfState->CalculateRegion2SpeedOfSound(testCase.Temperature, testCase.Pressure);

			Assert::AreEqual(testCase.ExpectedSpecificVolume, specificVolume, TOLERANCE, L"Specific volume mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificInternalEnergy, specificInternalEnergy, TOLERANCE, L"Specific internal energy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEntropy, specificEntropy, TOLERANCE, L"Specific entropy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEnthalpy, specificEnthalpy, TOLERANCE, L"Specific enthalpy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificIsobaricHeatCapacity, specificIsobaricHeatCapacity, TOLERANCE, L"Specific isobaric heat capacity mismatch.");
			Assert::AreEqual(testCase.ExpectedSpeedOfSound, speedOfSound, TOLERANCE, L"Speed of sound mismatch.");
		}

		static void VerifyRegion3Properties(const Region3TestCase& testCase)
		{
			const double pressure = WaterEquationOfState->CalculateRegion3Pressure(testCase.Temperature, testCase.Density);
			const double specificInternalEnergy = WaterEquationOfState->CalculateRegion3SpecificInternalEnergy(testCase.Temperature, testCase.Density);
			const double specificEntropy = WaterEquationOfState->CalculateRegion3SpecificEntropy(testCase.Temperature, testCase.Density);
			const double specificEnthalpy = WaterEquationOfState->CalculateRegion3SpecificEnthalpy(testCase.Temperature, testCase.Density);
			const double specificIsobaricHeatCapacity = WaterEquationOfState->CalculateRegion3SpecificIsobaricHeatCapacity(testCase.Temperature, testCase.Density);
			const double speedOfSound = WaterEquationOfState->CalculateRegion3SpeedOfSound(testCase.Temperature, testCase.Density);

			Assert::AreEqual(testCase.ExpectedPressure, pressure, TOLERANCE, L"Pressure mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificInternalEnergy, specificInternalEnergy, TOLERANCE, L"Specific internal energy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEntropy, specificEntropy, TOLERANCE, L"Specific entropy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEnthalpy, specificEnthalpy, TOLERANCE, L"Specific enthalpy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificIsobaricHeatCapacity, specificIsobaricHeatCapacity, TOLERANCE, L"Specific isobaric heat capacity mismatch.");
			Assert::AreEqual(testCase.ExpectedSpeedOfSound, speedOfSound, TOLERANCE, L"Speed of sound mismatch.");
		}

		static void VerifyRegion4PressureProperties(const Region4PressureTestCase& testCase)
		{
			const double saturationPressure = WaterEquationOfState->CalculateRegion4SaturationPressure(testCase.Temperature);
			Assert::AreEqual(testCase.ExpectedSaturationPressure, saturationPressure, TOLERANCE, L"Saturation pressure mismatch.");
		}

		static void VerifyRegion4TemperatureProperties(const Region4TemperatureTestCase& testCase)
		{
			const double saturationTemperature = WaterEquationOfState->CalculateRegion4SaturationTemperature(testCase.Pressure);
			Assert::AreEqual(testCase.ExpectedSaturationTemperature, saturationTemperature, 1e-2, L"Saturation pressure mismatch.");
		}

		static void VerifyRegion5Properties(const Region125TestCase& testCase)
		{
			const double specificVolume = WaterEquationOfState->CalculateRegion5SpecificVolume(testCase.Temperature, testCase.Pressure);
			const double specificInternalEnergy = WaterEquationOfState->CalculateRegion5SpecificInternalEnergy(testCase.Temperature, testCase.Pressure);
			const double specificEntropy = WaterEquationOfState->CalculateRegion5SpecificEntropy(testCase.Temperature, testCase.Pressure);
			const double specificEnthalpy = WaterEquationOfState->CalculateRegion5SpecificEnthalpy(testCase.Temperature, testCase.Pressure);
			const double specificIsobaricHeatCapacity = WaterEquationOfState->CalculateRegion5SpecificIsobaricHeatCapacity(testCase.Temperature, testCase.Pressure);
			const double speedOfSound = WaterEquationOfState->CalculateRegion5SpeedOfSound(testCase.Temperature, testCase.Pressure);

			Assert::AreEqual(testCase.ExpectedSpecificVolume, specificVolume, TOLERANCE, L"Specific volume mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificInternalEnergy, specificInternalEnergy, TOLERANCE, L"Specific internal energy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEntropy, specificEntropy, TOLERANCE, L"Specific entropy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificEnthalpy, specificEnthalpy, TOLERANCE, L"Specific enthalpy mismatch.");
			Assert::AreEqual(testCase.ExpectedSpecificIsobaricHeatCapacity, specificIsobaricHeatCapacity, TOLERANCE, L"Specific isobaric heat capacity mismatch.");
			Assert::AreEqual(testCase.ExpectedSpeedOfSound, speedOfSound, TOLERANCE, L"Speed of sound mismatch.");
		}
	};

	WaterSteamEquationOfState* WaterSteamEquationOfStateTestClass::WaterEquationOfState = nullptr;
}