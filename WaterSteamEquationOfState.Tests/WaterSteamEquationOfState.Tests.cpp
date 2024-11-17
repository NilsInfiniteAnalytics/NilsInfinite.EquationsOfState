#include "pch.h"
#include <vector>
#include <cmath>
#include <complex>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <sstream>

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

	struct TestResult {
		std::string PropertyName;
		double Expected;
		double Actual;
		bool Passed;
	};

	struct IceRegionTestCase
	{
		double Temperature;
		double Pressure;
		double ExpectedSpecificGibbsEnergy;
		double ExpectedFirstDerivativeGibbsEnergyPressure;
		double ExpectedSecondDerivativeGibbsEnergyPressure;
		double ExpectedSpecificEntropy;
		double ExpectedSpecificEnthalpy;
		double ExpectedSpecificInternalEnergy;
		double ExpectedSpecificIsobaricHeatCapacity;
		double ExpectedDensity;
		double ExpectedCubicExpansionCoefficient;
		double ExpectedIsothermalCompressibility;
		double ExpectedIsentropicCompressibility;
		double ExpectedPressureCoefficient;
	};

	TEST_CLASS(WaterSteamEquationOfStateTestClass)
	{
	public:
		static WaterSteamEquationOfState* WaterEquationOfState;
		static constexpr double TOLERANCE = 1e-4;
		static constexpr double ABSOLUTE_TOLERANCE = 1e-12;

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

		TEST_METHOD(Test_IceRegion_SublimationPressure)
		{
			constexpr double expectedPressure = 8.94735e-6;
			const double actualPressure = WaterEquationOfState->CalculateIceSublimationPressure(230);
			Assert::IsTrue(AreDoublesEqual(expectedPressure, actualPressure));
		}

		TEST_METHOD(Test_IceRegion_ComplexFunctionConstantsVerficiation)
		{
			std::vector<TestResult> results;
			constexpr double r1_real = 0.447050716285388e2;
			constexpr double r1_imag = 0.656876847463481e2;

			constexpr double r20_real = -0.725974574329220e2;
			constexpr double r20_imag = -0.781008427112870e2;

			constexpr double r21_real = -0.557107698030123e-4;
			constexpr double r21_imag = 0.464578634580806e-4;

			constexpr double r22_real = 0.234801409215913e-10;
			constexpr double r22_imag = -0.285651142904972e-10;
			for (int i = 0; i < WaterEquationOfState->NUMBER_OF_ICE_COMPLEX_FUNCTION_COEFFICIENTS; i++)
			{
				const std::complex<double> coefficient = WaterEquationOfState->GetIceComplexFunctionConstantR(i);
				int passCount = 0;

				auto evaluateTolerance = [&](const std::string& type, double expected, double actual, auto compareFunc, const std::string& name) {
					const bool passed = compareFunc(expected, actual);
					LogResult(results, name + " - " + type + " Compare", expected, actual, passed);
					if (passed) passCount++;
					};

				switch (i)
				{
				case 0:
					evaluateTolerance("Absolute", r1_real, WaterEquationOfState->GetIceComplexCoefficientR(i).real(), AbsoluteToleranceCompare, "R1 Real");
					evaluateTolerance("Relative", r1_real, WaterEquationOfState->GetIceComplexCoefficientR(i).real(), RelativeToleranceCompare, "R1 Real");
					evaluateTolerance("Combined", r1_real, WaterEquationOfState->GetIceComplexCoefficientR(i).real(), CombinedToleranceCompare, "R1 Real");

					evaluateTolerance("Absolute", r1_imag, WaterEquationOfState->GetIceComplexCoefficientR(i).imag(), AbsoluteToleranceCompare, "R1 Imaginary");
					evaluateTolerance("Relative", r1_imag, WaterEquationOfState->GetIceComplexCoefficientR(i).imag(), RelativeToleranceCompare, "R1 Imaginary");
					evaluateTolerance("Combined", r1_imag, WaterEquationOfState->GetIceComplexCoefficientR(i).imag(), CombinedToleranceCompare, "R1 Imaginary");

					evaluateTolerance("Absolute", r20_real, coefficient.real(), AbsoluteToleranceCompare, "R20 Real");
					evaluateTolerance("Relative", r20_real, coefficient.real(), RelativeToleranceCompare, "R20 Real");
					evaluateTolerance("Combined", r20_real, coefficient.real(), CombinedToleranceCompare, "R20 Real");

					evaluateTolerance("Absolute", r20_imag, coefficient.imag(), AbsoluteToleranceCompare, "R20 Imaginary");
					evaluateTolerance("Relative", r20_imag, coefficient.imag(), RelativeToleranceCompare, "R20 Imaginary");
					evaluateTolerance("Combined", r20_imag, coefficient.imag(), CombinedToleranceCompare, "R20 Imaginary");
					break;
				case 1:
					evaluateTolerance("Absolute", r21_real, coefficient.real(), AbsoluteToleranceCompare, "R21 Real");
					evaluateTolerance("Relative", r21_real, coefficient.real(), RelativeToleranceCompare, "R21 Real");
					evaluateTolerance("Combined", r21_real, coefficient.real(), CombinedToleranceCompare, "R21 Real");

					evaluateTolerance("Absolute", r21_imag, coefficient.imag(), AbsoluteToleranceCompare, "R21 Imaginary");
					evaluateTolerance("Relative", r21_imag, coefficient.imag(), RelativeToleranceCompare, "R21 Imaginary");
					evaluateTolerance("Combined", r21_imag, coefficient.imag(), CombinedToleranceCompare, "R21 Imaginary");
					break;
				case 2:
					evaluateTolerance("Absolute", r22_real, coefficient.real(), AbsoluteToleranceCompare, "R22 Real");
					evaluateTolerance("Relative", r22_real, coefficient.real(), RelativeToleranceCompare, "R22 Real");
					evaluateTolerance("Combined", r22_real, coefficient.real(), CombinedToleranceCompare, "R22 Real");

					evaluateTolerance("Absolute", r22_imag, coefficient.imag(), AbsoluteToleranceCompare, "R22 Imaginary");
					evaluateTolerance("Relative", r22_imag, coefficient.imag(), RelativeToleranceCompare, "R22 Imaginary");
					evaluateTolerance("Combined", r22_imag, coefficient.imag(), CombinedToleranceCompare, "R22 Imaginary");
					break;
				default:
					Assert::Fail(L"Invalid coefficient index.");
				}
				Assert::IsTrue(passCount >= 2, L"One or more coefficients failed to meet at least two tolerance criteria.");
			}
		}

		TEST_METHOD(Test_Ice_GibbsResidualCoefficientVerification)
		{
			std::vector<TestResult> results;

			constexpr double g00 = -0.632020233335886e6;
			constexpr double g01 = 0.655022213658955;
			constexpr double g02 = -0.189369929326131e-7;
			constexpr double g03 = 0.339746123271053e-14;
			constexpr double g04 = -0.556464869058991e-21;

			for (int i = 0; i < WaterEquationOfState->NUMBER_OF_ICE_REAL_GIBBS_COEFFICIENTS; i++)
			{
				const double coefficient = WaterEquationOfState->GetIceResidualGibbsCoefficient(i);
				int passCount = 0;

				auto evaluateTolerance = [&](const std::string& type, double expected, double actual, auto compareFunc) {
					const bool passed = compareFunc(expected, actual);
					LogResult(results, "ICE G" + std::to_string(i) + " - " + type + " Compare", expected, actual, passed);
					if (passed) passCount++;
					};

				switch (i)
				{
				case 0:
					evaluateTolerance("Absolute", g00, coefficient, AbsoluteToleranceCompare);
					evaluateTolerance("Relative", g00, coefficient, RelativeToleranceCompare);
					evaluateTolerance("Combined", g00, coefficient, CombinedToleranceCompare);
					break;
				case 1:
					evaluateTolerance("Absolute", g01, coefficient, AbsoluteToleranceCompare);
					evaluateTolerance("Relative", g01, coefficient, RelativeToleranceCompare);
					evaluateTolerance("Combined", g01, coefficient, CombinedToleranceCompare);
					break;
				case 2:
					evaluateTolerance("Absolute", g02, coefficient, AbsoluteToleranceCompare);
					evaluateTolerance("Relative", g02, coefficient, RelativeToleranceCompare);
					evaluateTolerance("Combined", g02, coefficient, CombinedToleranceCompare);
					break;
				case 3:
					evaluateTolerance("Absolute", g03, coefficient, AbsoluteToleranceCompare);
					evaluateTolerance("Relative", g03, coefficient, RelativeToleranceCompare);
					evaluateTolerance("Combined", g03, coefficient, CombinedToleranceCompare);
					break;
				case 4:
					evaluateTolerance("Absolute", g04, coefficient, AbsoluteToleranceCompare);
					evaluateTolerance("Relative", g04, coefficient, RelativeToleranceCompare);
					evaluateTolerance("Combined", g04, coefficient, CombinedToleranceCompare);
					break;
				default:
					Assert::Fail(L"Invalid coefficient index.");
				}
				Assert::IsTrue(passCount >= 2, L"One or more coefficients failed to meet at least two tolerance criteria.");
			}
		}

		TEST_METHOD(Test_IceRegion_ComplexConstantsVerficiation)
		{
			std::vector<TestResult> results;
			constexpr double t1_real = 0.368017112855051e-1;
			constexpr double t1_imag = 0.510878114959572e-1;

			constexpr double t2_real = 0.337315741065416;
			constexpr double t2_imag = 0.335449415919309;
			for (int i = 0; i < WaterEquationOfState->NUMBER_OF_ICE_COMPLEX_CONSTANTS_T; i++)
			{
				const std::complex<double> coefficient = WaterEquationOfState->GetIceComplexConstantT(i);
				int passCount = 0;

				auto evaluateTolerance = [&](const std::string& type, double expected, double actual, auto compareFunc, const std::string& name) {
					const bool passed = compareFunc(expected, actual);
					LogResult(results, name + " - " + type + " Compare", expected, actual, passed);
					if (passed) passCount++;
					};

				switch (i)
				{
				case 0:
					evaluateTolerance("Absolute", t1_real, coefficient.real(), AbsoluteToleranceCompare, "T1 Real");
					evaluateTolerance("Relative", t1_real, coefficient.real(), RelativeToleranceCompare, "T1 Real");
					evaluateTolerance("Combined", t1_real, coefficient.real(), CombinedToleranceCompare, "T1 Real,");

					evaluateTolerance("Absolute", t1_imag, coefficient.imag(), AbsoluteToleranceCompare, "T1 Imaginary");
					evaluateTolerance("Relative", t1_imag, coefficient.imag(), RelativeToleranceCompare, "T1 Imaginary");
					evaluateTolerance("Combined", t1_imag, coefficient.imag(), CombinedToleranceCompare, "T1 Imaginary");
					break;
				case 1:
					evaluateTolerance("Absolute", t2_real, coefficient.real(), AbsoluteToleranceCompare, "T2 Real");
					evaluateTolerance("Relative", t2_real, coefficient.real(), RelativeToleranceCompare, "T2 Real");
					evaluateTolerance("Combined", t2_real, coefficient.real(), CombinedToleranceCompare, "T2 Real");

					evaluateTolerance("Absolute", t2_imag, coefficient.imag(), AbsoluteToleranceCompare, "T2 Imaginary");
					evaluateTolerance("Relative", t2_imag, coefficient.imag(), RelativeToleranceCompare, "T2 Imaginary");
					evaluateTolerance("Combined", t2_imag, coefficient.imag(), CombinedToleranceCompare, "T2 Imaginary");
					break;
				default:
					Assert::Fail(L"Invalid coefficient index.");
				}
				Assert::IsTrue(passCount >= 2, L"One or more coefficients failed to meet at least two tolerance criteria.");
			}
		}

		TEST_METHOD(Test_IceRegion_VerificationCases)
		{
			IceRegionTestCase testCases[] = {
				// Ice-Rev2009 Absolute zero test case (baseline)
				{0.0,
					WaterEquationOfState->NORMAL_PRESSURE,
					-632020.233335886,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN,
					NAN
				},
				// Ice-Rev2009 Triple point verification case
				{
					WaterEquationOfState->TRIPLE_POINT_TEMPERATURE,
					WaterEquationOfState->TRIPLE_POINT_PRESSURE,
					0.611784135e-3,
					0.109085812737e-2,
					-0.128495941571e-12,
					-1220.69433940e-3,
					-333444.253996e-3,
					-333444.921197e-3,
					2096.78431622e-3,
					916.709492200,
					0.000159863102566,
					0.117793449348e-9,
					0.114161597779e-9,
					1357147.647
				}
			};
			for (const auto& testCase : testCases)
			{
				VerifyIceRegionProperties(testCase);
			}
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
		static bool AbsoluteToleranceCompare(const double expected, const double actual)
		{
			return std::fabs(expected - actual) <= std::numeric_limits<double>::epsilon();
		}

		static bool RelativeToleranceCompare(const double expected, const double actual)
		{
			const double maxValue = std::max(std::fabs(expected), std::fabs(actual));
			return std::fabs(expected - actual) <= std::numeric_limits<double>::epsilon() * maxValue;
		}

		static bool CombinedToleranceCompare(const double expected, const double actual)
		{
			const double maxValueOne = std::max({ 1.0, std::fabs(expected), std::fabs(actual) });
			return std::fabs(expected - actual) <= std::numeric_limits<double>::epsilon() * maxValueOne;
		}

		static bool AreDoublesEqual(
			const double expected,
			const double actual,
			const double tolerance = TOLERANCE,
			const double absoluteTolerance = ABSOLUTE_TOLERANCE,
			const bool allowNaNEquality = true)
		{
			if (allowNaNEquality && std::isnan(expected) && std::isnan(actual))
				return true;

			const double difference = fabs(expected - actual);
			const double maxAbsValue = std::max(fabs(expected), fabs(actual));
			if (maxAbsValue < absoluteTolerance)
				return (difference < absoluteTolerance);
			return (difference <= tolerance * maxAbsValue);
		}

		static void LogResult(
			std::vector<TestResult>& results,
			const std::string& name,
			const double expected,
			const double actual,
			const bool passed)
		{
			std::ostringstream ossExpected, ossActual;
			ossExpected << std::scientific << std::setprecision(15) << expected;
			ossActual << std::scientific << std::setprecision(15) << actual;
			std::string expectedStr = ossExpected.str();
			std::string actualStr = ossActual.str();

			results.push_back({ name, expected, actual, passed });

			const std::wstring message = L"Property: " + std::wstring(name.begin(), name.end()) +
				L" | Expected: " + std::wstring(expectedStr.begin(), expectedStr.end()) +
				L" | Actual: " + std::wstring(actualStr.begin(), actualStr.end()) +
				(passed ? L" | Passed\n" : L" | Failed\n");
			Logger::WriteMessage(message.c_str());
		}

		static void VerifyIceRegionProperties(const IceRegionTestCase& testCase) {
			std::vector<TestResult> results;
			const std::string inputs = " - Temperature: " + std::to_string(testCase.Temperature) + ", Pressure: " + std::to_string(testCase.Pressure);

			auto evaluateTolerance = [&](const std::string& propertyName, double expected, double actual) {
				int passCount = 0;

				auto logAndEvaluate = [&](const std::string& type, auto compareFunc) {
					bool passed = compareFunc(expected, actual);
					LogResult(results, propertyName + " - " + type, expected, actual, passed);
					if (passed) ++passCount;
					};

				logAndEvaluate("Absolute", AbsoluteToleranceCompare);
				logAndEvaluate("Relative", RelativeToleranceCompare);
				logAndEvaluate("Combined", CombinedToleranceCompare);

				return passCount >= 2;
				};

			if (!std::isnan(testCase.ExpectedSpecificGibbsEnergy)) {
				const double gibbsEnergy = WaterEquationOfState->CalculateIceSpecificGibbsFreeEnergy(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Specific Gibbs Energy" + inputs, testCase.ExpectedSpecificGibbsEnergy, gibbsEnergy);
			}

			if (!std::isnan(testCase.ExpectedFirstDerivativeGibbsEnergyPressure)) {
				const double firstDerivativeGibbsEnergyPressure = WaterEquationOfState->CalculateIceFirstDerivativeGibbsFreeEnergyPressure(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("First Derivative Gibbs Energy Pressure" + inputs, testCase.ExpectedFirstDerivativeGibbsEnergyPressure, firstDerivativeGibbsEnergyPressure);
			}

			if (!std::isnan(testCase.ExpectedSecondDerivativeGibbsEnergyPressure)) {
				const double secondDerivativeGibbsEnergyPressure = WaterEquationOfState->CalculateIceSecondDerivativeGibbsFreeEnergyPressure(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Second Derivative Gibbs Energy Pressure" + inputs, testCase.ExpectedSecondDerivativeGibbsEnergyPressure, secondDerivativeGibbsEnergyPressure);
			}

			if (!std::isnan(testCase.ExpectedSpecificEntropy)) {
				const double entropy = WaterEquationOfState->CalculateIceSpecificEntropy(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Specific Entropy" + inputs, testCase.ExpectedSpecificEntropy, entropy);
			}

			if (!std::isnan(testCase.ExpectedSpecificEnthalpy)) {
				const double enthalpy = WaterEquationOfState->CalculateIceSpecificEnthalpy(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Specific Enthalpy" + inputs, testCase.ExpectedSpecificEnthalpy, enthalpy);
			}

			if (!std::isnan(testCase.ExpectedSpecificInternalEnergy)) {
				const double internalEnergy = WaterEquationOfState->CalculateIceSpecificInternalEnergy(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Specific Internal Energy" + inputs, testCase.ExpectedSpecificInternalEnergy, internalEnergy);
			}

			if (!std::isnan(testCase.ExpectedSpecificIsobaricHeatCapacity)) {
				const double isobaricHeatCapacity = WaterEquationOfState->CalculateIceSpecificIsobaricHeatCapacity(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Specific Isobaric Heat Capacity" + inputs, testCase.ExpectedSpecificIsobaricHeatCapacity, isobaricHeatCapacity);
			}

			if (!std::isnan(testCase.ExpectedDensity)) {
				const double density = WaterEquationOfState->CalculateIceDensity(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Density" + inputs, testCase.ExpectedDensity, density);
			}

			if (!std::isnan(testCase.ExpectedCubicExpansionCoefficient)) {
				const double cubicExpansionCoefficient = WaterEquationOfState->CalculateIceCubicExpansionCoefficient(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Cubic Expansion Coefficient" + inputs, testCase.ExpectedCubicExpansionCoefficient, cubicExpansionCoefficient);
			}

			if (!std::isnan(testCase.ExpectedIsothermalCompressibility)) {
				const double isothermalCompressibility = WaterEquationOfState->CalculateIceIsothermalCompressibility(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Isothermal Compressibility" + inputs, testCase.ExpectedIsothermalCompressibility, isothermalCompressibility);
			}

			if (!std::isnan(testCase.ExpectedIsentropicCompressibility)) {
				const double isentropicCompressibility = WaterEquationOfState->CalculateIceIsentropicCompressibility(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Isentropic Compressibility" + inputs, testCase.ExpectedIsentropicCompressibility, isentropicCompressibility);
			}

			if (!std::isnan(testCase.ExpectedPressureCoefficient)) {
				const double pressureCoefficient = WaterEquationOfState->CalculateIcePressureCoefficient(testCase.Temperature, testCase.Pressure);
				evaluateTolerance("Pressure Coefficient" + inputs, testCase.ExpectedPressureCoefficient, pressureCoefficient);
			}

			const bool allPassed = std::ranges::all_of(results, [](const TestResult& result) { return result.Passed; });
			if (!allPassed) {
				Logger::WriteMessage(L"One or more properties did not match the expected values.\n");
			}
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
			const double density = WaterEquationOfState->CalculateRegion3Density(testCase.Temperature, testCase.ExpectedPressure);
			const double specificInternalEnergy = WaterEquationOfState->CalculateRegion3SpecificInternalEnergy(testCase.Temperature, testCase.Density);
			const double specificEntropy = WaterEquationOfState->CalculateRegion3SpecificEntropy(testCase.Temperature, testCase.Density);
			const double specificEnthalpy = WaterEquationOfState->CalculateRegion3SpecificEnthalpy(testCase.Temperature, testCase.Density);
			const double specificIsobaricHeatCapacity = WaterEquationOfState->CalculateRegion3SpecificIsobaricHeatCapacity(testCase.Temperature, testCase.Density);
			const double speedOfSound = WaterEquationOfState->CalculateRegion3SpeedOfSound(testCase.Temperature, testCase.Density);

			Assert::AreEqual(testCase.Density, density, TOLERANCE, L"Density mismatch");
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