# **NilsInfinite.EquationsOfState**

An implementation of the IAPWS-IF97 Equation of State for water and steam in C++, covering all five regions with plans for future extensions and improvements.

## **Table of Contents**

- [Introduction](#introduction)
- [Features](#features)
- [Implementation Details](#implementation-details)
- [Usage](#usage)
- [Dependencies](#dependencies)
- [Future Plans](#future-plans)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## **Introduction**

The **NilsInfinite.EquationsOfState** project provides a(nother) C++ implementation of the IAPWS-IF97 equation of state for water and steam.

This initial release includes the base implementation for all five regions defined in the IAPWS-IF97 formulation. Future extensions will include iterative solvers and enhanced validation through improved test cases.

## **Features**

- **Coverage of IF97 Regions:**
  - **Region 1:** Subcooled Water.
  - **Region 2:** Supercritical water/steam.
  - **Region 3:** Superheated steam.
  - **Region 4:** Saturation data.
  - **Region 5:** High-temperature steam.

- **Implemented Properties:**
  - Specific volume
  - Specific internal energy
  - Specific entropy
  - Specific enthalpy
  - Specific isobaric heat capacity
  - Specific isochoric heat capacity (where applicable)
  - Speed of sound
  - Saturation pressure and temperature calculations (Region 4)
  - Region boundary equations (e.g., Region 2/3 boundary)

- **Testing Suite:**
  - Unit tests covering key verification cases based on the IAPWS-IF97 tables.
  - Validation of property calculations across all regions.

## **Implementation Details**

### **Equation of State Formulation**

The implementation follows the official IAPWS-IF97 formulation. The equations and coefficients are based on the standards provided by the International Association for the Properties of Water and Steam (IAPWS).

For detailed information about the IAPWS-IF97 formulation, please refer to the [IAPWS-IF97 documentation](https://asmedigitalcollection.asme.org/gasturbinespower/article-abstract/122/1/150/461340/The-IAPWS-Industrial-Formulation-1997-for-the?redirectedFrom=fulltext).

### **Code Structure**

- **Regions:** Each region is implemented with its own set of equations and methods.
- **Coefficients:** The coefficients for the equations are stored in a SQLite database. The coefficients are loaded into memory during initialization. The database can be updated or modified as needed.

## **Usage**
- Requires the setup of a SQLite database containing the coefficients and correct table names. This will be adjusted in the future. Check my wrapper project for examples of setup in a C# project with P/Invoke.

## **Dependencies**
- **SQLite3:** The implementation uses the SQLite3 library for storing and retrieving the coefficients for the equations. The library is included in the project and does not require external installation. However, feel free to visit the [SQLite website](https://www.sqlite.org/index.html) for more information.

## **Future Plans**
- **Iterative Solvers:** Implement iterative solvers for regions with complex equations.
- **Enhanced Validation:** Expand the testing suite to cover more cases and edge conditions.
- **Backward Equations:** Implement backward equations for properties.
- **Issues:** Address known low tolerance issue with saturation temperature calculation.

## **License**
- **MIT License**

## **Acknowledgments**
- **IAPWS:** For providing the official formulation and standards for water and steam properties. [IAPWS Website](https://iapws.org/)]
