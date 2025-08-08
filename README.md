# aTBCM Code

## Description
This repository contains MATLAB code for the simulation of the Asymmetrical Triangular Body-Cover Model (aTBCM).
The aTBCM extends the TBCM framework by allowing asymmetries in muscle activation, viscoelastic properties, and geometry between the left and right vocal folds. This enables the study of vocal fold dynamics under pathological or experimentally induced imbalances.

The model is organized into folders according to its components, including asymmetry handling, subglottal tract models, muscle activation, and vocal tract.

## Project Structure
- `@BodyCoverModel`: Structure based on the BCM by Titze and Story
- `@MuscleActivation`: Muscle activation structure for the BCM (Titze's Rules)
- `@MuscleControlModel`: Structure for modeling muscle activation for the TBCM
- `@SubglottalTractModel`: Subglottal tract object introducing the WRA
- `@TriangularBodyCoverModel`: TBCM object (three masses)
- `@VocalTractModel`: Vocal tract object introducing the WRA
- `@TwoFolds`: Object for two Vfs interactions
- `+IntrinsicMuscles`: Constants for the Kelvin model of each laryngeal muscle
- `+measure`: Set of functions to calculate aerodynamic features
- `vfsolver`: Previous version implemented as individual scripts, uses `solveFlow`

## Example Codes

### `Simu_aTBCM.m`

This function-based script takes the following input values:
- Muscle activation
- Subglottal pressure
- Vowel
- Gender

Outputs include 50 ms of signals of interest:
- Geometry
- Mass positions
- Glottal area
- Subglottal pressure
- Collision pressure
- Output pressure (MIC)
- Glottal flow

### `Par_Asimetrico_CompTA.m`

This script uses the `Simu_aTBCM.m` function within for-loops to compute signals and acoustic parameters for multiple phonation configurations.

## References

- Parra, J. A., Calvache, C., Alzamendi, G. A., Ibarra, E. J., Soláque, L., Peterson, S. D., & Zañartu, M. (2024). Asymmetric triangular body-cover model of the vocal folds with bilateral intrinsic muscle activation. The Journal of the Acoustical Society of America, 156(2), 939-953. DOI: https://doi.org/10.1121/10.0028164

- Alzamendi, G. A., Peterson, S. D., Erath, B. D., Hillman, R. E., & Zañartu, M. (2022). Triangular body-cover model of the vocal folds with coordinated activation of the five intrinsic laryngeal muscles. The Journal of the Acoustical Society of America, 151(1), 17-30. DOI: https://doi.org/10.1121/10.0009169
